use std::collections::HashMap;
use std::collections::HashSet;
use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result};
use arrow::array::{ArrayRef, Int32Array, StringArray, UInt8Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use duckdb::Connection;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use tracing::info;

use crate::cli::FusionsArgs;
use crate::kmer::{kmer_iter, non_acgt_positions, render_layout_track, rescue_layout_track};

// ─── Index loading ────────────────────────────────────────────────────────────

struct FusionIndex {
    kmer_to_gene: HashMap<u64, u32>,
    gene_names: Vec<String>,
    gene_chroms: Vec<String>,
    /// Per-gene transcription strand ('+'/'-'/'.'); all '.' when `has_strand` is false.
    gene_strands: Vec<char>,
    /// Per-gene 0-based half-open transcript body bounds (0 when `has_strand` is false).
    gene_tx_start: Vec<i64>,
    gene_tx_end: Vec<i64>,
    /// True when the index carries strand + tx coords (meta `gene_strand=1`), enabling
    /// biological 5'->3' fusion-partner ordering. Legacy indexes are false → the scanner
    /// falls back to gene-index-order labels.
    has_strand: bool,
    /// Genome-unique (genome_copies==1) k-mer hashes, populated only when the unique-anchor
    /// filter is active (`--min-unique-anchor-reads > 0`) and the index carries copy counts.
    /// Used to tell a genuine "this read is from gene X" anchor from a multi-copy repeat hit.
    /// `None` when the filter is disabled, so the default hot path is unchanged.
    unique_kmers: Option<HashSet<u64>>,
    /// Per-gene count of genome-unique k-mers (parallel to `gene_names`); used to pick the
    /// higher-uniqueness partner of a pair as the side that must be unique-anchored. Populated
    /// alongside `unique_kmers`.
    gene_unique_total: Option<Vec<u32>>,
}

/// Read a value from the index `meta` table by key.
///
/// Returns `None` if the `meta` table does not exist (older indexes) or the key
/// is absent, so callers can gracefully skip checks against legacy indexes.
fn index_meta_value(conn: &Connection, key: &str) -> Result<Option<String>> {
    let has_meta: bool = conn
        .query_row(
            "SELECT count(*) FROM sqlite_master WHERE type='table' AND name='meta'",
            [],
            |r| r.get::<_, i64>(0),
        )
        .or_else(|_| {
            conn.query_row(
                "SELECT count(*) FROM duckdb_tables() WHERE table_name='meta'",
                [],
                |r| r.get::<_, i64>(0),
            )
        })
        .map(|n| n > 0)
        .unwrap_or(false);
    if !has_meta {
        return Ok(None);
    }
    Ok(conn
        .query_row(
            "SELECT value FROM meta WHERE key=?",
            duckdb::params![key],
            |r| r.get(0),
        )
        .ok())
}

fn index_kmer_size(conn: &Connection) -> Result<Option<usize>> {
    match index_meta_value(conn, "kmer_size")? {
        Some(s) => Ok(Some(s.parse::<usize>().with_context(|| {
            format!("index meta kmer_size is not an integer: {s:?}")
        })?)),
        None => Ok(None),
    }
}

fn index_geac_version(conn: &Connection) -> Result<Option<String>> {
    index_meta_value(conn, "geac_version")
}

/// Return true if the `kmers` table has a `genome_copies` column.
fn kmers_has_copies_column(conn: &Connection) -> Result<bool> {
    let mut stmt = conn.prepare("PRAGMA table_info('kmers')")?;
    let mut rows = stmt.query([])?;
    while let Some(row) = rows.next()? {
        // PRAGMA table_info columns: cid, name, type, notnull, dflt_value, pk
        let name: String = row.get(1)?;
        if name == "genome_copies" {
            return Ok(true);
        }
    }
    Ok(false)
}

fn load_index(
    index_path: &Path,
    k: usize,
    max_kmer_copies: Option<u32>,
    track_unique: bool,
    skip_version_check: bool,
) -> Result<FusionIndex> {
    let conn = Connection::open(index_path)
        .with_context(|| format!("failed to open fusion index: {}", index_path.display()))?;

    if let Some(index_k) = index_kmer_size(&conn)? {
        anyhow::ensure!(
            index_k == k,
            "--kmer-size {k} does not match the index (built with k={index_k}); \
             pass --kmer-size {index_k} or rebuild the index with k={k}"
        );
    }

    let running = env!("CARGO_PKG_VERSION");
    if let Some(index_v) = index_geac_version(&conn)? {
        if index_v != running {
            let msg = format!(
                "fusion index was built with GEAC v{index_v} but this binary is v{running}; \
                 rebuild the index with the current version or pass --skip-version-check to run anyway"
            );
            if skip_version_check {
                tracing::warn!("{msg}");
            } else {
                anyhow::bail!("{msg}");
            }
        }
    }

    // Indexes built with strand support carry a `gene_strand=1` meta flag and the
    // extra genes-table columns. Legacy indexes lack both — load chrom/name only and
    // leave strand unknown so the scanner falls back to gene-index-order labels.
    let has_strand = index_meta_value(&conn, "gene_strand")?.as_deref() == Some("1");
    let mut gene_names: Vec<String> = Vec::new();
    let mut gene_chroms: Vec<String> = Vec::new();
    let mut gene_strands: Vec<char> = Vec::new();
    let mut gene_tx_start: Vec<i64> = Vec::new();
    let mut gene_tx_end: Vec<i64> = Vec::new();

    if has_strand {
        let mut stmt = conn.prepare(
            "SELECT gene_index, gene_name, chrom, strand, tx_start, tx_end \
             FROM genes ORDER BY gene_index",
        )?;
        let mut rows = stmt.query([])?;
        while let Some(row) = rows.next()? {
            let _idx: u32 = row.get(0)?;
            gene_names.push(row.get(1)?);
            gene_chroms.push(row.get(2)?);
            let s: String = row.get(3)?;
            gene_strands.push(s.chars().next().unwrap_or('.'));
            gene_tx_start.push(row.get::<_, u64>(4)? as i64);
            gene_tx_end.push(row.get::<_, u64>(5)? as i64);
        }
    } else {
        let mut stmt =
            conn.prepare("SELECT gene_index, gene_name, chrom FROM genes ORDER BY gene_index")?;
        let mut rows = stmt.query([])?;
        while let Some(row) = rows.next()? {
            let _idx: u32 = row.get(0)?;
            gene_names.push(row.get(1)?);
            gene_chroms.push(row.get(2)?);
            gene_strands.push('.');
            gene_tx_start.push(0);
            gene_tx_end.push(0);
        }
    }
    info!(
        n_genes = gene_names.len(),
        has_strand, "loaded genes from index"
    );
    if !has_strand {
        // A legacy index (built before gene_strand was recorded) makes
        // `infer_five_prime` abstain on every call, so `partner_order` is `index`
        // cohort-wide. That silently reads downstream as "orientation is broken"
        // rather than "the index predates orientation support" — call it out once.
        tracing::warn!(
            "fusion index has no gene strand (gene_strand meta absent); 5'->3' \
             orientation is disabled and every partner_order will be `index` — \
             rebuild the index with a current geac to enable orientation"
        );
    }

    // The k-mer table dominates memory on large indexes (tens-to-hundreds of
    // millions of entries). Reserve up front from the row count so the map is
    // allocated once instead of repeatedly doubling — each doubling holds the old
    // and new buffers simultaneously, briefly spiking peak RSS to ~1.7x the final
    // size, which is enough to OOM a job that otherwise fits.
    let n_total: i64 = conn.query_row("SELECT count(*) FROM kmers", [], |r| r.get(0))?;

    // Read genome_copies when either a copy filter (--max-kmer-copies) or the
    // unique-anchor filter (--min-unique-anchor-reads) is active. Both need per-k-mer
    // copy counts; without one of them the default path stays identical (and compatible
    // with indexes built before the column existed).
    let need_copies = max_kmer_copies.is_some() || track_unique;
    if need_copies {
        anyhow::ensure!(
            kmers_has_copies_column(&conn)?,
            "this run needs per-k-mer genome copy counts (from --max-kmer-copies or \
             --min-unique-anchor-reads), but the index has no genome_copies column; \
             rebuild it with --check-genome-uniqueness or drop the flag"
        );
        // Fail fast on indexes whose genome_copies are all NULL (column present but
        // the genome pass never ran) before loading the whole kmers table into RAM.
        let n_with_copies: i64 =
            conn.query_row("SELECT count(genome_copies) FROM kmers", [], |r| r.get(0))?;
        anyhow::ensure!(
            n_total == 0 || n_with_copies > 0,
            "this run needs per-k-mer genome copy counts, but the index has none \
             (genome_copies is NULL for every k-mer); rebuild with --check-genome-uniqueness"
        );
    }

    let mut kmer_to_gene: HashMap<u64, u32> = HashMap::new();
    // Genome-unique k-mer set + per-gene unique totals, built only for the unique-anchor
    // filter. `unique_kmers` lets a per-read scan tell a genuine anchor from a repeat hit;
    // `gene_unique_total` picks the higher-uniqueness partner of a pair.
    let mut unique_kmers: Option<HashSet<u64>> = track_unique.then(HashSet::new);
    let mut gene_unique_total: Option<Vec<u32>> =
        track_unique.then(|| vec![0u32; gene_names.len()]);

    if need_copies {
        // Apply the optional copy filter at load time and store only k-mers that pass.
        // Excluded k-mers (copies > max, or NULL/unknown) never enter the map. When the
        // unique-anchor filter is on, also record which kept k-mers are genome-unique.
        let mut stmt = conn.prepare("SELECT kmer_hash, gene_index, genome_copies FROM kmers")?;
        let mut rows = stmt.query([])?;
        let mut n_excluded: u64 = 0;
        while let Some(row) = rows.next()? {
            let kmer_i64: i64 = row.get(0)?;
            let gene_idx: u32 = row.get(1)?;
            let copies: Option<i32> = row.get(2)?;
            // Keep the k-mer unless a max-copies filter excludes it (copies above the
            // ceiling, or copy count unknown when a ceiling is set).
            let keep = match (max_kmer_copies, copies) {
                (Some(max), Some(c)) => c >= 0 && (c as u32) <= max,
                (Some(_), None) => false,
                (None, _) => true,
            };
            if !keep {
                n_excluded += 1;
                continue;
            }
            let kmer = kmer_i64 as u64;
            kmer_to_gene.insert(kmer, gene_idx);
            if track_unique && copies == Some(1) {
                unique_kmers.as_mut().unwrap().insert(kmer);
                if let Some(v) = gene_unique_total.as_mut() {
                    if let Some(slot) = v.get_mut(gene_idx as usize) {
                        *slot += 1;
                    }
                }
            }
        }
        if n_excluded > 0 {
            info!(
                n_excluded,
                ?max_kmer_copies,
                "k-mers excluded by the genome-copy filter (too common or copy count unknown)"
            );
        }
    } else {
        kmer_to_gene.reserve(n_total.max(0) as usize);
        let mut stmt = conn.prepare("SELECT kmer_hash, gene_index FROM kmers")?;
        let mut rows = stmt.query([])?;
        while let Some(row) = rows.next()? {
            let kmer_i64: i64 = row.get(0)?;
            let gene_idx: u32 = row.get(1)?;
            kmer_to_gene.insert(kmer_i64 as u64, gene_idx);
        }
    }
    if let Some(uk) = unique_kmers.as_ref() {
        info!(
            n_kmers = kmer_to_gene.len(),
            n_unique_kmers = uk.len(),
            "loaded k-mer index (unique-anchor tracking on)"
        );
    } else {
        info!(n_kmers = kmer_to_gene.len(), "loaded k-mer index");
    }

    Ok(FusionIndex {
        kmer_to_gene,
        gene_names,
        gene_chroms,
        gene_strands,
        gene_tx_start,
        gene_tx_end,
        has_strand,
        unique_kmers,
        gene_unique_total,
    })
}

// ─── Read → gene assignment ───────────────────────────────────────────────────

struct ReadHit {
    gene_idx: u32,
    gene1_kmer_count: u32,
    gene1_kmer_min: u32,
    gene1_kmer_max: u32,
    /// Count of this read's gene1 k-mers that are genome-unique (genome_copies==1).
    /// 0 unless the unique-anchor filter is active (`index.unique_kmers.is_some()`).
    gene1_unique_count: u32,
    // Second-best gene on this read (if any k-mer hit a different gene).
    // Used by the junction-coherence filter to detect spanning reads.
    gene2_idx: Option<u32>,
    gene2_kmer_count: u32,
    gene2_kmer_min: u32,
    gene2_kmer_max: u32,
    /// Genome-unique k-mer count for gene2 (see `gene1_unique_count`).
    gene2_unique_count: u32,
    mapq: u8,
    // kmer_hits is populated only during the second BAM pass for kmer-hits-output;
    // empty during the first pass to avoid gigabytes of per-read detail on deep BAMs.
    kmer_hits: Vec<(u64, u32, usize)>, // (kmer_hash, matched_gene_idx, pos_in_read)
}

fn assign_gene(
    seq: &[u8],
    index: &FusionIndex,
    k: usize,
    min_hits: u32,
    collect_detail: bool,
    kmer_blacklist: Option<&std::collections::HashSet<u64>>,
) -> Option<ReadHit> {
    // gene_idx → (total_count, clean_count, unique_count, min_pos, max_pos)
    // total_count: all matching k-mers; clean_count: non-blacklisted k-mers only;
    // unique_count: matching k-mers that are genome-unique (only tracked when the
    // unique-anchor filter is active). Gene winner selection uses total_count for
    // stability; the min_hits gate uses clean_count so blacklisted k-mers alone
    // cannot qualify a read.
    let mut gene_hits: HashMap<u32, (u32, u32, u32, usize, usize)> = HashMap::new();
    let mut kmer_hits: Vec<(u64, u32, usize)> = Vec::new();
    let unique_kmers = index.unique_kmers.as_ref();
    // The genome-copy filter (--max-kmer-copies) is applied at index load time:
    // excluded k-mers are absent from kmer_to_gene, so a plain lookup is correct.
    for (kmer, pos_in_read) in kmer_iter(seq, k) {
        if let Some(&gene_idx) = index.kmer_to_gene.get(&kmer) {
            let e = gene_hits
                .entry(gene_idx)
                .or_insert((0, 0, 0, pos_in_read, pos_in_read));
            e.0 += 1;
            let is_clean = kmer_blacklist.is_none_or(|bl| !bl.contains(&kmer));
            if is_clean {
                e.1 += 1;
            }
            // Unique lookup only happens when the filter is on (unique_kmers is Some).
            if unique_kmers.is_some_and(|u| u.contains(&kmer)) {
                e.2 += 1;
            }
            e.3 = e.3.min(pos_in_read);
            e.4 = e.4.max(pos_in_read);
            // Retaining every matching k-mer for every assigned read costs
            // gigabytes on deep BAMs; only do it when the detail TSV is requested.
            if collect_detail {
                kmer_hits.push((kmer, gene_idx, pos_in_read));
            }
        }
    }

    // Sort by total_count desc, then gene_idx asc so the result is deterministic
    // (HashMap iteration order must not influence the winner choice).
    let mut sorted: Vec<(u32, u32, u32, u32, usize, usize)> = gene_hits
        .into_iter()
        .map(|(gi, (total, clean, uniq, mn, mx))| (gi, total, clean, uniq, mn, mx))
        .collect();
    sorted.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

    let (winning_gene, count, clean_count, g1_unique, g1_min, g1_max) = *sorted.first()?;
    // Gate on clean (non-blacklisted) k-mer count: blacklisted k-mers alone cannot
    // qualify a read, but they still contribute once clean evidence is sufficient.
    if clean_count < min_hits {
        return None;
    }

    // Expose the second gene whenever any k-mer matched it. Quality gating is
    // handled separately: coherence requires >= min_anchor_kmers from each side
    // (enforced inside read_coherence), and fusion pair selection via
    // fragment_top_pair still relies on gene_idx (the primary/winning gene) which
    // is gated on min_hits. Setting gene2_idx at the 1-k-mer level means
    // n_spanning_reads correctly reflects all reads that cross the junction,
    // even when the anchor on one side is short (asymmetric junction coverage).
    let (gene2_idx, gene2_kmer_count, g2_unique, g2_min, g2_max) = if sorted.len() >= 2 {
        let (g2, c2, _clean2, u2, mn2, mx2) = sorted[1];
        (Some(g2), c2, u2, mn2, mx2)
    } else {
        (None, 0, 0, 0, 0)
    };

    Some(ReadHit {
        gene_idx: winning_gene,
        gene1_kmer_count: count,
        gene1_kmer_min: g1_min as u32,
        gene1_kmer_max: g1_max as u32,
        gene1_unique_count: g1_unique,
        gene2_idx,
        gene2_kmer_count,
        gene2_kmer_min: g2_min as u32,
        gene2_kmer_max: g2_max as u32,
        gene2_unique_count: g2_unique,
        mapq: 0,
        kmer_hits,
    })
}

/// Select the fusion gene-pair for a fragment's read hits.
///
/// Primary votes: each read casts 1 vote for its `gene_idx`. If fewer than two
/// distinct genes are identified by primary votes alone, a supplemental round adds
/// `gene2_kmer_count` votes for each read's `gene2_idx` — weighting by k-mer count
/// so a genuine junction read (many k-mers on both sides) outweighs a noise hit
/// (1 k-mer). This handles single-read chimeric fragments (mate unmapped/filtered)
/// and multi-read fragments where the junction falls late in every read so all
/// primary assignments land on one gene.
///
/// The top two genes are chosen by vote count descending, breaking ties by
/// **ascending gene index** so the result is deterministic. The pair is normalized
/// to `(min_index, max_index)` so it is a stable key for the whole run.
fn fragment_top_pair(read_hits: &[ReadHit]) -> Option<(u32, u32)> {
    let mut gene_vote: HashMap<u32, usize> = HashMap::new();
    for rh in read_hits {
        *gene_vote.entry(rh.gene_idx).or_insert(0) += 1;
    }
    if gene_vote.len() >= 2 {
        // Primary votes identified two or more distinct genes — pick top two
        // by vote count, breaking ties by ascending gene index.
        let mut sorted: Vec<(u32, usize)> = gene_vote.into_iter().collect();
        sorted.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
        let ga = sorted[0].0;
        let gb = sorted[1].0;
        return Some(if ga <= gb { (ga, gb) } else { (gb, ga) });
    }
    // Primary votes identified only one gene. The primary gene must always be
    // part of the returned pair; supplemental evidence only selects its partner.
    // Adding supplemental k-mer-weighted votes into the same map would allow
    // two competing gene2 candidates to both outscore the primary gene and
    // exclude it from the pair — so we handle the two steps separately.
    let primary_gene = gene_vote.into_keys().next()?;
    let mut supp: HashMap<u32, usize> = HashMap::new();
    for rh in read_hits {
        if let Some(g2) = rh.gene2_idx {
            *supp.entry(g2).or_insert(0) += rh.gene2_kmer_count as usize;
        }
    }
    // Pick the partner with the highest total k-mer evidence; break ties by
    // ascending gene index so the choice is deterministic.
    let (partner, _) = supp
        .into_iter()
        .max_by(|a, b| a.1.cmp(&b.1).then_with(|| b.0.cmp(&a.0)))?;
    let (ga, gb) = (primary_gene, partner);
    Some(if ga <= gb { (ga, gb) } else { (gb, ga) })
}

/// Check whether a single read shows a coherent A-block→B-block k-mer partition for
/// the given fusion pair `(ga, gb)`.
///
/// Returns `(is_spanning, is_coherent)`:
/// - `is_spanning`: the read has k-mers assigned to both ga and gb.
/// - `is_coherent`: the read is spanning AND both gene blocks are non-overlapping
///   in `kmer_pos_in_read` AND each block has at least `min_anchor` k-mers.
///
/// A coherent read indicates a real junction splice (disjoint A-then-B or B-then-A
/// blocks). An interleaved read indicates homology / paralog artifact — the same
/// read bases match both genes, not a chimeric junction.
fn read_coherence(rh: &ReadHit, ga: u32, gb: u32, k: usize, min_anchor: u32) -> (bool, bool) {
    // Determine which side of the ReadHit corresponds to ga and which to gb.
    let (a_count, a_min, a_max, b_count, b_min, b_max) = if rh.gene_idx == ga {
        match rh.gene2_idx {
            Some(g2) if g2 == gb => (
                rh.gene1_kmer_count,
                rh.gene1_kmer_min,
                rh.gene1_kmer_max,
                rh.gene2_kmer_count,
                rh.gene2_kmer_min,
                rh.gene2_kmer_max,
            ),
            _ => return (false, false),
        }
    } else if rh.gene_idx == gb {
        match rh.gene2_idx {
            Some(g2) if g2 == ga => (
                rh.gene2_kmer_count,
                rh.gene2_kmer_min,
                rh.gene2_kmer_max,
                rh.gene1_kmer_count,
                rh.gene1_kmer_min,
                rh.gene1_kmer_max,
            ),
            _ => return (false, false),
        }
    } else {
        return (false, false);
    };

    // Spanning but without sufficient anchor on both sides → not coherent.
    if a_count < min_anchor || b_count < min_anchor {
        return (true, false);
    }

    // Two k-mer blocks [a_min, a_max+k) and [b_min, b_max+k) are non-overlapping iff
    // one ends before the other begins.
    let k32 = k as u32;
    let coherent = (a_max + k32 <= b_min) || (b_max + k32 <= a_min);
    (true, coherent)
}

/// Build the per-k-mer-window layout string for a read against fusion pair (ga, gb).
/// One character per k-mer start position (reference 5'→3', the BAM `SEQ` orientation):
/// `A`/`B` where the window's k-mer maps to gene ga/gb in the index, `N` for a window
/// containing a non-ACGT base (no k-mer is emitted there), and `.` for a k-mer matching
/// neither gene. The exact-match track from [`render_layout_track`] is then passed through
/// [`rescue_layout_track`], which upgrades `.`/`N` windows that are within edit distance 1
/// of a gene k-mer: lowercase `a`/`b` for SNP-rescued windows (bridging the `.` gap a single
/// SNP opens in an otherwise solid block), and CAPITAL `A`/`B` for single-`N` windows that
/// resolve uniquely. Exact `A`/`B` are never altered. Written to the `FL` BAM tag so each
/// evidence read carries its own A→B layout for IGV / samtools inspection. `diagnose-fusion`'s
/// `layout_5to3` applies the same render + rescue, so the two are byte-for-byte identical.
fn fusion_layout_track(
    seq: &[u8],
    k: usize,
    kmer_to_gene: &HashMap<u64, u32>,
    ga: u32,
    gb: u32,
) -> String {
    if seq.len() < k {
        return String::new();
    }
    let n_windows = seq.len() - k + 1;
    let mut a_positions = Vec::new();
    let mut b_positions = Vec::new();
    for (kmer, pos) in kmer_iter(seq, k) {
        if let Some(&g) = kmer_to_gene.get(&kmer) {
            if g == ga {
                a_positions.push(pos);
            } else if g == gb {
                b_positions.push(pos);
            }
        }
    }
    let mut track = render_layout_track(
        n_windows,
        &a_positions,
        &b_positions,
        &non_acgt_positions(seq),
        k,
    )
    .into_bytes();
    rescue_layout_track(&mut track, seq, k, kmer_to_gene, ga, gb);
    String::from_utf8(track).unwrap()
}

/// Per-fusion candidate statistics accumulated during fragment aggregation.
struct FusionStats {
    count: u32,
    min_mapq: u8,
    n_spanning_reads: u32,
    /// Number of fragments (qnames) where at least one read showed a coherent
    /// A-block→B-block k-mer partition. Incremented once per fragment, not per read.
    n_coherent_fragments: u32,
    /// Number of fragments with a genome-unique k-mer anchor on the higher-uniqueness
    /// partner (the unique-anchor specificity filter). 0 unless that filter is active.
    n_unique_anchored: u32,
}

/// Whether a read is anchored to `anchor_gene` by at least one genome-unique k-mer,
/// i.e. the read genuinely contains that gene's unique sequence rather than merely a
/// multi-copy repeat k-mer that also appears there. Drives the unique-anchor filter.
fn read_unique_anchored(rh: &ReadHit, anchor_gene: u32) -> bool {
    if rh.gene_idx == anchor_gene {
        rh.gene1_unique_count >= 1
    } else if rh.gene2_idx == Some(anchor_gene) {
        rh.gene2_unique_count >= 1
    } else {
        false
    }
}

/// Canonical `GENE_A::GENE_B` label for a normalized `(min_index, max_index)` pair.
/// Uses gene-index order so every output — Parquet, TSV, the `FX` BAM tag, and the
/// kmer-hits TSV — labels the same fusion identically. (The Parquet `gene_a`/`gene_b`
/// columns are filled from `gene_names[key.0]`/`[key.1]`, matching this order.)
fn fusion_pair_label(key: (u32, u32), gene_names: &[String]) -> String {
    format!(
        "{}::{}",
        gene_names[key.0 as usize], gene_names[key.1 as usize]
    )
}

// ─── Breakpoint accumulation types ───────────────────────────────────────────

struct SpanningReadData {
    chrom: String,
    pos: i64,
    last_a: usize,
    first_a: usize,
    last_b: usize,
    first_b: usize,
}

/// Positional data for a single-gene read approaching the fusion junction.
///
/// A read whose k-mers matched only one fusion partner — a concordant-pair mate.
/// Used as a breakpoint-estimation fallback when spanning reads are too sparse
/// (e.g. the breakpoint falls inside a repeat that is absent from the k-mer index,
/// so no single read can carry k-mers from both partners across the gap).
///
/// The junction-facing *alignment* edge estimates where this gene's sequence meets
/// the partner. We use the aligner's reference coordinates rather than k-mer offsets
/// because the aligner places the read through the repeat using its full sequence,
/// whereas the last/first matched k-mer stops at the edge of the indexed region and
/// scatters with where k-mers happen to fall in each read.
///   - gene-A reads: the alignment END   (gene-A sequence stops here)
///   - gene-B reads: the alignment START (gene-B sequence begins here)
struct SingleGeneReadData {
    chrom: String,
    aln_start: i64, // 0-based reference start
    aln_end: i64,   // 0-based exclusive reference end
}

struct BreakpointAccumulator {
    gene_a_idx: u32,
    gene_b_idx: u32,
    gene_a_chrom_votes: HashMap<String, u32>,
    gene_b_chrom_votes: HashMap<String, u32>,
    spanning_reads: Vec<SpanningReadData>,
    /// Single-gene gene-A reads: used as concordant-pair fallback when there
    /// are too few spanning reads to satisfy the chimera filter (e.g. when the
    /// breakpoint falls inside a repetitive element that is absent from the
    /// k-mer index).
    gene_a_single_reads: Vec<SingleGeneReadData>,
    /// Single-gene gene-B reads (see gene_a_single_reads).
    gene_b_single_reads: Vec<SingleGeneReadData>,
}

fn median_i64(values: &mut [i64]) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    values.sort_unstable();
    let n = values.len();
    if n.is_multiple_of(2) {
        Some((values[n / 2 - 1] + values[n / 2]) as f64 / 2.0)
    } else {
        Some(values[n / 2] as f64)
    }
}

fn std_dev_i64(values: &[i64]) -> Option<f64> {
    if values.len() < 2 {
        return None;
    }
    let mean = values.iter().sum::<i64>() as f64 / values.len() as f64;
    let variance = values
        .iter()
        .map(|&v| (v as f64 - mean).powi(2))
        .sum::<f64>()
        / (values.len() - 1) as f64;
    Some(variance.sqrt())
}

/// Half-width (bp) of the window used to isolate the dominant breakpoint cluster
/// before measuring consensus spread. A genuine junction concentrates read-edge
/// estimates within a read-length or so; concordant mates deep in the gene body and
/// paralog/chimera artifacts sit kilobases away. Trimming to this window around the
/// median lets a few distant reads be ignored without masking a truly scattered
/// (non-consensus) signal — that still fails the `--max-breakpoint-std` test, which
/// is far tighter than this window.
const BREAKPOINT_CLUSTER_WINDOW: i64 = 1000;

/// Keep only the estimates within `window` bp of the median — the dominant cluster.
/// Robust to a minority of far-off outliers because the median sits inside the
/// cluster whenever it holds the majority of reads.
fn cluster_around_median(estimates: &[i64], window: i64) -> Vec<i64> {
    if estimates.is_empty() {
        return Vec::new();
    }
    let mut sorted = estimates.to_vec();
    sorted.sort_unstable();
    let med = sorted[sorted.len() / 2];
    sorted
        .into_iter()
        .filter(|&e| (e - med).abs() <= window)
        .collect()
}

/// Per-fusion breakpoint consensus, derived from the spanning reads collected in a
/// `BreakpointAccumulator`. Drives both the breakpoints TSV and the
/// `--max-breakpoint-std` filter. `bp_*_std` is `None` when fewer than two reads
/// support that side (a standard deviation is undefined), which the consensus filter
/// treats as a failure.
struct BreakpointStats {
    gene_a_idx: u32,
    gene_b_idx: u32,
    chrom_a: String,
    chrom_b: String,
    bp_a: Option<f64>,
    bp_a_n: usize,
    bp_a_std: Option<f64>,
    bp_b: Option<f64>,
    bp_b_n: usize,
    bp_b_std: Option<f64>,
    n_spanning: usize,
    /// Biological 5'->3' orientation: Some(true) if gene_a (pair_key.0) is the 5'
    /// (donor) partner, Some(false) if gene_b is, None when undetermined (index lacks
    /// strand, unknown strands, ambiguous, or no read majority). Drives oriented
    /// partner ordering in the fusions/breakpoints output; None falls back to index order.
    gene_a_is_5prime: Option<bool>,
}

/// "Strong support" breakpoint tier. When BOTH partners are supported by this many
/// converging reads, the call is high-confidence enough to (a) tolerate a looser
/// breakpoint spread than `--max-breakpoint-std` — real high-depth junctions smear over
/// tens-to->100 bp from splice isoforms and k-mer transition-point noise — and (b) be
/// exempt from the same-locus filter. Single-locus paralog leakage is *spanning-read
/// dominated* at one locus and never accumulates this many INDEPENDENT (concordant)
/// reads on both partners, so this tier separates a real co-located junction (a true
/// translocation observed only through concordant pairs, e.g. a breakpoint buried in a
/// repeat) from leakage without weakening either artifact filter for low-support calls.
const STRONG_BREAKPOINT_READS: usize = 25;
const STRONG_BREAKPOINT_STD: f64 = 250.0;

impl BreakpointStats {
    /// Both sides have ≥ `STRONG_BREAKPOINT_READS` converging reads with spread within
    /// `STRONG_BREAKPOINT_STD`. Used to rescue real high-depth junctions from the
    /// `--max-breakpoint-std` ceiling. Counts ALL supporting reads (spanning included),
    /// because a real interchromosomal fusion is legitimately spanning-read dominated.
    fn strong_support(&self) -> bool {
        self.bp_a_n >= STRONG_BREAKPOINT_READS
            && self.bp_b_n >= STRONG_BREAKPOINT_READS
            && self.bp_a_std.is_some_and(|v| v <= STRONG_BREAKPOINT_STD)
            && self.bp_b_std.is_some_and(|v| v <= STRONG_BREAKPOINT_STD)
    }

    /// Asymmetric anchor tier (single-sample low-input rescue). The DOMINANT (better-covered)
    /// breakpoint must be a tight, well-supported junction (≥ `min_hi` reads, std ≤ `max_std`)
    /// and the partner breakpoint needs only ≥ 1 read. Recovers a real junction seen well on
    /// one side but by 1–2 reads on the other — where the symmetric tight tier (≥ N reads on
    /// BOTH sides) fails at low input. Requiring the *dominant* side to be tight (not just any
    /// side) rejects paralog-leakage artifacts, whose high-count side is smeared (loose std)
    /// while a minority side coincidentally converges. The caller additionally gates this on a
    /// MAPQ floor, since leakage reads are multi-mappers.
    fn asym_anchored(&self, min_hi: usize, max_std: f64) -> bool {
        let (hi_n, hi_std, lo_n) = if self.bp_a_n >= self.bp_b_n {
            (self.bp_a_n, self.bp_a_std, self.bp_b_n)
        } else {
            (self.bp_b_n, self.bp_b_std, self.bp_a_n)
        };
        hi_n >= min_hi && lo_n >= 1 && hi_std.is_some_and(|v| v <= max_std)
    }

    /// Like `strong_support`, but counts only INDEPENDENT (non-spanning) reads on each
    /// side. Single-locus paralog leakage is spanning-read dominated at one locus —
    /// the same molecule carries both partners' k-mers — so subtracting spanning reads
    /// collapses it, while a real junction observed through concordant pairs (e.g. a
    /// breakpoint buried in a repeat, where no read spans the gap) retains strong
    /// support. Used to exempt such real co-located junctions from the same-locus filter
    /// without weakening it for leakage.
    fn strong_concordant_support(&self) -> bool {
        let ca = self.bp_a_n.saturating_sub(self.n_spanning);
        let cb = self.bp_b_n.saturating_sub(self.n_spanning);
        ca >= STRONG_BREAKPOINT_READS
            && cb >= STRONG_BREAKPOINT_READS
            && self.bp_a_std.is_some_and(|v| v <= STRONG_BREAKPOINT_STD)
            && self.bp_b_std.is_some_and(|v| v <= STRONG_BREAKPOINT_STD)
    }
}

/// Infer which partner is the biological 5' (donor) gene from spanning-read geometry.
///
/// Returns `Some(true)` when gene_a (`acc.gene_a_idx`, i.e. `pair_key.0`) is the 5'
/// partner, `Some(false)` when gene_b is, and `None` when orientation cannot be
/// determined: the index carries no strand (`has_strand == false`), the relevant
/// strand is unknown (`'.'`), no spanning read is informative, or the vote ties.
///
/// For each spanning read we know `a_before_b` — gene_a's k-mer block precedes gene_b's
/// in the read `SEQ`, which is stored reference-forward, so "before" means lower genomic
/// coordinate at the aligned locus — and which gene the read aligns to (`chrom`/`pos`).
/// The junction is the aligned gene's block edge; with that gene's transcription strand:
///
///   aligned gene is the 5' partner  iff  (its block is first)  ==  (its strand is '+')
///
/// because a '+'-strand gene whose retained block sits at the lower-coordinate side keeps
/// its promoter-proximal 5' portion (donor), and the logic mirrors for '-' strand and for
/// the gene-b-aligned case. Votes are tallied across informative reads; the majority wins.
fn infer_five_prime(
    acc: &BreakpointAccumulator,
    index: &FusionIndex,
    chrom_a: &str,
    chrom_b: &str,
) -> Option<bool> {
    if !index.has_strand {
        return None;
    }
    let ga = acc.gene_a_idx as usize;
    let gb = acc.gene_b_idx as usize;
    let strand_a = index.gene_strands[ga];
    let strand_b = index.gene_strands[gb];
    let same_chrom = chrom_a == chrom_b;

    // Votes for "gene_a is the 5' partner" vs "gene_b is the 5' partner".
    let mut a5_votes: i32 = 0;
    let mut b5_votes: i32 = 0;

    for sd in &acc.spanning_reads {
        let a_before_b = sd.first_a < sd.first_b;
        // Which partner does this read align to?
        let aligns_a = if same_chrom {
            // Same chromosome: disambiguate by which gene body the alignment start sits in.
            let in_a = sd.pos >= index.gene_tx_start[ga] && sd.pos < index.gene_tx_end[ga];
            let in_b = sd.pos >= index.gene_tx_start[gb] && sd.pos < index.gene_tx_end[gb];
            match (in_a, in_b) {
                (true, false) => true,
                (false, true) => false,
                _ => continue, // ambiguous (overlapping or neither) → skip
            }
        } else if sd.chrom == chrom_a {
            true
        } else if sd.chrom == chrom_b {
            false
        } else {
            continue; // read on neither partner's chromosome → uninformative
        };

        if aligns_a {
            // gene_a is 5' iff (a_before_b == strand_a is '+')
            match strand_a {
                '+' => {
                    if a_before_b {
                        a5_votes += 1
                    } else {
                        b5_votes += 1
                    }
                }
                '-' => {
                    if a_before_b {
                        b5_votes += 1
                    } else {
                        a5_votes += 1
                    }
                }
                _ => continue, // unknown strand → uninformative
            }
        } else {
            // gene_b is 5' iff (!a_before_b == strand_b is '+')
            match strand_b {
                '+' => {
                    if a_before_b {
                        a5_votes += 1
                    } else {
                        b5_votes += 1
                    }
                }
                '-' => {
                    if a_before_b {
                        b5_votes += 1
                    } else {
                        a5_votes += 1
                    }
                }
                _ => continue,
            }
        }
    }

    // Confidence gate. Orientation is decided by majority vote, but two guards prevent a
    // misleading call:
    //   * Minimum support — a handful of scattered reads in a capture-limited background
    //     can vote the wrong way; require enough informative spanning reads that the
    //     majority is meaningful. (Validated: real calls carry >=13 such votes; the only
    //     misfires were 5-6-read chimera-filtered calls, which abstain here.)
    //   * Minimum fraction — guards against a ~50/50 coin flip. A real interchromosomal
    //     fusion legitimately yields a minority of reciprocal-derivative-junction reads
    //     (e.g. ~24% for EWSR1::FLI1), so this floor is intentionally permissive; the
    //     dominant orientation is the functional 5'->3' direction.
    // Below either guard, return None and the caller keeps gene-index order.
    let total = a5_votes + b5_votes;
    if total < ORIENT_MIN_VOTES {
        return None;
    }
    let winner = a5_votes.max(b5_votes);
    if (winner as f64) / (total as f64) < ORIENT_MIN_FRACTION {
        return None;
    }
    Some(a5_votes > b5_votes)
}

/// Minimum informative spanning-read votes before orientation is asserted.
const ORIENT_MIN_VOTES: i32 = 10;
/// Minimum share of votes the winning orientation must hold (else abstain). Permissive
/// because reciprocal-derivative junction reads form a legitimate minority.
const ORIENT_MIN_FRACTION: f64 = 0.6;

/// Estimate both breakpoints and their dispersion from one fusion's spanning reads.
/// Each spanning read votes a position for whichever side(s) its alignment chromosome
/// matches; the median is the breakpoint and the standard deviation is the consensus
/// tightness. A real junction yields a tight cluster (std ~ a few bp); scattered
/// estimates (std in the thousands+) indicate paralog/chimera artifacts.
fn compute_breakpoint_stats(
    acc: &BreakpointAccumulator,
    index: &FusionIndex,
    k: usize,
) -> BreakpointStats {
    // Most-voted chromosome, ties broken by lexicographically smallest name so the
    // choice is deterministic. Falls back to the gene's annotated chromosome when no
    // single-gene reads voted (e.g. every supporting read spans the junction).
    let modal_chrom = |votes: &HashMap<String, u32>, fallback: &str| -> String {
        votes
            .iter()
            .max_by(|a, b| a.1.cmp(b.1).then_with(|| b.0.cmp(a.0)))
            .map(|(c, _)| c.clone())
            .unwrap_or_else(|| fallback.to_string())
    };
    let chrom_a = modal_chrom(
        &acc.gene_a_chrom_votes,
        &index.gene_chroms[acc.gene_a_idx as usize],
    );
    let chrom_b = modal_chrom(
        &acc.gene_b_chrom_votes,
        &index.gene_chroms[acc.gene_b_idx as usize],
    );

    let mut bp_a_estimates: Vec<i64> = Vec::new();
    let mut bp_b_estimates: Vec<i64> = Vec::new();
    for sd in &acc.spanning_reads {
        let a_before_b = sd.first_a < sd.first_b;
        if sd.chrom == chrom_a {
            let est = if a_before_b {
                sd.pos + sd.last_a as i64 + k as i64
            } else {
                sd.pos + sd.first_a as i64
            };
            bp_a_estimates.push(est);
        }
        if sd.chrom == chrom_b {
            let est = if a_before_b {
                sd.pos + sd.first_b as i64
            } else {
                sd.pos + sd.last_b as i64 + k as i64
            };
            bp_b_estimates.push(est);
        }
    }

    // Concordant-pair fallback: when spanning reads are too few to pin a breakpoint
    // (e.g. it falls inside a repeat absent from the k-mer index, so no single read
    // carries k-mers from both partners across the gap), supplement with single-gene
    // mates. The junction-facing alignment edge estimates where each gene's sequence
    // meets the partner: gene-A reads stop at the breakpoint (alignment END), gene-B
    // reads start at it (alignment START). Off-modal-chromosome reads (a gene's native
    // locus, or paralogs) are excluded by the chrom filter; deep-body mates kilobases
    // away are trimmed by the clustering step below.
    for sg in &acc.gene_a_single_reads {
        if sg.chrom == chrom_a {
            bp_a_estimates.push(sg.aln_end);
        }
    }
    for sg in &acc.gene_b_single_reads {
        if sg.chrom == chrom_b {
            bp_b_estimates.push(sg.aln_start);
        }
    }

    // Reduce to the dominant breakpoint cluster before measuring spread, so a few
    // distant reads (deep-body concordant mates, paralog hits) cannot inflate the
    // std of an otherwise tight, well-supported junction. A genuinely scattered
    // (non-consensus) signal has no dominant cluster and still fails downstream.
    let bp_a_estimates = cluster_around_median(&bp_a_estimates, BREAKPOINT_CLUSTER_WINDOW);
    let bp_b_estimates = cluster_around_median(&bp_b_estimates, BREAKPOINT_CLUSTER_WINDOW);

    let n_spanning = acc.spanning_reads.len();
    let bp_a_std = std_dev_i64(&bp_a_estimates);
    let bp_b_std = std_dev_i64(&bp_b_estimates);
    let bp_a_n = bp_a_estimates.len();
    let bp_b_n = bp_b_estimates.len();
    let mut bp_a_sorted = bp_a_estimates;
    let mut bp_b_sorted = bp_b_estimates;
    let bp_a = median_i64(&mut bp_a_sorted);
    let bp_b = median_i64(&mut bp_b_sorted);

    let gene_a_is_5prime = infer_five_prime(acc, index, &chrom_a, &chrom_b);

    BreakpointStats {
        gene_a_idx: acc.gene_a_idx,
        gene_b_idx: acc.gene_b_idx,
        chrom_a,
        chrom_b,
        bp_a,
        bp_a_n,
        bp_a_std,
        bp_b,
        bp_b_n,
        bp_b_std,
        n_spanning,
        gene_a_is_5prime,
    }
}

fn decode_kmer(hash: u64, k: usize) -> String {
    const BASES: [u8; 4] = *b"ACGT";
    (0..k)
        .rev()
        .map(|i| BASES[((hash >> (2 * i)) & 3) as usize] as char)
        .collect()
}

// ─── TSV output ───────────────────────────────────────────────────────────────

fn write_fusion_tsv(records: &[FusionRecord], output: &Path) -> Result<()> {
    use std::io::{BufWriter, Write};
    let file = std::fs::File::create(output)
        .with_context(|| format!("failed to create TSV: {}", output.display()))?;
    let mut w = BufWriter::new(file);
    writeln!(w, "sample_id\tgene_a\tgene_b\tchrom_a\tchrom_b\tsupporting_reads\tmin_mapq\tn_spanning_reads\tn_coherent_fragments\tn_pon_samples\tpon_total_samples\tmax_pon_supporting_reads\tfilter\tpartner_order\tn_unique_anchored")?;
    for r in records {
        let max_pon = r
            .max_pon_supporting_reads
            .map(|v| v.to_string())
            .unwrap_or_else(|| "NA".to_string());
        let n_unique = r
            .n_unique_anchored
            .map(|v| v.to_string())
            .unwrap_or_else(|| "NA".to_string());
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.sample_id,
            r.gene_a,
            r.gene_b,
            r.chrom_a,
            r.chrom_b,
            r.supporting_reads,
            r.min_mapq,
            r.n_spanning_reads,
            r.n_coherent_fragments,
            r.n_pon_samples,
            r.pon_total_samples,
            max_pon,
            r.filter,
            r.partner_order,
            n_unique
        )?;
    }
    Ok(())
}

// ─── Parquet output ───────────────────────────────────────────────────────────

struct FusionRecord {
    // Not written to output — used to rebuild fusion_label after filtering.
    pair_key: (u32, u32),
    sample_id: String,
    gene_a: String,
    gene_b: String,
    chrom_a: String,
    chrom_b: String,
    supporting_reads: i32,
    min_mapq: u8,
    // Junction-coherence counts (from --min-coherent-fragments filter).
    // n_spanning_reads counts reads where gene2_idx matches the fusion partner,
    // including reads where the partner-side anchor is < min_anchor_kmers (those
    // reads are spanning but not coherent). Use n_coherent_fragments as the
    // quality-gated signal; n_spanning_reads is a raw sensitivity indicator.
    n_spanning_reads: i32,
    /// Fragments (qnames) with at least one coherent spanning read.
    n_coherent_fragments: i32,
    // Panel-of-Normals annotation. Without --fusion-pon, n_pon_samples and
    // pon_total_samples are 0 and max_pon_supporting_reads is None.
    n_pon_samples: i32,
    pon_total_samples: i32,
    max_pon_supporting_reads: Option<i32>,
    // VCF-style filter status: "PASS" by default, or "pon" when --max-pon-samples
    // is set and the call appears in more PoN samples than the threshold. Rows are
    // never dropped — the label lets downstream tools include or exclude them.
    filter: String,
    // How gene_a/gene_b are ordered: "5to3" when the biological 5'->3' orientation was
    // inferred from read geometry (gene_a is the 5'/donor partner), or "index" when it
    // could not be determined and partners fall back to gene-index order. Distinguishes a
    // biologically meaningful order from an arbitrary one for downstream/clinical readers.
    partner_order: String,
    // Supporting fragments anchored by a genome-unique (genome_copies==1) k-mer on the
    // higher-uniqueness partner — the raw evidence behind the unique-anchor filter. The
    // denominator is `supporting_reads`, so n_unique_anchored / supporting_reads is the
    // unique-anchored fraction; any threshold can be applied offline without re-running.
    // Some(_) only when uniqueness tracking ran (--min-unique-anchor-reads > 0 or
    // --emit-unique-anchor with a --check-genome-uniqueness index); None ("NA") otherwise,
    // so a strict/untracked run is not misread as "0 unique anchors".
    n_unique_anchored: Option<i32>,
}

fn write_fusion_parquet(records: &[FusionRecord], output: &Path) -> Result<()> {
    let schema = Arc::new(Schema::new(vec![
        Field::new("sample_id", DataType::Utf8, false),
        Field::new("gene_a", DataType::Utf8, false),
        Field::new("gene_b", DataType::Utf8, false),
        Field::new("chrom_a", DataType::Utf8, false),
        Field::new("chrom_b", DataType::Utf8, false),
        Field::new("supporting_reads", DataType::Int32, false),
        Field::new("min_mapq", DataType::UInt8, false),
        Field::new("n_spanning_reads", DataType::Int32, false),
        Field::new("n_coherent_fragments", DataType::Int32, false),
        Field::new("n_pon_samples", DataType::Int32, false),
        Field::new("pon_total_samples", DataType::Int32, false),
        Field::new("max_pon_supporting_reads", DataType::Int32, true),
        Field::new("filter", DataType::Utf8, false),
        Field::new("partner_order", DataType::Utf8, false),
        Field::new("n_unique_anchored", DataType::Int32, true),
    ]));

    let file = std::fs::File::create(output)
        .with_context(|| format!("failed to create output: {}", output.display()))?;
    let props = WriterProperties::builder()
        .set_compression(Compression::SNAPPY)
        .build();
    let mut writer = ArrowWriter::try_new(file, Arc::clone(&schema), Some(props))
        .context("failed to create Parquet writer")?;

    let batch = RecordBatch::try_new(
        Arc::clone(&schema),
        vec![
            Arc::new(StringArray::from(
                records
                    .iter()
                    .map(|r| r.sample_id.as_str())
                    .collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(StringArray::from(
                records
                    .iter()
                    .map(|r| r.gene_a.as_str())
                    .collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(StringArray::from(
                records
                    .iter()
                    .map(|r| r.gene_b.as_str())
                    .collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(StringArray::from(
                records
                    .iter()
                    .map(|r| r.chrom_a.as_str())
                    .collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(StringArray::from(
                records
                    .iter()
                    .map(|r| r.chrom_b.as_str())
                    .collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(Int32Array::from(
                records
                    .iter()
                    .map(|r| r.supporting_reads)
                    .collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(UInt8Array::from(
                records.iter().map(|r| r.min_mapq).collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(Int32Array::from(
                records
                    .iter()
                    .map(|r| r.n_spanning_reads)
                    .collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(Int32Array::from(
                records
                    .iter()
                    .map(|r| r.n_coherent_fragments)
                    .collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(Int32Array::from(
                records.iter().map(|r| r.n_pon_samples).collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(Int32Array::from(
                records
                    .iter()
                    .map(|r| r.pon_total_samples)
                    .collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(Int32Array::from(
                records
                    .iter()
                    .map(|r| r.max_pon_supporting_reads)
                    .collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(StringArray::from(
                records
                    .iter()
                    .map(|r| r.filter.as_str())
                    .collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(StringArray::from(
                records
                    .iter()
                    .map(|r| r.partner_order.as_str())
                    .collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(Int32Array::from(
                records
                    .iter()
                    .map(|r| r.n_unique_anchored)
                    .collect::<Vec<_>>(),
            )) as ArrayRef,
        ],
    )
    .context("failed to build record batch")?;

    writer.write(&batch).context("failed to write batch")?;
    writer.close().context("failed to finalize Parquet")?;
    Ok(())
}

// ─── Panel-of-Normals annotation ──────────────────────────────────────────────

/// Annotate fusion records against a fusion Panel-of-Normals DuckDB (a `geac merge`
/// of normal-sample `*.fusions.parquet` files, which carries a `fusions` table).
///
/// Matching is by the alphabetically-sorted gene-name pair, so a PoN built with one
/// k-mer index still matches calls made with another (gene-index order may differ).
fn annotate_fusion_pon(
    records: &mut [FusionRecord],
    pon_db: &Path,
    skip_version_check: bool,
) -> Result<()> {
    if !pon_db.exists() {
        anyhow::bail!("fusion PoN DuckDB not found: {}", pon_db.display());
    }
    let conn = Connection::open(pon_db)
        .with_context(|| format!("failed to open fusion PoN DuckDB: {}", pon_db.display()))?;

    conn.execute_batch("SELECT 1 FROM fusions LIMIT 0")
        .context(
            "fusion PoN DuckDB does not contain a fusions table — build it with \
         `geac merge` over normal-sample *.fusions.parquet files",
        )?;

    let running = env!("CARGO_PKG_VERSION");
    let pon_version: Option<String> = conn
        .query_row("SELECT geac_version FROM geac_metadata LIMIT 1", [], |r| {
            r.get(0)
        })
        .ok();
    if let Some(pon_v) = pon_version {
        if pon_v != running {
            let msg = format!(
                "fusion PoN was built with GEAC v{pon_v} but this binary is v{running}; \
                 rebuild the PoN with the current version or pass --skip-version-check to run anyway"
            );
            if skip_version_check {
                tracing::warn!("{msg}");
            } else {
                anyhow::bail!("{msg}");
            }
        }
    }

    let pon_total_samples: i32 =
        conn.query_row("SELECT COUNT(DISTINCT sample_id) FROM fusions", [], |r| {
            r.get(0)
        })?;

    // Aggregate per normalized (sorted) gene pair: how many distinct PoN samples carry
    // it and the highest supporting-read count seen.
    let mut stmt = conn.prepare(
        r#"
        SELECT
            LEAST(gene_a, gene_b)    AS g1,
            GREATEST(gene_a, gene_b) AS g2,
            COUNT(DISTINCT sample_id) AS n_samples,
            MAX(supporting_reads)     AS max_reads
        FROM fusions
        GROUP BY g1, g2
        "#,
    )?;
    let mut pon_agg: HashMap<(String, String), (i32, i32)> = HashMap::new();
    let mut rows = stmt.query([])?;
    while let Some(row) = rows.next()? {
        let g1: String = row.get(0)?;
        let g2: String = row.get(1)?;
        let n_samples: i32 = row.get(2)?;
        let max_reads: i32 = row.get(3)?;
        pon_agg.insert((g1, g2), (n_samples, max_reads));
    }

    for r in records.iter_mut() {
        r.pon_total_samples = pon_total_samples;
        let key = if r.gene_a <= r.gene_b {
            (r.gene_a.clone(), r.gene_b.clone())
        } else {
            (r.gene_b.clone(), r.gene_a.clone())
        };
        if let Some(&(n_samples, max_reads)) = pon_agg.get(&key) {
            r.n_pon_samples = n_samples;
            r.max_pon_supporting_reads = Some(max_reads);
        }
    }

    info!(
        pon_total_samples,
        n_pon_pairs = pon_agg.len(),
        "fusion PoN annotation complete"
    );
    Ok(())
}

// ─── Progress reporting ───────────────────────────────────────────────────────

/// Format a number of seconds as a compact human-readable duration (e.g. "1h05m",
/// "7m12s", "43s"). Returns "?" for non-finite or negative inputs.
fn format_duration(secs: f64) -> String {
    if !secs.is_finite() || secs < 0.0 {
        return "?".to_string();
    }
    let total = secs as u64;
    let h = total / 3600;
    let m = (total % 3600) / 60;
    let s = total % 60;
    if h > 0 {
        format!("{h}h{m:02}m")
    } else if m > 0 {
        format!("{m}m{s:02}s")
    } else {
        format!("{s}s")
    }
}

/// Emit a BAM scan progress line. `frac` is the estimated fraction of the file
/// scanned (0.0–1.0): for BAM it comes from the compressed byte offset, for CRAM
/// from the mapped genomic coordinate. When `frac` is `None` we report only
/// throughput.
fn log_scan_progress(
    reads_processed: u64,
    reads_assigned: u64,
    frac: Option<f64>,
    elapsed: std::time::Duration,
) {
    let secs = elapsed.as_secs_f64();
    let reads_per_sec = if secs > 0.0 {
        (reads_processed as f64 / secs) as u64
    } else {
        0
    };

    if let Some(frac) = frac {
        let frac = frac.clamp(0.0, 1.0);
        let eta_secs = if frac > 0.0 {
            secs * (1.0 - frac) / frac
        } else {
            f64::INFINITY
        };
        info!(
            reads_processed,
            reads_assigned,
            percent = format!("{:.1}", frac * 100.0),
            eta = %format_duration(eta_secs),
            reads_per_sec,
            "BAM scan progress"
        );
    } else {
        info!(
            reads_processed,
            reads_assigned, reads_per_sec, "BAM scan progress"
        );
    }
}

// ─── Main entry point ─────────────────────────────────────────────────────────

pub fn detect_fusions(args: &FusionsArgs) -> Result<()> {
    let k = args.kmer_size as usize;
    anyhow::ensure!(k > 0 && k <= 31, "--kmer-size must be between 1 and 31");
    anyhow::ensure!(
        args.max_pon_samples.is_none() || args.fusion_pon.is_some(),
        "--max-pon-samples requires --fusion-pon"
    );

    info!(index = %args.index.display(), "loading fusion k-mer index...");
    let index = load_index(
        &args.index,
        k,
        args.max_kmer_copies,
        args.min_unique_anchor_reads > 0 || args.emit_unique_anchor,
        args.skip_version_check,
    )?;

    let kmer_blacklist: Option<std::collections::HashSet<u64>> = if let Some(ref bl_path) =
        args.fusion_kmer_blacklist
    {
        let bl_conn =
            Connection::open_in_memory().context("failed to open DuckDB for k-mer blacklist")?;
        let escaped = bl_path.display().to_string().replace('\'', "''");
        let threshold = args.min_kmer_blacklist_samples;
        let mut stmt = bl_conn
                .prepare(&format!(
                    "SELECT kmer_hash FROM read_parquet('{escaped}') WHERE n_pon_samples >= {threshold}"
                ))
                .context("failed to query k-mer blacklist Parquet")?;
        let mut rows = stmt.query([]).context("failed to read k-mer blacklist")?;
        let mut set = std::collections::HashSet::new();
        while let Some(row) = rows.next()? {
            let h: i64 = row.get(0)?;
            set.insert(h as u64);
        }
        info!(
            n_blacklisted = set.len(),
            min_pon_samples = threshold,
            path = %bl_path.display(),
            "loaded k-mer blacklist"
        );
        Some(set)
    } else {
        None
    };

    let mut reader = bam::Reader::from_path(&args.bam)
        .with_context(|| format!("failed to open BAM/CRAM: {}", args.bam.display()))?;
    if let Some(ref ref_path) = args.reference {
        reader
            .set_reference(ref_path)
            .context("failed to set reference for CRAM decoding")?;
    }

    let sample_id = match &args.sample_id {
        Some(id) => id.clone(),
        None => {
            let header_text = std::str::from_utf8(reader.header().as_bytes())
                .context("BAM header is not valid UTF-8")?;
            let mut found = None;
            'outer: for line in header_text.lines() {
                if line.starts_with("@RG") {
                    for field in line.split('\t') {
                        if let Some(sm) = field.strip_prefix("SM:") {
                            found = Some(sm.to_string());
                            break 'outer;
                        }
                    }
                }
            }
            found.context("--sample-id not provided and no SM tag found in BAM/CRAM @RG header")?
        }
    };

    info!(
        sample_id = %sample_id,
        bam = %args.bam.display(),
        "scanning BAM for fusion evidence..."
    );

    // qname bytes → list of ReadHit for reads that hit a gene
    let mut qname_to_reads: HashMap<Vec<u8>, Vec<ReadHit>> = HashMap::new();
    let mut reads_processed: u64 = 0;
    let mut reads_assigned: u64 = 0;

    // Total compressed size, used to estimate how far through the file we are.
    let bytes_total = std::fs::metadata(&args.bam).map(|m| m.len()).unwrap_or(0);
    // reader.tell() returns a BGZF virtual offset whose high 48 bits are the
    // compressed byte position — meaningful for BAM but not for CRAM (whose file
    // pointer is a cram_fd, not a BGZF stream). For CRAM we instead estimate
    // progress from the mapped genomic coordinate against the total reference
    // length, which assumes coordinate-sorted input (the production norm).
    let is_cram = args
        .bam
        .extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case("cram"));
    // Prefix sums of reference lengths: ref_offsets[tid] is the genomic position
    // of the start of target `tid` in a single concatenated coordinate space.
    let (ref_offsets, total_ref_len) = {
        let header = reader.header();
        let n = header.target_count();
        let mut offsets = Vec::with_capacity(n as usize);
        let mut acc: u64 = 0;
        for tid in 0..n {
            offsets.push(acc);
            acc += header.target_len(tid).unwrap_or(0);
        }
        (offsets, acc)
    };
    let scan_start = std::time::Instant::now();

    let mut record = bam::Record::new();
    while let Some(result) = reader.read(&mut record) {
        result.context("error reading BAM record")?;
        reads_processed += 1;

        let flags = record.flags();
        // Skip duplicates only. Secondary, supplementary, and unmapped reads are
        // intentionally included: k-mer matching is alignment-independent, so reads
        // the aligner couldn't place or split chimerically (the most informative for
        // fusion detection) still carry useful k-mer signal. Per-fragment dedup is
        // by qname, so multiple records sharing a qname collapse into one fragment.
        if flags & 0x400 != 0 {
            continue;
        }
        // mapq filter applies only to mapped records; unmapped reads have mapq=0
        // by convention and would otherwise always be excluded.
        let is_unmapped = flags & 0x4 != 0;
        if !is_unmapped && record.mapq() < args.min_mapq {
            continue;
        }

        let seq = record.seq().as_bytes();
        if let Some(mut rh) = assign_gene(
            &seq,
            &index,
            k,
            args.min_kmer_hits,
            false,
            kmer_blacklist.as_ref(),
        ) {
            rh.mapq = record.mapq();
            qname_to_reads
                .entry(record.qname().to_vec())
                .or_default()
                .push(rh);
            reads_assigned += 1;
        }

        if reads_processed.is_multiple_of(5_000_000) {
            let frac = if is_cram {
                // Coordinate-based: unmapped reads (tid < 0) cluster at the end of
                // a coordinate-sorted file, so treat them as ~100% complete.
                if total_ref_len > 0 {
                    let tid = record.tid();
                    let pos = if tid >= 0 && (tid as usize) < ref_offsets.len() {
                        ref_offsets[tid as usize] + record.pos().max(0) as u64
                    } else {
                        total_ref_len
                    };
                    Some(pos as f64 / total_ref_len as f64)
                } else {
                    None
                }
            } else if bytes_total > 0 {
                // Byte-based: high 48 bits of the BGZF virtual offset.
                Some((reader.tell() >> 16) as f64 / bytes_total as f64)
            } else {
                None
            };
            log_scan_progress(reads_processed, reads_assigned, frac, scan_start.elapsed());
        }
    }

    info!(reads_processed, reads_assigned, "BAM scan complete");

    // Aggregate: for each fragment with reads assigned to 2+ different genes,
    // count it as a fusion candidate and accumulate junction-coherence statistics.
    // Also cache qname → pair_key so secondary BAM passes can build fusion_qnames
    // directly without re-running fragment_top_pair on every fragment.
    let mut fusion_counts: HashMap<(u32, u32), FusionStats> = HashMap::new();
    let mut qname_to_pair: HashMap<Vec<u8>, (u32, u32)> = HashMap::new();
    let min_anchor = args.min_anchor_kmers;

    for (qname, read_hits) in &qname_to_reads {
        let Some(key) = fragment_top_pair(read_hits) else {
            continue;
        };
        qname_to_pair.insert(qname.clone(), key);
        let (ga, gb) = key;
        let min_mq = read_hits.iter().map(|rh| rh.mapq).min().unwrap_or(0);

        let mut n_spanning = 0u32;
        let mut fragment_coherent = false;
        for rh in read_hits {
            let (spanning, coherent) = read_coherence(rh, ga, gb, k, min_anchor);
            if spanning {
                n_spanning += 1;
            }
            if coherent {
                fragment_coherent = true;
            }
        }

        // Unique-anchor filter: is this fragment anchored by a genome-unique k-mer on
        // the higher-uniqueness partner? The anchor side is the gene of the pair with
        // more genome-unique k-mers (a repeat partner cannot anchor, so we require it
        // only on the side that can). Only evaluated when the filter is active.
        let fragment_unique_anchored = index
            .gene_unique_total
            .as_ref()
            .map(|tot| {
                let anchor = if tot.get(ga as usize).copied().unwrap_or(0)
                    >= tot.get(gb as usize).copied().unwrap_or(0)
                {
                    ga
                } else {
                    gb
                };
                read_hits.iter().any(|rh| read_unique_anchored(rh, anchor))
            })
            .unwrap_or(false);

        let entry = fusion_counts.entry(key).or_insert(FusionStats {
            count: 0,
            min_mapq: 255,
            n_spanning_reads: 0,
            n_coherent_fragments: 0,
            n_unique_anchored: 0,
        });
        entry.count += 1;
        entry.min_mapq = entry.min_mapq.min(min_mq);
        entry.n_spanning_reads += n_spanning;
        if fragment_coherent {
            entry.n_coherent_fragments += 1;
        }
        if fragment_unique_anchored {
            entry.n_unique_anchored += 1;
        }
    }

    // Was unique-anchor tracking active this run? Drives whether n_unique_anchored is a
    // real count or NA (the index carries genome_copies and we built per-gene unique totals).
    let unique_tracked = index.gene_unique_total.is_some();

    let mut records: Vec<FusionRecord> = fusion_counts
        .into_iter()
        .filter(|(_, s)| s.count >= args.min_supporting_reads)
        .map(|((ga, gb), s)| FusionRecord {
            pair_key: (ga, gb),
            sample_id: sample_id.clone(),
            gene_a: index.gene_names[ga as usize].clone(),
            gene_b: index.gene_names[gb as usize].clone(),
            chrom_a: index.gene_chroms[ga as usize].clone(),
            chrom_b: index.gene_chroms[gb as usize].clone(),
            supporting_reads: s.count as i32,
            min_mapq: s.min_mapq,
            n_spanning_reads: s.n_spanning_reads as i32,
            n_coherent_fragments: s.n_coherent_fragments as i32,
            n_pon_samples: 0,
            pon_total_samples: 0,
            max_pon_supporting_reads: None,
            filter: "PASS".to_string(),
            // Set in the 5'->3' reorientation pass once breakpoint stats are known.
            partner_order: "index".to_string(),
            n_unique_anchored: unique_tracked.then_some(s.n_unique_anchored as i32),
        })
        .collect();

    // Sort by support descending; break ties on the normalized gene-pair key so the
    // output row order is deterministic (fusion_counts is a HashMap).
    records.sort_by(|a, b| {
        b.supporting_reads
            .cmp(&a.supporting_reads)
            .then_with(|| a.pair_key.cmp(&b.pair_key))
    });

    // Panel-of-Normals annotation and optional filtering, before any output is written.
    if let Some(ref pon_db) = args.fusion_pon {
        annotate_fusion_pon(&mut records, pon_db, args.skip_version_check)?;
    }
    // Tag (don't drop) fusions seen in more PoN samples than the threshold. The
    // rows stay in the output with filter="pon" so downstream tooling can include
    // or exclude them; nothing is removed here.
    if let Some(max_pon) = args.max_pon_samples {
        let mut flagged = 0usize;
        for r in records.iter_mut() {
            if r.n_pon_samples as u32 > max_pon {
                r.filter = "pon".to_string();
                flagged += 1;
            }
        }
        info!(
            flagged,
            max_pon_samples = max_pon,
            "flagged fusions present in PoN (filter=pon)"
        );
    }

    if args.min_coherent_fragments > 0 {
        let before = records.len();
        records.retain(|r| r.n_coherent_fragments >= args.min_coherent_fragments as i32);
        let dropped = before - records.len();
        info!(
            dropped,
            min_coherent_fragments = args.min_coherent_fragments,
            "filtered fusions lacking coherent fragments"
        );
    }

    // The main Parquet/TSV are written AFTER the optional second BAM pass below, so
    // that the --max-breakpoint-std consensus filter (which needs spanning-read
    // breakpoint estimates from that pass) can set the `filter` column first.

    // Build fusion_label from the post-filter survivors so secondary outputs
    // (reads BAM, kmer-hits TSV, breakpoints TSV) only include evidence for
    // fusions that appear in the Parquet output.
    let fusion_label: HashMap<(u32, u32), String> = records
        .iter()
        .map(|r| (r.pair_key, fusion_pair_label(r.pair_key, &index.gene_names)))
        .collect();

    // Build fusion_qnames once from the cached qname→pair mapping. Filtering by
    // fusion_label ensures only post-filter survivors appear. Both secondary BAM
    // passes share this map, avoiding re-running fragment_top_pair per qname.
    // qname → label, plus qname → (gene_a, gene_b) index pair (same gene-index order
    // as the FX label) used to draw the per-read k-mer layout tag in the evidence BAM.
    let mut fusion_qnames: HashMap<Vec<u8>, String> = HashMap::new();
    let mut fusion_pairs: HashMap<Vec<u8>, (u32, u32)> = HashMap::new();
    for (qname, key) in qname_to_pair {
        if let Some(label) = fusion_label.get(&key) {
            fusion_pairs.insert(qname.clone(), key);
            fusion_qnames.insert(qname, label.clone());
        }
    }

    // Optional second BAM pass: write all reads from fusion-supporting fragments.
    // NOTE: this pass intentionally writes every record for a matching qname —
    // including duplicates and low-MAPQ mates — so the output BAM contains the
    // complete fragment for IGV inspection. Only the kmer-hits/breakpoints pass
    // below applies the same read-quality filters as the first pass, because that
    // pass uses reads for quantitative analysis (position estimates, k-mer counts).
    if let Some(ref reads_output) = args.reads_output {
        info!(
            n_fragments = fusion_qnames.len(),
            output = %reads_output.display(),
            "writing fusion-supporting reads (second BAM pass)..."
        );

        let mut reader2 = bam::Reader::from_path(&args.bam)
            .with_context(|| format!("failed to re-open BAM: {}", args.bam.display()))?;
        if let Some(ref ref_path) = args.reference {
            reader2.set_reference(ref_path)?;
        }

        let header = bam::Header::from_template(reader2.header());
        let mut writer = bam::Writer::from_path(reads_output, &header, bam::Format::Bam)
            .with_context(|| format!("failed to create output BAM: {}", reads_output.display()))?;

        // Progress + ETA for the second pass. It re-reads the whole BAM, so the
        // total record count from pass 1 (`reads_processed`) is the denominator.
        let pass2_total = reads_processed;
        let pass2_start = std::time::Instant::now();
        let mut reads_written: u64 = 0;
        let mut records_scanned: u64 = 0;
        for result in reader2.records() {
            let mut record = result.context("error reading BAM record in second pass")?;
            records_scanned += 1;
            if records_scanned.is_multiple_of(5_000_000) {
                let secs = pass2_start.elapsed().as_secs_f64();
                let frac = if pass2_total > 0 {
                    (records_scanned as f64 / pass2_total as f64).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                let eta_secs = if frac > 0.0 {
                    secs * (1.0 - frac) / frac
                } else {
                    f64::INFINITY
                };
                let reads_per_sec = if secs > 0.0 {
                    (records_scanned as f64 / secs) as u64
                } else {
                    0
                };
                info!(
                    records_scanned,
                    reads_written,
                    percent = format!("{:.1}", frac * 100.0),
                    eta = %format_duration(eta_secs),
                    reads_per_sec,
                    "second BAM pass progress (writing fusion reads)"
                );
            }
            if let Some(label) = fusion_qnames.get(record.qname()) {
                // Tag the read with the fusion it supports, e.g. FX:Z:EWSR1::FLI1.
                // Remove any pre-existing FX so re-runs don't duplicate the tag.
                let _ = record.remove_aux(b"FX");
                record
                    .push_aux(b"FX", bam::record::Aux::String(label))
                    .context("failed to add FX tag to BAM record")?;
                // FL: per-k-mer-window layout of this read against the fusion pair —
                // A=gene-A k-mer, B=gene-B k-mer (A/B = first/second gene in FX),
                // N=window masked by a non-ACGT base, .=k-mer matching neither gene.
                // Lowercase a/b mark windows rescued to a gene at edit distance 1 (a SNP
                // under the window); a single-N window that resolves uniquely is promoted to
                // a capital A/B. See fusion_layout_track / kmer::rescue_layout_track.
                if let Some(&(ga, gb)) = fusion_pairs.get(record.qname()) {
                    let track = fusion_layout_track(
                        &record.seq().as_bytes(),
                        k,
                        &index.kmer_to_gene,
                        ga,
                        gb,
                    );
                    let _ = record.remove_aux(b"FL");
                    record
                        .push_aux(b"FL", bam::record::Aux::String(&track))
                        .context("failed to add FL tag to BAM record")?;
                }
                writer
                    .write(&record)
                    .context("failed to write BAM record")?;
                reads_written += 1;
            }
        }
        // Close before indexing.
        drop(writer);

        info!(reads_written, output = %reads_output.display(), "BAM written");

        bam::index::build(reads_output, None, bam::index::Type::Bai, 0)
            .with_context(|| format!("failed to index output BAM: {}", reads_output.display()))?;
        info!(output = %reads_output.display(), "BAI index written");
    }

    // Per-fusion breakpoint consensus stats, populated by the second BAM pass below
    // when breakpoints or the consensus filter are requested. Keyed by fusion label.
    let mut breakpoint_stats: HashMap<String, BreakpointStats> = HashMap::new();

    // Optional combined second BAM pass: produces kmer-hits TSV and/or breakpoints TSV,
    // and/or the spanning-read breakpoint estimates the --max-breakpoint-std filter needs.
    // All consumers share one BAM scan to avoid re-reading the file.
    if args.kmer_hits_output.is_some()
        || args.breakpoints_output.is_some()
        || args.max_breakpoint_std.is_some()
        || args.min_breakpoint_distance.is_some()
    {
        use std::io::{BufWriter, Write};

        info!(
            n_fragments = fusion_qnames.len(),
            "writing detail outputs (second BAM pass)..."
        );

        let mut reader2 = bam::Reader::from_path(&args.bam)
            .with_context(|| format!("failed to re-open BAM: {}", args.bam.display()))?;
        if let Some(ref ref_path) = args.reference {
            reader2.set_reference(ref_path)?;
        }

        let target_names: Vec<String> = {
            let hdr = reader2.header();
            (0..hdr.target_count())
                .map(|i| {
                    std::str::from_utf8(hdr.tid2name(i))
                        .unwrap_or("?")
                        .to_string()
                })
                .collect()
        };

        // kmer-hits TSV writer (optional).
        let mut kmer_hits_writer: Option<BufWriter<std::fs::File>> = if let Some(ref path) =
            args.kmer_hits_output
        {
            let file = std::fs::File::create(path)
                .with_context(|| format!("failed to create kmer hits TSV: {}", path.display()))?;
            let mut w = BufWriter::new(file);
            writeln!(w, "fusion\tsample_id\tread_name\tread_end\tchrom\tpos\tgene_matched\tkmer_pos_in_read\tkmer_hash\tkmer_seq")?;
            Some(w)
        } else {
            None
        };

        // Breakpoint accumulators (optional): one per passing fusion. Needed for the
        // breakpoints TSV and/or the --max-breakpoint-std consensus filter.
        let mut accumulators: Option<HashMap<String, BreakpointAccumulator>> =
            if args.breakpoints_output.is_some()
                || args.max_breakpoint_std.is_some()
                || args.min_breakpoint_distance.is_some()
            {
                let mut map: HashMap<String, BreakpointAccumulator> = HashMap::new();
                for (&(ga, gb), label) in &fusion_label {
                    map.insert(
                        label.clone(),
                        BreakpointAccumulator {
                            gene_a_idx: ga,
                            gene_b_idx: gb,
                            gene_a_chrom_votes: HashMap::new(),
                            gene_b_chrom_votes: HashMap::new(),
                            spanning_reads: Vec::new(),
                            gene_a_single_reads: Vec::new(),
                            gene_b_single_reads: Vec::new(),
                        },
                    );
                }
                Some(map)
            } else {
                None
            };

        let detail_total = reads_processed;
        let detail_start = std::time::Instant::now();
        let mut rows_written: u64 = 0;
        let mut detail_reads_processed: u64 = 0;
        for result in reader2.records() {
            let record = result.context("error reading BAM record in detail pass")?;
            detail_reads_processed += 1;
            if detail_reads_processed.is_multiple_of(5_000_000) {
                let secs = detail_start.elapsed().as_secs_f64();
                let frac = if detail_total > 0 {
                    (detail_reads_processed as f64 / detail_total as f64).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                let eta_secs = if frac > 0.0 {
                    secs * (1.0 - frac) / frac
                } else {
                    f64::INFINITY
                };
                info!(
                    detail_reads_processed,
                    rows_written,
                    percent = format!("{:.1}", frac * 100.0),
                    eta = %format_duration(eta_secs),
                    "detail pass progress"
                );
            }
            let Some(fusion_label_str) = fusion_qnames.get(record.qname()) else {
                continue;
            };
            // Apply the same filters as the first pass so kmer-hits and breakpoint
            // estimates are consistent with supporting_reads counts in the Parquet.
            let flags = record.flags();
            if flags & 0x400 != 0 {
                continue;
            }
            let is_unmapped = flags & 0x4 != 0;
            if !is_unmapped && record.mapq() < args.min_mapq {
                continue;
            }
            let seq = record.seq().as_bytes();
            let Some(rh) = assign_gene(
                &seq,
                &index,
                k,
                args.min_kmer_hits,
                true,
                kmer_blacklist.as_ref(),
            ) else {
                continue;
            };
            let read_end = if flags & 0x40 != 0 { "R1" } else { "R2" };
            let tid = record.tid();
            let chrom_str: &str = if tid >= 0 {
                target_names
                    .get(tid as usize)
                    .map(|s| s.as_str())
                    .unwrap_or("*")
            } else {
                "*"
            };
            let pos1 = record.pos() + 1; // 1-based alignment start

            // kmer-hits TSV rows.
            if let Some(ref mut w) = kmer_hits_writer {
                let read_name = std::str::from_utf8(record.qname()).unwrap_or("?");
                for &(kmer_hash, gene_idx, kmer_pos_in_read) in &rh.kmer_hits {
                    let gene_name = index
                        .gene_names
                        .get(gene_idx as usize)
                        .map(|s| s.as_str())
                        .unwrap_or("?");
                    let kmer_seq = decode_kmer(kmer_hash, k);
                    writeln!(
                        w,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        fusion_label_str,
                        sample_id,
                        read_name,
                        read_end,
                        chrom_str,
                        pos1,
                        gene_name,
                        kmer_pos_in_read,
                        kmer_hash as i64,
                        kmer_seq
                    )?;
                    rows_written += 1;
                }
            }

            // Breakpoint accumulation.
            if let Some(ref mut acc_map) = accumulators {
                if let Some(acc) = acc_map.get_mut(fusion_label_str) {
                    // Partition k-mer hit positions by gene.
                    let mut a_positions: Vec<usize> = Vec::new();
                    let mut b_positions: Vec<usize> = Vec::new();
                    for &(_, gene_idx, kmer_pos_in_read) in &rh.kmer_hits {
                        if gene_idx == acc.gene_a_idx {
                            a_positions.push(kmer_pos_in_read);
                        } else if gene_idx == acc.gene_b_idx {
                            b_positions.push(kmer_pos_in_read);
                        }
                    }

                    if chrom_str == "*" {
                        // Unmapped: skip for chromosome votes and breakpoint estimates.
                    } else if !a_positions.is_empty() && b_positions.is_empty() {
                        *acc.gene_a_chrom_votes
                            .entry(chrom_str.to_string())
                            .or_insert(0) += 1;
                        // Supplementary alignments of chimeric reads appear as separate
                        // BAM records with single-gene k-mers at off-junction positions.
                        // Exclude them so only concordant-pair primary reads contribute
                        // breakpoint estimates.
                        if flags & 0x800 == 0 {
                            acc.gene_a_single_reads.push(SingleGeneReadData {
                                chrom: chrom_str.to_string(),
                                aln_start: pos1 - 1, // convert to 0-based
                                aln_end: record.cigar().end_pos(),
                            });
                        }
                    } else if a_positions.is_empty() && !b_positions.is_empty() {
                        *acc.gene_b_chrom_votes
                            .entry(chrom_str.to_string())
                            .or_insert(0) += 1;
                        if flags & 0x800 == 0 {
                            acc.gene_b_single_reads.push(SingleGeneReadData {
                                chrom: chrom_str.to_string(),
                                aln_start: pos1 - 1, // convert to 0-based
                                aln_end: record.cigar().end_pos(),
                            });
                        }
                    } else if !a_positions.is_empty() && !b_positions.is_empty() {
                        // Spanning read.
                        let last_a = *a_positions.iter().max().unwrap();
                        let first_a = *a_positions.iter().min().unwrap();
                        let last_b = *b_positions.iter().max().unwrap();
                        let first_b = *b_positions.iter().min().unwrap();
                        acc.spanning_reads.push(SpanningReadData {
                            chrom: chrom_str.to_string(),
                            pos: pos1,
                            last_a,
                            first_a,
                            last_b,
                            first_b,
                        });
                    }
                }
            }
        }

        if let Some(ref path) = args.kmer_hits_output {
            info!(rows_written, output = %path.display(), "kmer hits TSV written");
        }

        // Reduce each fusion's accumulated spanning reads to breakpoint consensus
        // stats. Computed once here and reused by both the breakpoints TSV and the
        // --max-breakpoint-std filter applied after this pass.
        if let Some(ref acc_map) = accumulators {
            for (label, acc) in acc_map.iter() {
                breakpoint_stats.insert(label.clone(), compute_breakpoint_stats(acc, &index, k));
            }
        }

        // Write breakpoints TSV.
        if let Some(ref bp_path) = args.breakpoints_output {
            let file = std::fs::File::create(bp_path).with_context(|| {
                format!("failed to create breakpoints TSV: {}", bp_path.display())
            })?;
            let mut w = BufWriter::new(file);
            writeln!(w, "fusion\tgene_a\tchrom_a\tbreakpoint_a\tbp_a_n\tbp_a_std\tgene_b\tchrom_b\tbreakpoint_b\tbp_b_n\tbp_b_std\tn_spanning_reads")?;

            let mut fusions_sorted: Vec<&String> = breakpoint_stats.keys().collect();
            fusions_sorted.sort_unstable();

            let fmt_opt_f64 = |v: Option<f64>| {
                v.map(|x| format!("{:.1}", x))
                    .unwrap_or_else(|| "NA".to_string())
            };
            let fmt_opt_i64 = |v: Option<f64>| {
                v.map(|x| format!("{}", x as i64))
                    .unwrap_or_else(|| "NA".to_string())
            };

            for label in fusions_sorted {
                let s = &breakpoint_stats[label];
                let gene_a = &index.gene_names[s.gene_a_idx as usize];
                let gene_b = &index.gene_names[s.gene_b_idx as usize];
                // Present the A/B columns (and the fusion label) in biological 5'->3'
                // order when known, mirroring the fusions output; index order otherwise.
                let flip = s.gene_a_is_5prime == Some(false);
                let (g5, c5, b5, n5, s5, g3, c3, b3, n3, s3) = if flip {
                    (
                        gene_b, &s.chrom_b, s.bp_b, s.bp_b_n, s.bp_b_std, gene_a, &s.chrom_a,
                        s.bp_a, s.bp_a_n, s.bp_a_std,
                    )
                } else {
                    (
                        gene_a, &s.chrom_a, s.bp_a, s.bp_a_n, s.bp_a_std, gene_b, &s.chrom_b,
                        s.bp_b, s.bp_b_n, s.bp_b_std,
                    )
                };
                let oriented_label = if flip {
                    format!("{}::{}", gene_b, gene_a)
                } else {
                    label.clone()
                };
                writeln!(
                    w,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    oriented_label,
                    g5,
                    c5,
                    fmt_opt_i64(b5),
                    n5,
                    fmt_opt_f64(s5),
                    g3,
                    c3,
                    fmt_opt_i64(b3),
                    n3,
                    fmt_opt_f64(s3),
                    s.n_spanning,
                )?;
            }
            info!(output = %bp_path.display(), "breakpoints TSV written");
        }
    }

    // Breakpoint-consensus filter: a real fusion's spanning reads converge on one
    // junction base (tight std) with depth on both sides; paralog/PCR-chimera
    // artifacts splice at scattered positions. Tag failers filter="chimera" (kept,
    // not dropped) without clobbering a prior PoN tag.
    if let Some(max_std) = args.max_breakpoint_std {
        let min_n = args.min_breakpoint_reads as usize;
        let mut flagged = 0usize;
        for r in records.iter_mut() {
            if r.filter != "PASS" {
                continue;
            }
            let label = fusion_pair_label(r.pair_key, &index.gene_names);
            let rec_mapq = r.min_mapq;
            let passes = breakpoint_stats.get(&label).is_some_and(|s| {
                let tight = s.bp_a_n >= min_n
                    && s.bp_b_n >= min_n
                    && s.bp_a_std.is_some_and(|v| v <= max_std)
                    && s.bp_b_std.is_some_and(|v| v <= max_std);
                // A real high-depth junction can exceed `max_std` (splice-isoform +
                // transition-point spread); strong support on both sides rescues it
                // without admitting low-support artifacts (which never reach this tier).
                // Asymmetric anchor: one side a tight well-covered junction, the partner
                // seen by >=1 read — recovers low-input junctions captured well on only one
                // side (single-sample; gated on MAPQ to exclude multi-mapper paralog leakage).
                let asym = args.asym_anchor_reads.is_some_and(|n| {
                    rec_mapq >= args.asym_anchor_mapq
                        && s.asym_anchored(n as usize, args.asym_anchor_std)
                });
                tight || s.strong_support() || asym
            });
            if !passes {
                r.filter = "chimera".to_string();
                flagged += 1;
            }
        }
        info!(
            flagged,
            max_breakpoint_std = max_std,
            min_breakpoint_reads = args.min_breakpoint_reads,
            "tagged fusions lacking breakpoint consensus (filter=chimera)"
        );
    }

    // Same-locus / adjacency filter: both breakpoints on one chromosome within a tiny
    // window is the signature of single-locus paralog leakage (e.g. unindexed GNA12
    // reads split between GNA13/GNA11, both breakpoints landing at the GNA12 locus).
    // Such calls have a TIGHT, reproducible breakpoint, so the consensus filter passes
    // them; this catches them on geometry instead. Tag filter="samelocus" (kept, not
    // dropped); only touch calls still PASS so prior pon/chimera tags win.
    if let Some(min_dist) = args.min_breakpoint_distance {
        let mut flagged = 0usize;
        for r in records.iter_mut() {
            if r.filter != "PASS" {
                continue;
            }
            let label = fusion_pair_label(r.pair_key, &index.gene_names);
            let same_locus = breakpoint_stats.get(&label).is_some_and(|s| {
                // Single-locus leakage is spanning-read dominated at one locus and never
                // accrues strong INDEPENDENT (concordant) support on both partners. A
                // co-located call with strong concordant support is a real junction seen
                // through concordant pairs (e.g. a breakpoint buried in a repeat) — exempt
                // it. Total-support strength is NOT enough here: a high-volume leakage
                // locus can exceed the read threshold on spanning reads alone.
                if s.strong_concordant_support() {
                    return false;
                }
                match (s.bp_a, s.bp_b) {
                    (Some(a), Some(b)) => s.chrom_a == s.chrom_b && (a - b).abs() < min_dist as f64,
                    _ => false,
                }
            });
            if same_locus {
                r.filter = "samelocus".to_string();
                flagged += 1;
            }
        }
        info!(
            flagged,
            min_breakpoint_distance = min_dist,
            "tagged fusions whose breakpoints co-localize to one locus (filter=samelocus)"
        );
    }

    // Report partners in biological 5'->3' order. When orientation inference found that
    // gene_b (pair_key.1) is the 5' partner, swap the displayed gene/chrom so gene_a is 5'
    // and gene_b is 3'. Undetermined (legacy index, unknown strand, no read consensus)
    // leaves gene-index order. The internal pair_key and per-read FX/FL tags are unchanged
    // — those stay index-ordered read-level evidence.
    for r in records.iter_mut() {
        let label = fusion_pair_label(r.pair_key, &index.gene_names);
        match breakpoint_stats
            .get(&label)
            .and_then(|s| s.gene_a_is_5prime)
        {
            Some(true) => r.partner_order = "5to3".to_string(),
            Some(false) => {
                std::mem::swap(&mut r.gene_a, &mut r.gene_b);
                std::mem::swap(&mut r.chrom_a, &mut r.chrom_b);
                r.partner_order = "5to3".to_string();
            }
            None => {} // undetermined → leave gene-index order, partner_order = "index"
        }
    }

    // Unique-anchor specificity filter (applied last): a call stays filter="PASS" only if
    // enough of its supporting fragments are anchored by a genome-unique k-mer on the
    // higher-uniqueness partner. Repeat/segdup/pseudogene partners fabricate high-support
    // calls in nearly every sample, but those lack a unique-anchored junction. Only touches
    // calls still PASS (prior pon/chimera/samelocus tags win); tags, never drops.
    if args.min_unique_anchor_reads > 0 {
        let mut flagged = 0usize;
        for r in records.iter_mut() {
            if r.filter != "PASS" {
                continue;
            }
            // unique_tracked is guaranteed here (min_unique_anchor_reads > 0 enables it).
            let n = r.n_unique_anchored.unwrap_or(0).max(0) as u32;
            if n < args.min_unique_anchor_reads {
                r.filter = "no_unique_anchor".to_string();
                flagged += 1;
            }
        }
        info!(
            flagged,
            min_unique_anchor_reads = args.min_unique_anchor_reads,
            "tagged fusions lacking a genome-unique anchor (filter=no_unique_anchor)"
        );
    }

    info!(n_fusions = records.len(), "fusion candidates identified");
    write_fusion_parquet(&records, &args.output)?;
    info!(output = %args.output.display(), "fusion candidates written");

    if let Some(ref tsv_path) = args.tsv_output {
        write_fusion_tsv(&records, tsv_path)?;
        info!(output = %tsv_path.display(), "fusion TSV written");
    }

    // All outputs are flushed by this point. The in-memory index holds tens of
    // millions of k-mers; dropping it at scope exit can take several seconds, which
    // otherwise looks like a hang after the last "written" line. Release it
    // explicitly with bracketing logs so the pause is attributable.
    info!(
        n_kmers = index.kmer_to_gene.len(),
        "all outputs written; releasing in-memory index..."
    );
    drop(index);
    info!("fusion detection complete");

    Ok(())
}

// ─── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fusion_layout_track_blocks_n_and_classifies() {
        // k=4. Read: 4 bp of gene-A sequence, then an N, then 4 bp of gene-B.
        // Build a tiny index mapping the exact k-mers present to gene 0 (A) and 1 (B).
        let seq = b"ACGTNTTGG"; // positions: ACGT | N | TTGG
        let k = 4;
        let mut kmer_to_gene: HashMap<u64, u32> = HashMap::new();
        // gene 0 (A): canonical k-mer of "ACGT"
        for (kmer, _) in kmer_iter(b"ACGT", k) {
            kmer_to_gene.insert(kmer, 0);
        }
        // gene 1 (B): canonical k-mer of "TTGG"
        for (kmer, _) in kmer_iter(b"TTGG", k) {
            kmer_to_gene.insert(kmer, 1);
        }
        let track = fusion_layout_track(seq, k, &kmer_to_gene, 0, 1);
        // 9 bp, k=4 -> 6 windows (starts 0..5).
        // start0 ACGT -> A; starts 1..4 contain the N at index 4 -> N; start5 TTGG -> B.
        assert_eq!(track, "ANNNNB");
    }

    /// Build a ReadHit with the given primary/second gene assignments. Positions
    /// and counts default sensibly; tests override what they exercise.
    fn mk_read(
        gene_idx: u32,
        g1: (u32, u32, u32),                 // (count, min, max)
        gene2: Option<(u32, u32, u32, u32)>, // (idx, count, min, max)
    ) -> ReadHit {
        let (g2_idx, g2_count, g2_min, g2_max) = match gene2 {
            Some((i, c, mn, mx)) => (Some(i), c, mn, mx),
            None => (None, 0, 0, 0),
        };
        ReadHit {
            gene_idx,
            gene1_kmer_count: g1.0,
            gene1_kmer_min: g1.1,
            gene1_kmer_max: g1.2,
            gene1_unique_count: 0,
            gene2_idx: g2_idx,
            gene2_kmer_count: g2_count,
            gene2_kmer_min: g2_min,
            gene2_kmer_max: g2_max,
            gene2_unique_count: 0,
            mapq: 60,
            kmer_hits: Vec::new(),
        }
    }

    #[test]
    fn median_handles_parity_and_empty() {
        assert_eq!(median_i64(&mut []), None);
        assert_eq!(median_i64(&mut [5]), Some(5.0));
        assert_eq!(median_i64(&mut [3, 1, 2]), Some(2.0)); // odd, sorted internally
        assert_eq!(median_i64(&mut [4, 1, 3, 2]), Some(2.5)); // even → mean of middle two
    }

    #[test]
    fn std_dev_requires_two_values() {
        assert_eq!(std_dev_i64(&[]), None);
        assert_eq!(std_dev_i64(&[7]), None);
        // sample std dev of {2,4}: mean 3, variance (1+1)/(2-1)=2
        let s = std_dev_i64(&[2, 4]).unwrap();
        assert!((s - 2f64.sqrt()).abs() < 1e-9);
    }

    #[test]
    fn decode_kmer_roundtrips_known_values() {
        assert_eq!(decode_kmer(0, 4), "AAAA");
        assert_eq!(decode_kmer(255, 4), "TTTT");
        // Forward encoding of ACGT (A=0,C=1,G=2,T=3) = 0b00_01_10_11 = 27.
        assert_eq!(decode_kmer(27, 4), "ACGT");
    }

    #[test]
    fn fusion_pair_label_uses_index_order() {
        let names = vec!["A".to_string(), "B".to_string(), "C".to_string()];
        assert_eq!(fusion_pair_label((0, 2), &names), "A::C");
        assert_eq!(fusion_pair_label((1, 2), &names), "B::C");
    }

    // ── unique-anchor filter predicate ──────────────────────────────────────
    /// ReadHit with explicit unique counts on each side (gene1 = primary).
    fn mk_read_uniq(g1: u32, g1_uniq: u32, g2: Option<u32>, g2_uniq: u32) -> ReadHit {
        ReadHit {
            gene_idx: g1,
            gene1_kmer_count: 5,
            gene1_kmer_min: 0,
            gene1_kmer_max: 4,
            gene1_unique_count: g1_uniq,
            gene2_idx: g2,
            gene2_kmer_count: if g2.is_some() { 5 } else { 0 },
            gene2_kmer_min: 0,
            gene2_kmer_max: 0,
            gene2_unique_count: g2_uniq,
            mapq: 60,
            kmer_hits: Vec::new(),
        }
    }

    #[test]
    fn unique_anchor_true_on_primary_side() {
        // anchor gene 0 sits on gene1 with a unique k-mer → anchored.
        let rh = mk_read_uniq(0, 2, Some(1), 0);
        assert!(read_unique_anchored(&rh, 0));
    }

    #[test]
    fn unique_anchor_false_when_zero_unique_kmers() {
        // gene 0 is present but only via multi-copy (repeat) k-mers → not anchored.
        let rh = mk_read_uniq(0, 0, Some(1), 0);
        assert!(!read_unique_anchored(&rh, 0));
    }

    #[test]
    fn unique_anchor_checks_the_gene2_side() {
        // anchor gene 1 sits on the gene2 side with a unique k-mer → anchored.
        let rh = mk_read_uniq(0, 0, Some(1), 3);
        assert!(read_unique_anchored(&rh, 1));
        // ...but gene 0 (primary) has no unique k-mers → not anchored on that side.
        assert!(!read_unique_anchored(&rh, 0));
    }

    #[test]
    fn unique_anchor_false_when_anchor_gene_absent() {
        // anchor gene 7 is on neither side of the read.
        let rh = mk_read_uniq(0, 5, Some(1), 5);
        assert!(!read_unique_anchored(&rh, 7));
    }

    // ── 5'->3' orientation inference ─────────────────────────────────────────
    fn mk_orient_index(strands: &[char], chroms: &[&str], has_strand: bool) -> FusionIndex {
        let n = strands.len();
        FusionIndex {
            kmer_to_gene: HashMap::new(),
            gene_names: (0..n).map(|i| format!("G{i}")).collect(),
            gene_chroms: chroms.iter().map(|c| c.to_string()).collect(),
            gene_strands: strands.to_vec(),
            // Two well-separated bodies on each chrom for the same-chrom test.
            gene_tx_start: (0..n).map(|i| (i as i64) * 1000).collect(),
            gene_tx_end: (0..n).map(|i| (i as i64) * 1000 + 500).collect(),
            has_strand,
            unique_kmers: None,
            gene_unique_total: None,
        }
    }

    fn mk_acc(spans: Vec<SpanningReadData>) -> BreakpointAccumulator {
        BreakpointAccumulator {
            gene_a_idx: 0,
            gene_b_idx: 1,
            gene_a_chrom_votes: HashMap::new(),
            gene_b_chrom_votes: HashMap::new(),
            spanning_reads: spans,
            gene_a_single_reads: Vec::new(),
            gene_b_single_reads: Vec::new(),
        }
    }

    // a_before_b is set via first_a < first_b; chrom/pos pick the aligned gene.
    fn mk_span(chrom: &str, pos: i64, a_before_b: bool) -> SpanningReadData {
        let (first_a, first_b) = if a_before_b { (0, 30) } else { (30, 0) };
        SpanningReadData {
            chrom: chrom.to_string(),
            pos,
            first_a,
            last_a: first_a + 5,
            first_b,
            last_b: first_b + 5,
        }
    }

    // `n` identical spanning reads, to clear ORIENT_MIN_VOTES.
    fn mk_spans(chrom: &str, pos: i64, a_before_b: bool, n: usize) -> Vec<SpanningReadData> {
        (0..n).map(|_| mk_span(chrom, pos, a_before_b)).collect()
    }

    #[test]
    fn orient_none_without_strand() {
        let idx = mk_orient_index(&['+', '-'], &["chr1", "chr2"], false);
        let acc = mk_acc(mk_spans("chr1", 100, true, 12));
        assert_eq!(infer_five_prime(&acc, &idx, "chr1", "chr2"), None);
    }

    #[test]
    fn orient_gene_a_aligned_plus_block_first_is_5prime() {
        // gene_a on '+', its block first in SEQ → retains promoter-proximal 5' → gene_a is 5'.
        let idx = mk_orient_index(&['+', '-'], &["chr1", "chr2"], true);
        let acc = mk_acc(mk_spans("chr1", 100, true, 12));
        assert_eq!(infer_five_prime(&acc, &idx, "chr1", "chr2"), Some(true));
    }

    #[test]
    fn orient_gene_a_aligned_minus_block_first_is_3prime() {
        // gene_a on '-', block first → its 3' portion is retained → gene_a is the 3' partner.
        let idx = mk_orient_index(&['-', '+'], &["chr1", "chr2"], true);
        let acc = mk_acc(mk_spans("chr1", 100, true, 12));
        assert_eq!(infer_five_prime(&acc, &idx, "chr1", "chr2"), Some(false));
    }

    #[test]
    fn orient_gene_b_aligned_vote() {
        // Read aligns to gene_b (chr2), gene_b on '+', block NOT first (b before a in SEQ)
        // → gene_b retains its 5' → gene_b is 5' → gene_a is not 5'.
        let idx = mk_orient_index(&['+', '+'], &["chr1", "chr2"], true);
        let acc = mk_acc(mk_spans("chr2", 1000, false, 12));
        assert_eq!(infer_five_prime(&acc, &idx, "chr1", "chr2"), Some(false));
    }

    #[test]
    fn orient_majority_wins_with_minority_present() {
        // 10 gene_a-5' votes + 2 gene_b-5' votes (a reciprocal-junction-like minority).
        let idx = mk_orient_index(&['+', '-'], &["chr1", "chr2"], true);
        let mut spans = mk_spans("chr1", 100, true, 10); // gene_a 5'
        spans.extend(mk_spans("chr1", 100, false, 2)); // gene_b 5'
        assert_eq!(
            infer_five_prime(&mk_acc(spans), &idx, "chr1", "chr2"),
            Some(true)
        );
    }

    #[test]
    fn orient_abstains_below_min_votes() {
        // A decisive but low-support vote (capture-limited call) must abstain.
        let idx = mk_orient_index(&['+', '-'], &["chr1", "chr2"], true);
        let acc = mk_acc(mk_spans("chr1", 100, true, 6));
        assert_eq!(infer_five_prime(&acc, &idx, "chr1", "chr2"), None);
    }

    #[test]
    fn orient_abstains_on_coin_flip() {
        // Enough votes but split ~50/50 (below the fraction floor) → undetermined.
        let idx = mk_orient_index(&['+', '-'], &["chr1", "chr2"], true);
        let mut spans = mk_spans("chr1", 100, true, 6);
        spans.extend(mk_spans("chr1", 100, false, 6));
        assert_eq!(infer_five_prime(&mk_acc(spans), &idx, "chr1", "chr2"), None);
    }

    #[test]
    fn orient_same_chrom_disambiguates_by_position() {
        // Both genes on chr1; reads at pos in gene_a's body (0..500), gene_a '+', block first.
        let idx = mk_orient_index(&['+', '+'], &["chr1", "chr1"], true);
        let acc = mk_acc(mk_spans("chr1", 100, true, 12));
        assert_eq!(infer_five_prime(&acc, &idx, "chr1", "chr1"), Some(true));
    }

    #[test]
    fn top_pair_from_two_primary_genes() {
        let reads = vec![mk_read(0, (5, 0, 10), None), mk_read(3, (4, 0, 8), None)];
        assert_eq!(fragment_top_pair(&reads), Some((0, 3)));
    }

    #[test]
    fn top_pair_breaks_ties_by_ascending_index() {
        // Three genes, one primary vote each → top two by ascending index.
        let reads = vec![
            mk_read(5, (3, 0, 6), None),
            mk_read(1, (3, 0, 6), None),
            mk_read(3, (3, 0, 6), None),
        ];
        assert_eq!(fragment_top_pair(&reads), Some((1, 3)));
    }

    #[test]
    fn top_pair_uses_supplemental_votes_when_single_primary() {
        // Single primary gene (0); partner chosen from gene2 k-mer evidence.
        let reads = vec![
            mk_read(0, (8, 0, 20), Some((2, 6, 30, 50))),
            mk_read(0, (7, 0, 18), Some((2, 5, 28, 48))),
        ];
        assert_eq!(fragment_top_pair(&reads), Some((0, 2)));
    }

    #[test]
    fn top_pair_none_for_single_gene_no_partner() {
        let reads = vec![mk_read(0, (4, 0, 9), None)];
        assert_eq!(fragment_top_pair(&reads), None);
    }

    #[test]
    fn coherence_disjoint_blocks_are_coherent() {
        // gene0 k-mers in [0,10], gene1 k-mers in [40,50], k=20 → A ends (30) before B (40).
        let rh = mk_read(0, (5, 0, 10), Some((1, 5, 40, 50)));
        assert_eq!(read_coherence(&rh, 0, 1, 20, 3), (true, true));
    }

    #[test]
    fn coherence_interleaved_blocks_not_coherent() {
        // Overlapping position ranges → homology artifact, spanning but not coherent.
        let rh = mk_read(0, (5, 0, 50), Some((1, 5, 10, 60)));
        assert_eq!(read_coherence(&rh, 0, 1, 20, 3), (true, false));
    }

    #[test]
    fn coherence_requires_anchor_on_both_sides() {
        // gene1 has only 2 k-mers but min_anchor is 3 → spanning, not coherent.
        let rh = mk_read(0, (5, 0, 10), Some((1, 2, 40, 50)));
        assert_eq!(read_coherence(&rh, 0, 1, 20, 3), (true, false));
    }

    #[test]
    fn coherence_false_when_pair_absent() {
        // Read only hit gene0 → not spanning for the (0,1) pair.
        let rh = mk_read(0, (5, 0, 10), None);
        assert_eq!(read_coherence(&rh, 0, 1, 20, 3), (false, false));
        // Read's second gene is not the partner.
        let rh2 = mk_read(0, (5, 0, 10), Some((9, 5, 40, 50)));
        assert_eq!(read_coherence(&rh2, 0, 1, 20, 3), (false, false));
    }
}

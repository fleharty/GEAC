use std::collections::HashMap;
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
use crate::kmer::kmer_iter;

// ─── Index loading ────────────────────────────────────────────────────────────

struct FusionIndex {
    kmer_to_gene: HashMap<u64, u32>,
    gene_names: Vec<String>,
    gene_chroms: Vec<String>,
    /// Genome-wide copy number per k-mer. Populated only when a copy filter is
    /// requested and the index carries `genome_copies`; empty otherwise.
    kmer_to_copies: HashMap<u64, u8>,
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

fn load_index(index_path: &Path, max_kmer_copies: Option<u32>) -> Result<FusionIndex> {
    let conn = Connection::open(index_path)
        .with_context(|| format!("failed to open fusion index: {}", index_path.display()))?;

    let mut stmt =
        conn.prepare("SELECT gene_index, gene_name, chrom FROM genes ORDER BY gene_index")?;
    let mut gene_names: Vec<String> = Vec::new();
    let mut gene_chroms: Vec<String> = Vec::new();

    let mut rows = stmt.query([])?;
    while let Some(row) = rows.next()? {
        let _idx: u32 = row.get(0)?;
        gene_names.push(row.get(1)?);
        gene_chroms.push(row.get(2)?);
    }
    info!(n_genes = gene_names.len(), "loaded genes from index");

    // Only read genome_copies when a copy filter is active. This keeps the default
    // path identical (and compatible with indexes built before the column existed).
    let load_copies = max_kmer_copies.is_some();
    if load_copies {
        anyhow::ensure!(
            kmers_has_copies_column(&conn)?,
            "--max-kmer-copies requires an index built with --check-genome-uniqueness \
             (the index has no genome_copies column); rebuild the index or drop the flag"
        );
        // Fail fast on indexes whose genome_copies are all NULL (column present but
        // the genome pass never ran) before loading the whole kmers table into RAM.
        let n_total: i64 =
            conn.query_row("SELECT count(*) FROM kmers", [], |r| r.get(0))?;
        let n_with_copies: i64 =
            conn.query_row("SELECT count(genome_copies) FROM kmers", [], |r| r.get(0))?;
        anyhow::ensure!(
            n_total == 0 || n_with_copies > 0,
            "--max-kmer-copies requires per-k-mer genome copy counts, but the index has \
             none (genome_copies is NULL for every k-mer); rebuild with --check-genome-uniqueness"
        );
    }

    let mut kmer_to_gene: HashMap<u64, u32> = HashMap::new();
    let mut kmer_to_copies: HashMap<u64, u8> = HashMap::new();

    if load_copies {
        let mut stmt = conn.prepare("SELECT kmer_hash, gene_index, genome_copies FROM kmers")?;
        let mut rows = stmt.query([])?;
        let mut n_null: u64 = 0;
        while let Some(row) = rows.next()? {
            let kmer_i64: i64 = row.get(0)?;
            let gene_idx: u32 = row.get(1)?;
            let copies: Option<i32> = row.get(2)?;
            kmer_to_gene.insert(kmer_i64 as u64, gene_idx);
            match copies {
                Some(c) => {
                    kmer_to_copies.insert(kmer_i64 as u64, c.clamp(0, u8::MAX as i32) as u8);
                }
                None => n_null += 1,
            }
        }
        if n_null > 0 {
            info!(n_null, "some k-mers have NULL genome_copies and will be excluded by the copy filter");
        }
    } else {
        let mut stmt = conn.prepare("SELECT kmer_hash, gene_index FROM kmers")?;
        let mut rows = stmt.query([])?;
        while let Some(row) = rows.next()? {
            let kmer_i64: i64 = row.get(0)?;
            let gene_idx: u32 = row.get(1)?;
            kmer_to_gene.insert(kmer_i64 as u64, gene_idx);
        }
    }
    info!(n_kmers = kmer_to_gene.len(), "loaded k-mer index");

    Ok(FusionIndex {
        kmer_to_gene,
        gene_names,
        gene_chroms,
        kmer_to_copies,
    })
}

// ─── Read → gene assignment ───────────────────────────────────────────────────

struct ReadHit {
    gene_idx: u32,
    gene1_kmer_count: u32,
    gene1_kmer_min: u32,
    gene1_kmer_max: u32,
    // Second-best gene on this read (if any k-mer hit a different gene).
    // Used by the junction-coherence filter to detect spanning reads.
    gene2_idx: Option<u32>,
    gene2_kmer_count: u32,
    gene2_kmer_min: u32,
    gene2_kmer_max: u32,
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
    max_copies: Option<u32>,
) -> Option<ReadHit> {
    // gene_idx → (count, min_pos, max_pos)
    let mut gene_hits: HashMap<u32, (u32, usize, usize)> = HashMap::new();
    let mut kmer_hits: Vec<(u64, u32, usize)> = Vec::new();
    for (kmer, pos_in_read) in kmer_iter(seq, k) {
        if let Some(&gene_idx) = index.kmer_to_gene.get(&kmer) {
            // Optional genome-wide copy filter: skip k-mers that occur too often
            // (or whose copy number is unknown — NULL in the index).
            if let Some(max) = max_copies {
                match index.kmer_to_copies.get(&kmer) {
                    Some(&c) if c as u32 <= max => {}
                    _ => continue,
                }
            }
            let e = gene_hits.entry(gene_idx).or_insert((0, pos_in_read, pos_in_read));
            e.0 += 1;
            e.1 = e.1.min(pos_in_read);
            e.2 = e.2.max(pos_in_read);
            // Retaining every matching k-mer for every assigned read costs
            // gigabytes on deep BAMs; only do it when the detail TSV is requested.
            if collect_detail {
                kmer_hits.push((kmer, gene_idx, pos_in_read));
            }
        }
    }

    // Sort by count desc, then gene_idx asc so the result is deterministic
    // (HashMap iteration order must not influence the winner choice).
    let mut sorted: Vec<(u32, u32, usize, usize)> = gene_hits
        .into_iter()
        .map(|(gi, (cnt, mn, mx))| (gi, cnt, mn, mx))
        .collect();
    sorted.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

    let (winning_gene, count, g1_min, g1_max) = *sorted.first()?;
    if count < min_hits {
        return None;
    }

    // Expose the second gene whenever any k-mer matched it. Quality gating is
    // handled separately: coherence requires >= min_anchor_kmers from each side
    // (enforced inside read_coherence), and fusion pair selection via
    // fragment_top_pair still relies on gene_idx (the primary/winning gene) which
    // is gated on min_hits. Setting gene2_idx at the 1-k-mer level means
    // n_spanning_reads correctly reflects all reads that cross the junction,
    // even when the anchor on one side is short (asymmetric junction coverage).
    let (gene2_idx, gene2_kmer_count, g2_min, g2_max) = if sorted.len() >= 2 {
        let (g2, c2, mn2, mx2) = sorted[1];
        (Some(g2), c2, mn2, mx2)
    } else {
        (None, 0, 0, 0)
    };

    Some(ReadHit {
        gene_idx: winning_gene,
        gene1_kmer_count: count,
        gene1_kmer_min: g1_min as u32,
        gene1_kmer_max: g1_max as u32,
        gene2_idx,
        gene2_kmer_count,
        gene2_kmer_min: g2_min as u32,
        gene2_kmer_max: g2_max as u32,
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
                rh.gene1_kmer_count, rh.gene1_kmer_min, rh.gene1_kmer_max,
                rh.gene2_kmer_count, rh.gene2_kmer_min, rh.gene2_kmer_max,
            ),
            _ => return (false, false),
        }
    } else if rh.gene_idx == gb {
        match rh.gene2_idx {
            Some(g2) if g2 == ga => (
                rh.gene2_kmer_count, rh.gene2_kmer_min, rh.gene2_kmer_max,
                rh.gene1_kmer_count, rh.gene1_kmer_min, rh.gene1_kmer_max,
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

/// Per-fusion candidate statistics accumulated during fragment aggregation.
struct FusionStats {
    count: u32,
    min_mapq: u8,
    n_spanning_reads: u32,
    /// Number of fragments (qnames) where at least one read showed a coherent
    /// A-block→B-block k-mer partition. Incremented once per fragment, not per read.
    n_coherent_fragments: u32,
}

/// Canonical `GENE_A::GENE_B` label for a normalized `(min_index, max_index)` pair.
/// Uses gene-index order so every output — Parquet, TSV, the `FX` BAM tag, and the
/// kmer-hits TSV — labels the same fusion identically. (The Parquet `gene_a`/`gene_b`
/// columns are filled from `gene_names[key.0]`/`[key.1]`, matching this order.)
fn fusion_pair_label(key: (u32, u32), gene_names: &[String]) -> String {
    format!("{}::{}", gene_names[key.0 as usize], gene_names[key.1 as usize])
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

struct BreakpointAccumulator {
    gene_a_idx: u32,
    gene_b_idx: u32,
    gene_a_chrom_votes: HashMap<String, u32>,
    gene_b_chrom_votes: HashMap<String, u32>,
    spanning_reads: Vec<SpanningReadData>,
}

fn median_i64(values: &mut Vec<i64>) -> Option<f64> {
    if values.is_empty() {
        return None;
    }
    values.sort_unstable();
    let n = values.len();
    if n % 2 == 0 {
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
    let variance = values.iter().map(|&v| (v as f64 - mean).powi(2)).sum::<f64>()
        / (values.len() - 1) as f64;
    Some(variance.sqrt())
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
    writeln!(w, "sample_id\tgene_a\tgene_b\tchrom_a\tchrom_b\tsupporting_reads\tmin_mapq\tn_spanning_reads\tn_coherent_fragments\tn_pon_samples\tpon_total_samples\tmax_pon_supporting_reads")?;
    for r in records {
        let max_pon = r.max_pon_supporting_reads.map(|v| v.to_string()).unwrap_or_else(|| "NA".to_string());
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.sample_id, r.gene_a, r.gene_b, r.chrom_a, r.chrom_b, r.supporting_reads, r.min_mapq,
            r.n_spanning_reads, r.n_coherent_fragments,
            r.n_pon_samples, r.pon_total_samples, max_pon
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
                records.iter().map(|r| r.sample_id.as_str()).collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(StringArray::from(
                records.iter().map(|r| r.gene_a.as_str()).collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(StringArray::from(
                records.iter().map(|r| r.gene_b.as_str()).collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(StringArray::from(
                records.iter().map(|r| r.chrom_a.as_str()).collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(StringArray::from(
                records.iter().map(|r| r.chrom_b.as_str()).collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(Int32Array::from(
                records.iter().map(|r| r.supporting_reads).collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(UInt8Array::from(
                records.iter().map(|r| r.min_mapq).collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(Int32Array::from(
                records.iter().map(|r| r.n_spanning_reads).collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(Int32Array::from(
                records.iter().map(|r| r.n_coherent_fragments).collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(Int32Array::from(
                records.iter().map(|r| r.n_pon_samples).collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(Int32Array::from(
                records.iter().map(|r| r.pon_total_samples).collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(Int32Array::from(
                records.iter().map(|r| r.max_pon_supporting_reads).collect::<Vec<_>>(),
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
fn annotate_fusion_pon(records: &mut [FusionRecord], pon_db: &Path) -> Result<()> {
    if !pon_db.exists() {
        anyhow::bail!("fusion PoN DuckDB not found: {}", pon_db.display());
    }
    let conn = Connection::open(pon_db)
        .with_context(|| format!("failed to open fusion PoN DuckDB: {}", pon_db.display()))?;

    conn.execute_batch("SELECT 1 FROM fusions LIMIT 0").context(
        "fusion PoN DuckDB does not contain a fusions table — build it with \
         `geac merge` over normal-sample *.fusions.parquet files",
    )?;

    let pon_total_samples: i32 =
        conn.query_row("SELECT COUNT(DISTINCT sample_id) FROM fusions", [], |r| r.get(0))?;

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

    info!(pon_total_samples, n_pon_pairs = pon_agg.len(), "fusion PoN annotation complete");
    Ok(())
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
    let index = load_index(&args.index, args.max_kmer_copies)?;

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
            found.context(
                "--sample-id not provided and no SM tag found in BAM/CRAM @RG header",
            )?
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

    for result in reader.records() {
        let record = result.context("error reading BAM record")?;
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
        if let Some(mut rh) = assign_gene(&seq, &index, k, args.min_kmer_hits, false, args.max_kmer_copies) {
            rh.mapq = record.mapq();
            qname_to_reads
                .entry(record.qname().to_vec())
                .or_default()
                .push(rh);
            reads_assigned += 1;
        }

        if reads_processed % 5_000_000 == 0 {
            info!(reads_processed, reads_assigned, "BAM scan progress");
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
            if spanning { n_spanning += 1; }
            if coherent { fragment_coherent = true; }
        }

        let entry = fusion_counts.entry(key).or_insert(FusionStats {
            count: 0,
            min_mapq: 255,
            n_spanning_reads: 0,
            n_coherent_fragments: 0,
        });
        entry.count += 1;
        entry.min_mapq = entry.min_mapq.min(min_mq);
        entry.n_spanning_reads += n_spanning;
        if fragment_coherent { entry.n_coherent_fragments += 1; }
    }

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
        })
        .collect();

    records.sort_by(|a, b| b.supporting_reads.cmp(&a.supporting_reads));

    // Panel-of-Normals annotation and optional filtering, before any output is written.
    if let Some(ref pon_db) = args.fusion_pon {
        annotate_fusion_pon(&mut records, pon_db)?;
    }
    if let Some(max_pon) = args.max_pon_samples {
        let before = records.len();
        records.retain(|r| r.n_pon_samples as u32 <= max_pon);
        let dropped = before - records.len();
        info!(dropped, max_pon_samples = max_pon, "filtered fusions present in PoN");
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

    info!(n_fusions = records.len(), "fusion candidates identified");
    write_fusion_parquet(&records, &args.output)?;
    info!(output = %args.output.display(), "fusion candidates written");

    if let Some(ref tsv_path) = args.tsv_output {
        write_fusion_tsv(&records, tsv_path)?;
        info!(output = %tsv_path.display(), "fusion TSV written");
    }

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
    let fusion_qnames: HashMap<Vec<u8>, String> = qname_to_pair
        .into_iter()
        .filter_map(|(qname, key)| {
            let label = fusion_label.get(&key)?;
            Some((qname, label.clone()))
        })
        .collect();

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
        let mut writer =
            bam::Writer::from_path(reads_output, &header, bam::Format::Bam)
                .with_context(|| format!("failed to create output BAM: {}", reads_output.display()))?;

        let mut reads_written: u64 = 0;
        for result in reader2.records() {
            let mut record = result.context("error reading BAM record in second pass")?;
            if let Some(label) = fusion_qnames.get(record.qname()) {
                // Tag the read with the fusion it supports, e.g. FX:Z:EWSR1::FLI1.
                // Remove any pre-existing FX so re-runs don't duplicate the tag.
                let _ = record.remove_aux(b"FX");
                record
                    .push_aux(b"FX", bam::record::Aux::String(label))
                    .context("failed to add FX tag to BAM record")?;
                writer.write(&record).context("failed to write BAM record")?;
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

    // Optional combined second BAM pass: produces kmer-hits TSV and/or breakpoints TSV.
    // Both outputs share one BAM scan to avoid re-reading the file twice.
    if args.kmer_hits_output.is_some() || args.breakpoints_output.is_some() {
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
                .map(|i| std::str::from_utf8(hdr.tid2name(i)).unwrap_or("?").to_string())
                .collect()
        };

        // kmer-hits TSV writer (optional).
        let mut kmer_hits_writer: Option<BufWriter<std::fs::File>> =
            if let Some(ref path) = args.kmer_hits_output {
                let file = std::fs::File::create(path)
                    .with_context(|| format!("failed to create kmer hits TSV: {}", path.display()))?;
                let mut w = BufWriter::new(file);
                writeln!(w, "fusion\tsample_id\tread_name\tread_end\tchrom\tpos\tgene_matched\tkmer_pos_in_read\tkmer_hash\tkmer_seq")?;
                Some(w)
            } else {
                None
            };

        // Breakpoint accumulators (optional): one per passing fusion.
        let mut accumulators: Option<HashMap<String, BreakpointAccumulator>> =
            if args.breakpoints_output.is_some() {
                let mut map: HashMap<String, BreakpointAccumulator> = HashMap::new();
                for (&(ga, gb), label) in &fusion_label {
                    map.insert(label.clone(), BreakpointAccumulator {
                        gene_a_idx: ga,
                        gene_b_idx: gb,
                        gene_a_chrom_votes: HashMap::new(),
                        gene_b_chrom_votes: HashMap::new(),
                        spanning_reads: Vec::new(),
                    });
                }
                Some(map)
            } else {
                None
            };

        let mut rows_written: u64 = 0;
        for result in reader2.records() {
            let record = result.context("error reading BAM record in detail pass")?;
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
            let Some(rh) = assign_gene(&seq, &index, k, args.min_kmer_hits, true, args.max_kmer_copies) else {
                continue;
            };
            let read_end = if flags & 0x40 != 0 { "R1" } else { "R2" };
            let tid = record.tid();
            let chrom_str: &str = if tid >= 0 {
                target_names.get(tid as usize).map(|s| s.as_str()).unwrap_or("*")
            } else {
                "*"
            };
            let pos1 = record.pos() + 1; // 1-based alignment start

            // kmer-hits TSV rows.
            if let Some(ref mut w) = kmer_hits_writer {
                let read_name = std::str::from_utf8(record.qname()).unwrap_or("?");
                for &(kmer_hash, gene_idx, kmer_pos_in_read) in &rh.kmer_hits {
                    let gene_name = index.gene_names.get(gene_idx as usize).map(|s| s.as_str()).unwrap_or("?");
                    let kmer_seq = decode_kmer(kmer_hash, k);
                    writeln!(
                        w,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        fusion_label_str, sample_id, read_name, read_end, chrom_str, pos1,
                        gene_name, kmer_pos_in_read, kmer_hash as i64, kmer_seq
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
                        *acc.gene_a_chrom_votes.entry(chrom_str.to_string()).or_insert(0) += 1;
                    } else if a_positions.is_empty() && !b_positions.is_empty() {
                        *acc.gene_b_chrom_votes.entry(chrom_str.to_string()).or_insert(0) += 1;
                    } else if !a_positions.is_empty() && !b_positions.is_empty() {
                        // Spanning read.
                        let last_a = *a_positions.iter().max().unwrap();
                        let first_a = *a_positions.iter().min().unwrap();
                        let last_b = *b_positions.iter().max().unwrap();
                        let first_b = *b_positions.iter().min().unwrap();
                        acc.spanning_reads.push(SpanningReadData {
                            chrom: chrom_str.to_string(),
                            pos: pos1,
                            last_a, first_a, last_b, first_b,
                        });
                    }
                }
            }
        }

        if let Some(ref path) = args.kmer_hits_output {
            info!(rows_written, output = %path.display(), "kmer hits TSV written");
        }

        // Compute and write breakpoints TSV.
        if let Some(ref bp_path) = args.breakpoints_output {
            let file = std::fs::File::create(bp_path)
                .with_context(|| format!("failed to create breakpoints TSV: {}", bp_path.display()))?;
            let mut w = BufWriter::new(file);
            writeln!(w, "fusion\tgene_a\tchrom_a\tbreakpoint_a\tbp_a_n\tbp_a_std\tgene_b\tchrom_b\tbreakpoint_b\tbp_b_n\tbp_b_std\tn_spanning_reads")?;

            let acc_map = accumulators.as_mut().unwrap();
            let mut fusions_sorted: Vec<String> = acc_map.keys().cloned().collect();
            fusions_sorted.sort_unstable();

            for label in &fusions_sorted {
                let acc = acc_map.get_mut(label.as_str()).unwrap();
                let gene_a = &index.gene_names[acc.gene_a_idx as usize];
                let gene_b = &index.gene_names[acc.gene_b_idx as usize];

                let modal_chrom = |votes: &HashMap<String, u32>| -> String {
                    votes.iter()
                        .max_by_key(|(_, &v)| v)
                        .map(|(c, _)| c.clone())
                        .unwrap_or_else(|| "*".to_string())
                };
                let chrom_a = modal_chrom(&acc.gene_a_chrom_votes);
                let chrom_b = modal_chrom(&acc.gene_b_chrom_votes);

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

                let n_spanning = acc.spanning_reads.len();
                let bp_a = median_i64(&mut bp_a_estimates);
                let bp_a_std = std_dev_i64(&bp_a_estimates);
                let bp_a_n = bp_a_estimates.len();
                let bp_b = median_i64(&mut bp_b_estimates);
                let bp_b_std = std_dev_i64(&bp_b_estimates);
                let bp_b_n = bp_b_estimates.len();

                let fmt_opt_f64 = |v: Option<f64>| v.map(|x| format!("{:.1}", x)).unwrap_or_else(|| "NA".to_string());
                let fmt_opt_i64 = |v: Option<f64>| v.map(|x| format!("{}", x as i64)).unwrap_or_else(|| "NA".to_string());

                writeln!(
                    w,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    label, gene_a, chrom_a,
                    fmt_opt_i64(bp_a), bp_a_n, fmt_opt_f64(bp_a_std),
                    gene_b, chrom_b,
                    fmt_opt_i64(bp_b), bp_b_n, fmt_opt_f64(bp_b_std),
                    n_spanning,
                )?;
            }
            info!(output = %bp_path.display(), "breakpoints TSV written");
        }
    }

    Ok(())
}

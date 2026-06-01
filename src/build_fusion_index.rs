use std::collections::{BTreeSet, HashMap, HashSet};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use duckdb::Connection;
use flate2::read::GzDecoder;
use rust_htslib::faidx;
use tracing::info;

use crate::cli::BuildFusionIndexArgs;
use crate::kmer::kmer_iter;

// ─── Gene body ────────────────────────────────────────────────────────────────

struct GeneBody {
    gene_name: String,
    chrom: String,
    start: usize, // 0-based inclusive
    end: usize,   // 0-based exclusive (half-open)
}

fn parse_gene_bodies(gtf_path: &Path) -> Result<Vec<GeneBody>> {
    let file = std::fs::File::open(gtf_path)
        .with_context(|| format!("cannot open GTF: {}", gtf_path.display()))?;
    let reader: Box<dyn BufRead> =
        if gtf_path.extension().and_then(|e| e.to_str()) == Some("gz") {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

    let mut genes: Vec<GeneBody> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.splitn(9, '\t').collect();
        if fields.len() < 9 || fields[2] != "gene" {
            continue;
        }
        // GTF is 1-based inclusive → 0-based half-open
        let start: usize = match fields[3].parse::<usize>() {
            Ok(v) => v.saturating_sub(1),
            Err(_) => continue,
        };
        let end: usize = match fields[4].parse::<usize>() {
            Ok(v) => v,
            Err(_) => continue,
        };
        if end <= start {
            continue;
        }
        let attrs = fields[8];
        let gene_name = match extract_gtf_gene_name(attrs) {
            Some(g) => g,
            None => continue,
        };
        genes.push(GeneBody {
            gene_name,
            chrom: fields[0].to_string(),
            start,
            end,
        });
    }
    Ok(genes)
}

fn extract_gtf_gene_name(attrs: &str) -> Option<String> {
    for key in &["gene_name", "gene_id"] {
        let needle = format!("{} \"", key);
        if let Some(pos) = attrs.find(&needle) {
            let rest = &attrs[pos + needle.len()..];
            if let Some(end) = rest.find('"') {
                return Some(rest[..end].to_string());
            }
        }
    }
    None
}

// ─── Gene list filter ─────────────────────────────────────────────────────────

fn load_gene_list(path: &Path) -> Result<HashSet<String>> {
    let file = std::fs::File::open(path)
        .with_context(|| format!("cannot open gene list: {}", path.display()))?;
    let genes = BufReader::new(file)
        .lines()
        .map_while(Result::ok)
        .map(|l| l.trim().to_string())
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .collect::<HashSet<String>>();
    Ok(genes)
}

// ─── FASTA index parsing ──────────────────────────────────────────────────────

fn read_fai_sequences(fasta: &Path) -> Result<Vec<(String, usize)>> {
    let mut fai_str = fasta.to_string_lossy().into_owned();
    fai_str.push_str(".fai");
    let fai_path = std::path::PathBuf::from(&fai_str);

    let file = std::fs::File::open(&fai_path).with_context(|| {
        format!(
            "cannot open FASTA index {}; run 'samtools faidx' first",
            fai_path.display()
        )
    })?;

    let mut seqs = Vec::new();
    for line in BufReader::new(file).lines() {
        let line = line?;
        let mut fields = line.split('\t');
        let name = fields.next().unwrap_or("").to_string();
        let len: usize = fields.next().and_then(|s| s.parse().ok()).unwrap_or(0);
        if !name.is_empty() && len > 0 {
            seqs.push((name, len));
        }
    }
    Ok(seqs)
}

// ─── DuckDB output ────────────────────────────────────────────────────────────

fn write_index(
    kmer_to_gene: &HashMap<u64, u32>,
    kmer_to_pos: &HashMap<u64, (String, u32)>,
    genes: &[GeneBody],
    genome_counts: Option<&HashMap<u64, u8>>,
    k: usize,
    output: &Path,
) -> Result<()> {
    if output.exists() {
        std::fs::remove_file(output)
            .with_context(|| format!("failed to remove existing output: {}", output.display()))?;
    }

    let conn = Connection::open(output)
        .with_context(|| format!("failed to open DuckDB: {}", output.display()))?;

    // genome_copies is the genome-wide occurrence count of each k-mer, populated
    // only when --check-genome-uniqueness ran; NULL otherwise (count unknown).
    conn.execute_batch(
        "CREATE TABLE meta (key VARCHAR, value VARCHAR);
         CREATE TABLE genes (gene_index UINTEGER, gene_name VARCHAR, chrom VARCHAR);
         CREATE TABLE kmers (kmer_hash BIGINT, gene_index UINTEGER, genome_copies INTEGER);
         CREATE TABLE kmer_positions (kmer_hash BIGINT, chrom VARCHAR, pos INTEGER);",
    )?;

    // Record the k-mer size so `geac fusions` can verify the --kmer-size it was
    // given matches the value the index was built with. A mismatch otherwise
    // produces silently empty/garbage results (every hash fails to match).
    conn.execute(
        "INSERT INTO meta (key, value) VALUES ('kmer_size', ?)",
        duckdb::params![k.to_string()],
    )?;

    {
        let mut app = conn
            .appender("genes")
            .context("failed to create genes appender")?;
        for (i, gene) in genes.iter().enumerate() {
            app.append_row(duckdb::params![i as u32, &gene.gene_name, &gene.chrom])?;
        }
        app.flush()?;
    }

    {
        let mut app = conn
            .appender("kmers")
            .context("failed to create kmers appender")?;
        for (&kmer_hash, &gene_idx) in kmer_to_gene {
            let copies: Option<i32> = genome_counts
                .and_then(|gc| gc.get(&kmer_hash))
                .map(|&c| c as i32);
            app.append_row(duckdb::params![kmer_hash as i64, gene_idx, copies])?;
        }
        app.flush()?;
    }

    {
        let mut app = conn
            .appender("kmer_positions")
            .context("failed to create kmer_positions appender")?;
        for (&kmer_hash, (chrom, pos)) in kmer_to_pos {
            app.append_row(duckdb::params![kmer_hash as i64, chrom, *pos as i32])?;
        }
        app.flush()?;
    }

    conn.execute_batch(
        "CREATE INDEX kmers_hash_idx ON kmers(kmer_hash);
         CREATE INDEX kmer_positions_chrom_pos_idx ON kmer_positions(chrom, pos);",
    )?;

    info!(
        n_genes = genes.len(),
        n_kmers = kmer_to_gene.len(),
        "index written"
    );
    Ok(())
}

// ─── BED output ───────────────────────────────────────────────────────────────

/// Merge per-chromosome k-mer start positions into intervals [pos, pos+k) and
/// write them as BED rows, in `chrom_order` (FASTA index order) so the file is
/// already sorted for bedtools / IGV. Returns the number of intervals written.
fn write_merged_bed_intervals<W: std::io::Write>(
    w: &mut W,
    pos_by_chrom: &HashMap<String, Vec<u32>>,
    chrom_order: &[(String, usize)],
    k: usize,
) -> Result<usize> {
    let mut n_intervals: usize = 0;
    for (chrom_name, _) in chrom_order {
        let Some(positions) = pos_by_chrom.get(chrom_name.as_str()) else {
            continue;
        };

        let mut sorted = positions.clone();
        sorted.sort_unstable();

        // Merge: extend the current interval while the next k-mer starts before
        // the current interval ends (i.e. they are adjacent or overlapping).
        let mut iter = sorted.into_iter();
        let Some(first) = iter.next() else { continue };
        let mut interval_start = first;
        let mut interval_end = first + k as u32;

        for pos in iter {
            if pos <= interval_end {
                // Extend: this k-mer overlaps or is immediately adjacent.
                interval_end = interval_end.max(pos + k as u32);
            } else {
                writeln!(w, "{}\t{}\t{}", chrom_name, interval_start, interval_end)?;
                n_intervals += 1;
                interval_start = pos;
                interval_end = pos + k as u32;
            }
        }
        // Flush the last interval.
        writeln!(w, "{}\t{}\t{}", chrom_name, interval_start, interval_end)?;
        n_intervals += 1;
    }
    Ok(n_intervals)
}

/// Write a BED file of merged intervals covering all retained k-mer start positions.
fn write_kmer_bed(
    kmer_to_pos: &HashMap<u64, (String, u32)>,
    chrom_order: &[(String, usize)],
    k: usize,
    output: &Path,
) -> Result<()> {
    use std::io::BufWriter;

    // Group positions by chromosome.
    let mut pos_by_chrom: HashMap<String, Vec<u32>> = HashMap::new();
    for (chrom, pos) in kmer_to_pos.values() {
        pos_by_chrom.entry(chrom.clone()).or_default().push(*pos);
    }

    let file = std::fs::File::create(output)
        .with_context(|| format!("failed to create BED: {}", output.display()))?;
    let mut w = BufWriter::new(file);

    let n_intervals = write_merged_bed_intervals(&mut w, &pos_by_chrom, chrom_order, k)?;

    info!(
        n_intervals,
        output = %output.display(),
        "unique k-mer BED written"
    );
    Ok(())
}

/// Per-chromosome gene intervals over the *full* GTF (not just the panel), used to
/// annotate where a k-mer's genome occurrences fall. Intervals are sorted by start;
/// `max_len` bounds the backward scan during lookup.
struct ChromIntervals {
    intervals: Vec<(usize, usize, u32)>, // (start, end exclusive, name_idx), sorted by start
    max_len: usize,
}

struct GeneIntervals {
    by_chrom: HashMap<String, ChromIntervals>,
    names: Vec<String>,
}

impl GeneIntervals {
    fn build(genes: &[GeneBody]) -> Self {
        let mut name_to_idx: HashMap<&str, u32> = HashMap::new();
        let mut names: Vec<String> = Vec::new();
        let mut by_chrom: HashMap<String, Vec<(usize, usize, u32)>> = HashMap::new();
        for gene in genes {
            let idx = *name_to_idx.entry(gene.gene_name.as_str()).or_insert_with(|| {
                names.push(gene.gene_name.clone());
                (names.len() - 1) as u32
            });
            by_chrom
                .entry(gene.chrom.clone())
                .or_default()
                .push((gene.start, gene.end, idx));
        }
        let by_chrom = by_chrom
            .into_iter()
            .map(|(chrom, mut ivs)| {
                ivs.sort_unstable_by_key(|&(s, _, _)| s);
                let max_len = ivs.iter().map(|&(s, e, _)| e.saturating_sub(s)).max().unwrap_or(0);
                (chrom, ChromIntervals { intervals: ivs, max_len })
            })
            .collect();
        GeneIntervals { by_chrom, names }
    }

    /// Gene name indices of every gene whose body covers `pos` on `chrom`.
    fn genes_at(&self, chrom: &str, pos: usize) -> Vec<u32> {
        let Some(ci) = self.by_chrom.get(chrom) else {
            return Vec::new();
        };
        // All intervals with start <= pos; scan backward bounded by max_len.
        let upper = ci.intervals.partition_point(|&(s, _, _)| s <= pos);
        let lower = pos.saturating_sub(ci.max_len);
        let mut out = Vec::new();
        for &(s, e, idx) in ci.intervals[..upper].iter().rev() {
            if s < lower {
                break;
            }
            if pos < e {
                out.push(idx);
            }
        }
        out
    }
}

/// Per-k-mer annotation: the set of full-GTF genes its genome occurrences fall in,
/// plus whether any occurrence is intergenic (outside every gene).
#[derive(Default)]
struct GeneAnnot {
    genes: BTreeSet<u32>,
    intergenic: bool,
}

/// Format the BED column-4 label for a merged interval: a comma-separated, sorted
/// list of gene names (plus "intergenic" if any occurrence was outside a gene), or
/// just "intergenic" when no occurrence fell in a gene.
fn format_gene_label(genes: &BTreeSet<u32>, intergenic: bool, names: &[String]) -> String {
    let mut parts: Vec<&str> = genes.iter().map(|&i| names[i as usize].as_str()).collect();
    parts.sort_unstable();
    if intergenic {
        parts.push("intergenic");
    }
    if parts.is_empty() {
        "intergenic".to_string()
    } else {
        parts.join(",")
    }
}

/// Write one BED per genome-wide copy tier: `<prefix>.copies1.bed`,
/// `<prefix>.copies2.bed`, … up to `max_copies` (the retained range). Each file
/// holds intervals merged only among k-mers sharing that exact copy number, so the
/// strictly-unique (1×) and near-unique (2×, …) regions can be loaded and colored
/// separately in IGV, or used directly as target BEDs. Empty tiers are skipped.
///
/// For tiers ≥ 2, when `kmer_annot` is supplied a 4th column lists every full-GTF
/// gene the interval's k-mers occur in (or "intergenic"). The 1× tier is left as
/// plain 3-column BED — its only "match" is the k-mer's own gene.
fn write_kmer_bed_by_copies(
    kmer_to_pos: &HashMap<u64, (String, u32)>,
    genome_counts: &HashMap<u64, u8>,
    chrom_order: &[(String, usize)],
    k: usize,
    max_copies: u32,
    prefix: &Path,
    kmer_annot: Option<&HashMap<u64, GeneAnnot>>,
    gene_names: &[String],
) -> Result<()> {
    use std::io::{BufWriter, Write};

    // Group retained k-mers by (copy number → chromosome → [(pos, kmer_hash)]).
    // The k-mer hash is kept so the gene annotation can be looked up during merge.
    let mut by_copy: HashMap<u32, HashMap<String, Vec<(u32, u64)>>> = HashMap::new();
    for (kmer, (chrom, pos)) in kmer_to_pos {
        let copies = genome_counts.get(kmer).copied().unwrap_or(0) as u32;
        if copies >= 1 {
            by_copy
                .entry(copies)
                .or_default()
                .entry(chrom.clone())
                .or_default()
                .push((*pos, *kmer));
        }
    }

    for copies in 1..=max_copies {
        let Some(pos_by_chrom) = by_copy.get(&copies) else {
            continue;
        };
        let path = PathBuf::from(format!("{}.copies{}.bed", prefix.display(), copies));
        let file = std::fs::File::create(&path)
            .with_context(|| format!("failed to create BED: {}", path.display()))?;
        let mut w = BufWriter::new(file);

        // Annotate the gene column only for the multi-copy tiers.
        let annot = if copies >= 2 { kmer_annot } else { None };

        let mut n_intervals: usize = 0;
        for (chrom_name, _) in chrom_order {
            let Some(entries) = pos_by_chrom.get(chrom_name.as_str()) else {
                continue;
            };
            let mut sorted = entries.clone();
            sorted.sort_unstable_by_key(|&(p, _)| p);

            let mut iter = sorted.into_iter();
            let Some((first_pos, first_kmer)) = iter.next() else {
                continue;
            };
            let mut interval_start = first_pos;
            let mut interval_end = first_pos + k as u32;
            let mut cur_genes: BTreeSet<u32> = BTreeSet::new();
            let mut cur_intergenic = false;
            let accumulate = |kmer: u64, genes: &mut BTreeSet<u32>, inter: &mut bool| {
                if let Some(a) = annot.and_then(|m| m.get(&kmer)) {
                    genes.extend(a.genes.iter().copied());
                    *inter |= a.intergenic;
                }
            };
            accumulate(first_kmer, &mut cur_genes, &mut cur_intergenic);

            let flush = |w: &mut BufWriter<std::fs::File>,
                             start: u32,
                             end: u32,
                             genes: &BTreeSet<u32>,
                             inter: bool|
             -> Result<()> {
                if annot.is_some() {
                    let label = format_gene_label(genes, inter, gene_names);
                    writeln!(w, "{}\t{}\t{}\t{}", chrom_name, start, end, label)?;
                } else {
                    writeln!(w, "{}\t{}\t{}", chrom_name, start, end)?;
                }
                Ok(())
            };

            for (pos, kmer) in iter {
                if pos <= interval_end {
                    interval_end = interval_end.max(pos + k as u32);
                    accumulate(kmer, &mut cur_genes, &mut cur_intergenic);
                } else {
                    flush(&mut w, interval_start, interval_end, &cur_genes, cur_intergenic)?;
                    n_intervals += 1;
                    interval_start = pos;
                    interval_end = pos + k as u32;
                    cur_genes.clear();
                    cur_intergenic = false;
                    accumulate(kmer, &mut cur_genes, &mut cur_intergenic);
                }
            }
            flush(&mut w, interval_start, interval_end, &cur_genes, cur_intergenic)?;
            n_intervals += 1;
        }
        info!(
            copies,
            n_intervals,
            output = %path.display(),
            "per-copy-tier k-mer BED written"
        );
    }
    Ok(())
}

/// Merge k-mer start positions into intervals [pos, pos+k) and return total
/// covered bases. Adjacent/overlapping windows collapse, matching `write_kmer_bed`.
fn covered_bases_from_positions(positions: &[u32], k: usize) -> usize {
    if positions.is_empty() {
        return 0;
    }
    let mut sorted = positions.to_vec();
    sorted.sort_unstable();
    let mut total = 0usize;
    let mut iter = sorted.into_iter();
    let first = iter.next().unwrap();
    let mut cur_start = first;
    let mut cur_end = first + k as u32;
    for pos in iter {
        if pos <= cur_end {
            cur_end = cur_end.max(pos + k as u32);
        } else {
            total += (cur_end - cur_start) as usize;
            cur_start = pos;
            cur_end = pos + k as u32;
        }
    }
    total += (cur_end - cur_start) as usize;
    total
}

/// Write a per-gene uniqueness TSV.
///
/// Always reports the panel-wide tier (uniqueness among the indexed gene bodies,
/// i.e. what survives cross-gene dedup): `panel_unique_kmers`, `pct_unique_windows`
/// (panel-unique k-mers / k-mer windows) and `pct_unique_bases` (gene-body bases
/// covered by a panel-unique k-mer / body_len).
///
/// When `genome_counts` is provided (i.e. --check-genome-uniqueness ran), also
/// reports genome-wide copy-number buckets over the gene's panel-unique k-mers and
/// base coverage at the strict (==1) and relaxed (<= max_genome_copies) tiers.
///
/// Called before the genome-copy retain and the --min-gene-kmers prune so the maps
/// still hold the full panel-unique survivor set and every gene reports real values.
fn write_gene_stats(
    genes: &[GeneBody],
    kmer_to_gene: &HashMap<u64, u32>,
    kmer_to_pos: &HashMap<u64, (String, u32)>,
    k: usize,
    genome_counts: Option<&HashMap<u64, u8>>,
    max_genome_copies: u32,
    output: &Path,
) -> Result<()> {
    use std::io::{BufWriter, Write};

    // Group each gene's panel-unique k-mers as (position, genome_copies). When the
    // genome pass didn't run, copies is recorded as 0 (unknown) and unused.
    let mut by_gene: HashMap<u32, Vec<(u32, u8)>> = HashMap::new();
    for (kmer, &gene_idx) in kmer_to_gene {
        if let Some((_, pos)) = kmer_to_pos.get(kmer) {
            let copies = genome_counts
                .and_then(|gc| gc.get(kmer))
                .copied()
                .unwrap_or(0);
            by_gene.entry(gene_idx).or_default().push((*pos, copies));
        }
    }

    let file = std::fs::File::create(output)
        .with_context(|| format!("failed to create gene-stats TSV: {}", output.display()))?;
    let mut w = BufWriter::new(file);

    let have_genome = genome_counts.is_some();
    if have_genome {
        writeln!(
            w,
            "gene\tchrom\tstart\tend\tbody_len\tpanel_unique_kmers\tkmer_windows\t\
             pct_unique_windows\tcovered_bases\tpct_unique_bases\t\
             n_copies_1\tn_copies_2\tn_copies_3_10\tn_copies_gt10\t\
             genome_unique_pct_bases\trelaxed_pct_bases"
        )?;
    } else {
        writeln!(
            w,
            "gene\tchrom\tstart\tend\tbody_len\tpanel_unique_kmers\tkmer_windows\t\
             pct_unique_windows\tcovered_bases\tpct_unique_bases"
        )?;
    }

    for (gene_idx, gene) in genes.iter().enumerate() {
        let body_len = gene.end.saturating_sub(gene.start);
        // Number of k-mer windows in the body; 0 if the body is shorter than k.
        let kmer_windows = body_len.saturating_sub(k.saturating_sub(1));

        let entries = by_gene.get(&(gene_idx as u32));
        let panel_unique_kmers = entries.map_or(0, |e| e.len());

        let all_positions: Vec<u32> =
            entries.map_or_else(Vec::new, |e| e.iter().map(|(p, _)| *p).collect());
        let covered_bases = covered_bases_from_positions(&all_positions, k);

        let pct = |num: usize, den: usize| -> f64 {
            if den > 0 {
                num as f64 / den as f64 * 100.0
            } else {
                0.0
            }
        };
        let pct_windows = pct(panel_unique_kmers, kmer_windows);
        let pct_bases = pct(covered_bases, body_len);

        write!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{:.2}",
            gene.gene_name,
            gene.chrom,
            gene.start,
            gene.end,
            body_len,
            panel_unique_kmers,
            kmer_windows,
            pct_windows,
            covered_bases,
            pct_bases,
        )?;

        if have_genome {
            let mut n1 = 0usize;
            let mut n2 = 0usize;
            let mut n3_10 = 0usize;
            let mut ngt10 = 0usize;
            let mut pos_eq1: Vec<u32> = Vec::new();
            let mut pos_le_max: Vec<u32> = Vec::new();
            if let Some(e) = entries {
                for &(pos, copies) in e {
                    match copies {
                        1 => n1 += 1,
                        2 => n2 += 1,
                        3..=10 => n3_10 += 1,
                        _ => ngt10 += 1,
                    }
                    if copies == 1 {
                        pos_eq1.push(pos);
                    }
                    if copies as u32 <= max_genome_copies {
                        pos_le_max.push(pos);
                    }
                }
            }
            let genome_unique_pct = pct(covered_bases_from_positions(&pos_eq1, k), body_len);
            let relaxed_pct = pct(covered_bases_from_positions(&pos_le_max, k), body_len);
            write!(
                w,
                "\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}",
                n1, n2, n3_10, ngt10, genome_unique_pct, relaxed_pct
            )?;
        }
        writeln!(w)?;
    }

    info!(
        n_genes = genes.len(),
        output = %output.display(),
        "per-gene uniqueness TSV written"
    );
    Ok(())
}

/// Write a global histogram of genome-wide k-mer copy number over the whole
/// candidate set: two columns `copies` and `n_kmers`, sorted ascending by copies.
/// Copy counts saturate at 255 (the genome pass uses u8), so the last row is a
/// ">=255" bucket in practice.
fn write_copy_histogram(genome_counts: &HashMap<u64, u8>, output: &Path) -> Result<()> {
    use std::io::{BufWriter, Write};

    let mut hist: HashMap<u8, u64> = HashMap::new();
    for &c in genome_counts.values() {
        *hist.entry(c).or_insert(0) += 1;
    }
    let mut rows: Vec<(u8, u64)> = hist.into_iter().collect();
    rows.sort_unstable_by_key(|&(c, _)| c);

    let file = std::fs::File::create(output)
        .with_context(|| format!("failed to create copy-histogram TSV: {}", output.display()))?;
    let mut w = BufWriter::new(file);
    writeln!(w, "copies\tn_kmers")?;
    for (copies, n) in rows {
        writeln!(w, "{}\t{}", copies, n)?;
    }

    info!(output = %output.display(), "copy-number histogram written");
    Ok(())
}

// ─── Main entry point ─────────────────────────────────────────────────────────

const MULTI_GENE: u32 = u32::MAX;

pub fn build_fusion_index(args: &BuildFusionIndexArgs) -> Result<()> {
    let k = args.kmer_size as usize;
    anyhow::ensure!(k > 0 && k <= 31, "--kmer-size must be between 1 and 31");

    // --max-genome-copies and --copy-histogram-output both depend on the
    // genome-wide copy counts, which only exist when the genome pass runs.
    anyhow::ensure!(
        args.max_genome_copies >= 1,
        "--max-genome-copies must be >= 1"
    );
    if !args.check_genome_uniqueness {
        anyhow::ensure!(
            args.max_genome_copies == 1,
            "--max-genome-copies > 1 requires --check-genome-uniqueness (no genome-wide counts otherwise)"
        );
        anyhow::ensure!(
            args.copy_histogram_output.is_none(),
            "--copy-histogram-output requires --check-genome-uniqueness (the copy counts come from that pass)"
        );
        anyhow::ensure!(
            args.bed_output_by_copies.is_none(),
            "--bed-output-by-copies requires --check-genome-uniqueness (the copy counts come from that pass)"
        );
    }

    info!("parsing gene bodies from GTF...");
    let mut genes = parse_gene_bodies(&args.gtf)?;
    info!(n_genes = genes.len(), "gene bodies loaded");

    // Build the gene-annotation index from the FULL GTF (before the panel filter)
    // so per-copy BEDs can name paralogs that fall outside the --genes panel.
    let gene_annot: Option<GeneIntervals> = if args.bed_output_by_copies.is_some() {
        info!("building full-GTF gene annotation index for per-copy BED labels...");
        Some(GeneIntervals::build(&genes))
    } else {
        None
    };

    if let Some(ref genes_file) = args.genes {
        let allowed = load_gene_list(genes_file)?;
        let before = genes.len();
        genes.retain(|g| allowed.contains(g.gene_name.as_str()));
        info!(
            requested = allowed.len(),
            matched = genes.len(),
            skipped = before - genes.len(),
            "filtered to gene list"
        );
        if genes.is_empty() {
            anyhow::bail!(
                "no genes from --genes matched any gene_name in the GTF; \
                 check that names match exactly (e.g. 'BCR' not 'BCR-ABL1')"
            );
        }
    }

    let fai = faidx::Reader::from_path(&args.fasta)
        .with_context(|| format!("failed to open FASTA: {}", args.fasta.display()))?;

    // Step 1 — extract gene-body k-mers and track which gene each belongs to.
    // k-mers found in 2+ gene bodies are marked MULTI_GENE and discarded.
    //
    // We group genes by chromosome and fetch each chromosome sequence once,
    // then slice gene bodies from the in-memory sequence. This avoids 78K
    // individual faidx seeks (one per gene) and replaces them with ~25
    // sequential chromosome reads.
    info!("extracting gene-body k-mers (one chromosome fetch per sequence)...");

    // Build a per-chromosome index: chrom → [(start, end, gene_idx)]
    let mut genes_by_chrom: HashMap<String, Vec<(usize, usize, u32)>> = HashMap::new();
    for (gene_idx, gene) in genes.iter().enumerate() {
        genes_by_chrom
            .entry(gene.chrom.clone())
            .or_default()
            .push((gene.start, gene.end, gene_idx as u32));
    }

    // Process chromosomes in the order they appear in the FASTA index so that
    // progress is predictable and the standard chroms come first.
    let seq_list_for_genes = read_fai_sequences(&args.fasta)?;

    // Global maps accumulate survivors from all chromosomes.
    let mut kmer_to_gene: HashMap<u64, u32> = HashMap::new();
    // kmer_to_pos stores the 0-based chromosome position of the first occurrence
    // of each k-mer found during gene-body extraction.
    let mut kmer_to_pos: HashMap<u64, (String, u32)> = HashMap::new();
    let mut chroms_done = 0usize;

    for (chrom_name, chrom_len) in &seq_list_for_genes {
        let Some(chrom_genes) = genes_by_chrom.get(chrom_name.as_str()) else {
            continue;
        };

        info!(
            chrom = %chrom_name,
            n_genes = chrom_genes.len(),
            "fetching chromosome for gene k-mer extraction..."
        );

        let chrom_seq = match fai.fetch_seq(chrom_name, 0, chrom_len.saturating_sub(1)) {
            Ok(s) => s.to_vec(),
            Err(e) => {
                tracing::warn!(chrom = %chrom_name, "failed to fetch chromosome: {e}");
                continue;
            }
        };

        // Local map: scoped to this chromosome only. Stays small and cache-hot
        // regardless of how many chromosomes have already been processed.
        let mut local: HashMap<u64, u32> = HashMap::new();
        let mut local_pos: HashMap<u64, u32> = HashMap::new();
        for &(start, end, gene_idx) in chrom_genes {
            let end_clamped = end.min(chrom_seq.len());
            if start >= end_clamped {
                continue;
            }
            for (kmer, offset) in kmer_iter(&chrom_seq[start..end_clamped], k) {
                local
                    .entry(kmer)
                    .and_modify(|v| {
                        if *v != gene_idx {
                            *v = MULTI_GENE;
                        }
                    })
                    .or_insert(gene_idx);
                // Record the chromosome position of the first occurrence only.
                // offset is the true 0-based start position within the slice,
                // so (start + offset) gives the correct chromosome coordinate.
                local_pos.entry(kmer).or_insert((start + offset) as u32);
            }
        }

        // Prune within-chromosome cross-gene ambiguity, then merge survivors
        // into the global maps. Cross-chromosome duplicates are marked there.
        local.retain(|_, v| *v != MULTI_GENE);
        local_pos.retain(|kmer, _| local.contains_key(kmer));
        for (kmer, gene_idx) in local {
            let is_new = !kmer_to_gene.contains_key(&kmer);
            kmer_to_gene
                .entry(kmer)
                .and_modify(|v| {
                    *v = MULTI_GENE; // k-mer already seen on a different chromosome
                })
                .or_insert(gene_idx);
            if is_new {
                if let Some(&pos) = local_pos.get(&kmer) {
                    kmer_to_pos.insert(kmer, (chrom_name.clone(), pos));
                }
            }
        }

        chroms_done += 1;
        info!(
            chrom = %chrom_name,
            chroms_done,
            n_global = kmer_to_gene.len(),
            "chromosome done"
        );
    }

    // Single final prune for cross-chromosome duplicates.
    kmer_to_gene.retain(|_, v| *v != MULTI_GENE);
    kmer_to_pos.retain(|kmer, _| kmer_to_gene.contains_key(kmer));
    info!(
        n_candidate_kmers = kmer_to_gene.len(),
        "k-mers after cross-gene dedup"
    );

    // Step 2 — optional genome-wide uniqueness pass.
    // Scans the full FASTA and removes k-mers that appear more than once in the
    // genome (e.g. in repetitive intergenic regions). At k ≥ 19 the vast majority
    // of candidate k-mers that survived cross-gene dedup are already genome-unique,
    // so this pass is skipped by default. Enable with --check-genome-uniqueness when
    // you need the stricter guarantee and have sufficient RAM (the candidate set must
    // fit in a HashMap — can be tens of GB for whole-genome annotations).
    // Genome-wide copy counts, populated only when the genome pass runs. Kept
    // alive past the pass so the stats/histogram writers and the index can use it.
    let mut genome_copies: Option<HashMap<u64, u8>> = None;
    // Per-k-mer gene annotation (which full-GTF genes each occurrence falls in),
    // accumulated during the genome scan only when per-copy BEDs are requested.
    let mut kmer_annot: HashMap<u64, GeneAnnot> = HashMap::new();
    if args.check_genome_uniqueness {
        let seq_list = read_fai_sequences(&args.fasta)?;
        let total_bases: u64 = seq_list.iter().map(|(_, l)| *l as u64).sum();
        info!(
            n_sequences = seq_list.len(),
            total_bases,
            n_candidates = kmer_to_gene.len(),
            "genome-wide uniqueness pass enabled — this requires ~20 bytes × n_candidates RAM"
        );

        let mut genome_counts: HashMap<u64, u8> =
            kmer_to_gene.keys().map(|&kh| (kh, 0u8)).collect();

        const PROGRESS_INTERVAL: u64 = 100_000_000;
        let mut bases_done: u64 = 0;

        for (seq_name, seq_len) in &seq_list {
            info!(seq = %seq_name, len = seq_len, "scanning sequence...");
            let seq = match fai.fetch_seq(seq_name, 0, seq_len.saturating_sub(1)) {
                Ok(s) => s.to_vec(),
                Err(e) => {
                    tracing::warn!(seq = %seq_name, "failed to fetch sequence: {e}");
                    continue;
                }
            };

            let mut next_report = PROGRESS_INTERVAL;
            for (kmer, pos_in_seq) in kmer_iter(&seq, k) {
                if let Some(cnt) = genome_counts.get_mut(&kmer) {
                    *cnt = cnt.saturating_add(1);
                    // Record which full-GTF gene(s) this occurrence falls in (or
                    // intergenic) so multi-copy BED intervals can name the paralogs.
                    if let Some(ref annot) = gene_annot {
                        let entry = kmer_annot.entry(kmer).or_default();
                        let hits = annot.genes_at(seq_name, pos_in_seq as usize);
                        if hits.is_empty() {
                            entry.intergenic = true;
                        } else {
                            entry.genes.extend(hits);
                        }
                    }
                }
                let pos_in_seq = pos_in_seq as u64;
                if pos_in_seq >= next_report {
                    let total_done = bases_done + pos_in_seq;
                    let pct = (total_done * 100) / total_bases.max(1);
                    info!(
                        chrom = %seq_name,
                        position = pos_in_seq,
                        pct_complete = pct,
                        "genome-wide uniqueness pass"
                    );
                    next_report += PROGRESS_INTERVAL;
                }
            }

            bases_done += *seq_len as u64;
            info!(seq = %seq_name, bases_done, "sequence done");
        }

        // Defer the retain until after stats are written, so the stats see the
        // full panel-unique survivor set together with every genome-wide count.
        genome_copies = Some(genome_counts);
    }

    // Per-gene uniqueness stats and the global copy histogram are written here —
    // before the genome-copy retain and the --min-gene-kmers prune — so the maps
    // still hold the full panel-unique set and every gene reports real values.
    if let Some(ref stats_path) = args.gene_stats_output {
        info!(output = %stats_path.display(), "writing per-gene uniqueness TSV...");
        write_gene_stats(
            &genes,
            &kmer_to_gene,
            &kmer_to_pos,
            k,
            genome_copies.as_ref(),
            args.max_genome_copies,
            stats_path,
        )?;
    }
    if let (Some(hist_path), Some(counts)) = (&args.copy_histogram_output, &genome_copies) {
        info!(output = %hist_path.display(), "writing copy-number histogram...");
        write_copy_histogram(counts, hist_path)?;
    }

    // Apply the genome-wide copy filter now that stats are recorded: retain k-mers
    // occurring 1..=max_genome_copies times genome-wide.
    if let Some(ref counts) = genome_copies {
        let max = args.max_genome_copies;
        kmer_to_gene.retain(|kmer, _| {
            let c = counts.get(kmer).copied().unwrap_or(0) as u32;
            c >= 1 && c <= max
        });
        kmer_to_pos.retain(|kmer, _| kmer_to_gene.contains_key(kmer));
        info!(
            n_unique_kmers = kmer_to_gene.len(),
            max_genome_copies = max,
            "k-mers after genome-wide copy filter"
        );
    }

    // Step 3 — drop genes with too few unique k-mers.
    let min_kmers = args.min_gene_kmers;
    if min_kmers > 0 {
        let mut gene_kmer_counts = vec![0u32; genes.len()];
        for &gene_idx in kmer_to_gene.values() {
            gene_kmer_counts[gene_idx as usize] += 1;
        }
        let low_genes: HashSet<u32> = gene_kmer_counts
            .iter()
            .enumerate()
            .filter(|(_, &cnt)| cnt < min_kmers)
            .map(|(i, _)| i as u32)
            .collect();
        let n_removed = low_genes.len();
        kmer_to_gene.retain(|_, gene_idx| !low_genes.contains(gene_idx));
        kmer_to_pos.retain(|kmer, _| kmer_to_gene.contains_key(kmer));
        info!(
            n_genes_removed = n_removed,
            min_kmers,
            "genes removed due to insufficient unique k-mers"
        );
    }

    // Step 4 — write the index.
    info!(output = %args.output.display(), "writing fusion index...");
    write_index(&kmer_to_gene, &kmer_to_pos, &genes, genome_copies.as_ref(), k, &args.output)?;

    // Step 5 — optional BED of merged k-mer intervals (all retained k-mers).
    if let Some(ref bed_path) = args.bed_output {
        info!(output = %bed_path.display(), "writing k-mer coverage BED...");
        write_kmer_bed(&kmer_to_pos, &seq_list_for_genes, k, bed_path)?;
    }

    // Step 5b — optional per-copy-tier BEDs (requires the genome pass).
    if let Some(ref prefix) = args.bed_output_by_copies {
        if let Some(ref counts) = genome_copies {
            info!(prefix = %prefix.display(), "writing per-copy-tier k-mer BEDs...");
            let (annot, names): (Option<&HashMap<u64, GeneAnnot>>, &[String]) =
                match gene_annot {
                    Some(ref gi) => (Some(&kmer_annot), gi.names.as_slice()),
                    None => (None, &[]),
                };
            write_kmer_bed_by_copies(
                &kmer_to_pos,
                counts,
                &seq_list_for_genes,
                k,
                args.max_genome_copies,
                prefix,
                annot,
                names,
            )?;
        }
    }

    info!("done");
    Ok(())
}

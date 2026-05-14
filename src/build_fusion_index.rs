use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader};
use std::path::Path;

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
    output: &Path,
) -> Result<()> {
    if output.exists() {
        std::fs::remove_file(output)
            .with_context(|| format!("failed to remove existing output: {}", output.display()))?;
    }

    let conn = Connection::open(output)
        .with_context(|| format!("failed to open DuckDB: {}", output.display()))?;

    conn.execute_batch(
        "CREATE TABLE genes (gene_index UINTEGER, gene_name VARCHAR, chrom VARCHAR);
         CREATE TABLE kmers (kmer_hash BIGINT, gene_index UINTEGER);
         CREATE TABLE kmer_positions (kmer_hash BIGINT, chrom VARCHAR, pos INTEGER);",
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
            app.append_row(duckdb::params![kmer_hash as i64, gene_idx])?;
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

// ─── Main entry point ─────────────────────────────────────────────────────────

const MULTI_GENE: u32 = u32::MAX;

pub fn build_fusion_index(args: &BuildFusionIndexArgs) -> Result<()> {
    let k = args.kmer_size as usize;
    anyhow::ensure!(k > 0 && k <= 31, "--kmer-size must be between 1 and 31");

    info!("parsing gene bodies from GTF...");
    let mut genes = parse_gene_bodies(&args.gtf)?;
    info!(n_genes = genes.len(), "gene bodies loaded");

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
            for (offset, kmer) in kmer_iter(&chrom_seq[start..end_clamped], k).enumerate() {
                local
                    .entry(kmer)
                    .and_modify(|v| {
                        if *v != gene_idx {
                            *v = MULTI_GENE;
                        }
                    })
                    .or_insert(gene_idx);
                // Record the chromosome position of the first occurrence only.
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
            for (i, kmer) in kmer_iter(&seq, k).enumerate() {
                if let Some(cnt) = genome_counts.get_mut(&kmer) {
                    *cnt = cnt.saturating_add(1);
                }
                let pos_in_seq = i as u64;
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

        kmer_to_gene.retain(|kmer, _| genome_counts.get(kmer).copied().unwrap_or(0) == 1);
        kmer_to_pos.retain(|kmer, _| kmer_to_gene.contains_key(kmer));
        info!(
            n_unique_kmers = kmer_to_gene.len(),
            "k-mers after genome-wide uniqueness filter"
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
    write_index(&kmer_to_gene, &kmer_to_pos, &genes, &args.output)?;
    info!("done");
    Ok(())
}

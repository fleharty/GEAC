use std::collections::{HashMap, HashSet};
use std::path::Path;

use anyhow::{Context, Result};
use duckdb::Connection;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use tracing::info;

use crate::cli::ExtractGeneArgs;
use crate::kmer::kmer_iter;

// ─── Index loading ────────────────────────────────────────────────────────────

/// Load only the k-mers that belong to the requested target genes.
/// Returns (kmer_hash → gene_name, list of requested genes not found in the index).
fn load_target_kmers(
    index_path: &Path,
    target_genes: &[String],
) -> Result<(HashMap<u64, String>, Vec<String>)> {
    let conn = Connection::open(index_path)
        .with_context(|| format!("failed to open fusion index: {}", index_path.display()))?;

    let mut stmt = conn.prepare("SELECT gene_index, gene_name FROM genes")?;
    let mut gene_idx_to_name: HashMap<u32, String> = HashMap::new();
    let mut rows = stmt.query([])?;
    while let Some(row) = rows.next()? {
        let idx: u32 = row.get(0)?;
        let name: String = row.get(1)?;
        gene_idx_to_name.insert(idx, name);
    }

    let target_set: HashSet<&str> = target_genes.iter().map(|s| s.as_str()).collect();
    let target_indices: HashSet<u32> = gene_idx_to_name
        .iter()
        .filter(|(_, name)| target_set.contains(name.as_str()))
        .map(|(&idx, _)| idx)
        .collect();

    let unmatched: Vec<String> = target_genes
        .iter()
        .filter(|t| !gene_idx_to_name.values().any(|n| n == *t))
        .cloned()
        .collect();

    let mut stmt = conn.prepare("SELECT kmer_hash, gene_index FROM kmers")?;
    let mut kmer_to_gene: HashMap<u64, String> = HashMap::new();
    let mut rows = stmt.query([])?;
    while let Some(row) = rows.next()? {
        let kmer_i64: i64 = row.get(0)?;
        let gene_idx: u32 = row.get(1)?;
        if target_indices.contains(&gene_idx) {
            if let Some(name) = gene_idx_to_name.get(&gene_idx) {
                kmer_to_gene.insert(kmer_i64 as u64, name.clone());
            }
        }
    }

    info!(
        n_target_genes = target_indices.len(),
        n_target_kmers = kmer_to_gene.len(),
        "loaded target k-mers from index"
    );

    Ok((kmer_to_gene, unmatched))
}

// ─── K-mer helpers ────────────────────────────────────────────────────────────

fn decode_kmer(hash: u64, k: usize) -> String {
    const BASES: [u8; 4] = *b"ACGT";
    (0..k)
        .rev()
        .map(|i| BASES[((hash >> (2 * i)) & 3) as usize] as char)
        .collect()
}

// ─── Per-read hit detail (only collected when kmer_hits_output is requested) ──

struct ReadDetail {
    read_name: Vec<u8>,
    read_end: bool, // true = R1
    chrom: String,
    pos: i64, // 0-based
    kmer_hits: Vec<(u64, String)>, // (kmer_hash, gene_name)
}

// ─── TSV output ───────────────────────────────────────────────────────────────

fn write_kmer_hits_tsv(details: &[ReadDetail], sample_id: &str, k: usize, output: &Path) -> Result<()> {
    use std::io::{BufWriter, Write};
    let file = std::fs::File::create(output)
        .with_context(|| format!("failed to create kmer hits TSV: {}", output.display()))?;
    let mut w = BufWriter::new(file);
    writeln!(w, "gene_matched\tsample_id\tread_name\tread_end\tchrom\tpos\tkmer_hash\tkmer_seq")?;
    for d in details {
        let read_name = std::str::from_utf8(&d.read_name).unwrap_or("?");
        let read_end = if d.read_end { "R1" } else { "R2" };
        let pos1 = d.pos + 1;
        for (kmer_hash, gene_name) in &d.kmer_hits {
            let kmer_seq = decode_kmer(*kmer_hash, k);
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                gene_name,
                sample_id,
                read_name,
                read_end,
                d.chrom,
                pos1,
                *kmer_hash as i64,
                kmer_seq
            )?;
        }
    }
    Ok(())
}

// ─── Main entry point ─────────────────────────────────────────────────────────

pub fn extract_gene(args: &ExtractGeneArgs) -> Result<()> {
    let k = args.kmer_size as usize;
    anyhow::ensure!(k > 0 && k <= 31, "--kmer-size must be between 1 and 31");

    info!(index = %args.index.display(), genes = ?args.genes, "loading target k-mers...");
    let (kmer_to_gene, unmatched) = load_target_kmers(&args.index, &args.genes)?;

    if !unmatched.is_empty() {
        anyhow::bail!(
            "the following gene(s) were not found in the index: {}",
            unmatched.join(", ")
        );
    }
    if kmer_to_gene.is_empty() {
        anyhow::bail!("no k-mers found for the requested gene(s); check that the index was built with these genes");
    }

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

    // Pre-extract target names before entering the records loop.
    let target_names: Vec<String> = {
        let hdr = reader.header();
        (0..hdr.target_count())
            .map(|i| std::str::from_utf8(hdr.tid2name(i)).unwrap_or("?").to_string())
            .collect()
    };

    info!(
        sample_id = %sample_id,
        bam = %args.bam.display(),
        n_target_kmers = kmer_to_gene.len(),
        "scanning BAM for reads matching target gene(s)..."
    );

    let collect_detail = args.kmer_hits_output.is_some();
    let mut matching_qnames: HashSet<Vec<u8>> = HashSet::new();
    let mut read_details: Vec<ReadDetail> = Vec::new();
    let mut reads_processed: u64 = 0;
    let mut reads_matched: u64 = 0;

    for result in reader.records() {
        let record = result.context("error reading BAM record")?;
        reads_processed += 1;

        let flags = record.flags();
        // Skip secondary (0x100), supplementary (0x800), duplicate (0x400), unmapped (0x4)
        if flags & (0x100 | 0x800 | 0x400 | 0x4) != 0 {
            continue;
        }
        if record.mapq() < args.min_mapq {
            continue;
        }

        let seq = record.seq().as_bytes();
        let mut gene_hits: HashMap<&str, u32> = HashMap::new();
        let mut hit_list: Vec<(u64, String)> = Vec::new();

        for kmer in kmer_iter(&seq, k) {
            if let Some(gene_name) = kmer_to_gene.get(&kmer) {
                *gene_hits.entry(gene_name.as_str()).or_insert(0) += 1;
                if collect_detail {
                    hit_list.push((kmer, gene_name.clone()));
                }
            }
        }

        // A read matches if any target gene has at least min_kmer_hits.
        let any_match = gene_hits.values().any(|&cnt| cnt >= args.min_kmer_hits);
        if !any_match {
            continue;
        }

        matching_qnames.insert(record.qname().to_vec());
        reads_matched += 1;

        if collect_detail {
            let chrom = {
                let tid = record.tid();
                if tid >= 0 {
                    target_names.get(tid as usize).cloned().unwrap_or_else(|| "*".to_string())
                } else {
                    "*".to_string()
                }
            };
            read_details.push(ReadDetail {
                read_name: record.qname().to_vec(),
                read_end: flags & 0x40 != 0,
                chrom,
                pos: record.pos(),
                kmer_hits: hit_list,
            });
        }

        if reads_processed % 5_000_000 == 0 {
            info!(reads_processed, reads_matched, "BAM scan progress");
        }
    }

    info!(
        reads_processed,
        reads_matched,
        n_matching_fragments = matching_qnames.len(),
        "BAM scan complete"
    );

    // Write kmer hits TSV before the second BAM pass (reader is consumed).
    if let Some(ref tsv_path) = args.kmer_hits_output {
        // Sort by read_name for a stable output order.
        read_details.sort_by(|a, b| a.read_name.cmp(&b.read_name));
        write_kmer_hits_tsv(&read_details, &sample_id, k, tsv_path)?;
        info!(output = %tsv_path.display(), "kmer hits TSV written");
    }

    // Second BAM pass: write all reads from matching fragments (both R1 and R2).
    info!(
        n_fragments = matching_qnames.len(),
        output = %args.output.display(),
        "writing matching reads (second BAM pass)..."
    );

    let mut reader2 = bam::Reader::from_path(&args.bam)
        .with_context(|| format!("failed to re-open BAM: {}", args.bam.display()))?;
    if let Some(ref ref_path) = args.reference {
        reader2.set_reference(ref_path)?;
    }

    let header = bam::Header::from_template(reader2.header());
    let mut writer = bam::Writer::from_path(&args.output, &header, bam::Format::Bam)
        .with_context(|| format!("failed to create output BAM: {}", args.output.display()))?;

    let mut reads_written: u64 = 0;
    for result in reader2.records() {
        let record = result.context("error reading BAM record in second pass")?;
        if matching_qnames.contains(record.qname()) {
            writer.write(&record).context("failed to write BAM record")?;
            reads_written += 1;
        }
    }

    info!(reads_written, output = %args.output.display(), "done");
    Ok(())
}

use std::collections::{HashMap, HashSet};
use std::path::Path;

use anyhow::{Context, Result};
use duckdb::Connection;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use tracing::info;

use crate::cli::ExtractGeneArgs;
use crate::kmer::kmer_iter;

// ─── Region parsing ──────────────────────────────────────────────────────────

/// Parse a "chr:start-end" string (1-based inclusive) into
/// (chrom, start_0_inclusive, end_0_exclusive).
fn parse_region(s: &str) -> Result<(String, i64, i64)> {
    let (chrom, range) = s
        .rsplit_once(':')
        .with_context(|| format!("invalid region '{}': expected 'chr:start-end'", s))?;
    let (start_s, end_s) = range
        .split_once('-')
        .with_context(|| format!("invalid region '{}': expected 'chr:start-end'", s))?;
    let start_1: i64 = start_s
        .replace(',', "")
        .parse()
        .with_context(|| format!("invalid start in region '{}'", s))?;
    let end_1: i64 = end_s
        .replace(',', "")
        .parse()
        .with_context(|| format!("invalid end in region '{}'", s))?;
    anyhow::ensure!(
        start_1 >= 1 && end_1 >= start_1,
        "invalid region '{}': start must be ≥ 1 and end ≥ start",
        s
    );
    Ok((chrom.to_string(), start_1 - 1, end_1))
}

/// Short non-reversible identifier for a qname. Uses FNV-1a 64-bit and
/// renders as the low 8 hex digits — enough to correlate multiple records
/// from the same fragment in a log file but not to recover the original
/// qname or look it up in the BAM.
fn short_hash(bytes: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bytes {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{:08x}", h as u32)
}

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
    pos: i64,                      // 0-based
    kmer_hits: Vec<(u64, String)>, // (kmer_hash, gene_name)
}

// ─── TSV output ───────────────────────────────────────────────────────────────

fn write_kmer_hits_tsv(
    details: &[ReadDetail],
    sample_id: &str,
    k: usize,
    output: &Path,
) -> Result<()> {
    use std::io::{BufWriter, Write};
    let file = std::fs::File::create(output)
        .with_context(|| format!("failed to create kmer hits TSV: {}", output.display()))?;
    let mut w = BufWriter::new(file);
    writeln!(
        w,
        "gene_matched\tsample_id\tread_name\tread_end\tchrom\tpos\tkmer_hash\tkmer_seq"
    )?;
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
            found.context("--sample-id not provided and no SM tag found in BAM/CRAM @RG header")?
        }
    };

    // Pre-extract target names before entering the records loop.
    let target_names: Vec<String> = {
        let hdr = reader.header();
        (0..hdr.target_count())
            .map(|i| {
                std::str::from_utf8(hdr.tid2name(i))
                    .unwrap_or("?")
                    .to_string()
            })
            .collect()
    };

    info!(
        sample_id = %sample_id,
        bam = %args.bam.display(),
        n_target_kmers = kmer_to_gene.len(),
        "scanning BAM for reads matching target gene(s)..."
    );

    let collect_detail = args.kmer_hits_output.is_some();
    let debug_region = args.debug_region.as_deref().map(parse_region).transpose()?;
    if let Some((ref c, s, e)) = debug_region {
        info!(chrom = %c, start_0based = s, end_0based_excl = e, "debug-region active");
    }
    let mut matching_qnames: HashSet<Vec<u8>> = HashSet::new();
    let mut read_details: Vec<ReadDetail> = Vec::new();
    let mut reads_processed: u64 = 0;
    let mut reads_scanned: u64 = 0; // passed flag/mapq filter; eligible for k-mer matching
    let mut reads_matched: u64 = 0;
    let mut debug_reads_in_region_total: u64 = 0;
    let mut debug_reads_in_region_zero_kmers: u64 = 0;
    let mut debug_reads_in_region_zero_matches: u64 = 0;
    let mut debug_total_n_acgt: u64 = 0;
    let mut debug_total_n_n: u64 = 0;
    let mut debug_total_n_equals: u64 = 0;
    let mut debug_total_n_other: u64 = 0;

    for result in reader.records() {
        let record = result.context("error reading BAM record")?;
        reads_processed += 1;

        let flags = record.flags();
        // Skip secondary (0x100), duplicate (0x400), unmapped (0x4).
        // Supplementary alignments (0x800) are intentionally NOT filtered: for chimeric
        // reads whose primary alignment is outside the target gene, the supplementary
        // alignment may be the only record that covers a unique k-mer and we must scan
        // it to collect the qname.
        if flags & (0x100 | 0x400 | 0x4) != 0 {
            continue;
        }
        if record.mapq() < args.min_mapq {
            continue;
        }

        reads_scanned += 1;
        let seq = record.seq().as_bytes();
        let mut gene_hits: HashMap<&str, u32> = HashMap::new();
        let mut hit_list: Vec<(u64, String)> = Vec::new();
        let mut n_kmers_extracted: u32 = 0;
        let mut n_kmer_matches: u32 = 0;

        for (kmer, _) in kmer_iter(&seq, k) {
            n_kmers_extracted += 1;
            if let Some(gene_name) = kmer_to_gene.get(&kmer) {
                n_kmer_matches += 1;
                *gene_hits.entry(gene_name.as_str()).or_insert(0) += 1;
                if collect_detail {
                    hit_list.push((kmer, gene_name.clone()));
                }
            }
        }

        // Diagnostic: log reads overlapping the debug region.
        // Output is REDACTED: qnames are replaced by a short non-reversible hash
        // and read sequences are reported only as anonymized base-class counts
        // (n_acgt, n_n, n_equals, n_other). No raw sequence or read identifier
        // is logged, so the output can be shared without exposing BAM contents.
        if let Some((ref dr_chrom, dr_start, dr_end)) = debug_region {
            let tid = record.tid();
            let chrom_name: Option<&str> = if tid >= 0 {
                target_names.get(tid as usize).map(|s| s.as_str())
            } else {
                None
            };
            if chrom_name == Some(dr_chrom.as_str()) {
                let aln_start = record.pos();
                let aln_end = record.cigar().end_pos();
                if aln_start < dr_end && aln_end > dr_start {
                    debug_reads_in_region_total += 1;
                    if n_kmers_extracted == 0 {
                        debug_reads_in_region_zero_kmers += 1;
                    }
                    if n_kmer_matches == 0 {
                        debug_reads_in_region_zero_matches += 1;
                    }

                    // Anonymized base-class counts.
                    let (mut n_acgt, mut n_n, mut n_equals, mut n_other) = (0u32, 0u32, 0u32, 0u32);
                    for &b in seq.iter() {
                        match b | 0x20 {
                            b'a' | b'c' | b'g' | b't' => n_acgt += 1,
                            b'n' => n_n += 1,
                            b'=' => n_equals += 1,
                            _ => n_other += 1,
                        }
                    }
                    debug_total_n_acgt += n_acgt as u64;
                    debug_total_n_n += n_n as u64;
                    debug_total_n_equals += n_equals as u64;
                    debug_total_n_other += n_other as u64;

                    // Short non-reversible identifier so multiple records from
                    // the same fragment can be correlated in the log without
                    // revealing the original qname.
                    let qname_id = short_hash(record.qname());

                    info!(
                        qname_id = %qname_id,
                        flags = format!("0x{:x}", flags),
                        mapq = record.mapq(),
                        pos_1based = aln_start + 1,
                        aln_end_0based = aln_end,
                        seq_len = seq.len(),
                        n_kmers_extracted,
                        n_kmer_matches,
                        n_acgt,
                        n_n,
                        n_equals,
                        n_other,
                        "debug-region read"
                    );
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
                    target_names
                        .get(tid as usize)
                        .cloned()
                        .unwrap_or_else(|| "*".to_string())
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

        if reads_processed.is_multiple_of(5_000_000) {
            info!(reads_processed, reads_matched, "BAM scan progress");
        }
    }

    let reads_filtered = reads_processed - reads_scanned;
    let reads_no_kmer_hit = reads_scanned - reads_matched;
    info!(
        reads_processed,
        reads_filtered_flags_or_mapq = reads_filtered,
        reads_scanned_for_kmers = reads_scanned,
        reads_with_no_kmer_hit = reads_no_kmer_hit,
        reads_matched,
        n_matching_fragments = matching_qnames.len(),
        "BAM scan complete"
    );
    if debug_region.is_some() {
        info!(
            debug_reads_in_region_total,
            debug_reads_in_region_zero_kmers,
            debug_reads_in_region_zero_matches,
            total_bases_acgt = debug_total_n_acgt,
            total_bases_n = debug_total_n_n,
            total_bases_equals = debug_total_n_equals,
            total_bases_other = debug_total_n_other,
            "debug-region summary"
        );
    }

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

    let mut complement_writer: Option<bam::Writer> = args
        .complement_output
        .as_ref()
        .map(|p| {
            bam::Writer::from_path(p, &header, bam::Format::Bam)
                .with_context(|| format!("failed to create complement BAM: {}", p.display()))
        })
        .transpose()?;

    let mut reads_written: u64 = 0;
    let mut complement_written: u64 = 0;

    for result in reader2.records() {
        let record = result.context("error reading BAM record in second pass")?;
        if matching_qnames.contains(record.qname()) {
            writer
                .write(&record)
                .context("failed to write BAM record")?;
            reads_written += 1;
        } else if let Some(ref mut cw) = complement_writer {
            cw.write(&record)
                .context("failed to write complement BAM record")?;
            complement_written += 1;
        }
    }

    // Close writers before indexing — htslib requires files to be complete.
    drop(writer);
    info!(reads_written, output = %args.output.display(), "BAM written");
    bam::index::build(&args.output, None, bam::index::Type::Bai, 0)
        .with_context(|| format!("failed to index output BAM: {}", args.output.display()))?;
    info!(output = %args.output.display(), "BAI index written");

    if let Some(ref p) = args.complement_output {
        drop(complement_writer.take());
        bam::index::build(p, None, bam::index::Type::Bai, 0)
            .with_context(|| format!("failed to index complement BAM: {}", p.display()))?;
        info!(complement_written, output = %p.display(), "complement BAM written and indexed");
    }

    Ok(())
}

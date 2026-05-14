use std::collections::{HashMap, HashSet};
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
}

fn load_index(index_path: &Path) -> Result<FusionIndex> {
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

    let mut stmt = conn.prepare("SELECT kmer_hash, gene_index FROM kmers")?;
    let mut kmer_to_gene: HashMap<u64, u32> = HashMap::new();
    let mut rows = stmt.query([])?;
    while let Some(row) = rows.next()? {
        let kmer_i64: i64 = row.get(0)?;
        let gene_idx: u32 = row.get(1)?;
        kmer_to_gene.insert(kmer_i64 as u64, gene_idx);
    }
    info!(n_kmers = kmer_to_gene.len(), "loaded k-mer index");

    Ok(FusionIndex {
        kmer_to_gene,
        gene_names,
        gene_chroms,
    })
}

// ─── Read → gene assignment ───────────────────────────────────────────────────

struct ReadHit {
    gene_idx: u32,
    mapq: u8,
    chrom: String,
    pos: i64,   // 0-based leftmost position from BAM record
    is_read1: bool,
    kmer_hits: Vec<(u64, u32)>, // (kmer_hash, matched_gene_idx) for every kmer that hit any gene
}

fn assign_gene(seq: &[u8], index: &FusionIndex, k: usize, min_hits: u32) -> Option<ReadHit> {
    let mut gene_hits: HashMap<u32, u32> = HashMap::new();
    let mut kmer_hits: Vec<(u64, u32)> = Vec::new();
    for kmer in kmer_iter(seq, k) {
        if let Some(&gene_idx) = index.kmer_to_gene.get(&kmer) {
            *gene_hits.entry(gene_idx).or_insert(0) += 1;
            kmer_hits.push((kmer, gene_idx));
        }
    }
    let (winning_gene, count) = gene_hits
        .into_iter()
        .max_by_key(|&(_, count)| count)?;
    if count < min_hits {
        return None;
    }
    Some(ReadHit { gene_idx: winning_gene, mapq: 0, chrom: String::new(), pos: 0, is_read1: false, kmer_hits })
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
    writeln!(w, "sample_id\tgene_a\tgene_b\tchrom_a\tchrom_b\tsupporting_reads\tmin_mapq")?;
    for r in records {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.sample_id, r.gene_a, r.gene_b, r.chrom_a, r.chrom_b, r.supporting_reads, r.min_mapq
        )?;
    }
    Ok(())
}

// ─── K-mer hit detail TSV ─────────────────────────────────────────────────────

fn write_kmer_hits_tsv(
    fusion_reads: &[(&[u8], &Vec<ReadHit>, &str)], // (qname, hits_for_fragment, fusion_label)
    sample_id: &str,
    k: usize,
    gene_names: &[String],
    output: &Path,
) -> Result<()> {
    use std::io::{BufWriter, Write};
    let file = std::fs::File::create(output)
        .with_context(|| format!("failed to create kmer hits TSV: {}", output.display()))?;
    let mut w = BufWriter::new(file);
    writeln!(w, "fusion\tsample_id\tread_name\tread_end\tchrom\tpos\tgene_matched\tkmer_hash\tkmer_seq")?;
    for (qname, read_hits, fusion_label) in fusion_reads {
        let read_name = std::str::from_utf8(qname).unwrap_or("?");
        for rh in *read_hits {
            let read_end = if rh.is_read1 { "R1" } else { "R2" };
            let pos1 = rh.pos + 1; // convert 0-based to 1-based
            for &(kmer_hash, gene_idx) in &rh.kmer_hits {
                let gene_name = gene_names
                    .get(gene_idx as usize)
                    .map(|s| s.as_str())
                    .unwrap_or("?");
                let kmer_seq = decode_kmer(kmer_hash, k);
                writeln!(
                    w,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    fusion_label,
                    sample_id,
                    read_name,
                    read_end,
                    rh.chrom,
                    pos1,
                    gene_name,
                    kmer_hash as i64,
                    kmer_seq
                )?;
            }
        }
    }
    Ok(())
}

// ─── Parquet output ───────────────────────────────────────────────────────────

struct FusionRecord {
    sample_id: String,
    gene_a: String,
    gene_b: String,
    chrom_a: String,
    chrom_b: String,
    supporting_reads: i32,
    min_mapq: u8,
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
        ],
    )
    .context("failed to build record batch")?;

    writer.write(&batch).context("failed to write batch")?;
    writer.close().context("failed to finalize Parquet")?;
    Ok(())
}

// ─── Main entry point ─────────────────────────────────────────────────────────

pub fn detect_fusions(args: &FusionsArgs) -> Result<()> {
    let k = args.kmer_size as usize;
    anyhow::ensure!(k > 0 && k <= 31, "--kmer-size must be between 1 and 31");

    info!(index = %args.index.display(), "loading fusion k-mer index...");
    let index = load_index(&args.index)?;

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

    // Pre-extract target names so we can look up chrom during the records loop
    // without conflicting borrows on reader.
    let target_names: Vec<String> = {
        let hdr = reader.header();
        (0..hdr.target_count())
            .map(|i| std::str::from_utf8(hdr.tid2name(i)).unwrap_or("?").to_string())
            .collect()
    };

    // qname bytes → list of ReadHit for reads that hit a gene
    let mut qname_to_reads: HashMap<Vec<u8>, Vec<ReadHit>> = HashMap::new();
    let mut reads_processed: u64 = 0;
    let mut reads_assigned: u64 = 0;

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
        if let Some(mut rh) = assign_gene(&seq, &index, k, args.min_kmer_hits) {
            rh.mapq = record.mapq();
            rh.chrom = {
                let tid = record.tid();
                if tid >= 0 {
                    target_names.get(tid as usize).cloned().unwrap_or_else(|| "*".to_string())
                } else {
                    "*".to_string()
                }
            };
            rh.pos = record.pos();
            rh.is_read1 = flags & 0x40 != 0;
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
    // count it as a fusion candidate.
    let mut fusion_counts: HashMap<(u32, u32), (u32, u8)> = HashMap::new();

    for read_hits in qname_to_reads.values() {
        if read_hits.len() < 2 {
            continue;
        }

        let mut gene_vote: HashMap<u32, usize> = HashMap::new();
        let min_mq = read_hits.iter().map(|rh| rh.mapq).min().unwrap_or(0);
        for rh in read_hits {
            *gene_vote.entry(rh.gene_idx).or_insert(0) += 1;
        }
        if gene_vote.len() < 2 {
            continue;
        }

        let mut sorted: Vec<(u32, usize)> = gene_vote.into_iter().collect();
        sorted.sort_by(|a, b| b.1.cmp(&a.1));

        let ga = sorted[0].0;
        let gb = sorted[1].0;
        let key = if ga <= gb { (ga, gb) } else { (gb, ga) };

        let entry = fusion_counts.entry(key).or_insert((0, 255));
        entry.0 += 1;
        entry.1 = entry.1.min(min_mq);
    }

    // Determine which (gene_a, gene_b) pairs pass the supporting-read filter.
    let passing_pairs: HashSet<(u32, u32)> = fusion_counts
        .iter()
        .filter(|(_, (count, _))| *count >= args.min_supporting_reads)
        .map(|(pair, _)| *pair)
        .collect();

    let mut records: Vec<FusionRecord> = fusion_counts
        .into_iter()
        .filter(|(_, (count, _))| *count >= args.min_supporting_reads)
        .map(|((ga, gb), (count, min_mq))| FusionRecord {
            sample_id: sample_id.clone(),
            gene_a: index.gene_names[ga as usize].clone(),
            gene_b: index.gene_names[gb as usize].clone(),
            chrom_a: index.gene_chroms[ga as usize].clone(),
            chrom_b: index.gene_chroms[gb as usize].clone(),
            supporting_reads: count as i32,
            min_mapq: min_mq,
        })
        .collect();

    records.sort_by(|a, b| b.supporting_reads.cmp(&a.supporting_reads));

    info!(n_fusions = records.len(), "fusion candidates identified");
    write_fusion_parquet(&records, &args.output)?;
    info!(output = %args.output.display(), "fusion candidates written");

    if let Some(ref tsv_path) = args.tsv_output {
        write_fusion_tsv(&records, tsv_path)?;
        info!(output = %tsv_path.display(), "fusion TSV written");
    }

    // Build a helper: (min_gene_idx, max_gene_idx) → "GENE_A::GENE_B" label.
    let fusion_label: HashMap<(u32, u32), String> = passing_pairs
        .iter()
        .map(|&(ga, gb)| {
            let na = &index.gene_names[ga as usize];
            let nb = &index.gene_names[gb as usize];
            let label = if na <= nb {
                format!("{}::{}", na, nb)
            } else {
                format!("{}::{}", nb, na)
            };
            ((ga, gb), label)
        })
        .collect();

    // Helper closure: given per-fragment read hits, return the passing fusion key if any.
    let fragment_fusion_key = |read_hits: &Vec<ReadHit>| -> Option<(u32, u32)> {
        let mut gene_vote: HashMap<u32, usize> = HashMap::new();
        for rh in read_hits {
            *gene_vote.entry(rh.gene_idx).or_insert(0) += 1;
        }
        if gene_vote.len() < 2 {
            return None;
        }
        let mut sorted: Vec<(u32, usize)> = gene_vote.into_iter().collect();
        sorted.sort_by(|a, b| b.1.cmp(&a.1));
        let ga = sorted[0].0;
        let gb = sorted[1].0;
        let key = if ga <= gb { (ga, gb) } else { (gb, ga) };
        if passing_pairs.contains(&key) { Some(key) } else { None }
    };

    if let Some(ref kmer_hits_path) = args.kmer_hits_output {
        let mut fusion_reads: Vec<(&[u8], &Vec<ReadHit>, &str)> = Vec::new();
        for (qname, read_hits) in &qname_to_reads {
            if let Some(key) = fragment_fusion_key(read_hits) {
                if let Some(label) = fusion_label.get(&key) {
                    fusion_reads.push((qname.as_slice(), read_hits, label.as_str()));
                }
            }
        }
        fusion_reads.sort_by(|a, b| a.2.cmp(b.2).then(a.0.cmp(b.0)));
        write_kmer_hits_tsv(&fusion_reads, &sample_id, k, &index.gene_names, kmer_hits_path)?;
        info!(output = %kmer_hits_path.display(), "kmer hits TSV written");
    }

    // Optional second BAM pass: write all reads from fusion-supporting fragments.
    if let Some(ref reads_output) = args.reads_output {
        // Collect qnames of fragments that support a passing fusion.
        let fusion_qnames: HashSet<Vec<u8>> = qname_to_reads
            .iter()
            .filter(|(_, read_hits)| fragment_fusion_key(read_hits).is_some())
            .map(|(qname, _)| qname.clone())
            .collect();

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
            let record = result.context("error reading BAM record in second pass")?;
            if fusion_qnames.contains(record.qname()) {
                writer.write(&record).context("failed to write BAM record")?;
                reads_written += 1;
            }
        }

        info!(reads_written, output = %reads_output.display(), "fusion reads written");
    }

    Ok(())
}

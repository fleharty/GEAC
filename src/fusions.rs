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
    mapq: u8,
    // kmer_hits is populated only during the second BAM pass for kmer-hits-output;
    // empty during the first pass to avoid gigabytes of per-read detail on deep BAMs.
    kmer_hits: Vec<(u64, u32)>, // (kmer_hash, matched_gene_idx)
}

fn assign_gene(
    seq: &[u8],
    index: &FusionIndex,
    k: usize,
    min_hits: u32,
    collect_detail: bool,
    max_copies: Option<u32>,
) -> Option<ReadHit> {
    let mut gene_hits: HashMap<u32, u32> = HashMap::new();
    let mut kmer_hits: Vec<(u64, u32)> = Vec::new();
    for (kmer, _) in kmer_iter(seq, k) {
        if let Some(&gene_idx) = index.kmer_to_gene.get(&kmer) {
            // Optional genome-wide copy filter: skip k-mers that occur too often
            // (or whose copy number is unknown — NULL in the index).
            if let Some(max) = max_copies {
                match index.kmer_to_copies.get(&kmer) {
                    Some(&c) if c as u32 <= max => {}
                    _ => continue,
                }
            }
            *gene_hits.entry(gene_idx).or_insert(0) += 1;
            // Retaining every matching k-mer for every assigned read costs
            // gigabytes on deep BAMs; only do it when the detail TSV is requested.
            if collect_detail {
                kmer_hits.push((kmer, gene_idx));
            }
        }
    }
    // Highest k-mer count wins; ties broken by lowest gene index so the choice is
    // deterministic (HashMap iteration order must not influence the result).
    let (winning_gene, count) = gene_hits
        .into_iter()
        .max_by(|a, b| a.1.cmp(&b.1).then_with(|| b.0.cmp(&a.0)))?;
    if count < min_hits {
        return None;
    }
    Some(ReadHit { gene_idx: winning_gene, mapq: 0, kmer_hits })
}

/// Select the fusion gene-pair for a fragment's read hits.
///
/// Votes are gene → number of reads hitting it. Returns `None` for fragments hitting
/// fewer than two distinct genes. The top two genes are chosen by vote count
/// descending, breaking ties by **ascending gene index** so the result is
/// deterministic (HashMap iteration order must not influence it). The pair is
/// normalized to `(min_index, max_index)` so it is a stable key for the whole run.
fn fragment_top_pair(read_hits: &[ReadHit]) -> Option<(u32, u32)> {
    if read_hits.len() < 2 {
        return None;
    }
    let mut gene_vote: HashMap<u32, usize> = HashMap::new();
    for rh in read_hits {
        *gene_vote.entry(rh.gene_idx).or_insert(0) += 1;
    }
    if gene_vote.len() < 2 {
        return None;
    }
    let mut sorted: Vec<(u32, usize)> = gene_vote.into_iter().collect();
    // Vote count descending, then gene index ascending (deterministic tie-break).
    sorted.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    let ga = sorted[0].0;
    let gb = sorted[1].0;
    Some(if ga <= gb { (ga, gb) } else { (gb, ga) })
}

/// Canonical `GENE_A::GENE_B` label for a normalized `(min_index, max_index)` pair.
/// Uses gene-index order so every output — Parquet, TSV, the `FX` BAM tag, and the
/// kmer-hits TSV — labels the same fusion identically. (The Parquet `gene_a`/`gene_b`
/// columns are filled from `gene_names[key.0]`/`[key.1]`, matching this order.)
fn fusion_pair_label(key: (u32, u32), gene_names: &[String]) -> String {
    format!("{}::{}", gene_names[key.0 as usize], gene_names[key.1 as usize])
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
    // count it as a fusion candidate.
    let mut fusion_counts: HashMap<(u32, u32), (u32, u8)> = HashMap::new();

    for read_hits in qname_to_reads.values() {
        let Some(key) = fragment_top_pair(read_hits) else {
            continue;
        };
        let min_mq = read_hits.iter().map(|rh| rh.mapq).min().unwrap_or(0);
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

    // Canonical "GENE_A::GENE_B" label per passing pair — same gene-index order as
    // the Parquet/TSV gene_a/gene_b columns, so all outputs agree on the label.
    let fusion_label: HashMap<(u32, u32), String> = passing_pairs
        .iter()
        .map(|&key| (key, fusion_pair_label(key, &index.gene_names)))
        .collect();

    // Given a fragment's read hits, return its fusion key iff it passes the filter.
    let fragment_fusion_key = |read_hits: &[ReadHit]| -> Option<(u32, u32)> {
        let key = fragment_top_pair(read_hits)?;
        passing_pairs.contains(&key).then_some(key)
    };

    // Optional second BAM pass: write all reads from fusion-supporting fragments.
    if let Some(ref reads_output) = args.reads_output {
        // Collect qnames of fragments that support a passing fusion, mapping each
        // to its "GENE_A::GENE_B" label so we can tag every emitted record with FX.
        let fusion_qnames: HashMap<Vec<u8>, String> = qname_to_reads
            .iter()
            .filter_map(|(qname, read_hits)| {
                let key = fragment_fusion_key(read_hits)?;
                let label = fusion_label.get(&key)?;
                Some((qname.clone(), label.clone()))
            })
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

    // Optional second BAM pass: stream per-k-mer hit rows for fusion-supporting reads.
    // Doing this in a second pass (rather than accumulating kmer_hits during the first
    // pass) avoids storing gigabytes of per-read detail for reads that turn out not to
    // support any passing fusion.
    if let Some(ref kmer_hits_path) = args.kmer_hits_output {
        use std::io::{BufWriter, Write};

        let fusion_qnames: HashMap<Vec<u8>, String> = qname_to_reads
            .iter()
            .filter_map(|(qname, read_hits)| {
                let key = fragment_fusion_key(read_hits)?;
                let label = fusion_label.get(&key)?;
                Some((qname.clone(), label.clone()))
            })
            .collect();

        info!(
            n_fragments = fusion_qnames.len(),
            output = %kmer_hits_path.display(),
            "writing kmer hits TSV (second BAM pass)..."
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

        let file = std::fs::File::create(kmer_hits_path)
            .with_context(|| format!("failed to create kmer hits TSV: {}", kmer_hits_path.display()))?;
        let mut w = BufWriter::new(file);
        writeln!(w, "fusion\tsample_id\tread_name\tread_end\tchrom\tpos\tgene_matched\tkmer_hash\tkmer_seq")?;

        let mut rows_written: u64 = 0;
        for result in reader2.records() {
            let record = result.context("error reading BAM record in kmer hits pass")?;
            let Some(fusion_label_str) = fusion_qnames.get(record.qname()) else {
                continue;
            };
            let seq = record.seq().as_bytes();
            let Some(rh) = assign_gene(&seq, &index, k, args.min_kmer_hits, true, args.max_kmer_copies) else {
                continue;
            };
            let flags = record.flags();
            let read_end = if flags & 0x40 != 0 { "R1" } else { "R2" };
            let tid = record.tid();
            let chrom = if tid >= 0 {
                target_names.get(tid as usize).map(|s| s.as_str()).unwrap_or("*")
            } else {
                "*"
            };
            let pos1 = record.pos() + 1;
            let read_name = std::str::from_utf8(record.qname()).unwrap_or("?");
            for &(kmer_hash, gene_idx) in &rh.kmer_hits {
                let gene_name = index.gene_names.get(gene_idx as usize).map(|s| s.as_str()).unwrap_or("?");
                let kmer_seq = decode_kmer(kmer_hash, k);
                writeln!(
                    w,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    fusion_label_str, sample_id, read_name, read_end, chrom, pos1,
                    gene_name, kmer_hash as i64, kmer_seq
                )?;
                rows_written += 1;
            }
        }
        info!(rows_written, output = %kmer_hits_path.display(), "kmer hits TSV written");
    }

    Ok(())
}

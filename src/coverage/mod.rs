use std::collections::{HashMap, HashSet};
use std::path::Path;

use anyhow::{Context, Result};
use rust_htslib::bam::Read;
use tracing::info;

use crate::bam::{open_bam, read_group_sample_id, RefCache};
use crate::cli::CoverageArgs;
use crate::gene_annotations::GeneAnnotations;
use crate::record::CoverageRecord;
use crate::targets::TargetIntervals;
use crate::writer::parquet_coverage::CoverageWriter;

/// Process a BAM/CRAM file and return per-position coverage records.
///
/// When `--targets` is provided, zero-depth positions within target intervals are
/// included in the output (so dropout is captured). Without targets, only positions
/// covered by at least one read are emitted.
pub fn collect_coverage(
    args: &CoverageArgs,
    target_intervals: Option<&TargetIntervals>,
    gene_annots: Option<&GeneAnnotations>,
    writer: &mut CoverageWriter,
) -> Result<()> {
    let mut reader = open_bam(&args.input, &args.reference)?;
    let mut ref_cache = RefCache::new(&args.reference)?;

    let sample_id = match &args.sample_id {
        Some(id) => id.clone(),
        None => read_group_sample_id(reader.header())
            .context("--sample-id not provided and no SM tag in BAM/CRAM header")?,
    };

    let bam_contigs: Vec<(String, usize)> = {
        let header = reader.header();
        (0..header.target_count())
            .map(|tid| {
                let name = std::str::from_utf8(header.tid2name(tid))
                    .unwrap_or("unknown")
                    .to_string();
                let len = header.target_len(tid).unwrap_or(0) as usize;
                (name, len)
            })
            .collect()
    };

    if let Some(region) = &args.region {
        reader
            .fetch(region.as_str())
            .with_context(|| format!("failed to fetch region '{region}'"))?;
    }

    let has_targets = target_intervals.is_some();
    // When targets are provided, track which positions had reads so we can fill zeros.
    let mut seen: HashSet<(String, i64)> = HashSet::new();

    let mut positions_written: u64 = 0;
    let mut positions_seen: u64 = 0;
    let start = std::time::Instant::now();
    let progress_interval = std::time::Duration::from_secs(args.progress_interval);
    let mut last_progress = start;
    let mut acc = BinAccumulator::new(args.bin_size);

    for pileup in reader.pileup() {
        let pileup = pileup.context("error reading pileup")?;
        let tid = pileup.tid() as usize;
        let pos = pileup.pos() as i64;

        let (chrom, _) = ref_cache.get(&bam_contigs, tid, pos as usize)?;

        // Skip positions outside targets when targets are provided
        if let Some(ti) = target_intervals {
            if !ti.contains(&chrom, pos) {
                continue;
            }
        }

        let tally = tally_coverage(&pileup, args.min_map_qual);

        if tally.total_depth < args.min_depth {
            if has_targets {
                seen.insert((chrom.clone(), pos));
            }
            continue;
        }

        let gc = compute_gc_content(ref_cache.current_seq(), pos as usize, args.gc_window);
        let on_target = target_intervals.map(|ti| ti.contains(&chrom, pos));
        let gene = gene_annots.and_then(|g| g.get(&chrom, pos).map(str::to_owned));

        if has_targets {
            seen.insert((chrom.clone(), pos));
        }

        let record = build_record(&sample_id, &chrom, pos, &tally, gc, on_target, gene, args);
        if args.adaptive_depth_threshold.map_or(false, |t| tally.total_depth < t) {
            // Flush any in-progress bin, then emit this position at single-base resolution
            if let Some(bin) = acc.finish() {
                writer.push(bin)?;
                positions_written += 1;
            }
            acc = BinAccumulator::new(args.bin_size);
            writer.push(record)?;
            positions_written += 1;
        } else {
            if let Some(bin) = acc.push(record) {
                writer.push(bin)?;
                positions_written += 1;
            }
        }

        positions_seen += 1;

        // Check clock every 10_000 positions to limit overhead; skip if interval is 0
        if args.progress_interval > 0 && positions_seen % 10_000 == 0 {
            let now = std::time::Instant::now();
            if now.duration_since(last_progress) >= progress_interval {
                info!(
                    positions_written = %fmt_with_commas(positions_written),
                    locus = %format!("{}:{}", chrom, fmt_with_commas(pos as u64)),
                    elapsed = %fmt_elapsed(start.elapsed()),
                    "coverage progress"
                );
                last_progress = now;
            }
        }
    }

    // Flush any partial bin from the pileup phase before zero-depth fill
    if let Some(bin) = acc.finish() {
        writer.push(bin)?;
        positions_written += 1;
    }

    // Fill in zero-depth positions within targets that had no reads
    if let Some(ti) = target_intervals {
        let n_total = ti.total_bases();
        let n_with_reads = seen.len();
        info!(
            n_total,
            n_with_reads,
            n_zero = n_total - n_with_reads.min(n_total),
            "filling zero-depth target positions"
        );

        // Skip zero-depth positions when min_depth > 0 (they'd always be filtered)
        if args.min_depth == 0 {
            let mut zero_ref = ZeroDepthRefReader::open(&args.reference)?;

            let mut zero_acc = BinAccumulator::new(args.bin_size);
            let mut zero_err: Option<anyhow::Error> = None;
            ti.for_each_position(|chrom, pos| {
                if zero_err.is_some() { return; }
                if seen.contains(&(chrom.to_string(), pos)) {
                    return;
                }
                let gc = zero_ref.gc_content(chrom, pos, args.gc_window);
                let gene = gene_annots.and_then(|g| g.get(chrom, pos).map(str::to_owned));
                let record = build_zero_record(&sample_id, chrom, pos, gc, Some(true), gene, args);
                if let Some(bin) = zero_acc.push(record) {
                    if let Err(e) = writer.push(bin) {
                        zero_err = Some(e);
                    }
                }
            });
            if let Some(e) = zero_err {
                return Err(e);
            }
            if let Some(bin) = zero_acc.finish() {
                writer.push(bin)?;
            }
        }
    }

    info!(positions_written, "coverage collection complete");
    Ok(())
}

// ── Bin accumulator ───────────────────────────────────────────────────────────

struct BinAccumulator {
    bin_size: i64,
    n: u32,
    chrom: String,
    bin_start: i64,
    // depth
    sum_total: f64,
    min_total: i32,
    max_total: i32,
    sum_fwd: f64,
    sum_rev: f64,
    // dup
    sum_raw: f64,
    sum_frac_dup: f64,
    // overlap
    sum_overlap: f64,
    sum_frac_overlap: f64,
    // mapq
    sum_mean_mapq: f64,
    sum_frac_mapq0: f64,
    sum_frac_low_mapq: f64,
    // base qual
    sum_mean_bq: f64,
    min_bq: i32,
    max_bq: i32,
    sum_frac_low_bq: f64,
    // insert size
    sum_mean_ins: f64,
    min_ins: i32,
    max_ins: i32,
    n_ins_obs: i32,
    n_ins_pos: u32,
    // gc
    sum_gc: f64,
    // annotations
    on_target: Option<bool>,
    gene: Option<String>,
    // provenance
    sample_id: String,
    read_type: crate::record::ReadType,
    pipeline: crate::record::Pipeline,
    batch: Option<String>,
}

impl BinAccumulator {
    fn new(bin_size: i64) -> Self {
        Self {
            bin_size,
            n: 0,
            chrom: String::new(),
            bin_start: 0,
            sum_total: 0.0, min_total: 0, max_total: 0,
            sum_fwd: 0.0, sum_rev: 0.0,
            sum_raw: 0.0, sum_frac_dup: 0.0,
            sum_overlap: 0.0, sum_frac_overlap: 0.0,
            sum_mean_mapq: 0.0, sum_frac_mapq0: 0.0, sum_frac_low_mapq: 0.0,
            sum_mean_bq: 0.0, min_bq: 0, max_bq: 0, sum_frac_low_bq: 0.0,
            sum_mean_ins: 0.0, min_ins: 0, max_ins: 0, n_ins_obs: 0, n_ins_pos: 0,
            sum_gc: 0.0,
            on_target: None, gene: None,
            sample_id: String::new(),
            read_type: crate::record::ReadType::Duplex,
            pipeline: crate::record::Pipeline::Fgbio,
            batch: None,
        }
    }

    /// Push a single-position record. Returns a completed bin record if a boundary was crossed.
    fn push(&mut self, r: CoverageRecord) -> Option<CoverageRecord> {
        let new_bin_start = (r.pos / self.bin_size) * self.bin_size;
        let boundary = self.n > 0 && (r.chrom != self.chrom || new_bin_start != self.bin_start);
        let completed = if boundary { Some(self.build()) } else { None };
        if boundary || self.n == 0 {
            self.reset(&r, new_bin_start);
        }
        self.accumulate(&r);
        completed
    }

    /// Flush any remaining accumulated positions as a final bin record.
    fn finish(&mut self) -> Option<CoverageRecord> {
        if self.n > 0 { Some(self.build()) } else { None }
    }

    fn reset(&mut self, r: &CoverageRecord, bin_start: i64) {
        self.chrom = r.chrom.clone();
        self.bin_start = bin_start;
        self.n = 0;
        self.sum_total = 0.0; self.min_total = i32::MAX; self.max_total = i32::MIN;
        self.sum_fwd = 0.0; self.sum_rev = 0.0;
        self.sum_raw = 0.0; self.sum_frac_dup = 0.0;
        self.sum_overlap = 0.0; self.sum_frac_overlap = 0.0;
        self.sum_mean_mapq = 0.0; self.sum_frac_mapq0 = 0.0; self.sum_frac_low_mapq = 0.0;
        self.sum_mean_bq = 0.0; self.min_bq = i32::MAX; self.max_bq = i32::MIN; self.sum_frac_low_bq = 0.0;
        self.sum_mean_ins = 0.0; self.min_ins = i32::MAX; self.max_ins = i32::MIN;
        self.n_ins_obs = 0; self.n_ins_pos = 0;
        self.sum_gc = 0.0;
        self.on_target = None; self.gene = None;
        self.sample_id = r.sample_id.clone();
        self.read_type = r.read_type;
        self.pipeline = r.pipeline;
        self.batch = r.batch.clone();
    }

    fn accumulate(&mut self, r: &CoverageRecord) {
        self.n += 1;
        self.sum_total += r.total_depth as f64;
        self.min_total = self.min_total.min(r.total_depth);
        self.max_total = self.max_total.max(r.total_depth);
        self.sum_fwd += r.fwd_depth as f64;
        self.sum_rev += r.rev_depth as f64;
        self.sum_raw += r.raw_read_depth as f64;
        self.sum_frac_dup += r.frac_dup as f64;
        self.sum_overlap += r.overlap_depth as f64;
        self.sum_frac_overlap += r.frac_overlap as f64;
        self.sum_mean_mapq += r.mean_mapq as f64;
        self.sum_frac_mapq0 += r.frac_mapq0 as f64;
        self.sum_frac_low_mapq += r.frac_low_mapq as f64;
        self.sum_mean_bq += r.mean_base_qual as f64;
        self.min_bq = self.min_bq.min(r.min_base_qual_obs);
        self.max_bq = self.max_bq.max(r.max_base_qual_obs);
        self.sum_frac_low_bq += r.frac_low_bq as f64;
        if r.n_insert_size_obs > 0 {
            self.sum_mean_ins += r.mean_insert_size as f64;
            self.min_ins = self.min_ins.min(r.min_insert_size);
            self.max_ins = self.max_ins.max(r.max_insert_size);
            self.n_ins_obs += r.n_insert_size_obs;
            self.n_ins_pos += 1;
        }
        self.sum_gc += r.gc_content as f64;
        if self.gene.is_none() {
            self.gene = r.gene.clone();
        }
        match (self.on_target, r.on_target) {
            (_, Some(true)) => self.on_target = Some(true),
            (None, v) => self.on_target = v,
            _ => {}
        }
    }

    fn build(&self) -> CoverageRecord {
        let n = self.n as f64;
        let mean = |s: f64| (s / n) as f32;
        let mean_i = |s: f64| (s / n).round() as i32;
        let (mean_ins, min_ins, max_ins) = if self.n_ins_pos > 0 {
            (
                (self.sum_mean_ins / self.n_ins_pos as f64) as f32,
                self.min_ins,
                self.max_ins,
            )
        } else {
            (0.0, 0, 0)
        };
        CoverageRecord {
            sample_id: self.sample_id.clone(),
            chrom: self.chrom.clone(),
            pos: self.bin_start,
            end: self.bin_start + self.bin_size,
            bin_n: self.n as i32,
            total_depth: mean_i(self.sum_total),
            min_depth: self.min_total,
            max_depth: self.max_total,
            fwd_depth: mean_i(self.sum_fwd),
            rev_depth: mean_i(self.sum_rev),
            raw_read_depth: mean_i(self.sum_raw),
            frac_dup: mean(self.sum_frac_dup),
            overlap_depth: mean_i(self.sum_overlap),
            frac_overlap: mean(self.sum_frac_overlap),
            mean_mapq: mean(self.sum_mean_mapq),
            frac_mapq0: mean(self.sum_frac_mapq0),
            frac_low_mapq: mean(self.sum_frac_low_mapq),
            mean_base_qual: mean(self.sum_mean_bq),
            min_base_qual_obs: if self.min_bq == i32::MAX { 0 } else { self.min_bq },
            max_base_qual_obs: if self.max_bq == i32::MIN { 0 } else { self.max_bq },
            frac_low_bq: mean(self.sum_frac_low_bq),
            mean_insert_size: mean_ins,
            min_insert_size: min_ins,
            max_insert_size: max_ins,
            n_insert_size_obs: self.n_ins_obs,
            gc_content: mean(self.sum_gc),
            on_target: self.on_target,
            gene: self.gene.clone(),
            read_type: self.read_type,
            pipeline: self.pipeline,
            batch: self.batch.clone(),
        }
    }
}

// ── Pileup tallying ───────────────────────────────────────────────────────────

struct CoverageTally {
    raw_read_depth: i32,
    dup_count:      i32,

    /// Mapqs of non-dup reads (before mapq filter) — for mappability stats
    mapqs_non_dup: Vec<u8>,

    /// Base qualities from non-dup, passing-mapq reads
    base_quals: Vec<u8>,

    /// Insert sizes from properly-paired, non-dup, passing-mapq, R1 reads
    insert_sizes: Vec<i32>,

    // Fragment-level depth counts (collapsed by qname)
    total_depth:   i32,
    fwd_depth:     i32,
    rev_depth:     i32,
    overlap_depth: i32,
}

struct FragEntry {
    is_reverse: bool,
    is_first_in_pair: bool,
}

fn tally_coverage(
    pileup: &rust_htslib::bam::pileup::Pileup,
    min_map_qual: u8,
) -> CoverageTally {
    let mut raw_read_depth = 0i32;
    let mut dup_count = 0i32;
    let mut mapqs_non_dup: Vec<u8> = Vec::new();
    let mut base_quals: Vec<u8> = Vec::new();
    let mut insert_sizes: Vec<i32> = Vec::new();
    let mut by_qname: HashMap<Vec<u8>, Vec<FragEntry>> = HashMap::new();

    for alignment in pileup.alignments() {
        let record = alignment.record();

        // Secondary and supplementary reads double-count positions; exclude them.
        if record.is_secondary() || record.is_supplementary() {
            continue;
        }

        raw_read_depth += 1;

        if record.is_duplicate() {
            dup_count += 1;
            continue;
        }

        let mapq = record.mapq();
        mapqs_non_dup.push(mapq);

        if mapq < min_map_qual {
            continue;
        }

        // Skip deletions/ref-skips after mapq check — still counted in raw/mapq stats.
        if alignment.is_del() || alignment.is_refskip() {
            continue;
        }

        // Base quality at the query position
        if let Some(qpos) = alignment.qpos() {
            base_quals.push(record.qual()[qpos]);
        }

        // Insert size: properly-paired R1 reads only (avoid double-counting)
        if record.is_paired()
            && record.is_proper_pair()
            && record.is_first_in_template()
        {
            let tlen = record.insert_size().unsigned_abs() as i32;
            if tlen > 0 {
                insert_sizes.push(tlen);
            }
        }

        by_qname
            .entry(record.qname().to_vec())
            .or_default()
            .push(FragEntry {
                is_reverse:       record.is_reverse(),
                is_first_in_pair: record.flags() & 0x40 != 0,
            });
    }

    // Collapse by qname for fragment-level depth counts
    let mut total_depth   = 0i32;
    let mut fwd_depth     = 0i32;
    let mut rev_depth     = 0i32;
    let mut overlap_depth = 0i32;

    for reads in by_qname.values() {
        total_depth += 1;
        match reads.as_slice() {
            [r] => {
                if r.is_reverse { rev_depth += 1; } else { fwd_depth += 1; }
            }
            [r1, r2] => {
                overlap_depth += 1;
                let r1_is_rev = if r1.is_first_in_pair { r1.is_reverse } else { r2.is_reverse };
                if r1_is_rev { rev_depth += 1; } else { fwd_depth += 1; }
            }
            _ => {
                // >2 reads with same qname: treat as non-overlapping
                if let Some(r) = reads.first() {
                    if r.is_reverse { rev_depth += 1; } else { fwd_depth += 1; }
                }
            }
        }
    }

    CoverageTally {
        raw_read_depth,
        dup_count,
        mapqs_non_dup,
        base_quals,
        insert_sizes,
        total_depth,
        fwd_depth,
        rev_depth,
        overlap_depth,
    }
}

// ── GC content ────────────────────────────────────────────────────────────────

/// Fraction of G+C bases in a window centred on `pos`.
pub(crate) fn compute_gc_content(seq: &[u8], pos: usize, window: usize) -> f32 {
    if seq.is_empty() || window == 0 {
        return 0.0;
    }
    let half = window / 2;
    let start = pos.saturating_sub(half);
    let end = (pos + half + 1).min(seq.len());
    let slice = &seq[start..end];
    if slice.is_empty() {
        return 0.0;
    }
    let gc = slice.iter().filter(|&&b| b == b'G' || b == b'C').count();
    gc as f32 / slice.len() as f32
}

// ── Record builders ───────────────────────────────────────────────────────────

fn build_record(
    sample_id: &str,
    chrom: &str,
    pos: i64,
    t: &CoverageTally,
    gc: f32,
    on_target: Option<bool>,
    gene: Option<String>,
    args: &CoverageArgs,
) -> CoverageRecord {
    let (mean_mapq, frac_mapq0, frac_low_mapq) = mapq_stats(&t.mapqs_non_dup, args.min_map_qual);
    let (mean_bq, min_bq, max_bq, frac_low_bq) = base_qual_stats(&t.base_quals, args.min_base_qual);
    let (mean_ins, min_ins, max_ins, n_ins) = insert_size_stats(&t.insert_sizes);

    let frac_dup = if t.raw_read_depth > 0 {
        t.dup_count as f32 / t.raw_read_depth as f32
    } else {
        0.0
    };
    let frac_overlap = if t.total_depth > 0 {
        t.overlap_depth as f32 / t.total_depth as f32
    } else {
        0.0
    };

    CoverageRecord {
        sample_id: sample_id.to_string(),
        chrom: chrom.to_string(),
        pos,
        end: pos + 1,
        bin_n: 1,
        total_depth: t.total_depth,
        min_depth: t.total_depth,
        max_depth: t.total_depth,
        fwd_depth: t.fwd_depth,
        rev_depth: t.rev_depth,
        raw_read_depth: t.raw_read_depth,
        frac_dup,
        overlap_depth: t.overlap_depth,
        frac_overlap,
        mean_mapq,
        frac_mapq0,
        frac_low_mapq,
        mean_base_qual: mean_bq,
        min_base_qual_obs: min_bq,
        max_base_qual_obs: max_bq,
        frac_low_bq,
        mean_insert_size: mean_ins,
        min_insert_size: min_ins,
        max_insert_size: max_ins,
        n_insert_size_obs: n_ins,
        gc_content: gc,
        on_target,
        gene,
        read_type: args.read_type,
        pipeline: args.pipeline,
        batch: args.batch.clone(),
    }
}

fn build_zero_record(
    sample_id: &str,
    chrom: &str,
    pos: i64,
    gc: f32,
    on_target: Option<bool>,
    gene: Option<String>,
    args: &CoverageArgs,
) -> CoverageRecord {
    CoverageRecord {
        sample_id: sample_id.to_string(),
        chrom: chrom.to_string(),
        pos,
        end: pos + 1,
        bin_n: 1,
        total_depth: 0,
        min_depth: 0,
        max_depth: 0,
        fwd_depth: 0,
        rev_depth: 0,
        raw_read_depth: 0,
        frac_dup: 0.0,
        overlap_depth: 0,
        frac_overlap: 0.0,
        mean_mapq: 0.0,
        frac_mapq0: 0.0,
        frac_low_mapq: 0.0,
        mean_base_qual: 0.0,
        min_base_qual_obs: 0,
        max_base_qual_obs: 0,
        frac_low_bq: 0.0,
        mean_insert_size: 0.0,
        min_insert_size: 0,
        max_insert_size: 0,
        n_insert_size_obs: 0,
        gc_content: gc,
        on_target,
        gene,
        read_type: args.read_type,
        pipeline: args.pipeline,
        batch: args.batch.clone(),
    }
}

// ── Stat helpers ──────────────────────────────────────────────────────────────

fn mapq_stats(mapqs: &[u8], min_map_qual: u8) -> (f32, f32, f32) {
    if mapqs.is_empty() {
        return (0.0, 0.0, 0.0);
    }
    let n = mapqs.len() as f32;
    let mean = mapqs.iter().map(|&q| q as f32).sum::<f32>() / n;
    let frac_0   = mapqs.iter().filter(|&&q| q == 0).count() as f32 / n;
    let frac_low = mapqs.iter().filter(|&&q| q < min_map_qual).count() as f32 / n;
    (mean, frac_0, frac_low)
}

fn base_qual_stats(bqs: &[u8], min_bq: u8) -> (f32, i32, i32, f32) {
    if bqs.is_empty() {
        return (0.0, 0, 0, 0.0);
    }
    let n = bqs.len() as f32;
    let mean = bqs.iter().map(|&q| q as f32).sum::<f32>() / n;
    let min_obs = *bqs.iter().min().unwrap() as i32;
    let max_obs = *bqs.iter().max().unwrap() as i32;
    let frac_low = bqs.iter().filter(|&&q| q < min_bq).count() as f32 / n;
    (mean, min_obs, max_obs, frac_low)
}

fn insert_size_stats(sizes: &[i32]) -> (f32, i32, i32, i32) {
    if sizes.is_empty() {
        return (0.0, 0, 0, 0);
    }
    let n = sizes.len() as i32;
    let mean = sizes.iter().map(|&v| v as f32).sum::<f32>() / n as f32;
    let min_v = *sizes.iter().min().unwrap();
    let max_v = *sizes.iter().max().unwrap();
    (mean, min_v, max_v, n)
}

// ── Formatting helpers ────────────────────────────────────────────────────────

fn fmt_elapsed(d: std::time::Duration) -> String {
    let s = d.as_secs();
    if s < 60 {
        format!("{}s", s)
    } else if s < 3600 {
        format!("{}m{}s", s / 60, s % 60)
    } else {
        format!("{}h{}m{}s", s / 3600, (s % 3600) / 60, s % 60)
    }
}

fn fmt_with_commas(n: u64) -> String {
    let s = n.to_string();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    for (i, ch) in s.chars().rev().enumerate() {
        if i > 0 && i % 3 == 0 {
            out.push(',');
        }
        out.push(ch);
    }
    out.chars().rev().collect()
}

// ── Zero-depth reference reader ───────────────────────────────────────────────

/// Minimal reference reader for GC content at zero-depth positions.
struct ZeroDepthRefReader {
    fai:          rust_htslib::faidx::Reader,
    cached_chrom: String,
    cached_seq:   Vec<u8>,
}

impl ZeroDepthRefReader {
    fn open(reference: &Path) -> Result<Self> {
        let fai = rust_htslib::faidx::Reader::from_path(reference)
            .with_context(|| format!("failed to open reference: {}", reference.display()))?;
        Ok(Self { fai, cached_chrom: String::new(), cached_seq: Vec::new() })
    }

    fn gc_content(&mut self, chrom: &str, pos: i64, window: usize) -> f32 {
        if self.cached_chrom != chrom {
            match self.fai.fetch_seq(chrom, 0, usize::MAX.saturating_sub(1)) {
                Ok(seq) => {
                    let mut s = seq.to_vec();
                    s.make_ascii_uppercase();
                    self.cached_seq = s;
                }
                Err(_) => {
                    self.cached_seq.clear();
                }
            }
            self.cached_chrom = chrom.to_string();
        }
        compute_gc_content(&self.cached_seq, pos as usize, window)
    }
}

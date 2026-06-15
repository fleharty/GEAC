use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};

use anyhow::{bail, Result};
use rust_htslib::faidx;
use tracing::info;

use crate::build_fusion_index::read_fai_sequences;
use crate::cli::ComputeUniquenessMapArgs;
use crate::kmer::kmer_iter;

pub fn compute_uniqueness_map(args: &ComputeUniquenessMapArgs) -> Result<()> {
    if args.min_k == 0 || args.min_k > args.max_k {
        bail!("--min-k must be ≥ 1 and ≤ --max-k");
    }
    if args.max_k > 31 {
        bail!("--max-k must be ≤ 31 (limit of the 2-bit k-mer encoder)");
    }

    let seq_list = read_fai_sequences(&args.fasta)?;
    let fai = faidx::Reader::from_path(&args.fasta)?;

    // Load optional output-region filter: chrom → sorted Vec<(start, end)> (0-based half-open).
    let region_filter: Option<HashMap<String, Vec<(usize, usize)>>> = args
        .regions
        .as_ref()
        .map(|p| load_bed_regions(p))
        .transpose()?;

    // Allocate per-chromosome result buffers; 0 = unresolved.
    // u8 is sufficient because max_k ≤ 31 and sentinel is max_k+1 ≤ 32.
    let mut results: Vec<(String, Vec<u8>)> = seq_list
        .iter()
        .map(|(name, len)| (name.clone(), vec![0u8; *len]))
        .collect();

    let total_bases: u64 = seq_list.iter().map(|(_, l)| *l as u64).sum();
    info!(
        n_sequences = seq_list.len(),
        total_bases,
        min_k = args.min_k,
        max_k = args.max_k,
        "computing per-locus minimum unique k"
    );

    for k in args.min_k..=args.max_k {
        info!(k, "starting count pass");

        // ── Count pass: genome-wide k-mer occurrence counts (capped at 2) ────────
        let mut genome_counts: HashMap<u64, u8> = HashMap::new();
        for (seq_name, seq_len) in &seq_list {
            let seq = match fai.fetch_seq(seq_name, 0, seq_len.saturating_sub(1)) {
                Ok(s) => s.to_vec(),
                Err(e) => {
                    tracing::warn!(seq = %seq_name, "failed to fetch sequence: {e}");
                    continue;
                }
            };
            for (kmer, _pos) in kmer_iter(&seq, k) {
                genome_counts
                    .entry(kmer)
                    .and_modify(|c| {
                        if *c < 2 {
                            *c += 1;
                        }
                    })
                    .or_insert(1);
            }
        }

        // Count how many singletons exist before the query pass.
        let n_unique = genome_counts.values().filter(|&&c| c == 1).count();
        info!(k, n_unique, "count pass done; starting query pass");

        // ── Query pass: assign min_k to newly-resolved positions ─────────────────
        let mut newly_resolved: u64 = 0;
        for ((chrom, buf), (seq_name, seq_len)) in results.iter_mut().zip(seq_list.iter()) {
            debug_assert_eq!(chrom, seq_name);
            let seq = match fai.fetch_seq(seq_name, 0, seq_len.saturating_sub(1)) {
                Ok(s) => s.to_vec(),
                Err(e) => {
                    tracing::warn!(seq = %seq_name, "failed to fetch sequence: {e}");
                    continue;
                }
            };
            for (kmer, pos) in kmer_iter(&seq, k) {
                if buf[pos] == 0 {
                    if let Some(&1) = genome_counts.get(&kmer) {
                        buf[pos] = k as u8;
                        newly_resolved += 1;
                    }
                }
            }
        }

        info!(k, newly_resolved, "query pass done");
        drop(genome_counts);

        // Check whether all positions are resolved (early exit).
        let unresolved: u64 = results
            .iter()
            .map(|(_, buf)| buf.iter().filter(|&&v| v == 0).count() as u64)
            .sum();
        info!(k, unresolved, "positions still unresolved after this k");
        if unresolved == 0 {
            info!(k, "all positions resolved; stopping early");
            break;
        }
    }

    // Assign sentinel (max_k + 1) to positions that were never resolved, including
    // N-run positions that kmer_iter skips and highly repetitive loci.
    let sentinel = (args.max_k + 1) as u8;
    for (_, buf) in results.iter_mut() {
        for v in buf.iter_mut() {
            if *v == 0 {
                *v = sentinel;
            }
        }
    }

    // ── Write output ─────────────────────────────────────────────────────────────
    let mut writer = BufWriter::new(
        File::create(&args.output)
            .map_err(|e| anyhow::anyhow!("cannot create output {}: {e}", args.output.display()))?,
    );

    for (chrom, buf) in &results {
        // Determine which positions to output.
        let include: Box<dyn Fn(usize) -> bool> =
            match region_filter.as_ref().and_then(|rf| rf.get(chrom.as_str())) {
                Some(intervals) => Box::new(move |pos: usize| {
                    intervals
                        .binary_search_by(|&(s, e)| {
                            if e <= pos {
                                std::cmp::Ordering::Less
                            } else if s > pos {
                                std::cmp::Ordering::Greater
                            } else {
                                std::cmp::Ordering::Equal
                            }
                        })
                        .is_ok()
                }),
                None => Box::new(|_pos: usize| true),
            };

        if args.no_merge {
            for (pos, &val) in buf.iter().enumerate() {
                if include(pos) {
                    writeln!(writer, "{chrom}\t{pos}\t{}\t{val}", pos + 1)?;
                }
            }
        } else {
            write_merged_bedgraph(&mut writer, chrom, buf, &include)?;
        }
    }

    writer.flush()?;
    info!(output = %args.output.display(), "done");
    Ok(())
}

/// Emit one bedgraph line per run of equal values, skipping positions excluded
/// by `include`. Runs that span an excluded position are split at that boundary.
fn write_merged_bedgraph(
    writer: &mut impl Write,
    chrom: &str,
    buf: &[u8],
    include: &dyn Fn(usize) -> bool,
) -> Result<()> {
    let mut run_start: Option<usize> = None;
    let mut run_val: u8 = 0;

    for (pos, &val) in buf.iter().enumerate() {
        if !include(pos) {
            // Flush any open run before the excluded position.
            if let Some(start) = run_start.take() {
                writeln!(writer, "{chrom}\t{start}\t{pos}\t{run_val}")?;
            }
            continue;
        }
        match run_start {
            None => {
                run_start = Some(pos);
                run_val = val;
            }
            Some(start) => {
                if val != run_val {
                    writeln!(writer, "{chrom}\t{start}\t{pos}\t{run_val}")?;
                    run_start = Some(pos);
                    run_val = val;
                }
            }
        }
    }
    // Flush the final run.
    if let Some(start) = run_start {
        writeln!(writer, "{chrom}\t{start}\t{}\t{run_val}", buf.len())?;
    }
    Ok(())
}

/// Parse a BED file into a map of chrom → sorted, merged intervals (0-based half-open).
fn load_bed_regions(path: &std::path::Path) -> Result<HashMap<String, Vec<(usize, usize)>>> {
    use std::io::BufRead;
    let file = std::fs::File::open(path)
        .map_err(|e| anyhow::anyhow!("cannot open regions BED {}: {e}", path.display()))?;
    let mut map: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
    for line in std::io::BufReader::new(file).lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let mut fields = line.split('\t');
        let chrom = fields.next().unwrap_or("").to_string();
        let start: usize = fields.next().and_then(|s| s.parse().ok()).unwrap_or(0);
        let end: usize = fields.next().and_then(|s| s.parse().ok()).unwrap_or(0);
        if !chrom.is_empty() && end > start {
            map.entry(chrom).or_default().push((start, end));
        }
    }
    // Sort and merge each chromosome's intervals.
    for intervals in map.values_mut() {
        intervals.sort_unstable();
        let mut merged: Vec<(usize, usize)> = Vec::with_capacity(intervals.len());
        for &(s, e) in intervals.iter() {
            if let Some(last) = merged.last_mut() {
                if s <= last.1 {
                    last.1 = last.1.max(e);
                    continue;
                }
            }
            merged.push((s, e));
        }
        *intervals = merged;
    }
    Ok(map)
}

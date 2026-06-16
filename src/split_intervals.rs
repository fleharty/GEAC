//! `geac split-intervals` — split a BED / Picard interval list into N balanced
//! shards for scatter-gather collection.
//!
//! Shards are a true partition of the input target positions: intervals are first
//! merged (so the input is disjoint), then walked in order and packed greedily so
//! each shard holds roughly `total_bases / scatter_count` bases. Large intervals
//! (e.g. whole chromosomes in a WGS target) are subdivided at base boundaries so
//! shards stay balanced. Because every base lands in exactly one shard, running
//! `geac collect --region <shard>` on each shard and concatenating the results is
//! identical to an unsharded run (see the interval clip in `bam::collect_alt_bases`).

use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};

use crate::cli::SplitIntervalsArgs;
use crate::targets::TargetIntervals;

/// Split the input interval list and write shard BED files. Returns the number of
/// shards written.
pub fn run(args: &SplitIntervalsArgs) -> Result<usize> {
    let intervals = TargetIntervals::load(&args.input)?;
    let merged: Vec<(String, u32, u32)> = intervals
        .merged_intervals()
        .into_iter()
        .map(|(c, s, e)| (c.to_string(), s, e))
        .collect();
    let total_bases = intervals.total_bases();
    anyhow::ensure!(
        total_bases > 0,
        "no target positions to split in {}",
        args.input.display()
    );

    let shards = split(&merged, total_bases, args.scatter_count);

    std::fs::create_dir_all(&args.output_dir).with_context(|| {
        format!(
            "failed to create output directory {}",
            args.output_dir.display()
        )
    })?;
    for (i, shard) in shards.iter().enumerate() {
        let path = args
            .output_dir
            .join(format!("{}.{:04}.bed", args.output_prefix, i));
        write_bed(&path, shard)?;
    }
    Ok(shards.len())
}

/// Partition `merged` (disjoint, sorted intervals) into at most `scatter_count`
/// shards balanced by base count, subdividing intervals as needed.
fn split(
    merged: &[(String, u32, u32)],
    total_bases: usize,
    scatter_count: usize,
) -> Vec<Vec<(String, u32, u32)>> {
    // Can't have more shards than bases; always at least one.
    let n = scatter_count.clamp(1, total_bases);
    let target = total_bases.div_ceil(n);

    let mut shards: Vec<Vec<(String, u32, u32)>> = Vec::with_capacity(n);
    let mut cur: Vec<(String, u32, u32)> = Vec::new();
    let mut cur_bases = 0usize;

    for (chrom, start, end) in merged {
        let mut s = *start;
        while s < *end {
            // Close the current shard once it is full, unless we are already
            // filling the last shard (which absorbs all remaining bases).
            if cur_bases >= target && shards.len() + 1 < n {
                shards.push(std::mem::take(&mut cur));
                cur_bases = 0;
            }
            // The last shard takes whole remaining intervals; earlier shards take
            // only up to their remaining capacity, subdividing the interval.
            let want = if shards.len() + 1 < n {
                target - cur_bases
            } else {
                usize::MAX
            };
            let take = ((*end - s) as usize).min(want);
            let new_s = s + take as u32;
            cur.push((chrom.clone(), s, new_s));
            cur_bases += take;
            s = new_s;
        }
    }
    if !cur.is_empty() || shards.is_empty() {
        shards.push(cur);
    }
    shards
}

fn write_bed(path: &Path, shard: &[(String, u32, u32)]) -> Result<()> {
    let mut f = std::fs::File::create(path)
        .with_context(|| format!("failed to create shard file {}", path.display()))?;
    for (chrom, start, end) in shard {
        writeln!(f, "{chrom}\t{start}\t{end}")
            .with_context(|| format!("failed writing to {}", path.display()))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iv(chrom: &str, start: u32, end: u32) -> (String, u32, u32) {
        (chrom.to_string(), start, end)
    }

    fn total(merged: &[(String, u32, u32)]) -> usize {
        merged.iter().map(|(_, s, e)| (e - s) as usize).sum()
    }

    /// Flatten shards back to a sorted list of base positions for partition checks.
    fn positions(shards: &[Vec<(String, u32, u32)>]) -> Vec<(String, u32)> {
        let mut out = Vec::new();
        for shard in shards {
            for (chrom, s, e) in shard {
                for p in *s..*e {
                    out.push((chrom.clone(), p));
                }
            }
        }
        out.sort();
        out
    }

    fn assert_partition(merged: &[(String, u32, u32)], shards: &[Vec<(String, u32, u32)>]) {
        let got = positions(shards);
        let mut dedup = got.clone();
        dedup.dedup();
        assert_eq!(got.len(), dedup.len(), "shards must be disjoint (no dupes)");

        let mut want = Vec::new();
        for (chrom, s, e) in merged {
            for p in *s..*e {
                want.push((chrom.clone(), p));
            }
        }
        want.sort();
        assert_eq!(got, want, "shards must cover exactly the input positions");
    }

    #[test]
    fn subdivides_one_large_interval_into_balanced_shards() {
        let merged = vec![iv("chr1", 0, 1000)];
        let shards = split(&merged, total(&merged), 4);
        assert_eq!(shards.len(), 4);
        assert_partition(&merged, &shards);
        for shard in &shards {
            let bases: usize = shard.iter().map(|(_, s, e)| (e - s) as usize).sum();
            assert_eq!(bases, 250, "1000 bases / 4 shards = 250 each");
        }
    }

    #[test]
    fn packs_many_small_intervals_across_chroms() {
        let merged = vec![
            iv("chr1", 0, 100),
            iv("chr1", 200, 300),
            iv("chr2", 0, 100),
            iv("chr2", 500, 600),
        ];
        let shards = split(&merged, total(&merged), 2);
        assert_eq!(shards.len(), 2);
        assert_partition(&merged, &shards);
    }

    #[test]
    fn scatter_count_one_is_single_shard() {
        let merged = vec![iv("chr1", 0, 100), iv("chr2", 0, 100)];
        let shards = split(&merged, total(&merged), 1);
        assert_eq!(shards.len(), 1);
        assert_partition(&merged, &shards);
    }

    #[test]
    fn scatter_count_exceeding_bases_is_clamped() {
        let merged = vec![iv("chr1", 0, 3)];
        let shards = split(&merged, total(&merged), 100);
        assert!(shards.len() <= 3, "cannot exceed base count");
        assert_partition(&merged, &shards);
    }

    #[test]
    fn uneven_split_covers_all_and_stays_within_count() {
        let merged = vec![iv("chr1", 0, 997)];
        let shards = split(&merged, total(&merged), 4);
        assert!(shards.len() <= 4);
        assert_partition(&merged, &shards);
    }
}

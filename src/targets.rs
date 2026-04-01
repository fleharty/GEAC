use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};

/// Sorted, merged target intervals for one chromosome.
/// All coordinates are stored as 0-based half-open [start, end).
type Intervals = Vec<(u32, u32)>;

/// A single target interval with its optional name, in original (unmerged) form.
/// Coordinates are 0-based half-open [start, end).
#[derive(Debug, Clone)]
pub struct NamedInterval {
    pub chrom: String,
    pub start: u32,
    pub end:   u32,
    pub name:  Option<String>,
}

/// Pre-loaded target region lookup built from a BED file or a Picard interval list.
///
/// Auto-detects format: if the file contains lines starting with `@` it is
/// treated as a Picard interval list (1-based, end-inclusive); otherwise it
/// is treated as BED (0-based, end-exclusive).
///
/// Overlapping intervals within each chromosome are merged on load so that
/// binary search is always correct.
///
/// The original (unmerged) intervals with their names are preserved separately
/// for per-interval summary output.
pub struct TargetIntervals {
    /// Merged intervals per chromosome, used for fast containment queries.
    by_chrom: HashMap<String, Intervals>,
    /// Original intervals in sorted (chrom, start) order, with names.
    named: Vec<NamedInterval>,
    /// Maps chrom → (start, end) range into `named` for fast lookup.
    named_range: HashMap<String, (usize, usize)>,
}

impl TargetIntervals {
    pub fn load(path: &Path) -> Result<Self> {
        let content = std::fs::read_to_string(path)
            .with_context(|| format!("failed to read targets file: {}", path.display()))?;
        Self::parse(&content)
    }

    fn parse(content: &str) -> Result<Self> {
        // Detect format by the presence of @ header lines (Picard interval list).
        let is_interval_list = content.lines().any(|l| l.starts_with('@'));

        let mut raw: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
        let mut named_raw: Vec<NamedInterval> = Vec::new();

        for line in content.lines() {
            if line.starts_with('@') || line.trim().is_empty() {
                continue;
            }

            let fields: Vec<&str> = line.splitn(6, '\t').collect();
            if fields.len() < 3 {
                continue;
            }

            let chrom = fields[0].to_string();

            let (start, end, name) = if is_interval_list {
                // Picard: 1-based, end-inclusive → convert to 0-based half-open
                // col 4 (index 4) is the interval name
                let s: u32 = fields[1].parse().with_context(|| {
                    format!("invalid start coordinate in interval list: '{}'", fields[1])
                })?;
                let e: u32 = fields[2].parse().with_context(|| {
                    format!("invalid end coordinate in interval list: '{}'", fields[2])
                })?;
                let name = fields.get(4).map(|n| n.trim().to_string()).filter(|n| !n.is_empty());
                (s.saturating_sub(1), e, name)
            } else {
                // BED: 0-based, end-exclusive — use as-is
                // col 3 (index 3) is the interval name
                let s: u32 = fields[1].parse().with_context(|| {
                    format!("invalid start coordinate in BED: '{}'", fields[1])
                })?;
                let e: u32 = fields[2].parse().with_context(|| {
                    format!("invalid end coordinate in BED: '{}'", fields[2])
                })?;
                let name = fields.get(3).map(|n| n.trim().to_string()).filter(|n| !n.is_empty());
                (s, e, name)
            };

            raw.entry(chrom.clone()).or_default().push((start, end));
            named_raw.push(NamedInterval { chrom, start, end, name });
        }

        // Sort and merge overlapping intervals per chromosome for efficient lookup.
        let by_chrom = raw
            .into_iter()
            .map(|(chrom, mut ivs)| {
                ivs.sort_unstable_by_key(|&(s, _)| s);
                (chrom, merge_intervals(ivs))
            })
            .collect();

        // Sort named intervals by (chrom, start) and build per-chrom range index.
        named_raw.sort_unstable_by(|a, b| {
            a.chrom.cmp(&b.chrom).then(a.start.cmp(&b.start))
        });
        let mut named_range: HashMap<String, (usize, usize)> = HashMap::new();
        let mut i = 0;
        while i < named_raw.len() {
            let chrom = named_raw[i].chrom.clone();
            let start = i;
            while i < named_raw.len() && named_raw[i].chrom == chrom {
                i += 1;
            }
            named_range.insert(chrom, (start, i));
        }

        Ok(Self { by_chrom, named: named_raw, named_range })
    }

    /// Returns `true` if the 0-based position `pos` falls within any target interval.
    pub fn contains(&self, chrom: &str, pos: i64) -> bool {
        let Some(intervals) = self.by_chrom.get(chrom) else {
            return false;
        };
        let pos = pos as u32;
        // Find the last interval whose start <= pos, then check if it covers pos.
        // This is correct because intervals are sorted and non-overlapping after merging.
        let idx = intervals.partition_point(|&(s, _)| s <= pos);
        if idx == 0 {
            return false;
        }
        let (_, end) = intervals[idx - 1];
        end > pos
    }

    pub fn n_targets(&self) -> usize {
        self.by_chrom.values().map(|v| v.len()).sum()
    }

    /// Iterate every (chrom, pos) pair covered by any target interval (0-based).
    /// Chromosomes are visited in sorted order; positions within each chromosome
    /// are visited in ascending order.
    pub fn for_each_position<F: FnMut(&str, i64)>(&self, mut f: F) {
        let mut chroms: Vec<&String> = self.by_chrom.keys().collect();
        chroms.sort();
        for chrom in chroms {
            for &(start, end) in &self.by_chrom[chrom] {
                for pos in start..end {
                    f(chrom, pos as i64);
                }
            }
        }
    }

    /// Total number of base positions covered by all intervals.
    pub fn total_bases(&self) -> usize {
        self.by_chrom.values()
            .flat_map(|ivs| ivs.iter())
            .map(|(s, e)| (e - s) as usize)
            .sum()
    }

    /// All named intervals in sorted (chrom, start) order.
    pub fn named_intervals(&self) -> &[NamedInterval] {
        &self.named
    }

    /// Number of named (original, unmerged) intervals.
    pub fn n_named_intervals(&self) -> usize {
        self.named.len()
    }

    /// Returns the index of the named interval containing `pos` (0-based), or None.
    /// For non-overlapping intervals this is exact. For overlapping intervals,
    /// returns the last interval whose start ≤ pos and end > pos.
    pub fn interval_index(&self, chrom: &str, pos: i64) -> Option<usize> {
        let &(range_start, range_end) = self.named_range.get(chrom)?;
        let slice = &self.named[range_start..range_end];
        let pos = pos as u32;
        let idx = slice.partition_point(|iv| iv.start <= pos);
        if idx == 0 {
            return None;
        }
        let candidate = &slice[idx - 1];
        if candidate.end > pos {
            Some(range_start + idx - 1)
        } else {
            None
        }
    }

    /// Iterate every named interval's (index, chrom, pos) triple.
    /// Each position within each named interval is visited exactly once.
    pub fn for_each_named_position<F: FnMut(usize, &str, i64)>(&self, mut f: F) {
        for (i, iv) in self.named.iter().enumerate() {
            for pos in iv.start..iv.end {
                f(i, &iv.chrom, pos as i64);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(s: &str) -> TargetIntervals {
        TargetIntervals::parse(s).expect("parse failed")
    }

    // ── BED ───────────────────────────────────────────────────────────────────

    #[test]
    fn bed_inside_interval() {
        let t = parse("chr1\t100\t200\n");
        assert!(t.contains("chr1", 100));
        assert!(t.contains("chr1", 150));
        assert!(t.contains("chr1", 199));
    }

    #[test]
    fn bed_outside_interval() {
        let t = parse("chr1\t100\t200\n");
        assert!(!t.contains("chr1", 99));
        assert!(!t.contains("chr1", 200)); // end is exclusive
        assert!(!t.contains("chr1", 300));
    }

    #[test]
    fn bed_unknown_chrom() {
        let t = parse("chr1\t100\t200\n");
        assert!(!t.contains("chr2", 150));
        assert!(!t.contains("chrX", 150));
    }

    #[test]
    fn bed_multiple_chroms_independent() {
        let t = parse("chr1\t100\t200\nchr2\t500\t600\n");
        assert!( t.contains("chr1", 150));
        assert!(!t.contains("chr1", 550));
        assert!( t.contains("chr2", 550));
        assert!(!t.contains("chr2", 150));
    }

    #[test]
    fn bed_overlapping_intervals_merged() {
        // [100,250) and [200,300) should merge to [100,300)
        let t = parse("chr1\t100\t250\nchr1\t200\t300\n");
        assert_eq!(t.n_targets(), 1);
        assert!(t.contains("chr1", 100));
        assert!(t.contains("chr1", 249));
        assert!(t.contains("chr1", 250)); // inside the merged interval
        assert!(t.contains("chr1", 299));
        assert!(!t.contains("chr1", 300));
    }

    #[test]
    fn bed_abutting_intervals_merged() {
        // [100,200) and [200,300) share an edge and should merge
        let t = parse("chr1\t100\t200\nchr1\t200\t300\n");
        assert_eq!(t.n_targets(), 1);
        assert!(t.contains("chr1", 199));
        assert!(t.contains("chr1", 200));
        assert!(t.contains("chr1", 299));
    }

    #[test]
    fn bed_non_overlapping_not_merged() {
        let t = parse("chr1\t100\t200\nchr1\t300\t400\n");
        assert_eq!(t.n_targets(), 2);
        assert!( t.contains("chr1", 150));
        assert!(!t.contains("chr1", 250));
        assert!( t.contains("chr1", 350));
    }

    #[test]
    fn bed_skips_comment_and_blank_lines() {
        let t = parse("# comment\nchr1\t100\t200\n\nchr1\t300\t400\n");
        assert!(t.contains("chr1", 150));
        assert!(t.contains("chr1", 350));
    }

    // ── Picard interval list ───────────────────────────────────────────────────

    #[test]
    fn picard_converts_to_zero_based() {
        // Picard: 1-based inclusive [101, 200] → 0-based half-open [100, 200)
        let t = parse("@HD\tVN:1.6\nchr1\t101\t200\t+\ttarget\n");
        assert!( t.contains("chr1", 100)); // pos 100 is inside
        assert!( t.contains("chr1", 199)); // pos 199 is the last inside
        assert!(!t.contains("chr1", 99));  // pos 99 is before
        assert!(!t.contains("chr1", 200)); // pos 200 is outside (half-open)
    }

    #[test]
    fn picard_position_1_converts_correctly() {
        // Picard start=1 → 0-based start=0
        let t = parse("@HD\tVN:1.6\nchr1\t1\t100\t+\ttarget\n");
        assert!(t.contains("chr1", 0));
        assert!(t.contains("chr1", 99));
        assert!(!t.contains("chr1", 100));
    }

    #[test]
    fn picard_header_lines_skipped() {
        let t = parse("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:248956422\nchr1\t101\t200\t+\ttarget\n");
        assert!(t.contains("chr1", 100));
        assert_eq!(t.n_targets(), 1);
    }
}

/// Merge a sorted (by start) vec of intervals, combining any that overlap or abut.
fn merge_intervals(sorted: Vec<(u32, u32)>) -> Vec<(u32, u32)> {
    let mut merged: Vec<(u32, u32)> = Vec::with_capacity(sorted.len());
    for (start, end) in sorted {
        if let Some(last) = merged.last_mut() {
            if start <= last.1 {
                last.1 = last.1.max(end);
                continue;
            }
        }
        merged.push((start, end));
    }
    merged
}

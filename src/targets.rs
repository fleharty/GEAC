use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};

/// Sorted, merged target intervals for one chromosome.
/// All coordinates are stored as 0-based half-open [start, end).
type Intervals = Vec<(u32, u32)>;

/// Pre-loaded target region lookup built from a BED file or a Picard interval list.
///
/// Auto-detects format: if the file contains lines starting with `@` it is
/// treated as a Picard interval list (1-based, end-inclusive); otherwise it
/// is treated as BED (0-based, end-exclusive).
///
/// Overlapping intervals within each chromosome are merged on load so that
/// binary search is always correct.
pub struct TargetIntervals {
    by_chrom: HashMap<String, Intervals>,
}

impl TargetIntervals {
    pub fn load(path: &Path) -> Result<Self> {
        let content = std::fs::read_to_string(path)
            .with_context(|| format!("failed to read targets file: {}", path.display()))?;

        // Detect format by the presence of @ header lines (Picard interval list).
        let is_interval_list = content.lines().any(|l| l.starts_with('@'));

        let mut raw: HashMap<String, Vec<(u32, u32)>> = HashMap::new();

        for line in content.lines() {
            if line.starts_with('@') || line.trim().is_empty() {
                continue;
            }

            let fields: Vec<&str> = line.splitn(6, '\t').collect();
            if fields.len() < 3 {
                continue;
            }

            let (start, end) = if is_interval_list {
                // Picard: 1-based, end-inclusive → convert to 0-based half-open
                let s: u32 = fields[1].parse().with_context(|| {
                    format!("invalid start coordinate in interval list: '{}'", fields[1])
                })?;
                let e: u32 = fields[2].parse().with_context(|| {
                    format!("invalid end coordinate in interval list: '{}'", fields[2])
                })?;
                (s.saturating_sub(1), e)
            } else {
                // BED: 0-based, end-exclusive — use as-is
                let s: u32 = fields[1].parse().with_context(|| {
                    format!("invalid start coordinate in BED: '{}'", fields[1])
                })?;
                let e: u32 = fields[2].parse().with_context(|| {
                    format!("invalid end coordinate in BED: '{}'", fields[2])
                })?;
                (s, e)
            };

            raw.entry(fields[0].to_string()).or_default().push((start, end));
        }

        // Sort and merge overlapping intervals per chromosome for efficient lookup.
        let by_chrom = raw
            .into_iter()
            .map(|(chrom, mut ivs)| {
                ivs.sort_unstable_by_key(|&(s, _)| s);
                (chrom, merge_intervals(ivs))
            })
            .collect();

        Ok(Self { by_chrom })
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

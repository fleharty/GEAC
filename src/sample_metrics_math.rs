//! Shared depth statistics for sample metrics, computed from a covered-depth
//! histogram.
//!
//! Both the unsharded `geac collect` path and the sharded `geac aggregate-metrics`
//! path funnel through [`depth_stats`], so a sample collected in one pass and the
//! same sample collected across interval shards and aggregated produce identical
//! metrics — medians included, because the histogram is a lossless representation
//! of the covered depths.

use std::collections::BTreeMap;

/// Derived per-sample depth statistics.
pub struct DepthStats {
    pub n_target_positions_covered: i64,
    pub mean_target_depth_covered: Option<f32>,
    pub mean_target_depth_all: Option<f32>,
    pub median_target_depth_covered: Option<f32>,
    pub median_target_depth_all: Option<f32>,
    pub pct_fragment_bases_on_target: Option<f32>,
}

/// Compute depth statistics from the covered-depth histogram (`depth -> count`)
/// plus the scalar sufficient statistics. Target positions absent from `hist` are
/// treated as depth 0 for the `_all` variants. All covered depths are > 0 (the
/// collector skips zero-depth columns), so the implicit zeros sort first.
pub fn depth_stats(
    hist: &BTreeMap<i32, i64>,
    n_target_positions: i64,
    total_fragment_bases: i64,
    on_target_fragment_bases: i64,
) -> DepthStats {
    let n_covered: i64 = hist.values().sum();
    let sum_depths: i64 = hist.iter().map(|(&d, &c)| i64::from(d) * c).sum();
    let zero_count = (n_target_positions - n_covered).max(0);

    let mean_target_depth_covered = if n_covered > 0 {
        Some(sum_depths as f32 / n_covered as f32)
    } else {
        None
    };
    let mean_target_depth_all = if n_target_positions > 0 {
        Some(sum_depths as f32 / n_target_positions as f32)
    } else {
        None
    };

    DepthStats {
        n_target_positions_covered: n_covered,
        mean_target_depth_covered,
        mean_target_depth_all,
        median_target_depth_covered: median_from_histogram(hist, 0),
        median_target_depth_all: median_from_histogram(hist, zero_count),
        pct_fragment_bases_on_target: if total_fragment_bases > 0 {
            Some(on_target_fragment_bases as f32 / total_fragment_bases as f32)
        } else {
            None
        },
    }
}

/// Median over `zero_count` implicit zeros followed by the depths in `hist`
/// (ascending). Averages the two middle values for an even total, matching the
/// prior sorted-slice implementation. None when there are no values.
fn median_from_histogram(hist: &BTreeMap<i32, i64>, zero_count: i64) -> Option<f32> {
    let covered: i64 = hist.values().sum();
    let total = zero_count + covered;
    if total == 0 {
        return None;
    }
    let mid = total / 2;
    if total % 2 == 1 {
        Some(value_at(hist, zero_count, mid))
    } else {
        Some((value_at(hist, zero_count, mid - 1) + value_at(hist, zero_count, mid)) / 2.0)
    }
}

/// Value at sorted index `idx` in the sequence of `zero_count` zeros followed by
/// the histogram depths in ascending order.
fn value_at(hist: &BTreeMap<i32, i64>, zero_count: i64, idx: i64) -> f32 {
    if idx < zero_count {
        return 0.0;
    }
    let mut remaining = idx - zero_count;
    for (&depth, &count) in hist {
        if remaining < count {
            return depth as f32;
        }
        remaining -= count;
    }
    // Unreachable for idx < total; fall back to the largest depth.
    hist.keys().next_back().map(|&d| d as f32).unwrap_or(0.0)
}

/// Build a covered-depth histogram from a slice of depths (all expected > 0).
pub fn histogram_from_depths(depths: &[i32]) -> BTreeMap<i32, i64> {
    let mut hist = BTreeMap::new();
    for &d in depths {
        *hist.entry(d).or_insert(0) += 1;
    }
    hist
}

#[cfg(test)]
mod tests {
    use super::*;

    fn brute_median(values: &mut [f32]) -> Option<f32> {
        if values.is_empty() {
            return None;
        }
        values.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let n = values.len();
        let mid = n / 2;
        Some(if n % 2 == 1 {
            values[mid]
        } else {
            (values[mid - 1] + values[mid]) / 2.0
        })
    }

    #[test]
    fn histogram_median_matches_bruteforce_covered() {
        let depths = vec![3, 1, 4, 1, 5, 9, 2, 6];
        let hist = histogram_from_depths(&depths);
        let stats = depth_stats(&hist, depths.len() as i64, 0, 0);
        let mut vals: Vec<f32> = depths.iter().map(|&d| d as f32).collect();
        assert_eq!(stats.median_target_depth_covered, brute_median(&mut vals));
    }

    #[test]
    fn histogram_median_all_includes_zeros() {
        let depths = vec![10, 20, 30];
        let hist = histogram_from_depths(&depths);
        // 3 covered + 2 zero-depth = 5 target positions.
        let stats = depth_stats(&hist, 5, 0, 0);
        let mut vals = vec![0.0f32, 0.0, 10.0, 20.0, 30.0];
        assert_eq!(stats.median_target_depth_all, brute_median(&mut vals));
        assert_eq!(stats.median_target_depth_covered, Some(20.0));
        assert_eq!(stats.n_target_positions_covered, 3);
    }

    #[test]
    fn empty_histogram_yields_none() {
        let hist = BTreeMap::new();
        let stats = depth_stats(&hist, 0, 0, 0);
        assert_eq!(stats.median_target_depth_covered, None);
        assert_eq!(stats.mean_target_depth_all, None);
    }
}

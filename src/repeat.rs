/// Per-locus repetitiveness metrics computed from a reference sequence window.
pub struct RepeatMetrics {
    /// Length of the longest homopolymer run overlapping this position (always >= 1).
    pub homopolymer_len: i32,
    /// Period of the shortest tandem repeat unit whose tract includes this position
    /// (1 = homopolymer, 2 = dinucleotide, …). 0 if no repeat found.
    pub str_period: i32,
    /// Total length (bp) of the STR tract for `str_period`. 0 if no repeat found.
    pub str_len: i32,
}

/// Compute repeat metrics for a single 0-based `pos` on a chromosome.
///
/// `seq`    — the full (or cached) chromosome sequence, upper-cased ASCII bytes.
/// `pos`    — 0-based position of the locus.
/// `window` — number of bases to scan on each side of `pos`.
///
/// **Homopolymer:** longest uninterrupted run of the same base that includes `pos`.
///
/// **STR:** for periods k = 1 … 6, finds the maximal tandem repeat tract of period k
/// that includes `pos` and has ≥ 2 consecutive copies.  Returns the smallest k
/// (and corresponding tract length) for which this holds.  Homopolymers (k = 1)
/// are subsumed: if `homopolymer_len >= 2`, `str_period` will be 1.
pub fn compute_repeat_metrics(seq: &[u8], pos: usize, window: usize) -> RepeatMetrics {
    if pos >= seq.len() {
        return RepeatMetrics { homopolymer_len: 1, str_period: 0, str_len: 0 };
    }

    let start    = pos.saturating_sub(window);
    let end      = (pos + window + 1).min(seq.len());
    let w        = &seq[start..end];
    let pos_in_w = pos - start;

    // ── Homopolymer ───────────────────────────────────────────────────────────
    let base = w[pos_in_w];
    let mut left = 0usize;
    for i in (0..pos_in_w).rev() {
        if w[i] != base { break; }
        left += 1;
    }
    let mut right = 0usize;
    for i in (pos_in_w + 1)..w.len() {
        if w[i] != base { break; }
        right += 1;
    }
    let homopolymer_len = (left + right + 1) as i32;

    // ── STR ───────────────────────────────────────────────────────────────────
    // For each period k (ascending), slide a window through the local sequence
    // and find the maximal tract of consecutive equal k-mers.  The first k that
    // produces a tract of ≥ 2 copies covering pos_in_w wins.
    const MAX_PERIOD: usize = 6;
    let mut str_period = 0i32;
    let mut str_len    = 0i32;

    'outer: for k in 1..=MAX_PERIOD {
        let mut i = 0;
        while i + 2 * k <= w.len() {
            let unit = &w[i..i + k];
            let mut copies = 1usize;
            while i + (copies + 1) * k <= w.len()
                && &w[i + copies * k..i + (copies + 1) * k] == unit
            {
                copies += 1;
            }
            if copies >= 2 {
                let tract_end = i + copies * k;
                if i <= pos_in_w && pos_in_w < tract_end {
                    str_period = k as i32;
                    str_len    = (copies * k) as i32;
                    break 'outer;
                }
            }
            i += 1;
        }
    }

    RepeatMetrics { homopolymer_len, str_period, str_len }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn metrics(seq: &str, pos: usize) -> RepeatMetrics {
        compute_repeat_metrics(seq.as_bytes(), pos, 10)
    }

    #[test]
    fn single_base_no_repeat() {
        // No homopolymer or repeat around the G
        let m = metrics("ACGTCGAT", 2);
        assert_eq!(m.homopolymer_len, 1);
        assert_eq!(m.str_period, 0);
        assert_eq!(m.str_len, 0);
    }

    #[test]
    fn homopolymer() {
        let m = metrics("ACTTTTGC", 3); // middle T in TTTT
        assert_eq!(m.homopolymer_len, 4);
        assert_eq!(m.str_period, 1);
        assert_eq!(m.str_len, 4);
    }

    #[test]
    fn dinucleotide_str() {
        let m = metrics("ACACACACGT", 4); // middle of ACACACAC
        assert_eq!(m.str_period, 2);
        assert_eq!(m.str_len, 8);
    }

    #[test]
    fn homopolymer_trumps_dinucleotide() {
        // AAAAAA is both a homopolymer (k=1) and a dinucleotide repeat (k=2 "AA").
        // k=1 should win as it's the smallest period.
        let m = metrics("TTAAAAAATT", 5);
        assert_eq!(m.str_period, 1);
    }
}

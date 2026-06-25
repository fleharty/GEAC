use std::collections::HashMap;

/// Rolling canonical k-mer iterator.
///
/// Yields `(canonical_kmer, start_pos)` where `start_pos` is the 0-based index
/// of the first base of the k-mer within the input slice. Non-ACGT bytes reset
/// the sliding window; positions within reset runs are skipped in the output.
/// Encoding: A=0, C=1, G=2, T=3 (2-bit); k must be 1..=31.
pub struct KmerIter<'a> {
    seq: &'a [u8],
    k: usize,
    mask: u64,
    fwd: u64,
    rev: u64,
    valid: usize,
    pos: usize,
}

impl<'a> Iterator for KmerIter<'a> {
    type Item = (u64, usize); // (canonical_kmer, start_pos_in_slice)
    fn next(&mut self) -> Option<(u64, usize)> {
        while self.pos < self.seq.len() {
            let b = self.seq[self.pos];
            self.pos += 1;
            let enc: u64 = match b | 0x20 {
                b'a' => 0,
                b'c' => 1,
                b'g' => 2,
                b't' => 3,
                _ => {
                    self.valid = 0;
                    self.fwd = 0;
                    self.rev = 0;
                    continue;
                }
            };
            self.fwd = ((self.fwd << 2) | enc) & self.mask;
            // RC: complement (3^enc) shifted into the high end; old LSBs drop on >>2.
            self.rev = (self.rev >> 2) | ((3u64 ^ enc) << (2 * (self.k - 1)));
            self.valid += 1;
            if self.valid >= self.k {
                // pos has already been advanced past the last base, so the k-mer
                // starts at pos - k.
                return Some((self.fwd.min(self.rev), self.pos - self.k));
            }
        }
        None
    }
}

pub fn kmer_iter(seq: &[u8], k: usize) -> KmerIter<'_> {
    assert!(k > 0 && k <= 31, "k must be 1..=31");
    KmerIter {
        seq,
        k,
        mask: (1u64 << (2 * k)).wrapping_sub(1),
        fwd: 0,
        rev: 0,
        valid: 0,
        pos: 0,
    }
}

/// Positions of bases that [`kmer_iter`] treats as window breaks — anything other than
/// A/C/G/T (case-insensitive). Each such position masks the `k`-mer windows covering it.
/// Mirrors the validity test inside [`KmerIter::next`], so the two stay in sync.
pub fn non_acgt_positions(seq: &[u8]) -> Vec<usize> {
    seq.iter()
        .enumerate()
        .filter(|(_, &b)| !matches!(b | 0x20, b'a' | b'c' | b'g' | b't'))
        .map(|(i, _)| i)
        .collect()
}

/// Render a per-k-mer-window layout track of length `n_windows`, one character per k-mer
/// start position (in the sequence's stored orientation):
///
/// - `N` — the window covers a non-ACGT base (a position from `n_base_positions`), so
///   [`kmer_iter`] emits no k-mer there: a *masked* window, not a sequence gap.
/// - `A` / `B` — the window's k-mer matched gene A / gene B (its start position appears
///   in `a_positions` / `b_positions`, which are the positions yielded by `kmer_iter`).
/// - `.` — a k-mer was emitted but matched neither gene.
///
/// `N` windows are painted first. Because `kmer_iter` never emits a k-mer for a masked
/// window, `a_positions`/`b_positions` can't overlap them, so painting order is safe.
/// This renders the *exact* track; both the `FL` BAM tag and `diagnose-fusion`'s `layout_5to3`
/// column then run [`rescue_layout_track`] on top to mark edit-distance-1 windows (lowercase
/// `a`/`b`, and capital from a uniquely-resolving single-`N`).
pub fn render_layout_track(
    n_windows: usize,
    a_positions: &[usize],
    b_positions: &[usize],
    n_base_positions: &[usize],
    k: usize,
) -> String {
    let mut track = vec![b'.'; n_windows];
    for &p in n_base_positions {
        let start = p.saturating_sub(k - 1);
        let end = p.min(n_windows.saturating_sub(1));
        for slot in track.iter_mut().take(end + 1).skip(start) {
            *slot = b'N';
        }
    }
    for &p in a_positions {
        if p < n_windows {
            track[p] = b'A';
        }
    }
    for &p in b_positions {
        if p < n_windows {
            track[p] = b'B';
        }
    }
    String::from_utf8(track).unwrap()
}

/// Canonical k-mer of an all-ACGT window of length exactly `k`. Returns `None` if the
/// window contains a non-ACGT base (then [`kmer_iter`] yields nothing). Used by the
/// edit-distance-1 rescue below to probe substitution variants against the gene index.
fn window_canonical(window: &[u8], k: usize) -> Option<u64> {
    debug_assert_eq!(window.len(), k);
    kmer_iter(window, k).next().map(|(kmer, _)| kmer)
}

/// Resolve an ambiguous split window — one whose 1-edit neighbors reach *both* genes — by
/// the flanking exact-match block, reading the pre-rescue track `anchors` (which holds only
/// exact `A`/`B`/`N`/`.`). Returns the lowercase letter when the nearest exact anchors on
/// both sides agree, or when only one side has an anchor ("inferred, not observed"). Returns
/// `None` when the flanks disagree — the window sits at the junction between an A-block and a
/// B-block, genuinely undeterminable — or when no exact anchor exists at all.
fn flanking_letter(anchors: &[u8], i: usize) -> Option<u8> {
    let is_ab = |c: u8| c == b'A' || c == b'B';
    let left = anchors[..i].iter().rev().find(|&&c| is_ab(c)).copied();
    let right = anchors[i + 1..].iter().find(|&&c| is_ab(c)).copied();
    match (left, right) {
        (Some(l), Some(r)) if l == r => Some(l.to_ascii_lowercase()),
        (Some(l), None) => Some(l.to_ascii_lowercase()),
        (None, Some(r)) => Some(r.to_ascii_lowercase()),
        _ => None,
    }
}

/// Edit-distance-1 rescue layer for an exact-match layout track produced by
/// [`render_layout_track`]. Upgrades only `.` and `N` windows in place; exact `A`/`B`
/// are never touched, preserving the invariant "capital = unambiguous match".
///
/// - `.` (a k-mer was emitted but matched neither gene): probe the `3·k` single-substitution
///   neighbors of the window's k-mer. Unique gene → lowercase `a`/`b` (a SNP-rescue
///   *disagrees* at a known base, so it stays lowercase even when unique). No neighbor
///   matches → left `.`.
/// - `N` (window masked by a non-ACGT base): rescued only when the window contains *exactly
///   one* non-ACGT base — its uncertainty is then at a single known position, so the window
///   is within edit distance 1 by construction. Substitute that base with A/C/G/T and probe.
///   Unique gene → CAPITAL `A`/`B` (the `N` does not disagree with the gene, it is merely
///   unknown, so a unique resolution is as good as exact). Multi-`N` / no match → `N`.
/// - Split windows (a neighbor reaches *both* genes) are resolved by [`flanking_letter`]:
///   the surrounding exact block decides the gene and the window is rendered lowercase
///   ("inferred, not observed"); a window between an A-block and a B-block (the junction)
///   has no agreeing flank and is left `.`/`N`. Splits require cross-gene shared k-mers,
///   which the index `--edit-distance-filter`/`--max-genome-copies` passes suppress, so they
///   are rare on a well-built index.
///
/// `seq` is the read sequence in the same orientation the track was rendered from; window `i`
/// is `seq[i..i+k]`. The caller guarantees `track.len() == seq.len() - k + 1`.
pub fn rescue_layout_track(
    track: &mut [u8],
    seq: &[u8],
    k: usize,
    kmer_to_gene: &HashMap<u64, u32>,
    ga: u32,
    gb: u32,
) {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let is_acgt = |b: u8| matches!(b | 0x20, b'a' | b'c' | b'g' | b't');

    // Snapshot the exact (pre-rescue) track so flanking-block resolution always reads the
    // original A/B anchors, never partially-rescued output — keeping the result independent
    // of the left-to-right mutation order.
    let anchors = track.to_vec();

    for i in 0..track.len() {
        match track[i] {
            b'.' => {
                // k-mer matched neither gene exactly; the window is all-ACGT.
                let mut buf = seq[i..i + k].to_vec();
                let (mut hit_a, mut hit_b) = (false, false);
                'scan: for p in 0..k {
                    let orig = buf[p];
                    for &b in &BASES {
                        if b == orig {
                            continue;
                        }
                        buf[p] = b;
                        if let Some(&g) =
                            window_canonical(&buf, k).and_then(|km| kmer_to_gene.get(&km))
                        {
                            if g == ga {
                                hit_a = true;
                            } else if g == gb {
                                hit_b = true;
                            }
                            if hit_a && hit_b {
                                buf[p] = orig;
                                break 'scan; // split — flanking block decides
                            }
                        }
                    }
                    buf[p] = orig;
                }
                track[i] = match (hit_a, hit_b) {
                    (true, false) => b'a',
                    (false, true) => b'b',
                    (true, true) => flanking_letter(&anchors, i).unwrap_or(b'.'),
                    (false, false) => b'.',
                };
            }
            b'N' => {
                // Rescue only a single-unknown-base window.
                let window = &seq[i..i + k];
                let mut n_idx = None;
                let mut n_count = 0usize;
                for (j, &b) in window.iter().enumerate() {
                    if !is_acgt(b) {
                        n_count += 1;
                        n_idx = Some(j);
                    }
                }
                if n_count != 1 {
                    continue;
                }
                let j = n_idx.unwrap();
                let mut buf = window.to_vec();
                let (mut hit_a, mut hit_b) = (false, false);
                for &b in &BASES {
                    buf[j] = b;
                    if let Some(&g) = window_canonical(&buf, k).and_then(|km| kmer_to_gene.get(&km))
                    {
                        if g == ga {
                            hit_a = true;
                        } else if g == gb {
                            hit_b = true;
                        }
                    }
                }
                track[i] = match (hit_a, hit_b) {
                    (true, false) => b'A',
                    (false, true) => b'B',
                    (true, true) => flanking_letter(&anchors, i).unwrap_or(b'N'),
                    (false, false) => b'N',
                };
            }
            _ => {}
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_sequence() {
        assert_eq!(kmer_iter(b"", 4).count(), 0);
    }

    #[test]
    fn shorter_than_k() {
        assert_eq!(kmer_iter(b"ACG", 4).count(), 0);
    }

    #[test]
    fn single_kmer() {
        let hits: Vec<(u64, usize)> = kmer_iter(b"ACGT", 4).collect();
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].1, 0); // starts at position 0
    }

    #[test]
    fn positions_are_correct() {
        // ACGTACGT with k=4: 5 k-mers starting at positions 0,1,2,3,4
        let hits: Vec<(u64, usize)> = kmer_iter(b"ACGTACGT", 4).collect();
        assert_eq!(hits.len(), 5);
        for (i, &(_, pos)) in hits.iter().enumerate() {
            assert_eq!(pos, i);
        }
    }

    #[test]
    fn n_resets_window() {
        // N in the middle breaks the window; we get k-mers before and after.
        let a = kmer_iter(b"AAAA", 4).count();
        let b = kmer_iter(b"AAAANAAAA", 4).count();
        // 1 before N, 1 after N
        assert_eq!(a, 1);
        assert_eq!(b, 2);
    }

    #[test]
    fn position_after_n_run() {
        // Sequence: AAAA NNN CCCC (k=4)
        // k-mer 0: AAAA at pos 0
        // After NNN, next k-mer CCCC at pos 7
        let hits: Vec<(u64, usize)> = kmer_iter(b"AAAANNNCCCC", 4).collect();
        assert_eq!(hits.len(), 2);
        assert_eq!(hits[0].1, 0);
        assert_eq!(hits[1].1, 7); // 4 A's + 3 N's = position 7
    }

    #[test]
    fn canonical_is_palindrome_idempotent() {
        // A palindromic k-mer (its own RC) should yield the same value for
        // both the forward and RC sequences.
        // ACGT has RC ACGT — it's a palindrome for k=4.
        let fwd: Vec<(u64, usize)> = kmer_iter(b"ACGT", 4).collect();
        let rev: Vec<(u64, usize)> = kmer_iter(b"ACGT", 4).collect();
        assert_eq!(fwd, rev);
    }

    #[test]
    fn canonical_is_symmetric() {
        // kmer_iter on a sequence and its RC should yield the same canonical k-mers
        // (in reverse order since RC reverses the sequence).
        let seq = b"AACCTG";
        let rc = b"CAGGTT";
        let fwd_kmers: Vec<u64> = kmer_iter(seq, 4).map(|(k, _)| k).collect();
        let mut rc_kmers: Vec<u64> = kmer_iter(rc, 4).map(|(k, _)| k).collect();
        rc_kmers.reverse();
        assert_eq!(fwd_kmers, rc_kmers);
    }

    #[test]
    fn non_acgt_positions_flags_only_non_acgt() {
        assert_eq!(non_acgt_positions(b"ACGTacgt"), Vec::<usize>::new());
        assert_eq!(non_acgt_positions(b"ACNGT"), vec![2]);
        assert_eq!(non_acgt_positions(b"NN"), vec![0, 1]);
    }

    #[test]
    fn render_layout_track_paints_n_a_b_and_dot() {
        // k=11, 40 windows. N at base 5 masks windows whose [s, s+k) covers it: 0..=5.
        let track = render_layout_track(
            40,
            &(6..=14).collect::<Vec<_>>(),
            &(25..=39).collect::<Vec<_>>(),
            &[5],
            11,
        );
        assert_eq!(&track[0..6], "NNNNNN");
        assert_eq!(&track[6..15], "AAAAAAAAA");
        assert_eq!(&track[15..25], "..........");
        assert_eq!(&track[25..40], "BBBBBBBBBBBBBBB");
    }

    #[test]
    fn render_layout_track_clean_read_has_no_n() {
        let track = render_layout_track(
            40,
            &(0..=14).collect::<Vec<_>>(),
            &(25..=39).collect::<Vec<_>>(),
            &[],
            11,
        );
        assert!(!track.contains('N'));
        assert_eq!(&track[0..15], "AAAAAAAAAAAAAAA");
    }

    /// Index gene 0 with every k-mer of `gene_seq`; return (kmer_to_gene, k).
    fn index_one_gene(gene_seq: &[u8], k: usize, gene: u32) -> HashMap<u64, u32> {
        let mut m = HashMap::new();
        for (kmer, _) in kmer_iter(gene_seq, k) {
            m.insert(kmer, gene);
        }
        m
    }

    // Mirror of fusions::fusion_layout_track for unit-testing the rescue layer in isolation.
    fn layout_with_rescue(seq: &[u8], k: usize, m: &HashMap<u64, u32>, ga: u32, gb: u32) -> String {
        let n_windows = seq.len() - k + 1;
        let mut a = Vec::new();
        let mut b = Vec::new();
        for (kmer, pos) in kmer_iter(seq, k) {
            if let Some(&g) = m.get(&kmer) {
                if g == ga {
                    a.push(pos);
                } else if g == gb {
                    b.push(pos);
                }
            }
        }
        let mut track =
            render_layout_track(n_windows, &a, &b, &non_acgt_positions(seq), k).into_bytes();
        rescue_layout_track(&mut track, seq, k, m, ga, gb);
        String::from_utf8(track).unwrap()
    }

    #[test]
    fn rescue_snp_bridges_dot_run_with_lowercase() {
        // Gene A = the read's reference; a single mismatch in the middle knocks out the
        // exact match for every window covering it, but each is edit-distance 1 away.
        let gene = b"ACGATCGTAGCTAGGTCATGCA"; // 22 bp, low internal repetition
        let k = 6;
        let m = index_one_gene(gene, k, 0);
        let mut read = gene.to_vec();
        let snp = 11; // middle base
        read[snp] = if read[snp] == b'C' { b'A' } else { b'C' };
        let track = layout_with_rescue(&read, k, &m, 0, 1);

        // Windows covering the SNP (starts snp-k+1 ..= snp) miss exactly but rescue to 'a';
        // all other windows are exact 'A'. No '.' or 'N' should remain.
        assert!(!track.contains('.'), "no unbridged gaps: {track}");
        assert!(!track.contains('N'), "{track}");
        let lo = snp + 1 - k;
        for (i, c) in track.chars().enumerate() {
            if (lo..=snp).contains(&i) {
                assert_eq!(c, 'a', "window {i} should be SNP-rescued: {track}");
            } else {
                assert_eq!(c, 'A', "window {i} should be exact: {track}");
            }
        }
    }

    #[test]
    fn rescue_single_n_resolves_to_capital() {
        // One N in the read: every covering window is masked ('N'); substituting the unknown
        // base resolves uniquely to gene 0, so all are promoted to capital 'A'.
        let gene = b"ACGATCGTAGCTAGGTCATGCA";
        let k = 6;
        let m = index_one_gene(gene, k, 0);
        let mut read = gene.to_vec();
        read[11] = b'N';
        let track = layout_with_rescue(&read, k, &m, 0, 1);
        assert!(
            !track.contains('N'),
            "single-N windows should resolve: {track}"
        );
        assert!(track.chars().all(|c| c == 'A'), "{track}");
    }

    #[test]
    fn rescue_n_split_with_no_flank_stays_n() {
        // Two gene k-mers differ only at position 0: NTTT -> ATTT (gene A) and CTTT (gene B).
        // The single-N window is a split, and this 1-window read has no flanking exact block
        // to break the tie, so it stays 'N'.
        let k = 4;
        let mut m = HashMap::new();
        m.insert(kmer_iter(b"ATTT", k).next().unwrap().0, 0u32);
        m.insert(kmer_iter(b"CTTT", k).next().unwrap().0, 1u32);
        let track = layout_with_rescue(b"NTTT", k, &m, 0, 1);
        assert_eq!(track, "N");
    }

    #[test]
    fn rescue_snp_split_with_no_flank_stays_dot() {
        // GTTT matches neither gene exactly; substituting position 0 reaches gene A (ATTT)
        // and gene B (CTTT). Split with no flanking exact block -> left '.', not guessed.
        let k = 4;
        let mut m = HashMap::new();
        m.insert(kmer_iter(b"ATTT", k).next().unwrap().0, 0u32);
        m.insert(kmer_iter(b"CTTT", k).next().unwrap().0, 1u32);
        let track = layout_with_rescue(b"GTTT", k, &m, 0, 1);
        assert_eq!(track, ".");
    }

    #[test]
    fn rescue_leaves_truly_unmatched_dot() {
        // A window edit-distance >= 2 from any gene k-mer stays '.'.
        let k = 4;
        let m = index_one_gene(b"ATTT", k, 0);
        let track = layout_with_rescue(b"GGGG", k, &m, 0, 1);
        assert_eq!(track, ".");
    }

    #[test]
    fn flanking_letter_resolves_when_blocks_agree() {
        // Split window at index 2 surrounded by exact A on both sides -> lowercase 'a'.
        assert_eq!(flanking_letter(b"AA.AA", 2), Some(b'a'));
        // Only one flank present (read end) -> use it.
        assert_eq!(flanking_letter(b"..AAA", 0), Some(b'a'));
        assert_eq!(flanking_letter(b"BBB..", 4), Some(b'b'));
    }

    #[test]
    fn flanking_letter_abstains_at_junction_or_no_anchor() {
        // A-block left, B-block right: the junction is undeterminable -> None.
        assert_eq!(flanking_letter(b"AA.BB", 2), None);
        // No exact anchor anywhere -> None.
        assert_eq!(flanking_letter(b".....", 2), None);
        assert_eq!(flanking_letter(b"NNNNN", 2), None);
    }
}

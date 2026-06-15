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
/// This is the single source of truth for both `diagnose-fusion`'s `layout_5to3` column
/// and the `FL` BAM tag written by `fusions --reads-output`.
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
}

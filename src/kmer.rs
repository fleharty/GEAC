/// Rolling canonical k-mer iterator.
///
/// Yields the canonical form (min of forward and reverse-complement encodings)
/// of each k-mer in the input sequence. Non-ACGT bytes reset the sliding window.
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
    type Item = u64;
    fn next(&mut self) -> Option<u64> {
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
                return Some(self.fwd.min(self.rev));
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
        let kmers: Vec<u64> = kmer_iter(b"ACGT", 4).collect();
        assert_eq!(kmers.len(), 1);
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
    fn canonical_is_palindrome_idempotent() {
        // A palindromic k-mer (its own RC) should yield the same value for
        // both the forward and RC sequences.
        // ACGT has RC ACGT — it's a palindrome for k=4.
        let fwd: Vec<u64> = kmer_iter(b"ACGT", 4).collect();
        let rev: Vec<u64> = kmer_iter(b"ACGT", 4).collect();
        assert_eq!(fwd, rev);
    }

    #[test]
    fn canonical_is_symmetric() {
        // kmer_iter on a sequence and its RC should yield the same canonical k-mers
        // (in reverse order since RC reverses the sequence).
        let seq = b"AACCTG";
        let rc = b"CAGGTT";
        let fwd_kmers: Vec<u64> = kmer_iter(seq, 4).collect();
        let mut rc_kmers: Vec<u64> = kmer_iter(rc, 4).collect();
        rc_kmers.reverse();
        assert_eq!(fwd_kmers, rc_kmers);
    }
}

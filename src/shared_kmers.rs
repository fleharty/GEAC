use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::Path;

use anyhow::{bail, Context, Result};
use duckdb::Connection;
use rayon::prelude::*;
use rust_htslib::faidx;
use tracing::info;

use crate::build_fusion_index::{
    hamming_neighbors_canonical, parse_gene_bodies, read_fai_sequences, reverse_complement,
    GeneBody,
};
use crate::cli::SharedKmersArgs;
use crate::kmer::kmer_iter;

/// Decode a 2-bit canonical k-mer hash back to a DNA string. Encoding matches
/// `kmer_iter`: A=0 C=1 G=2 T=3, most-significant 2 bits are the leftmost base.
/// The reported sequence is the canonical (lexicographically smaller of the
/// forward/reverse-complement) orientation, not necessarily the genomic strand.
fn decode_kmer(hash: u64, k: usize) -> String {
    let mut bases = vec![0u8; k];
    for (i, slot) in bases.iter_mut().enumerate() {
        let shift = 2 * (k - 1 - i);
        *slot = match (hash >> shift) & 3 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        };
    }
    // Safe: every byte is one of ACGT.
    String::from_utf8(bases).unwrap()
}

/// Per-base (substitution) Hamming distance between two 2-bit k-mer encodings.
fn base_hamming(a: u64, b: u64, k: usize) -> u32 {
    let mut d = 0;
    for i in 0..k {
        let shift = 2 * i;
        if (a >> shift) & 3 != (b >> shift) & 3 {
            d += 1;
        }
    }
    d
}

/// Substitution distance between two canonical k-mers, minimised over relative
/// orientation. Each canonical form represents a k-mer that could have come from
/// either strand, so the biologically meaningful distance is the smaller of
/// comparing `a` to `b` and to `b`'s reverse complement.
fn canonical_distance(a: u64, b: u64, k: usize) -> u32 {
    base_hamming(a, b, k).min(base_hamming(a, reverse_complement(b, k), k))
}

/// Extract canonical k-mers from every body of a gene, returning a map from the
/// k-mer hash to the (chrom, 0-based position) of its first occurrence. A gene
/// name can map to several bodies (e.g. the same symbol on multiple contigs);
/// all are unioned, mirroring how build-fusion-index treats a gene.
fn gene_body_kmers(
    fai: &faidx::Reader,
    bodies: &[&GeneBody],
    k: usize,
) -> Result<HashMap<u64, (String, u32)>> {
    let mut kmers: HashMap<u64, (String, u32)> = HashMap::new();
    for body in bodies {
        if body.end <= body.start {
            continue;
        }
        // fetch_seq takes an inclusive end; body.end is half-open exclusive.
        let seq = fai
            .fetch_seq(&body.chrom, body.start, body.end - 1)
            .with_context(|| {
                format!(
                    "failed to fetch {}:{}-{} from FASTA",
                    body.chrom, body.start, body.end
                )
            })?
            .to_vec();
        for (kmer, offset) in kmer_iter(&seq, k) {
            // Record the first occurrence only; offset is the 0-based start within
            // the body slice, so body.start + offset is the chromosome coordinate.
            kmers
                .entry(kmer)
                .or_insert_with(|| (body.chrom.clone(), (body.start + offset) as u32));
        }
    }
    Ok(kmers)
}

/// Count genome-wide canonical-k-mer occurrences, but only for k-mers in `watch`.
/// Scans the full FASTA in parallel (one task per sequence, each opening its own
/// faidx reader since it is not Send) and keeps memory to O(|watch|).
fn genome_watch_counts(fasta: &Path, watch: &HashSet<u64>, k: usize) -> Result<HashMap<u64, u32>> {
    let mut seq_list = read_fai_sequences(fasta)?;
    // Longest sequences first so big chromosomes are handed out before threads go
    // idle on small alt contigs (longest-job-first scheduling).
    seq_list.sort_unstable_by_key(|b| std::cmp::Reverse(b.1));

    let per_seq: Vec<HashMap<u64, u32>> = seq_list
        .par_iter()
        .filter_map(|(name, len)| {
            let fai = match faidx::Reader::from_path(fasta) {
                Ok(r) => r,
                Err(e) => {
                    tracing::warn!(seq = %name, "failed to open FASTA in worker: {e}");
                    return None;
                }
            };
            let seq = match fai.fetch_seq(name, 0, len.saturating_sub(1)) {
                Ok(s) => s.to_vec(),
                Err(e) => {
                    tracing::warn!(seq = %name, "failed to fetch sequence: {e}");
                    return None;
                }
            };
            let mut local: HashMap<u64, u32> = HashMap::new();
            for (kmer, _pos) in kmer_iter(&seq, k) {
                if watch.contains(&kmer) {
                    *local.entry(kmer).or_insert(0) += 1;
                }
            }
            info!(seq = %name, "sequence scanned");
            Some(local)
        })
        .collect();

    let mut counts: HashMap<u64, u32> = HashMap::new();
    for local in per_seq {
        for (kmer, c) in local {
            *counts.entry(kmer).or_insert(0) += c;
        }
    }
    Ok(counts)
}

/// Load the set of canonical index k-mers assigned to `gene` from a built fusion
/// index, validating that the index k-mer size matches `k`. These are the k-mers
/// the fusion caller can actually vote on (a filtered subset of the gene body).
fn load_index_gene_kmers(index_path: &Path, gene: &str, k: usize) -> Result<HashSet<u64>> {
    let conn = Connection::open(index_path)
        .with_context(|| format!("failed to open fusion index: {}", index_path.display()))?;

    let index_k: Option<String> = conn
        .query_row("SELECT value FROM meta WHERE key = 'kmer_size'", [], |r| {
            r.get(0)
        })
        .ok();
    if let Some(s) = index_k {
        let ik: usize = s
            .parse()
            .with_context(|| format!("index meta kmer_size is not an integer: {s:?}"))?;
        anyhow::ensure!(
            ik == k,
            "--kmer-size {k} does not match the index (built with k={ik}); \
             pass --kmer-size {ik}"
        );
    }

    let mut stmt = conn.prepare(
        "SELECT k.kmer_hash FROM kmers k \
         JOIN genes g ON k.gene_index = g.gene_index \
         WHERE g.gene_name = ?",
    )?;
    let mut set = HashSet::new();
    let mut rows = stmt.query([gene])?;
    while let Some(row) = rows.next()? {
        let h: i64 = row.get(0)?;
        set.insert(h as u64);
    }
    Ok(set)
}

/// One A-k-mer ↔ B-k-mer match.
struct Match {
    a: u64,
    b: u64,
    dist: u32,
}

pub fn shared_kmers(args: &SharedKmersArgs) -> Result<()> {
    let k = args.kmer_size as usize;
    anyhow::ensure!((1..=31).contains(&k), "kmer-size must be between 1 and 31");

    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads as usize)
            .build_global()
            .unwrap_or_else(|e| tracing::warn!("could not set rayon thread count: {e}"));
    }

    let genes = parse_gene_bodies(&args.gene_annotation)?;
    let bodies_a: Vec<&GeneBody> = genes
        .iter()
        .filter(|g| g.gene_name == args.gene_a)
        .collect();
    let bodies_b: Vec<&GeneBody> = genes
        .iter()
        .filter(|g| g.gene_name == args.gene_b)
        .collect();
    if bodies_a.is_empty() {
        bail!("gene '{}' not found in annotation", args.gene_a);
    }
    if bodies_b.is_empty() {
        bail!("gene '{}' not found in annotation", args.gene_b);
    }

    let fai = faidx::Reader::from_path(&args.fasta)
        .with_context(|| format!("failed to open FASTA: {}", args.fasta.display()))?;

    let kmers_a = gene_body_kmers(&fai, &bodies_a, k)?;
    let kmers_b = gene_body_kmers(&fai, &bodies_b, k)?;

    // Optional: the actual index k-mer sets the fusion caller can vote on, used to
    // flag which matched k-mers could really drive a false call.
    let (index_a, index_b) = match &args.index {
        Some(p) => {
            let ia = load_index_gene_kmers(p, &args.gene_a, k)?;
            let ib = load_index_gene_kmers(p, &args.gene_b, k)?;
            info!(
                index = %p.display(),
                index_kmers_a = ia.len(),
                index_kmers_b = ib.len(),
                "loaded index k-mer sets for both genes"
            );
            (Some(ia), Some(ib))
        }
        None => (None, None),
    };

    // Match each gene-A k-mer to gene-B k-mers within `edit_distance` substitutions.
    // For each A k-mer the candidate set is itself (distance 0) plus its canonical
    // Hamming neighbours up to N; any candidate present in gene B is a match. At
    // N=0 this reduces to the exact shared-k-mer intersection.
    let n = args.edit_distance;
    let mut matches: Vec<Match> = Vec::new();
    for &a in kmers_a.keys() {
        if kmers_b.contains_key(&a) {
            matches.push(Match { a, b: a, dist: 0 });
        }
        if n > 0 {
            for b in hamming_neighbors_canonical(a, k, n) {
                if kmers_b.contains_key(&b) {
                    matches.push(Match {
                        a,
                        b,
                        dist: canonical_distance(a, b, k),
                    });
                }
            }
        }
    }
    // Deterministic order: by A hash, then B hash.
    matches.sort_unstable_by_key(|x| (x.a, x.b));

    // Optional genome scan: genome-wide copy number of each matched k-mer (both
    // the A side and the B side).
    let counts: HashMap<u64, u32> = if args.check_reference {
        let mut watch: HashSet<u64> = HashSet::new();
        for m in &matches {
            watch.insert(m.a);
            watch.insert(m.b);
        }
        info!(
            n_matches = matches.len(),
            watch_set = watch.len(),
            threads = rayon::current_num_threads(),
            "scanning reference for matched-k-mer copy numbers..."
        );
        genome_watch_counts(&args.fasta, &watch, k)?
    } else {
        HashMap::new()
    };

    let mut out: Box<dyn Write> = match &args.output {
        Some(p) => Box::new(std::io::BufWriter::new(
            std::fs::File::create(p)
                .with_context(|| format!("cannot create output: {}", p.display()))?,
        )),
        None => Box::new(std::io::BufWriter::new(std::io::stdout())),
    };

    // Header — optional columns appear only when their flag/input is present.
    let mut header =
        String::from("kmer_seq_a\tchrom_a\tpos_a\tkmer_seq_b\tchrom_b\tpos_b\tab_dist");
    if args.check_reference {
        header.push_str("\tref_copies_a\tref_copies_b");
    }
    if args.index.is_some() {
        header.push_str("\tin_index_a\tin_index_b");
    }
    writeln!(out, "{header}")?;

    for m in &matches {
        let (chrom_a, pos_a) = &kmers_a[&m.a];
        let (chrom_b, pos_b) = &kmers_b[&m.b];
        write!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            decode_kmer(m.a, k),
            chrom_a,
            pos_a,
            decode_kmer(m.b, k),
            chrom_b,
            pos_b,
            m.dist
        )?;
        if args.check_reference {
            write!(
                out,
                "\t{}\t{}",
                counts.get(&m.a).copied().unwrap_or(0),
                counts.get(&m.b).copied().unwrap_or(0)
            )?;
        }
        if let (Some(ia), Some(ib)) = (&index_a, &index_b) {
            write!(out, "\t{}\t{}", ia.contains(&m.a), ib.contains(&m.b))?;
        }
        writeln!(out)?;
    }
    out.flush()?;

    let n_exact = matches.iter().filter(|m| m.dist == 0).count();
    info!(
        gene_a = %args.gene_a,
        gene_b = %args.gene_b,
        kmers_a = kmers_a.len(),
        kmers_b = kmers_b.len(),
        edit_distance = n,
        matches = matches.len(),
        exact_matches = n_exact,
        "shared k-mer scan complete"
    );

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decode_kmer_roundtrips_known_values() {
        assert_eq!(decode_kmer(0, 4), "AAAA");
        assert_eq!(decode_kmer(255, 4), "TTTT");
        // 0b00_01_10_11 = ACGT
        assert_eq!(decode_kmer(0b00_01_10_11, 4), "ACGT");
    }

    #[test]
    fn canonical_distance_zero_for_identical() {
        // ACGT (k=4) vs itself.
        assert_eq!(canonical_distance(0b00_01_10_11, 0b00_01_10_11, 4), 0);
    }

    #[test]
    fn canonical_distance_counts_single_substitution() {
        // AAAA (0) vs AAAC (1): differ in one base. RC(AAAC)=GTTT, far from AAAA,
        // so the minimum-orientation distance is the direct one = 1.
        assert_eq!(canonical_distance(0b00_00_00_00, 0b00_00_00_01, 4), 1);
    }
}

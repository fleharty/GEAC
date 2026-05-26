use anyhow::{Context, Result};
use duckdb::Connection;
use tracing::info;

use crate::cli::LookupKmerArgs;

/// Encode a k-mer string into a (forward, reverse_complement) pair of u64s.
/// Same 2-bit scheme as `kmer.rs`: A=0, C=1, G=2, T=3.
fn encode_kmer(seq: &str, k: usize) -> Result<(u64, u64)> {
    let mask = (1u64 << (2 * k)).wrapping_sub(1);
    let mut fwd: u64 = 0;
    let mut rev: u64 = 0;
    for (i, b) in seq.bytes().enumerate() {
        let enc: u64 = match b | 0x20 {
            b'a' => 0,
            b'c' => 1,
            b'g' => 2,
            b't' => 3,
            _ => anyhow::bail!(
                "invalid base '{}' at position {} in k-mer '{}'",
                b as char,
                i,
                seq
            ),
        };
        fwd = ((fwd << 2) | enc) & mask;
        rev = (rev >> 2) | ((3u64 ^ enc) << (2 * (k - 1)));
    }
    Ok((fwd, rev))
}

pub fn lookup_kmer(args: &LookupKmerArgs) -> Result<()> {
    let kmer_seq = args.kmer.trim();
    let k = kmer_seq.len();
    anyhow::ensure!(k >= 1 && k <= 31, "k-mer length must be 1–31 (got {})", k);

    let (fwd, rev) = encode_kmer(kmer_seq, k)?;
    let canonical: u64 = fwd.min(rev);
    let canonical_signed: i64 = canonical as i64;

    info!(
        kmer = kmer_seq,
        k,
        canonical_u64 = canonical,
        canonical_signed = canonical_signed,
        is_palindrome = (fwd == rev),
        "computed canonical k-mer hash"
    );

    let conn = Connection::open(&args.index)
        .with_context(|| format!("failed to open fusion index: {}", args.index.display()))?;

    // gene_name and chrom from genes table; position (if any) from kmer_positions.
    let mut stmt = conn.prepare(
        "SELECT g.gene_name, g.chrom, kp.chrom AS pos_chrom, kp.pos
         FROM kmers k
         LEFT JOIN genes g ON k.gene_index = g.gene_index
         LEFT JOIN kmer_positions kp ON kp.kmer_hash = k.kmer_hash
         WHERE k.kmer_hash = ?",
    )?;
    let mut rows = stmt.query([canonical_signed])?;

    println!("kmer\t{}", kmer_seq);
    println!("k\t{}", k);
    println!("canonical_hash_u64\t{}", canonical);
    println!("canonical_hash_signed\t{}", canonical_signed);
    println!();
    println!("gene_name\tchrom\tpos_0based");

    let mut n_rows: usize = 0;
    while let Some(row) = rows.next()? {
        n_rows += 1;
        let gene_name: Option<String> = row.get(0).ok();
        let _gene_chrom: Option<String> = row.get(1).ok();
        let pos_chrom: Option<String> = row.get(2).ok();
        let pos: Option<i32> = row.get(3).ok();
        println!(
            "{}\t{}\t{}",
            gene_name.as_deref().unwrap_or("."),
            pos_chrom.as_deref().unwrap_or("."),
            pos.map(|p| p.to_string()).unwrap_or_else(|| ".".to_string()),
        );
    }

    if n_rows == 0 {
        println!(".\t.\t.");
        info!("k-mer not found in index");
    } else {
        info!(n_rows, "k-mer found");
    }

    Ok(())
}

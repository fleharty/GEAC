use anyhow::{Context, Result};
use duckdb::Connection;

use crate::cli::ScanReadArgs;
use crate::kmer::kmer_iter;

pub fn scan_read(args: &ScanReadArgs) -> Result<()> {
    let k = args.kmer_size as usize;
    anyhow::ensure!(k >= 1 && k <= 31, "--kmer-size must be 1–31");

    let conn = Connection::open(&args.index)
        .with_context(|| format!("failed to open fusion index: {}", args.index.display()))?;

    let mut stmt = conn.prepare(
        "SELECT g.gene_name, kp.chrom, kp.pos
         FROM kmers k
         LEFT JOIN genes g ON k.gene_index = g.gene_index
         LEFT JOIN kmer_positions kp ON kp.kmer_hash = k.kmer_hash
         WHERE k.kmer_hash = ?",
    )?;

    let seq = args.read.as_bytes();
    let total_kmers = seq.len().saturating_sub(k - 1);

    println!("offset\tkmer\tgene\tchrom\tpos_0based");

    let mut n_hits: usize = 0;
    for (canonical, offset) in kmer_iter(seq, k) {
        let canonical_signed = canonical as i64;
        let kmer_str = &args.read[offset..offset + k];

        let mut rows = stmt.query([canonical_signed])?;
        if let Some(row) = rows.next()? {
            n_hits += 1;
            let gene: Option<String> = row.get(0).ok();
            let chrom: Option<String> = row.get(1).ok();
            let pos: Option<i32> = row.get(2).ok();
            println!(
                "{}\t{}\t{}\t{}\t{}",
                offset,
                kmer_str,
                gene.as_deref().unwrap_or("."),
                chrom.as_deref().unwrap_or("."),
                pos.map(|p| p.to_string()).unwrap_or_else(|| ".".to_string()),
            );
        }
    }

    eprintln!(
        "scan complete: {}/{} k-mers matched (k={})",
        n_hits, total_kmers, k
    );

    Ok(())
}

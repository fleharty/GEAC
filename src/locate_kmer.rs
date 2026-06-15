use std::io::{self, BufRead, BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use rust_htslib::faidx;
use tracing::info;

use crate::cli::LocateKmerArgs;
use crate::gene_annotations::GeneAnnotations;

// ─── K-mer encoding ───────────────────────────────────────────────────────────

/// Encode a k-mer string into a (forward, reverse_complement) pair of u64s.
/// Uses the same 2-bit scheme as kmer_iter: A=0, C=1, G=2, T=3.
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

// ─── FAI helper ──────────────────────────────────────────────────────────────

fn read_fai_sequences(fasta: &Path) -> Result<Vec<(String, usize)>> {
    let mut fai_str = fasta.to_string_lossy().into_owned();
    fai_str.push_str(".fai");
    let fai_path = std::path::PathBuf::from(fai_str);

    let file = std::fs::File::open(&fai_path).with_context(|| {
        format!(
            "cannot open FASTA index {}; run 'samtools faidx' first",
            fai_path.display()
        )
    })?;

    let mut seqs = Vec::new();
    for line in io::BufReader::new(file).lines() {
        let line = line?;
        let mut fields = line.split('\t');
        let name = fields.next().unwrap_or("").to_string();
        let len: usize = fields.next().and_then(|s| s.parse().ok()).unwrap_or(0);
        if !name.is_empty() && len > 0 {
            seqs.push((name, len));
        }
    }
    Ok(seqs)
}

// ─── Main entry point ─────────────────────────────────────────────────────────

pub fn locate_kmer(args: &LocateKmerArgs) -> Result<()> {
    let kmer_seq = args.kmer.trim();
    let k = kmer_seq.len();
    anyhow::ensure!(
        (1..=31).contains(&k),
        "k-mer length must be 1–31 (got {})",
        k
    );

    let (target_fwd, target_rc) = encode_kmer(kmer_seq, k)?;
    let is_palindrome = target_fwd == target_rc;

    info!(
        kmer = kmer_seq,
        k, is_palindrome, "searching reference for k-mer occurrences..."
    );

    let gene_annots = args
        .gene_annotations
        .as_ref()
        .map(|p| {
            info!(gene_annotations = %p.display(), "loading gene annotations");
            GeneAnnotations::load(p)
        })
        .transpose()?;

    let fai = faidx::Reader::from_path(&args.fasta)
        .with_context(|| format!("failed to open FASTA: {}", args.fasta.display()))?;

    let seq_list = read_fai_sequences(&args.fasta)?;

    // Output: TSV with header, written to --output or stdout.
    let stdout = io::stdout();
    let mut w: Box<dyn Write> = match &args.output {
        Some(path) => {
            let file = std::fs::File::create(path)
                .with_context(|| format!("failed to create output: {}", path.display()))?;
            Box::new(BufWriter::new(file))
        }
        None => Box::new(BufWriter::new(stdout.lock())),
    };

    if gene_annots.is_some() {
        writeln!(w, "chrom\tstart\tend\tstrand\tgene\tregion")?;
    } else {
        writeln!(w, "chrom\tstart\tend\tstrand")?;
    }

    let mask = (1u64 << (2 * k)).wrapping_sub(1);
    let mut total_hits: u64 = 0;

    for (chrom_name, chrom_len) in &seq_list {
        let chrom_seq = match fai.fetch_seq(chrom_name, 0, chrom_len.saturating_sub(1)) {
            Ok(s) => s.to_vec(),
            Err(e) => {
                tracing::warn!(chrom = %chrom_name, "failed to fetch sequence: {e}");
                continue;
            }
        };

        let mut fwd: u64 = 0;
        let mut valid = 0usize;

        for (pos, &b) in chrom_seq.iter().enumerate() {
            let enc: u64 = match b | 0x20 {
                b'a' => 0,
                b'c' => 1,
                b'g' => 2,
                b't' => 3,
                _ => {
                    fwd = 0;
                    valid = 0;
                    continue;
                }
            };
            fwd = ((fwd << 2) | enc) & mask;
            valid += 1;

            if valid >= k {
                let start = pos + 1 - k; // 0-based start of this k-mer
                let end = start + k;

                let is_fwd_hit = fwd == target_fwd;
                let is_rc_hit = fwd == target_rc && !is_palindrome;

                if is_fwd_hit || is_rc_hit {
                    let annot_suffix = if let Some(ref ga) = gene_annots {
                        let ann = ga.get(chrom_name, start as i64);
                        let (gene, region) = match ann {
                            Some(a) => {
                                let region = match a.feature_type.as_deref() {
                                    Some(ft) => ft.to_string(),
                                    None => "intron".to_string(),
                                };
                                (a.gene, region)
                            }
                            None => (".".to_string(), "intergenic".to_string()),
                        };
                        format!("\t{}\t{}", gene, region)
                    } else {
                        String::new()
                    };

                    if is_fwd_hit {
                        writeln!(w, "{}\t{}\t{}\t+{}", chrom_name, start, end, annot_suffix)?;
                        total_hits += 1;
                    }
                    // When target_fwd == target_rc (palindrome) the k-mer is its own
                    // reverse complement; don't emit a duplicate "-" line.
                    if is_rc_hit {
                        writeln!(w, "{}\t{}\t{}\t-{}", chrom_name, start, end, annot_suffix)?;
                        total_hits += 1;
                    }
                }
            }
        }
    }

    info!(total_hits, kmer = kmer_seq, "search complete");
    Ok(())
}

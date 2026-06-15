use std::collections::HashSet;
use std::path::Path;

use anyhow::{Context, Result};
use rust_htslib::{bam, faidx};
use tracing::warn;

/// Caches one chromosome sequence at a time to avoid per-base faidx seeks.
/// When the chromosome changes, the new sequence is fetched and the old one dropped.
/// Contigs present in the BAM header but absent from the FASTA index are skipped
/// gracefully (returning 'N') rather than crashing via htslib.
pub(crate) struct RefCache {
    fai: faidx::Reader,
    /// Sequence names present in the FASTA index.
    fasta_seqs: HashSet<String>,
    current_tid: Option<usize>,
    current_chrom: String,
    /// Upper-cased sequence bytes for the cached chromosome (0-based)
    current_seq: Vec<u8>,
}

impl RefCache {
    pub(crate) fn new(reference: &Path) -> Result<Self> {
        let fai = faidx::Reader::from_path(reference)
            .with_context(|| format!("failed to open reference FASTA: {}", reference.display()))?;
        let fasta_seqs = fai
            .seq_names()
            .context("failed to read sequence names from FASTA index")?
            .into_iter()
            .collect();
        Ok(Self {
            fai,
            fasta_seqs,
            current_tid: None,
            current_chrom: String::new(),
            current_seq: Vec::new(),
        })
    }

    /// Returns (chrom_name, ref_base) for the given tid and 0-based position.
    /// Loads the chromosome sequence on first access or when the chromosome changes.
    /// `targets` is a pre-built slice of (chrom_name, chrom_len) indexed by tid.
    /// Returns `'N'` for contigs absent from the FASTA (e.g. decoy sequences like
    /// `hs37d5`) so the caller skips them without crashing.
    pub(crate) fn get(
        &mut self,
        targets: &[(String, usize)],
        tid: usize,
        pos: usize,
    ) -> Result<(String, char)> {
        if self.current_tid != Some(tid) {
            let (chrom, chrom_len) = targets
                .get(tid)
                .with_context(|| format!("tid {tid} not found in BAM header"))?;

            if !self.fasta_seqs.contains(chrom.as_str()) {
                warn!(
                    chrom,
                    "contig is in BAM header but not in FASTA index; skipping"
                );
                // Mark as loaded with an empty sequence so the warning fires once.
                self.current_tid = Some(tid);
                self.current_chrom = chrom.clone();
                self.current_seq = Vec::new();
                return Ok((chrom.clone(), 'N'));
            }

            self.current_seq = self
                .fai
                .fetch_seq(chrom, 0, chrom_len.saturating_sub(1))
                .with_context(|| format!("failed to fetch sequence for {chrom}"))?
                .to_vec();

            self.current_seq.make_ascii_uppercase();

            self.current_tid = Some(tid);
            self.current_chrom = chrom.clone();
        }

        let base = self.current_seq.get(pos).map(|&b| b as char).unwrap_or('N');

        Ok((self.current_chrom.clone(), base))
    }

    /// Returns the cached sequence for the current chromosome.
    /// Must be called after `get()` has loaded the chromosome.
    pub(crate) fn current_seq(&self) -> &[u8] {
        &self.current_seq
    }
}

/// Extract the SM (sample name) field from the first @RG line in the BAM header.
/// Returns an error if no @RG line exists or none has an SM tag.
pub fn read_group_sample_id(header: &bam::HeaderView) -> Result<String> {
    let header_text =
        std::str::from_utf8(header.as_bytes()).context("BAM header is not valid UTF-8")?;

    for line in header_text.lines() {
        if !line.starts_with("@RG") {
            continue;
        }
        for field in line.split('\t') {
            if let Some(sm) = field.strip_prefix("SM:") {
                return Ok(sm.to_string());
            }
        }
    }

    anyhow::bail!("no SM tag found in any @RG line of the BAM/CRAM header")
}

/// Compute GC fraction over a reference region `[start, start+len)` (0-based, half-open).
///
/// Clamps to the sequence bounds. Returns `None` if the region is empty.
pub(crate) fn gc_frac(seq: &[u8], start: usize, len: usize) -> Option<f32> {
    if seq.is_empty() || len == 0 {
        return None;
    }
    let end = (start + len).min(seq.len());
    if end <= start {
        return None;
    }
    let slice = &seq[start..end];
    let gc = slice.iter().filter(|&&b| b == b'G' || b == b'C').count();
    Some(gc as f32 / slice.len() as f32)
}

/// Convert a CLI `--max-pileup-depth` value to the cap to pass to
/// `rust_htslib::bam::pileup::Pileups::set_max_depth`. `0` means "unlimited"
/// from the user's perspective; htslib has no true unlimited, so we use the
/// largest signed value rust-htslib will accept.
pub(crate) fn resolve_max_pileup_depth(cli_value: u32) -> u32 {
    if cli_value == 0 {
        i32::MAX as u32
    } else {
        cli_value
    }
}

/// Open an indexed BAM/CRAM reader.
///
/// When `index` is `Some`, the reader is opened against that explicit index file,
/// which lets the index live under a different name or directory than `input`.
/// When `None`, htslib infers the index at the conventional location next to
/// `input` (e.g. `input.bam.bai` / `input.bam` → `input.bai`).
pub fn open_bam(
    input: &Path,
    reference: &Path,
    index: Option<&Path>,
) -> Result<bam::IndexedReader> {
    let mut reader = match index {
        Some(index_path) => bam::IndexedReader::from_path_and_index(input, index_path)
            .with_context(|| {
                format!(
                    "failed to open BAM/CRAM {} with index {}",
                    input.display(),
                    index_path.display()
                )
            })?,
        None => bam::IndexedReader::from_path(input)
            .with_context(|| format!("failed to open BAM/CRAM: {}", input.display()))?,
    };
    reader
        .set_reference(reference)
        .with_context(|| format!("failed to set reference: {}", reference.display()))?;
    Ok(reader)
}

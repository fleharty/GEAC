use std::path::Path;

use anyhow::{Context, Result};
use rust_htslib::{bam, faidx};

/// Caches one chromosome sequence at a time to avoid per-base faidx seeks.
/// When the chromosome changes, the new sequence is fetched and the old one dropped.
pub(crate) struct RefCache {
    fai: faidx::Reader,
    current_tid: Option<usize>,
    current_chrom: String,
    /// Upper-cased sequence bytes for the cached chromosome (0-based)
    current_seq: Vec<u8>,
}

impl RefCache {
    pub(crate) fn new(reference: &Path) -> Result<Self> {
        let fai = faidx::Reader::from_path(reference)
            .with_context(|| format!("failed to open reference FASTA: {}", reference.display()))?;
        Ok(Self {
            fai,
            current_tid: None,
            current_chrom: String::new(),
            current_seq: Vec::new(),
        })
    }

    /// Returns (chrom_name, ref_base) for the given tid and 0-based position.
    /// Loads the chromosome sequence on first access or when the chromosome changes.
    /// `targets` is a pre-built slice of (chrom_name, chrom_len) indexed by tid.
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

pub fn open_bam(input: &Path, reference: &Path) -> Result<bam::IndexedReader> {
    let mut reader = bam::IndexedReader::from_path(input)
        .with_context(|| format!("failed to open BAM/CRAM: {}", input.display()))?;
    reader
        .set_reference(reference)
        .with_context(|| format!("failed to set reference: {}", reference.display()))?;
    Ok(reader)
}

use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::Path;

use anyhow::{bail, Context, Result};
use duckdb::Connection;
use rust_htslib::bam::{self, Read};
use tracing::info;

use crate::build_fusion_index::hamming_neighbors_canonical;
use crate::cli::DiagnoseFusionArgs;
use crate::kmer::{kmer_iter, non_acgt_positions, render_layout_track};

/// Decode a 2-bit canonical k-mer hash to a DNA string (A=0 C=1 G=2 T=3, MSB =
/// leftmost base). The string is the canonical orientation.
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
    String::from_utf8(bases).unwrap()
}

/// Index k-mers for one gene: canonical hash → genome-wide copy count (None when
/// the index was built without --check-genome-uniqueness).
fn load_index_gene_kmers(
    index_path: &Path,
    gene: &str,
    k: usize,
) -> Result<HashMap<u64, Option<i32>>> {
    let conn = Connection::open(index_path)
        .with_context(|| format!("failed to open fusion index: {}", index_path.display()))?;

    if let Ok(s) = conn.query_row("SELECT value FROM meta WHERE key='kmer_size'", [], |r| {
        r.get::<_, String>(0)
    }) {
        let ik: usize = s
            .parse()
            .with_context(|| format!("index meta kmer_size is not an integer: {s:?}"))?;
        anyhow::ensure!(
            ik == k,
            "--kmer-size {k} does not match the index (built with k={ik}); pass --kmer-size {ik}"
        );
    }

    let n_in_genes: i64 = conn.query_row(
        "SELECT count(*) FROM genes WHERE gene_name = ?",
        [gene],
        |r| r.get(0),
    )?;
    if n_in_genes == 0 {
        bail!("gene '{gene}' is not present in the index");
    }

    let mut stmt = conn.prepare(
        "SELECT k.kmer_hash, k.genome_copies FROM kmers k \
         JOIN genes g ON k.gene_index = g.gene_index WHERE g.gene_name = ?",
    )?;
    let mut map = HashMap::new();
    let mut rows = stmt.query([gene])?;
    while let Some(row) = rows.next()? {
        let h: i64 = row.get(0)?;
        let copies: Option<i32> = row.get(1)?;
        map.insert(h as u64, copies);
    }
    Ok(map)
}

/// Per-read k-mer evidence for the (gene_a, gene_b) pair.
struct ReadEvidence {
    qname: String,
    chrom: String, // "unmapped" when unplaced
    pos: i64,      // 0-based; -1 when unmapped
    mapq: u8,
    a_positions: Vec<usize>,      // k-mer start positions matching gene A
    b_positions: Vec<usize>,      // k-mer start positions matching gene B
    a_hits: Vec<u64>,             // the actual gene-A k-mers hit (for error-path flagging)
    b_hits: Vec<u64>,             // the actual gene-B k-mers hit
    n_base_positions: Vec<usize>, // read positions of non-ACGT bases (N etc.)
    n_windows: usize,             // read_len - k + 1
}

impl ReadEvidence {
    fn a_count(&self) -> usize {
        self.a_positions.len()
    }
    fn b_count(&self) -> usize {
        self.b_positions.len()
    }
    fn spanning(&self) -> bool {
        self.a_count() > 0 && self.b_count() > 0
    }
    /// Spanning AND each gene anchored by >= min_anchor k-mers AND the two k-mer
    /// blocks are disjoint along the read (a real A→B or B→A junction). Mirrors
    /// `fusions::read_coherence`.
    fn coherent(&self, k: usize, min_anchor: usize) -> bool {
        if !self.spanning() || self.a_count() < min_anchor || self.b_count() < min_anchor {
            return false;
        }
        let a_min = *self.a_positions.iter().min().unwrap();
        let a_max = *self.a_positions.iter().max().unwrap();
        let b_min = *self.b_positions.iter().min().unwrap();
        let b_max = *self.b_positions.iter().max().unwrap();
        (a_max + k <= b_min) || (b_max + k <= a_min)
    }
    /// Per-window track of length n_windows showing where each gene's k-mers land:
    /// `A`/`B` = a gene-A/gene-B k-mer matched; `N` = the window contains a non-ACGT
    /// base so no k-mer could be emitted (masked); `.` = a k-mer was emitted but
    /// matched neither gene (chimeric junction or non-indexed sequence). Shared with
    /// the `fusions` `FL` BAM tag via [`render_layout_track`].
    fn layout(&self, k: usize) -> String {
        render_layout_track(
            self.n_windows,
            &self.a_positions,
            &self.b_positions,
            &self.n_base_positions,
            k,
        )
    }
}

fn median_usize(values: &mut [usize]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort_unstable();
    let n = values.len();
    if n % 2 == 1 {
        values[n / 2] as f64
    } else {
        (values[n / 2 - 1] + values[n / 2]) as f64 / 2.0
    }
}

pub fn diagnose_fusion(args: &DiagnoseFusionArgs) -> Result<()> {
    let k = args.kmer_size as usize;
    anyhow::ensure!((1..=31).contains(&k), "kmer-size must be between 1 and 31");
    let min_anchor = args.min_anchor as usize;

    let a_kmers = load_index_gene_kmers(&args.index, &args.gene_a, k)?;
    let b_kmers = load_index_gene_kmers(&args.index, &args.gene_b, k)?;
    info!(
        gene_a = %args.gene_a,
        gene_b = %args.gene_b,
        index_kmers_a = a_kmers.len(),
        index_kmers_b = b_kmers.len(),
        "loaded index k-mer sets"
    );

    let mut reader = bam::Reader::from_path(&args.reads)
        .with_context(|| format!("failed to open evidence BAM/CRAM: {}", args.reads.display()))?;
    if let Some(ref ref_path) = args.reference {
        reader
            .set_reference(ref_path)
            .context("failed to set reference for CRAM decoding")?;
    }
    let header = bam::Header::from_template(reader.header());
    let header_view = bam::HeaderView::from_header(&header);

    let want: HashSet<String> = [args.gene_a.clone(), args.gene_b.clone()]
        .into_iter()
        .collect();

    // Collect per-read evidence for fragments tagged with exactly this gene pair.
    let mut evidence: Vec<ReadEvidence> = Vec::new();
    let mut n_records = 0u64;
    let mut record = bam::Record::new();
    while let Some(r) = reader.read(&mut record) {
        r.context("error reading evidence BAM record")?;
        n_records += 1;

        // Filter to fragments whose FX tag names this exact pair (order-independent).
        let fx = match record.aux(b"FX") {
            Ok(bam::record::Aux::String(s)) => s,
            _ => continue,
        };
        let parts: HashSet<String> = fx.split("::").map(|s| s.to_string()).collect();
        if parts != want {
            continue;
        }

        let seq = record.seq().as_bytes();
        let n_windows = seq
            .len()
            .saturating_sub(k)
            .saturating_add(if seq.len() >= k { 1 } else { 0 });
        let mut a_positions = Vec::new();
        let mut b_positions = Vec::new();
        let mut a_hits = Vec::new();
        let mut b_hits = Vec::new();
        for (kmer, pos) in kmer_iter(&seq, k) {
            if a_kmers.contains_key(&kmer) {
                a_positions.push(pos);
                a_hits.push(kmer);
            } else if b_kmers.contains_key(&kmer) {
                b_positions.push(pos);
                b_hits.push(kmer);
            }
        }
        // Positions of non-ACGT bases (N etc.). kmer_iter treats any such base as a
        // window break, so these blank out k-mer windows in the layout track.
        let n_base_positions = non_acgt_positions(&seq);

        let (chrom, pos) = if record.tid() >= 0 {
            let name = std::str::from_utf8(header_view.tid2name(record.tid() as u32))
                .unwrap_or("?")
                .to_string();
            (name, record.pos())
        } else {
            ("unmapped".to_string(), -1)
        };

        evidence.push(ReadEvidence {
            qname: String::from_utf8_lossy(record.qname()).into_owned(),
            chrom,
            pos,
            mapq: record.mapq(),
            a_positions,
            b_positions,
            n_base_positions,
            a_hits,
            b_hits,
            n_windows,
        });
    }

    if evidence.is_empty() {
        bail!(
            "no reads in {} are tagged FX:Z for the pair {}::{} \
             (checked {} records) — is this the right evidence BAM and gene pair?",
            args.reads.display(),
            args.gene_a,
            args.gene_b,
            n_records
        );
    }

    // ── Build the report ──────────────────────────────────────────────────────
    let mut out: Box<dyn Write> = match &args.output {
        Some(p) => Box::new(std::io::BufWriter::new(
            std::fs::File::create(p)
                .with_context(|| format!("cannot create report: {}", p.display()))?,
        )),
        None => Box::new(std::io::BufWriter::new(std::io::stdout())),
    };

    let fragments: HashSet<&str> = evidence.iter().map(|e| e.qname.as_str()).collect();
    let spanning: Vec<&ReadEvidence> = evidence.iter().filter(|e| e.spanning()).collect();
    let n_coherent_reads = spanning
        .iter()
        .filter(|e| e.coherent(k, min_anchor))
        .count();
    let n_interleaved_reads = spanning.len() - n_coherent_reads;
    let n_no_hits = evidence
        .iter()
        .filter(|e| e.a_count() == 0 && e.b_count() == 0)
        .count();

    // Coherent fragments: a fragment is coherent if any of its reads is coherent.
    let mut coherent_frags: HashSet<&str> = HashSet::new();
    for e in &evidence {
        if e.coherent(k, min_anchor) {
            coherent_frags.insert(e.qname.as_str());
        }
    }

    let mut anchor_a: Vec<usize> = spanning.iter().map(|e| e.a_count()).collect();
    let mut anchor_b: Vec<usize> = spanning.iter().map(|e| e.b_count()).collect();
    let med_a = median_usize(&mut anchor_a);
    let med_b = median_usize(&mut anchor_b);

    writeln!(out, "diagnose-fusion: {}::{}", args.gene_a, args.gene_b)?;
    writeln!(out, "evidence BAM: {}", args.reads.display())?;
    writeln!(
        out,
        "index: {}   k={k}   min_anchor={min_anchor}",
        args.index.display()
    )?;
    writeln!(out)?;
    writeln!(out, "== Summary ==")?;
    writeln!(
        out,
        "fragments (FX={}::{}): {}",
        args.gene_a,
        args.gene_b,
        fragments.len()
    )?;
    writeln!(out, "evidence reads: {}", evidence.len())?;
    writeln!(
        out,
        "  reads with no {}/{} k-mers (mates etc.): {}",
        args.gene_a, args.gene_b, n_no_hits
    )?;
    writeln!(
        out,
        "spanning reads (k-mers from both genes): {}",
        spanning.len()
    )?;
    writeln!(
        out,
        "  coherent (disjoint A->B blocks = junction-like): {n_coherent_reads}"
    )?;
    writeln!(
        out,
        "  interleaved (overlapping/short anchor = homology-like): {n_interleaved_reads}"
    )?;
    writeln!(
        out,
        "coherent fragments (>=1 coherent read): {} / {}",
        coherent_frags.len(),
        fragments.len()
    )?;
    writeln!(
        out,
        "median anchor k-mers per spanning read: {}={med_a:.0}  {}={med_b:.0}",
        args.gene_a, args.gene_b
    )?;
    writeln!(
        out,
        "verdict: {}",
        verdict(
            spanning.len(),
            n_coherent_reads,
            med_a,
            med_b,
            &args.gene_a,
            &args.gene_b
        )
    )?;
    writeln!(out)?;

    // ── Section: original alignment map, grouped by which gene the read votes ──
    // Every read is bucketed so none are lost: a clear A/B winner, a tie (equal
    // k-mer support — typical of true junction reads), or no k-mer hits at all.
    writeln!(
        out,
        "== Original alignment of evidence reads (by voted gene) =="
    )?;
    let tie_label = format!("tie ({}=={})", args.gene_a, args.gene_b);
    let mut groups: Vec<(String, Vec<&ReadEvidence>)> = vec![
        (format!("{}-voting", args.gene_a), Vec::new()),
        (format!("{}-voting", args.gene_b), Vec::new()),
        (tie_label, Vec::new()),
        ("no k-mer hits".to_string(), Vec::new()),
    ];
    for e in &evidence {
        let idx = match (e.a_count(), e.b_count()) {
            (0, 0) => 3,
            (a, b) if a > b => 0,
            (a, b) if b > a => 1,
            _ => 2,
        };
        groups[idx].1.push(e);
    }
    for (label, group) in &groups {
        if group.is_empty() {
            continue;
        }
        writeln!(out, "{label} reads (n={}):", group.len())?;
        // chrom -> (count, min_pos, max_pos, mapqs)
        let mut by_chrom: HashMap<&str, (usize, i64, i64, Vec<u8>)> = HashMap::new();
        for e in group {
            let ent = by_chrom
                .entry(&e.chrom)
                .or_insert((0, i64::MAX, i64::MIN, Vec::new()));
            ent.0 += 1;
            if e.pos >= 0 {
                ent.1 = ent.1.min(e.pos);
                ent.2 = ent.2.max(e.pos);
            }
            ent.3.push(e.mapq);
        }
        let mut chroms: Vec<(&str, (usize, i64, i64, Vec<u8>))> = by_chrom.into_iter().collect();
        chroms.sort_by(|a, b| b.1 .0.cmp(&a.1 .0));
        for (chrom, (count, mn, mx, mut mapqs)) in chroms {
            let mut mq: Vec<usize> = mapqs.drain(..).map(|m| m as usize).collect();
            let med_mq = median_usize(&mut mq);
            if chrom == "unmapped" {
                writeln!(out, "  {chrom:<22} {count:>6}  median_mapq={med_mq:.0}")?;
            } else {
                writeln!(
                    out,
                    "  {chrom}:{mn}-{mx}   {count:>6}  median_mapq={med_mq:.0}"
                )?;
            }
        }
    }
    writeln!(out)?;

    // ── Section: suspicious k-mers from the minority (likely-spurious) side ────
    let total_a: usize = evidence.iter().map(|e| e.a_count()).sum();
    let total_b: usize = evidence.iter().map(|e| e.b_count()).sum();
    let (minority_gene, minority_b_side) = if total_b <= total_a {
        (&args.gene_b, true)
    } else {
        (&args.gene_a, false)
    };
    writeln!(
        out,
        "== Suspicious {minority_gene} k-mers (minority side: {} vs {} total k-mer hits) ==",
        total_a.min(total_b),
        total_a.max(total_b)
    )?;
    writeln!(out, "kmer_seq\tn_reads\t1edit_from_other\tref_copies")?;
    // The "other" gene's index set is the substitution target that manufactures
    // the false vote; the minority set supplies genome-copy annotation.
    let other_set: &HashMap<u64, Option<i32>> = if minority_b_side { &a_kmers } else { &b_kmers };
    let minority_set: &HashMap<u64, Option<i32>> =
        if minority_b_side { &b_kmers } else { &a_kmers };
    // Count, per minority k-mer, how many distinct reads carry it (cached hashes).
    let mut kmer_read_counts: HashMap<u64, usize> = HashMap::new();
    for e in &evidence {
        let hits = if minority_b_side {
            &e.b_hits
        } else {
            &e.a_hits
        };
        let seen: HashSet<u64> = hits.iter().copied().collect();
        for kmer in seen {
            *kmer_read_counts.entry(kmer).or_insert(0) += 1;
        }
    }

    let mut suspicious: Vec<(u64, usize)> = kmer_read_counts.into_iter().collect();
    suspicious.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    for (kmer, n_reads) in suspicious {
        // Is this minority k-mer within 1 substitution of an index k-mer of the
        // other gene? That is the sequencing-error path that manufactures the call.
        let one_edit = hamming_neighbors_canonical(kmer, k, 1)
            .into_iter()
            .find(|nb| other_set.contains_key(nb));
        let edit_str = match one_edit {
            Some(nb) => format!("yes:{}", decode_kmer(nb, k)),
            None => "no".to_string(),
        };
        let copies = minority_set
            .get(&kmer)
            .and_then(|c| *c)
            .map(|c| c.to_string())
            .unwrap_or_else(|| "NA".to_string());
        writeln!(
            out,
            "{}\t{}\t{}\t{}",
            decode_kmer(kmer, k),
            n_reads,
            edit_str,
            copies
        )?;
    }
    writeln!(out)?;

    // ── Section: per-read layout track ────────────────────────────────────────
    let per_read_header =
        "qname\tchrom\tpos\tmapq\ta_count\tb_count\tspanning\tcoherent\tlayout_5to3";
    if let Some(ref p) = args.per_read_output {
        let mut tsv = std::io::BufWriter::new(
            std::fs::File::create(p).with_context(|| format!("cannot create {}", p.display()))?,
        );
        writeln!(tsv, "{per_read_header}")?;
        for e in &evidence {
            write_read_row(&mut tsv, e, k, min_anchor)?;
        }
        writeln!(
            out,
            "== Per-read k-mer layout written to {} ==",
            p.display()
        )?;
    } else {
        writeln!(
            out,
            "== Per-read k-mer layout (first 15 reads; --per-read-output for all) =="
        )?;
        writeln!(out, "layout_5to3 legend: A=gene-A k-mer  B=gene-B k-mer  N=window masked by an N base  .=k-mer matched neither gene")?;
        writeln!(out, "{per_read_header}")?;
        for e in evidence.iter().take(15) {
            write_read_row(&mut out, e, k, min_anchor)?;
        }
    }
    out.flush()?;

    info!(
        fragments = fragments.len(),
        reads = evidence.len(),
        coherent_fragments = coherent_frags.len(),
        "fusion diagnosis complete"
    );
    Ok(())
}

fn write_read_row(w: &mut impl Write, e: &ReadEvidence, k: usize, min_anchor: usize) -> Result<()> {
    let pos = if e.pos >= 0 {
        e.pos.to_string()
    } else {
        "NA".to_string()
    };
    writeln!(
        w,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        e.qname,
        e.chrom,
        pos,
        e.mapq,
        e.a_count(),
        e.b_count(),
        e.spanning(),
        e.coherent(k, min_anchor),
        e.layout(k)
    )?;
    Ok(())
}

/// One-line heuristic interpretation of the coherence/anchor picture.
fn verdict(
    n_spanning: usize,
    n_coherent: usize,
    med_a: f64,
    med_b: f64,
    gene_a: &str,
    gene_b: &str,
) -> String {
    if n_spanning == 0 {
        return "no spanning reads — the pair is supported only by separate single-gene reads \
                (mate-based chimera); inspect the alignment map for where each mate sits"
            .to_string();
    }
    let coherent_frac = n_coherent as f64 / n_spanning as f64;
    let weak_anchor = med_a < 3.0 || med_b < 3.0;
    if coherent_frac < 0.1 {
        let weak = if weak_anchor {
            format!(
                " and one side is weakly anchored (median {gene_a}={med_a:.0}, {gene_b}={med_b:.0})"
            )
        } else {
            String::new()
        };
        format!(
            "{:.0}% of spanning reads are coherent{weak} → looks like a HOMOLOGY/PARALOG artifact \
             (interleaved k-mers, same bases matching both genes), not a chimeric junction",
            coherent_frac * 100.0
        )
    } else if coherent_frac > 0.5 {
        format!(
            "{:.0}% of spanning reads show disjoint A→B blocks → consistent with a REAL junction; \
             check the alignment map and breakpoints before dismissing",
            coherent_frac * 100.0
        )
    } else {
        format!(
            "{:.0}% of spanning reads are coherent → MIXED; inspect the per-read layout and the \
             suspicious-k-mer table to decide",
            coherent_frac * 100.0
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn evidence(
        a_positions: Vec<usize>,
        b_positions: Vec<usize>,
        n_base_positions: Vec<usize>,
        n_windows: usize,
    ) -> ReadEvidence {
        ReadEvidence {
            qname: "r".to_string(),
            chrom: "chr1".to_string(),
            pos: 0,
            mapq: 60,
            a_positions,
            b_positions,
            n_base_positions,
            a_hits: Vec::new(),
            b_hits: Vec::new(),
            n_windows,
        }
    }

    #[test]
    fn layout_marks_n_masked_windows() {
        // k=11, 50 bp read -> 40 windows. An N at read position 5 masks every window
        // whose [s, s+k) range contains position 5: starts 0..=5 (6 windows).
        let e = evidence((6..=14).collect(), (25..=39).collect(), vec![5], 40);
        let track = e.layout(11);
        assert_eq!(&track[0..6], "NNNNNN");
        assert_eq!(&track[6..15], "AAAAAAAAA");
        assert_eq!(&track[15..25], "..........");
        assert_eq!(&track[25..40], "BBBBBBBBBBBBBBB");
    }

    #[test]
    fn layout_clean_read_has_no_n() {
        let e = evidence((0..=14).collect(), (25..=39).collect(), vec![], 40);
        let track = e.layout(11);
        assert!(!track.contains('N'));
        assert_eq!(&track[0..15], "AAAAAAAAAAAAAAA");
    }
}

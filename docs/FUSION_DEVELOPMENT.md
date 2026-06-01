# GEAC — Fusion caller development & roadmap

> ⚠️ **Experimental.** This is a forward-looking design doc for the
> `geac experimental` fusion tooling. For the *current, shipped* command
> reference see [EXPERIMENTAL.md](EXPERIMENTAL.md); for completed-work history see
> [DEVELOPMENT_LOG.md](DEVELOPMENT_LOG.md). This file captures **where the fusion
> caller is going and why** — the niche it should own, and the prioritized work to
> get there.

## The niche: alignment-independent fusion detection from DNA / cfDNA

Mature fusion callers (STAR-Fusion, Arriba, FusionCatcher) are RNA-seq tools and are
better than GEAC for RNA. GEAC should **not** compete there. The open, underserved
gaps are in DNA/cfDNA:

1. **Alignment-independent DNA fusion detection.** GEAC assigns reads to genes by
  k-mer match, so it captures chimeric, soft-clipped, and unmapped reads that
   alignment-based SV callers (Manta, GRIDSS, Delly) place poorly — exactly the
   most informative reads for a junction. There is little open-source tooling for
   general fusion detection from DNA-seq BAMs that doesn't lean on the aligner.
2. **Tumor-informed longitudinal monitoring of a *known* fusion at ultra-low AF.**
  This is the differentiator with the clearest clinical value and essentially no
   open competitor. SNV-based MRD assays (NeXT Personal, Invitae PCM, GeneBits)
   detect down to ~0.008% AF but are **not** fusion-aware. Once a fusion's junction
   is known (from diagnosis), its junction-spanning k-mers are near-absent from the
   germline, so even 1–2 supporting reads in serial cfDNA is meaningful signal. A
   targeted, junction-specific index + sensitive scan would let GEAC track a fusion's
   burden across timepoints — a capability no general caller offers.

**Thesis: GEAC should be the best tool for (a) finding fusions in DNA/cfDNA without
re-alignment, and (b) monitoring a known fusion longitudinally at low allele
fraction.** Everything below is prioritized against that thesis.

## Current capabilities (shipped)

See [EXPERIMENTAL.md](EXPERIMENTAL.md) for full flag reference. In brief:

- `build-fusion-index` — gene-unique k-mer index from FASTA + GTF, with genome-wide
copy-number quantification and two-tier (panel vs genome) uniqueness.
- `fusions` — k-mer-vote read→gene assignment; fragments hitting two genes aggregated
into candidates with `--min-supporting-reads` / `--min-kmer-hits` / `--min-mapq`
thresholds. Outputs: Parquet, TSV, evidence BAM (`FX:Z:` tag), per-k-mer-hit TSV
(`kmer_pos_in_read`), breakpoint TSV, and Panel-of-Normals annotation columns.
- Junction-coherence filter: `--min-coherent-fragments` requires spanning reads to show
a disjoint A-block→B-block k-mer partition. `n_spanning_reads` and `n_coherent_reads`
are emitted as Parquet/TSV columns unconditionally.
- Fusion PoN: `--fusion-pon` annotates/filters against a `geac merge` of normal
`*.fusions.parquet` files.

## Why recurrent false positives happen (root-cause taxonomy)

Recurrent (every-sample) false calls are **structural** — caused by the
reference/annotation, not by random noise — so they're diagnosable and fixable at
the mechanism level rather than only blocklistable. In rough order of frequency:

1. **Paralog / homology.** Genes A and B share a sequence stretch. A single real,
  non-chimeric read from *one* locus carries k-mers the index treats as unique to
   *both* genes → looks like a junction. The PoN handles this worst, because the read
   is genuine signal, merely misattributed.
2. **Overlapping / adjacent genomic loci.** Nested, antisense, or abutting genes; a
  fragment spanning the boundary legitimately contains both genes' sequence.
3. **Processed pseudogenes / retrogenes / segmental duplications.** A duplicated copy
  of gene A inside gene B's locus: reads map to B but carry A's k-mers.
4. **Index built without `--check-genome-uniqueness`.** A k-mer unique among panel
  genes can still occur elsewhere genome-wide; reads from "elsewhere" misassign.

### The principled discriminator: k-mer co-linearity

`kmer_pos_in_read` (0-based k-mer offset within the read) distinguishes a real
junction from a homology artifact:

- **Real fusion:** on a spanning read, Gene-A k-mers occupy one contiguous block of
read positions and Gene-B k-mers occupy a *disjoint* block, with a single A→B
transition — two sequences spliced together.
- **Homology / paralog FP:** the two genes' k-mers are **interleaved or overlapping**
in `kmer_pos_in_read`, because the *same* read bases match both genes. No transition.

This is the highest-value signal we already emit, and it motivates the
**junction-coherence filter** below.

## Roadmap (prioritized against the niche)

### Tier 1 — Specificity (makes or breaks a clinical detector)

- **Junction-coherence / co-linearity filter.** ✅ Shipped. `--min-coherent-fragments`
requires spanning reads to show disjoint A/B k-mer blocks; `n_spanning_reads` and
`n_coherent_reads` emitted in all output. See `docs/DEVELOPMENT_LOG.md`.
- **Overlap / adjacency filter.** Reject pairs whose annotated gene bodies overlap
or sit within X bp, or whose breakpoints fall on the same chromosome within X bp
(causes #2/#3). Gene coordinates already live in the index.
- **Split-read vs discordant-pair separation.** Report `split_reads` (single read
spanning the junction; base-pair resolution) and `discordant_mates` (R1/R2 in
different genes) as distinct columns. Lets users filter on strong evidence and is a
major confidence signal. (Arriba's model.)
- **Fusion Panel-of-Normals.** Shipped — `--fusion-pon` / `--max-pon-samples`.
Treats the symptom; keep as a backstop after the mechanism-level filters above.

### Tier 2 — Clinical interpretation

- **Directionality (5′/3′ partner).** Recover which partner is 5′ from read
orientation + breakpoint geometry (currently dropped by pair normalization).
`EML4::ALK` ≠ `ALK::EML4` clinically.
- **Reading-frame awareness.** Approximate in-frame/out-of-frame prediction from
the GTF CDS; out-of-frame fusions are usually not oncogenic.
- **Known-fusion annotation.** Tag calls against COSMIC / Mitelman / ChimerDB so
known oncogenic fusions surface and novel ones are flagged.
- **VAF-like quantification.** Report supporting reads as a fraction of local
depth so a fusion's burden can be tracked across timepoints — the metric that turns
GEAC into an MRD tool. **Gates the longitudinal-monitoring niche.**

### Tier 3 — Precision & confidence

- **Base-pair breakpoint resolution.** Replace the k-mer-transition estimate with
the exact junction base from spanning reads (strand-correct). `reconstruct_fusions.sh`
does this via assembly today; fold a lightweight version in-binary.
- **Confidence score / tier.** Combine split-read count, anchor length, breakpoint
consensus tightness, and PoN status into a high/medium/low tier instead of hard
thresholds.
- **Anchor-length / complexity filter.** Require sufficient unique anchor on both
sides of the junction; reject low-complexity junctions.

### Tier 4 — The differentiated capability

- **Tumor-informed targeted monitoring mode.** Given a known junction (sequence or
breakpoint), build a tiny junction-specific index and scan serial cfDNA BAMs at high
sensitivity (allow 1–2 supporting reads), reporting per-timepoint supporting-read
counts and VAF. This is the niche no other open caller fills; depends on the VAF
work in Tier 2 and base-pair breakpoints in Tier 3.

## Benchmarking (gates everything above)

Without truth-set numbers, neither the tool nor any abstract has a credibility anchor,
and we can't tell whether a new filter actually helped.

- Assemble a fusion-positive reference set (e.g. SEQC2, EWSR1-FLI1 Ewing lines)
with cfDNA-like dilution series.
- Report sensitivity / specificity and breakpoint accuracy per release; re-run on
each Tier-1/Tier-2 change to confirm it improves rather than regresses.

## Diagnosing a specific recurrent FP (operational recipe)

For any recurrent `A::B`, pull its rows from the `--kmer-hits-output` TSV and check:

1. **Do A and B k-mers overlap in `kmer_pos_in_read`?** Yes → homology (#1).
2. **Where do supporting reads map (`chrom`, `pos`)?** Single cluster → pseudogene /
  seg-dup / overlap (#2/#3); cross-check the breakpoint TSV for `chrom_a == chrom_b`
   with small `|breakpoint_a − breakpoint_b|`.
3. `**genome_copies` of the driving k-mers?** Take `kmer_seq` values into
  `geac experimental lookup-kmer` / `locate-kmer`; multiple hits → non-unique k-mer
   leakage (#4), rebuild with `--check-genome-uniqueness` and call with
   `--max-kmer-copies 1`.


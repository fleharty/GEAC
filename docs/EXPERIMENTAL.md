# GEAC — Experimental commands

> ⚠️ **Experimental.** Everything under `geac experimental <...>` is pre-production.
> APIs, flags, output schemas, and the DuckDB index layout may change between
> releases without notice, and these commands are **not** covered by the same test
> guarantees as the stable `collect` / `coverage` / `cohort` surface. They are
> documented separately here so the main [README](../README.md) stays focused on
> tested functionality. Use them for exploration, not production pipelines.

These tools implement an alignment-independent, k-mer-based **gene-fusion detection**
workflow. The core idea: build an index of k-mers that are unique to each gene body,
then scan a BAM/CRAM and flag read fragments whose two ends match two *different*
genes — a signature of a fusion or chimeric read.

> 🧭 For the forward-looking design — the niche this caller aims to own, the
> false-positive root-cause taxonomy, and the prioritized improvement roadmap — see
> [FUSION_DEVELOPMENT.md](FUSION_DEVELOPMENT.md).

## Command overview

| Command | Purpose |
|---------|---------|
| [`build-fusion-index`](#build-fusion-index) | Build the per-gene unique-k-mer index (DuckDB) from a GTF + reference FASTA |
| [`fusions`](#fusions) | Scan a BAM/CRAM and report candidate gene fusions |
| [`build-fusion-kmer-blacklist`](#build-fusion-kmer-blacklist) | Aggregate per-sample kmer-hits TSVs from normal samples into a k-mer blacklist Parquet |
| [`extract-gene`](#extract-gene) | Pull all reads whose fragments hit a named gene's k-mers |
| [`locate-kmer`](#locate-kmer) | Find every occurrence of a single k-mer in the reference FASTA |
| [`lookup-kmer`](#lookup-kmer) | Look up which gene a single k-mer is assigned to in an index |
| [`scan-read`](#scan-read) | Show per-k-mer matches for one read sequence against an index |
| [`compute-uniqueness-map`](#compute-uniqueness-map) | For every genomic locus, compute the smallest k for which the k-mer at that position is genome-unique |
| [`shared-kmers`](#shared-kmers) | List the k-mers two genes' bodies share — the cross-gene k-mers that index building discards |
| [`diagnose-fusion`](#diagnose-fusion) | Diagnose a suspected false-positive fusion call from a `fusions` evidence BAM |

`lookup-kmer`, `scan-read`, `shared-kmers`, and `diagnose-fusion` are debugging aids for
inspecting the index and the calls it produces.

---

## Uniqueness model (read this first)

A k-mer's value for fusion detection depends on how unique it is. The index build
distinguishes two tiers:

- **Panel-wide unique** — the default. After extracting every gene body's k-mers, any
  k-mer occurring in **2+ of the indexed gene bodies** is discarded (cross-gene
  dedup). With `--genes <panel>`, "indexed gene bodies" means just your panel, so a
  k-mer can survive even if it is repetitive elsewhere in the genome.
- **Genome-wide unique** — enabled with `--check-genome-uniqueness`. Additionally
  scans the whole FASTA and records each surviving k-mer's genome-wide occurrence
  count, then (by default) keeps only those occurring exactly once.

Because no large gene is ever 100% genome-unique (paralogs, tandem repeats — e.g.
DUX4 sits near ~16% genome-unique), the genome tier also supports a **near-unique**
relaxation via `--max-genome-copies N`: keep k-mers occurring up to `N` times
genome-wide. The per-k-mer copy count is stored in the index so it can be re-tightened
at call time without rebuilding (see [`fusions --max-kmer-copies`](#fusions)).

---

## `build-fusion-index`

Build the DuckDB index of per-gene unique k-mers.

```bash
geac experimental build-fusion-index \
    --gene-annotation gencode.v47.annotation.gtf.gz \
    --fasta Homo_sapiens_assembly38.fasta \
    --genes my_genes.txt \
    --output panel_fusion_index.duckdb \
    --kmer-size 23 \
    --min-gene-kmers 1 \
    --check-genome-uniqueness \
    --max-genome-copies 2 \
    --bed-output panel_kmers.bed \
    --bed-output-by-copies panel_kmers \
    --gene-stats-output my_genes.stats.tsv \
    --copy-histogram-output my_genes.copy_hist.tsv
```

| Flag | Default | Description |
|------|---------|-------------|
| `--gene-annotation <PATH>` | — | Gene annotation file. Accepted formats: GTF (`.gtf`, `.gtf.gz`) — GENCODE, NCBI RefSeq, or UCSC GTF; GenePred/refFlat (`.txt`, `.txt.gz`) — UCSC genePredExt+bin (e.g. `ncbiRefSeq.txt.gz`), genePredExt without bin, or refFlat. Format is auto-detected from file extension and column layout. |
| `--fasta <PATH>` | — | Reference FASTA; requires a `.fai` (`samtools faidx`). |
| `--kmer-size <N>` | `23` | K-mer length (1–31). Must match the value used by `fusions`/`extract-gene`. |
| `--min-gene-kmers <N>` | `1` | Drop genes with fewer than `N` retained unique k-mers. Applied **after** the genome-copy filter. Default `1` keeps any gene with at least one usable k-mer; raise it to suppress repetitive genes. |
| `--output <PATH>` | — | Output DuckDB index. |
| `--genes <PATH>` | — | Restrict to the gene names in this file (one per line; `#` comments allowed). Builds a small targeted panel. |
| `--check-genome-uniqueness` | off | Scan the whole FASTA and record per-k-mer genome-wide copy counts. Required for `--max-genome-copies > 1`, `--copy-histogram-output`, and `--bed-output-by-copies`. Holds the candidate k-mer set in RAM (~20 bytes × n_candidates — tens of GB for a whole-genome build; cheap for a panel). |
| `--max-genome-copies <N>` | `1` | Retain k-mers occurring `1..=N` times genome-wide. `N>1` keeps near-unique k-mers (helps repetitive genes). Requires `--check-genome-uniqueness`. |
| `--bed-output <PATH>` | — | BED of merged intervals covering **all retained** k-mer positions. Sorted in FASTA order; usable directly in IGV / bedtools. |
| `--bed-output-by-copies <PREFIX>` | — | One BED **per genome-wide copy tier**: `<PREFIX>.copies1.bed`, `.copies2.bed`, … up to `--max-genome-copies`. Lets you load strictly-unique vs near-unique regions as separate tracks. Requires `--check-genome-uniqueness`. See [column 4](#per-copy-tier-bed-annotation) below. |
| `--gene-stats-output <PATH>` | — | Per-gene uniqueness TSV (see [columns](#gene-stats-tsv)). |
| `--copy-histogram-output <PATH>` | — | Genome-wide copy-number histogram TSV (`copies`, `n_kmers`). Requires `--check-genome-uniqueness`. |

### Index schema (DuckDB)

```
genes(gene_index UINTEGER, gene_name VARCHAR, chrom VARCHAR)
kmers(kmer_hash BIGINT, gene_index UINTEGER, genome_copies INTEGER)   -- genome_copies NULL unless --check-genome-uniqueness ran
kmer_positions(kmer_hash BIGINT, chrom VARCHAR, pos INTEGER)          -- first occurrence position
```

`genome_copies` is the genome-wide occurrence count per k-mer; it is `NULL` for
indexes built without `--check-genome-uniqueness`. The column is additive and does
not break older readers.

### gene-stats TSV

Always emitted (panel tier):

`gene, chrom, start, end, body_len, panel_unique_kmers, kmer_windows,
pct_unique_windows, covered_bases, pct_unique_bases`

- `pct_unique_windows` = `panel_unique_kmers / kmer_windows` (windows = `body_len − k + 1`).
- `pct_unique_bases` = gene-body bases covered by ≥1 panel-unique k-mer / `body_len`.

When `--check-genome-uniqueness` ran, these columns are appended:

`n_copies_1, n_copies_2, n_copies_3_10, n_copies_gt10, genome_unique_pct_bases,
relaxed_pct_bases`

- Copy buckets count the gene's panel-unique k-mers by genome-wide occurrence.
- `genome_unique_pct_bases` = coverage from strictly-unique (==1) k-mers.
- `relaxed_pct_bases` = coverage from k-mers with ≤ `--max-genome-copies` copies.

Stats are computed **before** the `--min-gene-kmers` drop, so even genes that get
excluded report their true uniqueness.

### Per-copy-tier BED annotation

For tiers ≥2 (`*.copies2.bed` and higher), each interval carries a 4th column listing
every **full-GTF** gene the interval's k-mers occur in (comma-separated, sorted), or
`intergenic`. This names the paralog(s) driving a region's non-uniqueness even when
they fall outside the `--genes` panel:

```
chr2  29450  29500  GENEB,GENEBP1
chr7  55020  55090  EGFR,intergenic
```

The `copies1.bed` tier is plain 3-column BED (a unique k-mer's only "match" is its own
gene). Adjacent tiers can overlap by up to `k−1` bases at boundaries — they are
independent k-mer sets with overlapping ±k windows, not a strict partition.

**Known limitations.** `--max-genome-copies` only relaxes *genome-wide* paralog
copies; it does not rescue *intra*-gene tandem repeats. The base-coverage metrics
slightly under-count tandem-repeat genes because only each k-mer's first position is
stored.

---

## `fusions`

Scan a BAM/CRAM and report candidate gene fusions. Reads are assigned to genes by
k-mer match (alignment-independent), then fragments whose ends map to two different
genes are aggregated. Secondary, supplementary, and unmapped reads are included
(only duplicates are skipped) since the most informative chimeric reads are often
poorly placed by the aligner.

```bash
geac experimental fusions \
    --bam sample.bam \
    --index panel_fusion_index.duckdb \
    --output sample.fusions.parquet \
    --tsv-output sample.fusions.tsv \
    --reads-output sample.fusion_reads.bam \
    --breakpoints-output sample.fusions.breakpoints.tsv \
    --max-kmer-copies 1
```

| Flag | Default | Description |
|------|---------|-------------|
| `--bam <PATH>` | — | Input BAM/CRAM. |
| `--index <PATH>` | — | Fusion index from `build-fusion-index`. |
| `--reference <PATH>` | — | Reference FASTA (required for CRAM). |
| `--sample-id <STR>` | `@RG SM` | Sample id for output; falls back to the read-group `SM` tag. |
| `--output <PATH>` | — | Fusion candidates Parquet. |
| `--kmer-size <N>` | `23` | Must match the index. |
| `--min-kmer-hits <N>` | `3` | Min unique k-mer hits for a read→gene assignment. |
| `--min-supporting-reads <N>` | `2` | Min fragments supporting a fusion to report it. |
| `--min-mapq <N>` | `0` | Mapping-quality floor (unmapped reads bypass it). |
| `--reads-output <PATH>` | — | BAM of all reads from fusion-supporting fragments. Each record is tagged `FX:Z:GENEA::GENEB` with the fusion it supports, and `FL:Z:<track>` with that read's per-k-mer-window layout against the pair (reference 5'→3', i.e. the BAM `SEQ` orientation; `A`/`B` = a gene-A/gene-B k-mer, where A/B is the first/second gene in `FX`; `N` = window masked by a non-ACGT base; `.` = k-mer matching neither gene — the same string as `diagnose-fusion`'s `layout_5to3`). A BAI is built automatically. The pair order in `FX` matches the `gene_a`/`gene_b` columns of the Parquet/TSV and the `fusion` column of `--kmer-hits-output` (all use gene-index order), so the tag joins cleanly back to the tables. |
| `--tsv-output <PATH>` | — | Human-readable fusion table. |
| `--kmer-hits-output <PATH>` | — | One row per k-mer hit from fusion reads. The `kmer_pos_in_read` column gives each k-mer's 0-based offset within the read, which locates the Gene A→Gene B transition on junction-spanning reads. |
| `--breakpoints-output <PATH>` | — | Per-fusion breakpoint TSV. Projects the k-mer transition on junction-spanning reads onto genomic coordinates; shares the same second BAM pass as `--kmer-hits-output` (only one extra scan when both are set). |
| `--max-kmer-copies <N>` | — | Ignore k-mers occurring > `N` times genome-wide (or with unknown copy count). Requires an index built with `--check-genome-uniqueness`; errors clearly otherwise. Lets you re-tighten or relax uniqueness at call time without rebuilding. |
| `--fusion-pon <PATH>` | — | Fusion Panel-of-Normals DuckDB (a `geac merge` of normal-sample `*.fusions.parquet` files). Annotates every call with `n_pon_samples`, `pon_total_samples`, `max_pon_supporting_reads`. Matching is by alphabetically-sorted gene-name pair, so PoN and call may use different indexes. |
| `--max-pon-samples <N>` | — | Tag fusions seen in strictly more than `N` PoN samples with `filter=pon` (rows are kept). Requires `--fusion-pon`. Default: annotate only, every row stays `filter=PASS`. |
| `--fusion-kmer-blacklist <PATH>` | — | K-mer blacklist Parquet from `build-fusion-kmer-blacklist`. A read must have at least `--min-kmer-hits` *non-blacklisted* k-mer matches to contribute to fusion evidence. Blacklisted k-mers still contribute once a read has passed on clean evidence — they just cannot be the sole basis for including a read. Complementary to `--fusion-pon`: the PoN filters at the gene-pair level; the k-mer blacklist filters at the read level, upstream. |
| `--min-kmer-blacklist-samples <N>` | `1` | Treat a k-mer as blacklisted only if it appears in at least `N` PoN samples in the blacklist Parquet (`n_pon_samples` column). Allows reusing one blacklist file at different stringency thresholds without rebuilding. |
| `--max-breakpoint-std <BP>` | — | Breakpoint-consensus filter (no PoN needed). Tags a fusion `filter=chimera` (rows kept) unless **both** breakpoints are supported by ≥ `--min-breakpoint-reads` spanning reads whose position estimates have a standard deviation ≤ `BP`. A real fusion's reads converge on one junction (std on the order of tens of bp — junction-adjacent splice isoforms plus k-mer transition-point estimation noise spread a high-depth real event beyond a single base); paralog/PCR-chimera artifacts splice at scattered positions (std in the thousands to millions of bp), so they fail. The real-vs-artifact gap spans orders of magnitude, so the exact cutoff is not delicate. **Recommended: `100`** (a tight 10 bp cutoff falsely rejects real high-depth fusions whose junction has any isoform/estimation spread — see DEVELOPMENT_LOG). When spanning reads are too sparse to pin a breakpoint — e.g. the junction falls inside an index-excluded repeat (an intronic Alu), so the indexed-k-mer gap is wider than a read and nothing spans it — the estimate is supplemented with **concordant single-gene mates**, using each read's junction-facing alignment edge; estimates are reduced to the dominant cluster (±1 kb of the median) before measuring spread, so a few deep-body mates can't inflate it. Triggers the same second BAM pass as `--breakpoints-output`. Default: disabled. |
| `--min-breakpoint-reads <N>` | `5` | Minimum spanning reads converging on **each** breakpoint under `--max-breakpoint-std`. Guards against tiny-`n` chance agreement (2 reads at the same coordinate give std 0 but mean nothing). Only consulted when `--max-breakpoint-std` is set. |
| `--min-breakpoint-distance <BP>` | — | Same-locus / adjacency filter. Tags a fusion `filter=samelocus` (rows kept) when **both** breakpoints fall on the same chromosome within `BP` bp. Both partners localizing to one spot is the signature of single-locus paralog leakage: reads from an *unindexed* paralog (e.g. GNA12, chr7) carry k-mers shared with two indexed cousins (GNA13 on chr17, GNA11 on chr19) and split between them, fabricating a fusion whose breakpoints both sit at the GNA12 locus (Δ ≈ a few bp). That fake breakpoint is tight and reproducible, so `--max-breakpoint-std` passes it — this catches it on geometry instead. A genuine fusion places partners on different chromosomes, or far apart on one. Recommended: `10000` (≫ a fragment, ≪ inter-gene distance). Triggers the same second BAM pass as `--breakpoints-output`. Default: disabled. |

Output Parquet/TSV columns: `sample_id, gene_a, gene_b, chrom_a, chrom_b,
supporting_reads, min_mapq, n_spanning_reads, n_coherent_fragments, n_pon_samples,
pon_total_samples, max_pon_supporting_reads, filter`. The three `*pon*` columns are
`0`/`0`/`NA` unless `--fusion-pon` is given. `filter` is `PASS` by default, `pon`
when `--max-pon-samples` flags a PoN-recurrent call, `chimera` when
`--max-breakpoint-std` rejects a call for lacking breakpoint consensus, or
`samelocus` when `--min-breakpoint-distance` rejects a call whose breakpoints
co-localize to one locus (single-locus paralog leakage).

`--kmer-hits-output` columns: `fusion, sample_id, read_name, read_end, chrom, pos,
gene_matched, kmer_pos_in_read, kmer_hash, kmer_seq`.

`--breakpoints-output` columns: `fusion, gene_a, chrom_a, breakpoint_a, bp_a_n,
bp_a_std, gene_b, chrom_b, breakpoint_b, bp_b_n, bp_b_std, n_spanning_reads`. The
`bp_*_n` columns count spanning reads contributing to each side and `bp_*_std` is
the spread in bp (low = tight consensus). Coordinates are 1-based; partners with no
spanning-read support report `NA`.

### Scan progress

The BAM/CRAM scan logs a progress line every 5 M reads with an estimate of how far
through the file it is and an ETA, so long runs are legible:

```
INFO geac::fusions: BAM scan progress reads_processed=5000000 reads_assigned=3434824 percent=23.4 eta=7m12s reads_per_sec=412000
```

- For **BAM**, `percent` comes from the compressed byte offset (BGZF virtual offset),
  so it tracks the file directly.
- For **CRAM**, the byte offset is not meaningful, so `percent` is estimated from the
  current read's mapped genomic coordinate against the total reference length. This
  assumes coordinate-sorted input (the production norm); during the trailing unmapped
  block it pins at ~100%. If the header carries no reference lengths, only
  `reads_per_sec` is reported.

---

## `build-fusion-kmer-blacklist`

Aggregate per-sample `--kmer-hits-output` TSV files from normal (PoN) samples into a
k-mer blacklist Parquet. Each row in the output records a k-mer hash and the number of
distinct normal samples in which that k-mer was observed supporting a fusion candidate.

The blacklist works at the read level: rather than suppressing entire gene pairs at call
time, individual reads are disqualified from contributing to fusion evidence if they do
not have enough *clean* (non-blacklisted) k-mer hits. This lets a genuine junction call
survive even when some of its supporting k-mers are noisy in normals, as long as the
read also carries clean junction-spanning k-mers.

```bash
# Step 1: run fusions on each normal with --kmer-hits-output.
for n in normal1 normal2; do
  geac experimental fusions --bam $n.bam --index panel.duckdb \
      --output $n.fusions.parquet \
      --kmer-hits-output $n.kmer_hits.tsv \
      --min-supporting-reads 1
done

# Step 2: aggregate into a k-mer blacklist.
geac experimental build-fusion-kmer-blacklist \
    --kmer-hits normal1.kmer_hits.tsv normal2.kmer_hits.tsv \
    --output normals.fusion_kmer_blacklist.parquet \
    --min-pon-samples 1
```

| Flag | Default | Description |
|------|---------|-------------|
| `--kmer-hits <PATH>...` | — | One or more `--kmer-hits-output` TSV files from normal-sample `fusions` runs. |
| `--output <PATH>` | — | Output Parquet. Schema: `kmer_hash BIGINT, n_pon_samples BIGINT`. |
| `--min-pon-samples <N>` | `1` | Only write k-mers seen in at least `N` distinct normal samples. `1` captures any k-mer observed in any normal; increase to require recurrence across multiple normals before blacklisting. |

### Output schema

```
kmer_hash    BIGINT   -- same encoding as the kmer_hash column in --kmer-hits-output
n_pon_samples BIGINT  -- number of distinct normal samples carrying this k-mer
```

### Using the blacklist at call time

Pass the blacklist Parquet to `fusions` with `--fusion-kmer-blacklist`. The
`--min-kmer-blacklist-samples` flag lets you re-apply a different threshold at call time
without rebuilding:

```bash
geac experimental fusions --bam tumor.bam --index panel.duckdb \
    --output tumor.fusions.parquet --tsv-output tumor.fusions.tsv \
    --fusion-pon fusion_pon.duckdb \
    --fusion-kmer-blacklist normals.fusion_kmer_blacklist.parquet \
    --min-kmer-blacklist-samples 1
```

The gene-pair PoN (`--fusion-pon`) and the k-mer blacklist (`--fusion-kmer-blacklist`)
are complementary and can be used together. The k-mer blacklist acts upstream (at read
assignment), the gene-pair PoN acts downstream (at call annotation/tagging).

---

## `extract-gene`

Pull every read whose fragment carries a k-mer hit to one or more named genes. Both
reads of a pair are written even if only one carries the hit.

```bash
geac experimental extract-gene \
    --bam sample.bam \
    --index panel_fusion_index.duckdb \
    --gene PAX3 --gene FOXO1 \
    --output pax3_foxo1.reads.bam \
    --complement-output other.reads.bam
```

| Flag | Default | Description |
|------|---------|-------------|
| `--bam <PATH>` | — | Input BAM/CRAM. |
| `--index <PATH>` | — | Fusion index. |
| `--reference <PATH>` | — | Reference (required for CRAM). |
| `--sample-id <STR>` | `@RG SM` | Sample id for TSV output. |
| `--gene <NAME>` | — | Target gene; repeatable. Required. |
| `--output <PATH>` | — | BAM of reads from fragments hitting any target gene. |
| `--kmer-size <N>` | `23` | Must match the index. |
| `--min-kmer-hits <N>` | `1` | Min unique k-mer hits to a target gene. |
| `--min-mapq <N>` | `0` | Mapping-quality floor. |
| `--kmer-hits-output <PATH>` | — | Per-k-mer-hit TSV for matching reads. |
| `--complement-output <PATH>` | — | BAM of reads from fragments with *no* target hit; every input record lands in exactly one of `--output`/`--complement-output`. BAI built automatically. |
| `--debug-region <chr:start-end>` | — | Log per-read k-mer-matching detail for reads overlapping this 1-based region (debugging "missing reads"). |

---

## `locate-kmer`

Find every occurrence of a single k-mer (and its reverse complement) in the reference
FASTA. Does not need an index.

```bash
geac experimental locate-kmer \
    --fasta Homo_sapiens_assembly38.fasta \
    --kmer ACGTACGTACGTACGTACGTACG \
    --gene-annotations gencode.v47.annotation.gtf.gz
```

| Flag | Default | Description |
|------|---------|-------------|
| `--fasta <PATH>` | — | Reference FASTA (`.fai` required). |
| `--kmer <SEQ>` | — | K-mer to search (length 1–31). Forward + reverse-complement reported; palindromes reported once as `+`. |
| `--output <PATH>` | stdout | TSV: `chrom, start, end, strand` (+ `gene, region` with annotations). |
| `--gene-annotations <PATH>` | — | GTF/GFF3/genePred to label each hit with gene name and region (exon/CDS/UTR/intron/intergenic). |

---

## `lookup-kmer`

Look up a single k-mer in an index and report the gene it is assigned to and its first
recorded position.

```bash
geac experimental lookup-kmer --index panel_fusion_index.duckdb --kmer ACGT...ACG
```

The canonical (min of forward/reverse-complement) form is queried. Length 1–31.

---

## `scan-read`

Scan one read sequence against an index and print every matching k-mer with its gene
and position — the read-level analog of what `fusions` does internally.

```bash
geac experimental scan-read --index panel_fusion_index.duckdb --read ACGT... --kmer-size 23
```

`N` bases reset the sliding window exactly as in `fusions`.

---

## `compute-uniqueness-map`

For every position in the genome, compute the smallest k such that the k-mer starting
at that locus appears **exactly once** genome-wide (in canonical / strand-collapsed
form). The result is a per-base signal: repetitive regions receive large values; regions
with high-complexity, unique sequence receive small values. It is useful for:

- Understanding which loci in a fusion gene panel are accessible at a given k-mer size.
- Choosing `--kmer-size` for `build-fusion-index` by inspecting the distribution over
  panel gene bodies.
- Loading alongside a `copies1.bed` track in IGV to see uniqueness dropouts in context.

```bash
# Whole-genome map (see RAM note below)
geac experimental compute-uniqueness-map \
    --fasta Homo_sapiens_assembly38.fasta \
    --output hg38_min_unique_k.bedgraph

# Targeted: write only positions overlapping a panel BED
geac experimental compute-uniqueness-map \
    --fasta Homo_sapiens_assembly38.fasta \
    --output panel_min_unique_k.bedgraph \
    --regions panel_fusion_kmers.bed
```

| Flag | Default | Description |
|------|---------|-------------|
| `--fasta <PATH>` | — | Reference FASTA; requires a `.fai` (`samtools faidx`). |
| `--output <PATH>` | — | Output bedgraph (see [output format](#output-format) below). |
| `--min-k <N>` | `15` | Smallest k to test. |
| `--max-k <N>` | `31` | Largest k to test (hard ceiling: 31, set by the 2-bit k-mer encoder). Positions with no unique k in `[min-k, max-k]` are assigned `max-k + 1`. |
| `--regions <PATH>` | — | Optional BED of output regions. The genome-wide count passes always scan the full FASTA (uniqueness is global), but only positions overlapping these regions appear in the output. |
| `--no-merge` | off | Write one line per base instead of merging adjacent equal-value positions. Produces a larger file; content is identical. |

### Output format

The output is a **bedgraph** (0-based half-open coordinates):

```
chrom   start   end   min_unique_k
chr1    0       1     18
chr1    1       2     22
chr1    2       3     17
chr1    950     952   32   ← repetitive; no unique k found in [15,31] → sentinel max_k+1
```

The value at position `p` is the smallest k such that the k-mer `[p, p+k)` appears
exactly once genome-wide. Adjacent positions with the same value are merged into a
single interval (run-length encoding). Because adjacent k-mer start positions almost
always differ in their min-unique-k, most intervals are 1 bp wide in practice. The
format is full base-pair resolution; IGV and bedtools load it correctly. Use
`--no-merge` if you need a strictly fixed-step layout.

Positions in N-runs are never yielded by the k-mer iterator and receive the sentinel
value `max-k + 1`.

### Algorithm

For each k from `min-k` to `max-k`:

1. **Count pass** — scan the full genome with the rolling canonical k-mer iterator,
   accumulating a `HashMap<u64, u8>` of genome-wide occurrence counts (capped at 2).
2. **Query pass** — scan again; for each position not yet resolved, check whether its
   k-mer has count == 1. If so, record `k` as that position's minimum unique k.
3. Drop the count map and advance k. Stop early if all positions are resolved.

### Resource requirements

| Genome | RAM (peak) | Typical runtime |
|--------|-----------|-----------------|
| Targeted panel (~few Mb output, genome-wide counting) | ~50–60 GB | 30–60 min |
| Full hg38 | ~55–65 GB | 3–8 hours |

The dominant cost is the `HashMap<u64, u8>` in the count pass. At k ≈ 23, nearly every
k-mer in hg38 is distinct (~3 billion entries × ~20 bytes ≈ 60 GB). A machine with
**≥ 64 GB RAM** is required for whole-genome use; a 128 GB server is comfortable.

`--regions` limits the output file size but does not reduce counting RAM — uniqueness
is a global property that requires scanning the whole genome regardless.

---

## Typical workflow

```bash
# 1. Build a targeted panel index with near-unique k-mers retained, plus diagnostics.
geac experimental build-fusion-index \
    --gene-annotation gencode.v47.annotation.gtf.gz --fasta ref.fasta \
    --genes panel.txt --output panel.duckdb \
    --min-gene-kmers 1 --check-genome-uniqueness --max-genome-copies 2 \
    --bed-output-by-copies panel_kmers \
    --gene-stats-output panel.stats.tsv --copy-histogram-output panel.copy_hist.tsv

# 2. Inspect: which genes have poor genome-unique coverage? (e.g. tandem repeats)
#    Load panel_kmers.copies1.bed and panel_kmers.copies2.bed as separate IGV tracks.

# 3. Call fusions, re-tightening to strictly-unique k-mers at call time.
#    --breakpoints-output localizes each junction from spanning reads.
geac experimental fusions --bam sample.bam --index panel.duckdb \
    --output sample.fusions.parquet --tsv-output sample.fusions.tsv \
    --reads-output sample.fusion_reads.bam \
    --breakpoints-output sample.fusions.breakpoints.tsv --max-kmer-copies 1

# 4. (Optional) reference-free reconstruction of the junction contig.
scripts/reconstruct_fusions.sh -b sample.fusion_reads.bam -r two_gene_miniref.fa
```

### Building a fusion Panel-of-Normals

GEAC provides two complementary PoN mechanisms:

- **Gene-pair PoN** (`--fusion-pon`): suppresses entire gene pairs that recur across
  normals at the call level. Simple but coarse — a real low-AF tumor call can be
  suppressed if the same pair appears in normals due to unrelated noise.
- **K-mer blacklist** (`--fusion-kmer-blacklist`): identifies the specific k-mers
  driving noise in normals. A read must have enough *clean* (non-blacklisted) k-mer
  support to count as evidence. Finer-grained: a tumor call can survive even if some
  of its k-mers appear in normals, as long as it also has junction-spanning k-mers
  that are clean.

Both can be built from the same normal runs and used together.

```bash
# Step 1: run each normal with sensitive settings and capture kmer-hits.
for n in normal1 normal2 normal3; do
  geac experimental fusions --bam $n.bam --index panel.duckdb \
      --output $n.fusions.parquet \
      --kmer-hits-output $n.kmer_hits.tsv \
      --min-supporting-reads 1
done

# Step 2a: build the gene-pair PoN DuckDB.
geac merge --output fusion_pon.duckdb \
    normal1.fusions.parquet normal2.fusions.parquet normal3.fusions.parquet

# Step 2b: build the k-mer blacklist.
geac experimental build-fusion-kmer-blacklist \
    --kmer-hits normal1.kmer_hits.tsv normal2.kmer_hits.tsv normal3.kmer_hits.tsv \
    --output fusion_kmer_blacklist.parquet \
    --min-pon-samples 1

# Step 3: call fusions on a tumor, applying both PoN layers.
geac experimental fusions --bam tumor.bam --index panel.duckdb \
    --output tumor.fusions.parquet --tsv-output tumor.fusions.tsv \
    --fusion-pon fusion_pon.duckdb \
    --fusion-kmer-blacklist fusion_kmer_blacklist.parquet
```

---

## `shared-kmers`

Pair up the k-mers that two genes' bodies have in common. At the default edit
distance of 0 these are the **exact** shared k-mers — precisely the cross-gene k-mers
that [`build-fusion-index`](#build-fusion-index) discards during panel-wide dedup (see
[Uniqueness model](#uniqueness-model-read-this-first)), so they never appear in a built
index. The usual trigger is a suspicious `GENEA::GENEB` call where you want to know
*which* sequence the two genes share and *where* it sits in each.

With `--edit-distance N` the match becomes **fuzzy**: a gene-A k-mer is paired with any
gene-B k-mer within `N` substitutions, even when the two sequences are not identical.
This catches **diverged paralogs** — gene-family members ~95% identical share few exact
23-mers but many near-identical ones, and those near-matches are exactly what can drive
cross-assignment even though the index's exact-dedup never removed them.

### Explaining false-positive fusion calls

Fuzzy matching plus `--index <built index>` is the tool for asking *"is this
`GENEA::GENEB` call a k-mer artifact?"*. A k-mer-driven false call happens when a read
that is truly from gene A picks up a sequencing error that turns one of its k-mers into
a k-mer the **index** has assigned to gene B — producing a spurious gene-B vote. So the
culprits are pairs where a gene-A k-mer and a gene-B k-mer are a single substitution
apart **and the voted gene's k-mer is actually in the index** (only index k-mers can be
voted; exact cross-gene k-mers were already deduped out).

Run with `--edit-distance 1 --index …` and look at the rows where `ab_dist > 0` and the
voted gene's `in_index` column is `true`: each is a concrete error path that can
manufacture the false call. (`--check-reference` composes here — a fragile index k-mer
that is also high-copy is doubly suspect.)

Gene bodies are parsed from the same annotation `build-fusion-index` consumes, k-mers
are canonicalized identically, and a gene name that maps to several bodies (e.g. the
same symbol on multiple contigs) has all of them unioned.

```bash
# Exact shared k-mers only
geac experimental shared-kmers \
    --gene-annotation gencode.v44.annotation.gtf \
    --fasta Homo_sapiens_assembly38.fasta \
    --gene-a BCR \
    --gene-b ABL1 \
    --kmer-size 23 \
    --output bcr_abl1_shared.tsv

# Fuzzy + index: pair k-mers within 1 substitution and flag which are real index
# k-mers — the rows that can drive a false BCR::ABL1 call
geac experimental shared-kmers \
    --gene-annotation gencode.v44.annotation.gtf \
    --fasta Homo_sapiens_assembly38.fasta \
    --gene-a BCR --gene-b ABL1 --kmer-size 23 \
    --edit-distance 1 --index hg38_fusion_index.duckdb \
    --check-reference --threads 8 \
    --output bcr_abl1_near.tsv
```

| Flag | Default | Description |
|------|---------|-------------|
| `--gene-annotation <PATH>` | — | GTF or GenePred annotation (same formats as `build-fusion-index`). |
| `--fasta <PATH>` | — | Reference FASTA; requires a `.fai` (`samtools faidx`). |
| `--gene-a <NAME>` | — | First gene name (must match a `gene_name` in the annotation). |
| `--gene-b <NAME>` | — | Second gene name. |
| `--kmer-size <N>` | `23` | K-mer length; use the same value as the index you are debugging. |
| `--edit-distance <N>` | `0` | Max substitution (Hamming) distance between a gene-A and gene-B k-mer for a match. `0` = exact shared only; `N > 0` also pairs near-identical k-mers. |
| `--output <PATH>` | stdout | Output TSV. |
| `--check-reference` | off | Scan the whole FASTA and add `ref_copies_a` / `ref_copies_b` columns: genome-wide occurrences of each matched k-mer. |
| `--index <PATH>` | — | Built fusion index (DuckDB). Adds `in_index_a` / `in_index_b` columns flagging whether each matched k-mer is a real index k-mer the caller can vote on. Must share the k-mer size. |
| `--threads <N>` | `0` (all cores) | Threads for the `--check-reference` genome scan. No effect unless that flag is set. |

### Output format

A TSV with one row per matched k-mer pair, sorted by k-mer. `pos_a` / `pos_b` are the
**0-based first-occurrence** chromosome coordinates in each gene's body; `ab_dist` is
the substitution distance between the two sequences (always `0` when
`--edit-distance 0`). A summary line (per-gene k-mer counts, match count, exact-match
count) is logged to stderr.

```
kmer_seq_a               chrom_a  pos_a    kmer_seq_b               chrom_b  pos_b    ab_dist
AAAACCCCGGGAACCCCGGGGTTT chr22    23290412 AAAACCCCGGGAACCCCGGGGTTT chr9     133710001 0
GGTTGAAAACCGCCCGGGGTTTAC chr22    23291004 GGTTCAAAACCGCCCGGGGTTTAC chr9     133710512 1
```

- `kmer_seq_a` / `kmer_seq_b` are the **canonical** (strand-collapsed) orientation. For
  exact matches (`ab_dist 0`) they are identical; for fuzzy matches they differ at up
  to `ab_dist` positions.
- `ref_copies_a` / `ref_copies_b` (only with `--check-reference`) count every
  genome-wide occurrence of each canonical k-mer. A value of `1` means the k-mer is
  unique to its gene body; higher values mean it also appears in other genes, paralogs,
  or repeats — i.e. it was repetitive sequence to begin with.
- `in_index_a` / `in_index_b` (only with `--index`) are `true` when that side's k-mer is
  an actual index k-mer the fusion caller can vote on. Exact shared k-mers (`ab_dist 0`)
  are typically `false` on both sides — they were removed by cross-gene dedup — whereas a
  near-pair with the voted gene's side `true` is a live false-positive path.

> A given gene-A k-mer can match several gene-B k-mers (and vice versa) under a non-zero
> edit distance, so each qualifying pair is emitted as its own row. The `--check-reference`
> scan reads the entire FASTA but only tracks the matched k-mers, so memory stays
> proportional to the match set, not the genome.

---

## `diagnose-fusion`

Take the **evidence BAM** written by [`fusions --reads-output`](#fusions) (reads tagged
`FX:Z:GENEA::GENEB`) plus the index, and explain *why* a suspected `GENEA::GENEB` call
fired. It re-derives each evidence read's gene-A / gene-B k-mer hits and reports four
things, each targeting a distinct failure mode.

```bash
geac experimental diagnose-fusion \
    --reads sample.fusion_reads.bam \
    --index hg38_fusion_index.duckdb \
    --gene-a BCR --gene-b ABL1 --kmer-size 23 \
    --per-read-output bcr_abl1_reads.tsv
```

| Flag | Default | Description |
|------|---------|-------------|
| `--reads <PATH>` | — | Evidence BAM/CRAM from `fusions --reads-output`. |
| `--index <PATH>` | — | The fusion index used for the call. Must share the k-mer size. |
| `--reference <PATH>` | — | Reference FASTA, required only if `--reads` is a CRAM. |
| `--gene-a` / `--gene-b` | — | The suspected pair (order-independent; matched against the `FX` tag). |
| `--kmer-size <N>` | `23` | K-mer length; must match the index and the original run. |
| `--min-anchor <N>` | `3` | K-mers per gene for a spanning read to count as anchored on both sides (matches `fusions --min-anchor-kmers`). |
| `--output <PATH>` | stdout | Human-readable report. |
| `--per-read-output <PATH>` | — | Per-read k-mer layout track as TSV (otherwise a short preview is in the report). |

### What it reports

1. **Homology vs junction summary** — per-fragment coherence (reusing the caller's own
   `read_coherence` logic): how many spanning reads have **disjoint A→B blocks**
   (junction-like) vs **interleaved** k-mers (homology/paralog artifact — the same read
   bases match both genes), plus median anchor k-mers per side and a one-line verdict.
2. **Original alignment map** — where the evidence reads actually align (their original
   BAM coordinates), bucketed by which gene they vote for. If reads voting gene B really
   sit at gene A's locus, they are mis-voting, not chimeric.
3. **Suspicious-k-mer table** — for the minority (likely-spurious) gene, each k-mer it
   contributed, how many reads carry it, whether it is **within one substitution of an
   index k-mer of the other gene** (the sequencing-error path that manufactures the
   call), and its genome-wide copy number.
4. **Per-read layout track** — a per-k-mer-window string (reference 5'→3', the BAM
   `SEQ` orientation) showing where each gene's k-mers land, so interleaving vs blocking
   is visible at a glance:

   | char | meaning |
   |------|---------|
   | `A` / `B` | a gene-A / gene-B k-mer matched at that window |
   | `N` | the window contains a non-ACGT base (an `N`), so no k-mer could be emitted — a *masked* gap, not a sequence gap |
   | `.` | a k-mer was emitted but matched neither gene (chimeric junction, or sequence not in the index) |

```
qname   chrom  pos  mapq  a_count  b_count  spanning  coherent  layout_5to3
read0   chr22  …    60    15       15       true      true      AAAAAAAAAAAAAAA..........BBBBBBBBBBBBBBB
read1   chr22  …    60     9       15       true      true      NNNNNNAAAAAAAAA..........BBBBBBBBBBBBBBB
read2   chr22  …    60     9        8       true      false     ABABABAABABBABAB...
```

The `N`/`.` distinction matters when interpreting a wide gap between the A and B blocks:
a single `N` masks roughly `k` windows (the k-mer iterator resets on any non-ACGT base
and must refill), so a gap of `N`s is read-quality masking (common in duplex/consensus
BAMs that write `N` at disagreement positions), whereas a gap of `.`s at the junction is
genuine chimeric sequence — a candidate breakpoint, possibly with a small insertion.

Because the track follows the BAM `SEQ` (reference-forward) orientation, a reverse-strand
read is laid out 3'→5' relative to the original fragment, so its A/B blocks appear in the
mirror order. The gene assignment is unaffected — k-mers are canonicalized — so counts,
spanning, and coherence are orientation-independent; only the left-right order of blocks
flips. (The `FL` tag written by `fusions --reads-output` uses this same orientation.)

A call where ~all spanning reads are **interleaved** (`ABAB`), the minority side is
**weakly anchored**, and its k-mers are **1 substitution from the other gene** is a
textbook k-mer artifact rather than a real fusion.

---

## Diagnosing false-positive fusion calls

`shared-kmers` and `diagnose-fusion` are complementary: they attack the same problem —
*why is GENEA::GENEB being called when it isn't real?* — from two different levels.

| Question | Tool | Inputs |
|----------|------|--------|
| Does the **index/reference** predispose this pair to confusion, for *any* sample? | [`shared-kmers`](#shared-kmers) | annotation + FASTA (+ index) — **no BAM** |
| Why did *this* call fire in *this* sample's **reads**? | [`diagnose-fusion`](#diagnose-fusion) | the sample's evidence BAM + index |

### The two mechanisms

A k-mer-based fusion caller produces a false `GENEA::GENEB` when reads pick up gene-A and
gene-B index k-mers without a real chimeric junction. Two distinct mechanisms cause this,
and the tools above measure each directly:

1. **Homology / paralog near-collisions.** Related genes share *near*-identical sequence.
   Exact shared k-mers are removed by cross-gene dedup at index-build time, but k-mers a
   single substitution apart survive in *both* genes' index entries. A read covering a
   conserved region then matches both — with the gene-A and gene-B k-mers **interleaved**
   along the read (the same bases support both), not partitioned into an A-block and a
   B-block.
2. **Sequencing-error paths.** A read truly from gene A picks up a base error that turns
   one of its k-mers into a gene-B *index* k-mer, casting a spurious gene-B vote.

`read_coherence` (used by both `fusions` and `diagnose-fusion`) is the discriminator:
real junctions give **disjoint, coherent** A→B blocks; both false mechanisms give
**interleaved or weakly-anchored** evidence.

### Recommended sequence

```bash
# A. Index level — is the pair intrinsically collision-prone? (no BAM needed)
#    Rows with ab_dist=1 and the voted gene's in_index=true are the index near-collisions.
geac experimental shared-kmers \
    --gene-annotation gencode.v44.gtf --fasta hg38.fa \
    --gene-a BCR --gene-b ABL1 --kmer-size 23 \
    --edit-distance 1 --index hg38_fusion_index.duckdb --check-reference \
    --output bcr_abl1_index.tsv

# B. Read level — re-run the original fusions call WITH the evidence BAM if you don't
#    already have one, then diagnose this sample's call.
geac experimental fusions --bam sample.bam --index hg38_fusion_index.duckdb \
    --kmer-size 23 --output sample.fusions.parquet \
    --reads-output sample.fusion_reads.bam

geac experimental diagnose-fusion \
    --reads sample.fusion_reads.bam --index hg38_fusion_index.duckdb \
    --gene-a BCR --gene-b ABL1 --kmer-size 23
```

### Reading the result, and what to do about it

- **`diagnose-fusion` verdict says "homology artifact"** (almost all spanning reads
  interleaved): a paralog collision. Confirm at the index level with `shared-kmers`
  (many `ab_dist=1`, `in_index=true` pairs). **Fix:** raise `--kmer-size`, or rebuild the
  index with `--check-genome-uniqueness` / a higher `--min-gene-kmers` so the colliding
  near-unique k-mers are dropped.
- **Suspicious-k-mer table shows `1edit_from_other = yes`** for the minority gene: an
  error-path artifact concentrated in a few fragile index k-mers. **Fix:** drop those
  k-mers via the [k-mer blacklist](#using-the-blacklist-at-call-time) (if they recur
  across normals) or rebuild with `--edit-distance-filter 1` so no index k-mer sits one
  substitution from a reference k-mer.
- **Alignment map shows the partner-voting reads sitting at the other gene's locus:**
  they are mis-voting, not chimeric — same conclusion, same fixes.
- **Verdict says "real junction"** (coherent blocks, both sides well-anchored, reads at
  two loci): treat the call as genuine; use `--breakpoints-output` and
  [`reconstruct_fusions.sh`](#reconstruct_fusionssh) to localize and assemble the
  junction.

Recurrent false pairs across many samples are better handled at the cohort level with a
[fusion Panel-of-Normals and k-mer blacklist](#building-a-fusion-panel-of-normals); the
two tools here are for understanding *why* a specific pair misbehaves so you can choose
the right cohort-level fix.

---

## Helper scripts

**Experimental.** These live in `scripts/` and complement the `geac experimental`
fusion commands; like the commands themselves, their interfaces may change.

### `reconstruct_fusions.sh`

Reference-free reconstruction of fusion junctions from a fusion *evidence BAM*
(records tagged `FX:Z:GENEA::GENEB` by `geac experimental fusions --reads-output`).
For each fusion it extracts that fusion's reads, assembles them de novo with CAP3,
and (optionally) aligns the contigs to a reference with minimap2 so the breakpoint
appears as a split/supplementary alignment. This is a sequence-level complement to
`--breakpoints-output`: the latter estimates coordinates from k-mer positions,
while this rebuilds the actual junction contig.

```bash
scripts/reconstruct_fusions.sh -b evidence.bam [-o outdir] [-r ref.fa] \
                               [-f "GENEA::GENEB"] [-t threads]
```

Requires `samtools` (≥1.12) and `cap3`; `minimap2` is needed only when `-r` is given.

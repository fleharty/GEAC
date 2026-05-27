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

## Command overview

| Command | Purpose |
|---------|---------|
| [`build-fusion-index`](#build-fusion-index) | Build the per-gene unique-k-mer index (DuckDB) from a GTF + reference FASTA |
| [`fusions`](#fusions) | Scan a BAM/CRAM and report candidate gene fusions |
| [`extract-gene`](#extract-gene) | Pull all reads whose fragments hit a named gene's k-mers |
| [`locate-kmer`](#locate-kmer) | Find every occurrence of a single k-mer in the reference FASTA |
| [`lookup-kmer`](#lookup-kmer) | Look up which gene a single k-mer is assigned to in an index |
| [`scan-read`](#scan-read) | Show per-k-mer matches for one read sequence against an index |

`lookup-kmer` and `scan-read` are debugging aids for inspecting the index.

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
    --gtf gencode.v47.annotation.gtf.gz \
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
| `--gtf <PATH>` | — | Gene annotation (`.gtf` / `.gtf.gz`). Only `gene` feature lines are read. |
| `--fasta <PATH>` | — | Reference FASTA; requires a `.fai` (`samtools faidx`). |
| `--kmer-size <N>` | `23` | K-mer length (1–31). Must match the value used by `fusions`/`extract-gene`. |
| `--min-gene-kmers <N>` | `100` | Drop genes with fewer than `N` retained unique k-mers. Applied **after** the genome-copy filter. |
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
| `--reads-output <PATH>` | — | BAM of all reads from fusion-supporting fragments. Each record is tagged `FX:Z:GENEA::GENEB` with the fusion it supports; a BAI is built automatically. |
| `--tsv-output <PATH>` | — | Human-readable fusion table. |
| `--kmer-hits-output <PATH>` | — | One row per k-mer hit from fusion reads (enables per-read detail collection). |
| `--max-kmer-copies <N>` | — | Ignore k-mers occurring > `N` times genome-wide (or with unknown copy count). Requires an index built with `--check-genome-uniqueness`; errors clearly otherwise. Lets you re-tighten or relax uniqueness at call time without rebuilding. |

Output Parquet/TSV columns: `sample_id, gene_a, gene_b, chrom_a, chrom_b,
supporting_reads, min_mapq`.

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

## Typical workflow

```bash
# 1. Build a targeted panel index with near-unique k-mers retained, plus diagnostics.
geac experimental build-fusion-index \
    --gtf gencode.v47.annotation.gtf.gz --fasta ref.fasta \
    --genes panel.txt --output panel.duckdb \
    --min-gene-kmers 1 --check-genome-uniqueness --max-genome-copies 2 \
    --bed-output-by-copies panel_kmers \
    --gene-stats-output panel.stats.tsv --copy-histogram-output panel.copy_hist.tsv

# 2. Inspect: which genes have poor genome-unique coverage? (e.g. tandem repeats)
#    Load panel_kmers.copies1.bed and panel_kmers.copies2.bed as separate IGV tracks.

# 3. Call fusions, re-tightening to strictly-unique k-mers at call time.
geac experimental fusions --bam sample.bam --index panel.duckdb \
    --output sample.fusions.parquet --tsv-output sample.fusions.tsv \
    --reads-output sample.fusion_reads.bam --max-kmer-copies 1
```

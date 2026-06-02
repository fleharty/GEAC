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
| `--reads-output <PATH>` | — | BAM of all reads from fusion-supporting fragments. Each record is tagged `FX:Z:GENEA::GENEB` with the fusion it supports; a BAI is built automatically. The pair order in `FX` matches the `gene_a`/`gene_b` columns of the Parquet/TSV and the `fusion` column of `--kmer-hits-output` (all use gene-index order), so the tag joins cleanly back to the tables. |
| `--tsv-output <PATH>` | — | Human-readable fusion table. |
| `--kmer-hits-output <PATH>` | — | One row per k-mer hit from fusion reads. The `kmer_pos_in_read` column gives each k-mer's 0-based offset within the read, which locates the Gene A→Gene B transition on junction-spanning reads. |
| `--breakpoints-output <PATH>` | — | Per-fusion breakpoint TSV. Projects the k-mer transition on junction-spanning reads onto genomic coordinates; shares the same second BAM pass as `--kmer-hits-output` (only one extra scan when both are set). |
| `--max-kmer-copies <N>` | — | Ignore k-mers occurring > `N` times genome-wide (or with unknown copy count). Requires an index built with `--check-genome-uniqueness`; errors clearly otherwise. Lets you re-tighten or relax uniqueness at call time without rebuilding. |
| `--fusion-pon <PATH>` | — | Fusion Panel-of-Normals DuckDB (a `geac merge` of normal-sample `*.fusions.parquet` files). Annotates every call with `n_pon_samples`, `pon_total_samples`, `max_pon_supporting_reads`. Matching is by alphabetically-sorted gene-name pair, so PoN and call may use different indexes. |
| `--max-pon-samples <N>` | — | Tag fusions seen in strictly more than `N` PoN samples with `filter=pon` (rows are kept). Requires `--fusion-pon`. Default: annotate only, every row stays `filter=PASS`. |
| `--fusion-kmer-blacklist <PATH>` | — | K-mer blacklist Parquet from `build-fusion-kmer-blacklist`. A read must have at least `--min-kmer-hits` *non-blacklisted* k-mer matches to contribute to fusion evidence. Blacklisted k-mers still contribute once a read has passed on clean evidence — they just cannot be the sole basis for including a read. Complementary to `--fusion-pon`: the PoN filters at the gene-pair level; the k-mer blacklist filters at the read level, upstream. |
| `--min-kmer-blacklist-samples <N>` | `1` | Treat a k-mer as blacklisted only if it appears in at least `N` PoN samples in the blacklist Parquet (`n_pon_samples` column). Allows reusing one blacklist file at different stringency thresholds without rebuilding. |

Output Parquet/TSV columns: `sample_id, gene_a, gene_b, chrom_a, chrom_b,
supporting_reads, min_mapq, n_pon_samples, pon_total_samples,
max_pon_supporting_reads`. The three `*pon*` columns are `0`/`0`/`NA` unless
`--fusion-pon` is given.

`--kmer-hits-output` columns: `fusion, sample_id, read_name, read_end, chrom, pos,
gene_matched, kmer_pos_in_read, kmer_hash, kmer_seq`.

`--breakpoints-output` columns: `fusion, gene_a, chrom_a, breakpoint_a, bp_a_n,
bp_a_std, gene_b, chrom_b, breakpoint_b, bp_b_n, bp_b_std, n_spanning_reads`. The
`bp_*_n` columns count spanning reads contributing to each side and `bp_*_std` is
the spread in bp (low = tight consensus). Coordinates are 1-based; partners with no
spanning-read support report `NA`.

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

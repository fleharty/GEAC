# GEAC CLI Reference

This page summarizes the stable GEAC commands. Run `geac <command> --help` for the
authoritative flag list for your installed binary.

## Commands

```text
collect          Process one BAM/CRAM into alt-base Parquet records
merge            Merge Parquet files or DuckDB databases into one cohort DuckDB
inspect          Check a merged cohort DuckDB for schema, metadata, and resource issues
qc               Print per-sample QC summaries from Parquet files
cohort           Summarize recurrent loci across per-sample Parquet files
sample-identity  Find samples in a merged cohort that look like the same individual
annotate-normal  Cross-check tumor loci against a paired normal BAM/CRAM
annotate-pon     Cross-check tumor loci against a Panel of Normals DuckDB
coverage         Emit per-position or binned coverage metrics
fragments        Emit per-fragment insert-size, GC, and end-motif metrics
experimental     K-mer / fusion tools; see docs/EXPERIMENTAL.md
```

## `geac collect`

Runs a pileup over one BAM/CRAM and writes one row per observed alt allele.

```bash
geac collect \
  --input sample.bam \
  --reference hg38.fa \
  --output SAMPLE_001.parquet \
  --read-type duplex \
  --pipeline fgbio
```

Important options:

| Flag | Default | Description |
|---|---:|---|
| `--sample-id` | BAM `SM` tag | Override sample identifier |
| `--subject-id` | none | Biological subject identifier shared by related samples |
| `--sample-type` | none | Free-text substrate type such as `cfDNA`, `tumor_tissue`, `normal_tissue` |
| `--batch` | none | Processing batch label |
| `--label1`, `--label2`, `--label3` | none | Free-text metadata labels |
| `--timepoint` | none | Longitudinal timepoint label |
| `--reads-output` | off | Also write `.reads.parquet` alongside `.locus.parquet` |
| `--input-checksum-sha256` | off | Hash input BAM/CRAM and store checksum in output provenance |
| `--read-type` | `duplex` | `raw`, `simplex`, or `duplex` |
| `--pipeline` | none | Free-text pipeline label stored in the `pipeline` column. `fgbio`/`dragen` (case-insensitive) also select built-in family-size tag schemes (`aD`/`bD`/`cD` and `XV`/`XW`); any other value records the label and reads no family size unless `--family-size-tags` is given. |
| `--family-size-tags` | none | Override family-size aux tags, as `ab=XX,ba=YY,total=ZZ[,fallback=sum\|none]` (e.g. `ab=aD,ba=bD,total=cD,fallback=sum`). Overrides the `--pipeline` preset. A tag absent from a read yields a null family size. |
| `--vcf` | none | VCF/BCF variant-call annotation |
| `--variants-tsv` | none | TSV variant annotation alternative to `--vcf` |
| `--gnomad` | none | bgzip+tabix gnomAD VCF/BCF for AF annotation |
| `--gnomad-af-field` | `AF` | gnomAD INFO field to read |
| `--gnomad-index`, `--gnomad-index-uri` | none | Path / canonical URI of the gnomAD index (`.tbi`/`.csi`), recorded as `gnomad_index_path` so IGV references it (`--gnomad-index-uri` wins). Note: geac's own AF annotation still requires the index co-located with `--gnomad` (htslib VCF reader limitation). |
| `--targets` | none | BED or Picard interval list for target annotation and sample metrics |
| `--targets-uri` | none | Canonical target path or URI to store in output metadata |
| `--gene-annotations` | none | GFF3, GTF, or UCSC genePred gene annotations |
| `--min-base-qual` | `1` | Minimum base quality |
| `--min-map-qual` | `0` | Minimum mapping quality |
| `--max-pileup-depth` | `0` | Max reads per pileup column; `0` disables htslib downsampling |
| `--include-duplicates` | off | Count duplicate reads |
| `--include-secondary` | off | Count secondary alignments |
| `--include-supplementary` | off | Count supplementary alignments |
| `--exclude-tag` | none | Drop reads whose aux tag equals a value, as `TAG:VALUE` (e.g. `RX:bad`). Exact string match across string/char/integer tags; reads lacking the tag are kept. Repeatable. |
| `--region` | whole input | Region string or BED/Picard interval list |
| `--repeat-window` | `10` | Bases on each side for homopolymer/STR metrics |
| `--index` | inferred next to `--input` | Path to the BAM/CRAM index (`.bai`/`.crai`) when it is not at the conventional location (different name or directory). Also recorded as `bai_path` for IGV unless `--bai-uri` is given. |
| `--bam-uri`, `--bai-uri` | local paths | Canonical BAM/index URIs for IGV and cloud provenance. `--bai-uri` takes precedence over `--index` for the stored `bai_path`. |
| `--variants-uri`, `--gnomad-uri` | local paths | Canonical annotation URIs for provenance |

When `--reads-output` is set, `sample.parquet` becomes:

- `sample.locus.parquet`
- `sample.reads.parquet`

When `--targets` is set, collect also writes `sample.sample_metrics.parquet`.

## `geac merge`

Builds a cohort DuckDB from Parquet files and/or existing DuckDB databases.

```bash
geac merge --output cohort.duckdb samples/*.parquet
```

Parquet routing is suffix-based:

| Suffix | DuckDB table |
|---|---|
| `.reads.parquet` | `alt_reads` |
| `.normal_evidence.parquet` | `normal_evidence` |
| `.pon_evidence.parquet` | `pon_evidence` |
| `.coverage.parquet` | `coverage` |
| `.coverage.intervals.parquet` | `coverage_intervals` |
| `.sample_metrics.parquet` | `sample_metrics` |
| `.fragments.parquet` | `fragments` |
| `.fusions.parquet` | `fusions` |
| anything else | `alt_bases` |

DuckDB inputs are attached and known tables are copied into the output. The `samples`
summary table is rebuilt from `alt_bases`.

Useful option:

| Flag | Description |
|---|---|
| `--on-target-tsv <PATH>` | Write a TSV of `on_target = true` alt-base calls |

## `geac inspect`

Checks a merged cohort DuckDB for common operational issues: missing provenance
tables, missing required `alt_bases` columns, empty core tables, optional feature
tables, conflicting sample metadata, missing embedded BAM/BAI/resource columns, and
local/non-URI resource paths that may not work in cloud or IGV contexts.

```bash
geac inspect --input cohort.duckdb
```

Useful option:

| Flag | Description |
|---|---|
| `--strict` | Exit non-zero when warnings are found, not only errors |

## `geac sample-identity`

Finds samples in a merged cohort whose common germline SNV fingerprints look like the
same individual.

```bash
geac sample-identity \
  --input cohort.duckdb \
  --output sample_identity_pairs.tsv
```

Tuning defaults:

| Flag | Default | Description |
|---|---:|---|
| `--min-depth` | `20` | Minimum marker depth |
| `--min-vaf` | `0.15` | Minimum non-reference marker VAF |
| `--hom-vaf` | `0.85` | Homozygous-alt genotype threshold |
| `--min-recurrence` | `2` | Marker must appear in at least this many samples |
| `--max-markers` | `5000` | Maximum marker panel size |
| `--max-samples` | `5000` | Pairwise comparison guardrail; `0` disables |
| `--min-jaccard` | `0.70` | Candidate-pair Jaccard threshold |
| `--min-concordance` | `0.95` | Genotype concordance threshold |
| `--common-gnomad-only` | off | Use common gnomAD markers when `gnomad_af` exists |
| `--gnomad-min-af` | `0.05` | Lower common-AF bound |
| `--gnomad-max-af` | `0.95` | Upper common-AF bound |
| `--all-pairs` | off | Write all pairs sharing at least one marker |

## `geac annotate-normal`

Cross-checks tumor alt loci against a paired normal BAM/CRAM.

```bash
geac annotate-normal \
  --tumor-parquet TUMOR.locus.parquet \
  --normal-bam NORMAL.bam \
  --reference hg38.fa \
  --output TUMOR.normal_evidence.parquet
```

Key options mirror pileup filters: `--normal-sample-id`, `--min-base-qual`,
`--min-map-qual`, `--max-pileup-depth`, `--include-duplicates`,
`--include-secondary`, and `--include-supplementary`.

## `geac annotate-pon`

Cross-checks tumor alt loci against a Panel of Normals DuckDB.

```bash
geac annotate-pon \
  --tumor-parquet TUMOR.locus.parquet \
  --pon-db pon.duckdb \
  --output TUMOR.pon_evidence.parquet
```

## `geac coverage`

Emits per-position or binned coverage records.

```bash
geac coverage \
  --input SAMPLE.bam \
  --reference hg38.fa \
  --output SAMPLE.coverage.parquet \
  --targets capture_targets.bed \
  --sample-id SAMPLE_001
```

Important options include `--subject-id`, `--sample-type`, `--batch`,
`--label1/2/3`, `--timepoint`, `--read-type`, `--pipeline`, `--targets`,
`--region`, `--gene-annotations`, `--min-map-qual`, `--max-pileup-depth`,
`--min-base-qual`, `--gc-window`, `--min-depth`, `--bin-size`,
`--intervals-output`, `--track NAME:FILE`, and `--fill-zeros`.

## `geac fragments`

Emits one row per proper-pair fragment with insert size, GC content, midpoint, and
4-mer end motifs.

```bash
geac fragments \
  --input SAMPLE.bam \
  --reference hg38.fa \
  --output SAMPLE.fragments.parquet \
  --read-type duplex \
  --pipeline fgbio
```

Important options include `--sample-id`, `--subject-id`, `--sample-type`,
`--batch`, `--label1/2/3`, `--timepoint`, `--read-type`, `--pipeline`,
`--region`, and `--min-map-qual`.

## `geac qc`

Prints per-sample QC summaries from one or more Parquet files.

```bash
geac qc --output qc.tsv --on-target-only samples/*.parquet
```

## `geac cohort`

Summarizes recurrent loci across per-sample Parquet files.

```bash
geac cohort --output recurrent.tsv --min-samples 3 samples/*.parquet
```

Use `--min-sample-fraction`, `--on-target-only`, and `--top-n` to tune reporting.

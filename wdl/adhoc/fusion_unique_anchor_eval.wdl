version 1.0

## fusion_unique_anchor_eval.wdl  [ADHOC — disposable]
##
## One-off analysis harness: re-run the fusion caller on the EVIDENCE BAMs from a prior
## `geac experimental fusions` cohort run (the small `*.fusions.bam` outputs, ~1 GB each)
## with `--emit-unique-anchor`, so each call gains the `n_unique_anchored` column, then
## merge the tiny per-sample Parquets into a single cohort DuckDB.
##
## Why this exists: the original cohort run was made with a build that predates the
## `n_unique_anchored` column. Rather than re-process the multi-GB *input* BAMs (expensive)
## or localize ~71 GB of evidence BAMs to a laptop, this re-runs only on the already-produced
## evidence BAMs (cheap, short tasks) and returns ONE small merged file. The unique-anchor
## N-sweep / precision-recall scoring against truth then runs offline on that single file.
##
## NOT a supported pipeline — kept in wdl/adhoc/ because it is tied to one investigation and
## likely short-lived. See wdl/adhoc/README.md.
##
## Inputs
##   evidence_bams   - Array of `*.fusions.bam` evidence BAMs from the prior run (gs:// paths).
##                     sample_id is derived from each filename stem.
##   fusion_index    - The SAME copy-labeled index used for the original run
##                     (built with --check-genome-uniqueness; has a genome_copies column).
##                     ⚠ MUST match the index that produced the evidence BAMs. Re-running with a
##                     different index changes k-mer matching and the spanning-read set, which makes
##                     the breakpoint-geometry tags (chimera/samelocus) diverge from the original run.
##                     `n_unique_anchored` itself is robust to this; the geometry filters are not.
##   docker_image    - geac image that includes the n_unique_anchored column (>= 0.4.60).
##   The filter inputs below default to the original investigation's report config; change them
##   only to keep the PASS / chimera / samelocus set faithful to the run you are evaluating.
##
## Output
##   cohort_duckdb   - one DuckDB with a `fusions` table for all samples, carrying
##                     n_unique_anchored + supporting_reads + filter (the N-sweep input).
##   per_sample_tsv  - the small per-sample fusions TSVs (convenience).

workflow FusionUniqueAnchorEval {
    input {
        Array[File] evidence_bams
        File        fusion_index
        String      docker_image

        File?   reference_fasta

        # Caller config — defaults match the investigation's report config. Keep these
        # aligned with the ORIGINAL run so the reproduced call set matches.
        Int     kmer_size               = 23
        Int     min_kmer_hits           = 1
        Int     min_supporting_reads    = 1
        Int     min_mapq                = 0
        Int     min_coherent_fragments  = 0
        Int     min_anchor_kmers        = 3
        Float   max_breakpoint_std      = 100.0
        Int     min_breakpoint_reads    = 5
        Int     min_breakpoint_distance = 10000
        Boolean skip_version_check      = true

        String  cohort_db_name          = "fusion_unique_anchor_cohort.duckdb"

        Int     memory_gb               = 16
        Int     disk_gb                 = 20
        Int     preemptible             = 2
    }

    scatter (bam in evidence_bams) {
        call ReEmitUniqueAnchor {
            input:
                evidence_bam            = bam,
                fusion_index            = fusion_index,
                reference_fasta         = reference_fasta,
                kmer_size               = kmer_size,
                min_kmer_hits           = min_kmer_hits,
                min_supporting_reads    = min_supporting_reads,
                min_mapq                = min_mapq,
                min_coherent_fragments  = min_coherent_fragments,
                min_anchor_kmers        = min_anchor_kmers,
                max_breakpoint_std      = max_breakpoint_std,
                min_breakpoint_reads    = min_breakpoint_reads,
                min_breakpoint_distance = min_breakpoint_distance,
                skip_version_check      = skip_version_check,
                docker_image            = docker_image,
                memory_gb               = memory_gb,
                disk_gb                 = disk_gb,
                preemptible             = preemptible
        }
    }

    call MergeCohort {
        input:
            parquets     = ReEmitUniqueAnchor.fusions_parquet,
            output_db    = cohort_db_name,
            docker_image = docker_image,
            preemptible  = preemptible
    }

    output {
        File        cohort_duckdb  = MergeCohort.cohort_duckdb
        Array[File] per_sample_tsv = ReEmitUniqueAnchor.fusions_tsv
    }
}

task ReEmitUniqueAnchor {
    input {
        File    evidence_bam
        File    fusion_index
        File?   reference_fasta
        Int     kmer_size
        Int     min_kmer_hits
        Int     min_supporting_reads
        Int     min_mapq
        Int     min_coherent_fragments
        Int     min_anchor_kmers
        Float   max_breakpoint_std
        Int     min_breakpoint_reads
        Int     min_breakpoint_distance
        Boolean skip_version_check
        String  docker_image
        Int     memory_gb
        Int     disk_gb
        Int     preemptible
    }

    # The evidence BAM stem is the sample id (e.g. "SAMPLE.fusions.bam" -> "SAMPLE").
    String stem            = sub(basename(evidence_bam), "\\.(fusions\\.)?bam$", "")
    String output_parquet  = stem + ".fusions.parquet"
    String output_tsv      = stem + ".fusions.tsv"

    # Deliberately NO --reads-output and NO --kmer-hits-output: those are the large
    # outputs (the evidence BAM is ~1 GB, kmer_hits ~30 GB). We only need the small
    # call table that now carries n_unique_anchored. --emit-unique-anchor turns on
    # uniqueness tracking without applying a filter, so any threshold is applied offline.
    command <<<
        set -euo pipefail

        geac experimental fusions \
            --bam                     ~{evidence_bam} \
            --index                   ~{fusion_index} \
            --output                  ~{output_parquet} \
            --tsv-output              ~{output_tsv} \
            --sample-id               ~{stem} \
            --kmer-size               ~{kmer_size} \
            --min-kmer-hits           ~{min_kmer_hits} \
            --min-supporting-reads    ~{min_supporting_reads} \
            --min-mapq                ~{min_mapq} \
            --min-coherent-fragments  ~{min_coherent_fragments} \
            --min-anchor-kmers        ~{min_anchor_kmers} \
            --max-breakpoint-std      ~{max_breakpoint_std} \
            --min-breakpoint-reads    ~{min_breakpoint_reads} \
            --min-breakpoint-distance ~{min_breakpoint_distance} \
            --emit-unique-anchor \
            ~{"--reference " + reference_fasta} \
            ~{true="--skip-version-check" false="" skip_version_check}
    >>>

    output {
        File fusions_parquet = output_parquet
        File fusions_tsv     = output_tsv
    }

    runtime {
        docker:      docker_image
        memory:      memory_gb + " GB"
        cpu:         1
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

task MergeCohort {
    input {
        Array[File] parquets
        String      output_db
        String      docker_image
        Int         preemptible
    }

    command <<<
        set -euo pipefail

        geac merge \
            --output ~{output_db} \
            ~{sep=" " parquets}
    >>>

    output {
        File cohort_duckdb = output_db
    }

    runtime {
        docker:      docker_image
        memory:      "8 GB"
        cpu:         1
        disks:       "local-disk 20 HDD"
        preemptible: preemptible
    }
}

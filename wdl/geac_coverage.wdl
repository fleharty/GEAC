version 1.0

## geac_coverage.wdl
##
## Scatter `geac coverage` across a cohort of BAM/CRAM files, then merge
## all per-sample coverage Parquet files into a single cohort DuckDB database.
##
## Intended for use on Terra. Provide a sample set where each member has
## bam / bai columns.
##
## Inputs (per-sample, parallel arrays — lengths must match):
##   input_bams            - BAM or CRAM files
##   input_bam_indices     - Corresponding .bai / .crai indices
##   sample_ids            - (optional) override sample IDs; defaults to SM tag per BAM
##   batches               - (optional) per-sample batch/group label stored as a column
##
## Inputs (shared across all samples):
##   reference_fasta       - Reference FASTA
##   reference_fasta_index - Corresponding .fai index
##   read_types            - (optional) per-sample array of duplex|simplex|raw; defaults to "duplex" for all
##   pipelines             - (optional) per-sample array of fgbio|dragen|raw; defaults to "fgbio" for all
##   targets               - (optional) BED or Picard interval list; when provided, zero-depth
##                           positions within targets are included in output (dropout captured)
##   gene_annotations      - (optional) GTF, GFF3, or UCSC genePred (.txt/.txt.gz)
##   region                - (optional) restrict all samples to a genomic region
##   min_map_qual          - minimum mapping quality for total_depth (default 0)
##   min_base_qual         - base quality threshold for frac_low_bq (default 20)
##   gc_window             - bp window centred on position for GC content (default 100)
##   min_depth             - suppress positions with total_depth below this (default 0)
##   bin_size              - aggregate consecutive positions into bins of this size (default 1 = per-position)
##   adaptive_depth_threshold - when set, positions below this depth are emitted at single-base resolution
##   cohort_name           - base name for the output DuckDB file (default: cohort)
##   docker_image          - geac Docker image, e.g. ghcr.io/fleharty/geac:latest
##
## Outputs:
##   coverage_parquets     - Per-sample coverage Parquet files from geac coverage
##   cohort_db             - Merged cohort DuckDB from geac merge

workflow GeacCoverage {

    input {
        # Per-sample parallel arrays
        Array[File]    input_bams
        Array[File]    input_bam_indices
        Array[String]? sample_ids    # optional; if provided must be same length as input_bams
        Array[String]? batches       # optional; per-sample batch/group label

        # Shared inputs
        File   reference_fasta
        File   reference_fasta_index
        Array[String]? read_types   # optional; defaults to "duplex" for all
        Array[String]? pipelines    # optional; defaults to "fgbio" for all

        File?   targets
        File?   gene_annotations
        String? region
        Int     min_map_qual  = 0
        Int     min_base_qual = 20
        Int     gc_window     = 100
        Int     min_depth     = 0
        Int     bin_size      = 1
        Int?    adaptive_depth_threshold

        String cohort_name = "cohort"

        # Resource settings
        String docker_image
        Int    coverage_memory_gb = 8
        Int    coverage_disk_gb   = 100
        Int    merge_memory_gb    = 16
        Int    merge_disk_gb      = 50
        Int    preemptible        = 2
    }

    scatter (i in range(length(input_bams))) {

        if (defined(sample_ids)) {
            String this_sample_id = select_first([sample_ids])[i]
        }
        if (defined(batches)) {
            String this_batch = select_first([batches])[i]
        }
        String this_read_type = if defined(read_types) then select_first([read_types])[i] else "duplex"
        String this_pipeline  = if defined(pipelines)  then select_first([pipelines])[i]  else "fgbio"

        call Coverage {
            input:
                input_bam             = input_bams[i],
                input_bam_index       = input_bam_indices[i],
                reference_fasta       = reference_fasta,
                reference_fasta_index = reference_fasta_index,
                read_type             = this_read_type,
                pipeline              = this_pipeline,
                sample_id             = this_sample_id,
                batch                 = this_batch,
                targets               = targets,
                gene_annotations      = gene_annotations,
                region                = region,
                min_map_qual          = min_map_qual,
                min_base_qual         = min_base_qual,
                gc_window             = gc_window,
                min_depth                 = min_depth,
                bin_size                  = bin_size,
                adaptive_depth_threshold  = adaptive_depth_threshold,
                docker_image              = docker_image,
                memory_gb             = coverage_memory_gb,
                disk_gb               = coverage_disk_gb,
                preemptible           = preemptible,
        }
    }

    call Merge {
        input:
            parquets     = Coverage.coverage_parquet,
            cohort_name  = cohort_name,
            docker_image = docker_image,
            memory_gb    = merge_memory_gb,
            disk_gb      = merge_disk_gb,
            preemptible  = preemptible,
    }

    output {
        Array[File] coverage_parquets = Coverage.coverage_parquet
        File        cohort_db         = Merge.cohort_db
    }
}

# ── Tasks ──────────────────────────────────────────────────────────────────────

task Coverage {

    input {
        File   input_bam
        File   input_bam_index
        File   reference_fasta
        File   reference_fasta_index
        String read_type
        String pipeline

        String? sample_id
        String? batch
        File?   targets
        File?   gene_annotations
        String? region
        Int     min_map_qual
        Int     min_base_qual
        Int     gc_window
        Int     min_depth
        Int     bin_size
        Int?    adaptive_depth_threshold

        String docker_image
        Int    memory_gb
        Int    disk_gb
        Int    preemptible
    }

    String stem   = sub(basename(input_bam), "\\.(bam|cram)$", "")
    String output_name = stem + ".coverage.parquet"

    command <<<
        set -euo pipefail

        geac coverage \
            --input            ~{input_bam} \
            --reference        ~{reference_fasta} \
            --output           ~{output_name} \
            --read-type        ~{read_type} \
            --pipeline         ~{pipeline} \
            --min-map-qual     ~{min_map_qual} \
            --min-base-qual    ~{min_base_qual} \
            --gc-window        ~{gc_window} \
            --min-depth        ~{min_depth} \
            --bin-size         ~{bin_size} \
            ~{"--adaptive-depth-threshold " + adaptive_depth_threshold} \
            ~{"--sample-id "        + sample_id} \
            ~{"--batch "            + batch} \
            ~{"--targets "          + targets} \
            ~{"--gene-annotations " + gene_annotations} \
            ~{"--region "           + region}
    >>>

    output {
        File coverage_parquet = output_name
    }

    runtime {
        docker:      docker_image
        memory:      memory_gb + " GB"
        cpu:         1
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

task Merge {

    input {
        Array[File] parquets
        String      cohort_name

        String docker_image
        Int    memory_gb
        Int    disk_gb
        Int    preemptible
    }

    String output_db = cohort_name + ".duckdb"

    command <<<
        set -euo pipefail

        geac merge \
            --output ~{output_db} \
            ~{sep=" " parquets}
    >>>

    output {
        File cohort_db = output_db
    }

    runtime {
        docker:      docker_image
        memory:      memory_gb + " GB"
        cpu:         2
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

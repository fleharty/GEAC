version 1.0

## geac_cohort_loci.wdl
##
## Run `geac cohort` on a set of per-sample locus Parquet files to identify
## recurrent alt-base loci across the cohort.  A locus is reported when it
## appears in at least `min_samples` samples (and optionally at least
## `min_sample_fraction` of all samples).
##
## This workflow wraps the `geac cohort` subcommand, which is distinct from
## geac_cohort.wdl (that workflow scatters `geac collect` across a sample set
## and is used to produce the per-sample Parquets in the first place).
##
## Inputs:
##   parquets             - Per-sample locus Parquet files (from geac collect)
##   output_name          - Base name for the output file (default: cohort_loci)
##   output_format        - "parquet" or "tsv" — controls output file extension (default: parquet)
##   min_samples          - Minimum number of samples a locus must appear in (default: 2)
##   min_sample_fraction  - Minimum fraction of samples (0.0–1.0, default: 0.0)
##   on_target_only       - Restrict to on-target loci only (requires on_target column)
##   top_n                - Number of top loci to print to stdout (default: 20)
##   docker_image         - geac Docker image, e.g. ghcr.io/fleharty/geac:latest
##   memory_gb            - Memory in GB (default 8)
##   disk_gb              - Disk space in GB (default 50)
##   preemptible          - Number of preemptible retries (default 2)
##
## Outputs:
##   cohort_loci          - Recurrent loci table (Parquet or TSV per output_format)

workflow GeacCohortLoci {

    input {
        Array[File] parquets
        String      output_name         = "cohort_loci"
        String      output_format       = "parquet"   # "parquet" or "tsv"
        Int         min_samples         = 2
        Float       min_sample_fraction = 0.0
        Boolean     on_target_only      = false
        Int         top_n               = 20

        String docker_image
        Int    memory_gb    = 8
        Int    disk_gb      = 50
        Int    preemptible  = 2
    }

    call CohortLoci {
        input:
            parquets             = parquets,
            output_name          = output_name,
            output_format        = output_format,
            min_samples          = min_samples,
            min_sample_fraction  = min_sample_fraction,
            on_target_only       = on_target_only,
            top_n                = top_n,
            docker_image         = docker_image,
            memory_gb            = memory_gb,
            disk_gb              = disk_gb,
            preemptible          = preemptible,
    }

    output {
        File cohort_loci = CohortLoci.cohort_loci
    }
}

# ── Tasks ──────────────────────────────────────────────────────────────────────

task CohortLoci {

    input {
        Array[File] parquets
        String      output_name
        String      output_format
        Int         min_samples
        Float       min_sample_fraction
        Boolean     on_target_only
        Int         top_n

        String docker_image
        Int    memory_gb
        Int    disk_gb
        Int    preemptible
    }

    String output_file = output_name + "." + output_format

    command <<<
        set -euo pipefail

        geac cohort \
            --output                ~{output_file} \
            --min-samples           ~{min_samples} \
            --min-sample-fraction   ~{min_sample_fraction} \
            --top-n                 ~{top_n} \
            ~{if on_target_only then "--on-target-only" else ""} \
            ~{sep=" " parquets}
    >>>

    output {
        File cohort_loci = output_file
    }

    runtime {
        docker:      docker_image
        memory:      memory_gb + " GB"
        cpu:         1
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

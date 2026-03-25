version 1.0

## geac_merge.wdl
##
## Standalone workflow that merges pre-existing per-sample Parquet files
## into a single cohort DuckDB database using `geac merge`.
##
## Use this when per-sample `geac collect` Parquets already exist and you
## only need to (re-)run the merge step, e.g. after adding new samples.
##
## Inputs:
##   parquets     - Per-sample Parquet files produced by `geac collect`
##   cohort_name  - Base name for the output DuckDB file (default: cohort)
##   docker_image - geac Docker image, e.g. ghcr.io/fleharty/geac:latest
##   memory_gb    - Memory in GB for the merge task (default: 16)
##   disk_gb      - Local disk in GB for the merge task (default: 50)
##   preemptible  - Number of preemptible retries (default: 2)
##
## Outputs:
##   cohort_db    - Merged cohort DuckDB database

workflow GeacMerge {

    input {
        Array[File] parquets

        String cohort_name = "cohort"

        String docker_image
        Int    memory_gb   = 16
        Int    disk_gb     = 50
        Int    preemptible = 2
    }

    call Merge {
        input:
            parquets     = parquets,
            cohort_name  = cohort_name,
            docker_image = docker_image,
            memory_gb    = memory_gb,
            disk_gb      = disk_gb,
            preemptible  = preemptible,
    }

    output {
        File cohort_db = Merge.cohort_db
    }
}

# ── Task ───────────────────────────────────────────────────────────────────────

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

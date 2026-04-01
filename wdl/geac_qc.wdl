version 1.0

## geac_qc.wdl
##
## Run `geac qc` on one or more per-sample Parquet files and produce a
## per-sample QC summary report.  Typically the inputs are the locus Parquets
## produced by geac_cohort.wdl or geac_collect.wdl.
##
## Inputs:
##   parquets          - One or more per-sample locus Parquet files
##   on_target_only    - Restrict QC to on-target loci only (requires on_target column)
##   output_name       - Base name for the output TSV (default: geac_qc)
##   docker_image      - geac Docker image, e.g. ghcr.io/fleharty/geac:latest
##   memory_gb         - Memory in GB (default 8)
##   disk_gb           - Disk space in GB (default 50)
##   preemptible       - Number of preemptible retries (default 2)
##
## Outputs:
##   qc_tsv            - Machine-readable TSV with one row per sample
##
## Notes:
##   `geac qc` also prints a human-readable report to stdout; that output is
##   captured in the task's stdout log and visible in Terra's job logs.

workflow GeacQc {

    input {
        Array[File] parquets
        Boolean     on_target_only = false
        String      output_name    = "geac_qc"

        String docker_image
        Int    memory_gb    = 8
        Int    disk_gb      = 50
        Int    preemptible  = 2
    }

    call Qc {
        input:
            parquets       = parquets,
            on_target_only = on_target_only,
            output_name    = output_name,
            docker_image   = docker_image,
            memory_gb      = memory_gb,
            disk_gb        = disk_gb,
            preemptible    = preemptible,
    }

    output {
        File qc_tsv = Qc.qc_tsv
    }
}

# ── Tasks ──────────────────────────────────────────────────────────────────────

task Qc {

    input {
        Array[File] parquets
        Boolean     on_target_only
        String      output_name

        String docker_image
        Int    memory_gb
        Int    disk_gb
        Int    preemptible
    }

    String output_tsv = output_name + ".tsv"

    command <<<
        set -euo pipefail

        geac qc \
            --output ~{output_tsv} \
            ~{if on_target_only then "--on-target-only" else ""} \
            ~{sep=" " parquets}
    >>>

    output {
        File qc_tsv = output_tsv
    }

    runtime {
        docker:      docker_image
        memory:      memory_gb + " GB"
        cpu:         1
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

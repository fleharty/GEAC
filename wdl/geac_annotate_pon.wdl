version 1.0

## geac_annotate_pon.wdl
##
## Run `geac annotate-pon` on a single tumor sample against a pre-built
## Panel of Normals (PoN) DuckDB.  The PoN DuckDB is produced once by running
## `geac collect` on each normal sample and then `geac merge` on the results.
##
## For each alt locus in the tumor Parquet, the subcommand queries how many
## PoN samples carry the same alt allele and at what VAF — all via DuckDB
## analytics, with no BAM pileup required.
##
## Inputs:
##   tumor_parquet  - Locus Parquet from `geac collect` for the tumor sample
##   pon_db         - PoN DuckDB from `geac merge` (must contain alt_bases table)
##   docker_image   - geac Docker image, e.g. ghcr.io/fleharty/geac:0.3.7
##
## Outputs:
##   pon_evidence_parquet - PoN evidence Parquet ({stem}.pon_evidence.parquet)
##                         Pass this to `geac merge` alongside locus Parquets
##                         to create a pon_evidence table in the cohort DuckDB.

workflow GeacAnnotatePon {

    input {
        File   tumor_parquet
        File   pon_db

        String docker_image
        Int    memory_gb   = 4
        Int    disk_gb     = 50
        Int    preemptible = 2
    }

    call AnnotatePon {
        input:
            tumor_parquet = tumor_parquet,
            pon_db        = pon_db,
            docker_image  = docker_image,
            memory_gb     = memory_gb,
            disk_gb       = disk_gb,
            preemptible   = preemptible,
    }

    output {
        File pon_evidence_parquet = AnnotatePon.pon_evidence_parquet
    }
}

# ── Task ────────────────────────────────────────────────────────────────────────

task AnnotatePon {

    input {
        File   tumor_parquet
        File   pon_db

        String docker_image
        Int    memory_gb
        Int    disk_gb
        Int    preemptible
    }

    # Derive output filename: strip .locus.parquet or .parquet, append .pon_evidence.parquet
    String stem        = sub(sub(basename(tumor_parquet), "\\.locus\\.parquet$", ""), "\\.parquet$", "")
    String output_name = stem + ".pon_evidence.parquet"

    command <<<
        set -euo pipefail

        geac annotate-pon \
            --tumor-parquet ~{tumor_parquet} \
            --pon-db        ~{pon_db} \
            --output        ~{output_name}
    >>>

    output {
        File pon_evidence_parquet = output_name
    }

    runtime {
        docker:      docker_image
        memory:      memory_gb + " GB"
        cpu:         1
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

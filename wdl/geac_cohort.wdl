version 1.0

## geac_cohort.wdl
##
## Scatter `geac collect` across a cohort of BAM/CRAM files, then merge
## all per-sample Parquet files into a single cohort DuckDB database.
##
## Intended for use on Terra. Provide a sample set where each member has
## bam / bai columns, plus optional per-sample variant TSV files.
##
## Inputs (per-sample, parallel arrays — lengths must match):
##   input_bams            - BAM or CRAM files
##   input_bam_indices     - Corresponding .bai / .crai indices
##   sample_ids            - (optional) override sample IDs; defaults to SM tag per BAM
##   variants_tsvs         - (optional) per-sample variant TSV files
##   vcfs                  - (optional) per-sample VCF/BCF files for variant annotation
##   vcf_indices           - (optional) per-sample .tbi / .csi indices; required when vcfs provided
##
## Inputs (shared across all samples):
##   reference_fasta       - Reference FASTA
##   reference_fasta_index - Corresponding .fai index
##   read_types            - (optional) per-sample array of duplex|simplex|raw; defaults to "duplex" for all
##   pipelines             - (optional) per-sample array of fgbio|dragen|raw; defaults to "fgbio" for all
##   targets               - (optional) BED or Picard interval list
##   gene_annotations      - (optional) GTF, GFF3, or UCSC genePred (.txt/.txt.gz)
##   region                - (optional) restrict all samples to a genomic region
##   repeat_window         - bases each side of locus for homopolymer/STR scan (default 10)
##   cohort_name           - Base name for the output DuckDB file (default: cohort)
##   docker_image          - geac Docker image, e.g. gcr.io/my-project/geac:latest
##
## Outputs:
##   parquets              - Per-sample Parquet files from geac collect
##   cohort_db             - Merged cohort DuckDB from geac merge

workflow GeacCohort {

    input {
        # Per-sample parallel arrays
        Array[File]    input_bams
        Array[File]    input_bam_indices
        Array[String]? sample_ids       # optional; if provided must be same length as input_bams
        Array[File]?   variants_tsvs    # optional; if provided must be same length as input_bams
        Array[File]?   vcfs             # optional; if provided must be same length as input_bams
        Array[File]?   vcf_indices      # optional; required when vcfs is provided

        # Shared inputs
        File   reference_fasta
        File   reference_fasta_index
        Array[String]? read_types      # optional; if provided must be same length as input_bams
        Array[String]? pipelines       # optional; if provided must be same length as input_bams

        File?   targets
        File?   gene_annotations
        String? region
        Int     repeat_window = 10

        Int min_base_qual = 1
        Int min_map_qual  = 0
        Int threads       = 4

        String cohort_name = "cohort"

        # Resource settings
        String docker_image
        Int    collect_memory_gb = 8
        Int    collect_disk_gb   = 100
        Int    merge_memory_gb   = 16
        Int    merge_disk_gb     = 50
        Int    preemptible       = 2
    }

    scatter (i in range(length(input_bams))) {

        # Safely index into optional per-sample arrays.
        # Variables declared inside `if` blocks are typed as optional (String?, File?)
        # in the enclosing scope, which is exactly what the Collect task expects.
        if (defined(sample_ids)) {
            String this_sample_id    = select_first([sample_ids])[i]
        }
        if (defined(variants_tsvs)) {
            File   this_variants_tsv = select_first([variants_tsvs])[i]
        }
        if (defined(vcfs)) {
            File   this_vcf          = select_first([vcfs])[i]
        }
        if (defined(vcf_indices)) {
            File   this_vcf_index    = select_first([vcf_indices])[i]
        }
        String this_read_type = if defined(read_types) then select_first([read_types])[i] else "duplex"
        String this_pipeline  = if defined(pipelines)  then select_first([pipelines])[i]  else "fgbio"

        call Collect {
            input:
                input_bam             = input_bams[i],
                input_bam_index       = input_bam_indices[i],
                reference_fasta       = reference_fasta,
                reference_fasta_index = reference_fasta_index,
                read_type             = this_read_type,
                pipeline              = this_pipeline,
                sample_id             = this_sample_id,
                variants_tsv          = this_variants_tsv,
                vcf                   = this_vcf,
                vcf_index             = this_vcf_index,
                targets               = targets,
                gene_annotations      = gene_annotations,
                region                = region,
                repeat_window         = repeat_window,
                min_base_qual         = min_base_qual,
                min_map_qual          = min_map_qual,
                threads               = threads,
                docker_image          = docker_image,
                memory_gb             = collect_memory_gb,
                disk_gb               = collect_disk_gb,
                preemptible           = preemptible,
        }
    }

    call Merge {
        input:
            parquets     = Collect.parquet,
            cohort_name  = cohort_name,
            docker_image = docker_image,
            memory_gb    = merge_memory_gb,
            disk_gb      = merge_disk_gb,
            preemptible  = preemptible,
    }

    output {
        Array[File] parquets  = Collect.parquet
        File        cohort_db = Merge.cohort_db
    }
}

# ── Tasks ──────────────────────────────────────────────────────────────────────

task Collect {

    input {
        File   input_bam
        File   input_bam_index
        File   reference_fasta
        File   reference_fasta_index
        String read_type
        String pipeline

        String? sample_id
        File?   variants_tsv
        File?   vcf
        File?   vcf_index
        File?   targets
        File?   gene_annotations
        String? region
        Int     repeat_window

        Int min_base_qual
        Int min_map_qual
        Int threads

        String docker_image
        Int    memory_gb
        Int    disk_gb
        Int    preemptible
    }

    String output_name = sub(basename(input_bam), "\\.(bam|cram)$", "") + ".parquet"

    command <<<
        set -euo pipefail

        geac collect \
            --input            ~{input_bam} \
            --reference        ~{reference_fasta} \
            --output           ~{output_name} \
            --read-type        ~{read_type} \
            --pipeline         ~{pipeline} \
            --min-base-qual    ~{min_base_qual} \
            --min-map-qual     ~{min_map_qual} \
            --threads          ~{threads} \
            ~{"--sample-id "        + sample_id} \
            ~{"--vcf "              + vcf} \
            ~{"--variants-tsv "     + variants_tsv} \
            ~{"--targets "          + targets} \
            ~{"--gene-annotations " + gene_annotations} \
            ~{"--region "           + region} \
            --repeat-window ~{repeat_window}
    >>>

    output {
        File parquet = output_name
    }

    runtime {
        docker:      docker_image
        memory:      memory_gb + " GB"
        cpu:         threads
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

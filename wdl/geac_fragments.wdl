version 1.0

## geac_fragments.wdl
##
## Scatter `geac fragments` across a cohort of BAM/CRAM files to produce
## per-sample fragment-level Parquet files (insert size, GC content, end motifs).
##
## Intended for use on Terra. Provide a sample set where each member has
## bam / bai columns.
##
## Fragment Parquets are large for WGS (hundreds of millions of rows per sample).
## Use --region or --targets to restrict to a region of interest for interactive
## exploratory analysis. For standard whole-genome collection, increase disk_gb.
##
## Inputs (per-sample parallel arrays — lengths must match):
##   input_bams            - BAM or CRAM files
##   input_bam_indices     - Corresponding .bai / .crai indices
##   sample_ids            - (optional) override sample IDs; defaults to SM tag per BAM
##   batches               - (optional) per-sample batch/group label
##
## Inputs (shared across all samples):
##   reference_fasta       - Reference FASTA
##   reference_fasta_index - Corresponding .fai index
##   read_types            - (optional) per-sample array of duplex|simplex|raw; defaults to "duplex"
##   pipelines             - (optional) per-sample array of fgbio|dragen|raw; defaults to "fgbio"
##   region                - (optional) restrict all samples to a genomic region string or BED file
##   min_map_qual          - minimum mapping quality (default 0)
##   docker_image          - geac Docker image, e.g. ghcr.io/fleharty/geac:latest
##
## Outputs:
##   fragments_parquets    - Per-sample fragment Parquet files

workflow GeacFragments {

    input {
        # Per-sample parallel arrays
        Array[File]    input_bams
        Array[File]    input_bam_indices
        Array[String]? sample_ids
        Array[String]? batches
        Array[String]? labels1
        Array[String]? labels2
        Array[String]? labels3
        Array[String]? timepoints

        # Shared inputs
        File   reference_fasta
        File   reference_fasta_index
        Array[String]? read_types
        Array[String]? pipelines

        String? region
        Int     min_map_qual = 0

        # Resource settings
        String docker_image
        Int    memory_gb    = 8
        Int    disk_gb      = 200
        Int    preemptible  = 2
    }

    scatter (i in range(length(input_bams))) {

        if (defined(sample_ids)) {
            String this_sample_id = select_first([sample_ids])[i]
        }
        if (defined(batches)) {
            String this_batch = select_first([batches])[i]
        }
        if (defined(labels1)) {
            String this_label1 = select_first([labels1])[i]
        }
        if (defined(labels2)) {
            String this_label2 = select_first([labels2])[i]
        }
        if (defined(labels3)) {
            String this_label3 = select_first([labels3])[i]
        }
        if (defined(timepoints)) {
            String this_timepoint = select_first([timepoints])[i]
        }
        String this_read_type = if defined(read_types) then select_first([read_types])[i] else "duplex"
        String this_pipeline  = if defined(pipelines)  then select_first([pipelines])[i]  else "fgbio"

        call Fragments {
            input:
                input_bam             = input_bams[i],
                input_bam_index       = input_bam_indices[i],
                reference_fasta       = reference_fasta,
                reference_fasta_index = reference_fasta_index,
                read_type             = this_read_type,
                pipeline              = this_pipeline,
                sample_id             = this_sample_id,
                batch                 = this_batch,
                label1                = this_label1,
                label2                = this_label2,
                label3                = this_label3,
                timepoint             = this_timepoint,
                region                = region,
                min_map_qual          = min_map_qual,
                docker_image          = docker_image,
                memory_gb             = memory_gb,
                disk_gb               = disk_gb,
                preemptible           = preemptible,
        }
    }

    output {
        Array[File] fragments_parquets = Fragments.fragments_parquet
    }
}

# ── Tasks ──────────────────────────────────────────────────────────────────────

task Fragments {

    input {
        File   input_bam
        File   input_bam_index
        File   reference_fasta
        File   reference_fasta_index
        String read_type
        String pipeline

        String? sample_id
        String? batch
        String? label1
        String? label2
        String? label3
        String? timepoint
        String? region
        Int     min_map_qual

        String docker_image
        Int    memory_gb
        Int    disk_gb
        Int    preemptible
    }

    String stem         = sub(basename(input_bam), "\\.(bam|cram)$", "")
    String output_name  = stem + ".fragments.parquet"

    command <<<
        set -euo pipefail

        geac fragments \
            --input            ~{input_bam} \
            --reference        ~{reference_fasta} \
            --output           ~{output_name} \
            --read-type        ~{read_type} \
            --pipeline         ~{pipeline} \
            --min-map-qual     ~{min_map_qual} \
            ~{"--sample-id "   + sample_id} \
            ~{"--batch "       + batch} \
            ~{"--label1 "      + label1} \
            ~{"--label2 "      + label2} \
            ~{"--label3 "      + label3} \
            ~{"--timepoint "   + timepoint} \
            ~{"--region "      + region}
    >>>

    output {
        File fragments_parquet = output_name
    }

    runtime {
        docker:      docker_image
        memory:      memory_gb + " GB"
        cpu:         1
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

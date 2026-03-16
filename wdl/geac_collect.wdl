version 1.0

## geac_collect.wdl
##
## Run `geac collect` on a single BAM/CRAM sample and produce a Parquet file
## of alt base metrics. Designed for use in Terra; scatter over a sample set
## to process a full cohort in parallel, then pass all Parquets to geac_merge.
##
## Inputs:
##   input_bam            - BAM or CRAM file
##   input_bam_index      - Corresponding .bai / .crai index
##   reference_fasta      - Reference FASTA
##   reference_fasta_index - Corresponding .fai index
##   read_type            - duplex | simplex | raw
##   pipeline             - fgbio | dragen | raw
##   sample_id            - (optional) override sample ID; defaults to BAM SM tag
##   vcf                  - (optional) VCF/BCF for variant call annotation
##   vcf_index            - (optional) Corresponding .tbi / .csi index
##   region               - (optional) restrict to a region, e.g. chr1:1-1000000
##   min_base_qual        - minimum base quality (default 1)
##   min_map_qual         - minimum mapping quality (default 20)
##   threads              - CPU threads for geac (default 4)
##   docker_image         - geac Docker image, e.g. gcr.io/my-project/geac:0.1.0
##   memory_gb            - memory in GB (default 8)
##   disk_gb              - disk space in GB (default 100)
##   preemptible          - number of preemptible retries (default 2)

workflow GeacCollect {

    input {
        File   input_bam
        File   input_bam_index
        File   reference_fasta
        File   reference_fasta_index
        String read_type
        String pipeline

        String? sample_id
        File?   vcf
        File?   vcf_index
        String? region

        Int min_base_qual = 1
        Int min_map_qual  = 20
        Int threads       = 4

        String docker_image
        Int    memory_gb    = 8
        Int    disk_gb      = 100
        Int    preemptible  = 2
    }

    call Collect {
        input:
            input_bam             = input_bam,
            input_bam_index       = input_bam_index,
            reference_fasta       = reference_fasta,
            reference_fasta_index = reference_fasta_index,
            read_type             = read_type,
            pipeline              = pipeline,
            sample_id             = sample_id,
            vcf                   = vcf,
            vcf_index             = vcf_index,
            region                = region,
            min_base_qual         = min_base_qual,
            min_map_qual          = min_map_qual,
            threads               = threads,
            docker_image          = docker_image,
            memory_gb             = memory_gb,
            disk_gb               = disk_gb,
            preemptible           = preemptible,
    }

    output {
        File parquet = Collect.parquet
    }
}

task Collect {

    input {
        File   input_bam
        File   input_bam_index
        File   reference_fasta
        File   reference_fasta_index
        String read_type
        String pipeline

        String? sample_id
        File?   vcf
        File?   vcf_index
        String? region

        Int min_base_qual
        Int min_map_qual
        Int threads

        String docker_image
        Int    memory_gb
        Int    disk_gb
        Int    preemptible
    }

    # Derive output filename from the input BAM/CRAM basename
    String output_name = sub(basename(input_bam), "\\.(bam|cram)$", "") + ".parquet"

    command <<<
        set -euo pipefail

        geac collect \
            --input ~{input_bam} \
            --reference ~{reference_fasta} \
            --output ~{output_name} \
            --read-type ~{read_type} \
            --pipeline ~{pipeline} \
            --min-base-qual ~{min_base_qual} \
            --min-map-qual ~{min_map_qual} \
            --threads ~{threads} \
            ~{"--sample-id " + sample_id} \
            ~{"--vcf " + vcf} \
            ~{"--region " + region}
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

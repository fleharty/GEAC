version 1.0

## geac_annotate_normal.wdl
##
## Run `geac annotate-normal` on a single tumor/normal pair.  Reads the tumor
## locus Parquet (produced by `geac collect`) and pileups the paired normal
## BAM/CRAM at each tumor alt-base locus, writing a `.normal_evidence.parquet`
## file that can be passed to `geac merge` alongside the locus Parquets.
##
## Inputs:
##   tumor_parquet         - Locus Parquet from `geac collect` for the tumor sample
##   normal_bam            - Normal BAM or CRAM file
##   normal_bam_index      - Corresponding .bai / .crai index
##   reference_fasta       - Reference FASTA
##   reference_fasta_index - Corresponding .fai index
##   normal_sample_id      - (optional) override normal sample ID; defaults to SM tag
##   min_base_qual         - minimum base quality (default 1)
##   min_map_qual          - minimum mapping quality (default 0)
##   include_duplicates    - include PCR/optical duplicate reads (FLAG 0x400); default false
##   include_secondary     - include secondary alignments (FLAG 0x100); default false
##   include_supplementary - include supplementary alignments (FLAG 0x800); default false
##   docker_image          - geac Docker image, e.g. ghcr.io/fleharty/geac:0.3.7
##
## Outputs:
##   normal_evidence_parquet - Normal evidence Parquet ({stem}.normal_evidence.parquet)

workflow GeacAnnotateNormal {

    input {
        File   tumor_parquet
        File   normal_bam
        File   normal_bam_index
        File   reference_fasta
        File   reference_fasta_index

        String? normal_sample_id

        Int     min_base_qual = 1
        Int     min_map_qual  = 0
        Boolean include_duplicates    = false
        Boolean include_secondary     = false
        Boolean include_supplementary = false

        String docker_image
        Int    memory_gb    = 8
        Int    disk_gb      = 100
        Int    preemptible  = 2
    }

    call AnnotateNormal {
        input:
            tumor_parquet         = tumor_parquet,
            normal_bam            = normal_bam,
            normal_bam_index      = normal_bam_index,
            reference_fasta       = reference_fasta,
            reference_fasta_index = reference_fasta_index,
            normal_sample_id      = normal_sample_id,
            min_base_qual         = min_base_qual,
            min_map_qual          = min_map_qual,
            include_duplicates    = include_duplicates,
            include_secondary     = include_secondary,
            include_supplementary = include_supplementary,
            docker_image          = docker_image,
            memory_gb             = memory_gb,
            disk_gb               = disk_gb,
            preemptible           = preemptible,
    }

    output {
        File normal_evidence_parquet = AnnotateNormal.normal_evidence_parquet
    }
}

# ── Task ────────────────────────────────────────────────────────────────────────

task AnnotateNormal {

    input {
        File   tumor_parquet
        File   normal_bam
        File   normal_bam_index
        File   reference_fasta
        File   reference_fasta_index

        String? normal_sample_id

        Int     min_base_qual
        Int     min_map_qual
        Boolean include_duplicates
        Boolean include_secondary
        Boolean include_supplementary

        String docker_image
        Int    memory_gb
        Int    disk_gb
        Int    preemptible
    }

    # Derive output filename: replace .parquet (or .locus.parquet) with .normal_evidence.parquet
    String stem   = sub(sub(basename(tumor_parquet), "\\.locus\\.parquet$", ""), "\\.parquet$", "")
    String output_name = stem + ".normal_evidence.parquet"

    command <<<
        set -euo pipefail

        geac annotate-normal \
            --tumor-parquet    ~{tumor_parquet} \
            --normal-bam       ~{normal_bam} \
            --reference        ~{reference_fasta} \
            --output           ~{output_name} \
            --min-base-qual    ~{min_base_qual} \
            --min-map-qual     ~{min_map_qual} \
            ~{"--normal-sample-id " + normal_sample_id} \
            ~{if include_duplicates    then "--include-duplicates"    else ""} \
            ~{if include_secondary     then "--include-secondary"     else ""} \
            ~{if include_supplementary then "--include-supplementary" else ""}
    >>>

    output {
        File normal_evidence_parquet = output_name
    }

    runtime {
        docker:      docker_image
        memory:      memory_gb + " GB"
        cpu:         1
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

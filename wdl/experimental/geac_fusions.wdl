version 1.0

## geac_fusions.wdl  [EXPERIMENTAL]
##
## Detect gene fusions from a single BAM/CRAM file using a pre-built k-mer index.
## Scatter over a sample set on Terra to process a cohort in parallel.
##
## The k-mer index is built once with geac_build_fusion_index.wdl and shared
## across all samples.
##
## NOTE: This workflow wraps an experimental GEAC subcommand whose API and outputs
## may change without notice. Do not rely on it in production pipelines.
##
## Inputs:
##   input_bam             - BAM or CRAM file
##   input_bam_index       - Corresponding .bai / .crai index
##   fusion_index          - DuckDB k-mer index from geac_build_fusion_index.wdl
##   reference_fasta       - (optional) Reference FASTA; required for CRAM input
##   reference_fasta_index - (optional) Corresponding .fai index
##   sample_id             - (optional) sample identifier; defaults to SM tag in BAM header
##   kmer_size             - K-mer length; must match the index (default: 23)
##   min_kmer_hits         - Minimum unique k-mer hits to a gene for a read to be assigned (default: 3)
##   min_supporting_reads  - Minimum fragment pairs to report a fusion (default: 2)
##   min_mapq              - Minimum mapping quality; unmapped reads bypass this filter (default: 0)
##   max_kmer_copies       - (optional) ignore k-mers seen more than this many times genome-wide;
##                           requires an index built with check_genome_uniqueness=true
##   docker_image          - geac Docker image, e.g. ghcr.io/fleharty/geac:latest
##   memory_gb             - Memory in GB (default: 32; increase to 64 for high-coverage WGS)
##   disk_gb               - Disk space in GB (default: 100)
##   preemptible           - Preemptible retries (default: 2)
##
## Outputs:
##   fusions_parquet       - Per-fusion detection results Parquet
##   reads_bam             - BAM of all reads from fusion-supporting fragments (empty BAM when no fusions)
##   fusions_tsv           - Human-readable fusion results TSV (header-only when no fusions)
##   kmer_hits_tsv         - Per-k-mer hits TSV for fusion-supporting reads (header-only when no fusions)

workflow GeacFusions {

    input {
        File  input_bam
        File  input_bam_index
        File  fusion_index
        File? reference_fasta
        File? reference_fasta_index

        String? sample_id
        Int     kmer_size            = 23
        Int     min_kmer_hits        = 3
        Int     min_supporting_reads = 2
        Int     min_mapq             = 0
        Int?    max_kmer_copies

        String docker_image
        Int    memory_gb   = 32
        Int    disk_gb     = 100
        Int    preemptible = 2
    }

    call Fusions {
        input:
            input_bam             = input_bam,
            input_bam_index       = input_bam_index,
            fusion_index          = fusion_index,
            reference_fasta       = reference_fasta,
            reference_fasta_index = reference_fasta_index,
            sample_id             = sample_id,
            kmer_size             = kmer_size,
            min_kmer_hits         = min_kmer_hits,
            min_supporting_reads  = min_supporting_reads,
            min_mapq              = min_mapq,
            max_kmer_copies       = max_kmer_copies,
            docker_image          = docker_image,
            memory_gb             = memory_gb,
            disk_gb               = disk_gb,
            preemptible           = preemptible,
    }

    output {
        File fusions_parquet = Fusions.fusions_parquet
        File reads_bam       = Fusions.reads_bam
        File fusions_tsv     = Fusions.fusions_tsv
        File kmer_hits_tsv   = Fusions.kmer_hits_tsv
    }
}

task Fusions {

    input {
        File  input_bam
        File  input_bam_index
        File  fusion_index
        File? reference_fasta
        File? reference_fasta_index

        String? sample_id
        Int     kmer_size
        Int     min_kmer_hits
        Int     min_supporting_reads
        Int     min_mapq
        Int?    max_kmer_copies

        String docker_image
        Int    memory_gb
        Int    disk_gb
        Int    preemptible
    }

    String stem             = sub(basename(input_bam), "\\.(bam|cram)$", "")
    String output_parquet   = stem + ".fusions.parquet"
    String output_reads_bam = stem + ".fusions.bam"
    String output_tsv       = stem + ".fusions.tsv"
    String output_kmer_hits = stem + ".fusions.kmer_hits.tsv"

    command <<<
        set -euo pipefail

        geac experimental fusions \
            --bam                  ~{input_bam} \
            --index                ~{fusion_index} \
            --output               ~{output_parquet} \
            --kmer-size            ~{kmer_size} \
            --min-kmer-hits        ~{min_kmer_hits} \
            --min-supporting-reads ~{min_supporting_reads} \
            --min-mapq             ~{min_mapq} \
            --reads-output         ~{output_reads_bam} \
            --tsv-output           ~{output_tsv} \
            --kmer-hits-output     ~{output_kmer_hits} \
            ~{"--reference "       + reference_fasta} \
            ~{"--sample-id "       + sample_id} \
            ~{"--max-kmer-copies " + max_kmer_copies}
    >>>

    output {
        File fusions_parquet = output_parquet
        File reads_bam       = output_reads_bam
        File fusions_tsv     = output_tsv
        File kmer_hits_tsv   = output_kmer_hits
    }

    runtime {
        docker:      docker_image
        memory:      memory_gb + " GB"
        cpu:         1
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

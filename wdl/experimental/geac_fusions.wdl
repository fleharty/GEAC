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
##   fusion_pon            - (optional) fusion Panel-of-Normals DuckDB (a geac merge of normal-sample
##                           *.fusions.parquet files); annotates each call with PoN columns
##   max_pon_samples       - (optional) drop fusions seen in > this many PoN samples; requires fusion_pon
##   min_coherent_fragments - require at least this many fragments (qnames) with a coherent A-block→B-block
##                           k-mer partition; 0 = disabled (default); set to 1+ to filter homology artifacts
##   min_anchor_kmers      - minimum k-mer hits from each gene for a spanning read to count as anchored
##                           (default: 3; only relevant when min_coherent_fragments > 0)
##   docker_image          - geac Docker image, e.g. ghcr.io/fleharty/geac:latest
##   memory_gb             - Memory in GB (default: 32; increase to 64 for high-coverage WGS)
##   disk_gb               - Disk space in GB (default: 100)
##   preemptible           - Preemptible retries (default: 2)
##
## Outputs:
##   fusions_parquet       - Per-fusion detection results Parquet
##   reads_bam             - BAM of all reads from fusion-supporting fragments (empty BAM when no fusions)
##   reads_bam_index       - Corresponding .bai index
##   fusions_tsv           - Human-readable fusion results TSV (header-only when no fusions)
##   kmer_hits_tsv         - Per-k-mer hits TSV for fusion-supporting reads (header-only when no fusions)
##   breakpoints_tsv       - Per-fusion breakpoint coordinates derived from junction-spanning reads
##   geac_version          - geac binary version string (e.g. "0.4.33") for provenance in the data table

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
        File?   fusion_pon
        Int?    max_pon_samples
        File?   fusion_kmer_blacklist
        Int?    min_kmer_blacklist_samples
        Int     min_coherent_fragments = 0
        Int     min_anchor_kmers       = 3

        String docker_image
        Int    memory_gb   = 32
        Int    disk_gb     = 100
        Int    preemptible = 2
    }

    call Fusions {
        input:
            input_bam                  = input_bam,
            input_bam_index            = input_bam_index,
            fusion_index               = fusion_index,
            reference_fasta            = reference_fasta,
            reference_fasta_index      = reference_fasta_index,
            sample_id                  = sample_id,
            kmer_size                  = kmer_size,
            min_kmer_hits              = min_kmer_hits,
            min_supporting_reads       = min_supporting_reads,
            min_mapq                   = min_mapq,
            max_kmer_copies            = max_kmer_copies,
            fusion_pon                 = fusion_pon,
            max_pon_samples            = max_pon_samples,
            fusion_kmer_blacklist      = fusion_kmer_blacklist,
            min_kmer_blacklist_samples = min_kmer_blacklist_samples,
            min_coherent_fragments     = min_coherent_fragments,
            min_anchor_kmers           = min_anchor_kmers,
            docker_image           = docker_image,
            memory_gb              = memory_gb,
            disk_gb                = disk_gb,
            preemptible            = preemptible,
    }

    output {
        File   fusions_parquet = Fusions.fusions_parquet
        File   reads_bam       = Fusions.reads_bam
        File   reads_bam_index = Fusions.reads_bam_index
        File   fusions_tsv     = Fusions.fusions_tsv
        File   kmer_hits_tsv   = Fusions.kmer_hits_tsv
        File   breakpoints_tsv = Fusions.breakpoints_tsv
        String geac_version    = Fusions.geac_version
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
        File?   fusion_pon
        Int?    max_pon_samples
        File?   fusion_kmer_blacklist
        Int?    min_kmer_blacklist_samples
        Int     min_coherent_fragments
        Int     min_anchor_kmers

        String docker_image
        Int    memory_gb
        Int    disk_gb
        Int    preemptible
    }

    String stem                  = sub(basename(input_bam), "\\.(bam|cram)$", "")
    String output_parquet        = stem + ".fusions.parquet"
    String output_reads_bam      = stem + ".fusions.bam"
    String output_reads_bam_bai  = stem + ".fusions.bam.bai"
    String output_tsv            = stem + ".fusions.tsv"
    String output_kmer_hits      = stem + ".fusions.kmer_hits.tsv"
    String output_breakpoints    = stem + ".fusions.breakpoints.tsv"

    command <<<
        set -euo pipefail

        # Record the geac version (e.g. "geac 0.4.33" -> "0.4.33") for the Terra
        # data table, so each run's output is traceable to the binary that produced it.
        geac --version | awk '{print $2}' > geac_version.txt

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
            --breakpoints-output   ~{output_breakpoints} \
            ~{"--reference "             + reference_fasta} \
            ~{"--sample-id "             + sample_id} \
            ~{"--max-kmer-copies "       + max_kmer_copies} \
            ~{"--fusion-pon "                  + fusion_pon} \
            ~{"--max-pon-samples "             + max_pon_samples} \
            ~{"--fusion-kmer-blacklist "       + fusion_kmer_blacklist} \
            ~{"--min-kmer-blacklist-samples "  + min_kmer_blacklist_samples} \
            --min-coherent-fragments ~{min_coherent_fragments} \
            --min-anchor-kmers       ~{min_anchor_kmers}
    >>>

    output {
        File   fusions_parquet  = output_parquet
        File   reads_bam        = output_reads_bam
        File   reads_bam_index  = output_reads_bam_bai
        File   fusions_tsv      = output_tsv
        File   kmer_hits_tsv    = output_kmer_hits
        File   breakpoints_tsv  = output_breakpoints
        String geac_version     = read_string("geac_version.txt")
    }

    runtime {
        docker:      docker_image
        memory:      memory_gb + " GB"
        cpu:         1
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

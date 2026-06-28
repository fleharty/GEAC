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
##   max_breakpoint_std    - (optional) tag fusions with filter="chimera" unless BOTH breakpoints are
##                           supported by >= min_breakpoint_reads reads whose position estimates have a
##                           standard deviation <= this many bp; rows are kept, not dropped. Recommended: 100.
##                           A "strong support" tier (>= 25 reads per side, std <= 250 bp) also PASSes: real
##                           high-depth junctions spread over tens-to->100 bp (splice isoforms + k-mer
##                           transition-point noise), and low-support artifacts never reach that tier.
##   min_breakpoint_reads  - minimum reads converging on EACH breakpoint under the
##                           max_breakpoint_std filter (default: 5; only consulted when max_breakpoint_std set)
##   min_breakpoint_distance - (optional) tag fusions with filter="samelocus" when both breakpoints fall on the
##                           same chromosome within this many bp (single-locus paralog leakage, e.g. unindexed
##                           GNA12 reads split between GNA13/GNA11); rows are kept. Recommended: 10000.
##                           Calls with strong INDEPENDENT (concordant, non-spanning) support on both partners
##                           are exempt — leakage is spanning-read dominated at one locus, whereas a real
##                           junction observed only through concordant pairs (breakpoint buried in a repeat) is not.
##   asym_anchor_reads     - (optional) single-sample low-input rescue: also keep filter="PASS" when the
##                           DOMINANT breakpoint has >= this many reads converging within asym_anchor_std bp,
##                           the partner breakpoint has >= 1 read, and min_mapq >= asym_anchor_mapq. Recovers a
##                           junction seen well on one side but by 1-2 reads on the other. Only consulted when
##                           max_breakpoint_std is set. Suggested: 10.
##   asym_anchor_std       - max breakpoint std (bp) on the dominant side for asym_anchor_reads (default: 25)
##   asym_anchor_mapq      - min call min_mapq for asym_anchor_reads; excludes multi-mapper leakage (default: 20)
##   min_unique_anchor_reads - (optional) copy-aware specificity filter: keep filter="PASS" only when at least
##                           this many supporting fragments are anchored by a genome-unique k-mer on the
##                           higher-uniqueness partner; else tag filter="no_unique_anchor". Suppresses
##                           repeat/segdup/pseudogene-partner artifacts. Requires an index built with
##                           --check-genome-uniqueness (genome_copies column). Default: 0 (disabled).
##   emit_unique_anchor    - (optional) emit the per-call n_unique_anchored count without filtering, so any
##                           unique-anchor threshold can be applied offline (n_unique_anchored /
##                           supporting_reads) from a single run. Implied when min_unique_anchor_reads > 0.
##                           Requires a --check-genome-uniqueness index. Default: false.
##   comment_text          - (optional) free-text note echoed straight to the data table (as the
##                           `comment` output); not processed by geac, never written to the Parquet/TSV
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
##   comment               - the free-text note passed in via comment_text (if any), as a data-table column

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
        File?    fusion_kmer_blacklist
        Int?     min_kmer_blacklist_samples
        Boolean? skip_version_check
        Int      min_coherent_fragments = 0
        Int      min_anchor_kmers       = 3
        Float?   max_breakpoint_std
        Int      min_breakpoint_reads   = 5
        Int?     min_breakpoint_distance
        Int?     asym_anchor_reads
        Float    asym_anchor_std        = 25.0
        Int      asym_anchor_mapq       = 20
        Int      min_unique_anchor_reads = 0
        Boolean  emit_unique_anchor      = false

        String?  comment_text

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
            skip_version_check         = skip_version_check,
            min_coherent_fragments     = min_coherent_fragments,
            min_anchor_kmers           = min_anchor_kmers,
            max_breakpoint_std         = max_breakpoint_std,
            min_breakpoint_reads       = min_breakpoint_reads,
            min_breakpoint_distance    = min_breakpoint_distance,
            asym_anchor_reads          = asym_anchor_reads,
            asym_anchor_std            = asym_anchor_std,
            asym_anchor_mapq           = asym_anchor_mapq,
            min_unique_anchor_reads    = min_unique_anchor_reads,
            emit_unique_anchor         = emit_unique_anchor,
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
        String? comment        = comment_text
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
        File?    fusion_kmer_blacklist
        Int?     min_kmer_blacklist_samples
        Boolean? skip_version_check
        Int      min_coherent_fragments
        Int      min_anchor_kmers
        Float?   max_breakpoint_std
        Int      min_breakpoint_reads
        Int?     min_breakpoint_distance
        Int?     asym_anchor_reads
        Float    asym_anchor_std
        Int      asym_anchor_mapq
        Int      min_unique_anchor_reads
        Boolean  emit_unique_anchor

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
            ~{true="--skip-version-check" false="" skip_version_check} \
            --min-coherent-fragments ~{min_coherent_fragments} \
            --min-anchor-kmers       ~{min_anchor_kmers} \
            ~{"--max-breakpoint-std "  + max_breakpoint_std} \
            --min-breakpoint-reads   ~{min_breakpoint_reads} \
            ~{"--min-breakpoint-distance " + min_breakpoint_distance} \
            ~{"--asym-anchor-reads " + asym_anchor_reads} \
            --asym-anchor-std        ~{asym_anchor_std} \
            --asym-anchor-mapq       ~{asym_anchor_mapq} \
            --min-unique-anchor-reads ~{min_unique_anchor_reads} \
            ~{true="--emit-unique-anchor" false="" emit_unique_anchor}
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

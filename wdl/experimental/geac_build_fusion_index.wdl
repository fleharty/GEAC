version 1.0

## geac_build_fusion_index.wdl  [EXPERIMENTAL]
##
## Build a gene-unique k-mer index from a reference genome and gene annotation.
## The resulting DuckDB index is used by geac_fusions.wdl to detect gene fusions.
## Run once per reference genome or gene panel; the index is shared across samples.
##
## NOTE: This workflow wraps an experimental GEAC subcommand whose API and outputs
## may change without notice. Do not rely on it in production pipelines.
##
## Inputs:
##   gene_annotation         - Gene annotation file. Accepted formats:
##                               GTF (.gtf, .gtf.gz)      — GENCODE, NCBI RefSeq, or UCSC GTF
##                               GenePred (.txt, .txt.gz) — UCSC genePredExt+bin (e.g. ncbiRefSeq.txt.gz),
##                                                          genePredExt without bin, or refFlat.
##                             Format is auto-detected from file extension and column layout.
##   fasta                   - Reference FASTA (must have a .fai index alongside it)
##   fasta_index             - Corresponding .fai index (localized alongside fasta by Cromwell)
##   index_name              - Base name for the output DuckDB index (default: fusion_index)
##   kmer_size               - K-mer length; must match geac_fusions.wdl (default: 23)
##   min_gene_kmers          - Drop genes with fewer unique k-mers than this (default: 1)
##   genes                   - (optional) text file listing gene names to index, one per line;
##                             blank lines and lines starting with '#' are ignored;
##                             when provided only those genes are indexed (e.g. a fusion panel)
##   check_genome_uniqueness - Scan full genome and remove k-mers that appear more than once
##                             genome-wide (default: false). Requires tens of GB of RAM.
##                             At k >= 19 cross-gene deduplication already eliminates nearly
##                             all non-unique k-mers; use this only for the strictest guarantee.
##   max_genome_copies       - Retain k-mers occurring up to this many times genome-wide
##                             (default: 1). Only meaningful with check_genome_uniqueness=true.
##   write_copy_histogram    - Write a TSV histogram of genome-wide k-mer copy numbers
##                             (default: false). Requires check_genome_uniqueness=true.
##   write_bed               - Write a BED of merged unique k-mer intervals (default: false)
##   write_bed_by_copies     - Write per-copy-tier BED files, one per copy count
##                             (default: false). Requires check_genome_uniqueness=true.
##   write_gene_stats        - Write a per-gene uniqueness TSV (default: false)
##   edit_distance_filter    - Discard candidates with any reference k-mer within Hamming
##                             distance N (default: 0 = disabled). Use 1 to keep only k-mers
##                             where no single-base sequencing error can produce a different
##                             reference k-mer. Triggers a full genome scan.
##   threads                 - CPU threads for the parallel genome scan (default: 8). Set to 0
##                             to use all logical CPUs on the machine. Has no effect when
##                             neither --check-genome-uniqueness nor --edit-distance-filter is set.
##   gene_padding            - Symmetrically extend each gene body by N bp before k-mer extraction
##                             (default: 0). Captures fusion breakpoints in flanking introns/UTRs
##                             just outside the annotated transcript; the genes table keeps true bounds.
##   docker_image            - geac Docker image, e.g. ghcr.io/fleharty/geac:latest
##   memory_gb               - Memory in GB (default: 64)
##   disk_gb                 - Disk space in GB (default: 100)
##   preemptible             - Preemptible retries (default: 0; long job — avoid mid-run eviction)
##
## Outputs:
##   fusion_index            - DuckDB k-mer index for use with geac_fusions.wdl
##   copy_histogram          - Copy-number histogram TSV (empty array when write_copy_histogram=false)
##   bed                     - Merged unique k-mer BED (empty array when write_bed=false)
##   beds_by_copies          - Per-copy-tier BED files, one per copy count (empty array when write_bed_by_copies=false)
##   gene_stats              - Per-gene uniqueness TSV (empty array when write_gene_stats=false)
##   geac_version            - geac binary version string (e.g. "0.4.33") for provenance in the data table

workflow GeacBuildFusionIndex {

    input {
        File   gene_annotation
        File   fasta
        File   fasta_index
        String index_name = "fusion_index"

        Int    kmer_size       = 23
        Int    min_gene_kmers  = 1
        File?  genes

        Boolean check_genome_uniqueness = false
        Int     max_genome_copies       = 1

        Boolean write_copy_histogram = false
        Boolean write_bed            = false
        Boolean write_bed_by_copies  = false
        Boolean write_gene_stats     = false
        Int     edit_distance_filter  = 0
        Int     threads               = 8
        Int     gene_padding          = 0

        String docker_image
        Int    memory_gb   = 64
        Int    disk_gb     = 100
        Int    preemptible = 0
    }

    call BuildFusionIndex {
        input:
            gene_annotation         = gene_annotation,
            fasta                   = fasta,
            fasta_index             = fasta_index,
            index_name              = index_name,
            kmer_size               = kmer_size,
            min_gene_kmers          = min_gene_kmers,
            genes                   = genes,
            check_genome_uniqueness = check_genome_uniqueness,
            max_genome_copies       = max_genome_copies,
            write_copy_histogram    = write_copy_histogram,
            write_bed               = write_bed,
            write_bed_by_copies     = write_bed_by_copies,
            write_gene_stats        = write_gene_stats,
            edit_distance_filter    = edit_distance_filter,
            threads                 = threads,
            gene_padding            = gene_padding,
            docker_image            = docker_image,
            memory_gb               = memory_gb,
            disk_gb                 = disk_gb,
            preemptible             = preemptible,
    }

    output {
        File        fusion_index   = BuildFusionIndex.fusion_index
        Array[File] copy_histogram = BuildFusionIndex.copy_histogram
        Array[File] bed            = BuildFusionIndex.bed
        Array[File] beds_by_copies = BuildFusionIndex.beds_by_copies
        Array[File] gene_stats     = BuildFusionIndex.gene_stats
        String      geac_version   = BuildFusionIndex.geac_version
    }
}

task BuildFusionIndex {

    input {
        File   gene_annotation
        File   fasta
        File   fasta_index
        String index_name

        Int    kmer_size
        Int    min_gene_kmers
        File?  genes

        Boolean check_genome_uniqueness
        Int     max_genome_copies

        Boolean write_copy_histogram
        Boolean write_bed
        Boolean write_bed_by_copies
        Boolean write_gene_stats
        Int     edit_distance_filter
        Int     threads
        Int     gene_padding

        String docker_image
        Int    memory_gb
        Int    disk_gb
        Int    preemptible
    }

    String output_db             = index_name + ".duckdb"
    String copy_histogram_tsv    = index_name + ".copy_histogram.tsv"
    String bed_file              = index_name + ".unique.bed"
    String bed_by_copies_prefix  = index_name
    String gene_stats_tsv        = index_name + ".gene_stats.tsv"

    command <<<
        set -euo pipefail

        # Record the geac version (e.g. "geac 0.4.33" -> "0.4.33") for the Terra
        # data table, so each index is traceable to the binary that built it.
        geac --version | awk '{print $2}' > geac_version.txt

        geac experimental build-fusion-index \
            --gene-annotation ~{gene_annotation} \
            --fasta           ~{fasta} \
            --output          ~{output_db} \
            --kmer-size       ~{kmer_size} \
            --min-gene-kmers  ~{min_gene_kmers} \
            ~{"--genes " + genes} \
            ~{if check_genome_uniqueness then "--check-genome-uniqueness" else ""} \
            ~{if check_genome_uniqueness then "--max-genome-copies " + max_genome_copies else ""} \
            ~{if write_copy_histogram then "--copy-histogram-output " + copy_histogram_tsv else ""} \
            ~{if write_bed            then "--bed-output "             + bed_file           else ""} \
            ~{if write_bed_by_copies  then "--bed-output-by-copies "   + bed_by_copies_prefix else ""} \
            ~{if write_gene_stats     then "--gene-stats-output "      + gene_stats_tsv     else ""} \
            ~{if edit_distance_filter > 0 then "--edit-distance-filter " + edit_distance_filter else ""} \
            ~{if gene_padding > 0 then "--gene-padding " + gene_padding else ""} \
            ~{if threads > 0 then "--threads " + threads else ""}
    >>>

    output {
        File        fusion_index   = output_db
        Array[File] copy_histogram = if write_copy_histogram then [copy_histogram_tsv] else []
        Array[File] bed            = if write_bed            then [bed_file]            else []
        Array[File] beds_by_copies = glob("*.copies*.bed")
        Array[File] gene_stats     = if write_gene_stats     then [gene_stats_tsv]      else []
        String      geac_version   = read_string("geac_version.txt")
    }

    runtime {
        docker:      docker_image
        memory:      memory_gb + " GB"
        cpu:         threads
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

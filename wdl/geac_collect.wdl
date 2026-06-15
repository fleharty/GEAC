version 1.0

## geac_collect.wdl
##
## Run `geac collect` on a single BAM/CRAM sample and produce a Parquet file
## of alt base metrics. Designed for use in Terra; scatter over a sample set
## to process a full cohort in parallel, then pass all Parquets to geac_merge.
##
## Inputs:
##   input_bam            - BAM or CRAM file
##   input_bam_index      - .bai / .crai index; passed to geac via --index, so it may live under a different name or directory than input_bam
##   bam_uri              - (optional) canonical BAM/CRAM URI stored in output metadata for IGV
##   bai_uri              - (optional) canonical index URI stored in output metadata for IGV
##   reference_fasta      - Reference FASTA
##   reference_fasta_index - Corresponding .fai index
##   read_type            - duplex | simplex | raw
##   pipeline             - free-text label; fgbio/dragen also select built-in family-size tag schemes
##   family_size_tags     - (optional) override family-size aux tags, e.g. "ab=aD,ba=bD,total=cD,fallback=sum"; overrides the pipeline preset
##   sample_id            - (optional) override sample ID; defaults to BAM SM tag
##   subject_id           - (optional) biological subject identifier (e.g. patient or animal)
##   sample_type          - (optional) sample substrate type (e.g. cfDNA, tumor_tissue)
##   batch                - (optional) batch/group label stored as a column in the output
##   label1               - (optional) free-text sample label 1 (e.g. tissue type)
##   label2               - (optional) free-text sample label 2 (e.g. library prep method)
##   label3               - (optional) free-text sample label 3 (e.g. sequencer type)
##   timepoint            - (optional) timepoint label for longitudinal studies (e.g. "T0", "week12")
##   vcf                  - (optional) VCF/BCF for variant call annotation
##   vcf_index            - (optional) Corresponding .tbi / .csi index
##   variants_tsv         - (optional) TSV variant list (chrom/pos_start/pos_end/ref/var, 0-based)
##                          Alternative to vcf; mutually exclusive.
##   gnomad               - (optional) bgzip+tabix-indexed gnomAD VCF/BCF for AF annotation
##   gnomad_index         - (optional) Corresponding .tbi / .csi index
##   gnomad_uri           - (optional) canonical gnomAD URI stored in output metadata for IGV
##   gnomad_index_uri     - (optional) canonical gnomAD index URI stored as gnomad_index_path for IGV; defaults to gnomad_index
##   gnomad_af_field      - INFO field to use as allele frequency (default "AF")
##   targets              - (optional) BED or Picard interval list; annotates on_target column
##   targets_uri          - (optional) canonical target-interval URI stored in output metadata for IGV
##   gene_annotations     - (optional) GFF3, GTF, or UCSC genePred; annotates gene column
##   region               - (optional) restrict to a region string, e.g. chr1:1-1000000
##   region_bed           - (optional) BED or Picard interval list; restricts pileup to those intervals
##   repeat_window        - bases each side of locus to scan for homopolymers/STRs (default 10)
##   min_base_qual        - minimum base quality (default 1)
##   min_map_qual         - minimum mapping quality (default 0)
##   max_pileup_depth     - max reads per pileup column; 0 = unlimited (default 0). htslib defaults to 8000 which silently downsamples high-coverage loci.
##   include_duplicates   - include PCR/optical duplicate reads (FLAG 0x400); default false
##   include_secondary    - include secondary alignments (FLAG 0x100); default false
##   include_supplementary - include supplementary alignments (FLAG 0x800); default false
##   exclude_tag          - drop reads whose aux tag equals a value, as "TAG:VALUE" (e.g. "RX:bad"); repeatable, default []
##   reads_output         - also write per-read detail Parquet (default false)
##   input_checksum_sha256 - compute SHA-256 for the input BAM/CRAM and store it in output Parquet provenance columns (default false)
##   threads              - VM CPU count (default 1); controls Terra VM sizing only — geac does not accept a --threads flag and uses its own internal parallelism
##   docker_image         - geac Docker image, e.g. ghcr.io/fleharty/geac:0.3.7
##   memory_gb            - memory in GB (default 32)
##   disk_gb              - disk space in GB (default 100)
##   preemptible          - number of preemptible retries (default 2)
##
## Outputs:
##   locus_parquet        - per-locus alt base Parquet ({stem}.locus.parquet or {stem}.parquet)
##   reads_parquets       - per-read detail Parquet array (one element when reads_output=true,
##                          empty array otherwise)
##   sample_metrics_parquets - sample-level target-depth metrics Parquet array (one element when
##                             targets are provided, empty array otherwise)

workflow GeacCollect {

    input {
        File   input_bam
        File   input_bam_index
        File   reference_fasta
        File   reference_fasta_index
        String read_type
        String pipeline

        String? family_size_tags
        String? sample_id
        String? subject_id
        String? sample_type
        String? batch
        String? label1
        String? label2
        String? label3
        String? timepoint
        File?   vcf
        File?   vcf_index
        File?   variants_tsv
        File?   gnomad
        File?   gnomad_index
        String? bam_uri
        String? bai_uri
        String? variants_uri
        String? gnomad_uri
        String? gnomad_index_uri
        String? targets_uri
        String  gnomad_af_field = "AF"
        File?   targets
        File?   gene_annotations
        String? region          # single region string, e.g. "chr1:1-1000000"
        File?   region_bed      # BED or Picard interval list; restricts processing to those intervals
        Int repeat_window = 10

        Int     min_base_qual  = 1
        Int     min_map_qual   = 0
        Int     max_pileup_depth = 0
        Boolean include_duplicates    = false
        Boolean include_secondary     = false
        Boolean include_supplementary = false
        Array[String] exclude_tag     = []
        Boolean reads_output   = false
        Boolean input_checksum_sha256 = false
        Int     threads        = 1

        String docker_image
        Int    memory_gb    = 32
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
            family_size_tags      = family_size_tags,
            sample_id             = sample_id,
            subject_id            = subject_id,
            sample_type           = sample_type,
            batch                 = batch,
            label1                = label1,
            label2                = label2,
            label3                = label3,
            timepoint             = timepoint,
            vcf                   = vcf,
            vcf_index             = vcf_index,
            variants_tsv          = variants_tsv,
            gnomad                = gnomad,
            gnomad_index          = gnomad_index,
            gnomad_af_field       = gnomad_af_field,
            targets               = targets,
            gene_annotations      = gene_annotations,
            region                = region,
            region_bed            = region_bed,
            repeat_window         = repeat_window,
            min_base_qual         = min_base_qual,
            min_map_qual          = min_map_qual,
            max_pileup_depth      = max_pileup_depth,
            include_duplicates    = include_duplicates,
            include_secondary     = include_secondary,
            include_supplementary = include_supplementary,
            exclude_tag           = exclude_tag,
            reads_output          = reads_output,
            input_checksum_sha256 = input_checksum_sha256,
            threads               = threads,
            bam_uri               = if defined(bam_uri) then select_first([bam_uri]) else input_bam,
            bai_uri               = if defined(bai_uri) then select_first([bai_uri]) else input_bam_index,
            variants_uri          = if defined(variants_uri) then variants_uri else if defined(vcf) then vcf else variants_tsv,
            gnomad_uri            = if defined(gnomad_uri) then gnomad_uri else gnomad,
            gnomad_index_uri      = if defined(gnomad_index_uri) then gnomad_index_uri else gnomad_index,
            targets_uri           = if defined(targets_uri) then targets_uri else targets,
            docker_image          = docker_image,
            memory_gb             = memory_gb,
            disk_gb               = disk_gb,
            preemptible           = preemptible,
    }

    output {
        File        locus_parquet  = Collect.locus_parquet
        Array[File] reads_parquets = Collect.reads_parquets
        Array[File] sample_metrics_parquets = Collect.sample_metrics_parquets
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

        String? family_size_tags
        String? sample_id
        String? subject_id
        String? sample_type
        String? batch
        String? label1
        String? label2
        String? label3
        String? timepoint
        File?   vcf
        File?   vcf_index
        File?   variants_tsv
        File?   gnomad
        File?   gnomad_index
        String  gnomad_af_field
        File?   targets
        File?   gene_annotations
        String? region
        File?   region_bed
        Int     repeat_window

        Int     min_base_qual
        Int     min_map_qual
        Int     max_pileup_depth
        Boolean include_duplicates
        Boolean include_secondary
        Boolean include_supplementary
        Array[String] exclude_tag
        Boolean reads_output
        Boolean input_checksum_sha256
        Int     threads

        String? bam_uri
        String? bai_uri
        String? variants_uri
        String? gnomad_uri
        String? gnomad_index_uri
        String? targets_uri

        String docker_image
        Int    memory_gb
        Int    disk_gb
        Int    preemptible
    }

    # When reads_output=true, geac derive two files: {stem}.locus.parquet and {stem}.reads.parquet.
    # When reads_output=false, geac writes a single {stem}.parquet.
    String stem        = sub(basename(input_bam), "\\.(bam|cram)$", "")
    String output_arg  = stem + ".parquet"
    String locus_name  = if reads_output then stem + ".locus.parquet" else stem + ".parquet"

    command <<<
        set -euo pipefail

        geac collect \
            --input            ~{input_bam} \
            --index            ~{input_bam_index} \
            --reference        ~{reference_fasta} \
            --output           ~{output_arg} \
            --read-type        ~{read_type} \
            --pipeline         ~{pipeline} \
            ~{"--family-size-tags " + family_size_tags} \
            --min-base-qual    ~{min_base_qual} \
            --min-map-qual     ~{min_map_qual} \
            --max-pileup-depth ~{max_pileup_depth} \
            ~{"--sample-id "        + sample_id} \
            ~{"--subject-id "       + subject_id} \
            ~{"--sample-type "      + sample_type} \
            ~{"--batch "            + batch} \
            ~{"--label1 "           + label1} \
            ~{"--label2 "           + label2} \
            ~{"--label3 "           + label3} \
            ~{"--timepoint "        + timepoint} \
            ~{"--vcf "              + vcf} \
            ~{"--variants-tsv "     + variants_tsv} \
            ~{"--gnomad "           + gnomad} \
            ~{if defined(gnomad) then "--gnomad-af-field " + gnomad_af_field else ""} \
            ~{"--targets "          + targets} \
            ~{"--gene-annotations " + gene_annotations} \
            ~{if defined(region) then "--region " + select_first([region]) else if defined(region_bed) then "--region " + select_first([region_bed]) else ""} \
            --repeat-window ~{repeat_window} \
            ~{if include_duplicates    then "--include-duplicates"    else ""} \
            ~{if include_secondary     then "--include-secondary"     else ""} \
            ~{if include_supplementary then "--include-supplementary" else ""} \
            ~{sep=" " prefix("--exclude-tag ", exclude_tag)} \
            ~{if input_checksum_sha256 then "--input-checksum-sha256" else ""} \
            ~{if reads_output then "--reads-output" else ""} \
            ~{"--bam-uri "      + bam_uri} \
            ~{"--bai-uri "      + bai_uri} \
            ~{"--variants-uri " + variants_uri} \
            ~{"--gnomad-uri "       + gnomad_uri} \
            ~{"--gnomad-index-uri " + gnomad_index_uri} \
            ~{"--targets-uri "  + targets_uri}
    >>>

    output {
        File        locus_parquet  = locus_name
        Array[File] reads_parquets = glob("*.reads.parquet")
        Array[File] sample_metrics_parquets = glob("*.sample_metrics.parquet")
    }

    runtime {
        docker:      docker_image
        memory:      memory_gb + " GB"
        cpu:         threads
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

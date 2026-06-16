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
##   input_bams              - BAM or CRAM files (localized by Cromwell for geac collect)
##   input_bam_indices       - .bai / .crai indices; passed to geac via --index, so each may live under a different name/directory than its BAM
##   bam_uris                - (optional) canonical BAM/CRAM URIs stored in output metadata for IGV
##   bai_uris                - (optional) canonical index URIs stored in output metadata for IGV
##   sample_ids              - (optional) override sample IDs; defaults to SM tag per BAM
##   subject_ids             - (optional) per-sample biological subject identifier
##   sample_types            - (optional) per-sample sample substrate type (e.g. cfDNA, tumor_tissue)
##   variants_tsvs           - (optional) per-sample variant TSV files
##   vcfs                    - (optional) per-sample VCF/BCF files for variant annotation
##   vcf_indices             - (optional) per-sample .tbi / .csi indices; required when vcfs provided
##
## Inputs (shared across all samples):
##   reference_fasta         - Reference FASTA
##   reference_fasta_index   - Corresponding .fai index
##   read_types              - (optional) per-sample array of duplex|simplex|raw; defaults to "duplex" for all
##   pipelines               - (optional) per-sample free-text pipeline label; fgbio/dragen select built-in family-size schemes; defaults to "fgbio" for all
##   family_size_tags        - (optional) one family-size aux-tag spec applied to ALL samples, e.g. "total=cD"; overrides the pipeline preset
##   batches                 - (optional) per-sample batch/group label stored as a column in each Parquet
##   labels1                 - (optional) per-sample free-text label 1 (e.g. tissue type)
##   labels2                 - (optional) per-sample free-text label 2 (e.g. library prep method)
##   labels3                 - (optional) per-sample free-text label 3 (e.g. sequencer type)
##   timepoints              - (optional) per-sample timepoint label for longitudinal studies (e.g. "T0", "week12")
##   gnomad                  - (optional) bgzip+tabix-indexed gnomAD VCF/BCF for AF annotation
##   gnomad_index            - (optional) Corresponding .tbi / .csi index
##   gnomad_uri              - (optional) canonical gnomAD URI stored in output metadata for IGV
##   gnomad_index_uri        - (optional) canonical gnomAD index URI stored as gnomad_index_path for IGV; defaults to gnomad_index
##   gnomad_af_field         - INFO field to use as allele frequency (default "AF")
##   targets                 - (optional) BED or Picard interval list
##   targets_uri             - (optional) canonical target-interval URI stored in output metadata for IGV
##   gene_annotations        - (optional) GTF, GFF3, or UCSC genePred (.txt/.txt.gz)
##   region                  - (optional) restrict all samples to a genomic region
##   scatter_interval_list    - (optional) BED/interval list; when set, each sample's collect is scattered across scatter_count shards, gathered, then merged (results identical to unsharded). Leave unset for one-collect-per-sample.
##   scatter_count            - number of interval shards per sample when scatter_interval_list is set (default 50)
##   repeat_window           - bases each side of locus for homopolymer/STR scan (default 10)
##   max_pileup_depth        - max reads per pileup column; 0 = unlimited (default 0). htslib defaults to 8000 which silently downsamples high-coverage loci.
##   include_duplicates      - include PCR/optical duplicate reads (FLAG 0x400); default false
##   include_secondary       - include secondary alignments (FLAG 0x100); default false
##   include_supplementary   - include supplementary alignments (FLAG 0x800); default false
##   exclude_tag             - drop reads whose aux tag equals a value, as "TAG:VALUE" (e.g. "cD:1"); repeatable, applied to all samples, default []
##   reads_output            - also write per-read detail Parquets and merge into alt_reads table (default false)
##   input_checksum_sha256   - compute SHA-256 for each input BAM/CRAM during collect (default false)
##   run_fragments           - also run geac fragments in parallel and ingest the fragment Parquets
##                             into the cohort DuckDB as a table (default false).
##                             Note: requires a second BAM pass per sample (parallel, not sequential).
##                             Fragment Parquets are large for WGS — increase fragments_disk_gb accordingly.
##   cohort_name             - Base name for the output DuckDB file (default: cohort)
##   docker_image            - geac Docker image, e.g. ghcr.io/fleharty/geac:latest
##
## Outputs:
##   locus_parquets          - Per-sample locus Parquet files from geac collect
##   reads_parquets          - Per-sample reads Parquet files (empty when reads_output=false)
##   sample_metrics_parquets - Per-sample sample_metrics Parquet files (empty when targets absent)
##   fragments_parquets      - Per-sample fragment Parquet files (empty when run_fragments=false)
##   cohort_db               - Merged cohort DuckDB from geac merge

workflow GeacCohort {

    input {
        # Per-sample parallel arrays
        Array[File]    input_bams
        Array[File]    input_bam_indices
        Array[String]? bam_uris           # optional; canonical BAM/CRAM URIs for embedded IGV paths
        Array[String]? bai_uris           # optional; canonical index URIs for embedded IGV paths
        Array[String]? sample_ids          # optional; if provided must be same length as input_bams
        Array[String]? subject_ids         # optional; if provided must be same length as input_bams
        Array[String]? sample_types        # optional; if provided must be same length as input_bams
        Array[File]?   variants_tsvs       # optional; if provided must be same length as input_bams
        Array[File]?   vcfs                # optional; if provided must be same length as input_bams
        Array[File]?   vcf_indices         # optional; required when vcfs is provided

        # Shared inputs
        File   reference_fasta
        File   reference_fasta_index
        Array[String]? read_types      # optional; if provided must be same length as input_bams
        Array[String]? pipelines       # optional; if provided must be same length as input_bams
        String? family_size_tags       # optional; one family-size tag spec applied to ALL samples (e.g. "total=cD")
        Array[String]? batches         # optional; per-sample batch/group label stored as a column
        Array[String]? labels1         # optional; per-sample free-text label 1
        Array[String]? labels2         # optional; per-sample free-text label 2
        Array[String]? labels3         # optional; per-sample free-text label 3
        Array[String]? timepoints      # optional; per-sample timepoint label for longitudinal studies

        File?   gnomad
        File?   gnomad_index
        String? gnomad_uri
        String? gnomad_index_uri
        String  gnomad_af_field = "AF"

        File?   targets
        String? targets_uri
        File?   gene_annotations
        String? region
        Int     repeat_window = 10

        # Interval-scatter: when scatter_interval_list is provided, each sample's
        # `geac collect` is split across scatter_count interval shards, run in
        # parallel, then gathered. Loci/reads Parquets concatenate; per-shard partial
        # sample-metrics are aggregated back into one exact row per sample. Leave
        # scatter_interval_list unset for the original one-Collect-per-sample behavior.
        File?   scatter_interval_list
        Int     scatter_count = 50

        Int     min_base_qual = 1
        Int     min_map_qual  = 0
        Int     max_pileup_depth = 0
        Boolean include_duplicates    = false
        Boolean include_secondary     = false
        Boolean include_supplementary = false
        Array[String] exclude_tag     = []
        Boolean reads_output  = false
        Boolean input_checksum_sha256 = true
        Boolean run_fragments = false
        Int     threads       = 1

        String cohort_name = "cohort"

        # Resource settings
        String docker_image
        Int    collect_memory_gb       = 8
        Int    collect_disk_gb         = 100
        Int    fragments_disk_gb       = 200
        Int    merge_memory_gb         = 16
        Int    merge_disk_gb           = 50
        Int    split_memory_gb         = 4
        Int    split_disk_gb           = 20
        Int    aggregate_memory_gb     = 4
        Int    aggregate_disk_gb       = 20
        Int    preemptible             = 2
    }

    # Split the interval list once; shards are shared across all samples.
    if (defined(scatter_interval_list)) {
        call SplitIntervals {
            input:
                interval_list = select_first([scatter_interval_list]),
                scatter_count = scatter_count,
                docker_image  = docker_image,
                memory_gb     = split_memory_gb,
                disk_gb       = split_disk_gb,
                preemptible   = preemptible,
        }
    }

    # Hoist the conditional task output into a concrete workflow-level Array[File]
    # (empty when not scattering). Referencing an optional call output from inside a
    # nested scatter trips Cromwell's sub-workflow input evaluation; a plain
    # Array[File] threads in cleanly.
    Array[File] shard_beds = select_first([SplitIntervals.shards, []])

    scatter (i in range(length(input_bams))) {

        # Safely index into optional per-sample arrays.
        # Variables declared inside `if` blocks are typed as optional (String?, File?)
        # in the enclosing scope, which is exactly what the Collect task expects.
        if (defined(sample_ids)) {
            String this_sample_id    = select_first([sample_ids])[i]
        }
        if (defined(subject_ids)) {
            String this_subject_id   = select_first([subject_ids])[i]
        }
        if (defined(sample_types)) {
            String this_sample_type  = select_first([sample_types])[i]
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
        if (defined(batches)) {
            String this_batch  = select_first([batches])[i]
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
        String this_bam_uri = if defined(bam_uris) then select_first([bam_uris])[i] else input_bams[i]
        String this_bai_uri = if defined(bai_uris) then select_first([bai_uris])[i] else input_bam_indices[i]

        Array[String] manifest_row = [
            this_bam_uri,
            this_bai_uri,
            select_first([this_sample_id, ""]),
            select_first([this_subject_id, ""]),
            select_first([this_sample_type, ""]),
            select_first([this_batch, ""]),
            this_read_type,
            this_pipeline,
            select_first([this_label1, ""]),
            select_first([this_label2, ""]),
            select_first([this_label3, ""]),
            select_first([this_timepoint, ""]),
            select_first([gnomad_uri, if defined(gnomad) then gnomad else ""]),
            select_first([targets_uri, if defined(targets) then targets else ""]),
        ]

        # ── Branch A: scatter this sample across interval shards ─────────────────
        # Each shard runs collect on its interval with --sample-metrics-partial; the
        # per-shard partial metrics are aggregated back into one exact row per sample.
        if (defined(scatter_interval_list)) {
            scatter (shard in shard_beds) {
                call Collect as CollectShard {
                    input:
                        input_bam             = input_bams[i],
                        input_bam_index       = input_bam_indices[i],
                        reference_fasta       = reference_fasta,
                        reference_fasta_index = reference_fasta_index,
                        read_type             = this_read_type,
                        pipeline              = this_pipeline,
                        family_size_tags      = family_size_tags,
                        batch                 = this_batch,
                        label1                = this_label1,
                        label2                = this_label2,
                        label3                = this_label3,
                        timepoint             = this_timepoint,
                        sample_id             = this_sample_id,
                        subject_id            = this_subject_id,
                        sample_type           = this_sample_type,
                        variants_tsv          = this_variants_tsv,
                        vcf                   = this_vcf,
                        vcf_index             = this_vcf_index,
                        gnomad                = gnomad,
                        gnomad_index          = gnomad_index,
                        gnomad_uri            = if defined(gnomad_uri) then select_first([gnomad_uri]) else gnomad,
                        gnomad_index_uri      = if defined(gnomad_index_uri) then select_first([gnomad_index_uri]) else gnomad_index,
                        gnomad_af_field       = gnomad_af_field,
                        targets               = targets,
                        targets_uri           = if defined(targets_uri) then select_first([targets_uri]) else targets,
                        gene_annotations      = gene_annotations,
                        region_bed            = shard,
                        sample_metrics_partial = defined(targets),
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
                        bam_uri               = this_bam_uri,
                        bai_uri               = this_bai_uri,
                        docker_image          = docker_image,
                        memory_gb             = collect_memory_gb,
                        disk_gb               = collect_disk_gb,
                        preemptible           = preemptible,
                }
            }

            if (defined(targets)) {
                call AggregateMetrics {
                    input:
                        partial_parquets = flatten(CollectShard.sample_metrics_partial_parquets),
                        output_basename  = sub(basename(input_bams[i]), "\\.(bam|cram)$", ""),
                        docker_image     = docker_image,
                        memory_gb        = aggregate_memory_gb,
                        disk_gb          = aggregate_disk_gb,
                        preemptible      = preemptible,
                }
            }

            Array[File] a_loci    = CollectShard.locus_parquet
            Array[File] a_reads   = flatten(CollectShard.reads_parquets)
            Array[File] a_metrics = if defined(AggregateMetrics.sample_metrics)
                                    then [select_first([AggregateMetrics.sample_metrics])]
                                    else []
        }

        # ── Branch B: original single-pass collect (no scatter list provided) ────
        if (!defined(scatter_interval_list)) {
            call Collect as CollectWhole {
                input:
                    input_bam             = input_bams[i],
                    input_bam_index       = input_bam_indices[i],
                    reference_fasta       = reference_fasta,
                    reference_fasta_index = reference_fasta_index,
                    read_type             = this_read_type,
                    pipeline              = this_pipeline,
                    family_size_tags      = family_size_tags,
                    batch                 = this_batch,
                    label1                = this_label1,
                    label2                = this_label2,
                    label3                = this_label3,
                    timepoint             = this_timepoint,
                    sample_id             = this_sample_id,
                    subject_id            = this_subject_id,
                    sample_type           = this_sample_type,
                    variants_tsv          = this_variants_tsv,
                    vcf                   = this_vcf,
                    vcf_index             = this_vcf_index,
                    gnomad                = gnomad,
                    gnomad_index          = gnomad_index,
                    gnomad_uri            = if defined(gnomad_uri) then select_first([gnomad_uri]) else gnomad,
                    gnomad_index_uri      = if defined(gnomad_index_uri) then select_first([gnomad_index_uri]) else gnomad_index,
                    gnomad_af_field       = gnomad_af_field,
                    targets               = targets,
                    targets_uri           = if defined(targets_uri) then select_first([targets_uri]) else targets,
                    gene_annotations      = gene_annotations,
                    region                = region,
                    sample_metrics_partial = false,
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
                    bam_uri               = this_bam_uri,
                    bai_uri               = this_bai_uri,
                    docker_image          = docker_image,
                    memory_gb             = collect_memory_gb,
                    disk_gb               = collect_disk_gb,
                    preemptible           = preemptible,
            }
            Array[File] b_loci    = [CollectWhole.locus_parquet]
            Array[File] b_reads   = CollectWhole.reads_parquets
            Array[File] b_metrics = CollectWhole.sample_metrics_parquets
        }

        # Unify the two mutually-exclusive branches into per-sample arrays. Exactly
        # one branch is defined; select_all drops the null branch and flatten yields
        # its contents (more robust than select_first, which throws if both are null).
        Array[File] sample_loci    = flatten(select_all([a_loci, b_loci]))
        Array[File] sample_reads   = flatten(select_all([a_reads, b_reads]))
        Array[File] sample_metrics = flatten(select_all([a_metrics, b_metrics]))

        if (run_fragments) {
            call Fragments {
                input:
                    input_bam             = input_bams[i],
                    input_bam_index       = input_bam_indices[i],
                    reference_fasta       = reference_fasta,
                    reference_fasta_index = reference_fasta_index,
                    read_type             = this_read_type,
                    pipeline              = this_pipeline,
                    sample_id             = this_sample_id,
                    subject_id            = this_subject_id,
                    sample_type           = this_sample_type,
                    batch                 = this_batch,
                    label1                = this_label1,
                    label2                = this_label2,
                    label3                = this_label3,
                    timepoint             = this_timepoint,
                    region                = region,
                    min_map_qual          = min_map_qual,
                    docker_image          = docker_image,
                    memory_gb             = collect_memory_gb,
                    disk_gb               = fragments_disk_gb,
                    preemptible           = preemptible,
            }
        }
    }

    # Gather per-sample outputs across the scatter. Each sample contributes an
    # Array[File] (one locus Parquet when unsharded, many when sharded); flatten to
    # one flat array. reads/metrics are empty arrays when not produced.
    Array[File] all_locus_parquets = flatten(sample_loci)
    Array[File] all_reads_parquets = flatten(sample_reads)
    Array[File] all_sample_metrics_parquets = flatten(sample_metrics)
    Array[File] all_fragments_parquets = select_all(Fragments.fragments_parquet)

    call Merge {
        input:
            parquets       = flatten([all_locus_parquets, all_reads_parquets, all_sample_metrics_parquets, all_fragments_parquets]),
            manifest_rows  = manifest_row,
            cohort_name    = cohort_name,
            docker_image   = docker_image,
            memory_gb      = merge_memory_gb,
            disk_gb        = merge_disk_gb,
            preemptible    = preemptible,
    }

    output {
        File cohort_db              = Merge.cohort_db
        File cohort_manifest        = Merge.cohort_manifest
        File cohort_on_target_tsv   = Merge.cohort_on_target_tsv
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

        String? family_size_tags
        String? sample_id
        String? subject_id
        String? sample_type
        String? batch
        String? label1
        String? label2
        String? label3
        String? timepoint
        File?   variants_tsv
        File?   vcf
        File?   vcf_index
        File?   gnomad
        File?   gnomad_index
        String  gnomad_af_field
        File?   targets
        File?   gene_annotations
        String? region
        # A BED/interval-list passed as a File so Cromwell localizes it into the
        # task (a File bound to a String input is NOT localized). Used for shard
        # scatter; takes precedence over the string `region` when set.
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
        Boolean sample_metrics_partial = false
        Boolean input_checksum_sha256
        Int     threads

        String? bam_uri
        String? bai_uri
        String? gnomad_uri
        String? gnomad_index_uri
        String? targets_uri

        String docker_image
        Int    memory_gb
        Int    disk_gb
        Int    preemptible
    }

    String stem       = sub(basename(input_bam), "\\.(bam|cram)$", "")
    String output_arg = stem + ".parquet"
    String locus_name = if reads_output then stem + ".locus.parquet" else stem + ".parquet"

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
            ~{"--gnomad-uri "       + gnomad_uri} \
            ~{"--gnomad-index-uri " + gnomad_index_uri} \
            ~{if defined(gnomad) then "--gnomad-af-field " + gnomad_af_field else ""} \
            ~{"--targets "          + targets} \
            ~{"--targets-uri "      + targets_uri} \
            ~{"--gene-annotations " + gene_annotations} \
            ~{if defined(region_bed) then "--region " + region_bed else if defined(region) then "--region " + region else ""} \
            --repeat-window ~{repeat_window} \
            ~{if include_duplicates    then "--include-duplicates"    else ""} \
            ~{if include_secondary     then "--include-secondary"     else ""} \
            ~{if include_supplementary then "--include-supplementary" else ""} \
            ~{sep=" " prefix("--exclude-tag ", exclude_tag)} \
            ~{if input_checksum_sha256 then "--input-checksum-sha256" else ""} \
            ~{if reads_output then "--reads-output" else ""} \
            ~{if sample_metrics_partial then "--sample-metrics-partial" else ""} \
            ~{"--bam-uri " + bam_uri} \
            ~{"--bai-uri " + bai_uri}
    >>>

    output {
        File        locus_parquet  = locus_name
        Array[File] reads_parquets = glob("*.reads.parquet")
        # Final sample_metrics (empty when --sample-metrics-partial); the partial
        # glob deliberately does not match the final ".sample_metrics.parquet".
        Array[File] sample_metrics_parquets = glob("*.sample_metrics.parquet")
        Array[File] sample_metrics_partial_parquets = glob("*.sample_metrics_partial.parquet")
    }

    runtime {
        docker:      docker_image
        memory:      memory_gb + " GB"
        cpu:         threads
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

task SplitIntervals {

    input {
        File   interval_list
        Int    scatter_count

        String docker_image
        Int    memory_gb
        Int    disk_gb
        Int    preemptible
    }

    command <<<
        set -euo pipefail
        mkdir -p shards
        geac split-intervals \
            --input ~{interval_list} \
            --scatter-count ~{scatter_count} \
            --output-dir shards
    >>>

    output {
        Array[File] shards = glob("shards/*.bed")
    }

    runtime {
        docker:      docker_image
        memory:      memory_gb + " GB"
        cpu:         1
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

task AggregateMetrics {

    input {
        Array[File] partial_parquets
        String      output_basename

        String docker_image
        Int    memory_gb
        Int    disk_gb
        Int    preemptible
    }

    command <<<
        set -euo pipefail
        geac aggregate-metrics \
            --output ~{output_basename}.sample_metrics.parquet \
            ~{sep=" " partial_parquets}
    >>>

    output {
        File sample_metrics = output_basename + ".sample_metrics.parquet"
    }

    runtime {
        docker:      docker_image
        memory:      memory_gb + " GB"
        cpu:         1
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

task Merge {

    input {
        Array[File]          parquets
        Array[Array[String]] manifest_rows
        String               cohort_name

        String docker_image
        Int    memory_gb
        Int    disk_gb
        Int    preemptible
    }

    # write_tsv is called inside the task so it writes to the worker VM's local disk,
    # avoiding the GCS-only filesystem restriction that applies at workflow scope.
    File manifest_data = write_tsv(manifest_rows)

    String output_db        = cohort_name + ".duckdb"
    String output_manifest  = cohort_name + ".manifest.tsv"
    String output_on_target = cohort_name + ".on_target.tsv"

    command <<<
        set -euo pipefail

        printf 'bam_path\tbai_path\tsample_id\tsubject_id\tsample_type\tbatch\tread_type\tpipeline\tlabel1\tlabel2\tlabel3\ttimepoint\tgnomad_path\ttargets_path\n' > ~{output_manifest}
        cat ~{manifest_data} >> ~{output_manifest}

        touch ~{output_on_target}

        geac merge \
            --output ~{output_db} \
            --on-target-tsv ~{output_on_target} \
            ~{sep=" " parquets}
    >>>

    output {
        File cohort_db            = output_db
        File cohort_manifest      = output_manifest
        File cohort_on_target_tsv = output_on_target
    }

    runtime {
        docker:      docker_image
        memory:      memory_gb + " GB"
        cpu:         2
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

task Fragments {

    input {
        File   input_bam
        File   input_bam_index
        File   reference_fasta
        File   reference_fasta_index
        String read_type
        String pipeline

        String? sample_id
        String? subject_id
        String? sample_type
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

    String stem        = sub(basename(input_bam), "\\.(bam|cram)$", "")
    String output_name = stem + ".fragments.parquet"

    command <<<
        set -euo pipefail

        geac fragments \
            --input            ~{input_bam} \
            --index            ~{input_bam_index} \
            --reference        ~{reference_fasta} \
            --output           ~{output_name} \
            --read-type        ~{read_type} \
            --pipeline         ~{pipeline} \
            --min-map-qual     ~{min_map_qual} \
            ~{"--sample-id "   + sample_id} \
            ~{"--subject-id "  + subject_id} \
            ~{"--sample-type " + sample_type} \
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

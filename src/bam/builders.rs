use crate::cli::CollectArgs;
use crate::record::{AltBase, AltRead, VariantType};
use crate::repeat::RepeatMetrics;

use super::indel::IndelCount;
use super::pileup::{true_cycle, BaseTally, ReadDetail};

pub(super) struct LocusContext {
    sample_id: String,
    chrom: String,
    pos: i64,
    ref_allele: String,
    total_depth: i32,
    ref_count: i32,
    fwd_depth: i32,
    rev_depth: i32,
    fwd_ref_count: i32,
    rev_ref_count: i32,
    overlap_depth: i32,
    overlap_ref_agree: i32,
    read_type: crate::record::ReadType,
    pipeline: crate::record::Pipeline,
    batch: Option<String>,
    label1: Option<String>,
    label2: Option<String>,
    label3: Option<String>,
    input_checksum_sha256: Option<String>,
    on_target: Option<bool>,
    gene: Option<String>,
    homopolymer_len: i32,
    str_period: i32,
    str_len: i32,
    trinuc_context: Option<String>,
}

impl LocusContext {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn new(
        args: &CollectArgs,
        sample_id: &str,
        chrom: &str,
        pos: i64,
        ref_base: char,
        total_depth: i32,
        fwd_depth: i32,
        rev_depth: i32,
        overlap_depth: i32,
        ref_tally: Option<&BaseTally>,
        on_target: Option<bool>,
        gene: Option<String>,
        repeat: &RepeatMetrics,
        trinuc_context: Option<String>,
        input_checksum_sha256: Option<String>,
    ) -> Self {
        Self {
            sample_id: sample_id.to_string(),
            chrom: chrom.to_string(),
            pos,
            ref_allele: ref_base.to_string(),
            total_depth,
            ref_count: ref_tally.map_or(0, |t| t.total),
            fwd_depth,
            rev_depth,
            fwd_ref_count: ref_tally.map_or(0, |t| t.fwd),
            rev_ref_count: ref_tally.map_or(0, |t| t.rev),
            overlap_depth,
            overlap_ref_agree: ref_tally.map_or(0, |t| t.overlap_alt_agree),
            read_type: args.read_type,
            pipeline: args.pipeline,
            batch: args.batch.clone(),
            label1: args.label1.clone(),
            label2: args.label2.clone(),
            label3: args.label3.clone(),
            input_checksum_sha256,
            on_target,
            gene,
            homopolymer_len: repeat.homopolymer_len,
            str_period: repeat.str_period,
            str_len: repeat.str_len,
            trinuc_context,
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn build_alt_base(
        &self,
        alt_allele: String,
        variant_type: VariantType,
        alt_count: i32,
        fwd_alt_count: i32,
        rev_alt_count: i32,
        overlap_alt_agree: i32,
        overlap_alt_disagree: i32,
        variant_called: Option<bool>,
        variant_filter: Option<String>,
        gnomad_af: Option<f32>,
        trinuc_context: Option<String>,
    ) -> AltBase {
        AltBase {
            sample_id: self.sample_id.clone(),
            chrom: self.chrom.clone(),
            pos: self.pos,
            ref_allele: self.ref_allele.clone(),
            alt_allele,
            variant_type,
            total_depth: self.total_depth,
            alt_count,
            ref_count: self.ref_count,
            fwd_depth: self.fwd_depth,
            rev_depth: self.rev_depth,
            fwd_alt_count,
            rev_alt_count,
            fwd_ref_count: self.fwd_ref_count,
            rev_ref_count: self.rev_ref_count,
            overlap_depth: self.overlap_depth,
            overlap_alt_agree,
            overlap_alt_disagree,
            overlap_ref_agree: self.overlap_ref_agree,
            read_type: self.read_type,
            pipeline: self.pipeline,
            batch: self.batch.clone(),
            label1: self.label1.clone(),
            label2: self.label2.clone(),
            label3: self.label3.clone(),
            input_checksum_sha256: self.input_checksum_sha256.clone(),
            variant_called,
            variant_filter,
            on_target: self.on_target,
            gene: self.gene.clone(),
            homopolymer_len: self.homopolymer_len,
            str_period: self.str_period,
            str_len: self.str_len,
            trinuc_context,
            gnomad_af,
        }
    }

    pub(super) fn build_alt_read(&self, alt_allele: &str, detail: &ReadDetail) -> AltRead {
        AltRead {
            sample_id: self.sample_id.clone(),
            chrom: self.chrom.clone(),
            pos: self.pos,
            alt_allele: alt_allele.to_string(),
            cycle: true_cycle(
                detail.qpos,
                detail.read_len,
                detail.is_reverse,
                detail.hard_clip_before,
            ),
            read_length: detail.read_len as i32,
            is_read1: detail.is_first_in_pair,
            ab_count: detail.ab_count,
            ba_count: detail.ba_count,
            family_size: detail.family_size,
            base_qual: detail.base_qual as i32,
            map_qual: detail.map_qual as i32,
            insert_size: detail.insert_size,
            input_checksum_sha256: self.input_checksum_sha256.clone(),
        }
    }

    pub(super) fn push_snv_record(
        &self,
        records: &mut Vec<AltBase>,
        base: char,
        tally: &BaseTally,
        variant_called: Option<bool>,
        variant_filter: Option<String>,
        gnomad_af: Option<f32>,
    ) {
        records.push(self.build_alt_base(
            base.to_string(),
            VariantType::Snv,
            tally.total,
            tally.fwd,
            tally.rev,
            tally.overlap_alt_agree,
            tally.overlap_alt_disagree,
            variant_called,
            variant_filter,
            gnomad_af,
            self.trinuc_context.clone(),
        ));
    }

    pub(super) fn push_indel_record(
        &self,
        records: &mut Vec<AltBase>,
        indel: &IndelCount,
        variant_called: Option<bool>,
        variant_filter: Option<String>,
        gnomad_af: Option<f32>,
    ) {
        records.push(self.build_alt_base(
            indel.alt_allele.clone(),
            indel.variant_type,
            indel.total,
            indel.fwd,
            indel.rev,
            indel.overlap_alt_agree,
            indel.overlap_alt_disagree,
            variant_called,
            variant_filter,
            gnomad_af,
            None,
        ));
    }
}

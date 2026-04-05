import pandas as pd

from pipeline_compare_helpers import (
    build_unique_pipeline_characterization_df,
    summarize_unique_pipeline_groups,
)


def test_build_unique_pipeline_characterization_df_normalizes_rows():
    pc_uniq = pd.DataFrame(
        [
            {
                "concordance": "only_a",
                "sample_id": "S1",
                "chrom": "chr1",
                "pos": 10,
                "alt_allele": "T",
                "variant_type": "SNV",
                "trinuc_context": "ACA",
                "ref_allele": "C",
                "vaf_a": 0.05,
                "vaf_b": None,
                "depth_a": 100,
                "depth_b": None,
                "alt_count_a": 5,
                "alt_count_b": None,
            },
            {
                "concordance": "only_b",
                "sample_id": "S2",
                "chrom": "chr2",
                "pos": 20,
                "alt_allele": "A",
                "variant_type": "SNV",
                "trinuc_context": "TGT",
                "ref_allele": "G",
                "vaf_a": None,
                "vaf_b": 0.12,
                "depth_a": None,
                "depth_b": 80,
                "alt_count_a": None,
                "alt_count_b": 10,
            },
            {
                "concordance": "shared",
                "sample_id": "S3",
                "chrom": "chr3",
                "pos": 30,
                "alt_allele": "C",
                "variant_type": "SNV",
                "trinuc_context": "ATA",
                "ref_allele": "T",
                "vaf_a": 0.2,
                "vaf_b": 0.19,
                "depth_a": 90,
                "depth_b": 88,
                "alt_count_a": 18,
                "alt_count_b": 17,
            },
        ]
    )

    out = build_unique_pipeline_characterization_df(pc_uniq, "PipeA", "PipeB")

    assert list(out["group"]) == ["Only PipeA", "Only PipeB"]
    assert list(out["pipeline"]) == ["PipeA", "PipeB"]
    assert list(out["vaf"]) == [0.05, 0.12]
    assert list(out["depth"]) == [100, 80]
    assert list(out["alt_count"]) == [5, 10]
    assert "shared" not in out.get("concordance", pd.Series(dtype=str)).tolist()


def test_summarize_unique_pipeline_groups_reports_pipeline_aware_metrics():
    uniq_cmp = pd.DataFrame(
        [
            {"group": "Only PipeA", "vaf": 0.03, "depth": 40, "alt_count": 2},
            {"group": "Only PipeA", "vaf": 0.05, "depth": 45, "alt_count": 3},
            {"group": "Only PipeB", "vaf": 0.10, "depth": 70, "alt_count": 7},
            {"group": "Only PipeB", "vaf": 0.12, "depth": 80, "alt_count": 8},
        ]
    )

    summary = summarize_unique_pipeline_groups(uniq_cmp, "PipeA", "PipeB")

    assert summary["metrics"]["Only PipeA"]["count"] == 2
    assert summary["metrics"]["Only PipeB"]["count"] == 2
    assert summary["metrics"]["Only PipeA"]["median_vaf"] == 0.04
    assert summary["metrics"]["Only PipeB"]["median_depth"] == 75.0
    assert "PipeA unique calls have lower median VAF" in summary["summary"]

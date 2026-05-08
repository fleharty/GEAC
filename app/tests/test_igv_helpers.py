"""Tests for IGV helper logic.

These tests cover pure-Python / DuckDB logic with no Streamlit dependency.
Run with: pytest app/tests/
"""

import sys
import os

import duckdb
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from igv_helpers import make_bed, query_distinct_samples, resolve_index_uri


# ── Fixtures ──────────────────────────────────────────────────────────────────

def _make_db(n_samples: int, rows_per_sample: int = 10) -> duckdb.DuckDBPyConnection:
    """Build an in-memory DuckDB with a synthetic alt_bases table."""
    con = duckdb.connect(":memory:")
    con.execute("""
        CREATE TABLE alt_bases (
            sample_id    VARCHAR,
            chrom        VARCHAR,
            pos          INTEGER,
            ref_allele   VARCHAR,
            alt_allele   VARCHAR,
            variant_type VARCHAR,
            alt_count    INTEGER,
            total_depth  INTEGER
        )
    """)
    rows = [
        (f"SAMPLE_{i:03d}", "chr1", pos * 1000, "A", "T", "SNV", 3, 100)
        for i in range(n_samples)
        for pos in range(rows_per_sample)
    ]
    con.executemany("INSERT INTO alt_bases VALUES (?,?,?,?,?,?,?,?)", rows)
    return con


# ── Tests ─────────────────────────────────────────────────────────────────────

class TestQueryDistinctSamples:

    def test_returns_all_samples(self):
        """Should return every sample regardless of row count."""
        con = _make_db(n_samples=10, rows_per_sample=10)
        result = query_distinct_samples(con, "alt_bases", "1=1")
        assert len(result) == 10

    def test_regression_display_limit_does_not_hide_samples(self):
        """Regression: IGV cap warning was suppressed because sample IDs were
        derived from a display-limited DataFrame (e.g. 100 rows) rather than
        queried from the full dataset.

        With 10 samples × 10 rows = 100 total rows, a display_limit=5 would
        return only 5 rows — potentially all from one sample — causing the
        cap warning to be silently skipped.

        query_distinct_samples must count samples from the full table, not from
        a pre-sliced DataFrame.
        """
        n_samples = 10
        con = _make_db(n_samples=n_samples, rows_per_sample=10)

        # Simulate what the old (buggy) code did: derive sample IDs from a
        # display-limited slice.  With LIMIT 5, you get at most 1 sample.
        display_df = con.execute(
            "SELECT * FROM alt_bases ORDER BY sample_id LIMIT 5"
        ).df()
        samples_from_display = display_df["sample_id"].unique().tolist()

        # The helper must return the full count, not the display-slice count.
        samples_from_query = query_distinct_samples(con, "alt_bases", "1=1")

        assert len(samples_from_query) == n_samples
        assert len(samples_from_display) < n_samples, (
            "Test setup error: display slice should contain fewer samples than the full dataset"
        )

    def test_where_filter_respected(self):
        """Extra WHERE conditions should narrow the sample list correctly."""
        con = _make_db(n_samples=10, rows_per_sample=2)
        result = query_distinct_samples(
            con, "alt_bases", "sample_id IN ('SAMPLE_000', 'SAMPLE_001')"
        )
        assert sorted(result) == ["SAMPLE_000", "SAMPLE_001"]

    def test_returns_sorted_list(self):
        """Result should be sorted for stable multiselect default ordering."""
        con = _make_db(n_samples=5, rows_per_sample=1)
        result = query_distinct_samples(con, "alt_bases", "1=1")
        assert result == sorted(result)

    def test_empty_result(self):
        """Should return an empty list when no rows match the filter."""
        con = _make_db(n_samples=3, rows_per_sample=1)
        result = query_distinct_samples(con, "alt_bases", "sample_id = 'DOES_NOT_EXIST'")
        assert result == []


class TestResolveIndexUri:

    def test_infers_bam_and_cram_indexes(self):
        assert resolve_index_uri("sample.bam", None) == "sample.bam.bai"
        assert resolve_index_uri("sample.cram", None) == "sample.cram.crai"

    def test_infers_variant_indexes(self):
        assert resolve_index_uri("gnomad.vcf.gz", None) == "gnomad.vcf.gz.tbi"
        assert resolve_index_uri("gnomad.vcf.bgz", None) == "gnomad.vcf.bgz.tbi"
        assert resolve_index_uri("gnomad.bcf", None) == "gnomad.bcf.csi"

    def test_prefers_explicit_index_when_provided(self):
        assert resolve_index_uri("gnomad.vcf.gz", "/tmp/custom.tbi") == "/tmp/custom.tbi"

    def test_infers_absolute_variant_index_when_track_path_is_absolute(self, tmp_path):
        track = str((tmp_path / "refs" / "gnomad.vcf.gz").resolve())
        assert resolve_index_uri(track, None) == track + ".tbi"


class TestMakeBed:
    """Verify BED coordinate generation, especially the deletion off-by-one.

    GEAC coordinate convention for deletions:
        pos       = anchor (0-based reference base before the deletion, VCF style)
        alt_allele = "-" + deleted_bases  (e.g. "-ACG" for a 3-base deletion)

    The deleted bases span pos+1 .. pos+del_len (0-based).
    The correct BED interval to cover anchor + all deleted bases is:
        [pos,  pos + len(alt_allele))
    which equals [pos, pos + 1 + del_len) — the '-' prefix in the alt allele
    length coincidentally supplies the +1 needed for the half-open BED end.

    BED start is always pos (the anchor), matching VCF POS convention.
    In IGV, the deletion gap in reads starts at pos+1 (one base after the anchor).
    """

    def _row(self, chrom, pos, alt_allele, variant_type="SNV"):
        import pandas as pd
        return pd.DataFrame([{
            "chrom": chrom,
            "pos": pos,
            "alt_allele": alt_allele,
            "variant_type": variant_type,
        }])

    def test_snv_single_base_interval(self):
        bed = make_bed(self._row("chr1", 100, "T", "SNV"))
        assert bed == "chr1\t100\t101\n"

    def test_deletion_3bp(self):
        # pos=100 (anchor), deleted bases ACG at positions 101-103 (0-based)
        # BED end = 100 + len("-ACG") = 100 + 4 = 104 (exclusive)
        bed = make_bed(self._row("chr1", 100, "-ACG", "deletion"))
        assert bed == "chr1\t100\t104\n"

    def test_deletion_1bp(self):
        # Single-base deletion: anchor at pos, one deleted base at pos+1
        # BED end = pos + len("-A") = pos + 2
        bed = make_bed(self._row("chr1", 500, "-A", "deletion"))
        assert bed == "chr1\t500\t502\n"

    def test_deletion_covers_anchor_and_all_deleted_bases(self):
        # Verify interval width: anchor (1) + deleted bases (del_len) = len(alt)
        import pandas as pd
        pos, alt = 200, "-TTTTTT"  # 6-base deletion
        df = self._row("chr1", pos, alt, "deletion")
        bed = make_bed(df)
        parts = bed.strip().split("\t")
        start, end = int(parts[1]), int(parts[2])
        del_len = len(alt) - 1  # exclude the '-' prefix
        assert start == pos
        assert end == pos + del_len + 1  # exclusive end covers through last deleted base
        assert end - start == len(alt)   # interval width = anchor + del_len

    def test_multi_locus_sorted_output(self):
        import pandas as pd
        df = pd.DataFrame([
            {"chrom": "chr1", "pos": 300, "alt_allele": "G",    "variant_type": "SNV"},
            {"chrom": "chr1", "pos": 100, "alt_allele": "-AC",  "variant_type": "deletion"},
            {"chrom": "chr2", "pos": 50,  "alt_allele": "T",    "variant_type": "SNV"},
        ])
        lines = make_bed(df).strip().split("\n")
        assert lines[0] == "chr1\t100\t103"  # deletion: 100 + len("-AC") = 103
        assert lines[1] == "chr1\t300\t301"  # SNV
        assert lines[2] == "chr2\t50\t51"    # SNV

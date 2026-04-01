"""Tests for IGV helper logic.

These tests cover pure-Python / DuckDB logic with no Streamlit dependency.
Run with: pytest app/tests/
"""

import sys
import os

import duckdb
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from igv_helpers import query_distinct_samples, resolve_index_uri


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

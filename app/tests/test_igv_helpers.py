"""Tests for IGV helper logic.

These tests cover pure-Python / DuckDB logic with no Streamlit dependency.
Run with: pytest app/tests/
"""

import sys
import os

import duckdb
import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from igv_helpers import make_vcf, query_distinct_samples, resolve_index_uri


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


_MISSING = "./.:.:.:."


class TestMakeVcf:
    """Verify VCF generation: coordinate conversion, allele encoding, multi-sample layout.

    GEAC coordinate convention:
        pos        = 0-based; VCF POS = pos + 1
        ref_allele = single anchor base
        alt_allele = plain base (SNV), "-" + deleted (deletion), "+" + inserted (insertion)
    """

    def _row(self, chrom, pos, ref_allele, alt_allele, variant_type="SNV",
             sample_id="S1", alt_count=3, ref_count=97, total_depth=100, vaf=0.03):
        return pd.DataFrame([{
            "chrom": chrom, "pos": pos, "ref_allele": ref_allele,
            "alt_allele": alt_allele, "variant_type": variant_type,
            "sample_id": sample_id, "alt_count": alt_count,
            "ref_count": ref_count, "total_depth": total_depth, "vaf": vaf,
        }])

    def _data_lines(self, vcf_str):
        return [l for l in vcf_str.splitlines() if l and not l.startswith("#")]

    def _col_header(self, vcf_str):
        for line in vcf_str.splitlines():
            if line.startswith("#CHROM"):
                return line
        return None

    def test_snv_pos_is_one_based(self):
        vcf = make_vcf(self._row("chr1", 100, "A", "T"))
        fields = self._data_lines(vcf)[0].split("\t")
        assert fields[0] == "chr1"
        assert fields[1] == "101"   # 0-based 100 → 1-based 101

    def test_snv_ref_alt_passthrough(self):
        vcf = make_vcf(self._row("chr1", 100, "A", "T"))
        fields = self._data_lines(vcf)[0].split("\t")
        assert fields[3] == "A"
        assert fields[4] == "T"

    def test_deletion_allele_expansion(self):
        # alt_allele="-ACG" → VCF REF=AACG, ALT=A
        vcf = make_vcf(self._row("chr1", 100, "A", "-ACG", variant_type="deletion"))
        fields = self._data_lines(vcf)[0].split("\t")
        assert fields[3] == "AACG"
        assert fields[4] == "A"

    def test_insertion_allele_expansion(self):
        # alt_allele="+ACGT" → VCF REF=A, ALT=AACGT
        vcf = make_vcf(self._row("chr1", 100, "A", "+ACGT", variant_type="insertion"))
        fields = self._data_lines(vcf)[0].split("\t")
        assert fields[3] == "A"
        assert fields[4] == "AACGT"

    def test_format_fields_present(self):
        vcf = make_vcf(self._row("chr1", 100, "A", "T",
                                  alt_count=3, ref_count=97, total_depth=100, vaf=0.03))
        fields = self._data_lines(vcf)[0].split("\t")
        assert fields[8] == "GT:DP:AD:VAF"
        assert fields[9] == "./.:100:97,3:0.03"

    def test_single_sample_column_in_header(self):
        vcf = make_vcf(self._row("chr1", 100, "A", "T", sample_id="MY_SAMPLE"))
        assert self._col_header(vcf).endswith("\tMY_SAMPLE")

    def test_multi_sample_absent_sample_gets_missing_sentinel(self):
        df = pd.DataFrame([
            {"chrom": "chr1", "pos": 100, "ref_allele": "A", "alt_allele": "T",
             "variant_type": "SNV", "sample_id": "S1",
             "alt_count": 5, "ref_count": 95, "total_depth": 100, "vaf": 0.05},
            {"chrom": "chr1", "pos": 200, "ref_allele": "C", "alt_allele": "G",
             "variant_type": "SNV", "sample_id": "S2",
             "alt_count": 2, "ref_count": 98, "total_depth": 100, "vaf": 0.02},
        ])
        vcf = make_vcf(df)
        col_hdr = self._col_header(vcf)
        assert "S1" in col_hdr and "S2" in col_hdr
        data = self._data_lines(vcf)
        assert len(data) == 2
        cols = col_hdr.split("\t")
        s1_col = cols.index("S1")
        s2_col = cols.index("S2")
        # First variant (pos 101): S1 present, S2 absent
        row1 = data[0].split("\t")
        assert row1[1] == "101"
        assert row1[s1_col] != _MISSING
        assert row1[s2_col] == _MISSING
        # Second variant (pos 201): S2 present, S1 absent
        row2 = data[1].split("\t")
        assert row2[1] == "201"
        assert row2[s2_col] != _MISSING
        assert row2[s1_col] == _MISSING

    def test_rows_sorted_by_chrom_then_pos(self):
        df = pd.DataFrame([
            {"chrom": "chr2", "pos": 50,  "ref_allele": "A", "alt_allele": "T",
             "variant_type": "SNV", "sample_id": "S1",
             "alt_count": 1, "ref_count": 9, "total_depth": 10, "vaf": 0.1},
            {"chrom": "chr1", "pos": 300, "ref_allele": "C", "alt_allele": "G",
             "variant_type": "SNV", "sample_id": "S1",
             "alt_count": 1, "ref_count": 9, "total_depth": 10, "vaf": 0.1},
            {"chrom": "chr1", "pos": 100, "ref_allele": "G", "alt_allele": "A",
             "variant_type": "SNV", "sample_id": "S1",
             "alt_count": 1, "ref_count": 9, "total_depth": 10, "vaf": 0.1},
        ])
        data = self._data_lines(make_vcf(df))
        positions = [(r.split("\t")[0], int(r.split("\t")[1])) for r in data]
        assert positions == [("chr1", 101), ("chr1", 301), ("chr2", 51)]

    def test_empty_dataframe_returns_header_only(self):
        vcf = make_vcf(pd.DataFrame())
        assert "##fileformat=VCFv4.2" in vcf
        assert self._data_lines(vcf) == []

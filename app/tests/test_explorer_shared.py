"""Tests for shared Explorer infrastructure."""

import os
import sys
import tempfile

import duckdb
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from explorer import COVERAGE_FILTER_STATE, MAIN_FILTER_STATE, load_schema_manifest
from explorer import schema as schema_module
from explorer.data_source import DataSource
from explorer.tabs import TAB_MODULES


def _write_alt_db(path: str) -> None:
    con = duckdb.connect(path)
    con.execute(
        """
        CREATE TABLE alt_bases (
            sample_id VARCHAR,
            chrom VARCHAR,
            pos BIGINT,
            ref_allele VARCHAR,
            alt_allele VARCHAR,
            variant_type VARCHAR,
            total_depth INTEGER,
            alt_count INTEGER,
            ref_count INTEGER,
            fwd_depth INTEGER,
            rev_depth INTEGER,
            fwd_alt_count INTEGER,
            rev_alt_count INTEGER,
            fwd_ref_count INTEGER,
            rev_ref_count INTEGER,
            overlap_depth INTEGER,
            overlap_alt_agree INTEGER,
            overlap_alt_disagree INTEGER,
            overlap_ref_agree INTEGER,
            read_type VARCHAR,
            pipeline VARCHAR,
            variant_called BOOLEAN
        );
        INSERT INTO alt_bases VALUES (
            'sample1', 'chr1', 10, 'A', 'T', 'SNV',
            100, 3, 97, 50, 50, 2, 1, 48, 49, 10, 1, 0, 9, 'raw', 'raw', TRUE
        );
        CREATE TABLE alt_reads (
            sample_id VARCHAR,
            chrom VARCHAR,
            pos BIGINT,
            alt_allele VARCHAR,
            cycle INTEGER,
            read_length INTEGER,
            is_read1 BOOLEAN,
            base_qual INTEGER,
            map_qual INTEGER
        );
        INSERT INTO alt_reads VALUES ('sample1', 'chr1', 10, 'T', 11, 150, TRUE, 40, 60);
        CREATE TABLE geac_metadata (
            schema_version VARCHAR,
            geac_version VARCHAR,
            created_at TIMESTAMPTZ,
            command_line VARCHAR,
            output_path VARCHAR,
            platform_os VARCHAR,
            platform_arch VARCHAR,
            platform_family VARCHAR,
            n_alt_bases_inputs BIGINT,
            n_alt_reads_inputs BIGINT,
            n_normal_evidence_inputs BIGINT,
            n_pon_evidence_inputs BIGINT,
            n_coverage_inputs BIGINT,
            n_duckdb_inputs BIGINT,
            n_samples BIGINT,
            alt_bases_rows BIGINT,
            alt_reads_rows BIGINT,
            normal_evidence_rows BIGINT,
            pon_evidence_rows BIGINT,
            coverage_rows BIGINT,
            samples_rows BIGINT
        );
        INSERT INTO geac_metadata VALUES (
            'duckdb-v2', '0.3.17', current_timestamp, 'geac merge', '/tmp/cohort.duckdb',
            'macos', 'arm64', 'unix',
            1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1
        );
        """
    )
    con.close()


def _write_coverage_db(path: str) -> None:
    con = duckdb.connect(path)
    con.execute(
        """
        CREATE TABLE coverage (
            sample_id VARCHAR,
            chrom VARCHAR,
            pos BIGINT,
            "end" BIGINT,
            total_depth INTEGER,
            min_depth INTEGER,
            max_depth INTEGER,
            fwd_depth INTEGER,
            rev_depth INTEGER,
            raw_read_depth INTEGER,
            frac_dup DOUBLE,
            overlap_depth INTEGER,
            frac_overlap DOUBLE,
            mean_mapq DOUBLE,
            frac_mapq0 DOUBLE,
            frac_low_mapq DOUBLE,
            mean_base_qual DOUBLE,
            min_base_qual_obs INTEGER,
            max_base_qual_obs INTEGER,
            frac_low_bq DOUBLE,
            mean_insert_size DOUBLE,
            min_insert_size INTEGER,
            max_insert_size INTEGER,
            n_insert_size_obs INTEGER,
            gc_content DOUBLE,
            read_type VARCHAR,
            pipeline VARCHAR,
            bin_n INTEGER,
            gene VARCHAR,
            on_target BOOLEAN
        );
        INSERT INTO coverage VALUES (
            'sample1', 'chr1', 10, 11, 100, 100, 100, 50, 50, 100,
            0.0, 10, 0.1, 60.0, 0.0, 0.0, 35.0, 30, 40, 0.0,
            200.0, 200, 200, 1, 0.5, 'raw', 'raw', 1, 'TP53', TRUE
        );
        CREATE TABLE geac_metadata (
            schema_version VARCHAR,
            geac_version VARCHAR,
            created_at TIMESTAMPTZ,
            command_line VARCHAR,
            output_path VARCHAR,
            platform_os VARCHAR,
            platform_arch VARCHAR,
            platform_family VARCHAR,
            n_alt_bases_inputs BIGINT,
            n_alt_reads_inputs BIGINT,
            n_normal_evidence_inputs BIGINT,
            n_pon_evidence_inputs BIGINT,
            n_coverage_inputs BIGINT,
            n_duckdb_inputs BIGINT,
            n_samples BIGINT,
            alt_bases_rows BIGINT,
            alt_reads_rows BIGINT,
            normal_evidence_rows BIGINT,
            pon_evidence_rows BIGINT,
            coverage_rows BIGINT,
            samples_rows BIGINT
        );
        INSERT INTO geac_metadata VALUES (
            'duckdb-v2', '0.3.17', current_timestamp, 'geac merge', '/tmp/coverage.duckdb',
            'macos', 'arm64', 'unix',
            0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0
        );
        """
    )
    con.close()


class TestSchemaManifest:
    def test_manifest_lists_feature_tables(self):
        manifest = load_schema_manifest()
        assert "alt_reads" in manifest.feature_tables
        assert "coverage" in manifest.tables

    def test_schema_path_supports_packaged_libexec_layout(self, tmp_path, monkeypatch):
        libexec = tmp_path / "libexec"
        explorer_dir = libexec / "explorer"
        schema_dir = libexec / "schema"
        explorer_dir.mkdir(parents=True)
        schema_dir.mkdir()
        packaged_schema = schema_dir / "geac_schema.json"
        packaged_schema.write_text('{"feature_tables":[],"tables":{}}', encoding="utf-8")

        monkeypatch.setattr(schema_module, "__file__", str(explorer_dir / "schema.py"))

        assert schema_module._schema_path() == packaged_schema


class TestDataSource:
    def test_open_alt_bases_detects_optional_tables(self, tmp_path):
        db = tmp_path / "cohort.duckdb"
        _write_alt_db(str(db))

        ds = DataSource.open_alt_bases(str(db))

        assert ds.is_duckdb is True
        assert ds.has_optional_table("alt_reads") is True
        assert ds.has_optional_table("pon_evidence") is False
        assert ds.required_columns_missing() == []
        assert ds.distinct_values("sample_id") == ["sample1"]
        assert ds.db_version == "0.3.17"
        assert ds.db_schema_version == "duckdb-v2"
        assert list(ds.metadata_header()["schema_version"]) == ["duckdb-v2"]
        assert list(ds.metadata_inputs().columns) == []

    def test_open_coverage_uses_same_metadata_and_contract_logic(self, tmp_path):
        db = tmp_path / "coverage.duckdb"
        _write_coverage_db(str(db))

        ds = DataSource.open_coverage(str(db))

        assert ds.is_duckdb is True
        assert ds.required_columns_missing() == []
        assert ds.has_column("bin_n") is True
        assert ds.has_non_null("gene") is True
        assert ds.db_version == "0.3.17"
        assert ds.db_schema_version == "duckdb-v2"
        assert list(ds.metadata_header()["schema_version"]) == ["duckdb-v2"]


class TestFilterState:
    def test_main_filter_reset_applies_overrides(self):
        state = {}

        MAIN_FILTER_STATE.reset(
            state,
            overrides={"sample_recurrence": (1, 5), "family_size_range": (0, 7)},
        )

        assert state["chrom_sel"] == "All"
        assert state["variant_sel"] == ["SNV", "insertion", "deletion"]
        assert state["sample_recurrence"] == (1, 5)
        assert state["family_size_range"] == (0, 7)

    def test_clear_removes_only_known_keys(self):
        state = {"sample_sel": ["s1"], "chrom_sel": "chr1", "keep_me": True}

        COVERAGE_FILTER_STATE.clear(state)

        assert "sample_sel" not in state
        assert "chrom_sel" not in state
        assert state["keep_me"] is True


class TestTabModules:
    def test_each_tab_module_has_non_empty_label(self):
        for module in TAB_MODULES:
            assert hasattr(module, "LABEL"), f"{module.__name__} missing LABEL"
            assert isinstance(module.LABEL, str) and module.LABEL.strip(), (
                f"{module.__name__}.LABEL must be a non-empty string"
            )

    def test_tab_module_labels_are_unique(self):
        labels = [module.LABEL for module in TAB_MODULES]
        assert len(labels) == len(set(labels)), "Duplicate LABEL values in TAB_MODULES"


def _write_alt_parquet(path: str) -> None:
    table = pa.table(
        {
            "sample_id":             ["sample1"],
            "chrom":                 ["chr1"],
            "pos":                   pa.array([10], type=pa.int64()),
            "ref_allele":            ["A"],
            "alt_allele":            ["T"],
            "variant_type":          ["SNV"],
            "total_depth":           pa.array([100], type=pa.int32()),
            "alt_count":             pa.array([3], type=pa.int32()),
            "ref_count":             pa.array([97], type=pa.int32()),
            "fwd_depth":             pa.array([50], type=pa.int32()),
            "rev_depth":             pa.array([50], type=pa.int32()),
            "fwd_alt_count":         pa.array([2], type=pa.int32()),
            "rev_alt_count":         pa.array([1], type=pa.int32()),
            "fwd_ref_count":         pa.array([48], type=pa.int32()),
            "rev_ref_count":         pa.array([49], type=pa.int32()),
            "overlap_depth":         pa.array([10], type=pa.int32()),
            "overlap_alt_agree":     pa.array([1], type=pa.int32()),
            "overlap_alt_disagree":  pa.array([0], type=pa.int32()),
            "overlap_ref_agree":     pa.array([9], type=pa.int32()),
            "read_type":             ["raw"],
            "pipeline":              ["raw"],
            "variant_called":        pa.array([True], type=pa.bool_()),
        }
    )
    pq.write_table(table, path)


class TestDataSourceParquet:
    def test_parquet_mode_detects_columns_and_queries(self, tmp_path):
        parquet_path = str(tmp_path / "cohort.alt_bases.parquet")
        _write_alt_parquet(parquet_path)

        ds = DataSource.open_alt_bases(parquet_path)

        assert ds.is_duckdb is False
        assert ds.has_column("variant_type") is True
        assert ds.has_column("nonexistent_col") is False
        assert ds.distinct_values("sample_id") == ["sample1"]
        assert ds.db_version is None  # no metadata in parquet mode
        assert ds.db_schema_version is None
        assert ds.required_columns_missing() == []

    def test_parquet_has_non_null_batch(self, tmp_path):
        parquet_path = str(tmp_path / "cohort.alt_bases.parquet")
        _write_alt_parquet(parquet_path)

        ds = DataSource.open_alt_bases(parquet_path)
        result = ds.has_non_null_batch(["variant_called", "nonexistent_col"])

        assert result["variant_called"] is True
        assert result["nonexistent_col"] is False

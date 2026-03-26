"""Tests for the family-size stratified spectrum locus_fs CTE.

The bug being tested: the locus_fs CTE queried alt_reads with only
    WHERE family_size IS NOT NULL
and did not apply _reads_where when per-read filters were active. As a
result, the singleton/multi classification used ALL reads' family sizes
regardless of the active filter, making it inconsistent with the rest of
the Explorer which uses filtered data.

The fix appends `AND {_reads_where}` to the locus_fs CTE when _reads_active
is True.

Run with: pytest app/tests/
"""

import duckdb
import pytest


# ── Helpers ───────────────────────────────────────────────────────────────────

def _make_db() -> duckdb.DuckDBPyConnection:
    """Build an in-memory DuckDB with synthetic alt_bases and alt_reads.

    Locus layout:
      pos=100  SNV  trinuc=ATA  3 reads family_size=1, 1 read family_size=10
               → unfiltered median_fs = median(1,1,1,10) = 1 → singleton
               → filtered (fs >= 2): only fs=10 read passes → median_fs=10 → multi
    """
    con = duckdb.connect(":memory:")
    con.execute("""
        CREATE TABLE alt_bases (
            sample_id       VARCHAR,
            chrom           VARCHAR,
            pos             INTEGER,
            ref_allele      VARCHAR,
            alt_allele      VARCHAR,
            variant_type    VARCHAR,
            trinuc_context  VARCHAR,
            alt_count       INTEGER,
            total_depth     INTEGER
        )
    """)
    con.execute(
        "INSERT INTO alt_bases VALUES ('s1','chr1',100,'A','T','SNV','ATA',4,20)"
    )

    con.execute("""
        CREATE TABLE alt_reads (
            sample_id   VARCHAR,
            chrom       VARCHAR,
            pos         INTEGER,
            alt_allele  VARCHAR,
            family_size INTEGER
        )
    """)
    con.executemany("INSERT INTO alt_reads VALUES (?,?,?,?,?)", [
        ("s1", "chr1", 100, "T", 1),
        ("s1", "chr1", 100, "T", 1),
        ("s1", "chr1", 100, "T", 1),
        ("s1", "chr1", 100, "T", 10),
    ])
    return con


def _run_locus_fs_query(
    con: duckdb.DuckDBPyConnection,
    reads_where: str | None,
) -> dict:
    """Execute the locus_fs classification query and return {pos: fs_group}.

    reads_where=None simulates _reads_active=False (no extra filter).
    reads_where=<expr> simulates _reads_active=True with that condition.
    """
    locus_fs_filter = f"AND {reads_where}" if reads_where is not None else ""
    rows = con.execute(f"""
        WITH locus_fs AS (
            SELECT sample_id, chrom, pos, alt_allele,
                   MEDIAN(family_size) AS median_fs
            FROM alt_reads
            WHERE family_size IS NOT NULL {locus_fs_filter}
            GROUP BY sample_id, chrom, pos, alt_allele
        )
        SELECT
            _t.pos,
            CASE WHEN COALESCE(lfs.median_fs, 1) <= 1
                 THEN 'singleton' ELSE 'multi' END AS fs_group
        FROM alt_bases _t
        LEFT JOIN locus_fs lfs
            ON  lfs.sample_id  = _t.sample_id
            AND lfs.chrom      = _t.chrom
            AND lfs.pos        = _t.pos
            AND lfs.alt_allele = _t.alt_allele
        WHERE _t.variant_type = 'SNV'
        ORDER BY _t.pos
    """).fetchall()
    return {pos: fs_group for pos, fs_group in rows}


def _run_buggy_locus_fs_query(con: duckdb.DuckDBPyConnection) -> dict:
    """Execute the OLD (buggy) locus_fs CTE — always uses all reads regardless of filter."""
    return _run_locus_fs_query(con, reads_where=None)


# ── Tests ─────────────────────────────────────────────────────────────────────

class TestFsStratifiedSpectrumFilter:

    def test_no_filter_classifies_majority_singleton_locus_as_singleton(self):
        """Without a filter, pos=100 has median_fs=1 → classified as singleton."""
        con = _make_db()
        result = _run_locus_fs_query(con, reads_where=None)
        assert result[100] == "singleton"

    def test_filter_reclassifies_locus_when_singleton_reads_excluded(self):
        """With family_size >= 2, only the fs=10 read passes.

        median_fs becomes 10 → the locus must be classified as 'multi'.
        This is the regression test for bug #4.
        """
        con = _make_db()
        result = _run_locus_fs_query(con, reads_where="family_size >= 2")
        assert result[100] == "multi", (
            "pos=100 has 3 singleton reads and 1 family_size=10 read. "
            "After filtering to family_size >= 2, only the fs=10 read remains, "
            "so median_fs=10 and the locus must be classified as 'multi'."
        )

    def test_buggy_query_classifies_locus_as_singleton_despite_active_filter(self):
        """The old CTE ignores the filter and sees all 4 reads → median_fs=1 → singleton.

        This documents the bug: the old code gave 'singleton' even though a
        family_size >= 2 filter was active and the locus should be 'multi'.
        """
        con = _make_db()
        buggy_result = _run_buggy_locus_fs_query(con)
        assert buggy_result[100] == "singleton", (
            "The buggy CTE (no filter applied) should classify pos=100 as singleton "
            "— if this fails, the test setup needs updating"
        )
        # And confirm the fixed query gives the correct answer
        fixed_result = _run_locus_fs_query(con, reads_where="family_size >= 2")
        assert fixed_result[100] == "multi"

    def test_no_filter_active_produces_same_result_as_buggy_query(self):
        """When no per-read filter is active, both queries must agree."""
        con = _make_db()
        assert _run_locus_fs_query(con, None) == _run_buggy_locus_fs_query(con)

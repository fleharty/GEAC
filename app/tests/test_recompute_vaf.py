"""Tests for the re-aggregation (recompute_vaf) table_expr query.

The fix being tested: when recompute_vaf=True and a SNV locus has reads in
alt_reads but NONE pass the per-read filter, the locus should show alt_count=0,
not the original alt_count.

Root cause of the bug: the old query used
    COUNT(*) ... WHERE <filter>   (returns no row when nothing passes)
followed by COALESCE(filtered_alt_count, ab.alt_count).  Since "no row" and
"no reads at all" (indels) both produce NULL from the LEFT JOIN, the COALESCE
treated them identically and fell back to the original alt_count for both.

Fix: use COUNT(*) FILTER (WHERE ...) with a separate TRUE AS has_reads flag so
the CASE can distinguish the two NULL sources.

Run with: pytest app/tests/
"""

import duckdb
import pytest


# ── Helpers ───────────────────────────────────────────────────────────────────

def _make_db() -> duckdb.DuckDBPyConnection:
    """Build an in-memory DuckDB with synthetic alt_bases and alt_reads tables.

    Locus layout:
      pos=100  SNV   alt_count=5  reads: 3× family_size=1, 2× family_size=3
      pos=200  SNV   alt_count=4  reads: 4× family_size=1  (all fail filter)
      pos=300  ins   alt_count=3  no rows in alt_reads      (indel, pass-through)
    """
    con = duckdb.connect(":memory:")

    con.execute("""
        CREATE TABLE alt_bases (
            sample_id   VARCHAR,
            chrom       VARCHAR,
            pos         INTEGER,
            ref_allele  VARCHAR,
            alt_allele  VARCHAR,
            variant_type VARCHAR,
            alt_count   INTEGER,
            total_depth INTEGER
        )
    """)
    con.executemany("INSERT INTO alt_bases VALUES (?,?,?,?,?,?,?,?)", [
        ("s1", "chr1", 100, "A", "T",    "SNV",       5, 20),
        ("s1", "chr1", 200, "C", "G",    "SNV",       4, 20),
        ("s1", "chr1", 300, "A", "+ATG", "insertion", 3, 20),
    ])

    con.execute("""
        CREATE TABLE alt_reads (
            sample_id   VARCHAR,
            chrom       VARCHAR,
            pos         INTEGER,
            alt_allele  VARCHAR,
            family_size INTEGER
        )
    """)
    # pos=100: 3 singleton reads + 2 family-size-3 reads
    con.executemany("INSERT INTO alt_reads VALUES (?,?,?,?,?)", [
        ("s1", "chr1", 100, "T", 1),
        ("s1", "chr1", 100, "T", 1),
        ("s1", "chr1", 100, "T", 1),
        ("s1", "chr1", 100, "T", 3),
        ("s1", "chr1", 100, "T", 3),
        # pos=200: all singletons — none will pass family_size >= 2
        ("s1", "chr1", 200, "G", 1),
        ("s1", "chr1", 200, "G", 1),
        ("s1", "chr1", 200, "G", 1),
        ("s1", "chr1", 200, "G", 1),
        # pos=300 (insertion) intentionally has no rows
    ])
    return con


def _run_fixed_query(con: duckdb.DuckDBPyConnection, reads_where: str) -> dict:
    """Execute the fixed re-aggregation table_expr and return {pos: alt_count}."""
    rows = con.execute(f"""
        SELECT ab.pos, ab.alt_allele,
            CASE
                WHEN ar_agg.has_reads IS NULL THEN ab.alt_count
                ELSE COALESCE(ar_agg.filtered_alt_count, 0)
            END AS alt_count
        FROM alt_bases ab
        LEFT JOIN (
            SELECT sample_id, chrom, pos, alt_allele,
                   COUNT(*) FILTER (WHERE {reads_where}) AS filtered_alt_count,
                   TRUE AS has_reads
            FROM alt_reads
            GROUP BY sample_id, chrom, pos, alt_allele
        ) ar_agg ON ab.sample_id = ar_agg.sample_id
                 AND ab.chrom    = ar_agg.chrom
                 AND ab.pos      = ar_agg.pos
                 AND ab.alt_allele = ar_agg.alt_allele
        ORDER BY ab.pos
    """).fetchall()
    return {pos: alt_count for pos, _, alt_count in rows}


def _run_buggy_query(con: duckdb.DuckDBPyConnection, reads_where: str) -> dict:
    """Execute the OLD (buggy) re-aggregation query and return {pos: alt_count}.

    The bug: COUNT(*) WHERE filter returns no row when nothing passes, so
    COALESCE(NULL, ab.alt_count) silently falls back to the original count.
    """
    rows = con.execute(f"""
        SELECT ab.pos, ab.alt_allele,
            COALESCE(ar_agg.filtered_alt_count, ab.alt_count) AS alt_count
        FROM alt_bases ab
        LEFT JOIN (
            SELECT sample_id, chrom, pos, alt_allele,
                   COUNT(*) AS filtered_alt_count
            FROM alt_reads
            WHERE {reads_where}
            GROUP BY sample_id, chrom, pos, alt_allele
        ) ar_agg ON ab.sample_id = ar_agg.sample_id
                 AND ab.chrom    = ar_agg.chrom
                 AND ab.pos      = ar_agg.pos
                 AND ab.alt_allele = ar_agg.alt_allele
        ORDER BY ab.pos
    """).fetchall()
    return {pos: alt_count for pos, _, alt_count in rows}


# ── Tests ─────────────────────────────────────────────────────────────────────

class TestRecomputeVafCoalesce:

    def test_snv_with_passing_reads_shows_filtered_count(self):
        """SNV locus where some reads pass: alt_count = number of passing reads."""
        con = _make_db()
        result = _run_fixed_query(con, "family_size >= 2")
        assert result[100] == 2, (
            "pos=100 has 2 reads with family_size=3; filtered count should be 2"
        )

    def test_snv_with_no_passing_reads_shows_zero(self):
        """SNV locus where ALL reads fail the filter: alt_count must be 0, not original.

        This is the regression test for the COALESCE bug. Before the fix,
        COALESCE(NULL, ab.alt_count) would return 4 (the original count).
        """
        con = _make_db()
        result = _run_fixed_query(con, "family_size >= 2")
        assert result[200] == 0, (
            "pos=200 has 4 reads, all with family_size=1; none pass family_size>=2, "
            "so filtered alt_count must be 0"
        )

    def test_indel_with_no_alt_reads_rows_preserves_original_count(self):
        """Indel locus with no alt_reads rows: original alt_count is preserved."""
        con = _make_db()
        result = _run_fixed_query(con, "family_size >= 2")
        assert result[300] == 3, (
            "pos=300 is an insertion with no alt_reads rows; "
            "original alt_count=3 should be preserved"
        )

    def test_buggy_query_would_have_failed(self):
        """Confirm that the old query returns the wrong answer for pos=200.

        This documents why the fix was necessary: the buggy COALESCE gives
        the original alt_count (4) instead of 0 when all reads fail the filter.
        """
        con = _make_db()
        buggy_result = _run_buggy_query(con, "family_size >= 2")
        assert buggy_result[200] == 4, (
            "The buggy query should return the original alt_count (4) for pos=200 "
            "— if this fails, the buggy query has been corrected and this test "
            "needs updating"
        )
        # And confirm the fixed query gives 0 for the same input
        fixed_result = _run_fixed_query(con, "family_size >= 2")
        assert fixed_result[200] == 0

    def test_no_filter_active_all_loci_use_original_count(self):
        """With a trivially true filter, every locus reflects its original alt_count."""
        con = _make_db()
        result = _run_fixed_query(con, "TRUE")
        assert result[100] == 5
        assert result[200] == 4
        assert result[300] == 3

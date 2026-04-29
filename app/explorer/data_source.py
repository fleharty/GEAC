from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Iterable, Optional

import duckdb
import pandas as pd

from .schema import SchemaManifest, load_schema_manifest


def _sql_str(value: str) -> str:
    return value.replace("'", "''")


@dataclass
class GenomicRegion:
    """Parsed result of a combined region-or-gene string.

    Exactly one of (chrom, gene) will be non-None, or both are None for 'All'.
    """
    chrom: Optional[str]  # None means no chrom filter
    start: Optional[int]  # 0-based; None means no position filter
    end: Optional[int]    # 0-based inclusive; None means no end bound
    gene: Optional[str]   # None means no gene filter; set for exact gene match


def parse_region(text: str, known_chroms: list[str]) -> GenomicRegion:
    """Parse a region-or-gene string into a GenomicRegion.

    Accepted formats:
      - "" or "All"              → no filter
      - "1" or "chr1"            → chromosome only
      - "1:1000" or "chr1:1,000" → single position (1-based, converted to 0-based)
      - "1:1000-2000"            → range (1-based inclusive, converted to 0-based)
      - "TP53", "BRCA1", …       → exact gene match (anything not resolving as a chrom)

    Normalises the chromosome name against known_chroms, trying with/without
    the "chr" prefix so "1" matches "chr1" and vice versa.
    """
    text = text.strip()
    if not text or text.lower() == "all":
        return GenomicRegion(None, None, None, None)

    if ":" in text:
        chrom_part, pos_part = text.split(":", 1)
        pos_part = pos_part.strip().replace(",", "")
        m = re.match(r"^(\d+)(?:[:\-](\d+))?$", pos_part)
        # User input is 1-based (IGV/VCF convention); DB stores 0-based positions.
        start = int(m.group(1)) - 1 if m else None
        end   = int(m.group(2)) - 1 if m and m.group(2) else None
        chrom = _resolve_chrom(chrom_part.strip(), known_chroms)
        return GenomicRegion(chrom, start, end, None)

    # No colon — try to resolve as a chromosome name first.
    resolved = _try_resolve_chrom(text, known_chroms)
    if resolved is not None:
        return GenomicRegion(resolved, None, None, None)

    # Not a known chromosome — treat as an exact gene name.
    return GenomicRegion(None, None, None, text)


def _try_resolve_chrom(name: str, known_chroms: list[str]) -> Optional[str]:
    """Return the matching chrom from known_chroms, or None if not found."""
    if name in known_chroms:
        return name
    with_chr = "chr" + name
    if with_chr in known_chroms:
        return with_chr
    without_chr = name.removeprefix("chr")
    if without_chr in known_chroms:
        return without_chr
    return None


def _resolve_chrom(name: str, known_chroms: list[str]) -> str:
    """Match name against known_chroms, trying with/without 'chr' prefix.

    Falls back to returning the name as-is (for region strings with a colon
    where the chrom part should pass through even if unknown).
    """
    return _try_resolve_chrom(name, known_chroms) or name


def sort_chroms(chroms: list[str]) -> list[str]:
    """Sort chromosome names in natural genomic order.

    Handles both bare names (1, 2, X) and prefixed names (chr1, chr2, chrX).
    Order: numeric contigs by number, then X, Y, M/MT, then anything else
    lexicographically.
    """
    _SPECIAL = {"x": 100, "y": 101, "m": 102, "mt": 102}

    def _key(name: str) -> tuple:
        bare = name.lower().removeprefix("chr")
        if bare.isdigit():
            return (0, int(bare), name)
        special = _SPECIAL.get(bare)
        if special is not None:
            return (1, special, name)
        return (2, 0, name)

    return sorted(chroms, key=_key)


@dataclass
class DataSource:
    path: str
    table_name: str
    manifest: SchemaManifest = field(default_factory=load_schema_manifest)
    con: duckdb.DuckDBPyConnection = field(init=False)
    table_expr: str = field(init=False)
    is_duckdb: bool = field(init=False)
    schema_cols: set[str] = field(init=False)
    available_tables: set[str] = field(init=False, default_factory=set)
    db_version: str | None = field(init=False, default=None)
    db_created: object | None = field(init=False, default=None)
    db_schema_version: str | None = field(init=False, default=None)

    @classmethod
    def open_alt_bases(cls, path: str) -> "DataSource":
        return cls(path=path, table_name="alt_bases")

    @classmethod
    def open_coverage(cls, path: str) -> "DataSource":
        return cls(path=path, table_name="coverage")

    def __post_init__(self) -> None:
        self.is_duckdb = self.path.endswith(".duckdb")
        if self.is_duckdb:
            self.con = duckdb.connect(self.path, read_only=True)
            self.table_expr = self.table_name
            self.available_tables = {
                row[0]
                for row in self.con.execute(
                    "SELECT table_name FROM information_schema.tables WHERE table_schema = 'main'"
                    " UNION ALL "
                    "SELECT table_name FROM information_schema.views WHERE table_schema = 'main'"
                ).fetchall()
            }
            if self.table_name not in self.available_tables:
                raise ValueError(
                    f"No `{self.table_name}` table found in this DuckDB."
                )
            self._load_metadata()
        else:
            escaped = _sql_str(self.path)
            self.con = duckdb.connect()
            self.table_expr = f"read_parquet('{escaped}', union_by_name=true)"

        self.schema_cols = set(
            self.con.execute(f"DESCRIBE SELECT * FROM {self.table_expr} LIMIT 0")
            .df()["column_name"]
            .tolist()
        )

    def _load_metadata(self) -> None:
        if "geac_metadata" not in self.available_tables:
            return
        metadata_cols = set(
            self.con.execute("DESCRIBE SELECT * FROM geac_metadata LIMIT 0")
            .df()["column_name"]
            .tolist()
        )
        select_cols = ["geac_version", "created_at"]
        if "schema_version" in metadata_cols:
            select_cols.append("schema_version")
        row = self.con.execute(
            f"SELECT {', '.join(select_cols)} FROM geac_metadata LIMIT 1"
        ).fetchone()
        if row:
            self.db_version = row[0]
            self.db_created = row[1]
            if len(row) > 2:
                self.db_schema_version = row[2]

    def table_exists(self, table: str) -> bool:
        return table in self.available_tables

    def has_optional_table(self, table: str) -> bool:
        return table in self.manifest.feature_tables and self.table_exists(table)

    def has_column(self, column: str) -> bool:
        return column in self.schema_cols

    def has_non_null(self, column: str) -> bool:
        if column not in self.schema_cols:
            return False
        count = self.con.execute(
            f"SELECT COUNT(*) FROM {self.table_expr} WHERE {column} IS NOT NULL"
        ).fetchone()[0]
        return count > 0

    def has_non_null_batch(self, columns: Iterable[str]) -> dict[str, bool]:
        """Check multiple columns for non-null values in a single query.

        More efficient than calling has_non_null() repeatedly when many optional
        columns need to be probed at startup.
        """
        columns = list(columns)
        present = [c for c in columns if c in self.schema_cols]
        result: dict[str, bool] = {c: False for c in columns}
        if not present:
            return result
        parts = [
            f"COUNT(*) FILTER (WHERE {c} IS NOT NULL) > 0 AS \"{c}\""
            for c in present
        ]
        row = self.con.execute(
            f"SELECT {', '.join(parts)} FROM {self.table_expr}"
        ).fetchone()
        for col, val in zip(present, row):
            result[col] = bool(val)
        return result

    def distinct_values(
        self,
        column: str,
        *,
        not_null: bool = True,
        extra_where: Iterable[str] | None = None,
    ) -> list[object]:
        clauses = []
        if not_null:
            clauses.append(f"{column} IS NOT NULL")
        if extra_where:
            clauses.extend(extra_where)
        where = f"WHERE {' AND '.join(clauses)}" if clauses else ""
        return (
            self.con.execute(
                f"SELECT DISTINCT {column} FROM {self.table_expr} {where} ORDER BY {column}"
            )
            .df()[column]
            .tolist()
        )

    def required_columns_missing(self) -> list[str]:
        contract = self.manifest.table(self.table_name)
        return [
            column for column in contract.required_columns if column not in self.schema_cols
        ]

    def summary_stats(self) -> pd.DataFrame:
        return self.con.execute(
            f"""
            SELECT
                COUNT(*) AS n_records,
                COUNT(DISTINCT sample_id) AS n_samples,
                SUM(alt_count) AS total_alt_bases,
                ROUND(AVG(alt_count * 1.0 / total_depth), 4) AS mean_vaf,
                ROUND(AVG(total_depth), 1) AS mean_depth,
                COUNT(*) FILTER (WHERE variant_called IS NOT NULL) AS n_annotated,
                COUNT(*) FILTER (WHERE variant_called = true) AS n_called
            FROM {self.table_expr}
            """
        ).df()

    def metadata_header(self) -> pd.DataFrame:
        if not self.is_duckdb or "geac_metadata" not in self.available_tables:
            return pd.DataFrame()
        return self.con.execute("SELECT * FROM geac_metadata LIMIT 1").df()

    def metadata_inputs(self) -> pd.DataFrame:
        if not self.is_duckdb or "geac_inputs" not in self.available_tables:
            return pd.DataFrame()
        return self.con.execute(
            """
            SELECT *
            FROM geac_inputs
            ORDER BY input_kind, input_path
            """
        ).df()

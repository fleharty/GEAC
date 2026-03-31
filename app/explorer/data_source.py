from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

import duckdb
import pandas as pd

from .schema import SchemaManifest, load_schema_manifest


def _sql_str(value: str) -> str:
    return value.replace("'", "''")


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
        row = self.con.execute(
            "SELECT geac_version, created_at FROM geac_metadata LIMIT 1"
        ).fetchone()
        if row:
            self.db_version, self.db_created = row

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

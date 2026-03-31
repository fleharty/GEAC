import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path


GEAC_VERSION = "0.3.16"


@dataclass(frozen=True)
class TableSchema:
    name: str
    required_columns: tuple[str, ...]
    optional_columns: tuple[str, ...]

    @property
    def known_columns(self) -> tuple[str, ...]:
        return self.required_columns + self.optional_columns


@dataclass(frozen=True)
class SchemaManifest:
    feature_tables: tuple[str, ...]
    tables: dict[str, TableSchema]

    def table(self, name: str) -> TableSchema:
        return self.tables[name]


def _schema_path() -> Path:
    return Path(__file__).resolve().parents[2] / "schema" / "geac_schema.json"


@lru_cache(maxsize=1)
def load_schema_manifest() -> SchemaManifest:
    raw = json.loads(_schema_path().read_text())
    tables = {
        name: TableSchema(
            name=name,
            required_columns=tuple(spec["required_columns"]),
            optional_columns=tuple(spec["optional_columns"]),
        )
        for name, spec in raw["tables"].items()
    }
    return SchemaManifest(
        feature_tables=tuple(raw["feature_tables"]),
        tables=tables,
    )

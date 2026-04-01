import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path


GEAC_VERSION = "0.3.18"


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
    here = Path(__file__).resolve()
    candidates = [
        here.parents[2] / "schema" / "geac_schema.json",  # source checkout layout
        here.parents[1] / "schema" / "geac_schema.json",  # packaged libexec layout
        here.parent / "geac_schema.json",                 # fallback for direct bundling
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        "Could not find geac_schema.json. Looked in: "
        + ", ".join(str(path) for path in candidates)
    )


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

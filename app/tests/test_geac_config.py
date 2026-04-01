"""Tests for geac.toml config loading and path normalization."""

from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import geac_config


def _write_config(path: Path, body: str) -> None:
    path.write_text(body, encoding="utf-8")


def _load_from_config(monkeypatch, config_path: Path, cwd: Path) -> dict:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "streamlit",
            "run",
            "app/geac_explorer.py",
            "--",
            "--config",
            str(config_path),
        ],
    )
    monkeypatch.chdir(cwd)
    return geac_config.load()


def test_relative_paths_are_resolved_against_config_directory(monkeypatch, tmp_path):
    config_dir = tmp_path / "config"
    config_dir.mkdir()
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    config_path = config_dir / "geac.toml"
    _write_config(
        config_path,
        """
data = "data/cohort.duckdb"
manifest = "manifests/manifest.tsv"
cosmic = "refs/COSMIC.txt"
target_regions = "refs/targets.bed"
gnomad_track = "refs/gnomad.vcf.gz"
gnomad_track_index = "refs/gnomad.vcf.gz.tbi"
genome_build = "hg38"
auto_launch_igv = true
""".strip(),
    )

    cfg = _load_from_config(monkeypatch, config_path, work_dir)

    assert cfg["data"] == str((config_dir / "data" / "cohort.duckdb").resolve())
    assert cfg["manifest"] == str((config_dir / "manifests" / "manifest.tsv").resolve())
    assert cfg["cosmic"] == str((config_dir / "refs" / "COSMIC.txt").resolve())
    assert cfg["target_regions"] == str((config_dir / "refs" / "targets.bed").resolve())
    assert cfg["gnomad_track"] == str((config_dir / "refs" / "gnomad.vcf.gz").resolve())
    assert cfg["gnomad_track_index"] == str(
        (config_dir / "refs" / "gnomad.vcf.gz.tbi").resolve()
    )
    assert cfg["genome_build"] == "hg38"
    assert cfg["auto_launch_igv"] is True


def test_absolute_and_uri_paths_are_preserved(monkeypatch, tmp_path):
    config_dir = tmp_path / "config"
    config_dir.mkdir()
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    config_path = config_dir / "geac.toml"
    absolute_data = (tmp_path / "data" / "cohort.duckdb").resolve()
    _write_config(
        config_path,
        f"""
data = "{absolute_data}"
manifest = "gs://bucket/project/manifest.tsv"
cosmic = "https://example.com/COSMIC.txt"
target_regions = "http://example.com/targets.bed"
gnomad_track = "gs://bucket/refs/gnomad.vcf.gz"
gnomad_track_index = "https://example.com/gnomad.vcf.gz.tbi"
""".strip(),
    )

    cfg = _load_from_config(monkeypatch, config_path, work_dir)

    assert cfg["data"] == str(absolute_data)
    assert cfg["manifest"] == "gs://bucket/project/manifest.tsv"
    assert cfg["cosmic"] == "https://example.com/COSMIC.txt"
    assert cfg["target_regions"] == "http://example.com/targets.bed"
    assert cfg["gnomad_track"] == "gs://bucket/refs/gnomad.vcf.gz"
    assert cfg["gnomad_track_index"] == "https://example.com/gnomad.vcf.gz.tbi"


def test_config_path_beats_current_working_directory(monkeypatch, tmp_path):
    config_dir = tmp_path / "project"
    config_dir.mkdir()
    other_dir = tmp_path / "elsewhere"
    other_dir.mkdir()
    (other_dir / "geac.toml").write_text('data = "wrong/place.duckdb"\n', encoding="utf-8")

    config_path = config_dir / "geac.toml"
    _write_config(config_path, 'data = "correct/place.duckdb"\n')

    cfg = _load_from_config(monkeypatch, config_path, other_dir)

    assert cfg["data"] == str((config_dir / "correct" / "place.duckdb").resolve())

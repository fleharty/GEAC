from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from typing import Any, Mapping, MutableMapping


@dataclass(frozen=True)
class FilterState:
    keys: tuple[str, ...]
    defaults: Mapping[str, Any] = field(default_factory=dict)

    def reset(
        self,
        session_state: MutableMapping[str, Any],
        overrides: Mapping[str, Any] | None = None,
    ) -> None:
        values = {key: deepcopy(value) for key, value in self.defaults.items()}
        if overrides:
            values.update({key: deepcopy(value) for key, value in overrides.items()})
        for key, value in values.items():
            session_state[key] = value

    def clear(self, session_state: MutableMapping[str, Any]) -> None:
        for key in self.keys:
            session_state.pop(key, None)


MAIN_FILTER_KEYS = (
    "chrom_sel",
    "sample_sel",
    "sample_recurrence",
    "batch_sel",
    "label1_sel",
    "label2_sel",
    "label3_sel",
    "gene_text",
    "variant_sel",
    "vaf_range",
    "min_alt",
    "min_fwd_alt",
    "min_rev_alt",
    "min_overlap_agree",
    "min_overlap_disagree",
    "variant_called_sel",
    "variant_filter_sel",
    "on_target_sel",
    "gnomad_af_range",
    "gnomad_include_null",
    "homopolymer_range",
    "str_len_range",
    "min_depth",
    "max_depth",
    "table_limit_sel",
    "recompute_vaf",
    "family_size_range",
    "cycle_range",
    "map_qual_range",
    "insert_size_range",
    "fs_exclude_mode",
    "cycle_exclude_mode",
    "mq_exclude_mode",
    "is_exclude_mode",
    "read_strand_sel",
)

MAIN_TAB_UI_KEYS = (
    "sb_scale",
    "sb_color_by",
    "sb_show_all",
    "dfe_color_by",
    "dfe_y_mode",
    "dfe_show_r1r2",
    "bq_color_by",
    "bq_show_r1r2",
    "ins_color_by",
    "ins_y_mode",
    "af_ins_color_by",
    "af_ins_y_mode",
    "fs_color_by",
    "fs_y_mode",
    "fs_x_range",
    "sbs_y_mode",
    "top_n_sig",
    "cmp_top_n",
)

MAIN_FILTER_STATE = FilterState(
    keys=MAIN_FILTER_KEYS,
    defaults={
        "chrom_sel": "All",
        "sample_sel": [],
        "batch_sel": [],
        "label1_sel": [],
        "label2_sel": [],
        "label3_sel": [],
        "gene_text": "",
        "variant_sel": ["SNV", "insertion", "deletion"],
        "vaf_range": (0.0, 1.0),
        "min_alt": 1,
        "max_alt": 0,
        "min_fwd_alt": 0,
        "min_rev_alt": 0,
        "min_overlap_agree": 0,
        "min_overlap_disagree": 0,
        "variant_called_sel": "All",
        "variant_filter_sel": [],
        "on_target_sel": "All",
        "gnomad_af_range": ("0", "1.0"),
        "gnomad_include_null": True,
        "homopolymer_range": (0, 20),
        "str_len_range": (0, 50),
        "min_depth": 0,
        "max_depth": 0,
        "table_limit_sel": 500,
        "recompute_vaf": False,
        "fs_exclude_mode": False,
        "cycle_exclude_mode": False,
        "mq_exclude_mode": False,
        "is_exclude_mode": False,
        "read_strand_sel": "All",
    },
)

COVERAGE_FILTER_KEYS = (
    "sample_sel",
    "chrom_sel",
    "gene_text",
    "on_target_sel",
)

COVERAGE_FILTER_STATE = FilterState(keys=COVERAGE_FILTER_KEYS)

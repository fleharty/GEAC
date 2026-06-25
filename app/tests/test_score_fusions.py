"""Tests for the fusion benchmarking scorer (experimental/fusion_benchmark/score_fusions.py).

Pure-stdlib logic; covers the TP/FN/FP scoring, the expected-level semantics, negative-control
FPs, orientation accuracy, and the recall-by-dilution breakdown. Run by CI via pytest app/tests/.
"""

from __future__ import annotations

import os
import sys

sys.path.insert(
    0, os.path.join(os.path.dirname(__file__), "..", "..", "experimental", "fusion_benchmark")
)
import score_fusions as sf


def T(sid, a, b, exp, dil=""):
    return sf.Truth({"sample_id": sid, "gene_a": a, "gene_b": b, "expected": exp, "dilution": dil})


def C(a, b, filt, po="index", supp=10):
    return sf.Call(
        {"gene_a": a, "gene_b": b, "filter": filt, "partner_order": po, "supporting_reads": str(supp)}
    )


def test_self_test_passes():
    # The module's own synthetic check should hold.
    sf.self_test()


def test_pass_requires_PASS_filter():
    truth = {"s": [T("s", "EWSR1", "FLI1", "pass")]}
    # a tagged (chimera) match is detected but does NOT satisfy expected=pass
    _, r = sf.score(truth, {"s": [C("EWSR1", "FLI1", "chimera")]})
    assert r["tp_pass"] == 0 and r["fn_pass"] == 1


def test_detect_accepts_tagged_match():
    truth = {"s": [T("s", "PAX3", "FOXO1", "detect")]}
    _, r = sf.score(truth, {"s": [C("PAX3", "FOXO1", "samelocus")]})
    assert r["tp_detect"] == 1 and r["fn_detect"] == 0


def test_detection_match_is_unordered_and_case_insensitive():
    truth = {"s": [T("s", "EWSR1", "FLI1", "pass")]}
    _, r = sf.score(truth, {"s": [C("fli1", "ewsr1", "PASS")]})
    assert r["tp_pass"] == 1


def test_negative_control_pass_call_is_fp():
    truth = {"s": [T("s", "NONE", "NONE", "negative")]}
    _, r = sf.score(truth, {"s": [C("BCR", "ABL1", "PASS")]})
    assert r["fp_pass"] == 1


def test_offtarget_pass_on_real_sample_is_fp_not_fn():
    truth = {"s": [T("s", "EWSR1", "FLI1", "pass")]}
    # one correct PASS (TP) + one off-target PASS (FP)
    calls = {"s": [C("EWSR1", "FLI1", "PASS"), C("NCOA1", "TERT", "PASS")]}
    _, r = sf.score(truth, calls)
    assert r["tp_pass"] == 1 and r["fp_pass"] == 1 and r["fn_pass"] == 0


def test_orientation_accuracy_only_counts_5to3():
    truth = {
        "ok": [T("ok", "EWSR1", "FLI1", "pass")],
        "rev": [T("rev", "EWSR1", "FLI1", "pass")],
        "abstain": [T("abstain", "EWSR1", "FLI1", "pass")],
    }
    calls = {
        "ok": [C("EWSR1", "FLI1", "PASS", po="5to3")],       # correct order
        "rev": [C("FLI1", "EWSR1", "PASS", po="5to3")],      # asserted but reversed -> miss
        "abstain": [C("FLI1", "EWSR1", "PASS", po="index")],  # abstained -> not counted
    }
    _, r = sf.score(truth, calls)
    assert r["orient_total"] == 2 and r["orient_ok"] == 1


def test_recall_by_dilution():
    truth = {
        "hi": [T("hi", "EWSR1", "FLI1", "pass", "100")],
        "lo": [T("lo", "EWSR1", "FLI1", "pass", "0.1")],
    }
    calls = {"hi": [C("EWSR1", "FLI1", "PASS")], "lo": [C("EWSR1", "FLI1", "chimera")]}
    _, r = sf.score(truth, calls)
    assert r["by_dilution"]["100"] == [1, 0]
    assert r["by_dilution"]["0.1"] == [0, 1]


def test_present_miss_not_counted_as_fn():
    truth = {"s": [T("s", "CIC", "DUX4", "present")]}
    _, r = sf.score(truth, {"s": []})  # not detected
    assert r["fn_pass"] == 0 and r["fn_detect"] == 0  # present-miss is not a headline FN

#!/usr/bin/env python3
"""Score GEAC fusion calls against a known-truth manifest.

Phase 1 of the fusion benchmarking harness: given a truth manifest and a directory of
per-sample ``*.fusions.tsv`` files (from ``geac experimental fusions``), report TP / FN / FP
per sample and a cohort rollup, stratified by confidence tier (the ``filter`` column) and by
dilution, plus 5'->3' orientation accuracy. Pure stdlib; no third-party deps.

NOTE: keep real cohort truth manifests (which carry sample/cell-line identifiers) OUT of this
public repo. This script is generic — point it at a local manifest.

Truth manifest (TSV, header row):
    sample_id  gene_a  gene_b  expected  [dilution]  [notes]
  expected ∈ {pass, detect, present, negative}
    pass     -> the matching call must have filter==PASS (headline recall)
    detect   -> any matching call (PASS or tagged) counts as detected
    present  -> biologically present; a miss is reported but not a headline FN
    negative -> no fusion expected; any PASS call on the sample is an FP
  Use gene_a=gene_b=NONE for negative rows.

Detection match is the UNORDERED, case-insensitive gene pair. Orientation correctness is
scored separately: for a matched call with partner_order==5to3, does (gene_a,gene_b) equal
the truth 5'->3' order? FP = filter==PASS calls matching no truth pair for that sample.

Usage:
    score_fusions.py --truth truth.tsv --results results/        # dir of *.fusions.tsv
    score_fusions.py --truth truth.tsv --results results/ --tsv-out scored.tsv
    score_fusions.py --self-test
"""

import argparse
import csv
import glob
import os
import sys
from collections import defaultdict

EXPECTED_LEVELS = {"pass", "detect", "present", "negative"}


def norm_gene(g):
    return (g or "").strip().upper()


def pair_key(a, b):
    """Unordered, case-insensitive gene-pair key."""
    return frozenset((norm_gene(a), norm_gene(b)))


class Truth:
    __slots__ = ("sample_id", "gene_a", "gene_b", "expected", "dilution", "notes", "key")

    def __init__(self, row):
        self.sample_id = row["sample_id"].strip()
        self.gene_a = norm_gene(row.get("gene_a"))
        self.gene_b = norm_gene(row.get("gene_b"))
        self.expected = row["expected"].strip().lower()
        self.dilution = (row.get("dilution") or "").strip()
        self.notes = (row.get("notes") or "").strip()
        self.key = pair_key(self.gene_a, self.gene_b)
        if self.expected not in EXPECTED_LEVELS:
            raise ValueError(
                f"{self.sample_id}: expected='{self.expected}' not in {sorted(EXPECTED_LEVELS)}"
            )


class Call:
    __slots__ = ("gene_a", "gene_b", "filt", "partner_order", "supporting", "key")

    def __init__(self, row):
        self.gene_a = norm_gene(row.get("gene_a"))
        self.gene_b = norm_gene(row.get("gene_b"))
        self.filt = (row.get("filter") or "").strip()
        self.partner_order = (row.get("partner_order") or "").strip()
        try:
            self.supporting = int(row.get("supporting_reads") or 0)
        except ValueError:
            self.supporting = 0
        self.key = pair_key(self.gene_a, self.gene_b)

    @property
    def is_pass(self):
        return self.filt == "PASS"


def load_truth(path):
    """sample_id -> list[Truth]. Errors on unknown expected levels / missing columns."""
    by_sample = defaultdict(list)
    with open(path, newline="") as f:
        rdr = csv.DictReader(f, delimiter="\t")
        required = {"sample_id", "gene_a", "gene_b", "expected"}
        missing = required - set(rdr.fieldnames or [])
        if missing:
            sys.exit(f"truth manifest is missing required columns: {sorted(missing)}")
        for row in rdr:
            if not (row.get("sample_id") or "").strip():
                continue
            t = Truth(row)
            by_sample[t.sample_id].append(t)
    return by_sample


def load_results(results_dir):
    """sample_id -> list[Call], read from every *.fusions.tsv under results_dir.

    Keying is by the sample_id column INSIDE each file (robust to filenames), so a
    per-sample subdirectory layout or a flat directory both work.
    """
    by_sample = defaultdict(list)
    files = glob.glob(os.path.join(results_dir, "**", "*.fusions.tsv"), recursive=True)
    if not files:
        sys.exit(f"no *.fusions.tsv found under {results_dir}")
    for fp in files:
        with open(fp, newline="") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                sid = (row.get("sample_id") or "").strip()
                if sid:
                    by_sample[sid].append(Call(row))
    return by_sample, len(files)


def score(truth_by_sample, calls_by_sample):
    """Return (per_sample_rows, rollup). Pure; no I/O."""
    per_sample = []
    tp = defaultdict(int)  # expected level -> matched
    fn = defaultdict(int)  # expected level -> missed
    fp_total = 0
    orient_ok = orient_total = 0
    by_dil = defaultdict(lambda: [0, 0])  # dilution -> [tp_pass, fn_pass]

    all_samples = sorted(set(truth_by_sample) | set(calls_by_sample))
    for sid in all_samples:
        truths = truth_by_sample.get(sid, [])
        calls = calls_by_sample.get(sid, [])
        if not truths:
            per_sample.append(
                {"sample_id": sid, "status": "UNSCORED (no truth row)",
                 "expected": "", "fusion": "", "called_filter": "",
                 "partner_order": "", "n_pass_calls": sum(c.is_pass for c in calls)}
            )
            continue

        truth_keys = {t.key for t in truths if t.expected != "negative"}

        for t in truths:
            if t.expected == "negative":
                continue
            matches = [c for c in calls if c.key == t.key]
            pass_match = next((c for c in matches if c.is_pass), None)
            any_match = matches[0] if matches else None

            if t.expected == "pass":
                hit = pass_match is not None
            else:  # detect / present
                hit = any_match is not None

            status = "TP" if hit else ("FN" if t.expected in ("pass", "detect") else "MISS(present)")
            if hit:
                tp[t.expected] += 1
            elif t.expected in ("pass", "detect"):
                fn[t.expected] += 1

            if t.expected == "pass" and t.dilution:
                by_dil[t.dilution][0 if hit else 1] += 1

            chosen = pass_match or any_match
            po = chosen.partner_order if chosen else ""
            if chosen and chosen.partner_order == "5to3":
                orient_total += 1
                if (chosen.gene_a, chosen.gene_b) == (t.gene_a, t.gene_b):
                    orient_ok += 1

            per_sample.append(
                {"sample_id": sid, "status": status, "expected": t.expected,
                 "fusion": f"{t.gene_a}::{t.gene_b}",
                 "called_filter": (chosen.filt if chosen else "—"),
                 "partner_order": po,
                 "n_pass_calls": sum(c.is_pass for c in calls)}
            )

        # False positives: PASS calls whose pair is in no (non-negative) truth for this sample.
        for c in calls:
            if c.is_pass and c.key not in truth_keys:
                fp_total += 1
                per_sample.append(
                    {"sample_id": sid, "status": "FP", "expected": "",
                     "fusion": f"{c.gene_a}::{c.gene_b}", "called_filter": c.filt,
                     "partner_order": c.partner_order, "n_pass_calls": ""}
                )

    rollup = {
        "tp_pass": tp["pass"], "fn_pass": fn["pass"],
        "tp_detect": tp["detect"], "fn_detect": fn["detect"],
        "fp_pass": fp_total,
        "orient_ok": orient_ok, "orient_total": orient_total,
        "by_dilution": dict(by_dil),
    }
    return per_sample, rollup


def _num(s):
    try:
        return float(s)
    except ValueError:
        return float("inf")


def fmt_rollup(rollup):
    def rate(n, d):
        return f"{n}/{d} = {100*n/d:.1f}%" if d else f"{n}/0 = n/a"

    tp_p, fn_p = rollup["tp_pass"], rollup["fn_pass"]
    tp_d, fn_d = rollup["tp_detect"], rollup["fn_detect"]
    fp = rollup["fp_pass"]
    lines = ["=" * 60, "COHORT ROLLUP", "=" * 60]
    lines.append(f"PASS recall   (expected=pass):   {rate(tp_p, tp_p + fn_p)}")
    lines.append(f"Detect recall (expected=detect): {rate(tp_d, tp_d + fn_d)}")
    lines.append(f"PASS precision:                  {rate(tp_p, tp_p + fp)}   (FP={fp})")
    if rollup["orient_total"]:
        lines.append(
            f"Orientation accuracy (5to3 calls): {rate(rollup['orient_ok'], rollup['orient_total'])}"
        )
    if rollup["by_dilution"]:
        lines.append("\nPASS recall by dilution:")
        for dil in sorted(rollup["by_dilution"], key=_num):
            hit, miss = rollup["by_dilution"][dil]
            lines.append(f"  {dil:>8}: {rate(hit, hit + miss)}")
    return "\n".join(lines)


def write_tsv(per_sample, path):
    cols = ["sample_id", "status", "expected", "fusion", "called_filter",
            "partner_order", "n_pass_calls"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in per_sample:
            w.writerow(r)


def run(truth_path, results_dir, tsv_out):
    truth = load_truth(truth_path)
    calls, n_files = load_results(results_dir)
    per_sample, rollup = score(truth, calls)
    print(f"loaded {sum(len(v) for v in truth.values())} truth rows over {len(truth)} samples; "
          f"{n_files} result files over {len(calls)} samples\n")
    missing = sorted(set(truth) - set(calls))
    unscored = sorted(set(calls) - set(truth))
    if missing:
        print(f"WARNING: {len(missing)} truth samples have no results: {missing}")
    if unscored:
        print(f"WARNING: {len(unscored)} result samples have no truth row: {unscored}")
    print()
    for r in per_sample:
        if r["status"] != "TP":
            print(f"  [{r['status']:>14}] {r['sample_id']:<14} {r['fusion']:<16} "
                  f"filter={r['called_filter']} {r['partner_order']}")
    print("\n" + fmt_rollup(rollup))
    if tsv_out:
        write_tsv(per_sample, tsv_out)
        print(f"\nper-call scored table -> {tsv_out}")


def self_test():
    """Synthetic check of the scoring logic (no files needed). Gene symbols are public; the
    sample_ids are synthetic placeholders."""
    def T(sid, a, b, exp, dil=""):
        return Truth({"sample_id": sid, "gene_a": a, "gene_b": b, "expected": exp, "dilution": dil})

    def C(a, b, filt, po, supp):
        return Call({"gene_a": a, "gene_b": b, "filter": filt, "partner_order": po, "supporting_reads": str(supp)})

    truth = {
        "S_pass_hit": [T("S_pass_hit", "EWSR1", "FLI1", "pass", "100")],
        "S_pass_miss": [T("S_pass_miss", "EWSR1", "ERG", "pass", "0.1")],
        "S_detect": [T("S_detect", "PAX3", "FOXO1", "detect")],
        "S_orient_ok": [T("S_orient_ok", "EWSR1", "FLI1", "pass", "100")],
        "S_neg": [T("S_neg", "NONE", "NONE", "negative")],
    }
    calls = {
        # reversed order + asserted 5to3 -> TP detection but orientation MISS
        "S_pass_hit": [C("FLI1", "EWSR1", "PASS", "5to3", 50)],
        # only a tagged call -> FN for expected=pass
        "S_pass_miss": [C("EWSR1", "ERG", "chimera", "index", 3)],
        # tagged match counts for expected=detect
        "S_detect": [C("PAX3", "FOXO1", "samelocus", "index", 5)],
        # correct 5to3 order -> orientation hit
        "S_orient_ok": [C("EWSR1", "FLI1", "PASS", "5to3", 80)],
        # a PASS call on a negative control -> FP
        "S_neg": [C("BCR", "ABL1", "PASS", "index", 9)],
    }
    _, r = score(truth, calls)
    assert r["tp_pass"] == 2 and r["fn_pass"] == 1, r
    assert r["tp_detect"] == 1 and r["fn_detect"] == 0, r
    assert r["fp_pass"] == 1, r
    assert r["orient_total"] == 2 and r["orient_ok"] == 1, r
    assert r["by_dilution"]["100"] == [2, 0] and r["by_dilution"]["0.1"] == [0, 1], r
    print("self-test OK")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", help="truth manifest TSV")
    ap.add_argument("--results", help="directory containing *.fusions.tsv (searched recursively)")
    ap.add_argument("--tsv-out", help="optional path to write the per-call scored table")
    ap.add_argument("--self-test", action="store_true", help="run the built-in logic check and exit")
    args = ap.parse_args()
    if args.self_test:
        self_test()
        return
    if not (args.truth and args.results):
        ap.error("--truth and --results are required (or use --self-test)")
    run(args.truth, args.results, args.tsv_out)


if __name__ == "__main__":
    main()

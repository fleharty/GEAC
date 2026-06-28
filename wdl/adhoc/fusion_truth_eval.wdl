version 1.0

## fusion_truth_eval.wdl  [ADHOC — disposable]
##
## Score fusion calls against per-sample truth and sweep the unique-anchor threshold, producing a
## cohort sensitivity / precision curve. Pairs with fusion_unique_anchor_eval.wdl: that re-emits the
## per-sample calls carrying `n_unique_anchored`; this consumes those TSVs plus the truth label.
##
## Truth flows in as a Terra column: scatter a fusion_sample_set and Terra populates the parallel
## arrays from `this.fusion_samples.fusions_tsv`, `this.fusion_samples.truth_fusion`, etc. The arrays
## MUST be index-aligned (they are when sourced from one set). truth_fusion is "GENE_A::GENE_B" for a
## positive, or "negative" for a control.
##
## Pure-stdlib Python scoring (no duckdb/pandas) on the small fusions TSVs — runs on a plain python
## image, independent of the geac image. NOT a supported pipeline; see wdl/adhoc/README.md.
##
## At each threshold N a call is "kept" when filter == "PASS" and n_unique_anchored >= N. For a
## positive sample, detection = a kept call whose gene pair equals the truth pair (order-insensitive);
## every other kept call is a false positive. For a negative sample, every kept call is a false
## positive. Outputs cohort_sweep.tsv (per-N sensitivity/precision) and cohort_detail.tsv (per sample).

workflow FusionTruthEval {
    input {
        # Parallel, index-aligned arrays (Terra fills these from a fusion_sample_set).
        Array[File]   fusions_tsvs       # per-sample fusions TSV carrying n_unique_anchored
        Array[String] sample_ids         # this.fusion_samples.fusion_sample_id
        Array[String] truth_fusions      # this.fusion_samples.truth_fusion
        Array[String] dilutions          # this.fusion_samples.truth_dilution (free text; "" ok)

        Array[Int]    thresholds   = [0, 1, 2, 5, 10, 15, 20, 25, 50, 100]
        String        python_docker = "python:3.11-slim"
    }

    scatter (i in range(length(fusions_tsvs))) {
        call ScorePerSample {
            input:
                fusions_tsv  = fusions_tsvs[i],
                sample_id    = sample_ids[i],
                truth_fusion = truth_fusions[i],
                dilution     = dilutions[i],
                thresholds   = thresholds,
                python_docker = python_docker
        }
    }

    call AggregateSweep {
        input:
            per_sample_sweeps = ScorePerSample.sweep_tsv,
            python_docker     = python_docker
    }

    output {
        File cohort_sweep  = AggregateSweep.cohort_sweep   # threshold -> sensitivity / precision
        File cohort_detail = AggregateSweep.cohort_detail  # per sample x threshold
    }
}

task ScorePerSample {
    input {
        File          fusions_tsv
        String        sample_id
        String        truth_fusion
        String        dilution
        Array[Int]    thresholds
        String        python_docker
    }

    command <<<
        set -euo pipefail
        cat > score.py <<'EOF'
import csv, sys
tsv, sample_id, truth_fusion, dilution, thr_csv, out = sys.argv[1:7]
thresholds = [int(x) for x in thr_csv.split(",")]
truth = None if truth_fusion == "negative" else frozenset(truth_fusion.split("::"))
rows = list(csv.DictReader(open(tsv), delimiter="\t"))

def nua(r):
    v = r.get("n_unique_anchored", "NA")
    try:
        return int(v)
    except (TypeError, ValueError):
        return 0  # NA / untracked -> treated as 0 unique anchors

with open(out, "w") as o:
    o.write("sample_id\ttruth_fusion\tdilution\tthreshold\tis_positive\tdetected\tn_fp\tn_pass\n")
    for N in thresholds:
        passing = [r for r in rows if r.get("filter") == "PASS" and nua(r) >= N]
        detected, n_fp = 0, 0
        for r in passing:
            pair = frozenset([r.get("gene_a"), r.get("gene_b")])
            if truth is not None and pair == truth:
                detected = 1
            else:
                n_fp += 1
        is_pos = 1 if truth is not None else 0
        o.write(f"{sample_id}\t{truth_fusion}\t{dilution}\t{N}\t{is_pos}\t{detected}\t{n_fp}\t{len(passing)}\n")
EOF
        python3 score.py "~{fusions_tsv}" "~{sample_id}" "~{truth_fusion}" "~{dilution}" "~{sep=',' thresholds}" "~{sample_id}.sweep.tsv"
    >>>

    output {
        File sweep_tsv = "~{sample_id}.sweep.tsv"
    }

    runtime {
        docker:      python_docker
        memory:      "2 GB"
        cpu:         1
        disks:       "local-disk 10 HDD"
        preemptible: 2
    }
}

task AggregateSweep {
    input {
        Array[File] per_sample_sweeps
        String      python_docker
    }

    command <<<
        set -euo pipefail
        cat > agg.py <<'EOF'
import csv, sys
from collections import defaultdict
files = sys.argv[1:]
rows = []
for f in files:
    rows += list(csv.DictReader(open(f), delimiter="\t"))

with open("cohort_detail.tsv", "w") as o:
    w = csv.DictWriter(o, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(rows)

by_n = defaultdict(list)
for r in rows:
    by_n[int(r["threshold"])].append(r)

with open("cohort_sweep.tsv", "w") as o:
    o.write("threshold\tn_positive\tn_detected\tsensitivity\tn_negative\t"
            "samples_with_fp\ttotal_fp\tprecision\n")
    for N in sorted(by_n):
        rs = by_n[N]
        pos = [r for r in rs if r["is_positive"] == "1"]
        neg = [r for r in rs if r["is_positive"] == "0"]
        n_det = sum(int(r["detected"]) for r in pos)
        total_fp = sum(int(r["n_fp"]) for r in rs)
        samples_with_fp = sum(1 for r in rs if int(r["n_fp"]) > 0)
        sens = n_det / len(pos) if pos else 0.0
        prec = n_det / (n_det + total_fp) if (n_det + total_fp) > 0 else 0.0
        o.write(f"{N}\t{len(pos)}\t{n_det}\t{sens:.4f}\t{len(neg)}\t"
                f"{samples_with_fp}\t{total_fp}\t{prec:.4f}\n")
EOF
        python3 agg.py ~{sep=' ' per_sample_sweeps}
    >>>

    output {
        File cohort_sweep  = "cohort_sweep.tsv"
        File cohort_detail = "cohort_detail.tsv"
    }

    runtime {
        docker:      python_docker
        memory:      "2 GB"
        cpu:         1
        disks:       "local-disk 10 HDD"
        preemptible: 2
    }
}

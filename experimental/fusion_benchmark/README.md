# Fusion benchmarking harness

Score `geac experimental fusions` calls against a known-truth manifest, to make "did this
change help or hurt?" a one-command answer instead of an eyeball over per-sample TSVs.

This is **phase 1** (local scoring) of the harness. Planned phases: (2) a Terra WDL that runs
the cohort and scores as a workflow step; (3) a GitHub Action that triggers it and posts the
rollup as a CI summary. See the Benchmarking-harness item in [TODO.md](../../TODO.md).

> **Containment:** keep real truth manifests (which carry sample / cell-line identifiers) OUT
> of this public repo. `score_fusions.py` is generic — point it at a local manifest.

## Usage

```bash
# score a directory of per-sample *.fusions.tsv against a truth manifest
python3 experimental/fusion_benchmark/score_fusions.py --truth truth.tsv --results results/ --tsv-out scored.tsv

# logic self-check (synthetic; no files needed)
python3 experimental/fusion_benchmark/score_fusions.py --self-test
```

`--results` is searched recursively for `*.fusions.tsv`; scoring keys on the `sample_id`
column *inside* each file, so a flat directory or per-sample subdirectories both work.

## Truth manifest (TSV, header row)

| column | required | meaning |
|---|---|---|
| `sample_id` | yes | must match the `sample_id` value inside the per-sample `fusions.tsv` |
| `gene_a` | yes | expected **5′** partner (donor); `NONE` for a negative control |
| `gene_b` | yes | expected **3′** partner; `NONE` for a negative control |
| `expected` | yes | `pass` \| `detect` \| `present` \| `negative` (below) |
| `dilution` | no | numeric spike level / VAF, for the recall-by-dilution breakdown |
| `notes` | no | free text (cell line, known complication) |

`expected` levels:
- **`pass`** — the matching call must reach `filter==PASS` (drives the headline PASS recall).
- **`detect`** — any matching call (PASS *or* tagged) counts as detected.
- **`present`** — biologically present but not required to be seen (e.g. trace dilution); a
  miss is reported, not counted as a headline FN.
- **`negative`** — no fusion expected; any `PASS` call on the sample is a false positive.

One row per expected fusion (a sample with two known fusions gets two rows).

## Scoring rules

- **Detection match** is the unordered, case-insensitive gene pair — truth `EWSR1::FLI1`
  matches a called `FLI1::EWSR1`.
- **Orientation accuracy** is scored separately: for a matched call that asserted an order
  (`partner_order==5to3`), does `(gene_a, gene_b)` equal the truth 5′→3′ order? (Calls that
  abstain — `partner_order==index` — are not counted for or against orientation.)
- **False positives** are `filter==PASS` calls whose pair matches no truth fusion for that
  sample. Tagged calls (`chimera`/`samelocus`/`pon`/…) are reported as context, not FPs.

## Output

A per-sample list of the actionable rows (FN / FP / present-misses) followed by a cohort
rollup: PASS recall, detect recall, PASS precision, orientation accuracy, and PASS recall
stratified by dilution (the limit-of-detection curve). `--tsv-out` writes the full per-call
scored table.

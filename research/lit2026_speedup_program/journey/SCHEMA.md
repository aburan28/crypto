# Journey entry schema (`lit2026-journey-entry/v1`)

One JSON file per measured step: `journey/entries/<YYYYMMDD>-<workstream>-<slug>.json`.
Each step adds a new file. Never edit another step's file; a correction is a new entry
whose `supersedes` names the old one.

All values are copied from the step's run report. The renderer
(`render_journey.py`) plots these recorded values and computes nothing.

```json
{
  "schema": "lit2026-journey-entry/v1",
  "workstream": "W4",
  "slug": "f4-early-abort",
  "date": "2026-10-12",
  "pr": "aburan28/crypto#NNNN",
  "papers": ["arXiv 2603.24890"],
  "kind": "speedup | no-change | regression | diagnostic | negative",
  "workload": {
    "id": "gf2-f4-descended-s3-n19",
    "description": "one-line description, including the exact curve id and subgroup",
    "curve_ids": ["..."]
  },
  "metric": "median wall time per decomposition attempt",
  "unit": "ms | Ir | S | counted ops | rows",
  "isolation": "L0 | L1 | L2",
  "samples": 24,
  "arms": [
    {"role": "original", "label": "program baseline", "commit": "4b6e09b78", "value": 0.0},
    {"role": "previous", "label": "best before this step", "commit": "...", "value": 0.0},
    {"role": "candidate", "label": "this step", "commit": "...", "value": 0.0}
  ],
  "ratio_original_over_candidate": {"value": 0.0, "interval": [0.0, 0.0], "interval_kind": "paired bootstrap 95% | min-max of rounds | none"},
  "ratio_previous_over_candidate": {"value": 0.0, "interval": [0.0, 0.0], "interval_kind": "..."},
  "outputs_identical": true,
  "run_report": "research/<study>/reports/run-<UTC>/report.pdf",
  "note": "scope and limits in one or two sentences"
}
```

Rules:

- **`original` arm.** It is the program baseline binary, built at commit `4b6e09b78`
  (main when the program started) and re-measured in the same interleaved session as
  the candidate. Ratios are never chained across sessions. If the original commit
  cannot run the workload, because the feature did not exist, set the `original`
  arm's value to `null` and use `previous`.
- **Ratios.** Record exactly what the harness or run report computed, with its
  interval.
- **Non-speedup kinds.** A `diagnostic`, `negative` or `no-change` step still gets an
  entry. Its arms may be empty, and the `note` carries the result.
- **Wall time.** Admit wall time only at the isolation level the harness grants.
  Deterministic counters (rows, Ir, counted ops) are preferred.

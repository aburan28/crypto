#!/usr/bin/env python3
import json

import validate as v

CONFIGS = {
    "batch24-cache6-min3": 24,
    "batch32-cache8-min3": 32,
}

reference_builtin = v.OUT / "selected-builtin-test.log"
builtin_rows = []
for label in CONFIGS:
    binary = v.ROOT / "build" / f"g7-{label}"
    code, elapsed, _ = v.run([binary, "--test"],
                             v.OUT / f"{label}-builtin-test.log", 600)
    row = {"label": label, "returncode": code, "elapsed_seconds": elapsed,
           "log_sha256": v.sha(v.OUT / f"{label}-builtin-test.log")}
    builtin_rows.append(row)
assert all(row["returncode"] == 0 for row in builtin_rows)
assert all(row["log_sha256"] == v.sha(reference_builtin) for row in builtin_rows)

state_rows = []
for runid in (0, 139):
    reference = v.client(v.SELECTED, 16, "followup-selected", runid)
    for label, batch in CONFIGS.items():
        candidate = v.client(v.ROOT / "build" / f"g7-{label}", batch,
                             label, runid)
        assert candidate["checkpoint_sha256"] == reference["checkpoint_sha256"]
        assert candidate["dp_multiset_sha256"] == reference["dp_multiset_sha256"]
        assert candidate["dp_records"] == reference["dp_records"]
        state_rows.append({"label": label, "run_id": runid,
                           "reference": reference, "candidate": candidate})
result = {"builtin_tests": builtin_rows, "state_checks": state_rows}
(v.OUT / "followup-validation.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({"builtins": len(builtin_rows),
                  "state_checks": len(state_rows)}, indent=2))

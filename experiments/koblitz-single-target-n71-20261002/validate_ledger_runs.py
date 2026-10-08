#!/usr/bin/env python3
"""Independent Python arithmetic replay for the n=71 ledger runs (R1-R3)."""
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
REPLAY = REPO / "research/sat_factor_base_review_20260908/autolab_orbit_extract_20260924/independent_replay.py"
spec = importlib.util.spec_from_file_location("independent_replay", REPLAY)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

header = json.loads((HERE / "base_n71_K600.jsonl").read_text().splitlines()[0])
curve = module.Curve(module.Field(header["n"], header["field_modulus_low_terms"]), header["a"])
r = header["subgroup_order"]
x_codes = set()
points_by_x = {}
for point in header["factor_base_point_coordinates"]:
    if point is not None:
        x, y = (int(v) for v in point)
        points_by_x.setdefault(x, []).append((x, y))
        x_codes.add(x)
answer = int((HERE / "ledger-freeze" / "known-answer.txt").read_text().strip())
results = []
for run in ("ledger-R1", "ledger-R2", "ledger-R3"):
    records = [json.loads(line) for line in (HERE / "runs" / run / "ic.jsonl").read_text().splitlines() if line.strip()]
    assert len(records) == 1, (run, len(records))
    record = records[0]
    assert [int(v) for v in record["target"]] == [int(v) for v in record["published_q"]]
    record["published_fixture_scalar"] = answer
    # Normalize the decimal-string wire format to ints for replay.
    record["generator"] = [int(v) for v in record["generator"]]
    record["target"] = [int(v) for v in record["target"]]
    record["published_q"] = [int(v) for v in record["published_q"]]
    record["x_codes"] = [int(v) for v in record["x_codes"]]
    record["pinned_intermediates"] = [int(v) for v in record["pinned_intermediates"]]
    checks = module.replay_record(record, header, curve, r, x_codes, points_by_x)
    rho_rows = [json.loads(line) for line in (HERE / "runs" / run / "rho.jsonl").read_text().splitlines() if line.strip()]
    rho = next(x for x in rho_rows if x.get("kind") == "rho_public_fixture")
    checks["rho_same_public_point"] = [int(v) for v in rho["published_q"]] == [int(v) for v in record["published_q"]]
    checks["rho_scalar_matches_sidecar"] = int(rho["recovered_fixture_scalar"]) == answer
    checks["rho_verified"] = rho["verified"] is True
    results.append({"run_id": run, "checks": checks, "pass": all(checks.values())})
assert all(result["pass"] for result in results), results
report = {
    "method": "Python GF(2^n) and affine-curve replay independent of the Rust producers",
    "runs": results,
    "all_pass": True,
}
(HERE / "independent-replay-ledger.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))

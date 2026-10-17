"""Audit source binding and paired stage accounting for cubic-cut variants."""

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


path = HERE / "type_i_fast_cubic_probe.json"
receipt = json.loads(path.read_text())
assert receipt["status"] == "matched_prime_field_cubic_reduction_probe"
assert receipt["source_sha256"] == sha(HERE / "line_quotient_cubic_fast_probe.py")
assert receipt["frozen_sha256"] == sha(HERE / "type_i_two_pivot_reduced64.json")
assert receipt["baseline_source_sha256"] == sha(HERE / "line_quotient_cubic_probe.py")
assert receipt["plane_sha256"] == sha(HERE / "type_i_plane_quartic.json")
runtime = HERE / "type_i_fast_cubic_runtime_info.json"
assert receipt["runtime_info_sha256"] == sha(runtime)
assert json.loads(runtime.read_text())["status"] == "verified"
assert receipt["classes"] == len(receipt["pairs"]) == 16
assert receipt["residual_cuts"] == 32
log = (HERE / "type_i_fast_cubic_probe.log").read_text()
assert [int(line.split()[1]) for line in log.splitlines()
        if line.startswith("FAST_CUBIC_PAIR_PASS ")] == list(range(16))
assert "FAST_CUBIC_PASS 16 " in log
totals = {method: 0 for method in ("fast", "hybrid", "original")}
for index, pair in enumerate(receipt["pairs"]):
    assert pair["index"] == index
    assert set(pair["order"]) == set(totals)
    for method in totals:
        for side in ("positive_input_residual", "negative_input_residual"):
            elapsed = pair[f"{side}_{method}_ns"]
            assert isinstance(elapsed, int) and elapsed > 0
            totals[method] += elapsed
for method, elapsed in totals.items():
    assert receipt[f"{method}_total_ns_unisolated"] == elapsed
output = {
    "status": "verified_matched_cubic_probe_receipt",
    "classes": 16, "residual_cuts": 32,
    "totals_ns_unisolated": totals,
    "producer_sha256": sha(path), "source_sha256": sha(Path(__file__)),
    "boundary": "Source, input, runtime, progression, rotated order, and exclusive stage totals audited. The Sage producer checks exact equality to frozen cubic and residual coefficients.",
}
(HERE / "type_i_fast_cubic_probe_verification.json").write_text(
    json.dumps(output, indent=2, sort_keys=True) + "\n")
print("FAST_CUBIC_VERIFICATION_PASS", output["residual_cuts"])

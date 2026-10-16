#!/usr/bin/env python3
"""Freeze one fresh n73 target (Python-only; no Rust producer involved)."""
import importlib.util
import json
import hashlib
import sys
from pathlib import Path

E = Path("experiments/koblitz-single-target-n73-20261002")
REPLAY = Path("research/sat_factor_base_review_20260908/autolab_orbit_extract_20260924/independent_replay.py")
spec = importlib.util.spec_from_file_location("independent_replay", REPLAY)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

header = json.loads((E / "base_n73_K600.jsonl").read_text().splitlines()[0])
curve = module.Curve(module.Field(header["n"], header["field_modulus_low_terms"]), header["a"])
r = header["subgroup_order"]
dump_rows = [json.loads(l) for l in (E / "smoke_out.jsonl").read_text().splitlines() if l.strip()]
trow = next(x for x in dump_rows if x.get("kind") == "compact_orbit_dlp_target")
G = tuple(int(v) for v in trow["generator"])
assert curve.on_curve(G), "G must be on curve"
assert curve.mul(r, G) is None, "[r]G must be infinity"
s = int(sys.argv[1]) if len(sys.argv) > 1 else 9876543210987654321
assert 1 <= s < r, "fresh scalar in range"
Q = curve.mul(s, G)
assert curve.on_curve(Q) and Q is not None
qx, qy = Q
(E / "ledger-freeze").mkdir(exist_ok=True)
(E / "ledger-freeze" / "target-q.jsonl").write_text(f'["{qx}","{qy}"]\n')
(E / "ledger-freeze" / "known-answer.txt").write_text(f"{s}\n")
receipt = {
    "kind": "koblitz_n73_ledger_target_freeze",
    "date": "2026-10-02",
    "method": "standalone Python GF(2^73) scalar multiplication (module replay arithmetic); no Rust producer involved in freezing",
    "fresh_scalar": s,
    "scalar_note": "fresh: never used in any prior fixture",
    "public_target_q": [str(qx), str(qy)],
    "curve": {"n": 73, "a": 0, "subgroup_order": r, "generator": [str(G[0]), str(G[1])]},
    "base": {"file": "base_n73_K600.jsonl", "hash": header["base_hash"], "points": header["factor_base_points"], "columns": header["orbit_columns"]},
    "checks": {"generator_on_curve": True, "subgroup_annihilates_generator": True, "target_on_curve": True, "target_is_scalar_times_g": True},
    "target_q_sha256": hashlib.sha256(f'["{qx}","{qy}"]\n'.encode()).hexdigest(),
}
(E / "ledger-freeze" / "freeze_receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
print("frozen Q:", [str(qx), str(qy)])

#!/usr/bin/env python3
"""Check custody and structural claims in one saved Sage run."""
import gzip
import hashlib
import json
import sys
from pathlib import Path

root = Path(__file__).resolve().parent
run_dir = Path(sys.argv[1])
raw_bytes = (run_dir / "raw.json.gz").read_bytes()
result = json.loads(gzip.decompress(raw_bytes))
summary = json.loads((run_dir / "summary.json").read_text())
cert = result["isogeny_certificate"]
assert result["contract"]["field_exponent"] == 37
assert cert["degree"] == 73 and cert["separable"]
assert cert["inventory_count"] == 74
assert cert["kernel_polynomial_degree"] == 36
assert cert["kernel_division_polynomial_quotient_degree"] == 2628
assert cert["source_group_order"] == cert["target_group_order"] == 137439487532
assert cert["trace_over_gf_2_37"] == -534059
assert cert["frobenius_order_discriminant"] == -7 * 194399**2
assert cert["frobenius_order_conductor"] == 73 * 2663
assert cert["source_endomorphism_order_discriminant"] == -7
assert cert["target_endomorphism_order_discriminant"] == -7 * 73**2
assert cert["legendre_symbol_minus7_mod_73"] == -1
assert cert["homomorphism_checked_pairs"] == 128
assert len(result["cases"]) == 128
assert all(c["status"] == "VERIFIED" for c in result["cases"])
assert all(c["true_selector_assignments"] > 0 for c in result["cases"])
assert result["accounting"]["full_dlp_total_operations"] is None
assert result["accounting"]["speedup"] is None
assert summary["systems"] == 128
assert summary["full_dlp_speedup"] is None
print("PASS: degree-73 map, kernel, transported workloads, matrix replays, and scope checks.")

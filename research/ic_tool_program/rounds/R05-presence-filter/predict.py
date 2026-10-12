#!/usr/bin/env python3
"""R05's model prediction, before any R05 run (PROTOCOL.md, "Prediction").

    python3 predict.py > prediction.json

From frozen files only: R04's stage costs and filter check, and v0's
profile from R01. For each suite size it gives the presence filter's
pass rate for an absent key, as built now and as the candidate builds
it, and what that does to collection and to cold time if nothing else
changes.

The pass rate is the blocked filter's: a key's bits all lie in one
64-bit word, and the word holds a Poisson number of the stored keys,
`λ = 64 · stored pairs / filter bits` on average. With `k` bits a key,

    pass(k) = Σ_j Poisson(j; λ) · (1 − (63/64)^(k·j))^k.

With `k = 1` and four bits a key it reproduces R04's measured admitted
fraction to within 2% from `2^36.6` up (`filter_check.json`).

The model then scales R04's admitted stage by the admitted keys a
summand, true hits kept, and leaves every other stage as R04 measured
it. It leaves out the larger filter's own cost, so it is an upper
estimate of the gain where the filter's misses matter.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
R04 = HERE.parent / "R04-scan-probes"
R01 = HERE.parent / "R01-baseline-v0"


def filter_bits(keys: int, per_key: int) -> int:
    """`filter_bits_for` in `koblitz_index_calculus.rs`, for `per_key`
    bits a key: the bit length of `keys · per_key`, within 6..32."""
    return min(max((max(keys, 1) * per_key).bit_length(), 6), 32)


def pass_rate(keys: int, bits: int, probes: int) -> float:
    lam = 64 * keys / 2 ** bits
    total, p, j = 0.0, math.exp(-lam), 0
    while j < 64 * 8 or p > 1e-18:
        if j:
            p *= lam / j
        total += p * (1 - (63 / 64) ** (probes * j)) ** probes
        j += 1
    return total


def main() -> None:
    stages = {s["slug"]: s for s in json.loads((R04 / "analysis.json").read_text())["sizes"]}
    check = {s["slug"]: s for s in json.loads((R04 / "filter_check.json").read_text())["sizes"]}
    profile = {(p["a"], p["n"]): p for p in json.loads((R01 / "analysis.json").read_text())["profile"]}
    rows = []
    for slug, c in check.items():
        s = stages[slug]
        p = profile[(s["a"], s["n"])]
        keys = c["stored_pairs"]
        old_bits, new_bits = filter_bits(keys, 4), filter_bits(keys, 8)
        old, new = pass_rate(keys, old_bits, 1), pass_rate(keys, new_bits, 3)
        true = c["true_hits_per_summand"]
        measured = c["measured_admitted_per_summand"]
        admitted_new = true + new * (1 - true)
        st = s["stages"]
        collect_ns = st["scan_ns_per_summand"] + st["trial"]["ns_per_summand"]
        saved = st["admitted"]["ns_per_summand"] * (1 - admitted_new / measured)
        collect_ratio = collect_ns / (collect_ns - saved)
        # Collection's share of cold time in v0's profile.
        share = p["setup_phase_shares"]["collect"] * p["setup_ms_median"] / p["cold_ms_median"]
        rows.append({
            "slug": slug, "log2_r": s["log2_r"], "stored_pairs": keys,
            "filter_bits_now": old_bits, "filter_bits_candidate": new_bits,
            "filter_mib_now": 2 ** old_bits / 8 / 2 ** 20, "filter_mib_candidate": 2 ** new_bits / 8 / 2 ** 20,
            "pass_rate_now": round(old, 4), "pass_rate_candidate": round(new, 4),
            "admitted_per_summand_measured": measured, "admitted_per_summand_candidate": round(admitted_new, 4),
            "collection_share_of_cold_v0": round(share, 3),
            "model_collection_ratio": round(collect_ratio, 3),
            "model_cold_ratio": round(1 / (1 - share * (1 - 1 / collect_ratio)), 3),
        })
    rows.sort(key=lambda r: r["log2_r"])
    print(json.dumps({
        "what_this_is": "R05's model prediction, made before any R05 run: the presence filter's pass rate "
                        "now and as the candidate builds it, and the collection and cold-time ratios that "
                        "follow if every other stage stays as R04 measured it",
        "sources": ["../R04-scan-probes/analysis.json", "../R04-scan-probes/filter_check.json",
                    "../R01-baseline-v0/analysis.json"],
        "sizes": rows,
    }, indent=1))


if __name__ == "__main__":
    main()

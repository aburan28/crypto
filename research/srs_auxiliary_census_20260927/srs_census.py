"""Cheon exposure of deployed powers-of-tau setups.

For each published setup: the largest exponent q with [tau^q]G1 published,
the largest divisor d of r-1 with d <= q (the p-1 variant needs [tau^d]G),
the largest divisor d' of r+1 with 2d' <= q (the p+1 variant needs 2d'
powers), and the resulting generic cost against Pollard rho on the same
group.  Factorisations reuse ../torsion_auxiliary_inputs_20260918/
curve_divisors.py unchanged, so every d is exact over the split part of
r -/+ 1 and a lower bound only where a cofactor is left unsplit.

Costs are in group operations, constant 1, both baby-step/giant-step
stages charged:
    p-1:  sqrt(r/d) + sqrt(d)
    p+1:  sqrt(r/d) + d
    rho:  sqrt(pi*r/4)            (negation map)
They are the generic floors the algorithms sit on, not measurements at
this size; the constant the parent round measured at 24-48 bits (~20x
the floor with fixed-base tables) is reported separately and is an
extrapolation.

    python3 srs_census.py            # writes results/srs_census.{md,json}
"""

from __future__ import annotations

import json
import math
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "torsion_auxiliary_inputs_20260918"))

from curve_divisors import CURVES, factor_bounded, largest_divisor_below  # noqa: E402

R = {"BLS12-381": CURVES["BLS12-381 (r)"], "BN254": CURVES["BN254 (r)"]}

# (setup, curve, q = largest published G1 exponent, source, published prior figure)
SETUPS = [
    ("EIP-4844 KZG (mainnet, n=12)", "BLS12-381", 2**12 - 1,
     "ethereum/kzg-ceremony-specs README + docs/participant/participant.md "
     "(4 sub-ceremonies, independent secrets; mainnet blobs use FIELD_ELEMENTS_PER_BLOB = 4096)",
     None),
    ("EIP-4844 KZG (largest, n=15)", "BLS12-381", 2**15 - 1,
     "ethereum/kzg-ceremony-specs (separate secret; not the mainnet transcript)",
     None),
    ("Zcash Sapling powers of tau", "BLS12-381", 2**22 - 2,
     "powersoftau src/lib.rs: TAU_POWERS_G1_LENGTH = (1<<21 << 1) - 1",
     {"d": 2**21, "cost": "2^117.2 exponentiations", "by": "ethresear.ch #6692 (daira)"}),
    ("Filecoin powers of tau", "BLS12-381", 2**28 - 2,
     "Filecoin 'Trusted setup' posts (64x Zcash, 2^27 tau powers); G1 = 2*2^27 - 1",
     {"d": 2**27, "cost": "~2^114", "by": "ethresear.ch #6692"}),
    ("Aztec Ignition", "BN254", 100_800_000,
     "AztecProtocol/ignition-verification Transcript_spec.md: x^1..x^100,800,000 in G1",
     {"d": 3 * 2**25, "cost": "~2^114", "by": "ethresear.ch #6692"}),
    ("Perpetual Powers of Tau", "BN254", 2**29 - 2,
     "privacy-ethereum/perpetualpowersoftau: 'up to 536870911 powers' = 2^29 - 1 points",
     {"d": 2**28, "cost": "~2^114", "by": "ethresear.ch #6692"}),
]

MEASURED_CONSTANT = 20.0  # parent round, p-1 case, fixed-base tables, 24-48 bits


def lg(x: float) -> float:
    return math.log2(x)


def main() -> None:
    rng = random.Random(1)
    fac = {}
    for name, r in R.items():
        fm, um = factor_bounded(r - 1, rng)
        fp, up = factor_bounded(r + 1, rng)
        fac[name] = {"r": r, "minus": fm, "minus_unsplit": um, "plus": fp, "plus_unsplit": up}

    rows = []
    for setup, curve, q, source, prior in SETUPS:
        f = fac[curve]
        r = f["r"]
        d1 = largest_divisor_below(f["minus"], q)
        d2 = largest_divisor_below(f["plus"], q // 2)
        c1 = math.sqrt(r / d1) + math.sqrt(d1)
        c2 = math.sqrt(r / d2) + d2
        rho = math.sqrt(math.pi * r / 4)
        best_side, best = ("p-1", c1) if c1 <= c2 else ("p+1", c2)
        row = {
            "setup": setup, "curve": curve, "q": q, "log2_q": round(lg(q), 3),
            "source": source,
            "d_pminus1": d1, "log2_d_pminus1": round(lg(d1), 3),
            "d_pplus1": d2, "log2_d_pplus1": round(lg(d2), 3),
            "log2_cost_pminus1": round(lg(c1), 2), "log2_cost_pplus1": round(lg(c2), 2),
            "best_side": best_side, "log2_cost_best": round(lg(best), 2),
            "log2_rho": round(lg(rho), 2),
            "bits_lost_vs_rho": round(lg(rho) - lg(best), 2),
            "log2_cost_at_measured_constant": round(lg(MEASURED_CONSTANT * best), 2),
            "prior": prior,
        }
        if prior:
            row["prior_d_divides"] = (f["r"] - 1) % prior["d"] == 0
            row["prior_d_leq_q"] = prior["d"] <= q
            row["d_gain_over_prior_bits"] = round(lg(d1) - lg(prior["d"]), 3)
        rows.append(row)

    out = HERE / "results"
    (out / "srs_census.json").write_text(json.dumps(
        {"factorisations": {k: {"r": str(v["r"]),
                                "r_minus_1": {str(p): e for p, e in v["minus"].items()},
                                "r_minus_1_unsplit": [str(u) for u in v["minus_unsplit"]],
                                "r_plus_1": {str(p): e for p, e in v["plus"].items()},
                                "r_plus_1_unsplit": [str(u) for u in v["plus_unsplit"]]}
                            for k, v in fac.items()},
         "rows": rows}, indent=2))

    L = []
    L.append("# Cheon exposure of deployed powers-of-tau setups\n")
    L.append("Generated by `srs_census.py`. Costs are generic floors in group operations "
             "(constant 1), not measurements at this size.\n")
    L.append("| setup | curve | top G1 exponent `q` | best `d` (p−1) | `log₂ d` | Cheon floor | rho | bits lost |")
    L.append("|:--|:--|--:|--:|--:|--:|--:|--:|")
    for w in rows:
        L.append(f"| {w['setup']} | {w['curve']} | {w['q']:,} | {w['d_pminus1']:,} | {w['log2_d_pminus1']:.2f} | "
                 f"`2^{w['log2_cost_best']:.2f}` | `2^{w['log2_rho']:.2f}` | **{w['bits_lost_vs_rho']:.2f}** |")
    L.append("\n## Cross-check against the published analysis (ethresear.ch #6692)\n")
    L.append("| setup | published `d` | divides `r−1` | `≤ q` | this census `d` | gain over published `d` |")
    L.append("|:--|--:|:--:|:--:|--:|--:|")
    for w in rows:
        if w["prior"]:
            L.append(f"| {w['setup']} | {w['prior']['d']:,} | {'yes' if w['prior_d_divides'] else '**no**'} | "
                     f"{'yes' if w['prior_d_leq_q'] else '**no**'} | {w['d_pminus1']:,} | "
                     f"{w['d_gain_over_prior_bits']:+.3f} bits of `d` ({w['d_gain_over_prior_bits']/2:+.3f} of security) |")
    L.append("\n## p+1 side (needs `2d` powers, `d | r+1`)\n")
    L.append("| setup | best `d` (p+1) | Cheon p+1 floor | vs p−1 floor |")
    L.append("|:--|--:|--:|--:|")
    for w in rows:
        L.append(f"| {w['setup']} | {w['d_pplus1']:,} | `2^{w['log2_cost_pplus1']:.2f}` | "
                 f"{w['log2_cost_pplus1'] - w['log2_cost_pminus1']:+.2f} bits |")
    L.append("\n## Factorisations used\n")
    for k, v in fac.items():
        fm = " · ".join(f"{p}" + (f"^{e}" if e > 1 else "") for p, e in sorted(v["minus"].items()))
        um = ", ".join(f"c{len(str(u))}" for u in v["minus_unsplit"]) or "—"
        L.append(f"- **{k}** `r−1` = {fm}; unsplit: {um}")
    (out / "srs_census.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()

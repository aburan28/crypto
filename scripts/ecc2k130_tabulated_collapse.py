#!/usr/bin/env python3
"""The tabulated rows of the ECC2K-130 decomposition note, with the orbit collapse.

`research/notes/ecc2k130/RESEARCH_ECC2K130_DECOMPOSITION.md` §6 establishes that
any materialised factor base can be a union of Frobenius orbits and take the
`n`-fold collapse.  Its tabulated rows -- every `split >= 2` cell, including
the memory-capped ones -- materialise the base by construction, yet are built
with `cost_cell(..., frobenius=False)`.  And `cost_cell`'s own collapse path
folds the relation count and the linear algebra but not the subset-sum table.

This script re-prices those rows with both folds, using the note's own
`cost_cell` unchanged, and sets the result beside two independent numbers for
the same method family: the tournament's family model (`docs/ic/README.md`,
"Round 3") and its measured IC/rho crossover carried to `r = 2^129`
(`research/ic_candidate_tournament_20260915/campaign_20260916/RESULTS.md`,
round 0022).  It is additive: the note's frozen artefact is not regenerated.

Everything here is derived.  Nothing is a run.

    python3 scripts/ecc2k130_tabulated_collapse.py
"""

import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ecc2k130_point_decomposition as note  # noqa: E402

N = 131
LOG2_N = math.log2(N)
LOG2_R = 129.0000               # r, the order of <G>, from the note's §0
LOG2_E = 131.0                  # #E = 4r
RHO = 60.8090                   # rho with <-1> x <pi>, the note's reference

# Round 0022 of the IC candidate tournament: IC/rho on Koblitz curves with a
# balanced Frobenius-orbit base, 64 fixtures a cell, both arms complete.
TOURNAMENT_CELLS = [
    ("n23a1", 4_196_903, 0.831),
    ("n37a0", 230_603_167, 1.396),
    ("n43a1", 4_644_189_029, 2.207),
]
TOURNAMENT_RATES = {"low": 0.075, "step1_low": 0.130, "derived": 1 / 6, "high": 0.188}


def cell(m, l, split, frobenius, fold_table):
    """The note's `cost_cell`, optionally folding the subset-sum table by `n`.

    Folding stores one representative per sigma-orbit of subset sums and
    canonicalises the lookup key; it is the tournament's `t = n`.  The stored
    base itself stays `2^l` points.
    """
    c = note.cost_cell(m, l, N, split, frobenius=frobenius)
    if fold_table and frobenius and c["log2_table_entries"] is not None:
        table = c["log2_table_entries"]
        if split > 1 and table > l:
            table = max(table - LOG2_N, l)
        terms = [c["log2_collection"], c["log2_linear_algebra"], table]
        c = dict(c, log2_table_entries=round(table, 2),
                 log2_total=round(note.log2_add(*terms), 2))
    return c


def best(m, split, frobenius, fold_table, cap=None):
    out = None
    for step in range(32, 4 * 70):
        c = cell(m, step / 4.0, split, frobenius, fold_table)
        if cap is not None and (c["log2_table_entries"] or 0) > cap:
            continue
        if out is None or c["log2_total"] < out["log2_total"]:
            out = c
    return out


def best_capped(cap, corrected):
    out = None
    for m in range(2, 9):
        for split in range(1, m + 1):
            c = best(m, split, corrected, corrected, cap)
            if c is not None and (out is None or c["log2_total"] < out["log2_total"]):
                out = dict(c, m=m, split=split)
    return out


def family_model(t, k):
    """`ops = 0.75 (#E t/k)^{2/3} / t` at `F* = (#E t/k)^{1/3}`, in log2."""
    inner = LOG2_E + math.log2(t) - math.log2(k)
    return {
        "t": t, "k": k,
        "log2_base": round(inner / 3, 2),
        "log2_table_entries": round(2 * inner / 3 - 2 - math.log2(t), 2),
        "log2_total": round(math.log2(0.75) + 2 * inner / 3 - math.log2(t), 2),
    }


def main():
    reconciliation = []
    for m in (2, 3):
        for label, frob, fold in [
            ("as committed (no collapse)", False, False),
            ("+ orbit collapse (the note's own frobenius path)", True, False),
            ("+ collapse + folded table", True, True),
        ]:
            c = best(m, 2, frob, fold)
            reconciliation.append({
                "m": m, "variant": label, "log2_base": c["l"],
                "log2_table_entries": c["log2_table_entries"],
                "log2_total": c["log2_total"],
                "log2_ratio_to_rho": round(c["log2_total"] - RHO, 2),
            })

    models = [family_model(1, 1), family_model(N, N)]
    for mdl in models:
        mdl["log2_ratio_to_rho"] = round(mdl["log2_total"] - RHO, 2)

    name, r, ratio = TOURNAMENT_CELLS[-1]
    anchor = math.log2(r)
    extrapolation = {
        "anchor_cell": name, "anchor_log2_r": round(anchor, 2),
        "anchor_log2_ratio": round(math.log2(ratio), 2),
        "log2_ratio_at_2^129": {
            key: round(math.log2(ratio) + rate * (LOG2_R - anchor), 2)
            for key, rate in TOURNAMENT_RATES.items()
        },
        "status": "extrapolation over 97 bits from three cells whose rate ranges "
                  "overlap; the lane itself reports no resolvable curvature",
    }

    memory = []
    for cap in (30, 40, 50, 60, 70):
        before = best_capped(cap, False)
        after = best_capped(cap, True)
        memory.append({
            "log2_entries_cap": cap,
            "as_committed": {"log2_total": before["log2_total"], "m": before["m"],
                             "split": before["split"]},
            "corrected": {"log2_total": after["log2_total"], "m": after["m"],
                          "split": after["split"], "log2_base": after["l"]},
            "bsgs_as_priced_by_the_note": round(max(cap, LOG2_R - cap), 2),
            "rho": RHO,
            "corrected_minus_rho": round(after["log2_total"] - RHO, 2),
        })

    # The tournament's own family -- m = 2, a folded pair table -- under a cap,
    # as against the best cell of any shape above.
    pair_family_capped = []
    for cap in (40, 50, 60, 70):
        c = best(2, 2, True, True, cap)
        pair_family_capped.append({
            "log2_entries_cap": cap, "log2_base": c["l"],
            "log2_total": c["log2_total"],
            "log2_ratio_to_rho": round(c["log2_total"] - RHO, 2),
        })

    report = {
        "curve": "K0 over F_2^131 (ECC2K-130)",
        "class": "accounting (AGENTS.md section 3): the note's own collapse, applied "
                 "where its own section 6 says it applies; no algorithm changed",
        "rho_log2": RHO,
        "reconciliation": reconciliation,
        "tournament_family_model": models,
        "tournament_measured_rate_extrapolated": extrapolation,
        "memory_capped": memory,
        "pair_family_m2_capped": pair_family_capped,
    }
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                       "experiments", "ecc2k130_tabulated_collapse.json")
    with open(out, "w") as fh:
        json.dump(report, fh, indent=2)
        fh.write("\n")

    for row in reconciliation:
        print(f"m={row['m']} {row['variant']:48s} 2^{row['log2_total']:6.2f}"
              f"  = 2^{row['log2_ratio_to_rho']:+6.2f} x rho")
    for mdl in models:
        print(f"family model t=k={mdl['t']:<3d}                              "
              f"       2^{mdl['log2_total']:6.2f}  = 2^{mdl['log2_ratio_to_rho']:+6.2f} x rho")
    for key, v in extrapolation["log2_ratio_at_2^129"].items():
        print(f"measured rate ({key}) carried to r = 2^129:                     "
              f"= 2^{v:+6.2f} x rho")
    for row in memory:
        print(f"cap 2^{row['log2_entries_cap']}: committed 2^{row['as_committed']['log2_total']:.2f}"
              f"  corrected 2^{row['corrected']['log2_total']:.2f}"
              f"  bsgs(as priced) 2^{row['bsgs_as_priced_by_the_note']:.2f}"
              f"  rho 2^{RHO:.2f}")
    for row in pair_family_capped:
        print(f"m=2 pair family, cap 2^{row['log2_entries_cap']}: 2^{row['log2_total']:.2f}"
              f" = 2^{row['log2_ratio_to_rho']:+.2f} x rho")


if __name__ == "__main__":
    main()

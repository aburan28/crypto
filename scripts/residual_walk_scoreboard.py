#!/usr/bin/env python3
"""Optimisation scoreboard for the residual-walk relation generators.

Reads the JSON written by `residual_walk_bench --baseline --json FILE`
(a list of strategy reports; the `recovery` list of a `--panel` JSON is
accepted too) and, for every (bits, B, dp, strategy) cell, decomposes the
cost of recovering d into the three factors the research note argues
about:

    S  = total_ops / sqrt(n)               the score (lower is better)
    kappa = samples / sqrt(n)              the collision *count* factor
    c  = walk_ops / samples                group operations per residual
    overheads = setup, replay, verify      as fractions of total_ops

together with the generic floors:

    kappa_floor = sqrt(2 (B + 1))   for any residual-collision method (birthday)
    kappa_rho   = sqrt(pi / 2)      for plain rho (first collision)
    S_floor     = kappa_floor * 1   one group operation per residual, no overhead

Usage:
    scripts/residual_walk_scoreboard.py RUN.json
    scripts/residual_walk_scoreboard.py RUN.json --baseline experiments/20_residual_walk_baseline.json
    scripts/residual_walk_scoreboard.py RUN.json --baseline BASE.json --fail-on-regression [--tolerance 0.10]

With --baseline, every cell prints the improvement factor
baseline_S / run_S (>1 is better) and the kappa ratio; a kappa ratio
below 1 - tolerance on a strategy would mean the *count* invariant was
beaten, which is the research target.  --fail-on-regression exits 1 if
any cell's S worsens by more than the tolerance.
"""
import argparse
import json
import math
import sys
from collections import defaultdict

TAG_ORDER = ["A", "B", "C1", "C2", "R"]


def load(path):
    with open(path) as f:
        data = json.load(f)
    if isinstance(data, dict):
        data = data["recovery"]
    return data


def mean(xs):
    xs = list(xs)
    return sum(xs) / len(xs) if xs else float("nan")


def variant(r):
    """Which generic levers a report ran with (absent fields = off)."""
    parts = []
    tag = r.get("tag")
    if r.get("negation_map") and r.get("dp_bits", 0) == 0:
        parts.append("neg")
    if r.get("diff_table") and tag == "B":
        parts.append("diff")
    if r.get("segment_len") and tag in ("C1", "R"):
        parts.append(f"seg{r['segment_len']}")
    if r.get("continue_after_collision") and tag == "B":
        parts.append("cont")
    if r.get("fold_order", 1) == 6:
        parts.append("aut6")
    if r.get("seed_pairs") or r.get("seeded_points"):
        parts.append("seed2")
    if r.get("j_zero"):
        parts.append("j0")
    if r.get("s3_oracle"):
        parts.append("s3")
    if r.get("s4_oracle"):
        parts.append("s4")
    if r.get("mitm_neighbours"):
        parts.append("mitm3")
    return "+".join(parts) or "plain"


def cells(rows):
    groups = defaultdict(list)
    for r in rows:
        groups[(r["bits"], r["factor_base"], r["dp_bits"], r["tag"], variant(r))].append(r)
    out = {}
    for key, g in groups.items():
        bits, fb, dp, tag, var = key
        sq = [math.sqrt(r["n"]) for r in g]
        n_ops = [r["total_ops"] for r in g]
        kappa_floor = math.sqrt(math.pi / 2) if tag == "R" else math.sqrt(2 * (fb + 1))
        # Folding by a group of order gamma (2 for the negation map, 6
        # with the j = 0 automorphism) divides the space every residual
        # is compared against, so the birthday floor drops by sqrt(gamma).
        gamma = g[0].get("fold_order") or (2 if g[0].get("negation_map") else 1)
        kappa_floor /= math.sqrt(gamma)
        out[key] = {
            "bits": bits,
            "B": fb,
            "dp": dp,
            "tag": tag,
            "variant": var,
            "seeds": len(g),
            "kappa": mean(r["samples"] / s for r, s in zip(g, sq)),
            # Seeded residuals have a known decomposition and collide like
            # walked ones; the honest count includes them.
            "kappa_total": mean((r["samples"] + r.get("seeded_points", 0)) / s for r, s in zip(g, sq)),
            "seeded_over_sqrt_n": mean(r.get("seeded_points", 0) / s for r, s in zip(g, sq)),
            "kappa_floor": kappa_floor,
            "c": mean(r["walk_ops"] / r["samples"] for r in g),
            "setup_frac": mean(r["setup_ops"] / t for r, t in zip(g, n_ops)),
            "replay_frac": mean(r["replay_ops"] / t for r, t in zip(g, n_ops)),
            "verify_frac": mean(r["verify_ops"] / t for r, t in zip(g, n_ops)),
            "oracle_frac": mean(r.get("oracle_ops", 0) / t for r, t in zip(g, n_ops)),
            "oracle_hits": mean(r.get("oracle_hits", 0) for r in g),
            "S": mean(t / s for t, s in zip(n_ops, sq)),
            "S_floor": kappa_floor,
            "ops_per_relation_over_sqrt_n": mean(
                t / max(r["relations_independent"], 1) / s for r, t, s in zip(g, n_ops, sq)
            ),
            "stored_over_sqrt_n": mean(r["table_entries"] / s for r, s in zip(g, sq)),
            "trivial": mean(r["collisions_trivial"] for r in g),
            "correct": all(r["correct"] for r in g),
        }
    return out


def sort_key(k):
    bits, fb, dp, tag, var = k
    return (bits, fb, dp, TAG_ORDER.index(tag) if tag in TAG_ORDER else 99, var)


def fmt(v, d=2):
    if isinstance(v, bool):
        return "yes" if v else "NO"
    if isinstance(v, float):
        if math.isnan(v):
            return "-"
        return f"{v:,.{d}f}"
    return str(v)


def print_scoreboard(c):
    print("| bits | B | dp | tag | variant | seeds | κ = samples/√n | seeded/√n | κ_total | κ floor | κ_total/floor | c ops/residual | setup | replay | verify | oracle | S = ops/√n | S floor | S/floor | ops/rel/√n | stored/√n | trivial | correct |")
    print("|---:|---:|---:|:--|:--|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:--|")
    for k in sorted(c, key=sort_key):
        x = c[k]
        print("| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |".format(
            x["bits"], x["B"], x["dp"], x["tag"], x["variant"], x["seeds"],
            fmt(x["kappa"]), fmt(x["seeded_over_sqrt_n"]), fmt(x["kappa_total"]),
            fmt(x["kappa_floor"]), fmt(x["kappa_total"] / x["kappa_floor"]),
            fmt(x["c"], 1),
            fmt(100 * x["setup_frac"], 1) + "%", fmt(100 * x["replay_frac"], 1) + "%", fmt(100 * x["verify_frac"], 1) + "%",
            fmt(100 * x["oracle_frac"], 1) + "%",
            fmt(x["S"], 1), fmt(x["S_floor"], 2), fmt(x["S"] / x["S_floor"], 2),
            fmt(x["ops_per_relation_over_sqrt_n"], 3), fmt(x["stored_over_sqrt_n"], 3),
            fmt(x["trivial"], 0), fmt(x["correct"])))
    print()


def print_ledger(c):
    """Step-by-step gains at each (bits, B) for the dp = 0 cells."""
    print("Optimisation ledger (dp = 0 cells; factors are ratios of S, > 1 means cheaper):\n")
    print("| bits | B | A → B (mutation) | A → C1 (r-adding) | A → C2 (fresh hash) | B → R (drop the base) | C1 → R | B / S floor | C1 / S floor | R / rho floor |")
    print("|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    keys = sorted({(k[0], k[1], k[4]) for k in c if k[2] == 0})
    for bits, fb, var in keys:
        def s(tag):
            return c.get((bits, fb, 0, tag, var), {}).get("S", float("nan"))
        floor = c.get((bits, fb, 0, "B", var), c.get((bits, fb, 0, "A", var), {"S_floor": math.sqrt(2 * (fb + 1))}))["S_floor"]
        print("| {} | {} ({}) | {}× | {}× | {}× | {}× | {}× | {}× | {}× | {}× |".format(
            bits, fb, var,
            fmt(s("A") / s("B")), fmt(s("A") / s("C1")), fmt(s("A") / s("C2")),
            fmt(s("B") / s("R")), fmt(s("C1") / s("R")),
            fmt(s("B") / floor), fmt(s("C1") / floor),
            fmt(s("R") / c.get((bits, fb, 0, "R", var), {"S_floor": math.sqrt(math.pi / 2)})["S_floor"])))
    print()
    dp_keys = sorted({(k[0], k[1], k[3], k[4]) for k in c if k[2] > 0})
    if dp_keys:
        print("Distinguished points (memory bought per operation spent):\n")
        print("| bits | B | tag | dp | stored/√n dp=0 | stored/√n dp | memory ÷ | S dp=0 | S dp | ops × |")
        print("|---:|---:|:--|---:|---:|---:|---:|---:|---:|---:|")
        for bits, fb, tag, var in dp_keys:
            base = c.get((bits, fb, 0, tag, var)) or c.get((bits, fb, 0, tag, "plain"))
            for k in sorted(c, key=sort_key):
                if k[:2] == (bits, fb) and k[3] == tag and k[4] == var and k[2] > 0 and base:
                    x = c[k]
                    print("| {} | {} | {} | {} | {} | {} | {}× | {} | {} | {}× |".format(
                        bits, fb, tag, k[2], fmt(base["stored_over_sqrt_n"], 3), fmt(x["stored_over_sqrt_n"], 4),
                        fmt(base["stored_over_sqrt_n"] / x["stored_over_sqrt_n"], 0),
                        fmt(base["S"], 1), fmt(x["S"], 1), fmt(x["S"] / base["S"])))
        print()


def compare(run, base, tolerance):
    print(f"Comparison against baseline (tolerance {tolerance:.0%}):\n")
    print("A run cell is matched to the baseline cell with the same (bits, B, dp, tag); "
          "the baseline variant is used when the run's variant is absent from it. "
          "κ_total (walked + seeded residuals) is compared as κ_total/κ_floor, with the floor "
          "adjusted for the fold order, so neither folding nor precomputation is mistaken for a beaten count.\n")
    print("| bits | B | dp | tag | variant | S baseline | S run | improvement (base/run) | κ_total/floor baseline | κ_total/floor run | ratio | count invariant beaten? | correct | verdict |")
    print("|---:|---:|---:|:--|:--|---:|---:|---:|---:|---:|---:|:--|:--|:--|")
    regressions = 0
    for k in sorted(run, key=sort_key):
        x = run[k]
        b = base.get(k) or base.get(k[:4] + ("plain",))
        if b is None:
            # Fall back to any baseline cell with the same (bits, B, dp, tag);
            # prefer the one with the fewest levers.
            candidates = [v for kk, v in base.items() if kk[:4] == k[:4]]
            b = min(candidates, key=lambda v: len(v["variant"]), default=None)
        if b is None:
            print("| {} | {} | {} | {} | {} | - | {} | new cell | - | {} | - | - | {} | - |".format(
                x["bits"], x["B"], x["dp"], x["tag"], x["variant"], fmt(x["S"], 1),
                fmt(x["kappa_total"] / x["kappa_floor"]), fmt(x["correct"])))
            continue
        imp = b["S"] / x["S"]
        kr = (x["kappa_total"] / x["kappa_floor"]) / (b["kappa_total"] / b["kappa_floor"])
        # Plain rho is a single-collision process whose first-collision
        # time has a wide spread; its kappa is the reference, not a target.
        # An oracle run tests each residual against virtual points it never
        # paid for, so its walked count is not comparable to the floor: it
        # is scored on S, and the count flag is withheld.
        oracle = x["oracle_frac"] > 0
        beaten = x["tag"] != "R" and kr < 1 - tolerance and x["correct"] and not oracle
        if not x["correct"]:
            verdict = "WRONG ANSWER"
            regressions += 1
        elif imp < 1 - tolerance:
            verdict = "regression"
            regressions += 1
        elif imp > 1 + tolerance:
            verdict = "improvement"
        else:
            verdict = "unchanged"
        print("| {} | {} | {} | {} | {} | {} | {} | {}× | {} | {} | {} | {} | {} | {} |".format(
            x["bits"], x["B"], x["dp"], x["tag"], x["variant"], fmt(b["S"], 1), fmt(x["S"], 1), fmt(imp),
            fmt(b["kappa_total"] / b["kappa_floor"]), fmt(x["kappa_total"] / x["kappa_floor"]), fmt(kr),
            "n/a (oracle)" if oracle else ("YES" if beaten else "no"), fmt(x["correct"]), verdict))
    print()
    return regressions


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--baseline")
    ap.add_argument("--tolerance", type=float, default=0.10)
    ap.add_argument("--fail-on-regression", action="store_true")
    ap.add_argument("--write-cells", help="write the per-cell scoreboard as JSON")
    args = ap.parse_args()

    run = cells(load(args.run))
    print_scoreboard(run)
    print_ledger(run)
    if args.write_cells:
        with open(args.write_cells, "w") as f:
            json.dump({"|".join(map(str, k)): v for k, v in sorted(run.items(), key=lambda kv: sort_key(kv[0]))}, f, indent=1)
    if args.baseline:
        base = cells(load(args.baseline))
        regressions = compare(run, base, args.tolerance)
        if regressions and args.fail_on_regression:
            print(f"{regressions} regression(s)", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()

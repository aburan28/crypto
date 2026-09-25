#!/usr/bin/env python3
"""Gates and tables for the registered runs of RESEARCH_GROEBNER_STAGE_INSTRUCTIONS.md.

    python3 research/groebner_instructions_20260925/table.py > research/groebner_instructions_20260925/tables.md

G1  identity: every counter (every field but *_ns) of every stage rung equals
    the registered run of the same arm (A1, A2: round 2; A0: round 1's
    reference, on the keys both harness versions record); whole-log trials,
    relations and oracle word operations of A1 and A2 equal round 2's.
G2  determinism: each rerun's stage instructions within 1e-6 of the run's.
G3  every whole logarithm and every rho run verified.
Then: the stage table (word operations and instructions, all arms, rounds 1
and 2 in both units), the phase shares of the default, and the whole-log S
from process instructions against rho measured the same way.  Exits non-zero
if a gate fails.
"""
import json, math, pathlib, statistics, sys

HERE = pathlib.Path(__file__).resolve().parent
R1 = HERE.parent / "chain_split_order_20260924"
R2 = HERE.parent / "support_local_multipliers_20260924"
SUITES = ["frozen", "chain", "chain-holdout", "chain-holdout-2", "r2-holdout"]
ARMS = ["A0", "A1", "A2"]
REGISTERED = {"A0": (R1, "reference"), "A1": (R2, "reference"), "A2": (R2, "candidate")}
PHASES = {
    "system construction": ["crypto_lib::cryptanalysis::koblitz_groebner::build_decomposition_system"],
    "root builds": ["crypto_lib::cryptanalysis::inherited_f4::ReducedBasis::from_system_with"],
    "specialisation": ["crypto_lib::cryptanalysis::koblitz_groebner::InheritedBases::specialise_owned"],
    "of which insert": ["crypto_lib::cryptanalysis::inherited_f4::ReducedBasis::insert"],
    "linear elimination": ["crypto_lib::cryptanalysis::koblitz_groebner::eliminate_linear_generators"],
    "counted elimination kernels": [
        "crypto_lib::cryptanalysis::koblitz_groebner::echelon_f2_suffix_counted",
        "crypto_lib::cryptanalysis::koblitz_groebner::echelon_f2_m4ri_counted",
        "crypto_lib::cryptanalysis::koblitz_groebner::echelon_f2_counted",
    ],
}
CELLS = {
    "K0_2^13": list(range(201, 211)) + list(range(301, 306)),
    "K0_2^9": list(range(201, 206)),
}
bad = []


def load(p):
    return json.loads(p.read_text())


def counters(row):
    return {k: v for k, v in row.items() if not k.endswith("_ns")}


def stage_rows():
    """suite -> list of (rung index, {arm: summary})"""
    out = {}
    for suite in SUITES:
        rungs = sorted((HERE / "stage" / suite / "A2").glob("rung*/summary.json"), key=lambda p: int(p.parent.name[4:]))
        rows = []
        for p in rungs:
            if load(p).get("skipped"):
                continue
            i = int(p.parent.name[4:])
            rows.append((i, {a: load(HERE / "stage" / suite / a / f"rung{i}" / "summary.json") for a in ARMS}))
        out[suite] = rows
    return out


def gate_identity(rows):
    print("## G1 — identity with the registered counters\n")
    for suite in SUITES:
        for arm in ARMS:
            root, name = REGISTERED[arm]
            reg_path = root / suite / name / "rep1" / "stage.json"
            if not reg_path.exists():
                print(f"- {suite}, {arm}: no registered run (reported without the gate)")
                continue
            reg = load(reg_path)["rows"]
            # rungs in registered order = admissible rungs in ladder order
            mine = [s[arm]["row"] for _, s in rows[suite]]
            if len(mine) != len(reg):
                bad.append(f"G1 {suite} {arm}: {len(mine)} rungs vs {len(reg)} registered")
                print(f"- {suite}, {arm}: **rung count differs** ({len(mine)} vs {len(reg)})")
                continue
            diffs, keys = 0, None
            for a, b in zip(mine, reg):
                ca, cb = counters(a), counters(b)
                common = sorted(set(ca) & set(cb))
                keys = common
                if any(ca[k] != cb[k] for k in common):
                    diffs += 1
                    bad.append(f"G1 {suite} {arm} {a['curve']} m={a['m']}")
            print(f"- {suite}, {arm}: {len(mine)} rungs, {diffs} differing, on {len(keys)} fields")
    for arm, name in (("A1", "reference"), ("A2", "candidate")):
        diffs = 0
        for cell, seeds in CELLS.items():
            for s in seeds:
                mine = load(HERE / "e2e" / arm / f"{cell}_seed{s}" / "summary.json")["counts"]
                reg = load(R2 / "e2e" / name / f"{cell}_seed{s}.json")["counts"]
                if any(mine[k] != reg[k] for k in ("trials", "relations", "f4_word_ops")):
                    diffs += 1
                    bad.append(f"G1 e2e {arm} {cell} {s}")
        print(f"- whole logarithms, {arm}: 20 seeds, {diffs} differing from round 2 in trials, relations, oracle word ops")


def gate_determinism():
    print("\n## G2 — determinism\n")
    for p in sorted((HERE / "determinism").glob("*/*/rung*/summary.json")):
        suite, arm, rung = p.parts[-4], p.parts[-3], p.parts[-2]
        a = load(p)["ir"]
        b = load(HERE / "stage" / suite / arm / rung / "summary.json")["ir"]
        rel = abs(a - b) / b
        if rel >= 1e-6:
            bad.append(f"G2 {suite} {arm} {rung}")
        print(f"- {suite} {rung} {arm}: {b:,} vs {a:,} (relative {rel:.1e})")


def stage_table(rows):
    print("\n## Stage ladders: word operations and instructions\n")
    print("Instructions are those executed inside `groebner_decompose`; ratios are old / new.\n")
    print("| suite | rung | word ops A0 / A1 / A2 | instructions A0 / A1 / A2 | instr. per word op (A2) | round 1 (A0/A1) words, instr. | round 2 (A1/A2) words, instr. |")
    print("|:--|:--|:--|:--|--:|:--|:--|")
    for suite in SUITES:
        tw = {a: 0 for a in ARMS}
        ti = {a: 0 for a in ARMS}
        for _, s in rows[suite]:
            w = {a: s[a]["row"]["word_ops"] for a in ARMS}
            ir = {a: s[a]["ir"] for a in ARMS}
            for a in ARMS:
                tw[a] += w[a]
                ti[a] += ir[a]
            r = s["A2"]["row"]
            name = f"`{r['curve']}` m={r['m']}" + (f" div {r['factor_index']}" if r.get("factor_index") else "") + (
                f" (targets {r['first_target']}…)" if r.get("first_target") else "")
            print(
                f"| {suite} | {name} | {w['A0']:,} / {w['A1']:,} / {w['A2']:,} | {ir['A0']:,} / {ir['A1']:,} / {ir['A2']:,} "
                f"| {ir['A2'] / w['A2']:.0f} | {w['A0'] / w['A1']:.2f}×, {ir['A0'] / ir['A1']:.2f}× | {w['A1'] / w['A2']:.2f}×, {ir['A1'] / ir['A2']:.2f}× |"
            )
        print(
            f"| **{suite} total** | | {tw['A0']:,} / {tw['A1']:,} / {tw['A2']:,} | {ti['A0']:,} / {ti['A1']:,} / {ti['A2']:,} "
            f"| {ti['A2'] / tw['A2']:.0f} | **{tw['A0'] / tw['A1']:.2f}×, {ti['A0'] / ti['A1']:.2f}×** | **{tw['A1'] / tw['A2']:.2f}×, {ti['A1'] / ti['A2']:.2f}×** |"
        )


def phase_table(rows):
    print("\n## Where the default's stage instructions go (A2)\n")
    print("| suite | " + " | ".join(PHASES) + " | instructions |")
    print("|:--|" + "--:|" * (len(PHASES) + 1))
    grand = {k: 0 for k in PHASES}
    gt = 0
    for suite in SUITES:
        acc = {k: 0 for k in PHASES}
        t = 0
        for _, s in rows[suite]:
            inc = s["A2"]["inclusive_ir"]
            t += s["A2"]["ir"]
            for k, fs in PHASES.items():
                acc[k] += sum(inc.get(f, 0) for f in fs)
        for k in PHASES:
            grand[k] += acc[k]
        gt += t
        print(f"| {suite} | " + " | ".join(f"{100 * acc[k] / t:.0f}%" for k in PHASES) + f" | {t:,} |")
    print("| **all rungs** | " + " | ".join(f"**{100 * grand[k] / gt:.0f}%**" for k in PHASES) + f" | {gt:,} |")


def whole_logs():
    print("\n## Whole logarithms: `S` from process instructions\n")
    per_add = {cell: load(HERE / "add" / cell / "summary.json")["ir_per_add"] for cell in CELLS}
    print("Conversion: " + ", ".join(f"{c}: {v:.1f} instructions per addition" for c, v in per_add.items()) + ".\n")
    print("| cell | arm | runs verified | trials | process instructions (sum) | `S` per run, mean [min–max] | ratio to rho (means) |")
    print("|:--|:--|--:|--:|--:|:--|--:|")
    for cell, seeds in CELLS.items():
        rho = [load(HERE / "rho" / f"{cell}_seed{s}" / "summary.json") for s in seeds]
        for x in rho:
            if not (x["result"]["verified"]):
                bad.append(f"G3 rho {cell} {x['seed']}")
        rho_s = [x["ir"] / per_add[cell] / math.sqrt(x["r"]) for x in rho]
        rho_mean = statistics.fmean(rho_s)
        counted = statistics.fmean(x["result"]["s"] for x in rho)
        print(f"| `{cell}` | rho | {sum(x['result']['verified'] for x in rho)} / {len(rho)} | — | {sum(x['ir'] for x in rho):,} "
              f"| {rho_mean:.1f} [{min(rho_s):.1f}–{max(rho_s):.1f}] (counted: {counted:.1f}) | 1 |")
        for arm in ARMS:
            runs = [load(HERE / "e2e" / arm / f"{cell}_seed{s}" / "summary.json") for s in seeds]
            for s, x in zip(seeds, runs):
                if not x["verified"]:
                    bad.append(f"G3 e2e {arm} {cell} {s}")
            ss = [x["ir"] / per_add[cell] / math.sqrt(x["r"]) for x in runs]
            print(f"| `{cell}` | {arm} | {sum(x['verified'] for x in runs)} / {len(runs)} | {sum(x['counts']['trials'] for x in runs)} "
                  f"| {sum(x['ir'] for x in runs):,} | {statistics.fmean(ss):,.0f} [{min(ss):,.0f}–{max(ss):,.0f}] | {statistics.fmean(ss) / rho_mean:,.0f}× |")
    print("\nRound ratios in process instructions (sum over the cell's seeds; old / new):\n")
    for cell, seeds in CELLS.items():
        t = {a: sum(load(HERE / "e2e" / a / f"{cell}_seed{s}" / "summary.json")["ir"] for s in seeds) for a in ARMS}
        st = {a: sum(load(HERE / "e2e" / a / f"{cell}_seed{s}" / "summary.json")["inclusive_ir"].get(
            "crypto_lib::cryptanalysis::koblitz_index_calculus::groebner_decompose", 0) for s in seeds) for a in ARMS}
        print(f"- `{cell}`: round 1 {t['A0'] / t['A1']:.2f}×, round 2 {t['A1'] / t['A2']:.2f}×, both {t['A0'] / t['A2']:.2f}×; "
              f"the oracle's share of the process (A2) {100 * st['A2'] / t['A2']:.0f}%")


def main():
    rows = stage_rows()
    gate_identity(rows)
    gate_determinism()
    stage_table(rows)
    phase_table(rows)
    whole_logs()
    print("\n## Gates\n")
    print("All gates pass." if not bad else "**Failures:**\n" + "\n".join(f"- {b}" for b in bad))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()

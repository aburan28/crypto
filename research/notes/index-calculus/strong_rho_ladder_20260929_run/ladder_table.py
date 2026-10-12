#!/usr/bin/env python3
"""Build the ladder table from the raw logs (no hand-typed numbers).

Sources:
  callgrind_rung{0..3}.stderr.log   'I   refs:' line = whole-process Ir
  callgrind_ic.stderr.log           same, for the IC arm on the current source
  ../matched_rho_orbit_dlp_20260928_run/callgrind_ks_matched.annotate.txt
                                    PR #955's R0 (its own binary), PROGRAM TOTALS
  rung{0..3}_native.jsonl           steps, verification, in-process wall
Prints a markdown table and the pre-registered decision.
"""
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
PR955 = os.path.join(
    HERE, "..", "matched_rho_orbit_dlp_20260928_run", "callgrind_ks_matched.annotate.txt"
)


def ir_from_log(path):
    for line in open(path):
        if "I   refs:" in line:
            return int(line.split("refs:")[1].strip().replace(",", ""))
    raise SystemExit(f"no Ir in {path}")


def ir_from_annotate(path):
    for line in open(path):
        m = re.match(r"\s*([\d,]+)\s+\(100\.0%\)\s+PROGRAM TOTALS", line)
        if m:
            return int(m.group(1).replace(",", ""))
    raise SystemExit(f"no PROGRAM TOTALS in {path}")


def native(rung):
    fixtures, summary = 0, None
    for line in open(os.path.join(HERE, f"rung{rung}_native.jsonl")):
        obj = json.loads(line)
        if obj["kind"] == "rho_ks_batch_fixture":
            fixtures += obj["verified"] and (
                obj["recovered_fixture_scalar"] == obj["published_fixture_scalar"]
            )
        else:
            summary = obj
    return fixtures, summary


def main():
    ic = ir_from_log(os.path.join(HERE, "callgrind_ic.stderr.log"))
    rows = [("R0 as merged in PR #955 (its own binary)", ir_from_annotate(PR955), None)]
    labels = {
        0: "R0 regression: this binary, rho-local software field",
        1: "R1 library Gf2 field (hardware clmul, table reduction)",
        2: "R2 + normal-coordinate canonicalization, fast mulmod",
        3: "R3 + 32 lockstep lanes, Gf2::batch_inv",
    }
    passing = {}
    for r in range(4):
        ir = ir_from_log(os.path.join(HERE, f"callgrind_rung{r}.stderr.log"))
        ok, summary = native(r)
        rows.append((labels[r], ir, (ok, summary)))
        if ok == 1024 and summary["all_verified"]:
            passing[r] = ir
    print("| reference rho | Ir (whole process) | IC / rho | steps | verified |")
    print("|:--|--:|--:|--:|--:|")
    for label, ir, meta in rows:
        steps = f"{meta[1]['total_walk_steps']:,}" if meta else "19,103,507"
        ver = f"{meta[0]:,}/1,024" if meta else "1,024/1,024"
        print(f"| {label} | {ir:,} | {ic / ir:.3f} | {steps} | {ver} |")
    print(f"| IC `koblitz_orbit_dlp_fast` (current source) | {ic:,} | — | — | 1,024/1,024 |")
    eligible = {r: ir for r, ir in passing.items() if r in (1, 2, 3)}
    best_r = min(eligible, key=eligible.get)
    rho_star = ic / eligible[best_r]
    verdict = "survives" if rho_star < 0.8 else ("dies" if rho_star >= 1.0 else "inconclusive")
    print()
    print(f"rho_best = R{best_r} ({eligible[best_r]:,} Ir); rho* = IC/rho_best = {rho_star:.3f} -> lead {verdict}")


if __name__ == "__main__":
    main()

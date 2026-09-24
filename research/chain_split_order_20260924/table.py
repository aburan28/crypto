#!/usr/bin/env python3
"""Tabulate the registered runs of RESEARCH_CHAIN_SPLIT_ORDER.md.

    python3 research/chain_split_order_20260924/table.py [ROOT] > tables.md

One table per suite: rows are rungs, columns are the arms of the O x L x D
factorial, cells are median word operations over the three repetitions with
the ratio reference / arm.  Beside them, per rung, the reference and
candidate node counts, targets decomposed and exhausted budgets, and a flag
for a rung whose counters differ between repetitions (which would void it).
Unit: 64-bit word operations of the Groebner stage; a stage diagnostic.
"""
import json, pathlib, statistics, sys

ARMS = ["reference", "D", "L", "LD", "O", "OD", "candidate", "OLD"]
SUITES = ["frozen", "chain", "chain-holdout", "chain-holdout-2"]


def load(root, suite, arm):
    reps = sorted((root / suite / arm).glob("rep*/stage.json"))
    return [json.loads(p.read_text()) for p in reps]


def main():
    root = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else pathlib.Path(__file__).parent)
    for suite in SUITES:
        docs = {arm: load(root, suite, arm) for arm in ARMS}
        docs = {a: d for a, d in docs.items() if d}
        if "reference" not in docs:
            continue
        ref = docs["reference"]
        print(f"\n### {suite}\n")
        head = "| rung | targets | decomposed | nodes ref → cand | " + " | ".join(docs) + " |"
        print(head)
        print("|:--|--:|--:|:--|" + "--:|" * len(docs))
        totals = {a: 0 for a in docs}
        for i, row in enumerate(ref[0]["rows"]):
            name = f"`{row['curve']}` m={row['m']}" + (
                f" divisor {row['factor_index']}" if row.get("factor_index") else ""
            ) + (
                f" (targets {row.get('first_target', 0)}…)" if row.get("first_target") else ""
            )
            cells = []
            flags = []
            for arm, d in docs.items():
                rows = [x["rows"][i] for x in d]
                assert all((r["curve"], r["m"]) == (row["curve"], row["m"]) for r in rows)
                ops = [r["word_ops"] for r in rows]
                if len(set(ops)) != 1:
                    flags.append(f"{arm}: word ops differ across reps")
                med = statistics.median(ops)
                totals[arm] += med
                if arm == "reference":
                    cells.append(f"{med:,.0f}")
                else:
                    ratio = statistics.median(x["rows"][i]["word_ops"] for x in ref) / med
                    cells.append(f"{med:,.0f} ({ratio:.2f}×)")
                if any(r["exhausted"] for r in rows):
                    flags.append(f"{arm}: {rows[0]['exhausted']} exhausted")
                if rows[0]["decomposed"] != row["decomposed"]:
                    flags.append(f"{arm}: decomposed {rows[0]['decomposed']}")
            cand = docs.get("candidate", ref)[0]["rows"][i]
            nodes = f"{row['reductions']:,} → {cand['reductions']:,}"
            print(f"| {name} | {row['targets']} | {row['decomposed']} | {nodes} | " + " | ".join(cells) + " |")
            for f in flags:
                print(f"|  ⚠ {f} | | | |" + " |" * len(docs))
        tot = ["**{:,.0f}**".format(totals["reference"])] + [
            "**{:,.0f} ({:.2f}×)**".format(totals[a], totals["reference"] / totals[a]) for a in docs if a != "reference"
        ]
        print("| **total** | | | | " + " | ".join(tot) + " |")


if __name__ == "__main__":
    main()

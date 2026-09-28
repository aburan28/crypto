#!/usr/bin/env python3
"""Summarise one callgrind job of run.sh into summary.json beside it.

    python3 summarise.py {stage|e2e|rho|add} DIR

Reads DIR/callgrind.out (the collected instructions, `summary:` line) and the
job's own output, and keeps the inclusive instructions of every function at
or above 0.2% of the total, so the phase shares can be audited without the
raw profile.  A stage job whose rung is inadmissible (no rows) is recorded
as skipped.
"""
import json, pathlib, re, subprocess, sys


def total(d):
    for line in (d / "callgrind.out").read_text(errors="replace").splitlines():
        if line.startswith("summary:"):
            return int(line.split()[1])
    raise SystemExit(f"{d}: no summary line")


def inclusive(d, t):
    out = subprocess.run(
        ["callgrind_annotate", "--inclusive=yes", "--threshold=100", str(d / "callgrind.out")],
        capture_output=True, text=True, check=True,
    ).stdout
    funcs = {}
    for line in out.splitlines():
        m = re.match(r"\s*([\d,]+)\s+\(\s*[\d.]+%\)\s+(\S.*)$", line)
        if not m or m.group(2).startswith("PROGRAM TOTALS"):
            continue
        name = re.sub(r"\s*\[.*\]$", "", m.group(2)).replace("???:", "")
        ir = int(m.group(1).replace(",", ""))
        if ir >= 0.002 * t:
            funcs[name] = max(funcs.get(name, 0), ir)
    return funcs


def main():
    kind, d = sys.argv[1], pathlib.Path(sys.argv[2])
    doc = {"kind": kind}
    if kind == "stage":
        rows = json.loads((d / "stage" / "stage.json").read_text())["rows"]
        if not rows:
            (d / "summary.json").write_text(json.dumps({"kind": kind, "skipped": True}) + "\n")
            (d / "callgrind.out").unlink(missing_ok=True)
            return
        t = total(d)
        doc.update(row=rows[0], ir=t, inclusive_ir=inclusive(d, t))
    elif kind == "e2e":
        run = json.loads((d / "run.json").read_text())
        t = total(d)
        doc.update(
            ir=t,
            verified=run["result"]["verified"] and run["result"]["expected"] == run["result"]["recovered"],
            r=int(run["parameters"]["subgroup_order"]),
            counts={k: run["counts"][k] for k in ("trials", "relations", "f4_word_ops", "f4_reductions")},
            elapsed_seconds=run["elapsed_seconds"],
            inclusive_ir=inclusive(d, t),
        )
    elif kind == "rho":
        run = json.loads((d / "run.json").read_text())
        doc.update(ir=total(d), r=run["r"], seed=run["seed"], result=run["result"])
    elif kind == "add":
        run = json.loads((d / "run.json").read_text())
        t = total(d)
        doc.update(ir=t, samples=run["samples"], r=run["r"], ir_per_add=t / run["samples"])
    (d / "summary.json").write_text(json.dumps(doc, indent=1) + "\n")


if __name__ == "__main__":
    main()

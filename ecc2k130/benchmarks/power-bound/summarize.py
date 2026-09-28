#!/usr/bin/env python3
"""Freeze a power-bound job run into summary.json (POWER-BOUND.md).

    python3 benchmarks/power-bound/summarize.py /tmp/power-bound benchmarks/power-bound/rtx-pro-6000

For each bench the job wrote the rate (`finished:` line) and a 100 ms power
trace taken while the binary ran.  The trace covers setup as well as the timed
launches, so the loaded window is taken as the samples drawing at least half
the trace's peak power, less the first five of them (clock ramp); mean power
over that window divided by the rate is the energy per update.  Cites; does
not decide.
"""
import csv
import json
import pathlib
import re
import shutil
import statistics
import sys


def window(path):
    rows = []
    with open(path) as f:
        for r in csv.reader(f):
            try: rows.append((float(r[1]), float(r[2]), float(r[3])))
            except (ValueError, IndexError): pass
    if not rows: return None
    peak = max(p for p, _, _ in rows)
    loaded = [r for r in rows if r[0] >= 0.5 * peak][5:]
    if len(loaded) < 5: return None
    return {"samples": len(loaded), "meanPowerW": statistics.mean(p for p, _, _ in loaded),
            "meanClockMHz": statistics.mean(c for _, c, _ in loaded), "maxTempC": max(t for _, _, t in loaded)}


def main():
    run = pathlib.Path(sys.argv[1]).resolve(); target = pathlib.Path(sys.argv[2]).resolve()
    res = run / "results"
    out = {}
    for line in (res / "bench.txt").read_text().splitlines():
        m = re.match(r"(uncapped|capped) (\S+) rep (\d+):\s+finished: ([0-9.]+) M it/s", line)
        if not m: continue
        tag, name, rep, rate = m.group(1), m.group(2), int(m.group(3)), float(m.group(4)) / 1000
        w = window(res / ("power-%s-%s-%d.csv" % (tag, name, rep)))
        e = {"rep": rep, "rateBps": rate, **(w or {})}
        if w: e["updatesPerJouleM"] = rate * 1e3 / w["meanPowerW"]
        out.setdefault(tag, {}).setdefault(name, []).append(e)
    summary = {"note": "POWER-BOUND.md", "host": (res / "host.txt").read_text().strip(),
               "powercap": (res / "powercap.txt").read_text().strip() if (res / "powercap.txt").exists() else None,
               "launch": json.loads((run / "launch.json").read_text()), "passes": {}}
    for tag, variants in out.items():
        rows = {}
        for name, reps in variants.items():
            rows[name] = {"reps": reps, "medianBps": statistics.median(r["rateBps"] for r in reps),
                          "medianPowerW": statistics.median(r["meanPowerW"] for r in reps if "meanPowerW" in r) if any("meanPowerW" in r for r in reps) else None,
                          "medianClockMHz": statistics.median(r["meanClockMHz"] for r in reps if "meanClockMHz" in r) if any("meanClockMHz" in r for r in reps) else None,
                          "medianUpdatesPerJouleM": statistics.median(r["updatesPerJouleM"] for r in reps if "updatesPerJouleM" in r) if any("updatesPerJouleM" in r for r in reps) else None}
            v = (res / ("verify-%s.txt" % name))
            rows[name]["verify"] = v.read_text().strip() if v.exists() else None
        summary["passes"][tag] = rows
    target.mkdir(parents=True, exist_ok=True)
    (target / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    raw = target / "raw"; raw.mkdir(exist_ok=True)
    for f in res.iterdir():
        if f.suffix in (".txt", ".log", ".csv") and f.stat().st_size < 4 << 20: shutil.copy(f, raw / f.name)
    shutil.copy(run / "launch.json", raw / "launch.json")
    for tag, rows in summary["passes"].items():
        print("== %s" % tag)
        for n, r in rows.items():
            print("  %-15s %.3f B/s  %s W  %s MHz  %s M upd/J" % (n, r["medianBps"], "%.1f" % r["medianPowerW"] if r["medianPowerW"] else "-",
                  "%.0f" % r["medianClockMHz"] if r["medianClockMHz"] else "-", "%.2f" % r["medianUpdatesPerJouleM"] if r["medianUpdatesPerJouleM"] else "-"))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Freeze one autosweep run into summary.json (AUTOSWEEP.md).

    python3 benchmarks/autosweep/summarize.py /tmp/autosweep-b200 benchmarks/autosweep/b200

Reads the job's sweep.json (arms, builds, screening rates, star ranking,
greedy trace), the verify-*.txt / dp-identity.txt of the finalists and the
final alternating bench from job.log, and writes summary.json plus the raw
text files into the target directory.  Cites; does not decide.
"""
import json
import pathlib
import re
import shutil
import statistics
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "two-chains"))
import summarize as S  # noqa: E402


def main():
    run = pathlib.Path(sys.argv[1]).resolve(); target = pathlib.Path(sys.argv[2]).resolve()
    results = run / "results"
    launch = json.loads((run / "launch.json").read_text())
    sweep = json.loads((results / "sweep.json").read_text())
    knobs = {"base": ""}
    for arm in sweep["arms"]:
        name, _, k = arm.partition("|"); knobs[name] = k
    for line in sweep["greedy"]:
        m = re.match(r"greedy (combo\d+) = (\S+) \+ (\S+):", line)
        if m: knobs[m.group(1)] = (knobs.get(m.group(2), "") + " " + knobs[m.group(3)]).strip()
    star = []
    base_screen = None
    for line in sweep["star"]:
        if line.startswith("base "): base_screen = float(line.split()[1]); continue
        ratio, name, med = line.split(); star.append({"arm": name, "knobs": knobs.get(name, ""), "screenedBps": float(med),
                                                      "ratioToBase": float(ratio), **sweep["build"].get(name, {})})
    greedy = []
    for line in sweep["greedy"]:
        m = re.match(r"greedy (combo\d+) = (\S+) \+ (\S+): ([0-9.]+) B/s vs ([0-9.]+)", line)
        if m: greedy.append({"name": m.group(1), "from": m.group(2), "added": m.group(3), "screenedBps": float(m.group(4)),
                             "against": float(m.group(5)), "knobs": knobs[m.group(1)], **sweep["build"].get(m.group(1), {})})
        elif line.strip() in ("keep", "drop") and greedy: greedy[-1]["decision"] = line.strip()
    result_line = [l for l in sweep["greedy"] if l.startswith("greedy result")]
    finalists = [re.match(r"=== bench (\S+) rep 1", l).group(1) for l in (results / "job.log").read_text().splitlines() if re.match(r"=== bench (\S+) rep 1", l)]
    S.VARIANTS = finalists
    bench = S.parse_bench((results / "job.log").read_text())
    identity = S.parse_identity(results / "dp-identity.txt")
    final = {}
    ref = finalists[0] if finalists else "base"
    for name in finalists:
        samples = bench.get(name, []); rates = [s[0] for s in samples]
        e = {"knobs": knobs.get(name, ""), "samplesBps": rates, "medianBps": statistics.median(rates) if rates else None,
             "smClockMHz": [s[1] for s in samples], "powerW": [s[2] for s in samples], "tempC": [s[3] for s in samples],
             "verify": S.parse_verify(results / ("verify-%s.txt" % name)), "dpIdentity": identity.get(name),
             **sweep["build"].get(name, {})}
        if name != ref and rates and bench.get(ref):
            pairs = [r / bench[ref][i][0] for i, r in enumerate(rates) if i < len(bench[ref])]
            e.update(pairedRatioToBase=pairs, pairedRatioMedian=statistics.median(pairs), pairedRatioMin=min(pairs), pairedRatioMax=max(pairs))
        final[name] = e
    summary = {"note": "AUTOSWEEP.md", "unit": "billions of complete scalar updates per second, `finished:` line",
               "host": (results / "host.txt").read_text().strip(),
               "launch": {k: launch.get(k) for k in ("provider", "gpu", "cudaImage", "gitRev", "gitDirty", "startedAt", "finishedAt", "exitCode", "seconds", "deviceLine", "token")},
               "baseKnobs": [l for l in (results / "host.txt").read_text().splitlines() if l.startswith("base:")],
               "screenCommand": "--bench --steps 1024 --launches 32, two passes alternating, median",
               "finalCommand": "--bench --steps 1024 --launches 64, alternating, five repetitions",
               "baseScreenedBps": base_screen, "star": star, "greedy": greedy, "greedyResult": result_line,
               "final": final}
    target.mkdir(parents=True, exist_ok=True)
    (target / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    raw = target / "raw"; raw.mkdir(exist_ok=True)
    for f in results.iterdir():
        if f.suffix in (".txt", ".log", ".json") and f.stat().st_size < 4 << 20: shutil.copy(f, raw / f.name)
    shutil.copy(run / "launch.json", raw / "launch.json")
    print("base screened %.3f" % base_screen)
    for s in star: print("  %-22s %.4f  %.3f  regs %s spill %s  %s" % (s["arm"], s["ratioToBase"], s["screenedBps"], s.get("registers"), s.get("spillBytes"), s["knobs"]))
    for g in greedy: print("  %s = %s + %s: %.3f vs %.3f -> %s" % (g["name"], g["from"], g["added"], g["screenedBps"], g["against"], g.get("decision")))
    for n, e in final.items():
        print("final %-14s median %s paired %s verified %s/%s dropped %s DP %s" % (
            n, "%.3f" % e["medianBps"] if e["medianBps"] else "-", "%.4f" % e["pairedRatioMedian"] if e.get("pairedRatioMedian") else "-",
            e["verify"]["verified"], e["verify"]["distinguishedPoints"], e["verify"]["dropped"],
            "identical" if (e["dpIdentity"] or {}).get("identicalToRef") else e["dpIdentity"]))


if __name__ == "__main__":
    main()

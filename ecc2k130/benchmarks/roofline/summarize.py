#!/usr/bin/env python3
"""Freeze one run of benchmarks/roofline/gpujob.sh into summary.json (ROOFLINE.md section 6).

    python3 benchmarks/roofline/summarize.py /tmp/roofline-run1 benchmarks/roofline

Reads the launcher's launch.json and the job's results/: the pipe probe
(pipes.jsonl) with the SASS audit of every timed loop, each binary's roofline
prediction (roofline-<arm>.json), verification, distinguished-point identity,
the alternating bench, and Nsight Compute output where the host allowed it.
Writes summary.json and the raw text files into the target directory, and for
every arm the measured rate beside the rate the roofline predicted for it (the
reference's measured rate times the ratio of the binding pipe's work, under
each pipe model of ROOFLINE.md section 3).  It cites; it does not decide.
"""
import importlib.util
import json
import pathlib
import shutil
import statistics
import sys

HERE = pathlib.Path(__file__).resolve().parent
# The two-chains job writes the same bench, verify and identity files.
_spec = importlib.util.spec_from_file_location("two_chains_summarize", HERE.parent / "two-chains" / "summarize.py")
two_chains = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(two_chains)
parse_bench, parse_identity, parse_verify = two_chains.parse_bench, two_chains.parse_identity, two_chains.parse_verify

VARIANTS = ["ref", "sqtab", "invpoly", "invpoly2", "both", "both2", "fused2"]
MODELS = {  # lanes per SM-clock, ROOFLINE.md section 3
    "M4 alu 64 / fma 64 / clmad 2.0": {"alu": 64.0, "fma": 64.0, "clmad": 2.0},
    "M5 int (alu+fma) 69 / clmad 2.0": {"int": 69.0, "clmad": 2.0},
}


def pipe_work(per_update, pipe):
    if pipe == "int":
        return per_update.get("alu", 0.0) + per_update.get("fma", 0.0)
    return per_update.get(pipe, 0.0)


def floor_clocks(per_update, model):
    return max(pipe_work(per_update, p) / rate for p, rate in model.items())


def main():
    run, target = pathlib.Path(sys.argv[1]).resolve(), pathlib.Path(sys.argv[2]).resolve()
    results = run / "results"
    launch = json.loads((run / "launch.json").read_text()) if (run / "launch.json").exists() else {}
    two_chains.VARIANTS = VARIANTS
    src = results / "job.log" if (results / "job.log").exists() else results / "bench.txt"
    bench = parse_bench(src.read_text()) if src.exists() else {}
    identity = parse_identity(results / "dp-identity.txt")
    roof = {}
    for name in VARIANTS:
        p = results / ("roofline-%s.json" % name)
        if p.exists():
            d = json.loads(p.read_text())
            roof[name] = {"perUpdate": d["perUpdate"], "registers": d["registers"], "spills": d["spills"],
                          "expectedClmad": d["expectedClmad"]}
    variants = {}
    for name in VARIANTS:
        rates = [s[0] for s in bench.get(name, [])]
        e = {"samplesBps": rates, "medianBps": statistics.median(rates) if rates else None,
             "smClockMHz": [s[1] for s in bench.get(name, [])], "powerW": [s[2] for s in bench.get(name, [])],
             "verify": parse_verify(results / ("verify-%s.txt" % name)), "dpIdentity": identity.get(name),
             "roofline": roof.get(name)}
        if name != "ref" and rates and bench.get("ref"):
            pairs = [r / bench["ref"][i][0] for i, r in enumerate(rates) if i < len(bench["ref"])]
            e.update(pairedRatioToRef=pairs, pairedRatioMedian=statistics.median(pairs))
        if name != "ref" and name in roof and "ref" in roof:
            # predicted ratio: binding-pipe floor of the reference over this arm's
            e["predictedRatio"] = {m: floor_clocks(roof["ref"]["perUpdate"], mod) /
                                   floor_clocks(roof[name]["perUpdate"], mod) for m, mod in MODELS.items()}
        variants[name] = e
    pipes, audit = [], {}
    if (results / "pipes-sass-loops.jsonl").exists():
        for line in (results / "pipes-sass-loops.jsonl").read_text().splitlines():
            d = json.loads(line)
            audit[d["op"]] = d
    if (results / "pipes.jsonl").exists():
        k = 0
        for line in (results / "pipes.jsonl").read_text().splitlines():
            if line.startswith("{\"stream\""):
                d = json.loads(line)
                d["sassLoop"] = audit.get(k)
                pipes.append(d)
                k += 1
    summary = {
        "note": "ROOFLINE.md",
        "unit": "billions of complete scalar updates per second, `finished:` line",
        "host": (results / "host.txt").read_text().strip() if (results / "host.txt").exists() else None,
        "launch": {k: launch.get(k) for k in ("provider", "gpu", "cudaImage", "gitRev", "gitDirty", "startedAt",
                                              "finishedAt", "exitCode", "seconds", "deviceLine")},
        "pipeModels": MODELS,
        "variants": variants,
        "pipes": pipes,
        "ncu": {f.name: f.read_text()[-20000:] for f in sorted(results.glob("ncu-*.txt"))
                if f.name != "ncu-query.txt"},
    }
    target.mkdir(parents=True, exist_ok=True)
    (target / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    raw = target / "raw"
    raw.mkdir(exist_ok=True)
    for f in results.iterdir():
        if f.suffix in (".txt", ".log", ".jsonl", ".json") and f.stat().st_size < 4 << 20 and f.name != "ncu-query.txt":
            shutil.copy(f, raw / f.name)
    if (run / "launch.json").exists():
        shutil.copy(run / "launch.json", raw / "launch.json")
    for name, e in variants.items():
        print("%-8s median %s  paired %s  predicted %s  verified %s/%s  DP %s" % (
            name, "%.3f" % e["medianBps"] if e["medianBps"] else "-",
            "%.4f" % e["pairedRatioMedian"] if e.get("pairedRatioMedian") else "-",
            {m.split()[0]: round(v, 4) for m, v in e.get("predictedRatio", {}).items()} or "-",
            e["verify"]["verified"], e["verify"]["distinguishedPoints"],
            "identical" if (e["dpIdentity"] or {}).get("identicalToRef") else e["dpIdentity"]))
    for d in pipes:
        print("%-62s %7.3f lanes/SM-clk" % (d["stream"][:62], d["lanesPerSmClock"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())

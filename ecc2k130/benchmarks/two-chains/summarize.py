#!/usr/bin/env python3
"""Freeze one GPU job run into summary.json (TWO-CHAINS.md sections 5 and 6).

    python3 benchmarks/two-chains/summarize.py /tmp/two-chains-run1 benchmarks/two-chains
    python3 benchmarks/two-chains/summarize.py /tmp/fast-clmad-b200 benchmarks/fast-clmad \
        --variants ref,clsq,topclmad,topclmad-clsq,onbinv,c2-256x32 --profiles ref-prof,topclmad-clsq-prof

Reads the launcher's launch.json and the job's results/ (job.log, verify-*.txt,
dp-identity.txt, profile-*.txt, build-*.txt, host.txt), computes per-variant
medians and the paired ratios against the reference within each repetition,
and writes summary.json plus the raw text files into the target directory.
It cites; it does not decide: the class and verdict columns are filled by the
note against the target it declared.
"""
import json
import pathlib
import re
import shutil
import statistics
import sys

VARIANTS = ["ref", "c2-256x32", "c2-384x32", "c2-512x16", "ref-alusqr", "c2-256x32-alusqr"]


def parse_bench(text):
    """{variant: [(rate_B_per_s, sm_mhz, watts, temp_c), ...]} in repetition order."""
    out = {v: [] for v in VARIANTS}
    current = None
    for line in text.splitlines():
        m = re.match(r"=== bench (\S+) rep (\d+)", line)
        if m:
            current = m.group(1); continue
        if line.startswith("==="):
            current = None; continue
        m = re.search(r"finished: ([0-9.]+) M it/s", line)
        if m and current:
            out.setdefault(current, []).append([float(m.group(1)) / 1000.0, None, None, None])
        # The SM clock / power / temperature sample follows the rate, on the
        # same line (the job joins them) or the next.
        m = re.search(r"(\d+) MHz, ([0-9.]+) W, (\d+)\s*$", line)
        if m and current and out[current] and out[current][-1][1] is None:
            out[current][-1][1:] = [int(m.group(1)), float(m.group(2)), int(m.group(3))]
    return out


def parse_verify(path):
    text = path.read_text() if path.exists() else ""
    m = re.search(r"finished: [0-9.]+ M it/s, (\d+) distinguished points \((\d+) verified against the reference, (\d+) dropped\)", text)
    mism = len(re.findall(r"MISMATCH", text))
    reg = re.search(r"(\d+) registers/thread", text)
    res = re.search(r"(\d+) block\(s\) of (\d+) packed threads resident per SM", text)
    return {"distinguishedPoints": int(m.group(1)) if m else None,
            "verified": int(m.group(2)) if m else None,
            "dropped": int(m.group(3)) if m else None,
            "mismatches": mism,
            "registers": int(reg.group(1)) if reg else None,
            "residentBlocks": int(res.group(1)) if res else None,
            "blockThreads": int(res.group(2)) if res else None}


def parse_identity(path):
    out = {}
    if not path.exists(): return out
    for line in path.read_text().splitlines():
        m = re.match(r"(\S+)\s+(\d+) records\s+sha256 (\w+)\s+(IDENTICAL to \w+|DIFFERS from \w+)", line)
        if m: out[m.group(1)] = {"records": int(m.group(2)), "sha256_16": m.group(3), "identicalToRef": m.group(4).startswith("IDENTICAL")}
        m = re.match(r"(\S+)\s+missing", line)
        if m: out[m.group(1)] = {"missing": True}
    return out


def parse_profile(path):
    if not path.exists(): return None
    text = path.read_text()
    m = re.search(r"per warp-step (\d+) cycles = forward (\d+) \(([0-9.]+)%\) \+ inversion (\d+) \(([0-9.]+)%\) \+ reverse (\d+) \(([0-9.]+)%\); per update ([0-9.]+) warp-cycles", text)
    rate = re.search(r"finished: ([0-9.]+) M it/s", text)
    if not m: return {"raw": text.strip()}
    return {"cyclesPerWarpStep": int(m.group(1)), "forward": int(m.group(2)), "forwardPct": float(m.group(3)),
            "inversion": int(m.group(4)), "inversionPct": float(m.group(5)), "reverse": int(m.group(6)),
            "reversePct": float(m.group(7)), "warpCyclesPerUpdate": float(m.group(8)),
            "rateBps": float(rate.group(1)) / 1000.0 if rate else None}


def main():
    global VARIANTS
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("run"); ap.add_argument("target")
    ap.add_argument("--variants", default=",".join(VARIANTS), help="bench variant names, reference first")
    ap.add_argument("--profiles", default="ref-prof,c2-256x32-prof", help="PHASE_PROFILE binaries")
    a = ap.parse_args()
    VARIANTS = a.variants.split(",")
    run = pathlib.Path(a.run).resolve()
    target = pathlib.Path(a.target).resolve()
    results = run / "results"
    launch = json.loads((run / "launch.json").read_text())
    # job.log carries the "=== bench <variant> rep <n>" headers with the samples;
    # bench.txt is the tee of the samples alone.
    bench_src = results / "job.log" if (results / "job.log").exists() else results / "bench.txt"
    bench = parse_bench(bench_src.read_text()) if bench_src.exists() else {}
    identity = parse_identity(results / "dp-identity.txt")
    variants = {}
    reps = min((len(v) for v in bench.values() if v), default=0)
    for name in VARIANTS:
        samples = bench.get(name, [])
        rates = [s[0] for s in samples]
        entry = {"samplesBps": rates,
                 "medianBps": statistics.median(rates) if rates else None,
                 "smClockMHz": [s[1] for s in samples], "powerW": [s[2] for s in samples], "tempC": [s[3] for s in samples],
                 "verify": parse_verify(results / ("verify-%s.txt" % name)),
                 "dpIdentity": identity.get(name),
                 "build": (results / ("build-%s.txt" % name)).read_text().strip() if (results / ("build-%s.txt" % name)).exists() else None}
        if name != "ref" and rates and bench.get("ref"):
            pairs = [r / bench["ref"][i][0] for i, r in enumerate(rates) if i < len(bench["ref"])]
            entry["pairedRatioToRef"] = pairs
            entry["pairedRatioMedian"] = statistics.median(pairs) if pairs else None
            entry["pairedRatioMin"] = min(pairs) if pairs else None
            entry["pairedRatioMax"] = max(pairs) if pairs else None
        variants[name] = entry
    summary = {
        "note": "TWO-CHAINS.md",
        "unit": "billions of complete scalar updates per second, `finished:` line",
        "command": "--curve 131 --packed --bench --steps 1024 --launches 64 --verify 0, automatic worker count, alternating binaries",
        "verifyCommand": "--curve 131 --packed --threads <forced common count> --dp-weight 48 --dp-cap 262144 --steps 96 --launches 6 --verify 300 --run-id 7 --dp-file",
        "repetitions": reps,
        "host": (results / "host.txt").read_text().strip() if (results / "host.txt").exists() else None,
        "launch": {k: launch.get(k) for k in ("provider", "gpu", "cudaImage", "gitRev", "gitDirty", "startedAt", "finishedAt", "exitCode", "seconds", "deviceLine")},
        "floorBps": {"carrylessUnit33p1Clmad": 22.3, "carrylessUnit31p85Clmad": 23.1, "fiveProductsOnly": 24.6,
                     "basis": "1.62 lane-CLMADs per SM-clock, 188 SMs, 2.42 GHz (ONE-BLOCK-GEOMETRY.md section 1)"},
        "referenceBps": {"oneBlockGeometry": 20.078, "thisSession": variants["ref"]["medianBps"]},
        "variants": variants,
        "phaseProfile": {name: parse_profile(results / ("profile-%s.txt" % name)) for name in a.profiles.split(",")},
    }
    target.mkdir(parents=True, exist_ok=True)
    (target / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    raw = target / "raw"; raw.mkdir(exist_ok=True)
    for f in results.iterdir():
        if f.suffix in (".txt", ".log") and f.stat().st_size < 4 << 20:
            shutil.copy(f, raw / f.name)
    shutil.copy(run / "launch.json", raw / "launch.json")
    for name, entry in variants.items():
        print("%-20s median %s  paired %s  verified %s/%s dropped %s  DP %s" % (
            name, "%.3f" % entry["medianBps"] if entry["medianBps"] else "-",
            "%.4f" % entry["pairedRatioMedian"] if entry.get("pairedRatioMedian") else "-",
            entry["verify"]["verified"], entry["verify"]["distinguishedPoints"], entry["verify"]["dropped"],
            "identical" if (entry["dpIdentity"] or {}).get("identicalToRef") else entry["dpIdentity"]))
    for name, prof in summary["phaseProfile"].items():
        if prof: print("profile %-20s %s" % (name, {k: v for k, v in prof.items() if k in ("cyclesPerWarpStep", "forward", "inversion", "reverse", "warpCyclesPerUpdate", "rateBps")}))
    return 0


if __name__ == "__main__":
    sys.exit(main())

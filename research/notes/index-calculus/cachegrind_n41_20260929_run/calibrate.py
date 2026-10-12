#!/usr/bin/env python3
"""Host calibration for RESEARCH_CACHEGRIND_N41_20260929.md (run before the arms are measured).

Writes calibration.json: effective clock, load-to-use latency by working-set size,
memory-level parallelism, branch-mispredict cost. Nothing here touches the IC or rho arms.
usage: calibrate.py <path to mem_calib binary> <run index>
"""
import json
import statistics
import subprocess
import sys

BIN = sys.argv[1]
RUN = sys.argv[2]


def run(*args):
    out = subprocess.run([BIN, *map(str, args)], capture_output=True, text=True, check=True).stdout
    return json.loads(out.strip().splitlines()[-1])


KIB, MIB = 1024, 1024 * 1024
result = {"bin": BIN}
freqs = [run("freq")["f_hz"] for _ in range(3)]
result["f_hz_runs"] = freqs
result["f_hz"] = statistics.median(freqs)

sizes = [16 * KIB, 256 * KIB, 512 * KIB, MIB, 2 * MIB, 4 * MIB, 6 * MIB, 8 * MIB, 12 * MIB, 16 * MIB,
         32 * MIB, 64 * MIB, 256 * MIB, 1024 * MIB]
lat = {}
for w in sizes:
    steps = 20_000_000 if w <= MIB else 3_000_000 if w <= 64 * MIB else 1_500_000
    rows = [run("chase", w, 1, steps, 5, steps // 10) for _ in range(3)]
    lat[str(w)] = {
        "ns_per_load_median_of_runs": statistics.median(r["ns_per_step_median"] for r in rows),
        "runs": [r["ns_per_step_median"] for r in rows],
    }
result["latency_ns_by_bytes"] = lat

mlp = {}
base = None
for k in (1, 2, 4, 6, 8, 10, 12, 16):
    rows = [run("chase", 256 * MIB, k, 1_000_000, 5, 100_000) for _ in range(3)]
    ns_group = statistics.median(r["ns_per_step_median"] for r in rows)
    if k == 1:
        base = ns_group
    mlp[str(k)] = {"ns_per_group_step": ns_group, "ns_per_load": ns_group / k, "mlp": base / (ns_group / k)}
result["mlp_256MiB"] = mlp
result["mlp_max"] = max(v["mlp"] for v in mlp.values())

br = [run("branch") for _ in range(3)]
result["branch_ns_per_mispredict_runs"] = [b["ns_per_mispredict"] for b in br]
result["branch_ns_per_mispredict"] = statistics.median(result["branch_ns_per_mispredict_runs"])

json.dump(result, open(f"calibration_run{RUN}.json", "w"), indent=2)
print(f"run {RUN}:", {k: round(v["ns_per_load_median_of_runs"], 2) for k, v in lat.items()})
print(f"run {RUN} mlp:", {k: round(v["mlp"], 2) for k, v in mlp.items()}, "f_GHz", round(result["f_hz"] / 1e9, 3),
      "branch_ns", result["branch_ns_per_mispredict"])

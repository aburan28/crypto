#!/usr/bin/env python3
"""One cell of RESEARCH_PRECOMPUTED_START_RHO_20260929.md (macOS arm64).

usage: measure_cell.py <n> <L> <K> <corpus> [repeats=5]

setup   rho R3 once, uninstrumented; writes scalars.txt (the corpus' planted
        scalars, fed to IC). Not a measured repeat.
repeats IC, R3, R4 alternating, each under /usr/bin/time -l; the metric is
        "instructions retired" for the whole process. Before every process the
        1-minute load average must be below 14 (polled every 60 s).

Everything lands in cell_n<n>_L<L>_K<K>/ next to this script; results.json
holds per-run records, gates and medians. No verdict is written here.
"""
import json
import os
import re
import statistics
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
RHO = os.path.join(ROOT, "target", "release", "examples", "koblitz_rho_batch_ks_strong_ps")
IC = os.path.join(ROOT, "target", "release", "examples", "koblitz_orbit_dlp_fast")
BATCH_SEED = "531310"
LOAD_LIMIT = 14.0


def load1():
    out = subprocess.run(["sysctl", "-n", "vm.loadavg"], capture_output=True, text=True).stdout
    return float(out.strip("{} \n").split()[0])


def wait_quiet(log):
    while True:
        l1 = load1()
        if l1 < LOAD_LIMIT:
            return l1
        log.write(f"{time.strftime('%FT%TZ', time.gmtime())} wait load1={l1}\n")
        log.flush()
        time.sleep(60)


def rho_env(corpus, rung):
    env = dict(os.environ)
    env.update(KIC_RHO_RUNG=str(rung), KIC_RHO_LANES="32", KIC_RHO_DP_BITS="4",
               KIC_RHO_BATCH_CORPUS=corpus)
    return env


def timed(cmd, stdout_path, stderr_path, env=None):
    t0 = time.perf_counter()
    with open(stdout_path, "w") as out, open(stderr_path, "w") as err:
        rc = subprocess.run(["/usr/bin/time", "-l"] + cmd, stdout=out, stderr=err, env=env).returncode
    wall = time.perf_counter() - t0
    text = open(stderr_path).read()

    def grab(label):
        m = re.findall(r"^\s*(\d+)\s+" + re.escape(label) + r"\s*$", text, re.M)
        return int(m[-1]) if m else None

    rus = re.findall(r"([\d.]+) real\s+([\d.]+) user\s+([\d.]+) sys", text)
    return {
        "rc": rc,
        "wall_s": wall,
        "instructions_retired": grab("instructions retired"),
        "cycles_elapsed": grab("cycles elapsed"),
        "max_rss_bytes": grab("maximum resident set size"),
        "user_s": float(rus[-1][1]) if rus else None,
        "sys_s": float(rus[-1][2]) if rus else None,
    }


def parse_rho(path):
    fixtures, summary = [], None
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        obj = json.loads(line)
        if obj.get("kind") == "rho_ks_batch_fixture":
            fixtures.append(obj)
        elif obj.get("kind") == "rho_ks_batch_summary":
            summary = obj
    return fixtures, summary


def rho_run(n, L, corpus, rung, tag):
    rec = timed([RHO, str(n), "0", "signed_frobenius", str(L), BATCH_SEED],
                f"{tag}.jsonl", f"{tag}.stderr.log", rho_env(corpus, rung))
    fixtures, summary = parse_rho(f"{tag}.jsonl")
    rec["G1"] = (rec["rc"] == 0 and len(fixtures) == L and summary is not None
                 and summary.get("all_verified") is True
                 and all(f["verified"] and f["recovered_fixture_scalar"] == f["published_fixture_scalar"]
                         for f in fixtures))
    rec["total_walk_steps"] = summary and summary.get("total_walk_steps")
    rec["producer_version"] = summary and summary.get("producer_version")
    rec["precomputed_starts"] = summary and summary.get("precomputed_starts")
    return rec, fixtures


def ic_run(n, L, K, tag):
    rec = timed([IC, f"construct:{n}:0:{K}", "scalars.txt", "7", f"{tag}.jsonl"],
                f"{tag}.summary.json", f"{tag}.stderr.log")
    try:
        ic = json.load(open(f"{tag}.summary.json"))
    except Exception:  # noqa: BLE001 - reported as a failed gate
        ic = {}
    rec["G1"] = rec["rc"] == 0 and ic.get("targets_solved") == L and ic.get("targets_failed") == 0
    for k in ("targets_solved", "targets_failed", "rank", "orbit_columns", "rank_failures"):
        rec[k] = ic.get(k)
    return rec


def main():
    n, L, K, corpus = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    repeats = int(sys.argv[5]) if len(sys.argv) > 5 else 5
    cell = os.path.join(HERE, f"cell_n{n}_L{L}_K{K}")
    os.makedirs(cell, exist_ok=True)
    os.chdir(cell)
    log = open("measure.log", "a")
    res = {"n": n, "L": L, "K": K, "corpus": corpus, "batch_seed": BATCH_SEED, "repeats": repeats,
           "host": dict(zip(("sysname", "nodename", "release", "version", "machine"), os.uname())), "load_limit_1min": LOAD_LIMIT, "runs": []}

    res["setup_load1"] = wait_quiet(log)
    setup, fixtures = rho_run(n, L, corpus, 3, "setup_rho_r3")
    with open("scalars.txt", "w") as fh:
        for f in fixtures:
            fh.write(f"{f['published_fixture_scalar']}\n")
    res["setup"] = setup
    if not setup["G1"]:
        res["stopped"] = "setup R3 failed G1; no repeats run"
        json.dump(res, open("results.json", "w"), indent=1)
        return 1

    for r in range(repeats):
        for arm in ("ic", "r3", "r4"):
            l1 = wait_quiet(log)
            tag = f"rep{r}_{arm}"
            rec = ic_run(n, L, K, tag) if arm == "ic" else rho_run(n, L, corpus, 3 if arm == "r3" else 4, tag)[0]
            rec.update(arm=arm, repeat=r, load1_before=l1, load1_after=load1())
            res["runs"].append(rec)
            log.write(f"{time.strftime('%FT%TZ', time.gmtime())} {tag} rc={rec['rc']} "
                      f"ir={rec['instructions_retired']} G1={rec['G1']}\n")
            log.flush()
            json.dump(res, open("results.json", "w"), indent=1)

    med = {}
    for arm in ("ic", "r3", "r4"):
        rs = [x for x in res["runs"] if x["arm"] == arm]
        med[arm] = {
            "all_G1": all(x["G1"] for x in rs),
            "instructions_retired_median": statistics.median(x["instructions_retired"] for x in rs)
            if all(x["instructions_retired"] for x in rs) else None,
            "wall_s_median": statistics.median(x["wall_s"] for x in rs),
            "max_rss_bytes_median": statistics.median(x["max_rss_bytes"] or 0 for x in rs),
        }
        if arm != "ic":
            med[arm]["total_walk_steps_median"] = statistics.median(x["total_walk_steps"] or 0 for x in rs)
    s3, s4 = med["r3"]["total_walk_steps_median"], med["r4"]["total_walk_steps_median"]
    res["G3_r4_over_r3_steps"] = (s4 / s3) if s3 else None
    res["G3_pass"] = bool(s3) and 0.85 <= s4 / s3 <= 1.25
    res["medians"] = med
    json.dump(res, open("results.json", "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""The timing protocol of note section 15.3 (AGENTS.md section 10).

    python3 research/pkm_tower_round6_20260930/bench.py run BASELINE_BIN CANDIDATE_BIN
    python3 research/pkm_tower_round6_20260930/bench.py summary
    # the second measurement of note section 15.8: fresh runs, 15 pairs each
    python3 research/pkm_tower_round6_20260930/bench.py run BASELINE_BIN CANDIDATE_BIN --rounds 15 --dir bench2
    python3 research/pkm_tower_round6_20260930/bench.py summary --rounds 15 --dir bench2

`run` times three systems, target 0 each, at one thread
(`RAYON_NUM_THREADS=1`), every run through `tools/isolated_bench.py` on CPU 3:
- A/A: the baseline against a byte-identical copy of itself, interleaved over
  five rounds, for the noise floor;
- A/B: the baseline against the candidate, interleaved over five rounds.

Each run's isolated-bench record goes to `bench/records.jsonl`, its row to
`bench/rows/<label>.jsonl` and its stderr to `bench/logs/<label>.log`, where the
label is `<system>-<phase><round>-<build>`. A run the tool refuses (a busy
machine) is retried after a pause, and the refusal is kept in
`bench/refusals.txt`.

`summary` prints, per system: the median and minimum wall time of each build,
the A/A spread, the A/B ratio, the contended runs, and the multiply-adds of
each build, which must not vary between repeats or between the two copies of
the baseline. It also prints, per system and phase, the paired ratios of the
interleaved rounds (A2/A, B/A), their mean and its 95% confidence interval
(Student's t): the test of note section 15.8.
"""

import hashlib
import json
import os
import shutil
import statistics
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
OUT = os.path.join(HERE, "bench")
TOOL = os.path.join(ROOT, "tools", "isolated_bench.py")
P1 = "2013265921"
COMMON = ["--engine", "tower", "--p", P1, "--kinds", "kummer", "--controls", "tower",
          "--planted", "0", "--random", "1", "--ladder-t", "none", "--budget", "3600",
          "--max-nnz", "2000000000"]
SYSTEMS = [
    ("m4-N12", ["--m", "4", "--t-min", "3", "--max-t-m4", "3"]),
    ("m3-N12", ["--m", "3", "--t-min", "4", "--max-t-m3", "4"]),
    ("m2-N16", ["--m", "2", "--t-min", "8", "--max-t", "8"]),
]
ROUNDS = 5
CPU = "3"
# Two-sided 95% Student's t quantiles, by degrees of freedom.
T95 = {4: 2.776, 9: 2.262, 14: 2.145, 19: 2.093, 29: 2.045}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def options(argv):
    """`--rounds N` and `--dir NAME` (under this directory)."""
    global OUT, ROUNDS
    if "--rounds" in argv:
        ROUNDS = int(argv[argv.index("--rounds") + 1])
    if "--dir" in argv:
        OUT = os.path.join(HERE, argv[argv.index("--dir") + 1])


def paired(xs, ys):
    """Mean of ys[i] / xs[i] and its 95% confidence interval."""
    ratios = [y / x for x, y in zip(xs, ys)]
    m = statistics.mean(ratios)
    h = T95[len(ratios) - 1] * statistics.stdev(ratios) / len(ratios) ** 0.5
    return ratios, m, m - h, m + h


def one(label, binary, flags):
    rows = os.path.join(OUT, "rows", label + ".jsonl")
    log = os.path.join(OUT, "logs", label + ".log")
    env = dict(os.environ, RAYON_NUM_THREADS="1")
    cmd = [sys.executable, TOOL, "run", "--cpus", CPU, "--wait", "--label", label,
           "--out", os.path.join(OUT, "records.jsonl"), "--", binary, *COMMON, *flags, "--out", rows]
    for attempt in range(20):
        with open(log, "w") as err:
            code = subprocess.call(cmd, env=env, stdout=subprocess.DEVNULL, stderr=err)
        text = open(log).read()
        if "machine is busy" not in text:
            if code != 0:
                raise SystemExit(f"{label}: exit {code}; see {log}")
            return
        with open(os.path.join(OUT, "refusals.txt"), "a") as f:
            f.write(f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {label} attempt {attempt}: "
                    f"{text.strip().splitlines()[-1]}\n")
        time.sleep(30)
    raise SystemExit(f"{label}: the machine stayed busy")


def run(base, cand):
    for d in ("rows", "logs"):
        os.makedirs(os.path.join(OUT, d), exist_ok=True)
    tmp = tempfile.mkdtemp()
    copy = os.path.join(tmp, "baseline-copy")
    shutil.copy2(base, copy)
    builds = {"A": base, "A2": copy, "B": cand}
    with open(os.path.join(OUT, "builds.json"), "w") as f:
        json.dump({k: {"path": v, "sha256": sha256(v)} for k, v in builds.items()}, f, indent=1)
    for system, flags in SYSTEMS:
        for r in range(ROUNDS):
            for b in ("A", "A2"):
                one(f"{system}-aa{r}-{b}", builds[b], flags)
        for r in range(ROUNDS):
            for b in ("A", "B"):
                one(f"{system}-ab{r}-{b}", builds[b], flags)
    shutil.rmtree(tmp)


def summary():
    recs = {}
    with open(os.path.join(OUT, "records.jsonl")) as f:
        for line in f:
            r = json.loads(line)
            recs[r["label"]] = r
    print("| system | phase | build | runs | contended | median s | min s | multiply-adds |")
    print("|:--|:--|:--|--:|--:|--:|--:|--:|")
    bad = 0
    ratios = []
    for system, _ in SYSTEMS:
        med = {}
        for phase, builds in (("aa", ("A", "A2")), ("ab", ("A", "B"))):
            for b in builds:
                labels = [f"{system}-{phase}{r}-{b}" for r in range(ROUNDS)]
                runs = [recs[l]["run"] for l in labels if l in recs]
                clean = [x["wall_seconds"] for x in runs if not x["contended"]]
                contended = sum(x["contended"] for x in runs)
                mul = set()
                for l in labels:
                    path = os.path.join(OUT, "rows", l + ".jsonl")
                    if os.path.exists(path):
                        for line in open(path):
                            row = json.loads(line)
                            if "N" in row:
                                mul.add(row["muladds"])
                if len(mul) > 1:
                    bad += 1
                med[(phase, b)] = (statistics.median(clean) if clean else None,
                                   min(clean) if clean else None, mul)
                m, lo, _ = med[(phase, b)]
                print(f"| {system} | {phase.upper()} | {b} | {len(runs)} | {contended} | "
                      f"{m if m is None else f'{m:.3f}'} | {lo if lo is None else f'{lo:.3f}'} | "
                      f"{', '.join(str(x) for x in sorted(mul))} |")
        a, a2 = med[("aa", "A")], med[("aa", "A2")]
        b0, b = med[("ab", "A")], med[("ab", "B")]
        if a[2] != a2[2]:
            bad += 1
        if a[0] and a2[0] and b0[0] and b[0]:
            spread = abs(a2[0] - a[0]) / a[0]
            ratios.append(f"- {system}: A/A spread of the medians {100 * spread:.2f}%; "
                          f"A/B median ratio B/A {b[0] / b0[0]:.4f}, minimum ratio "
                          f"{b[1] / b0[1]:.4f}; multiply-adds B/A {min(b[2]) / min(b0[2]):.4f}")
    print()
    print("\n".join(ratios))
    print()
    for system, _ in SYSTEMS:
        for phase, (x, y) in (("aa", ("A", "A2")), ("ab", ("A", "B"))):
            pairs = [(f"{system}-{phase}{r}-{x}", f"{system}-{phase}{r}-{y}") for r in range(ROUNDS)]
            pairs = [(a, b) for a, b in pairs if a in recs and b in recs
                     and not recs[a]["run"]["contended"] and not recs[b]["run"]["contended"]]
            if len(pairs) < 2:
                continue
            _, m, lo, hi = paired([recs[a]["run"]["wall_seconds"] for a, _ in pairs],
                                  [recs[b]["run"]["wall_seconds"] for _, b in pairs])
            verdict = "slower" if lo > 1 else "faster" if hi < 1 else "no difference resolved"
            print(f"- {system} {phase.upper()} ({y}/{x}, {len(pairs)} uncontended pairs): mean {m:.4f}, "
                  f"95% interval [{lo:.4f}, {hi:.4f}]: {verdict}")
    if bad:
        print(f"{bad} build(s) whose multiply-adds varied between runs or copies")
    return 1 if bad else 0


if __name__ == "__main__":
    options(sys.argv)
    if sys.argv[1:2] == ["run"]:
        run(sys.argv[2], sys.argv[3])
    elif sys.argv[1:2] == ["summary"]:
        sys.exit(summary())
    else:
        raise SystemExit(__doc__)

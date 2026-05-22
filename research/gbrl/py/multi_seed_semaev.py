# Multi-seed PPO on Semaev systems (arith_ops reward), evaluating each seed
# greedily across train + held-out instances. Reports per-instance mean ± std
# and whether the policies strictly beat sugar.
#
#     python3 py/multi_seed_semaev.py --seeds 5

import argparse
import json
import os
import re
import subprocess
import time

import numpy as np


def parse_per_instance_best(log_path):
    """Walk the training log and track the best untruncated greedy seen for
    each instance — best by arith_ops, ties broken by steps."""
    pattern = re.compile(
        r"\[eval\]\s+OK\s+(\S+):\s+steps=(\d+)\s+basis=(\d+)\s+arith_ops=(\d+)"
    )
    best = {}  # inst -> (steps, basis, ops)
    if not os.path.exists(log_path):
        return best
    with open(log_path) as f:
        for line in f:
            m = pattern.search(line)
            if not m:
                continue
            inst = m.group(1)
            steps = int(m.group(2))
            basis = int(m.group(3))
            ops = int(m.group(4))
            cur = best.get(inst)
            if cur is None or ops < cur[2] or (ops == cur[2] and steps < cur[0]):
                best[inst] = (steps, basis, ops)
    return best


def train_one_seed(args, seed):
    save = f"/tmp/semaev_ms_seed{seed}.pt"
    log = f"/tmp/semaev_ms_seed{seed}.log"
    cmd = [
        "python3", "py/train.py",
        "--instances", *args.train_instances,
        "--seed", str(seed),
        "--updates", str(args.updates),
        "--batch-episodes", str(args.batch_episodes),
        "--parallel", str(args.parallel),
        "--max-steps", str(args.max_steps),
        "--step-timeout", "5",
        "--eval-max-wall", "30",
        "--lr", str(args.lr),
        "--entropy", str(args.entropy),
        "--reward-mode", args.reward_mode,
        "--init", args.init,
        "--save", save,
    ]
    t0 = time.time()
    with open(log, "w") as f:
        subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT,
                       env={**os.environ, "PYTHONUNBUFFERED": "1"})
    wall = time.time() - t0
    best = parse_per_instance_best(log)
    return {"seed": seed, "wall_s": wall, "log": log, "save": save, "best": best}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train-instances", nargs="+",
                    default=["semaev2-4", "semaev2-5", "semaev2-6",
                             "semaev2-7", "semaev2-8"])
    ap.add_argument("--init", default="/tmp/semaev_il.pt")
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--updates", type=int, default=30)
    ap.add_argument("--batch-episodes", type=int, default=32)
    ap.add_argument("--parallel", type=int, default=4)
    ap.add_argument("--max-steps", type=int, default=200)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--entropy", type=float, default=0.01)
    ap.add_argument("--reward-mode", default="arith_ops")
    args = ap.parse_args()

    # Sugar / normal baselines (hard-coded from earlier baseline.py output).
    sugar = {
        "semaev2-4":  (23, 130),
        "semaev2-5":  (34, 243),
        "semaev2-6":  (43, 360),
        "semaev2-7":  (51, 478),
        "semaev2-8":  (59, 612),
        "semaev2-9":  (64, 768),
        "semaev2-10": (72, 959),
        "semaev2-11": (81, 1196),
        "semaev2-12": (90, 1471),
    }

    runs = []
    for s in range(args.seeds):
        print(f"\n=== seed {s} ===", flush=True)
        r = train_one_seed(args, s)
        runs.append(r)
        print(f"  wall {r['wall_s']:.0f}s; best per instance:", flush=True)
        for inst in sorted(r["best"].keys()):
            st, ba, op = r["best"][inst]
            su = sugar.get(inst, (None, None))
            tag = ""
            if su[0] is not None:
                if st < su[0] and op < su[1]:
                    tag = "  ★ beats sugar BOTH"
                elif st <= su[0] and op <= su[1] and (st < su[0] or op < su[1]):
                    tag = "  ☆ beats sugar BOTH (ties)"
            print(f"    {inst:>12}: steps={st:4d}  arith_ops={op:>6}{tag}",
                  flush=True)

    print("\n=== multi-seed summary (best per (seed, instance)) ===\n")
    all_instances = sorted(sugar.keys())
    print(f"{'instance':>12} {'sugar_s':>8} {'sugar_op':>9}  "
          f"{'mean_s':>7} {'std_s':>5} {'min_s':>5}  "
          f"{'mean_op':>8} {'std_op':>6} {'min_op':>6}  "
          f"{'beat_both':>10}")
    for inst in all_instances:
        per_seed = [r["best"].get(inst) for r in runs if r["best"].get(inst)]
        if not per_seed:
            continue
        steps_arr = np.array([s for s, _, _ in per_seed])
        ops_arr = np.array([o for _, _, o in per_seed])
        sug_s, sug_op = sugar[inst]
        n_beat = int(((steps_arr < sug_s) & (ops_arr < sug_op)).sum())
        print(f"{inst:>12} {sug_s:>8} {sug_op:>9}  "
              f"{steps_arr.mean():>7.1f} {steps_arr.std():>5.1f} {steps_arr.min():>5}  "
              f"{ops_arr.mean():>8.0f} {ops_arr.std():>6.0f} {ops_arr.min():>6}  "
              f"{n_beat:>4}/{len(per_seed)}")
    print()

    with open("/tmp/semaev_multiseed.json", "w") as f:
        json.dump(runs, f, indent=2, default=str)
    print("summary dumped to /tmp/semaev_multiseed.json")


if __name__ == "__main__":
    main()

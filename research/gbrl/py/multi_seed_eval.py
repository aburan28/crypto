# Train K seeds with the best hyperparameters and report mean ± std of greedy
# (steps, arith_ops) at the end. Decides whether PPO actually beats sugar.
#
#     python3 py/multi_seed_eval.py --seeds 5 --updates 30 --lr 3e-4 --entropy 0.01

import argparse
import json
import os
import re
import subprocess
import time

import numpy as np


def parse_best(log_path):
    ops = steps = None
    if not os.path.exists(log_path):
        return ops, steps
    with open(log_path) as f:
        for line in f:
            if "best untruncated greedy" in line:
                m = re.search(r"arith_ops=(\d+)", line)
                if m:
                    ops = int(m.group(1))
                m = re.search(r"steps=(\d+)", line)
                if m:
                    steps = int(m.group(1))
    return ops, steps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--instance", default="cyclic-6")
    ap.add_argument("--init", default="/tmp/cyc56_il.pt")
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--updates", type=int, default=30)
    ap.add_argument("--batch-episodes", type=int, default=32)
    ap.add_argument("--parallel", type=int, default=4)
    ap.add_argument("--max-steps", type=int, default=800)
    ap.add_argument("--step-timeout", type=float, default=5.0)
    ap.add_argument("--eval-max-wall", type=float, default=120.0)
    ap.add_argument("--lr", type=float, default=3e-4)
    ap.add_argument("--entropy", type=float, default=0.01)
    args = ap.parse_args()

    runs = []
    for seed in range(args.seeds):
        save = f"/tmp/cyc6_ms_seed{seed}.pt"
        log = f"/tmp/cyc6_ms_seed{seed}.log"
        cmd = [
            "python3", "py/train.py",
            "--instance", args.instance,
            "--seed", str(seed),
            "--updates", str(args.updates),
            "--batch-episodes", str(args.batch_episodes),
            "--parallel", str(args.parallel),
            "--max-steps", str(args.max_steps),
            "--step-timeout", str(args.step_timeout),
            "--eval-max-wall", str(args.eval_max_wall),
            "--lr", str(args.lr),
            "--entropy", str(args.entropy),
            "--init", args.init,
            "--save", save,
        ]
        print(f"\n=== seed {seed} ===", flush=True)
        t0 = time.time()
        with open(log, "w") as f:
            subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT,
                            env={**os.environ, "PYTHONUNBUFFERED": "1"})
        wall = time.time() - t0
        ops, steps = parse_best(log)
        print(f"  best: steps={steps}  arith_ops={ops}  wall={wall:.0f}s",
              flush=True)
        runs.append({"seed": seed, "best_steps": steps, "best_ops": ops,
                     "wall_s": wall, "best_path": save + ".best", "log": log})

    print("\n=== multi-seed summary ===")
    valid = [r for r in runs if r["best_ops"] is not None]
    if not valid:
        print("ALL SEEDS TRUNCATED — PPO failed to find an untruncated greedy "
              "policy in this budget.")
        print("sugar baseline (cyclic-6): steps=379  arith_ops=14864")
        return
    ops = np.array([r["best_ops"] for r in valid])
    steps = np.array([r["best_steps"] for r in valid])
    print(f"valid seeds: {len(valid)}/{args.seeds}")
    print(f"  steps     : mean={steps.mean():7.1f}  std={steps.std():6.1f}  "
          f"min={steps.min()}  max={steps.max()}")
    print(f"  arith_ops : mean={ops.mean():7.0f}  std={ops.std():6.0f}  "
          f"min={ops.min()}  max={ops.max()}")
    print("\nbaselines (cyclic-6, post-GM):")
    print("  sugar  : steps=379  arith_ops=14864")
    print("  normal : steps=672  arith_ops=25162")
    n_beat_both = int(((steps < 379) & (ops < 14864)).sum())
    n_beat_steps = int((steps < 379).sum())
    n_beat_ops = int((ops < 14864).sum())
    print(f"\nseeds beating sugar by steps      : {n_beat_steps}/{len(valid)}")
    print(f"seeds beating sugar by arith_ops  : {n_beat_ops}/{len(valid)}")
    print(f"seeds beating sugar by BOTH       : {n_beat_both}/{len(valid)}")

    summary = {"runs": runs,
               "steps_mean": float(steps.mean()), "steps_std": float(steps.std()),
               "ops_mean": float(ops.mean()), "ops_std": float(ops.std()),
               "n_valid": len(valid), "n_seeds": args.seeds,
               "n_beat_both": n_beat_both}
    with open("/tmp/cyc6_multiseed.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nsummary dumped to /tmp/cyc6_multiseed.json")


if __name__ == "__main__":
    main()

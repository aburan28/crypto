# Small hyperparameter sweep for cyclic-6 PPO. For each config, run a short
# training and record the best untruncated greedy (steps, arith_ops).
#
#     python3 py/sweep_cyc6.py --updates 12 --init /tmp/cyc4.pt
#
# Configs to try (per the user's request):
#   lr      ∈ {1e-4, 3e-4, 1e-3}
#   entropy ∈ {0.00, 0.01, 0.03}
#   parallel ∈ {4, 8}
#   batch_episodes ∈ {32, 64}
# Full grid is 36 configs — too expensive. The default sweep below probes the
# corners and centers; extend manually if a corner looks promising.

import argparse
import json
import os
import re
import subprocess
import time


CONFIGS = [
    # (lr, entropy, parallel, batch_episodes)
    (3e-4, 0.01, 4, 32),
    (3e-4, 0.00, 4, 32),
    (3e-4, 0.03, 4, 32),
    (1e-4, 0.01, 4, 32),
    (1e-3, 0.01, 4, 32),
    (3e-4, 0.01, 8, 64),
]


def parse_best(log_path):
    """Pull the 'best untruncated greedy' line written at the end of training."""
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


def run_one(args, lr, ec, par, be):
    label = f"lr{lr:.0e}_ec{ec}_p{par}_b{be}".replace("+0", "").replace("-0", "-")
    save = f"/tmp/cyc6_sweep_{label}.pt"
    log = f"/tmp/cyc6_sweep_{label}.log"
    cmd = [
        "python3", "py/train.py",
        "--instance", args.instance,
        "--updates", str(args.updates),
        "--batch-episodes", str(be),
        "--parallel", str(par),
        "--max-steps", str(args.max_steps),
        "--step-timeout", str(args.step_timeout),
        "--eval-max-wall", str(args.eval_max_wall),
        "--lr", str(lr),
        "--entropy", str(ec),
        "--init", args.init,
        "--save", save,
        "--seed", str(args.seed),
    ]
    print(f"\n=== {label} ===  ({' '.join(cmd[2:])})", flush=True)
    t0 = time.time()
    with open(log, "w") as f:
        subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, env={**os.environ, "PYTHONUNBUFFERED": "1"})
    wall = time.time() - t0
    ops, steps = parse_best(log)
    print(f"  best: steps={steps}  arith_ops={ops}  wall={wall:.0f}s   log={log}",
          flush=True)
    return {"label": label, "lr": lr, "entropy": ec, "parallel": par,
            "batch_episodes": be, "best_steps": steps, "best_ops": ops,
            "wall_s": wall, "log": log, "save": save}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--instance", default="cyclic-6")
    ap.add_argument("--init", default="/tmp/cyc56_il.pt")
    ap.add_argument("--updates", type=int, default=12)
    ap.add_argument("--max-steps", type=int, default=800)
    ap.add_argument("--step-timeout", type=float, default=5.0)
    ap.add_argument("--eval-max-wall", type=float, default=120.0)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    results = []
    for lr, ec, par, be in CONFIGS:
        results.append(run_one(args, lr, ec, par, be))

    print("\n=== sweep summary (sorted by arith_ops, untruncated only) ===")
    valid = [r for r in results if r["best_ops"] is not None]
    valid.sort(key=lambda r: r["best_ops"])
    print(f"{'label':<35} {'steps':>6} {'arith_ops':>10} {'wall_s':>7}")
    for r in valid:
        print(f"{r['label']:<35} {r['best_steps']:>6} {r['best_ops']:>10} {r['wall_s']:>7.0f}")
    truncated = [r for r in results if r["best_ops"] is None]
    if truncated:
        print(f"\nNo untruncated greedy in {len(truncated)} configs: "
              f"{[r['label'] for r in truncated]}")

    # Also dump JSON for downstream multi-seed eval.
    out_json = "/tmp/cyc6_sweep_results.json"
    with open(out_json, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nresults dumped to {out_json}")
    print(f"sugar baseline (cyclic-6): steps=379  arith_ops=14864")


if __name__ == "__main__":
    main()

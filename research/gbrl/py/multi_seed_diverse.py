# Multi-seed PPO + per-instance eval across train / held-out splits of the
# diverse Semaev family. Reports whether the learned policies beat sugar
# on instances they did and did NOT train on.
#
#     python3 py/multi_seed_diverse.py --seeds 5

import argparse
import json
import os
import re
import subprocess
import time

import numpy as np
import torch

from env import GbrlEnv, featurize
from model import PointerPolicy


def train_one_seed(args, seed, train_instances):
    save = f"/tmp/semaev_d_ms_seed{seed}.pt"
    log = f"/tmp/semaev_d_ms_seed{seed}.log"
    cmd = [
        "python3", "py/train.py",
        "--instances", *train_instances,
        "--seed", str(seed),
        "--updates", str(args.updates),
        "--batch-episodes", str(args.batch_episodes),
        "--parallel", str(args.parallel),
        "--max-steps", str(args.max_steps),
        "--step-timeout", "5",
        "--eval-max-wall", "30",
        "--eval-every", "9999",  # skip per-update eval; we eval separately
        "--lr", str(args.lr),
        "--entropy", str(args.entropy),
        "--reward-mode", "arith_ops",
        "--init", args.init,
        "--save", save,
    ]
    t0 = time.time()
    with open(log, "w") as f:
        subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT,
                       env={**os.environ, "PYTHONUNBUFFERED": "1"})
    return {"seed": seed, "wall_s": time.time() - t0, "log": log, "save": save}


def greedy(env, policy, instance, max_steps=2000, max_wall=30.0):
    import time as _t
    obs = env.reset(instance)
    steps = 0
    t0 = _t.time()
    while True:
        f = featurize(obs)
        if f is None:
            break
        if _t.time() - t0 > max_wall:
            return None
        pf = torch.from_numpy(f[0]); gf = torch.from_numpy(f[1])
        with torch.no_grad():
            logits, _ = policy(pf, gf)
            a = int(torch.argmax(logits).item())
        r = env.step(a)
        steps += 1
        obs = r.obs
        if r.done or steps >= max_steps:
            break
    return steps, obs.get("basis_size", 0), obs.get("arith_ops", 0)


def get_sugar(env, instance):
    """Sugar greedy: at each state pick min (sugar, lcm_degree)."""
    obs = env.reset(instance)
    steps = 0
    while True:
        pairs = obs["pairs"]
        if not pairs:
            break
        idx = min(range(len(pairs)),
                  key=lambda k: (pairs[k]["sugar"], pairs[k]["lcm_degree"]))
        obs = env.step(idx).obs
        steps += 1
        if obs["done"]:
            break
    return steps, obs.get("basis_size", 0), obs.get("arith_ops", 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train-seeds", type=int, default=50)
    ap.add_argument("--test-seeds", type=int, default=50)
    ap.add_argument("--seeds", type=int, default=5, help="training random seeds")
    ap.add_argument("--init", default="/tmp/semaev_d_il.pt")
    ap.add_argument("--updates", type=int, default=30)
    ap.add_argument("--batch-episodes", type=int, default=32)
    ap.add_argument("--parallel", type=int, default=4)
    ap.add_argument("--max-steps", type=int, default=400)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--entropy", type=float, default=0.01)
    args = ap.parse_args()

    train_inst = [f"semaev2d-{i}" for i in range(args.train_seeds)]
    test_inst = [f"semaev2d-{i}" for i in range(args.train_seeds,
                                                 args.train_seeds + args.test_seeds)]

    # Train all seeds first.
    runs = []
    for s in range(args.seeds):
        print(f"\n=== seed {s}: training ===", flush=True)
        r = train_one_seed(args, s, train_inst)
        runs.append(r)
        print(f"  wall {r['wall_s']:.0f}s", flush=True)

    # Eval each policy on train + test, plus sugar baseline.
    binary = os.path.abspath(os.path.join(os.path.dirname(__file__), "..",
                                           "target", "release", "rl_server"))
    env = GbrlEnv(binary)
    try:
        print("\nrunning sugar baselines...", flush=True)
        sugar = {}
        for inst in train_inst + test_inst:
            sugar[inst] = get_sugar(env, inst)

        per_seed_results = []
        for r in runs:
            policy = PointerPolicy()
            policy.load_state_dict(torch.load(r["save"], map_location="cpu"))
            policy.eval()
            inst_results = {}
            for inst in train_inst + test_inst:
                inst_results[inst] = greedy(env, policy, inst)
            per_seed_results.append({"seed": r["seed"], "results": inst_results})
            print(f"  seed {r['seed']}: evaluated {len(inst_results)} instances",
                  flush=True)

        # Aggregate per-instance stats, separated train/test.
        def report(label, instances):
            print(f"\n=== {label} ({len(instances)} instances) ===")
            beat_both_total = 0
            beat_steps_total = 0
            beat_ops_total = 0
            n_total = 0
            sugar_steps_sum = 0
            sugar_ops_sum = 0
            ppo_steps_sum = 0
            ppo_ops_sum = 0
            for inst in instances:
                sug_s, _, sug_op = sugar[inst]
                # Best PPO across seeds (lowest arith_ops, ties broken by steps)
                vals = [psr["results"][inst] for psr in per_seed_results
                        if psr["results"][inst]]
                if not vals:
                    continue
                vals.sort(key=lambda x: (x[2], x[0]))
                best_s, _, best_op = vals[0]
                n_total += 1
                sugar_steps_sum += sug_s
                sugar_ops_sum += sug_op
                ppo_steps_sum += best_s
                ppo_ops_sum += best_op
                if best_s < sug_s and best_op < sug_op: beat_both_total += 1
                if best_s < sug_s: beat_steps_total += 1
                if best_op < sug_op: beat_ops_total += 1
            print(f"  PPO best-of-{args.seeds}-seed vs sugar (per-instance):")
            print(f"    beat sugar on BOTH metrics:     {beat_both_total}/{n_total}")
            print(f"    beat sugar on STEPS:            {beat_steps_total}/{n_total}")
            print(f"    beat sugar on ARITH_OPS:        {beat_ops_total}/{n_total}")
            if n_total:
                print(f"  total step count : sugar {sugar_steps_sum}  PPO {ppo_steps_sum}  "
                      f"({100*(sugar_steps_sum-ppo_steps_sum)/sugar_steps_sum:+.1f}%)")
                print(f"  total arith_ops  : sugar {sugar_ops_sum}  PPO {ppo_ops_sum}  "
                      f"({100*(sugar_ops_sum-ppo_ops_sum)/sugar_ops_sum:+.1f}%)")

        report("TRAIN", train_inst)
        report("TEST (held-out)", test_inst)

        out = {"sugar": {k: v for k, v in sugar.items()},
               "per_seed": per_seed_results}
        with open("/tmp/semaev_diverse_eval.json", "w") as f:
            json.dump(out, f, indent=2, default=str)
        print("\ndumped to /tmp/semaev_diverse_eval.json", flush=True)
    finally:
        env.close()


if __name__ == "__main__":
    main()

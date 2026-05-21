# Side-by-side comparison: built-in heuristics vs trained policies, across
# multiple instances. Drives the final results table.
#
#     python3 compare.py --weights label1=path1 label2=path2 \
#                        --instances cyclic-3 cyclic-4 cyclic-5

import argparse
import os

import torch

from baseline import run as run_baseline
from env import GbrlEnv, featurize
from model import PointerPolicy


def greedy(env, policy, instance, max_steps=50_000):
    obs = env.reset(instance)
    steps = 0
    while True:
        f = featurize(obs)
        if f is None:
            break
        pf = torch.from_numpy(f[0])
        gf = torch.from_numpy(f[1])
        with torch.no_grad():
            logits, _ = policy(pf, gf)
            a = int(torch.argmax(logits).item())
        r = env.step(a)
        steps += 1
        obs = r.obs
        if r.done or steps >= max_steps:
            break
    return steps, obs["basis_size"], obs["arith_ops"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", default=os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "target", "release", "rl_server")))
    ap.add_argument("--instances", nargs="+", required=True)
    ap.add_argument("--weights", nargs="+", default=[],
                    help='label=path entries, e.g. ppo=/tmp/cyc4.pt')
    args = ap.parse_args()

    labeled_policies = []
    for spec in args.weights:
        if "=" not in spec:
            raise SystemExit(f"--weights entries must be label=path, got: {spec}")
        label, path = spec.split("=", 1)
        p = PointerPolicy()
        p.load_state_dict(torch.load(path, map_location="cpu"))
        p.eval()
        labeled_policies.append((label, p))

    print(f"{'strategy':>14} {'instance':>10} {'steps':>7} {'basis':>6} {'arith_ops':>10}")
    env = GbrlEnv(args.binary)
    try:
        for inst in args.instances:
            for strat in ("sugar", "normal"):
                steps, basis, ops = run_baseline(args.binary, inst, strat)
                print(f"{strat:>14} {inst:>10} {steps:>7} {basis:>6} {ops:>10}")
            for label, policy in labeled_policies:
                steps, basis, ops = greedy(env, policy, inst)
                print(f"{label:>14} {inst:>10} {steps:>7} {basis:>6} {ops:>10}")
            print()
    finally:
        env.close()


if __name__ == "__main__":
    main()

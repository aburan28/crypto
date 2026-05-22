# Evaluate a saved policy greedily across instances, alongside sugar/normal/random.
#
#     python3 eval.py --weights /tmp/cyc4.pt --instances cyclic-4 cyclic-5

import argparse
import os

import numpy as np
import torch

from baseline import run as run_baseline
from env import GbrlEnv, featurize
from model import PointerPolicy


def greedy(env: GbrlEnv, policy: PointerPolicy, instance: str, max_steps: int = 50_000):
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
    ap.add_argument("--weights", required=True)
    ap.add_argument("--instances", nargs="+", default=["cyclic-4", "cyclic-5"])
    args = ap.parse_args()

    policy = PointerPolicy()
    policy.load_state_dict(torch.load(args.weights, map_location="cpu"))
    policy.eval()

    print(f"{'strategy':>10} {'instance':>10} {'steps':>7} {'basis':>6} {'arith_ops':>10}")
    env = GbrlEnv(args.binary)
    try:
        for inst in args.instances:
            for strat in ("sugar", "normal"):
                steps, basis, ops = run_baseline(args.binary, inst, strat)
                print(f"{strat:>10} {inst:>10} {steps:>7} {basis:>6} {ops:>10}")
            steps, basis, ops = greedy(env, policy, inst)
            print(f"{'learned':>10} {inst:>10} {steps:>7} {basis:>6} {ops:>10}")
    finally:
        env.close()


if __name__ == "__main__":
    main()

# Run greedy rollouts on N random inits — establishes whether "greedy=35"
# in training was actually learned, or just luck from the small logit init.

import os
import sys

import numpy as np
import torch

from env import GbrlEnv, featurize
from model import PointerPolicy


def greedy(env, policy, instance, max_steps=10_000):
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
    return steps


def main():
    binary = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "target", "release", "rl_server"))
    instance = sys.argv[1] if len(sys.argv) > 1 else "cyclic-4"
    env = GbrlEnv(binary)
    try:
        results = []
        for seed in range(20):
            torch.manual_seed(seed)
            policy = PointerPolicy()
            policy.eval()
            s = greedy(env, policy, instance)
            results.append(s)
            print(f"seed={seed:2d}  greedy {instance}: steps={s}")
        results = np.array(results)
        print(f"\nrandom-init greedy {instance}: "
              f"mean={results.mean():.1f}  min={results.min()}  "
              f"max={results.max()}  median={int(np.median(results))}")
    finally:
        env.close()


if __name__ == "__main__":
    main()

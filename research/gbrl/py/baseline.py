# Run built-in strategies through the same JSON-RPC env as the trainer.
# Use these numbers as the target for the learned policy to beat.
#
#     python3 baseline.py --instances cyclic-3 cyclic-4 cyclic-5

import argparse
import json
import os
import subprocess

import numpy as np


def _spawn(binary: str) -> subprocess.Popen:
    return subprocess.Popen(
        [binary],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    )


def _call(proc: subprocess.Popen, req: dict) -> dict:
    proc.stdin.write(json.dumps(req) + "\n")
    proc.stdin.flush()
    return json.loads(proc.stdout.readline())


def run(binary: str, instance: str, strategy: str, seed: int = 0, max_steps: int = 200_000):
    proc = _spawn(binary)
    rng = np.random.default_rng(seed)
    try:
        obs = _call(proc, {"method": "Reset", "instance": instance})["obs"]
        while not obs["done"]:
            pairs = obs["pairs"]
            if not pairs:
                break
            if strategy == "sugar":
                idx = min(range(len(pairs)), key=lambda k: (pairs[k]["sugar"], pairs[k]["lcm_degree"]))
            elif strategy == "normal":
                idx = min(range(len(pairs)), key=lambda k: pairs[k]["lcm_degree"])
            elif strategy == "random":
                idx = int(rng.integers(0, len(pairs)))
            else:
                raise ValueError(f"unknown strategy: {strategy}")
            obs = _call(proc, {"method": "Step", "pair_idx": idx})["obs"]
            if obs["step_count"] >= max_steps:
                break
        return obs["step_count"], obs["basis_size"], obs["arith_ops"]
    finally:
        try:
            proc.stdin.close()
            proc.wait(timeout=2)
        except Exception:
            proc.kill()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", default=os.path.join(os.path.dirname(__file__), "..", "target", "release", "rl_server"))
    ap.add_argument("--instances", nargs="+", default=["cyclic-3", "cyclic-4"])
    ap.add_argument("--strategies", nargs="+", default=["sugar", "normal", "random"])
    args = ap.parse_args()

    binary = os.path.abspath(args.binary)
    print(f"{'strategy':>9} {'instance':>10} {'steps':>7} {'basis':>6} {'arith_ops':>10}")
    for inst in args.instances:
        for strat in args.strategies:
            steps, basis, ops = run(binary, inst, strat)
            print(f"{strat:>9} {inst:>10} {steps:>7} {basis:>6} {ops:>10}")

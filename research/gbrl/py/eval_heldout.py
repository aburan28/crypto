# Greedy eval each multi-seed checkpoint on held-out Semaev instances and
# compare to sugar.

import argparse
import os

import numpy as np
import torch

from env import GbrlEnv, featurize
from model import PointerPolicy


def greedy(env, policy, instance, max_steps=2000, max_wall=60.0):
    import time as _t
    obs = env.reset(instance)
    steps = 0
    t0 = _t.time()
    while True:
        f = featurize(obs)
        if f is None:
            break
        if _t.time() - t0 > max_wall:
            return None  # truncated
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
    return steps, obs.get("basis_size", 0), obs.get("arith_ops", 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", default=os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "target", "release", "rl_server")))
    ap.add_argument("--instances", nargs="+",
                    default=["semaev2-9", "semaev2-10", "semaev2-11", "semaev2-12"])
    ap.add_argument("--checkpoints", nargs="+", required=True,
                    help="paths to PointerPolicy state dicts to eval")
    args = ap.parse_args()

    sugar = {
        "semaev2-9":  (64, 768),
        "semaev2-10": (72, 959),
        "semaev2-11": (81, 1196),
        "semaev2-12": (90, 1471),
    }

    env = GbrlEnv(args.binary)
    try:
        # Greedy each checkpoint × each instance.
        results = {}  # (ckpt, inst) -> (steps, basis, ops)
        for ckpt in args.checkpoints:
            policy = PointerPolicy()
            policy.load_state_dict(torch.load(ckpt, map_location="cpu"))
            policy.eval()
            for inst in args.instances:
                out = greedy(env, policy, inst)
                results[(ckpt, inst)] = out
                if out:
                    s, b, o = out
                    sug = sugar.get(inst, (None, None))
                    tag = ""
                    if sug[0] is not None:
                        if s < sug[0] and o < sug[1]: tag = "  ★"
                        elif s <= sug[0] and o <= sug[1] and (s < sug[0] or o < sug[1]): tag = "  ☆"
                    print(f"  {os.path.basename(ckpt):<35} {inst:>12}: "
                          f"steps={s:>4}  arith_ops={o:>6}{tag}", flush=True)
                else:
                    print(f"  {os.path.basename(ckpt):<35} {inst:>12}: TRUNCATED",
                          flush=True)

        # Aggregate across checkpoints per instance.
        print("\n=== held-out summary across checkpoints ===")
        print(f"{'instance':>12} {'sugar_s':>8} {'sugar_op':>9}  "
              f"{'mean_s':>7} {'std_s':>5} {'min_s':>5}  "
              f"{'mean_op':>8} {'std_op':>6} {'min_op':>6}  "
              f"{'beat_both':>10}")
        for inst in args.instances:
            vals = [results[(c, inst)] for c in args.checkpoints
                    if results[(c, inst)] is not None]
            if not vals:
                continue
            steps_arr = np.array([s for s, _, _ in vals])
            ops_arr = np.array([o for _, _, o in vals])
            sug_s, sug_op = sugar[inst]
            n_beat = int(((steps_arr < sug_s) & (ops_arr < sug_op)).sum())
            print(f"{inst:>12} {sug_s:>8} {sug_op:>9}  "
                  f"{steps_arr.mean():>7.1f} {steps_arr.std():>5.1f} {steps_arr.min():>5}  "
                  f"{ops_arr.mean():>8.0f} {ops_arr.std():>6.0f} {ops_arr.min():>6}  "
                  f"{n_beat:>4}/{len(vals)}")
    finally:
        env.close()


if __name__ == "__main__":
    main()

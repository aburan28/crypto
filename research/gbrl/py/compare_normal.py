# Compare best PPO seed vs normal (F4-style lowest-lcm-degree selection)
# across the diverse Semaev test set. Normal is the closest single-pair
# analog of F4's degree-by-degree selection — beating it is the strongest
# "vs F4 default ordering" claim our framework can make.

import argparse
import os

import torch

from baseline import run as baseline_run
from env import GbrlEnv, featurize
from model import PointerPolicy


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", default=os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "target", "release", "rl_server")))
    ap.add_argument("--checkpoints", nargs="+", required=True)
    ap.add_argument("--instances", nargs="+", required=True)
    args = ap.parse_args()

    env = GbrlEnv(args.binary)
    try:
        # Per instance: collect sugar / normal / each PPO seed
        rows = []
        for inst in args.instances:
            sug = baseline_run(args.binary, inst, "sugar")
            nrm = baseline_run(args.binary, inst, "normal")
            best_ppo = None
            for ckpt in args.checkpoints:
                policy = PointerPolicy()
                policy.load_state_dict(torch.load(ckpt, map_location="cpu"))
                policy.eval()
                r = greedy(env, policy, inst)
                if r and (best_ppo is None or r[2] < best_ppo[2] or
                          (r[2] == best_ppo[2] and r[0] < best_ppo[0])):
                    best_ppo = r
            rows.append((inst, sug, nrm, best_ppo))

        # Print
        print(f"{'instance':>14} | {'sugar':>14} | {'normal':>14} | {'PPO-best':>14}"
              f" | vs sugar | vs normal")
        for inst, sug, nrm, ppo in rows:
            sug_s, _, sug_op = sug
            nrm_s, _, nrm_op = nrm
            if ppo is None:
                print(f"{inst:>14} | PPO TRUNCATED")
                continue
            p_s, _, p_op = ppo
            vs_sugar = ""
            if p_s < sug_s and p_op < sug_op: vs_sugar = "★ both"
            elif p_s < sug_s: vs_sugar = "steps"
            elif p_op < sug_op: vs_sugar = "ops"
            vs_normal = ""
            if p_s < nrm_s and p_op < nrm_op: vs_normal = "★ both"
            elif p_s < nrm_s: vs_normal = "steps"
            elif p_op < nrm_op: vs_normal = "ops"
            print(f"{inst:>14} | {sug_s:>5}/{sug_op:>7} | {nrm_s:>5}/{nrm_op:>7} | "
                  f"{p_s:>5}/{p_op:>7} | {vs_sugar:>8} | {vs_normal:>9}")

        # Aggregate
        n = len(rows)
        ppo_beat_sugar_both = sum(1 for _, s, _, p in rows
                                  if p and p[0] < s[0] and p[2] < s[2])
        ppo_beat_normal_both = sum(1 for _, _, nrm, p in rows
                                   if p and p[0] < nrm[0] and p[2] < nrm[2])
        ppo_beat_sugar_steps = sum(1 for _, s, _, p in rows if p and p[0] < s[0])
        ppo_beat_normal_steps = sum(1 for _, _, nrm, p in rows if p and p[0] < nrm[0])
        ppo_beat_sugar_ops = sum(1 for _, s, _, p in rows if p and p[2] < s[2])
        ppo_beat_normal_ops = sum(1 for _, _, nrm, p in rows if p and p[2] < nrm[2])
        print(f"\nacross {n} instances, best-of-{len(args.checkpoints)} PPO:")
        print(f"  vs sugar  : steps {ppo_beat_sugar_steps}/{n}  "
              f"ops {ppo_beat_sugar_ops}/{n}  both {ppo_beat_sugar_both}/{n}")
        print(f"  vs normal : steps {ppo_beat_normal_steps}/{n}  "
              f"ops {ppo_beat_normal_ops}/{n}  both {ppo_beat_normal_both}/{n}")
    finally:
        env.close()


if __name__ == "__main__":
    main()

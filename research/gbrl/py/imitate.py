# Imitation-learning warm start. Run the sugar strategy across many rollouts,
# record (state, sugar's-action) at every step, train the pointer policy with
# cross-entropy. Output weights can be loaded by train.py for PPO fine-tuning.
#
# Why bother: random-init PPO on cyclic-6 spends many updates just discovering
# the sugar heuristic. With an IL warm start, the policy already matches sugar
# greedily; PPO then only has to find improvements over sugar.
#
#     python3 imitate.py --instances cyclic-5 cyclic-6 --rollouts 30 --save /tmp/cyc56_il.pt

import argparse
import os
import random
import time

import numpy as np
import torch
import torch.nn.functional as F

from env import GbrlEnv, featurize
from model import PointerPolicy


def sugar_action(obs):
    pairs = obs["pairs"]
    return min(range(len(pairs)),
               key=lambda k: (pairs[k]["sugar"], pairs[k]["lcm_degree"]))


def collect(env, instances, n_rollouts, max_steps_per_rollout=3000):
    """Run sugar on `n_rollouts` instances (sampled from `instances`),
    recording (pair_feats, global_feats, sugar_action_idx) at every state."""
    data = []
    rollout_lens = []
    for ri in range(n_rollouts):
        inst = random.choice(instances)
        obs = env.reset(inst)
        steps = 0
        while True:
            f = featurize(obs)
            if f is None:
                break
            pair_feats, global_feats = f
            a = sugar_action(obs)
            data.append((pair_feats, global_feats, a))
            r = env.step(a)
            steps += 1
            obs = r.obs
            if r.done or steps >= max_steps_per_rollout:
                break
        rollout_lens.append(steps)
        if (ri + 1) % 10 == 0:
            print(f"  rollout {ri+1}/{n_rollouts}: last={steps} avg={np.mean(rollout_lens):.0f}",
                  flush=True)
    return data, rollout_lens


def train_supervised(policy, optimizer, data, epochs, batch_size, device):
    n = len(data)
    feat_dim = data[0][0].shape[1]
    global_dim = data[0][1].shape[0]
    for ep in range(epochs):
        perm = np.random.permutation(n)
        losses = []
        accs = []
        t0 = time.time()
        for start in range(0, n, batch_size):
            idx = perm[start:start + batch_size]
            B = len(idx)
            P_max = max(data[i][0].shape[0] for i in idx)
            pair_feats_pad = np.zeros((B, P_max, feat_dim), dtype=np.float32)
            pair_mask = np.zeros((B, P_max), dtype=bool)
            global_feats_stack = np.zeros((B, global_dim), dtype=np.float32)
            actions = np.zeros(B, dtype=np.int64)
            for i, di in enumerate(idx):
                pf, gf, a = data[di]
                pair_feats_pad[i, :pf.shape[0]] = pf
                pair_mask[i, :pf.shape[0]] = True
                global_feats_stack[i] = gf
                actions[i] = a
            pf_t = torch.from_numpy(pair_feats_pad).to(device)
            gf_t = torch.from_numpy(global_feats_stack).to(device)
            mask_t = torch.from_numpy(pair_mask).to(device)
            actions_t = torch.from_numpy(actions).to(device)

            logits, _ = policy.forward_batched(pf_t, gf_t, mask_t)
            # Cross-entropy over the masked logit distribution.
            # Logits at padded positions are already -inf; CE handles correctly.
            loss = F.cross_entropy(logits, actions_t)

            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(policy.parameters(), 1.0)
            optimizer.step()

            losses.append(loss.item())
            with torch.no_grad():
                preds = logits.argmax(dim=-1)
                accs.append((preds == actions_t).float().mean().item())
        dt = time.time() - t0
        print(f"  epoch {ep+1:2d}/{epochs}: loss={np.mean(losses):.4f}  "
              f"acc={np.mean(accs):.3f}  ({dt:.1f}s)", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", default=os.path.join(
        os.path.dirname(__file__), "..", "target", "release", "rl_server"))
    ap.add_argument("--instances", nargs="+", required=True)
    ap.add_argument("--rollouts", type=int, default=30)
    ap.add_argument("--rollout-max-steps", type=int, default=3000)
    ap.add_argument("--epochs", type=int, default=10)
    ap.add_argument("--batch-size", type=int, default=64)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--save", required=True)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    device = torch.device(args.device)

    binary = os.path.abspath(args.binary)
    env = GbrlEnv(binary)
    try:
        print(f"collecting sugar rollouts: {args.rollouts} across {args.instances} "
              f"(max {args.rollout_max_steps} steps)...", flush=True)
        data, lens = collect(env, args.instances, args.rollouts, args.rollout_max_steps)
        print(f"collected {len(data)} (state, action) pairs from {len(lens)} rollouts "
              f"(mean len {np.mean(lens):.0f})", flush=True)
    finally:
        env.close()

    policy = PointerPolicy().to(device)
    optimizer = torch.optim.Adam(policy.parameters(), lr=args.lr)
    print(f"training {args.epochs} epochs, batch {args.batch_size}, lr {args.lr}...", flush=True)
    train_supervised(policy, optimizer, data, args.epochs, args.batch_size, device)

    torch.save(policy.state_dict(), args.save)
    print(f"saved to {args.save}", flush=True)


if __name__ == "__main__":
    main()

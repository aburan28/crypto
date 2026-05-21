# REINFORCE with learned baseline. Simpler than PPO: no clipping, no value
# regression replay, no advantage normalization heroics — just one update per
# trajectory with (R_t - V(s_t)) as the policy gradient signal.
#
# Why have this alongside PPO: smaller surface area, easier to diagnose, and
# the canonical algorithm Peifer et al used for Buchberger pair selection.
#
#     python3 train_reinforce.py --instances cyclic-4 --episodes 400

import argparse
import os
import random
from collections import deque

import numpy as np
import torch
import torch.nn.functional as F
from torch.distributions import Categorical

from env import GbrlEnv, featurize
from model import PointerPolicy


def rollout_with_grad(env, policy, instance, device, max_steps):
    obs = env.reset(instance)
    logps, values, rewards, entropies = [], [], [], []
    steps = 0
    while True:
        f = featurize(obs)
        if f is None:
            break
        pair_feats, global_feats = f
        pf = torch.from_numpy(pair_feats).to(device)
        gf = torch.from_numpy(global_feats).to(device)
        logits, value = policy(pf, gf)
        dist = Categorical(logits=logits)
        action = dist.sample()
        logps.append(dist.log_prob(action))
        values.append(value)
        entropies.append(dist.entropy())
        r = env.step(int(action.item()))
        rewards.append(r.reward)
        steps += 1
        obs = r.obs
        if r.done or steps >= max_steps:
            break
    return logps, values, rewards, entropies, steps


def discounted_returns(rewards, gamma: float):
    R = 0.0
    out = [0.0] * len(rewards)
    for t in reversed(range(len(rewards))):
        R = rewards[t] + gamma * R
        out[t] = R
    return out


def greedy_rollout(env, policy, instance, device, max_steps):
    obs = env.reset(instance)
    steps = 0
    while True:
        f = featurize(obs)
        if f is None:
            break
        pf = torch.from_numpy(f[0]).to(device)
        gf = torch.from_numpy(f[1]).to(device)
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
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", default=os.path.join(
        os.path.dirname(__file__), "..", "target", "release", "rl_server"))
    ap.add_argument("--instances", nargs="+", default=["cyclic-4"])
    ap.add_argument("--episodes", type=int, default=400)
    ap.add_argument("--max-steps", type=int, default=2000)
    ap.add_argument("--lr", type=float, default=3e-4)
    ap.add_argument("--gamma", type=float, default=0.99)
    ap.add_argument("--vc", type=float, default=0.5)
    ap.add_argument("--ec", type=float, default=0.0)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--eval-every", type=int, default=20)
    ap.add_argument("--save", default=None)
    args = ap.parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    device = torch.device(args.device)

    binary = os.path.abspath(args.binary)
    env = GbrlEnv(binary)
    policy = PointerPolicy().to(device)
    optimizer = torch.optim.Adam(policy.parameters(), lr=args.lr)

    history = {inst: deque(maxlen=40) for inst in args.instances}
    try:
        for ep in range(args.episodes):
            instance = random.choice(args.instances)
            logps, values, rewards, entropies, steps = rollout_with_grad(
                env, policy, instance, device, args.max_steps)
            if not logps:
                continue

            returns_raw = torch.tensor(discounted_returns(rewards, args.gamma),
                                        dtype=torch.float32, device=device)
            values_t = torch.stack(values)
            logps_t = torch.stack(logps)
            entropies_t = torch.stack(entropies)

            # Normalize the advantage signal — without this REINFORCE on negative
            # unbounded rewards (-1 per step) produces huge, unstable gradients.
            adv_raw = (returns_raw - values_t).detach()
            advantage = (adv_raw - adv_raw.mean()) / (adv_raw.std() + 1e-8)
            policy_loss = -(logps_t * advantage).mean()
            # Train value head on raw (un-normalized) returns so it stays calibrated.
            value_loss = F.mse_loss(values_t, returns_raw)
            entropy = entropies_t.mean()

            loss = policy_loss + args.vc * value_loss - args.ec * entropy
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(policy.parameters(), 0.5)
            optimizer.step()

            history[instance].append(steps)
            if (ep + 1) % 10 == 0:
                avg = {k: f"{np.mean(v):5.1f}" for k, v in history.items() if v}
                print(f"ep={ep+1:4d} inst={instance:9s} steps={steps:4d} "
                      f"pl={policy_loss.item():+.3f} vl={value_loss.item():6.2f} "
                      f"ent={entropy.item():.2f} avg={avg}", flush=True)
            if (ep + 1) % args.eval_every == 0:
                for inst in args.instances:
                    g = greedy_rollout(env, policy, inst, device, args.max_steps)
                    print(f"  [eval] greedy {inst}: steps={g}", flush=True)
    finally:
        if args.save:
            torch.save(policy.state_dict(), args.save)
            print(f"saved weights to {args.save}", flush=True)
        env.close()


if __name__ == "__main__":
    main()

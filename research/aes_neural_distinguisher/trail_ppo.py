"""PPO trail proposer with per-round reward shaping.

Improvements over trail_rl.py (REINFORCE):
  - Actor-critic with learned value function for credit assignment
  - Per-step reward (each S-box visit costs -1)
  - Round-end shaping: at the boundary, extra penalty of -active_count_in_next_round
  - GAE advantage estimation
  - Clipped surrogate objective (PPO)
"""

from __future__ import annotations

import argparse
import time
from typing import Dict, List

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from trail_rl import TrailEnv, parse_delta


class ActorCritic(nn.Module):
    def __init__(self, obs_dim: int = 49, hidden: int = 256, n_actions: int = 256):
        super().__init__()
        self.body = nn.Sequential(
            nn.Linear(obs_dim, hidden), nn.ReLU(inplace=True),
            nn.Linear(hidden, hidden),  nn.ReLU(inplace=True),
            nn.Linear(hidden, hidden),  nn.ReLU(inplace=True),
        )
        self.policy_head = nn.Linear(hidden, n_actions)
        self.value_head  = nn.Linear(hidden, 1)

    def forward(self, x):
        h = self.body(x)
        return self.policy_head(h), self.value_head(h).squeeze(-1)


def load_actor_critic_checkpoint(model: ActorCritic, path: str, device: str):
    """Load either a native ActorCritic checkpoint or an older PolicyMLP
    imitation checkpoint into the actor side of PPO.

    Older DAgger runs saved trail_rl.PolicyMLP as net.{0,2,4,6}.  The PPO model
    has the same three hidden layers plus a separate policy head, so those
    weights can warm-start the actor while leaving the value head freshly
    initialized.
    """
    ckpt = torch.load(path, map_location=device, weights_only=False)
    state = ckpt.get("model", ckpt)
    own = model.state_dict()

    if "body.0.weight" in state:
        copied = 0
        partial = 0
        for key, dst_t in own.items():
            if key not in state:
                continue
            src_t = state[key]
            if src_t.shape == dst_t.shape:
                dst_t.copy_(src_t)
                copied += 1
            elif src_t.ndim == 2 and dst_t.ndim == 2:
                rows = min(src_t.shape[0], dst_t.shape[0])
                cols = min(src_t.shape[1], dst_t.shape[1])
                dst_t.zero_()
                dst_t[:rows, :cols].copy_(src_t[:rows, :cols])
                partial += 1
        model.load_state_dict(own)
        return f"actor-critic tensors={copied} partial={partial}"

    policy_map = {
        "net.0.weight": "body.0.weight",
        "net.0.bias": "body.0.bias",
        "net.2.weight": "body.2.weight",
        "net.2.bias": "body.2.bias",
        "net.4.weight": "body.4.weight",
        "net.4.bias": "body.4.bias",
        "net.6.weight": "policy_head.weight",
        "net.6.bias": "policy_head.bias",
    }
    copied = 0
    for src, dst in policy_map.items():
        if src not in state or dst not in own:
            continue
        src_t = state[src]
        dst_t = own[dst]
        if src_t.shape == dst_t.shape:
            dst_t.copy_(src_t)
            copied += 1
        elif src.endswith("weight") and src_t.ndim == 2 and dst_t.ndim == 2:
            rows = min(src_t.shape[0], dst_t.shape[0])
            cols = min(src_t.shape[1], dst_t.shape[1])
            dst_t.zero_()
            dst_t[:rows, :cols].copy_(src_t[:rows, :cols])
            copied += 1
    model.load_state_dict(own)
    return f"policy-mlp actor warm-start tensors={copied}"


def collect_rollouts(env: TrailEnv, model: ActorCritic, delta: np.ndarray,
                     n_rollouts: int, device: str, shape_coef: float) -> Dict:
    obs_list, act_list, mask_list = [], [], []
    rew_list, done_list, val_list, lp_list = [], [], [], []
    actives, ep_lengths = [], []

    for _ in range(n_rollouts):
        obs = env.reset(delta)
        ep_len = 0
        done = False
        while not done:
            mask_np = env.valid_action_mask().copy()
            x = torch.from_numpy(obs).to(device)
            mask_t = torch.from_numpy(mask_np).to(device)
            with torch.no_grad():
                logits, value = model(x)
                logits = logits.masked_fill(~mask_t, float("-inf"))
                dist = torch.distributions.Categorical(logits=logits)
                a = dist.sample()
                log_prob = dist.log_prob(a)
            action = int(a.item())

            obs_list.append(obs.copy())
            act_list.append(action)
            mask_list.append(mask_np)
            lp_list.append(float(log_prob.item()))
            val_list.append(float(value.item()))

            prev_round = env.round
            obs, _sbox_logp, done, info = env.step(action)

            r = -1.0
            if (env.round > prev_round) and not done:
                # Round just ended and a new one begins; penalize active count of the new round.
                r -= shape_coef * len(env.active_in_round)
            rew_list.append(r)
            done_list.append(done)
            ep_len += 1

        actives.append(info["active_total"])
        ep_lengths.append(ep_len)

    return {
        "obs":       np.stack(obs_list).astype(np.float32),
        "actions":   np.array(act_list, dtype=np.int64),
        "masks":     np.stack(mask_list),
        "log_probs": np.array(lp_list, dtype=np.float32),
        "values":    np.array(val_list, dtype=np.float32),
        "rewards":   np.array(rew_list, dtype=np.float32),
        "dones":     np.array(done_list, dtype=np.bool_),
        "actives":   actives,
        "ep_lengths": ep_lengths,
    }


def compute_gae(rewards: np.ndarray, values: np.ndarray, dones: np.ndarray,
                gamma: float = 0.99, lam: float = 0.95):
    T = len(rewards)
    advantages = np.zeros(T, dtype=np.float32)
    a_next = 0.0
    v_next = 0.0
    for t in reversed(range(T)):
        if dones[t]:
            v_next = 0.0
            a_next = 0.0
        delta = rewards[t] + gamma * v_next - values[t]
        a_next = delta + gamma * lam * a_next
        advantages[t] = a_next
        v_next = values[t]
    returns = advantages + values
    return advantages, returns


def ppo_update(model: ActorCritic, opt: torch.optim.Optimizer, data: Dict,
               advantages: np.ndarray, returns: np.ndarray,
               clip: float, epochs: int, minibatch: int, ent_coef: float,
               vf_coef: float, device: str):
    obs       = torch.from_numpy(data["obs"]).to(device)
    actions   = torch.from_numpy(data["actions"]).to(device)
    masks     = torch.from_numpy(data["masks"]).to(device)
    lp_old    = torch.from_numpy(data["log_probs"]).to(device)
    adv_t     = torch.from_numpy(advantages).to(device)
    ret_t     = torch.from_numpy(returns).to(device)
    adv_t     = (adv_t - adv_t.mean()) / (adv_t.std() + 1e-8)

    N = obs.shape[0]
    for _ in range(epochs):
        perm = torch.randperm(N)
        for s in range(0, N, minibatch):
            idx = perm[s:s + minibatch]
            logits, values = model(obs[idx])
            logits = logits.masked_fill(~masks[idx], float("-inf"))
            dist = torch.distributions.Categorical(logits=logits)
            log_probs = dist.log_prob(actions[idx])
            ratio = (log_probs - lp_old[idx]).exp()
            adv = adv_t[idx]
            surr1 = ratio * adv
            surr2 = torch.clamp(ratio, 1 - clip, 1 + clip) * adv
            policy_loss = -torch.min(surr1, surr2).mean()
            value_loss  = F.mse_loss(values, ret_t[idx])
            entropy     = dist.entropy().mean()
            loss = policy_loss + vf_coef * value_loss - ent_coef * entropy
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 0.5)
            opt.step()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rounds", type=int, default=3)
    ap.add_argument("--delta",  type=str, default="00008000000000808000000000800000",
                    help="Input difference. Default is MILP-optimal 3-round AES.")
    ap.add_argument("--updates",    type=int, default=200)
    ap.add_argument("--rollouts",   type=int, default=64)
    ap.add_argument("--lr",         type=float, default=3e-4)
    ap.add_argument("--clip",       type=float, default=0.2)
    ap.add_argument("--ppo-epochs", type=int, default=4)
    ap.add_argument("--minibatch",  type=int, default=256)
    ap.add_argument("--shape-coef", type=float, default=1.0)
    ap.add_argument("--entropy",    type=float, default=0.02)
    ap.add_argument("--vf-coef",    type=float, default=0.5)
    ap.add_argument("--gamma",      type=float, default=0.99)
    ap.add_argument("--lam",        type=float, default=0.95)
    ap.add_argument("--device",     type=str, default="cpu")
    ap.add_argument("--log-every",  type=int, default=5)
    ap.add_argument("--seed",       type=int, default=0)
    ap.add_argument("--init-from",  type=str, default=None,
                    help="Path to a PPO checkpoint to warm-start from (for curriculum).")
    ap.add_argument("--save",       type=str, default=None,
                    help="Path to save the final policy.")
    args = ap.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    delta = parse_delta(args.delta)
    env = TrailEnv(args.rounds)
    model = ActorCritic().to(args.device)
    if args.init_from:
        msg = load_actor_critic_checkpoint(model, args.init_from, args.device)
        print(f"warm-start from {args.init_from}  {msg}")
    opt = torch.optim.Adam(model.parameters(), lr=args.lr, eps=1e-5)

    print(f"device={args.device}  rounds={args.rounds}  delta={args.delta}  "
          f"input_active={int((delta != 0).sum())}")
    print(f"updates={args.updates}  rollouts/update={args.rollouts}  "
          f"shape_coef={args.shape_coef}  entropy={args.entropy}")

    best_active = 10**9
    t0 = time.time()
    for upd in range(args.updates):
        data = collect_rollouts(env, model, delta, args.rollouts, args.device, args.shape_coef)
        advantages, returns = compute_gae(data["rewards"], data["values"], data["dones"],
                                          gamma=args.gamma, lam=args.lam)
        ppo_update(model, opt, data, advantages, returns,
                   clip=args.clip, epochs=args.ppo_epochs, minibatch=args.minibatch,
                   ent_coef=args.entropy, vf_coef=args.vf_coef, device=args.device)

        mean_act = float(np.mean(data["actives"]))
        cur_best = int(min(data["actives"]))
        best_active = min(best_active, cur_best)

        if (upd + 1) % args.log_every == 0:
            elapsed = time.time() - t0
            print(f"upd {upd+1:4d}/{args.updates}  mean_active={mean_act:5.2f}  "
                  f"batch_best={cur_best:3d}  alltime_best={best_active:3d}  ({elapsed:.0f}s)")

    print(f"\nFinal: alltime_best_active = {best_active}")
    if args.save:
        torch.save({"model": model.state_dict(), "args": vars(args), "best_active": best_active}, args.save)
        print(f"saved -> {args.save}")


if __name__ == "__main__":
    main()

"""DAgger-style imitation learning: train a policy to reproduce the oracle's
optimal value-level trail for 3-round AES.

The oracle (trail_oracle.find_best_3round_value_trail) returns the highest-
probability 9-active S-box trail for a given diagonal-input value 0x?? in
{1..255}. We generate a teacher dataset by varying the input value, recording
(observation, optimal_action) pairs as the env walks the trail, and train the
trail_rl.PolicyMLP via supervised cross-entropy. We then evaluate by rolling
the trained policy out *greedily* on held-out input values and counting active
S-boxes — a successful student should reproduce 9 active even on values it
didn't train on.

The trained checkpoint can be passed to trail_ppo.py via --init-from for further
RL fine-tuning.
"""

from __future__ import annotations

import argparse
import time
from typing import List, Tuple

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from trail_rl import TrailEnv
from trail_ppo import ActorCritic
from trail_oracle import find_best_3round_value_trail, assemble_input_delta


def parse_values(spec: str) -> List[int]:
    values: List[int] = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            lo_s, hi_s = part.split("-", 1)
            lo, hi = int(lo_s, 0), int(hi_s, 0)
            values.extend(range(lo, hi + 1))
        else:
            values.append(int(part, 0))
    out = sorted(set(values))
    if any(v < 1 or v > 255 for v in out):
        raise ValueError("teacher values must be in 1..255")
    return out


def oracle_trajectory_for_value(v: int):
    """Return (action_sequence, log2p, active_count) for the highest-probability
    9-active trail with diagonal input value v, optimized over cancel_row."""
    best = None
    for row in range(4):
        b = find_best_3round_value_trail(v, cancel_row=row)
        if b is not None and (best is None or b["log2p"] > best["log2p"]):
            best = b
            best["row"] = row
    if best is None:
        return None
    actions = list(best["r1_v"]) + [best["r2_dout"]] + list(best["r3_douts"])
    return actions, best["log2p"], best["active"]


def collect_teacher_pairs(values: List[int], rounds: int = 3):
    """For each diagonal-input value v in `values`, drive the env through the
    oracle's optimal action sequence and record (obs, action) pairs."""
    env = TrailEnv(rounds)
    pairs: List[Tuple[np.ndarray, int]] = []
    skipped = 0
    for v in values:
        traj = oracle_trajectory_for_value(v)
        if traj is None:
            skipped += 1
            continue
        actions, _logp, _active = traj
        delta_bytes = np.frombuffer(assemble_input_delta(v), dtype=np.uint8).copy()
        obs = env.reset(delta_bytes)
        for action in actions:
            pairs.append((obs.copy(), int(action)))
            obs, _, done, _ = env.step(action)
            if done:
                break
    return pairs, skipped


def policy_logits(model: ActorCritic, x: torch.Tensor) -> torch.Tensor:
    logits, _value = model(x)
    return logits


def evaluate_greedy(model: ActorCritic, values: List[int], rounds: int, device: str):
    env = TrailEnv(rounds)
    model.eval()
    actives, exactly_optimal = [], 0
    for v in values:
        delta_bytes = np.frombuffer(assemble_input_delta(v), dtype=np.uint8).copy()
        obs = env.reset(delta_bytes)
        done = False
        while not done:
            x = torch.from_numpy(obs).to(device)
            mask = torch.from_numpy(env.valid_action_mask()).to(device)
            with torch.no_grad():
                logits = policy_logits(model, x).masked_fill(~mask, float("-inf"))
                action = int(torch.argmax(logits).item())
            obs, _, done, info = env.step(action)
        actives.append(info["active_total"])
        if info["active_total"] == 9:
            exactly_optimal += 1
    model.train()
    return actives, exactly_optimal


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rounds", type=int, default=3)
    ap.add_argument("--train-values", type=int, default=200,
                    help="Use input values 1..N for training (N <= 255).")
    ap.add_argument("--values", type=str, default=None,
                    help="Explicit teacher values, e.g. '0x80' or '1-32,0x80'. Overrides --train-values.")
    ap.add_argument("--epochs", type=int, default=80)
    ap.add_argument("--batch", type=int, default=128)
    ap.add_argument("--lr", type=float, default=3e-4)
    ap.add_argument("--device", type=str, default="cpu")
    ap.add_argument("--save", type=str, default="ckpt_dagger.pt")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    if args.values:
        train_values = parse_values(args.values)
        train_set = set(train_values)
        test_values = [v for v in range(1, 256) if v not in train_set]
    else:
        train_values = list(range(1, args.train_values + 1))
        test_values  = list(range(args.train_values + 1, 256))

    print(f"Building teacher dataset via oracle...", flush=True)
    train_pairs, train_skipped = collect_teacher_pairs(train_values, args.rounds)
    print(f"  train: {len(train_pairs)} (obs, action) pairs from "
          f"{len(train_values) - train_skipped} values ({train_skipped} skipped)", flush=True)

    model = ActorCritic().to(args.device)
    opt = torch.optim.Adam(model.parameters(), lr=args.lr)

    obs_t   = torch.from_numpy(np.stack([p[0] for p in train_pairs])).to(args.device)
    acts_t  = torch.tensor([p[1] for p in train_pairs], dtype=torch.long, device=args.device)
    N = obs_t.shape[0]

    print(f"\nTraining DAgger student ({args.epochs} epochs, lr={args.lr})...", flush=True)
    t0 = time.time()
    for ep in range(args.epochs):
        perm = torch.randperm(N)
        running, n = 0.0, 0
        for s in range(0, N, args.batch):
            idx = perm[s : s + args.batch]
            logits = policy_logits(model, obs_t[idx])
            loss = F.cross_entropy(logits, acts_t[idx])
            opt.zero_grad()
            loss.backward()
            opt.step()
            running += loss.item() * idx.numel()
            n += idx.numel()
        if (ep + 1) % 10 == 0:
            print(f"  epoch {ep+1:3d}/{args.epochs}  loss={running/n:.4f}  ({time.time()-t0:.0f}s)", flush=True)

    print()
    print("Evaluating greedily...")
    tr_act, tr_opt = evaluate_greedy(model, train_values, args.rounds, args.device)
    te_act, te_opt = evaluate_greedy(model, test_values,  args.rounds, args.device)
    print(f"  train (n={len(train_values)}): mean_active={np.mean(tr_act):5.2f}  "
          f"min={min(tr_act)}  max={max(tr_act)}  pct@9={100*tr_opt/len(tr_act):.1f}%")
    print(f"  test  (n={len(test_values)}): mean_active={np.mean(te_act):5.2f}  "
          f"min={min(te_act)}  max={max(te_act)}  pct@9={100*te_opt/len(te_act):.1f}%")

    if args.save:
        torch.save({
            "model": model.state_dict(),
            "args": vars(args),
            "format": "actor_critic_from_oracle_imitation",
        }, args.save)
        print(f"\nsaved -> {args.save}")


if __name__ == "__main__":
    main()

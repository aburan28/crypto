"""REINFORCE trail proposer for round-reduced AES.

The agent walks through the cipher round-by-round, picking an output difference
delta_out for each active S-box. After all S-boxes in a round are decided, the
round's SR + MC are applied deterministically and we move to the next round.

State (per step):
    [ delta_in for this round (16 bytes / 255)
    , post-SB state filled so far (16 bytes / 255)
    , round_idx / total_rounds
    , position one-hot (16) ]

Action: index in {0..255} for delta_out. Invalid choices (DDT[delta_in, .] == 0)
are masked to -inf.

Reward: episode-level R = -(total active S-boxes visited) + (small log-prob bonus
to break ties between trails of equal active-S-box count). Maximizing R drives
the agent toward truncated-differential-optimal trails.
"""

from __future__ import annotations

import argparse
import math
import time
from collections import deque
from typing import List, Tuple

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from aes import shift_rows, mix_columns, bytes_to_state, state_to_bytes
from trail_search import DDT


# --- environment -------------------------------------------------------------

class TrailEnv:
    """Stateful trail-building environment.

    reset(delta_in)  -> initial observation
    step(action)     -> (obs, reward, done, info)
    """

    def __init__(self, rounds: int):
        self.rounds = rounds

    def reset(self, delta_in_bytes: np.ndarray):
        self.diff = bytes_to_state(delta_in_bytes[None])[0].astype(np.uint8)  # pre-SB for current round
        self.post_sb = np.zeros((4, 4), dtype=np.uint8)
        self.round = 0
        self.active_in_round = self._collect_active(self.diff)
        self.active_idx = 0
        self.active_total = 0
        self.delta0 = self.diff.copy()
        return self._obs()

    @staticmethod
    def _collect_active(state: np.ndarray) -> List[Tuple[int, int]]:
        return [(i, j) for i in range(4) for j in range(4) if state[i, j] != 0]

    def current_position(self) -> Tuple[int, int]:
        return self.active_in_round[self.active_idx]

    def current_delta_in(self) -> int:
        i, j = self.current_position()
        return int(self.diff[i, j])

    def valid_action_mask(self) -> np.ndarray:
        d_in = self.current_delta_in()
        return DDT[d_in] > 0                          # (256,)

    def _obs(self):
        i, j = self.current_position()
        pos = np.zeros(16, dtype=np.float32)
        pos[4 * i + j] = 1.0
        return np.concatenate([
            (self.diff.flatten().astype(np.float32) / 255.0),
            (self.post_sb.flatten().astype(np.float32) / 255.0),
            np.array([self.round / max(1, self.rounds)], dtype=np.float32),
            pos,
        ])

    def step(self, action: int):
        i, j = self.current_position()
        d_in = self.current_delta_in()
        # log_prob_increment is computed by the caller from DDT; here we just record
        self.post_sb[i, j] = action
        self.active_total += 1
        self.active_idx += 1
        done = False
        sbox_logp = math.log2(DDT[d_in, action] / 256.0)   # negative

        if self.active_idx >= len(self.active_in_round):
            # End of round: apply SR (always) and MC (if not last round)
            sr = shift_rows(self.post_sb[None])[0]
            if self.round < self.rounds - 1:
                self.diff = mix_columns(sr[None])[0]
            else:
                self.diff = sr
            self.round += 1
            self.post_sb = np.zeros((4, 4), dtype=np.uint8)
            if self.round >= self.rounds:
                done = True
            else:
                self.active_in_round = self._collect_active(self.diff)
                self.active_idx = 0
                if not self.active_in_round:
                    # Whole state collapsed to zero — trail terminates early (probability stays at current value)
                    self.round = self.rounds
                    done = True

        return (None if done else self._obs()), sbox_logp, done, {
            "active_total": self.active_total,
        }


# --- policy -----------------------------------------------------------------

class PolicyMLP(nn.Module):
    def __init__(self, obs_dim: int = 49, hidden: int = 256, n_actions: int = 256):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(obs_dim, hidden), nn.ReLU(inplace=True),
            nn.Linear(hidden, hidden),  nn.ReLU(inplace=True),
            nn.Linear(hidden, hidden),  nn.ReLU(inplace=True),
            nn.Linear(hidden, n_actions),
        )

    def forward(self, x):
        return self.net(x)


# --- rollouts ---------------------------------------------------------------

def rollout(policy: PolicyMLP, env: TrailEnv, delta_in: np.ndarray, device: str,
            greedy: bool = False):
    obs = env.reset(delta_in)
    log_probs: List[torch.Tensor] = []
    entropies: List[torch.Tensor] = []
    obs_hist: List[np.ndarray] = []
    action_hist: List[int] = []
    mask_hist: List[np.ndarray] = []
    total_logp_bits = 0.0
    done = False
    while not done:
        mask_np = env.valid_action_mask().copy()
        x = torch.from_numpy(obs).to(device)
        logits = policy(x)
        mask = torch.from_numpy(mask_np).to(device)
        logits = logits.masked_fill(~mask, float("-inf"))
        dist = torch.distributions.Categorical(logits=logits)
        if greedy:
            a = torch.argmax(logits)
        else:
            a = dist.sample()
        action = int(a.item())
        log_probs.append(dist.log_prob(a))
        entropies.append(dist.entropy())
        obs_hist.append(obs.copy())
        action_hist.append(action)
        mask_hist.append(mask_np)
        obs, sbox_logp_bits, done, info = env.step(action)
        total_logp_bits += sbox_logp_bits
    return {
        "log_probs": log_probs,
        "entropies": entropies,
        "obs": obs_hist,
        "actions": action_hist,
        "masks": mask_hist,
        "logp_bits": total_logp_bits,
        "active": info["active_total"],
    }


def replay_log_prob(policy: PolicyMLP, obs_list, action_list, mask_list, device: str) -> torch.Tensor:
    """Recompute the sum of log-probabilities of a stored trajectory under the current policy."""
    obs_t = torch.from_numpy(np.stack(obs_list)).to(device)
    masks_t = torch.from_numpy(np.stack(mask_list)).to(device)
    actions_t = torch.tensor(action_list, dtype=torch.long, device=device)
    logits = policy(obs_t)
    logits = logits.masked_fill(~masks_t, float("-inf"))
    log_probs = F.log_softmax(logits, dim=-1)
    return log_probs.gather(1, actions_t.unsqueeze(1)).squeeze(1).sum()


# --- training loop ----------------------------------------------------------

def parse_delta(spec: str) -> np.ndarray:
    spec = spec.replace("_", "").replace(" ", "")
    if len(spec) != 32:
        raise ValueError(f"delta must be 32 hex chars, got {len(spec)}")
    return np.frombuffer(bytes.fromhex(spec), dtype=np.uint8)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--delta", type=str, default="80000000008000000000800000000080",
                    help="Input difference; default is the AES diagonal pattern.")
    ap.add_argument("--rounds", type=int, default=3)
    ap.add_argument("--episodes", type=int, default=20000)
    ap.add_argument("--batch", type=int, default=64,
                    help="Rollouts per gradient update.")
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--entropy", type=float, default=0.01,
                    help="Entropy bonus coefficient.")
    ap.add_argument("--alpha", type=float, default=1.0,
                    help="Weight on active-S-box term in reward (vs trail-logprob tie-break).")
    ap.add_argument("--device", type=str, default="cpu")
    ap.add_argument("--log-every", type=int, default=20)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    device = args.device

    delta = parse_delta(args.delta)
    env = TrailEnv(args.rounds)
    policy = PolicyMLP().to(device)
    opt = torch.optim.Adam(policy.parameters(), lr=args.lr)
    baseline = 0.0
    baseline_count = 0
    history = deque(maxlen=200)

    best_active = 10**9
    best_logp = -10**9
    best_choices = None

    print(f"device={device}  rounds={args.rounds}  delta={args.delta}")
    print(f"input active bytes={int((delta != 0).sum())}  episodes={args.episodes}  batch={args.batch}")

    # Replay buffer: stores stateless data (obs, actions, masks, reward) so we can
    # recompute current-policy log_probs on each replay and get a valid gradient.
    replay = []  # list of dicts
    replay_cap = 16

    t0 = time.time()
    n_updates = args.episodes // args.batch
    for upd in range(n_updates):
        traj: List[dict] = []
        for _ in range(args.batch):
            tr = rollout(policy, env, delta, device, greedy=False)
            tr["reward"] = -args.alpha * tr["active"] + 0.01 * tr["logp_bits"]
            traj.append(tr)
            history.append(tr["reward"])
            if tr["active"] < best_active or (tr["active"] == best_active and tr["logp_bits"] > best_logp):
                best_active = tr["active"]
                best_logp = tr["logp_bits"]

            # Replay: stateless snapshot.
            snapshot = {
                "obs": tr["obs"], "actions": tr["actions"], "masks": tr["masks"],
                "reward": tr["reward"],
            }
            replay.append(snapshot)
            replay.sort(key=lambda d: -d["reward"])
            replay = replay[:replay_cap]

        # REINFORCE with running-mean baseline + entropy bonus + replay.
        baseline = sum(history) / len(history)
        loss = torch.tensor(0.0, device=device)
        ent_term = torch.tensor(0.0, device=device)
        for tr in traj:
            advantage = tr["reward"] - baseline
            ep_lp = torch.stack(tr["log_probs"]).sum()
            loss = loss - advantage * ep_lp
            ent_term = ent_term + torch.stack(tr["entropies"]).mean()

        # Replay-based gradient: recompute log_probs under current policy.
        for snap in replay:
            advantage = snap["reward"] - baseline
            if advantage > 0:
                ep_lp = replay_log_prob(policy, snap["obs"], snap["actions"], snap["masks"], device)
                loss = loss - 0.25 * advantage * ep_lp

        loss = loss / args.batch - args.entropy * ent_term / args.batch
        opt.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(policy.parameters(), 5.0)
        opt.step()

        if (upd + 1) % args.log_every == 0:
            mean_act = sum(t["active"] for t in traj) / len(traj)
            mean_lp = sum(t["logp_bits"] for t in traj) / len(traj)
            elapsed = time.time() - t0
            print(f"upd {upd+1:5d}/{n_updates}  ep={((upd+1)*args.batch):6d}  "
                  f"mean_active={mean_act:5.2f}  mean_log2p={mean_lp:7.2f}  "
                  f"best_active={best_active}  best_log2p={best_logp:.2f}  ({elapsed:.0f}s)")

    # Greedy eval
    policy.eval()
    with torch.no_grad():
        tr = rollout(policy, env, delta, device, greedy=True)
    print(f"\nGreedy eval: active={tr['active']}  log2(prob)={tr['logp_bits']:.2f}")
    print(f"Best found  : active={best_active}  log2(prob)={best_logp:.2f}")


if __name__ == "__main__":
    main()

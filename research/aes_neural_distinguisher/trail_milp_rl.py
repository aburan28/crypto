"""MILP-RL hybrid trail search.

Stage 1: MILP solves the truncated-differential structure (which bytes are
active per round).
Stage 2: a PPO agent picks specific byte values at each active position, with
the environment enforcing the truncated pattern as a hard constraint: any
round whose post-MC active pattern doesn't match the MILP prescription is
terminated with a large penalty.

Effect: the policy's exploration is collapsed from the full ~127^k value
combinations down to just the ones that satisfy the structural cancellation,
which is exactly where pure RL choked on 3-round AES.
"""

from __future__ import annotations

import argparse
import time
from typing import Dict, List

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from aes import shift_rows, mix_columns, bytes_to_state
from trail_search import DDT
from trail_rl import TrailEnv, parse_delta
from trail_milp import solve_truncated_trail, pattern_to_delta_bytes
from trail_ppo import ActorCritic, compute_gae, ppo_update, load_actor_critic_checkpoint


def sr_pattern(pattern: np.ndarray) -> np.ndarray:
    out = np.zeros_like(pattern)
    for i in range(4):
        out[i] = np.roll(pattern[i], -i)
    return out


def oracle_milp_patterns(input_pattern: np.ndarray, cancel_row: int) -> List[np.ndarray]:
    """Return a concrete 3-round 4->1->4 MILP-optimal activity pattern.

    The unconstrained MILP may choose any symmetry-equivalent cancellation row.
    For DAgger/PPO warm starts we sometimes need the exact symmetry used by the
    oracle teacher, so this constructs the same wide-trail structure directly.
    """
    sr0 = sr_pattern(input_pattern)
    active_cols = np.flatnonzero(sr0.any(axis=0))
    if len(active_cols) != 1:
        raise ValueError("--oracle-pattern-row requires an input pattern that ShiftRows maps into one column")
    c0 = int(active_cols[0])

    p1 = np.zeros((4, 4), dtype=np.uint8)
    p1[cancel_row, c0] = 1

    sr1 = sr_pattern(p1)
    active_cols = np.flatnonzero(sr1.any(axis=0))
    c1 = int(active_cols[0])
    p2 = np.zeros((4, 4), dtype=np.uint8)
    p2[:, c1] = 1

    p3 = sr_pattern(p2)
    return [input_pattern.astype(np.uint8), p1, p2, p3]


class GuidedTrailEnv(TrailEnv):
    """TrailEnv that terminates with a penalty whenever a round's post-MC
    activity pattern diverges from the MILP-prescribed truncated trail."""

    def __init__(self, rounds: int, target_patterns: List[np.ndarray]):
        super().__init__(rounds)
        self.target_patterns = target_patterns  # length rounds+1

    def _obs(self):
        base = super()._obs()
        next_round = min(self.round + 1, len(self.target_patterns) - 1)
        target = self.target_patterns[next_round].flatten().astype(np.float32)
        return np.concatenate([base, target])

    def step(self, action: int):
        prev_round = self.round
        result = super().step(action)
        obs, sbox_lp, done, info = result
        info = dict(info)
        info["pattern_violation"] = False
        if self.round > prev_round:
            # A round just transitioned. Compare env activity to MILP target.
            actual = (self.diff != 0).astype(np.uint8)
            target = self.target_patterns[min(self.round, len(self.target_patterns) - 1)]
            if not np.array_equal(actual, target):
                info["pattern_violation"] = True
                done = True
            else:
                info["pattern_match"] = True
        return obs, sbox_lp, done, info


def collect_rollouts_guided(env: GuidedTrailEnv, model: ActorCritic, delta: np.ndarray,
                            n_rollouts: int, device: str, violation_penalty: float = 50.0,
                            match_bonus: float = 3.0):
    obs_list, act_list, mask_list = [], [], []
    rew_list, done_list, val_list, lp_list = [], [], [], []
    actives, scores, ep_lengths, violations, valid = [], [], [], 0, 0

    for _ in range(n_rollouts):
        obs = env.reset(delta)
        ep_len = 0
        done = False
        while not done:
            mask_np = env.valid_action_mask().copy()
            x = torch.from_numpy(obs).to(device)
            mt = torch.from_numpy(mask_np).to(device)
            with torch.no_grad():
                logits, value = model(x)
                logits = logits.masked_fill(~mt, float("-inf"))
                dist = torch.distributions.Categorical(logits=logits)
                a = dist.sample()
                lp = dist.log_prob(a)
            action = int(a.item())

            obs_list.append(obs.copy())
            act_list.append(action)
            mask_list.append(mask_np)
            lp_list.append(float(lp.item()))
            val_list.append(float(value.item()))

            prev_round = env.round
            obs, _sbox_lp, done, info = env.step(action)

            r = -1.0
            if info.get("pattern_violation", False):
                r -= violation_penalty
                violations += 1
            elif info.get("pattern_match", False):
                r += match_bonus
            rew_list.append(r)
            done_list.append(done)
            ep_len += 1

        actives.append(info["active_total"])
        is_valid = not info.get("pattern_violation", False)
        valid += int(is_valid)
        scores.append(info["active_total"] if is_valid else info["active_total"] + violation_penalty)
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
        "scores":    scores,
        "ep_lengths": ep_lengths,
        "violations": violations,
        "valid":     valid,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rounds",     type=int, default=3)
    ap.add_argument("--updates",    type=int, default=200)
    ap.add_argument("--rollouts",   type=int, default=64)
    ap.add_argument("--lr",         type=float, default=3e-4)
    ap.add_argument("--clip",       type=float, default=0.2)
    ap.add_argument("--ppo-epochs", type=int, default=4)
    ap.add_argument("--minibatch",  type=int, default=256)
    ap.add_argument("--entropy",    type=float, default=0.05)
    ap.add_argument("--vf-coef",    type=float, default=0.5)
    ap.add_argument("--gamma",      type=float, default=0.99)
    ap.add_argument("--lam",        type=float, default=0.95)
    ap.add_argument("--violation",  type=float, default=50.0,
                    help="Reward penalty for a pattern-violating action.")
    ap.add_argument("--match-bonus", type=float, default=3.0,
                    help="Reward bonus whenever a completed round matches the target pattern.")
    ap.add_argument("--device",     type=str, default="cpu")
    ap.add_argument("--log-every",  type=int, default=10)
    ap.add_argument("--delta",      type=str, default=None,
                    help="Optional fixed input difference. When set, MILP is solved with this input activity pattern.")
    ap.add_argument("--oracle-pattern-row", type=int, default=None,
                    help="For 3-round SR-aligned inputs, force the symmetry-equivalent 4->1->4 pattern to this cancellation row.")
    ap.add_argument("--init-from",  type=str, default=None,
                    help="Warm-start actor from an ActorCritic/PPO checkpoint or an older PolicyMLP DAgger checkpoint.")
    ap.add_argument("--save",       type=str, default=None,
                    help="Path to save the final hybrid policy.")
    ap.add_argument("--seed",       type=int, default=0)
    args = ap.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    print(f"Stage 1: MILP truncated trail for {args.rounds} rounds...")
    fixed_delta = None
    input_pattern = None
    if args.delta:
        fixed_delta = parse_delta(args.delta)
        input_pattern = (bytes_to_state(fixed_delta[None])[0] != 0).astype(np.uint8)
    sol = solve_truncated_trail(args.rounds, input_pattern=input_pattern)
    if sol is None:
        print("MILP failed."); return
    if args.oracle_pattern_row is not None:
        if args.rounds != 3 or input_pattern is None:
            raise ValueError("--oracle-pattern-row currently requires --rounds 3 and --delta")
        if not 0 <= args.oracle_pattern_row <= 3:
            raise ValueError("--oracle-pattern-row must be in 0..3")
        sol["patterns"] = oracle_milp_patterns(input_pattern, args.oracle_pattern_row)
        sol["objective"] = sum(int(p.sum()) for p in sol["patterns"][:-1])
        sol["status"] = "Optimal (fixed oracle symmetry)"
    print(f"  truncated optimum: {sol['objective']} active S-boxes")
    target_patterns = sol["patterns"]
    if fixed_delta is None:
        delta = np.frombuffer(pattern_to_delta_bytes(target_patterns[0]), dtype=np.uint8).copy()
    else:
        delta = fixed_delta
    print(f"  input delta (hex): {delta.tobytes().hex()}")
    print(f"  input active     : {int((delta != 0).sum())}")
    print()
    print("Stage 2: PPO with truncated-pattern enforcement...")

    env = GuidedTrailEnv(args.rounds, target_patterns)
    model = ActorCritic(obs_dim=65).to(args.device)
    if args.init_from:
        msg = load_actor_critic_checkpoint(model, args.init_from, args.device)
        print(f"warm-start from {args.init_from}  {msg}")
    opt = torch.optim.Adam(model.parameters(), lr=args.lr, eps=1e-5)

    best_score = 10**9
    successful_episodes = 0
    t0 = time.time()
    for upd in range(args.updates):
        data = collect_rollouts_guided(env, model, delta, args.rollouts, args.device,
                                        violation_penalty=args.violation,
                                        match_bonus=args.match_bonus)
        advantages, returns = compute_gae(data["rewards"], data["values"], data["dones"],
                                          gamma=args.gamma, lam=args.lam)
        ppo_update(model, opt, data, advantages, returns,
                   clip=args.clip, epochs=args.ppo_epochs, minibatch=args.minibatch,
                   ent_coef=args.entropy, vf_coef=args.vf_coef, device=args.device)

        mean_score = float(np.mean(data["scores"]))
        cur_best = int(min(data["scores"]))
        best_score = min(best_score, cur_best)
        # An episode is "successful" if it reached the full trail (no truncated violation)
        # AND the active count matches the MILP optimum (e.g., 9 for 3-round).
        n_optimal = sum(1 for a in data["actives"] if a == sol["objective"])
        successful_episodes += n_optimal

        if (upd + 1) % args.log_every == 0:
            elapsed = time.time() - t0
            print(f"upd {upd+1:4d}/{args.updates}  mean_score={mean_score:5.2f}  "
                  f"batch_best={cur_best:3d}  alltime_best={best_score:3d}  "
                  f"violations/rollout={data['violations']/args.rollouts:.2f}  "
                  f"valid={data['valid']}/{args.rollouts}  "
                  f"optimal_in_batch={n_optimal}/{args.rollouts}  ({elapsed:.0f}s)")

    print(f"\nFinal: alltime_best_score = {best_score}  "
          f"target = {sol['objective']}  "
          f"total_optimal_episodes = {successful_episodes}")
    if args.save:
        torch.save({"model": model.state_dict(), "args": vars(args), "best_score": best_score}, args.save)
        print(f"saved -> {args.save}")


if __name__ == "__main__":
    main()

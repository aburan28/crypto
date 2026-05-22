"""Differential trail search for round-reduced AES (classical baseline).

A trail specifies, at every active S-box of the cipher, an input/output
difference pair. All non-SubBytes operations are linear or permutational on
differences, so the only probabilistic step is SubBytes. For a single S-box,
P(delta_in -> delta_out) = DDT[delta_in][delta_out] / 256.

This file computes:
  - the AES S-box DDT
  - forward trail propagation (deterministic given per-S-box output choices)
  - a beam search over per-round S-box choices that maximizes trail probability

The output is a high-probability trail plus the *input difference* implied by
the search start, which can be fed back into the distinguisher as a better
choice than an arbitrary single-byte delta.
"""

from __future__ import annotations

import argparse
import heapq
from dataclasses import dataclass, field
from typing import List, Tuple

import numpy as np

from aes import SBOX, bytes_to_state, state_to_bytes, shift_rows, mix_columns


# --- DDT ---------------------------------------------------------------------

def build_ddt() -> np.ndarray:
    """DDT[a, b] = #{x : S(x) ^ S(x ^ a) = b}. Shape (256, 256), values in [0, 256]."""
    ddt = np.zeros((256, 256), dtype=np.int32)
    xs = np.arange(256, dtype=np.uint8)
    for a in range(256):
        a_u8 = np.uint8(a)
        out = SBOX[xs] ^ SBOX[xs ^ a_u8]
        # Count occurrences of each output difference
        counts = np.bincount(out, minlength=256)
        ddt[a] = counts
    return ddt


DDT = build_ddt()
# Cache, per input difference, the list of (out_diff, prob) with non-zero prob,
# sorted by descending probability.
DDT_ROWS: List[List[Tuple[int, float]]] = []
for a in range(256):
    row = [(b, DDT[a, b] / 256.0) for b in range(256) if DDT[a, b] > 0]
    row.sort(key=lambda t: -t[1])
    DDT_ROWS.append(row)


# --- difference propagation (linear ops are deterministic) ------------------

def propagate_linear(diff_state: np.ndarray, apply_mc: bool) -> np.ndarray:
    """Apply ShiftRows (always) and MixColumns (if not last round) to a difference state.

    Input/output shape: (4, 4) uint8.
    """
    s = diff_state[None]                   # (1, 4, 4)
    s = shift_rows(s)
    if apply_mc:
        s = mix_columns(s)
    return s[0]


# --- beam search over trails -------------------------------------------------

@dataclass(order=True)
class TrailNode:
    neg_logprob: float                      # priority key: smaller = higher probability
    diff: np.ndarray = field(compare=False) # current (post-MC, pre-SubBytes-of-next-round) difference, shape (4, 4)
    choices: List[Tuple[int, int, int]] = field(default_factory=list, compare=False)
    # `choices` records, per active S-box visited: (round, byte_index_in_state, delta_out)


def _expand_round(diff_state: np.ndarray, max_branch_per_sbox: int):
    """Yield (new_state_after_SB, logprob_inc) candidates for one round's SubBytes step.

    We iterate active S-boxes one at a time, beam-pruning inside the round so the
    expansion stays bounded. Returns a list of (post_SB_state, logprob_inc, choices_in_round).
    """
    state_flat = diff_state.flatten()                          # 16 active S-boxes max
    active = [i for i, d in enumerate(state_flat) if d != 0]
    if not active:
        return [(diff_state.copy(), 0.0, [])]                  # all-zero difference; trail trivially has prob 1
    partial = [(np.zeros(16, dtype=np.uint8), 0.0, [])]        # (post_SB_so_far, lp, choices_in_round)
    for byte_idx in active:
        delta_in = int(state_flat[byte_idx])
        row = DDT_ROWS[delta_in][:max_branch_per_sbox]
        new_partial = []
        for post, lp, ch in partial:
            for delta_out, prob in row:
                np_post = post.copy()
                np_post[byte_idx] = delta_out
                new_partial.append((np_post, lp + float(np.log2(prob)), ch + [(byte_idx, delta_out)]))
        # Prune within the round to stay bounded.
        new_partial.sort(key=lambda t: -t[1])
        partial = new_partial[: max_branch_per_sbox * max_branch_per_sbox]  # rough cap
    # Reshape post_SB back to (4, 4)
    return [(p[0].reshape(4, 4), p[1], p[2]) for p in partial]


def beam_search(delta_in_bytes: np.ndarray, rounds: int, beam: int = 32, branch_per_sbox: int = 8):
    """Find the top-`beam` highest-probability trails over `rounds` AES rounds starting from `delta_in_bytes`.

    Returns a list of (log2_prob, trail_choices, terminal_diff) sorted by descending probability.
    """
    init_state = bytes_to_state(delta_in_bytes[None])[0]      # (4, 4)
    frontier: List[TrailNode] = [TrailNode(neg_logprob=0.0, diff=init_state, choices=[])]

    for r in range(rounds):
        is_last = (r == rounds - 1)
        new_frontier: List[TrailNode] = []
        for node in frontier:
            for post_sb, lp_inc, ch in _expand_round(node.diff, branch_per_sbox):
                post_lin = propagate_linear(post_sb, apply_mc=not is_last)
                new_choices = node.choices + [(r, idx, delta_out) for idx, delta_out in ch]
                new_frontier.append(TrailNode(
                    neg_logprob=node.neg_logprob - lp_inc,
                    diff=post_lin,
                    choices=new_choices,
                ))
        new_frontier.sort()
        frontier = new_frontier[:beam]

    return [(-n.neg_logprob, n.choices, n.diff) for n in frontier]


# --- helpers ----------------------------------------------------------------

def parse_delta(spec: str) -> np.ndarray:
    spec = spec.replace("_", "").replace(" ", "")
    if len(spec) != 32:
        raise ValueError(f"delta must be 32 hex chars, got {len(spec)}")
    return np.frombuffer(bytes.fromhex(spec), dtype=np.uint8)


def hex_of_state(s: np.ndarray) -> str:
    return state_to_bytes(s[None])[0].tobytes().hex()


def active_byte_count(s: np.ndarray) -> int:
    return int((s != 0).sum())


# --- entry point -------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--delta", type=str, default="80000000000000000000000000000000",
                    help="32-hex-char input difference.")
    ap.add_argument("--rounds", type=int, default=3)
    ap.add_argument("--beam", type=int, default=64,
                    help="Beam width across rounds.")
    ap.add_argument("--branch", type=int, default=8,
                    help="Top-K differences kept at each S-box within a round.")
    ap.add_argument("--show", type=int, default=5, help="Number of trails to print.")
    args = ap.parse_args()

    delta = parse_delta(args.delta)
    print(f"input delta : {args.delta}")
    print(f"rounds      : {args.rounds}")
    print(f"DDT max     : {DDT[1:, :].max()} / 256 = 2^{np.log2(DDT[1:, :].max()/256):.2f}")
    print(f"active S-boxes at round 1 input: {int((delta != 0).sum())}")
    print()

    trails = beam_search(delta, args.rounds, beam=args.beam, branch_per_sbox=args.branch)
    print(f"found {len(trails)} trails; top {min(args.show, len(trails))}:")
    for i, (lp, choices, term) in enumerate(trails[: args.show]):
        print(f"  #{i:2d}  log2(prob) = {lp:.2f}   active S-boxes used = {len(choices)}   "
              f"terminal diff hex = {hex_of_state(term)}   terminal active bytes = {active_byte_count(term)}")


if __name__ == "__main__":
    main()

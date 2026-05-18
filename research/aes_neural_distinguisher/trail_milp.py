"""Truncated-differential MILP for round-reduced AES (classical SOTA baseline).

Models only the *activity pattern* (which bytes are non-zero) of each state, not
the byte values. Minimizes total active S-box count under the wide-trail
constraint (MC branch number >= 5 per column).

This is the standard formulation used in MILP-based cryptanalysis (Mouha et al.,
Sun et al.). It gives the ground-truth minimum that ML approaches should match.

State variables:
    a[r, i, j] in {0,1}  — active flag for state byte (i, j) at "round boundary" r.
                            r = 0 is the input difference;
                            r = R is the final state (after the last round).
                            For r in [0, R-1], a[r] is the input to SubBytes of
                            round r+1, so |a[r]| = number of active S-boxes in
                            round r+1.

Round transitions:
    For rounds with MixColumns (r = 0..R-2):
        post_SR[i, j] = a[r, i, (j + i) mod 4]   # ShiftRows on rows
        For each column j:
            d_col[r, j] in {0,1}                  # is column j active in this transition?
            sum_i (post_SR[i, j] + a[r+1, i, j]) >= 5 * d_col[r, j]
            d_col[r, j] >= each byte in post_SR[:, j] and a[r+1, :, j]
    For the final round (no MC), the SR permutation is applied directly:
        a[R, i, j] = a[R-1, i, (j + i) mod 4]

Objective:
    minimize sum over r in [0..R-1] of |a[r]|
"""

from __future__ import annotations

import argparse

import numpy as np
import pulp


def solve_truncated_trail(rounds: int, input_pattern=None, min_input_active: int = 1,
                          time_limit: int = 60, verbose: bool = False):
    """Find the minimum active-S-box truncated trail for round-reduced AES.

    rounds         : number of AES rounds (last round has no MC).
    input_pattern  : optional (4, 4) ndarray of 0/1 fixing the input activity.
    min_input_active: minimum number of active input bytes (only used if input_pattern is None).
    """
    prob = pulp.LpProblem("aes_trunc_trail", pulp.LpMinimize)

    # Active byte flags.
    a = {(r, i, j): pulp.LpVariable(f"a_{r}_{i}_{j}", cat="Binary")
         for r in range(rounds + 1) for i in range(4) for j in range(4)}

    # Column-activity flags for MC rounds (r = 0..rounds-2 transitions).
    d = {(r, k): pulp.LpVariable(f"d_{r}_{k}", cat="Binary")
         for r in range(rounds - 1) for k in range(4)}

    # Input constraints.
    if input_pattern is not None:
        for i in range(4):
            for j in range(4):
                prob += a[0, i, j] == int(input_pattern[i, j])
    else:
        prob += pulp.lpSum(a[0, i, j] for i in range(4) for j in range(4)) >= min_input_active

    # MC transitions for r = 0..rounds-2.
    for r in range(rounds - 1):
        for k in range(4):
            # post-SR of round (r+1) feeds MC, but SR happens AFTER SB in round r+1.
            # By our convention, a[r] is the state right BEFORE SB of round r+1.
            # SubBytes preserves activity, so post-SB equals a[r]. Then SR maps
            # (i, j) -> (i, (j - i) mod 4). The byte landing at output col k row i
            # came from input col (k + i) mod 4.
            pre_MC = [a[r, i, (k + i) % 4] for i in range(4)]
            post_MC = [a[r + 1, i, k] for i in range(4)]
            col_sum = pulp.lpSum(pre_MC + post_MC)

            for x in pre_MC + post_MC:
                prob += d[r, k] >= x
            prob += d[r, k] <= col_sum
            prob += col_sum >= 5 * d[r, k]

    # Last round: SR only, no MC.
    for i in range(4):
        for j in range(4):
            prob += a[rounds, i, j] == a[rounds - 1, i, (j + i) % 4]

    # Objective: total active S-boxes (one S-box per active byte at a[0..R-1]).
    prob += pulp.lpSum(a[r, i, j] for r in range(rounds) for i in range(4) for j in range(4))

    solver = pulp.PULP_CBC_CMD(msg=verbose, timeLimit=time_limit)
    status = prob.solve(solver)

    if pulp.LpStatus[status] not in ("Optimal", "Not Solved"):
        return None

    patterns = []
    for r in range(rounds + 1):
        p = np.zeros((4, 4), dtype=np.uint8)
        for i in range(4):
            for j in range(4):
                p[i, j] = int(round(pulp.value(a[r, i, j])))
        patterns.append(p)

    return {
        "status": pulp.LpStatus[status],
        "objective": int(round(pulp.value(prob.objective))),
        "patterns": patterns,
    }


def pattern_to_delta_bytes(pattern: np.ndarray, active_byte_value: int = 0x80) -> bytes:
    """Convert a (4, 4) activity pattern back into a 16-byte difference (column-major)."""
    out = bytearray(16)
    for i in range(4):
        for j in range(4):
            if pattern[i, j]:
                out[4 * j + i] = active_byte_value
    return bytes(out)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rounds", type=int, default=3)
    ap.add_argument("--time-limit", type=int, default=60)
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--show-patterns", action="store_true")
    args = ap.parse_args()

    res = solve_truncated_trail(args.rounds, time_limit=args.time_limit, verbose=args.verbose)
    if res is None:
        print("MILP failed")
        return

    print(f"status     : {res['status']}")
    print(f"rounds     : {args.rounds}")
    print(f"min active S-boxes (truncated optimum): {res['objective']}")
    print()
    if args.show_patterns:
        for r, p in enumerate(res["patterns"]):
            tag = "input" if r == 0 else (f"after round {r}" if r < args.rounds else "final state")
            n = int(p.sum())
            print(f"  state {r} ({tag}, active={n}):")
            for row in p:
                print("    " + " ".join(str(int(b)) for b in row))
        print()

    delta = pattern_to_delta_bytes(res["patterns"][0])
    print(f"input active pattern (hex, value 0x80 per active byte): {delta.hex()}")
    print()
    print("Active S-box breakdown by round:")
    for r in range(args.rounds):
        n = int(res["patterns"][r].sum())
        print(f"  round {r+1}: {n} active S-boxes")


if __name__ == "__main__":
    main()

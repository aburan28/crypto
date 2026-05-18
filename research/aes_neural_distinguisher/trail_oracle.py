"""Value-level oracle for AES truncated trails.

Given a truncated trail from the MILP (which bytes are active per round), find
specific byte-value differences that realize the trail, maximizing the trail
probability over the choice of values within each active byte position.

For 3-round AES with the wide-trail input (4 active bytes that SR-align to a
single column), the key constraint is:
    MC * (delta_R1_out_col) = (one non-zero target, 0, 0, 0)        (cancellation)
Since MC is MDS (invertible), this gives a 1-parameter family
    delta_R1_out_col = MC_inverse[:, 0] * target
for any non-zero target byte. We then check, per target, that each entry of
delta_R1_out_col lies in the DDT row of its corresponding delta_in (i.e. the
SubBytes transition exists), and compute the resulting trail probability.

The output is the highest-probability value-level realization plus the implied
sequence of (delta_in, delta_out) choices at each active S-box -- exactly the
'teacher trajectory' a DAgger student would imitate.
"""

from __future__ import annotations

import argparse
from typing import List, Tuple

import numpy as np

from aes import SBOX, bytes_to_state, state_to_bytes, shift_rows, mix_columns, _xtime
from trail_search import DDT, build_ddt


# GF(2^8) primitives reused from aes.py via _xtime.
def gf_mul(a: int, b: int) -> int:
    """GF(2^8) multiplication, modulo AES polynomial x^8 + x^4 + x^3 + x + 1."""
    r = 0
    a &= 0xff
    b &= 0xff
    for _ in range(8):
        if b & 1:
            r ^= a
        hi = a & 0x80
        a = (a << 1) & 0xff
        if hi:
            a ^= 0x1b
        b >>= 1
    return r


# AES MixColumns matrix and its inverse over GF(2^8).
MC_M = np.array([
    [0x02, 0x03, 0x01, 0x01],
    [0x01, 0x02, 0x03, 0x01],
    [0x01, 0x01, 0x02, 0x03],
    [0x03, 0x01, 0x01, 0x02],
], dtype=np.uint8)

MC_INV = np.array([
    [0x0e, 0x0b, 0x0d, 0x09],
    [0x09, 0x0e, 0x0b, 0x0d],
    [0x0d, 0x09, 0x0e, 0x0b],
    [0x0b, 0x0d, 0x09, 0x0e],
], dtype=np.uint8)


def mc_column(col: np.ndarray) -> np.ndarray:
    """Apply MC to a single column (4-byte vector) over GF(2^8)."""
    out = np.zeros(4, dtype=np.uint8)
    for i in range(4):
        v = 0
        for j in range(4):
            v ^= gf_mul(int(MC_M[i, j]), int(col[j]))
        out[i] = v
    return out


def mc_inv_column(col: np.ndarray) -> np.ndarray:
    out = np.zeros(4, dtype=np.uint8)
    for i in range(4):
        v = 0
        for j in range(4):
            v ^= gf_mul(int(MC_INV[i, j]), int(col[j]))
        out[i] = v
    return out


def find_best_3round_value_trail(delta_in_value: int = 0x80,
                                 cancel_row: int = 0,
                                 verbose: bool = False):
    """For 3-round AES with diagonal input (all 4 active bytes have value
    `delta_in_value`), find the value-level trail of maximum probability whose
    MC R1 output keeps only `cancel_row`-th byte non-zero.

    Returns a dict with the trail probability (log2), per-S-box (delta_in,
    delta_out, ddt_count) records, and the implied input difference (column-major hex).
    """
    # Step 1: the 4 active input bytes are at positions (i, i) (the diagonal).
    # After SR row i shifts left by i, so byte (i, i) goes to (i, (i - i) mod 4) = (i, 0).
    # So post-SR all 4 active bytes lie in column 0 of round 1's SR output.
    # Then MC R1 acts on column 0. We want MC * v = e_{cancel_row} * target with
    # target != 0, where v is the post-SB-R1 byte values at row i of column 0.

    # Step 2: enumerate target in [1, 255]; compute v = MC_inv * (target at cancel_row, 0s elsewhere)
    best = None
    rows_for_dout = DDT[delta_in_value]    # DDT row: counts of each delta_out given delta_in
    for target in range(1, 256):
        rhs = np.zeros(4, dtype=np.uint8)
        rhs[cancel_row] = target
        v = mc_inv_column(rhs)              # the 4 post-SB-R1 delta_outs
        # Check each v[i] is a valid DDT transition from delta_in_value
        counts = [int(rows_for_dout[int(v[i])]) for i in range(4)]
        if any(c == 0 for c in counts):
            continue
        # R1 contributes 4 active S-boxes. Trail prob factor = prod(counts) / 256^4
        # We'll also need R2 (1 S-box: input = target, after SR R2 it sits in some position)
        # and R3 (4 S-boxes: inputs are MC R2 output bytes).
        # For maximum trail probability, R2 and R3 pick max-DDT outputs.

        # R2: a single active byte with input = target. We need a valid output.
        r2_row = DDT[target]
        r2_max = int(r2_row.max())
        if r2_max == 0:
            continue
        r2_dout = int(np.argmax(r2_row))

        # After R2 SB: state has 1 active byte (at the position where R1 MC landed,
        # then permuted by R2 SR). Then MC R2 of that column produces 4 active bytes.
        # The actual byte values: column = [a, b, c, d]^T where exactly one entry is r2_dout
        # and the others are 0. Say r2_dout is at row k; then MC * e_k * r2_dout has 4 entries
        # = (MC[0,k]*r2_dout, MC[1,k]*r2_dout, MC[2,k]*r2_dout, MC[3,k]*r2_dout).
        # Each of these is the delta_in for an R3 S-box.

        # We need to know `k` = row index within R2's MC input column. That's determined
        # by which byte position the R1 MC output landed in, then through SR R2.
        # R1 MC output was at (cancel_row, 0). After SR R2 row r shifts left by r:
        # (cancel_row, 0) -> (cancel_row, -cancel_row mod 4) = (cancel_row, (4 - cancel_row) % 4)
        # That's the position entering R2 MC. It sits in column ((4 - cancel_row) % 4),
        # at row `cancel_row`. So within that column, k = cancel_row.
        k = cancel_row

        r3_inputs = [gf_mul(int(MC_M[i, k]), r2_dout) for i in range(4)]
        # R3: 4 S-boxes, each with the above inputs. Pick max-DDT output for each.
        r3_counts = []
        r3_douts = []
        ok = True
        for d_in in r3_inputs:
            row = DDT[d_in]
            m = int(row.max())
            if m == 0:
                ok = False; break
            r3_counts.append(m)
            r3_douts.append(int(np.argmax(row)))
        if not ok:
            continue

        # Trail probability (log2): sum of log2(count / 256) per active S-box.
        all_counts = counts + [r2_max] + r3_counts
        log2p = sum(np.log2(c / 256.0) for c in all_counts)
        active = len(all_counts)

        if best is None or log2p > best["log2p"]:
            best = {
                "log2p": float(log2p),
                "active": active,
                "target": target,
                "r1_v": [int(x) for x in v],
                "r1_counts": counts,
                "r2_dout": r2_dout, "r2_count": r2_max,
                "r3_inputs": r3_inputs,
                "r3_douts": r3_douts, "r3_counts": r3_counts,
            }

    return best


def assemble_input_delta(delta_in_value: int = 0x80) -> bytes:
    """The 4-byte diagonal input difference, value-level: bytes (0,0), (1,1), (2,2), (3,3) set."""
    out = bytearray(16)
    for i in range(4):
        out[4 * i + i] = delta_in_value  # column-major: state[i, i] -> index 4*i + i
    return bytes(out)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--delta-in", type=int, default=0x80,
                    help="Common value for the 4 diagonal input difference bytes.")
    args = ap.parse_args()

    print(f"Input difference: {assemble_input_delta(args.delta_in).hex()}")
    print(f"All 4 diagonal bytes set to 0x{args.delta_in:02x}")
    print()

    best_over_rows = None
    for row in range(4):
        b = find_best_3round_value_trail(args.delta_in, cancel_row=row)
        if b is None:
            print(f"cancel_row={row}: no valid value-level realization")
            continue
        print(f"cancel_row={row}: best log2(p) = {b['log2p']:.2f}, active = {b['active']}")
        if best_over_rows is None or b["log2p"] > best_over_rows["log2p"]:
            best_over_rows = b
            best_over_rows["row"] = row

    if best_over_rows is None:
        print("\nNo valid trail found.")
        return

    print()
    print(f"=== Best 3-round value-level trail (input value 0x{args.delta_in:02x}) ===")
    b = best_over_rows
    print(f"  log2(probability)   : {b['log2p']:.2f}")
    print(f"  active S-boxes      : {b['active']}")
    print(f"  cancellation target : 0x{b['target']:02x} at row {b['row']} of column 0")
    print(f"  R1 delta_outs       : {[hex(x) for x in b['r1_v']]}")
    print(f"  R1 DDT counts       : {b['r1_counts']}")
    print(f"  R2 delta_out        : 0x{b['r2_dout']:02x}  (DDT count {b['r2_count']})")
    print(f"  R3 delta_ins        : {[hex(x) for x in b['r3_inputs']]}")
    print(f"  R3 delta_outs       : {[hex(x) for x in b['r3_douts']]}")
    print(f"  R3 DDT counts       : {b['r3_counts']}")


if __name__ == "__main__":
    main()

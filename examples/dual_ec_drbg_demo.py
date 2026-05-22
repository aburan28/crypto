#!/usr/bin/env python3
"""
Dual EC DRBG backdoor demo (Shumow & Ferguson, 2007).

An attacker who chose Q with Q = d*P (knowing d) can recover the
generator's internal state from a single observed output block, and
predict all subsequent output.

Run:
  python3 examples/dual_ec_drbg_demo.py           # 8-bit truncation (~seconds)
  python3 examples/dual_ec_drbg_demo.py --full    # 16-bit truncation (slow)

NIST SP 800-90A truncates the top 16 bits of each 256-bit x-coordinate
output. The brute force is 2**bits_truncated curve-mul attempts, which is
expensive in pure Python -- the Rust/Sage demos in this directory do the
real 16-bit attack.
"""

import argparse
import secrets
import sys
import time

# NIST P-256
p  = 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
a  = -3 % p
b  = 0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b
n  = 0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
Gx = 0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296
Gy = 0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5
G  = (Gx, Gy)


def ec_add(P, Q):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0:
            return None
        m = (3 * x1 * x1 + a) * pow(2 * y1, -1, p) % p
    else:
        m = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (m * m - x1 - x2) % p
    y3 = (m * (x1 - x3) - y1) % p
    return (x3, y3)


def ec_mul(k, P):
    k %= n
    R, Q = None, P
    while k:
        if k & 1:
            R = ec_add(R, Q)
        Q = ec_add(Q, Q)
        k >>= 1
    return R


def sqrt_mod_p(c):
    # p ≡ 3 (mod 4): sqrt = c^((p+1)/4)
    r = pow(c, (p + 1) // 4, p)
    return r if (r * r) % p == c else None


def drbg_step(state, P, Q, mask):
    """One DRBG iteration. Returns (new_state, output_block)."""
    s_new = ec_mul(state, P)[0]
    r = ec_mul(s_new, Q)[0]
    return s_new, r & mask


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true",
                    help="use real 16-bit truncation (slow in Python)")
    ap.add_argument("--bits", type=int, default=None,
                    help="bits truncated from top of x-coord output")
    args = ap.parse_args()

    bits_truncated = args.bits if args.bits is not None else (16 if args.full else 8)
    output_bits = 256 - bits_truncated
    mask = (1 << output_bits) - 1
    candidates = 1 << bits_truncated

    print(f"=== Dual EC DRBG backdoor demo ===")
    print(f"truncating top {bits_truncated} bits -> {output_bits}-bit output blocks")
    print(f"brute-force search space: 2^{bits_truncated} = {candidates}\n")

    # 1. Attacker chooses backdoor d, publishes Q = d*P.
    d = secrets.randbelow(n - 1) + 1
    P = G
    Q = ec_mul(d, P)
    e = pow(d, -1, n)   # so P = e*Q

    # 2. Victim seeds DRBG with a secret state and emits two blocks.
    s0 = secrets.randbelow(n - 1) + 1
    s1, out1 = drbg_step(s0, P, Q, mask)
    s2, out2 = drbg_step(s1, P, Q, mask)
    print(f"victim seed (secret): {hex(s0)}")
    print(f"observed block 1:     {hex(out1)}")
    print(f"observed block 2:     {hex(out2)}\n")

    # 3. Attack. For each guess of the truncated top bits, reconstruct the
    #    full x-coordinate, decompress to a point R = s1*Q, then compute
    #    e*R = s1*(e*Q) = s1*P -- whose x-coord IS the next state s2.
    print("attacker brute-forcing top bits...")
    t0 = time.time()
    found = None
    tried = 0
    for top in range(candidates):
        x_full = (top << output_bits) | out1
        if x_full >= p:
            continue
        rhs = (x_full * x_full * x_full + a * x_full + b) % p
        y = sqrt_mod_p(rhs)
        if y is None:
            continue
        for cand_y in (y, p - y):
            tried += 1
            R = (x_full, cand_y)
            s2_guess = ec_mul(e, R)[0]      # = s1 * P  →  next state
            # Verify against block 2.
            pred_out2 = ec_mul(s2_guess, Q)[0] & mask
            if pred_out2 == out2:
                found = (top, s2_guess)
                break
        if found:
            break

    dt = time.time() - t0
    if found is None:
        print(f"[!] no candidate matched after {tried} tries", file=sys.stderr)
        sys.exit(1)

    top, s2_recovered = found
    print(f"[+] recovered after {tried} tries ({dt:.1f}s)")
    print(f"    top {bits_truncated} bits: {top:#x}")
    print(f"    recovered state s2: {hex(s2_recovered)}")
    print(f"    actual state s2:    {hex(s2)}")
    assert s2_recovered == s2

    # 4. Predict the next block from the recovered state.
    _, out3_actual = drbg_step(s2, P, Q, mask)
    _, out3_pred   = drbg_step(s2_recovered, P, Q, mask)
    print(f"\nactual block 3:     {hex(out3_actual)}")
    print(f"predicted block 3:  {hex(out3_pred)}")
    assert out3_actual == out3_pred
    print("\n[+] future output fully predictable from one observed block.")


if __name__ == "__main__":
    main()

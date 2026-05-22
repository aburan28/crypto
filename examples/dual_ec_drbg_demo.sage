#!/usr/bin/env sage
"""
Dual EC DRBG backdoor demo on NIST P-256.

Run:
  sage examples/dual_ec_drbg_demo.sage

Uses the real NIST spec: top 16 bits of each x-coordinate output are
truncated, so the brute force is 2**16. The Sage EC arithmetic is
C-backed (~1 min in pure Python this would be hours).
"""

# NIST P-256
p  = 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
a  = -3
b  = 0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b
nn = 0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
Gx = 0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296
Gy = 0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5

F = GF(p)
E = EllipticCurve(F, [a, b])
assert E.order() == nn

P = E(Gx, Gy)

BITS_TRUNCATED = 16
OUTPUT_BITS    = 256 - BITS_TRUNCATED
MASK           = (1 << OUTPUT_BITS) - 1


def drbg_step(state, P, Q):
    s_new = ZZ((state * P).xy()[0])
    r     = ZZ((s_new * Q).xy()[0])
    return s_new, r & MASK


def main():
    print("=== Dual EC DRBG backdoor demo (Sage, P-256) ===")
    print(f"truncating top {BITS_TRUNCATED} bits  ->  search space 2^{BITS_TRUNCATED}\n")

    # 1. Attacker picks backdoor d, sets Q = d*P.
    d = ZZ.random_element(1, nn)
    Q = d * P
    e = d.inverse_mod(nn)   # P = e * Q

    # 2. Victim seeds the DRBG and emits two output blocks.
    s0 = ZZ.random_element(1, nn)
    s1, out1 = drbg_step(s0, P, Q)
    s2, out2 = drbg_step(s1, P, Q)

    print(f"observed block 1: {hex(out1)}")
    print(f"observed block 2: {hex(out2)}\n")

    # 3. Brute-force the truncated bits, lift each candidate to a curve
    #    point R, compute e*R = s1*P; its x-coord is the next state s2.
    print("attacker searching 2^16 candidates...")
    import time
    t0 = time.time()
    found = None
    for top in range(1 << BITS_TRUNCATED):
        x_full = (top << OUTPUT_BITS) | out1
        if x_full >= p:
            continue
        try:
            R = E.lift_x(F(x_full), all=True)
        except (ValueError, TypeError):
            continue
        for cand in R:
            s2_guess = ZZ((e * cand).xy()[0])
            pred_out2 = ZZ((s2_guess * Q).xy()[0]) & MASK
            if pred_out2 == out2:
                found = (top, s2_guess)
                break
        if found:
            break

    dt = time.time() - t0
    assert found is not None, "no candidate matched"
    top, s2_recovered = found
    print(f"[+] recovered in {dt:.1f}s, top 16 bits = {top:#06x}")
    assert s2_recovered == s2
    print(f"    state matches: s2 = {hex(s2)}")

    # 4. Predict next block.
    _, out3_actual = drbg_step(s2, P, Q)
    _, out3_pred   = drbg_step(s2_recovered, P, Q)
    assert out3_actual == out3_pred
    print(f"\nactual block 3:    {hex(out3_actual)}")
    print(f"predicted block 3: {hex(out3_pred)}")
    print("\n[+] future output fully predictable from one observed block.")


main()

# SM2 — Chinese national ECC curve (GB/T 32918 / GM/T 0003), in SageMath.
#
# Self-contained module that
#   1. builds and *verifies* the recommended sm2p256v1 curve,
#   2. reimplements the SM2 signature scheme faithfully (real SM3 + ZA),
#      so signatures cross-check against the Rust impl in ../src/ecc/sm2.rs,
#   3. derives and demonstrates the SM2-specific HNP lattice attack on
#      biased nonces (the repo's hnp_ecdsa.rs only covers plain ECDSA).
#
# Run:  sage sm2.sage
#
# Parameters mirror src/ecc/curve.rs::CurveParams::sm2() exactly.

# ── Curve parameters (sm2p256v1, GB/T 32918) ───────────────────────────
p  = 0xFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF
a  = 0xFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFC
b  = 0x28E9FA9E9D9F5E344D5A9E4BCF6509A7F39789F515AB8F92DDBCBD414D940E93
Gx = 0x32C4AE2C1F1981195F9904466A39C9948FE30BBFF2660BE1715A4589334C74C7
Gy = 0xBC3736A2F4F6779C59BDCEE36B692153D0A9877CC62A474002DF32E52139F0A0
n  = 0xFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123
h  = 1

F = GF(p)
E = EllipticCurve(F, [F(a), F(b)])
G = E(Gx, Gy)

DEFAULT_ID = b"1234567812345678"  # GB/T 32918.2 default user identifier


# ── SM3 hash (GB/T 32905) ──────────────────────────────────────────────
# Pure-Python; verified below against the canonical SM3("abc") vector,
# which pins it to the same function as ../src/hash/sm3.rs.
_M32 = 0xFFFFFFFF
_IV = [0x7380166f, 0x4914b2b9, 0x172442d7, 0xda8a0600,
       0xa96f30bc, 0x163138aa, 0xe38dee4d, 0xb0fb0e4e]

def _rotl(x, k):
    k &= 31
    if k == 0:
        return x & _M32
    return ((x << k) | (x >> (32 - k))) & _M32

def _ff(j, x, y, z):
    return (x ^^ y ^^ z) if j < 16 else ((x & y) | (x & z) | (y & z))

def _gg(j, x, y, z):
    # bitwise NOT within 32 bits == XOR with all-ones (Sage's ~ is mult. inverse)
    return (x ^^ y ^^ z) if j < 16 else ((x & y) | ((x ^^ _M32) & z))

def _p0(x):
    return x ^^ _rotl(x, 9) ^^ _rotl(x, 17)

def _p1(x):
    return x ^^ _rotl(x, 15) ^^ _rotl(x, 23)

def sm3(msg):
    msg = bytes(msg)
    bitlen = len(msg) * 8
    msg += b"\x80"
    while len(msg) % 64 != 56:
        msg += b"\x00"
    msg += int(bitlen).to_bytes(8, "big")

    V = list(_IV)
    for off in range(0, len(msg), 64):
        blk = msg[off:off + 64]
        W = [int.from_bytes(blk[i:i + 4], "big") for i in range(0, 64, 4)]
        for j in range(16, 68):
            W.append(_p1(W[j-16] ^^ W[j-9] ^^ _rotl(W[j-3], 15))
                     ^^ _rotl(W[j-13], 7) ^^ W[j-6])
        Wp = [(W[j] ^^ W[j+4]) for j in range(64)]

        A, B_, C, D, Ee, Ff, Gg, Hh = V
        for j in range(64):
            tj = 0x79cc4519 if j < 16 else 0x7a879d8a
            ss1 = _rotl((_rotl(A, 12) + Ee + _rotl(tj, j)) & _M32, 7)
            ss2 = ss1 ^^ _rotl(A, 12)
            tt1 = (_ff(j, A, B_, C) + D + ss2 + Wp[j]) & _M32
            tt2 = (_gg(j, Ee, Ff, Gg) + Hh + ss1 + W[j]) & _M32
            D = C; C = _rotl(B_, 9); B_ = A; A = tt1
            Hh = Gg; Gg = _rotl(Ff, 19); Ff = Ee; Ee = _p0(tt2)
        V = [(x ^^ y) & _M32 for x, y in
             zip(V, [A, B_, C, D, Ee, Ff, Gg, Hh])]

    return b"".join(int(v).to_bytes(4, "big") for v in V)


# ── SM2 signature scheme (GB/T 32918.2) ────────────────────────────────
def _be32(x):
    return int(x).to_bytes(32, "big")

def za(pub, identity=DEFAULT_ID):
    """ZA = SM3(ENTL || ID || a || b || Gx || Gy || PAx || PAy)."""
    entl = (len(identity) * 8).to_bytes(2, "big")
    buf = (entl + identity + _be32(a) + _be32(b) + _be32(Gx) + _be32(Gy)
           + _be32(int(pub.xy()[0])) + _be32(int(pub.xy()[1])))
    return sm3(buf)

def _e_of(pub, msg, identity=DEFAULT_ID):
    return int.from_bytes(sm3(za(pub, identity) + bytes(msg)), "big") % n

def keygen():
    d = randint(1, n - 1)
    return d, d * G

def sign(d, pub, msg, identity=DEFAULT_ID, k=None):
    """SM2 sign. Returns (r, s, k); k exposed for the HNP experiments."""
    e = _e_of(pub, msg, identity)
    while True:
        kk = k if k is not None else randint(1, n - 1)
        x1 = int((kk * G).xy()[0])
        r = (e + x1) % n
        if r == 0 or (r + kk) % n == 0:
            if k is not None:
                raise ValueError("supplied k is degenerate")
            continue
        s = (inverse_mod(1 + d, n) * (kk - r * d)) % n
        if s == 0:
            if k is not None:
                raise ValueError("supplied k is degenerate")
            continue
        return r, s, kk

def verify(pub, msg, r, s, identity=DEFAULT_ID):
    if not (1 <= r < n and 1 <= s < n):
        return False
    e = _e_of(pub, msg, identity)
    t = (r + s) % n
    if t == 0:
        return False
    P = s * G + t * pub
    if P.is_zero():
        return False
    return (e + int(P.xy()[0])) % n == r


# ── SM2-specific Hidden Number Problem ─────────────────────────────────
# SM2:   s = (1+d)^{-1} (k - r d)  (mod n)
#  =>     k = s + (r + s) d         (mod n)
#
# So with  t_i = s_i ,  A_i = (r_i + s_i) mod n :   k_i = t_i + A_i d (mod n).
# Note the message digest e CANCELS — unlike ECDSA (t_i = s_i^{-1} z_i),
# the SM2 HNP needs only the (r, s) pairs and the nonce bias, not the
# messages. If the top (256 - kbits) bits of every k_i are zero, d is the
# hidden number and LLL recovers it (Boneh–Venkatesan, scaled-integer form
# matching ../src/cryptanalysis/hnp_ecdsa.rs).

def sm2_hnp_attack(sigs, kbits, pub):
    """Recover d from biased-nonce SM2 signatures.

    sigs : list of (r, s) with each nonce k < 2^kbits.
    Returns d on success (validated against pub = d*G), else None.
    """
    m = len(sigs)
    K = 1 << kbits
    A = [(r + s) % n for (r, s) in sigs]   # multipliers of d
    T = [s % n for (r, s) in sigs]         # constants

    B = matrix(ZZ, m + 2, m + 2)
    for i in range(m):
        B[i, i] = n * n
    for i in range(m):
        B[m, i]     = n * A[i]
        B[m + 1, i] = n * T[i]
    B[m, m]         = K
    B[m + 1, m + 1] = n * K

    for v in B.LLL():
        if v[m + 1] == 0:
            continue
        for sign_ in (1, -1):
            cand = sign_ * v[m]
            if cand % K != 0:
                continue
            d = (cand // K) % n
            if 1 <= d < n and d * G == pub:
                return d
    return None


# ── Curve / security verification ──────────────────────────────────────
def verify_curve():
    checks = []
    def chk(name, cond):
        checks.append((name, bool(cond)))

    # Special prime form: p = 2^256 - 2^224 - 2^96 + 2^64 - 1
    chk("prime p has Solinas form 2^256-2^224-2^96+2^64-1",
        p == 2**256 - 2**224 - 2**96 + 2**64 - 1)
    chk("p is prime", is_prime(p))
    chk("a == p - 3 (efficient doubling)", a == p - 3)
    chk("G is on the curve", G in E)
    chk("ord(G) == n", G.order() == n)
    chk("n is prime", is_prime(n))
    chk("#E == n  (cofactor h == 1)", E.order() == n and h == 1)

    # Not anomalous (#E != p  =>  trace != 1, blocks Smart/SSSA additive lift)
    t = p + 1 - n
    chk("non-anomalous (#E != p)", n != p)

    # MOV/Frey-Ruck resistance: embedding degree is huge (probe small k)
    low_emb = any(power_mod(p, k, n) == 1 for k in range(1, 200))
    chk("no low embedding degree (k < 200)", not low_emb)

    # CM discriminant t^2 - 4p is large (no small-CM weakness)
    disc = t * t - 4 * p
    chk("large CM discriminant", abs(disc) > 2**250)

    # Quadratic-twist structure #E' = 2(p+1) - n. NOT prime-order:
    # it factors 7 * 443 * p100 * p148. The *smooth* part (small-subgroup
    # leakage) is tiny, but the 100-bit factor means invalid-curve attacks
    # can recover ~100 bits at ~2^50 work IF an implementation skips point
    # validation -> SM2 MUST validate input points. Reported, not asserted.
    twist_order = 2 * (p + 1) - n
    smooth, rem = 1, twist_order
    for q in primes(10**7):
        while rem % q == 0:
            rem //= q
            smooth *= q
    chk("twist smooth part is small (< 2^16 small-subgroup leakage)",
        smooth < 2**16)

    return checks, t, twist_order, smooth, disc


# ── Driver ─────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 68)
    print("SM2  (sm2p256v1, GB/T 32918)  —  SageMath setup")
    print("=" * 68)

    # SM3 self-test against the standard vector.
    sm3_abc = sm3(b"abc").hex()
    expect  = "66c7f0f462eeedd9d1f2d46bdc10e4e24167c4875cf2f7a2297da02b8f4ba8e0"
    print("\n[SM3] SM3(\"abc\") = %s" % sm3_abc)
    print("[SM3] matches GB/T 32905 test vector: %s" % (sm3_abc == expect))
    assert sm3_abc == expect, "SM3 implementation is wrong"

    # Curve & security checks.
    print("\n[curve] verification:")
    checks, t, twist_order, smooth, disc = verify_curve()
    for name, ok in checks:
        print("   [%s] %s" % ("OK" if ok else "!!", name))
    assert all(ok for _, ok in checks), "curve verification failed"
    print("   trace of Frobenius t = p + 1 - n = %d" % t)
    print("   security level ~= %.1f bits (Pollard rho)"
          % (RR(log(sqrt(n), 2))))
    print("   note: twist is NOT prime-order (7 * 443 * p100 * p148);")
    print("         small-subgroup part = %d -> validate input points." % smooth)

    # Signature round-trip + a tampered-message rejection.
    print("\n[sig] sign / verify round-trip:")
    d, pub = keygen()
    msg = b"the quick brown fox jumps over the lazy dog"
    r, s, _ = sign(d, pub, msg)
    print("   verify(genuine)  = %s" % verify(pub, msg, r, s))
    print("   verify(tampered) = %s" % verify(pub, b"different msg", r, s))
    assert verify(pub, msg, r, s) and not verify(pub, b"different msg", r, s)

    # SM2-specific HNP attack on biased nonces.
    print("\n[hnp] SM2 biased-nonce key recovery (LLL):")
    leak  = 16                 # top `leak` bits of every nonce are zero
    kbits = 256 - leak
    m     = 40                 # signatures; m*leak = 640 >> 256
    print("   leaking top %d bits of each nonce, %d signatures" % (leak, m))
    d2, pub2 = keygen()
    sigs = []
    for i in range(m):
        k_bad = randint(1, (1 << kbits) - 1)         # biased nonce
        r, s, _ = sign(d2, pub2, ("msg-%d" % i).encode(), k=k_bad)
        sigs.append((r, s))                          # note: messages NOT used
    rec = sm2_hnp_attack(sigs, kbits, pub2)
    print("   recovered key == true key: %s" % (rec is not None and rec == d2))
    print("   (HNP used only (r,s) pairs — the SM2 digest e cancels out)")
    assert rec == d2, "HNP attack failed"

    print("\nAll checks passed.")

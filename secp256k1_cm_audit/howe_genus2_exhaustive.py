#!/usr/bin/env python3
"""
howe_genus2_exhaustive.py  --  Thread 22 (autolab 2026-07-29)

Companion to howe_c3_family_sweep.py.

That script found that 62 of 67 Howe-qualifying pairs of j=0 sextic twists are
realised by a genus-2 curve inside the zeta_3-equivariant family
y^2 = x^6 + a x^3 + c.  All 5 misses are *quadratic-twist* pairs (E_i, E_i^t),
whose Weil polynomial has c1 = 0:

    P(T) = (T^2 - tT + p)(T^2 + tT + p) = T^4 + (2p - t^2) T^2 + p^2

This is exactly the shape of the secp256k1 target pair (0,3) of Thread 18.
So the misses have to be explained: is the zeta_3 family simply too small, or
does the isogeny class of E x E^t contain no Jacobian at all?

This script settles that by searching ALL genus-2 curves over F_p.

COMPLETENESS OF THE MODEL
-------------------------
For p >= 7 every genus-2 curve C/F_p has a model

    y^2 = lc * (x^6 + b4 x^4 + b3 x^3 + b2 x^2 + b1 x + b0),   lc in {1, nu}

Reason: C is y^2 = F(x,z), F a binary sextic.  P^1(F_p) has p+1 >= 8 points and
F has at most 6 roots, so some F_p-rational point of P^1 is a non-root; move it
to infinity to get a degree-6 model with non-zero leading coefficient.  Scaling
y puts that coefficient in {1, nu} (nu a fixed non-residue), and translating x
kills the x^5 coefficient (needs 6 invertible, true for p >= 7).  Hence the
p^5 monic forms above, times the two twists, hit every isomorphism class.

SOUNDNESS OF A NEGATIVE ANSWER
------------------------------
Singular f (non-squarefree) are *not* filtered out: computing a degree-6
discriminant for ~10^8 quintuples is more expensive than the sweep itself.
This is safe in the direction that matters.  Singular members can only ADD
spurious (c1,c2) matches, never remove a genuine one, so "zero matches" is a
conclusive negative.  Any positive match is re-verified individually with PARI
hyperellcharpoly, which is defined only for nonsingular models.

METHOD
------
Same Weil-coefficient recovery as howe_c3_family_sweep.py (see its docstring):
c1 = -1 - S1 and c2 = (c1^2 + 1 + S2)/2, with S1, S2 the quadratic-character
sums over F_p and F_p^2.  Two accelerations:

 * The character argument is linear in (b4,b3,b2,b1,b0), so evaluating f over
   the whole field is a matrix product B @ V handled by BLAS in float64
   (all intermediates are < 6 p^2 << 2^53, so the arithmetic is exact).
 * Two-stage filter: S1 costs p evaluations per curve and S2 costs p^2.  We
   run the cheap S1 pass over all p^5 curves first, keep only those with the
   target c1 (a ~1/(4 sqrt p) fraction), and run the S2 pass on survivors only.

USAGE
-----
    python3 howe_genus2_exhaustive.py -p 7
    python3 howe_genus2_exhaustive.py -p 19 --pairs 2,5
    python3 howe_genus2_exhaustive.py -p 7,19 --all-qualifying
"""

import argparse
import subprocess
import sys

import numpy as np

from howe_c3_family_sweep import (
    howe_pairs,
    legendre_table,
    non_residue,
    primitive_root,
    sextic_twists,
)


# --------------------------------------------------------------------------
# field-evaluation tables
# --------------------------------------------------------------------------


def fp_tables(p):
    """Powers x^6 and [x^4,x^3,x^2,x,1] for x in F_p."""
    x = np.arange(p, dtype=np.int64)
    x2 = x * x % p
    x3 = x2 * x % p
    x4 = x3 * x % p
    x6 = x3 * x3 % p
    V = np.stack([x4, x3, x2, x, np.ones(p, dtype=np.int64)]).astype(np.float64)
    return x6.astype(np.int64), V


def fp2_tables(p, r):
    """Powers of X over F_p^2 = F_p[t]/(t^2-r), as (real, imag) tables.

    Returns (X6r, X6i, Vr, Vi) with Vr/Vi the (5, p^2) tables of
    [X^4, X^3, X^2, X, 1] split into real and imaginary parts.
    """
    a = np.repeat(np.arange(p, dtype=np.int64), p)  # real part
    b = np.tile(np.arange(p, dtype=np.int64), p)  # imag part

    def mul(u, v):
        ur, ui = u
        vr, vi = v
        return ((ur * vr + r * (ui * vi % p)) % p, (ur * vi + ui * vr) % p)

    X1 = (a, b)
    X2 = mul(X1, X1)
    X3 = mul(X2, X1)
    X4 = mul(X3, X1)
    X6 = mul(X3, X3)
    one_r = np.ones(p * p, dtype=np.int64)
    one_i = np.zeros(p * p, dtype=np.int64)
    Vr = np.stack([X4[0], X3[0], X2[0], X1[0], one_r]).astype(np.float64)
    Vi = np.stack([X4[1], X3[1], X2[1], X1[1], one_i]).astype(np.float64)
    return X6[0], X6[1], Vr, Vi


# --------------------------------------------------------------------------
# the two passes
# --------------------------------------------------------------------------


def index_to_coeffs(p, idx, family):
    """Decode a flat index into the coefficient matrix (n, 5) = (b4,b3,b2,b1,b0).

    family="full": idx ranges over p^5, all monic sextics with no x^5 term.
    family="even": idx ranges over p^3, the even sextics x^6+b4x^4+b2x^2+b0
                   (b3 = b1 = 0).  These are the bielliptic ones: u = x^2 maps
                   C onto the elliptic curve y^2 = u^3 + b4 u^2 + b2 u + b0.
    """
    B = np.zeros((idx.size, 5), dtype=np.float64)
    tmp = idx.copy()
    slots = [0, 1, 2, 3, 4] if family == "full" else [0, 2, 4]
    for k in reversed(slots):
        B[:, k] = tmp % p
        tmp //= p
    return B


def family_size(p, family):
    return p**5 if family == "full" else p**3


def s1_pass(p, target_c1, family="full", chunk_rows=200_000, progress=True):
    """All curves in the family whose monic model has c1 == +-target_c1.

    Also returns the leading coefficient realising the sign: a curve matches
    the target as lc=1 if c1 == target_c1 and as lc=nu if c1 == -target_c1.
    Both are collected (with a flag) since c2 is twist-invariant.
    """
    chi = legendre_table(p).astype(np.int64)
    x6, V = fp_tables(p)
    total = family_size(p, family)

    keep_idx = []
    keep_lc = []
    for start in range(0, total, chunk_rows):
        end = min(start + chunk_rows, total)
        idx = np.arange(start, end, dtype=np.int64)
        B = index_to_coeffs(p, idx, family)
        vals = (B @ V) % p  # (nchunk, p)
        vals = (vals.astype(np.int64) + x6[None, :]) % p
        S1 = chi[vals].sum(axis=1)
        c1 = -1 - S1
        m_plus = c1 == target_c1
        m_minus = c1 == -target_c1
        if m_plus.any():
            keep_idx.append(idx[m_plus])
            keep_lc.append(np.ones(int(m_plus.sum()), dtype=np.int64))
        if target_c1 != 0 and m_minus.any():
            keep_idx.append(idx[m_minus])
            keep_lc.append(np.full(int(m_minus.sum()), non_residue(p), dtype=np.int64))
        if progress and (start // chunk_rows) % 25 == 0:
            done = end / total
            print(f"    S1 pass {done:6.1%}", end="\r", flush=True)
    if progress:
        print("    S1 pass 100.0%  ", flush=True)

    if not keep_idx:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
    return np.concatenate(keep_idx), np.concatenate(keep_lc)


def s2_pass(p, idx, target_c1, family="full", budget=4_000_000, progress=True):
    """c2 of each survivor from the S1 pass (c2 is twist-invariant).

    The working arrays are (chunk_rows, p^2), so chunk_rows is chosen from an
    element budget rather than fixed: at p=127 a fixed 20k rows would be
    20000*127^2 = 3.2e8 entries (~2.6 GB per temporary).
    """
    chunk_rows = max(1, budget // (p * p))
    chi = legendre_table(p).astype(np.int64)
    r = non_residue(p)
    X6r, X6i, Vr, Vi = fp2_tables(p, r)

    out = np.empty(idx.size, dtype=np.int64)
    for start in range(0, idx.size, chunk_rows):
        end = min(start + chunk_rows, idx.size)
        sub = idx[start:end]
        B = index_to_coeffs(p, sub, family)
        R = ((B @ Vr) % p).astype(np.int64)
        I = ((B @ Vi) % p).astype(np.int64)
        R = (R + X6r[None, :]) % p
        I = (I + X6i[None, :]) % p
        nrm = (R * R - r * (I * I % p)) % p
        S2 = chi[nrm].sum(axis=1)
        out[start:end] = (target_c1 * target_c1 + 1 + S2) // 2
        if progress and (start // chunk_rows) % 50 == 0:
            print(f"    S2 pass {end / max(idx.size,1):6.1%}", end="\r", flush=True)
    if progress:
        print("    S2 pass 100.0%  ", flush=True)
    return out


def decode(p, i, family="full"):
    """Index -> (b4,b3,b2,b1,b0)."""
    B = index_to_coeffs(p, np.array([i], dtype=np.int64), family)
    return tuple(int(v) for v in B[0])


def pari_verify(p, lc, coeffs, want_c1, want_c2):
    """Re-check one hit with PARI hyperellcharpoly (rejects singular models)."""
    b4, b3, b2, b1, b0 = coeffs
    f = f"{lc}*(x^6+{b4}*x^4+{b3}*x^3+{b2}*x^2+{b1}*x+{b0})"
    script = (
        f"f = Mod(1,{p})*({f});"
        f"if(poldisc(f)==0, print(\"SINGULAR\"),"
        f'  cp=hyperellcharpoly(f); print(-polcoef(cp,3)," ",polcoef(cp,2)));'
    )
    out = subprocess.run(
        ["gp", "-q"], input=script, capture_output=True, text=True, timeout=300
    )
    lines = [ln for ln in out.stdout.strip().splitlines() if ln.strip()]
    if not lines:
        return "PARI-ERROR", None
    last = lines[-1].strip()
    if "SINGULAR" in last:
        return "SINGULAR", None
    got = tuple(int(v) for v in last.split())
    return ("MATCH" if got == (want_c1, want_c2) else "MISMATCH"), got


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------


def run(p, wanted_pairs, max_verify=5, family="full"):
    tw = sextic_twists(p)
    pairs = {tuple(r["pair"]): r for r in howe_pairs(tw, p)}

    label = "EXHAUSTIVE (all genus-2)" if family == "full" else "EVEN/bielliptic family"
    print(f"===== {label} search, p = {p} "
          f"({family_size(p, family)} monic sextics x 2 twists) =====")
    for e in tw:
        typ = "II(split)" if e["tors2"] == (1, 1, 1) else "I (irred)"
        print(f"    k={e['k']}  b={e['b']:4d}  t={e['t']:5d}  #E={e['n']:5d}  "
              f"2-tors {typ}")

    for pr in wanted_pairs:
        row = pairs[pr]
        tc1, tc2 = row["c1"], row["c2"]
        print(f"\n  --- pair {pr}: Howe={row['ok']}  target (c1,c2)=({tc1},{tc2}) ---")
        idx, lc = s1_pass(p, tc1, family=family)
        print(f"    S1 survivors (c1 = +-{tc1}): {idx.size}")
        neg = ("RESULT: NO genus-2 curve over F_p has this Weil polynomial."
               if family == "full" else
               "RESULT: none in this family (search NOT exhaustive - "
               "no conclusion about other genus-2 curves).")
        if idx.size == 0:
            print(f"    {neg}")
            continue
        c2 = s2_pass(p, idx, tc1, family=family)
        hit = np.nonzero(c2 == tc2)[0]
        print(f"    full matches (c1 and c2): {hit.size}")
        if hit.size == 0:
            print(f"    {neg}")
            continue
        print(f"    RESULT: {hit.size} model(s) match; verifying with PARI:")
        nverified = 0
        for h in hit:
            if nverified >= max_verify:
                break
            co = decode(p, int(idx[h]), family)
            status, got = pari_verify(p, int(lc[h]), co, tc1, tc2)
            print(f"      lc={int(lc[h])} (b4,b3,b2,b1,b0)={co}  -> {status} {got or ''}")
            if status != "SINGULAR":
                nverified += 1
        if nverified == 0:
            print("    NOTE: every inspected match was a SINGULAR model; "
                  "inspect more (raise --max-verify) before concluding.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-p", "--primes", required=True)
    ap.add_argument("--pairs", default=None,
                    help="semicolon-separated pairs, e.g. '2,5;1,4'")
    ap.add_argument("--all-qualifying", action="store_true")
    ap.add_argument("--max-verify", type=int, default=5)
    ap.add_argument("--family", choices=["full", "even"], default="full",
                    help="'full' = all p^5 genus-2 models (exhaustive); "
                         "'even' = the p^3 bielliptic even sextics (cheap, "
                         "not exhaustive)")
    args = ap.parse_args()

    for p in [int(s) for s in args.primes.split(",")]:
        if p < 7 or p % 6 != 1:
            print(f"skip p={p}: need p >= 7 and p = 1 mod 6")
            continue
        tw = sextic_twists(p)
        rows = howe_pairs(tw, p)
        if args.pairs:
            want = [tuple(int(v) for v in s.split(",")) for s in args.pairs.split(";")]
        elif args.all_qualifying:
            want = [tuple(r["pair"]) for r in rows if r["ok"]]
        else:
            # default: the quadratic-twist pairs among the qualifying ones
            want = [
                tuple(r["pair"])
                for r in rows
                if r["ok"] and r["c1"] == 0
            ]
        if not want:
            print(f"p={p}: no pairs selected")
            continue
        run(p, want, max_verify=args.max_verify, family=args.family)
        print()


if __name__ == "__main__":
    sys.exit(main())

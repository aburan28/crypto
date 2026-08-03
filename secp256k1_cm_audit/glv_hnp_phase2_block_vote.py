"""
Thread 23, exploit follow-up: per-block corruption and majority-vote readout.

glv_hnp_phase2_exact_cvp.py produced an apparent contradiction: at eff=0.05,
12/50 instances have a lattice point STRICTLY CLOSER to the target than the
true nonce offset, yet exact CVP still recovered d on 50/50.

Resolution: the rival points differ from the true offset by a lambda-block
vector w supported on a single column pair (i, m+i).  Adding w corrupts
(k1_i, k2_i) but leaves every other block untouched.  The readout used so far
reads d from signature 0 only, so it survives iff block 0 happens to be clean.

Two consequences to test here:

  F1  How many of the m blocks does the returned point actually corrupt?
      If the answer is "one or two", the failures are near-misses, not noise.

  F2  Majority-vote readout: compute
          d_i = (k1_i + lam*k2_i - A_i) * B_i^{-1}   (mod n)
      for EVERY i and take the modal value.  This recovers d whenever a
      plurality of blocks are clean, instead of requiring block 0 to be clean.

Note on the baseline: the original arm-O lattice carries d in its own column,
which a lambda-block vector never touches, so arm O is already immune to
block corruption.  That is why arm O matched exact CVP in the previous
experiment.  Majority voting is what gives the d-free formulations the same
immunity -- and, if corruption is sparse enough, more.

Run: python3 glv_hnp_phase2_block_vote.py
"""

import math
import random
from collections import Counter

import sympy
from fpylll import CVP, IntegerMatrix, LLL

import glv_hnp_phase2_projected as T23
from glv_hnp_phase2_projected import (
    scales, gen_signatures, norm, modinv, find_generator, build_curve,
    eisenstein_decompose, j0_traces, glv_roots,
    build_lattice_C, target_C, babai_nearest_plane, lll_rows,
)
from glv_hnp_phase2_exact_cvp import (
    l0_basis, exact_cvp_offset, search_curves,
)

SEEDS = [42, 1234, 9999, 555, 31337]
M_SIGS = 10


def d_candidates(offset, sigs, n, lam, k1b, k2b):
    """d_i from every block i; None for blocks whose entries are not on the
    scale lattice (those are corrupted beyond a clean (k1_i,k2_i) reading)."""
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1b, k2b)
    out = []
    for i in range(m):
        a, c = offset[i], offset[m + i]
        if a % S_K1 or c % S_K2:
            out.append(None)
            continue
        k1_i, k2_i = a // S_K1, c // S_K2
        try:
            Binv = modinv(sigs[i]['B'] % n, n)
        except ValueError:
            out.append(None)
            continue
        out.append((k1_i + lam * k2_i - sigs[i]['A']) % n * Binv % n)
    return out


def modal_d(cands):
    vals = [c for c in cands if c is not None and c != 0]
    if not vals:
        return None, 0
    v, cnt = Counter(vals).most_common(1)[0]
    return v, cnt


print("=" * 78)
print("Thread 23 -- per-block corruption and majority-vote readout")
print("=" * 78)

CURVES = search_curves(2**16, 2**17, want=10)
print(f"\n{len(CURVES)} fresh 17-bit j=0 GLV curves, m={M_SIGS}, "
      f"L0 dimension {2*M_SIGS}")

print("\n" + "-" * 78)
print("F1/F2 -- clean-block count and majority vote, exact-CVP offsets")
print("-" * 78)
print("clean  = blocks whose d_i equals the true d")
print("O      = historical arm (LLL, d read from its own column)")
print("X@0    = exact CVP, d read from block 0 only  (previous experiment)")
print("X@maj  = exact CVP, modal d over all m blocks (this experiment)")
print()
print(f"{'eff':>5} | {'O':>7} {'X@0':>7} {'X@maj':>7} | "
      f"{'mean clean':>10} {'median clean':>12} {'all-m clean':>12}")
print("-" * 78)

rows = []
for eff in (0.05, 0.10, 0.15, 0.20, 0.25, 0.35, 0.50):
    wo = wx0 = wxm = tot = 0
    cleans = []
    allclean = 0
    for (p, b, n, lam, G) in CURVES:
        k2b = math.isqrt(n) + 1
        k1b = max(2, int(eff * n / k2b))
        for seed in SEEDS:
            d_trial = random.Random(seed + 7777).randint(1, n - 1)
            sigs = gen_signatures(G, d_trial, M_SIGS, n, lam, p, k1b, k2b, seed)
            if len(sigs) < M_SIGS:
                continue
            tot += 1

            res = T23.trial((p, b, n, lam, G), M_SIGS, d_trial, k1b, seed,
                            arms=('O',))
            wo += bool(res['O'][0])

            basis = l0_basis(sigs, n, lam, k1b, k2b)
            tgt = target_C(sigs, n, k1b, k2b)
            off = exact_cvp_offset(basis, tgt)

            best = None
            for sgn in (-1, 1):
                cands = d_candidates([sgn * x for x in off], sigs, n, lam,
                                     k1b, k2b)
                nclean = sum(1 for c in cands if c == d_trial)
                mv, _ = modal_d(cands)
                if best is None or nclean > best[0]:
                    best = (nclean, cands[0] == d_trial, mv == d_trial)
            nclean, ok0, okm = best
            cleans.append(nclean)
            allclean += (nclean == M_SIGS)
            wx0 += ok0
            wxm += okm

    cleans.sort()
    mean_c = sum(cleans) / len(cleans) if cleans else float('nan')
    med_c = cleans[len(cleans) // 2] if cleans else float('nan')
    rows.append((eff, wo, wx0, wxm, tot, mean_c, med_c, allclean))
    print(f"{eff:>5.2f} | {wo:>3}/{tot:<3} {wx0:>3}/{tot:<3} {wxm:>3}/{tot:<3} | "
          f"{mean_c:>10.2f} {med_c:>12} {allclean:>7}/{tot:<4}")

print("\n" + "-" * 78)
print("F3 -- verdict")
print("-" * 78)
to = sum(r[1] for r in rows); t0 = sum(r[2] for r in rows)
tm = sum(r[3] for r in rows); tt = sum(r[4] for r in rows)
print(f"  totals over {tt} instances:")
print(f"    arm O   (LLL, direct d column) : {to}/{tt}")
print(f"    arm X@0 (exact CVP, block 0)   : {t0}/{tt}")
print(f"    arm X@maj (exact CVP, modal d) : {tm}/{tt}")
print()
print("  If X@maj > X@0 but X@maj ~= O, majority voting only restores the")
print("  immunity arm O already had, and Phase 2 gains nothing.")
print("  If X@maj > O, the block structure is genuinely exploitable and the")
print("  historical readout was leaving recoveries on the table.")

print("\n" + "=" * 78)
print("done")
print("=" * 78)

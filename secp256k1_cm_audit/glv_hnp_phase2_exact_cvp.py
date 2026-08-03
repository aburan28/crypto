"""
Thread 23, decisive experiment: is the GLV-HNP Phase-2 wall a
reduction-quality wall or an information-theoretic one?

Established earlier in this session (glv_hnp_phase2_projected*.py):

  * Projecting out the d-column removes the trivial shortest vector
    n*S_D*e_m but changes recovery on 0 of 300 instances -- the two arms
    succeed on exactly the same instance set.
  * Whenever the planted vector is NOT lambda_1 of the projected lattice, the
    actual lambda_1 has Kannan coordinate 0 (33/33 such instances).  Those
    rivals are filtered out at readout by the |last| == S_KANNAN test, which is
    why removing them cannot help.

So "make the planted vector lambda_1" is the wrong target: recovery is a coset
condition, i.e. a CVP/BDD instance in the dim-2m sublattice
    L0 = { v in L' : v_kannan = 0 },
with target t = (-A_i*S_K1 | 0) and true offset (k1_i*S_K1 | k2_i*S_K2).

This script replaces every approximation with the exact answer:

  arm O   original lattice, LLL              (historical baseline)
  arm C   L0 + Babai nearest-plane           (approximate CVP)
  arm X   L0 + fpylll enumeration            (EXACT closest vector)

Decision rule:
  * arm X recovers where O/C fail   -> the wall is reduction quality; better
    reduction (BKZ, sieving) buys real ground and Phase 2 is not at its ceiling.
  * arm X also fails, AND the returned lattice point is strictly closer to the
    target than the true offset -> the wall is INFORMATION-THEORETIC: the true
    nonce pair is not the closest point, so no lattice algorithm whatsoever
    recovers it from this formulation.

Run: python3 glv_hnp_phase2_exact_cvp.py
"""

import math
import random

import sympy
from fpylll import CVP, IntegerMatrix, LLL

import glv_hnp_phase2_projected as T23
from glv_hnp_phase2_projected import (
    scales, gen_signatures, norm, lam_star, find_generator, build_curve,
    eisenstein_decompose, j0_traces, glv_roots,
    build_lattice_C, target_C, babai_nearest_plane, d_from_k,
    lll_rows,
)

SEEDS = [42, 1234, 9999, 555, 31337]
M_SIGS = 10          # dim 2m = 20 for the enumeration


def l0_basis(sigs, n, lam, k1b, k2b):
    """Full-rank LLL-reduced basis of L0 (the rank-2m Kannan-free sublattice)."""
    m = len(sigs)
    rows = build_lattice_C(sigs, n, lam, k1b, k2b)   # 2m+1 gens, rank 2m
    red = lll_rows(rows, 2 * m)
    red = [r for r in red if any(r)]
    assert len(red) == 2 * m, f"rank {len(red)} != {2*m}"
    return red


def exact_cvp_offset(basis, target):
    B = IntegerMatrix.from_matrix(basis)
    LLL.reduction(B)
    v = CVP.closest_vector(B, tuple(target))
    return [t - x for t, x in zip(target, v)]


def check_offset(offset, sigs, n, lam, k1b, k2b, d_secret):
    m = len(sigs)
    S_K1, S_D, S_K2, S_KAN = scales(n, k1b, k2b)
    for sgn in (-1, 1):
        k1_0, k2_0 = sgn * offset[0], sgn * offset[m]
        if k1_0 % S_K1 or k2_0 % S_K2:
            continue
        d_cand = d_from_k(k1_0 // S_K1, k2_0 // S_K2, sigs[0], n, lam)
        if d_cand and d_cand == d_secret:
            return True
    return False


def search_curves(lo, hi, want=10):
    out = []
    p = int(sympy.nextprime(lo))
    while p < hi and len(out) < want:
        if p % 3 == 1:
            eis = eisenstein_decompose(p)
            if eis is not None:
                a_e, b_e = eis
                for t in j0_traces(a_e, b_e):
                    nc = p + 1 - t
                    if nc < 2 or nc % 3 != 1 or not sympy.isprime(nc):
                        continue
                    if glv_roots(nc) is None:
                        continue
                    cur = build_curve(p, nc)
                    if cur is None:
                        continue
                    out.append((p, cur[1], nc, glv_roots(nc)[0], cur[4]))
                    break
        p = int(sympy.nextprime(p))
    return out



def main():
    print("=" * 78)
    print("Thread 23 -- exact CVP: reduction-quality wall or information-theoretic?")
    print("=" * 78)


    CURVES = search_curves(2**16, 2**17, want=10)
    print(f"\n{len(CURVES)} fresh 17-bit j=0 GLV curves, m={M_SIGS} signatures, "
          f"L0 dimension {2*M_SIGS}")

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("E1 -- recovery by arm, and WHY exact CVP fails when it fails")
    print("-" * 78)
    print("d_true / d_cvp : distance from target to the true offset vs. to the")
    print("exact closest lattice point.  ratio < 1 means the true nonce pair is")
    print("NOT the closest point -> information-theoretic failure.")
    print()
    print(f"{'eff':>5} | {'arm O':>8} {'arm C':>8} {'arm X':>8} | "
          f"{'X closer-than-true':>19} {'mean d_cvp/d_true':>18}")
    print("-" * 78)

    summary = []
    for eff in (0.05, 0.10, 0.15, 0.20, 0.25, 0.35, 0.50):
        wo = wc = wx = tot = 0
        strictly_closer = 0
        ratios = []
        for (p, b, n, lam, G) in CURVES:
            k2b = math.isqrt(n) + 1
            k1b = max(2, int(eff * n / k2b))
            S_K1, S_D, S_K2, S_KAN = scales(n, k1b, k2b)
            for seed in SEEDS:
                d_trial = random.Random(seed + 7777).randint(1, n - 1)
                sigs = gen_signatures(G, d_trial, M_SIGS, n, lam, p, k1b, k2b, seed)
                if len(sigs) < M_SIGS:
                    continue
                tot += 1
                m = M_SIGS

                res = T23.trial((p, b, n, lam, G), M_SIGS, d_trial, k1b, seed,
                                arms=('O',))
                wo += bool(res['O'][0])

                basis = l0_basis(sigs, n, lam, k1b, k2b)
                tgt = target_C(sigs, n, k1b, k2b)

                off_b = babai_nearest_plane(basis, tgt)
                wc += check_offset(off_b, sigs, n, lam, k1b, k2b, d_trial)

                off_x = exact_cvp_offset(basis, tgt)
                wx += check_offset(off_x, sigs, n, lam, k1b, k2b, d_trial)

                true_off = [0] * (2 * m)
                for i in range(m):
                    true_off[i] = sigs[i]['k1'] * S_K1
                    true_off[m + i] = sigs[i]['k2'] * S_K2
                d_true, d_cvp = norm(true_off), norm(off_x)
                ratios.append(d_cvp / d_true if d_true else float('nan'))
                if d_cvp < d_true * (1 - 1e-12):
                    strictly_closer += 1

        mr = sum(ratios) / len(ratios) if ratios else float('nan')
        summary.append((eff, wo, wc, wx, tot, strictly_closer, mr))
        print(f"{eff:>5.2f} | {wo:>4}/{tot:<3} {wc:>4}/{tot:<3} {wx:>4}/{tot:<3} | "
              f"{strictly_closer:>10}/{tot:<8} {mr:>18.4f}")

    # ---------------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("E2 -- verdict")
    print("-" * 78)
    gained = [(e, wx - wo) for e, wo, wc, wx, t, sc, mr in summary]
    print("  exact-CVP gain over the historical LLL arm, per eff:")
    for e, g in gained:
        print(f"    eff={e:.2f}: {g:+d}")
    tot_true = sum(t for _, _, _, _, t, _, _ in summary)
    tot_sc = sum(sc for _, _, _, _, _, sc, _ in summary)
    tot_wx = sum(wx for _, _, _, wx, _, _, _ in summary)
    tot_wo = sum(wo for _, wo, _, _, _, _, _ in summary)
    print()
    print(f"  totals: arm O {tot_wo}/{tot_true}, arm X {tot_wx}/{tot_true}")
    print(f"  instances where a strictly closer lattice point than the true nonce")
    print(f"  pair exists: {tot_sc}/{tot_true} "
          f"({100.0 * tot_sc / tot_true:.1f}%)")
    print()
    print("  Reading: an exact-CVP failure with a strictly closer point is an")
    print("  information-theoretic failure of THIS lattice formulation -- no")
    print("  amount of BKZ/sieving can recover those instances.  An exact-CVP")
    print("  failure WITHOUT a closer point would be a bug, not a wall.")

    print("\n" + "=" * 78)
    print("done")
    print("=" * 78)


if __name__ == "__main__":
    main()

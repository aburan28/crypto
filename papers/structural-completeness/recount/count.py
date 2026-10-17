#!/usr/bin/env python3
"""Recount the instance totals the paper quotes in Section 4 (Remark after
Corollary B5 universality) and Remark rem:split-toy (Thread 12 sweep).

Uses the cited scripts' own point counter and primality test (imported from
secp256k1_cm_audit/) and the same enumeration as
exp_zprime2_generalized_lemma.residue_breakdown (l in {3,5,7}, 5 <= p < 500,
p != 0 mod l, b in 1..p-1, l | #E) and
exp_zdoubleprime_lemma_extended.residue_breakdown (l in {11,13}, p = 1 mod 3,
p != l).  Each instance is split by

  * curve type: y^2 = x^3 + b has j = 0, so it is ordinary iff p = 1 (mod 3)
    and supersingular iff p = 2 (mod 3) (p >= 5);
  * residue class: r = p mod l inside the lemma's hypothesis (r != l-1) or
    the excluded class r = l-1.

It also counts the primes visited by thread12_cm73_sweep.gp, whose loop is
forprime(p=7,5000, ...) with an early return when p % 6 != 1.

Run from the repository root:
    python3 -I -B papers/structural-completeness/recount/count.py
"""

import os
import sys
from collections import Counter

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "secp256k1_cm_audit"))

from exp_zprime2_generalized_lemma import count_points_j0, is_prime  # noqa: E402


def enumerate_instances(ells, max_p, ordinary_only, skip_p_eq_ell):
    c = Counter()
    for p in range(5, max_p):
        if not is_prime(p):
            continue
        ordinary = (p % 3 == 1)
        if ordinary_only and not ordinary:
            continue
        orders = None
        for ell in ells:
            if skip_p_eq_ell and p == ell:
                continue
            r = p % ell
            if r == 0:
                continue
            if orders is None:
                orders = [count_points_j0(b, p) for b in range(1, p)]
            for n_E in orders:
                if n_E % ell:
                    continue
                n_Et = 2 * (p + 1) - n_E
                cls = "excluded" if r == ell - 1 else "hypothesis"
                kind = "ordinary" if ordinary else "supersingular"
                c[(ell, kind, cls, "l|#E^t" if n_Et % ell == 0 else "l∤#E^t")] += 1
    return c


def report(title, c):
    print(title)
    for k in sorted(c):
        print("  l=%-2d %-13s %-10s %-7s %6d" % (k[0], k[1], k[2], k[3], c[k]))
    tot = sum(c.values())
    ordn = sum(v for k, v in c.items() if k[1] == "ordinary")
    ss = tot - ordn
    hyp_ord = sum(v for k, v in c.items() if k[1] == "ordinary" and k[2] == "hypothesis")
    hyp_ss = sum(v for k, v in c.items() if k[1] == "supersingular" and k[2] == "hypothesis")
    spurious = sum(v for k, v in c.items() if k[2] == "hypothesis" and k[3] == "l|#E^t")
    excl_fail = sum(v for k, v in c.items() if k[2] == "excluded" and k[3] == "l|#E^t")
    excl = sum(v for k, v in c.items() if k[2] == "excluded")
    print("  total instances                       %6d" % tot)
    print("    ordinary (p = 1 mod 3)              %6d" % ordn)
    print("    supersingular (p = 2 mod 3)         %6d" % ss)
    print("  in lemma hypothesis class (r != l-1)  %6d  (ordinary %d, supersingular %d)"
          % (hyp_ord + hyp_ss, hyp_ord, hyp_ss))
    print("    violations (l | #E^t)               %6d" % spurious)
    print("  in excluded class r = l-1             %6d  (l | #E^t in %d)" % (excl, excl_fail))
    print()
    return tot, ordn, hyp_ord, spurious


def main():
    a = enumerate_instances([3, 5, 7], 500, ordinary_only=False, skip_p_eq_ell=False)
    ta, oa, ha, sa = report("[1] exp_zprime2_generalized_lemma.py enumeration: "
                            "l in {3,5,7}, 5 <= p < 500, all j=0 curves", a)
    b = enumerate_instances([11, 13], 500, ordinary_only=True, skip_p_eq_ell=True)
    tb, ob, hb, sb = report("[2] exp_zdoubleprime_lemma_extended.py enumeration: "
                            "l in {11,13}, 5 <= p < 500, p = 1 mod 3", b)
    print("[3] combined l in {3,5,7,11,13}")
    print("  total instances                       %6d" % (ta + tb))
    print("  ordinary instances                    %6d" % (oa + ob))
    print("  ordinary, in lemma hypothesis class   %6d" % (ha + hb))
    print("  violations in hypothesis class        %6d" % (sa + sb))
    print()
    primes = [p for p in range(7, 5001) if is_prime(p)]
    swept = [p for p in primes if p % 6 == 1]
    print("[4] thread12_cm73_sweep.gp loop: forprime(p=7,5000) with p % 6 == 1")
    print("  primes visited by forprime            %6d" % len(primes))
    print("  primes p = 1 mod 6 processed          %6d  (first %d, last %d)"
          % (len(swept), swept[0], swept[-1]))
    print("  all primes p <= 5000                  %6d"
          % len([p for p in range(2, 5001) if is_prime(p)]))


if __name__ == "__main__":
    main()

#!/usr/bin/env sage -python
"""Independent check of every curve in results/*.txt with Sage.

For each file: q, N, ell prime; t^2 - 4q = -19 ell^2; every curve has
order exactly N (SEA, not the prototype's [N]P test); the crater has
j = -884736 and every floor curve does not; Frobenius discriminant is
-19 ell^2 on all of them; embedding degree exceeds 20.
Run:  sage -python verify_sage.py results/*.txt
"""
import re, sys
from sage.all import GF, EllipticCurve, is_prime, power_mod

ok_all = True
for fn in sys.argv[1:]:
    txt = open(fn).read()
    q = int(re.search(r"q = (\d+)", txt).group(1))
    t = int(re.search(r" t = (\d+)", txt).group(1))
    ell = int(re.search(r"ell = (\d+)", txt).group(1))
    N = int(re.search(r"N = #E = (\d+)", txt).group(1))
    F = GF(q)
    head = [is_prime(q), is_prime(N), is_prime(ell), t * t - 4 * q == -19 * ell * ell]
    print(f"{fn}: q={q} ({q.bit_length()} bits) ell={ell} N={N}  q,N,ell prime & disc: {head}")
    ok_all &= all(head)
    for m in re.finditer(r"(crater|floor) E(\d+): a = (\d+) b = (\d+) j = (\d+)", txt):
        level = m.group(1)
        E = EllipticCurve(F, [int(m.group(3)), int(m.group(4))])
        checks = {
            "order==N": E.order() == N,
            "j==reported": E.j_invariant() == F(int(m.group(5))),
            "j==-884736": (E.j_invariant() == F(-884736)) == (level == "crater"),
            "frob_disc==-19ell^2": E.frobenius_order().discriminant() == -19 * ell * ell,
            "embedding_degree>20": all(power_mod(q, k, N) != 1 for k in range(1, 21)),
        }
        ok_all &= all(checks.values())
        print(f"  {level} E{m.group(2)}: " + "  ".join(f"{k}={v}" for k, v in checks.items()))
print("ALL OK" if ok_all else "FAILURES PRESENT")
sys.exit(0 if ok_all else 1)

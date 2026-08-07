"""
GLV-HNP Phase 2, Thread 25c: the constrained product law.

Thread 25b fit an unconstrained exponent on the two coordinates and got

    recover  iff  NU * nu_hat^0.950  <  0.963

with the exponent landing on 0.950, i.e. indistinguishable from 1.  This
script tests the constrained form P = NU * nu_hat directly, brackets it the
way Thread 23/24 bracketed NU, and measures the correlation of each candidate
coordinate with the bias strength eff = K1*K2/n.

That last check is the point.  Thread 24's W5 nominated mu = lambda_1(L2) as
the second coordinate on the strength of per-eff-stratum AUCs.  If mu is
itself correlated with eff, then holding eff fixed is exactly the condition
under which mu and nu_hat cannot be told apart, and W5 could not have seen
the difference.  Measured here:

    rho(nu_hat, eff) = +0.005     rho(mu, eff) = -0.654    rho(NU, eff) = +0.721

so nu_hat is orthogonal to the bias axis and mu is not.  NU carries the
bias-strength axis; nu_hat is the independent curve-geometry axis.

Run: python3 glv_hnp_phase2_nuband_product.py [table.json]
"""

import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from glv_hnp_phase2_gsprofile import auc, spearman

if __name__ == "__main__":
    path = (sys.argv[1] if len(sys.argv) > 1
            else os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "glv_hnp_phase2_nuband_table.json"))
    rows = json.load(open(path))
    y = [r['ok'] for r in rows]

    print("Thread 25c — constrained product law  P = NU * nu_hat")
    print("=" * 70)
    for name, f in (("NU*nu_hat^1.000", lambda r: r['NU'] * r['nuhat']),
                    ("NU*nu_hat^0.950", lambda r: r['NU'] * r['nuhat'] ** 0.950),
                    ("NU*nu_hat^0.500", lambda r: r['NU'] * r['nuhat'] ** 0.5),
                    ("NU alone", lambda r: r['NU']),
                    ("nu_hat alone", lambda r: r['nuhat'])):
        P = [f(r) for r in rows]
        a = auc([p for p, t in zip(P, y) if t],
                [p for p, t in zip(P, y) if not t])
        print(f"  {name:<18} AUC {a:.4f}")

    P = [r['NU'] * r['nuhat'] for r in rows]
    acc, c = max((sum(1 for p, t in zip(P, y) if (p < cc) == t) / len(y), cc)
                 for cc in sorted(P))
    tp = sum(1 for p, t in zip(P, y) if p < c and t)
    fp = sum(1 for p, t in zip(P, y) if p < c and not t)
    fn = sum(1 for p, t in zip(P, y) if p >= c and t)
    tn = sum(1 for p, t in zip(P, y) if p >= c and not t)
    print(f"\n  best single threshold on NU*nu_hat : {c:.5f}  acc {acc:.3f}")
    print(f"  TP {tp} FP {fp} FN {fn} TN {tn}")

    pos = [p for p, t in zip(P, y) if t]
    neg = [p for p, t in zip(P, y) if not t]
    print(f"  P | success : min {min(pos):.4f} "
          f"median {sorted(pos)[len(pos)//2]:.4f} max {max(pos):.4f}")
    print(f"  P | failure : min {min(neg):.4f} "
          f"median {sorted(neg)[len(neg)//2]:.4f} max {max(neg):.4f}")
    print(f"  bracket     : sufficient P < {min(neg):.4f} , "
          f"necessary P > {max(pos):.4f}")

    print(f"\n  Spearman(nu_hat, eff) = "
          f"{spearman([r['nuhat'] for r in rows], [r['eff'] for r in rows]):+.4f}")
    print(f"  Spearman(NU, eff)     = "
          f"{spearman([r['NU'] for r in rows], [r['eff'] for r in rows]):+.4f}")
    print(f"  Spearman(mu, eff)     = "
          f"{spearman([r['mu'] for r in rows], [r['eff'] for r in rows]):+.4f}")
    print(f"  Spearman(NU, nu_hat)  = "
          f"{spearman([r['NU'] for r in rows], [r['nuhat'] for r in rows]):+.4f}")

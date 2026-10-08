# Amendment 1 rescoring of E1: the as-run exp.sage built the predicted relation
# with constant Tr(x3) instead of the registered Tr(b/x3^2). The linear part
# (sum_j Tr(v_j)(y1j + y2j)) was built correctly and is stored in predicted_L,
# and span(F) ∩ R_{<=1} is stored as deg_le1 (an echelon basis). This recomputes
# the registered constant from the stored x3, b and re-tests membership.
#   sage rescore_e1.sage exp/e1_*.jsonl
import sys, json
out = []
for path in sys.argv[1:]:
    for line in open(path):
        r = json.loads(line); n = r["n"]
        K = GF(2**n, 'z', modulus='minimal_weight'); z = K.gen()
        x3 = K(sage_eval(r["x3"], locals={'z': z})); b = K(sage_eval(r["b"], locals={'z': z}))
        c = int((b / x3**2).trace()) if x3 != 0 else None
        N = 2 * r["nprime"]
        B = BooleanPolynomialRing(N, 'x', order='deglex'); X = B.gens()
        loc = {('x%d' % i): X[i] for i in range(N)}
        lin = B(sage_eval(r["predicted_L"], locals=loc)) if r["predicted_L"] not in ("0", "1") else B(0)
        # strip the as-run constant, keep the linear part, add the registered constant
        lin = lin + B(lin.constant_coefficient()) if lin != 0 else B(0)
        L = lin + B(c)
        basis = [B(sage_eval(p, locals=loc)) for p in r["deg_le1"]]
        def rank(ps):
            ps = [p for p in ps if p != 0]
            if not ps: return 0
            mons = sorted({m for p in ps for m in p.monomials()}, reverse=True)
            ix = {m: i for i, m in enumerate(mons)}
            M = matrix(GF(2), len(ps), len(mons))
            for i, p in enumerate(ps):
                for m in p.monomials(): M[i, ix[m]] = 1
            return M.rank()
        in_span = (L == 0) or rank(basis + [L]) == rank(basis)
        one_in_span = rank(basis + [B(1)]) == rank(basis) and basis != []
        out.append(dict(n=n, family=r["family"], seed=r["seed"], trace_on_V_zero=r["trace_on_V_zero"],
                        c_registered=c, L_registered=str(L), L_in_span=bool(in_span),
                        one_in_span=bool(one_in_span), unsat=r["unsat"], ffd=r["ffd"],
                        d_last=r["d_last"], deg_le1_dim=len(basis)))
for o in out: print(json.dumps(o))

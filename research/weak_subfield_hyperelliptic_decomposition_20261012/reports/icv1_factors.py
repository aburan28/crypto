"""ICV1 identities of the elliptic factors E1, E2 cited by the E4/E4x records.

Usage: python3 icv1_factors.py <results.jsonl> <out.json>

The lab records E1 : v^2 = u^3 + a u^2 + b u + c and E2 : v^2 = c u^3 + b u^2 + a u + 1
over F_q (coefficients `abc_r`), measured over K = F_q[t]/(t^3 + m2 t^2 + m1 t + m0),
`modulus_low = [m0, m1, m2]`. ICV1 (docs/curves/ICV1.md) names the short Weierstrass
model over K, so each curve is moved to y^2 = x^3 + A x + B by an F_q-isomorphism
(E2 first scaled to monic: (U, V) = (c u, c v), V^2 = U^3 + b U^2 + a c U + c^2),
which preserves #E(K). The order is recomputed here by an independent point count
and must equal the lab's E1_K / E2_K. Nothing is computed from the lab's output other
than a, b, c, q, n, the modulus and those two counts.
"""
import json, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts"))
import curve_id as cid  # noqa: E402


def short(q, a2, a4, a6):
    """y^2 = x^3 + a2 x^2 + a4 x + a6 over F_q (q prime > 3) to y^2 = x^3 + A x + B."""
    i3, i27 = pow(3, -1, q), pow(27, -1, q)
    return (a4 - a2 * a2 * i3) % q, (a6 - a2 * a4 * i3 + 2 * a2 ** 3 * i27) % q


def main(src, out):
    rows = [json.loads(l) for l in open(src)]
    res, seen = [], set()
    for r in rows:
        if r["exp"] not in ("E4", "E4x"):
            continue
        q, n, mod = r["q"], r["n"], r["modulus_low"]
        a, b, c, _ = r["abc_r"]
        for name, (a2, a4, a6), order in (("E1", (a, b, c), r["E1_K"]),
                                          ("E2", (b, a * c, c * c), r["E2_K"])):
            A, B = short(q, a2, a4, a6)
            ea, eb = [A] + [0] * (n - 1), [B] + [0] * (n - 1)
            counted = cid._count_extension(q, mod, ea, eb)
            assert counted == order, (r["exp"], q, name, counted, order)
            ident = cid.extension_id(q, n, mod, ea, eb, order)
            res.append({"exp": r["exp"], "q": q, "n": n, "factor": name, "weierstrass_over_Fq": [a2, a4, a6],
                        "short": {"p": q, "k": n, "modulus_low": mod, "a": ea, "b": eb}, "order_K": order,
                        "order_recounted": counted, "slug": ident["slug"], "icv1": ident["icv1"],
                        "model_sha256": ident["model_sha256"]})
            seen.add(ident["slug"])
    Path(out).write_text(json.dumps(res, indent=1) + "\n")
    print(len(res), "factors,", len(seen), "distinct slugs")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

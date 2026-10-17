"""Screen every selected-point quartic line for residual multiplicities.

For a pivot P and line direction (1,m), precompute the Taylor coefficients
of F(P+t(1,m)) as short polynomials in m. Batch-invert all nonvertical
abscissa differences to obtain slopes with one field inverse per pivot.
Dividing by the two known point roots leaves a residual quadratic.
"""

import hashlib
import json
import time
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent
INPUTS = ("type_i_plane_quartic.json", "type_i_line_base_points.tsv",
          "type_i_line_new_points.tsv")
sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
plane = json.loads((HERE / INPUTS[0]).read_text())
p = int(plane["p"])
old = [tuple(map(int, row.split()[1:])) for row in
       (HERE / INPUTS[1]).read_text().splitlines()]
new = [tuple(map(int, row.split())) for row in
       (HERE / INPUTS[2]).read_text().splitlines()]
points = old + new
n = len(points)
assert len(old) == 141 and len(new) == 9784 and len(set(points)) == n == 9925
terms = [(tuple(map(int, key.split(","))), int(c) % p)
         for key, c in plane["quartic_coefficients"].items()]
assert all(sum(exponents) == 4 for exponents, _ in terms)


def taylor_at(x, y):
    """Return c_k(m) with F(x+t,y+mt,1) = sum_k c_k(m)t^k."""
    xp = [pow(x, i, p) for i in range(5)]
    yp = [pow(y, i, p) for i in range(5)]
    out = [[0] * (k + 1) for k in range(5)]
    for (a, b, _), value in terms:
        for ix in range(a + 1):
            for iy in range(b + 1):
                k = ix + iy
                out[k][iy] = (out[k][iy] + value * comb(a, ix)
                               * comb(b, iy) * xp[a - ix] * yp[b - iy]) % p
    assert out[0][0] == 0
    return out


def eval_poly(coeffs, m):
    z = 0
    for c in reversed(coeffs):
        z = (z * m + c) % p
    return z


counts = {"pairs": 0, "vertical_pairs": 0, "degree_below_four": 0,
          "residual_double_root": 0, "tangent_at_first": 0,
          "tangent_at_second": 0}
exceptions = []
exception_hash = hashlib.sha256()
start = time.monotonic_ns()


def inspect(i, j, d, c1, c2, c3, c4, vertical=False):
    counts["pairs"] += 1
    counts["vertical_pairs"] += vertical
    A = c4
    B = (c3 + d * A) % p
    C = (c2 + d * B) % p
    assert (c1 + d * C) % p == 0
    kinds = []
    if A == 0:
        counts["degree_below_four"] += 1
        kinds.append("degree_below_four")
    elif (B * B - 4 * A * C) % p == 0:
        counts["residual_double_root"] += 1
        kinds.append("residual_double_root")
    if C == 0:
        counts["tangent_at_first"] += 1
        kinds.append("tangent_at_first")
    if (A * d * d + B * d + C) % p == 0:
        counts["tangent_at_second"] += 1
        kinds.append("tangent_at_second")
    if kinds:
        row = [i, j, *kinds]
        exception_hash.update((json.dumps(row, separators=(",", ":")) + "\n").encode())
        if len(exceptions) < 1000:
            exceptions.append(row)


for i, (x, y) in enumerate(points):
    t = taylor_at(x, y)
    js, prefix, vertical = [], [], []
    product = 1
    for j in range(i + 1, n):
        dx = (points[j][0] - x) % p
        if dx:
            product = product * dx % p
            js.append(j)
            prefix.append(product)
        else:
            vertical.append(j)
    if js:
        inverse_product = pow(product, -1, p)
        for k in range(len(js) - 1, -1, -1):
            j = js[k]
            u, v = points[j]
            d = (u - x) % p
            inverse_dx = inverse_product * (prefix[k - 1] if k else 1) % p
            inverse_product = inverse_product * d % p
            m = (v - y) * inverse_dx % p
            inspect(i, j, d, eval_poly(t[1], m), eval_poly(t[2], m),
                    eval_poly(t[3], m), eval_poly(t[4], m))
        assert inverse_product == 1
    for j in vertical:
        d = (points[j][1] - y) % p
        assert d
        inspect(i, j, d, t[1][1], t[2][2], t[3][3], t[4][4], True)
    if (i + 1) % 250 == 0:
        print("FULL_RESIDUAL_DISCRIMINANT_PROGRESS", i + 1, n,
              counts, flush=True)

assert counts["pairs"] == n * (n - 1) // 2
out = {
    "status": "FULL_9925_POINT_RESIDUAL_DISCRIMINANT_CENSUS",
    "field_prime": str(p), "point_count": n, "counts": counts,
    "exception_witnesses_first_1000": exceptions,
    "exception_stream_sha256": exception_hash.hexdigest(),
    "source_sha256": sha(__file__),
    "input_sha256": {name: sha(HERE / name) for name in INPUTS},
    "elapsed_ns_unisolated": time.monotonic_ns() - start,
    "boundary": "Complete degree and residual-discriminant test for lines through two selected affine quartic points. If all exception counts vanish, every such cut has four distinct affine intersections; relation rank and target overlap remain separate.",
}
(HERE / "type_i_line_full_residual_discriminant.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n")
print("FULL_RESIDUAL_DISCRIMINANT_PASS", counts, flush=True)

#!/usr/bin/env python3
"""Complete full-point Boolean addition edge for y²+xy=x³+a*x²+b.

This keeps the frozen source-curve DAG untouched and returns its Relation
interface, so a later factor-domain chain can bind the same point roles.
"""
from __future__ import annotations

from pathlib import Path
import sys

PARENT = Path(__file__).resolve().parents[1] / "symbolic_dag_fullpoint_20260925"
sys.path.insert(0, str(PARENT))
from dag import Dag, Field, Relation  # noqa: E402


def build_relation(n: int, modulus: int, a: int, b: int) -> Relation:
    """Build an exact edge on a nonsingular ordinary binary curve.

    The branch cover is unchanged from the source-curve proof: for finite
    points with the same x, the curve equations imply y1+y2 is either 0 or
    x. The latter is the inverse branch; the former is doubling unless x=0,
    where the unique finite point is self-inverse. No field inversion appears
    in the Boolean circuit.
    """
    d = Dag()
    f = Field(d, n, modulus)
    curve_a = f.const(a)
    curve_b = f.const(b)
    if b == 0:
        raise ValueError("b=0 is a singular binary curve")

    points = []
    for label in ("p", "q", "r"):
        points.append((d.var(f"{label}_o"), f.var(f"{label}_x"), f.var(f"{label}_y")))
    lam = f.var("lambda")
    (po, px, py), (qo, qx, qy), (ro, rx, ry) = points

    def point_valid(point: tuple[int, tuple[int, ...], tuple[int, ...]]) -> int:
        o, x, y = point
        canonical_infinity = d.implies(o, d.and_(f.zero(x), f.zero(y)))
        lhs = f.add(f.square(y), f.mul(x, y))
        x_squared = f.square(x)
        rhs = f.add(f.add(f.mul(x_squared, x), f.mul(curve_a, x_squared)), curve_b)
        affine_on_curve = d.implies(d.not_(o), f.equal(lhs, rhs))
        return d.and_(canonical_infinity, affine_on_curve)

    validity = d.all_([point_valid(point) for point in points])
    finite = d.and_(d.not_(po), d.not_(qo))
    same_x = f.equal(px, qx)
    same_y = f.equal(py, qy)
    inverse = f.equal(f.add(py, qy), px)
    nonzero_x = d.not_(f.zero(px))
    branches = {
        "copy_q": po,
        "copy_p": d.and_(d.not_(po), qo),
        "inverse": d.all_([finite, same_x, inverse]),
        "double": d.all_([finite, same_x, d.not_(inverse), same_y, nonzero_x]),
        "generic": d.and_(finite, d.not_(same_x)),
    }
    cover = d.any_(list(branches.values()))

    def equal_points(left, right) -> int:
        return d.all_([
            d.eq(left[0], right[0]),
            f.equal(left[1], right[1]),
            f.equal(left[2], right[2]),
        ])

    copy_q = d.implies(branches["copy_q"], equal_points(points[2], points[1]))
    copy_p = d.implies(branches["copy_p"], equal_points(points[2], points[0]))
    inv_result = d.implies(branches["inverse"], d.and_(ro, d.and_(f.zero(rx), f.zero(ry))))

    double_slope = f.equal(f.mul(lam, px), f.add(f.square(px), py))
    double_x = f.add(f.add(f.square(lam), lam), curve_a)
    double_y = f.add(f.square(px), f.mul(f.add(lam, f.const(1)), rx))
    double_result = d.implies(branches["double"], d.all_([
        d.not_(ro), double_slope, f.equal(rx, double_x), f.equal(ry, double_y),
    ]))

    generic_slope = f.equal(f.mul(lam, f.add(px, qx)), f.add(py, qy))
    generic_x = f.add(f.add(f.add(f.square(lam), lam), f.add(px, qx)), curve_a)
    generic_y = f.add(f.add(f.mul(lam, f.add(px, rx)), rx), py)
    generic_result = d.implies(branches["generic"], d.all_([
        d.not_(ro), generic_slope, f.equal(rx, generic_x), f.equal(ry, generic_y),
    ]))

    output = d.all_([
        validity, cover, copy_q, copy_p, inv_result, double_result, generic_result,
    ])
    return Relation(d, f, output, branches)

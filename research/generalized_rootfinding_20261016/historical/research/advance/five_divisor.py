"""Exact five-source divisor support profiler, deliberately not a solver.

For E:y²+xy=x³+A*x²+B and T=-Q=(z,yT), a pole-six function is
F=P(x)+(h*x+k)y with P=x³+r*x²+s*x+t.  Set
t=z³+r*z²+s*z+(h*z+k)*yT, so F(T)=0.  Its norm, divided by X+z,
is monic H of degree five.  If its source abscissae are in a+V,
the elementary coefficients of H(X+a) lie in product spans W1,...,W5.

At fixed h, the remaining coefficients contain k*r and k*s; they are
not jointly binary-affine in r,s,k.  At fixed h,k the r,s constraints
are exactly binary-affine.  On n=83,t=18 prefix domains, W1..W5 have
dimensions 18,35,52,69,83.  There are at most 93 constraints on 166
conditional bits, giving nullity >=73 for any consistent family.
This backend exposes that obstruction instead of enumerating it.
"""

import json
from pathlib import Path

from rootfinder.binaryfield import BinaryField
from rootfinder import polynomial as poly
try:
    from .divisor_triple import prepare_model, _basis, _residual, _solve, enumerate_family
except ImportError:
    from divisor_triple import prepare_model, _basis, _residual, _solve, enumerate_family


def prepare_five_model(field, curve_a, curve_b, target, offset, basis):
    """Measure exact support ranks without enumerating field or domain."""
    model = prepare_model(field, curve_a, curve_b, target, offset, basis)
    spans = [model["basis"]]
    pivots = [model["domain_pivots"]]
    for _ in range(4):
        span, echelon = _basis(field.mul(left, right)
                               for left in spans[-1] for right in model["basis"])
        spans.append(span)
        pivots.append(echelon)
    model["product_spans"] = tuple(spans)
    model["product_pivots"] = tuple(pivots)
    model["support_dimensions"] = tuple(map(len, spans))
    model["conditional_constraint_ceiling"] = sum(field.degree - len(span)
                                                    for span in spans[1:])
    model["conditional_nullity_floor"] = max(
        0, 2 * field.degree - model["conditional_constraint_ceiling"])
    return model


def five_coefficients_formula(model, h, k, r, s):
    """Exact elementary coefficients, including curve A and T=-Q sign."""
    f = model["field"]
    h, k, r, s = map(f.validate, (h, k, r, s))
    z, yt = model["negative_target"]
    a = model["curve_a"]
    h2, k2, r2, s2, z2 = map(f.square, (h, k, r, s, z))
    z3, z4 = f.mul(z2, z), f.square(z2)
    e1 = h2 ^ h ^ z
    e2 = r2 ^ f.mul(h, r) ^ k ^ f.mul(h2, a) ^ f.mul(z, e1)
    e3 = f.mul(h, s) ^ f.mul(k, r) ^ k2 ^ f.mul(z, e2)
    e4 = (s2 ^ f.mul(k, s) ^ f.mul(z, f.mul(k, r)) ^ f.mul(z2, r2)
          ^ f.mul(k2, a ^ z) ^ f.mul(z2, k) ^ f.mul(h2, f.square(yt))
          ^ f.mul(f.mul(h, k), yt) ^ z4)
    e5 = (f.mul(z, s2) ^ f.mul(z3, r2) ^ f.mul(k2, yt ^ f.mul(z, a) ^ z2)
          ^ f.mul(z, f.mul(h2, f.square(yt))) ^ f.mul(z4, z))
    return e5, e4, e3, e2, e1, 1


def five_coefficients_norm(model, h, k, r, s):
    """Independent literal norm/division; cross-check the simplified formulas."""
    f = model["field"]
    h, k, r, s = map(f.validate, (h, k, r, s))
    z, yt = model["negative_target"]
    t = (f.mul(f.square(z), z) ^ f.mul(r, f.square(z)) ^ f.mul(s, z)
         ^ f.mul(f.mul(h, z) ^ k, yt))
    p, q = (t, s, r, 1), (k, h)
    norm = poly.add(poly.mul(p, p, f), poly.mul(poly.mul((0, 1), q, f), p, f), f)
    curve = (model["curve_b"], 0, model["curve_a"], 1)
    norm = poly.add(norm, poly.mul(poly.mul(q, q, f), curve, f), f)
    quotient, remainder = poly.divmod_poly(norm, (z, 1), f)
    if remainder or quotient != five_coefficients_formula(model, h, k, r, s):
        raise ArithmeticError("five-source coefficient formulas disagree with literal norm")
    return quotient


def shift_source_polynomial(model, coefficients):
    """H(X+a) has the source displacement roots in V, including offsets."""
    f, a = model["field"], model["offset"]
    result = ()
    for coefficient in reversed(coefficients):
        result = poly.add(poly.mul(result, (a, 1), f), (coefficient,), f)
    return result


def support_conditions(model, coefficients):
    shifted = shift_source_polynomial(model, coefficients)
    return tuple(_residual(shifted[5 - j], model["product_pivots"][j - 1]) == 0
                 for j in range(1, 6))


def conditional_five_family(model, h, k):
    """Exact r,s family after fixing h,k; never enumerate its cardinality."""
    f = model["field"]
    h, k = map(f.validate, (h, k))
    base = five_coefficients_formula(model, h, k, 0, 0)
    if not support_conditions(model, base)[0]:
        return {"consistent": False, "outer_admitted": False,
                "particular": None, "basis": (), "rank": 0,
                "nullity": 0, "cardinality": 0}

    def packed(coefficients):
        shifted, value = shift_source_polynomial(model, coefficients), 0
        for j in range(2, 6):
            value |= _residual(shifted[5 - j], model["product_pivots"][j - 1]) << (
                (j - 2) * f.degree)
        return value

    rhs = packed(base)
    columns = tuple(packed(five_coefficients_formula(model, h, k, unit, 0)) ^ rhs
                    for unit in model["units"])
    columns += tuple(packed(five_coefficients_formula(model, h, k, 0, unit)) ^ rhs
                     for unit in model["units"])
    result = _solve(columns, rhs)
    result.update({"outer_admitted": True, "variables": 2 * f.degree,
                   "constraint_rank_ceiling": model["conditional_constraint_ceiling"],
                   "nullity_floor": model["conditional_nullity_floor"]})
    return result


def profile_instance(path, outer_samples=2, d0_samples=(0, 1)):
    """Read an actual neutral instance and return untimed exact rank evidence."""
    data = json.loads(Path(path).read_text())
    spec, curve = data["extension_field"], data["curve"]
    f = BinaryField(spec["degree"], int(spec["modulus_hex"], 16))
    target = (int(data["target"]["x_hex"], 16), int(data["target"]["y_hex"], 16))
    basis = tuple(int(value, 16) for value in data["source_basis_hex"])
    offset = int(data.get("source_offset_hex", "0x0"), 16)
    model = prepare_five_model(f, int(curve["a_hex"], 16), int(curve["b_hex"], 16),
                               target, offset, basis)
    rows = []
    for index, h in enumerate(enumerate_family(model["outer_family"])):
        if index >= outer_samples:
            break
        for k in d0_samples:
            family = conditional_five_family(model, h, k)
            # Replay one admissible parameter from a consistent family;
            # this checks the norm and every support constraint, not roots.
            if family["consistent"]:
                r = family["particular"] & (f.order - 1)
                s = family["particular"] >> f.degree
                polynomial = five_coefficients_norm(model, h, k, r, s)
                if not all(support_conditions(model, polynomial)):
                    raise ArithmeticError("conditional solution failed support replay")
            rows.append({"d1_hex": hex(h), "d0_hex": hex(k),
                         "consistent": family["consistent"],
                         "rank": family["rank"], "nullity": family["nullity"],
                         "relaxed_cardinality_exact": family["cardinality"]})
    return {"schema": "five-source-divisor-support-profiler-v1",
            "instance_id": data["instance_id"], "field_degree": f.degree,
            "field_modulus_hex": hex(f.modulus), "source_dimension": len(basis),
            "target": [hex(value) for value in target],
            "negative_target": [hex(value) for value in model["negative_target"]],
            "support_dimensions": list(model["support_dimensions"]),
            "outer_rank": model["outer_family"]["rank"],
            "outer_nullity": model["outer_family"]["nullity"],
            "outer_cardinality_exact": model["outer_family"]["cardinality"],
            "conditional_variable_bits": 2 * f.degree,
            "conditional_constraint_rank_ceiling": model["conditional_constraint_ceiling"],
            "conditional_nullity_floor": model["conditional_nullity_floor"],
            "nonlinear_products_at_fixed_d1": ["d0*c2", "d0*c1"],
            "profiles": rows, "solver_complete": False,
            "claim_scope": "Exact necessary-condition ranks only; no root solving, timing, speed or novelty claim."}


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("instance")
    parser.add_argument("--outer-samples", type=int, default=2)
    arguments = parser.parse_args()
    print(json.dumps(profile_instance(arguments.instance, arguments.outer_samples), indent=2))

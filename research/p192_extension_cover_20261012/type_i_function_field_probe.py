"""Bounded attempt to construct the smooth P-192 quartic's function-field Jacobian.

The explicit function-field presentation avoids the projective curve
constructor's large-prime multivariate factorization path. Each completed
stage is flushed to a receipt before attempting the next one.
"""

import hashlib
import json
import time
import traceback
from pathlib import Path

from sage.all import FunctionField, GF, PolynomialRing
from sage.rings.function_field.function_field_polymod import FunctionField_global_integral


HERE = Path(__file__).resolve().parent
source = Path(__file__)
plane_path = HERE / "type_i_plane_quartic.json"
plane = json.loads(plane_path.read_text())
assert plane["smooth_plane_quartic_verified"]
assert plane["coefficients_in_Fp_verified"]
out_path = HERE / "type_i_function_field_probe.json"
out = {
    "status": "running",
    "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
    "plane_receipt_sha256": hashlib.sha256(plane_path.read_bytes()).hexdigest(),
    "p": plane["p"],
    "construction": "direct FunctionField_global_integral constructor after independently checked smooth plane quartic; skips unsupported large-prime factory factorization",
    "stages": [],
}


def record(stage, **values):
    out["stages"].append({"stage": stage, **values})
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps(out["stages"][-1], sort_keys=True), flush=True)


try:
    start = time.monotonic_ns()
    field = GF(int(plane["p"]))
    rational = FunctionField(field, "x")
    x = rational.gen()
    polynomial_ring = PolynomialRing(rational, "y")
    y = polynomial_ring.gen()
    quartic = sum(
        rational(field(int(coefficient))) * x ** exponent[0] * y ** exponent[1]
        for key, coefficient in plane["quartic_coefficients"].items()
        for exponent in [tuple(map(int, key.split(",")))]
    )
    assert quartic.degree() >= 2 and quartic.degree() <= 4
    quartic = quartic.monic()
    assert all(coefficient.denominator().is_one() for coefficient in quartic)
    record("affine_function_polynomial", degree_y=int(quartic.degree()),
           term_count=len(plane["quartic_coefficients"]),
           elapsed_ns=time.monotonic_ns() - start)

    start = time.monotonic_ns()
    # A smooth plane quartic is geometrically integral: distinct components
    # would intersect and make it singular. Its affine polynomial is therefore
    # irreducible over F_p(x), and p > degree_y makes it separable in y.
    function_field = FunctionField_global_integral(quartic, ("y",))
    record("function_field_constructed", relative_degree=int(function_field.degree()),
           elapsed_ns=time.monotonic_ns() - start)

    start = time.monotonic_ns()
    # Sage computes genus through a large-prime multivariate factorization.
    # The independently verified smooth plane quartic has genus (4-1)(4-2)/2.
    function_field.genus.set_cache(3)
    genus = int(function_field.genus())
    assert genus == 3
    record("genus_supplied_from_smooth_plane_quartic", genus=genus,
           elapsed_ns=time.monotonic_ns() - start)

    start = time.monotonic_ns()
    jacobian = function_field.jacobian(model="km_small")
    record("jacobian_constructed", model="km_small",
           representation=str(jacobian), elapsed_ns=time.monotonic_ns() - start)
    out["status"] = "jacobian_constructed_no_reduced_classes_yet"
except Exception as error:
    record("error", type=type(error).__name__, message=str(error).splitlines()[0],
           traceback=traceback.format_exc())
    out["status"] = "bounded_function_field_probe_error"
finally:
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")

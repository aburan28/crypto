#!/usr/bin/env python3
"""Frozen, explicitly opt-in degree-7 leaf-bridge certificate producer.

This file is syntax/hash checked in CI but never executed there. Run only
after an independent review of its exact commit and FROZEN.json.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import resource
import signal
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "research/ecc2k130_relations"))
sys.path.insert(0, str(REPO / "research/ecc2k130_oriented_transport_20260924"))
from relations import Koblitz, challenge_curve  # noqa: E402
from oriented_velu import BinaryVeluMap  # noqa: E402
from cubic_velu import CubicVeluMap  # noqa: E402


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frozen_inputs(spec: dict) -> None:
    for rel, expected in spec["input_sha256"].items():
        actual = digest(REPO / rel)
        if actual != expected:
            raise AssertionError(f"input hash mismatch: {rel}: {actual}")
    for name, expected in spec["implementation_sha256"].items():
        actual = digest(HERE / name)
        if actual != expected:
            raise AssertionError(f"implementation hash mismatch: {name}: {actual}")


def point_hex(point):
    return None if point is None else [hex(point[0]), hex(point[1])]


def line_apply(A, v, p=263):
    return [(A[0][0] * v[0] + A[0][1] * v[1]) % p,
            (A[1][0] * v[0] + A[1][1] * v[1]) % p]


def normalize(v, p=263):
    if v[0]:
        return [1, v[1] * pow(v[0], -1, p) % p]
    return [0, 1]


def matmul(A, B, p=263):
    return [[sum(A[i][k] * B[k][j] for k in range(2)) % p
             for j in range(2)] for i in range(2)]


def matrix_preflight(spec):
    saved = spec["saved_kernel"]
    col0, col1 = saved["twist_frobenius_matrix_columns_mod263"]
    M = [[col0[0], col1[0]], [col0[1], col1[1]]]
    A = [[(int(i == j) - 2 * M[i][j]) % 263 for j in range(2)]
         for i in range(2)]
    assert A == saved["source_norm7_matrix_rows_mod263"]
    assert matmul(A, A) == [[256, 0], [0, 256]]
    assert normalize(line_apply(A, saved["source_line"])) == saved["norm7_image_line"]
    bad = [[(int(i == j) + 2 * M[i][j]) % 263 for j in range(2)]
           for i in range(2)]
    assert normalize(line_apply(bad, saved["source_line"])) != saved["norm7_image_line"]
    v = saved["other_orbit_line"]
    hits = []
    for k in range(131):
        if normalize(v) == saved["norm7_image_line"]:
            hits.append(k)
        v = line_apply(M, v)
    assert hits == [saved["other_orbit_frobenius_exponent"]] == [23]
    return {"A_rows": A, "A_squared_rows": matmul(A, A),
            "conjugacy_exponents": hits, "wrong_sign_rejected": True}


def alpha(curve: Koblitz, point):
    return curve.add(point, curve.mul(curve.frobenius(point), 2))


def public_cases(E, P, Q, r, spec):
    cases = [
        ("P", P, None, None),
        ("Q", Q, None, None),
        ("P+Q", E.add(P, Q), None, None),
        ("2P", E.mul(P, 2), None, None),
    ]
    labels = spec["public_target_labels"]
    for i in range(labels["count"]):
        prefix = labels["prefix"]
        u = int.from_bytes(hashlib.sha256(f"{prefix}|{i}|u".encode()).digest(), "big") % r
        v = 1 + int.from_bytes(hashlib.sha256(f"{prefix}|{i}|v".encode()).digest(), "big") % (r - 1)
        point = E.add(E.mul(P, u), E.mul(Q, v))
        assert point is not None
        cases.append((f"target-{i}", point, u, v))
    return cases


def checked_x_isomorphism(source: Koblitz, target: Koblitz):
    """Normalized characteristic-two model change (x,y)->(x,y+s*x)."""
    assert (source.F.deg, source.F.irr) == (target.F.deg, target.F.irr)
    assert source.b == target.b
    delta = source.a ^ target.a
    if delta == 0:
        return [0, 1]
    s = source.F.solve_artin_schreier(delta)
    assert s is not None and source.F.sqr(s) ^ s == delta
    return [s, s ^ 1]


def apply_x_isomorphism(F, point, s):
    return None if point is None else (point[0], point[1] ^ F.mul(s, point[0]))


def build_sage_field(spec):
    from sage.all import GF, PolynomialRing  # noqa: PLC0415

    b = GF(2)
    R = PolynomialRing(b, "z0")
    bits = int(spec["field"]["modulus_hex"], 16)
    modulus = R([b((bits >> i) & 1) for i in range(bits.bit_length())])
    assert bits.bit_length() - 1 == spec["field"]["degree"]
    F = GF(2 ** spec["field"]["degree"], name="z", modulus=modulus)
    assert F.modulus() == modulus
    return F


def sage_element(F, value: int):
    z = F.gen()
    out = F.zero()
    for i in range(F.degree() - 1, -1, -1):
        out = out * z + F((value >> i) & 1)
    return out


def bit_integer(value):
    return sum(int(coef) << i for i, coef in enumerate(value.polynomial().list()))


def sage_curve(F, source):
    from sage.all import EllipticCurve  # noqa: PLC0415
    return EllipticCurve(F, [F.one(), sage_element(F, source.a),
                             F.zero(), F.zero(), sage_element(F, source.b)])


def sage_point(curve, point):
    if point is None:
        return curve(0)
    F = curve.base_ring()
    return curve([sage_element(F, point[0]), sage_element(F, point[1])])


def native_point(point):
    if point.is_zero():
        return None
    return (bit_integer(point[0]), bit_integer(point[1]))


def choose_sage_iso(psi, target_curve, witness, expected):
    choices = psi.codomain().isomorphisms(target_curve)
    match = [iso for iso in choices if iso(psi(witness)) == expected]
    assert len(match) == 1, f"expected one oriented codomain isomorphism; found {len(match)}"
    return match[0], len(choices)


def encode_cubic(value):
    """Canonical 3 Fq coefficients of a class in Fq[U]/h(U)."""
    coeffs = value.lift().list()
    assert len(coeffs) <= 3
    return [hex(bit_integer(coeffs[i])) if i < len(coeffs) else "0x0"
            for i in range(3)]


def encode_sextic(value):
    """Canonical 2 by 3 Fq coefficients for Fq[U,V]/(h,V²+V+1)."""
    coeffs = value.lift().list()
    assert len(coeffs) <= 2
    zero = value.parent().base_ring().zero()
    return [encode_cubic(coeffs[i] if i < len(coeffs) else zero)
            for i in range(2)]


def exact_kernel_points(curve, isogeny, kernel_poly, q, expected_count):
    roots = kernel_poly.roots()
    assert len(roots) == 3 and all(int(multiplicity) == 1 for _, multiplicity in roots)
    root_records = []
    point_records = []
    seen = set()
    eigen_ok = 0
    for x, _ in roots:
        root_records.append(encode_sextic(x))
        lifts = curve.lift_x(x, all=True)
        assert len(lifts) == 2
        for point in lifts:
            assert not point.is_zero() and (7 * point).is_zero()
            image5 = 5 * point
            assert (point[0] ** q, point[1] ** q) == (image5[0], image5[1])
            assert isogeny(point).is_zero()
            record = {"x": encode_sextic(point[0]), "y": encode_sextic(point[1])}
            key = json.dumps(record, sort_keys=True)
            assert key not in seen
            seen.add(key)
            point_records.append(record)
            eigen_ok += 1
    assert len(seen) == expected_count == 6
    assert int(kernel_poly.degree()) == 3
    return {"root_count": len(roots), "point_count": len(seen),
            "roots": sorted(root_records, key=lambda row: json.dumps(row)),
            "points": sorted(point_records, key=lambda row: json.dumps(row, sort_keys=True)),
            "kernel_polynomial_low_to_high":
                [encode_sextic(kernel_poly[i]) for i in range(4)],
            "exact_order_seven": True, "q_frobenius_eigenvalue": 5,
            "eigenvalue_checks": eigen_ok, "all_mapped_to_infinity": True}


def clock():
    usage = resource.getrusage(resource.RUSAGE_SELF)
    # macOS reports bytes; Linux reports KiB.
    peak = int(usage.ru_maxrss)
    if sys.platform.startswith("linux"):
        peak *= 1024
    return {"wall_seconds": time.monotonic(), "cpu_seconds": time.process_time(),
            "peak_rss_bytes": peak}


def timeout_handler(_signum, _frame):
    raise TimeoutError("frozen 300-second wall cap reached")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    if args.out.exists():
        raise SystemExit("refusing to overwrite an existing receipt")
    spec = json.loads((HERE / "FROZEN.json").read_text())
    started = clock()
    phase = "initialization"
    receipt = {"schema": "ecc2k130-leaf7-bridge-certificate-v1",
               "status": "STOP", "structural_gate_status": "UNRESOLVED",
               "freeze_sha256": digest(HERE / "FROZEN.json"),
               "source_main_commit": spec["source_main_commit"],
               "parent_note_commit": spec["parent_note_commit"],
               "release_main_head": spec["release_main_head"],
               "release_gate": spec["release_gate"],
               "input_sha256": spec["input_sha256"],
               "implementation_sha256": spec["implementation_sha256"],
               "phase": phase}
    try:
        assert spec["status"] == "protocol_only_no_outcome"
        if spec["release_main_head"] is None:
            raise RuntimeError("parent note unmerged: post-merge re-freeze and review required before outcome")
        signal.signal(signal.SIGALRM, timeout_handler)
        signal.setitimer(signal.ITIMER_REAL, spec["caps"]["child_wall_seconds"])
        try:
            resource.setrlimit(resource.RLIMIT_AS,
                               (spec["caps"]["child_peak_rss_bytes"],
                                spec["caps"]["child_peak_rss_bytes"]))
        except (AttributeError, OSError, ValueError) as exc:
            receipt["memory_limit_setup"] = f"STOP: unsupported {type(exc).__name__}"
            raise RuntimeError("frozen 2-GiB memory cap could not be enforced") from exc
        receipt["memory_limit_setup"] = "enforced: RLIMIT_AS"
        phase = "frozen_hashes_and_matrix"
        frozen_inputs(spec)
        receipt["matrix"] = matrix_preflight(spec)

        phase = "source_and_saved_263_maps"
        Fpy, E, P, Q, r = challenge_curve()
        assert r == int(spec["field"]["subgroup_order"]) and E.mul(P, r) is None
        assert E.mul(Q, r) is None and P != Q
        assert alpha(E, alpha(E, P)) == E.neg(E.mul(P, 7))
        assert alpha(E, alpha(E, Q)) == E.neg(E.mul(Q, 7))
        saved = json.loads((REPO / "research/ecc2k130_direction_review_20260924/"
                            "twist_torsion_results.json").read_text())
        row = next(row for row in saved["runs"]
                   if row["seed"] == spec["saved_kernel"]["torsion_seed"])
        G, H = [tuple(int(coordinate, 16) for coordinate in point)
                for point in row["basis_on_twist"]]
        twist = Koblitz(Fpy, a=1, b=1)
        assert twist.on_curve(G) and twist.on_curve(H)
        Gprime = twist.add(G, twist.neg(twist.mul(twist.frobenius(G), 2)))
        assert (twist.mul(Gprime, pow(155, -1, 263))
                == twist.add(G, twist.mul(H, 74)))
        phi0 = BinaryVeluMap.from_generator(E, twist, G, 263)
        phi1 = BinaryVeluMap.from_generator(E, twist, Gprime, 263)
        other_generator = twist.add(G, twist.mul(H, 4))
        phi_other = BinaryVeluMap.from_generator(E, twist, other_generator, 263)
        coefficient_hits = [k for k in range(131)
                            if Fpy.frobenius(phi_other.codomain.b, k) == phi1.codomain.b]
        assert coefficient_hits == [spec["saved_kernel"]["other_orbit_frobenius_exponent"]]
        assert phi0.codomain.b != phi1.codomain.b
        receipt["other_orbit_b"] = hex(phi_other.codomain.b)
        receipt["model_conjugacy_exponents"] = coefficient_hits
        receipt["other_orbit_map_role"] = "structural control only; excluded from direct second-leaf cost arm"
        cases = public_cases(E, P, Q, r, spec)
        case_data = []
        for label, T, u, v in cases:
            AT = alpha(E, T)
            assert E.mul(AT, r) is None
            case_data.append({"label": label, "u": None if u is None else str(u),
                              "v": None if v is None else str(v),
                              "source": point_hex(T), "alpha_source": point_hex(AT),
                              "leaf0": point_hex(phi0(T)),
                              "direct_leaf1": point_hex(phi1(AT))})
        receipt["cases"] = case_data
        receipt["leaf0_b"] = hex(phi0.codomain.b)
        receipt["leaf1_b"] = hex(phi1.codomain.b)
        receipt["saved_generator"] = point_hex(G)
        receipt["derived_generator"] = point_hex(Gprime)

        phase = "sage_7_division_and_cubic"
        import sage.all as sage  # noqa: PLC0415
        from sage.version import version as sage_version  # noqa: PLC0415
        assert str(sage_version).startswith(spec["sage_version"])
        Fs = build_sage_field(spec)
        E0s = sage_curve(Fs, phi0.codomain)
        div7 = E0s.division_polynomial(7)
        factors = [(factor, int(power)) for factor, power in div7.factor()]
        degrees = sorted((int(factor.degree()), power) for factor, power in factors)
        assert [[degree, power] for degree, power in degrees] == spec["gate"]["predicted_7_division_factor_degrees"]
        h = next(factor.monic() for factor, power in factors
                 if int(factor.degree()) == 3 and power == 1)
        assert h.is_irreducible()
        cubic_coeffs = tuple(bit_integer(h[i]) for i in range(3))
        cubic = CubicVeluMap(phi0.codomain, cubic_coeffs)
        assert cubic.codomain.b == phi1.codomain.b
        receipt["division_factor_degrees"] = [[d, power] for d, power in degrees]
        receipt["kernel_polynomial_low_to_high"] = [hex(x) for x in cubic_coeffs] + ["0x1"]
        receipt["cubic_codomain_b"] = hex(cubic.codomain.b)
        psi = E0s.isogeny(h)
        assert int(psi.degree()) == 7 and psi.is_separable()
        dual = psi.dual()
        Ecs = sage_curve(Fs, cubic.codomain)
        E1s = sage_curve(Fs, phi1.codomain)
        witness = sage_point(E0s, tuple(int(x, 16) for x in case_data[0]["leaf0"]))
        expected_cubic = sage_point(Ecs, cubic(tuple(int(x, 16) for x in case_data[0]["leaf0"])))
        iso_c, count_c = choose_sage_iso(psi, Ecs, witness, expected_cubic)
        target_candidates = checked_x_isomorphism(cubic.codomain, phi1.codomain)
        oriented = [s for s in target_candidates
                    if apply_x_isomorphism(Fpy,
                                           cubic(tuple(int(x, 16) for x in case_data[0]["leaf0"])),
                                           s)
                    == tuple(int(x, 16) for x in case_data[0]["direct_leaf1"])]
        assert len(oriented) == 1
        s = oriented[0]
        iso_t, count_t = choose_sage_iso(
            psi, E1s, witness,
            sage_point(E1s, tuple(int(x, 16) for x in case_data[0]["direct_leaf1"])))
        assert count_c >= 1 and count_t >= 1
        receipt["codomain_model_isomorphism_s"] = hex(s)
        receipt["sage_codomain_ainvs"] = [hex(bit_integer(v)) for v in psi.codomain().a_invariants()]
        receipt["sage_to_cubic_isomorphism_count"] = count_c
        receipt["sage_to_target_isomorphism_count"] = count_t

        phase = "full_point_commuting_square"
        for case in case_data:
            T0 = tuple(int(x, 16) for x in case["leaf0"])
            direct = tuple(int(x, 16) for x in case["direct_leaf1"])
            cubic_image = cubic(T0)
            assert apply_x_isomorphism(Fpy, cubic_image, s) == direct
            sage_image = psi(sage_point(E0s, T0))
            assert iso_c(sage_image) == sage_point(Ecs, cubic_image)
            assert iso_t(sage_image) == sage_point(E1s, direct)
            assert dual(sage_image) == 7 * sage_point(E0s, T0)
            case["cubic_leaf1"] = point_hex(cubic_image)
            case["sage_oriented_leaf1"] = point_hex(native_point(iso_t(sage_image)))
        assert cubic(None) is None and psi(E0s(0)).is_zero()
        assert len(case_data) == spec["gate"]["predicted_dual_composition_controls"]
        receipt["dual_composition_controls"] = len(case_data)
        try:
            cubic((0, 0))
        except ValueError:
            receipt["off_curve_rejected"] = True
        else:
            raise AssertionError("off-curve input was accepted")

        phase = "forward_and_dual_exceptional_points"
        # The cubic abscissae live over q^3; their point lifts live over q^6.
        F3 = Fs.extension(h, "u")
        # [F3:F2]=393 is odd, so Tr_F3/F2(1)=1 and V²+V+1 is irreducible.
        R3 = sage.PolynomialRing(F3, "w")
        w = R3.gen()
        F6 = F3.extension(w ** 2 + w + 1, "v")
        E0e = E0s.base_extend(F6)
        he = h.change_ring(F6)
        psie = E0e.isogeny(he)
        duale = psie.dual()
        q = 2 ** spec["field"]["degree"]
        receipt["forward_kernel"] = exact_kernel_points(
            E0e, psie, he, q, spec["gate"]["predicted_forward_kernel_points"])
        receipt["reverse_kernel"] = exact_kernel_points(
            psie.codomain(), duale, duale.kernel_polynomial(), q,
            spec["gate"]["predicted_reverse_kernel_points"])
        assert psie(E0e(0)).is_zero()
        assert duale(psie.codomain()(0)).is_zero()
        receipt["infinity_mapped_to_infinity"] = True

        phase = "complete"
        receipt["status"] = "PRODUCER_PASS"
        receipt["structural_gate_status"] = "UNRESOLVED_PENDING_REPLAY"
        receipt["cost"] = {
            "shared_263_torsion_reused_for_direct_second_leaf": True,
            "direct_second_leaf_rederived_from_saved_basis": True,
            "incremental_mul_equivalent_ratio": None,
            "cold_mul_equivalent_ratio": None,
            "ratio_unset_reason": "Sage kernel discovery and extension operations are not metered in a comparable Fq unit",
            "promotion": False}
    except BaseException as exc:
        receipt["status"] = "STOP"
        receipt["error_type"] = type(exc).__name__
        receipt["error"] = str(exc)
        receipt["traceback"] = traceback.format_exc(limit=12)
    finally:
        end = clock()
        receipt["phase"] = phase
        receipt["resources"] = {
            "wall_seconds": end["wall_seconds"] - started["wall_seconds"],
            "cpu_seconds": end["cpu_seconds"] - started["cpu_seconds"],
            "peak_rss_bytes": end["peak_rss_bytes"]}
        if receipt["resources"]["wall_seconds"] > spec["caps"]["child_wall_seconds"]:
            receipt["status"] = "STOP"
            receipt.setdefault("error", "wall time exceeded frozen cap")
        if end["peak_rss_bytes"] > spec["caps"]["child_peak_rss_bytes"]:
            receipt["status"] = "STOP"
            receipt.setdefault("error", "peak RSS exceeded frozen cap")
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
        signal.setitimer(signal.ITIMER_REAL, 0)
    return 0 if receipt["status"] == "PRODUCER_PASS" else 1


if __name__ == "__main__":
    sys.exit(main())

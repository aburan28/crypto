#!/usr/bin/env python3
"""Independent Sage check of public fixtures at the eleven exact NIST sizes.

Run with Sage's Python. This validates constants and known point sums only;
it performs no factor-base search, relation collection or logarithm recovery.
"""

import argparse
import hashlib
import json
import re
from pathlib import Path

from sage.all import EllipticCurve, GF, Integer, PolynomialRing
from sage.version import version as sage_version

ROOT = Path(__file__).resolve().parents[1]
NAMES = (
    "p256", "sect163k1", "sect163r2", "sect233k1", "sect233r1",
    "sect283k1", "sect283r1", "sect409k1", "sect409r1",
    "sect571k1", "sect571r1",
)


def body(source, name):
    match = re.search(r"pub fn " + re.escape(name) + r"\(\) -> Self \{(.*?)\n    \}", source, re.S)
    if not match:
        raise ValueError(f"constructor {name} missing")
    return match.group(1)


def parameters(name, prime, binary, fields):
    if name == "p256":
        constructor = body(prime, name)
        values = {}
        for key in ("p", "a", "b", "gx", "gy", "n"):
            match = re.search(key + r': BigUint::parse_bytes\(\s*b"([0-9a-fA-F]+)"', constructor)
            if not match:
                raise ValueError(f"{name}: cannot read {key}")
            values[key] = Integer(match.group(1), 16)
        h = Integer(re.search(r"h:\s*(\d+)", constructor).group(1))
        field = GF(values["p"])
        record = {
            "model": "short-weierstrass", "encoding": "hex-integer",
            "field": {"characteristic": format(int(values["p"]), "x"), "degree": 1},
        }
        return field, values["a"], values["b"], values["gx"], values["gy"], values["n"], h, record
    constructor = body(binary, name)
    m = int(re.search(r"sect(\d+)", name).group(1))
    modulus_constructor = body(fields, f"deg_{m}")
    low = [int(v) for v in re.search(r"low_terms: vec!\[(.*?)\]", modulus_constructor).group(1).split(",")]
    terms = sorted(low + [m])
    ring = PolynomialRing(GF(2), "z")
    z = ring.gen()
    modulus = sum(z**i for i in terms)
    if not modulus.is_irreducible():
        raise ValueError(f"{name}: reducible field modulus")
    field = GF(2**m, name="z", modulus=modulus)
    if "Self::from_hex_parts" in constructor:
        raw = re.findall(r'"([0-9a-fA-F]+)"', constructor)
        if len(raw) != 5:
            raise ValueError(f"{name}: expected five parameter strings")
        a, b, gx, gy, n = (Integer(s, 16) for s in raw)
        h = Integer(re.search(r'"[0-9a-fA-F]+",\s*(\d+),\s*\)', constructor).group(1))
    else:
        def element(key):
            if re.search(r"let " + key + r" = F2mElement::one\(m\)", constructor):
                return Integer(1)
            return Integer(re.search(r"let " + key + r' = F2mElement::from_hex\("([0-9a-fA-F]+)"', constructor).group(1), 16)
        a, b, gx, gy = (element(key) for key in ("a", "b", "gx", "gy"))
        n = Integer(re.search(r'BigUint::parse_bytes\(b"([0-9a-fA-F]+)"', constructor).group(1), 16)
        h = Integer(re.search(r"cofactor: BigUint::from\((\d+)u32\)", constructor).group(1))
    record = {
        "model": "binary-weierstrass", "encoding": "hex-polynomial-bits",
        "field": {"characteristic": "2", "degree": m, "basis": "polynomial",
                  "modulus_exponents": terms},
    }
    return field, a, b, gx, gy, n, h, record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-revision", help="Reviewed source revision, when known")
    args = parser.parse_args()
    paths = ("src/ecc/curve.rs", "src/binary_ecc/curve.rs", "src/binary_ecc/f2m.rs")
    sources = [(ROOT / p).read_text() for p in paths]
    rows = []
    for name in NAMES:
        field, a, b, gx, gy, n, h, record = parameters(name, *sources)
        lift = field if name == "p256" else field.from_integer
        coefficients = [0, 0, 0, lift(a), lift(b)] if name == "p256" else [1, lift(a), 0, 0, lift(b)]
        curve = EllipticCurve(field, coefficients)
        generator = curve(lift(gx), lift(gy))
        q = field.cardinality()
        checks = {
            "generator_on_curve": True,
            "generator_subgroup": bool(n * generator == curve(0)),
            "subgroup_order_prime": bool(n.is_prime(proof=True)),
            "positive_cofactor": bool(h > 0),
            "hasse_bound": bool((h*n - (q+1))**2 <= 4*q),
        }
        record.update(a=format(int(a), "x"), b=format(int(b), "x"),
                      generator={"x": format(int(gx), "x"), "y": format(int(gy), "x")},
                      subgroup_order=format(int(n), "x"), cofactor=format(int(h), "x"))
        encode = (lambda el: format(int(el), "x")) if name == "p256" else (lambda el: format(int(el.to_integer()), "x"))
        fixtures = []
        for u, v in ((Integer(1), Integer(2)), (n-2, n-3)):
            w = (u + v) % n
            p, point_q, target = u * generator, v * generator, w * generator
            wrong = ((w+1) % n) * generator
            x1, x2, x3 = p[0], point_q[0], target[0]
            def s3(x):
                if name == "p256":
                    return (x1-x2)**2*x**2 - 2*((x1+x2)*(x1*x2+lift(a))+2*lift(b))*x + (x1*x2-lift(a))**2 - 4*lift(b)*(x1+x2)
                return (x1+x2)**2*x**2 + x1*x2*x + (x1*x2)**2 + lift(b)
            fixture = {
                "u": format(int(u), "x"), "v": format(int(v), "x"), "w": format(int(w), "x"),
                "target": {"x": encode(target[0]), "y": encode(target[1])},
                "point_sum_verified": bool(p + point_q == target), "s3_zero": bool(s3(x3) == 0),
                "wrong_target_rejected": bool(p + point_q != wrong),
                "wrong_s3_rejected": bool(s3(wrong[0]) != 0),
            }
            fixtures.append(fixture)
        verified = all(checks.values()) and all(all(f[k] for k in (
            "point_sum_verified", "s3_zero", "wrong_target_rejected", "wrong_s3_rejected"
        )) for f in fixtures)
        rows.append({"curve": name, "exact_parameters": record, "checks": checks,
                     "fixtures": fixtures, "verified": verified})
        print(f"{name}: {'PASS' if verified else 'FAIL'}", flush=True)
    report = {
        "schema_version": 1, "scope": "independent public-parameter/S3 diagnostic; no DLP recovery",
        "source_revision": args.source_revision,
        "sage_version": sage_version,
        "source_sha256": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        "curves": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if not all(row["verified"] for row in rows):
        raise SystemExit("parameter diagnostic failed; inspect the saved receipt")


if __name__ == "__main__":
    main()

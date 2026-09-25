#!/usr/bin/env python3
"""Frozen width refusal and exhaustive n13 K0 <-> split-mu4 point-map audit."""

from __future__ import annotations

import argparse
import hashlib
import json
import resource
import subprocess
import time
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
O = None


def canonical(obj: object) -> bytes:
    return (json.dumps(obj, sort_keys=True, separators=(",", ":")) + "\n").encode()


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def preflight() -> tuple[dict, dict]:
    frozen = json.loads((HERE / "FROZEN.json").read_text())
    for name, digest in frozen["files"].items():
        path = ROOT / name
        assert path.is_file(), f"frozen file missing: {name}"
        assert sha(path.read_bytes()) == digest, f"frozen hash drift: {name}"
    inp = json.loads((HERE / "INPUT.json").read_text())
    assert inp["schema"] == "symbolic-oaware-width-gate-v1"
    assert inp["toy"] == {"n": 13, "modulus": 8219, "curve_a": 0, "curve_b": 1}
    assert [a["id"] for a in inp["arms"]] == [
        "unequal-m9", "unequal-m10", "balanced-m10-control"
    ]
    assert len(inp["required_for_n131_admission"]) == 5
    return frozen, inp


def mul(a: int, b: int, n: int, mod: int) -> int:
    mask = (1 << n) - 1
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & (1 << n):
            a ^= mod
    return result & mask


def sq(a: int, n: int, mod: int) -> int:
    return mul(a, a, n, mod)


def power(a: int, exponent: int, n: int, mod: int) -> int:
    z = 1
    while exponent:
        if exponent & 1:
            z = mul(z, a, n, mod)
        a = sq(a, n, mod)
        exponent >>= 1
    return z


def inv(a: int, n: int, mod: int) -> int:
    assert a != 0
    out = power(a, (1 << n) - 2, n, mod)
    assert mul(a, out, n, mod) == 1
    return out


def trace(a: int, n: int, mod: int) -> int:
    t = a
    z = 0
    for _ in range(n):
        z ^= t
        t = sq(t, n, mod)
    assert z in (0, 1) and t == a
    return z


def halftrace(a: int, n: int, mod: int) -> int:
    assert n % 2 == 1
    z = 0
    t = a
    for _ in range((n + 1) // 2):
        z ^= t
        t = sq(sq(t, n, mod), n, mod)
    assert sq(z, n, mod) ^ z == a
    return z


def on_k0(point: tuple[int, int] | None, n: int, mod: int) -> bool:
    if point is O:
        return True
    x, y = point
    return (sq(y, n, mod) ^ mul(x, y, n, mod)) == (mul(sq(x, n, mod), x, n, mod) ^ 1)


def add(p: tuple[int, int] | None, q: tuple[int, int] | None, n: int, mod: int):
    assert on_k0(p, n, mod) and on_k0(q, n, mod)
    if p is O:
        return q
    if q is O:
        return p
    x1, y1 = p
    x2, y2 = q
    if x1 == x2:
        if y1 ^ y2 == x1:
            return O
        assert y1 == y2 and x1 != 0
        lam = x1 ^ mul(y1, inv(x1, n, mod), n, mod)
        xr = sq(lam, n, mod) ^ lam
        yr = sq(x1, n, mod) ^ mul(lam ^ 1, xr, n, mod)
    else:
        lam = mul(y1 ^ y2, inv(x1 ^ x2, n, mod), n, mod)
        xr = sq(lam, n, mod) ^ lam ^ x1 ^ x2
        yr = mul(lam, x1 ^ xr, n, mod) ^ xr ^ y1
    out = (xr, yr)
    assert on_k0(out, n, mod)
    return out


def image(point: tuple[int, int] | None, n: int, mod: int) -> tuple[int, int, int, int]:
    if point is O:
        return (1, 1, 0, 1)
    x, y = point
    x2 = sq(x, n, mod)
    return (x2, x2 ^ y, 1, x2 ^ x ^ y)


def on_mu4(v: tuple[int, int, int, int], n: int, mod: int) -> bool:
    a, b, c, d = v
    if not any(v):
        return False
    return sq(a ^ c, n, mod) == mul(b, d, n, mod) and sq(b ^ d, n, mod) == mul(a, c, n, mod)


def inverse_image(v: tuple[int, int, int, int], n: int, mod: int):
    assert on_mu4(v, n, mod)
    a, b, c, d = v
    if c == 0:
        assert a == b == d != 0
        return O
    c_inv = inv(c, n, mod)
    out = (mul(b ^ d, c_inv, n, mod), mul(a ^ b, c_inv, n, mod))
    assert on_k0(out, n, mod)
    return out


def source_points(n: int, mod: int) -> set[tuple[int, int]]:
    points = {(0, 1)}
    for x in range(1, 1 << n):
        xi = inv(x, n, mod)
        rhs = x ^ sq(xi, n, mod)
        if trace(rhs, n, mod):
            continue
        z = halftrace(rhs, n, mod)
        y = mul(x, z, n, mod)
        points.update(((x, y), (x, y ^ x)))
    return points


def mu4_chart_points(n: int, mod: int) -> set[tuple[int, int, int, int]]:
    """Enumerate the X2=1 chart from its quadrics, not from source y roots."""
    points = set()
    for a in range(1 << n):
        x = power(a, 1 << (n - 1), n, mod)  # unique square root
        assert sq(x, n, mod) == a
        if x == 0:
            points.add((0, 1, 1, 1))
            continue
        t = mul(sq(a ^ 1, n, mod), inv(a, n, mod), n, mod)
        if trace(t, n, mod):
            continue
        z = halftrace(t, n, mod)
        b = mul(x, z, n, mod)
        points.update(((a, b, 1, b ^ x), (a, b ^ x, 1, b)))
    assert all(on_mu4(v, n, mod) for v in points)
    return points


def width_decision(inp: dict) -> dict:
    source = (ROOT / "src/cryptanalysis/pq_descent_symbolic.rs").read_text()
    groebner = (ROOT / "src/cryptanalysis/pq_groebner_f2.rs").read_text()
    sat = (ROOT / "src/cryptanalysis/semaev_sat.rs").read_text()
    gf = (ROOT / "src/cryptanalysis/semaev_decomp.rs").read_text()
    assert "pub const MAX_VARS: u32 = 64;" in source
    assert "HashMap<u64, u64>" in source
    assert "pub mask: u64" in groebner
    assert 'assert!(k < 64, "monomial cap is 64 variables")' in groebner
    assert 'assert!(n_vars <= 64, "problem variables are packed into a u64")' in sat
    assert "pub fn model_assignment(&self) -> u64" in sat
    assert 'assert!(irr.degree <= 63, "Gf2 handles n ≤ 63")' in gf
    arms = []
    for arm in inp["arms"]:
        n, slots = arm["n"], arm["slots"]
        raw = sum(slots)
        chain = max(0, len(slots) - 2) * n
        arms.append({"id": arm["id"], "n": n, "m": len(slots),
                     "factor_bits": raw, "intermediate_x_bits": chain,
                     "affine_chain_floor_bits": raw + chain,
                     "current_u64_symbolic_admitted": False})
        assert raw > 64 and n > 63
    assert [(a["factor_bits"], a["affine_chain_floor_bits"]) for a in arms] == [
        (131, 1048), (131, 1179), (130, 1178)
    ]
    return {"status": "NOT_ADMITTED", "arms": arms,
            "blocking_interfaces": ["Gf2_degree_le_63", "FieldBoolPoly_u64_mask",
                                    "F2BoolMono_u64_mask", "SAT_wrapper_u64_model",
                                    "no_full_point_symbolic_circuit_or_n131_receipt"],
            "missing_required_artifacts": inp["required_for_n131_admission"]}


def toy_map(inp: dict) -> dict:
    n, mod = inp["toy"]["n"], inp["toy"]["modulus"]
    src = source_points(n, mod)
    images = set()
    for p in sorted(src):
        assert on_k0(p, n, mod)
        v = image(p, n, mod)
        assert on_mu4(v, n, mod)
        assert inverse_image(v, n, mod) == p
        scale = (p[0] ^ p[1]) or 1
        scaled = tuple(mul(c, scale, n, mod) for c in v)
        assert on_mu4(scaled, n, mod)
        assert inverse_image(scaled, n, mod) == p
        assert v not in images
        images.add(v)
    chart = mu4_chart_points(n, mod)
    assert images == chart
    assert inverse_image(image(O, n, mod), n, mod) is O
    assert not on_mu4((1, 0, 0, 1), n, mod)
    damaged = list(image((0, 1), n, mod))
    damaged[1] ^= 1
    assert not on_mu4(tuple(damaged), n, mod)
    try:
        c = image(O, n, mod)[2]
        assert c != 0, "O has no affine inverse chart"
    except AssertionError as e:
        assert str(e) == "O has no affine inverse chart"
    else:
        raise AssertionError("O chart guard failed")
    t = (1, 0)
    chain = []
    p = t
    for _ in range(4):
        p = add(p, t, n, mod)
        chain.append(p)
    assert chain == [(0, 1), (1, 1), O, t]
    return {"status": "PASS", "source_affine_points": len(src),
            "mu4_affine_chart_points": len(chart), "projective_points": len(src) + 1,
            "source_points_sha256": sha(canonical(sorted(src))),
            "mu4_chart_sha256": sha(canonical(sorted(chart))),
            "controls": {name: "PASS" for name in inp["toy_controls"]},
            "order4_chain": [[0, 1], [1, 1], "O", [1, 0]]}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    assert args.preflight != args.run
    frozen, inp = preflight()
    if args.preflight:
        print("FROZEN_HASHES_PASS")
        return
    assert args.output is not None and not args.output.exists()
    start = time.monotonic()
    cpu_start = time.process_time()
    gate = width_decision(inp)
    toy = toy_map(inp)
    wall = time.monotonic() - start
    cpu = time.process_time() - cpu_start
    usage = resource.getrusage(resource.RUSAGE_SELF)
    rss = usage.ru_maxrss if subprocess.run(["uname", "-s"], capture_output=True,
                                             text=True, check=True).stdout.strip() == "Darwin" else usage.ru_maxrss * 1024
    receipt = {"schema": "symbolic-oaware-width-gate-receipt-v1",
               "frozen_sha256": sha((HERE / "FROZEN.json").read_bytes()),
               "source_commit": frozen["parent_commit"],
               "width": gate, "toy_mu4_map": toy,
               "cost": {"wall_seconds": wall, "cpu_seconds": cpu, "peak_rss_bytes": rss}}
    assert wall <= inp["caps"]["wall_seconds"]
    assert rss <= inp["caps"]["rss_bytes"]
    blob = canonical(receipt)
    assert len(blob) <= inp["caps"]["receipt_bytes"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(blob)
    print(json.dumps({"status": "PASS", "receipt_sha256": sha(blob),
                      "receipt_bytes": len(blob), "width": gate["status"],
                      "toy": toy["status"]}, sort_keys=True))


if __name__ == "__main__":
    main()

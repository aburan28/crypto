#!/usr/bin/env python3
"""Independent saved-evidence check for the degree-7 paired base pilot."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "research/ecc2k130_relations"))
sys.path.insert(0, str(ROOT / "research/ecc2k130_oriented_transport_20260924"))
from fastfield import FastGF2m  # noqa: E402
from relations import Koblitz  # noqa: E402
from oriented_velu import BinaryVeluMap  # noqa: E402


def point(value):
    return None if value is None else tuple(value)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def independent_rank(rows, modulus):
    """RREF from scratch, with no reuse of the producer's incremental tracker."""
    a = [list(row) for row in rows]
    width = len(a[0]) if a else 17
    h = 0
    for col in range(width):
        pivot = next((i for i in range(h, len(a)) if a[i][col] % modulus), None)
        if pivot is None:
            continue
        a[h], a[pivot] = a[pivot], a[h]
        inv = pow(a[h][col], -1, modulus)
        a[h] = [(v * inv) % modulus for v in a[h]]
        for i in range(len(a)):
            if i != h and a[i][col] % modulus:
                f = a[i][col]
                a[i] = [(x - f * y) % modulus for x, y in zip(a[i], a[h])]
        h += 1
        if h == width:
            break
    return h


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("result", type=Path)
    args = parser.parse_args()
    d = json.loads(args.result.read_text())
    assert d["status"] == "PASS"
    src = d["source_sha256"]
    names = {
        "protocol": Path(__file__).with_name("PROTOCOL.md"),
        "runner": Path(__file__).with_name("run.py"),
        "fastfield": ROOT / "research/ecc2k130_relations/fastfield.py",
        "relations": ROOT / "research/ecc2k130_relations/relations.py",
        "oriented_velu": ROOT / "research/ecc2k130_oriented_transport_20260924/oriented_velu.py",
    }
    assert all(sha(path) == src[key] for key, path in names.items())
    p = d["parameters"]
    assert (p["degree"], p["irreducible"], p["isogeny_degree"]) == (21, "0x200005", 7)
    assert p["source_order"] == 2099948 and p["twist_order"] == 2094358
    assert p["subgroup_order"] == 421 and p["cofactor"] == 4988
    F = FastGF2m(21, 0x200005)
    E0, Tw = Koblitz(F, 0, 1), Koblitz(F, 1, 1)
    geom = d["geometry"]
    assert geom["line_count"] == 8 and geom["self_j_lines"] == 1
    assert tuple(geom["selected_kernel_x"]) == min(
        tuple(v["kernel_x"]) for v in geom["all_lines"]
        if v["codomain_b"] != 1)
    H = point(geom["selected_generator"])
    assert Tw.on_curve(H) and Tw.mul(H, 7) is None and H is not None
    half = [Tw.mul(H, i) for i in range(1, 4)]
    assert sorted(P[0] for P in half) == geom["selected_kernel_x"]
    phi = BinaryVeluMap.from_generator(E0, Tw, H, 7)
    E1 = phi.codomain
    assert E1.b == geom["selected_codomain_b"]
    G, Q = point(d["challenge"]["generator"]), point(d["challenge"]["target"])
    k = d["challenge"]["secret_audit_only"]
    assert E0.on_curve(G) and E0.mul(G, 421) is None and G is not None
    assert Q == E0.mul(G, k)
    G1, Q1 = phi(G), phi(Q)
    assert Q1 == E1.mul(G1, k)
    bases = {name: [point(P) for P in d["bases"][name]]
             for name in ["original", "transported", "descendant_native", "pullback"]}
    assert all(len(b) == len(set(b)) == 16 and None not in b
               for b in bases.values())
    assert bases["transported"] == [phi(P) for P in bases["original"]]
    assert bases["descendant_native"] == [phi(P) for P in bases["pullback"]]
    for name, base in bases.items():
        E = E1 if name in ("transported", "descendant_native") else E0
        assert all(E.on_curve(P) and E.mul(P, 421) is None for P in base)
    outcomes = {}
    for name, result in d["variants"].items():
        E = E1 if name in ("transported", "descendant_native") else E0
        gen, target = (G1, Q1) if E is E1 else (G, Q)
        base = bases[name]
        lookup = {}
        for a in range(16):
            for b in range(a, 16):
                lookup.setdefault(E.add(base[a], base[b]), []).append((a, b))
        rows, hits, first_rank = [], 0, None
        assert len(result["cases"]) == 512
        for i, case in enumerate(result["cases"]):
            assert case["i"] == i
            u, v = case["u"], case["v"]
            T = E.add(E.mul(gen, u), E.mul(target, v))
            assert point(case["target"]) == T
            if E is E1:
                assert T == phi(E0.add(E0.mul(G, u), E0.mul(Q, v)))
            matches = lookup.get(T, [])
            assert len(matches) == case["witness_count"]
            assert (list(matches[0]) if matches else None) == case["first_witness"]
            if matches:
                hits += 1
            if matches and first_rank is None:
                a, b = matches[0]
                row = [0] * 17
                row[a] += 1
                row[b] += 1
                row[-1] = -v
                rows.append(row)
                rank = independent_rank(rows, 421)
                assert case["independent_before_rank_stop"] == (
                    rank > independent_rank(rows[:-1], 421))
                if rank == 17:
                    first_rank = i + 1
            else:
                assert not case["independent_before_rank_stop"]
        assert result["hits"] == hits and result["rank"] == 17
        assert result["first_full_rank_attempt"] == first_rank
        assert result["solution"][-1] == k
        assert all(E.mul(gen, value) == P
                   for value, P in zip(result["solution"][:-1], base))
        outcomes[name] = {"hits": hits, "first_full_rank_attempt": first_rank}
    assert outcomes["original"] == outcomes["transported"]
    assert outcomes["descendant_native"] == outcomes["pullback"]
    print(json.dumps({"status": "PASS", "cases_checked": 4 * 512,
                      "outcomes": outcomes}, indent=2))


if __name__ == "__main__":
    main()

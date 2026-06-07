#!/usr/bin/env python3
"""
Phase 1 (Mode A / Mode D) - shallow 2-isogeny crawl around P-256.

Builds the 2-isogeny graph by BFS using exact Velu steps (so every node is an
honest same-order neighbour of P-256: #E'(F_p) = n by construction). Each edge
is a rational 2-isogeny; each node records its j-invariant, model (a,b), the
number of rational 2-torsion points (= 2-isogeny degree in the graph), and the
local volcano signal v_2(Frobenius disc).

Because v_2(Delta_pi)=0 for P-256, the 2-volcano has height 0: P-256 is on the
crater, all 2-isogenies are horizontal, and the crawl walks the crater cycle.

Outputs (Phase 12 artifacts):
    curves.jsonl   one node per line
    edges.jsonl    one directed edge per line
"""

from __future__ import annotations
import json
import hashlib
import sys
from collections import deque

from ecfp import Curve, INF, two_torsion_x, velu_2_isogeny, order_is
import params as PP


def node_id(j: int) -> str:
    return hashlib.sha256(str(j).encode()).hexdigest()[:16]


def v_ell(m: int, ell: int) -> int:
    v, m = 0, abs(m)
    while m and m % ell == 0:
        m //= ell
        v += 1
    return v


def crawl(max_depth: int = 4, verify_order_first: int = 6):
    p, n = PP.P, PP.N
    t = p + 1 - n
    disc = t * t - 4 * p
    v2 = v_ell(disc, 2)

    E0 = Curve(PP.A, PP.B, p)
    j0 = E0.j_invariant()

    nodes = {}        # j -> record
    edges = []        # list of edge dicts
    seen_curves = {}  # j -> Curve (a representative model)

    start = node_id(j0)
    q = deque()
    q.append((j0, E0, 0, "root", []))
    order_checks = []

    while q:
        j, E, depth, path, ellpath = q.popleft()
        nid = node_id(j)
        roots = two_torsion_x(E)
        if nid not in nodes:
            # verify the same-order invariant on the first few nodes
            verified = None
            if len(order_checks) < verify_order_first:
                verified = order_is(E, n, p, trials=6, seed=12345 + depth)
                order_checks.append((j, verified))
            nodes[nid] = {
                "node_id": nid,
                "j": j,
                "depth": depth,
                "path": path,
                "ell_path": ellpath,
                "curve_model": [E.a, E.b],
                "num_2_torsion": len(roots),
                "two_isogeny_degree": len(roots),
                "order_verified": verified,
                "volcano_signal": {
                    "ell": 2,
                    "v_ell_discriminant": v2,
                    "direction": "horizontal" if v2 == 0 else "unknown",
                },
            }
            seen_curves[j] = E

        if depth >= max_depth:
            continue

        for x0 in roots:
            E2 = velu_2_isogeny(E, x0)
            j2 = E2.j_invariant()
            edges.append({
                "from": nid,
                "to": node_id(j2),
                "ell": 2,
                "kernel_x": x0,
            })
            if node_id(j2) not in nodes:
                q.append((j2, E2, depth + 1, path + f"->2", ellpath + [2]))

    return nodes, edges, order_checks, v2, j0


def main():
    max_depth = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    nodes, edges, order_checks, v2, j0 = crawl(max_depth)

    # dedup edges (undirected multiplicity) for a clean summary
    undirected = set()
    for e in edges:
        undirected.add(frozenset((e["from"], e["to"])) if e["from"] != e["to"]
                       else (e["from"],))

    degree_hist = {}
    for nd in nodes.values():
        d = nd["two_isogeny_degree"]
        degree_hist[d] = degree_hist.get(d, 0) + 1

    summary = {
        "root_j": j0,
        "max_depth": max_depth,
        "ell": 2,
        "v2_frobenius_disc": v2,
        "volcano_height_ell2": v2 // 2,  # height = v_ell(conductor) <= v_ell(disc)/2
        "num_nodes": len(nodes),
        "num_directed_edges": len(edges),
        "num_undirected_edges": len(undirected),
        "two_isogeny_degree_histogram": {str(k): v for k, v in sorted(degree_hist.items())},
        "order_verification_sample": [
            {"j": j, "same_order_as_p256": ok} for j, ok in order_checks
        ],
        "all_sampled_nodes_same_order": all(ok for _, ok in order_checks),
    }
    print(json.dumps(summary, indent=2))

    with open("curves.jsonl", "w") as f:
        for nd in nodes.values():
            f.write(json.dumps(nd) + "\n")
    with open("edges.jsonl", "w") as f:
        for e in edges:
            f.write(json.dumps(e) + "\n")

    print(f"\nwrote {len(nodes)} nodes -> curves.jsonl", file=sys.stderr)
    print(f"wrote {len(edges)} edges -> edges.jsonl", file=sys.stderr)


if __name__ == "__main__":
    main()

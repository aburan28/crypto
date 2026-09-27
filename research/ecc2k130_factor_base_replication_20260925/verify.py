#!/usr/bin/env python3
"""Independent full-point, witness, rank and scalar replay of the frozen cells."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import random
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "research/ecc2k130_relations"))
sys.path.insert(0, str(ROOT / "research/ecc2k130_oriented_transport_20260924"))
sys.path.insert(0, str(ROOT / "research/ecc2k130_dual_transport_20260925"))
from fastfield import FastGF2m  # noqa: E402
from relations import Koblitz  # noqa: E402
from oriented_velu import BinaryVeluMap  # noqa: E402
from dual_transport import DualTransport  # noqa: E402

R, N, COFACTOR, SIZE, ATTEMPTS = 421, 21, 4988, 16, 512
VARIANTS = ("original", "transported", "descendant_native", "pullback")


def point(value):
    return None if value is None else tuple(value)


def sha_int(value):
    return int.from_bytes(hashlib.sha256(value.encode()).digest(), "big")


def canonical(E, P):
    if P is None:
        return None
    orbit = []
    q = P
    for _ in range(N):
        orbit.extend((q, E.neg(q)))
        q = (E.F.sqr(q[0]), E.F.sqr(q[1]))
    assert q == P
    return min(orbit)


def independent_orbits(E, G):
    groups = {}
    for k in range(1, R):
        P = E.mul(G, k)
        groups.setdefault(canonical(E, P), set()).add(P)
    assert len(groups) == 10 and all(len(v) == 42 for v in groups.values())
    return sorted(groups)


def complement_on_twist(Twist, chosen_x):
    rng = random.Random(2026092407)
    for _ in range(4096):
        for P in Twist.points_over(rng.getrandbits(N)):
            H = Twist.mul(P, 2094358 // 343)
            if H is None:
                continue
            while Twist.mul(H, 7) is not None:
                H = Twist.mul(H, 7)
            xs = tuple(sorted(Twist.mul(H, j)[0] for j in range(1, 4)))
            if xs != chosen_x:
                return H
    raise AssertionError("complement not found")


def replay_base(E, seed, role, source, D, quota):
    out, back, used_x, unique = [], [], set(), set()
    counts = Counter()
    for trial in range(4096):
        x = sha_int(f"degree7-base-v1|{seed}|{role}|{trial}") & ((1 << N) - 1)
        if x in used_x:
            continue
        used_x.add(x)
        for P in E.points_over(x):
            Q = E.mul(P, COFACTOR)
            if Q is None or Q in unique:
                continue
            unique.add(Q)
            B = Q if role == "source" else source.mul(D.dual(Q), pow(7, -1, R))
            rep = canonical(source, B)
            if rep in quota and counts[rep] < quota[rep]:
                out.append(Q)
                back.append(B)
                counts[rep] += 1
                if len(out) == SIZE:
                    assert all(counts[k] == quota[k] for k in quota)
                    return out, back, trial + 1
    raise AssertionError("incomplete independent base replay")


def independent_rank(rows, modulus=R, augmented=False):
    width = SIZE + 1 + int(augmented)
    matrix = [list(row[:width]) for row in rows]
    pivots = []
    cursor = 0
    for column in range(width):
        found = next((i for i in range(cursor, len(matrix))
                      if matrix[i][column] % modulus), None)
        if found is None:
            continue
        matrix[cursor], matrix[found] = matrix[found], matrix[cursor]
        inv = pow(matrix[cursor][column], -1, modulus)
        matrix[cursor] = [x * inv % modulus for x in matrix[cursor]]
        for i in range(len(matrix)):
            if i == cursor:
                continue
            coeff = matrix[i][column] % modulus
            if coeff:
                matrix[i] = [(x - coeff * y) % modulus
                             for x, y in zip(matrix[i], matrix[cursor])]
        pivots.append(column)
        cursor += 1
        if cursor == len(matrix):
            break
    return len(pivots), matrix, pivots


def independent_solution(rows):
    rank, matrix, pivots = independent_rank(rows, augmented=True)
    assert rank == SIZE + 1 and pivots == list(range(SIZE + 1))
    return [matrix[i][-1] for i in range(SIZE + 1)]


def pair_table(E, base):
    sums = {}
    for i in range(SIZE):
        for j in range(i, SIZE):
            sums.setdefault(E.add(base[i], base[j]), []).append([i, j])
    return sums


def replay_stream(E, G, Q, allowed, label):
    accepted = []
    draw = 0
    for draw in range(1, 8193):
        trial = draw - 1
        u = sha_int(f"degree7-target-v1|{label}|{trial}|u") % R
        v = 1 + sha_int(f"degree7-target-v1|{label}|{trial}|v") % (R - 1)
        T = E.add(E.mul(G, u), E.mul(Q, v))
        if T is not None and canonical(E, T) in allowed:
            accepted.append((u, v, T))
            if len(accepted) == ATTEMPTS:
                return accepted, draw
    raise AssertionError("target stream incomplete")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    assert not args.out.exists()
    raw = json.loads(args.input.read_text())
    for name, relative in {
        "protocol": "research/ecc2k130_factor_base_replication_20260925/PROTOCOL.md",
        "runner": "research/ecc2k130_factor_base_replication_20260925/run.py",
        "pilot": "research/ecc2k130_factor_base_pilot_20260924/run.py",
        "pilot_result": "research/ecc2k130_factor_base_pilot_20260924/results_final.json",
        "dual": "research/ecc2k130_dual_transport_20260925/dual_transport.py",
        "field": "research/ecc2k130_relations/fastfield.py",
        "curve": "research/ecc2k130_relations/relations.py",
        "velu": "research/ecc2k130_oriented_transport_20260924/oriented_velu.py",
    }.items():
        assert hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() == raw["source_sha256"][name]
    F = FastGF2m(N, 0x200005)
    E0, Twist = Koblitz(F, 0, 1), Koblitz(F, 1, 1)
    G, Q = point(raw["challenge"]["G"]), point(raw["challenge"]["Q"])
    secret = sha_int("degree7-replication-secret-v1") % R
    assert 0 < secret < R and raw["challenge"]["secret_audit_only"] == secret
    assert E0.mul(G, R) is None and E0.mul(G, secret) == Q
    old = json.loads((ROOT / "research/ecc2k130_factor_base_pilot_20260924/results_final.json").read_text())
    H = point(old["geometry"]["selected_generator"])
    selected_x = tuple(old["geometry"]["selected_kernel_x"])
    assert tuple(raw["geometry"]["selected_kernel_x"]) == selected_x
    phi = BinaryVeluMap.from_generator(E0, Twist, H, 7)
    D = DualTransport(E0, Twist, H, complement_on_twist(Twist, selected_x), 7, G)
    E1 = Koblitz(F, phi.codomain.a, phi.codomain.b)
    assert E1.b == raw["geometry"]["selected_codomain_b"]
    assert D.compose(G) == E0.mul(G, 7)
    reps = independent_orbits(E0, G)
    assert [list(p) for p in reps] == raw["geometry"]["orbit_representatives"]
    quota = {k: 4 for k in reps[:4]}
    assert set(reps[:5]).isdisjoint(reps[5:])
    replay_bases = {}
    for seed in (2026092511, 2026092512):
        B0, _, trials0 = replay_base(E0, seed, "source", E0, D, quota)
        B1, back, trials1 = replay_base(E1, seed, "leaf", E0, D, quota)
        actual = {"original": B0, "transported": [phi(P) for P in B0],
                  "descendant_native": B1, "pullback": back}
        saved = {name: [point(p) for p in raw["bases"][str(seed)][name]]
                 for name in VARIANTS}
        assert actual == saved
        assert trials0 == raw["base_meta"][str(seed)]["source"]["candidate_trials"]
        assert trials1 == raw["base_meta"][str(seed)]["native"]["candidate_trials"]
        for a, b in zip(back, B1):
            assert phi(a) == b and E0.mul(a, R) is None and E1.mul(b, R) is None
        replay_bases[str(seed)] = actual
    assert replay_bases["2026092511"]["original"] != replay_bases["2026092512"]["original"]
    G1, Q1 = phi(G), phi(Q)
    streams = {}
    for label, allowed in (("A", set(reps[:5])), ("B", set(reps[5:]))):
        targets, draw = replay_stream(E0, G, Q, allowed, label)
        assert draw == raw["target_streams"][label]["meta"]["draws"]
        assert [[u, v] for u, v, _ in targets] == raw["target_streams"][label]["coefficients"]
        assert [reps.index(canonical(E0, T)) for _, _, T in targets] == raw["target_streams"][label]["point_orbits"]
        streams[label] = targets
    cells = {}
    for seed in replay_bases:
        cells[seed] = {}
        for label in ("A", "B"):
            cells[seed][label] = {}
            for name in VARIANTS:
                leaf = name in ("transported", "descendant_native")
                E = E1 if leaf else E0
                base = replay_bases[seed][name]
                table = pair_table(E, base)
                saved = raw["variants"][seed][label][name]
                independent_rows, all_rows = [], []
                first_rank, hits, dependent = None, 0, 0
                assert len(saved["cases"]) == ATTEMPTS
                for index, ((u, v, T0), case) in enumerate(zip(streams[label], saved["cases"])):
                    T = E.add(E.mul(G1, u), E.mul(Q1, v)) if leaf else T0
                    if leaf:
                        assert phi(T0) == T
                    witnesses = table.get(T, [])
                    assert len(witnesses) == case[0]
                    assert (witnesses[0] if witnesses else None) == case[1]
                    if not witnesses:
                        assert case[2] is False
                        continue
                    hits += 1
                    expected_independent = False
                    if first_rank is None:
                        a, b = witnesses[0]
                        row = [0] * (SIZE + 1)
                        row[a] += 1
                        row[b] += 1
                        row[-1] = -v % R
                        augmented = row + [u]
                        old_rank = independent_rank(all_rows)[0]
                        all_rows.append(augmented)
                        new_rank = independent_rank(all_rows)[0]
                        assert independent_rank(all_rows, augmented=True)[0] == new_rank
                        expected_independent = new_rank > old_rank
                        if expected_independent:
                            independent_rows.append(augmented)
                        else:
                            dependent += 1
                        if new_rank == SIZE + 1:
                            first_rank = index + 1
                    assert case[2] is expected_independent
                assert hits == saved["hits"]
                assert first_rank == saved["first_full_rank_attempt"]
                assert dependent == saved["dependent_before_stop"]
                assert independent_rank(all_rows)[0] == saved["rank"]
                assert len(table) == saved["distinct_pair_sums"]
                assert sum(map(len, table.values())) == saved["pair_entries"] == 136
                if first_rank is not None:
                    logs = independent_solution(independent_rows)
                    assert logs == saved["solution"]
                    assert logs[-1] == secret and E.mul(G1 if leaf else G, logs[-1]) == (Q1 if leaf else Q)
                    assert all(E.mul(G1 if leaf else G, d) == P for d, P in zip(logs[:-1], base))
                cells[seed][label][name] = {"rank": saved["rank"], "first_rank": first_rank,
                                            "hits": hits, "witnesses_replayed": ATTEMPTS}
    # Reassemble the cold prefix from disjoint measured phases independently of
    # the producer's stored total. Group-add counts are diagnostics overlapping
    # field arithmetic; they are never added as a second cost unit.
    def sum_costs(chunks):
        total = Counter()
        for chunk in chunks:
            assert all(type(value) is int and value >= 0 for value in chunk.values())
            total.update(chunk)
        return dict(sorted(total.items()))
    ledger = raw["phase_costs"]
    for seed in raw["variants"]:
        for label in ("A", "B"):
            stream = raw["target_streams"][label]
            assert stream["meta"]["draws"] == 512 + sum(stream["meta"]["rejections"].values())
            assert len(raw["target_prefix_costs"][label]) == 512
            assert len(raw["codomain_target_prefix_costs"][label]) == 512
            for name in VARIANTS:
                saved = raw["variants"][seed][label][name]
                limit = saved["first_full_rank_attempt"] or 512
                costs = saved["costs"]
                components = [ledger["field_setup"], ledger["generator"],
                              ledger["challenge"], ledger["source_orbit_universe"],
                              *raw["target_prefix_costs"][label][:limit],
                              costs["pair_table"], costs["scan_to_rank"],
                              costs["recovery_verify"]]
                if name == "original":
                    components += [ledger[f"source_base_{seed}"]]
                elif name == "transported":
                    components += [ledger["kernel_search"], ledger["forward_setup"],
                                   ledger[f"source_base_{seed}"],
                                   ledger[f"transport_base_{seed}"],
                                   ledger["codomain_generator_challenge"],
                                   *raw["codomain_target_prefix_costs"][label][:limit]]
                elif name == "descendant_native":
                    components += [ledger["kernel_search"], ledger["dual_setup"],
                                   ledger[f"native_base_{seed}"],
                                   ledger["codomain_generator_challenge"],
                                   *raw["codomain_target_prefix_costs"][label][:limit]]
                else:
                    components += [ledger["kernel_search"], ledger["dual_setup"],
                                   ledger[f"native_base_{seed}"]]
                assert sum_costs(components) == raw["cold_cost_to_rank_or_512"][seed][label][name]
    output = {"schema": "degree7-disjoint-orbit-independent-replay-v1", "status": "PASS",
              "input_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
              "verifier_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "cells": cells, "total_cases": 2 * 2 * 4 * ATTEMPTS,
              "base_replays": 4, "target_candidate_replays": sum(
                  raw["target_streams"][h]["meta"]["draws"] for h in ("A", "B"))}
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": output["status"], "cases": output["total_cases"],
                      "target_candidate_replays": output["target_candidate_replays"],
                      "output": str(args.out)}, indent=2))


if __name__ == "__main__":
    main()

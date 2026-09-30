#!/usr/bin/env python3
"""Fresh-process brute-triple, rank, scalar and bit-polynomial replay."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import resource
import signal
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "research/ecc2k130_factor_base_pilot_20260924"))
sys.path.insert(0, str(ROOT / "research/ecc2k130_dual_transport_20260925"))
import run as pilot  # noqa: E402
from dual_transport import DualTransport  # noqa: E402

CONFIG = HERE / "CONFIG.json"
VARIANTS = ("original", "transported", "descendant_native", "pullback")
R = 421


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def hash_int(label: str) -> int:
    return int.from_bytes(hashlib.sha256(label.encode()).digest(), "big")


def point(value):
    return None if value is None else tuple(value)


def save(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def peak_rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def checked_config() -> dict:
    config = json.loads(CONFIG.read_text())
    assert config["schema"] == "ecc2k130-degree7-m3-four-policy-v1"
    assert config["field_degree"] == 21 and config["base_size"] == 8
    for relative, expected in config["inputs_sha256"].items():
        assert sha(ROOT / relative) == expected, relative
    lock = json.loads((HERE / "FROZEN.json").read_text())
    assert lock["schema"] == "ecc2k130-degree7-m3-four-policy-source-lock-v1"
    for relative, expected in lock["sha256"].items():
        assert sha(ROOT / relative) == expected, relative
    return config


class ReferenceField:
    """Direct polynomial arithmetic, independent of FastGF2m's tables."""

    degree = 21
    modulus = 0x200005
    mask = (1 << degree) - 1

    def mul(self, a: int, b: int) -> int:
        result = 0
        while b:
            if b & 1:
                result ^= a
            b >>= 1
            a <<= 1
            if a & (1 << self.degree):
                a ^= self.modulus
        return result & self.mask

    def sqr(self, a: int) -> int:
        return self.mul(a, a)

    def inv(self, a: int) -> int:
        assert a != 0
        value, power, exponent = 1, a, (1 << self.degree) - 2
        while exponent:
            if exponent & 1:
                value = self.mul(value, power)
            power = self.sqr(power)
            exponent >>= 1
        assert self.mul(a, value) == 1
        return value


class ReferenceCurve:
    """Separate bit-polynomial curve law for selected full-point controls."""

    def __init__(self, b: int):
        self.F, self.b = ReferenceField(), b

    def on_curve(self, P) -> bool:
        if P is None:
            return True
        x, y = P
        F = self.F
        return (F.sqr(y) ^ F.mul(x, y)) == (
            F.mul(x, F.sqr(x)) ^ self.b)

    @staticmethod
    def neg(P):
        return None if P is None else (P[0], P[1] ^ P[0])

    def add(self, P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        F = self.F
        x1, y1 = P
        x2, y2 = Q
        if x1 == x2:
            if y2 == (y1 ^ x1) or x1 == 0:
                return None
            slope = x1 ^ F.mul(y1, F.inv(x1))
            x3 = F.sqr(slope) ^ slope
            return x3, F.sqr(x1) ^ F.mul(slope ^ 1, x3)
        slope = F.mul(y1 ^ y2, F.inv(x1 ^ x2))
        x3 = F.sqr(slope) ^ slope ^ x1 ^ x2
        return x3, F.mul(slope, x1 ^ x3) ^ x3 ^ y1

    def mul(self, P, k: int):
        assert k >= 0
        acc, multiple = None, P
        while k:
            if k & 1:
                acc = self.add(acc, multiple)
            multiple = self.add(multiple, multiple)
            k >>= 1
        return acc


def canonical(E, P):
    assert P is not None
    seen, current = [], P
    for _ in range(21):
        seen.extend((current, E.neg(current)))
        current = (E.F.sqr(current[0]), E.F.sqr(current[1]))
    assert current == P
    return min(seen)


def independent_orbits(E, G):
    by_point, members = {}, defaultdict(set)
    for k in range(1, R):
        P = E.mul(G, k)
        rep = canonical(E, P)
        by_point[P] = rep
        members[rep].add(P)
    representatives = sorted(members)
    assert len(representatives) == 10
    assert all(len(members[rep]) == 42 for rep in representatives)
    return representatives, by_point, members


def replay_base(E, seed: int, role: str, source, dual, by_point,
                quotas: dict, config: dict):
    out, back, used_x, unique = [], [], set(), set()
    counts, meta = Counter(), Counter()
    for trial in range(config["max_base_x_trials_per_curve_seed"]):
        x = hash_int(f"degree7-m3-four-policy-base-v1|{seed}|{role}|{trial}") & (
            (1 << 21) - 1)
        if x in used_x:
            meta["repeated_x"] += 1
            continue
        used_x.add(x)
        meta["x_scanned"] += 1
        for P in E.points_over(x):
            meta["raw_projected"] += 1
            Q = E.mul(P, pilot.COFACTOR)
            if Q is None:
                meta["infinity"] += 1
                continue
            if Q in unique:
                meta["duplicate"] += 1
                continue
            unique.add(Q)
            B = Q if role == "source" else source.mul(
                dual.dual(Q), pow(pilot.ELL, -1, R))
            assert B is not None
            if role == "leaf":
                meta["dual_pullbacks"] += 1
            orbit = by_point.get(B)
            assert orbit is not None
            meta["classified"] += 1
            if orbit in quotas and counts[orbit] < quotas[orbit]:
                out.append(Q)
                back.append(B)
                counts[orbit] += 1
                if len(out) == config["base_size"]:
                    meta["candidate_trials"] = trial + 1
                    assert all(counts[key] == quotas[key] for key in quotas)
                    return out, back, dict(meta), dict(counts)
            else:
                meta["out_of_quota"] += 1
    raise AssertionError(f"independent {role} base incomplete")


def target_stream(E, G, Q, by_point, allowed, label, config):
    targets, rejections = [], Counter()
    for trial in range(config["max_target_draws_per_holdout"]):
        u = hash_int(f"degree7-m3-four-policy-target-v1|{label}|{trial}|u") % R
        v = 1 + hash_int(f"degree7-m3-four-policy-target-v1|{label}|{trial}|v") % (R - 1)
        T = E.add(E.mul(G, u), E.mul(Q, v))
        orbit = by_point.get(T)
        if orbit is None:
            rejections["infinity"] += 1
        elif orbit not in allowed:
            rejections["other_orbit"] += 1
        else:
            targets.append((u, v, T))
            if len(targets) == config["accepted_targets_per_holdout"]:
                return targets, {"draws": trial + 1, "accepted": len(targets),
                                 "rejections": dict(rejections)}
    raise AssertionError(f"independent holdout {label} incomplete")


def brute_triples(E, base):
    triples = defaultdict(list)
    pair_sums = set()
    for a in range(len(base)):
        for b in range(a, len(base)):
            pair_sums.add(E.add(base[a], base[b]))
    for k in range(len(base)):
        for a in range(len(base)):
            for b in range(a, len(base)):
                T = E.add(E.add(base[a], base[b]), base[k])
                triples[T].append((k, a, b))
    for witnesses in triples.values():
        witnesses.sort()
    assert len(triples) <= 120
    return triples, len(pair_sums)


def independent_linear_system(rows: list[tuple[list[int], int]], width: int):
    matrix = [[value % R for value in coeffs] + [rhs % R]
              for coeffs, rhs in rows]
    pivots = []
    cursor = 0
    for col in range(width):
        row = next((j for j in range(cursor, len(matrix)) if matrix[j][col]), None)
        if row is None:
            continue
        matrix[cursor], matrix[row] = matrix[row], matrix[cursor]
        scale = pow(matrix[cursor][col], -1, R)
        matrix[cursor] = [x * scale % R for x in matrix[cursor]]
        for j in range(len(matrix)):
            if j != cursor and matrix[j][col]:
                factor = matrix[j][col]
                matrix[j] = [(a - factor * b) % R
                             for a, b in zip(matrix[j], matrix[cursor])]
        pivots.append(col)
        cursor += 1
        if cursor == len(matrix):
            break
    assert all(any(row[:width]) or row[-1] == 0 for row in matrix)
    solution = None
    if len(pivots) == width:
        solution = [0] * width
        for i, col in enumerate(pivots):
            solution[col] = matrix[i][-1]
        assert all(sum(a * b for a, b in zip(coeffs, solution)) % R == rhs % R
                   for coeffs, rhs in rows)
    return len(pivots), solution


def sum_costs(parts):
    total = Counter()
    for part in parts:
        assert all(type(value) is int and value >= 0 for value in part.values())
        total.update(part)
    return dict(sorted(total.items()))


def replay(folder: Path) -> dict:
    config = checked_config()
    result = json.loads((folder / "result.json").read_text())
    assert result["schema"] == "ecc2k130-degree7-m3-four-policy-result-v1"
    assert result["status"] in ("PASS_PANEL", "NEGATIVE_PANEL")
    assert result["config_sha256"] == sha(CONFIG)
    assert result["frozen_sha256"] == sha(HERE / "FROZEN.json")
    F = pilot.FastGF2m(21, pilot.IRR)
    E0 = pilot.Koblitz(F, 0, 1)
    twist = pilot.Koblitz(F, 1, 1)
    lines, selected, line_trials = pilot.order_seven_lines(F, twist)
    saved = json.loads((ROOT / "research/ecc2k130_factor_base_pilot_20260924/results_final.json").read_text())
    assert list(selected) == saved["geometry"]["selected_kernel_x"]
    kernel = lines[selected]
    phi = pilot.BinaryVeluMap.from_generator(E0, twist, kernel["generator"], 7)
    complement = next(lines[key]["generator"] for key in sorted(lines) if key != selected)
    G, generator_trials = pilot.choose_generator(E0)
    secret = 1 + hash_int(config["labels"]["secret"]) % (R - 1)
    Q = E0.mul(G, secret)
    D = DualTransport(E0, twist, kernel["generator"], complement, 7, G)
    E1 = pilot.Koblitz(F, phi.codomain.a, phi.codomain.b)
    assert E1.a == 0, "reference point law handles the a=0 codomain"
    G1, Q1 = phi(G), phi(Q)
    assert D.compose(G) == E0.mul(G, 7)
    assert E1.mul(G1, secret) == Q1
    reps, by_point, members = independent_orbits(E0, G)
    assert result["geometry"] == {
        "selected_kernel_x": list(selected), "codomain_b": E1.b,
        "line_trials": line_trials, "generator_trials": generator_trials,
        "orbit_representatives": [list(rep) for rep in reps],
        "orbit_sizes": [len(members[rep]) for rep in reps]}
    assert result["challenge"] == {"G": list(G), "Q": list(Q),
                                   "secret_audit_only": secret}
    quotas = {rep: config["base_points_per_source_signed_orbit"]
              for rep in reps[:4]}
    bases = {}
    reference_checks = 0
    refs = {False: ReferenceCurve(1), True: ReferenceCurve(E1.b)}
    for curve, points in ((refs[False], (G, Q)), (refs[True], (G1, Q1))):
        for P in points:
            assert curve.on_curve(P) and curve.mul(P, R) is None
            reference_checks += 1
    for seed in config["base_seeds"]:
        key = str(seed)
        B0, _, meta0, occ0 = replay_base(E0, seed, "source", E0, D,
                                          by_point, quotas, config)
        B1, back, meta1, occ1 = replay_base(E1, seed, "leaf", E0, D,
                                             by_point, quotas, config)
        bases[key] = {"original": B0, "transported": [phi(P) for P in B0],
                      "descendant_native": B1, "pullback": back}
        assert result["base_meta"][key] == {"source": meta0, "native": meta1}
        assert result["occupancy"][key] == {
            "source": {str(rep): occ0[rep] for rep in quotas},
            "native": {str(rep): occ1[rep] for rep in quotas}}
        assert all([list(P) for P in values] == result["bases"][key][name]
                   for name, values in bases[key].items())
        assert all(phi(P) == Q for P, Q in zip(back, B1))
        for name, points in bases[key].items():
            ref = refs[name in ("transported", "descendant_native")]
            for P in points:
                assert ref.on_curve(P) and ref.mul(P, R) is None
                reference_checks += 1
    assert bases[str(config["base_seeds"][0])]["original"] != bases[
        str(config["base_seeds"][1])]["original"]
    targets = {}
    for label, allowed in (("A", set(reps[:5])), ("B", set(reps[5:]))):
        source, meta = target_stream(E0, G, Q, by_point, allowed, label, config)
        leaf = [(u, v, E1.add(E1.mul(G1, u), E1.mul(Q1, v)))
                for u, v, _ in source]
        assert all(phi(T[2]) == U[2] for T, U in zip(source, leaf))
        assert result["target_streams"][label] == {
            "coefficients": [[u, v] for u, v, _ in source], "meta": meta,
            "point_orbits": [reps.index(by_point[T]) for _, _, T in source]}
        targets[label], targets[f"{label}_leaf"] = source, leaf
        for cases, ref, gen, challenge in ((source, refs[False], G, Q),
                                           (leaf, refs[True], G1, Q1)):
            for u, v, T in cases[:2] + cases[-2:]:
                assert ref.on_curve(T)
                assert ref.add(ref.mul(gen, u), ref.mul(challenge, v)) == T
                reference_checks += 1
    cells, total_cases = {}, 0
    for seed in config["base_seeds"]:
        key = str(seed)
        cells[key] = {}
        for label in config["holdouts"]:
            cells[key][label] = {}
            trajectories = {}
            for name in VARIANTS:
                leaf = name in ("transported", "descendant_native")
                E, gen, challenge = (E1, G1, Q1) if leaf else (E0, G, Q)
                target_cases = targets[f"{label}_leaf"] if leaf else targets[label]
                base = bases[key][name]
                path = folder / "cells" / key / label / f"{name}.json"
                assert sha(path) == result["cell_sha256"][key][label][name]
                case_result = json.loads(path.read_text())
                assert {k: v for k, v in case_result.items() if k != "cases"} == (
                    result["variants"][key][label][name])
                assert len(case_result["cases"]) == len(target_cases) == 512
                triples, pair_sums = brute_triples(E, base)
                assert result["triple_support_cardinality"][key][label][name] == len(triples)
                assert case_result["pair_entries"] == 36
                assert case_result["distinct_pair_sums"] == pair_sums
                rows = []
                rank = hits = dependent = 0
                first_rank = solution = None
                for index, ((u, v, T), record) in enumerate(zip(
                        target_cases, case_result["cases"])):
                    witnesses = triples.get(T, ())
                    first = witnesses[0] if witnesses else None
                    witness = [first[1], first[2], first[0]] if first else None
                    multiplicity = (sum(1 for item in witnesses if item[0] == first[0])
                                    if first else 0)
                    tested = first[0] + 1 if first else len(base)
                    assert record["i"] == index and record["u"] == u and record["v"] == v
                    assert record["target"] == list(T)
                    assert record["first_witness"] == witness
                    assert record["selected_pair_multiplicity"] == multiplicity
                    assert record["third_candidates_tested"] == tested
                    assert all(type(value) is int and value >= 0
                               for value in record["scan_cost"].values())
                    hits += bool(first)
                    independent = False
                    if first and first_rank is None:
                        coeffs = [0] * (len(base) + 1)
                        for j in witness:
                            coeffs[j] += 1
                        coeffs[-1] = -v % R
                        rows.append((coeffs, u))
                        new_rank, candidate = independent_linear_system(rows, len(base) + 1)
                        independent = new_rank > rank
                        dependent += not independent
                        rank = new_rank
                        if rank == len(base) + 1:
                            first_rank, solution = index + 1, candidate
                    assert record["independent_before_rank_stop"] is independent
                    assert record["rank_after"] == rank
                    total_cases += 1
                assert case_result["hits"] == hits
                assert case_result["dependent_before_stop"] == dependent
                assert case_result["rank"] == rank
                assert case_result["first_full_rank_attempt"] == first_rank
                assert case_result["solution"] == solution
                assert case_result["verified"] == (solution is not None)
                limit = first_rank or 512
                for metric in ("group_add", "mul", "sqr", "inv"):
                    assert sum(record["scan_cost"].get(metric, 0)
                               for record in case_result["cases"][:limit]) == (
                        case_result["costs"]["scan_to_rank"].get(metric, 0))
                if solution is not None:
                    assert solution[-1] == secret
                    assert E.mul(gen, solution[-1]) == challenge
                    assert all(E.mul(gen, d) == P for d, P in zip(solution[:-1], base))
                    ref = refs[leaf]
                    assert ref.mul(gen, solution[-1]) == challenge
                    assert all(ref.mul(gen, d) == P for d, P in zip(solution[:-1], base))
                    reference_checks += len(base) + 1
                cells[key][label][name] = {"hits": hits, "rank": rank,
                                           "first_full_rank_attempt": first_rank}
                trajectories[name] = [
                    (case["first_witness"], case["rank_after"])
                    for case in case_result["cases"]]
            original, transported, native, pullback = (
                cells[key][label][name] for name in VARIANTS)
            assert original == transported and native == pullback
            assert trajectories["original"] == trajectories["transported"]
            assert trajectories["descendant_native"] == trajectories["pullback"]
            expected_controls = {
                "forward_hit_parity": True, "dual_hit_parity": True,
                "forward_rank_parity": True, "dual_rank_parity": True,
                "forward_trajectory_parity": True, "dual_trajectory_parity": True,
                "all_verified": all(result["variants"][key][label][name]["verified"]
                                    for name in VARIANTS)}
            assert result["controls"][key][label] == expected_controls
    assert total_cases == 16 * 512
    expected_status = "PASS_PANEL" if all(
        result["variants"][key][label][name]["verified"]
        for key in result["variants"] for label in result["variants"][key]
        for name in VARIANTS) else "NEGATIVE_PANEL"
    assert result["status"] == expected_status
    ledger = result["phase_costs"]
    assert len(result["target_prefix_costs"]["A"]) == 512
    assert len(result["target_prefix_costs"]["B"]) == 512
    for key in result["variants"]:
        for label in config["holdouts"]:
            for name in VARIANTS:
                cell = result["variants"][key][label][name]
                limit = cell["first_full_rank_attempt"] or 512
                shared = [ledger["field_setup"], ledger["generator"],
                          ledger["challenge"], ledger["source_orbit_universe"],
                          *result["target_prefix_costs"][label][:limit]]
                chunks = [cell["costs"]["pair_table"],
                          cell["costs"]["scan_to_rank"],
                          cell["costs"]["recovery_verify"]]
                if name == "original":
                    chunks.append(ledger[f"source_base_{key}"])
                elif name == "transported":
                    chunks += [ledger["kernel_search"], ledger["forward_setup"],
                               ledger[f"source_base_{key}"],
                               ledger[f"transport_base_{key}"],
                               ledger["codomain_generator_challenge"],
                               *result["codomain_target_prefix_costs"][label][:limit]]
                elif name == "descendant_native":
                    chunks += [ledger["kernel_search"], ledger["dual_setup"],
                               ledger[f"native_base_{key}"],
                               ledger["codomain_generator_challenge"],
                               *result["codomain_target_prefix_costs"][label][:limit]]
                else:
                    chunks += [ledger["kernel_search"], ledger["dual_setup"],
                               ledger[f"native_base_{key}"]]
                assert sum_costs(shared + chunks) == (
                    result["cold_cost_to_rank_or_512"][key][label][name])
    return {"schema": "ecc2k130-degree7-m3-four-policy-replay-v1",
            "status": "PASS", "panel_status": expected_status,
            "result_sha256": sha(folder / "result.json"),
            "config_sha256": sha(CONFIG),
            "frozen_sha256": sha(HERE / "FROZEN.json"),
            "cells_replayed": 16, "cases_replayed": total_cases,
            "bit_polynomial_point_checks": reference_checks,
            "cells": cells}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite a replay receipt"
    config = json.loads(CONFIG.read_text())
    started = time.perf_counter()
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(
        TimeoutError("frozen verifier wall cap")))
    signal.alarm(config["verifier_wall_cap_seconds"])
    try:
        receipt = replay(args.evidence.resolve())
        assert peak_rss_bytes() <= config["process_rss_cap_bytes"]
        save(args.out, receipt)
        print(json.dumps({"status": receipt["status"],
                          "panel_status": receipt["panel_status"],
                          "cases_replayed": receipt["cases_replayed"]}, sort_keys=True))
    except BaseException as error:
        save(args.out, {"status": "CENSORED" if isinstance(error, (
            TimeoutError, MemoryError)) else "FAIL", "error_type": type(error).__name__,
            "error": str(error), "traceback": traceback.format_exc(),
            "elapsed_seconds": time.perf_counter() - started,
            "peak_rss_bytes": peak_rss_bytes()})
        raise
    finally:
        signal.alarm(0)


if __name__ == "__main__":
    main()

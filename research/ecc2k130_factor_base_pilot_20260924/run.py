#!/usr/bin/env python3
"""Frozen, paired degree-7 factor-base pilot; see PROTOCOL.md before running."""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import random
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "research/ecc2k130_relations"))
sys.path.insert(0, str(ROOT / "research/ecc2k130_oriented_transport_20260924"))
from fastfield import FastGF2m  # noqa: E402
from relations import Koblitz  # noqa: E402
import oriented_velu  # noqa: E402
from oriented_velu import BinaryVeluMap  # noqa: E402

SEED = 2026092407
N = 21
IRR = 0x200005
ELL = 7
ORDER = 2099948
TWIST_ORDER = 2094358
R = 421
COFACTOR = ORDER // R
BASE_SIZE = 16
ATTEMPTS = 512


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode_point(P):
    return None if P is None else [P[0], P[1]]


class Meter:
    def __init__(self):
        self.calls = Counter()
        self.group_adds = 0

    def snapshot(self):
        return dict(self.calls), self.group_adds, time.process_time_ns()

    @staticmethod
    def delta(before, after):
        a, g, t = before
        b, h, u = after
        return {**{key: b.get(key, 0) - a.get(key, 0)
                   for key in sorted(set(a) | set(b))},
                "group_add": h - g, "cpu_ns": u - t}


class CountedField(FastGF2m):
    def __init__(self, meter: Meter, deg=N, irr=IRR):
        self.meter = meter
        super().__init__(deg, irr)

    def _sqr_slow(self, a):
        self.meter.calls["setup_sqr_slow"] += 1
        return super()._sqr_slow(a)

    def mul(self, a, b):
        self.meter.calls["mul"] += 1
        return super().mul(a, b)

    def sqr(self, a):
        self.meter.calls["sqr"] += 1
        return super().sqr(a)

    def inv(self, a):
        self.meter.calls["inv"] += 1
        return super().inv(a)

    def trace(self, a):
        self.meter.calls["trace"] += 1
        return super().trace(a)

    def half_trace(self, a):
        self.meter.calls["half_trace"] += 1
        return super().half_trace(a)


class CountedCurve(Koblitz):
    def __init__(self, meter, field, a=0, b=1):
        self.meter = meter
        super().__init__(field, a=a, b=b)

    def add(self, P, Q):
        self.meter.group_adds += 1
        return super().add(P, Q)


class InternalCountedCurve(Koblitz):
    """Charge group additions inside the map's copied curves."""

    def __init__(self, field, a=0, b=1):
        self.meter = field.meter
        super().__init__(field, a=a, b=b)

    def add(self, P, Q):
        self.meter.group_adds += 1
        return super().add(P, Q)


oriented_velu.Koblitz = InternalCountedCurve


def phase(meter: Meter, ledger: dict, name: str, fn):
    before = meter.snapshot()
    value = fn()
    ledger[name] = meter.delta(before, meter.snapshot())
    return value


def order_seven_lines(F, twist, *, limit=4096):
    rng = random.Random(SEED)
    lines = {}
    for trial in range(limit):
        x = rng.getrandbits(N)
        points = twist.points_over(x)
        if not points:
            continue
        H = twist.mul(points[0], TWIST_ORDER // (ELL ** 3))
        if H is None:
            continue
        while twist.mul(H, ELL) is not None:
            H = twist.mul(H, ELL)
        kernel = [twist.mul(H, j) for j in range(1, (ELL + 1) // 2)]
        assert all(P is not None and twist.on_curve(P) for P in kernel)
        key = tuple(sorted(P[0] for P in kernel))
        if key in lines:
            continue
        t = 0
        for u in key:
            t ^= u
        b = 1 ^ t ^ F.sqr(t)
        lines[key] = {"generator": H, "codomain_b": b, "trial": trial,
                      "kernel_x": list(key)}
        if len(lines) == ELL + 1:
            break
    assert len(lines) == ELL + 1, f"only {len(lines)} kernel lines"
    assert all(twist.mul(v["generator"], ELL) is None for v in lines.values())
    selected_key = next(key for key in sorted(lines) if lines[key]["codomain_b"] != 1)
    return lines, selected_key, trial + 1


def choose_generator(E):
    rng = random.Random(SEED)
    for trial in range(4096):
        x = rng.getrandbits(N)
        for P in E.points_over(x):
            G = E.mul(P, COFACTOR)
            if G is not None:
                assert E.mul(G, R) is None and E.on_curve(G)
                return G, trial + 1
    raise AssertionError("no order-421 generator")


def candidate_base(E):
    base, seen, projected, scanned = [], set(), 0, 0
    for x in range(64):
        scanned += 1
        for P in E.points_over(x):
            projected += 1
            Q = E.mul(P, COFACTOR)
            if Q is None or Q in seen:
                continue
            assert E.mul(Q, R) is None and E.on_curve(Q)
            seen.add(Q)
            base.append(Q)
            if len(base) == BASE_SIZE:
                return base, {"x_scanned": scanned, "raw_projected": projected,
                              "discarded": projected - len(base)}
    raise AssertionError(f"only {len(base)} unique useful base points")


def target_stream(E, G, meter):
    before = meter.snapshot()
    secret = 1 + int.from_bytes(
        hashlib.sha256(f"degree7-secret|{SEED}".encode()).digest(), "big") % (R - 1)
    Q = E.mul(G, secret)
    assert Q is not None and E.mul(Q, R) is None
    challenge_cost = meter.delta(before, meter.snapshot())
    rng = random.Random(SEED)
    targets, per_attempt = [], []
    for i in range(ATTEMPTS):
        before = meter.snapshot()
        u, v = rng.randrange(R), rng.randrange(1, R)
        T = E.add(E.mul(G, u), E.mul(Q, v))
        targets.append((u, v, T))
        per_attempt.append(meter.delta(before, meter.snapshot()))
    return secret, Q, targets, challenge_cost, per_attempt


def transport_stream(phi, G, Q, targets, meter):
    before = meter.snapshot()
    G1, Q1 = phi(G), phi(Q)
    setup = meter.delta(before, meter.snapshot())
    output, per_attempt = [], []
    for u, v, T in targets:
        before = meter.snapshot()
        output.append((u, v, phi(T)))
        per_attempt.append(meter.delta(before, meter.snapshot()))
    return G1, Q1, output, setup, per_attempt


def direct_codomain_targets(curve, G, Q, paired_targets, meter):
    output, per_attempt = [], []
    for u, v, _ in paired_targets:
        before = meter.snapshot()
        T = curve.add(curve.mul(G, u), curve.mul(Q, v))
        output.append((u, v, T))
        per_attempt.append(meter.delta(before, meter.snapshot()))
    return output, per_attempt


def add_costs(*costs):
    total = Counter()
    for cost in costs:
        total.update(cost)
    return dict(sorted(total.items()))


def pair_table(E, base):
    table = {}
    for i, P in enumerate(base):
        for j in range(i, len(base)):
            T = E.add(P, base[j])
            table.setdefault(T, []).append((i, j))
    return table


class RankTracker:
    def __init__(self, width=BASE_SIZE + 1, modulus=R):
        self.width, self.modulus = width, modulus
        self.pivots = {}
        self.ops = Counter()

    def add(self, coefficients, rhs):
        row, q = list(coefficients), self.modulus
        rhs %= q
        for pivot, (lead, const) in sorted(self.pivots.items()):
            f = row[pivot]
            if not f:
                continue
            for j in range(pivot, self.width):
                row[j] = (row[j] - f * lead[j]) % q
                self.ops["mul_mod_r"] += 1
                self.ops["add_mod_r"] += 1
            rhs = (rhs - f * const) % q
            self.ops["mul_mod_r"] += 1
            self.ops["add_mod_r"] += 1
        pivot = next((i for i, a in enumerate(row) if a), None)
        if pivot is None:
            assert rhs == 0, "inconsistent group relation"
            return False
        inv = pow(row[pivot], -1, q)
        self.ops["inv_mod_r"] += 1
        for j in range(pivot, self.width):
            row[j] = row[j] * inv % q
            self.ops["mul_mod_r"] += 1
        rhs = rhs * inv % q
        self.ops["mul_mod_r"] += 1
        self.pivots[pivot] = (row, rhs)
        return True

    def solve(self):
        assert len(self.pivots) == self.width
        x = [0] * self.width
        for p in sorted(self.pivots, reverse=True):
            row, rhs = self.pivots[p]
            for j in range(p + 1, self.width):
                rhs = (rhs - row[j] * x[j]) % self.modulus
                self.ops["mul_mod_r"] += 1
                self.ops["add_mod_r"] += 1
            x[p] = rhs
        return x


def relation_run(E, G, Q, targets, base, *, expected_secret, meter):
    before_pair = meter.snapshot()
    table = pair_table(E, base)
    after_pair = meter.snapshot()
    tracker = RankTracker()
    cases, first_full_rank, hits, dependent = [], None, 0, 0
    at_rank = None
    for i, (u, v, T) in enumerate(targets):
        witnesses = table.get(T, [])
        if witnesses:
            hits += 1
        independent = False
        if witnesses and first_full_rank is None:
            a, b = witnesses[0]
            row = [0] * (BASE_SIZE + 1)
            row[a] += 1
            row[b] += 1
            row[-1] = -v % R
            independent = tracker.add(row, u)
            if not independent:
                dependent += 1
            if len(tracker.pivots) == tracker.width:
                first_full_rank = i + 1
                at_rank = meter.snapshot()
        cases.append({"i": i, "u": u, "v": v, "target": encode_point(T),
                      "witness_count": len(witnesses),
                      "first_witness": list(witnesses[0]) if witnesses else None,
                      "independent_before_rank_stop": independent})
    after_scan = meter.snapshot()
    solution = tracker.solve() if first_full_rank is not None else None
    verified = False
    if solution is not None:
        verified = solution[-1] == expected_secret and E.mul(G, solution[-1]) == Q
        verified &= all(E.mul(G, d) == P for d, P in zip(solution[:-1], base))
    after_verify = meter.snapshot()
    costs = {
        "pair_table": meter.delta(before_pair, after_pair),
        "scan_to_rank": meter.delta(after_pair, at_rank or after_scan),
        "audit_tail": meter.delta(at_rank or after_scan, after_scan),
        "recovery_verify": meter.delta(after_scan, after_verify),
    }
    return {"attempts": ATTEMPTS, "hits": hits, "dependent_before_stop": dependent,
            "costs": costs,
            "first_full_rank_attempt": first_full_rank, "rank": len(tracker.pivots),
            "mod_r_ops": dict(tracker.ops), "solution": solution,
            "verified": bool(verified), "cases": cases,
            "pair_entries": sum(map(len, table.values())),
            "distinct_pair_sums": len(table)}


def direct_paired_velu(F, P, half_kernel):
    """Independent F_(q²) paired-sum reference for ordinary source points."""
    if P is None:
        return None
    x, y = P

    def eadd(a, b):
        return a[0] ^ b[0], a[1] ^ b[1]

    def emul(a, b):
        ac, bd = F.mul(a[0], b[0]), F.mul(a[1], b[1])
        middle = F.mul(a[0] ^ a[1], b[0] ^ b[1]) ^ ac
        return ac ^ bd, middle

    def esqr(a):
        a2, b2 = F.sqr(a[0]), F.sqr(a[1])
        return a2 ^ b2, b2

    X, Y, t = (x, 0), (y, 0), 0
    for u, v in half_kernel:
        assert x != u, "a rational source point cannot share twist kernel x"
        t ^= u
        inv = F.inv(x ^ u)
        pair = []
        for w in (v, v ^ u):
            lam = F.mul(y ^ w, inv), F.mul(u, inv)
            x3 = eadd(eadd(esqr(lam), lam), (x ^ u, 0))
            y3 = eadd(eadd(emul(lam, eadd((x, 0), x3)), x3), (y, 0))
            pair.append((x3, y3))
        X = eadd(X, eadd(pair[0][0], pair[1][0]))
        Y = eadd(Y, eadd(eadd(pair[0][1], pair[1][1]), (u, 0)))
    assert X[1] == 0 and Y[1] == 0
    return X[0], Y[0] ^ t


def comparison_digest(cases):
    return hashlib.sha256(json.dumps(
        [(x["witness_count"], x["independent_before_rank_stop"]) for x in cases],
        separators=(",", ":")).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "preserve prior results; choose a new output"
    meter, ledger = Meter(), {}
    start = meter.snapshot()
    F = CountedField(meter)
    F.frobenius(1, N - 1)  # charge the shared sqrt table before either base scan
    ledger["field_setup"] = meter.delta(start, meter.snapshot())
    E0, Twist = CountedCurve(meter, F, 0, 1), CountedCurve(meter, F, 1, 1)
    assert ORDER % ELL and TWIST_ORDER % (ELL ** 3) == 0
    lines, selected, trials = phase(
        meter, ledger, "kernel_line_search",
        lambda: order_seven_lines(F, Twist))
    record = lines[selected]
    phi = phase(
        meter, ledger, "selected_map_construction",
        lambda: BinaryVeluMap.from_generator(
            E0, Twist, record["generator"], ELL))
    E1 = CountedCurve(meter, F, phi.codomain.a, phi.codomain.b)
    assert phi.codomain.b == record["codomain_b"]
    assert phi(None) is None
    G, g_trials = phase(meter, ledger, "subgroup_generator",
                        lambda: choose_generator(E0))
    secret, Q, targets0, challenge_cost, target_attempt_costs = phase(
        meter, ledger, "target_stream", lambda: target_stream(E0, G, meter))
    B0, B0_meta = phase(meter, ledger, "original_base",
                        lambda: candidate_base(E0))
    B1, B1_meta = phase(meter, ledger, "native_base",
                        lambda: candidate_base(E1))
    G1, Q1, targets1, transport_setup_cost, transport_attempt_costs = phase(
        meter, ledger, "target_transport",
        lambda: transport_stream(phi, G, Q, targets0, meter))
    assert G1 is not None and E1.mul(G1, R) is None
    assert Q1 == E1.mul(G1, secret)
    assert all(T is None or E1.on_curve(T) for _, _, T in targets1)
    targets1_direct, codomain_attempt_costs = phase(
        meter, ledger, "codomain_target_generation",
        lambda: direct_codomain_targets(E1, G1, Q1, targets0, meter))
    assert targets1_direct == targets1, "paired map/direct target disagreement"
    B_transport = phase(meter, ledger, "base_transport",
                        lambda: [phi(P) for P in B0])
    assert len(set(B_transport)) == BASE_SIZE
    assert all(P is not None and E1.on_curve(P) and E1.mul(P, R) is None
               for P in B_transport)

    def inverse_subgroup_table():
        table = {}
        P = None
        for _ in range(R):
            image = phi(P)
            assert image not in table
            table[image] = P
            P = E0.add(P, G)
        assert P is None and len(table) == R
        return table

    inverse = phase(meter, ledger, "toy_inverse_table", inverse_subgroup_table)
    B_pullback = phase(meter, ledger, "base_pullback",
                       lambda: [inverse[P] for P in B1])
    assert len(set(B_pullback)) == BASE_SIZE
    assert all(phi(P) == Q for P, Q in zip(B_pullback, B1))
    assert all(P is not None and E0.mul(P, R) is None for P in B_pullback)

    half_kernel = [Twist.mul(record["generator"], j)
                   for j in range(1, (ELL + 1) // 2)]

    def independent_lift_checks():
        points = [G, Q, *B0, *[T for _, _, T in targets0[:32] if T is not None]]
        differences = []
        for P in points:
            direct = direct_paired_velu(F, P, half_kernel)
            mapped = phi(P)
            if direct != mapped:
                differences.append({"source": encode_point(P),
                                    "direct": encode_point(direct),
                                    "mapped": encode_point(mapped)})
            assert E1.on_curve(direct)
        for P, Qp in zip(points[:16], points[1:17]):
            assert phi(E0.add(P, Qp)) == E1.add(phi(P), phi(Qp))
        return {"points": len(points), "additivity_pairs": 16,
                "differences": differences}

    lift = phase(meter, ledger, "independent_full_lift_audit",
                 independent_lift_checks)
    assert not lift["differences"]

    variants = {}
    for name, curve, gen, challenge, targets, base in [
        ("original", E0, G, Q, targets0, B0),
        ("transported", E1, G1, Q1, targets1_direct, B_transport),
        ("descendant_native", E1, G1, Q1, targets1_direct, B1),
        ("pullback", E0, G, Q, targets0, B_pullback),
    ]:
        variants[name] = phase(
            meter, ledger, f"relation_{name}",
            lambda curve=curve, gen=gen, challenge=challenge,
                   targets=targets, base=base: relation_run(
                       curve, gen, challenge, targets, base,
                       expected_secret=secret, meter=meter))
        assert variants[name]["verified"], f"{name} failed exact recovery"
    for left, right in [("original", "transported"),
                        ("descendant_native", "pullback")]:
        a, b = variants[left], variants[right]
        assert [c["witness_count"] for c in a["cases"]] == [
            c["witness_count"] for c in b["cases"]], (left, right)
        assert a["first_full_rank_attempt"] == b["first_full_rank_attempt"]
        assert a["rank"] == b["rank"] == BASE_SIZE + 1

    cold_to_rank = {}
    for name, v in variants.items():
        limit = v["first_full_rank_attempt"]
        assert limit is not None
        shared = [ledger["field_setup"], ledger["subgroup_generator"],
                  challenge_cost]
        shared.extend((codomain_attempt_costs if name in
                       ("transported", "descendant_native") else
                       target_attempt_costs)[:limit])
        relation = v["costs"]
        specific = [relation["pair_table"], relation["scan_to_rank"],
                    relation["recovery_verify"]]
        if name == "original":
            specific.append(ledger["original_base"])
        elif name == "transported":
            specific.extend([ledger["kernel_line_search"],
                             ledger["selected_map_construction"],
                             ledger["original_base"], ledger["base_transport"],
                             transport_setup_cost])
        elif name == "descendant_native":
            specific.extend([ledger["kernel_line_search"],
                             ledger["selected_map_construction"],
                             ledger["native_base"], transport_setup_cost])
        else:
            specific.extend([ledger["kernel_line_search"],
                             ledger["selected_map_construction"],
                             ledger["native_base"], ledger["toy_inverse_table"],
                             ledger["base_pullback"]])
        cold_to_rank[name] = add_costs(*shared, *specific)

    raw = {
        "schema": "ecc2k130-degree7-paired-factor-base-pilot-v1",
        "status": "PASS",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "parameters": {"seed": SEED, "degree": N, "irreducible": hex(IRR),
                       "isogeny_degree": ELL, "source_order": ORDER,
                       "twist_order": TWIST_ORDER, "subgroup_order": R,
                       "cofactor": COFACTOR, "base_size": BASE_SIZE,
                       "attempts": ATTEMPTS},
        "source_sha256": {
            "protocol": digest(Path(__file__).with_name("PROTOCOL.md")),
            "runner": digest(Path(__file__)),
            "fastfield": digest(ROOT / "research/ecc2k130_relations/fastfield.py"),
            "relations": digest(ROOT / "research/ecc2k130_relations/relations.py"),
            "oriented_velu": digest(
                ROOT / "research/ecc2k130_oriented_transport_20260924/oriented_velu.py"),
        },
        "geometry": {"kernel_sampling_trials": trials,
                     "line_count": len(lines), "self_j_lines": sum(
                         v["codomain_b"] == 1 for v in lines.values()),
                     "selected_kernel_x": list(selected),
                     "selected_generator": encode_point(record["generator"]),
                     "selected_codomain_b": record["codomain_b"],
                     "all_lines": [
                         {"kernel_x": list(key),
                          "codomain_b": v["codomain_b"],
                          "first_sample_trial": v["trial"]}
                         for key, v in sorted(lines.items())]},
        "challenge": {"generator": encode_point(G), "generator_x_trials": g_trials,
                      "target": encode_point(Q), "secret_audit_only": secret},
        "bases": {"original": [encode_point(P) for P in B0],
                  "transported": [encode_point(P) for P in B_transport],
                  "descendant_native": [encode_point(P) for P in B1],
                  "pullback": [encode_point(P) for P in B_pullback],
                  "original_scan": B0_meta, "native_scan": B1_meta},
        "independent_lifts": lift,
        "phase_costs": ledger,
        "prefix_costs": {"challenge_setup": challenge_cost,
                         "target_generation": target_attempt_costs,
                         "transport_setup": transport_setup_cost,
                         "target_transport_audit": transport_attempt_costs,
                         "codomain_generation": codomain_attempt_costs},
        "cold_cost_to_verified_rank": cold_to_rank,
        "variants": variants,
        "comparison_controls": {
            "original_transport_hit_digest":
                comparison_digest(variants["original"]["cases"]),
            "native_pullback_hit_digest":
                comparison_digest(variants["descendant_native"]["cases"]),
            "all_full_lifts_verified": True,
            "all_scalars_recovered": True,
        },
        "limitations": [
            "toy field and m=2 pair table; no ECC2K-130 PDP solver",
            "degree 7 is ramified in CM, unlike split degree 263",
            "finite subgroup inverse table is an audit device, not scalable",
            "native counters are exclusive by phase but not calibrated to curve additions",
            "single fixed target stream; no generalizable speed claim",
            "per-target phi(T) is audit-only; codomain policies generate paired targets directly",
        ],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(raw, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": raw["status"], "line_count": len(lines),
                      "selected_codomain_b": hex(record["codomain_b"]),
                      "hits": {name: v["hits"] for name, v in variants.items()},
                      "first_full_rank_attempt": {
                          name: v["first_full_rank_attempt"]
                          for name, v in variants.items()},
                      "output": str(args.out)}, indent=2))


if __name__ == "__main__":
    main()

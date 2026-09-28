#!/usr/bin/env python3
"""Complete n19 rotated six-sum census for a preregistered normal generator."""
from __future__ import annotations

import argparse
import json
import resource
import signal
import struct
import sys
import tarfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
CORPUS = HERE.parent / "rotated_pdp_corpus_20260925"
PARENT = HERE.parent / "rotated_subspace_support_20260925"
sys.path.insert(0, str(CORPUS))
sys.path.insert(0, str(PARENT))
import gate  # noqa: E402
from producer import lifts, pjson, key, save  # noqa: E402 - #767 exact lift

REFERENCE_ARCHIVE = CORPUS / "evidence" / "raw.tar.gz"
Q = 130873
H = (385982, 301867)
BETA_REFERENCE_SUPPORT = 62389
CAP_SECONDS = 300
CAP_RSS = 512 * 1024 * 1024
Point = tuple[int, int] | None


def rss() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def point(value) -> Point:
    return None if value is None else tuple(value)


def reference():
    with tarfile.open(REFERENCE_ARCHIVE, "r:gz") as archive:
        def load(name: str):
            stream = archive.extractfile("raw/n19-m6/" + name)
            assert stream is not None
            return stream.read().decode("utf8")
        projected = {}
        for line in load("projected_histogram.jsonl").splitlines():
            row = json.loads(line)
            projected[point(row["point"])] = row["count"]
        targets = json.loads(load("targets.json"))
    assert len(projected) == BETA_REFERENCE_SUPPORT and len(targets) == 8
    assert sum(projected.values()) == 7 ** 6
    return projected, targets


def run_arm(beta: int, out: Path) -> None:
    started, cpu_started = time.perf_counter(), time.process_time()
    def expired(_signal, _frame):
        raise TimeoutError(f"beta {beta}: {CAP_SECONDS}s producer cap")
    signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, CAP_SECONDS)
    try:
        selection = json.loads((HERE / "selection.json").read_text())
        assert selection["status"] == "success"
        assert beta in [row["beta"] for row in selection["selected"]]
        reference_projected, frozen_targets = reference()
        field = gate.Field(19, [0, 1, 2, 5])
        field.rabin_prime_degree()
        assert gate.source_group_order(19) == 4 * Q and gate.is_prime_by_trial(Q)
        curve = gate.Curve(field)
        assert curve.on_curve(H) and curve.scalar(H, Q) is None
        assert curve.tau(H) == curve.scalar(H, 41811)
        conjugates = gate.normal_conjugates(field, beta)
        bases = gate.subspace_basis(conjugates, 6, 2)
        factors = [sorted(p for x in gate.all_x(basis) for p in lifts(curve, x))
                   for basis in bases]
        assert [len(factor) for factor in factors] == [7] * 6
        assert all(factors[i + 1] == sorted(curve.tau(p) for p in factors[i])
                   for i in range(5))
        columns = [{curve.scalar(p, 4) for p in factor} for factor in factors]
        column_signed = [len(column) for column in columns]
        column_quotient = [len({None if p is None else min(p, curve.neg(p))
                                for p in column}) for column in columns]
        selected = next(row for row in selection["selected"] if row["beta"] == beta)
        assert column_signed[0] == selected["projected_signed"]
        assert column_quotient[0] == selected["projected_sign_quotient"]
        setup_wall = time.perf_counter() - started
        setup_ops = {"field": dict(field.operations), "curve": dict(curve.operations)}

        full, witnesses = gate.full_histogram(curve, factors)
        assert sum(full.values()) == 7 ** 6
        histogram_wall = time.perf_counter() - started - setup_wall
        projected, projected_witness = {}, {}
        for s, count in full.items():
            r = curve.scalar(s, 4)
            projected[r] = projected.get(r, 0) + count
            projected_witness.setdefault(r, (s, witnesses[s]))
        assert sum(projected.values()) == 7 ** 6
        projection_wall = time.perf_counter() - started - setup_wall - histogram_wall

        four_h = curve.scalar(H, 4)
        assert four_h is not None and curve.scalar(four_h, Q) is None
        current = None
        counts = []
        for _ in range(Q):
            counts.append(projected.get(current, 0))
            current = curve.add(current, four_h)
        assert current is None and sum(counts) == 7 ** 6
        assert sum(c > 0 for c in counts) == len(projected)
        reference_set, current_set = set(reference_projected), set(projected)
        energy = sum(count * count for count in projected.values())
        reference_energy = sum(count * count for count in reference_projected.values())
        assert (energy - 7 ** 6) % 2 == 0
        overlap = len(reference_set & current_set)
        union = len(reference_set | current_set)
        torsion = [None, (0, 1), (1, 0), (1, 1)]
        frozen_rows = []
        for i, row in enumerate(frozen_targets):
            qpoint, r = point(row["Q"]), point(row["R"])
            assert curve.scalar(qpoint, 4) == r
            cosets = [curve.add(qpoint, t) for t in torsion]
            multiplicities = [full.get(p, 0) for p in cosets]
            assert sum(multiplicities) == projected.get(r, 0)
            frozen_rows.append({"index": i, "reference_class": row["class"],
                                "Q": row["Q"], "R": row["R"],
                                "projected_multiplicity": projected.get(r, 0),
                                "coset_multiplicities": multiplicities,
                                "coset_witness_indices": [list(witnesses[p]) if p in witnesses else None
                                                          for p in cosets]})
        comparison_wall = time.perf_counter() - started - setup_wall - histogram_wall - projection_wall
        save(out / "factors.json", [[pjson(p) for p in factor] for factor in factors])
        save(out / "fixed_targets.json", frozen_rows)
        (out / "target_counts.u32le").write_bytes(b"".join(struct.pack("<I", c) for c in counts))
        with (out / "full_histogram.jsonl").open("w") as stream:
            for p in sorted(full, key=key):
                stream.write(json.dumps({"point": pjson(p), "count": full[p],
                                         "witness_indices": list(witnesses[p])},
                                        sort_keys=True, separators=(",", ":")) + "\n")
        with (out / "projected_histogram.jsonl").open("w") as stream:
            for p in sorted(projected, key=key):
                s, witness = projected_witness[p]
                stream.write(json.dumps({"point": pjson(p), "count": projected[p],
                                         "full_sum": pjson(s),
                                         "witness_indices": list(witness)},
                                        sort_keys=True, separators=(",", ":")) + "\n")
        summary = {"beta": beta, "q": Q, "m": 6, "d": 2,
                   "factor_sizes": [len(factor) for factor in factors],
                   "projected_column_signed": column_signed,
                   "projected_column_sign_quotient": column_quotient,
                   "labelled_tuples": 7 ** 6,
                   "distinct_full_sums": len(full),
                   "distinct_projected_sums": len(projected),
                   "full_collisions": 7 ** 6 - len(full),
                   "projected_collisions": 7 ** 6 - len(projected),
                   "projected_energy": energy,
                   "projected_colliding_pairs": (energy - 7 ** 6) // 2,
                   "effective_support_floor_num": (7 ** 6) ** 2,
                   "effective_support_floor_den": energy,
                   "reference_projected_energy": reference_energy,
                   "reference_effective_support_floor_num": (7 ** 6) ** 2,
                   "reference_effective_support_floor_den": reference_energy,
                   "exact_misses": Q - len(projected),
                   "counting_minimum_misses": Q - 7 ** 6,
                   "reference_support": len(reference_set),
                   "support_intersection": overlap,
                   "support_union": union,
                   "candidate_only": len(current_set - reference_set),
                   "reference_only": len(reference_set - current_set),
                   "common_misses": Q - union,
                   "fixed_reference_targets_positive": sum(row["projected_multiplicity"] > 0 for row in frozen_rows),
                   "setup_wall_seconds": setup_wall,
                   "histogram_wall_seconds": histogram_wall,
                   "projection_wall_seconds": projection_wall,
                   "comparison_wall_seconds": comparison_wall,
                   "total_wall_seconds": time.perf_counter() - started,
                   "total_cpu_seconds": time.process_time() - cpu_started,
                   "peak_rss_bytes": rss(),
                   "setup_operations": setup_ops,
                   "total_operations": {"field": dict(field.operations),
                                        "curve": dict(curve.operations)}}
        save(out / "summary.json", summary)
        assert summary["total_wall_seconds"] <= CAP_SECONDS and rss() <= CAP_RSS
    except Exception as error:
        save(out / "producer_failure.json", {"beta": beta, "error": repr(error),
                                               "wall_seconds": time.perf_counter() - started,
                                               "peak_rss_bytes": rss()})
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--beta", type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    run_arm(args.beta, args.out)


if __name__ == "__main__":
    main()

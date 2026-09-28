"""Tiny-factor-base full-Frobenius-orbit MITM relation-stage control."""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

import run as base

HERE = Path(__file__).resolve().parent
PRIOR = HERE / "results" / "run_20260927"


def record(**kwargs):
    print(json.dumps(kwargs, sort_keys=True), flush=True)


def child(seed):
    resource.setrlimit(resource.RLIMIT_AS, (1 << 30, 1 << 30))
    started = time.perf_counter()
    cutoff = started + 67
    counts = collections.Counter()
    phase = ["setup"]
    setup = base.Setup(83, 2, seed, "ternary_inline", base.MODULUS)
    curve = setup.curve
    original_add = curve.add

    def counted_add(p, q):
        counts[phase[0]] += 1
        return original_add(p, q)

    curve.add = counted_add
    g = base.field_generator(curve)
    root = base.sqrt_mod(-7 % curve.r, curve.r)
    roots = [((-1 + sign * root) * pow(2, -1, curve.r)) % curve.r for sign in (1, -1)]
    lam, = [a for a in roots if curve.scale(g, a) == curve.frob(g)]
    assert pow(lam, 83, curve.r) == 1 and curve.order == 4 * curve.r
    planted, natural, used = base.make_inputs(setup, seed, g)
    assert len(planted) == 2
    inputs = planted + natural
    input_hashes = []
    for index, item in enumerate(inputs):
        saved = json.loads((PRIOR / f"fixed-s{seed}-i{index}-fes.json").read_text())
        prior = json.loads(saved["stdout"])
        assert prior["status"] == "complete"
        assert base.digest(item) == prior["input_sha256"]
        assert list(item["target"]) == prior["input"]["target"]
        input_hashes.append(prior["input_sha256"])
    record(stage="setup", seed=seed, normal_generator=setup.normal.beta,
           modulus=hex(base.MODULUS), generator=g, subgroup_order=curve.r,
           eigenvalue=lam, usable_lifts=len(used), input_sha256=input_hashes,
           group_additions=counts["setup"], seconds=time.perf_counter() - started)

    phase[0] = "orbit"
    reps = sorted({min(point, curve.neg(point)) for _, point in used})
    orbit = {}
    collisions = 0
    for column, rep in enumerate(reps):
        point, power = rep, 1
        for frob_phase in range(83):
            for sign, signed in ((1, point), (-1, curve.neg(point))):
                assert curve.valid(signed) and curve.scale(signed, curve.r) is None
                description = (column, power * sign % curve.r, frob_phase, sign)
                if signed in orbit:
                    collisions += 1
                else:
                    orbit[signed] = description
            point = curve.frob(point)
            power = power * lam % curve.r
        assert point == rep and power == 1
    points = sorted(orbit)
    orbit_hash = base.digest([[point, orbit[point]] for point in points])
    record(stage="orbit", seed=seed, representatives=reps, signed_orbit_size=len(points),
           collisions=collisions, orbit_sha256=orbit_hash,
           group_additions=counts["orbit"], seconds=time.perf_counter() - started)

    phase[0] = "pairs"
    pair_table = collections.defaultdict(list)
    pair_count = 0
    for i, left in enumerate(points):
        if time.perf_counter() > cutoff:
            raise TimeoutError(f"pair table expired at i={i}, pairs={pair_count}")
        for j in range(i, len(points)):
            result = curve.add(left, points[j])
            pair_count += 1
            if result is not None:
                pair_table[result].append((i, j))
    record(stage="pairs", seed=seed, attempts=pair_count, unique_sums=len(pair_table),
           group_additions=counts["pairs"], seconds=time.perf_counter() - started,
           peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024)

    natural_matrix = {}
    for index, item in enumerate(inputs):
        phase[0] = "search"
        target = item["target"]
        candidates = set()
        opposite_rejections = 0
        for third in points:
            if time.perf_counter() > cutoff:
                raise TimeoutError(f"target search expired at i={index}")
            needed = curve.add(target, curve.neg(third))
            for a, b in pair_table.get(needed, []):
                key = tuple(sorted((points[a], points[b], third)))
                if key in candidates:
                    continue
                triples = (points[a], points[b], third)
                if not base.c.proper(curve, triples):
                    opposite_rejections += 1
                    continue
                candidates.add(key)
        phase[0] = "verify"
        row_set = set()
        tautologies = 0
        for triple in sorted(candidates):
            assert all(curve.valid(p) and curve.scale(p, curve.r) is None for p in triple)
            assert base.c.group_sum(curve, triple) == target
            row = [0] * len(reps)
            for point in triple:
                col, coeff, _, _ = orbit[point]
                row[col] = (row[col] + coeff) % curve.r
            row = tuple(row)
            row_set.add(row)
            if target in orbit:
                col, coeff, _, _ = orbit[target]
                single = [0] * len(reps)
                single[col] = coeff
                tautologies += row == tuple(single)
            if item["kind"] == "natural":
                base.c.rank_add(natural_matrix, list(row), curve.r)
        record(stage="target", seed=seed, index=index, kind=item["kind"], target=target,
               input_sha256=input_hashes[index], scalar=item.get("scalar"),
               verified_triples=len(candidates), distinct_rows=len(row_set),
               rows=sorted(row_set), tautological_triples=tautologies,
               opposite_rejections=opposite_rejections, natural_rank=len(natural_matrix),
               group_additions_search=counts["search"],
               group_additions_verify=counts["verify"],
               seconds=time.perf_counter() - started)
    record(stage="final", seed=seed, status="complete", natural_rank=len(natural_matrix),
           representative_count=len(reps), natural_targets=len(natural),
           group_additions=dict(counts), seconds=time.perf_counter() - started,
           peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
           full_dlp_verified=False, matched_rho_operations=None, full_cost_S=None)


def parent(out):
    out.mkdir(parents=True, exist_ok=False)
    files = sorted((HERE / "source").rglob("*.py")) + [HERE / "run.py", HERE / "run_orbit.py",
            HERE / "PROTOCOL_ORBIT_RELATIONS.md"]
    (out / "manifest.json").write_text(json.dumps({
        "seeds": base.SEEDS, "input_results": str(PRIOR.relative_to(HERE)),
        "limits": {"address_space_bytes": 1 << 30, "wall_seconds": 70},
        "source_sha256": {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in files}}, sort_keys=True, indent=2) + "\n")
    for seed in base.SEEDS:
        cmd = [sys.executable, str(HERE / "run_orbit.py"), "--child", "--seed", str(seed)]
        base.run_child(out, f"orbit-s{seed}", cmd)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path)
    p.add_argument("--child", action="store_true")
    p.add_argument("--seed", type=int, choices=base.SEEDS)
    a = p.parse_args()
    if a.child:
        try:
            child(a.seed)
        except (TimeoutError, MemoryError) as error:
            record(stage="final", seed=a.seed,
                   status="timeout" if isinstance(error, TimeoutError) else "oom",
                   message=str(error))
    else:
        parent(a.out)

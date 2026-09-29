#!/usr/bin/env python3
"""Replay fixed public n=53 fixtures and verify their relations in Python.

This verifier uses polynomial convolution/long division, exponentiation for
inversion, and direct affine group addition. It does not import Rust arithmetic
or trust the producer's group-valid flag. It is restricted to the retained
public synthetic base and the twelve archived fixture labels.
"""

import argparse
import hashlib
import gzip
import itertools
import json
import os
from pathlib import Path
import subprocess
import time


REPO = Path(__file__).resolve().parents[4]
STUDY = REPO / "research/sat_factor_base_review_20260908"
ARCHIVE = STUDY / "autolab_orbit_extract_20260924"
BASE = ARCHIVE / "independent_replay_20260924_codex/base_header.jsonl.gz"
EXE = REPO / "target/release/examples/koblitz_s5_sat_instance"
SOURCE = REPO / "examples/koblitz_s5_sat_instance.rs"
EXPECTED_BASE_SHA256 = "d859319015ea405fd18aee41b51396ce4edcab64ef66265d8edcdeb5e040eb71"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Curve:
    def __init__(self, header):
        assert header["n"] == 53 and header["a"] == 0
        self.n = header["n"]
        self.modulus = (1 << self.n) | sum(1 << bit for bit in header["field_modulus_low_terms"])
        self.order = header["subgroup_order"]
        self.generator = tuple(header["generator"])

    def multiply(self, a, b):
        # Full carryless convolution, then long reduction, unlike the Rust
        # fast-word field kernel and the shift-and-reduce Python panel replay.
        product = 0
        for bit in range(b.bit_length()):
            if b >> bit & 1:
                product ^= a << bit
        while product.bit_length() > self.n:
            product ^= self.modulus << (product.bit_length() - self.n - 1)
        return product

    def power(self, a, exponent):
        result = 1
        while exponent:
            if exponent & 1:
                result = self.multiply(result, a)
            a = self.multiply(a, a)
            exponent >>= 1
        return result

    def inverse(self, a):
        assert a
        value = self.power(a, (1 << self.n) - 2)
        assert self.multiply(a, value) == 1
        return value

    def on_curve(self, point):
        if point is None:
            return True
        x, y = point
        if not (0 <= x < 1 << self.n and 0 <= y < 1 << self.n):
            return False
        return self.multiply(y, y) ^ self.multiply(x, y) == self.multiply(self.multiply(x, x), x) ^ 1

    def add(self, p, q):
        if p is None:
            return q
        if q is None:
            return p
        x, y = p
        z, w = q
        if x == z and (y ^ w) == x:
            return None
        if x == z:
            if x == 0:
                return None
            lam = x ^ self.multiply(y, self.inverse(x))
            out_x = self.multiply(lam, lam) ^ lam
            out_y = self.multiply(x, x) ^ self.multiply(lam ^ 1, out_x)
        else:
            lam = self.multiply(y ^ w, self.inverse(x ^ z))
            out_x = self.multiply(lam, lam) ^ lam ^ x ^ z
            out_y = self.multiply(lam, x ^ out_x) ^ out_x ^ y
        result = out_x, out_y
        assert self.on_curve(result)
        return result

    def scalar(self, p, k):
        acc = None
        for bit in bin(k)[2:]:
            acc = self.add(acc, acc)
            if bit == "1":
                acc = self.add(acc, p)
        return acc

    def semaev3(self, u, v, w):
        cross = self.multiply(u, v) ^ self.multiply(u, w) ^ self.multiply(v, w)
        return self.multiply(cross, cross) ^ self.multiply(self.multiply(u, v), w) ^ 1


def check_witness(curve, points_by_x, target, codes, intermediates):
    assert len(codes) == 4 and len(intermediates) == 2
    assert all(code in points_by_x for code in codes)
    u, v = intermediates
    assert curve.semaev3(codes[0], codes[1], u) == 0
    assert curve.semaev3(codes[2], codes[3], v) == 0
    assert curve.semaev3(u, v, target[0]) == 0
    for choice in itertools.product(*(points_by_x[x] for x in codes)):
        first = curve.add(choice[0], choice[1])
        second = curve.add(choice[2], choice[3])
        if (first is not None and second is not None and first[0] == u
                and second[0] == v and curve.add(first, second) == target):
            return [list(point) for point in choice]
    raise AssertionError("no sign choice matches both pinned pair sums and the public target")


def run(out, variant):
    out.mkdir(parents=True, exist_ok=False)
    original = json.loads((ARCHIVE / "results.json").read_text())
    claim = json.loads((ARCHIVE / "claim_relation_yield.json").read_text())
    archived = original["retained_base_panel"]["runs"]
    assert len(archived) == 12
    assert [(r["target_class"], r["seed"]) for r in archived] == (
        [("planted", k) for k in range(1, 5)] + [("natural", k) for k in range(1, 9)]
    )
    assert original["status"] == claim["status"] == "PENDING_INDEPENDENT_VALIDATION"
    assert claim["base_hash"] == EXPECTED_BASE_SHA256
    original_sha = sha256(ARCHIVE / "results.json")
    source_sha = sha256(SOURCE)
    exe_sha = sha256(EXE)
    assert source_sha != claim["executable_or_source_hash"]["source_sha256"], "candidate must differ from archive"
    with gzip.open(BASE, "rb") as stream:
        header_bytes = stream.readline()
        assert not stream.read(), "expected exactly one bundled header"
    base_path = out / "base_header.jsonl"
    base_path.write_bytes(header_bytes)
    header = json.loads(header_bytes)
    assert header["base_hash"] == EXPECTED_BASE_SHA256
    assert header["subgroup_order"] == 21044858204113
    assert header["field_modulus_low_terms"] == [0, 1, 2, 6]
    assert len(header["factor_base_point_coordinates"]) == 23320
    assert len(header["factor_base_representatives"]) == 220
    curve = Curve(header)
    assert curve.on_curve(curve.generator)
    assert curve.scalar(curve.generator, curve.order) is None

    points = [tuple(p) for p in header["factor_base_point_coordinates"]]
    assert len(set(points)) == len(points)
    by_x = {}
    for point in points:
        assert curve.on_curve(point)
        by_x.setdefault(point[0], []).append(point)
    assert len(by_x) == 11660 and all(len(group) == 2 for group in by_x.values())
    assert all(group[0][1] ^ group[1][1] == x for x, group in by_x.items())

    env = {key: val for key, val in os.environ.items() if not key.startswith("KIC_")}
    env.update({
        "KIC_ALGEBRA_ENCODING": "orbit_factorized",
        "KIC_ORBIT_LAZY_RELATIVE_SUPPORT": "1",
        "KIC_ORBIT_BRANCH_ORDER": "pair_then_pair",
        "KIC_ORBIT_REP_ENCODING": "one_hot",
        "KIC_FACTOR_BASE_JSONL": str(base_path),
        "KIC_TASK_ID": "TASK-IC-ORBIT-INDEPENDENT-REPLAY-20260924",
    })
    rows = []
    for prior in archived:
        kind, seed = prior["target_class"], prior["seed"]
        command = [str(EXE), "53", "0", "1", "10", kind, str(seed), "2000", "1", "internal"]
        started = time.perf_counter()
        process = subprocess.run(command, capture_output=True, text=True, env=env, timeout=30)
        elapsed = time.perf_counter() - started
        stdout_path = out / f"{kind}_{seed:02d}.stdout.jsonl"
        stderr_path = out / f"{kind}_{seed:02d}.stderr.txt"
        stdout_path.write_text(process.stdout)
        stderr_path.write_text(process.stderr)
        assert process.returncode == 0, (kind, seed, process.stderr[-1000:])
        assert len(process.stdout.splitlines()) == 1
        current = json.loads(process.stdout)
        assert current["factor_base_input_hash"] == EXPECTED_BASE_SHA256
        # The archived BLAKE3 covers its original whole JSONL container,
        # including appended observations. This replay uses a different
        # retained JSONL container with the same certified point-set hash.
        assert current["factor_base_input_blake3"] != claim["fixture_hash"]
        assert current["n"] == 53 and current["a"] == 0
        assert current["target_class"] == kind and current["seed"] == seed
        assert current["pair_table_entries"] == current["pair_selector_variables"] == 0
        assert current["compact_orbit_extraction"]["pair_table_entries"] == 0
        assert current["compact_orbit_extraction"]["edge_selectors"] == 0
        assert current["decomposition_verdict"] == "SAT"
        assert current["valid_x_tuples"] >= 1 and current["invalid_group_lifts"] == 0
        assert current["conflicts"] == 0
        if kind == "natural":
            label = current["published_scalar_validator_label"]
            assert isinstance(label, int) and 1 <= label < curve.order
            target = curve.scalar(curve.generator, label)
        else:
            assert current["published_scalar_validator_label"] is None
            planted = current["planted_point_indices_validator_only"]
            assert len(planted) == 4 and all(0 <= i < len(points) for i in planted)
            target = None
            for index in planted:
                target = curve.add(target, points[index])
        assert target is not None and curve.on_curve(target)
        new_relation = current["compact_orbit_extraction"]
        assert new_relation["enabled"] and new_relation["group_valid"]
        replay_lift = check_witness(curve, by_x, target, new_relation["x_codes"], new_relation["pinned_intermediates"])
        archived_lift = check_witness(curve, by_x, target, prior["x_codes"], prior["pinned_intermediates"])
        for point in set(map(tuple, replay_lift + archived_lift)):
            assert curve.scalar(point, curve.order) is None
        row = {
            "target_class": kind, "seed": seed, "published_scalar": current["published_scalar_validator_label"],
            "target": list(target), "independent_replay_lift": replay_lift,
            "independent_archived_lift": archived_lift,
            "replay_x_codes": new_relation["x_codes"], "archived_x_codes": prior["x_codes"],
            "replay_pinned_intermediates": new_relation["pinned_intermediates"],
            "archived_pinned_intermediates": prior["pinned_intermediates"],
            "replay_conflicts": current["conflicts"], "pair_table_entries": current["pair_table_entries"],
            "edge_selectors": new_relation["edge_selectors"], "process_wall_seconds": elapsed,
            "producer_stdout_sha256": hashlib.sha256(process.stdout.encode()).hexdigest(),
            "producer_stderr_sha256": hashlib.sha256(process.stderr.encode()).hexdigest(),
            "producer_stdout": stdout_path.name, "producer_stderr": stderr_path.name,
        }
        rows.append(row)
        print(json.dumps({k: row[k] for k in ("target_class", "seed", "process_wall_seconds", "pair_table_entries")}))

    assert sha256(SOURCE) == source_sha and sha256(EXE) == exe_sha
    assert sha256(ARCHIVE / "results.json") == original_sha
    base_path.unlink()
    report = {
        "schema_version": "2.0", "accepted": True,
        "classification": f"INDEPENDENT_PUBLIC_GROUP_REPLAY_OF_12_{variant.upper()}_COMPACT_ORBIT_RELATIONS",
        "scope": "Fixed public synthetic n=53 retained base; verifier recomputes curve arithmetic and S3 identities, including archived and newly replayed witnesses. No rho or end-to-end IC comparison.",
        "independent_verifier_sha256": sha256(Path(__file__)),
        "producer_source_sha256": source_sha, "producer_executable_sha256": exe_sha,
        "bundled_base_gzip_sha256": sha256(BASE),
        "archive_results_sha256": original_sha,
        "base_header_sha256": hashlib.sha256(header_bytes).hexdigest(),
        "archived_base_file_blake3": claim["fixture_hash"],
        "replay_base_file_blake3": current["factor_base_input_blake3"],
        "base_hash": EXPECTED_BASE_SHA256,
        "field_modulus": curve.modulus, "generator": list(curve.generator), "subgroup_order": curve.order,
        "base_points_checked_on_curve": len(points),
        "archived_and_replay_pair_intermediates_exact": True,
        "replayed_targets": len(rows), "independently_checked_archived_relations": len(rows),
        "natural_targets": sum(row["target_class"] == "natural" for row in rows),
        "all_producer_conflicts_zero": all(row["replay_conflicts"] == 0 for row in rows),
        "pair_table_entries": 0, "edge_selectors": 0,
        "rows": rows,
        "official_ledger_promotion": False,
        "limitation": "The replay binary regenerates public fixture scalars and planted labels. The bundled header has the same certified point-set hash but a different whole-file BLAKE3 from the archived file. This positive replay does not certify exhaustive negative coverage, full-rank collection, or end-to-end IC cost.",
    }
    (out / "validation.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: report[k] for k in ("accepted", "classification", "replayed_targets", "independently_checked_archived_relations", "natural_targets")}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True, help="new output directory")
    parser.add_argument("--variant", choices=("deterministic", "unsorted"), default="deterministic")
    args = parser.parse_args()
    run(args.out, args.variant)

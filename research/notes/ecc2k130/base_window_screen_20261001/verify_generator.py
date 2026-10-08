#!/usr/bin/env python3
"""Independent scan and signed-Frobenius replay for point-defined base windows.

The arithmetic below comes from the repository's pure-Python independent
replayer. This tool reads no Q, scalar fixture, rank row, or timing outcome.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import sys
import tarfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(
    0,
    str(ROOT / "research/sat_factor_base_review_20260908/autolab_orbit_extract_20260924"),
)
from independent_replay import Curve, Field, points_with_x  # noqa: E402

ARCHIVE = ROOT / (
    "research/notes/ecc2k130/disjoint_cold_v2_outcome_20261001/"
    "evidence_run_36803331080/raw/n41_L1024.tar.gz"
)
ARCHIVE_SHA256 = "261b5f3e7b209c2163d807dc78554b96b80873eda052405156b1178b29dce4a4"
ARCHIVE_BASE = "disjoint-cold-v2-n41_L1024/n41_L1024/b00_ic_a.base.jsonl"
ARCHIVE_MATERIALIZATION = "disjoint-cold-v2-n41_L1024/materialization.json"
FROZEN_COMPACT_SHA256 = "702a0a05709bc14bc10bafbb08edb2a2f5f86794e970d84c317da3bc65cdcf38"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def orbit_key(field: Field, x: int) -> int:
    values = []
    for _ in range(field.n):
        values.append(x)
        x = field.mul(x, x)
    assert x == values[0]
    return min(values)


def sqrt_mod(a: int, p: int) -> int:
    """Tonelli-Shanks, used to derive the Frobenius eigenvalue independently."""
    a %= p
    assert p > 2 and p % 2 and pow(a, (p - 1) // 2, p) == 1
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, x = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, power = 0, t
        while power != 1:
            power = power * power % p
            i += 1
            assert i < m
        b = pow(c, 1 << (m - i - 1), p)
        x, t, c, m = x * b % p, t * b * b % p, b * b % p, i
    assert x * x % p == a
    return x


def frobenius_lambda(curve: Curve, rep: tuple[int, int], r: int) -> int:
    root = sqrt_mod(-7, r)
    trace = -1 if curve.a == 0 else 1
    inv_two = pow(2, -1, r)
    candidates = {((trace + sign * root) * inv_two) % r for sign in (-1, 1)}
    f = curve.f
    image = (f.mul(rep[0], rep[0]), f.mul(rep[1], rep[1]))
    matched = [candidate for candidate in candidates if curve.mul(candidate, rep) == image]
    assert len(matched) == 1
    return matched[0]


def replay_scan(scan_path: Path, receipt: dict, curve: Curve) -> tuple[list[int], list[int], list[tuple[int, int]]]:
    lines = [json.loads(line) for line in scan_path.read_text().splitlines()]
    assert len(lines) == receipt["trace_decisions"]
    assert Counter(line["status"] for line in lines) == receipt["trace_status_counts"]
    by_x = defaultdict(list)
    for line in lines:
        by_x[line["raw_x"]].append(line)
    assert sorted(by_x) == list(range(1, receipt["raw_x_scanned"] + 1))
    assert receipt["raw_x_scanned"] <= receipt["raw_x_cap"]
    seen: set[int] = set()
    accepted: list[int] = []
    selected_keys: list[int] = []
    selected_reps: list[tuple[int, int]] = []
    start, end = receipt["accepted_start"], receipt["accepted_end"]
    for raw_x, records in sorted(by_x.items()):
        expected_lifts = set(points_with_x(curve, raw_x))
        if not expected_lifts:
            assert records == [{"raw_x": raw_x, "status": "no_lift"}]
            continue
        observed_lifts: set[tuple[int, int]] = set()
        for index, row in enumerate(records):
            assert row["lift_index"] == index
            lifted = tuple(row["lifted"])
            assert lifted[0] == raw_x and lifted in expected_lifts
            assert lifted not in observed_lifts and curve.on_curve(lifted)
            observed_lifts.add(lifted)
            projected = curve.mul(receipt["cofactor"], lifted)
            if projected is None:
                assert row["status"] == "projected_infinity"
                assert "projected" not in row
                continue
            assert tuple(row["projected"]) == projected
            if projected[0] <= 1:
                assert row["status"] == "small_projected_x"
                continue
            key = orbit_key(curve.f, projected[0])
            assert row["orbit_key"] == key
            if key in seen:
                assert row["status"] == "duplicate_orbit"
                continue
            assert row["status"] == "accepted"
            assert row["accepted_ordinal"] == len(accepted)
            selected = start <= len(accepted) < end
            assert row["selected"] is selected
            seen.add(key)
            accepted.append(key)
            if selected:
                selected_keys.append(key)
                selected_reps.append(projected)
        if raw_x != receipt["raw_x_scanned"] or receipt["status"] == "raw_x_cap":
            assert observed_lifts == expected_lifts
        else:
            assert observed_lifts <= expected_lifts
    assert accepted == receipt["accepted_orbit_keys"]
    assert selected_keys == receipt["selected_orbit_keys"]
    assert len(accepted) == receipt["accepted_orbits"]
    assert len(selected_keys) == receipt["selected_orbits"]
    if receipt["status"] == "complete":
        assert len(accepted) == end and len(selected_reps) == receipt["columns"]
        assert lines[-1]["status"] == "accepted"
        assert lines[-1]["accepted_ordinal"] == end - 1
    else:
        assert receipt["status"] == "raw_x_cap"
        assert receipt["raw_x_scanned"] == receipt["raw_x_cap"]
        assert len(accepted) < end
    return accepted, selected_keys, selected_reps


def verify_header(header_path: Path, receipt: dict, curve: Curve,
                  reps: list[tuple[int, int]]) -> dict:
    header_bytes = header_path.read_bytes()
    assert header_bytes.endswith(b"\n") and header_bytes.count(b"\n") == 1
    base = json.loads(header_bytes)
    n, a, r, columns = (base[key] for key in ("n", "a", "subgroup_order", "orbit_columns"))
    assert base["kind"] == "point_defined_factor_base"
    assert (n, a, columns) == (receipt["n"], receipt["a"], receipt["columns"])
    assert base["field_modulus_low_terms"] == [
        i for i in range(n) if curve.f.modulus & (1 << i)
    ]
    assert [tuple(rep) for rep in base["factor_base_representatives"]] == reps
    assert len(reps) == columns
    assert len(base["factor_base_point_coordinates"]) == len(base["factor_base_point_labels"]) == 2 * n * columns
    assert base["factor_base_points"] == 2 * n * columns
    assert len(base["base_hash"]) == 64
    assert all(c in "0123456789abcdef" for c in base["base_hash"])
    lam = frobenius_lambda(curve, reps[0], r)
    f = curve.f
    verified_points = 0
    for column, rep in enumerate(reps):
        assert curve.on_curve(rep) and curve.mul(r, rep) is None
        assert orbit_key(f, rep[0]) == receipt["selected_orbit_keys"][column]
        x, y = rep
        coefficient = 1
        for step in range(n):
            at = 2 * (column * n + step)
            assert base["factor_base_point_coordinates"][at] == [x, y]
            assert base["factor_base_point_coordinates"][at + 1] == [x, x ^ y]
            assert base["factor_base_point_labels"][at] == [column, coefficient]
            assert base["factor_base_point_labels"][at + 1] == [column, (-coefficient) % r]
            assert curve.on_curve((x, y)) and curve.on_curve((x, x ^ y))
            verified_points += 2
            x, y = f.mul(x, x), f.mul(y, y)
            coefficient = coefficient * lam % r
        assert (x, y) == rep and coefficient == 1
    return {
        "header_sha256": hashlib.sha256(header_bytes).hexdigest(),
        "base_hash": base["base_hash"],
        "frobenius_lambda": lam,
        "representatives_order_checked": columns,
        "representatives_subgroup_checked": columns,
        "orbit_members_checked": verified_points,
    }


def verify(header_path: Path | None, scan_path: Path, receipt_path: Path,
           archived_control: bool = False, protocol: bool = False) -> dict:
    receipt = json.loads(receipt_path.read_text())
    assert receipt["schema"] == "koblitz-base-window-generator-receipt-v1"
    n, a = receipt["n"], receipt["a"]
    assert n % 2 == 1 and a in (0, 1)
    if protocol:
        config = json.loads((HERE / "CONFIG.json").read_text())
        assert (n, a, receipt["columns"]) == (
            config["n"], config["a"], config["useful_orbit_columns"]
        )
        assert receipt["window"] in range(len(config["window_starts"]))
        assert receipt["accepted_start"] == config["window_starts"][receipt["window"]]
        assert receipt["raw_x_cap"] == config["raw_x_trial_cap"]
        assert receipt["accepted_end"] <= config["accepted_orbits_to_scan"]
    assert receipt["accepted_start"] == receipt["window"] * receipt["columns"]
    assert receipt["accepted_end"] == receipt["accepted_start"] + receipt["columns"]
    # The point law is independent of the Rust generator; n41 has b=1.
    modulus = receipt["field_modulus_low_terms"]
    if protocol:
        assert modulus == [0, 3] and receipt["cofactor"] == config["cofactor"]
        assert receipt["subgroup_order"] == config["subgroup_order"]
    field = Field(n, modulus)
    curve = Curve(field, a)
    assert isinstance(receipt["cofactor"], int) and receipt["cofactor"] > 0
    accepted, keys, reps = replay_scan(scan_path, receipt, curve)
    if receipt["status"] == "complete":
        assert header_path and header_path.is_file()
        base_checks = verify_header(header_path, receipt, curve, reps)
    else:
        assert header_path is None or not header_path.exists()
        base_checks = {}
    if archived_control:
        assert (n, a, receipt["columns"], receipt["window"]) == (41, 0, 255, 0)
        assert receipt["status"] == "complete"
        assert sha(ARCHIVE) == ARCHIVE_SHA256
        with tarfile.open(ARCHIVE, "r:gz") as archive:
            archived_bytes = archive.extractfile(ARCHIVE_BASE).read()
            materialization = json.loads(archive.extractfile(ARCHIVE_MATERIALIZATION).read())
        assert materialization["schema"] == "compact-frozen-source-materialization-v1"
        assert materialization["pinned_files"] == 20
        assert materialization["replayed_from_snapshot"][
            "examples/koblitz_orbit_dlp_s3_batch.rs"
        ]["frozen_sha256"] == FROZEN_COMPACT_SHA256
        assert header_path.read_bytes() == archived_bytes
        base_checks["archived_v2_header_sha256"] = hashlib.sha256(archived_bytes).hexdigest()
        base_checks["archived_v2_compact_source_sha256"] = FROZEN_COMPACT_SHA256
    return {
        "schema": "koblitz-base-window-independent-replay-v1",
        "status": "PASS",
        "generator_status": receipt["status"],
        "n": n,
        "a": a,
        "window": receipt["window"],
        "accepted_orbits_replayed": len(accepted),
        "selected_orbits_replayed": len(keys),
        "scan_decisions_replayed": receipt["trace_decisions"],
        "scan_sha256": sha(scan_path),
        "generator_receipt_sha256": sha(receipt_path),
        "archived_control": archived_control,
        **base_checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--header", type=Path)
    parser.add_argument("--scan", required=True, type=Path)
    parser.add_argument("--receipt", required=True, type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--archived-control", action="store_true")
    parser.add_argument("--protocol", action="store_true")
    args = parser.parse_args()
    result = verify(args.header, args.scan, args.receipt, args.archived_control, args.protocol)
    if args.out:
        args.out.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()

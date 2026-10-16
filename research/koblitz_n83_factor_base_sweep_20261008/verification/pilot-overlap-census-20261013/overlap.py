#!/usr/bin/env python3
"""Independently replay all primary N83 pilot archives and their orbit overlaps."""

from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import sys

import blake3


HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
REPO = STUDY.parents[1]
sys.path.insert(0, str(HERE.parent / "paired-capacity-screen-20261013"))
import screen  # noqa: E402
PANEL = STUDY / "pilot-01"
PILOT = HERE.parent / "cross-stream-pilot-import-source-20261013/census.json"
FIRST_AUDIT = HERE.parent / "cross-stream-first-gate-source-20261013/audit.json"
OUTPUT = HERE / "overlap.json"
INPUT_SHA256 = {
    "manifest.json": "5bd862987844f709e03beeab9fd265e3d2024011bcadcb82f5af9c72018df6a6",
    "replay.json": "fa5e545c27dd64c32a47ddae05a2a0b3a031545d700897f3e0e34f59abc0d72e",
    "upload-receipt.json": "0eae40ffa0e938470671cf796633c1e47b61ff676df0807e7a7e8592e56bb828",
    "pilot-census.json": "d1a490e2e73552d96a5ee327e34367c5bb32fbc6f6d591eb9c0e5f029118a285",
    "first-audit.json": "dfc96f11e4a399a42559847dd39e94b753f37911a54d806e0e54e830e7735478",
    "koblitz_large_prime.rs": "640598c7e43229812c21d8fb84948104f1ef7e650a9d91f5f7ece764d95e98a5",
    "koblitz_n83_mixed_orbit_worker.rs": "5f310f7495989a94004a6545bed5682d54d3bc4b4f41db543ed24408a1fa1b7b",
    "paired-capacity-screen.py": "116617e9de7552e48b548e632e80121508515e7079878305de0b975c8a924c5d",
}
N = 83
SIZE = 1 << N
LOW_POLY = (1 << 45) | 7
POLY = SIZE | LOW_POLY
ORDER = 2417851639230796216685689
LAMBDA = 254512724090651164922414
CURVE = "icv1-f2m83-tm6151469093347-debefd74"
GENERATOR = ["5402229379304697578719232", "3600185575789834070382019"]
SQUARE_BYTE = tuple(sum(((byte >> bit) & 1) << (2 * bit)
                        for bit in range(8)) for byte in range(256))


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(value: dict) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2) + "\n").encode()


def square(value: int) -> int:
    require(0 <= value < SIZE, "field element outside degree 83")
    result = 0
    for byte in range(11):
        result ^= SQUARE_BYTE[(value >> (8 * byte)) & 255] << (16 * byte)
    while result.bit_length() > N:
        result ^= POLY << (result.bit_length() - N - 1)
    return result


def multiply(a: int, b: int) -> int:
    require(0 <= a < SIZE and 0 <= b < SIZE,
            "field multiplier outside degree 83")
    result = 0
    for _ in range(N):
        if b & 1:
            result ^= a
        b >>= 1
        carry = a >> (N - 1)
        a = (a << 1) & (SIZE - 1)
        if carry:
            a ^= LOW_POLY
    return result


def on_curve(x: int, y: int) -> bool:
    return square(y) ^ multiply(x, y) == multiply(square(x), x) ^ 1


def checked_inputs() -> tuple[dict, dict, dict, dict, dict]:
    paths = {
        "manifest.json": PANEL / "manifest.json",
        "replay.json": PANEL / "replay.json",
        "upload-receipt.json": PANEL / "upload-receipt.json",
        "pilot-census.json": PILOT,
        "first-audit.json": FIRST_AUDIT,
        "koblitz_large_prime.rs": REPO / "src/cryptanalysis/koblitz_large_prime.rs",
        "koblitz_n83_mixed_orbit_worker.rs":
            REPO / "examples/koblitz_n83_mixed_orbit_worker.rs",
        "paired-capacity-screen.py": Path(screen.__file__),
    }
    for label, path in paths.items():
        require(sha256(path) == INPUT_SHA256[label],
                f"frozen input changed: {label}")
    manifest, replay, upload, pilot, audit = (
        json.loads(paths[label].read_text()) for label in (
            "manifest.json", "replay.json", "upload-receipt.json",
            "pilot-census.json", "first-audit.json")
    )
    require(manifest["schema"] == "n83.factor-base-panel/v1"
            and manifest["completed_base_count"] == 54
            and len(manifest["bases"]) == 54
            and replay["status"] == upload["status"] == "PASS"
            and len(replay["checks"]) == len(upload["objects"]) == 54
            and pilot["primary_rows"] == 27
            and pilot["ordered_cross_stream_pairs"] == 648
            and audit["status"] == "PASS_receipts_with_censored_search",
            "pinned panel evidence changed")
    return manifest, replay, upload, pilot, audit


def load_archive(row: dict, objects_root: Path) -> dict:
    digest = row["compressed_blake3"]
    expected_name = f"{digest}.jsonl.gz"
    require(row["object"] == f"objects/{expected_name}"
            and row["points"] == row["columns"] * 166,
            "manifest object shape mismatch")
    path = objects_root / expected_name
    compressed = path.read_bytes()
    require(len(compressed) == row["bytes"]
            and blake3.blake3(compressed).hexdigest() == digest,
            f"compressed archive mismatch: {digest}")
    stream_hash = blake3.blake3()
    orbit_keys: set[int] = set()
    x_domain: set[int] = set()
    with gzip.open(path, "rb") as stream:
        header_line = stream.readline()
        stream_hash.update(header_line)
        header = json.loads(header_line)
        require(header["schema"] == "n83.factor-base/v1"
                and header["study"] == STUDY.name
                and header["curve_slug"] == CURVE
                and header["n"] == N
                and header["curve_a"] == 0
                and header["curve_b"] == "1"
                and header["modulus_low_terms"] == [0, 1, 2, 45]
                and header["subgroup_order"] == str(ORDER)
                and header["frobenius_lambda"] == str(LAMBDA)
                and header["cofactor"] == "4"
                and header["generator"] == GENERATOR
                and header["closure"] == "signed_frobenius"
                and header["policy"] == row["policy"]
                and header["seed"] == row["seed"]
                and header["orbit_columns"] == row["columns"]
                and header["point_count"] == row["points"],
                f"archive header mismatch: {digest}")
        reps = header["representatives"]
        require(len(reps) == row["columns"]
                and len(header["source_points_before_cofactor_clearing"])
                    == row["columns"]
                and len(header["accepted_candidate_indices"])
                    == row["columns"],
                f"archive column metadata mismatch: {digest}")
        for column, point in enumerate(reps):
            x, y = map(int, point)
            require(0 <= x < SIZE and 0 <= y < SIZE and on_curve(x, y),
                    f"representative off curve: {digest}:{column}")
            start = (x, y)
            phase_x: set[int] = set()
            coefficient = 1
            for phase in range(N):
                require(x not in phase_x and x not in x_domain,
                        f"short or duplicate signed orbit: {digest}:{column}")
                phase_x.add(x)
                x_domain.add(x)
                for negative in (False, True):
                    line = stream.readline()
                    require(bool(line), f"truncated point rows: {digest}")
                    stream_hash.update(line)
                    archived = json.loads(line)
                    expected_point = [str(x), str(y ^ (x if negative else 0))]
                    expected_coefficient = ORDER - coefficient if negative else coefficient
                    require(archived["point"] == expected_point
                            and archived["column"] == column
                            and archived["phase"] == phase
                            and archived["negative"] is negative
                            and archived["coefficient"] == str(expected_coefficient),
                            f"point or label differs: {digest}:{column}:{phase}")
                x, y = square(x), square(y)
                coefficient = coefficient * LAMBDA % ORDER
            require((x, y) == start and len(phase_x) == N,
                    f"signed orbit did not close: {digest}:{column}")
            key = min(phase_x)
            require(key not in orbit_keys,
                    f"duplicate signed orbit key: {digest}:{column}")
            orbit_keys.add(key)
        require(stream.readline() == b"", f"trailing archive rows: {digest}")
    require(stream_hash.hexdigest() == row["plain_blake3"]
            and len(orbit_keys) == row["columns"]
            and len(x_domain) == row["columns"] * N,
            f"archive plain hash or orbit count mismatch: {digest}")
    return {"orbit_keys": orbit_keys, "x_domain": x_domain,
            "plain_blake3": stream_hash.hexdigest()}


def build(objects_root: Path) -> dict:
    manifest, replay, upload, pilot, audit = checked_inputs()
    replay_map = {row["object"]: row for row in replay["checks"]}
    upload_map = {row["compressed_blake3"]: row for row in upload["objects"]}
    pilot_map = {(row["policy"], row["seed"], row["columns"]): row
                 for row in pilot["primary_objects"]}
    rows = sorted((row for row in manifest["bases"] if row["a"] == 0),
                  key=lambda row: (row["policy"], row["seed"], row["columns"]))
    require(len(rows) == len(pilot_map) == 27,
            "primary pilot population changed")
    archives = {}
    for row in rows:
        identity = (row["policy"], row["seed"], row["columns"])
        replay_row = replay_map[row["object"]]
        upload_row = upload_map[row["compressed_blake3"]]
        require(pilot_map[identity]["compressed_blake3"]
                == row["compressed_blake3"]
                and pilot_map[identity]["s3_uri"] == row["s3_uri"]
                and replay_row["status"] == "PASS"
                and replay_row["points_checked"] == row["points"]
                and replay_row["representatives_checked"] == row["columns"]
                and upload_row["s3_uri"] == row["s3_uri"]
                and upload_row["downloaded_hash_matches"] is True,
                "pilot replay or S3 round-trip receipt mismatch")
        archives[row["compressed_blake3"]] = load_archive(row, objects_root)
    pairs = []
    overlap_histogram: Counter[int] = Counter()
    for base, pool in itertools.permutations(rows, 2):
        if (base["policy"], base["seed"]) == (pool["policy"], pool["seed"]):
            continue
        base_id, pool_id = base["compressed_blake3"], pool["compressed_blake3"]
        shared = len(archives[base_id]["orbit_keys"]
                     & archives[pool_id]["orbit_keys"])
        shared_x = len(archives[base_id]["x_domain"]
                       & archives[pool_id]["x_domain"])
        require(shared_x == N * shared,
                f"orbit-key and exact x-domain overlap differ: {base_id}:{pool_id}")
        residual = pool["columns"] - shared
        require(0 <= shared <= min(base["columns"], pool["columns"])
                and 0 <= residual <= pool["columns"],
                "pair overlap outside set bounds")
        overlap_histogram[shared] += 1
        model_shapes = []
        if residual >= 2:
            for summands, residual_summands in screen.SHAPES:
                model = screen.requirements(
                    base["columns"], base["columns"] + residual,
                    summands, residual_summands)
                require(model["current_worker_cap_pass"],
                        "exact pilot model exceeds frozen worker caps")
                model_shapes.append({
                    "summands": summands,
                    "residual_summands": residual_summands,
                    "first_pattern_total_variables": model["total_variables"],
                    "first_pattern_domain_clauses": model["domain_clauses"],
                    "first_pattern_xor_rows": model["xor_rows"],
                })
        pairs.append({
            "base_object_blake3": base_id,
            "pool_object_blake3": pool_id,
            "base_policy": base["policy"],
            "pool_policy": pool["policy"],
            "base_seed": base["seed"],
            "pool_seed": pool["seed"],
            "base_orbits": base["columns"],
            "pool_orbits": pool["columns"],
            "shared_signed_orbits": shared,
            "residual_signed_orbits": residual,
            "shared_x_values": shared_x,
            "passes_current_min_two_residual_screen": residual >= 2,
            "source_derived_exact_model_shapes": model_shapes,
        })
    require(len(pairs) == 648 and sum(overlap_histogram.values()) == 648,
            "ordered cross-stream pilot population changed")
    first_matches = [row for row in pairs
                     if row["base_object_blake3"] == audit["base_object_blake3"]
                     and row["pool_object_blake3"] == audit["pool_object_blake3"]]
    require(len(first_matches) == 1
            and first_matches[0]["shared_signed_orbits"]
                == audit["shared_signed_orbits"] == 0
            and first_matches[0]["residual_signed_orbits"]
                == audit["residual_signed_orbits"] == 600,
            "retained first-pair overlap does not replay")
    policy_groups = {}
    for (base_policy, pool_policy), group in itertools.groupby(
            sorted(pairs, key=lambda row: (row["base_policy"], row["pool_policy"])),
            key=lambda row: (row["base_policy"], row["pool_policy"])):
        items = list(group)
        policy_groups[f"{base_policy}->{pool_policy}"] = {
            "ordered_pairs": len(items),
            "disjoint": sum(row["shared_signed_orbits"] == 0 for row in items),
            "empty_residual": sum(row["residual_signed_orbits"] == 0
                                  for row in items),
            "singleton_residual": sum(row["residual_signed_orbits"] == 1
                                      for row in items),
            "minimum_shared_orbits": min(row["shared_signed_orbits"]
                                          for row in items),
            "maximum_shared_orbits": max(row["shared_signed_orbits"]
                                          for row in items),
        }
    shape_summaries = []
    for summands, residual_summands in screen.SHAPES:
        models = [model for pair in pairs
                  for model in pair["source_derived_exact_model_shapes"]
                  if (model["summands"], model["residual_summands"])
                     == (summands, residual_summands)]
        require(len(models) == 558,
                "exact model-shape population differs from minimum-size screen")
        shape_summaries.append({
            "summands": summands,
            "residual_summands": residual_summands,
            "pair_count": len(models),
            "minimum_first_pattern_total_variables": min(
                model["first_pattern_total_variables"] for model in models),
            "maximum_first_pattern_total_variables": max(
                model["first_pattern_total_variables"] for model in models),
            "minimum_first_pattern_domain_clauses": min(
                model["first_pattern_domain_clauses"] for model in models),
            "maximum_first_pattern_domain_clauses": max(
                model["first_pattern_domain_clauses"] for model in models),
        })
    return {
        "schema": "n83.primary-pilot-exact-overlap-census/v1",
        "status": "PASS_archive_backed_exact_orbit_overlap",
        "input_sha256": INPUT_SHA256,
        "curve_slug": CURVE,
        "primary_archive_count": len(rows),
        "ordered_cross_stream_pair_count": len(pairs),
        "overlap_histogram": {str(k): v for k, v in sorted(overlap_histogram.items())},
        "disjoint_pair_count": overlap_histogram[0],
        "fully_absorbed_pool_count": sum(row["residual_signed_orbits"] == 0
                                          for row in pairs),
        "singleton_residual_count": sum(row["residual_signed_orbits"] == 1
                                        for row in pairs),
        "passes_current_min_two_residual_screen_count": sum(
            row["passes_current_min_two_residual_screen"] for row in pairs),
        "policy_pair_summaries": policy_groups,
        "source_derived_exact_model_shape_summaries": shape_summaries,
        "first_retained_pair": first_matches[0],
        "archive_objects": [{
            "policy": row["policy"], "seed": row["seed"],
            "orbit_columns": row["columns"],
            "compressed_blake3": row["compressed_blake3"],
            "plain_blake3": row["plain_blake3"],
            "s3_uri": row["s3_uri"],
            "independent_full_point_row_replay": True,
        } for row in rows],
        "ordered_pairs": pairs,
        "new_sat_search_executed": False,
        "selected_best_total_runtime": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    parser.add_argument("--objects-root", type=Path, required=True,
                        help="local directory of the manifest-pinned .jsonl.gz objects")
    args = parser.parse_args()
    rendered = canonical(build(args.objects_root))
    if args.write:
        OUTPUT.write_bytes(rendered)
        print(f"wrote {OUTPUT}")
    else:
        require(OUTPUT.read_bytes() == rendered,
                "saved exact-overlap census changed")
        print("PASS: pinned full-archive primary pilot overlap census")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Replay every public-Q selection and rejection from a frozen input manifest."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(
    0,
    str(ROOT / "research/sat_factor_base_review_20260908/autolab_orbit_extract_20260924"),
)
from independent_replay import Curve, Field  # noqa: E402

FROZEN = HERE / "FROZEN.json"
R = 549756390943
G = (2056947637384, 1635505394702)
N = 41


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha(path: Path) -> str:
    return digest(path.read_bytes())


def lines(path: Path) -> list:
    return [json.loads(line) for line in path.read_bytes().splitlines() if line.strip()]


def orbit_key(x: int, field: Field) -> int:
    smallest = x
    current = x
    for _ in range(N):
        smallest = min(smallest, current)
        current = field.mul(current, current)
    assert current == x
    return smallest


def prior_inventory(curve: Curve, field: Field, commit: str) -> tuple[list[dict], str, set[int], int]:
    names = [name for name in subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", commit], cwd=ROOT, text=True
    ).splitlines() if name.endswith(".points.jsonl")]
    assert all(not name.startswith(str(HERE.relative_to(ROOT)) + "/") for name in names)
    assert names == sorted(names) and names
    records, seen = [], set()
    n41_rows = 0
    for name in names:
        path = ROOT / name
        match = re.match(r"n(\d+)_", path.name)
        assert match, name
        degree = int(match.group(1))
        raw = subprocess.check_output(["git", "show", f"{commit}:{name}"], cwd=ROOT)
        points = [json.loads(line) for line in raw.splitlines() if line.strip()]
        assert points, name
        for point in points:
            assert (isinstance(point, list) and len(point) == 2 and
                    all(isinstance(word, int) and 0 <= word < (1 << degree)
                        for word in point))
            if degree == N:
                assert curve.on_curve(tuple(point))
                seen.add(orbit_key(point[0], field))
                n41_rows += 1
        records.append({"path": name, "sha256": digest(raw), "rows": len(points)})
    inventory_bytes = (json.dumps(
        records, sort_keys=True, separators=(",", ":")
    ) + "\n").encode("utf-8")
    return records, digest(inventory_bytes), seen, n41_rows


def verify() -> dict:
    frozen = json.loads(FROZEN.read_text())
    config = json.loads((HERE / "CONFIG.json").read_text())
    assert frozen["schema"] == "ecc2k130-base-window-input-freeze-v1"
    assert frozen["curve_slug"] == config["curve_slug"]
    assert (frozen["n"], frozen["a"], frozen["subgroup_order"],
            frozen["cofactor"]) == (N, 0, R, 4)
    assert frozen["field_modulus_low_terms"] == [0, 3]
    assert frozen["generator"] == list(G)
    assert frozen["protocol_config_sha256"] == sha(HERE / "CONFIG.json")
    assert frozen["protocol_sha256"] == sha(HERE / "PROTOCOL.md")
    assert frozen["prepare_sha256"] == sha(HERE / "prepare_inputs.py")
    source = frozen["source_lock"]
    commit = source["merge_commit"]
    assert len(commit) == 40 and all(c in "0123456789abcdef" for c in commit)
    subprocess.run(
        ["git", "merge-base", "--is-ancestor", commit, "HEAD"], cwd=ROOT, check=True
    )
    assert len(source["sha256"]) >= 10
    for name, expected in source["sha256"].items():
        assert sha(ROOT / name) == expected, name
        committed = subprocess.check_output(["git", "show", f"{commit}:{name}"], cwd=ROOT)
        assert digest(committed) == expected, name
    field = Field(N, [0, 3])
    curve = Curve(field, 0)
    assert curve.on_curve(G) and curve.mul(R, G) is None
    prior_tree = frozen["prior_point_tree_commit"]
    assert len(prior_tree) == 40 and all(c in "0123456789abcdef" for c in prior_tree)
    subprocess.run(
        ["git", "merge-base", "--is-ancestor", commit, prior_tree], cwd=ROOT, check=True
    )
    subprocess.run(
        ["git", "merge-base", "--is-ancestor", prior_tree, "HEAD"], cwd=ROOT, check=True
    )
    records, inventory_digest, prior, n41_rows = prior_inventory(curve, field, prior_tree)
    assert records == frozen["prior_point_inventory"]
    assert inventory_digest == frozen["prior_inventory_digest"]
    assert n41_rows == frozen["prior_n41_rows"]
    assert len(prior) == frozen["prior_n41_orbits"]
    assert frozen["candidate_hash_domain"] == config["target_domain"]
    assert len(frozen["blocks"]) == config["blocks"] == 5
    seen: set[int] = set()
    block_checks = []
    rejection_totals = {"zero_scalar": 0, "prior_orbit": 0, "new_orbit": 0}
    for block, item in enumerate(frozen["blocks"]):
        assert item["block"] == block
        assert item["rho_seed"] == config["rho_seed_base"] + block
        assert item["corpus"] == f"base-window-screen-n41-L1024-b{block:02d}-20261001-v1"
        prefix = f"fixtures/n41_L1024_b{block:02d}"
        assert item["fixture_file"] == prefix + ".fixture.jsonl"
        assert item["points_file"] == prefix + ".points.jsonl"
        fixture_path, point_path = HERE / item["fixture_file"], HERE / item["points_file"]
        assert sha(fixture_path) == item["fixture_sha256"]
        assert sha(point_path) == item["points_sha256"]
        fixtures, points = lines(fixture_path), lines(point_path)
        assert len(fixtures) == len(points) == config["public_targets_per_block"]
        attempts = item["candidate_attempts"]
        assert config["public_targets_per_block"] <= attempts <= config["target_candidate_cap_per_block"]
        rejected = {"zero_scalar": 0, "prior_orbit": 0, "new_orbit": 0}
        accepted = 0
        for index in range(attempts):
            message = f"{config['target_domain']}|{block}|{index}".encode("utf-8")
            hash_bytes = hashlib.sha256(message).digest()
            scalar = int.from_bytes(hash_bytes, "big") % R
            if scalar == 0:
                rejected["zero_scalar"] += 1
                continue
            q = curve.mul(scalar, G)
            assert q is not None and curve.on_curve(q) and curve.mul(R, q) is None
            key = orbit_key(q[0], field)
            if key in prior:
                rejected["prior_orbit"] += 1
                continue
            if key in seen:
                rejected["new_orbit"] += 1
                continue
            seen.add(key)
            assert accepted < config["public_targets_per_block"]
            assert fixtures[accepted] == {
                "kind": "base_window_public_fixture",
                "block": block,
                "fixture_index": accepted,
                "candidate_index": index,
                "candidate_sha256": hash_bytes.hex(),
                "published_fixture_scalar": scalar,
                "published_q": list(q),
                "orbit_key": key,
            }
            assert points[accepted] == list(q)
            accepted += 1
        assert accepted == config["public_targets_per_block"]
        assert fixtures[-1]["candidate_index"] == attempts - 1
        assert rejected == item["candidate_rejections"]
        for reason in rejection_totals:
            rejection_totals[reason] += rejected[reason]
        block_checks.append({
            "block": block,
            "candidate_attempts": attempts,
            "rejections": rejected,
            "points_sha256": item["points_sha256"],
            "fixture_sha256": item["fixture_sha256"],
        })
    assert len(seen) == frozen["new_n41_orbits"] == 5120
    return {
        "schema": "ecc2k130-base-window-input-replay-v1",
        "status": "PASS",
        "frozen_sha256": sha(FROZEN),
        "prior_inventory_digest": inventory_digest,
        "prior_n41_rows": n41_rows,
        "prior_n41_orbits": len(prior),
        "new_n41_orbits": len(seen),
        "point_equations_verified": 5120,
        "candidate_rejections": rejection_totals,
        "blocks": block_checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite an input replay receipt"
    result = verify()
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "blocks"}, sort_keys=True))


if __name__ == "__main__":
    main()

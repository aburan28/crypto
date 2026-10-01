#!/usr/bin/env python3
"""Freeze new n41 public-Q blocks after the base-window source lock merges.

Only this preparer handles known-answer scalars. The timed runner opens the
point-only JSONL files and never imports this module or reads fixture files.
"""
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

OLD_FREEZE = ROOT / "research/notes/ecc2k130/disjoint_cold_v2_20261001/FROZEN.json"
OLD_FREEZE_SHA256 = "da958a3f1117dd1b88703fed2055c5c9f64255516e1e920f35d193bf53c6a88d"
SOURCE_PATHS = (
    ".github/workflows/ecc2k130-base-window-source-lock.yml",
    "examples/koblitz_base_window.rs",
    "research/notes/ecc2k130/base_window_screen_20261001/CONFIG.json",
    "research/notes/ecc2k130/base_window_screen_20261001/PROTOCOL.md",
    "research/notes/ecc2k130/base_window_screen_20261001/SOURCE_LOCK.md",
    "research/notes/ecc2k130/base_window_screen_20261001/build_frozen.py",
    "research/notes/ecc2k130/base_window_screen_20261001/prepare_inputs.py",
    "research/notes/ecc2k130/base_window_screen_20261001/verify_inputs.py",
    "research/notes/ecc2k130/base_window_screen_20261001/run_screen.py",
    "research/notes/ecc2k130/base_window_screen_20261001/verify_screen.py",
    "research/notes/ecc2k130/base_window_screen_20261001/analyze_screen.py",
    "research/notes/ecc2k130/base_window_screen_20261001/test_source_lock.py",
    "research/notes/ecc2k130/base_window_screen_20261001/verify_generator.py",
    "research/notes/ecc2k130/base_window_screen_20261001/test_generator.py",
    "research/notes/ecc2k130/compact_ir_ledger_20260930/run_panel.py",
    "research/notes/ecc2k130/compact_orbit_rank_evidence_20260929/verify_rank.py",
    "research/notes/ecc2k130/disjoint_cold_v2_20261001/check_pairing.py",
    "research/sat_factor_base_review_20260908/autolab_orbit_extract_20260924/independent_replay.py",
    "tools/isolated_bench.py",
)
FIELD_DEGREE = 41
CURVE_A = 0
FIELD_LOW_TERMS = [0, 3]
GENERATOR = (2056947637384, 1635505394702)
ORDER = 549756390943
COFACTOR = 4


def sha_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha(path: Path) -> str:
    return sha_bytes(path.read_bytes())


def canon_bytes(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")


def config() -> dict:
    result = json.loads((HERE / "CONFIG.json").read_text())
    assert result["schema"] == "ecc2k130-base-window-screen-protocol-v1"
    assert (result["n"], result["a"], result["subgroup_order"], result["cofactor"]) == (
        FIELD_DEGREE, CURVE_A, ORDER, COFACTOR
    )
    assert (result["blocks"], result["public_targets_per_block"]) == (5, 1024)
    assert result["target_candidate_cap_per_block"] == 20000
    return result


def fixed_curve() -> tuple[Curve, Field]:
    assert sha(OLD_FREEZE) == OLD_FREEZE_SHA256
    old = json.loads(OLD_FREEZE.read_text())
    old_spec = old["specs"]["n41_L1024"]
    assert (old_spec["n"], old_spec["a"], old_spec["subgroup_order"]) == (
        FIELD_DEGREE, CURVE_A, ORDER
    )
    assert old_spec["field_modulus_low_terms"] == FIELD_LOW_TERMS
    assert old_spec["generator"] == list(GENERATOR)
    field = Field(FIELD_DEGREE, FIELD_LOW_TERMS)
    curve = Curve(field, CURVE_A)
    assert curve.on_curve(GENERATOR) and curve.mul(ORDER, GENERATOR) is None
    return curve, field


def orbit_key(x: int, field: Field) -> int:
    values = []
    for _ in range(FIELD_DEGREE):
        values.append(x)
        x = field.mul(x, x)
    assert x == values[0]
    return min(values)


def candidate(domain: str, block: int, index: int) -> tuple[int, str]:
    assert block >= 0 and index >= 0
    digest = hashlib.sha256(f"{domain}|{block}|{index}".encode("utf-8")).digest()
    return int.from_bytes(digest, "big") % ORDER, digest.hex()


def inventory() -> tuple[list[dict], str, set[int], int, str]:
    """Hash every point corpus in the preparation HEAD's immutable Git tree."""
    tree_commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    names = [name for name in subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", tree_commit], cwd=ROOT, text=True
    ).splitlines() if name.endswith(".points.jsonl")]
    assert names == sorted(names) and names
    assert all(not name.startswith(str(HERE.relative_to(ROOT)) + "/") for name in names)
    curve, field = fixed_curve()
    rows = []
    n41_orbits: set[int] = set()
    n41_points = 0
    for name in names:
        path = ROOT / name
        assert path.is_file(), f"tracked point corpus is absent: {name}"
        match = re.match(r"n(\d+)_", path.name)
        assert match, f"cannot classify point corpus: {name}"
        degree = int(match.group(1))
        raw = subprocess.check_output(["git", "show", f"{tree_commit}:{name}"], cwd=ROOT)
        assert raw == path.read_bytes(), f"working point corpus differs from HEAD: {name}"
        lines = [json.loads(line) for line in raw.splitlines() if line.strip()]
        assert lines, name
        for q in lines:
            assert isinstance(q, list) and len(q) == 2
            assert all(isinstance(word, int) and 0 <= word < (1 << degree) for word in q)
            if degree == FIELD_DEGREE:
                assert curve.on_curve(tuple(q)), name
                n41_orbits.add(orbit_key(q[0], field))
                n41_points += 1
        rows.append({"path": name, "sha256": sha_bytes(raw), "rows": len(lines)})
    return rows, sha_bytes(canon_bytes(rows)), n41_orbits, n41_points, tree_commit


def source_lock(commit: str) -> dict:
    assert len(commit) == 40 and all(c in "0123456789abcdef" for c in commit)
    subprocess.run(
        ["git", "merge-base", "--is-ancestor", commit, "HEAD"],
        cwd=ROOT, check=True,
    )
    paths = {name: sha(ROOT / name) for name in SOURCE_PATHS}
    for name, expected in paths.items():
        committed = subprocess.check_output(
            ["git", "show", f"{commit}:{name}"], cwd=ROOT
        )
        assert sha_bytes(committed) == expected, name
    return {"merge_commit": commit, "sha256": paths}


def generate(source_commit: str) -> dict:
    cfg = config()
    curve, field = fixed_curve()
    prior_rows, prior_digest, prior, prior_points, prior_tree_commit = inventory()
    source = source_lock(source_commit)
    assert not (HERE / "FROZEN.json").exists()
    assert not (HERE / "fixtures").exists()
    seen: set[int] = set()
    pending: list[tuple[str, bytes]] = []
    blocks = []
    for block in range(cfg["blocks"]):
        accepted, point_rows = [], []
        rejected = {"zero_scalar": 0, "prior_orbit": 0, "new_orbit": 0}
        for index in range(cfg["target_candidate_cap_per_block"]):
            scalar, digest = candidate(cfg["target_domain"], block, index)
            if scalar == 0:
                rejected["zero_scalar"] += 1
                continue
            q = curve.mul(scalar, GENERATOR)
            assert q is not None and curve.on_curve(q)
            key = orbit_key(q[0], field)
            if key in prior:
                rejected["prior_orbit"] += 1
                continue
            if key in seen:
                rejected["new_orbit"] += 1
                continue
            seen.add(key)
            accepted.append({
                "kind": "base_window_public_fixture",
                "block": block,
                "fixture_index": len(accepted),
                "candidate_index": index,
                "candidate_sha256": digest,
                "published_fixture_scalar": scalar,
                "published_q": list(q),
                "orbit_key": key,
            })
            point_rows.append(list(q))
            if len(accepted) == cfg["public_targets_per_block"]:
                break
        assert len(accepted) == cfg["public_targets_per_block"], (
            block, "candidate cap exhausted"
        )
        prefix = f"fixtures/n41_L1024_b{block:02d}"
        fixture_file, points_file = prefix + ".fixture.jsonl", prefix + ".points.jsonl"
        fixture_bytes = b"".join(canon_bytes(row) for row in accepted)
        points_bytes = b"".join(canon_bytes(row) for row in point_rows)
        pending.extend(((fixture_file, fixture_bytes), (points_file, points_bytes)))
        blocks.append({
            "block": block,
            "corpus": f"base-window-screen-n41-L1024-b{block:02d}-20261001-v1",
            "rho_seed": cfg["rho_seed_base"] + block,
            "candidate_attempts": accepted[-1]["candidate_index"] + 1,
            "candidate_rejections": rejected,
            "fixture_file": fixture_file,
            "fixture_sha256": sha_bytes(fixture_bytes),
            "points_file": points_file,
            "points_sha256": sha_bytes(points_bytes),
        })
    frozen = {
        "schema": "ecc2k130-base-window-input-freeze-v1",
        "curve_slug": cfg["curve_slug"],
        "n": FIELD_DEGREE,
        "a": CURVE_A,
        "field_modulus_low_terms": FIELD_LOW_TERMS,
        "subgroup_order": ORDER,
        "cofactor": COFACTOR,
        "generator": list(GENERATOR),
        "source_lock": source,
        "prior_point_tree_commit": prior_tree_commit,
        "protocol_config_sha256": sha(HERE / "CONFIG.json"),
        "protocol_sha256": sha(HERE / "PROTOCOL.md"),
        "prepare_sha256": sha(Path(__file__)),
        "prior_point_inventory": prior_rows,
        "prior_inventory_digest": prior_digest,
        "prior_n41_rows": prior_points,
        "prior_n41_orbits": len(prior),
        "new_n41_orbits": len(seen),
        "candidate_hash_domain": cfg["target_domain"],
        "blocks": blocks,
    }
    (HERE / "fixtures").mkdir()
    for relative, data in pending:
        (HERE / relative).write_bytes(data)
    (HERE / "FROZEN.json").write_text(json.dumps(frozen, sort_keys=True, indent=2) + "\n")
    return frozen


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-merge-commit", required=True)
    args = parser.parse_args()
    result = generate(args.source_merge_commit)
    print(json.dumps({
        "status": "PASS",
        "new_n41_orbits": result["new_n41_orbits"],
        "prior_inventory_digest": result["prior_inventory_digest"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()

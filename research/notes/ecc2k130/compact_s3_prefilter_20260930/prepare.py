#!/usr/bin/env python3
"""Generate disjoint known-answer fixtures and separate point-only inputs."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

SEED = 623051
N_VALUES = (41, 53)
LENGTH = 1024
ROOT = Path(__file__).resolve().parents[4]
PRIOR_PANELS = (
    "compact_orbit_point_panel_20260929",
    "compact_orbit_strong_rho_20260929",
    "compact_k_boundary_20260929",
    "compact_swap_quotient_20260929",
    "compact_s3_batch_20260929/pilot_622935",
    "compact_s3_batch_20260929",
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def corpus(n: int) -> str:
    return f"compact-s3-prefilter-n{n}-L1024-eval-20260930-v1"


def generate(rho_v2: Path, out: Path) -> dict:
    assert not out.exists(), "never overwrite frozen point files"
    out.mkdir(parents=True)
    specs = {}
    for n in N_VALUES:
        name = f"n{n}_L1024_eval"
        env = dict(os.environ)
        for key in ("KIC_RHO_POINT_INPUT", "KIC_RHO_CANON_BACKEND"):
            env.pop(key, None)
        env.update(KIC_RHO_BATCH_CORPUS=corpus(n), KIC_RHO_GENERATE_ONLY="1")
        command = [str(rho_v2.resolve()), str(n), "0", "signed_frobenius",
                   str(LENGTH), str(SEED)]
        result = subprocess.run(command, env=env, capture_output=True, text=True,
                                check=True)
        assert not result.stderr.strip(), result.stderr
        rows = [json.loads(line) for line in result.stdout.splitlines()]
        assert len(rows) == LENGTH
        new_points = {tuple(row["published_q"]) for row in rows}
        assert len(new_points) == LENGTH
        prior_hashes = {}
        prior_points = set()
        for panel in PRIOR_PANELS:
            folder = ROOT / "research/notes/ecc2k130" / panel / "fixtures"
            files = sorted(folder.glob(f"n{n}_*.points.jsonl"))
            assert files, (panel, n)
            for path in files:
                prior_hashes[str(path.relative_to(ROOT))] = sha(path)
                prior_points.update(tuple(json.loads(line))
                                    for line in path.read_text().splitlines())
        assert not (new_points & prior_points), "new Q overlap an earlier point panel"
        assert all(row["kind"] == "rho_ks_public_fixture"
                   and row["fixture_index"] == index
                   and row["n"] == n and row["a"] == 0
                   and row["batch_seed"] == SEED
                   and row["corpus"] == corpus(n)
                   for index, row in enumerate(rows))
        fixture = out / f"{name}.fixture.jsonl"
        points = out / f"{name}.points.jsonl"
        fixture.write_text(result.stdout)
        points.write_text("".join(json.dumps(row["published_q"],
                                            separators=(",", ":")) + "\n"
                                  for row in rows))
        specs[name] = {
            "n": n, "a": 0, "L": LENGTH, "corpus": corpus(n), "seed": SEED,
            "fixture_file": f"fixtures/{fixture.name}",
            "fixture_sha256": sha(fixture),
            "points_file": f"fixtures/{points.name}",
            "points_sha256": sha(points),
            "generator": rows[0]["generator"],
            "subgroup_order": rows[0]["subgroup_order"],
            "field_modulus_low_terms": rows[0]["field_modulus_low_terms"],
            "automorphism_size": rows[0]["automorphism_size"],
            "prior_point_sha256": prior_hashes,
            "prior_points_checked": len(prior_points),
            "prior_overlap_count": 0,
        }
    return specs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rho-v2", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(generate(args.rho_v2, args.out), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

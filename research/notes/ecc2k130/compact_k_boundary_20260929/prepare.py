#!/usr/bin/env python3
"""Generate disjoint known-answer fixtures and separate point-only inputs."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

SEED = 622933
N_VALUES = (41, 53)
LENGTH = 1024


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def corpus(n: int) -> str:
    return f"compact-k-n{n}-L1024-eval-20260929-v1"


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
        assert len({tuple(row["published_q"]) for row in rows}) == LENGTH
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

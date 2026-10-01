#!/usr/bin/env python3
"""Generate frozen known-answer fixtures, then release only point JSONL to the arms."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess


SPECS = [(n, length, "eval")
         for n in (37, 41, 53)
         for length in (1, 1024)]
SEED = 622901


def corpus(n: int, length: int, kind: str) -> str:
    return f"strong-rho-n{n}-L{length}-{kind}-20260929-v1"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(rho: Path, destination: Path) -> dict:
    destination.mkdir(parents=True, exist_ok=True)
    frozen = {}
    for n, length, kind in SPECS:
        name = f"n{n}_L{length}_{kind}"
        fixture_path = destination / f"{name}.fixture.jsonl"
        points_path = destination / f"{name}.points.jsonl"
        assert not fixture_path.exists() and not points_path.exists(), name
        name_corpus = corpus(n, length, kind)
        env = dict(__import__("os").environ)
        env.pop("KIC_RHO_POINT_INPUT", None)
        env["KIC_RHO_BATCH_CORPUS"] = name_corpus
        env["KIC_RHO_GENERATE_ONLY"] = "1"
        command = [str(rho.resolve()), str(n), "0", "signed_frobenius",
                   str(length), str(SEED)]
        result = subprocess.run(command, env=env, capture_output=True, text=True,
                                check=True)
        assert not result.stderr.strip(), result.stderr
        records = [json.loads(line) for line in result.stdout.splitlines()]
        assert len(records) == length
        assert all(record["kind"] == "rho_ks_public_fixture"
                   and record["fixture_index"] == index
                   and record["n"] == n and record["a"] == 0
                   and record["corpus"] == name_corpus
                   for index, record in enumerate(records))
        fixture_path.write_text(result.stdout)
        points_path.write_text("".join(json.dumps(record["published_q"],
                                                  separators=(",", ":")) + "\n"
                                       for record in records))
        frozen[name] = {
            "n": n, "a": 0, "L": length, "kind": kind,
            "corpus": name_corpus, "seed": SEED,
            "fixture_file": str(fixture_path.relative_to(destination.parent)),
            "fixture_sha256": sha(fixture_path),
            "points_file": str(points_path.relative_to(destination.parent)),
            "points_sha256": sha(points_path),
            "subgroup_order": records[0]["subgroup_order"],
            "automorphism_size": records[0]["automorphism_size"],
            "generator": records[0]["generator"],
            "field_modulus_low_terms": records[0]["field_modulus_low_terms"],
        }
    return frozen


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rho", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.rho, args.out), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

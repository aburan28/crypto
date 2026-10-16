#!/usr/bin/env sage -python
"""Freeze twelve one-target public points after the S3 source freeze."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing


HERE = Path(__file__).resolve().parent
OUT = HERE / "holdout"
COUNT = 12
SEED = b"s3-paired-inversion-holdout-20261007-v1"


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")


def write_immutable(path: Path, value: dict | list) -> None:
    data = json.dumps(value, indent=2, sort_keys=True) + "\n"
    if path.exists():
        if path.read_text() != data:
            raise RuntimeError(f"frozen artifact would change: {path}")
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(data)


def main() -> None:
    primary = json.loads((HERE / "workload.json").read_text())
    baseline = json.loads((HERE / "baseline_manifest.json").read_text())
    candidate = json.loads((HERE / "candidate_manifest.json").read_text())
    freeze = json.loads((HERE / "freeze_receipt.json").read_text())
    for name, manifest in (("baseline", baseline), ("candidate", candidate)):
        path = HERE / f"{name}_manifest.json"
        assert hashlib.sha256(path.read_bytes()).hexdigest() == (
            freeze["variants"][name]["manifest_sha256"]
        )
        source = Path(manifest["candidate_freeze"]["source_path"])
        binary = Path(manifest["candidate_freeze"]["binary_path"])
        assert hashlib.sha256(source.read_bytes()).hexdigest() == (
            manifest["candidate_freeze"]["source_sha256"]
        )
        assert hashlib.sha256(binary.read_bytes()).hexdigest() == (
            manifest["candidate_freeze"]["binary_sha256"]
        )

    curve = candidate["identity_record"]["curve"]
    n = candidate["identity_record"]["field"]["degree"]
    r = curve["subgroup_order"]
    assert n == 53
    R = PolynomialRing(GF(2), "t")
    t = R.gen()
    F = GF(2**n, name="z", modulus=t**53 + t**6 + t**2 + t + 1)
    powers = [F.gen() ** i for i in range(n)]

    def from_word(word: int):
        value = F(0)
        while word:
            bit = (word & -word).bit_length() - 1
            value += powers[bit]
            word &= word - 1
        return value

    def to_word(value) -> int:
        poly = value.polynomial()
        result = sum(int(c) << i for i, c in enumerate(poly.list()))
        assert 0 <= result < 1 << n and from_word(result) == value
        return result

    E = EllipticCurve(F, [F(1), F(0), F(0), F(0), F(1)])
    G = E(*(from_word(word) for word in curve["generator"]))
    assert r * G == E(0)
    fixtures = []
    seen = set()
    for index in range(1, COUNT + 1):
        digest = hashlib.sha256(SEED + index.to_bytes(4, "big")).digest()
        scalar = 1 + (int.from_bytes(digest, "big") % (r - 1))
        assert scalar not in seen
        seen.add(scalar)
        point = scalar * G
        words = [to_word(point[0]), to_word(point[1])]
        assert E(*(from_word(word) for word in words)) == point
        label = f"T{index:02d}"
        target_path = OUT / label / "public_target.json"
        write_immutable(target_path, words)
        workload = copy.deepcopy(primary["record"])
        workload["targets"] = [words]
        workload["target_generation_seed"] = f"{SEED.decode()}:T{index:02d}"
        workload["input_law"] = (
            "one deterministic pseudorandom public subgroup point supplied "
            "as coordinates; generation scalar is outside both timed intervals "
            "and is absent from the solver input"
        )
        workload["comparison_question"] = (
            "fresh-target one-target correctness and latency variability for "
            "paired baseline versus shared-inversion S3 solver"
        )
        workload["rho_reference_id"] = None
        workload["rho_fixture_or_walk_seeds"] = []
        workload["cold_or_warm_target_count"]["rho_queries"] = 0
        identity = hashlib.sha256(canonical(workload)).hexdigest()
        workload_id = identity[:12]
        write_immutable(OUT / label / "workload.json", {
            "identity_sha256": identity,
            "record": workload,
            "workload_id": workload_id,
        })
        fixtures.append({
            "label": label,
            "fixture_scalar": scalar,
            "public_point": words,
            "public_target_sha256": hashlib.sha256(target_path.read_bytes()).hexdigest(),
            "scalar_derivation_sha256": digest.hex(),
            "workload_id": workload_id,
        })
    manifest = {
        "kind": "s3_pair_fresh_target_fixture_generation",
        "question": "Does paired S3 inversion preserve correctness across fresh one-target workloads?",
        "count": COUNT,
        "seed_label": SEED.decode(),
        "source_candidate_ids": {
            "baseline": baseline["candidate_id"],
            "candidate": candidate["candidate_id"],
        },
        "fixture_scalar_visibility": (
            "stored only in this generator output and independent replay; "
            "each solver reads only its public_target.json"
        ),
        "fixtures": fixtures,
    }
    write_immutable(OUT / "fixtures.json", manifest)
    print(json.dumps({"count": COUNT, "workload_ids": [f["workload_id"] for f in fixtures]}))


if __name__ == "__main__":
    main()

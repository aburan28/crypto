"""Freeze matched fresh N41/N53 public targets through checked Sage.

Run with /Volumes/SSD990/cryptanalysis/sage -python after saving
the launcher's --runtime-info output. No solver reads fixture scalars.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
COUNT = 200
SEED_LABEL = "s3-crosscurve-rank243-20261009-v1"
SEED_PREFIX = b"s3-pair-target-v1:"
RELATION_SEED = 20260928
WORKERS = 14
CURVES = {
    41: {"r": 549756390943, "cofactor": 4,
         "modulus_powers": [41, 3, 0],
         "generator": [2056947637384, 1635505394702],
         "curve_id": "EC1N41Ckb1hcbf4a37f6e66",
         "factor_base_points": 20008,
         "reference_candidate_id": "IC1N41Ckb1fb20008PDP4rootRCsampleLAgaussTDdirectISO0h3f64f5604376",
         "rank243_candidate_id": "IC1N41Ckb1fb20008PDP4rootRCsampleLAgaussTDdirectISO0h0f38cf1ef60a"},
    53: {"r": 21044858204113, "cofactor": 428,
         "modulus_powers": [53, 6, 2, 1, 0],
         "generator": [198217578752339, 7929897206038174],
         "curve_id": "EC1N53Ckb1h738a572b91f2",
         "factor_base_points": 25864,
         "reference_candidate_id": "IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0h5ed3e1d8954b",
         "rank243_candidate_id": "IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0hbd278c9cf111"},
}


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def digest(value: object) -> str:
    return hashlib.sha256(canonical(value)).hexdigest()


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_immutable(path: Path, content: str) -> None:
    if path.exists():
        assert path.read_text() == content, f"frozen artifact changed: {path}"
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)


def main() -> None:
    prior = json.loads((ROOT / "experiments/koblitz-s3-partial-prepared-20261009/rank243/FREEZE.json").read_text())
    assert prior["rank_policies"]["deep"] == 243
    assert prior["relation_seed"] == RELATION_SEED and prior["workers"] == WORKERS
    reference_path = ROOT / "experiments/koblitz-s3-pair-query-20261007-v4/candidate.rs"
    partial_path = ROOT / "examples/s3_partial_prepared_targetspan.rs"
    assert sha(reference_path) == prior["reference_source_sha256"]
    assert sha(partial_path) == prior["partial_source_sha256"]
    polynomial = PolynomialRing(GF(2), "t")
    t = polynomial.gen()
    targets = []
    for n, cfg in CURVES.items():
        field = GF(2**n, name="z", modulus=sum(t**p for p in cfg["modulus_powers"]))
        powers = [field.gen() ** bit for bit in range(n)]

        def from_word(word: int):
            assert 0 <= word < 1 << n
            result = field(0)
            while word:
                bit = (word & -word).bit_length() - 1
                result += powers[bit]
                word &= word - 1
            return result

        def to_word(value) -> int:
            result = sum(int(bit) << index for index, bit in enumerate(value.polynomial().list()))
            assert from_word(result) == value
            return result

        curve = EllipticCurve(field, [field(1), field(0), field(0), field(0), field(1)])
        generator = curve(*(from_word(word) for word in cfg["generator"]))
        assert cfg["r"] * generator == curve(0)
        fixtures = []
        seen_scalars = set()
        seen_points = set()
        prior_points = set()
        prior_root = ROOT / "experiments/koblitz-s3-partial-prepared-20261009"
        for relative in (f"holdout/n{n}/fixtures.json",
                         f"rank242/holdout/n{n}/fixtures.json",
                         f"rank243/holdout/n{n}/fixtures.json"):
            prior_path = prior_root / relative
            if prior_path.exists():
                prior_points.update(tuple(item["public_point"])
                                    for item in json.loads(prior_path.read_text())["fixtures"])
        for index in range(1, COUNT + 1):
            label = f"T{index:03}"
            raw_seed = (SEED_PREFIX + SEED_LABEL.encode() + b":" +
                        n.to_bytes(4, "big") + index.to_bytes(4, "big"))
            scalar_hash = hashlib.sha256(raw_seed).digest()
            scalar = 1 + int.from_bytes(scalar_hash, "big") % (cfg["r"] - 1)
            assert scalar not in seen_scalars
            seen_scalars.add(scalar)
            point = scalar * generator
            words = [to_word(point[0]), to_word(point[1])]
            assert tuple(words) not in seen_points and tuple(words) not in prior_points
            seen_points.add(tuple(words))
            assert curve(*(from_word(word) for word in words)) == point
            relative = f"holdout/n{n}/{label}/public_target.json"
            public_path = HERE / relative
            write_immutable(public_path, json.dumps(words, separators=(",", ":")) + "\n")
            fixture = {"n": n, "label": label, "family": "fresh_matched_curve",
                       "point_file_relative": relative,
                       "point_sha256": sha(public_path), "public_point": words,
                       "fixture_scalar": scalar,
                       "scalar_derivation_sha256": scalar_hash.hex()}
            workload = {
                "curve_id": cfg["curve_id"], "subgroup_order_r": cfg["r"],
                "target_count": 1, "target_point": words, "target_label": label,
                "target_input_law": "known_scalar_sha256_public_point_fixture_v1",
                "target_seed_label": SEED_LABEL, "target_index": index,
                "scalar_derivation_sha256": scalar_hash.hex(),
                "relation_seed": RELATION_SEED,
                "preparation_state_at_online_start": "candidate_policy_reusable_setup_ready"}
            workload_id = digest(workload)[:12]
            fixture["workload_id"] = workload_id
            fixture["workload_identity_sha256"] = digest(workload)
            write_immutable(HERE / "workloads" / f"W{workload_id}.json",
                            json.dumps({"workload_id": workload_id, "record": workload,
                                        "identity_sha256": digest(workload)},
                                       indent=2, sort_keys=True) + "\n")
            fixtures.append(fixture)
            targets.append(fixture)
        write_immutable(HERE / f"holdout/n{n}/fixtures.json",
                        json.dumps({"kind": "s3_crosscurve_fresh_public_fixtures_v1",
                                    "seed_label": SEED_LABEL, "count": COUNT,
                                    "curve_id": cfg["curve_id"], "fixtures": fixtures},
                                   indent=2, sort_keys=True) + "\n")
    plan = {
        "kind": "s3_rank243_matched_curve_fresh_panel_v1",
        "seed_label": SEED_LABEL,
        "input_law": "1 + int_be(SHA256('s3-pair-target-v1:' || seed_label || ':' || u32be(n) || u32be(index))) mod (r-1)",
        "target_count_per_curve": COUNT,
        "initial_pilot_count_per_curve": 24,
        "one_public_target_per_process": True,
        "relation_seed": RELATION_SEED,
        "workers": WORKERS,
        "timeout_seconds_per_run": 180,
        "repetitions_per_arm_per_target": 2,
        "paired_order": "ABBA_or_BAAB_rotated_by_target_and_curve; one additional prepared A/A run",
        "curves": {str(n): cfg for n, cfg in CURVES.items()},
        "reference_source_sha256": sha(reference_path),
        "partial_source_sha256": sha(partial_path),
        "reference_binary_sha256_local": prior["reference_binary_sha256"],
        "partial_binary_sha256_local": prior["partial_binary_sha256"],
        "parent_rank243_plan_sha256": sha(ROOT / "experiments/koblitz-s3-partial-prepared-20261009/rank243/FREEZE.json"),
        "controlled_speedup": None,
        "targets": targets,
    }
    write_immutable(HERE / "FREEZE.json", canonical(plan).decode() + "\n")
    print(json.dumps({"plan_sha256": sha(HERE / "FREEZE.json"),
                      "targets_per_curve": COUNT, "workloads": len(targets)}))


if __name__ == "__main__":
    main()

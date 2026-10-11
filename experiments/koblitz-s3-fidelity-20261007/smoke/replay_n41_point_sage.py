"""Independent Sage point replay for the N41 S3 feasibility smoke.

This checks the public curve, target and recovered scalar. It does not replay
the factor-base relation witnesses or certify the unisolated timing ratio.
Launch through /Volumes/SSD990/cryptanalysis/sage -python.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing
from sage.env import SAGE_VERSION


HERE = Path(__file__).resolve().parent
N = 41
R = 549756390943
GENERATOR = [2056947637384, 1635505394702]
FIXTURE_SCALAR = 123212651130


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    receipt_path = HERE / "independent_point_replay.json"
    if receipt_path.exists():
        raise FileExistsError(f"refusing to overwrite {receipt_path}")
    runtime_path = HERE / "sage_runtime_info.json"
    assert runtime_path.is_file(), "save checked Sage --runtime-info first"
    runtime = json.loads(runtime_path.read_text())
    assert runtime["status"] == "verified"
    assert runtime["sage_version"] == SAGE_VERSION

    polynomial = PolynomialRing(GF(2), "t")
    t = polynomial.gen()
    field = GF(2**N, name="z", modulus=t**41 + t**3 + 1)
    powers = [field.gen() ** i for i in range(N)]

    def from_word(word: int):
        assert isinstance(word, int) and 0 <= word < 1 << N
        value = field(0)
        while word:
            bit = (word & -word).bit_length() - 1
            value += powers[bit]
            word &= word - 1
        return value

    curve = EllipticCurve(field, [field(1), field(0), field(0), field(0), field(1)])

    def point(words: list[int]):
        assert len(words) == 2
        return curve(from_word(words[0]), from_word(words[1]))

    generator = point(GENERATOR)
    assert generator != curve(0) and R * generator == curve(0)
    target_path = HERE / "public_target.json"
    target_words = json.loads(target_path.read_text())
    target = point(target_words)
    assert target != curve(0) and R * target == curve(0)
    assert FIXTURE_SCALAR * generator == target

    reports = {}
    for arm in ("baseline", "candidate"):
        path = HERE / f"{arm}.jsonl"
        report = json.loads(path.read_text())
        assert report["n"] == N and report["a"] == 0
        assert report["subgroup_order"] == R
        assert report["target"] == target_words
        assert report["group_verified"] is True
        recovered = int(report["recovered_scalar"])
        assert recovered == FIXTURE_SCALAR
        assert recovered * generator == target
        reports[arm] = report

    baseline, candidate = reports["baseline"], reports["candidate"]
    equal_fields = [
        "factor_base_points", "factor_base_digest", "orbit_columns",
        "rank", "rank_attempts", "rank_attempts_completed", "rank_new_rows",
        "rank_dependent_rows", "rank_failures", "target_relation_indices",
        "target_span_stop_relation_prefix", "recovered_scalar",
    ]
    assert all(baseline[key] == candidate[key] for key in equal_fields)
    assert len(baseline["rank_relation_witnesses"]) == 252
    assert len(candidate["rank_relation_witnesses"]) == 252
    for left, right in zip(baseline["rank_relation_witnesses"],
                           candidate["rank_relation_witnesses"]):
        for key in ("point_indices", "rank_gain", "relation_scalar", "x_codes"):
            assert left[key] == right[key]

    receipt = {
        "kind": "s3_n41_smoke_independent_sage_point_replay_v1",
        "scope": "curve, subgroup, fixture target, recovered scalar, and baseline/candidate witness-field agreement; relation witnesses not independently replayed",
        "sage_version": SAGE_VERSION,
        "sage_runtime_info_sha256": sha(runtime_path),
        "script_sha256": sha(Path(__file__)),
        "public_target_sha256": sha(target_path),
        "raw_sha256": {arm: sha(HERE / f"{arm}.jsonl") for arm in reports},
        "curve": "y^2 + xy = x^3 + 1 over GF(2^41), polynomial t^41+t^3+1",
        "subgroup_order": R,
        "generator": GENERATOR,
        "target": target_words,
        "recovered_scalar": FIXTURE_SCALAR,
        "point_replay_verified": True,
        "witness_fields_equal": True,
        "relation_witnesses_independently_verified": False,
        "cpu_timing_status": "exploratory_unisolated",
    }
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"point_replay_verified": True,
                      "witness_fields_equal": True,
                      "relation_witnesses": 252}))


if __name__ == "__main__":
    main()

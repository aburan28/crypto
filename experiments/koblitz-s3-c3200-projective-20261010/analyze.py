"""Audit the frozen T025 projective scalar-control archive."""

from __future__ import annotations

import gzip
import hashlib
import json
import math
from pathlib import Path
import statistics


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ARCHIVE = HERE / "local/verified_t025"
REFERENCE = ROOT / "experiments/koblitz-s3-c3200-exact-20261010/local/verified_t025/report.json.gz"
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def audit() -> dict:
    freeze_raw = (HERE / "CONTROL_FREEZE.json").read_bytes()
    freeze = json.loads(freeze_raw)
    receipt = json.loads((ARCHIVE / "receipt.json").read_bytes())
    assert receipt["freeze_sha256"] == digest(freeze_raw)
    assert receipt["source_sha256"] == freeze["source_sha256"]
    assert receipt["binary_sha256"] == freeze["binary_sha256"]
    assert all(receipt[key] is True for key in (
        "scalar_verified", "diagnostic_verified", "full_rank_before_target",
        "index_layout_ok", "memory_within_envelope", "phase_accounting_ok",
        "stdout_matches_report"))
    assert receipt["unexpected_reference_differences"] == []
    assert receipt["peak_rss_bytes"] <= freeze["memory_envelope_bytes"]

    archived = {}
    for name, expected in receipt["archive"].items():
        packed = (ARCHIVE / (name + ".gz")).read_bytes()
        raw = gzip.decompress(packed)
        assert digest(packed) == expected["gzip_sha256"]
        assert digest(raw) == expected["raw_sha256"]
        assert len(packed) == expected["gzip_bytes"]
        assert len(raw) == expected["raw_bytes"]
        archived[name] = raw
    assert archived["report.json"] == archived["stdout.bin"]
    assert archived["stderr.txt"] == b""
    report = json.loads(archived["report.json"])
    reference = json.loads(gzip.decompress(REFERENCE.read_bytes()))
    assert digest(REFERENCE.read_bytes()) == freeze["reference_report_gzip_sha256"]
    assert report["group_verified"] and report["recovered_scalar"] == freeze["fixture_scalar"]
    assert report["target_relation_indices"] == reference["target_relation_indices"]
    assert report["target_s3_calls"] == reference["target_s3_calls"]
    assert report["prepared_rank"] == freeze["prepare_rank"]
    assert report["factor_base_points"] == 339200
    assert receipt["factor_base_digest"] == reference["factor_base_digest"]
    for source in (report, reference):
        timing = source["timing_ms"]
        assert abs(sum(timing[key] for key in PHASES) -
                   timing["target_online_after_reusable_setup"]) < .02

    rows = [json.loads(line) for line in archived["diagnostic.jsonl"].splitlines()]
    assert len(rows) == freeze["diagnostic_pairs"] * 2
    metrics = {}
    for phase in ("validation", "replay"):
        for clock in ("wall", "thread_cpu"):
            key = f"{phase}_{clock}_ns"
            samples = []
            for pair in range(freeze["diagnostic_pairs"]):
                first, second = rows[2 * pair:2 * pair + 2]
                by_mode = {row["mode_name"]: row for row in (first, second)}
                assert first["pair"] == second["pair"] == pair
                assert set(by_mode) == {"affine", "projective"}
                assert first["mode_name"] == ("affine" if pair % 2 == 0 else "projective")
                for row in (first, second):
                    assert row["subgroup_verified"] and row["scalar_replay_verified"]
                    assert row["target"] == freeze["public_point"]
                    assert row["recovered_scalar"] == freeze["fixture_scalar"]
                affine, projective = by_mode["affine"][key], by_mode["projective"][key]
                if affine is not None and projective is not None:
                    assert isinstance(affine, int) and isinstance(projective, int)
                    assert affine > 0 and projective > 0
                    samples.append((affine, projective))
            ratios = sorted(affine / projective for affine, projective in samples)
            metrics[f"{phase}_{clock}"] = {
                "complete_pairs": len(samples),
                "affine_median_ns": statistics.median(a for a, _ in samples),
                "projective_median_ns": statistics.median(b for _, b in samples),
                "paired_ratio_median": statistics.median(ratios),
                "paired_ratio_geometric_mean": math.exp(statistics.mean(map(math.log, ratios))),
                "paired_ratio_p10": ratios[int(.1 * (len(ratios) - 1))],
                "paired_ratio_p90": ratios[int(.9 * (len(ratios) - 1))],
            }

    manifest_path, = (HERE / "identities/candidates").glob("*.json")
    manifest_raw = manifest_path.read_bytes()
    manifest = json.loads(manifest_raw)
    candidate_id = manifest.pop("candidate_id")
    canonical = json.dumps(manifest, sort_keys=True, separators=(",", ":"),
                           ensure_ascii=False).encode()
    assert candidate_id.endswith("h" + digest(canonical)[:12])
    assert manifest["implementation"]["example_source_sha256"] == freeze["source_sha256"]
    replay = json.loads((ARCHIVE / "independent_sage_replay.json").read_bytes())
    assert replay["candidate_id"] == candidate_id
    assert replay["candidate_manifest_sha256"] == digest(manifest_raw)
    assert replay["report_raw_sha256"] == digest(archived["report.json"])
    assert replay["generator_order_verified"]
    assert replay["scalar_times_generator_equals_public_target"]

    return {
        "kind": "s3_n53_c3200_projective_control_analysis_v1",
        "candidate_id": candidate_id,
        "workload_id": freeze["workload_id"],
        "freeze_sha256": digest(freeze_raw),
        "source_sha256": freeze["source_sha256"],
        "binary_sha256": freeze["binary_sha256"],
        "verified_target_scalar": report["recovered_scalar"],
        "same_first_relation_witness": True,
        "same_target_s3_calls": report["target_s3_calls"],
        "prepared_rank": report["prepared_rank"],
        "actual_usable_points_B": report["factor_base_points"],
        "peak_rss_bytes": receipt["peak_rss_bytes"],
        "online_phase_ms": {"reference": {key: reference["timing_ms"][key] for key in PHASES},
                            "projective": {key: report["timing_ms"][key] for key in PHASES}},
        "online_ms": {"reference": reference["timing_ms"]["target_online_after_reusable_setup"],
                      "projective": report["timing_ms"]["target_online_after_reusable_setup"]},
        "online_ratio_exploratory": (reference["timing_ms"]["target_online_after_reusable_setup"] /
                                     report["timing_ms"]["target_online_after_reusable_setup"]),
        "diagnostic_scope": "256 balanced check pairs on one already-verified public point",
        "diagnostic": metrics,
        "independent_sage_replay_verified": True,
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }


if __name__ == "__main__":
    result = audit()
    (HERE / "analysis.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"candidate_id": result["candidate_id"],
                      "online_ratio_exploratory": result["online_ratio_exploratory"],
                      "diagnostic": result["diagnostic"]}, sort_keys=True))

#!/usr/bin/env python3
"""Apply the preregistered admission rule to four sealed N83 screens."""

import argparse
import hashlib
import json
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
SOURCE_SHA256 = "f2745227d95e36ef45c73c5104a23ada6cebb0794405368581d2be17044845dc"
TARGET_SHA256 = "cfa943ae538fcc41d2af8bb8fd1c8fc0f29adc340b55644e90382da25b213424"
ANSWER_SHA256 = "1dbd4250796de4434badc2458089d7ad9fe63cc59d7b9d3f6825ad431476a13b"
PANELS = (
    ("direct", "n83-best-e2e-38768544528", 38768544528,
     "ae6a10b641e17629bd282d8205787ede129dbf98",
     "n83.best-e2e-analysis/v1", "n83_best_e2e_results_20261015", 4000),
    ("legacy", "n83-legacy-challenger-38770004923", 38770004923,
     "d6293e278a037ef34a6de0101844fcfbf97703e2",
     "n83.legacy-challenger-analysis/v1",
     "n83_legacy_challenger_results_20261015", 1000),
    ("root4_bloom", "n83-root4-bloom-38772262240", 38772262240,
     "afea7179c40c9e00f8701a32c596e47cd8129f62",
     "n83.root4-bloom-analysis/v1", "n83_root4_bloom_results_20261015", 1000),
    ("root11_highk", "n83-root11-highk-38773883285", 38773883285,
     "e66475a416dead205ec797ca16f025abcea34d16",
     "n83.root11-highk-analysis/v1", "n83_root11_highk_results_20261015", 1000),
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def method(row):
    return (row["k"], row["variant"], row["bloom_profile"],
            row["base_hash"], row["root_table_entries"], SOURCE_SHA256)


def eligible_pipeline(analysis, pipeline):
    rows = [row for row in analysis["arms"] if row.get("pipeline") == pipeline]
    if len(rows) != 2 or {row["arm"] for row in rows} != {
            pipeline + "_a", pipeline + "_b"}:
        return None
    rows.sort(key=lambda row: row["arm"])
    a, b = rows
    for row in rows:
        if (row["status"] != "PASS_INDEPENDENT_REPLAY"
                or not row["resource_envelope_status"].startswith("WITHIN_")
                or row["cpu_throttled_periods_delta"] != 0
                or row["oom_kills_delta"] != 0
                or row["actual_factor_base_points"] != 166 * row["k"]
                or row["effective_columns"] != row["k"]
                or len(row["targets"]) != analysis["target_count"]
                or any(target["scalar_replay"] != "PASS"
                       or target["relation_replay"] != "PASS"
                       for target in row["targets"])):
            return None
    if method(a) != method(b):
        return None
    paired = [(x["online_ms_exploratory"] + y["online_ms_exploratory"]) / 2
              for x, y in zip(a["targets"], b["targets"])]
    return {"pipeline": pipeline, "k": a["k"], "variant": a["variant"],
            "bloom_profile": a["bloom_profile"],
            "actual_usable_points": a["actual_factor_base_points"],
            "base_hash": a["base_hash"],
            "root_table_entries": a["root_table_entries"],
            "screen_mean_online_ms": statistics.mean(paired),
            "screen_peak_rss_bytes": max(a["peak_rss_bytes"], b["peak_rss_bytes"]),
            "arm_names": [a["arm"], b["arm"]]}


def screen(spec, artifact_dir):
    label, archive_name, run_id, commit, schema, result_dir, count = spec
    zip_path = artifact_dir / (archive_name + ".zip")
    root = artifact_dir / archive_name / "results"
    selection_path = REPO / "research" / result_dir / "selection.json"
    if not (zip_path.is_file() and (root / "analysis.json").is_file()
            and selection_path.is_file()):
        raise ValueError(f"{label}: source artifact or sealed result missing")
    selection = read(selection_path)
    analysis = read(root / "analysis.json")
    source = read(root / "source.json")
    host = read(root / "host.json")
    artifact_sha = digest(zip_path)
    analysis_sha = digest(root / "analysis.json")
    if (artifact_sha != selection["artifact_sha256"]
            or analysis_sha != selection["analysis_sha256"]
            or selection["source_commit"] != commit
            or host["git_sha"] != commit
            or source["hashes"]["examples/koblitz_orbit_dlp_fast.rs"] != SOURCE_SHA256
            or analysis["schema"] != schema
            or analysis["hardware_tier"] != "cloud_instance_exploratory"
            or analysis["target_count"] != count
            or analysis["arm_order"] != source["arm_order"]):
        raise ValueError(f"{label}: source or sealed evidence identity differs")
    pipelines = sorted({row.get("pipeline") for row in analysis["arms"]
                        if row.get("pipeline") is not None})
    candidates = [candidate for pipeline in pipelines
                  if (candidate := eligible_pipeline(analysis, pipeline)) is not None]
    winner = min(candidates,
                 key=lambda row: (row["screen_mean_online_ms"], row["pipeline"])) if candidates else None
    if winner is not None:
        path = REPO / "research" / result_dir / (
            f"candidate_K{winner['k']}_{winner['variant']}_{winner['bloom_profile']}.json")
        manifest = read(path)
        if (manifest["factor_base"]["actual_usable_points"] != winner["actual_usable_points"]
                or manifest["factor_base"]["source_base_hash_blake3"] != winner["base_hash"]
                or manifest["implementation"]["source_sha256"]["ic"] != SOURCE_SHA256):
            raise ValueError(f"{label}: candidate manifest differs from replayed base")
        winner.update({"candidate_id": manifest["candidate_id"],
                       "candidate_manifest_sha256": digest(path)})
    return {"screen": label, "run_id": run_id, "source_commit": commit,
            "artifact_sha256": artifact_sha, "analysis_sha256": analysis_sha,
            "source_status": analysis["status"],
            "source_selected_pipeline": (analysis["balanced_comparison"] or {}).get(
                "selected_pipeline"),
            "eligible_pipeline_count": len(candidates), "representative": winner}, candidates


def admit(artifact_dir):
    if (digest(HERE / "targets.jsonl") != TARGET_SHA256
            or digest(HERE / "answers.jsonl") != ANSWER_SHA256):
        raise ValueError("frozen final public points or scalar sidecar changed")
    screens = []
    by_screen = {}
    for spec in PANELS:
        record, candidates = screen(spec, artifact_dir)
        screens.append(record)
        by_screen[spec[0]] = {row["pipeline"]: row for row in candidates}
    controls = (by_screen["direct"].get("wide3"),
                by_screen["direct"].get("root4k3596"))
    if any(row is None for row in controls):
        raise ValueError("standing controls lack two eligible direct-panel runs")
    order = [("control", row) for row in controls]
    order.extend((record["screen"], record["representative"])
                 for record in screens if record["representative"] is not None)
    admitted = []
    seen = set()
    for origin, row in order:
        key = method(row)
        if key in seen:
            continue
        seen.add(key)
        admitted.append({"position": len(admitted), "origin": origin, **row,
                         "source_sha256": SOURCE_SHA256})
    if not 2 <= len(admitted) <= 6:
        raise ValueError("final candidate count is outside the frozen bounds")
    return {"schema": "n83.final-tournament-admission/v1",
            "status": "FROZEN_ADMISSION", "source_screens": screens,
            "candidate_order": admitted,
            "mirrored_arm_order": [row["pipeline"] for row in admitted] +
                                  [row["pipeline"] for row in reversed(admitted)],
            "target_sha256": TARGET_SHA256,
            "answer_sidecar_sha256": ANSWER_SHA256}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact-dir", required=True, type=Path)
    args = parser.parse_args()
    record = admit(args.artifact_dir.resolve())
    path = HERE / "admission.json"
    blob = json.dumps(record, indent=2, sort_keys=True) + "\n"
    if path.exists():
        if path.read_text() != blob:
            raise ValueError("frozen candidate admission changed")
    else:
        path.write_text(blob)
    print(json.dumps({"candidate_order": record["mirrored_arm_order"],
                      "status": record["status"]}, sort_keys=True))


if __name__ == "__main__":
    main()

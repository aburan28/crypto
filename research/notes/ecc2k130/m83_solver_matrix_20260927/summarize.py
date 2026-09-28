"""Validate and summarize immutable m83 run receipts; fail on any model mismatch."""
from __future__ import annotations

import collections
import json
from pathlib import Path
import statistics

HERE = Path(__file__).resolve().parent
ROOT = HERE / "results"


def receipt(path):
    wrapper = json.loads(path.read_text())
    lines = [json.loads(s) for s in wrapper["stdout"].splitlines() if s.strip()]
    assert lines, (path, wrapper["stderr"])
    return lines[-1], wrapper


def group(paths, oracle_dir=None):
    entries = collections.defaultdict(list)
    accepted = {"complete", "timeout", "oom", "setup_failed"}
    for path in paths:
        record, wrapper = receipt(path)
        assert record["status"] in accepted, (path, record, wrapper["stderr"])
        if record["status"] == "complete":
            assert wrapper["returncode"] == 0
            if oracle_dir is not None:
                oracle_name = (f"fixed-s{record['seed']}-i{record['index']}-fes.json"
                               if oracle_dir.name == "run_20260927" else
                               path.name.replace(record["engine"] + ".json", "fes.json"))
                oracle_path = oracle_dir / oracle_name
                assert oracle_path.exists(), (path, oracle_path)
                reference, _ = receipt(oracle_path)
                assert reference["status"] == "complete"
                for key in ("input_sha256", "equation_sha256", "algebraic_roots"):
                    assert record[key] == reference[key], (path, key)
                assert record["group_verification"] == reference["group_verification"], path
        kind = record.get("input", {}).get("kind", "planted" if record["index"] == 0 else "natural")
        key = (record.get("factor_base_policy", "quotient"),
               record.get("payload_bits", 2), record["engine"])
        entries[key].append((kind, record))
    summary = []
    for (policy, s, engine), rows in sorted(entries.items()):
        successes = [r for _, r in rows if r["status"] == "complete"]
        times = [r["timing"]["solver_including_setup"]["wall_s"] for r in successes]
        summary.append({"policy": policy, "payload_bits": s, "engine": engine,
                        "planned": len(rows), "completed": len(successes),
                        "timeout": sum(r["status"] == "timeout" for _, r in rows),
                        "setup_failed": sum(r["status"] == "setup_failed" for _, r in rows),
                        "natural_verified_queries": sum(kind == "natural" and
                            bool(r["group_verification"]["verified"]) for kind, r in rows
                            if r["status"] == "complete"),
                        "algebraic_false_positives": sum(len(r["group_verification"]["false_positive"])
                            for r in successes),
                        "median_solver_stage_wall_s": statistics.median(times) if times else None,
                        "max_solver_stage_wall_s": max(times) if times else None})
    return summary


def main():
    original = ROOT / "run_20260927"
    extension = ROOT / "extension_20260927"
    assert len(list(original.glob("fixed*.json"))) == 48
    assert len(list(original.glob("phase*.json"))) == 4
    assert len(list((extension / "quotient_s2").glob("*.json"))) == 12
    assert len(list((extension / "plain_subspace").glob("*.json"))) == 60
    replay = json.loads((extension / "quotient_preflight_replay.json").read_text())
    assert len(replay) == 4 and all(x["status"] == "setup_failed" for x in replay)
    base = group(original.glob("fixed*.json"), original)
    native = group((extension / "quotient_s2").glob("*.json"), original)
    plain = group((extension / "plain_subspace").glob("*.json"), extension / "plain_subspace")
    phase_rows = []
    for p in sorted(original.glob("phase*.json")):
        wrapper = json.loads(p.read_text())
        stages = [json.loads(line) for line in wrapper["stdout"].splitlines() if line.strip()]
        phase_rows.append({"case": p.stem, "returncode": wrapper["returncode"],
                           "stages": [{k: v for k, v in row.items()
                                       if k in ("stage", "seconds", "status", "reason",
                                                "nvars", "boolean_terms", "peak_rss_bytes")}
                                      for row in stages],
                           "total_wall_s": wrapper["timing"]["wall_s"],
                           "failed_at": "SAT_setup" if wrapper["returncode"] and
                              stages[-1]["stage"] == "equations" else None})
    orbit_rows = []
    for seed in (260938, 260939):
        path = ROOT / "orbit_20260927" / f"orbit-s{seed}.json"
        wrapper = json.loads(path.read_text())
        assert wrapper["returncode"] == 0, path
        stages = [json.loads(line) for line in wrapper["stdout"].splitlines() if line.strip()]
        assert [s["stage"] for s in stages] == ["setup", "orbit", "pairs"] + ["target"] * 6 + ["final"]
        assert stages[-1]["status"] == "complete" and stages[-1]["natural_rank"] == 0
        for index, target in enumerate(stages[3:9]):
            oracle, _ = receipt(original / f"fixed-s{seed}-i{index}-fes.json")
            assert target["input_sha256"] == oracle["input_sha256"]
            assert target["target"] == oracle["input"]["target"]
            assert target["verified_triples"] >= 1 if index < 2 else target["verified_triples"] == 0
        orbit_rows.append({"seed": seed, "signed_orbit_size": stages[1]["signed_orbit_size"],
                           "pair_attempts": stages[2]["attempts"],
                           "natural_verified_queries": sum(bool(x["verified_triples"]) for x in stages[5:9]),
                           "natural_rank": stages[-1]["natural_rank"],
                           "group_additions": stages[-1]["group_additions"],
                           "wall_seconds": wrapper["timing"]["wall_s"]})
    output = {"fixed_phase_quotient": base + native, "fixed_phase_normal_subspace": plain,
              "phase_choice": phase_rows, "full_orbit_relation_control": orbit_rows,
              "quotient_size_preflight": replay,
              "full_dlp_cost": None, "matched_rho_cost": None,
              "cross_engine_operation_ratio": None}
    (ROOT / "validated_summary.json").write_text(json.dumps(output, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"fixed_phase_quotient": output["fixed_phase_quotient"],
                      "fixed_phase_normal_subspace": output["fixed_phase_normal_subspace"]}, indent=2))


if __name__ == "__main__":
    main()

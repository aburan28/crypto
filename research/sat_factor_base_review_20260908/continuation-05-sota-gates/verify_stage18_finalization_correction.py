#!/usr/bin/env python3
"""Fail-closed verifier for the additive Stage 18 finalization correction."""

from __future__ import annotations

import argparse
import copy
from datetime import datetime
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Any

import verify_stage18_degree23_panel as stage18


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PANEL = HERE / "stage-18-degree23-panel-lock-corrected-20260909"
AMENDMENT = HERE / "stage-18-amendment-02-finalization.json"
CORRECTED_SUMMARY = HERE / "stage-18-finalization-corrected-summary.json"
CERTIFICATE = HERE / "stage-18-finalization-certificate.json"
RESULTS = HERE / "STAGE18_RESULTS.md"
BUNDLE_MANIFEST = HERE / "stage-18-finalization-correction-manifest.json"
EXECUTION_COMMIT = "5d069cf659fe6f6f9badb817ef446d64f3ae2782"
AMENDMENT_FILE_SHA256 = "3197b0c34e56585c6f24ea32ce0d5aab2180ca9c578b1fec639c7054eda747c9"
AMENDMENT_CANONICAL_SHA256 = "95a08cf7ab3be6519119fdc74af3119b10120dfb93f452459d08af35171a2285"
PANEL_TREE_SHA256 = "f3b4da6dcf36366663c4c851b9b1deda06c43114376a14f9f55314bb4aa9dea5"
CORRECTED_STATUS = "complete_verified_scientific_task_panel_with_additive_finalization_correction"
RUNNER_RELATIVE = Path("research/sat_factor_base_review_20260908/continuation-05-sota-gates/run_stage18_degree23_panel.py")
VERIFIER_RELATIVE = Path("research/sat_factor_base_review_20260908/continuation-05-sota-gates/verify_stage18_degree23_panel.py")


class CorrectionError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise CorrectionError(message)


def _object_no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    value: dict[str, Any] = {}
    for key, item in pairs:
        if key in value:
            raise CorrectionError(f"duplicate JSON key: {key}")
        value[key] = item
    return value


def read_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(), object_pairs_hook=_object_no_duplicates)
    except CorrectionError:
        raise
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise CorrectionError(f"cannot read JSON {path}: {error}") from error
    require(isinstance(value, dict), f"expected JSON object: {path}")
    return value


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def output_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def panel_inventory(panel: Path) -> tuple[list[dict], int, str]:
    require(panel.is_dir() and not panel.is_symlink(), "immutable panel root is missing or non-regular")
    rows = []
    total_bytes = 0
    for path in sorted(panel.rglob("*")):
        if path.is_dir():
            require(not path.is_symlink(), f"symlink directory in immutable panel: {path}")
            continue
        require(path.is_file() and not path.is_symlink(), f"non-regular immutable panel artifact: {path}")
        relative = str(path.relative_to(panel))
        size = path.stat().st_size
        rows.append({"path": relative, "bytes": size, "sha256": sha256_file(path)})
        total_bytes += size
    return rows, total_bytes, hashlib.sha256(canonical_bytes(rows)).hexdigest()


def validate_amendment(amendment: dict, amendment_path: Path) -> None:
    require(sha256_file(amendment_path) == AMENDMENT_FILE_SHA256, "Amendment 02 bytes changed")
    require(canonical_sha256(amendment) == AMENDMENT_CANONICAL_SHA256, "Amendment 02 canonical hash changed")
    require(
        amendment.get("schema") == "koblitz_degree23_finalization_amendment.v1"
        and amendment.get("amendment_id") == "STAGE18-AMENDMENT-02-ADDITIVE-FINALIZATION"
        and amendment.get("execution_commit") == EXECUTION_COMMIT
        and amendment.get("corrected_status") == CORRECTED_STATUS,
        "Amendment 02 identity changed",
    )
    require(amendment.get("immutable_panel", {}).get("tree_sha256") == PANEL_TREE_SHA256, "panel tree binding changed")
    require(amendment.get("additive_outputs", {}).get("existing_panel_mutation_forbidden") is True, "panel mutation boundary changed")


def validate_panel_tree(panel: Path, amendment: dict) -> list[dict]:
    rows, total_bytes, tree_hash = panel_inventory(panel)
    frozen = amendment["immutable_panel"]
    require(len(rows) == frozen["file_count"] == 77, "immutable panel file count changed")
    require(total_bytes == frozen["total_bytes"] == 214914, "immutable panel byte count changed")
    require(tree_hash == frozen["tree_sha256"] == PANEL_TREE_SHA256, "immutable panel tree changed")
    critical = {
        "run.json": "run_sha256",
        "summary.json": "provisional_summary_sha256",
        "artifact-manifest.json": "artifact_manifest_sha256",
        "verification.json": "provisional_verification_sha256",
        "outer-attempts/0001/receipt.json": "outer_receipt_sha256",
        "outer-attempts/0001/stderr.txt": "outer_stderr_sha256",
    }
    indexed = {row["path"]: row for row in rows}
    for relative, field in critical.items():
        require(indexed[relative]["sha256"] == frozen[field], f"critical panel artifact changed: {relative}")
    return rows


def git_bytes(commit: str, relative: Path) -> bytes:
    completed = subprocess.run(
        ["git", "show", f"{commit}:{relative}"], cwd=REPO,
        check=False, capture_output=True,
    )
    require(completed.returncode == 0, f"cannot read historical Git file: {relative}")
    return completed.stdout


def git_blob_oid(commit: str, relative: Path) -> str:
    completed = subprocess.run(
        ["git", "rev-parse", f"{commit}:{relative}"], cwd=REPO,
        check=False, text=True, capture_output=True,
    )
    require(completed.returncode == 0, f"cannot read historical Git object: {relative}")
    return completed.stdout.strip()


def validate_historical_implementation(amendment: dict, run: dict) -> dict:
    observed = {}
    for key, relative in (("runner", RUNNER_RELATIVE), ("verifier", VERIFIER_RELATIVE)):
        data = git_bytes(EXECUTION_COMMIT, relative)
        identity = {
            "path": str(relative), "bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest(),
            "git_blob_oid": git_blob_oid(EXECUTION_COMMIT, relative),
        }
        require(identity == amendment["historical_implementation"][key], f"historical {key} amendment binding changed")
        recorded = run["implementation"][key]
        require(
            recorded.get("bytes") == identity["bytes"]
            and recorded.get("sha256") == identity["sha256"],
            f"run did not use the Git-authenticated historical {key}",
        )
        observed[key] = identity
    return observed


def historical_provisional_summary(panel: Path) -> dict:
    with tempfile.TemporaryDirectory(prefix="stage18-finalization-historical-") as temporary:
        clone = Path(temporary) / "repo"
        cloned = subprocess.run(
            ["git", "clone", "--quiet", "--no-hardlinks", str(REPO), str(clone)],
            check=False, text=True, capture_output=True,
        )
        require(cloned.returncode == 0, f"cannot clone historical checkout: {cloned.stderr.strip()}")
        checked = subprocess.run(
            ["git", "-C", str(clone), "checkout", "--quiet", "--detach", EXECUTION_COMMIT],
            check=False, text=True, capture_output=True,
        )
        require(checked.returncode == 0, f"cannot check out historical implementation: {checked.stderr.strip()}")
        historical_here = clone / RUNNER_RELATIVE.parent
        command = [
            str(Path(sys.executable).resolve()), str(historical_here / VERIFIER_RELATIVE.name),
            "--protocol", str(historical_here / "stage-18-degree23-panel-protocol.json"),
            "--panel", str(panel.resolve()), "--allow-incomplete",
        ]
        completed = subprocess.run(
            command, cwd=clone, check=False, text=True, capture_output=True,
            env={**dict(__import__("os").environ), "PYTHONDONTWRITEBYTECODE": "1"},
        )
        require(completed.returncode == 0, f"historical allow-incomplete verifier failed: {completed.stderr.strip()}")
        try:
            summary = json.loads(completed.stdout, object_pairs_hook=_object_no_duplicates)
        except (CorrectionError, json.JSONDecodeError) as error:
            raise CorrectionError(f"historical verifier emitted invalid JSON: {error}") from error
        require(isinstance(summary, dict), "historical verifier did not emit an object")
        return summary


def parse_timestamp(value: Any, label: str) -> datetime:
    require(isinstance(value, str) and value.endswith("Z"), f"invalid chronology timestamp: {label}")
    try:
        return datetime.fromisoformat(value[:-1] + "+00:00")
    except ValueError as error:
        raise CorrectionError(f"invalid chronology timestamp: {label}") from error


def validate_chronology_values(intervals: list[dict], outer_started_at: str,
                               outer_finished_at: str, post_execution_present: bool) -> dict:
    require(len(intervals) == 12, "chronology does not contain twelve task intervals")
    expected_order = [
        "discovery-a0", "discovery-a1",
        *(f"row-{index:02d}-{mode}" for index in range(1, 6) for mode in ("ic", "rho-auto")),
    ]
    require([row.get("task_id") for row in intervals] == expected_order, "task chronology order changed")
    outer_start = parse_timestamp(outer_started_at, "outer start")
    outer_finish = parse_timestamp(outer_finished_at, "outer finish")
    require(outer_start < outer_finish, "outer chronology is reversed")
    previous_finish = None
    parsed_finishes = []
    for position, interval in enumerate(intervals):
        require(
            isinstance(interval, dict)
            and set(interval) == {"task_id", "started_at", "finished_at"},
            f"invalid task chronology row: {position}",
        )
        started = parse_timestamp(interval["started_at"], f"{interval['task_id']} start")
        finished = parse_timestamp(interval["finished_at"], f"{interval['task_id']} finish")
        require(started < finished, f"task interval is reversed: {interval['task_id']}")
        require(outer_start < started and finished < outer_finish, f"task interval escapes outer interval: {interval['task_id']}")
        if previous_finish is not None:
            require(previous_finish < started, f"task intervals overlap or are out of order: {interval['task_id']}")
        previous_finish = finished
        parsed_finishes.append(finished)
    require(post_execution_present, "post-execution custody is absent after the task intervals")
    require(intervals[-1]["task_id"] == "row-05-rho-auto", "last task changed")
    require(
        intervals[-1]["finished_at"] == "2026-09-10T06:18:16.105512Z",
        "last task finish changed",
    )
    require(outer_started_at == "2026-09-10T06:05:40.384729Z", "outer start changed")
    require(outer_finished_at == "2026-09-10T06:18:16.907775Z", "outer finish changed")
    require(max(parsed_finishes) < outer_finish, "maximum task finish does not precede outer finish")
    return {
        "task_count": len(intervals),
        "task_intervals": intervals,
        "strictly_sequential": True,
        "outer_started_at": outer_started_at,
        "outer_finished_at": outer_finished_at,
        "last_task_id": intervals[-1]["task_id"],
        "last_task_finished_at": intervals[-1]["finished_at"],
        "maximum_task_finish_precedes_outer_finish": True,
        "post_execution_present": True,
        "chronology_sha256": canonical_sha256({
            "task_intervals": intervals,
            "outer_started_at": outer_started_at,
            "outer_finished_at": outer_finished_at,
            "post_execution_present": True,
        }),
    }


def validate_chronology(panel: Path, amendment: dict, run: dict) -> dict:
    expected = amendment.get("chronology", {})
    expected_intervals = expected.get("task_intervals")
    require(isinstance(expected_intervals, list), "Amendment 02 chronology is missing")
    actual = []
    for expected_row in expected_intervals:
        task_id = expected_row.get("task_id")
        invocation = read_json(panel / "tasks" / str(task_id) / "invocation.json")
        actual.append({
            "task_id": task_id,
            "started_at": invocation.get("started_at"),
            "finished_at": invocation.get("finished_at"),
        })
    require(actual == expected_intervals, "task invocation chronology differs from Amendment 02")
    outer = read_json(panel / "outer-attempts/0001/invocation.json")
    require(outer.get("started_at") == expected.get("outer_started_at"), "outer start differs from Amendment 02")
    require(outer.get("finished_at") == expected.get("outer_finished_at"), "outer finish differs from Amendment 02")
    require(expected.get("last_task_finished_at") == actual[-1]["finished_at"], "amendment last task finish changed")
    return validate_chronology_values(
        actual, outer["started_at"], outer["finished_at"],
        isinstance(run.get("post_execution"), dict),
    )


def validate_scientific_completion(panel: Path, amendment: dict, provisional: dict) -> dict:
    archived = read_json(panel / "summary.json")
    require(archived == provisional, "historical allow-incomplete recomputation differs from provisional summary")
    require(
        provisional.get("status") == "provisional_complete_pending_final_controls"
        and provisional.get("task_panel_complete") is True
        and provisional.get("expected_tasks") == provisional.get("attempted_tasks")
        == provisional.get("verified_tasks") == 12
        and provisional.get("status_counts") == {"verified": 12},
        "provisional summary does not contain a complete verified task panel",
    )
    run = read_json(panel / "run.json")
    require(run.get("source_revision", {}).get("commit") == EXECUTION_COMMIT, "execution commit changed")
    validate_historical_implementation(amendment, run)
    expected_task_ids = {
        "discovery-a0", "discovery-a1",
        *(f"row-{index:02d}-{mode}" for index in range(1, 6) for mode in ("ic", "rho-auto")),
    }
    require(set(run.get("tasks", {})) == expected_task_ids, "run task index changed")
    for task_id, entry in run["tasks"].items():
        receipt = panel / "tasks" / task_id / "receipt.json"
        require(entry.get("status") == "verified", f"task is not verified: {task_id}")
        require(entry.get("receipt_sha256") == sha256_file(receipt), f"task receipt hash changed: {task_id}")
    chronology = validate_chronology(panel, amendment, run)
    require(run.get("sources") == run.get("post_execution", {}).get("sources"), "source post-custody changed")
    require(run.get("tools") == run.get("post_execution", {}).get("binaries"), "binary post-custody changed")
    require(run.get("post_execution", {}).get("root_lock_sha256") == stage18.LOCK_SHA256, "post-execution lock changed")
    require(run.get("build_receipt_sha256") == sha256_file(panel / "build/receipt.json"), "build receipt hash changed")

    rows = provisional.get("rows")
    require(isinstance(rows, list) and len(rows) == 5, "five-row panel changed")
    for row, (index, secret, seed) in zip(rows, stage18.EXPECTED_RUNS):
        require(
            row.get("index") == index and row.get("secret") == secret
            and row.get("seed") == str(seed) and row.get("status") == "verified_pair"
            and row.get("target") == stage18.EXPECTED_TARGETS[index],
            f"frozen matched target row changed: {index}",
        )
    require(provisional.get("comparison") == "rho_faster_all_pairs", "pairwise comparison changed")
    accounting = provisional.get("accounting", {})
    require(
        accounting.get("online_ic_over_rho_core_ratio") is not None
        and accounting.get("strict_setup_charged_ic_over_rho_core_ratio") is not None,
        "algorithm ratios are missing",
    )

    outer_meter = read_json(panel / "outer-attempts/0001/metrics.json")
    outer_receipt = read_json(panel / "outer-attempts/0001/receipt.json")
    frozen_outer = amendment["observed_finalization_failure"]
    require(
        outer_meter.get("returncode") == outer_receipt.get("returncode") == frozen_outer["outer_returncode"] == 1
        and outer_meter.get("timed_out") is False and outer_meter.get("orphan_group_terminated") is False,
        "outer finalization failure terminal changed",
    )
    metrics = outer_meter.get("metrics", {})
    require(
        metrics.get("wall_seconds") == frozen_outer["outer_wall_seconds"]
        and metrics.get("total_core_seconds") == frozen_outer["outer_total_core_seconds"]
        and metrics.get("peak_rss_bytes") == frozen_outer["outer_peak_rss_bytes"],
        "outer finalization envelope changed",
    )
    require(
        (panel / "outer-attempts/0001/stderr.txt").read_text().strip()
        == frozen_outer["inner_error"],
        "outer finalization failure cause changed",
    )
    require(provisional.get("outer_driver_receipts", {}).get("excluded_from_algorithm_ratios") is True, "outer failure entered algorithm ratios")
    return {"run": run, "chronology": chronology}


def corrected_summary(amendment: dict, provisional: dict, chronology: dict) -> dict:
    corrected = copy.deepcopy(provisional)
    original_claim = corrected["claim_boundary"]
    corrected["schema"] = "koblitz_degree23_finalization_corrected_summary.v1"
    corrected["status"] = CORRECTED_STATUS
    corrected["original_provisional_status"] = provisional["status"]
    corrected["original_provisional_summary_sha256"] = amendment["immutable_panel"]["provisional_summary_sha256"]
    corrected["amendment_02_sha256"] = AMENDMENT_CANONICAL_SHA256
    corrected["scientific_claim_boundary"] = original_claim
    corrected["operational_finalization_failure"] = {
        **amendment["observed_finalization_failure"],
        "does_not_change_task_receipts": True,
        "does_not_change_algorithm_accounting": True,
    }
    corrected["verified_chronology"] = chronology
    corrected["claim_boundary"] = amendment["claim_boundary"]
    require(corrected["accounting"] == provisional["accounting"], "corrected summary changed algorithm accounting")
    return corrected


def certificate(amendment: dict, corrected: dict, implementation: dict,
                chronology: dict) -> dict:
    return {
        "schema": "koblitz_degree23_finalization_correction_certificate.v1",
        "status": "verified",
        "classification": CORRECTED_STATUS,
        "amendment_file_sha256": AMENDMENT_FILE_SHA256,
        "amendment_canonical_sha256": AMENDMENT_CANONICAL_SHA256,
        "execution_commit": EXECUTION_COMMIT,
        "immutable_panel": amendment["immutable_panel"],
        "historical_implementation": implementation,
        "verified_chronology": chronology,
        "corrected_summary_sha256": hashlib.sha256(output_bytes(corrected)).hexdigest(),
        "corrected_summary_canonical_sha256": canonical_sha256(corrected),
        "checks": [
            "77-file immutable panel tree",
            "historical runner and verifier from Git execution commit",
            "historical allow-incomplete provisional summary recomputation",
            "12 exact verified task receipts",
            "five exact matched target and scalar rows",
            "source, locked build, binary and post-execution custody",
            "outer returncode-1 finalization failure and full envelope",
            "outer envelope excluded from unchanged IC/rho ratios",
            "twelve exact sequential task intervals inside the outer interval",
            "post-execution custody present after the final task",
        ],
        "scientific_tasks_added_by_correction": 0,
        "existing_panel_files_modified_by_correction": 0,
        "outer_operational_failure": amendment["observed_finalization_failure"],
        "claim_boundary": amendment["claim_boundary"],
    }


def results_markdown(amendment: dict, corrected: dict) -> bytes:
    accounting = corrected["accounting"]
    outer = amendment["observed_finalization_failure"]
    chronology = corrected["verified_chronology"]
    text = f"""# Stage 18: five-target degree-23 panel

The scientific task panel is complete and verified through the additive finalization certificate. Both public discovery processes and all ten target processes have exact verified receipts. Every one of the five frozen scalar-derived targets was recovered independently by the index-calculus arm and its same-target signed-Frobenius rho control.

The original panel controls remain byte-identical and retain status `provisional_complete_pending_final_controls`. The inner driver reached final summary composition only after all twelve tasks and post-execution custody had completed, but it inspected `outer-attempts` while the active outer directory was still named `.staging-0001`. That raised `invalid or non-sequential outer attempt directory`; the outer meter retained return code 1 and then atomically renamed the directory to `0001` before writing the provisional controls.

Amendment 02 and its standalone verifier authenticate the immutable {amendment['immutable_panel']['file_count']}-file panel tree, load the execution-era runner and verifier from Git commit `{EXECUTION_COMMIT}`, reproduce the original summary with the historical verifier's `--allow-incomplete` path, and revalidate all task, source, build, lock, binary, post-execution, chronology, and outer receipts. The additive status is `{CORRECTED_STATUS}`. No panel file was changed and no scientific task was added or rerun.

The twelve task intervals are strictly sequential inside the outer interval. The outer interval ran from `{chronology['outer_started_at']}` through `{chronology['outer_finished_at']}`. The last task, `{chronology['last_task_id']}`, finished at `{chronology['last_task_finished_at']}`, before outer finalization failed.

The five index-calculus processes used {accounting['index_calculus_online']['total_core_seconds_sum']:.6f} core-seconds, while the five matched rho controls used {accounting['automorphism_rho']['total_core_seconds_sum']:.6f} core-seconds. The online aggregate IC/rho core-time ratio is {accounting['online_ic_over_rho_core_ratio']:.6f}. Charging the two public discovery processes once raises the strict setup-charged ratio to {accounting['strict_setup_charged_ic_over_rho_core_ratio']:.6f}. Rho was faster in every matched pair. These ratios use the unchanged task receipts.

The failed outer envelope remains charged separately: {outer['outer_total_core_seconds']:.6f} core-seconds, {outer['outer_wall_seconds']:.6f} wall-seconds, and {outer['outer_peak_rss_bytes']} peak RSS bytes. It includes nested work and is excluded from the IC/rho algorithm ratios; it must not be added to the separately summed task costs.

This is repeated scalar-blind evidence at one toy degree and one public factor base. It does not establish an asymptotic crossover, a cryptographic-size attack, novelty, independent external reproduction, or a state-of-the-art result.
"""
    return text.encode()


def _manifest_row(path: Path, data: bytes) -> dict:
    return {
        "path": str(path.resolve().relative_to(REPO.resolve())),
        "bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
    }


def bundle_manifest(amendment: dict, corrected: dict, cert: dict,
                    results: bytes) -> dict:
    component_data = {
        AMENDMENT: AMENDMENT.read_bytes(),
        HERE / "compose_stage18_finalization_correction.py": (HERE / "compose_stage18_finalization_correction.py").read_bytes(),
        HERE / "verify_stage18_finalization_correction.py": (HERE / "verify_stage18_finalization_correction.py").read_bytes(),
        HERE / "verify_stage18_degree23_panel.py": (HERE / "verify_stage18_degree23_panel.py").read_bytes(),
        CORRECTED_SUMMARY: output_bytes(corrected),
        CERTIFICATE: output_bytes(cert),
        RESULTS: results,
    }
    rows = sorted((_manifest_row(path, data) for path, data in component_data.items()), key=lambda row: row["path"])
    return {
        "schema": "koblitz_degree23_finalization_correction_bundle_manifest.v1",
        "status": "verified",
        "classification": CORRECTED_STATUS,
        "self_excluded": True,
        "manifest_path": str(BUNDLE_MANIFEST.resolve().relative_to(REPO.resolve())),
        "immutable_panel": amendment["immutable_panel"],
        "files": rows,
        "file_count": len(rows),
        "scientific_tasks_added_by_bundle": 0,
        "immutable_panel_files_modified_by_bundle": 0,
    }


def compose_expected(panel: Path = PANEL, amendment_path: Path = AMENDMENT) -> tuple[dict, dict]:
    amendment = read_json(amendment_path)
    validate_amendment(amendment, amendment_path)
    validate_panel_tree(panel, amendment)
    provisional = historical_provisional_summary(panel)
    completion = validate_scientific_completion(panel, amendment, provisional)
    run = completion["run"]
    chronology = completion["chronology"]
    implementation = validate_historical_implementation(amendment, run)
    summary = corrected_summary(amendment, provisional, chronology)
    return summary, certificate(amendment, summary, implementation, chronology)


def compose_bundle_expected(panel: Path = PANEL, amendment_path: Path = AMENDMENT) -> tuple[dict, dict, bytes, dict]:
    amendment = read_json(amendment_path)
    summary, cert = compose_expected(panel, amendment_path)
    results = results_markdown(amendment, summary)
    manifest = bundle_manifest(amendment, summary, cert, results)
    return summary, cert, results, manifest


def verify_outputs(panel: Path = PANEL, amendment_path: Path = AMENDMENT,
                   summary_path: Path = CORRECTED_SUMMARY,
                   certificate_path: Path = CERTIFICATE,
                   results_path: Path = RESULTS,
                   manifest_path: Path = BUNDLE_MANIFEST) -> dict:
    expected_summary, expected_certificate, expected_results, expected_manifest = compose_bundle_expected(panel, amendment_path)
    require(read_json(summary_path) == expected_summary, "additive corrected summary differs from recomputation")
    require(read_json(certificate_path) == expected_certificate, "additive finalization certificate differs from recomputation")
    require(results_path.read_bytes() == expected_results, "authenticated Stage 18 results Markdown differs from recomputation")
    require(read_json(manifest_path) == expected_manifest, "correction bundle manifest differs from recomputation")
    require(
        sha256_file(summary_path) == expected_certificate["corrected_summary_sha256"],
        "corrected summary byte hash changed",
    )
    return expected_certificate


def self_test() -> dict:
    summary, cert, results, manifest = compose_bundle_expected()
    require(summary["status"] == CORRECTED_STATUS and summary["verified_tasks"] == 12, "corrected completion failed")
    checks = 8
    changed = copy.deepcopy(read_json(AMENDMENT))
    changed["immutable_panel"]["tree_sha256"] = "0" * 64
    try:
        validate_amendment(changed, AMENDMENT)
    except CorrectionError:
        checks += 1
    else:
        raise AssertionError("mutated Amendment 02 was accepted")
    provisional = read_json(PANEL / "summary.json")
    bad_provisional = copy.deepcopy(provisional)
    bad_provisional["verified_tasks"] = 11
    try:
        validate_scientific_completion(PANEL, read_json(AMENDMENT), bad_provisional)
    except CorrectionError:
        checks += 1
    else:
        raise AssertionError("incomplete provisional task panel was accepted")
    bad_run = read_json(PANEL / "run.json")
    bad_run["implementation"]["verifier"]["sha256"] = "0" * 64
    try:
        validate_historical_implementation(read_json(AMENDMENT), bad_run)
    except CorrectionError:
        checks += 1
    else:
        raise AssertionError("wrong historical verifier identity was accepted")
    bad_intervals = copy.deepcopy(read_json(AMENDMENT)["chronology"]["task_intervals"])
    bad_intervals[1]["started_at"] = bad_intervals[0]["finished_at"]
    try:
        validate_chronology_values(
            bad_intervals, "2026-09-10T06:05:40.384729Z",
            "2026-09-10T06:18:16.907775Z", True,
        )
    except CorrectionError:
        checks += 1
    else:
        raise AssertionError("out-of-order task chronology was accepted")
    with tempfile.TemporaryDirectory(prefix="stage18-finalization-tree-control-") as temporary:
        copied_panel = Path(temporary) / "panel"
        shutil.copytree(PANEL, copied_panel)
        receipt = copied_panel / "tasks/row-01-ic/receipt.json"
        receipt.write_bytes(receipt.read_bytes() + b" ")
        try:
            validate_panel_tree(copied_panel, read_json(AMENDMENT))
        except CorrectionError:
            checks += 1
        else:
            raise AssertionError("mutated task receipt tree was accepted")
    require(results == results_markdown(read_json(AMENDMENT), summary), "deterministic results rendering changed")
    require(manifest["file_count"] == 7 and manifest["self_excluded"] is True, "bundle manifest inventory changed")
    checks += 2
    if all(path.is_file() for path in (CORRECTED_SUMMARY, CERTIFICATE, RESULTS, BUNDLE_MANIFEST)):
        verify_outputs()
        checks += 1
    require(not CORRECTED_SUMMARY.resolve().is_relative_to(PANEL.resolve()), "corrected summary is inside immutable panel")
    require(not CERTIFICATE.resolve().is_relative_to(PANEL.resolve()), "certificate is inside immutable panel")
    checks += 2
    return {
        "self_test": "pass", "checks": checks,
        "status": summary["status"], "verified_tasks": summary["verified_tasks"],
        "matched_pairs": len(summary["rows"]), "panel_tree_sha256": PANEL_TREE_SHA256,
        "outer_returncode": cert["outer_operational_failure"]["outer_returncode"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, default=PANEL)
    parser.add_argument("--amendment", type=Path, default=AMENDMENT)
    parser.add_argument("--summary", type=Path, default=CORRECTED_SUMMARY)
    parser.add_argument("--certificate", type=Path, default=CERTIFICATE)
    parser.add_argument("--results", type=Path, default=RESULTS)
    parser.add_argument("--manifest", type=Path, default=BUNDLE_MANIFEST)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    try:
        result = self_test() if args.self_test else verify_outputs(
            args.panel.resolve(), args.amendment.resolve(),
            args.summary.resolve(), args.certificate.resolve(),
            args.results.resolve(), args.manifest.resolve(),
        )
        print(json.dumps(result, indent=2, sort_keys=True))
    except (CorrectionError, OSError, subprocess.SubprocessError) as error:
        parser.exit(1, f"Stage 18 finalization verification failed: {error}\n")


if __name__ == "__main__":
    main()

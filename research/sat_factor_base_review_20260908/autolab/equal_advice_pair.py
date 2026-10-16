"""Equal-advice BL archive reuse + optional compact-orbit IC pairing.

Honest advance toward schema-complete vs_rho for
``koblitz.compact_orbit.n61_bl_reference``: reuse a hash-verified BL
L>=1024 archive, measure IC online on the same public points, and emit
only measured fields. Does not invent manifests, replay certificates,
rho_policy, or same_resource_envelope when arms are not co-measured.
"""
from __future__ import annotations

import json
import math
import os
import shutil
import statistics
import subprocess
from pathlib import Path
from typing import Any

IC_PHASE_FIELDS = {
    "T_target_query_ms": "target_query_ms",
    "T_target_PDP_ms": "target_pdp_ms",
    "T_target_relation_check_ms": "target_relation_check_ms",
    "T_target_descent_ms": "target_descent_ms",
    "T_target_recovery_check_ms": "target_recovery_check_ms",
}


def resolve_bl_archive(path: Path) -> Path:
    """Accept the archive root or the nested run directory."""
    path = path.resolve()
    if (path / "logs" / "bl.jsonl").is_file() and (path / "inputs" / "points.txt").is_file():
        return path
    nested = sorted(
        p for p in path.iterdir() if p.is_dir() and (p / "logs" / "bl.jsonl").is_file()
    ) if path.is_dir() else []
    if len(nested) == 1:
        return nested[0]
    raise ValueError(
        f"BL archive must contain logs/bl.jsonl and inputs/points.txt: {path}"
    )


def parse_points_file(path: Path) -> list[list[int]]:
    points: list[list[int]] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        value = json.loads(line)
        if not (isinstance(value, list) and len(value) == 2 and all(type(v) is int for v in value)):
            raise ValueError(f"point line must be [x,y] ints: {line[:80]}")
        points.append(value)
    return points


def parse_scalars_file(path: Path) -> list[int]:
    scalars: list[int] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if line:
            scalars.append(int(line))
    return scalars


def materialize_bl_archive(
    *,
    archive: Path,
    run: Path,
    lab: Any,
    targets: int,
    table_bytes: int,
) -> dict[str, Any]:
    """Copy archived BL inputs/logs into the run and hash-verify counts."""
    archive = resolve_bl_archive(archive)
    for name in ("artifacts", "inputs", "logs", "receipts"):
        (run / name).mkdir(parents=True, exist_ok=True)

    copies = {
        "inputs/points.txt": archive / "inputs" / "points.txt",
        "inputs/scalars.txt": archive / "inputs" / "scalars.txt",
        "logs/bl.jsonl": archive / "logs" / "bl.jsonl",
    }
    for stdout_name in ("bl.stdout", "bl.stderr", "targets.stdout", "targets.stderr"):
        src = archive / "logs" / stdout_name
        if src.is_file():
            copies[f"logs/{stdout_name}"] = src
    for relative, src in copies.items():
        require_file = src.is_file()
        if not require_file:
            raise FileNotFoundError(f"missing archive file: {src}")
        dest = run / relative
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)

    points = parse_points_file(run / "inputs/points.txt")
    scalars = parse_scalars_file(run / "inputs/scalars.txt")
    if len(points) != targets:
        raise ValueError(f"archive points {len(points)} != requested targets {targets}")
    if len(scalars) != targets:
        raise ValueError(f"archive scalars {len(scalars)} != requested targets {targets}")

    records = lab.parse_json_lines((run / "logs/bl.jsonl").read_text())
    precompute = next((r for r in records if r.get("kind") == "bl_precompute"), None)
    summary = next(
        (r for r in records if r.get("kind") in ("bl_summary", "bl_level_summary")),
        None,
    )
    online_rows = [r for r in records if r.get("kind") == "bl_online_target"]
    primary = [
        r
        for r in online_rows
        if int(r.get("table_bytes") or 0) == int(table_bytes)
    ]
    if len(primary) != targets:
        raise ValueError(
            f"archive primary-table online rows {len(primary)} != targets {targets} "
            f"(table_bytes={table_bytes})"
        )
    pin = {
        "kind": "reused_bl_archive",
        "archive_run": str(archive),
        "archive_run_id": archive.name,
        "points_sha256": lab.sha256(run / "inputs/points.txt"),
        "scalars_sha256": lab.sha256(run / "inputs/scalars.txt"),
        "bl_jsonl_sha256": lab.sha256(run / "logs/bl.jsonl"),
        "targets": targets,
        "table_bytes": table_bytes,
        "bl_online_rows": len(online_rows),
        "bl_primary_online_rows": len(primary),
        "bl_precompute": precompute,
        "bl_summary": summary,
        "remeasured": False,
        "note": (
            "BL arm reused byte-for-byte from archive; not co-measured with the "
            "IC arm in this process, so same_resource_envelope stays unset."
        ),
    }
    lab.write_json(run / "artifacts/bl_archive_pin.json", pin)
    return pin


def select_bl_primary_rows(
    records: list[dict[str, Any]], *, table_bytes: int
) -> list[dict[str, Any]]:
    rows = [
        r
        for r in records
        if r.get("kind") == "bl_online_target" and int(r.get("table_bytes") or 0) == int(table_bytes)
    ]
    rows.sort(key=lambda r: int(r["target_index"]))
    return rows


def write_ic_points_subset(points: list[list[int]], path: Path, count: int) -> None:
    path.write_text("".join(f"[{p[0]},{p[1]}]\n" for p in points[:count]))


def run_ic_arm(
    *,
    lab: Any,
    run: Path,
    binaries: dict[str, str],
    n: int,
    a: int,
    k: int,
    rank_seed: int,
    points: list[list[int]],
    ic_target_count: int,
) -> dict[str, Any]:
    count = int(ic_target_count)
    if count < 1 or count > len(points):
        raise ValueError(f"ic_target_count must be in 1..{len(points)}")
    points_path = run / "inputs/ic_points.txt"
    write_ic_points_subset(points, points_path, count)
    out_path = run / "logs/ic.jsonl"
    stdout_path = run / "logs/ic.stdout"
    stderr_path = run / "logs/ic.stderr"
    cmd = [
        binaries["ic"],
        f"construct:{n}:{a}:{k}",
        str(points_path),
        str(rank_seed),
        str(out_path),
    ]
    env = os.environ.copy()
    env["RAYON_NUM_THREADS"] = "1"
    # Stream to files: capture_output can OOM / appear hung on L>=1024 panels.
    with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
        completed = subprocess.run(cmd, cwd=lab.REPO, stdout=stdout, stderr=stderr, env=env)
    if completed.returncode != 0:
        err_tail = stderr_path.read_text()[-4000:] if stderr_path.is_file() else ""
        raise RuntimeError(f"IC producer failed:\n{err_tail}")
    if not out_path.is_file() or out_path.stat().st_size == 0:
        raise RuntimeError("IC producer did not write jsonl")
    rows = lab.parse_json_lines(out_path.read_text())
    targets = [r for r in rows if r.get("kind") == "compact_orbit_dlp_target"]
    summaries = [r for r in rows if r.get("kind") == "compact_orbit_dlp_summary"]
    if len(targets) != count:
        raise RuntimeError(f"IC online rows {len(targets)} != ic_target_count {count}")
    return {
        "command": cmd,
        "exit_code": completed.returncode,
        "ic_target_count": count,
        "k": k,
        "rank_seed": rank_seed,
        "points_sha256": lab.sha256(points_path),
        "ic_jsonl_sha256": lab.sha256(out_path),
        "target_rows": targets,
        "summary": summaries[-1] if summaries else None,
        "stdout_tail": stdout_path.read_text()[-2000:] if stdout_path.is_file() else "",
    }


def verify_scalar(recovered: Any, truth: int) -> bool:
    try:
        return int(recovered) == int(truth)
    except (TypeError, ValueError):
        return False


def draft_paired_vs_rho_slice(
    *,
    lab: Any,
    beat: dict[str, Any],
    beat_id: str,
    run: Path,
    pin: dict[str, Any],
    ic_obs: dict[str, Any],
    bl_records: list[dict[str, Any]],
    points: list[list[int]],
    scalars: list[int],
    table_bytes: int,
    binaries: dict[str, str],
    git_head: str,
) -> dict[str, Any]:
    """Build a one-target schema-facing claim from measured IC+BL fields only."""
    bl_primary = select_bl_primary_rows(bl_records, table_bytes=table_bytes)
    ic_row = ic_obs["target_rows"][0]
    bl_row = bl_primary[0]
    point = points[0]
    truth = scalars[0]
    ic_q = ic_row.get("published_q") or ic_row.get("target")
    if not (isinstance(ic_q, list) and len(ic_q) == 2):
        ic_q = point
    same_point = (
        isinstance(ic_q, list)
        and [int(ic_q[0]), int(ic_q[1])] == [int(point[0]), int(point[1])]
    )
    ic_ms = lab.positive_cost(ic_row.get("online_ms"))
    bl_ms = lab.positive_cost(bl_row.get("online_ms"))
    ic_phases = {}
    for schema_key, field in IC_PHASE_FIELDS.items():
        value = ic_row.get(field)
        if type(value) in (int, float) and math.isfinite(value) and value >= 0:
            ic_phases[schema_key] = float(value)
    ic_phase_sum_ok = (
        len(ic_phases) == len(IC_PHASE_FIELDS)
        and ic_ms is not None
        and math.isclose(math.fsum(ic_phases.values()), float(ic_ms), rel_tol=1e-8, abs_tol=0.02)
    )
    # BL emits one online interval without exclusive sub-phases; record that
    # single measured interval rather than inventing walk/recovery splits.
    bl_phases = {"bl_online": float(bl_ms)} if bl_ms is not None else None
    ic_verified = bool(
        ic_row.get("group_verified") is True
        and verify_scalar(ic_row.get("recovered_scalar"), truth)
        and same_point
        and ic_phase_sum_ok
    )
    bl_verified = bool(
        (bl_row.get("solved") or bl_row.get("verified") or bl_row.get("verified_vs_truth"))
        and verify_scalar(bl_row.get("recovered_scalar"), truth)
    )
    speedup = (
        float(bl_ms) / float(ic_ms)
        if ic_ms is not None and bl_ms is not None and ic_verified and bl_verified
        else None
    )
    ic_hash = lab.public_target_hash(beat, point)
    probes = ic_row.get("probes")
    steps = bl_row.get("online_steps")
    operation_accounting = None
    if all(
        isinstance(value, (int, float)) and not isinstance(value, bool) and value > 0
        for value in (probes, steps)
    ):
        operation_accounting = {
            "schema_version": "1.0",
            "unit_assumption": (
                "IC root-index probes and BL online walk steps are distinct "
                "uncalibrated counters"
            ),
            "comparison_status": "native_counters_only",
            "operation_units": {"ic": "root-index probes", "rho": "BL online walk steps"},
            "ic_online_operations": probes,
            "ic_online_operations_basis": "probes of the one online target",
            "rho_online_operations": steps,
            "rho_online_operations_basis": "measured BL online_steps at primary table_bytes",
            "rho_per_ic_native_counter": float(steps) / float(probes),
            "ops_speedup_online": None,
        }
    missing = [
        "candidate_id / candidate_manifest_sha256 (IC1 identity manifests not minted)",
        "workload_id / workload_manifest_sha256",
        "run_id (IC1W…R… form)",
        "rho_policy (BL equal-advice is not KS batch rho; policy object not minted)",
        "same_resource_envelope (BL reused from archive; not co-measured)",
        "independent_validation / replay certificates",
        "ic_resource_envelope / rho_resource_envelope matched pair",
        "rho_online_phase_ms declared KS phases (BL emits one online interval only)",
        "timing_class single_target_online primary (beat goal is single_target_online_wall)",
        "record_class from vs_rho enum",
    ]
    claim = {
        "schema_version": 2,
        "task_id": lab.TASK_ID,
        "beat_id": beat_id,
        "stage": "vs_rho",
        "status": "SCHEMA_INCOMPLETE",
        "diagnostic_only": True,
        "reference_arm": "bernstein_lange_equal_advice",
        "regime": beat.get("regime", "koblitz"),
        "n": beat["n"],
        "a": beat["a"],
        "n_or_bits": beat["n"],
        "timing_class": beat.get("timing_class_goal", "single_target_online_wall"),
        "record_class": None,
        "candidate_id": None,
        "candidate_manifest_sha256": None,
        "workload_id": None,
        "workload_manifest_sha256": None,
        "run_id": None,
        "autolab_run_id": run.name,
        "target_count": 1,
        "panel_target_count": len(points),
        "paired_ic_target_count": ic_obs["ic_target_count"],
        "ic_K": ic_obs["k"],
        "ic_target_hash": ic_hash,
        "rho_target_hash": ic_hash if same_point else None,
        "ic_online_wall_ms": ic_ms,
        "rho_online_wall_ms": bl_ms,
        "ic_online_ms": ic_ms,
        "rho_online_ms": bl_ms,
        "ic_cost": ic_ms,
        "rho_cost": bl_ms,
        "online_speedup": speedup,
        "online_interval": {
            "ic_start_event": ic_row.get("online_start_event"),
            "ic_stop_event": ic_row.get("online_stop_event"),
            "ic_start_ns": ic_row.get("online_start_ns"),
            "ic_stop_ns": ic_row.get("online_stop_ns"),
            "bl_note": (
                "BL producer records online_ms without start/stop event markers; "
                "archived arm reused, not co-measured"
            ),
        },
        "same_resource_envelope": None,
        "independent_validation": False,
        "ic_replay_certificate_sha256": None,
        "rho_replay_certificate_sha256": None,
        "ic_resource_envelope": None,
        "rho_resource_envelope": None,
        "rho_policy": None,
        "ic_scalar_verified": ic_verified,
        "rho_scalar_verified": bl_verified,
        "ic_verified": ic_verified,
        "rho_verified": bl_verified,
        "automorphism_discount": lab.automorphism_discount(int(beat["n"])),
        "all_stages_charged_same_series": False,
        "paired_target": {
            "ic_public_q": [int(point[0]), int(point[1])],
            "rho_public_q": [int(point[0]), int(point[1])],
            "same_public_point": same_point,
            "ic_fixture_scalar": truth,
            "rho_fixture_scalar": truth,
            "same_fixture_scalar": True,
            "target_index": 0,
            "table_bytes": table_bytes,
        },
        "ic_online_phase_ms": ic_phases or None,
        "rho_online_phase_ms": bl_phases,
        "ic_online_interval": (
            f"{ic_row.get('online_start_event')} -> {ic_row.get('online_stop_event')}"
            if ic_row.get("online_start_event")
            else None
        ),
        "rho_online_interval": (
            "BL online_ms (walk through recovery); no exclusive sub-phase markers"
        ),
        "operation_accounting": operation_accounting,
        "bl_archive_pin": {
            "archive_run_id": pin.get("archive_run_id"),
            "points_sha256": pin.get("points_sha256"),
            "bl_jsonl_sha256": pin.get("bl_jsonl_sha256"),
            "remeasured": pin.get("remeasured"),
        },
        "bl_online_row": {
            "online_ms": bl_row.get("online_ms"),
            "online_steps": bl_row.get("online_steps"),
            "table_bytes": bl_row.get("table_bytes"),
            "solved": bl_row.get("solved"),
            "recovered_scalar": bl_row.get("recovered_scalar"),
        },
        "ic_online_row": {
            "online_ms": ic_row.get("online_ms"),
            "probes": ic_row.get("probes"),
            "group_verified": ic_row.get("group_verified"),
            "recovered_scalar": ic_row.get("recovered_scalar"),
            "K": ic_obs["k"],
        },
        "fixture_hash": ic_hash if same_point else None,
        "executable_or_source_hash": {
            "ic": lab.sha256(binaries["ic"]),
            "bl_source": lab.sha256(lab.REPO / beat["panel_producers"]["bl"]["source"]),
            "ic_source": lab.sha256(lab.REPO / beat["panel_producers"]["ic"]["source"]),
            "git_head": git_head,
        },
        "host_id": lab.host_record(),
        "resource_caps": {"rayon_threads": 1, "common_cap_bytes": beat.get("resource_cap_bytes")},
        "seeds": {
            "ic_rank_seed": ic_obs["rank_seed"],
            "target_label": beat.get("target_label_default"),
            "archive_run_id": pin.get("archive_run_id"),
        },
        "claim_boundary": (
            "Diagnostic equal-advice BL reference with paired compact-orbit IC online "
            "on the same public synthetic points. BL arm reused from archive (not "
            "co-measured). Fail-closed for vs_rho ledger promotion; not KS batch rho; "
            "not asymptotic sub-rho; PENDING_INDEPENDENT_VALIDATION at most."
        ),
        "claim_boundary_non_claims": [
            "not key recovery",
            "not asymptotic sub-sqrt",
            "not imported/external points",
            "not a KS v2 batch rho comparator",
            "not same_resource_envelope (archived BL + later IC)",
            "not ledger promotion until schema-complete co-measured pair and independent-host replay",
        ],
        "missing_vs_rho_fields_note": "; ".join(missing),
        "remaining_schema_gaps": missing,
        "independent_replay_pointer": None,
        "verdict": "DIAGNOSTIC_PAIRED_IC_SLICE_SCHEMA_INCOMPLETE",
        "existing_study": beat.get("existing_study"),
    }
    return claim


def panel_pairing_summary(
    *,
    ic_obs: dict[str, Any],
    bl_records: list[dict[str, Any]],
    points: list[list[int]],
    scalars: list[int],
    table_bytes: int,
) -> dict[str, Any]:
    bl_primary = select_bl_primary_rows(bl_records, table_bytes=table_bytes)
    count = min(len(ic_obs["target_rows"]), len(bl_primary), len(points), len(scalars))
    ratios: list[float] = []
    ic_ms_list: list[float] = []
    bl_ms_list: list[float] = []
    verified = 0
    for i in range(count):
        ic_row = ic_obs["target_rows"][i]
        bl_row = bl_primary[i]
        truth = scalars[i]
        ic_ok = ic_row.get("group_verified") is True and verify_scalar(
            ic_row.get("recovered_scalar"), truth
        )
        bl_ok = bool(bl_row.get("solved")) and verify_scalar(bl_row.get("recovered_scalar"), truth)
        ic_ms = ic_row.get("online_ms")
        bl_ms = bl_row.get("online_ms")
        if isinstance(ic_ms, (int, float)) and isinstance(bl_ms, (int, float)) and ic_ms > 0:
            ic_ms_list.append(float(ic_ms))
            bl_ms_list.append(float(bl_ms))
            if ic_ok and bl_ok:
                verified += 1
                ratios.append(float(bl_ms) / float(ic_ms))
    def describe(values: list[float]) -> dict[str, Any] | None:
        if not values:
            return None
        return {
            "count": len(values),
            "median": statistics.median(values),
            "mean": statistics.fmean(values),
            "min": min(values),
            "max": max(values),
        }

    return {
        "kind": "equal_advice_ic_bl_panel_pairing",
        "paired_targets": count,
        "verified_pairs": verified,
        "table_bytes": table_bytes,
        "ic_K": ic_obs["k"],
        "ic_online_ms": describe(ic_ms_list),
        "bl_online_ms": describe(bl_ms_list),
        "online_speedup_bl_over_ic": describe(ratios),
        "note": (
            "Secondary multi-target pairing stats; primary vs_rho schema still "
            "requires target_count=1 with co-measured envelopes and identity manifests."
        ),
    }

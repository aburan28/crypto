#!/usr/bin/env python3
"""Compose the exact-commit selected-default replay."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


STAGE = Path(__file__).resolve().parent
EXPECTED_EQUATIONS = "02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def receipt(path: Path) -> dict:
    return {"path": str(path.relative_to(STAGE)), "bytes": path.stat().st_size, "sha256": sha256(path)}


def process_artifact(root: Path) -> dict:
    metrics = root / "metrics.json"
    stdout = root / "stdout.txt"
    stderr = root / "stderr.txt"
    return {
        "process": json.loads(metrics.read_text()),
        "artifacts": {
            "metrics": receipt(metrics),
            "stdout": receipt(stdout),
            "stderr": receipt(stderr),
        },
    }


failed_build = process_artifact(STAGE / "development" / "failed-build-locked")
build = process_artifact(STAGE / "development" / "build")
tests = {}
for name, expected in (
    ("default_f4", "10 passed; 0 failed"),
    ("quadratic_f4", "10 passed; 0 failed"),
    ("default_backend", "3 passed; 0 failed"),
    ("quadratic_backend", "3 passed; 0 failed"),
):
    root = STAGE / "development" / f"test-{name.replace('_', '-')}"
    item = process_artifact(root)
    item["passed"] = (
        item["process"]["returncode"] == 0
        and not item["process"]["timed_out"]
        and expected in (root / "stdout.txt").read_text()
    )
    tests[name] = item
if not all(item["passed"] for item in tests.values()):
    raise SystemExit("post-selection validation failed")

replay_root = STAGE / "development" / "replay"


def load_replay(backend: str) -> dict:
    root = replay_root / backend
    metrics_path = root / "metrics.json"
    stdout_path = root / "stdout.json"
    stderr_path = root / "stderr.txt"
    return {
        "process": json.loads(metrics_path.read_text()),
        "report": json.loads(stdout_path.read_text()),
        "artifacts": {
            "metrics": receipt(metrics_path),
            "stdout": receipt(stdout_path),
            "stderr": receipt(stderr_path),
        },
    }


native = load_replay("native-f4")
direct = load_replay("direct-mitm")
extra = native["report"]["cost"]["extra"]
native_correct = (
    native["process"]["returncode"] == 0
    and not native["process"]["timed_out"]
    and native["report"]["status"] == "unsat"
    and native["report"]["exhaustive"] is True
    and native["report"]["source_instance_verified"] is True
    and native["report"]["regenerated_source_exact"] is True
    and native["report"]["solver_equations_blake3"] == EXPECTED_EQUATIONS
    and native["report"]["cost"]["ops"] == 319_313_687_585
    and extra["word_xors_performed"] == 147_794_583_858
    and extra["pair_dense_select_calls"] == 1_011_275
    and extra["pair_quadratic_select_calls"] == 0
    and extra["pair_candidate_visits"] == 1_137_001_812
    and extra["pair_lcm_groups"] == 594_604_504
    and extra["pair_cover_lookups"] == 698_657_372
    and extra["pair_dense_scratch_bytes_max"] == 4_472_832
    and extra["full_m4ri_matrices"] == 0
    and native["report"]["single_thread_requested"] is False
    and native["report"]["conflicts"] is None
)
direct_correct = (
    direct["process"]["returncode"] == 0
    and not direct["process"]["timed_out"]
    and direct["report"]["status"] == "unsat"
    and direct["report"]["exhaustive"] is True
    and direct["report"]["source_instance_verified"] is True
)
if not native_correct or not direct_correct:
    raise SystemExit("selected replay correctness failure")

all_processes = [
    failed_build["process"],
    build["process"],
    *(item["process"] for item in tests.values()),
    native["process"],
    direct["process"],
]
charge = {
    "components": len(all_processes),
    "wall_seconds_sum": sum(item["metrics"]["wall_seconds"] for item in all_processes),
    "total_core_seconds_sum": sum(item["metrics"]["total_core_seconds"] for item in all_processes),
    "peak_rss_bytes_max": max(item["metrics"]["peak_rss_bytes"] for item in all_processes),
}
nm = native["process"]["metrics"]
dm = direct["process"]["metrics"]
result = {
    "schema": "koblitz_stage185_dense_pair_default_replay.v1",
    "claim_boundary": "Exact selected-default implementation replay on one opened n=59 target; not a full index-calculus run or SOTA.",
    "source": {
        "selection_commit": "8014149a2b55cd2cca202237a302643e39a50f6e",
        "binary_sha256": "3aefd5cd0daf73602bf8a38e7fdcb363e630ffd57bfebe47a37123b1d03429b1",
        "source_checkout": "/Volumes/SSD990/crypto-kic-stage185-source",
        "supplied_lock_sha256": "4f17b356fa7bac392b6d801d1c74fb9e36b6517f9465c8ebc19bb9a2792a84c5",
        "lock_tracked_at_selection_commit": False,
    },
    "failed_locked_build": failed_build,
    "build": build,
    "tests": tests,
    "replay": {
        "default_native_f4": native,
        "direct_mitm": direct,
        "default_native_correct": native_correct,
        "direct_correct": direct_correct,
        "valid_single_core_seconds": None,
        "single_core_reason": "The backend requested 12 Rayon workers; process-meter single_core_seconds is not a valid single-core row for this replay.",
    },
    "same_binary_native_over_direct": {
        "wall": nm["wall_seconds"] / dm["wall_seconds"],
        "core": nm["total_core_seconds"] / dm["total_core_seconds"],
        "rss": nm["peak_rss_bytes"] / dm["peak_rss_bytes"],
    },
    "decision": {
        "status": "SELECTED_DEFAULT_REPLAY_PASS",
        "dense_default_verified": True,
        "quadratic_control_verified": True,
        "lockfile_restore_required": True,
        "reason": "Default routing, control routing, tests, and exact target replay all passed; restore the supplied lockfile to make clean --locked builds self-contained.",
    },
    "campaign_charge": charge,
    "conflicts": None,
    "wall_interpretation": "The host was heavily descheduled during the replay; this run validates routing and correctness, while Stage 184 supplies paired performance evidence.",
}
(STAGE / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

lines = [
    "# Stage 185: selected dense-pair default replay",
    "",
    "An exact detached checkout of selection commit `8014149a2` passed ten F4 tests and three backend tests in both default-dense and explicit quadratic-control modes. With no selector environment variable, the target replay reported 1,011,275 dense selections, zero quadratic selections, zero full-M4RI matrices, and exact exhaustive UNSAT.",
    "",
    f"The load-affected default replay took {nm['wall_seconds']:.6f} wall seconds, {nm['total_core_seconds']:.6f} core-seconds, and {nm['peak_rss_bytes']} bytes RSS. Same-binary direct MITM took {dm['wall_seconds']:.6f} wall seconds, {dm['total_core_seconds']:.6f} core-seconds, and {dm['peak_rss_bytes']} bytes RSS. The resulting stage ratios are {result['same_binary_native_over_direct']['wall']:.2f}x wall, {result['same_binary_native_over_direct']['core']:.2f}x CPU, and {result['same_binary_native_over_direct']['rss']:.2f}x RSS; they are not a full-method comparison.",
    "",
    "The first exact `--locked` build failed because the selection commit did not track `Cargo.lock`. The successful retry supplied lock SHA-256 `4f17b356fa7bac392b6d801d1c74fb9e36b6517f9465c8ebc19bb9a2792a84c5`; restoring that file to the PR is required.",
    "",
    f"The failed build, successful build, four validation commands, and two replay processes charge {charge['wall_seconds_sum']:.6f} wall seconds, {charge['total_core_seconds_sum']:.6f} core-seconds, and {charge['peak_rss_bytes_max']} bytes maximum RSS across {charge['components']} components.",
    "",
    "The dense pair default is validated on this one target. The result remains implementation engineering, not a full attack or SOTA evidence.",
    "",
]
(STAGE / "RESULTS.md").write_text("\n".join(lines))
print(json.dumps({"decision": result["decision"], "same_binary_native_over_direct": result["same_binary_native_over_direct"], "campaign_charge": charge}, indent=2))

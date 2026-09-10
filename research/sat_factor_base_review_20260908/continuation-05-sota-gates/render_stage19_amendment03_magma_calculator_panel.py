#!/usr/bin/env python3
"""Render the nine-task Stage 19 Amendment 03 Magma Calculator successor."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import stat
import tempfile
from typing import Any


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PROTOCOL = HERE / "stage-19-amendment-03-magma-calculator-panel-protocol.json"
DEFAULT_ARTIFACT = HERE / "stage-19-magma-calculator-panel-amendment-03-20260910"
STAGE13_PROTOCOL = HERE / "stage-13-pdp-panel-protocol.json"
STAGE13_CUSTODY = HERE / "stage-13-custody.json"
STAGE15_VERIFIER = HERE / "verify_stage15_magma_calculator.py"
STAGE16_VERIFIER = HERE / "verify_stage16_magma_compact.py"
STAGE17_VERIFIER = HERE / "verify_stage17_magma_n41.py"
STAGE15_SUMMARY = HERE / "stage-15-magma-calculator-20260909" / "stage-15-summary.json"
STAGE16_SUMMARY = HERE / "stage-16-magma-compact-20260909" / "stage-16-summary.json"
STAGE17_SUMMARY = HERE / "stage-17-magma-n41-20260909" / "stage-17-summary.json"
AMENDMENT03_VERIFIER = HERE / "verify_stage19_amendment03.py"
AMENDMENT03_SUMMARY = HERE / "stage-19-amendment-03-summary.json"
PREDECESSOR_ARTIFACTS = (
    HERE / "stage-19-magma-calculator-panel-20260909",
    HERE / "stage-19-magma-calculator-panel-amendment-01-20260910",
    HERE / "stage-19-magma-calculator-panel-amendment-02-20260910",
)
PLAN_SCHEMA = "koblitz_magma_calculator_stage19_amendment03_plan.v1"


class RenderError(RuntimeError):
    """The frozen source inventory cannot produce the declared panel."""


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RenderError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


STAGE17 = load_module("stage17_for_stage19_renderer", STAGE17_VERIFIER)
STAGE16 = STAGE17.STAGE16
STAGE15 = STAGE16.STAGE15
AMENDMENT03 = load_module("amendment03_for_stage19_renderer", AMENDMENT03_VERIFIER)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def read_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise RenderError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise RenderError(f"expected JSON object in {path}")
    return value


def task_key(seed: int, cell_id: str) -> str:
    return f"seed-{seed}/{cell_id}"


def input_name(ordinal: int, seed: int, cell_id: str) -> str:
    return f"{ordinal:02d}-seed-{seed}-{cell_id}.magma"


def source_paths(seed: int, cell_id: str) -> tuple[Path, Path, Path]:
    anf, magma = STAGE16.task_paths(seed, cell_id)
    binding = (
        HERE
        / "stage-13-panel-20260909"
        / "tasks"
        / f"seed-{seed}"
        / cell_id
        / "source-binding.json"
    )
    return anf, magma, binding


def relative_record(path: Path) -> dict:
    data = path.read_bytes()
    return {
        "path": str(path.relative_to(REPO)),
        "bytes": len(data),
        "sha256": sha256_bytes(data),
    }


def validate_protocol(protocol: dict) -> None:
    if protocol.get("schema") != "koblitz_magma_calculator_stage19_amendment03_protocol.v1":
        raise RenderError("unexpected Stage 19 Amendment 03 protocol schema")
    if protocol.get("stage") != 19 or protocol.get("amendment") != 3:
        raise RenderError("Stage 19 Amendment 03 protocol number changed")
    service = protocol.get("service")
    policy = protocol.get("execution_policy")
    admission = protocol.get("admission")
    if not all(isinstance(value, dict) for value in (service, policy, admission)):
        raise RenderError("Stage 19 Amendment 03 protocol lacks service, policy, or admission")
    if service != {
        "endpoint": "https://magma.maths.usyd.edu.au/xml/calculator.xml",
        "form_field": "input",
        "published_calculator_page": "https://magma.maths.usyd.edu.au/calc/",
        "expected_version": "2.29-10",
        "max_input_bytes": 50000,
        "max_service_seconds": 60,
        "client_timeout_seconds": 75,
        "max_response_bytes": 1048576,
        "user_agent": "aburan28-crypto-koblitz-stage19-amendment03/1",
    }:
        raise RenderError("Stage 19 Amendment 03 calculator service contract changed")
    if policy != {
        "task_order": "remaining seed-major frozen order",
        "parallel_workers": 1,
        "attempts_per_task": 1,
        "retries_per_task": 0,
        "inter_request_delay_seconds": 2,
        "halt_after_nonclean_receipt": True,
        "resume_skips_every_started_task": True,
        "halted_state_is_terminal": True,
        "predecessor_ledgers_are_terminal_and_never_resumed": True,
        "predecessor_started_tasks_excluded": True,
        "explicit_execute_flag_required": True,
    }:
        raise RenderError("Stage 19 Amendment 03 one-attempt sequential policy changed")
    tasks = protocol.get("tasks")
    if not isinstance(tasks, list) or len(tasks) != 9:
        raise RenderError("Stage 19 Amendment 03 protocol must contain exactly nine tasks")
    if [row.get("ordinal") for row in tasks if isinstance(row, dict)] != list(range(1, 10)):
        raise RenderError("Stage 19 Amendment 03 task ordinals changed")
    if [row.get("original_stage19_ordinal") for row in tasks] != list(range(2, 11)):
        raise RenderError("Stage 19 Amendment 03 original ordinals changed")
    if len({(row.get("seed"), row.get("cell_id")) for row in tasks}) != 9:
        raise RenderError("Stage 19 Amendment 03 task inventory contains a duplicate")
    prior = protocol.get("retained_prior_f4_terminals")
    if not isinstance(prior, list) or len(prior) != 5:
        raise RenderError("Stage 19 Amendment 03 prior retained-terminal inventory changed")
    if protocol.get("stage13_archive_commit") != "5dbe54dd468e914b5872d10d06e7c44d01b0a236":
        raise RenderError("Stage 13 archive pin changed")
    if protocol.get("prepared_from_revision") != AMENDMENT03.PREDECESSOR_COMMIT:
        raise RenderError("Stage 19 Amendment 03 predecessor commit changed")
    admission = protocol["admission"]
    if (
        admission.get("generic_backslash_unfolding_permitted") is not False
        or "high-half marker" not in str(admission.get("source_identity"))
    ):
        raise RenderError("Stage 19 Amendment 03 split-identity policy changed")


def predecessor_attempt_starts(protocol: dict) -> list[dict]:
    declared = protocol.get("predecessors")
    expected_names = [path.name for path in PREDECESSOR_ARTIFACTS]
    if not isinstance(declared, list) or [row.get("artifact") for row in declared] != expected_names:
        raise RenderError("Stage 19 Amendment 03 predecessor inventory changed")
    records = []
    for root, row in zip(PREDECESSOR_ARTIFACTS, declared, strict=True):
        paths = sorted(root.rglob("attempt-start.json"))
        if len(paths) != 1:
            raise RenderError(f"{root.name}: expected exactly one predecessor attempt-start")
        start = read_json(paths[0])
        allowed = row.get("allowed_attempted_task_ids")
        if allowed != [start.get("task_id")]:
            raise RenderError(f"{root.name}: predecessor attempted-task declaration changed")
        records.append(
            {
                "artifact": root.name,
                "task_id": start["task_id"],
                "attempt_ordinal": start["attempt_ordinal"],
                "attempt_start": relative_record(paths[0]),
            }
        )
    return records


def verified_recovered_terminal(protocol: dict) -> tuple[dict, dict]:
    summary = AMENDMENT03.summarize()
    if AMENDMENT03_SUMMARY.read_bytes() != canonical_bytes(summary):
        raise RenderError("Amendment 03 adjudication summary no longer regenerates byte-for-byte")
    recovered = summary["recovered_terminal"]
    declared = protocol.get("recovered_terminal")
    if (
        declared != {
            "task_id": recovered["task_id"],
            "amendment_summary": AMENDMENT03_SUMMARY.name,
            "classification": recovered["classification"],
        }
        or recovered.get("source_equivalent") is not True
    ):
        raise RenderError("protocol recovered-terminal declaration changed")
    receipt = {
        "task_id": recovered["task_id"],
        "seed": 2026091301,
        "cell_id": "n31-l5-m3-ggmp-a0-f0",
        "stage": 19,
        "amendment": 3,
        "source": "additive_exact_fold_adjudication",
        "response_sha256": recovered["response"]["sha256"],
    }
    return receipt, relative_record(AMENDMENT03_SUMMARY)


def verified_prior_terminals(protocol: dict) -> tuple[list[dict], list[dict]]:
    stage15 = STAGE15.summarize()
    stage16 = STAGE16.summarize()
    stage17 = STAGE17.summarize()
    summaries = [
        (15, STAGE15_SUMMARY, stage15),
        (16, STAGE16_SUMMARY, stage16),
        (17, STAGE17_SUMMARY, stage17),
    ]
    summary_records = []
    for stage, path, summary in summaries:
        rendered = canonical_bytes(summary)
        if path.is_symlink() or not path.is_file() or path.read_bytes() != rendered:
            raise RenderError(f"Stage {stage} expected summary no longer regenerates byte-for-byte")
        summary_records.append({"stage": stage, **relative_record(path)})
    observed = [
        *(
            {
                "seed": int(case["seed"]),
                "cell_id": str(case["cell_id"]),
                "stage": 15,
                "response_sha256": case["adapted_submission"]["response"]["response_sha256"],
            }
            for case in stage15["cases"]
            if case["adapted_submission"]["classification"]
            == "sat_basis_certificate_unverified_model"
        ),
        *(
            {
                "seed": int(case["seed"]),
                "cell_id": str(case["cell_id"]),
                "stage": 16,
                "response_sha256": case["response_sha256"],
            }
            for case in stage16["cases"]
            if case["classification"] == "sat_basis_certificate_unverified_model"
        ),
        {
            "seed": int(stage17["source_equivalence"]["seed"]),
            "cell_id": str(stage17["source_equivalence"]["cell_id"]),
            "stage": 17,
            "response_sha256": stage17["response_sha256"],
        },
    ]
    declared = [
        {"seed": row["seed"], "cell_id": row["cell_id"], "stage": row["stage"]}
        for row in protocol["retained_prior_f4_terminals"]
    ]
    if [
        {"seed": row["seed"], "cell_id": row["cell_id"], "stage": row["stage"]}
        for row in observed
    ] != declared:
        raise RenderError("protocol prior coverage differs from reverified Stages 15-17")
    return observed, summary_records


def render_bound_named(
    task_id: str, source_sha256: str, variables: int, equations: list[list[int]]
) -> bytes:
    if not source_sha256 or len(source_sha256) != 64 or set(source_sha256) - set("0123456789abcdef"):
        raise RenderError(f"{task_id}: source SHA-256 is malformed")
    if any(character not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_/" for character in task_id):
        raise RenderError(f"unsafe task identity marker {task_id!r}")
    base = STAGE17.render_named(variables, equations)
    marker = b'printf "KOBLITZ_MAGMA_SCHEMA=koblitz_magma_f4_terminal.v1\\n";'
    if base.count(marker) != 1:
        raise RenderError("named renderer terminal marker changed")
    identity = (
        f'printf "KOBLITZ_MAGMA_TASK_ID={task_id}\\n";'
        f'printf "KOBLITZ_MAGMA_SOURCE_SHA256_HI={source_sha256[:32]}\\n";'
        f'printf "KOBLITZ_MAGMA_SOURCE_SHA256_LO={source_sha256[32:]}\\n";'
    ).encode()
    for line in (
        f"KOBLITZ_MAGMA_TASK_ID={task_id}",
        f"KOBLITZ_MAGMA_SOURCE_SHA256_HI={source_sha256[:32]}",
        f"KOBLITZ_MAGMA_SOURCE_SHA256_LO={source_sha256[32:]}",
    ):
        if len(line.encode()) > 79 or "\\" in line:
            raise RenderError(f"{task_id}: identity output could be folded")
    return base.replace(marker, identity + marker)


def build_plan() -> tuple[dict, dict[str, bytes]]:
    protocol = read_json(PROTOCOL)
    validate_protocol(protocol)
    equivalence = STAGE16.complete_equivalence_inventory()
    if (
        equivalence.get("tasks") != 20
        or equivalence.get("all_equation_monomial_lists_equal") is not True
    ):
        raise RenderError("all-twenty Stage 13 source equivalence is not established")
    named = STAGE17.named_inventory(equivalence)
    eligible = [task_key(row["seed"], row["cell_id"]) for row in named["eligible_tasks"]]
    prior_receipts, prior_summary_records = verified_prior_terminals(protocol)
    recovered_receipt, recovered_summary_record = verified_recovered_terminal(protocol)
    attempt_starts = predecessor_attempt_starts(protocol)
    prior = {task_key(row["seed"], row["cell_id"]) for row in prior_receipts}
    expected = [
        task_key(int(row["seed"]), str(row["cell_id"])) for row in protocol["tasks"]
    ]
    remaining = [key for key in eligible if key not in prior]
    recovered_id = recovered_receipt["task_id"]
    predecessor_ids = {row["task_id"] for row in attempt_starts}
    if (
        len(eligible) != 15
        or remaining != [recovered_id, *expected]
        or prior - set(eligible)
        or predecessor_ids != {recovered_id}
        or predecessor_ids & set(expected)
    ):
        raise RenderError("successor tasks are not exactly the nine never-started eligible tasks")

    by_key = {
        task_key(int(row["seed"]), str(row["cell_id"])): row
        for row in equivalence["rows"]
    }
    named_by_key = {
        task_key(int(row["seed"]), str(row["cell_id"])): row
        for row in named["rows"]
    }
    inputs: dict[str, bytes] = {}
    tasks = []
    for declared in protocol["tasks"]:
        ordinal = int(declared["ordinal"])
        original_ordinal = int(declared["original_stage19_ordinal"])
        seed = int(declared["seed"])
        cell_id = str(declared["cell_id"])
        key = task_key(seed, cell_id)
        source = by_key[key]
        named_source = named_by_key[key]
        anf_path, magma_path, binding_path = source_paths(seed, cell_id)
        variables, equations = STAGE16.parse_anf(anf_path)
        base_rendered = STAGE17.render_named(variables, equations)
        filename = input_name(ordinal, seed, cell_id)
        if (
            len(base_rendered) != named_source["bytes"]
            or sha256_bytes(base_rendered) != named_source["sha256"]
        ):
            raise RenderError(f"{key}: named renderer or eligibility changed")
        binding = read_json(binding_path)
        source_sha256 = binding.get("source_instance_sha256")
        rendered = render_bound_named(key, source_sha256, variables, equations)
        if len(rendered) > int(protocol["service"]["max_input_bytes"]):
            raise RenderError(f"{key}: identity-bound named input exceeds calculator cap")
        if source["anf_sha256"] != sha256_file(anf_path):
            raise RenderError(f"{key}: ANF hash changed after equivalence check")
        if source["verbose_magma_sha256"] != sha256_file(magma_path):
            raise RenderError(f"{key}: verbose Magma hash changed after equivalence check")
        inputs[filename] = rendered
        tasks.append(
            {
                "ordinal": ordinal,
                "original_stage19_ordinal": original_ordinal,
                "id": key,
                "seed": seed,
                "cell_id": cell_id,
                "variables": variables,
                "equations": len(equations),
                "source_instance_sha256": source_sha256,
                "producer_source_instance_blake3": binding.get(
                    "producer_source_instance_blake3"
                ),
                "source_anf": relative_record(anf_path),
                "source_verbose_magma": relative_record(magma_path),
                "source_binding": relative_record(binding_path),
                "named_input": {
                    "path": f"inputs/{filename}",
                    "bytes": len(rendered),
                    "sha256": sha256_bytes(rendered),
                },
                "equation_monomial_lists_equal": True,
                "response_identity_markers": {
                    "task_id": key,
                    "source_instance_sha256_hi": source_sha256[:32],
                    "source_instance_sha256_lo": source_sha256[32:],
                    "generic_backslash_unfolding_permitted": False,
                },
                "controls": {
                    "boolean_ring_order": "grevlex",
                    "gpu_disabled": True,
                    "deterministic_seed": 1,
                    "algorithm": "direct-f4-sparse",
                    "dense": False,
                    "threads": 1,
                },
            }
        )
    dependencies = [
        STAGE13_PROTOCOL,
        STAGE13_CUSTODY,
        STAGE15_VERIFIER,
        STAGE16_VERIFIER,
        STAGE17_VERIFIER,
        STAGE15_SUMMARY,
        STAGE16_SUMMARY,
        STAGE17_SUMMARY,
        AMENDMENT03_VERIFIER,
        AMENDMENT03_SUMMARY,
    ]
    plan = {
        "schema": PLAN_SCHEMA,
        "protocol": relative_record(PROTOCOL),
        "source_archive": {
            "commit": protocol["stage13_archive_commit"],
            "source_revision": protocol["stage13_source_revision"],
            "verified_tasks": 20,
            "all_equation_monomial_lists_equal": True,
        },
        "dependencies": [relative_record(path) for path in dependencies],
        "representation": protocol["admission"]["representation"],
        "named_eligible_tasks": 15,
        "retained_prior_f4_terminals": 6,
        "retained_prior_receipts": [*prior_receipts, recovered_receipt],
        "retained_prior_expected_summaries": [
            *prior_summary_records,
            {"stage": 19, "amendment": 3, **recovered_summary_record},
        ],
        "predecessor_attempt_starts": attempt_starts,
        "predecessor_attempted_task_ids": sorted(predecessor_ids),
        "successor_tasks_disjoint_from_predecessor_attempts": True,
        "planned_requests": len(tasks),
        "task_order": protocol["execution_policy"]["task_order"],
        "tasks": tasks,
        "execution_policy": protocol["execution_policy"],
        "service": protocol["service"],
        "claim_boundary": protocol["claim_boundary"],
    }
    return plan, inputs


def ensure_single_link_regular(path: Path, label: str) -> None:
    info = path.lstat()
    if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
        raise RenderError(f"{label} is not a single-link regular file: {path}")


def atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.parent.is_symlink() or not path.parent.is_dir():
        raise RenderError(f"atomic-write parent is not a real directory: {path.parent}")
    if path.exists() or path.is_symlink():
        ensure_single_link_regular(path, "atomic-write destination")
    descriptor, temporary_name = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.tmp-")
    temporary = Path(temporary_name)
    try:
        info = os.fstat(descriptor)
        if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
            raise RenderError("atomic temporary is not a single-link regular file")
        with os.fdopen(descriptor, "wb", closefd=True) as handle:
            descriptor = -1
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        ensure_single_link_regular(temporary, "atomic temporary")
        os.replace(temporary, path)
        directory_fd = os.open(path.parent, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def write_new_or_identical(path: Path, data: bytes, refresh: bool = False) -> None:
    if path.exists() or path.is_symlink():
        if path.is_symlink() or not path.is_file():
            raise RenderError(f"refusing to replace changed artifact {path}")
        ensure_single_link_regular(path, "prepared artifact")
        if path.read_bytes() == data:
            return
        if not refresh:
            raise RenderError(f"refusing to replace changed artifact {path}")
    atomic_write(path, data)


def materialize(artifact: Path, refresh: bool = False) -> dict:
    if refresh and any(
        (artifact / name).exists() or (artifact / name).is_symlink()
        for name in ("execution-manifest.json", "run.json", "attempts", "summary.json")
    ):
        raise RenderError("prepared artifacts cannot refresh after execution binding or attempts")
    plan, inputs = build_plan()
    write_new_or_identical(artifact / "plan.json", canonical_bytes(plan), refresh)
    for filename, data in inputs.items():
        write_new_or_identical(artifact / "inputs" / filename, data, refresh)
    discovered = sorted(path.name for path in (artifact / "inputs").glob("*.magma"))
    if discovered != sorted(inputs):
        raise RenderError("prepared artifact contains an unexpected named input")
    return plan


def self_test() -> dict:
    plan, inputs = build_plan()
    if plan["planned_requests"] != 9 or len(inputs) != 9:
        raise AssertionError("wrong Stage 19 Amendment 03 request inventory")
    if plan["tasks"][0]["id"] != "seed-2026091302/n31-l5-m3-standard-a1-f0":
        raise AssertionError("wrong first Stage 19 Amendment 03 task")
    if plan["tasks"][-1]["id"] != "seed-2026091305/n41-l5-m3-standard-a1-f0":
        raise AssertionError("wrong last Stage 19 Amendment 03 task")
    if any(task["named_input"]["bytes"] > 50_000 for task in plan["tasks"]):
        raise AssertionError("Stage 19 Amendment 03 includes an over-cap input")
    return {
        "self_test": "pass",
        "planned_requests": len(inputs),
        "predecessor_attempted_tasks": 1,
        "split_source_identity_markers": True,
        "plan_sha256": sha256_bytes(canonical_bytes(plan)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_ARTIFACT)
    parser.add_argument("--materialize", action="store_true")
    parser.add_argument("--refresh-prepared", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        print(json.dumps(self_test(), indent=2, sort_keys=True))
        return
    if args.materialize and args.refresh_prepared:
        parser.error("choose --materialize or --refresh-prepared")
    plan = (
        materialize(args.artifact.resolve(), refresh=args.refresh_prepared)
        if args.materialize or args.refresh_prepared
        else build_plan()[0]
    )
    print(canonical_bytes(plan).decode(), end="")


if __name__ == "__main__":
    main()

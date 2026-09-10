#!/usr/bin/env python3
"""Verify the bounded Stage 15 Magma Calculator receipts fail closed."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import sys
import xml.etree.ElementTree as ET


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PANEL = HERE / "stage-13-panel-20260909"
ARTIFACT = HERE / "stage-15-magma-calculator-20260909"
RUNNER = REPO / "scripts" / "run_koblitz_pdp_matrix.py"
CELL = "n31-l5-m3-standard-a1-f0"
SEEDS = (2026091301, 2026091303)
REMOVED_CALCULATOR_LINE = b"SetNthreads(1);\n"
ENDPOINT = "https://magma.maths.usyd.edu.au/xml/calculator.xml"
OBSERVED_VERSION = "2.29-10"
OBSERVED_MAX_TIME = 60
OBSERVED_MAX_INPUT = 50_000
WITNESS_SUFFIX = b'''cpu_start := Cputime();
wall_start := Realtime();
sat, S := SAT(F);
cpu_seconds := Cputime(cpu_start);
wall_seconds := Realtime(wall_start);
printf "KOBLITZ_MAGMA_WITNESS_SCHEMA=koblitz_magma_sat_witness.v1\\n";
if sat then
  printf "KOBLITZ_MAGMA_WITNESS_STATUS=SAT\\n";
  printf "KOBLITZ_MAGMA_WITNESS_ASSIGNMENT=%o\\n", S;
else
  printf "KOBLITZ_MAGMA_WITNESS_STATUS=UNSAT\\n";
  printf "KOBLITZ_MAGMA_WITNESS_ASSIGNMENT=[]\\n";
end if;
printf "KOBLITZ_MAGMA_WITNESS_CPU_SECONDS=%o\\n", cpu_seconds;
printf "KOBLITZ_MAGMA_WITNESS_WALL_SECONDS=%o\\n", wall_seconds;
quit;
'''

sys.path.insert(0, str(REPO))
from scripts.run_koblitz_pdp_matrix import parse_magma_terminal  # noqa: E402


class VerificationError(RuntimeError):
    """A retained service receipt violates its source or claim boundary."""


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def read_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise VerificationError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise VerificationError(f"expected JSON object in {path}")
    return value


def exact_instance(seed: int) -> Path:
    return (
        PANEL
        / "tasks"
        / f"seed-{seed}"
        / CELL
        / "matrix"
        / CELL
        / "instance.magma"
    )


def source_binding(seed: int) -> Path:
    return PANEL / "tasks" / f"seed-{seed}" / CELL / "source-binding.json"


def parse_memory_mb(value: str) -> float:
    match = re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)MB", value)
    if match is None:
        raise VerificationError(f"unexpected calculator memory field {value!r}")
    return float(match.group(1))


def parse_calculator_xml(path: Path) -> dict:
    try:
        root = ET.fromstring(path.read_bytes())
    except (OSError, ET.ParseError) as error:
        raise VerificationError(f"cannot parse calculator response {path}: {error}") from error
    if root.tag != "calculator":
        raise VerificationError(f"{path}: unexpected XML root")
    root_children = list(root)
    if [child.tag for child in root_children] != ["headers", "results"]:
        raise VerificationError(f"{path}: missing calculator headers or results")
    headers_node, results_node = root_children
    header_nodes = list(headers_node)
    header_names = [child.tag for child in header_nodes]
    if len(set(header_names)) != len(header_names):
        raise VerificationError(f"{path}: duplicate calculator header")
    headers = {child.tag: child.text or "" for child in header_nodes}
    required = {"max_time", "max_input", "seed", "version", "time", "memory"}
    if not required.issubset(headers):
        raise VerificationError(f"{path}: incomplete calculator headers")
    try:
        maximum_time = int(headers["max_time"])
        maximum_input = int(headers["max_input"])
        server_seed = int(headers["seed"])
        service_time = float(headers["time"])
    except ValueError as error:
        raise VerificationError(f"{path}: malformed numeric calculator header") from error
    if (
        maximum_time != OBSERVED_MAX_TIME
        or maximum_input != OBSERVED_MAX_INPUT
        or server_seed < 0
        or not math.isfinite(service_time)
        or service_time < 0
    ):
        raise VerificationError(f"{path}: unexpected calculator resource header")
    if headers["version"] != OBSERVED_VERSION:
        raise VerificationError(f"{path}: unexpected Magma version {headers['version']!r}")
    result_nodes = list(results_node)
    if any(line.tag != "line" or len(line) != 0 for line in result_nodes):
        raise VerificationError(f"{path}: calculator results are not flat line records")
    lines = [line.text or "" for line in result_nodes]
    output = "\n".join(lines) + "\n"
    return {
        "response_path": path.name,
        "response_bytes": path.stat().st_size,
        "response_sha256": sha256_file(path),
        "service": {
            "endpoint": ENDPOINT,
            "version": headers["version"],
            "max_time_seconds": maximum_time,
            "max_input_bytes": maximum_input,
            "server_seed": server_seed,
            "service_reported_time_seconds": service_time,
            "service_reported_memory_text": headers["memory"],
            "service_reported_memory_mb": parse_memory_mb(headers["memory"]),
            "warning": headers.get("warning"),
        },
        "output": output,
    }


def verify_source(seed: int) -> dict:
    input_path = exact_instance(seed)
    if not input_path.is_file() or input_path.is_symlink():
        raise VerificationError(f"seed {seed}: exact Magma input is missing or a symlink")
    binding = read_json(source_binding(seed))
    export = binding.get("exports", {}).get("magma_boolean_f4")
    if not isinstance(export, dict):
        raise VerificationError(f"seed {seed}: source binding lacks Magma export custody")
    exact = input_path.read_bytes()
    if export.get("path") != "instance.magma":
        raise VerificationError(f"seed {seed}: unexpected Magma export path")
    if export.get("bytes") != len(exact) or export.get("sha256") != sha256_bytes(exact):
        raise VerificationError(f"seed {seed}: exact input differs from Stage 13 custody")
    if not exact.startswith(REMOVED_CALCULATOR_LINE):
        raise VerificationError(f"seed {seed}: exact input lacks frozen global-thread request")
    required = (
        b"SetGPU(false);",
        b"BooleanPolynomialRing(",
        b'GroebnerBasis(I : Al := "Direct", Faugere := true, Dense := false, Nthreads := 1)',
    )
    if any(marker not in exact for marker in required):
        raise VerificationError(f"seed {seed}: exact input lacks a required F4 control")
    return {
        "path": str(input_path.relative_to(REPO)),
        "bytes": len(exact),
        "sha256": sha256_bytes(exact),
        "source_instance_sha256": binding.get("source_instance_sha256"),
        "producer_source_instance_blake3": binding.get("producer_source_instance_blake3"),
    }


def panel_inventory() -> dict:
    protocol = read_json(HERE / "stage-13-pdp-panel-protocol.json")
    seeds = protocol.get("replicate_seeds")
    cells = protocol.get("cells")
    if not isinstance(seeds, list) or not isinstance(cells, list):
        raise VerificationError("Stage 13 protocol lacks its frozen task inventory")
    paths = []
    for seed in seeds:
        for cell in cells:
            cell_id = cell.get("id")
            path = (
                PANEL
                / "tasks"
                / f"seed-{seed}"
                / str(cell_id)
                / "matrix"
                / str(cell_id)
                / "instance.magma"
            )
            if not path.is_file() or path.is_symlink():
                raise VerificationError(f"missing frozen Magma input {path}")
            paths.append(path)
    discovered = sorted(PANEL.glob("tasks/seed-*/*/matrix/*/instance.magma"))
    if sorted(paths) != discovered or len(paths) != 20:
        raise VerificationError("Stage 13 Magma inventory is not the exact frozen 20 tasks")
    eligible = [path for path in paths if path.stat().st_size < OBSERVED_MAX_INPUT]
    expected_eligible = [exact_instance(seed) for seed in SEEDS]
    if sorted(eligible) != sorted(expected_eligible):
        raise VerificationError("calculator-eligible frozen input set changed")
    ineligible_sizes = [path.stat().st_size for path in paths if path not in eligible]
    return {
        "tasks": len(paths),
        "eligible": len(eligible),
        "without_retained_response": len(paths) - len(eligible),
        "ineligible_min_bytes": min(ineligible_sizes),
        "ineligible_max_bytes": max(ineligible_sizes),
    }


def terminal_identity(terminal: dict) -> dict:
    return {
        key: terminal[key]
        for key in (
            "schema",
            "algorithm",
            "terminal_status",
            "f4_step_degrees",
            "basis_size",
            "single_thread_requested",
            "gpu_disabled",
        )
    }


def require_clean_f4_output(output: str, seed: int) -> None:
    numeric = r"[0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?"
    pattern = (
        r"KOBLITZ_MAGMA_SCHEMA=koblitz_magma_f4_terminal\.v1\n"
        r"KOBLITZ_MAGMA_ALGORITHM=direct-f4-sparse\n"
        r"KOBLITZ_MAGMA_STATUS=(?:SAT|UNSAT)\n"
        r"KOBLITZ_MAGMA_F4_DEGREES=\[\s*(?:[0-9]+(?:\s*,\s*[0-9]+)*)?\s*\]\n"
        r"KOBLITZ_MAGMA_BASIS_SIZE=[0-9]+\n"
        rf"KOBLITZ_MAGMA_CPU_SECONDS={numeric}\n"
        rf"KOBLITZ_MAGMA_WALL_SECONDS={numeric}\n\n"
    )
    if re.fullmatch(pattern, output) is None:
        raise VerificationError(f"seed {seed}: adapted response contains nonterminal output")


def verify_case(seed: int) -> dict:
    source = verify_source(seed)
    exact_bytes = exact_instance(seed).read_bytes()
    adapted_path = ARTIFACT / f"seed-{seed}-adapted-input.magma"
    adapted = adapted_path.read_bytes()
    if adapted != exact_bytes[len(REMOVED_CALCULATOR_LINE) :]:
        raise VerificationError(
            f"seed {seed}: calculator adaptation changed more than SetNthreads(1)"
        )
    if len(adapted) >= 50_000:
        raise VerificationError(f"seed {seed}: adapted input exceeds the calculator cap")

    exact_response = parse_calculator_xml(ARTIFACT / f"seed-{seed}-exact-response.xml")
    exact_terminal = parse_magma_terminal(exact_response["output"])
    if exact_terminal is None:
        raise VerificationError(f"seed {seed}: exact response lacks a complete F4 marker block")
    exact_warning = exact_response["service"]["warning"]
    if (
        exact_warning is None
        or "SetNthreads" not in exact_response["output"]
        or "Illegal operation" not in exact_response["output"]
    ):
        raise VerificationError(f"seed {seed}: exact response lacks its retained setup error")

    adapted_response = parse_calculator_xml(ARTIFACT / f"seed-{seed}-adapted-response.xml")
    require_clean_f4_output(adapted_response["output"], seed)
    adapted_terminal = parse_magma_terminal(adapted_response["output"])
    if adapted_terminal is None or adapted_response["service"]["warning"] is not None:
        raise VerificationError(f"seed {seed}: adapted response is not a clean F4 terminal")
    if terminal_identity(exact_terminal) != terminal_identity(adapted_terminal):
        raise VerificationError(f"seed {seed}: exact and adapted F4 terminals disagree")
    if adapted_terminal["terminal_status"] != "sat":
        raise VerificationError(f"seed {seed}: planted instance did not produce a proper ideal")

    witness_input = (ARTIFACT / f"seed-{seed}-sat-witness-input.magma").read_bytes()
    source_prefix = adapted.split(b"I := ideal<R | F>;\n", 1)[0]
    if witness_input != source_prefix + WITNESS_SUFFIX:
        raise VerificationError(f"seed {seed}: witness attempt changed its source or operation")
    witness_response = parse_calculator_xml(
        ARTIFACT / f"seed-{seed}-sat-witness-response.xml"
    )
    if (
        witness_response["service"]["warning"] is None
        or "GetTempDir" not in witness_response["output"]
        or "KOBLITZ_MAGMA_WITNESS_ASSIGNMENT=" in witness_response["output"]
    ):
        raise VerificationError(f"seed {seed}: witness failure boundary changed")

    return {
        "seed": seed,
        "cell_id": CELL,
        "source": source,
        "exact_submission": {
            "classification": "f4_markers_after_calculator_setup_error",
            "scientific_terminal_admitted": False,
            "terminal": exact_terminal,
            "response": exact_response,
        },
        "calculator_adaptation": {
            "rule": "remove exactly the unsupported global SetNthreads(1) line; retain Nthreads := 1 on GroebnerBasis",
            "path": adapted_path.name,
            "bytes": len(adapted),
            "sha256": sha256_bytes(adapted),
            "source_system_byte_identical": True,
            "f4_controls_retained": True,
        },
        "adapted_submission": {
            "classification": "sat_basis_certificate_unverified_model",
            "terminal": adapted_terminal,
            "response": adapted_response,
            "process_scoped_resources_complete": False,
            "point_witness_validated": False,
        },
        "witness_attempt": {
            "classification": "calculator_temp_directory_unavailable_operational",
            "path": f"seed-{seed}-sat-witness-input.magma",
            "bytes": len(witness_input),
            "sha256": sha256_bytes(witness_input),
            "source_system_byte_identical": True,
            "response": witness_response,
            "assignment_obtained": False,
        },
    }


def verify_timeout_attempts() -> dict:
    attempts = read_json(ARTIFACT / "timeout-attempts.json")
    if attempts.get("schema") != "koblitz_magma_calculator_timeout_attempts.v1":
        raise VerificationError("unexpected timeout-attempt schema")
    rows = attempts.get("attempts")
    if not isinstance(rows, list) or len(rows) != 2:
        raise VerificationError("expected the retained variety and basis-export attempts")
    expected = {
        "variety": ARTIFACT / "seed-2026091301-variety-input.magma",
        "basis_export": ARTIFACT / "seed-2026091301-basis-export-input.magma",
    }
    if {row.get("kind") for row in rows if isinstance(row, dict)} != set(expected):
        raise VerificationError("timeout attempt kinds are missing or duplicated")
    if attempts.get("recorded_at") != "2026-09-09":
        raise VerificationError("unexpected timeout-attempt record date")
    adapted = (ARTIFACT / "seed-2026091301-adapted-input.magma").read_bytes()
    if not adapted.endswith(b"quit;\n"):
        raise VerificationError("adapted input lacks its terminal quit")
    common_prefix = adapted[: -len(b"quit;\n")]
    expected_suffix = {
        "variety": b'''variety_cpu_start := Cputime();
variety_wall_start := Realtime();
V := VarietySequence(I);
variety_cpu_seconds := Cputime(variety_cpu_start);
variety_wall_seconds := Realtime(variety_wall_start);
printf "KOBLITZ_MAGMA_VARIETY_SCHEMA=koblitz_magma_variety_witness.v1\\n";
printf "KOBLITZ_MAGMA_VARIETY_SIZE=%o\\n", #V;
if #V gt 0 then
  printf "KOBLITZ_MAGMA_VARIETY_ASSIGNMENT=%o\\n", V[1];
else
  printf "KOBLITZ_MAGMA_VARIETY_ASSIGNMENT=[]\\n";
end if;
printf "KOBLITZ_MAGMA_VARIETY_CPU_SECONDS=%o\\n", variety_cpu_seconds;
printf "KOBLITZ_MAGMA_VARIETY_WALL_SECONDS=%o\\n", variety_wall_seconds;
quit;
''',
        "basis_export": b'''printf "KOBLITZ_MAGMA_BASIS_BEGIN\\n";
print G: Magma;
printf "KOBLITZ_MAGMA_BASIS_END\\n";
quit;
''',
    }
    for row in rows:
        path = expected.get(row.get("kind"))
        if path is None:
            raise VerificationError("unknown timeout attempt kind")
        data = path.read_bytes()
        if (
            row.get("input_path") != path.name
            or row.get("input_bytes") != len(data)
            or row.get("input_sha256") != sha256_bytes(data)
            or row.get("client_watchdog_seconds") != 75
            or row.get("curl_exit_code") != 28
            or row.get("response_file_created") is not False
            or row.get("classification") != "transport_timeout_operational"
            or row.get("endpoint") != ENDPOINT
            or row.get("asserts_nothing_about")
            != "Magma variety or basis-export correctness, runtime, memory, or terminal status"
        ):
            raise VerificationError(f"invalid {row.get('kind')} timeout receipt")
        if data != common_prefix + expected_suffix[row["kind"]]:
            raise VerificationError(f"{row['kind']} timeout input changed its source or operation")
    return attempts


def summarize() -> dict:
    inventory = panel_inventory()
    cases = [verify_case(seed) for seed in SEEDS]
    timeouts = verify_timeout_attempts()
    expected_responses = {
        ARTIFACT / f"seed-{seed}-{kind}-response.xml"
        for seed in SEEDS
        for kind in ("exact", "adapted", "sat-witness")
    }
    if set(ARTIFACT.glob("*-response.xml")) != expected_responses:
        raise VerificationError("retained calculator response inventory changed")
    return {
        "schema": "koblitz_magma_calculator_stage15_summary.v1",
        "evidence_class": "exploratory_external_service_receipt",
        "service": {
            "endpoint": ENDPOINT,
            "observed_version": OBSERVED_VERSION,
            "observed_max_time_seconds": OBSERVED_MAX_TIME,
            "observed_max_input_bytes": OBSERVED_MAX_INPUT,
        },
        "frozen_panel_tasks": inventory["tasks"],
        "exact_archived_inputs_under_service_cap": inventory["eligible"],
        "frozen_inputs_without_retained_response": inventory["without_retained_response"],
        "ineligible_input_size_range_bytes": [
            inventory["ineligible_min_bytes"],
            inventory["ineligible_max_bytes"],
        ],
        "exact_archived_submissions": 2,
        "clean_source_equivalent_f4_terminals": 2,
        "validated_magma_point_witnesses": 0,
        "magma_process_resource_receipts_complete": False,
        "full_magma_matrix_executed": False,
        "full_solver_matrix_gate_passed": False,
        "cases": cases,
        "timeout_attempts": timeouts,
        "claim": (
            "Two target-matched n=31 standard systems have supporting Magma F4 proper-ideal "
            "certificates. They lack Magma-derived point witnesses and process-scoped resource "
            "receipts; the other eighteen frozen inputs have no retained calculator responses "
            "and were outside the observed input cap. This is partial external-service evidence, "
            "not a complete Magma benchmark or SOTA result."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--expected", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        assert parse_memory_mb("32.09MB") == 32.09
        assert terminal_identity(
            {
                "schema": "s",
                "algorithm": "a",
                "terminal_status": "sat",
                "f4_step_degrees": [2, 3],
                "basis_size": 2,
                "single_thread_requested": True,
                "gpu_disabled": True,
                "cpu_seconds": 1.0,
                "wall_seconds": 2.0,
            }
        )["basis_size"] == 2
        print(json.dumps({"self_test": "pass"}, indent=2))
        return
    summary = summarize()
    rendered = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    if args.expected is not None and args.expected.read_text() != rendered:
        raise VerificationError("recomputed Stage 15 summary differs from expected summary")
    if args.output is not None:
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()

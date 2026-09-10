#!/usr/bin/env python3
"""Run matched binary-Koblitz PDP instances through available solver backends.

Every subprocess is wrapped by macOS ``/usr/bin/time -lp``.  Missing tools and
watchdog expirations are recorded as operational outcomes; they are never
converted to UNSAT or scientific evidence.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time


def run_timed(command: list[str], timeout: float, cwd: Path) -> dict:
    meter = Path(__file__).with_name("process_meter.py")
    with tempfile.TemporaryDirectory() as directory:
        scratch = Path(directory)
        stdout_path = scratch / "stdout"
        stderr_path = scratch / "stderr"
        metrics_path = scratch / "metrics.json"
        subprocess.run(
            [
                sys.executable,
                str(meter),
                "--cwd",
                str(cwd),
                "--timeout",
                str(timeout),
                "--stdout",
                str(stdout_path),
                "--stderr",
                str(stderr_path),
                "--metrics",
                str(metrics_path),
                "--",
                *command,
            ],
            check=True,
        )
        record = json.loads(metrics_path.read_text())
        record["stdout"] = stdout_path.read_text()
        record["stderr"] = stderr_path.read_text()
        return record


def parse_wdsat_model(stdout: str, n_vars: int) -> list[bool] | None:
    for line in stdout.splitlines():
        value = line.strip()
        if len(value) >= n_vars and set(value) <= {"0", "1"}:
            return [bit == "1" for bit in value[:n_vars]]
    return None


def validate_wdsat_anf(path: Path, model: list[bool]) -> bool:
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    for line in lines[1:]:
        tokens = line.split()[1:-1]
        monomials: list[list[int]] = []
        has_t = False
        i = 0
        while i < len(tokens):
            token = tokens[i]
            if token == "T":
                has_t = True
                i += 1
            elif token.startswith("."):
                degree = int(token[1:])
                monomials.append([int(value) - 1 for value in tokens[i + 1 : i + 1 + degree]])
                i += degree + 1
            else:
                monomials.append([int(token) - 1])
                i += 1
        # Export convention: T is present exactly when the original
        # polynomial has no constant term.
        parity = not has_t
        for monomial in monomials:
            parity ^= all(model[index] for index in monomial)
        if parity:
            return False
    return True


def parse_cms_model(stdout: str, max_var: int) -> list[bool] | None:
    values: dict[int, bool] = {}
    for line in stdout.splitlines():
        if not line.startswith("v "):
            continue
        for token in line.split()[1:]:
            literal = int(token)
            if literal:
                values[abs(literal)] = literal > 0
    if not all(index in values for index in range(1, max_var + 1)):
        return None
    return [values[index] for index in range(1, max_var + 1)]


def validate_xor_dimacs(path: Path, model: list[bool]) -> bool:
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line[0] in {"c", "p"}:
            continue
        xor = line.startswith("x")
        tokens = line.split()[1:-1] if xor else line.split()[:-1]
        literals = [int(value) for value in tokens]
        values = [model[abs(literal) - 1] if literal > 0 else not model[abs(literal) - 1] for literal in literals]
        if xor:
            if sum(values) % 2 != 1:
                return False
        elif not any(values):
            return False
    return True


def solver_status(run: dict, solver: str, instance: Path, manifest: dict) -> dict:
    if run["timed_out"]:
        status = "timeout_inconclusive"
        conflicts = None
        model_valid = None
    elif solver == "wdsat":
        lines = [line.strip() for line in run["stdout"].splitlines() if line.strip()]
        model = parse_wdsat_model(run["stdout"], manifest["source_variables"])
        if model is not None:
            status = "sat"
            model_valid = validate_wdsat_anf(instance / "instance.anf", model)
        elif any("UNSAT" in line for line in lines):
            status = "unsat"
            model_valid = None
        else:
            status = "solver_error"
            model_valid = None
        conflicts = next((int(line) for line in reversed(lines) if line.isdigit()), None)
    elif solver == "cryptominisat":
        if "s SATISFIABLE" in run["stdout"]:
            status = "sat"
            max_var = manifest["exports"]["cryptominisat_xor_dimacs"]["variables"]
            model = parse_cms_model(run["stdout"], max_var)
            model_valid = None if model is None else validate_xor_dimacs(instance / "instance.xor.cnf", model)
        elif "s UNSATISFIABLE" in run["stdout"]:
            status = "unsat"
            model_valid = None
        elif run["returncode"] == 0:
            status = "unknown_inconclusive"
            model_valid = None
        else:
            status = "solver_error"
            model_valid = None
        matches = re.findall(r"(?:conflicts|Conflicts)\s*[:=]?\s*(\d+)", run["stdout"] + run["stderr"])
        conflicts = int(matches[-1]) if matches else None
    else:
        text = run["stdout"] + run["stderr"]
        status = "completed" if run["returncode"] == 0 else "solver_error"
        conflicts = None
        model_valid = None
        if "Runtime error" in text or "Error" in text:
            status = "solver_error"
    return {
        "solver": solver,
        "status": status,
        "conflicts": conflicts,
        "source_model_valid": model_valid,
        "returncode": run["returncode"],
        "timed_out": run["timed_out"],
        "metrics": run["metrics"],
        "command": run["command"],
    }


def version(binary: str | None) -> dict:
    if binary is None:
        return {"available": False, "path": None, "version": None}
    path = shutil.which(binary) if "/" not in binary else binary
    if not path or not Path(path).exists():
        return {"available": False, "path": path, "version": None}
    for option in ([path, "--version"], [path, "-h"]):
        try:
            result = subprocess.run(option, capture_output=True, text=True, timeout=10)
            text = (result.stdout + result.stderr).strip()
            if text:
                return {"available": True, "path": str(Path(path).resolve()), "version": text[:1000]}
        except (OSError, subprocess.SubprocessError):
            pass
    return {"available": True, "path": str(Path(path).resolve()), "version": "unreported"}


def build_wdsat(source: Path, instance: Path, manifest: dict, timeout: float) -> tuple[dict, str | None]:
    """Build an instance-sized WDSat binary and return its charged receipt."""
    build_root = instance / "wdsat-build"
    build_root.mkdir()
    copy_receipt = run_timed(
        ["cp", "-R", str((source / "src").resolve()), str(build_root / "src")],
        timeout,
        instance,
    )
    cms = manifest["exports"]["cryptominisat_xor_dimacs"]
    n_source = int(manifest["source_variables"])
    max_degree = int(manifest["source_max_degree"])
    max_id = int(cms["variables"])
    max_eq = int(cms["cnf_clauses"])
    max_xeq = int(cms["xor_rows"])
    max_terms = int(manifest["source_max_monomials_per_equation"])
    config = "\n".join(
        [
            "#define __XG_ENHANCED__",
            f"#define __MAX_ANF_ID__ {n_source + 1}",
            f"#define __MAX_DEGREE__ {max_degree + 1}",
            f"#define __MAX_ID__ {max_id}",
            f"#define __MAX_BUFFER_SIZE__ {max(5000, max_terms * 8, max_id * 16)}",
            f"#define __MAX_EQ__ {max(64, max_eq + 16)}",
            f"#define __MAX_EQ_SIZE__ {max_degree + 2}",
            f"#define __MAX_XEQ__ {max(8, max_xeq + 2)}",
            f"#define __MAX_XEQ_SIZE__ {max(max_id + 1, max_terms + 2)}",
            "",
        ]
    )
    config_start = time.perf_counter()
    (build_root / "src" / "config.h").write_text(config)
    config_wall = time.perf_counter() - config_start
    run = run_timed(["make", "-C", "src"], timeout, build_root)
    (build_root / "build.stdout").write_text(run["stdout"])
    (build_root / "build.stderr").write_text(run["stderr"])
    binary = build_root / "wdsat_solver"
    receipt = {
        "status": "completed" if run["returncode"] == 0 and binary.exists() else (
            "timeout_inconclusive" if run["timed_out"] else "failed_operational"
        ),
        "source_copy": {
            "returncode": copy_receipt["returncode"],
            "timed_out": copy_receipt["timed_out"],
            "metrics": copy_receipt["metrics"],
        },
        "configuration_wall_seconds": config_wall,
        "metrics": run["metrics"],
        "command": run["command"],
        "config": config,
    }
    return receipt, str(binary.resolve()) if binary.exists() else None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="new result directory")
    parser.add_argument("--exporter", type=Path, required=True, help="built koblitz_pdp_export binary")
    parser.add_argument("--wdsat", help="WDSat binary")
    parser.add_argument("--wdsat-source", type=Path, help="WDSat checkout to right-size and build per cell")
    parser.add_argument("--cryptominisat", help="CryptoMiniSat binary")
    parser.add_argument("--magma", help="Magma binary")
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--conflicts", type=int, default=100_000)
    parser.add_argument(
        "--configs",
        default="15:5:standard,31:5:standard,31:5:ggmp,41:5:standard,59:9:standard,67:9:standard",
        help="comma-separated n:ell:basis cells",
    )
    parser.add_argument("--seed", type=int, default=20260909)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output must be a new path")
    args.output.mkdir(parents=True)
    if not args.exporter.exists():
        parser.error("exporter does not exist; build the release example first")

    tools = {
        "wdsat": (
            {
                "available": args.wdsat_source is not None and (args.wdsat_source / "src" / "makefile").exists(),
                "path": str(args.wdsat_source.resolve()) if args.wdsat_source else None,
                "version": "per-instance source build",
            }
            if args.wdsat_source
            else version(args.wdsat)
        ),
        "cryptominisat": version(args.cryptominisat),
        "magma": version(args.magma),
    }
    report = {
        "schema": "koblitz_pdp_matched_matrix.v1",
        "scope": "PDP construction and solver matrix; not an end-to-end ECDLP result",
        "host": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": platform.python_version(),
            "logical_cpus": os.cpu_count(),
        },
        "policy": {
            "unknown_is_unsat": False,
            "all_process_resources_charged": True,
            "single_thread_requested": True,
            "missing_tool_is_negative_evidence": False,
        },
        "tools": tools,
        "instances": [],
    }

    for cell_index, cell in enumerate(args.configs.split(",")):
        n_text, ell_text, basis = cell.split(":")
        n, ell = int(n_text), int(ell_text)
        instance = args.output / f"n{n}-l{ell}-m3-{basis}"
        command = [
            str(args.exporter.resolve()),
            str(n),
            str(ell),
            basis,
            str(args.seed + cell_index),
            str(args.conflicts),
            str(instance),
        ]
        generated = run_timed(command, args.timeout, args.exporter.parent)
        (args.output / f"n{n}-l{ell}-m3-{basis}.export.stdout").write_text(generated["stdout"])
        (args.output / f"n{n}-l{ell}-m3-{basis}.export.stderr").write_text(generated["stderr"])
        entry = {
            "cell": {"n": n, "ell": ell, "m": 3, "basis": basis},
            "generator": {
                "status": "timeout_inconclusive" if generated["timed_out"] else (
                    "completed" if generated["returncode"] == 0 else "failed_operational"
                ),
                "returncode": generated["returncode"],
                "metrics": generated["metrics"],
                "command": command,
            },
            "solvers": [],
        }
        manifest_path = instance / "manifest.json"
        if not manifest_path.exists():
            report["instances"].append(entry)
            continue
        manifest = json.loads(manifest_path.read_text())
        entry["manifest"] = manifest

        wdsat_binary = tools["wdsat"]["path"]
        if args.wdsat_source and tools["wdsat"]["available"]:
            entry["wdsat_build"], wdsat_binary = build_wdsat(
                args.wdsat_source, instance, manifest, args.timeout
            )
        if tools["wdsat"]["available"] and wdsat_binary:
            branch_variables = ",".join(str(index) for index in range(1, 3 * ell + 1))
            run = run_timed(
                [wdsat_binary, "-i", str(instance / "instance.anf"), "-g", branch_variables],
                args.timeout,
                instance,
            )
            (instance / "wdsat.stdout").write_text(run["stdout"])
            (instance / "wdsat.stderr").write_text(run["stderr"])
            entry["solvers"].append(solver_status(run, "wdsat", instance, manifest))
        else:
            entry["solvers"].append({"solver": "wdsat", "status": "unavailable_operational"})

        if tools["cryptominisat"]["available"]:
            run = run_timed(
                [tools["cryptominisat"]["path"], "--verb", "1", "--threads", "1", str(instance / "instance.xor.cnf")],
                args.timeout,
                instance,
            )
            (instance / "cryptominisat.stdout").write_text(run["stdout"])
            (instance / "cryptominisat.stderr").write_text(run["stderr"])
            entry["solvers"].append(solver_status(run, "cryptominisat", instance, manifest))
        else:
            entry["solvers"].append({"solver": "cryptominisat", "status": "unavailable_operational"})

        if tools["magma"]["available"]:
            run = run_timed(
                [tools["magma"]["path"], "-b", str(instance / "instance.magma")],
                args.timeout,
                instance,
            )
            (instance / "magma.stdout").write_text(run["stdout"])
            (instance / "magma.stderr").write_text(run["stderr"])
            entry["solvers"].append(solver_status(run, "magma-f4", instance, manifest))
        else:
            entry["solvers"].append({"solver": "magma-f4", "status": "unavailable_operational"})
        report["instances"].append(entry)
        (args.output / "progress.json").write_text(json.dumps(report, indent=2) + "\n")

    (args.output / "result.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

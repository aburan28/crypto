#!/usr/bin/env sage -python
"""Run the frozen matched radical pair-image solver-stage suite."""

from __future__ import annotations

import argparse
import hashlib
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

from sage.env import SAGE_VERSION

from radical_pair import BinaryCell, PrimeCell, system_hash


ROUND_RE = re.compile(r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+x\s+(\d+)\s+")


def sha256_bytes(blob):
    return hashlib.sha256(blob).hexdigest()


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def command_output(command):
    return subprocess.run(command, capture_output=True, text=True, check=True).stdout.strip()


def msolve_input(polynomials):
    ring = polynomials[0].parent()
    lines = [",".join(ring.variable_names()), str(ring.base_ring().characteristic())]
    lines.append(",\n".join(p._repr_().replace(" ", "") for p in polynomials))
    return "\n".join(lines) + "\n"


def parse_log(stdout):
    rounds = []
    for line in stdout.splitlines():
        match = ROUND_RE.match(line)
        if match:
            degree, selected, pairs, rows, columns = map(int, match.groups())
            rounds.append({"degree": degree, "selected": selected, "pairs": pairs,
                           "rows": rows, "columns": columns,
                           "matrix_area": rows * columns})
    max_matrix = re.search(r"max\. matrix data\s+(\d+)\s+x\s+(\d+)", stdout)
    reported = re.search(r"overall\(elapsed\)\s+([0-9.]+) sec", stdout)
    return {
        "rounds": rounds,
        "max_f4_degree": max((r["degree"] for r in rounds), default=None),
        "max_matrix_rows": int(max_matrix.group(1)) if max_matrix else None,
        "max_matrix_columns": int(max_matrix.group(2)) if max_matrix else None,
        "peak_round_matrix_area": max((r["matrix_area"] for r in rounds), default=0),
        "reported_elapsed_seconds": float(reported.group(1)) if reported else None,
    }


def run_process(msolve, input_path, stdout_path, stderr_path, timeout):
    command = [msolve, "-f", str(input_path), "-g", "1", "-v", "2",
               "-t", "1", "--random-seed", "0"]
    started = time.perf_counter()
    try:
        completed = subprocess.run(command, capture_output=True, text=True, timeout=timeout)
        status = "COMPLETED" if completed.returncode == 0 else "ERROR"
        stdout, stderr, returncode = completed.stdout, completed.stderr, completed.returncode
    except subprocess.TimeoutExpired as error:
        status = "TIMEOUT"
        stdout = error.stdout or ""
        stderr = error.stderr or ""
        if isinstance(stdout, bytes):
            stdout = stdout.decode(errors="replace")
        if isinstance(stderr, bytes):
            stderr = stderr.decode(errors="replace")
        returncode = None
    wall = time.perf_counter() - started
    stdout_path.write_text(stdout)
    stderr_path.write_text(stderr)
    parsed = parse_log(stdout)
    return {
        "status": status,
        "returncode": returncode,
        "wall_seconds": wall,
        "command": command,
        "stdout": str(stdout_path),
        "stderr": str(stderr_path),
        "stdout_sha256": sha256_bytes(stdout.encode()),
        "stderr_sha256": sha256_bytes(stderr.encode()),
        **parsed,
    }


def source_hashes(root):
    files = [root / "radical_pair.py", root / "run.py", root / "summarize.py",
             root / "audit.py", root / "test_radical_pair.py", root / "contract.json"]
    return {path.name: sha256_file(path) for path in files if path.exists()}


def target_values(cell, indices=None):
    if indices is None:
        return list(cell.roots)
    values = list(cell.field)
    return [values[index] for index in indices]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=120.0)
    args = parser.parse_args()

    output = Path(args.output).resolve()
    if output.exists():
        raise SystemExit(f"refusing to overwrite {output}")
    raw_dir, input_dir = output / "raw", output / "inputs"
    raw_dir.mkdir(parents=True)
    input_dir.mkdir()
    root = Path(__file__).resolve().parent
    msolve = shutil.which("msolve")
    if not msolve:
        raise SystemExit("msolve is required")

    suite = [
        ("prime-main-b3", lambda: PrimeCell(17, 2, 2, 3), [0, 1, 3]),
        ("prime-main-b5", lambda: PrimeCell(17, 2, 2, 5), [0, 1, 3]),
        ("prime-main-b7", lambda: PrimeCell(17, 2, 2, 7), [0, 1, 3]),
        ("prime-main-b9", lambda: PrimeCell(17, 2, 2, 9), [0, 1, 3]),
        ("prime-holdout-b5", lambda: PrimeCell(29, 1, 1, 5), [0, 1, 6]),
        ("prime-holdout-b7", lambda: PrimeCell(29, 1, 1, 7), [0, 1, 6]),
        ("prime-holdout-b9", lambda: PrimeCell(29, 1, 1, 9), [0, 1, 6]),
        ("prime-holdout-b11", lambda: PrimeCell(29, 1, 1, 11), [0, 1, 6]),
        ("binary-k2", lambda: BinaryCell(2), None),
        ("binary-k3", lambda: BinaryCell(3), None),
    ]

    processes = []
    certificates = {}
    cells = []
    for cell_index, (label, factory, indices) in enumerate(suite):
        setup_started = time.perf_counter()
        cell = factory()
        setup_seconds = time.perf_counter() - setup_started
        targets = target_values(cell, indices)
        variants = ["direct", "radical_symmetric"]
        if isinstance(cell, PrimeCell):
            variants.insert(1, "norm_symmetric")
        cell_record = {
            "label": label,
            "cell_id": cell.cell_id,
            "field": str(cell.field),
            "factor_base": [str(x) for x in cell.roots],
            "factor_base_size": len(cell.roots),
            "targets": [str(x) for x in targets],
            "setup_seconds": setup_seconds,
            "radical_preprocessing": cell.preprocessing,
            "variants": variants,
        }
        cells.append(cell_record)

        for target_index, target in enumerate(targets):
            verification_started = time.perf_counter()
            certificate = cell.certificate(target)
            verification_seconds = time.perf_counter() - verification_started
            certificate["verification_seconds"] = verification_seconds
            certificate["target"] = str(target)
            certificates[f"{label}/target-{target_index}"] = certificate

            encoding_started = time.perf_counter()
            systems = cell.systems(target)
            encoding_seconds = time.perf_counter() - encoding_started
            for variant, polynomials in systems.items():
                content = msolve_input(polynomials)
                input_path = input_dir / f"{label}-t{target_index}-{variant}.msolve"
                input_path.write_text(content)
                input_sha = sha256_bytes(content.encode())
                variant_encoding = encoding_seconds / len(systems)
                # Alternate execution order by rotating below; inputs are frozen first.
                for repetition in range(args.repetitions):
                    processes.append({
                        "cell_index": cell_index,
                        "cell": label,
                        "target_index": target_index,
                        "target": str(target),
                        "variant": variant,
                        "repetition": repetition,
                        "input": str(input_path.relative_to(output)),
                        "input_sha256": input_sha,
                        "system_sha256": system_hash(polynomials),
                        "variables": polynomials[0].parent().ngens(),
                        "equations": len(polynomials),
                        "input_degree": max(int(f.total_degree()) for f in polynomials),
                        "encoding_seconds_share": variant_encoding,
                        "certificate_key": f"{label}/target-{target_index}",
                    })

    # Rotate variant order inside each repetition/target to avoid a fixed warm-order bias.
    processes.sort(key=lambda item: (
        item["cell_index"], item["target_index"], item["repetition"],
        (item["cell_index"] + item["target_index"] + item["repetition"]
         + ["direct", "norm_symmetric", "radical_symmetric"].index(item["variant"])) % 3,
    ))

    with (output / "processes.jsonl").open("w") as process_file:
        for ordinal, record in enumerate(processes):
            stem = (f"{ordinal:04d}-{record['cell']}-t{record['target_index']}-"
                    f"{record['variant']}-r{record['repetition']}")
            result = run_process(
                msolve,
                output / record["input"],
                raw_dir / f"{stem}.stdout",
                raw_dir / f"{stem}.stderr",
                args.timeout,
            )
            result["stdout"] = str(Path(result["stdout"]).relative_to(output))
            result["stderr"] = str(Path(result["stderr"]).relative_to(output))
            record.update(result)
            record["ordinal"] = ordinal
            process_file.write(json.dumps(record, sort_keys=True) + "\n")
            process_file.flush()
            print(record["cell"], record["target_index"], record["variant"],
                  record["repetition"], record["status"], record["max_f4_degree"],
                  f"{record['wall_seconds']:.4f}", flush=True)

    (output / "cells.json").write_text(json.dumps(cells, indent=2, sort_keys=True) + "\n")
    (output / "certificates.json").write_text(
        json.dumps(certificates, indent=2, sort_keys=True) + "\n")
    hardware = {
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "python": sys.version,
        "cpu_count": os.cpu_count(),
    }
    try:
        hardware["cpu_brand"] = command_output(["sysctl", "-n", "machdep.cpu.brand_string"])
    except Exception:
        hardware["cpu_brand"] = None
    manifest = {
        "status": "complete",
        "base_git_commit": command_output(["git", "rev-parse", "HEAD"]),
        "source_hashes": source_hashes(root),
        "sage": SAGE_VERSION,
        "msolve_version": command_output([msolve, "--version"]).splitlines()[0],
        "msolve_path": msolve,
        "repetitions": args.repetitions,
        "timeout_seconds": args.timeout,
        "msolve_options": ["-g", "1", "-v", "2", "-t", "1", "--random-seed", "0"],
        "processes": len(processes),
        "hardware": hardware,
        "evidence_type": "measured_local_solver_stage",
        "end_to_end_metrics": None,
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()

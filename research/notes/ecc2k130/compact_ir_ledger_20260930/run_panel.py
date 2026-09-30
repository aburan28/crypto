#!/usr/bin/env python3
"""Whole-process Callgrind Ir panel on frozen compact/rho source and public Q."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
CONFIG = HERE / "CONFIG.json"


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_cell(cell_id: str) -> tuple[dict, dict, dict, Path, Path]:
    config = json.loads(CONFIG.read_text())
    assert config["schema"] == "ecc2k130-compact-ir-ledger-v1"
    cells = [cell for cell in config["cells"] if cell["id"] == cell_id]
    assert len(cells) == 1, cell_id
    cell = cells[0]
    assert (cell["n"], cell["L"]) in ((37, 1), (37, 1024), (41, 1),
                                       (41, 1024), (53, 1), (53, 1024))
    assert cell["arms"] in (["off", "rho"], ["off", "blocked", "rho"])
    assert cell["arms"] == (["off", "blocked", "rho"] if cell["L"] == 1024
                             and cell["n"] in (41, 53) else ["off", "rho"])
    assert cell.get("repeat_control", []) == (["off", "rho"] if cell_id == "n37_L1" else [])
    source_path = ROOT / config["source_freeze"]
    point_path = ROOT / config["point_panel_freeze"]
    assert sha(source_path) == config["source_freeze_sha256"]
    assert sha(point_path) == config["point_panel_freeze_sha256"]
    source = json.loads(source_path.read_text())
    assert source["schema"] == "compact-s3-prefilter-evaluation-freeze-v1"
    frozen = (json.loads(point_path.read_text()) if cell["freeze_kind"] == "point_panel"
              else source)
    spec = frozen["specs"][cell["freeze_key"]]
    assert (spec["n"], spec["L"], spec["a"]) == (cell["n"], cell["L"], 0)
    assert spec["subgroup_order"] > 0 and spec["seed"] > 0
    directory = ROOT / "research/notes/ecc2k130" / (
        "compact_orbit_point_panel_20260929" if cell["freeze_kind"] == "point_panel"
        else "compact_s3_prefilter_20260930")
    points = directory / spec["points_file"]
    fixture = directory / spec["fixture_file"]
    assert sha(points) == cell["points_sha256"] == spec["points_sha256"]
    assert sha(fixture) == cell["fixture_sha256"] == spec["fixture_sha256"]
    assert sum(bool(line.strip()) for line in points.read_text().splitlines()) == cell["L"]
    if cell["L"] == 1024:
        selected = (source["k_by_n"] if cell["freeze_kind"] == "prefilter"
                    else frozen["batch_k"])
        assert cell["k"] == selected[str(cell["n"])]
    else:
        assert cell["k"] in frozen["single_k_grid"][str(cell["n"])]
    assert (config["rank_seed"], config["s3_batch_window"],
            config["rho_dp_bits"], config["rho_canonicalization_backend"],
            config["rho_quotient_mode"]) == (7, 64, 4, "normal_basis", "signed_frobenius")
    return config, cell, spec, points, fixture


def parse_ir(path: Path) -> int:
    """Use Callgrind's whole-run summary, not a possibly smaller line total."""
    events: list[str] | None = None
    summary: int | None = None
    totals: int | None = None
    with path.open(errors="replace") as stream:
        for line in stream:
            if line.startswith("events:"):
                events = line.split()[1:]
                assert events.count("Ir") == 1
            if line.startswith(("summary:", "totals:")):
                assert events is not None
                values = [int(value) for value in line.split()[1:]]
                assert len(values) == len(events)
                value = values[events.index("Ir")]
                if line.startswith("summary:"):
                    summary = value
                else:
                    totals = value
    result = summary if summary is not None else totals
    assert result is not None and result > 0, path
    if summary is not None and totals is not None:
        assert 0 < totals <= summary, (summary, totals)
    return result


def descendant_pids(pid: int) -> set[int]:
    remaining = [pid]
    found: set[int] = set()
    while remaining:
        child = remaining.pop()
        if child in found:
            continue
        found.add(child)
        path = Path(f"/proc/{child}/task/{child}/children")
        try:
            remaining.extend(int(value) for value in path.read_text().split())
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            pass
    return found


def rss_bytes(pid: int) -> int:
    total = 0
    for child in descendant_pids(pid):
        try:
            for line in Path(f"/proc/{child}/status").read_text().splitlines():
                if line.startswith("VmRSS:"):
                    total += int(line.split()[1]) * 1024
                    break
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            pass
    return total


def run_child(command: list[str], env: dict[str, str], stdout: Path,
              stderr: Path, timeout: int, rss_limit: int) -> dict:
    started = time.monotonic()
    peak_rss = 0
    stopped_for: str | None = None
    with stdout.open("wb") as out, stderr.open("wb") as err:
        process = subprocess.Popen(command, stdout=out, stderr=err, env=env,
                                   start_new_session=True)
        while process.poll() is None:
            peak_rss = max(peak_rss, rss_bytes(process.pid))
            elapsed = time.monotonic() - started
            if elapsed > timeout:
                stopped_for = "timeout"
            elif peak_rss > rss_limit:
                stopped_for = "rss_limit"
            if stopped_for:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                break
            time.sleep(0.2)
        exit_code = process.wait()
    return {"exit_code": exit_code, "stopped_for": stopped_for,
            "elapsed_under_backend_seconds_not_a_cost": time.monotonic() - started,
            "observed_peak_rss_bytes": peak_rss}


def write_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def run(cell_id: str, frozen_root: Path, batch: Path, rho: Path, materialization: Path,
        output: Path, backend: str, smoke_arm: str | None) -> dict:
    assert sys.platform == "linux" or backend == "native", "Callgrind requires Linux"
    assert backend in ("callgrind", "native")
    assert smoke_arm is None or backend == "native"
    config, cell, spec, points, _fixture = load_cell(cell_id)
    source = json.loads((ROOT / config["source_freeze"]).read_text())
    receipt = json.loads(materialization.read_text())
    assert receipt["schema"] == "compact-frozen-source-materialization-v1"
    assert receipt["pinned_files"] == len(source["source_sha256"])
    for filename, expected in source["source_sha256"].items():
        assert sha(frozen_root / filename) == expected, filename
    assert batch.is_file() and rho.is_file()
    assert not output.exists(), "never overwrite a run directory"
    output.mkdir(parents=True)
    if backend == "callgrind":
        version = subprocess.check_output(["valgrind", "--version"], text=True).strip()
        assert version.startswith("valgrind-")
    else:
        version = None
    sequence = list(cell["arms"])
    sequence.extend(f"{arm}_repeat" for arm in cell.get("repeat_control", []))
    if smoke_arm is not None:
        assert smoke_arm in sequence
        sequence = [smoke_arm]
    report = {"schema": "ecc2k130-compact-ir-run-v1", "status": "RUNNING",
              "backend": backend, "cell": cell, "spec": spec,
              "config_sha256": sha(CONFIG), "source_freeze_sha256": sha(ROOT / config["source_freeze"]),
              "materialization_sha256": sha(materialization),
              "materialization": receipt,
              "host": {"platform": platform.platform(), "machine": platform.machine(),
                       "python": sys.version, "valgrind": version,
                       "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"],
                                                           cwd=ROOT, text=True).strip()},
              "binaries": {"compact_sha256": sha(batch), "rho_sha256": sha(rho)},
              "points_sha256": sha(points), "sequence": sequence, "runs": []}
    write_json(output / "run.json", report)
    for name in sequence:
        policy = name.removesuffix("_repeat")
        prefix = output / name
        environment = {key: value for key, value in os.environ.items()
                       if not key.startswith("KIC_")}
        environment.update({"RAYON_NUM_THREADS": "1", "LC_ALL": "C"})
        if policy in ("off", "blocked"):
            environment.update({
                "KIC_DUMP_BASE": str(prefix.with_suffix(".base.jsonl")),
                "KIC_DUMP_RANK": str(prefix.with_suffix(".rank.jsonl")),
                "KIC_S3_BATCH_WINDOW": str(config["s3_batch_window"]),
                "KIC_S3_PREFILTER": policy,
            })
            program = [str(batch), f"construct:{cell['n']}:0:{cell['k']}",
                       str(points), str(config["rank_seed"]),
                       str(prefix.with_suffix(".target.jsonl"))]
        else:
            assert policy == "rho"
            environment.update({
                "KIC_RHO_POINT_INPUT": str(points),
                "KIC_RHO_BATCH_CORPUS": spec["corpus"],
                "KIC_RHO_DP_BITS": str(config["rho_dp_bits"]),
                "KIC_RHO_CANON_BACKEND": config["rho_canonicalization_backend"],
            })
            program = [str(rho), str(cell["n"]), "0", config["rho_quotient_mode"],
                       str(cell["L"]), str(spec["seed"])]
        callgrind = prefix.with_suffix(".callgrind.out")
        command = ([
            "valgrind", "--tool=callgrind", "--instr-atstart=yes",
            "--collect-atstart=yes",
            "--cache-sim=no", "--branch-sim=no", "--error-exitcode=97",
            f"--callgrind-out-file={callgrind}", *program,
        ] if backend == "callgrind" else program)
        assert not any(".fixture.jsonl" in value for value in
                       command + list(environment.values()))
        stdout, stderr = prefix.with_suffix(".stdout.jsonl"), prefix.with_suffix(".stderr.txt")
        result = run_child(command, environment, stdout, stderr,
                           config["timeout_seconds_per_arm"],
                           config["rss_limit_bytes_per_arm"])
        ir = None
        if backend == "callgrind" and result["exit_code"] == 0 and result["stopped_for"] is None:
            try:
                ir = parse_ir(callgrind)
            except (AssertionError, ValueError, FileNotFoundError) as error:
                result["parse_error"] = f"{type(error).__name__}: {error}"
        record = {"arm": name, "policy": policy, "command": command,
                  "environment": {key: environment[key] for key in environment
                                  if key.startswith("KIC_") or key in ("RAYON_NUM_THREADS", "LC_ALL")},
                  "stdout": stdout.name, "stdout_sha256": sha(stdout),
                  "stderr": stderr.name, "stderr_sha256": sha(stderr),
                  "callgrind": callgrind.name if callgrind.exists() else None,
                  "callgrind_sha256": sha(callgrind) if callgrind.exists() else None,
                  "callgrind_bytes": callgrind.stat().st_size if callgrind.exists() else None,
                  "Ir": ir, **result}
        report["runs"].append(record)
        write_json(output / "run.json", report)
        print(json.dumps({"cell": cell_id, "arm": name, "Ir": ir,
                          "exit_code": result["exit_code"],
                          "stopped_for": result["stopped_for"]}), flush=True)
    required = len(report["runs"]) == len(sequence) and all(
        item["exit_code"] == 0 and item["stopped_for"] is None
        and (backend == "native" or item["Ir"] is not None) for item in report["runs"])
    if backend == "callgrind" and cell_id == "n37_L1" and required:
        by_name = {item["arm"]: item["Ir"] for item in report["runs"]}
        report["deterministic_control"] = {
            arm: abs(by_name[arm] - by_name[f"{arm}_repeat"]) / by_name[arm]
            for arm in ("off", "rho")}
        required = all(value <= 0.001 for value in report["deterministic_control"].values())
    report["status"] = ("SMOKE_PASS" if required and backend == "native"
                        else "PASS" if required else "FAIL")
    write_json(output / "run.json", report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cell", required=True)
    parser.add_argument("--frozen-root", type=Path, required=True)
    parser.add_argument("--batch", type=Path, required=True)
    parser.add_argument("--rho", type=Path, required=True)
    parser.add_argument("--materialization", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--backend", choices=("callgrind", "native"), default="callgrind")
    parser.add_argument("--smoke-arm")
    args = parser.parse_args()
    result = run(args.cell, args.frozen_root.resolve(), args.batch.resolve(),
                 args.rho.resolve(), args.materialization.resolve(),
                 args.out.resolve(), args.backend, args.smoke_arm)
    if result["status"] not in ("PASS", "SMOKE_PASS"):
        raise SystemExit(1)


if __name__ == "__main__":
    main()

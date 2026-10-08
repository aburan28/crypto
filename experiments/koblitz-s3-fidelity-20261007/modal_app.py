"""Modal CPU diagnostic for the frozen N41/N53 paired-S3 experiment.

Run from this checkout:
  modal run experiments/koblitz-s3-fidelity-20261007/modal_app.py::probe
  modal run experiments/koblitz-s3-fidelity-20261007/modal_app.py::smoke
  modal run experiments/koblitz-s3-fidelity-20261007/modal_app.py::pilot

Modal container evidence is exploratory unless a separate auditable host-level
isolation receipt satisfies the repository's CPU gate. No result from this
module assigns an IC1 candidate ID or declares a controlled speedup.
"""

from __future__ import annotations

import base64
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import statistics
import subprocess
import tempfile
import time

import modal


REMOTE_ROOT = Path("/root/crypto")
LOCAL_ROOT = Path(__file__).resolve().parents[2] if modal.is_local() else REMOTE_ROOT
PILOT = REMOTE_ROOT / "experiments/koblitz-s3-fidelity-20261007/pilot"
WARMUP = REMOTE_ROOT / "experiments/koblitz-s3-fidelity-20261007/warmup"
BIN = REMOTE_ROOT / "target/release/examples"
SOURCE_COMMIT = "833ef86191947ba8e4def56a78ee6c51a4ea0ff2"
WORKERS = 14
COLUMNS = 244
RELATION_SEED = 20260928
CPU_REQUEST = 16.0  # Modal describes this as physical cores.
MEMORY_MIB = 8192


def local_file(path: str, remote: str):
    return (Path(LOCAL_ROOT) / path, str(REMOTE_ROOT / remote))


image = (
    modal.Image.from_registry("rust:1.98-bookworm", add_python="3.12")
    .apt_install("linux-perf", "numactl", "util-linux", "procps")
    .add_local_file(*local_file("Cargo.toml", "Cargo.toml"), copy=True)
    .add_local_file(*local_file("Cargo.lock", "Cargo.lock"), copy=True)
    .add_local_dir(LOCAL_ROOT / "src", remote_path=str(REMOTE_ROOT / "src"), copy=True)
    .add_local_dir(LOCAL_ROOT / "docs", remote_path=str(REMOTE_ROOT / "docs"), copy=True)
    .add_local_dir(LOCAL_ROOT / "tests", remote_path=str(REMOTE_ROOT / "tests"), copy=True)
    .add_local_dir(LOCAL_ROOT / "benches", remote_path=str(REMOTE_ROOT / "benches"), copy=True)
    .add_local_file(*local_file("experiments/koblitz-s3-pair-query-20261007-v4/baseline.rs", "experiments/koblitz-s3-pair-query-20261007-v4/baseline.rs"), copy=True)
    .add_local_file(*local_file("experiments/koblitz-s3-pair-query-20261007-v4/candidate.rs", "experiments/koblitz-s3-pair-query-20261007-v4/candidate.rs"), copy=True)
    .add_local_dir(
        LOCAL_ROOT / "experiments/koblitz-s3-fidelity-20261007/pilot",
        remote_path=str(PILOT), copy=True,
    )
    .add_local_dir(
        LOCAL_ROOT / "experiments/koblitz-s3-fidelity-20261007/warmup",
        remote_path=str(WARMUP), copy=True,
    )
    .run_commands(
        f"cd {REMOTE_ROOT} && cargo build --release --locked "
        "--example s3_pair_baseline_frozen --example s3_pair_candidate_frozen"
    )
)
app = modal.App("s3-pair-two-curve-modal-pilot")


def _read(path: str, limit: int = 400_000) -> str | None:
    try:
        return Path(path).read_text(errors="replace")[:limit]
    except OSError:
        return None


def _command(args: list[str], timeout: int = 10) -> dict:
    try:
        run = subprocess.run(args, capture_output=True, text=True, timeout=timeout)
        return {"argv": args, "exit_code": run.returncode,
                "stdout": run.stdout[:30_000], "stderr": run.stderr[:10_000]}
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"argv": args, "error": str(exc)}


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _cpu_list(spec: str | None) -> list[int]:
    out = []
    for token in (spec or "").strip().split(","):
        if not token:
            continue
        bounds = token.split("-")
        out.extend(range(int(bounds[0]), int(bounds[-1]) + 1))
    return sorted(set(out))


def _node(cpu: int) -> int | None:
    found = sorted(Path(f"/sys/devices/system/cpu/cpu{cpu}").glob("node[0-9]*"))
    return int(found[0].name[4:]) if found else None


def _cpu_selection() -> dict:
    allowed = sorted(os.sched_getaffinity(0))
    by_node: dict[int | None, list[int]] = {}
    for cpu in allowed:
        by_node.setdefault(_node(cpu), []).append(cpu)
    choices = []
    for node, cpus in by_node.items():
        seen = set()
        physical = []
        topology_known = True
        for cpu in cpus:
            sibling_text = _read(f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list")
            if sibling_text is None:
                topology_known = False
            siblings = tuple(_cpu_list(sibling_text)) if sibling_text else (cpu,)
            if siblings not in seen:
                seen.add(siblings)
                physical.append(cpu)
        choices.append((len(physical), node, physical, topology_known))
    choices.sort(key=lambda x: x[0], reverse=True)
    _, node, physical, topology_known = choices[0]
    selected = physical[: WORKERS + 1]
    spare = [cpu for cpu in allowed if cpu not in selected]
    affinity_probe = _command(["taskset", "-c", ",".join(map(str, selected)), "true"])
    return {"allowed": allowed, "selected": selected, "sampler_cpu": spare[0] if spare else None,
            "node": node, "distinct_physical_cores": topology_known and len(selected) == WORKERS + 1,
            "topology_known": topology_known, "node_physical_cores_available": len(physical),
            "affinity_probe": affinity_probe, "affinity_supported": affinity_probe.get("exit_code") == 0,
            "nodes": by_node}


def _psi() -> dict:
    return {name: _read(f"/proc/pressure/{name}") for name in ("cpu", "memory", "io")}


def _cgroup() -> dict:
    names = (
        "cpuset.cpus", "cpuset.cpus.effective", "cpuset.cpus.exclusive.effective",
        "cpuset.cpus.partition", "cpuset.mems", "cpuset.mems.effective",
        "cpu.max", "cpu.stat", "memory.max", "memory.peak", "memory.events",
        "memory.numa_stat", "cpu.pressure", "memory.pressure", "io.pressure",
    )
    return {name: _read(f"/sys/fs/cgroup/{name}") for name in names}


def _interrupt_counts(cpus: list[int]) -> dict:
    raw = _read("/proc/interrupts", 2_000_000)
    if not raw:
        return {"available": False}
    lines = raw.splitlines()
    headers = [int(x[3:]) for x in re.findall(r"CPU\d+", lines[0])]
    indexes = {cpu: headers.index(cpu) for cpu in cpus if cpu in headers}
    totals = {str(cpu): 0 for cpu in indexes}
    for line in lines[1:]:
        words = line.split()
        for cpu, index in indexes.items():
            if index + 1 < len(words) and words[index + 1].isdigit():
                totals[str(cpu)] += int(words[index + 1])
    return {"available": bool(indexes), "per_cpu_total": totals,
            "raw_sha256": hashlib.sha256(raw.encode()).hexdigest()}


def _steal_ticks(cpus: list[int]) -> dict:
    raw = _read("/proc/stat") or ""
    wanted = {f"cpu{cpu}" for cpu in cpus}
    return {line.split()[0]: int(line.split()[8]) for line in raw.splitlines()
            if line.split() and line.split()[0] in wanted and len(line.split()) > 8}


def _counter_snapshot(cpus: list[int]) -> dict:
    return {"monotonic_ns": time.monotonic_ns(), "psi": _psi(), "cgroup": _cgroup(),
            "interrupts": _interrupt_counts(cpus), "steal_ticks": _steal_ticks(cpus)}


def _host() -> dict:
    selection = _cpu_selection()
    return {
        "captured_unix_ns": time.time_ns(), "hostname": platform.node(),
        "platform": platform.platform(), "uname": list(os.uname()),
        "modal_image_id": os.environ.get("MODAL_IMAGE_ID"),
        "modal_cloud_provider": os.environ.get("MODAL_CLOUD_PROVIDER"),
        "modal_region": os.environ.get("MODAL_REGION"),
        "selection": selection, "lscpu": _command(["lscpu", "-J"]),
        "numactl_hardware": _command(["numactl", "--hardware"]),
        "cpuinfo": _read("/proc/cpuinfo", 300_000),
        "kernel_command_line": _read("/proc/cmdline"),
        "self_cgroup": _read("/proc/self/cgroup"),
        "perf_event_paranoid": _read("/proc/sys/kernel/perf_event_paranoid"),
        "perf_version": _command(["perf", "--version"]),
        "cgroup": _cgroup(), "psi": _psi(),
        "interrupts": _interrupt_counts(selection["selected"]),
        "steal_ticks": _steal_ticks(selection["selected"]),
        "rustc": _command(["rustc", "--version", "--verbose"]),
        "cargo": _command(["cargo", "--version"]),
        "source_commit": SOURCE_COMMIT,
        "source_sha256": {
            "baseline": _hash(REMOTE_ROOT / "experiments/koblitz-s3-pair-query-20261007-v4/baseline.rs"),
            "candidate": _hash(REMOTE_ROOT / "experiments/koblitz-s3-pair-query-20261007-v4/candidate.rs"),
        },
        "binary_sha256": {
            arm: _hash(BIN / f"s3_pair_{arm}_frozen")
            for arm in ("baseline", "candidate")
        },
    }


def _phase_check(report: dict) -> bool:
    timing = report["timing_ms"]
    keys = ("target_query", "target_pdp_charged", "target_relation_check",
            "target_descent", "target_recovery_check")
    return all(isinstance(timing.get(key), (int, float)) and timing[key] >= 0 for key in keys) and abs(
        sum(timing[key] for key in keys) - timing["target_online_after_reusable_setup"]
    ) < 0.02


def _thread_masks(pid: int) -> dict:
    masks = set()
    count = 0
    for status in Path(f"/proc/{pid}/task").glob("*/status"):
        text = _read(str(status), 15_000)
        if text:
            count += 1
            found = re.search(r"^Cpus_allowed_list:\s*(.+)$", text, re.M)
            if found:
                masks.add(found.group(1).strip())
    return {"count": count, "masks": sorted(masks)}


def _run_solver(arm: str, n: int, target_file: Path, scalar: int,
                selection: dict, membind: bool, tag: str) -> dict:
    binary = BIN / f"s3_pair_{arm}_frozen"
    binary_before = _hash(binary)
    point = json.loads(target_file.read_text())
    with tempfile.TemporaryDirectory(prefix="s3-modal-") as tmp:
        output = Path(tmp) / "solver.jsonl"
        stdout = Path(tmp) / "stdout"
        stderr = Path(tmp) / "stderr"
        argv = (["taskset", "-c", ",".join(map(str, selection["selected"]))]
                if selection["affinity_supported"] else [])
        if membind:
            argv += ["numactl", f"--membind={selection['node']}"]
        argv += [str(binary), str(n), "0", str(COLUMNS), str(RELATION_SEED),
                 str(target_file), str(output), str(WORKERS)]
        before = _counter_snapshot(selection["selected"])
        started = time.time_ns()
        max_threads = 0
        masks = set()
        with stdout.open("w") as out, stderr.open("w") as err:
            process = subprocess.Popen(argv, stdout=out, stderr=err)
            while process.poll() is None:
                sample = _thread_masks(process.pid)
                max_threads = max(max_threads, sample["count"])
                masks.update(sample["masks"])
                time.sleep(0.05)
            exit_code = process.wait()
        finished = time.time_ns()
        after = _counter_snapshot(selection["selected"])
        raw = output.read_text() if output.exists() else None
        stdout_text = stdout.read_text()
        report = json.loads(raw) if raw else None
        binary_after = _hash(binary)
        stdout_matches_raw = stdout_text == raw
        valid = bool(
            exit_code == 0 and binary_before == binary_after and stdout_matches_raw
            and report and report.get("n") == n and report.get("a") == 0
            and report.get("target") == point and report.get("target_count") == 1
            and report.get("relation_seed") == RELATION_SEED
            and report.get("orbit_columns") == COLUMNS
            and report.get("relation_collection_workers") == WORKERS
            and report.get("group_verified") is True
            and report.get("recovered_scalar") == scalar and _phase_check(report)
        )
        return {
            "tag": tag, "arm": arm, "n": n, "public_point": point,
            "target_file_sha256": _hash(target_file),
            "binary_sha256_before": binary_before,
            "binary_sha256_after": binary_after, "argv": argv,
            "started_unix_ns": started, "finished_unix_ns": finished,
            "outer_wall_ns": finished - started, "exit_code": exit_code,
            "verified_native_and_fixture": valid, "thread_masks_seen": sorted(masks),
            "max_threads_seen": max_threads, "before": before, "after": after,
            "stdout_sha256": _hash(stdout), "stdout_matches_raw": stdout_matches_raw,
            "stdout_text_if_different": None if stdout_matches_raw else stdout_text,
            "stderr": stderr.read_text()[:20_000],
            "raw_sha256": hashlib.sha256(raw.encode()).hexdigest() if raw else None,
            "raw_text": raw,
            "report": report,
        }


def _perf_probe(binary: Path, n: int, target_file: Path, selection: dict) -> dict:
    with tempfile.TemporaryDirectory(prefix="s3-perf-") as tmp:
        output = Path(tmp) / "solver.jsonl"
        argv = (["taskset", "-c", ",".join(map(str, selection["selected"]))]
                if selection["affinity_supported"] else [])
        argv += ["perf", "stat", "-x", ",", "-e",
                "task-clock,cycles,instructions,context-switches,cpu-migrations,page-faults,branches,cache-misses",
                "--", str(binary), str(n), "0", str(COLUMNS), str(RELATION_SEED),
                str(target_file), str(output), str(WORKERS)]
        result = _command(argv, timeout=180)
        result["scope"] = "whole process including reusable setup; separate diagnostic pass"
        result["solver_raw_sha256"] = _hash(output) if output.exists() else None
        return result


def _semantic_witness(report: dict) -> tuple:
    return (
        report["factor_base_digest"], report["factor_base_points"], report["orbit_columns"],
        report["rank"], report["rank_attempts"], report["rank_new_rows"],
        report["target_relation_indices"], report["target_span_stop_relation_prefix"],
        report["recovered_scalar"],
        [(w["point_indices"], w["rank_gain"], w["relation_scalar"], w["x_codes"])
         for w in report["rank_relation_witnesses"]],
    )


def _p95(values: list[float]) -> float:
    return sorted(values)[math.ceil(0.95 * len(values)) - 1]


def _panel(count: int, aa: bool) -> dict:
    host = _host()
    selection = host["selection"]
    if len(selection["selected"]) < WORKERS + 1:
        raise RuntimeError("Modal container did not expose 15 logical CPUs")
    if selection["sampler_cpu"] is not None:
        try:
            os.sched_setaffinity(0, {selection["sampler_cpu"]})
            host["sampler_affinity_set"] = True
        except OSError as exc:
            host["sampler_affinity_set"] = False
            host["sampler_affinity_error"] = str(exc)
    bind_probe = (_command(["numactl", f"--membind={selection['node']}", "true"])
                  if selection["node"] is not None else {"error": "NUMA node unknown"})
    membind = bind_probe.get("exit_code") == 0
    runs = []
    blocks = []
    # Warmups use a separately derived public point outside the measured pilot.
    for n in (41, 53):
        target = WARMUP / f"n{n}/T001/public_target.json"
        fixture = json.loads((WARMUP / f"n{n}/fixtures.json").read_text())["fixtures"][0]
        runs.append(_run_solver("baseline", n, target, int(fixture["fixture_scalar"]),
                                selection, membind, f"warmup-n{n}"))
    for index in range(1, count + 1):
        for n in ((41, 53) if index % 2 else (53, 41)):
            fixture_doc = json.loads((PILOT / f"n{n}/fixtures.json").read_text())
            fixture = fixture_doc["fixtures"][index - 1]
            target = PILOT / f"n{n}/T{index:03}/public_target.json"
            scalar = int(fixture["fixture_scalar"])
            arms = ("baseline", "candidate", "candidate", "baseline") if index % 2 else (
                "candidate", "baseline", "baseline", "candidate")
            block_runs = []
            for position, arm in enumerate(arms):
                result = _run_solver(arm, n, target, scalar, selection, membind,
                                     f"n{n}-T{index:03}-pair-{position + 1}")
                runs.append(result)
                block_runs.append(result)
            if aa:
                for position in range(2):
                    result = _run_solver("baseline", n, target, scalar, selection, membind,
                                         f"n{n}-T{index:03}-aa-{position + 1}")
                    runs.append(result)
                    block_runs.append(result)
            good = all(run["verified_native_and_fixture"] for run in block_runs)
            semantic_equal = good and _semantic_witness(block_runs[0]["report"]) == _semantic_witness(
                block_runs[1]["report"])
            if semantic_equal:
                semantic_equal = all(_semantic_witness(run["report"]) == _semantic_witness(block_runs[0]["report"])
                                     for run in block_runs[:4])
            baseline = [run["report"]["timing_ms"]["target_online_after_reusable_setup"]
                        for run in block_runs[:4] if run["arm"] == "baseline" and good]
            candidate = [run["report"]["timing_ms"]["target_online_after_reusable_setup"]
                         for run in block_runs[:4] if run["arm"] == "candidate" and good]
            ratio = math.sqrt(baseline[0] * baseline[1] / (candidate[0] * candidate[1])) if good else None
            aa_noise = (abs(math.log(block_runs[4]["report"]["timing_ms"]["target_online_after_reusable_setup"] /
                                     block_runs[5]["report"]["timing_ms"]["target_online_after_reusable_setup"]))
                        if aa and good else None)
            blocks.append({"n": n, "label": fixture["label"], "point": fixture["public_point"],
                           "scalar_fixture": scalar, "order": arms,
                           "all_native_verified": good, "semantic_equal": semantic_equal,
                           "exploratory_online_ratio": ratio, "aa_abs_log_ratio": aa_noise,
                           "run_tags": [run["tag"] for run in block_runs]})
            print(f"n{n} {fixture['label']} verified={good} semantic={semantic_equal} ratio={ratio}", flush=True)
    perf = []
    for n in (41, 53):
        for arm in ("baseline", "candidate"):
            perf.append({"n": n, "arm": arm, "result": _perf_probe(
                BIN / f"s3_pair_{arm}_frozen", n,
                PILOT / f"n{n}/T001/public_target.json", selection)})
    summary = {}
    for n in (41, 53):
        curve_blocks = [b for b in blocks if b["n"] == n]
        ratios = [b["exploratory_online_ratio"] for b in curve_blocks if b["exploratory_online_ratio"]]
        aa_values = [b["aa_abs_log_ratio"] for b in curve_blocks if b["aa_abs_log_ratio"] is not None]
        summary[str(n)] = {
            "targets": len(curve_blocks), "all_verified": all(b["all_native_verified"] for b in curve_blocks),
            "all_semantic_equal": all(b["semantic_equal"] for b in curve_blocks),
            "geometric_mean_ratio": math.exp(statistics.mean(map(math.log, ratios))) if ratios else None,
            "sample_sd_log_ratio": statistics.stdev(map(math.log, ratios)) if len(ratios) > 1 else None,
            "aa_p95_abs_log_ratio": _p95(aa_values) if aa_values else None,
            "aa_noise_gate_pass": _p95(aa_values) < math.log(1.05) if aa_values else None,
        }
    return {
        "kind": "modal_s3_two_curve_container_diagnostic_v1", "source_commit": SOURCE_COMMIT,
        "runner_sha256": _hash(Path(__file__)),
        "modal_cpu_request_physical_cores": CPU_REQUEST, "modal_memory_request_mib": MEMORY_MIB,
        "strict_host_isolation_receipt": None, "controlled_speedup": None,
        "claim_status": "exploratory: Modal container counters do not prove host-wide exclusive CPUs, IRQs, and NUMA",
        "host": host, "numa_bind_probe": bind_probe, "numa_bind_used": membind,
        "runs": runs, "blocks": blocks, "perf_whole_process": perf, "summary": summary,
        "finished_unix_ns": time.time_ns(),
    }


def _encode(doc: dict) -> str:
    return base64.b64encode(gzip.compress(json.dumps(doc, separators=(",", ":")).encode())).decode()


def _save(encoded: str, output: str) -> None:
    path = Path(output)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = base64.b64decode(encoded)
    path.write_bytes(raw)
    doc = json.loads(gzip.decompress(raw))
    print(json.dumps({"artifact": str(path), "sha256": _hash(path),
                      "kind": doc["kind"], "summary": doc.get("summary")}, sort_keys=True))


@app.function(image=image, cpu=CPU_REQUEST, memory=MEMORY_MIB, timeout=120)
def remote_probe() -> str:
    return _encode({"kind": "modal_s3_host_probe_v1", "host": _host(),
                    "strict_host_isolation_receipt": None})


@app.function(image=image, cpu=CPU_REQUEST, memory=MEMORY_MIB, timeout=3600)
def remote_panel(count: int, aa: bool) -> str:
    if count not in (1, 12):
        raise ValueError("count must be 1 (smoke) or 12 (frozen pilot)")
    return _encode(_panel(count, aa))


@app.local_entrypoint()
def probe(output: str = "experiments/koblitz-s3-fidelity-20261007/modal/probe.json.gz"):
    _save(remote_probe.remote(), output)


@app.local_entrypoint()
def smoke(output: str = "experiments/koblitz-s3-fidelity-20261007/modal/smoke.json.gz"):
    _save(remote_panel.remote(1, False), output)


@app.local_entrypoint()
def pilot(output: str = "experiments/koblitz-s3-fidelity-20261007/modal/pilot.json.gz"):
    _save(remote_panel.remote(12, True), output)

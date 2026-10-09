"""Paired S3 prepared-log timing and correctness panel on a Modal VM.

The VM records guest topology, affinity, NUMA, PSI, cgroup, interrupt, steal,
and optional whole-process perf diagnostics. Host isolation remains a separate
gate. Run with ``python3 modal_prepared_panel.py --count 1|12 --output FILE``.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import uuid

import modal


ROOT = Path(__file__).resolve().parents[2]
REMOTE = Path("/root/crypto")
SOURCE_COMMIT = "833ef86191947ba8e4def56a78ee6c51a4ea0ff2"


def image() -> modal.Image:
    item = modal.Image.from_registry("rust:1.98-bookworm", add_python="3.12")
    item = item.apt_install("linux-perf", "numactl", "util-linux", "procps")
    for name in ("Cargo.toml", "Cargo.lock"):
        item = item.add_local_file(str(ROOT / name), str(REMOTE / name), copy=True)
    for name in ("src", "docs/ic", "docs/ecbench", "docs/curves"):
        item = item.add_local_dir(str(ROOT / name), remote_path=str(REMOTE / name), copy=True)
    for name in (
        "examples/s3_prepared_targetspan.rs",
        "experiments/koblitz-s3-pair-query-20261007-v4/baseline.rs",
        "experiments/koblitz-s3-pair-query-20261007-v4/candidate.rs",
    ):
        item = item.add_local_file(str(ROOT / name), str(REMOTE / name), copy=True)
    for name in (
        "experiments/koblitz-s3-fidelity-20261007/pilot",
        "experiments/koblitz-s3-prepared-20261009/warmup",
    ):
        item = item.add_local_dir(str(ROOT / name), remote_path=str(REMOTE / name), copy=True)
    return item.run_commands(
        f"cd {REMOTE} && cargo build --release --locked "
        "--example s3_pair_candidate_frozen --example s3_prepared_targetspan"
    )


REMOTE_CODE = r'''
import gzip, hashlib, json, math, os, platform, re, statistics, subprocess, sys, tempfile, time
from pathlib import Path

ROOT = Path("/root/crypto")
PILOT = ROOT / "experiments/koblitz-s3-fidelity-20261007/pilot"
WARMUP = ROOT / "experiments/koblitz-s3-prepared-20261009/warmup"
BIN = ROOT / "target/release/examples"
WORKERS = 14
PHASES = ("target_query", "target_pdp_charged", "target_relation_check",
          "target_descent", "target_recovery_check")

def read(path, limit=500000):
    try:
        return Path(path).read_text(errors="replace")[:limit]
    except OSError:
        return None

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def command(args, timeout=10):
    try:
        run = subprocess.run(args, capture_output=True, text=True, timeout=timeout)
        return {"argv": args, "exit_code": run.returncode,
                "stdout": run.stdout[:100000], "stderr": run.stderr[:20000]}
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"argv": args, "error": str(exc)}

def cpu_list(spec):
    result = []
    for token in (spec or "").strip().split(","):
        if token:
            bounds = token.split("-")
            result.extend(range(int(bounds[0]), int(bounds[-1]) + 1))
    return sorted(set(result))

def node_for(cpu):
    nodes = sorted(Path(f"/sys/devices/system/cpu/cpu{cpu}").glob("node[0-9]*"))
    return int(nodes[0].name[4:]) if nodes else None

def select_cpus():
    allowed = sorted(os.sched_getaffinity(0))
    by_node = {}
    for cpu in allowed:
        by_node.setdefault(node_for(cpu), []).append(cpu)
    choices = []
    for node, cpus in by_node.items():
        seen, physical = set(), []
        for cpu in cpus:
            siblings = tuple(cpu_list(read(f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list"))) or (cpu,)
            if siblings not in seen and cpu != 0:
                physical.append(cpu)
                seen.add(siblings)
        choices.append((len(physical), node, physical))
    choices.sort(key=lambda row: row[0], reverse=True)
    _, node, physical = choices[0]
    selected = physical[:WORKERS + 1]
    if len(selected) != WORKERS + 1:
        raise RuntimeError("VM did not expose 15 physical execution CPUs on one NUMA node")
    sampler = next((cpu for cpu in allowed if cpu not in selected), None)
    probe = command(["taskset", "-c", ",".join(map(str, selected)), "true"])
    return {"allowed": allowed, "selected": selected, "node": node,
            "sampler_cpu": sampler, "by_node": by_node, "affinity_probe": probe,
            "distinct_physical_cpus": len(selected)}

def interrupt_counts(cpus):
    raw = read("/proc/interrupts", 3000000)
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
    return {"available": bool(indexes), "totals": totals,
            "raw_sha256": hashlib.sha256(raw.encode()).hexdigest()}

def steal_ticks(cpus):
    wanted = {f"cpu{cpu}" for cpu in cpus}
    return {words[0]: int(words[8]) for line in (read("/proc/stat") or "").splitlines()
            if (words := line.split()) and words[0] in wanted and len(words) > 8}

def snapshot(cpus):
    return {"monotonic_ns": time.monotonic_ns(),
            "psi": {name: read(f"/proc/pressure/{name}") for name in ("cpu", "memory", "io")},
            "cgroup": {name: read(f"/sys/fs/cgroup/{name}") for name in
                       ("cpu.stat", "cpu.pressure", "memory.pressure", "io.pressure",
                        "memory.peak", "memory.events", "memory.numa_stat")},
            "interrupts": interrupt_counts(cpus), "steal_ticks": steal_ticks(cpus)}

def host(selection):
    cpu = selection["selected"][0]
    return {"platform": platform.platform(), "uname": list(os.uname()),
            "hostname": platform.node(), "selection": selection,
            "modal_image_id": os.environ.get("MODAL_IMAGE_ID"),
            "modal_cloud_provider": os.environ.get("MODAL_CLOUD_PROVIDER"),
            "modal_region": os.environ.get("MODAL_REGION"),
            "lscpu": command(["lscpu", "-J"]), "numactl_hardware": command(["numactl", "--hardware"]),
            "cpuinfo": read("/proc/cpuinfo", 300000), "kernel_command_line": read("/proc/cmdline"),
            "self_cgroup": read("/proc/self/cgroup"),
            "perf_event_paranoid": read("/proc/sys/kernel/perf_event_paranoid"),
            "perf_version": command(["perf", "--version"]),
            "cpufreq_governor": read(f"/sys/devices/system/cpu/cpu{cpu}/cpufreq/scaling_governor"),
            "cpufreq_min": read(f"/sys/devices/system/cpu/cpu{cpu}/cpufreq/scaling_min_freq"),
            "cpufreq_max": read(f"/sys/devices/system/cpu/cpu{cpu}/cpufreq/scaling_max_freq"),
            "cgroup": {name: read(f"/sys/fs/cgroup/{name}") for name in
                       ("cpuset.cpus", "cpuset.cpus.effective", "cpuset.cpus.exclusive.effective",
                        "cpuset.cpus.partition", "cpuset.mems", "cpuset.mems.effective",
                        "cpu.max", "cpu.stat", "memory.max", "memory.peak")},
            "initial_counters": snapshot(selection["selected"]),
            "rustc": command(["rustc", "--version", "--verbose"]),
            "cargo": command(["cargo", "--version"])}

def thread_masks(pid):
    masks = set()
    count = 0
    for status in Path(f"/proc/{pid}/task").glob("*/status"):
        value = read(status, 15000)
        if value:
            count += 1
            match = re.search(r"^Cpus_allowed_list:\s*(.+)$", value, re.M)
            if match:
                masks.add(match.group(1).strip())
    return count, masks

def run_solver(arm, n, target, scalar, selection, membind, tag):
    binary = BIN / ("s3_prepared_targetspan" if arm == "prepared" else "s3_pair_candidate_frozen")
    point = json.loads(target.read_text())
    with tempfile.TemporaryDirectory(prefix="s3-prepared-") as directory:
        output = Path(directory) / "raw.jsonl"
        stdout = Path(directory) / "stdout"
        stderr = Path(directory) / "stderr"
        argv = ["taskset", "-c", ",".join(map(str, selection["selected"]))]
        if membind:
            argv += ["numactl", f"--membind={selection['node']}"]
        argv += [str(binary), str(n), "0", "244", "20260928", str(target), str(output), str(WORKERS)]
        before = snapshot(selection["selected"])
        started = time.time_ns()
        masks, max_threads, error = set(), 0, None
        with stdout.open("w") as out, stderr.open("w") as err:
            process = subprocess.Popen(argv, stdout=out, stderr=err)
            while process.poll() is None:
                count, current = thread_masks(process.pid)
                max_threads = max(max_threads, count)
                masks.update(current)
                if time.time_ns() - started > 180_000_000_000:
                    process.kill()
                    error = "timeout_180s"
                    break
                time.sleep(0.02)
            exit_code = process.wait()
        finished = time.time_ns()
        after = snapshot(selection["selected"])
        raw = output.read_text() if output.exists() else None
        stdout_text = stdout.read_text()
        try:
            report = json.loads(raw) if raw else None
        except json.JSONDecodeError:
            report = None
            error = "invalid_json"
        times = report.get("timing_ms", {}) if report else {}
        phase_ok = bool(report and all(isinstance(times.get(key), (int, float)) and times[key] >= 0
                                       for key in PHASES) and
                        abs(sum(times[key] for key in PHASES) -
                            times["target_online_after_reusable_setup"]) < 0.02)
        verified = bool(exit_code == 0 and report and report.get("group_verified") is True
                        and report.get("recovered_scalar") == scalar and report.get("target") == point
                        and report.get("n") == n and report.get("a") == 0
                        and report.get("target_count") == 1 and report.get("relation_seed") == 20260928
                        and report.get("factor_base_points") is not None
                        and report.get("relation_collection_workers") == WORKERS
                        and (arm != "prepared" or report.get("rank_complete") is True)
                        and phase_ok and stdout_text == raw)
        return {"tag": tag, "n": n, "arm": arm, "point": point, "fixture_scalar": scalar,
                "target_sha256": sha(target), "binary_sha256": sha(binary), "argv": argv,
                "started_unix_ns": started, "finished_unix_ns": finished,
                "outer_wall_ns": finished - started, "exit_code": exit_code, "error": error,
                "verified_native_and_fixture": verified, "phase_ok": phase_ok,
                "max_threads_seen": max_threads, "thread_masks_seen": sorted(masks),
                "before": before, "after": after, "stdout_sha256": sha(stdout),
                "stdout_matches_raw": stdout_text == raw,
                "stderr": stderr.read_text()[:20000],
                "raw_sha256": hashlib.sha256(raw.encode()).hexdigest() if raw else None,
                "raw_text": raw}

def perf_probe(arm, n, target, selection, membind):
    binary = BIN / ("s3_prepared_targetspan" if arm == "prepared" else "s3_pair_candidate_frozen")
    with tempfile.TemporaryDirectory(prefix="s3-prepared-perf-") as directory:
        output = Path(directory) / "raw.jsonl"
        argv = ["taskset", "-c", ",".join(map(str, selection["selected"]))]
        if membind:
            argv += ["numactl", f"--membind={selection['node']}"]
        argv += ["perf", "stat", "-x", ",", "-e",
                 "cycles,instructions,task-clock,context-switches,cpu-migrations,page-faults,branches,cache-misses",
                 "--", str(binary), str(n), "0", "244", "20260928", str(target), str(output), str(WORKERS)]
        result = command(argv, timeout=180)
        result["scope"] = "whole process including reusable preparation and one target"
        result["raw_sha256"] = sha(output) if output.exists() else None
        return result

def main(count):
    selection = select_cpus()
    if selection["sampler_cpu"] is not None:
        os.sched_setaffinity(0, {selection["sampler_cpu"]})
    bind_probe = command(["numactl", f"--membind={selection['node']}", "true"])
    membind = bind_probe.get("exit_code") == 0
    runs, blocks = [], []
    for n in (41, 53):
        fixture = json.loads((WARMUP / f"n{n}/fixtures.json").read_text())["fixtures"][0]
        runs.append(run_solver("prepared", n, WARMUP / f"n{n}/T001/public_target.json",
                               int(fixture["fixture_scalar"]), selection, membind, f"warmup-n{n}"))
    for index in range(1, count + 1):
        for n in ((41, 53) if index % 2 else (53, 41)):
            fixture = json.loads((PILOT / f"n{n}/fixtures.json").read_text())["fixtures"][index - 1]
            point = PILOT / f"n{n}/T{index:03}/public_target.json"
            scalar = int(fixture["fixture_scalar"])
            arms = ("candidate", "prepared", "prepared", "candidate") if index % 2 else (
                   "prepared", "candidate", "candidate", "prepared")
            block = []
            for position, arm in enumerate(arms + ("prepared", "prepared")):
                row = run_solver(arm, n, point, scalar, selection, membind,
                                 f"n{n}-T{index:03}-{position + 1}")
                runs.append(row)
                block.append(row)
            good = all(row["verified_native_and_fixture"] for row in block)
            reports = [json.loads(row["raw_text"]) for row in block] if good else []
            by_arm = {arm: [report["timing_ms"]["target_online_after_reusable_setup"]
                            for row, report in zip(block[:4], reports[:4]) if row["arm"] == arm]
                      for arm in ("candidate", "prepared")} if good else {}
            ratio = math.sqrt(math.prod(by_arm["candidate"]) / math.prod(by_arm["prepared"])) if good else None
            aa = abs(math.log(reports[4]["timing_ms"]["target_online_after_reusable_setup"] /
                              reports[5]["timing_ms"]["target_online_after_reusable_setup"])) if good else None
            blocks.append({"n": n, "label": fixture["label"], "point": fixture["public_point"],
                           "order": arms, "run_tags": [row["tag"] for row in block],
                           "all_verified": good, "exploratory_online_ratio": ratio,
                           "prepared_aa_abs_log_ratio": aa})
            print(f"n{n} {fixture['label']} verified={good} ratio={ratio}", file=sys.stderr, flush=True)
    summary = {}
    for n in (41, 53):
        rows = [row for row in blocks if row["n"] == n]
        ratios = [row["exploratory_online_ratio"] for row in rows if row["exploratory_online_ratio"] is not None]
        aa = [row["prepared_aa_abs_log_ratio"] for row in rows if row["prepared_aa_abs_log_ratio"] is not None]
        summary[str(n)] = {"targets": len(rows), "all_verified": all(row["all_verified"] for row in rows),
                           "geometric_mean_ratio": math.exp(statistics.mean(map(math.log, ratios))) if ratios else None,
                           "min_target_ratio": min(ratios) if ratios else None,
                           "sample_sd_log_ratio": statistics.stdev(map(math.log, ratios)) if len(ratios) > 1 else None,
                           "prepared_aa_p95_abs_log_ratio": sorted(aa)[math.ceil(.95 * len(aa)) - 1] if aa else None}
    perf = [{"n": n, "arm": arm,
             "result": perf_probe(arm, n, PILOT / f"n{n}/T001/public_target.json", selection, membind)}
            for n in (41, 53) for arm in ("candidate", "prepared")]
    doc = {"kind": "modal_vm_s3_prepared_single_target_panel_v1", "source_commit": "833ef86191947ba8e4def56a78ee6c51a4ea0ff2",
           "source_sha256": {"candidate": sha(ROOT / "experiments/koblitz-s3-pair-query-20261007-v4/candidate.rs"),
                             "prepared": sha(ROOT / "examples/s3_prepared_targetspan.rs")},
           "binary_sha256": {"candidate": sha(BIN / "s3_pair_candidate_frozen"),
                             "prepared": sha(BIN / "s3_prepared_targetspan")},
           "host": host(selection), "numa_bind_probe": bind_probe, "numa_bind_used": membind,
           "strict_host_isolation_receipt": None, "controlled_speedup": None,
           "runs": runs, "blocks": blocks, "summary": summary, "perf_whole_process": perf,
           "finished_unix_ns": time.time_ns()}
    compressed = gzip.compress(json.dumps(doc, separators=(",", ":")).encode())
    Path("/results/panel.json.gz").write_bytes(compressed)
    print(json.dumps({"bytes": len(compressed), "sha256": hashlib.sha256(compressed).hexdigest(),
                      "runs": len(runs), "blocks": len(blocks), "summary": summary}), flush=True)

main(int(sys.argv[1]))
'''


def run(count: int, output: Path, region: str | None) -> None:
    if output.exists():
        raise FileExistsError(output)
    application = modal.App("s3-prepared-targetspan-vm-panel")
    volume_name = f"s3-prepared-panel-{uuid.uuid4().hex[:20]}"
    volume = modal.Volume.from_name(volume_name, create_if_missing=True)
    print(json.dumps({"volume_name": volume_name,
                      "remote_runner_sha256": hashlib.sha256(REMOTE_CODE.encode()).hexdigest()}), flush=True)
    with application.run():
        sandbox = modal.Sandbox.create(
            "python3", "-c", REMOTE_CODE, str(count), app=application,
            image=image(), cpu=16, memory=8192, timeout=1800, region=region,
            volumes={"/results": volume}, experimental_options={"vm_runtime": True},
        )
        print(json.dumps({"sandbox_id": sandbox.object_id}), flush=True)
        try:
            sandbox.wait()
            stdout = sandbox.stdout.read().strip()
            stderr = sandbox.stderr.read()
            if not stdout:
                raise RuntimeError("VM emitted no summary: " + stderr[-6000:])
            expected = json.loads(stdout)
        finally:
            sandbox.terminate()
        raw = b"".join(volume.read_file("panel.json.gz"))
        digest = hashlib.sha256(raw).hexdigest()
        if len(raw) != expected["bytes"] or digest != expected["sha256"]:
            raise RuntimeError("Modal Volume transfer failed size or SHA-256 validation")
        doc = json.loads(gzip.decompress(raw))
        if len(doc["runs"]) != 2 + count * 2 * 6 or len(doc["blocks"]) != count * 2:
            raise RuntimeError("incomplete raw run ledger")
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_bytes(raw)
        print(json.dumps({"output": str(output), "sha256": digest,
                          "runs": len(doc["runs"]), "blocks": len(doc["blocks"]),
                          "summary": doc["summary"]}, sort_keys=True), flush=True)
        volume.remove_file("panel.json.gz")
    modal.Volume.delete(volume_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--count", type=int, choices=(1, 12), default=1)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--region", default=None)
    args = parser.parse_args()
    run(args.count, args.output, args.region)

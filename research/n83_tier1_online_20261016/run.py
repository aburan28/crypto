#!/usr/bin/env python3
"""Five fresh same-point IC/rho pairs on one bounded N83 cloud allocation."""

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import urllib.request

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
DYNAMIC_RUNNER = ROOT / "research/n83_root11_dynamic_20261016/run.py"
SOURCE_SHA256 = "3ee6bc11399b382e72789b152043e5a67abc5b783c4e998ad8c534d60baebc03"
TARGETS_SHA256 = "4976679ac6c9e3858cd26312ec2e50b8eb30971041ea38286c9582dd9217a1c1"
ANSWERS_SHA256 = "41004f94d7257c26a434a291767a06b87d41d590d952e5551edca1aa62c5545d"
CALIBRATION_SEED = "2026101601"
RHO_SEED_BASE = 2026101700
IC_CAP_SECONDS = 900
RHO_CAP_SECONDS = 120
EXPECTED_AZ = "us-east-1b"
EXPECTED_TYPE = "g5g.metal"

spec = importlib.util.spec_from_file_location("n83_tier1_dynamic_runner", DYNAMIC_RUNNER)
dynamic = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dynamic)
cold = dynamic.cold
parent = dynamic.parent
dynamic.SOURCE_SHA256 = SOURCE_SHA256
parent.EXPECTED[parent.SOURCE] = SOURCE_SHA256


def one_json(path):
    return json.loads(path.read_text())


def fixture_preflight():
    hashes = dynamic.frozen_inputs()
    manifest = one_json(HERE / "fixture_manifest.json")
    if (parent.digest(HERE / "targets.jsonl") != TARGETS_SHA256
            or parent.digest(HERE / "answers.jsonl") != ANSWERS_SHA256
            or manifest["targets_sha256"] != TARGETS_SHA256
            or manifest["answers_sha256"] != ANSWERS_SHA256
            or manifest["target_count"] != 5
            or manifest["domain"] != "n83-tier1-online-five-fresh-20261016-v1"
            or parent.digest(dynamic.SOURCE) != SOURCE_SHA256):
        raise ValueError("frozen Tier 1 fixture or merged IC source differs")
    points = [json.loads(line) for line in
              (HERE / "targets.jsonl").read_text().splitlines()]
    if len(points) != 5 or len({tuple(point) for point in points}) != 5:
        raise ValueError("five distinct fresh public points required")
    hashes.update({
        str((HERE / "targets.jsonl").relative_to(ROOT)): TARGETS_SHA256,
        str((HERE / "answers.jsonl").relative_to(ROOT)): ANSWERS_SHA256,
        str((HERE / "fixture_manifest.json").relative_to(ROOT)):
            parent.digest(HERE / "fixture_manifest.json"),
    })
    return hashes, points


def selection_preflight():
    selection_path = HERE / "selection.json"
    if not selection_path.exists():
        raise ValueError("versioned candidate selection is required before panel launch")
    selection = one_json(selection_path)
    k = selection.get("k")
    layout = selection.get("layout")
    candidate_path = ROOT / selection.get("candidate_manifest", "")
    if (selection.get("schema") != "n83.tier1-online-selection/v1"
            or selection.get("candidate_id") is None
            or not isinstance(k, int) or k <= 0
            or layout not in ("root22_dynamic", "root11_dynamic", "root11_fixed")
            or selection.get("source_sha256") != SOURCE_SHA256
            or selection.get("targets_sha256") != TARGETS_SHA256
            or not candidate_path.is_file()
            or parent.digest(candidate_path) != selection.get("candidate_manifest_sha256")):
        raise ValueError("candidate identity, source or selection record differs")
    candidate = one_json(candidate_path)
    body = {key: value for key, value in candidate.items()
            if key != "candidate_id"}
    digest = hashlib.sha256(json.dumps(body, sort_keys=True,
                                       separators=(",", ":"),
                                       ensure_ascii=False).encode()).hexdigest()[:12]
    expected_id = (f"IC1N83Ckb1a1fb{166 * k}PDP4rootRCguidedLAgaussTDdirectISO0h"
                   + digest)
    flags = candidate.get("implementation", {}).get("algorithm_flags", {})
    if (candidate.get("candidate_id") != selection["candidate_id"]
            or candidate["candidate_id"] != expected_id
            or candidate.get("factor_base", {}).get("actual_usable_points") != 166 * k
            or candidate["factor_base"].get("effective_columns") != k
            or candidate["implementation"]["source_sha256"]["ic"] != SOURCE_SHA256
            or flags.get("root_only22") != (layout == "root22_dynamic")
            or flags.get("root_only11") != (layout != "root22_dynamic")
            or flags.get("indexed5_load75_request") != (layout != "root11_fixed")
            or flags.get("indexed5_fixed30_request") != (layout == "root11_fixed")
            or flags.get("root_bloom") is not False):
        raise ValueError("selected candidate manifest does not bind the base")
    return selection, candidate_path


def imds_identity():
    token_request = urllib.request.Request(
        "http://169.254.169.254/latest/api/token", data=b"", method="PUT",
        headers={"X-aws-ec2-metadata-token-ttl-seconds": "21600"})
    with urllib.request.urlopen(token_request, timeout=3) as response:
        token = response.read().decode()
    request = urllib.request.Request(
        "http://169.254.169.254/latest/dynamic/instance-identity/document",
        headers={"X-aws-ec2-metadata-token": token})
    with urllib.request.urlopen(request, timeout=3) as response:
        return json.load(response)


def github_identity():
    if (os.environ.get("GITHUB_ACTIONS") != "true"
            or os.environ.get("GITHUB_REPOSITORY") != "aburan28/crypto"
            or os.environ.get("RUNNER_OS") != "Linux"):
        raise ValueError("GitHub-hosted Tier 1 identity is unavailable")
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"],
                                     cwd=ROOT, text=True).strip()
    run_id = os.environ.get("GITHUB_RUN_ID", "")
    attempt = os.environ.get("GITHUB_RUN_ATTEMPT", "")
    if (not run_id.isdecimal() or not attempt.isdecimal()
            or os.environ.get("GITHUB_SHA") != commit
            or not os.environ.get("RUNNER_NAME")):
        raise ValueError("GitHub run, attempt, commit or runner differs")
    return {"provider": "github_actions", "repository": "aburan28/crypto",
            "runId": int(run_id), "runAttempt": int(attempt),
            "sha": commit, "runnerName": os.environ["RUNNER_NAME"],
            "runnerOS": "Linux",
            "instanceId": f"github-actions-{run_id}-{attempt}"}


def github_seal_preflight(selection, identity):
    path = HERE / "github_artifact_seal.json"
    if not path.exists():
        raise ValueError("GitHub selector artifact is not sealed")
    seal = one_json(path)
    if (seal.get("schema") != "n83.tier1-github-selector-seal/v1"
            or seal.get("artifact_id", 0) <= 0
            or seal.get("artifact_digest") is None
            or seal.get("selection_sha256") != parent.digest(
                HERE / "selection.json")
            or seal.get("candidate_manifest_sha256") != selection[
                "candidate_manifest_sha256"]
            or seal.get("selector_host_sha256") != selection[
                "screen_host_sha256"]
            or len(seal.get("selector_cpus", [])) != 2
            or seal.get("run_id") != identity["runId"]
            or seal.get("run_attempt") != identity["runAttempt"]
            or seal.get("head_sha") != identity["sha"]
            or seal.get("repository") != identity["repository"]):
        raise ValueError("GitHub selector artifact receipt differs")
    return seal


def cpu_list(value):
    result = set()
    for part in value.strip().split(","):
        if not part:
            continue
        ends = part.split("-", 1)
        result.update(range(int(ends[0]), int(ends[-1]) + 1))
    return result


def choose_cpus(count=4):
    available = set(os.sched_getaffinity(0))
    nodes = sorted(Path("/sys/devices/system/node").glob("node[0-9]*"))
    for node in nodes:
        cpulist = parent.read_text(node / "cpulist")
        if not cpulist:
            continue
        candidates = sorted(available & cpu_list(cpulist))
        chosen = []
        occupied_siblings = set()
        for cpu in candidates:
            if cpu in occupied_siblings:
                continue
            siblings = parent.read_text(Path(
                f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list"))
            chosen.append(cpu)
            occupied_siblings.update(cpu_list(siblings or str(cpu)))
            if len(chosen) == count:
                return chosen, node.name
    raise ValueError(f"{count} distinct physical cores on one visible NUMA node required")


def host_telemetry(cpus):
    stat = parent.read_text(Path("/proc/stat")) or ""
    lines = {line.split()[0]: line.split()[1:] for line in stat.splitlines()
             if line.startswith("cpu") or line.startswith("intr ")}
    steal = sum(int(lines[f"cpu{cpu}"][7]) for cpu in cpus)
    interrupts = int(lines.get("intr", ["0"])[0])
    freq = {}
    for cpu in cpus:
        path = Path(f"/sys/devices/system/cpu/cpu{cpu}/cpufreq/scaling_cur_freq")
        value = parent.read_text(path)
        freq[str(cpu)] = int(value.strip()) if value and value.strip().isdigit() else None
    return {"recorded_utc": parent.now(), "steal_ticks": steal,
            "interrupt_count": interrupts, "frequency_khz": freq}


def one_run(output, name, command, environment, cap, group, cpus):
    folder = output / "runs" / name
    before = host_telemetry(cpus)
    usage = parent.one_run(folder, command, environment, cap, group)
    after = host_telemetry(cpus)
    ticks = os.sysconf("SC_CLK_TCK")
    steal_seconds = (after["steal_ticks"] - before["steal_ticks"]) / ticks
    telemetry = {"before": before, "after": after,
                 "steal_seconds": steal_seconds,
                 "steal_fraction_of_wall": steal_seconds / usage["wall_seconds"],
                 "steal_fraction_of_cpu_envelope":
                     steal_seconds / (usage["wall_seconds"] * len(cpus)),
                 "interrupt_delta": after["interrupt_count"] - before["interrupt_count"],
                 "selected_cpus": cpus}
    parent.write_json(folder / "telemetry.json", telemetry)
    print(f"{name}: {usage['status']} after {usage['wall_seconds']:.3f}s",
          flush=True)
    return usage


def make_inputs(output, points):
    directory = output / "inputs"
    directory.mkdir()
    result = []
    for index, point in enumerate(points, start=1):
        path = directory / f"point_{index}.jsonl"
        path.write_text(json.dumps(point, separators=(",", ":")) + "\n")
        result.append(path)
    return result


def run(binary, rho_binary, output, github_host=False):
    if output.exists():
        raise ValueError("output directory already exists")
    if github_host:
        identity = github_identity()
        selection, candidate_path = selection_preflight()
        github_seal_preflight(selection, identity)
        hashes, points = fixture_preflight()
    else:
        identity = imds_identity()
        if (identity.get("availabilityZone") != EXPECTED_AZ
                or identity.get("instanceType") != EXPECTED_TYPE):
            raise ValueError("AWS allocation type or zone differs")
        hashes, points = fixture_preflight()
        selection, candidate_path = selection_preflight()
    cpus, node = choose_cpus(2 if github_host else 4)
    if github_host and cpus != github_seal_preflight(selection, identity)[
            "selector_cpus"]:
        raise ValueError("selector and panel CPU allocation differs")
    os.sched_setaffinity(0, cpus)
    output.mkdir(parents=True)
    group = parent.cgroup_root()
    host = parent.host_receipt(group)
    topology = subprocess.run(["lscpu", "-J"], capture_output=True, text=True,
                              check=False)
    selected_topology = {}
    for cpu in cpus:
        root = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        selected_topology[str(cpu)] = {
            name: parent.read_text(root / name) for name in
            ("core_id", "physical_package_id", "thread_siblings_list")}
    provider = "github_actions" if github_host else "aws_ec2_spot"
    host.update({"schema": "n83.tier1-online-host/v1", "provider": provider,
                 "instance_identity": identity, "instance_id": identity["instanceId"],
                 "cpu_affinity": cpus, "numa_node": node,
                 "selected_cpu_topology": selected_topology,
                 "lscpu_json": json.loads(topology.stdout) if topology.returncode == 0
                 else None,
                 "lscpu_stderr": topology.stderr,
                 "cpuinfo_sha256": parent.digest(Path("/proc/cpuinfo")),
                 "clock_ticks_per_second": os.sysconf("SC_CLK_TCK"),
                 "recorded_start_utc": parent.now()})
    host["run_url"] = (f"https://github.com/aburan28/crypto/actions/runs/"
                       f"{identity['runId']}" if github_host else None)
    parent.write_json(output / "host.json", host)
    memory_max = (host["cgroup_before"]["memory.max"] or "").strip()
    if (host["mem_total_kib"] is None
            or host["mem_total_kib"] * 1024 < cold.RSS_CAP_BYTES + 256 * 1024**2
            or host["disk_free_bytes"] < 2 * 1024**3
            or (memory_max.isdigit() and int(memory_max) < cold.RSS_CAP_BYTES)):
        parent.write_json(output / "preflight.json", {"status": "REJECTED_HOST_CAPACITY"})
        return 2
    parent.write_json(output / "preflight.json", {
        "status": "PASS_CAPACITY_PREFLIGHT", "rss_cap_bytes": cold.RSS_CAP_BYTES,
        "selected_cpus": cpus})
    inputs = make_inputs(output, points)
    schedule = [{"point_index": i, "pair_order": ["ic", "rho"] if i % 2 else
                 ["rho", "ic"]} for i in range(1, 6)]
    parent.write_json(output / "source.json", {
        "schema": "n83.tier1-online-source/v1", "hashes": hashes,
        "source_sha256": SOURCE_SHA256,
        "binary_sha256": parent.digest(binary),
        "rho_binary_sha256": parent.digest(rho_binary),
        "runner_sha256": parent.digest(Path(__file__)),
        "selection_sha256": parent.digest(HERE / "selection.json"),
        "candidate_manifest_sha256": parent.digest(candidate_path),
        "candidate_id": selection["candidate_id"],
        "point_input_sha256": {str(i): parent.digest(path)
                               for i, path in enumerate(inputs, start=1)},
        "schedule": schedule, "rho_seed_base": RHO_SEED_BASE,
        "calibration_seed": CALIBRATION_SEED,
    })
    for row in schedule:
        i = row["point_index"]
        calibrate = lambda suffix: one_run(
            output, f"p{i}_cal_{suffix}",
            [str(rho_binary), str(cold.TARGETS), CALIBRATION_SEED,
             str(cold.RHO_CONFIG),
             str(output / "runs" / f"p{i}_cal_{suffix}" / "receipt.json")],
            {"PATH": os.environ.get("PATH", "/usr/bin:/bin")},
            RHO_CAP_SECONDS, group, cpus)
        calibrate("before")
        for arm in row["pair_order"]:
            name = f"p{i}_{arm}"
            folder = output / "runs" / name
            if arm == "ic":
                command = [str(binary), f"construct:83:1:{selection['k']}",
                           str(inputs[i - 1]), "29", str(folder / "targets_result.jsonl")]
                environment = dynamic.environment(folder, selection["layout"])
                cap = IC_CAP_SECONDS
            else:
                command = [str(rho_binary), str(inputs[i - 1]),
                           str(RHO_SEED_BASE + i), str(cold.RHO_CONFIG),
                           str(folder / "receipt.json")]
                environment = {"PATH": os.environ.get("PATH", "/usr/bin:/bin")}
                cap = RHO_CAP_SECONDS
            one_run(output, name, command, environment, cap, group, cpus)
        calibrate("after")
    parent.write_json(output / "finish.json", {"finished_utc": parent.now(),
                                               "instance_id": identity["instanceId"]})
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fixture-preflight", action="store_true")
    parser.add_argument("--github-host", action="store_true")
    parser.add_argument("--binary", type=Path)
    parser.add_argument("--rho-binary", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.fixture_preflight:
        hashes, points = fixture_preflight()
        print(json.dumps({"hashes": hashes, "point_count": len(points)},
                         sort_keys=True))
        return 0
    if not args.binary or not args.rho_binary or not args.output:
        parser.error("--binary, --rho-binary and --output are required")
    return run(args.binary.resolve(), args.rho_binary.resolve(),
               args.output.resolve(), args.github_host)


if __name__ == "__main__":
    raise SystemExit(main())

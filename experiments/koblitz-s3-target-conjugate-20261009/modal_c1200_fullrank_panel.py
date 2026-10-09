"""Prepare or launch the frozen C1200 single-target panel on a Modal VM.

The default command only verifies and prints the local upload preview. The
--launch option creates one bounded VM sandbox and transfers the exact files
in MODAL_UPLOAD_MANIFEST.json. It must be used only after upload approval.
"""

from __future__ import annotations

import argparse
import ast
import gzip
import hashlib
import json
from pathlib import Path
import uuid


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MANIFEST = HERE / "MODAL_UPLOAD_MANIFEST.json"
TELEMETRY_SOURCE = ROOT / "experiments/koblitz-s3-prepared-20261009/modal_prepared_panel.py"
REMOTE = Path("/root/crypto")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def remote_code() -> str:
    """Reuse the previously measured guest telemetry functions unchanged."""
    tree = ast.parse(TELEMETRY_SOURCE.read_text())
    matches = [node for node in tree.body if isinstance(node, ast.Assign)
               and any(isinstance(target, ast.Name) and target.id == "REMOTE_CODE"
                       for target in node.targets)]
    assert len(matches) == 1
    previous = ast.literal_eval(matches[0].value)
    telemetry, delimiter, _ = previous.partition("\ndef run_solver(")
    assert delimiter and "def select_cpus():" in telemetry
    assert "def snapshot(cpus):" in telemetry and "def host(selection):" in telemetry
    return telemetry + REMOTE_PANEL


REMOTE_PANEL = r'''

HOLDOUT = ROOT / "experiments/koblitz-s3-target-conjugate-20261009/MODAL_HOLDOUT_FREEZE.json"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009"
MANIFEST = ROOT / "experiments/koblitz-s3-target-conjugate-20261009/MODAL_UPLOAD_MANIFEST.json"
BINARY = {"baseline": "s3_partial_prepared_parallel_targetspan",
          "expanded": "s3_partial_prepared_fullrank_targetspan"}

def verify_bundle():
    manifest = json.loads(MANIFEST.read_text())
    for entry in manifest["files"]:
        path = ROOT / entry["path"]
        if not path.is_file() or path.stat().st_size != entry["bytes"] or sha(path) != entry["sha256"]:
            raise RuntimeError("uploaded file mismatch: " + entry["path"])
    if sha(HOLDOUT) != manifest["holdout_freeze_sha256"]:
        raise RuntimeError("holdout freeze changed")
    plan = json.loads(HOLDOUT.read_text())
    parent = json.loads((PARENT / "FREEZE.json").read_text())
    assert plan["parent_plan_sha256"] == sha(PARENT / "FREEZE.json")
    assert plan["modal_cpu_request"] == 16 and plan["modal_memory_mib_request"] == 32768
    assert plan["repetitions_per_arm_per_target"] == 2 and plan["target_workers"] == 1
    assert plan["relation_workers"] == 14 and len(plan["targets"]) == 24
    assert [row["label"] for row in plan["targets"]] == [f"T{i:03}" for i in range(33, 57)]
    targets = {(row["n"], row["label"]): row for row in parent["targets"]}
    for arm, name in BINARY.items():
        if sha(ROOT / "examples" / (name + ".rs")) != plan["arms"][arm]["source_sha256"]:
            raise RuntimeError("source hash changed: " + arm)
    return plan, targets

def run_solver(arm, target, frozen, selection, membind, tag):
    binary = BIN / BINARY[arm]
    point_file = PARENT / target["point_file_relative"]
    point = json.loads(point_file.read_text())
    assert sha(point_file) == target["point_sha256"] and point == target["public_point"]
    with tempfile.TemporaryDirectory(prefix="s3-fullrank-") as directory:
        report_path = Path(directory) / "report.json"
        stdout_path = Path(directory) / "stdout"
        stderr_path = Path(directory) / "stderr"
        argv = ["taskset", "-c", ",".join(map(str, selection["selected"]))]
        if membind:
            argv += ["numactl", f"--membind={selection['node']}"]
        argv += [str(binary), "53", "0", str(frozen["columns"]), "20260928",
                 str(point_file), str(report_path), "14",
                 str(frozen["requested_pre_target_rank"]), "1"]
        before = snapshot(selection["selected"])
        started = time.time_ns()
        masks, max_threads, max_vm_rss_kb, error = set(), 0, 0, None
        with stdout_path.open("w") as out, stderr_path.open("w") as err:
            process = subprocess.Popen(argv, stdout=out, stderr=err)
            while process.poll() is None:
                count, current = thread_masks(process.pid)
                max_threads = max(max_threads, count)
                masks.update(current)
                status = read(f"/proc/{process.pid}/status", 20000) or ""
                rss = re.search(r"^VmHWM:\s+(\d+) kB$", status, re.M)
                if rss:
                    max_vm_rss_kb = max(max_vm_rss_kb, int(rss.group(1)))
                if time.time_ns() - started > 180_000_000_000:
                    process.kill()
                    error = "timeout_180s"
                    break
                time.sleep(.02)
            exit_code = process.wait()
        finished = time.time_ns()
        after = snapshot(selection["selected"])
        raw = report_path.read_text() if report_path.exists() else None
        stdout = stdout_path.read_text()
        try:
            report = json.loads(raw) if raw else None
        except json.JSONDecodeError:
            report, error = None, "invalid_json"
        timing = report.get("timing_ms", {}) if report else {}
        phase_ok = bool(report and all(isinstance(timing.get(key), (int, float))
                                       and timing[key] >= 0 for key in PHASES)
                        and abs(sum(timing[key] for key in PHASES) -
                                timing["target_online_after_reusable_setup"]) < .02)
        verified = bool(exit_code == 0 and report and phase_ok and stdout == raw
                        and report.get("group_verified") is True
                        and report.get("target_count") == 1
                        and report.get("target") == point
                        and report.get("recovered_scalar") == target["fixture_scalar"]
                        and report.get("factor_base_points") == frozen["factor_base_points"]
                        and report.get("factor_base_digest") == frozen["factor_base_digest"]
                        and report.get("orbit_columns") == frozen["columns"]
                        and report.get("root_table_slots") == frozen["root_table_slots"]
                        and report.get("prepared_rank") >= frozen["requested_pre_target_rank"]
                        and report.get("target_query_workers") == 1
                        and (arm == "baseline" or (
                            report.get("prepared_rank") == 1200
                            and report.get("target_full_rank_direct_recovery") is True
                            and report.get("target_span_checks") == 0
                            and report.get("target_conjugate_count") == 53
                            and report.get("target_query_batch_width") == 8)))
        return {"kind": "modal_s3_fullrank_run_v1", "tag": tag, "arm": arm,
                "label": target["label"], "candidate_id": frozen["candidate_id"],
                "workload_id": target["workload_id"], "point_sha256": target["point_sha256"],
                "source_sha256": frozen["source_sha256"], "binary_sha256": sha(binary),
                "argv": argv, "started_unix_ns": started, "finished_unix_ns": finished,
                "outer_wall_ns": finished - started, "exit_code": exit_code,
                "error": error, "verified_native_and_fixture": verified,
                "phase_ok": phase_ok, "max_threads_seen": max_threads,
                "max_vm_rss_kb_seen": max_vm_rss_kb,
                "thread_masks_seen": sorted(masks), "before": before, "after": after,
                "stdout_sha256": sha(stdout_path), "stdout_matches_raw": stdout == raw,
                "stderr": stderr_path.read_text()[:20000],
                "raw_sha256": hashlib.sha256(raw.encode()).hexdigest() if raw else None,
                "report": report}

def perf_probe(arm, target, frozen, selection, membind):
    binary = BIN / BINARY[arm]
    point_file = PARENT / target["point_file_relative"]
    with tempfile.TemporaryDirectory(prefix="s3-fullrank-perf-") as directory:
        report_path = Path(directory) / "report.json"
        argv = ["taskset", "-c", ",".join(map(str, selection["selected"]))]
        if membind:
            argv += ["numactl", f"--membind={selection['node']}"]
        argv += ["perf", "stat", "-x", ",", "-e",
                 "cycles,instructions,task-clock,context-switches,cpu-migrations,page-faults,branches,cache-misses",
                 "--", str(binary), "53", "0", str(frozen["columns"]), "20260928",
                 str(point_file), str(report_path), "14",
                 str(frozen["requested_pre_target_rank"]), "1"]
        result = command(argv, timeout=180)
        result["scope"] = "whole process including reusable preparation and one target"
        result["report_sha256"] = sha(report_path) if report_path.exists() else None
        return result

def main():
    plan, targets = verify_bundle()
    selection = select_cpus()
    os.sched_setaffinity(0, {selection["sampler_cpu"]})
    bind_probe = command(["numactl", f"--membind={selection['node']}", "true"])
    membind = bind_probe.get("exit_code") == 0
    machine = host(selection)
    warmup = targets[(53, "T025")]
    runs = [run_solver(arm, warmup, plan["arms"][arm], selection, membind,
                       "warmup-" + arm) for arm in ("baseline", "expanded")]
    blocks = []
    for index, target in enumerate(plan["targets"]):
        arms = ("baseline", "expanded", "expanded", "baseline") if index % 2 else (
            "expanded", "baseline", "baseline", "expanded")
        block = []
        for slot, arm in enumerate(arms, 1):
            row = run_solver(arm, target, plan["arms"][arm], selection, membind,
                             f"{target['label']}-{slot}")
            row["slot"] = slot
            row["repetition"] = 1 + sum(old["arm"] == arm for old in block)
            row["run_id"] = (f"{row['candidate_id']}W{row['workload_id']}"
                             f"R{row['repetition']}")
            runs.append(row)
            block.append(row)
        good = all(row["verified_native_and_fixture"] for row in block)
        if good:
            by_arm = {arm: [row["report"]["timing_ms"]["target_online_after_reusable_setup"]
                            for row in block if row["arm"] == arm]
                      for arm in BINARY}
            ratio = math.sqrt(math.prod(by_arm["baseline"]) / math.prod(by_arm["expanded"]))
        else:
            ratio = None
        blocks.append({"label": target["label"], "order": arms, "all_verified": good,
                       "exploratory_online_ratio": ratio})
        print(f"{target['label']} verified={good} ratio={ratio}", file=sys.stderr, flush=True)
    perf = [{"arm": arm, "result": perf_probe(arm, warmup, plan["arms"][arm], selection, membind)}
            for arm in BINARY]
    doc = {"kind": "modal_s3_fullrank_panel_v1", "holdout_freeze_sha256": sha(HOLDOUT),
           "upload_manifest_sha256": sha(MANIFEST), "host": machine,
           "numa_bind_probe": bind_probe, "numa_bind_used": membind,
           "strict_host_isolation_receipt": None, "controlled_speedup": None,
           "runs": runs, "blocks": blocks, "perf_whole_process": perf,
           "finished_unix_ns": time.time_ns()}
    raw = gzip.compress(json.dumps(doc, separators=(",", ":")).encode())
    Path("/results/panel.json.gz").write_bytes(raw)
    print(json.dumps({"bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest(),
                      "runs": len(runs), "blocks": len(blocks),
                      "all_verified": all(row["verified_native_and_fixture"] for row in runs)}),
          flush=True)

main()
'''


def verify_upload_manifest() -> dict:
    manifest = json.loads(MANIFEST.read_text())
    assert manifest["kind"] == "s3_fullrank_modal_upload_manifest_v1"
    assert manifest["holdout_freeze_sha256"] == sha((HERE / "MODAL_HOLDOUT_FREEZE.json").read_bytes())
    for entry in manifest["files"]:
        path = ROOT / entry["path"]
        assert path.is_file() and path.stat().st_size == entry["bytes"]
        assert sha(path.read_bytes()) == entry["sha256"]
    return manifest


def image(manifest: dict):
    import modal
    item = modal.Image.from_registry("rust:1.98-bookworm", add_python="3.12")
    item = item.apt_install("linux-perf", "numactl", "util-linux", "procps")
    item = item.add_local_dir(str(ROOT / "src"), remote_path=str(REMOTE / "src"), copy=True)
    for entry in manifest["files"]:
        name = entry["path"]
        if name.startswith("src/"):
            continue
        item = item.add_local_file(str(ROOT / name), str(REMOTE / name), copy=True)
    item = item.add_local_file(str(MANIFEST), str(REMOTE / MANIFEST.relative_to(ROOT)), copy=True)
    return item.run_commands(
        f"cd {REMOTE} && cargo build --release --locked "
        "--example s3_partial_prepared_parallel_targetspan "
        "--example s3_partial_prepared_fullrank_targetspan")


def launch(output: Path, region: str | None, manifest: dict, code: str) -> None:
    import modal
    if output.exists():
        raise FileExistsError(output)
    app = modal.App("s3-fullrank-c1200-vm-panel")
    volume_name = f"s3-fullrank-c1200-{uuid.uuid4().hex[:20]}"
    volume = modal.Volume.from_name(volume_name, create_if_missing=True)
    print(json.dumps({"volume_name": volume_name, "remote_runner_sha256": sha(code.encode())}), flush=True)
    with app.run():
        sandbox = modal.Sandbox.create(
            "python3", "-c", code, app=app, image=image(manifest),
            cpu=16, memory=32768, timeout=3600, region=region,
            volumes={"/results": volume},
            experimental_options={"vm_runtime": True})
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
        if len(raw) != expected["bytes"] or sha(raw) != expected["sha256"]:
            raise RuntimeError("Modal result transfer failed size or SHA-256 validation")
        doc = json.loads(gzip.decompress(raw))
        if len(doc["runs"]) != 98 or len(doc["blocks"]) != 24:
            raise RuntimeError("incomplete raw run ledger")
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_bytes(raw)
        print(json.dumps({"output": str(output), "sha256": sha(raw),
                          "runs": len(doc["runs"]), "blocks": len(doc["blocks"]),
                          "all_verified": expected["all_verified"]}, sort_keys=True), flush=True)
        volume.remove_file("panel.json.gz")
    modal.Volume.delete(volume_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--region", default=None)
    parser.add_argument("--launch", action="store_true")
    args = parser.parse_args()
    files = verify_upload_manifest()
    code = remote_code()
    print(json.dumps({"mode": "launch" if args.launch else "local_preview",
                      "upload_manifest_sha256": sha(MANIFEST.read_bytes()),
                      "files": len(files["files"]),
                      "bytes": sum(row["bytes"] for row in files["files"]),
                      "holdout_targets": files["target_count"],
                      "remote_runner_sha256": sha(code.encode()),
                      "output": str(args.output)}, sort_keys=True), flush=True)
    if args.launch:
        launch(args.output, args.region, files, code)

"""Prepare a Modal VM ABBA panel for cooperative S3 target queries.

--preview-upload stays local. A cloud invocation transfers the listed source
and public fixtures; it requires explicit upload approval before execution.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import uuid

import modal


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
REMOTE = Path("/root/crypto")
PLAN = HERE / "FREEZE.json"
GROUPLIFT = ROOT / "experiments/koblitz-s3-grouplift-20261009/FREEZE.json"
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def upload_files() -> list[Path]:
    directories = [ROOT / name for name in ("src", "docs/ic", "docs/ecbench", "docs/curves")]
    files = [ROOT / name for name in (
        "Cargo.toml", "Cargo.lock",
        "examples/s3_partial_prepared_grouplift_targetspan.rs",
        "examples/s3_partial_prepared_parallel_targetspan.rs",
        "experiments/koblitz-s3-grouplift-20261009/FREEZE.json",
        "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json",
        "experiments/koblitz-s3-parallel-target-20261009/FREEZE.json",
    )]
    for directory in directories:
        files.extend(path for path in directory.rglob("*") if path.is_file())
    plan = json.loads(PLAN.read_text())
    files.extend(ROOT / "experiments/koblitz-s3-crosscurve-20261009/holdout" /
                 f"n{target['n']}" / target["label"] / "public_target.json"
                 for target in plan["targets"])
    assert len(files) == len(set(files)) and all(path.is_file() for path in files)
    return sorted(files)


def preview_upload() -> dict:
    plan = json.loads(PLAN.read_text())
    files = upload_files()
    assert len(plan["targets"]) == 48
    old = json.loads(GROUPLIFT.read_text())
    assert sha(GROUPLIFT) == plan["grouplift_plan_sha256"]
    assert sha(PARENT) == plan["parent_plan_sha256"]
    assert sha(ROOT / "examples/s3_partial_prepared_grouplift_targetspan.rs") == old["source_sha256"]
    assert sha(ROOT / "examples/s3_partial_prepared_parallel_targetspan.rs") == plan["source_sha256"]
    for target in plan["targets"]:
        point_file = (ROOT / "experiments/koblitz-s3-crosscurve-20261009/holdout" /
                      f"n{target['n']}" / target["label"] / "public_target.json")
        assert point_file in files and sha(point_file) == target["point_sha256"]
    return {"frozen_plan_sha256": sha(PLAN), "file_count": len(files),
            "total_bytes": sum(path.stat().st_size for path in files),
            "public_points_per_curve": 24,
            "source_sha256": {"grouplift": old["source_sha256"],
                              "parallel": plan["source_sha256"]},
            "included_roots": ["Cargo.toml", "Cargo.lock", "src/", "docs/ic/",
                               "docs/ecbench/", "docs/curves/", "two S3 examples",
                               "three frozen plans and 48 public target points"]}


def image() -> modal.Image:
    item = modal.Image.from_registry("rust:1.98-bookworm", add_python="3.12")
    item = item.apt_install("linux-perf", "numactl", "util-linux", "procps")
    directory_names = ("src", "docs/ic", "docs/ecbench", "docs/curves")
    for name in directory_names:
        item = item.add_local_dir(str(ROOT / name), remote_path=str(REMOTE / name), copy=True)
    directory_roots = tuple(ROOT / name for name in directory_names)
    for path in upload_files():
        if any(root in path.parents for root in directory_roots):
            continue
        item = item.add_local_file(str(path), str(REMOTE / path.relative_to(ROOT)), copy=True)
    return item.run_commands(
        f"cd {REMOTE} && cargo build --release --locked "
        "--example s3_partial_prepared_grouplift_targetspan "
        "--example s3_partial_prepared_parallel_targetspan"
    )


def remote_code() -> str:
    prior_path = ROOT / "experiments/koblitz-s3-prepared-20261009/modal_prepared_panel.py"
    spec = importlib.util.spec_from_file_location("prior_modal_s3", prior_path)
    prior = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prior)
    common, marker, _ = prior.REMOTE_CODE.partition("def run_solver(")
    assert marker and "def select_cpus():" in common
    return common + REMOTE_PANEL_CODE


REMOTE_PANEL_CODE = r'''
PLAN = ROOT / "experiments/koblitz-s3-parallel-target-20261009/FREEZE.json"
GROUPLIFT_PLAN = ROOT / "experiments/koblitz-s3-grouplift-20261009/FREEZE.json"
PARENT_PLAN = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
POINTS = ROOT / "experiments/koblitz-s3-crosscurve-20261009/holdout"

def target_order(n, label):
    index = int(label[1:])
    return (("grouplift", "parallel", "parallel", "grouplift")
            if (index + (n == 53)) % 2 else
            ("parallel", "grouplift", "grouplift", "parallel"))

def run_solver(arm, target, repetition, selection, membind, tag, plan):
    n = target["n"]
    binary = BIN / ("s3_partial_prepared_grouplift_targetspan" if arm == "grouplift"
                    else "s3_partial_prepared_parallel_targetspan")
    point_file = POINTS / f"n{n}" / target["label"] / "public_target.json"
    assert sha(point_file) == target["point_sha256"]
    point = json.loads(point_file.read_text())
    candidate_id = plan["curves"][str(n)][
        "grouplift_candidate_id" if arm == "grouplift" else "parallel_candidate_id"]
    run_id = f"{candidate_id}W{target['workload_id']}R{repetition}"
    with tempfile.TemporaryDirectory(prefix="s3-crosscurve-") as directory:
        output = Path(directory) / "raw.jsonl"
        stdout = Path(directory) / "stdout"
        stderr = Path(directory) / "stderr"
        argv = ["taskset", "-c", ",".join(map(str, selection["selected"]))]
        if membind:
            argv += ["numactl", f"--membind={selection['node']}"]
        argv += [str(binary), str(n), "0", "244", str(plan["relation_seed"]),
                 str(point_file), str(output), str(WORKERS)]
        argv.append("243")
        if arm == "parallel":
            argv.append(str(plan["target_workers_by_curve"][str(n)]))
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
            report, error = None, "invalid_json"
        times = report.get("timing_ms", {}) if report else {}
        phase_ok = bool(report and all(isinstance(times.get(key), (int, float))
                                       and times[key] >= 0 for key in PHASES)
                        and abs(sum(times[key] for key in PHASES) -
                                times["target_online_after_reusable_setup"]) < 0.02)
        verified = bool(exit_code == 0 and report and report.get("group_verified") is True
                        and report.get("recovered_scalar") == target["fixture_scalar"]
                        and report.get("target") == point and report.get("target_count") == 1
                        and report.get("n") == n and report.get("a") == 0
                        and report.get("relation_seed") == plan["relation_seed"]
                        and report.get("factor_base_points") ==
                            plan["factor_base_points_by_curve"][str(n)]
                        and report.get("relation_collection_workers") == WORKERS
                        and report.get("prepared_rank_reached") is True
                        and report.get("prepare_rank_requested") == 243
                        and (arm == "grouplift" or report.get("target_query_workers") ==
                             plan["target_workers_by_curve"][str(n)])
                        and phase_ok and stdout_text == raw)
        return {"tag": tag, "n": n, "family": target["family"],
                "label": target["label"], "arm": arm, "run_id": run_id,
                "candidate_id": candidate_id, "workload_id": target["workload_id"],
                "point": point, "target_sha256": sha(point_file),
                "binary_sha256": sha(binary), "argv": argv,
                "started_unix_ns": started, "finished_unix_ns": finished,
                "outer_wall_ns": finished - started, "exit_code": exit_code,
                "error": error, "verified_native_and_fixture": verified,
                "phase_ok": phase_ok, "max_threads_seen": max_threads,
                "thread_masks_seen": sorted(masks), "before": before, "after": after,
                "stdout_sha256": sha(stdout), "stdout_matches_raw": stdout_text == raw,
                "stderr": stderr.read_text()[:20000],
                "raw_sha256": hashlib.sha256(raw.encode()).hexdigest() if raw else None,
                "raw_text": raw}

def perf_probe(arm, target, selection, membind, plan):
    n = target["n"]
    binary = BIN / ("s3_partial_prepared_grouplift_targetspan" if arm == "grouplift"
                    else "s3_partial_prepared_parallel_targetspan")
    point_file = POINTS / f"n{n}" / target["label"] / "public_target.json"
    with tempfile.TemporaryDirectory(prefix="s3-crosscurve-perf-") as directory:
        output = Path(directory) / "raw.jsonl"
        argv = ["taskset", "-c", ",".join(map(str, selection["selected"]))]
        if membind:
            argv += ["numactl", f"--membind={selection['node']}"]
        argv += ["perf", "stat", "-x", ",", "-e",
                 "cycles,instructions,task-clock,context-switches,cpu-migrations,page-faults,branches,cache-misses",
                 "--", str(binary), str(n), "0", "244", str(plan["relation_seed"]),
                 str(point_file), str(output), str(WORKERS)]
        argv.append("243")
        if arm == "parallel":
            argv.append(str(plan["target_workers_by_curve"][str(n)]))
        result = command(argv, timeout=180)
        result["scope"] = "whole process including reusable preparation and one target"
        result["raw_sha256"] = sha(output) if output.exists() else None
        return result

def main(first, count):
    plan = json.loads(PLAN.read_text())
    old = json.loads(GROUPLIFT_PLAN.read_text())
    parent = json.loads(PARENT_PLAN.read_text())
    assert plan["target_count_per_curve"] == 24 and plan["workers"] == WORKERS
    assert sha(GROUPLIFT_PLAN) == plan["grouplift_plan_sha256"]
    assert sha(PARENT_PLAN) == plan["parent_plan_sha256"]
    plan["factor_base_points_by_curve"] = {
        str(n): parent["curves"][str(n)]["factor_base_points"] for n in (41, 53)}
    assert sha(ROOT / "examples/s3_partial_prepared_grouplift_targetspan.rs") == old["source_sha256"]
    assert sha(ROOT / "examples/s3_partial_prepared_parallel_targetspan.rs") == plan["source_sha256"]
    parent_targets = {(row["n"], row["label"]): row for row in parent["targets"]}
    by_target = {}
    for row in plan["targets"]:
        key = (row["n"], row["label"])
        selected = parent_targets[key]
        assert selected["workload_id"] == row["workload_id"]
        assert selected["point_sha256"] == row["point_sha256"]
        by_target[key] = selected
    assert len(by_target) == 48 and 1 <= first and 1 <= count <= 24
    assert first + count <= 25
    selection = select_cpus()
    if selection["affinity_probe"].get("exit_code") != 0:
        raise RuntimeError("selected CPU affinity was rejected")
    os.sched_setaffinity(0, {selection["sampler_cpu"]})
    bind_probe = command(["numactl", f"--membind={selection['node']}", "true"])
    membind = bind_probe.get("exit_code") == 0
    runs, blocks = [], []
    for index in range(first, first + count):
        for n in ((41, 53) if index % 2 else (53, 41)):
            target = by_target[(n, f"T{index:03}")]
            order = target_order(n, target["label"])
            repetitions = {"grouplift": 0, "parallel": 0}
            block = []
            for position, arm in enumerate(order):
                repetitions[arm] += 1
                row = run_solver(arm, target, repetitions[arm], selection, membind,
                                 f"n{n}-{target['label']}-{position + 1}", plan)
                runs.append(row)
                block.append(row)
            good = all(row["verified_native_and_fixture"] for row in block)
            times = [(row["arm"], json.loads(row["raw_text"])["timing_ms"][
                "target_online_after_reusable_setup"]) for row in block] if good else []
            by_arm = {arm: [value for name, value in times if name == arm]
                      for arm in ("grouplift", "parallel")} if good else {}
            ratio = (math.sqrt(math.prod(by_arm["grouplift"]) /
                               math.prod(by_arm["parallel"])) if good else None)
            aa = ({arm: abs(math.log(values[0] / values[1]))
                   for arm, values in by_arm.items()} if good else None)
            blocks.append({"n": n, "label": target["label"],
                           "workload_id": target["workload_id"],
                           "point": json.loads((POINTS / f"n{n}" / target["label"] /
                                                "public_target.json").read_text()),
                           "order": order,
                           "run_tags": [row["tag"] for row in block],
                           "all_verified": good, "exploratory_online_ratio": ratio,
                           "aa_abs_log_ratio_by_arm": aa})
            print(f"n{n} {target['label']} verified={good} ratio={ratio}", file=sys.stderr, flush=True)
    summary = {}
    for n in (41, 53):
        selected = [block for block in blocks if block["n"] == n]
        ratios = [block["exploratory_online_ratio"] for block in selected
                  if block["exploratory_online_ratio"] is not None]
        summary[str(n)] = {"targets": len(selected),
                           "verified": len(ratios),
                           "minimum_online_ratio": min(ratios) if ratios else None,
                           "geometric_mean_online_ratio":
                               math.exp(statistics.mean(map(math.log, ratios))) if ratios else None,
                           "at_least_2_count": sum(ratio >= 2 for ratio in ratios)}
    first_targets = {n: by_target[(n, f"T{first:03}")] for n in (41, 53)}
    perf = [{"n": n, "arm": arm,
             "result": perf_probe(arm, first_targets[n], selection, membind, plan)}
            for n in (41, 53) for arm in ("grouplift", "parallel")]
    doc = {"kind": "modal_vm_s3_parallel_target_panel_v1",
           "first_target_index": first, "target_count_per_curve": count,
           "frozen_plan_sha256": sha(PLAN),
           "source_sha256": {"grouplift": sha(ROOT / "examples/s3_partial_prepared_grouplift_targetspan.rs"),
                             "parallel": sha(ROOT / "examples/s3_partial_prepared_parallel_targetspan.rs")},
           "binary_sha256": {"grouplift": sha(BIN / "s3_partial_prepared_grouplift_targetspan"),
                             "parallel": sha(BIN / "s3_partial_prepared_parallel_targetspan")},
           "host": host(selection), "numa_bind_probe": bind_probe,
           "numa_bind_used": membind,
           "strict_host_isolation_receipt": None, "controlled_speedup": None,
           "runs": runs, "blocks": blocks, "summary": summary,
           "perf_whole_process": perf, "finished_unix_ns": time.time_ns()}
    compressed = gzip.compress(json.dumps(doc, separators=(",", ":")).encode())
    Path("/results/panel.json.gz").write_bytes(compressed)
    print(json.dumps({"bytes": len(compressed),
                      "sha256": hashlib.sha256(compressed).hexdigest(),
                      "runs": len(runs), "blocks": len(blocks),
                      "summary": summary}), flush=True)

main(int(sys.argv[1]), int(sys.argv[2]))
'''


def run(first: int, count: int, output: Path, region: str | None) -> None:
    if not (1 <= first <= 24 and 1 <= count <= 24 and first + count <= 25):
        raise ValueError("select 1..24 consecutive indices within the frozen 1..24 range")
    if output.exists():
        raise FileExistsError(output)
    preview = preview_upload()
    application = modal.App("s3-parallel-target-matched-curve-vm-panel")
    volume_name = f"s3-parallel-{uuid.uuid4().hex[:20]}"
    volume = modal.Volume.from_name(volume_name, create_if_missing=True)
    code = remote_code()
    print(json.dumps({"upload": preview, "volume_name": volume_name,
                      "remote_runner_sha256": hashlib.sha256(code.encode()).hexdigest()}),
          flush=True)
    with application.run():
        sandbox = modal.Sandbox.create(
            "python3", "-c", code, str(first), str(count), app=application,
            image=image(), cpu=16, memory=16384, timeout=1800, region=region,
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
        if len(doc["runs"]) != count * 2 * 4 or len(doc["blocks"]) != count * 2:
            raise RuntimeError("incomplete raw run ledger")
        if doc["frozen_plan_sha256"] != preview["frozen_plan_sha256"]:
            raise RuntimeError("frozen plan hash mismatch")
        if doc["source_sha256"]["grouplift"] != preview["source_sha256"]["grouplift"]:
            raise RuntimeError("group-lift source hash mismatch")
        if doc["source_sha256"]["parallel"] != preview["source_sha256"]["parallel"]:
            raise RuntimeError("parallel source hash mismatch")
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_bytes(raw)
        print(json.dumps({"output": str(output), "sha256": digest,
                          "runs": len(doc["runs"]), "blocks": len(doc["blocks"]),
                          "summary": doc["summary"]}, sort_keys=True), flush=True)
        volume.remove_file("panel.json.gz")
    modal.Volume.delete(volume_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first", type=int, default=1)
    parser.add_argument("--count", type=int, choices=range(1, 25), default=1)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--region", default=None)
    parser.add_argument("--preview-upload", action="store_true")
    args = parser.parse_args()
    if args.preview_upload:
        print(json.dumps(preview_upload(), sort_keys=True))
    else:
        if args.output is None:
            parser.error("--output is required for a Modal run")
        run(args.first, args.count, args.output, args.region)

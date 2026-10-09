"""Prepare a Modal VM panel for frozen N53 rank-243 S3 one-target pairs.

``--preview-upload`` is local and read-only. The actual Modal invocation uses
``--count 1`` for a fresh-target/full-rank smoke or ``--count 24`` for the
complete frozen panel. Guest counters are diagnostics; controlled CPU timing
still requires an auditable host-wide isolation receipt.
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
SOURCE_COMMIT = "833ef86191947ba8e4def56a78ee6c51a4ea0ff2"
PLAN = HERE / "rank243/FREEZE.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def upload_files() -> list[Path]:
    directories = [ROOT / name for name in ("src", "docs/ic", "docs/ecbench", "docs/curves")]
    files = [ROOT / name for name in (
        "Cargo.toml", "Cargo.lock",
        "examples/s3_partial_prepared_targetspan.rs",
        "experiments/koblitz-s3-pair-query-20261007-v4/baseline.rs",
        "experiments/koblitz-s3-pair-query-20261007-v4/candidate.rs",
        "experiments/koblitz-s3-partial-prepared-20261009/rank243/FREEZE.json",
        "experiments/koblitz-s3-partial-prepared-20261009/rank243/holdout/n53/fixtures.json",
    )]
    for directory in directories:
        files.extend(path for path in directory.rglob("*") if path.is_file())
    fresh_dir = HERE / "rank243/holdout/n53"
    files.extend(sorted(fresh_dir.glob("T*/public_target.json")))
    files.extend(HERE / f"holdout/n53/T{index:03}/public_target.json"
                 for index in (67, 85, 96))
    assert len(files) == len(set(files))
    assert all(path.is_file() for path in files)
    return sorted(files)


def preview_upload() -> dict:
    files = upload_files()
    plan = json.loads(PLAN.read_text())
    assert sha(ROOT / "examples/s3_partial_prepared_targetspan.rs") == plan["partial_source_sha256"]
    assert sha(ROOT / "experiments/koblitz-s3-pair-query-20261007-v4/candidate.rs") == plan["reference_source_sha256"]
    assert len(plan["targets"]) == 27
    for target in plan["targets"]:
        point = HERE / target["point_file_relative"]
        assert point in files and sha(point) == target["point_sha256"]
    return {"source_commit": SOURCE_COMMIT, "frozen_plan_sha256": sha(PLAN),
            "file_count": len(files), "total_bytes": sum(path.stat().st_size for path in files),
            "source_sha256": {"reference": plan["reference_source_sha256"],
                              "rank243": plan["partial_source_sha256"]},
            "fresh_public_points": 24, "stress_public_points": 3,
            "included_roots": ["Cargo.toml", "Cargo.lock", "src/", "docs/ic/",
                               "docs/ecbench/", "docs/curves/", "two S3 examples",
                               "frozen plan and 27 public target points"]}


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
        relative = path.relative_to(ROOT)
        item = item.add_local_file(str(path), str(REMOTE / relative), copy=True)
    return item.run_commands(
        f"cd {REMOTE} && cargo build --release --locked "
        "--example s3_pair_candidate_frozen --example s3_partial_prepared_targetspan"
    )


def remote_code() -> str:
    prior_path = ROOT / "experiments/koblitz-s3-prepared-20261009/modal_prepared_panel.py"
    spec = importlib.util.spec_from_file_location("prior_modal_s3", prior_path)
    prior = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prior)
    common, marker, _ = prior.REMOTE_CODE.partition("def run_solver(")
    assert marker and "def thread_masks(pid):" in common
    return common + REMOTE_PANEL_CODE


REMOTE_PANEL_CODE = r'''
PLAN = ROOT / "experiments/koblitz-s3-partial-prepared-20261009/rank243/FREEZE.json"
EXPERIMENT = ROOT / "experiments/koblitz-s3-partial-prepared-20261009"
ARMS = ("reference", "shallow", "deep")

def run_solver(arm, target, selection, membind, tag, plan):
    binary = BIN / ("s3_pair_candidate_frozen" if arm == "reference"
                    else "s3_partial_prepared_targetspan")
    point_file = EXPERIMENT / target["point_file_relative"]
    assert sha(point_file) == target["point_sha256"]
    point = json.loads(point_file.read_text())
    with tempfile.TemporaryDirectory(prefix="s3-partial-") as directory:
        output = Path(directory) / "raw.jsonl"
        stdout = Path(directory) / "stdout"
        stderr = Path(directory) / "stderr"
        argv = ["taskset", "-c", ",".join(map(str, selection["selected"]))]
        if membind:
            argv += ["numactl", f"--membind={selection['node']}"]
        argv += [str(binary), "53", "0", "244", str(plan["relation_seed"]),
                 str(point_file), str(output), str(WORKERS)]
        if arm != "reference":
            argv.append(str(plan["rank_policies"][arm]))
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
        phase_ok = bool(report and all(isinstance(times.get(key), (int, float))
                                       and times[key] >= 0 for key in PHASES)
                        and abs(sum(times[key] for key in PHASES) -
                                times["target_online_after_reusable_setup"]) < 0.02)
        expected_rank = plan["rank_policies"].get(arm)
        verified = bool(exit_code == 0 and report and report.get("group_verified") is True
                        and report.get("recovered_scalar") == target["fixture_scalar"]
                        and report.get("target") == point and report.get("target_count") == 1
                        and report.get("n") == 53 and report.get("a") == 0
                        and report.get("relation_seed") == plan["relation_seed"]
                        and report.get("factor_base_points") == 25864
                        and report.get("relation_collection_workers") == WORKERS
                        and (expected_rank is None or
                             (report.get("prepared_rank_reached") is True
                              and report.get("prepare_rank_requested") == expected_rank))
                        and phase_ok and stdout_text == raw)
        return {"tag": tag, "family": target["family"], "label": target["label"],
                "arm": arm, "point": point,
                "target_sha256": sha(point_file), "binary_sha256": sha(binary),
                "argv": argv, "started_unix_ns": started,
                "finished_unix_ns": finished, "outer_wall_ns": finished - started,
                "exit_code": exit_code, "error": error,
                "verified_native_and_fixture": verified, "phase_ok": phase_ok,
                "max_threads_seen": max_threads, "thread_masks_seen": sorted(masks),
                "before": before, "after": after, "stdout_sha256": sha(stdout),
                "stdout_matches_raw": stdout_text == raw,
                "stderr": stderr.read_text()[:20000],
                "raw_sha256": hashlib.sha256(raw.encode()).hexdigest() if raw else None,
                "raw_text": raw}

def perf_probe(arm, target, selection, membind, plan):
    binary = BIN / ("s3_pair_candidate_frozen" if arm == "reference"
                    else "s3_partial_prepared_targetspan")
    point_file = EXPERIMENT / target["point_file_relative"]
    with tempfile.TemporaryDirectory(prefix="s3-partial-perf-") as directory:
        output = Path(directory) / "raw.jsonl"
        argv = ["taskset", "-c", ",".join(map(str, selection["selected"]))]
        if membind:
            argv += ["numactl", f"--membind={selection['node']}"]
        argv += ["perf", "stat", "-x", ",", "-e",
                 "cycles,instructions,task-clock,context-switches,cpu-migrations,page-faults,branches,cache-misses",
                 "--", str(binary), "53", "0", "244", str(plan["relation_seed"]),
                 str(point_file), str(output), str(WORKERS)]
        if arm != "reference":
            argv.append(str(plan["rank_policies"][arm]))
        result = command(argv, timeout=180)
        result["scope"] = "whole process including reusable preparation and one target"
        result["raw_sha256"] = sha(output) if output.exists() else None
        return result

def main(count):
    plan = json.loads(PLAN.read_text())
    assert plan["rank_policies"] == {"shallow": 242, "deep": 243}
    assert plan["fresh_target_count"] == 24 and plan["stress_target_count"] == 3
    targets = ([target for target in plan["targets"]
                if target["family"] == "fresh"][:count] +
               [target for target in plan["targets"]
                if target["family"] == "full_rank_stress_control"])
    selection = select_cpus()
    if selection["affinity_probe"].get("exit_code") != 0:
        raise RuntimeError("selected CPU affinity was rejected")
    if selection["sampler_cpu"] is not None:
        os.sched_setaffinity(0, {selection["sampler_cpu"]})
    bind_probe = command(["numactl", f"--membind={selection['node']}", "true"])
    membind = bind_probe.get("exit_code") == 0
    runs, blocks = [], []
    for index, target in enumerate(targets):
        if index % 2:
            order = ("deep", "reference", "reference", "deep")
        else:
            order = ("reference", "deep", "deep", "reference")
        block = []
        for position, arm in enumerate(order + ("shallow",)):
            row = run_solver(arm, target, selection, membind,
                             f"{target['label']}-{position + 1}", plan)
            runs.append(row)
            block.append(row)
        good = all(row["verified_native_and_fixture"] for row in block)
        by_arm = {arm: [json.loads(row["raw_text"])["timing_ms"]["target_online_after_reusable_setup"]
                       for row in block if row["arm"] == arm]
                  for arm in ARMS} if good else {}
        ratio = (math.sqrt(math.prod(by_arm["reference"]) /
                           math.prod(by_arm["deep"])) if good else None)
        aa = (abs(math.log(by_arm["deep"][0] / by_arm["deep"][1])) if good else None)
        blocks.append({"label": target["label"], "family": target["family"],
                       "point": target["public_point"], "order": order + ("shallow",),
                       "run_tags": [row["tag"] for row in block],
                       "all_verified": good, "exploratory_online_ratio": ratio,
                       "deep_aa_abs_log_ratio": aa})
        print(f"{target['label']} verified={good} ratio={ratio}", file=sys.stderr, flush=True)
    summary = {}
    for family in ("fresh", "full_rank_stress_control"):
        selected = [row for row in blocks if row["family"] == family]
        ratios = [row["exploratory_online_ratio"] for row in selected
                  if row["exploratory_online_ratio"] is not None]
        aa = [row["deep_aa_abs_log_ratio"] for row in selected
              if row["deep_aa_abs_log_ratio"] is not None]
        summary[family] = {"targets": len(selected),
                           "all_verified": all(row["all_verified"] for row in selected),
                           "geometric_mean_ratio": math.exp(statistics.mean(map(math.log, ratios))) if ratios else None,
                           "min_target_ratio": min(ratios) if ratios else None,
                           "pairs_at_least_2": sum(ratio >= 2 for ratio in ratios),
                           "sample_sd_log_ratio": statistics.stdev(map(math.log, ratios)) if len(ratios) > 1 else None,
                           "deep_aa_p95_abs_log_ratio": sorted(aa)[math.ceil(.95 * len(aa)) - 1] if aa else None}
    perf = [{"arm": arm,
             "result": perf_probe(arm, targets[0], selection, membind, plan)}
            for arm in ARMS]
    doc = {"kind": "modal_vm_s3_rank243_single_target_panel_v1",
           "source_commit": "833ef86191947ba8e4def56a78ee6c51a4ea0ff2",
           "frozen_plan_sha256": sha(PLAN),
           "source_sha256": {"reference": sha(ROOT / "experiments/koblitz-s3-pair-query-20261007-v4/candidate.rs"),
                             "partial": sha(ROOT / "examples/s3_partial_prepared_targetspan.rs")},
           "binary_sha256": {"reference": sha(BIN / "s3_pair_candidate_frozen"),
                             "partial": sha(BIN / "s3_partial_prepared_targetspan")},
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

main(int(sys.argv[1]))
'''


def run(count: int, output: Path, region: str | None) -> None:
    if output.exists():
        raise FileExistsError(output)
    preview = preview_upload()
    application = modal.App("s3-rank243-targetspan-vm-panel")
    volume_name = f"s3-rank243-panel-{uuid.uuid4().hex[:20]}"
    volume = modal.Volume.from_name(volume_name, create_if_missing=True)
    code = remote_code()
    print(json.dumps({"upload": preview, "volume_name": volume_name,
                      "remote_runner_sha256": hashlib.sha256(code.encode()).hexdigest()}),
          flush=True)
    with application.run():
        sandbox = modal.Sandbox.create(
            "python3", "-c", code, str(count), app=application,
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
        expected_blocks = count + 3
        if len(doc["runs"]) != expected_blocks * 5 or len(doc["blocks"]) != expected_blocks:
            raise RuntimeError("incomplete raw run ledger")
        if doc["frozen_plan_sha256"] != preview["frozen_plan_sha256"]:
            raise RuntimeError("Modal frozen-plan hash mismatch")
        if doc["source_sha256"]["reference"] != preview["source_sha256"]["reference"]:
            raise RuntimeError("Modal reference source hash mismatch")
        if doc["source_sha256"]["partial"] != preview["source_sha256"]["rank243"]:
            raise RuntimeError("Modal rank-243 source hash mismatch")
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_bytes(raw)
        print(json.dumps({"output": str(output), "sha256": digest,
                          "runs": len(doc["runs"]), "blocks": len(doc["blocks"]),
                          "summary": doc["summary"]}, sort_keys=True), flush=True)
        volume.remove_file("panel.json.gz")
    modal.Volume.delete(volume_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--count", type=int, choices=(1, 24), default=1)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--region", default=None)
    parser.add_argument("--preview-upload", action="store_true")
    args = parser.parse_args()
    if args.preview_upload:
        print(json.dumps(preview_upload(), sort_keys=True))
    else:
        if args.output is None:
            parser.error("--output is required for a Modal run")
        run(args.count, args.output, args.region)

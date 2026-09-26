#!/usr/bin/env python3
"""Seal six clean-source paired runs into a deterministic, replayable archive."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import shutil
import statistics
import tarfile
from pathlib import Path

RUNS = tuple((n, i) for n in (37, 41) for i in range(3))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    runs = args.runs.resolve()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    stage = out / "contents"
    stage.mkdir()
    summary = {"classification": "clean main-derived source; paired full-rank control", "runs": [], "aggregate": {}}
    for n, index in RUNS:
        name = f"n{n}_seed{index}"
        source = runs / name
        target = stage / name
        target.mkdir()
        for filename in ("ic.jsonl", "rho.jsonl", "ic.stderr.txt", "rho.stderr.txt", "manifest.json", "replay.json", "replay.stderr.txt"):
            shutil.copyfile(source / filename, target / filename)
        manifest = json.loads((target / "manifest.json").read_text())
        replay = json.loads((target / "replay.json").read_text())
        assert manifest["clean_checkout"] and manifest["includes_pinned_main_ref"], name
        assert manifest["replay_verdict"] == replay["verdict"] == "PASS", name
        assert len(manifest["arms"]) == 2
        arms = {item["arm"]: item for item in manifest["arms"]}
        assert set(arms) == {"ic", "rho"}
        assert all(item["returncode"] == 0 and not item["timed_out"] for item in arms.values())
        assert replay["ic_stdout_sha256"] == sha(target / "ic.jsonl")
        assert replay["rho_stdout_sha256"] == sha(target / "rho.jsonl")
        ic, rho = arms["ic"], arms["rho"]
        summary["runs"].append({
            "n": n,
            "seed_index": index,
            "seed": replay["public_hash_seed"],
            "q": replay["q"],
            "recovered_scalar": replay["recovered_scalar"],
            "relations_replayed": replay["relations_replayed"],
            "base_labels_replayed": replay["base_points_orbit_labels_replayed"],
            "terminal_rank": replay["terminal_rank"],
            "rho_automorphism_size": replay["rho_automorphism_size"],
            "rho_steps_over_ideal": replay["rho_steps_over_ideal"],
            "ic_wall_ms": ic["wall_ms"],
            "rho_wall_ms": rho["wall_ms"],
            "rho_over_ic_wall": rho["wall_ms"] / ic["wall_ms"],
            "ic_cpu_s": ic["cpu_s"],
            "rho_cpu_s": rho["cpu_s"],
            "rho_over_ic_cpu": rho["cpu_s"] / ic["cpu_s"],
            "ic_peak_rss_bytes": ic["peak_rss_bytes"],
            "rho_peak_rss_bytes": rho["peak_rss_bytes"],
            "ic_stdout_sha256": ic["stdout_sha256"],
            "rho_stdout_sha256": rho["stdout_sha256"],
            "ic_binary_sha256": ic["binary_sha256"],
            "rho_binary_sha256": rho["binary_sha256"],
            "checkout_head": manifest["checkout_head"],
            "source_sha256": manifest["source_sha256"],
        })
    for n in (37, 41):
        samples = [item for item in summary["runs"] if item["n"] == n]
        summary["aggregate"][str(n)] = {
            "rho_over_ic_wall_median": statistics.median(item["rho_over_ic_wall"] for item in samples),
            "rho_over_ic_cpu_median": statistics.median(item["rho_over_ic_cpu"] for item in samples),
            "relations_replayed_total": sum(item["relations_replayed"] for item in samples),
            "max_ic_rss_bytes": max(item["ic_peak_rss_bytes"] for item in samples),
            "all_replay_pass": True,
        }
    summary_path = out / "clean_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    shutil.copyfile(summary_path, stage / "clean_summary.json")
    paths = sorted(p for p in stage.rglob("*") if p.is_file())
    (stage / "SHA256SUMS").write_text(
        "".join(f"{sha(path)}  {path.relative_to(stage)}\n" for path in paths)
    )
    archive = out / "clean_evidence.tar.gz"
    with archive.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=9, mtime=0) as zipped:
            with tarfile.open(fileobj=zipped, mode="w") as tar:
                for path in sorted(p for p in stage.rglob("*") if p.is_file()):
                    data = path.read_bytes()
                    info = tarfile.TarInfo(str(path.relative_to(stage)))
                    info.size = len(data)
                    info.mode = 0o644
                    info.mtime = 0
                    info.uid = info.gid = 0
                    info.uname = info.gname = ""
                    tar.addfile(info, io.BytesIO(data))
    (out / "clean_archive_manifest.json").write_text(json.dumps({
        "archive_sha256": sha(archive),
        "archive_bytes": archive.stat().st_size,
        "reassembly": "tar -xzf clean_evidence.tar.gz && shasum -a 256 -c SHA256SUMS",
        "summary_sha256": sha(summary_path),
        "classification": summary["classification"],
    }, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary["aggregate"], sort_keys=True))


if __name__ == "__main__":
    main()

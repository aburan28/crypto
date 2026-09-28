#!/usr/bin/env python3
"""Count unused affine rank on discovery traces after the timing run finishes."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
PRIMARY = HERE.parent / "run_01"

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert (PRIMARY / "manifest.json").exists(), "Wait for the uninstrumented run."
    expected = {}
    keys = ["outcome", "model", "trace", "nodes", "decisions", "forced", "kernel_calls", "specialized_terms", "source_rows", "source_columns", "max_depth"]
    for n in [16,24]:
        for line in (PRIMARY / f"raw-n{n}.jsonl").read_text().splitlines():
            row = json.loads(line)
            if row["type"] == "sample" and row["split"] == "discovery" and row["variant"] == "tail_wide" and row["rep"] == 0:
                expected[row["cell"]] = {k:row[k] for k in keys}
    assert len(expected) == 12
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    for path in HERE.iterdir():
        if path.is_file(): shutil.copyfile(path,out/path.name)
    metadata = {
        "scope": "Discrete diagnostic counts on twelve discovery fixtures. Instrumented timings are not performance evidence; no holdouts are sampled.",
        "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "primary_manifest_sha256": sha(PRIMARY/"manifest.json"),
        "source_hashes": {p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()},
        "complete": False,
    }
    with tempfile.TemporaryDirectory(prefix="boolean-affine-opportunity-") as temp:
        executable = Path(temp)/"probe"
        command = [shutil.which("rustc"),"--edition","2021","-O",str(out/"worker.rs"),"-o",str(executable)]
        result = subprocess.run(command,capture_output=True,text=True,timeout=60)
        (out/"compile.txt").write_text(result.stdout+result.stderr)
        metadata.update(compiler_command=command,compiler_exit_code=result.returncode)
        (out/"receipt.json").write_text(json.dumps(metadata,indent=2)+"\n")
        assert result.returncode == 0
        metadata["executable_sha256"] = sha(executable)
        result = subprocess.run([str(executable)],capture_output=True,text=True,timeout=120)
        (out/"raw.jsonl").write_text(result.stdout)
        (out/"stderr.txt").write_text(result.stderr)
        metadata.update(command=[str(executable)],exit_code=result.returncode)
        (out/"receipt.json").write_text(json.dumps(metadata,indent=2)+"\n")
        assert result.returncode == 0
    records = [json.loads(s) for s in result.stdout.splitlines()]
    assert len(records)==12
    for row in records:
        cell = f"n{row['n']}-discovery-{row['seed']}-{row['family']}"
        assert {k:row[k] for k in keys} == expected[cell]
        assert row["affine_opportunity_calls"] <= row["decisions"]
        assert row["affine_opportunity_calls"] <= row["affine_opportunity_rank_sum"]
        assert row["affine_opportunity_rank_sum"] <= row["nonunit_affine_rows"]
    summary = {
        "scope": metadata["scope"],
        "trace_model_and_all_logical_counters_match": True,
        "cases": records,
        "total_decisions": sum(r["decisions"] for r in records),
        "opportunity_calls": sum(r["affine_opportunity_calls"] for r in records),
        "rank_sum": sum(r["affine_opportunity_rank_sum"] for r in records),
        "interpretation": "An opportunity is a mixed nonlinear, consistent state with no remaining unit force but at least one nonunit affine row. Rank sums repeat constraints across search states; they are not independent collected relations or a measured speedup.",
        "substitution_solver_cost": None,
    }
    (out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    metadata.update(complete=True,ended_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    (out/"receipt.json").write_text(json.dumps(metadata,indent=2)+"\n")
    (out/"manifest.json").write_text(json.dumps({"files":{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()},"scope":metadata["scope"]},indent=2)+"\n")
    print(json.dumps({k:summary[k] for k in ["total_decisions","opportunity_calls","rank_sum","trace_model_and_all_logical_counters_match"]}))

if __name__ == "__main__":
    main()

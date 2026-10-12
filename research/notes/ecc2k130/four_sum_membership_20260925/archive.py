#!/usr/bin/env python3
"""Archive six frozen producer runs and three independent validations losslessly."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import shutil

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
ARMS = [(37, 3), (41, 8), (41, 12)]
ATTEMPTS = {
    "extractor-n37-R3-host-contended": "valid target/correctness; local wall overlaps #758 input generation 10:19:42–10:19:50 UTC",
    "oracle-n41-R12-full-file-prefix": "valid first128 oracle answers; superseded to ensure byte-identical literal input file",
    "extractor-n41-R12-512-invalid": "invalid 512-target input file for preregistered 128-target R12 arm",
}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def frozen_gzip(data: bytes) -> bytes:
    buffer = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buffer, mtime=0,
                       compresslevel=9) as stream:
        stream.write(data)
    return buffer.getvalue()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs-root", type=Path, required=True)
    parser.add_argument("--validation-root", type=Path, required=True)
    parser.add_argument("--old-runner", type=Path, required=True)
    parser.add_argument("--old-protocol", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=HERE / "evidence")
    args = parser.parse_args()
    evidence = args.out.resolve()
    assert not evidence.exists(), "immutable archive already exists"
    (evidence / "runs").mkdir(parents=True)
    (evidence / "source").mkdir()
    (evidence / "validations").mkdir()
    (evidence / "attempts").mkdir()
    meta = {"schema_version": "1.0", "run_names": [], "source_sha256": {},
            "stdout_sha256": {}, "validation_sha256": {}, "attempts": ATTEMPTS,
            "method_source_sha256": {}, "historical_source_sha256": {},
            "original_source_input_freeze_commit": "12d2e489de77dc633f02df47de2debf865a3f301",
            "budget_amendment_commit": "4d1b48634d559e38369e85936f5c7a2190e24323",
            "final_pre_outcome_format_commit": "089052161aaf8a3c37fc7346418596ac1e650552"}
    old_sources = ((args.old_runner, "run_before_prefix.py",
                    "3d90a2d5656db2ac44681a44d517a74e62eb8feb1e31c9cb7a59ceb914c27335"),
                   (args.old_protocol, "PROTOCOL_before_prefix.md",
                    "9d6f222efca7dd149aa556212e3fef1cfa48d2b87ef87f6973cbb753b5ddec79"))
    for source, filename, expected in old_sources:
        data = source.read_bytes()
        assert sha(data) == expected
        (evidence / "source" / filename).write_bytes(data)
        meta["historical_source_sha256"][filename] = expected
    for filename in ("PROTOCOL.md", "make_inputs.py", "run.py", "verify.py",
                     "archive.py", "replay_all.py", "input_manifest.json",
                     "input_amendment.json", "target_points_n37.jsonl",
                     "target_points_n41.jsonl", "target_points_n41_R12.jsonl"):
        meta["method_source_sha256"][filename] = sha((HERE / filename).read_bytes())
    for n, r in ARMS:
        for mode in ("oracle", "extractor"):
            name = f"{mode}-n{n}-R{r}"
            source = args.runs_root / name
            dest = evidence / "runs" / name
            dest.mkdir()
            manifest = json.loads((source / "manifest.json").read_bytes())
            receipt = json.loads((source / "receipt.json").read_bytes())
            assert (manifest["mode"], manifest["n"], manifest["R"]) == (mode, n, r)
            assert receipt["returncode"] == 0 and not receipt["timed_out"] and not receipt["sampled_rss_stop"]
            for filename in ("manifest.json", "receipt.json", "producer.stderr.txt"):
                shutil.copyfile(source / filename, dest / filename)
            stdout = (source / "producer.stdout.jsonl").read_bytes()
            assert sha(stdout) == receipt["stdout_sha256"]
            (dest / "producer.stdout.jsonl.gz").write_bytes(frozen_gzip(stdout))
            meta["run_names"].append(name)
            meta["stdout_sha256"][name] = sha(stdout)
            rel_source = manifest["source_path"]
            raw_source = (REPO / rel_source).read_bytes()
            assert sha(raw_source) == manifest["source_sha256"]
            snap = evidence / "source" / (Path(rel_source).name + ".gz")
            compressed = frozen_gzip(raw_source)
            if snap.exists():
                assert snap.read_bytes() == compressed
            else:
                snap.write_bytes(compressed)
            meta["source_sha256"][rel_source] = sha(raw_source)
        validation = args.validation_root / f"independent-n{n}-R{r}.json"
        data = validation.read_bytes()
        json.loads(data)
        filename = validation.name
        (evidence / "validations" / filename).write_bytes(data)
        meta["validation_sha256"][filename] = sha(data)
    for name in ATTEMPTS:
        source = args.runs_root / "attempts" / name
        dest = evidence / "attempts" / name
        dest.mkdir()
        manifest = json.loads((source / "manifest.json").read_bytes())
        receipt = json.loads((source / "receipt.json").read_bytes())
        for filename in ("manifest.json", "receipt.json", "producer.stderr.txt"):
            shutil.copyfile(source / filename, dest / filename)
        stdout = (source / "producer.stdout.jsonl").read_bytes()
        assert sha(stdout) == receipt["stdout_sha256"]
        (dest / "producer.stdout.jsonl.gz").write_bytes(frozen_gzip(stdout))
        if name == "extractor-n41-R12-512-invalid":
            observed = json.loads(stdout)
            assert manifest["count"] == 128
            assert observed["compact_orbit_point_batch"]["targets_requested"] == 512
            amendment = json.loads((HERE / "input_amendment.json").read_bytes())
            assert sha((source / "manifest.json").read_bytes()) == amendment["failed_512_attempt_manifest_sha256"]
            assert sha((source / "receipt.json").read_bytes()) == amendment["failed_512_attempt_receipt_sha256"]
            assert sha(stdout) == amendment["failed_512_attempt_stdout_sha256"]
        elif name == "oracle-n41-R12-full-file-prefix":
            observed = [json.loads(line) for line in stdout.splitlines()]
            assert manifest["count"] == 128 and len(observed) == 129
        else:
            assert name == "extractor-n37-R3-host-contended"
            assert manifest["count"] == 512
    (evidence / "archive_manifest.json").write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    files = sorted(path for path in evidence.rglob("*") if path.is_file())
    sums = "".join(f"{sha(path.read_bytes())}  {path.relative_to(evidence)}\n" for path in files)
    (evidence / "SHA256SUMS").write_text(sums)
    print(json.dumps({"archive_manifest_sha256": sha((evidence / "archive_manifest.json").read_bytes()),
                      "SHA256SUMS_sha256": sha((evidence / "SHA256SUMS").read_bytes()),
                      "run_names": meta["run_names"]}, sort_keys=True))


if __name__ == "__main__":
    main()

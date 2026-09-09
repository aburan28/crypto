"""Validate and measure the packed multiplier on one Modal GPU allocation.

ECC_GPU=RTX-PRO-6000 ECC_PACKED_SINGLE_PRODUCT=1 modal run packed_audit.py \
    --output packed-audit.json
"""
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import time

import modal
# Modal re-imports this entry point as /root/packed_audit.py; dependencies live
# in the source tree copied into client.image.
if not modal.is_local():
    sys.path.insert(0, "/root/ecc2k130")
import modal_app as client

app = modal.App("ecc2k130-packed-audit")


@app.function(image=client.image, gpu=client.DEFAULT_GPU, timeout=1800,
              volumes={"/data": client.volume})
def runAudit(minBlocks=4, repeats=3):
    result = dict(valid=False, minBlocks=minBlocks, repeats=repeats,
                  batch=32, blockThreads=128, steps=1024, launches=32)
    try:
        if minBlocks <= 0 or repeats <= 0:
            raise ValueError("min-blocks and repeats must be positive")
        ok, build = client.buildFor(32, 128, 0, minBlocks=minBlocks)
        result["build"] = build
        if not ok:
            raise RuntimeError("CUDA build failed")
        result["identity"] = client.benchmarkIdentity(packed=True)

        def run(command, timeout):
            p = subprocess.run(command, cwd=client.REMOTE, capture_output=True,
                               text=True, timeout=timeout)
            row = dict(command=command, returncode=p.returncode,
                       output=p.stdout + p.stderr)
            print(json.dumps(row), flush=True)
            return row

        result["integration"] = run(
            ["python3", "codegen/testpackedclient.py", "./ecc2k130"], 600)
        if result["integration"]["returncode"]:
            raise RuntimeError("packed GPU integration failed")

        result["benchmark"] = client.measureBench(1024, 32, 0, False, repeats, packed=True)
        if not result["benchmark"]["valid"]:
            raise RuntimeError("throughput benchmark failed")

        # Time real collection with CPU trail replay disabled only after the
        # integration test has independently replayed reports. Each sample
        # starts with its own corpus, so no earlier reports enter its count.
        samples = result["collection"] = []
        for repeat in range(repeats):
            with tempfile.TemporaryDirectory() as directory:
                corpus = Path(directory) / "points.bin"
                row = run(["./ecc2k130", "--packed", "--curve", "131",
                           "--dp-weight", "34", "--steps", "1024", "--launches", "32",
                           "--verify", "0", "--dp-file", str(corpus)], 180)
                sample = client.benchResult(row["command"], row["returncode"], row["output"])
                records = corpus.stat().st_size if corpus.exists() else 0
                final = re.findall(r"finished:.*?, (\d+) distinguished points "
                                   r"\(0 verified against the reference, (\d+) dropped\)",
                                   row["output"])
                sample["corpusBytes"] = records
                sample["corpusRecords"] = records // 32
                sample["valid"] = (sample["valid"] and len(final) == 1
                                   and int(final[0][0]) > 0 and int(final[0][1]) == 0
                                   and records == int(final[0][0]) * 32)
                if not sample["valid"]:
                    sample["rate"] = 0.0
                samples.append(sample)
                if not sample["valid"]:
                    raise RuntimeError("collection failed or report counts disagree")
        result["collectionSummary"] = client.summarizeSamples(samples)
        result["valid"] = True
    except Exception as exc:
        result["error"] = str(exc)
    artifact = Path("/data/packed-audit") / (str(time.time_ns()) + ".json")
    artifact.parent.mkdir(parents=True, exist_ok=True)
    result["remoteArtifact"] = str(artifact)
    artifact.write_text(json.dumps(result, indent=2) + "\n")
    client.volume.commit()
    return result


@app.local_entrypoint()
def main(output: str = "packed-audit.json", min_blocks: int = 4, repeats: int = 3):
    result = runAudit.remote(minBlocks=min_blocks, repeats=repeats)
    Path(output).write_text(json.dumps(result, indent=2) + "\n")
    print(f"Audit saved to {output}")
    if not result["valid"]:
        raise RuntimeError(result.get("error", "packed audit failed"))
    print(f"Benchmark median: {result['benchmark']['rate']:.3f} M iterations/s")
    print(f"Collection median: {result['collectionSummary']['rate']:.3f} M iterations/s")

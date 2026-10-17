"""Resumable checked-Sage reduction of all 9,784 new quotient classes."""

import argparse
import concurrent.futures
import hashlib
import json
import subprocess
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
SAGE = Path("/Volumes/SSD990/cryptanalysis/sage")
SOURCE = HERE / "line_new_quotient_reduction.py"
INPUT = HERE / "type_i_line_full_new_quotient_points.json"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_chunk(job):
    start, count, directory = job
    name = f"{start:05d}-{start + count - 1:05d}"
    output = directory / f"{name}.json"
    log = directory / f"{name}.log"
    error = directory / f"{name}.stderr"
    if output.exists() and log.exists():
        try:
            prior = json.loads(output.read_text())
            if (prior["start_index"] == start and prior["input_classes"] ==
                    prior["reduced_classes"] == count and
                    prior["source_sha256"] == sha(SOURCE) and
                    prior["input_sha256"] == sha(INPUT) and
                    "NEW_QUOTIENT_REDUCTION_PASS" in log.read_text()):
                return {"start": start, "count": count,
                        "receipt_sha256": sha(output), "resumed": True}
        except (KeyError, ValueError):
            pass
        raise RuntimeError(f"stale or incomplete existing chunk: {name}")
    command = [str(SAGE), "-python", str(SOURCE), "--start", str(start),
               "--count", str(count), "--out", str(output)]
    with log.open("w") as out, error.open("w") as err:
        result = subprocess.run(command, cwd=HERE.parents[1], stdout=out, stderr=err)
    if result.returncode:
        raise RuntimeError(f"chunk {name} exit {result.returncode}; inspect {error}")
    receipt = json.loads(output.read_text())
    assert receipt["start_index"] == start
    assert receipt["input_classes"] == receipt["reduced_classes"] == count
    assert receipt["unreduced_classes"] == 0
    assert len(receipt["rows"]) == len(receipt["representatives"]) == count
    return {"start": start, "count": count,
            "receipt_sha256": sha(output), "resumed": False}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--chunk-size", type=int, default=512)
    parser.add_argument("--directory", type=Path,
                        default=HERE / "type_i_line_new_all_chunks")
    parser.add_argument("--out", type=Path,
                        default=HERE / "type_i_line_new_all_reduction.json")
    args = parser.parse_args()
    assert SAGE.is_file() and SOURCE.is_file()
    assert 1 <= args.workers <= 8 and 1 <= args.chunk_size <= 9784
    points = json.loads(INPUT.read_text())
    assert len(points) == 9784
    args.directory.mkdir(parents=True, exist_ok=True)
    jobs = [(start, min(args.chunk_size, len(points) - start), args.directory)
            for start in range(0, len(points), args.chunk_size)]
    start_ns = time.monotonic_ns()
    chunks = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_chunk, job): job for job in jobs}
        for future in concurrent.futures.as_completed(futures):
            chunk = future.result()
            chunks.append(chunk)
            print("NEW_ALL_CHUNK_PASS", chunk["start"], chunk["count"],
                  len(chunks), len(jobs), flush=True)
    chunks.sort(key=lambda item: item["start"])
    assert [item["start"] for item in chunks] == [job[0] for job in jobs]
    assert sum(item["count"] for item in chunks) == len(points)
    receipt = {
        "status": "all_9784_new_quotient_classes_reduced",
        "field_prime": str(int("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF", 16)),
        "input_classes": len(points), "reduced_classes": len(points),
        "chunk_size": args.chunk_size, "workers": args.workers,
        "chunks": chunks,
        "elapsed_ns_unisolated": time.monotonic_ns() - start_ns,
        "source_sha256": sha(Path(__file__)),
        "worker_source_sha256": sha(SOURCE),
        "input_sha256": sha(INPUT),
        "runtime_info_sha256": sha(HERE / "type_i_line_new64_runtime_info.json"),
        "boundary": "All chunked fixed-origin reductions passed the checked Sage constructor; independent native cubic-cut replay is a separate gate.",
    }
    args.out.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print("NEW_ALL_REDUCTION_PASS", json.dumps({
        "classes": len(points), "chunks": len(chunks),
        "output_sha256": sha(args.out),
    }, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()

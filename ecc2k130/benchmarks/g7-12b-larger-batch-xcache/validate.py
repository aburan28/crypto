#!/usr/bin/env python3
from pathlib import Path
import array
import hashlib
import json
import os
import struct
import subprocess
import time

from build import CONFIGS, ROOT, OUT

SELECTED = ROOT / "build/ecc2k130-local-packed"
LOGICAL_WALKS = 12288

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def run(command, log, timeout=300):
    assert not log.exists()
    start = time.monotonic()
    with log.open("x") as stream:
        process = subprocess.run([str(x) for x in command], stdout=stream,
                                 stderr=subprocess.STDOUT, timeout=timeout,
                                 env=dict(os.environ, OMP_NUM_THREADS="8",
                                          CUDA_DISABLE_PTX_JIT="1"))
    return process.returncode, time.monotonic() - start, log.read_text()

def normalized_checkpoint(path):
    data = path.read_bytes()
    magic, version, degree, threads, batch, lanes, runid, iters = struct.unpack_from(
        "<8s6IQ", data)
    assert (magic, version, degree, lanes) == (b"ECC2K130", 2, 131, 1)
    count = threads * batch
    assert count == LOGICAL_WALKS and len(data) == 40 + count * 60
    if batch == 16:
        return data
    reference_workers = LOGICAL_WALKS // 16
    result = bytearray(struct.pack("<8s6IQ", magic, version, degree,
                                   reference_workers, 16, lanes, runid, iters))
    offset = 40
    for _ in range(2):
        source = array.array("I")
        source.frombytes(data[offset:offset + count * 20])
        dest = array.array("I", [0]) * (count * 5)
        for i in range(count):
            for word in range(5):
                dest[((i // reference_workers) * 5 + word) * reference_workers
                     + i % reference_workers] = source[
                         ((i // threads) * 5 + word) * threads + i % threads]
        result.extend(dest.tobytes())
        offset += count * 20
    result.extend(data[offset:])
    return bytes(result)

def client(binary, batch, label, runid):
    workers = LOGICAL_WALKS // batch
    checkpoint = ROOT / "build" / f"g7-xcache-{label}-{runid}.ckpt"
    corpus = ROOT / "build" / f"g7-xcache-{label}-{runid}.dp"
    checkpoint.unlink(missing_ok=True)
    corpus.unlink(missing_ok=True)
    command = [binary, "--packed", "--threads", workers, "--steps", 64,
               "--launches", 4, "--run-id", runid, "--checkpoint", checkpoint,
               "--dp-weight", 52, "--verify", 16, "--dp-file", corpus]
    code, elapsed, text = run(command, OUT / f"{label}-{runid}.log")
    assert code == 0 and "MISMATCH" not in text and "0 dropped" in text
    assert "16 verified against the reference" in text
    dp = corpus.read_bytes()
    assert len(dp) % 32 == 0
    return {"checkpoint_sha256": hashlib.sha256(normalized_checkpoint(checkpoint)).hexdigest(),
            "dp_multiset_sha256": hashlib.sha256(b"".join(sorted(
                dp[i:i + 32] for i in range(0, len(dp), 32)))).hexdigest(),
            "dp_records": len(dp) // 32, "elapsed_seconds": elapsed}

def main():
    builtin_rows = []
    for label in ["selected", *CONFIGS]:
        binary = SELECTED if label == "selected" else ROOT / "build" / f"g7-{label}"
        code, elapsed, _ = run([binary, "--test"], OUT / f"{label}-builtin-test.log", 600)
        builtin_rows.append({"label": label, "returncode": code,
                             "elapsed_seconds": elapsed,
                             "log_sha256": sha(OUT / f"{label}-builtin-test.log")})
    assert all(row["returncode"] == 0 for row in builtin_rows)
    assert len({row["log_sha256"] for row in builtin_rows}) == 1
    (OUT / "builtin-tests.json").write_text(json.dumps(builtin_rows, indent=2) + "\n")

    state_rows = []
    for runid in (0, 139):
        reference = client(SELECTED, 16, "selected", runid)
        for label, (batch, _, _) in CONFIGS.items():
            candidate = client(ROOT / "build" / f"g7-{label}", batch, label, runid)
            assert candidate["checkpoint_sha256"] == reference["checkpoint_sha256"]
            assert candidate["dp_multiset_sha256"] == reference["dp_multiset_sha256"]
            assert candidate["dp_records"] == reference["dp_records"]
            state_rows.append({"label": label, "run_id": runid,
                               "reference": reference, "candidate": candidate})
    (OUT / "state-validation.json").write_text(json.dumps(state_rows, indent=2) + "\n")
    print(json.dumps({"builtins": builtin_rows, "state_checks": len(state_rows)}, indent=2))


if __name__ == "__main__":
    main()

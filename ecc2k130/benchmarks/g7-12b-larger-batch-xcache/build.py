#!/usr/bin/env python3
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib
import json
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
SRC = ROOT / "build/g7-larger-batch-xcache-src"
NVCC = ROOT / "build/assembler134-toolchain/bin/nvcc"
FLAGS = {
    "SINGLE_PRODUCT": 0, "CACHE_DENOM": 1, "BY_VALUE": 1,
    "PERM_SIGMA": 3, "POLY_CHAIN": 1, "UNROLL_INV": 1,
    "PAIR_PRODUCTS": 1, "POLY_STATE": 1, "DIRECT_REDUCE": 1,
    "GENERATED_PRODUCT": 1, "CLMAD": 1, "WEIGHTED_PREFIX": 2,
    "COMPACT_STATE": 1, "SHARED_SIGMA": 1, "STATE_TILE": 256,
    "INLINE": 3, "BLOCK_INVERSE": 1, "BATCH_SPLIT": 2,
    "PARTIAL_SIGMA": 1, "LAST_SLOT_CACHE": 1, "TAIL_LAYOUT": 1,
    "SIGMA_ORDER": 1, "BYTE_SIGMA": 1, "FAST_CONVERT": 1,
    "LOGICAL_PAIR_INVERSE": 0,
}
CONFIGS = {
    "batch24-cache4-min2": (24, 4, 2),
    "batch24-cache6-min2": (24, 6, 2),
    "batch32-cache4-min2": (32, 4, 2),
    "batch32-cache8-min2": (32, 8, 2),
}

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def command(label, binary):
    batch, slots, minblocks = CONFIGS[label]
    flags = FLAGS | {"SHARED_X_SLOTS": slots}
    return [str(NVCC), "-O3", "-std=c++17", "-gencode",
            "arch=compute_120,code=sm_120", f"-DECC_BATCH={batch}",
            "-DECC_THREADS=256", f"-DECC_MINBLOCKS={minblocks}",
            *[f"-DECC_PACKED_{key}={value}" for key, value in flags.items()],
            "-Xptxas", "-v", "-lineinfo", "-Xcompiler", "-O3",
            "-Xcompiler", "-fopenmp", str(SRC / "src/main.cu"),
            "-o", str(binary), "-lgomp"]

def build(label):
    binary = ROOT / "build" / f"g7-{label}"
    log = OUT / f"{label}-build.log"
    assert not binary.exists() and not log.exists()
    cmd = command(label, binary)
    start = time.monotonic()
    with log.open("x") as stream:
        process = subprocess.run(cmd, stdout=stream, stderr=subprocess.STDOUT,
                                 timeout=900)
    row = {"label": label, "batch": CONFIGS[label][0],
           "shared_x_slots": CONFIGS[label][1],
           "minimum_blocks_per_sm": CONFIGS[label][2], "command": cmd,
           "returncode": process.returncode,
           "elapsed_seconds": time.monotonic() - start,
           "log_sha256": sha(log),
           "binary": str(binary.relative_to(ROOT)),
           "binary_sha256": sha(binary) if binary.exists() else None}
    (OUT / f"{label}-build.json").write_text(json.dumps(row, indent=2) + "\n")
    print(json.dumps(row), flush=True)
    return row

if __name__ == "__main__":
    with ThreadPoolExecutor(max_workers=2) as pool:
        rows = list(pool.map(build, CONFIGS))
    assert all(row["returncode"] == 0 for row in rows)
    (OUT / "builds.json").write_text(json.dumps(rows, indent=2) + "\n")

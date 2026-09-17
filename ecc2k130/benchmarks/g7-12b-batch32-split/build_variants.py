#!/usr/bin/env python3
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib, json, subprocess, time

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
SRC = ROOT / "build/g7-batch32-split-src"
NVCC = ROOT / "build/assembler134-toolchain/bin/nvcc"
FLAGS = {
    "SINGLE_PRODUCT": 0, "CACHE_DENOM": 1, "BY_VALUE": 1,
    "PERM_SIGMA": 3, "POLY_CHAIN": 1, "UNROLL_INV": 1,
    "PAIR_PRODUCTS": 1, "POLY_STATE": 1, "DIRECT_REDUCE": 1,
    "GENERATED_PRODUCT": 1, "CLMAD": 1, "WEIGHTED_PREFIX": 2,
    "COMPACT_STATE": 1, "SHARED_SIGMA": 1, "STATE_TILE": 256,
    "INLINE": 3, "BLOCK_INVERSE": 1, "BATCH_SPLIT": 2,
    "PARTIAL_SIGMA": 1, "LAST_SLOT_CACHE": 1, "TAIL_LAYOUT": 1,
    "SIGMA_ORDER": 1, "BYTE_SIGMA": 1, "SHARED_X_SLOTS": 4,
    "FAST_CONVERT": 1, "LOGICAL_PAIR_INVERSE": 0,
}

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def command(batch, binary):
    return [str(NVCC), "-O3", "-std=c++17", "-gencode",
            "arch=compute_120,code=sm_120", f"-DECC_BATCH={batch}",
            "-DECC_THREADS=256", "-DECC_MINBLOCKS=3",
            *[f"-DECC_PACKED_{key}={value}" for key, value in FLAGS.items()],
            "-Xptxas", "-v", "-lineinfo", "-Xcompiler", "-O3",
            "-Xcompiler", "-fopenmp", str(SRC / "src/main.cu"),
            "-o", str(binary), "-lgomp"]

def build(batch):
    label=f"batch{batch}"
    binary=ROOT / "build" / f"g7-split-{label}"
    log=OUT / f"{label}-build.log"
    receipt=OUT / f"{label}-build.json"
    assert not binary.exists() and not log.exists() and not receipt.exists()
    cmd=command(batch,binary); start=time.monotonic()
    with log.open("x") as stream:
        process=subprocess.run(cmd,stdout=stream,stderr=subprocess.STDOUT,timeout=600)
    row={"batch":batch,"command":cmd,"returncode":process.returncode,
         "elapsed_seconds":time.monotonic()-start,"log_sha256":sha(log),
         "binary":str(binary.relative_to(ROOT)),
         "binary_sha256":sha(binary) if binary.exists() else None}
    receipt.write_text(json.dumps(row,indent=2)+"\n")
    return row

with ThreadPoolExecutor(max_workers=2) as pool:
    rows=list(pool.map(build,[16,24,32]))
print(json.dumps(rows,indent=2))
assert all(row["returncode"]==0 for row in rows)

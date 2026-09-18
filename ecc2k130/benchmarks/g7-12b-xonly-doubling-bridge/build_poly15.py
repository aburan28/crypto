#!/usr/bin/env python3
"""Build mod72 control, arithmetic-only diagnostic, and poly12 complete-map binaries."""
from pathlib import Path
import hashlib, json, subprocess, time

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
NVCC = ROOT / "build/assembler134-toolchain/bin/nvcc"
SRC = ROOT / "src/main.cu"
BASE = {
    "SINGLE_PRODUCT": 0, "CACHE_DENOM": 1, "BY_VALUE": 1, "PERM_SIGMA": 3,
    "POLY_CHAIN": 1, "UNROLL_INV": 1, "PAIR_PRODUCTS": 1, "POLY_STATE": 1,
    "DIRECT_REDUCE": 1, "GENERATED_PRODUCT": 1, "CLMAD": 1, "WEIGHTED_PREFIX": 2,
    "COMPACT_STATE": 1, "SHARED_SIGMA": 1, "STATE_TILE": 256, "INLINE": 3,
    "BLOCK_INVERSE": 1, "BATCH_SPLIT": 2, "PARTIAL_SIGMA": 1, "LAST_SLOT_CACHE": 1,
    "TAIL_LAYOUT": 1, "SIGMA_ORDER": 1, "BYTE_SIGMA": 1, "FAST_CONVERT": 1,
    "SHARED_X_SLOTS": 4, "XONLY_23": 1, "XONLY_DOUBLE_ONLY": 1,
    "XONLY_BRIDGE1_COMMON": 1, "XONLY_SKIP_EMPTY_BRIDGE": 0,
}

VARIANTS = {
    "control": {
        "XONLY_BRIDGE3": 1, "XONLY_BRIDGE_MOD72": 1,
        "XONLY_ARITHMETIC_ONLY": 0, "XONLY_POLY_SELECT": 0, "XONLY_POLY_DP_CONVERT": 1,
    },
    "arith": {
        "XONLY_BRIDGE3": 0, "XONLY_BRIDGE_MOD72": 0,
        "XONLY_ARITHMETIC_ONLY": 1, "XONLY_POLY_SELECT": 0, "XONLY_POLY_DP_CONVERT": 1,
    },
    "poly12": {
        "XONLY_BRIDGE3": 1, "XONLY_BRIDGE_MOD72": 0,
        "XONLY_ARITHMETIC_ONLY": 0, "XONLY_POLY_SELECT": 1, "XONLY_POLY_DP_CONVERT": 0,
        "XONLY_POLY_SELECT_MASK": 4095, "XONLY_POLY_SELECT_TARGET": 14,
    },
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def build(label):
    binary = ROOT / "build" / f"g7-xonly-15b-{label}"
    log = OUT / f"poly15-{label}-build.log"
    flags = dict(BASE, **VARIANTS[label])
    cmd = [
        str(NVCC), "-O3", "-std=c++17",
        "-gencode", "arch=compute_120,code=sm_120",
        "-DECC_BATCH=16", "-DECC_THREADS=256", "-DECC_MINBLOCKS=3",
        *[f"-DECC_PACKED_{k}={v}" for k, v in flags.items()],
        "-Xptxas", "-v", "-lineinfo",
        "-Xcompiler", "-O3", "-Xcompiler", "-fopenmp",
        str(SRC), "-o", str(binary), "-lgomp",
    ]
    started = time.monotonic()
    with log.open("w") as stream:
        result = subprocess.run(cmd, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT, timeout=900)
    row = {
        "label": label,
        "flags": flags,
        "command": cmd,
        "returncode": result.returncode,
        "elapsed_seconds": time.monotonic() - started,
        "binary": str(binary.relative_to(ROOT)),
        "binary_sha256": sha(binary) if binary.exists() else None,
        "log_sha256": sha(log),
        "source_sha256": sha(SRC),
        "kernel_sha256": sha(ROOT / "include/packedxonly23.cuh"),
    }
    (OUT / f"poly15-{label}-build.json").write_text(json.dumps(row, indent=2) + "\n")
    print(json.dumps(row, indent=2), flush=True)
    if result.returncode:
        raise SystemExit(result.returncode)
    return row


if __name__ == "__main__":
    (ROOT / "build").mkdir(exist_ok=True)
    for label in ("control", "arith", "poly12"):
        build(label)

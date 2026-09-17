#!/usr/bin/env python3
from pathlib import Path
import difflib
import shutil

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
BASE = ROOT / "build/g7-batch32-split-src"
SRC = ROOT / "build/g7-larger-batch-xcache-src"

assert BASE.is_dir()
assert not SRC.exists()
shutil.copytree(BASE, SRC, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
path = SRC / "include/packedkernels.cuh"
before = path.read_text()
old = '''#if ECC_PACKED_SHARED_X_SLOTS != 0 && ECC_PACKED_SHARED_X_SLOTS != 2 && ECC_PACKED_SHARED_X_SLOTS != 4
#error "SHARED_X_SLOTS must be 0, 2 or 4"
'''
new = '''#if ECC_PACKED_SHARED_X_SLOTS < 0 || ECC_PACKED_SHARED_X_SLOTS > 8
#error "SHARED_X_SLOTS must be from 0 through 8"
'''
assert before.count(old) == 1
after = before.replace(old, new)
path.write_text(after)
diff = difflib.unified_diff(
    before.splitlines(True), after.splitlines(True),
    fromfile="base/include/packedkernels.cuh",
    tofile="candidate/include/packedkernels.cuh")
(OUT / "candidate.patch").write_text("".join(diff))
print(SRC)

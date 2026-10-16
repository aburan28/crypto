#!/usr/bin/env python3
"""Check walk_step_problem.Model bit for bit against include/ref.h.

    python3 ecc2k130/kernelagent/check_reference.py [--n 4096]

Builds check_reference.cpp with g++ into a temporary directory, feeds it the
same inputs as the PyTorch reference (random elements plus zero, one and
single-bit edge cases) and fails on the first differing limb.  Run it after
any change to walk_step_problem.py: it is the only thing that ties the
problem KernelAgent optimizes to the walk the GPU client runs.
"""
import argparse
import pathlib
import subprocess
import sys
import tempfile

import torch

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import walk_step_problem as wsp  # noqa: E402

MASK64 = (1 << 64) - 1


def edge_cases() -> torch.Tensor:
    rows = [[0, 0, 0], [-1, -1, 7]]  # zero and one (all-ones normal-basis vector)
    for b in (0, 1, 63, 64, 127, 128, 130):
        limb = [0, 0, 0]
        limb[b // 64] = 1 << (b % 64)
        if limb[b // 64] >= 1 << 63:
            limb[b // 64] -= 1 << 64
        rows.append(limb)
    return torch.tensor(rows, dtype=torch.long)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=4096)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    g = torch.Generator().manual_seed(args.seed)
    e = edge_cases()
    xs = torch.cat([e.repeat_interleave(len(e), 0), wsp._random_elements(args.n, g)])
    ys = torch.cat([e.repeat(len(e), 1), wsp._random_elements(args.n, g)])

    got = wsp.Model()(xs, ys).reshape(len(xs), 6)

    with tempfile.TemporaryDirectory() as tmp:
        exe = pathlib.Path(tmp) / "check_reference"
        subprocess.run(["g++", "-O2", "-std=c++17", "-o", str(exe), str(HERE / "check_reference.cpp")],
                       check=True)
        lines = "".join("%x %x %x %x %x %x\n" % tuple(v & MASK64 for v in row)
                        for row in torch.cat([xs, ys], 1).tolist())
        out = subprocess.run([str(exe)], input=lines, capture_output=True, text=True, check=True).stdout

    want = [[int(t, 16) for t in line.split()] for line in out.splitlines()]
    if len(want) != len(xs):
        print("oracle returned %d rows for %d inputs" % (len(want), len(xs)), file=sys.stderr)
        return 1
    for i, (w, gr) in enumerate(zip(want, got.tolist())):
        if w != [v & MASK64 for v in gr]:
            print("mismatch at row %d: x=%s y=%s ref.h=%s torch=%s"
                  % (i, xs[i].tolist(), ys[i].tolist(), [hex(v) for v in w],
                     [hex(v & MASK64) for v in gr]), file=sys.stderr)
            return 1
    print("walk_step_problem matches include/ref.h on %d inputs (%d edge cases)"
          % (len(xs), len(e) ** 2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

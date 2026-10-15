#!/usr/bin/env python3
"""Time a KernelAgent kernel for the ECC2K-130 walk step against the reference.

    python3 bench_candidate.py --kernel /results/kernel.py --out /results [--reps 30]

Before any timing the kernel must reproduce walk_step_problem.Model exactly
on fresh inputs (a seed the agent's test never used) at every batch size
timed; a mismatch stops the run with no throughput recorded.  Each batch
size is then timed with CUDA events, kernel and reference interleaved per
repetition, and reported as median walk steps per second and the median of
the per-repetition speedup ratios.

The reference is a readable PyTorch program, so a speedup over it says
nothing yet about the production CUDA client: that runs a different
iteration (the x-only doubling map with a sigma^3 bridge, see ../README.md)
and is measured in complete scalar updates/s.  A kernel worth keeping is
ported and measured against the production kernel with the repository's
paired protocol before any claim.
"""
import argparse
import importlib.util
import json
import pathlib
import statistics
import sys

import torch

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import walk_step_problem as wsp  # noqa: E402


def load_kernel(path: pathlib.Path):
    spec = importlib.util.spec_from_file_location("candidate_kernel", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    if hasattr(mod, "kernel_function"):
        return mod.kernel_function
    for name in ("ModelNew", "Model"):
        if hasattr(mod, name):
            return getattr(mod, name)().cuda()
    raise SystemExit("%s defines no kernel_function, ModelNew or Model" % path)


def time_ms(fn, x, y) -> float:
    start, end = torch.cuda.Event(enable_timing=True), torch.cuda.Event(enable_timing=True)
    start.record()
    fn(x, y)
    end.record()
    torch.cuda.synchronize()
    return start.elapsed_time(end)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kernel", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--reps", type=int, default=30)
    ap.add_argument("--batches", default="65536,262144,1048576")
    ap.add_argument("--seed", type=int, default=20261014)
    args = ap.parse_args()

    kernel = load_kernel(pathlib.Path(args.kernel))
    ref = wsp.Model().cuda()
    g = torch.Generator().manual_seed(args.seed)
    report = {"gpu": torch.cuda.get_device_name(), "torch": torch.__version__,
              "seed": args.seed, "reps": args.reps, "batches": []}

    for n in (int(b) for b in args.batches.split(",")):
        x = wsp._random_elements(n, g).cuda()
        y = wsp._random_elements(n, g).cuda()
        want, got = ref(x, y), kernel(x, y)
        if not (isinstance(got, torch.Tensor) and got.shape == want.shape and torch.equal(got, want)):
            report["verdict"] = "mismatch at batch %d" % n
            print(report["verdict"], file=sys.stderr)
            break
        for _ in range(3):
            ref(x, y)
            kernel(x, y)
        ref_ms, ker_ms = [], []
        for _ in range(args.reps):
            ref_ms.append(time_ms(ref, x, y))
            ker_ms.append(time_ms(kernel, x, y))
        ratios = sorted(r / k for r, k in zip(ref_ms, ker_ms))
        row = {"batch": n,
               "reference_ms": statistics.median(ref_ms), "kernel_ms": statistics.median(ker_ms),
               "reference_steps_per_s": n / statistics.median(ref_ms) * 1e3,
               "kernel_steps_per_s": n / statistics.median(ker_ms) * 1e3,
               "speedup_median": statistics.median(ratios),
               "speedup_min": ratios[0], "speedup_max": ratios[-1]}
        report["batches"].append(row)
        print(json.dumps(row), flush=True)
    else:
        report["verdict"] = "exact on all batches"

    pathlib.Path(args.out, "bench.json").write_text(json.dumps(report, indent=2))
    return 0 if report["verdict"] == "exact on all batches" else 1


if __name__ == "__main__":
    raise SystemExit(main())

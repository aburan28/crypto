"""Run the secp256k1 GPU kernels on Modal, on an actual GPU.

    ECC_GPU=H100 modal run modal_app.py::selftest
    ECC_GPU=H100 modal run modal_app.py::bench
    ECC_GPU=H100 modal run modal_app.py::tune

The reason this file exists is `selftest`.  `FP_PTX=1` replaces seven
multiprecision primitives in `fp256.cuh` with hand-written carry chains worth
55% of a field multiply, and it is off by default because that assembly has
never executed.  Every asm block is guarded on `__CUDA_ARCH__`, so the host
always compiles the portable path -- which is exactly what makes
`./bench selftest` built with `-DFP_PTX=1` a differential test of the
assembly against the portable code, and what means no CPU-only test can
substitute for it.  `ptx_asm_check.py` narrows that gap without hardware;
this closes it.

`selftest` builds both configurations and runs each against the host
reference.  It is the one command standing between `FP_PTX` and being the
default, and it exits non-zero if either configuration disagrees.

`bench` then answers the question the instruction counts only predict: 210
PTX instructions against 472 is a 55% reduction on paper, but PTX is not
SASS and the kernels are occupancy-bound, so the throughput ratio is an
open question until measured.

`tune` sweeps RHO_MIN_BLOCKS, which `ptxas` says trades registers for
resident warps (224 regs at 12.5% occupancy unset, 128 at 25% at four
blocks/SM).  Which point wins is not decidable offline.

The GPU comes from the ECC_GPU environment variable, read at import time and
baked into the function definitions, so it works on every Modal version.
Newer clients also accept --gpu on the entry points.
"""

import os
import pathlib
import re
import subprocess

import modal

CUDA_VERSION = "12.8.1"

DEFAULT_GPU = os.environ.get("ECC_GPU", "H100")

# Compute capability per Modal GPU type.  sm_120 is the Blackwell workstation
# part (RTX PRO 6000), sm_100 is B200/B300, sm_90 is H100/H200, sm_89 is
# L40S/L4, sm_86 is A10, sm_80 is A100, sm_75 is T4.
GPU_ARCH = {
    "T4": "75", "L4": "89", "L40S": "89", "A10": "86", "A10G": "86",
    "A100": "80", "A100-80GB": "80", "H100": "90", "H200": "90",
    "B200": "100", "B300": "100", "RTX-PRO-6000": "120",
}

REMOTE = "/root/ecc"
LOCAL = pathlib.Path(__file__).parent
HOUR = 60 * 60

image = (
    modal.Image.from_registry(
        f"nvidia/cuda:{CUDA_VERSION}-devel-ubuntu24.04", add_python="3.12"
    )
    .entrypoint([])
    .apt_install("build-essential")
    .add_local_dir(
        LOCAL,
        remote_path=REMOTE,
        copy=True,
        # Built binaries only -- do not glob `test_*`, that would drop
        # test_cpu.cpp and the image's `make test` with it.
        ignore=[
            "bench", "bench_0", "bench_1", "bench_tune",
            "test_secp_fast", "test_secp_mont", "test_toy_mont",
            "__pycache__", "*.pyc", "*.o",
        ],
    )
    # Generate the curve/vector headers once at image build.  They are pure
    # Python and identical on every architecture, so this need not be redone
    # per call.  The CPU suites and ptx_asm_check run here too: if they fail,
    # the image is broken and there is no point paying for a GPU to find out.
    .run_commands(
        f"cd {REMOTE} && make gen",
        f"cd {REMOTE} && make test",
    )
)

app = modal.App("ecc-secp256k1")


def onGpu(fn, gpu):
    """Point a function at a GPU type, if the client supports it."""
    if not gpu or gpu == DEFAULT_GPU:
        return fn
    if hasattr(fn, "with_options"):
        return fn.with_options(gpu=gpu)
    raise SystemExit(
        "This Modal client cannot change the GPU per call "
        "(Function.with_options was added in 0.72).\n"
        "Select the GPU through the environment instead:\n\n"
        f"    ECC_GPU={gpu} modal run modal_app.py::<entrypoint>\n"
    )


# ---------------------------------------------------------------------------


def sh(cmd, cwd=REMOTE, timeout=None):
    r = subprocess.run(
        cmd, shell=True, cwd=cwd, capture_output=True, text=True, timeout=timeout
    )
    return r.returncode, r.stdout + r.stderr


def gpuInfo():
    rc, out = sh(
        "nvidia-smi --query-gpu=name,compute_cap,clocks.max.sm,memory.total "
        "--format=csv,noheader"
    )
    if rc != 0 or not out.strip():
        return "unknown", "90"
    row = [f.strip() for f in out.strip().splitlines()[0].split(",")]
    name = row[0]
    cap = row[1].replace(".", "") if len(row) > 1 else "90"
    if len(row) > 3:
        name = f"{name} ({row[2]}, {row[3]})"
    return name, cap


def build(ptx, arch, minBlocks=0, out="bench"):
    """Build the bench binary for one configuration.  Returns (ok, log)."""
    rc, log = sh(
        f"make -B bench ARCH=sm_{arch} FP_PTX={int(ptx)} "
        f"MIN_BLOCKS={minBlocks} && mv bench {out}",
        timeout=30 * 60,
    )
    return rc == 0, log


def confName(ptx, minBlocks=0):
    s = "FP_PTX=1" if ptx else "portable"
    return s if not minBlocks else f"{s} minBlocks={minBlocks}"


# ---------------------------------------------------------------------------


@app.function(image=image, gpu=DEFAULT_GPU, timeout=2 * HOUR)
def runSelftest():
    """Build both configurations and check each against the host reference."""
    name, cap = gpuInfo()
    print(f"GPU: {name}, sm_{cap}\n")

    results = {}
    for ptx in (False, True):
        conf = confName(ptx)
        print(f"=== building {conf} ===")
        ok, log = build(ptx, cap, out=f"bench_{int(ptx)}")
        if not ok:
            print(log)
            results[conf] = ("BUILD FAILED", log)
            continue
        print(f"=== ./bench selftest ({conf}) ===")
        rc, out = sh(f"./bench_{int(ptx)} selftest", timeout=30 * 60)
        print(out)
        results[conf] = ("PASS" if rc == 0 else "FAIL", out)

    print("\n" + "=" * 60)
    for conf, (verdict, _) in results.items():
        print(f"  {conf:<24} {verdict}")

    bad = [c for c, (v, _) in results.items() if v != "PASS"]
    if bad:
        print(f"\n{', '.join(bad)} did not pass -- FP_PTX must stay off.")
        return {"ok": False, "gpu": name, "results": {k: v for k, (v, _) in results.items()}}

    print(
        "\nBoth configurations agree with the host reference on real hardware.\n"
        "That is the check OPTIMIZATION_BLACKWELL.md section 2 calls for: the\n"
        "inline assembly is now executed, not merely assembled and interpreted.\n"
        "FP_PTX=1 can be made the default on this architecture."
    )
    return {"ok": True, "gpu": name, "results": {k: v for k, (v, _) in results.items()}}


RATE = re.compile(r"([\d.]+)\s*([MG])\s*(?:it|step)s?/s", re.I)


def parseRate(text):
    """Best (largest) throughput figure in a bench output, in Mit/s."""
    best = 0.0
    for value, unit in RATE.findall(text):
        v = float(value) * (1000.0 if unit.upper() == "G" else 1.0)
        best = max(best, v)
    return best


@app.function(image=image, gpu=DEFAULT_GPU, timeout=2 * HOUR)
def runBench(walks: int = 0, iters: int = 0, w: int = 0, variant: str = ""):
    """Measure whether the 55% instruction reduction is a 55% throughput win."""
    name, cap = gpuInfo()
    print(f"GPU: {name}, sm_{cap}\n")

    opts = ""
    if walks:
        opts += f" --walks {walks}"
    if iters:
        opts += f" --iters {iters}"
    if w:
        opts += f" --w {w}"
    if variant:
        opts += f" --variant {variant}"

    rates = {}
    for ptx in (False, True):
        conf = confName(ptx)
        ok, log = build(ptx, cap, out=f"bench_{int(ptx)}")
        if not ok:
            print(f"{conf}: BUILD FAILED\n{log}")
            continue
        binary = f"./bench_{int(ptx)}"
        for cmd in ("field", "mul", f"rho{opts}"):
            print(f"=== {conf}: ./bench {cmd} ===")
            rc, out = sh(f"{binary} {cmd}", timeout=45 * 60)
            print(out)
            if rc == 0 and cmd.startswith("rho"):
                rates[conf] = parseRate(out)

    print("\n" + "=" * 60)
    for conf, r in rates.items():
        print(f"  {conf:<24} {r:10.1f} Mstep/s")
    if len(rates) == 2:
        a = rates.get("portable", 0.0)
        b = rates.get("FP_PTX=1", 0.0)
        if a > 0:
            print(f"\n  FP_PTX speedup: {b / a:.2f}x")
            print(
                "  Static PTX counts predict 472 -> 210 instructions (2.25x).\n"
                "  A smaller ratio here means the kernel is not instruction-issue\n"
                "  bound -- occupancy or memory, which `tune` probes."
            )
    return {"gpu": name, "rates": rates}


@app.function(image=image, gpu=DEFAULT_GPU, timeout=4 * HOUR)
def runTune(minBlocks: str = "0,3,4", ptxBoth: bool = True,
            walks: int = 0, iters: int = 0):
    """Sweep RHO_MIN_BLOCKS: registers against resident warps.

    ptxas offline says 224 registers / 12.5% occupancy unset, 168 / 18.8% at
    three blocks per SM, 128 / 25% at four.  Which is fastest depends on
    whether the kernel is issue-bound or latency-bound, and only a device
    settles that."""
    name, cap = gpuInfo()
    print(f"GPU: {name}, sm_{cap}\n")

    opts = ""
    if walks:
        opts += f" --walks {walks}"
    if iters:
        opts += f" --iters {iters}"

    ptxValues = (False, True) if ptxBoth else (True,)
    rows = []
    for ptx in ptxValues:
        for mb in [int(x) for x in minBlocks.split(",") if x.strip()]:
            conf = confName(ptx, mb)
            ok, log = build(ptx, cap, minBlocks=mb, out="bench_tune")
            if not ok:
                print(f"{conf}: BUILD FAILED\n{log}")
                rows.append((conf, 0.0))
                continue
            rc, out = sh(f"./bench_tune rho{opts}", timeout=45 * 60)
            rate = parseRate(out) if rc == 0 else 0.0
            print(f"=== {conf}: {rate:.1f} Mstep/s ===")
            print(out)
            rows.append((conf, rate))

    print("\n" + "=" * 60)
    rows.sort(key=lambda r: -r[1])
    for conf, rate in rows:
        print(f"  {conf:<34} {rate:10.1f} Mstep/s")
    if rows and rows[0][1] > 0:
        print(f"\n  best: {rows[0][0]}")
        print(
            "  Update the RHO_MIN_BLOCKS recommendation in\n"
            "  OPTIMIZATION_BLACKWELL.md section 3 with this, per architecture:\n"
            "  the offline register figures do not decide it."
        )
    return {"gpu": name, "rows": rows}


# ---------------------------------------------------------------------------


@app.local_entrypoint()
def selftest(gpu: str = ""):
    onGpu(runSelftest, gpu).remote()


@app.local_entrypoint()
def bench(gpu: str = "", walks: int = 0, iters: int = 0, w: int = 0,
          variant: str = ""):
    onGpu(runBench, gpu).remote(walks=walks, iters=iters, w=w, variant=variant)


@app.local_entrypoint()
def tune(gpu: str = "", min_blocks: str = "0,3,4", ptx_both: bool = True,
         walks: int = 0, iters: int = 0):
    onGpu(runTune, gpu).remote(
        minBlocks=min_blocks, ptxBoth=ptx_both, walks=walks, iters=iters
    )

"""Run the ECC2K-130 client on Modal GPUs.

    ECC_GPU=RTX-PRO-6000 modal run modal_app.py::validate
    ECC_GPU=RTX-PRO-6000 modal run modal_app.py::bench
    ECC_GPU=RTX-PRO-6000 modal run modal_app.py::autolab    # no GPU, picks candidates
    ECC_GPU=RTX-PRO-6000 modal run modal_app.py::autotune   # measures them
    ECC_GPU=RTX-PRO-6000 modal run modal_app.py::campaign   # both, in one call
    ECC_GPU=RTX-PRO-6000 modal run modal_app.py::profile    # Nsight Compute
    ECC_GPU=RTX-PRO-6000 modal run modal_app.py::search --curve 97 --hours 4
    ECC_GPU=RTX-PRO-6000 modal run modal_app.py::fanout --curve 97 --count 8 --hours 4

Curve 97 is ECC2K-95, a 2^44 iteration problem: feasible in GPU-hours, and its
answer has been public since Harley's group solved it in 1998, so a recovered
logarithm can be checked rather than merely believed.

The GPU comes from the ECC_GPU environment variable, which is read when this
file is imported and baked into the function definitions.  That works on every
Modal version.  Newer clients also accept --gpu on the entry points, which
overrides it per call via Function.with_options.

The image builds a fat binary covering Ampere through Blackwell so any GPU type
works without a rebuild.  The autotuner rebuilds for the local architecture
only, which takes seconds rather than minutes.

Distinguished points and results land in a Modal Volume, so a search can be
stopped and resumed and several containers can contribute to one corpus.
"""

import hashlib
import json
import os
import pathlib
import re
import signal
import struct
import subprocess
import sys
import time

import modal
# Modal imports this module as /root/modal_app.py, separately from the source
# tree copied into the image. Resolve helpers from that tree on the container.
REMOTE = "/root/ecc2k130"
helperRoot = pathlib.Path(__file__).parent if modal.is_local() else pathlib.Path(REMOTE)
sys.path.insert(0, str(helperRoot))
from codegen.benchreport import benchResult, bestResult, parseRate, reportsVerified, summarizeSamples

CUDA_VERSION = os.environ.get("ECC_CUDA_VERSION", "12.8.1")
if not re.fullmatch(r'\d+\.\d+\.\d+', CUDA_VERSION):
    raise ValueError("ECC_CUDA_VERSION must be a toolkit version such as 13.0.2")

# Valid values include T4, L4, A10, L40S, A100, A100-80GB, RTX-PRO-6000, H100,
# H200, B200 and B300; append ":n" for several of them.
DEFAULT_GPU = os.environ.get("ECC_GPU", "H100")

# Compute capability per Modal GPU type.  sm_120 is the Blackwell workstation
# part (RTX PRO 6000), sm_100 is B200/B300, sm_90 is H100/H200, sm_89 is
# L40S/L4, sm_86 is A10, sm_80 is A100, sm_75 is T4.
GPU_ARCH = {
    "T4": "75", "L4": "89", "L40S": "89", "A10": "86", "A10G": "86",
    "A100": "80", "A100-80GB": "80", "H100": "90", "H200": "90",
    "B200": "100", "B300": "100", "RTX-PRO-6000": "120",
}
ALL_ARCHES = ("80", "89", "90", "100", "120")


def archesFor(gpu):
    """Architectures to bake into the image.

    Every architecture is a separate full ptxas pass over a kernel that takes
    minutes to compile, so a fat binary covering five of them costs about four
    minutes on every image rebuild -- and add_local_dir invalidates the layer on
    any source change, so that is every iteration.  Build only what the chosen
    GPU can run.  An unrecognised name falls back to the full set rather than
    guessing, since a missing architecture is a runtime failure, not a slow
    build."""
    arch = GPU_ARCH.get(gpu.split(":")[0].strip())
    return (arch,) if arch else ALL_ARCHES


BAKED_ARCHES = archesFor(DEFAULT_GPU)
GENCODE = " ".join(
    "-gencode arch=compute_%s,code=sm_%s" % (a, a) for a in BAKED_ARCHES
)
LOCAL = pathlib.Path(__file__).parent

# What the image already contains.  A request for exactly this on a matching
# architecture needs no rebuild, and the default ::bench is exactly this -- it
# was spending minutes recompiling a binary it already had, in silence.
BAKED = {"batch": 32, "threads": 128, "leaf": 0, "minBlocks": 2}

# Cleared the first time buildFor actually builds.  `make -B gpu` replaces
# ./ecc2k130 in place, so once anything has rebuilt, the image's baked binary is
# gone and a later request for BAKED has to build it again.  Without this a
# sweep that rebuilt for one point and then reached BAKED would time the
# previous point's binary and report it under BAKED's knobs -- the same trap the
# generator comment in buildFor warns about, one level up.
bakedIntact = [True]

image = (
    modal.Image.from_registry(
        f"nvidia/cuda:{CUDA_VERSION}-devel-ubuntu24.04", add_python="3.12"
    )
    .entrypoint([])
    # Containers re-import this module; preserve the settings that selected
    # their image and baked architecture rather than reverting to defaults.
    .env({"ECC_CUDA_VERSION": CUDA_VERSION, "ECC_GPU": DEFAULT_GPU})
    .apt_install("build-essential")
    .add_local_dir(
        LOCAL,
        remote_path=REMOTE,
        copy=True,
        ignore=["ecc2k130-cpu", "ecc2k130", "build/*", "__pycache__", "*.pyc"],
    )
    .run_commands(
        # Regenerate before building so the baked binary's leaf is known to be
        # BAKED["leaf"] rather than whatever headers the local checkout held.
        f"cd {REMOTE}/codegen && python3 gen.py --out ../generated "
        f"--leaf {BAKED['leaf']}",
        # x86-64-v3 keeps the host binary runnable on any Modal machine; the
        # GPU client picks its own word width on the device.
        f"cd {REMOTE} && make cpu MARCH=x86-64-v3",
        f'cd {REMOTE} && make gpu ARCH="{GENCODE}" BATCH={BAKED["batch"]} '
        f'THREADS={BAKED["threads"]} MINBLOCKS={BAKED["minBlocks"]}',
    )
)

# Nsight Compute lives in its own image.  It is about two gigabytes and only the
# profiler wants it, so putting it in the main image would slow every bench and
# search rebuild for a tool most runs never invoke.  It layers on top, so the
# baked binary and the source are already present.
profileImage = image.apt_install("nsight-compute")

volume = modal.Volume.from_name("ecc2k130", create_if_missing=True)
app = modal.App("ecc2k130")

HOUR = 60 * 60

def onGpu(fn, gpu):
    """Point a function at a GPU type.

    Modal clients from 0.72 on can retarget a call with Function.with_options;
    older ones cannot, so fall back to the ECC_GPU environment variable, which
    is baked in at import time and therefore always works."""
    if not gpu or gpu == DEFAULT_GPU:
        return fn
    if hasattr(fn, "with_options"):
        return fn.with_options(gpu=gpu)
    raise SystemExit(
        "This Modal client cannot change the GPU per call "
        "(Function.with_options was added in 0.72).\n"
        "Either upgrade with `pip install -U modal`, or select the GPU through "
        "the environment:\n\n"
        f"    ECC_GPU={gpu} modal run modal_app.py::<entrypoint>\n"
    )


# ---------------------------------------------------------------------------
def sh(cmd, cwd=REMOTE, timeout=None):
    r = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True, timeout=timeout)
    return r.returncode, r.stdout + r.stderr


def shStream(cmd, cwd=REMOTE, timeout=None, prefix=""):
    """Run a command, echoing each line as it arrives and returning them too.

    Modal forwards a container's stdout to the caller's terminal live, but only
    if something writes to it.  Capturing the whole of a multi-minute compile and
    printing it at the end looks identical to a hang, which is what running
    ::bench used to look like."""
    p = subprocess.Popen(cmd, shell=True, cwd=cwd, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, text=True, bufsize=1)
    lines = []
    try:
        for line in p.stdout:
            lines.append(line)
            print(prefix + line.rstrip(), flush=True)
        p.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        p.kill()
        return 1, "".join(lines) + "\ntimeout"
    return p.returncode, "".join(lines)


def computeCapability():
    rc, out = sh("nvidia-smi --query-gpu=compute_cap --format=csv,noheader")
    if rc != 0 or not out.strip():
        return "90"
    return out.strip().splitlines()[0].strip().replace(".", "")


def gpuName():
    rc, out = sh("nvidia-smi --query-gpu=name --format=csv,noheader")
    return out.strip().splitlines()[0].strip() if rc == 0 and out.strip() else "unknown"


def buildFor(batch, threads, leaf, arch=None, minBlocks=2,
             streamKarat=False, smemSpill=False, globalCg=False):
    """Rebuild the client for one architecture and one set of knobs."""
    arch = arch or computeCapability()
    # leaf 0 means "let the generator choose by register budget", which is what
    # gen.py does with no --leaf: the largest halving-chain size whose
    # straight-line multiply still fits 255 live values.  On a 255-register GPU
    # that is 66 words at m=131; it is emphatically not the right answer for the
    # 16-register host, so the CPU reference build and the device build do not
    # want the same leaf, and only a GPU settles which one this is.
    # Regenerate unconditionally, passing leaf 0 straight through: gen.py reads
    # 0 as "choose by register budget".  Skipping the generator instead would
    # compile whatever headers the container already holds, which during an
    # autotune sweep is the previous point's leaf -- so the sweep would score
    # that build twice and label one of them the register-budget choice.
    want = {"batch": batch, "threads": threads, "leaf": leaf, "minBlocks": minBlocks}
    if smemSpill and int(CUDA_VERSION.split('.')[0]) < 13:
        return False, "--smem-spill requires ECC_CUDA_VERSION=13.x.y (CUDA 13 or newer)"
    experimental = streamKarat or smemSpill or globalCg
    if want == BAKED and not experimental and arch in BAKED_ARCHES and bakedIntact[0]:
        print(f"sm_{arch}: batch {batch}, threads {threads}, leaf {leaf}, "
              f"minBlocks {minBlocks} is what the image already holds; not rebuilding",
              flush=True)
        return True, "baked into the image"
    print(f"building for sm_{arch}: batch {batch}, threads {threads}, leaf {leaf}, "
          f"minBlocks {minBlocks} -- ptxas takes a few minutes on this kernel",
          flush=True)
    bakedIntact[0] = False
    rc, out = sh(f"cd codegen && python3 gen.py --out ../generated --leaf {leaf}")
    if rc != 0:
        return False, out
    gencode = f"-gencode arch=compute_{arch},code=sm_{arch}"
    rc, out = shStream(
        f'make -B gpu ARCH="{gencode}" BATCH={batch} THREADS={threads} '
        f"MINBLOCKS={minBlocks} STREAM_KARAT={int(streamKarat)} "
        f"SMEM_SPILL={int(smemSpill)} GLOBAL_CG={int(globalCg)}",
        timeout=1800,
        prefix="  build| ",
    )
    return rc == 0, out


def benchmarkIdentity():
    """Identity travels with the result even though the image has no .git."""
    root = pathlib.Path(REMOTE)
    source = hashlib.sha256()
    paths = [root / 'Makefile', root / 'modal_app.py']
    for folder in ('include', 'src', 'generated', 'codegen'):
        paths += sorted(p for p in (root / folder).rglob('*')
                        if p.suffix in ('.h', '.cu', '.cpp', '.py'))
    for path in paths:
        source.update(str(path.relative_to(root)).encode() + b'\0' + path.read_bytes() + b'\0')
    rc, compiler = sh('nvcc --version')
    gpuRc, gpu = sh('nvidia-smi --query-gpu=name,uuid,driver_version,pstate,clocks.current.sm,'
                   'clocks.current.memory,power.limit,temperature.gpu --format=csv')
    leaf = re.search(r'LEAF = (\d+)', (root / 'generated/eccF131.h').read_text())
    return dict(sourceSha256=source.hexdigest(),
                binarySha256=hashlib.sha256((root / 'ecc2k130').read_bytes()).hexdigest(),
                actualLeaf=int(leaf.group(1)), compiler=compiler, compilerReturncode=rc,
                gpuState=gpu, gpuStateReturncode=gpuRc, cudaImageVersion=CUDA_VERSION)


def measureBench(steps, launches, workers, preferL1, repeats):
    if steps <= 0 or launches <= 0 or repeats <= 0 or workers < 0:
        raise ValueError('steps, launches and repeats must be positive; workers must be nonnegative')
    cmd = f'./ecc2k130 --curve 131 --bench --steps {steps} --launches {launches} --verify 0'
    if workers:
        cmd += f' --threads {workers}'
    if preferL1:
        cmd += ' --prefer-l1'
    samples = []
    for repeat in range(repeats):
        try:
            rc, out = sh(cmd, timeout=1800)
        except subprocess.TimeoutExpired as exc:
            rc = 124
            out = exc.stdout or ''
            if isinstance(out, bytes):
                out = out.decode(errors='replace')
            out += '\nbenchmark timed out'
        sample = benchResult(cmd, rc, out)
        samples.append(sample)
        print(f'  repeat {repeat + 1}/{repeats}: {sample["rate"]:.3f} M it/s '
              f'({"complete" if sample["valid"] else "FAILED"})', flush=True)
        if not sample['valid']:
            break
    return summarizeSamples(samples)


# ---------------------------------------------------------------------------
@app.function(image=image, gpu=DEFAULT_GPU, timeout=2 * HOUR, volumes={"/data": volume})
def runValidate(batch=32, threads=128, leaf=0, minBlocks=2,
                streamKarat=False, smemSpill=False, globalCg=False, preferL1=False):
    """Field arithmetic, orbit invariants, solver, and end-to-end discrete
    logarithms recovered on the GPU itself."""
    cc = computeCapability()
    out = ["device: " + gpuName(), "compute capability: " + cc, ""]
    # Validate the requested candidate on the actual architecture. buildFor
    # reuses the image only when its untouched binary matches the request.
    ok, log = buildFor(batch, threads, leaf, arch=cc, minBlocks=minBlocks,
                       streamKarat=streamKarat, smemSpill=smemSpill, globalCg=globalCg)
    if not ok:
        raise RuntimeError("\n".join(out + ["rebuild failed", log]))
    out.append(json.dumps(benchmarkIdentity(), indent=2))
    cacheFlag = ' --prefer-l1' if preferL1 else ''
    rc, t = sh("./ecc2k130-cpu --test")
    out.append(t.strip())
    if rc != 0:
        raise RuntimeError("\n".join(out) + "\nHOST VALIDATION FAILED")

    out.append("\n--- end-to-end on the GPU ---")
    ok = True
    # A small curve reaches a collision almost immediately, so keep the walk
    # count modest: a million walks would overrun the report buffer on the
    # first launch and throw most of the points away.  Curves 19 and 13 have no
    # normal basis, so they exercise the polynomial-basis backend that ECC2K-95
    # depends on; curve 41 is solved through both backends.
    for curve, instances, workers, steps in (("23", 4, 256, 8), ("19", 4, 256, 8),
                                             ("13", 2, 128, 4), ("41", 4, 2048, 32),
                                             ("41 --poly-basis", 4, 2048, 32)):
        for i in range(instances):
            rc, t = sh(
                f"./ecc2k130 --curve {curve} --instance {i} --threads {workers} "
                f"--steps {steps} --dp-cap 262144 --verify 4{cacheFlag}"
            )
            line = [l for l in t.splitlines() if "planted" in l or "MISMATCH" in l]
            good = rc == 0 and 'MISMATCH' not in t and any("yes" in l for l in line)
            ok = ok and good
            out.append(f"curve {curve} instance {i}: " + ("; ".join(line) if line else t.strip()[-200:]))
    out.append("GPU END TO END: " + ("all instances solved" if ok else "FAILED"))

    # ECC2K-95 itself: no collision in a short run, but the reports have to be
    # reproducible from their seeds, which is what the server depends on.
    #
    # dp weight 36 is theta = 1/139 at m=97, far looser than any real search
    # would use: the point is to make a short run report at all, so there is
    # something to recompute.  Size the run to that.  4096 threads for 4
    # launches walked 268M iterations to verify 16 reports and took ten minutes
    # -- most of a validate run -- while overrunning the report buffer, so the
    # gate people are told to run first was the slowest thing here.  512 threads
    # for 2 launches still produces about 120k reports from half a million
    # parallel walks, stays inside --dp-cap, and takes well under a minute.
    # Staying inside the cap matters beyond speed: it means a "dropped" count in
    # validate output is a real signal rather than the expected state.
    for curve, weight, workers in ((97, 36, 512), (131, 50, 128)):
        rc, t = sh(f"./ecc2k130 --curve {curve} --dp-weight {weight} --threads {workers} "
                   f"--steps 16 --launches 2 --dp-cap 262144 --verify 16{cacheFlag}", timeout=1800)
        out.append(f"\n--- GF(2^{curve}) reporting path ---")
        out.append(t.strip())
        passed = reportsVerified(rc, t)
        out.append('REPORT REPLAY: ' + ('PASS' if passed else 'FAILED'))
        ok = ok and passed
    if not ok:
        raise RuntimeError("\n".join(out))
    return "\n".join(out)


@app.function(image=image, gpu=DEFAULT_GPU, timeout=2 * HOUR)
def runBench(batch=32, threads=128, leaf=0, minBlocks=2, steps=64, launches=20,
             workers=0, rebuild=True, repeats=3, streamKarat=False,
             smemSpill=False, globalCg=False, preferL1=False):
    """Completed-run median on the challenge curve, with reproducible identity."""
    info = dict(gpu=gpuName(), cc=computeCapability(), batch=batch,
                threads=threads, leaf=leaf, minBlocks=minBlocks, workers=workers,
                steps=steps, launches=launches, repeats=repeats, streamKarat=streamKarat,
                smemSpill=smemSpill, globalCg=globalCg, preferL1=preferL1)
    want = dict(batch=batch, threads=threads, leaf=leaf, minBlocks=minBlocks)
    if not rebuild and (streamKarat or smemSpill or globalCg or not bakedIntact[0]
                        or want != BAKED or info['cc'] not in BAKED_ARCHES):
        raise ValueError('rebuild=False requires the untouched matching baked binary')
    if rebuild:
        ok, log = buildFor(batch, threads, leaf, minBlocks=minBlocks,
                           streamKarat=streamKarat, smemSpill=smemSpill, globalCg=globalCg)
        info['buildLog'] = log
        if not ok:
            return dict(info, valid=False, rate=0.0, error=log)
    info['identity'] = benchmarkIdentity()
    info.update(measureBench(steps, launches, workers, preferL1, repeats))
    return info


def autotuneConfigs(batches, threadCounts, leaves, minBlocksList, configs):
    """The (leaf, batch, threads, minBlocks) builds to measure.

    `configs` names them one per entry as leaf:batch:threads:minBlocks, which is
    what autolab.py emits.  Its Pareto front is a handful of points, and
    expanding the distinct values of those points back into a cross product
    measures two or three times as many builds as the front has rows.  Without
    it the four lists are crossed, as before."""
    if configs.strip():
        out = []
        for entry in configs.split(","):
            if not entry.strip():
                continue
            leaf, batch, threads, mb = (int(x) for x in entry.split(":"))
            out.append((leaf, batch, threads, mb))
        return out
    ints = lambda t: [int(x) for x in t.split(",") if x]
    return [(leaf, batch, threads, mb)
            for leaf in ints(leaves)
            for threads in ints(threadCounts)
            for batch in ints(batches)
            for mb in ints(minBlocksList)]


@app.function(image=image, gpu=DEFAULT_GPU, timeout=4 * HOUR, volumes={"/data": volume})
def runAutotune(batches="8,16,32,64", threadCounts="64,128,256", leaves="0,17,33,66",
                minBlocksList="2,4,8", steps=64, launches=12, configs="",
                workers=0, repeats=3, streamKarat=False, smemSpill=False,
                globalCg=False, preferL1=False):
    """Rank completed repetitions by median; retain failed candidates as errors."""
    results = []
    arch = computeCapability()
    name = gpuName()
    plan = autotuneConfigs(batches, threadCounts, leaves, minBlocksList, configs)
    print(f"{len(plan)} builds to measure on {name} (sm_{arch})")
    for leaf, batch, threads, mb in plan:
        t0 = time.time()
        ok, log = buildFor(batch, threads, leaf, arch, mb,
                           streamKarat=streamKarat, smemSpill=smemSpill, globalCg=globalCg)
        cfg = dict(batch=batch, threads=threads, leaf=leaf, minBlocks=mb,
                   workers=workers, repeats=repeats, steps=steps, launches=launches,
                   streamKarat=streamKarat, smemSpill=smemSpill,
                   globalCg=globalCg, preferL1=preferL1,
                   buildSeconds=round(time.time() - t0, 1), buildLog=log)
        if not ok:
            results.append(dict(cfg, valid=False, rate=0.0, error=log))
            continue
        cfg['identity'] = benchmarkIdentity()
        cfg.update(measureBench(steps, launches, workers, preferL1, repeats))
        results.append(cfg)
        print(f"leaf {leaf} threads {threads} batch {batch} "
              f"minBlocks {mb}: {cfg['rate']:.3f} M it/s median")
    results.sort(key=lambda r: -r.get("rate", 0.0))
    report = {"gpu": name, "cc": arch, "results": results,
              "best": bestResult(results)}
    os.makedirs("/data/autotune", exist_ok=True)
    # Each sweep is evidence, including experiments that lose. Keep the legacy
    # latest-result path too for existing consumers.
    stem = name.replace(' ', '_')
    path = f"/data/autotune/{stem}-{time.time_ns()}.json"
    report['reportPath'] = path
    for dest in (path, f"/data/autotune/{stem}.json"):
        with open(dest, "w") as fh:
            json.dump(report, fh, indent=2)
    volume.commit()
    return report


@app.function(image=image, timeout=4 * HOUR, volumes={"/data": volume})
def runAutolab(batches="4,8,16,32", threadCounts="64,128,256",
               minBlocksList="1,2,3,4", leaves="0,17,33,66,131", arch="", top=6):
    """Search the build space offline, on a CPU container, and return a shortlist.

    ptxas is deterministic and needs no device, so registers, spill traffic and
    the occupancy that follows from them are all measurable without renting a
    GPU -- and this image already carries nvcc and ptxas for the real build.
    What it cannot measure is time, so it produces the few configurations worth
    running through ::autotune rather than a verdict.

    The metric cache lives in the volume keyed by architecture, so a second run
    only compiles what the first one did not."""
    arch = arch or (BAKED_ARCHES[0] if len(BAKED_ARCHES) == 1 else "90")
    os.makedirs("/data/autolab", exist_ok=True)
    out = f"/data/autolab/sm{arch}.json"
    rc, log = sh(
        f"python3 autolab.py --compiler nvcc --cuda-path=/usr/local/cuda "
        f"--ptxas=/usr/local/cuda/bin/ptxas --arch={arch} --batch {batches} "
        f"--threads {threadCounts} --min-blocks {minBlocksList} --leaf {leaves} "
        f"--top {top} --out {out}",
        cwd=f"{REMOTE}/codegen",
        timeout=4 * HOUR - 600,
    )
    volume.commit()
    # The front is the product, so read it back rather than scraping the log.
    front, configs = [], ""
    if os.path.exists(out):
        with open(out) as fh:
            saved = json.load(fh)
        front, configs = saved.get("front", []), saved.get("configs", "")
    return {"arch": arch, "returncode": rc, "cache": out, "front": front,
            "configs": configs, "log": log}


@app.function(image=profileImage, gpu=DEFAULT_GPU, timeout=2 * HOUR)
def runProfile(batch=32, threads=128, leaf=0, minBlocks=2, steps=4, launches=1,
               section="", metrics="", workers=0, streamKarat=False,
               smemSpill=False, globalCg=False, preferL1=False):
    """Profile the walk kernel with Nsight Compute, or say precisely why not.

    Whether this works at all is a property of the host, not of this code.
    NVIDIA has restricted the GPU performance counters to administrators since
    driver 418.43, and a container gets at them only if the host loaded the
    driver with NVreg_RestrictProfilingToAdminUsers=0 or the container was given
    CAP_SYS_ADMIN.  Neither is ours to set on Modal.  So this runs ncu and, if
    the counters are refused, reports that as the answer rather than as an
    error -- the run is still informative, because it settles the question.

    `steps` and `launches` default low: ncu serialises and replays each kernel
    to collect counters, so a profiled launch is orders of magnitude slower than
    a real one.  Four steps is plenty to characterise a kernel whose every step
    is identical."""
    arch = computeCapability()
    name = gpuName()
    rc, ver = sh("ncu --version")
    if rc != 0:
        return {"gpu": name, "available": False,
                "why": "ncu is not on PATH in this image", "log": ver[-2000:]}

    ok, log = buildFor(batch, threads, leaf, arch, minBlocks,
                       streamKarat=streamKarat, smemSpill=smemSpill, globalCg=globalCg)
    if not ok:
        return {"gpu": name, "available": False, "why": "build failed",
                "log": log[-2000:]}

    what = ""
    if metrics:
        what = "--metrics %s" % metrics
    elif section:
        what = "--section %s" % section
    else:
        # Enough to answer where the walk kernel's time goes without asking for
        # the full set, which replays the kernel many more times.
        what = ("--section SpeedOfLight --section MemoryWorkloadAnalysis "
                "--section LaunchStats --section Occupancy "
                "--section WarpStateStats")
    runFlags = f" --threads {workers}" if workers else ""
    if preferL1:
        runFlags += " --prefer-l1"
    identity = benchmarkIdentity()
    rc, out = shStream(
        f"ncu --target-processes all --kernel-name eccWalkKernel "
        f"--launch-count 1 {what} "
        f"./ecc2k130 --curve 131 --bench --steps {steps} --launches {launches} "
        f"--verify 0{runFlags}",
        timeout=1 * HOUR,
        prefix="  ncu| ",
    )

    denied = ("ERR_NVGPU_DEBUG_PERF_COUNTER_ACCESS_DENIED" in out
              or "The user does not have permission" in out
              or "insufficient permissions" in out.lower())
    if denied:
        return {
            "gpu": name, "available": False,
            "why": ("the GPU performance counters are restricted to "
                    "administrators on this host, which is a driver and "
                    "container-capability setting Modal controls, not this app"),
            "remedy": ("ask Modal whether profiling can be enabled for this GPU "
                       "class.  Failing that, the multiplier's share can be "
                       "measured without any counters by adding K extra "
                       "multiplications by a runtime-supplied identity to each "
                       "step and fitting the slope of rate against K -- the "
                       "marginal cost of one multiply, straight off the card"),
            "log": out[-4000:],
        }
    return {"gpu": name, "cc": arch, "available": rc == 0,
            "batch": batch, "threads": threads, "leaf": leaf,
            "minBlocks": minBlocks, "workers": workers, "identity": identity,
            "streamKarat": streamKarat, "smemSpill": smemSpill,
            "globalCg": globalCg, "preferL1": preferL1,
            "report": out[-2000:] if rc else out}


# Expected rho iterations, and the weight cutoff that makes walks short enough
# that most of them actually report within the run.
CURVE_FACTS = {
    # curve: (field size m, log2 of expected iterations)
    131: (131, 60.9),
    97: (97, 44.0),
    83: (83, 37.1),
    41: (41, 16.6),
    23: (23, 8.1),
}


def recommendedWeight(curve, walks, plannedIters=0):
    """Pick the distinguished-point cutoff for a given amount of parallelism.

    A walk reports after about 1/theta steps, and each of `walks` parallel walks
    only gets budget/walks steps, so the cutoff has to satisfy 1/theta well below
    that or most walks never report at all.  Aim for a quarter of the budget.

    The budget is the whole expected run only when the run can actually finish.
    ECC2K-130 needs 2^60.9 iterations -- decades of GPU time -- so sizing its
    cutoff against that yields a weight no walk reaches in any session anyone
    will ever run: at the full-run choice a four-hour pass on a fast GPU reports
    about two hundred points in total.  Pass plannedIters and the cutoff is
    sized for the run being done instead, which stores more points per useful
    iteration but is the difference between collecting a corpus and collecting
    nothing.  Storing more is never wrong -- collision probability depends on
    iterations walked, not on how often walks report -- it only costs disk."""
    if curve not in CURVE_FACTS:
        return -1
    m, logIters = CURVE_FACTS[curve]
    budget = 2.0 ** logIters
    if plannedIters > 0:
        budget = min(budget, float(plannedIters))
    budget /= max(1.0, float(walks))
    target = max(64.0, budget / 4.0)
    total = 0
    for k in range(0, m + 1):
        c = 1
        for i in range(k):
            c = c * (m - i) // (i + 1)
        total += c
        p = float(total) / (2.0 ** m)
        if p > 0 and 1.0 / p <= target:
            return k
    return m


DP_RECORD = struct.Struct("<Q3Q")  # seed, then the canonical orbit hash


def corpusFiles(curve, root="/data/dp"):
    """Every distinguished-point file in the volume for one curve."""
    if not os.path.isdir(root):
        return []
    return [os.path.join(root, fn) for fn in sorted(os.listdir(root))
            if fn.startswith(f"curve{curve}-") and fn.endswith(".bin")]


def corpusCount(path):
    """Records in a corpus file.  Fixed-width, so this is a stat, not a scan."""
    if not os.path.exists(path):
        return 0
    return os.path.getsize(path) // DP_RECORD.size


# The client prints a progress line every couple of seconds:
#   "     12.0 s     842.135 M it/s   1234567 iterations    890 dp    889 stored"
PROGRESS_RE = re.compile(
    r"([\d.]+)\s+s\s+([\d.]+)\s+M it/s\s+(\d+)\s+iterations\s+(\d+)\s+dp\s+(\d+)\s+stored"
    r"(?:\s+(\d+)\s+dropped)?")


def parseProgress(line):
    """Pull the numbers out of one client progress line, or None."""
    m = PROGRESS_RE.search(line)
    if not m:
        return None
    return {"seconds": float(m.group(1)), "rate": float(m.group(2)),
            "iters": int(m.group(3)), "dp": int(m.group(4)), "stored": int(m.group(5)),
            "dropped": int(m.group(6)) if m.group(6) else 0}


def humanRate(r):
    """M it/s spans four orders of magnitude between a debug run on a CPU and a
    real one on a GPU, so let the small end keep its digits."""
    return ("%.3f" % r) if r < 10 else ("%.1f" % r)


def humanCount(n):
    """Counts here run to the trillions, where digits stop being readable."""
    for scale, suffix in ((1e12, "T"), (1e9, "G"), (1e6, "M"), (1e3, "k")):
        if n >= scale:
            return "%.2f%s" % (n / scale, suffix)
    return str(int(n))


def humanTime(sec):
    if sec < 0:
        return "0s"
    if sec < 3600:
        return "%dm%02ds" % (sec // 60, sec % 60)
    return "%dh%02dm" % (sec // 3600, (sec % 3600) // 60)


def humanBytes(n):
    for scale, suffix in ((1 << 30, "GB"), (1 << 20, "MB"), (1 << 10, "kB")):
        if n >= scale:
            return "%.1f %s" % (float(n) / scale, suffix)
    return "%d B" % n


@app.function(image=image, gpu=DEFAULT_GPU, timeout=24 * HOUR, volumes={"/data": volume})
def runSearch(hours=1.0, curve=97, batch=8, threads=128, leaf=0, dpWeight=-1,
              runId=1, steps=256, workers=0, rebuild=True, walksTarget=4000000,
              checkpointEvery=300, resume=True, loadMax=50000000):
    """Collect distinguished points into the volume until the time budget runs
    out.  Records are 32 bytes of (seed, canonical orbit hash); a collision is
    resolved by recomputing both walks from their seeds.

    Nothing here is throwaway.  The container dies at the deadline, but the run
    does not: the client checkpoints its live walks to the volume, reloads the
    corpus at startup so old points still collide with new ones, and is stopped
    with SIGTERM rather than killed so it writes a final checkpoint first.  At
    any instant about a quarter of a run's iterations sit in walks that have not
    yet reported, so a hard kill would throw away 25% of the work done.

    Resuming needs the same shape it saved: curve, runId, worker thread count
    and the build-time batch size all appear in the checkpoint header, and a
    mismatch makes the client start fresh rather than misread the file.  So pass
    the same batch/workers/walksTarget you passed the first time.

    curve=97 is ECC2K-95, which Harley's group solved in 1998 after about
    2.16e13 iterations; the published answer is baked into the generated header
    so a recovered logarithm can be checked against it."""
    if rebuild:
        ok, log = buildFor(batch, threads, leaf)
        if not ok:
            return {"error": log[-2000:]}
    name = gpuName()
    os.makedirs("/data/dp", exist_ok=True)
    os.makedirs("/data/ckpt", exist_ok=True)
    dpFile = f"/data/dp/curve{curve}-run{runId}.bin"
    ckFile = f"/data/ckpt/curve{curve}-run{runId}.ck"
    rc, out = sh(f"./ecc2k130 --curve {curve} --bench --steps 8 --launches 4 --verify 0")
    rate = parseRate(out) or 1.0
    # The number of reports comes out at roughly four times the number of
    # parallel walks, because each walk has only total/walks steps to spend and
    # the cutoff is set so it reports a few times within that.  So the walk
    # count, not the GPU, decides how much storage the run needs: cap it.
    workerThreads = workers
    if not workerThreads and walksTarget:
        perThread = batch * 32
        workerThreads = max(1024, int(walksTarget) // perThread)
    walks = workerThreads * batch * 32
    # The bench above measured this GPU, so the length of the pass about to run
    # is known rather than guessed; size the cutoff against that.
    plannedIters = rate * 1e6 * hours * HOUR
    if dpWeight < 0:
        dpWeight = recommendedWeight(curve, walks, plannedIters)
    print(f"{walks} parallel walks, distinguished-point weight {dpWeight}, "
          f"{plannedIters:.3g} iterations planned this pass")
    cmd = (f"./ecc2k130 --curve {curve} --steps {steps} --run-id {runId} "
           f"--dp-file {dpFile} --verify 4 --launches 0 --threads {workerThreads} "
           f"--checkpoint-every {int(checkpointEvery)} --load-max {int(loadMax)}")
    if resume:
        cmd += f" --checkpoint {ckFile}"
    if dpWeight >= 0:
        cmd += f" --dp-weight {dpWeight}"
    # Every other worker's corpus counts too: a collision between this run and a
    # sibling's is just as good as one within a single run, and finding it here
    # beats waiting for an offline merge.  The client's own dp file reloads
    # itself, so it is not named twice.  It reads newest-first and stops at
    # loadMax, because a collection run that never finishes -- ECC2K-130 is
    # decades of GPU time -- grows a corpus no container can hold in memory.
    for other in sorted(corpusFiles(curve)):
        if other != dpFile:
            cmd += f" --load {other}"
    deadline = time.time() + hours * HOUR
    print(f"{name}: collecting into {dpFile} for {hours} h at ~{rate:.2f} M it/s")
    proc = subprocess.Popen(cmd, shell=True, cwd=REMOTE, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True)
    lines = []
    solved = None
    started = time.time()
    lastCommit = started
    lastReport = started
    last = None
    stopped = ""
    # A container's output is the only window into a run that will outlive the
    # terminal that started it, so summarise on a fixed clock rather than
    # relaying the client's own line every two seconds.
    expected = 2.0 ** CURVE_FACTS[curve][1] if curve in CURVE_FACTS else 0.0
    print(f"progress every 60 s; corpus {dpFile}"
          + (f", checkpoint {ckFile}" if resume else ""), flush=True)
    try:
        while proc.poll() is None:
            line = proc.stdout.readline()
            if not line:
                break
            lines.append(line.rstrip())
            if "k = " in line:
                solved = line.strip()
            prog = parseProgress(line)
            if prog:
                last = prog
            # Anything that is not a progress line is an event -- a collision, a
            # verification failure, a checkpoint warning -- and is worth showing
            # as it happens rather than only in the tail at the end.
            elif line.strip():
                print("  " + line.rstrip(), flush=True)
            now = time.time()
            if now - lastReport >= 60:
                lastReport = now
                if last:
                    frac = (" (%.3f%% of 2^%.1f)" % (100.0 * last["iters"] / expected,
                                                     CURVE_FACTS[curve][1])) if expected else ""
                    drop = ("  %s DROPPED" % humanCount(last["dropped"])) if last["dropped"] else ""
                    print("[%s] %s M it/s  %s iters%s  %s dp  %s distinct  "
                          "corpus %s  %s left%s"
                          % (humanTime(now - started), humanRate(last["rate"]),
                             humanCount(last["iters"]), frac,
                             humanCount(last["dp"]), humanCount(last["stored"]),
                             humanBytes(os.path.getsize(dpFile) if os.path.exists(dpFile) else 0),
                             humanTime(deadline - now), drop), flush=True)
                else:
                    print("[%s] no progress line yet (still starting up?)"
                          % humanTime(now - started), flush=True)
            # Commit on a clock, not on a line count: the client's output rate
            # depends on the launch size, so counting lines would space the
            # commits arbitrarily far apart on a quiet run.
            if now - lastCommit > 60:
                volume.commit()
                lastCommit = now
            if now > deadline:
                stopped = "deadline"
                break
    finally:
        if proc.poll() is None:
            # SIGTERM, not kill: the client finishes the launch in flight, then
            # writes its checkpoint and flushes the corpus.  Draining stdout
            # while it does keeps the pipe from filling and wedging the exit.
            proc.send_signal(signal.SIGTERM)
            graceful = time.time() + 600
            while proc.poll() is None and time.time() < graceful:
                line = proc.stdout.readline()
                if not line:
                    break
                lines.append(line.rstrip())
            try:
                proc.wait(timeout=max(1, graceful - time.time()))
            except Exception:
                stopped = "killed before it could checkpoint"
                proc.kill()
                proc.wait(timeout=60)
    # Commit last, so the checkpoint the client just wrote is part of the
    # snapshot rather than the one before it.
    volume.commit()
    return {"gpu": name, "distinguishedPoints": corpusCount(dpFile), "file": dpFile,
            "checkpoint": ckFile if os.path.exists(ckFile) else None,
            "checkpointBytes": os.path.getsize(ckFile) if os.path.exists(ckFile) else 0,
            "iterations": last["iters"] if last else 0,
            "rate": last["rate"] if last else 0.0,
            "elapsed": round(time.time() - started, 1),
            "stopped": stopped, "solved": solved, "tail": lines[-25:]}


@app.function(image=image, timeout=2 * HOUR, volumes={"/data": volume})
def mergeCorpus(curve=131):
    """Merge every distinguished-point file in the volume and report duplicate
    hashes, which are the candidate collisions.

    Records are the client's 32-byte binary format, so a partial trailing record
    (a container that died mid-write) is ignored rather than misparsed."""
    seen = {}
    dup = []
    total = 0
    short = 0
    files = corpusFiles(curve)
    if not files:
        return {"error": "no distinguished points yet"}
    for path in files:
        with open(path, "rb") as fh:
            while True:
                rec = fh.read(DP_RECORD.size)
                if len(rec) < DP_RECORD.size:
                    short += len(rec)
                    break
                seed, h0, h1, h2 = DP_RECORD.unpack(rec)
                total += 1
                h = (h0, h1, h2)
                if h in seen and seen[h] != seed:
                    dup.append(("%016x" % seen[h], "%016x" % seed,
                                "%016x%016x%016x" % (h2, h1, h0)))
                else:
                    seen[h] = seed
    return {"files": len(files), "records": total, "distinct": len(seen),
            "truncatedBytes": short, "collisions": dup[:50],
            "collisionCount": len(dup)}


@app.function(image=image, timeout=4 * HOUR, volumes={"/data": volume})
def solveCorpus(curve=131, loadMax=0):
    """Recover the logarithm from a collision spanning several corpora.

    mergeCorpus reports that two seeds reached the same orbit; it does not say
    what k is, and a campaign that ends with collisionCount 1 and no answer has
    not finished.  The client already knows how: --load rebuilds the store from
    the corpora, and a collision found during that reload is rewalked from both
    seeds and solved on the spot.  Nothing here needs a GPU -- rewalking two
    walks is cheap -- so this runs the host build and the whole thing costs a
    container without a card.

    --load-max defaults to off here, which is the opposite of the search, and
    the reason is that the two runs want opposite things.  The cap keeps a long
    collection run inside memory; this wants the whole corpus resident at once,
    because the colliding pair is exactly what a cap might drop.  Which risk
    matters depends on whether the run can finish: a complete ECC2K-95 corpus
    is about 17M points, near 1 GB resident, so uncapped is free.  An ECC2K-130
    collection run does not finish and grows by roughly that much every pass,
    so there uncapped eventually means out of memory -- hence the parameter,
    and hence reporting the size before asking for it rather than after."""
    files = corpusFiles(curve)
    if not files:
        return {"error": "no distinguished points yet"}
    records = 0
    for f in files:
        records += corpusCount(f)
    resident = records * 64          # key, entry and table overhead, roughly
    cmd = ("./ecc2k130-cpu --curve %d --threads 1 --steps 1 --launches 1 "
           "--verify 0 --run-id 65535" % curve)
    if loadMax:
        cmd += " --load-max %d" % int(loadMax)
    for f in files:
        cmd += " --load %s" % f
    print("solving from %d file(s), %s points, roughly %s resident%s"
          % (len(files), humanCount(records), humanBytes(resident),
             "" if not loadMax else " (capped at %s)" % humanCount(loadMax)),
          flush=True)
    rc, out = sh(cmd, timeout=4 * HOUR - 60)
    k = None
    verified = False
    matches = None
    for line in out.splitlines():
        t = line.strip()
        if t.startswith("k = "):
            k = t[4:].strip()
        if "verified [k]P == Q" in t:
            verified = True
        if "matches the published solution:" in t:
            matches = t.split(":")[-1].strip()
    return {"files": len(files), "records": records,
            "residentBytes": resident, "k": k, "verified": verified,
            "matchesPublished": matches, "tail": out.strip()[-1500:]}


# ---------------------------------------------------------------------------
@app.function(image=image, timeout=1 * HOUR)
def runCompileCheck(arch="120", streamKarat=False, smemSpill=False, globalCg=False):
    """Compile experimental device code on a CPU container, without renting a GPU."""
    if not re.fullmatch(r'\d+', arch):
        raise ValueError('arch must be a numeric compute capability')
    ok, log = buildFor(32, 128, 0, arch=arch, streamKarat=streamKarat,
                       smemSpill=smemSpill, globalCg=globalCg)
    if not ok:
        raise RuntimeError(log)
    return dict(arch=arch, cudaImageVersion=CUDA_VERSION, streamKarat=streamKarat,
                smemSpill=smemSpill, globalCg=globalCg, buildLog=log,
                binarySha256=hashlib.sha256(pathlib.Path(REMOTE, 'ecc2k130').read_bytes()).hexdigest())


@app.local_entrypoint()
def compile_check(arch: str = "120", stream_karat: bool = False,
                  smem_spill: bool = False, global_cg: bool = False):
    print(json.dumps(runCompileCheck.remote(arch=arch, streamKarat=stream_karat,
                                           smemSpill=smem_spill, globalCg=global_cg), indent=2))


@app.local_entrypoint()
def validate(gpu: str = "", batch: int = 32, threads: int = 128, leaf: int = 0,
             min_blocks: int = 2, stream_karat: bool = False, smem_spill: bool = False,
             global_cg: bool = False, prefer_l1: bool = False):
    print(onGpu(runValidate, gpu).remote(batch=batch, threads=threads, leaf=leaf,
          minBlocks=min_blocks, streamKarat=stream_karat, smemSpill=smem_spill,
          globalCg=global_cg, preferL1=prefer_l1))


@app.local_entrypoint()
def bench(gpu: str = "", batch: int = 32, threads: int = 128, leaf: int = 0,
          min_blocks: int = 2, steps: int = 64, launches: int = 20,
          workers: int = 0, repeats: int = 3, stream_karat: bool = False,
          smem_spill: bool = False, global_cg: bool = False, prefer_l1: bool = False):
    r = onGpu(runBench, gpu).remote(batch=batch, threads=threads, leaf=leaf,
                                    minBlocks=min_blocks, steps=steps, launches=launches,
                                    workers=workers, repeats=repeats, streamKarat=stream_karat,
                                    smemSpill=smem_spill, globalCg=global_cg, preferL1=prefer_l1)
    print(json.dumps(r, indent=2))
    if not r.get('valid'):
        raise RuntimeError('benchmark did not complete successfully')


@app.local_entrypoint()
def autotune(gpu: str = "", batches: str = "8,16,32,64",
             thread_counts: str = "64,128,256", leaves: str = "0,17,33,66",
             min_blocks_list: str = "2,4,8", configs: str = "",
             steps: int = 64, launches: int = 12, workers: int = 0, repeats: int = 3,
             stream_karat: bool = False, smem_spill: bool = False,
             global_cg: bool = False, prefer_l1: bool = False):
    """--configs takes leaf:batch:threads:minBlocks entries, as ::autolab
    prints them, and measures exactly those instead of a cross product."""
    r = onGpu(runAutotune, gpu).remote(batches=batches, threadCounts=thread_counts,
                                       leaves=leaves, minBlocksList=min_blocks_list,
                                       configs=configs, steps=steps, launches=launches,
                                       workers=workers, repeats=repeats, streamKarat=stream_karat,
                                       smemSpill=smem_spill, globalCg=global_cg, preferL1=prefer_l1)
    print(json.dumps(r, indent=2))
    if r['best'] is None:
        raise RuntimeError('no benchmark candidate completed successfully')


@app.local_entrypoint()
def autolab(batches: str = "4,8,16,32", thread_counts: str = "64,128,256",
            min_blocks_list: str = "1,2,3,4", leaves: str = "0,17,33,66,131",
            arch: str = "", top: int = 6):
    """Offline build-space search.  No GPU is rented; ptxas does not need one."""
    r = runAutolab.remote(batches=batches, threadCounts=thread_counts,
                          minBlocksList=min_blocks_list, leaves=leaves,
                          arch=arch, top=top)
    print(r["log"])
    print("cache: %s (in the ecc2k130 volume)" % r["cache"])


@app.local_entrypoint()
def campaign(gpu: str = "", batches: str = "4,8,16,32",
             thread_counts: str = "64,128,256", min_blocks_list: str = "1,2,3,4",
             leaves: str = "0,17,33,66,131", top: int = 6, arch: str = ""):
    """Offline search, then measure its Pareto front on the card.

    The two halves are the point: the CPU container maps the whole space for
    the price of compile time, and the GPU only ever runs the handful of builds
    that survived.  Nothing here decides a winner offline -- the ranking that
    comes out of the search is a proxy, and the rates that come out of the
    measurement are the result."""
    r = runAutolab.remote(batches=batches, threadCounts=thread_counts,
                          minBlocksList=min_blocks_list, leaves=leaves,
                          arch=arch, top=top)
    print(r["log"])
    if not r["configs"]:
        print("the offline search produced no front; not renting a GPU")
        return
    print("measuring %d builds on a GPU: %s\n" % (len(r["front"]), r["configs"]))
    m = onGpu(runAutotune, gpu).remote(configs=r["configs"])
    print(json.dumps(m, indent=2))


@app.local_entrypoint()
def profile(gpu: str = "", batch: int = 32, threads: int = 128, leaf: int = 0,
            min_blocks: int = 2, steps: int = 4, section: str = "",
            metrics: str = "", workers: int = 0, stream_karat: bool = False,
            smem_spill: bool = False, global_cg: bool = False, prefer_l1: bool = False):
    """Nsight Compute on the walk kernel.

    Whether the counters are readable is a host setting rather than anything
    this app controls, so a refusal is itself the answer worth having."""
    r = onGpu(runProfile, gpu).remote(batch=batch, threads=threads, leaf=leaf,
                                      minBlocks=min_blocks, steps=steps,
                                      section=section, metrics=metrics, workers=workers,
                                      streamKarat=stream_karat, smemSpill=smem_spill,
                                      globalCg=global_cg, preferL1=prefer_l1)
    if not r.get("available"):
        print("Nsight Compute did not run: %s" % r.get("why"))
        if r.get("remedy"):
            print("  %s" % r["remedy"])
        print(r.get("log", "")[-2000:])
        return
    print(json.dumps({k: v for k, v in r.items() if k != 'report'}, indent=2))
    print(r["report"])


@app.local_entrypoint()
def search(gpu: str = "", hours: float = 1.0, curve: int = 97, batch: int = 8,
           threads: int = 128, leaf: int = 0, dp_weight: int = -1, run_id: int = 1,
           walks: int = 4000000, load_max: int = 50000000):
    r = onGpu(runSearch, gpu).remote(hours=hours, curve=curve, batch=batch,
                                     threads=threads, leaf=leaf, dpWeight=dp_weight,
                                     runId=run_id, walksTarget=walks, loadMax=load_max)
    print(json.dumps(r, indent=2))


@app.local_entrypoint()
def fanout(gpu: str = "", count: int = 4, hours: float = 1.0, curve: int = 97,
           batch: int = 8, threads: int = 128, leaf: int = 0, dp_weight: int = -1,
           walks: int = 4000000, load_max: int = 50000000):
    """Run `count` independent searchers, each with its own run id so their
    seeds never collide, then merge what they produced."""
    fn = onGpu(runSearch, gpu)
    calls = [fn.spawn(hours=hours, curve=curve, batch=batch, threads=threads, leaf=leaf,
                      dpWeight=dp_weight, runId=i + 1, walksTarget=walks, loadMax=load_max)
             for i in range(count)]
    for c in calls:
        print(json.dumps(c.get(), indent=2))
    print(json.dumps(mergeCorpus.remote(curve=curve), indent=2))


@app.local_entrypoint()
def merge(curve: int = 131, solve: bool = True, load_max: int = 0):
    """Scan the corpora for collisions, and recover k when one is there."""
    r = mergeCorpus.remote(curve=curve)
    print(json.dumps(r, indent=2))
    if solve and r.get("collisionCount"):
        print("\n%d collision(s); recovering the logarithm" % r["collisionCount"])
        print(json.dumps(solveCorpus.remote(curve=curve, loadMax=load_max), indent=2))

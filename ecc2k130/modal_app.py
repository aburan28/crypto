"""Run the ECC2K-130 client on Modal GPUs.

    ECC_GPU=RTX-PRO-6000 modal run modal_app.py::validate
    ECC_GPU=RTX-PRO-6000 modal run modal_app.py::bench
    ECC_GPU=RTX-PRO-6000 modal run modal_app.py::autotune
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

import json
import os
import pathlib
import re
import signal
import struct
import subprocess
import time

import modal

CUDA_VERSION = "12.8.1"
# sm_120 is the Blackwell workstation part (RTX PRO 6000), sm_100 is B200,
# sm_90 is H100/H200, sm_89 is L40S, sm_80 is A100.  Every architecture adds a
# full ptxas pass to the image build, so trim this list if build time matters.
GENCODE = " ".join(
    "-gencode arch=compute_%s,code=sm_%s" % (a, a) for a in ("80", "89", "90", "100", "120")
)
REMOTE = "/root/ecc2k130"
LOCAL = pathlib.Path(__file__).parent

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
        ignore=["ecc2k130-cpu", "ecc2k130", "build/*", "__pycache__", "*.pyc"],
    )
    .run_commands(
        # x86-64-v3 keeps the host binary runnable on any Modal machine; the
        # GPU client picks its own word width on the device.
        f"cd {REMOTE} && make cpu MARCH=x86-64-v3",
        f'cd {REMOTE} && make gpu ARCH="{GENCODE}" BATCH=32 THREADS=128 MINBLOCKS=2',
    )
)

volume = modal.Volume.from_name("ecc2k130", create_if_missing=True)
app = modal.App("ecc2k130")

HOUR = 60 * 60

# Valid values include T4, L4, A10, L40S, A100, A100-80GB, RTX-PRO-6000, H100,
# H200, B200 and B300; append ":n" for several of them.
DEFAULT_GPU = os.environ.get("ECC_GPU", "H100")


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


def computeCapability():
    rc, out = sh("nvidia-smi --query-gpu=compute_cap --format=csv,noheader")
    if rc != 0 or not out.strip():
        return "90"
    return out.strip().splitlines()[0].strip().replace(".", "")


def gpuName():
    rc, out = sh("nvidia-smi --query-gpu=name --format=csv,noheader")
    return out.strip().splitlines()[0].strip() if rc == 0 and out.strip() else "unknown"


def buildFor(batch, threads, leaf, arch=None, minBlocks=2):
    """Rebuild the client for one architecture and one set of knobs."""
    arch = arch or computeCapability()
    # leaf 0 means "let the generator choose by register budget", which is what
    # gen.py does with no --leaf: the largest halving-chain size whose
    # straight-line multiply still fits 255 live values.  On a 255-register GPU
    # that is 66 words at m=131; it is emphatically not the right answer for the
    # 16-register host, so the CPU reference build and the device build do not
    # want the same leaf, and only a GPU settles which one this is.
    if leaf:
        rc, out = sh(f"cd codegen && python3 gen.py --out ../generated --leaf {leaf}")
        if rc != 0:
            return False, out
    gencode = f"-gencode arch=compute_{arch},code=sm_{arch}"
    rc, out = sh(
        f'make -B gpu ARCH="{gencode}" BATCH={batch} THREADS={threads} '
        f"MINBLOCKS={minBlocks}",
        timeout=1800,
    )
    return rc == 0, out


def parseRate(text):
    """Iterations per second, in millions, from the client's own report.  The
    "finished" line carries the average over the whole run, so prefer it."""
    best = 0.0
    for line in text.splitlines():
        if "M it/s" not in line:
            continue
        head = line.split("M it/s")[0].split()
        if not head:
            continue
        try:
            rate = float(head[-1])
        except ValueError:
            continue
        if line.strip().startswith("finished"):
            return rate
        best = max(best, rate)
    return best


# ---------------------------------------------------------------------------
@app.function(image=image, gpu=DEFAULT_GPU, timeout=2 * HOUR, volumes={"/data": volume})
def runValidate():
    """Field arithmetic, orbit invariants, solver, and end-to-end discrete
    logarithms recovered on the GPU itself."""
    out = ["device: " + gpuName(), "compute capability: " + computeCapability(), ""]
    rc, t = sh("./ecc2k130-cpu --test")
    out.append(t.strip())
    if rc != 0:
        return "\n".join(out) + "\nHOST VALIDATION FAILED"

    out.append("\n--- end-to-end on the GPU ---")
    ok = True
    # A small curve reaches a collision almost immediately, so keep the walk
    # count modest: a million walks would overrun the report buffer on the
    # first launch and throw most of the points away.  Curves 19 and 13 have no
    # normal basis, so they exercise the polynomial-basis backend that ECC2K-95
    # depends on; curve 41 is solved through both backends.
    for curve, instances, threads, steps in (("23", 4, 256, 8), ("19", 4, 256, 8),
                                             ("13", 2, 128, 4), ("41", 4, 2048, 32),
                                             ("41 --poly-basis", 4, 2048, 32)):
        for i in range(instances):
            rc, t = sh(
                f"./ecc2k130 --curve {curve} --instance {i} --threads {threads} "
                f"--steps {steps} --dp-cap 262144 --verify 4"
            )
            line = [l for l in t.splitlines() if "planted" in l or "MISMATCH" in l]
            good = any("yes" in l for l in line)
            ok = ok and good
            out.append(f"curve {curve} instance {i}: " + ("; ".join(line) if line else t.strip()[-200:]))
    out.append("GPU END TO END: " + ("all instances solved" if ok else "FAILED"))

    # ECC2K-95 itself: no collision in a short run, but the reports have to be
    # reproducible from their seeds, which is what the server depends on.
    rc, t = sh("./ecc2k130 --curve 97 --dp-weight 36 --threads 4096 --steps 16 "
               "--launches 4 --dp-cap 262144 --verify 16", timeout=1800)
    out.append("\n--- ECC2K-95 reporting path ---")
    out.append("\n".join(t.strip().splitlines()[-3:]))
    if "MISMATCH" in t:
        out.append("ECC2K-95 REPORTS DID NOT REPRODUCE")
        ok = False
    return "\n".join(out)


@app.function(image=image, gpu=DEFAULT_GPU, timeout=1 * HOUR)
def runBench(batch=32, threads=128, leaf=0, minBlocks=2, steps=64, launches=20,
             workers=0, rebuild=True):
    """Throughput on the challenge curve."""
    info = {"gpu": gpuName(), "cc": computeCapability(), "batch": batch,
            "threads": threads, "leaf": leaf, "minBlocks": minBlocks}
    if rebuild:
        ok, log = buildFor(batch, threads, leaf, minBlocks=minBlocks)
        if not ok:
            info["error"] = log[-2000:]
            return info
    cmd = f"./ecc2k130 --curve 131 --bench --steps {steps} --launches {launches} --verify 0"
    if workers:
        cmd += f" --threads {workers}"
    rc, out = sh(cmd, timeout=1800)
    info["rate"] = parseRate(out)
    info["raw"] = out.strip()[-1200:]
    return info


@app.function(image=image, gpu=DEFAULT_GPU, timeout=4 * HOUR, volumes={"/data": volume})
def runAutotune(batches="8,16,32,64", threadCounts="64,128,256", leaves="0,17,33,66",
                minBlocksList="2,4,8", steps=64, launches=12):
    """Sweep the build-time knobs on the real device and report the best.

    minBlocks is the interesting one: asking ptxas for more resident blocks per
    SM trades registers for occupancy.  Offline, 2 is free (255 registers, no
    extra spills) while 8 cuts registers to 64 and nearly doubles spill traffic,
    so which side wins is exactly what this measures."""
    results = []
    arch = computeCapability()
    name = gpuName()
    for leaf in [int(x) for x in leaves.split(",") if x]:
        for threads in [int(x) for x in threadCounts.split(",") if x]:
            for batch in [int(x) for x in batches.split(",") if x]:
                for mb in [int(x) for x in minBlocksList.split(",") if x]:
                    t0 = time.time()
                    ok, log = buildFor(batch, threads, leaf, arch, mb)
                    cfg = {"batch": batch, "threads": threads, "leaf": leaf,
                           "minBlocks": mb}
                    if not ok:
                        results.append(dict(cfg, rate=0.0, error=log[-400:]))
                        continue
                    rc, out = sh(
                        f"./ecc2k130 --curve 131 --bench --steps {steps} "
                        f"--launches {launches} --verify 0",
                        timeout=1800,
                    )
                    rate = parseRate(out)
                    results.append(dict(cfg, rate=rate,
                                        buildSeconds=round(time.time() - t0, 1)))
                    print(f"leaf {leaf} threads {threads} batch {batch} "
                          f"minBlocks {mb}: {rate:.3f} M it/s")
    results.sort(key=lambda r: -r.get("rate", 0.0))
    report = {"gpu": name, "cc": arch, "results": results, "best": results[0] if results else None}
    os.makedirs("/data/autotune", exist_ok=True)
    path = f"/data/autotune/{name.replace(' ', '_')}.json"
    with open(path, "w") as fh:
        json.dump(report, fh, indent=2)
    volume.commit()
    return report


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
    r"([\d.]+)\s+s\s+([\d.]+)\s+M it/s\s+(\d+)\s+iterations\s+(\d+)\s+dp\s+(\d+)\s+stored")


def parseProgress(line):
    """Pull the numbers out of one client progress line, or None."""
    m = PROGRESS_RE.search(line)
    if not m:
        return None
    return {"seconds": float(m.group(1)), "rate": float(m.group(2)),
            "iters": int(m.group(3)), "dp": int(m.group(4)), "stored": int(m.group(5))}


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
                    print("[%s] %s M it/s  %s iters%s  %s dp  %s distinct  "
                          "corpus %s  %s left"
                          % (humanTime(now - started), humanRate(last["rate"]),
                             humanCount(last["iters"]), frac,
                             humanCount(last["dp"]), humanCount(last["stored"]),
                             humanBytes(os.path.getsize(dpFile) if os.path.exists(dpFile) else 0),
                             humanTime(deadline - now)), flush=True)
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


# ---------------------------------------------------------------------------
@app.local_entrypoint()
def validate(gpu: str = ""):
    print(onGpu(runValidate, gpu).remote())


@app.local_entrypoint()
def bench(gpu: str = "", batch: int = 32, threads: int = 128, leaf: int = 0,
          minBlocks: int = 2, steps: int = 64, launches: int = 20):
    r = onGpu(runBench, gpu).remote(batch=batch, threads=threads, leaf=leaf,
                                    minBlocks=minBlocks, steps=steps, launches=launches)
    print(json.dumps({k: v for k, v in r.items() if k != "raw"}, indent=2))
    if "raw" in r:
        print(r["raw"])


@app.local_entrypoint()
def autotune(gpu: str = "", batches: str = "8,16,32,64",
             threadCounts: str = "64,128,256", leaves: str = "0,17,33,66",
             minBlocksList: str = "2,4,8"):
    r = onGpu(runAutotune, gpu).remote(batches=batches, threadCounts=threadCounts,
                                       leaves=leaves, minBlocksList=minBlocksList)
    print(json.dumps(r, indent=2))


@app.local_entrypoint()
def search(gpu: str = "", hours: float = 1.0, curve: int = 97, batch: int = 8,
           threads: int = 128, leaf: int = 0, dpWeight: int = -1, runId: int = 1,
           walks: int = 4000000, loadMax: int = 50000000):
    r = onGpu(runSearch, gpu).remote(hours=hours, curve=curve, batch=batch,
                                     threads=threads, leaf=leaf, dpWeight=dpWeight,
                                     runId=runId, walksTarget=walks, loadMax=loadMax)
    print(json.dumps(r, indent=2))


@app.local_entrypoint()
def fanout(gpu: str = "", count: int = 4, hours: float = 1.0, curve: int = 97,
           batch: int = 8, threads: int = 128, leaf: int = 0, dpWeight: int = -1,
           walks: int = 4000000, loadMax: int = 50000000):
    """Run `count` independent searchers, each with its own run id so their
    seeds never collide, then merge what they produced."""
    fn = onGpu(runSearch, gpu)
    calls = [fn.spawn(hours=hours, curve=curve, batch=batch, threads=threads, leaf=leaf,
                      dpWeight=dpWeight, runId=i + 1, walksTarget=walks, loadMax=loadMax)
             for i in range(count)]
    for c in calls:
        print(json.dumps(c.get(), indent=2))
    print(json.dumps(mergeCorpus.remote(curve=curve), indent=2))


@app.local_entrypoint()
def merge(curve: int = 131):
    print(json.dumps(mergeCorpus.remote(curve=curve), indent=2))

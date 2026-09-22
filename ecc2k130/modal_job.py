"""Run one benchmark job script on a Modal GPU and bring the results back.

    modal run modal_job.py --job benchmarks/two-chains/gpujob.sh --out /tmp/two-chains
    modal run modal_job.py --job benchmarks/two-chains/gpujob.sh --out /tmp/x --gpu RTX-PRO-6000 --env REPS=3

The job script is any script in this tree written for the container contract
that benchmarks/*/gpujob.sh use: the tree is at /work, results go to
$RESULTS (/results), nvcc 13.3 and nvidia-smi are on PATH.  It runs inside
nvidia/cuda:13.3.1-devel-ubuntu24.04 on the requested GPU (default the RTX
PRO 6000, the campaign's SKU), and everything the job wrote under /results
comes back as results.tgz, unpacked into --out beside a launch.json receipt.

Compared with modal_app.py this image bakes no binary: the job builds what it
measures, so the same script gives the same receipts here, on RunPod
(runpod_job.py) and on EC2 (aws/bench_job.py).  Credentials: MODAL_TOKEN_ID
and MODAL_TOKEN_SECRET in the environment, or `modal token set`.
"""
import io
import json
import os
import pathlib
import subprocess
import tarfile
import time

import modal

LOCAL = pathlib.Path(__file__).parent
REMOTE = "/work"
CUDA_VERSION = os.environ.get("ECC_CUDA_VERSION", "13.3.1")
DEFAULT_GPU = os.environ.get("ECC_GPU", "RTX-PRO-6000")
HOUR = 60 * 60

image = (
    modal.Image.from_registry(f"nvidia/cuda:{CUDA_VERSION}-devel-ubuntu24.04", add_python="3.12")
    .entrypoint([])
    .apt_install("build-essential", "make", "g++", "python3")
    .add_local_dir(
        LOCAL, remote_path=REMOTE, copy=True,
        ignore=["ecc2k130-cpu", "ecc2k130-cpu-v3", "ecc2k130-cpu-v4", "ecc2k130", "ecc2k130-*",
                "build/*", "__pycache__", "*.pyc", "*.bin", "*.tgz"],
    )
)
app = modal.App("ecc2k130-job")


@app.function(image=image, gpu=DEFAULT_GPU, timeout=3 * HOUR)
def run(job: str, env: dict) -> dict:
    """Run /work/<job> with the environment, streaming its output; return the results tarball."""
    results = "/results"
    os.makedirs(results, exist_ok=True)
    full = dict(os.environ, RESULTS=results, **env)
    started = time.time()
    with open(os.path.join(results, "job.log"), "w") as log:
        p = subprocess.Popen(["bash", os.path.join(REMOTE, job)], cwd=REMOTE, env=full,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
        for line in p.stdout:
            print(line.rstrip(), flush=True)
            log.write(line)
        rc = p.wait()
    with open(os.path.join(results, "exit-code"), "w") as f:
        f.write("%d\n" % rc)
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as archive:
        archive.add(results, arcname=".")
    gpu = subprocess.run("nvidia-smi --query-gpu=name,driver_version --format=csv,noheader",
                         shell=True, capture_output=True, text=True).stdout.strip()
    return {"exitCode": rc, "seconds": time.time() - started, "gpu": gpu, "results": buf.getvalue()}


@app.local_entrypoint()
def main(job: str = "benchmarks/two-chains/gpujob.sh", out: str = "/tmp/ecc2k130-job",
         gpu: str = "", env: str = ""):
    """--env takes comma-separated KEY=VALUE pairs passed to the job (e.g. REPS=3)."""
    if not (LOCAL / job).is_file():
        raise SystemExit("job script not found: %s" % (LOCAL / job))
    rev = subprocess.run(["git", "rev-parse", "HEAD"], cwd=LOCAL, capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain", "--", "."], cwd=LOCAL,
                                capture_output=True, text=True).stdout.strip())
    environment = {"SOURCE_REV": rev + (" (uncommitted changes)" if dirty else "")}
    for item in filter(None, env.split(",")):
        key, _, value = item.partition("=")
        environment[key] = value
    fn = run if not gpu or gpu == DEFAULT_GPU else run.with_options(gpu=gpu)
    outdir = pathlib.Path(out)
    outdir.mkdir(parents=True, exist_ok=False)
    receipt = {"job": job, "gpu": gpu or DEFAULT_GPU, "cudaImage": f"nvidia/cuda:{CUDA_VERSION}-devel-ubuntu24.04",
               "gitRev": rev, "gitDirty": dirty, "env": environment, "provider": "modal",
               "startedAt": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    result = fn.remote(job, environment)
    (outdir / "results.tgz").write_bytes(result["results"])
    with tarfile.open(outdir / "results.tgz") as archive:
        archive.extractall(outdir / "results", filter="data")
    receipt.update(exitCode=result["exitCode"], seconds=result["seconds"], deviceLine=result["gpu"],
                   finishedAt=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    (outdir / "launch.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print("job exit code", result["exitCode"], "on", result["gpu"], "- results in", outdir / "results")

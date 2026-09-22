"""Run one benchmark job script on a Modal GPU and bring the results back.

    modal run --detach modal_job.py --job benchmarks/two-chains/gpujob.sh --out /tmp/two-chains
    modal run --detach modal_job.py --job benchmarks/two-chains/gpujob.sh --out /tmp/x --gpu B200 --env REPS=3
    modal run modal_job.py --fetch <token> --out /tmp/x      # collect a run whose client died

The job script is any script in this tree written for the container contract
that benchmarks/*/gpujob.sh use: the tree is at /work, results go to
$RESULTS (/results), nvcc 13.3 and nvidia-smi are on PATH.  It runs inside
nvidia/cuda:13.3.1-devel-ubuntu24.04 on the requested GPU (default the RTX
PRO 6000, the campaign's SKU), and everything the job wrote under /results
is tarred into the `ecc2k130-jobs` Volume under the run's token, then
downloaded and unpacked into --out beside a launch.json receipt.

Results go through the Volume rather than the function's return value so that
a client that loses its connection mid-run (a gRPC heartbeat timeout stopped
an ephemeral app 35 minutes into a B200 run) loses nothing: with --detach the
app keeps running, and `--fetch <token>` collects the results afterwards.

Compared with modal_app.py this image bakes no binary: the job builds what it
measures, so the same script gives the same receipts here, on RunPod
(runpod_job.py) and on EC2 (aws/bench_job.py).  Credentials: MODAL_TOKEN_ID
and MODAL_TOKEN_SECRET in the environment, or `modal token set`.
"""
import json
import os
import pathlib
import subprocess
import tarfile
import time
import uuid

import modal

LOCAL = pathlib.Path(__file__).parent
REMOTE = "/work"
DATA = "/data"
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
volume = modal.Volume.from_name("ecc2k130-jobs", create_if_missing=True)
app = modal.App("ecc2k130-job")


# The verification phase of a job re-walks device reports and files a million
# distinguished points on the host, so the container gets a few cores.
@app.function(image=image, gpu=DEFAULT_GPU, cpu=8, timeout=3 * HOUR, volumes={DATA: volume})
def run(job: str, env: dict, token: str) -> dict:
    """Run /work/<job> with the environment; results tarball to the Volume under <token>."""
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
    gpu = subprocess.run("nvidia-smi --query-gpu=name,driver_version --format=csv,noheader",
                         shell=True, capture_output=True, text=True).stdout.strip()
    meta = {"exitCode": rc, "seconds": time.time() - started, "gpu": gpu, "token": token}
    out = os.path.join(DATA, token)
    os.makedirs(out, exist_ok=True)
    with tarfile.open(os.path.join(out, "results.tgz"), "w:gz") as archive:
        archive.add(results, arcname=".")
    with open(os.path.join(out, "done.json"), "w") as f:
        json.dump(meta, f)
    volume.commit()
    return meta


def fetch(token: str, outdir: pathlib.Path):
    """Download <token>/results.tgz from the Volume into outdir/results; return done.json or None."""
    try:
        meta = json.loads(b"".join(volume.read_file(f"{token}/done.json")).decode())
    except Exception:
        return None
    data = b"".join(volume.read_file(f"{token}/results.tgz"))
    (outdir / "results.tgz").write_bytes(data)
    with tarfile.open(outdir / "results.tgz") as archive:
        archive.extractall(outdir / "results", filter="data")
    return meta


@app.local_entrypoint()
def main(job: str = "benchmarks/two-chains/gpujob.sh", out: str = "/tmp/ecc2k130-job",
         gpu: str = "", env: str = "", fetch_token: str = ""):
    """--env takes comma-separated KEY=VALUE pairs passed to the job (e.g. REPS=3)."""
    outdir = pathlib.Path(out)
    if fetch_token:
        outdir.mkdir(parents=True, exist_ok=True)
        meta = fetch(fetch_token, outdir)
        if meta is None:
            raise SystemExit("no finished run under token %s yet" % fetch_token)
        receipt = json.loads((outdir / "launch.json").read_text()) if (outdir / "launch.json").exists() else {"token": fetch_token}
        receipt.update(exitCode=meta["exitCode"], seconds=meta["seconds"], deviceLine=meta["gpu"],
                       fetchedAt=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
        (outdir / "launch.json").write_text(json.dumps(receipt, indent=2) + "\n")
        print("job exit code", meta["exitCode"], "on", meta["gpu"], "- results in", outdir / "results")
        if meta["exitCode"] != 0:
            raise SystemExit(1)
        return
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
    outdir.mkdir(parents=True, exist_ok=False)
    token = uuid.uuid4().hex
    receipt = {"job": job, "gpu": gpu or DEFAULT_GPU, "cudaImage": f"nvidia/cuda:{CUDA_VERSION}-devel-ubuntu24.04",
               "gitRev": rev, "gitDirty": dirty, "env": environment, "provider": "modal", "token": token,
               "volume": "ecc2k130-jobs", "startedAt": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    (outdir / "launch.json").write_text(json.dumps(receipt, indent=2) + "\n")
    call = fn.spawn(job, environment, token)
    receipt["functionCallId"] = call.object_id
    (outdir / "launch.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print("spawned", call.object_id, "token", token, "- if this client dies, collect with:",
          "modal run modal_job.py --fetch-token", token, "--out", out, flush=True)
    while True:
        try:
            meta = call.get(timeout=60)
            break
        except TimeoutError:
            continue
    fetched = fetch(token, outdir)
    meta = fetched or meta
    receipt.update(exitCode=meta["exitCode"], seconds=meta["seconds"], deviceLine=meta["gpu"],
                   finishedAt=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    (outdir / "launch.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print("job exit code", meta["exitCode"], "on", meta["gpu"], "- results in", outdir / "results")
    if meta["exitCode"] != 0:
        raise SystemExit(1)

"""Run KernelAgent on the ECC2K-130 walk step on a Modal GPU, then time the result.

    export FIREWORKS_API_KEY=...          # read here, handed to Modal as a secret
    modal run --detach modal_kernelagent.py --out /tmp/ka-walk
    modal run --detach modal_kernelagent.py --out /tmp/ka-walk --model accounts/fireworks/models/... --gpu B200
    modal run modal_kernelagent.py --fetch-token <token> --out /tmp/ka-walk   # collect a run whose client died

The container clones the KernelAgent fork at a pinned commit, runs
run_kernelagent.py (Fireworks model, strict exact-match test) and then
bench_candidate.py on the kernel it produced.  Everything under /results
(launch.json, kernelagent.json, kernel.py, test_walk_step.py, bench.json,
job.log, KernelAgent's own run directories) is tarred into the
`ecc2k130-jobs` Volume under the run's token and unpacked into --out.

The key never reaches a file: it is read from FIREWORKS_API_KEY in the
launching shell and passed with modal.Secret.from_local_environ.  Modal
credentials: MODAL_TOKEN_ID and MODAL_TOKEN_SECRET, or `modal token set`.
"""
import json
import os
import pathlib
import shutil
import subprocess
import tarfile
import time
import uuid

import modal

LOCAL = pathlib.Path(__file__).parent
REMOTE = "/work/kernelagent-harness"
KA_DIR = "/opt/kernelagent"
KA_REPO = "https://github.com/aburan28/KernelAgent"
KA_COMMIT = os.environ.get("KERNELAGENT_COMMIT", "23a1fb4ebd611e2faeb7effc87eb5b23732039cf")
DATA = "/data"
DEFAULT_GPU = os.environ.get("ECC_GPU", "RTX-PRO-6000")
HOUR = 60 * 60

image = (
    modal.Image.from_registry("nvidia/cuda:12.8.1-devel-ubuntu24.04", add_python="3.12")
    .entrypoint([])
    .apt_install("git", "build-essential")
    .pip_install("torch==2.8.0", index_url="https://download.pytorch.org/whl/cu128")
    .run_commands(
        f"git clone {KA_REPO} {KA_DIR} && git -C {KA_DIR} checkout {KA_COMMIT}",
        f"pip install -e {KA_DIR}",
    )
    .add_local_dir(LOCAL, remote_path=REMOTE, copy=True, ignore=["__pycache__", "*.pyc"])
)
volume = modal.Volume.from_name("ecc2k130-jobs", create_if_missing=True)
app = modal.App("ecc2k130-kernelagent")


def sh(cmd, log, env, cwd):
    p = subprocess.Popen(cmd, cwd=cwd, env=env, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, text=True, bufsize=1)
    for line in p.stdout:
        print(line.rstrip(), flush=True)
        log.write(line)
    return p.wait()


@app.function(image=image, gpu=DEFAULT_GPU, cpu=8, timeout=6 * HOUR, volumes={DATA: volume},
              secrets=[modal.Secret.from_local_environ(["FIREWORKS_API_KEY"])])
def run(model: str, workers: int, rounds: int, token: str) -> dict:
    results = "/results"
    os.makedirs(results, exist_ok=True)
    env = dict(os.environ, PYTHONPATH=KA_DIR)
    started = time.time()
    with open(os.path.join(results, "job.log"), "w") as log:
        cmd = ["python3", f"{REMOTE}/run_kernelagent.py", "--out", results,
               "--workers", str(workers), "--rounds", str(rounds)]
        if model:
            cmd += ["--model", model]
        rc = sh(cmd, log, env, cwd=results)
        bench_rc = None
        if os.path.exists(os.path.join(results, "kernel.py")):
            bench_rc = sh(["python3", f"{REMOTE}/bench_candidate.py",
                           "--kernel", os.path.join(results, "kernel.py"), "--out", results],
                          log, env, cwd=results)
    for name in (".fuse", ".optimize"):
        src = os.path.join(results, name)
        if os.path.isdir(src):
            shutil.move(src, os.path.join(results, name.lstrip(".") + "-runs"))
    gpu = subprocess.run("nvidia-smi --query-gpu=name,driver_version --format=csv,noheader",
                         shell=True, capture_output=True, text=True).stdout.strip()
    meta = {"exitCode": rc, "benchExitCode": bench_rc, "seconds": time.time() - started,
            "gpu": gpu, "token": token, "kernelagentCommit": KA_COMMIT}
    with open(os.path.join(results, "meta.json"), "w") as f:
        json.dump(meta, f, indent=2)
    out = os.path.join(DATA, token)
    os.makedirs(out, exist_ok=True)
    with tarfile.open(os.path.join(out, "results.tgz"), "w:gz") as archive:
        archive.add(results, arcname=".")
    volume.commit()
    return meta


def fetch(token: str, out: str) -> None:
    dest = pathlib.Path(out)
    dest.mkdir(parents=True, exist_ok=True)
    tgz = dest / "results.tgz"
    with open(tgz, "wb") as f:
        for chunk in volume.read_file(f"{token}/results.tgz"):
            f.write(chunk)
    with tarfile.open(tgz) as archive:
        archive.extractall(dest, filter="data")
    tgz.unlink()
    print("results in", dest)


@app.local_entrypoint()
def main(out: str, model: str = "", gpu: str = "", workers: int = 4, rounds: int = 10,
         fetch_token: str = ""):
    if fetch_token:
        fetch(fetch_token, out)
        return
    if not os.environ.get("FIREWORKS_API_KEY"):
        raise SystemExit("set FIREWORKS_API_KEY in this shell first")
    token = time.strftime("%Y%m%d-%H%M%S-") + uuid.uuid4().hex[:8]
    pathlib.Path(out).mkdir(parents=True, exist_ok=True)
    pathlib.Path(out, "launch.json").write_text(json.dumps(
        {"token": token, "gpu": gpu or DEFAULT_GPU, "model": model or "(auto: GLM 5.3 Flash)",
         "workers": workers, "rounds": rounds, "kernelagentCommit": KA_COMMIT}, indent=2))
    print("token", token, "- if this client dies: modal run modal_kernelagent.py --fetch-token",
          token, "--out", out)
    fn = run.with_options(gpu=gpu) if gpu else run
    meta = fn.remote(model, workers, rounds, token)
    print(json.dumps(meta, indent=2))
    fetch(token, out)

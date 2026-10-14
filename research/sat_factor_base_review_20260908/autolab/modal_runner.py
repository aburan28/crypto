"""Modal runner for the IC boundary autolab (independent-host replay).

Runs boundary_autolab.py beats on a clean x86_64 Linux container with Rust
toolchain, git-clones the crypto repo at a pinned commit, pins every producer
to one CPU via taskset, and uploads the run directory to a Modal Volume.
"""
import modal

REPO = "https://github.com/aburan28/crypto.git"
COMMIT = "e3090cbc9"  # main: autolab comparator fix (includes all merged PRs)

image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("git", "build-essential", "util-linux", "curl", "ca-certificates", "pkg-config", "libssl-dev")
    .run_commands(
        "curl -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain stable --profile minimal",
    )
    .env({"PATH": "/root/.cargo/bin:$PATH"})
    .run_commands(f"git clone {REPO} /crypto && cd /crypto && git checkout {COMMIT}")
    .pip_install()
)

app = modal.App("ic-boundary-autolab")
vol = modal.Volume.from_name("ic-autolab-runs", create_if_missing=True)


@app.function(
    image=image,
    cpu=8,
    memory=16 * 1024,
    timeout=24 * 60 * 60,
    volumes={"/runs": vol},
)
def run_beat(beat, mode="single", cpu_pin=2, extra=None):
    import os, subprocess, tarfile, datetime

    os.chdir("/crypto")
    git_head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    print(f"repo head: {git_head}")
    nproc = int(subprocess.run(["nproc"], capture_output=True, text=True).stdout.strip())
    print(f"nproc: {nproc}")

    cmd = [
        "python3", "research/sat_factor_base_review_20260908/autolab/boundary_autolab.py",
        "launch-single" if mode == "single" else "launch",
        "--beat", beat,
    ]
    if mode == "single" and cpu_pin is not None and cpu_pin < nproc:
        cmd += ["--cpu", str(cpu_pin)]
    if extra:
        cmd += extra
    print("cmd:", " ".join(cmd))
    r = subprocess.run(cmd)
    print(f"autolab exit: {r.returncode}")

    # locate newest run dir and tar it into the volume
    runs_root = "research/sat_factor_base_review_20260908/autolab/runs"
    import pathlib
    run_dirs = sorted(pathlib.Path(runs_root).glob("2026*"), key=lambda p: p.stat().st_mtime)
    if not run_dirs:
        print("no run dirs found")
        return {"exit": r.returncode, "run_id": None}
    run_dir = run_dirs[-1]
    stamp = datetime.datetime.now().strftime("%Y%m%dT%H%M%SZ")
    tar_path = f"/runs/{stamp}-{run_dir.name}.tar.gz"
    with tarfile.open(tar_path, "w:gz") as tf:
        tf.add(run_dir, arcname=run_dir.name)
    print(f"archived to volume: {tar_path}")
    return {"exit": r.returncode, "run_id": run_dir.name, "archive": tar_path}


@app.local_entrypoint()
def main(beat: str = "koblitz.compact_orbit.n61_single_target", mode: str = "single", cpu: int = 2):
    print(run_beat.remote(beat=beat, mode=mode, cpu_pin=cpu))

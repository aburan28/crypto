"""Run the rho walk leads harness on Modal CPU containers and collect its JSON lines.

    modal run modal_run.py --plan full --out runs/<tag>

Launcher only: it builds harness/ (Rust, no dependencies) in the image, fans
the plan's shards out over containers and writes each command's stdout,
unchanged, to <out>/<exp>-<bits>.jsonl beside launch.json. Every number comes
from the Rust binary. Credentials: MODAL_TOKEN_ID / MODAL_TOKEN_SECRET.
"""
import json
import pathlib
import subprocess
import time

import modal

HERE = pathlib.Path(__file__).parent
CPU = 16

image = (
    modal.Image.from_registry("rust:1.88-slim-bookworm", add_python="3.12")
    .add_local_dir(HERE / "harness", remote_path="/h", copy=True, ignore=["target"])
    .run_commands("cd /h && cargo build --release && cp target/release/rho_walk_leads /usr/local/bin/")
)
app = modal.App("rho-walk-leads")


@app.function(image=image, cpu=CPU, memory=16384, timeout=3 * 60 * 60)
def run(argv: list) -> dict:
    t0 = time.time()
    host = subprocess.run(["sh", "-c", "nproc; grep -m1 'model name' /proc/cpuinfo"], capture_output=True, text=True).stdout
    out = subprocess.run(["rho_walk_leads", *argv], capture_output=True, text=True, check=True).stdout
    return {"argv": argv, "stdout": out, "wall_s": time.time() - t0, "host": host.strip()}


def shards(cmd, bits, trials, per, extra=()):
    exe = "e1" if cmd.startswith("e1") and cmd != "e1dp" else cmd
    return [(f"{cmd}-{bits}", [exe, "--bits", str(bits), "--t0", str(t0), "--trials", str(min(per, trials - t0)), *extra])
            for t0 in range(0, trials, per)]


PLANS = {
    "smoke": lambda: shards("e1", 24, 64, 32),
    "full": lambda: (
        shards("e1", 32, 2000, 250)
        + shards("e1", 36, 2000, 125)
        + shards("e1", 40, 500, 32)
        + shards("e1shift", 36, 2000, 125, ("--escape", "shift"))
        + [("e1noesc-32", ["e1", "--bits", "32", "--trials", "200", "--escape", "off"])]
        + shards("e1dp", 36, 1000, 125)
        + shards("rf", 32, 4000, 1000)
        + shards("rf", 36, 2000, 250)
        + shards("e2", 36, 1000, 125)
        + shards("e4", 36, 1000, 125)
        + [("bench", ["bench", "--walks", "1024", "--steps", "4000", "--reps", "21"])]
    ),
}


@app.local_entrypoint()
def main(plan: str = "smoke", out: str = "runs/smoke"):
    jobs = PLANS[plan]()
    od = HERE / out
    od.mkdir(parents=True, exist_ok=True)
    started = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, cwd=HERE).stdout.strip()
    results = list(run.map([a for _, a in jobs], order_outputs=True))
    files = {}
    for (name, _), r in zip(jobs, results):
        files.setdefault(name, []).append(r)
        with open(od / f"{name}.jsonl", "a") as f:
            f.write(r["stdout"])
    receipt = {
        "plan": plan, "started_utc": started, "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "commit": commit, "image": "rust:1.88-slim-bookworm", "cpu_per_container": CPU,
        "jobs": [{"file": n, "argv": r["argv"], "wall_s": round(r["wall_s"], 1), "host": r["host"]} for (n, _), r in zip(jobs, results)],
    }
    (od / "launch.json").write_text(json.dumps(receipt, indent=1))
    print(f"{len(jobs)} jobs -> {od}")

#!/bin/sh
# Offline Linux arm64 build; capped preparation, no factor-base import/search.
set -eu
python3 - "$@" <<'PY'
from decimal import Decimal
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import uuid

if len(sys.argv) != 3:
    raise SystemExit("usage: build-linux.sh /absolute/fresh/build-dir /absolute/authorization.json")
checkout = Path(subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip()).resolve(strict=True)
output = Path(sys.argv[1])
auth_path = Path(sys.argv[2]).resolve(strict=True)
image = "sha256:7c4ae649a84014c467d79319bbf17ce2632ae8b8be123ac2fb2ea5be46823f31"
if not output.is_absolute() or output.resolve().is_relative_to(checkout) or output.exists():
    raise SystemExit("build output must be a fresh absolute path outside the checkout")
output = output.resolve()
auth_bytes = auth_path.read_bytes()
auth = json.loads(auth_bytes)
if (auth.get("schema") != "n83.additional-local-budget/v1"
        or auth.get("study") != "koblitz_n83_factor_base_sweep_20261008"
        or auth.get("purpose") != "structured_orbit_import_and_one_m5_trial"
        or auth.get("prior_pilot_charged_seconds") != "7200"
        or auth.get("build_dir") != str(output.resolve())
        or not isinstance(auth.get("user_approval_reference"), str)
        or not auth["user_approval_reference"].strip()
        or not isinstance(auth.get("additional_authorized_seconds"), str)
        or not auth["additional_authorized_seconds"].isdecimal()
        or Decimal(auth["additional_authorized_seconds"]) < 1800):
    raise SystemExit("build requires a run-specific 1800-second additional authorization")
if subprocess.check_output(["git", "-C", str(checkout), "status", "--porcelain", "--untracked-files=normal"]):
    raise SystemExit("build requires a clean committed checkout")
commit = subprocess.check_output(["git", "-C", str(checkout), "rev-parse", "HEAD"], text=True).strip()
registry = Path(os.environ.get("CARGO_HOME", str(Path.home() / ".cargo"))) / "registry"
if not registry.is_dir():
    raise SystemExit("offline Cargo registry is unavailable")
inspected = subprocess.run(["docker", "image", "inspect", image,
                            "--format", "{{.Id}} {{.Architecture}} {{.Os}}"],
                           capture_output=True, text=True, timeout=15, check=True)
if inspected.stdout.strip() != f"{image} arm64 linux":
    raise SystemExit("pinned Linux arm64 build image is unavailable")

started = time.monotonic()
output.mkdir(parents=True, exist_ok=False)
target = output / "target"
target.mkdir()
name = f"n83-orbit-build-{uuid.uuid4().hex[:16]}"
command = [
    "docker", "run", "--pull", "never", "--platform", "linux/arm64",
    "--name", name, "--network", "none", "--cpus", "4",
    "--memory", "8192m", "--memory-swap", "8192m", "--pids-limit", "256",
    "--cap-drop", "ALL", "--security-opt", "no-new-privileges",
    "--user", f"{os.getuid()}:{os.getgid()}", "--tmpfs", "/cargo:rw,mode=1777",
    "--mount", f"type=bind,src={registry},dst=/cargo/registry,readonly",
    "--mount", f"type=bind,src={checkout},dst=/src,readonly",
    "--mount", f"type=bind,src={target},dst=/target",
    "--env", "CARGO_HOME=/cargo", "--env", "CARGO_TARGET_DIR=/target",
    "--workdir", "/src", image, "cargo", "build", "--release", "--offline",
    "--locked", "--no-default-features", "--example", "koblitz_n83_structured_orbit_worker",
]
timed_out = False
launch_error = None
cleanup_error = None
exit_code = None
with (output / "stdout.log").open("xb") as stdout, (output / "stderr.log").open("xb") as stderr:
    try:
        child = subprocess.Popen(command, stdout=stdout, stderr=stderr, start_new_session=True)
        try:
            exit_code = child.wait(timeout=450)
        except subprocess.TimeoutExpired:
            timed_out = True
            try:
                killed = subprocess.run(["docker", "kill", name], capture_output=True,
                                        text=True, timeout=15, check=False)
                if killed.returncode != 0:
                    cleanup_error = killed.stderr.strip() or "docker kill failed"
            except (OSError, subprocess.SubprocessError) as error:
                cleanup_error = f"docker kill failed: {error}"
            try:
                exit_code = child.wait(timeout=15)
            except subprocess.TimeoutExpired:
                child.kill()
                exit_code = child.wait()
    except (OSError, subprocess.SubprocessError) as error:
        launch_error = f"{type(error).__name__}: {error}"

def docker(*args):
    try:
        return subprocess.run(["docker", *args], capture_output=True, text=True,
                              timeout=15, check=False)
    except (OSError, subprocess.SubprocessError):
        return None

state_call = docker("inspect", "--format", "{{json .State}}", name)
state = None
if state_call and state_call.returncode == 0:
    try:
        state = json.loads(state_call.stdout)
    except json.JSONDecodeError:
        pass
removed = docker("rm", "-f", name) if state else None
remaining = docker("ps", "-aq", "--filter", f"name=^{name}$")
cleanup_ok = bool(remaining and remaining.returncode == 0 and not remaining.stdout.strip()
                  and (removed is None or removed.returncode == 0) and not cleanup_error)
source_clean = not subprocess.check_output(["git", "-C", str(checkout), "status", "--porcelain", "--untracked-files=normal"])
source_same = subprocess.check_output(["git", "-C", str(checkout), "rev-parse", "HEAD"], text=True).strip() == commit
binary = output / "worker"
artifact = target / "release" / "examples" / "koblitz_n83_structured_orbit_worker"
status = "PRODUCER_FAILURE_build"
if not cleanup_ok:
    status = "PRODUCER_FAILURE_cleanup"
elif timed_out:
    status = "UNKNOWN_build_wall_cap"
elif not source_clean or not source_same:
    status = "PRODUCER_FAILURE_source_changed"
elif (launch_error or exit_code != 0 or not isinstance(state, dict)
      or state.get("ExitCode") != 0 or state.get("OOMKilled") is not False):
    status = "PRODUCER_FAILURE_build"
elif not artifact.is_file():
    status = "PRODUCER_FAILURE_binary"
else:
    with artifact.open("rb") as stream:
        magic = stream.read(4)
    if magic == b"\x7fELF":
        shutil.copyfile(artifact, binary)
        binary.chmod(0o755)
        status = "PASS_build_only"
    else:
        status = "PRODUCER_FAILURE_binary"
elapsed = time.monotonic() - started
if elapsed > 600:
    status = "PRODUCER_FAILURE_budget_overshoot"
receipt = {
    "schema": "n83.structured-orbit-linux-build/v1", "status": status,
    "source_commit": commit,
    "binary_sha256": sha256(binary.read_bytes()).hexdigest() if binary.is_file() else None,
    "authorization_sha256": sha256(auth_bytes).hexdigest(),
    "container_image_id": image, "binary": str(binary) if binary.is_file() else None,
    "command": command, "container_state": state, "container_cleanup_verified": cleanup_ok,
    "exit_code": exit_code, "timed_out": timed_out, "launch_error": launch_error,
    "cleanup_error": cleanup_error, "source_clean": source_clean, "source_same": source_same,
    "elapsed_seconds_L0": elapsed, "conservative_charge_seconds": 600,
    "solver_search_executed": False,
}
with (output / "build.json").open("x") as stream:
    json.dump(receipt, stream, sort_keys=True, indent=2)
    stream.write("\n")
print(json.dumps({"status": status, "build_receipt": str(output / "build.json")}))
raise SystemExit(0 if status == "PASS_build_only" else 2)
PY

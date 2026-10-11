"""One bounded, source-matched offline Linux build for the N83 comparison."""

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import subprocess
import time


IMAGE = "sha256:7c4ae649a84014c467d79319bbf17ce2632ae8b8be123ac2fb2ea5be46823f31"
NAME = "n83-ordered-pair-build-20261011"
WALL_CAP = 300
EXAMPLE = "koblitz_n83_factor_base_export"


def checked_output(command: list[str]) -> str:
    return subprocess.check_output(command, text=True, timeout=15).strip()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    checkout = Path(checked_output(["git", "rev-parse", "--show-toplevel"])).resolve(strict=True)
    target = args.target.resolve(strict=True)
    output = args.output.resolve()
    if not target.is_dir() or not output.is_absolute() or output.exists():
        raise SystemExit("target must exist; output must be a fresh absolute directory")
    if output.is_relative_to(checkout) or target.is_relative_to(checkout):
        raise SystemExit("build artifacts must be outside the clean checkout")
    if checked_output(["git", "status", "--porcelain"]):
        raise SystemExit("build requires a clean committed checkout")
    commit = checked_output(["git", "rev-parse", "HEAD"])
    registry = Path.home() / ".cargo" / "registry"
    if not registry.is_dir():
        raise SystemExit("offline Cargo registry missing")
    inspected = checked_output(["docker", "image", "inspect", IMAGE,
                                "--format", "{{.Id}} {{.Architecture}} {{.Os}}"])
    if inspected != f"{IMAGE} arm64 linux":
        raise SystemExit("pinned Linux image missing or changed")
    if checked_output(["docker", "ps", "-aq", "--filter", f"name=^{NAME}$"]):
        raise SystemExit("named build container already exists")

    output.mkdir(parents=True)
    command = [
        "docker", "run", "--rm", "--pull", "never", "--platform", "linux/arm64",
        "--name", NAME, "--network", "none", "--cpus", "4", "--memory", "8192m",
        "--memory-swap", "8192m", "--pids-limit", "256", "--cap-drop", "ALL",
        "--security-opt", "no-new-privileges", "--user", f"{os.getuid()}:{os.getgid()}",
        "--tmpfs", "/cargo:rw,mode=1777", "--tmpfs", "/tmp:rw,mode=1777",
        "--mount", f"type=bind,src={registry},dst=/cargo/registry,readonly",
        "--mount", f"type=bind,src={checkout},dst=/src,readonly",
        "--mount", f"type=bind,src={target},dst=/target",
        "--env", "CARGO_HOME=/cargo", "--env", "CARGO_TARGET_DIR=/target",
        "--workdir", "/src", IMAGE,
        "cargo", "build", "--release", "--offline", "--locked",
        "--no-default-features", "--example", EXAMPLE,
    ]
    started = time.monotonic()
    timed_out = False
    exit_code = None
    launch_error = None
    with (output / "stdout.log").open("xb") as stdout, (output / "stderr.log").open("xb") as stderr:
        try:
            child = subprocess.Popen(command, stdout=stdout, stderr=stderr, start_new_session=True)
            try:
                exit_code = child.wait(timeout=WALL_CAP)
            except subprocess.TimeoutExpired:
                timed_out = True
                try:
                    subprocess.run(["docker", "kill", NAME], capture_output=True, timeout=15, check=False)
                except (OSError, subprocess.SubprocessError):
                    pass
                try:
                    exit_code = child.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    child.kill()
                    exit_code = child.wait()
        except (OSError, subprocess.SubprocessError) as error:
            launch_error = f"{type(error).__name__}: {error}"
    cleanup_error = None
    try:
        subprocess.run(["docker", "rm", "-f", NAME], capture_output=True, timeout=15, check=False)
        remaining = checked_output(["docker", "ps", "-aq", "--filter", f"name=^{NAME}$"])
        cleanup_ok = not remaining
    except (OSError, subprocess.SubprocessError) as error:
        cleanup_ok = False
        cleanup_error = f"{type(error).__name__}: {error}"
    source_clean = not checked_output(["git", "status", "--porcelain"])
    source_same = checked_output(["git", "rev-parse", "HEAD"]) == commit
    binary = target / "release" / "examples" / EXAMPLE
    binary_bytes = binary.read_bytes() if binary.is_file() else None
    status = "PASS_build_only"
    if not cleanup_ok:
        status = "PRODUCER_FAILURE_cleanup"
    elif timed_out:
        status = "UNKNOWN_build_wall_cap"
    elif launch_error or exit_code != 0 or binary_bytes is None or binary_bytes[:4] != b"\x7fELF":
        status = "PRODUCER_FAILURE_build"
    elif not source_clean or not source_same:
        status = "PRODUCER_FAILURE_source_changed"
    retained_binary = output / "worker"
    if status == "PASS_build_only":
        shutil.copy2(binary, retained_binary)
        retained_binary.chmod(0o755)
    receipt = {
        "schema": "n83.ordered-first-pair-linux-build/v1", "status": status,
        "source_commit": commit, "command": command, "image_id": IMAGE,
        "binary": str(retained_binary) if status == "PASS_build_only" else None,
        "binary_sha256": sha256(binary_bytes).hexdigest() if binary_bytes else None,
        "binary_bytes": len(binary_bytes) if binary_bytes else None,
        "exit_code": exit_code, "timed_out": timed_out, "launch_error": launch_error,
        "cleanup_verified": cleanup_ok, "cleanup_error": cleanup_error,
        "source_clean": source_clean, "source_same": source_same,
        "wall_cap_seconds": WALL_CAP, "wall_elapsed_seconds_L0": time.monotonic() - started,
        "conservative_charge_seconds": WALL_CAP, "solver_search_executed": False,
    }
    (output / "build.json").write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"status": status, "receipt": str(output / "build.json")}, sort_keys=True))
    return 0 if status == "PASS_build_only" else 2


if __name__ == "__main__":
    raise SystemExit(main())

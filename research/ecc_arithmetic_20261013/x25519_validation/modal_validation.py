"""Run the frozen X25519 correctness drivers on an ephemeral x86-64 worker."""

import base64
import gzip
import hashlib
import json
import pathlib
import subprocess
import time

import modal


HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[2]
app = modal.App("ecc-x25519-validation-20261013")
image = modal.Image.from_registry("rust:1.93-slim-bookworm", add_python="3.11")
for local, remote in (
    (ROOT / "src/ecc/f25519.rs", "/work/src/ecc/f25519.rs"),
    (ROOT / "src/ecc/x25519.rs", "/work/src/ecc/x25519.rs"),
    (HERE / "ecc.rs", "/work/research/ecc_arithmetic_20261013/x25519_validation/ecc.rs"),
    (HERE / "ecc_f25519_check.rs", "/work/research/ecc_arithmetic_20261013/x25519_validation/ecc_f25519_check.rs"),
    (HERE / "ecc_x25519_check.rs", "/work/research/ecc_arithmetic_20261013/x25519_validation/ecc_x25519_check.rs"),
    (HERE / "ecc_f25519_check.py", "/work/research/ecc_arithmetic_20261013/x25519_validation/ecc_f25519_check.py"),
    (HERE / "ecc_x25519_check.py", "/work/research/ecc_arithmetic_20261013/x25519_validation/ecc_x25519_check.py"),
):
    image = image.add_local_file(str(local), remote_path=remote)


@app.function(image=image, cpu=2, memory=2048, timeout=900)
def validate():
    import os
    import platform
    import resource

    work = pathlib.Path("/work/research/ecc_arithmetic_20261013/x25519_validation")
    cpuinfo = pathlib.Path("/proc/cpuinfo").read_text()
    model = next(
        (line.split(":", 1)[1].strip() for line in cpuinfo.splitlines() if line.startswith("model name")),
        "unknown",
    )
    flags = next(
        (line.split(":", 1)[1].strip() for line in cpuinfo.splitlines() if line.startswith("flags")),
        "",
    )
    start = time.time()
    commands = []

    def run(command, output=None):
        began = time.time()
        if output is None:
            process = subprocess.run(command, cwd=work, capture_output=True, text=True, timeout=300)
        else:
            with (work / output).open("wb") as stream:
                process = subprocess.run(command, cwd=work, stdout=stream, stderr=subprocess.PIPE, timeout=300)
        receipt = {
            "command": command,
            "started_unix": began,
            "ended_unix": time.time(),
            "exit": process.returncode,
            "stdout": process.stdout if output is None else None,
            "stderr": process.stderr.decode() if output is not None else process.stderr,
            "output_path": output,
        }
        commands.append(receipt)
        return process.returncode

    field_source = "ecc_f25519_check.rs"
    x25519_source = "ecc_x25519_check.rs"
    phase_specs = (
        (["rustc", "--edition=2021", "-C", "opt-level=2", field_source, "-o", "/tmp/ecc_f25519_check"], None),
        (["/tmp/ecc_f25519_check"], "field_rows.txt"),
        (["python3", "ecc_f25519_check.py"], None),
        (["rustc", "--edition=2021", "-C", "opt-level=2", x25519_source, "-o", "/tmp/ecc_x25519_check"], None),
        (["/tmp/ecc_x25519_check"], "x25519_rows.txt"),
        (["python3", "ecc_x25519_check.py"], None),
    )
    for command, output in phase_specs:
        if run(command, output) != 0:
            break

    files = {}
    for name in ("field_rows.txt", "x25519_rows.txt"):
        path = work / name
        if path.exists():
            raw = path.read_bytes()
            files[name] = {
                "sha256": hashlib.sha256(raw).hexdigest(),
                "bytes": len(raw),
                "gzip_base64": base64.b64encode(gzip.compress(raw, mtime=0)).decode(),
            }
    sources = {}
    for name in ("f25519.rs", "x25519.rs"):
        raw = (pathlib.Path("/work/src/ecc") / name).read_bytes()
        sources[name] = hashlib.sha256(raw).hexdigest()
    return {
        "started_unix": start,
        "ended_unix": time.time(),
        "platform": platform.platform(),
        "cpu_model": model,
        "cpu_flags_relevant": [
            flag for flag in flags.split() if flag in {"adx", "avx2", "avx512f", "bmi2", "pclmulqdq", "vpclmulqdq"}
        ],
        "affinity": sorted(os.sched_getaffinity(0)),
        "rustc_version": subprocess.run(["rustc", "--version"], capture_output=True, text=True, check=True).stdout.strip(),
        "python_version": platform.python_version(),
        "max_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "sources": sources,
        "commands": commands,
        "files": files,
    }


@app.local_entrypoint()
def main(label: str = "run1"):
    if label not in {"run1", "run2"}:
        raise ValueError("label must be run1 or run2")
    receipt_path = HERE / f"modal_validation_{label}.json"
    if receipt_path.exists():
        raise FileExistsError(receipt_path)
    result = validate.remote()
    receipt_path.write_text(json.dumps(result, indent=2) + "\n")
    for name, receipt in result["files"].items():
        raw = gzip.decompress(base64.b64decode(receipt["gzip_base64"]))
        assert hashlib.sha256(raw).hexdigest() == receipt["sha256"]
        (HERE / f"modal_{label}_{name}").write_bytes(raw)
    print(
        f"label={label} phases={len(result['commands'])}/6 "
        f"all_ok={all(item['exit'] == 0 for item in result['commands'])} "
        f"cpu={result['cpu_model']}"
    )

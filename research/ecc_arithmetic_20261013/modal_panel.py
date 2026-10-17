"""Run the exact paired Rust panel on a bounded Modal x86-64 CPU worker."""

import hashlib
import json
import pathlib
import subprocess

import modal


STUDY = pathlib.Path(__file__).resolve().parent
app = modal.App("ecc-arithmetic-panel-20261013")
image = (
    modal.Image.from_registry("rust:1.93-slim-bookworm", add_python="3.11")
    .add_local_file(str(STUDY / "Cargo.toml"), remote_path="/work/study/Cargo.toml")
    .add_local_file(str(STUDY / "Cargo.lock"), remote_path="/work/study/Cargo.lock")
    .add_local_file(str(STUDY / "src/main.rs"), remote_path="/work/study/src/main.rs")
    .add_local_file(str(STUDY / "src/bin/compare.rs"), remote_path="/work/study/src/bin/compare.rs")
    .add_local_file(str(STUDY / "baseline_src/f2m.rs"), remote_path="/work/study/baseline_src/f2m.rs")
    .add_local_file(str(STUDY / "candidate_src/f2m.rs"), remote_path="/work/src/binary_ecc/f2m.rs")
)


@app.function(image=image, cpu=2, memory=2048, timeout=1200)
def run_panel():
    import os
    import platform

    cpuinfo = pathlib.Path("/proc/cpuinfo").read_text()
    model = next((line.split(":", 1)[1].strip() for line in cpuinfo.splitlines() if line.startswith("model name")), "unknown")
    flags = next((line.split(":", 1)[1].strip() for line in cpuinfo.splitlines() if line.startswith("flags")), "")
    manifest = "/work/study/Cargo.toml"
    base = ["cargo", "--version"]
    cargo_version = subprocess.run(base, capture_output=True, text=True, check=True).stdout.strip()
    rustc_version = subprocess.run(["rustc", "--version"], capture_output=True, text=True, check=True).stdout.strip()
    test = subprocess.run(
        ["cargo", "test", "--release", "--manifest-path", manifest, "--bin", "compare"],
        cwd="/work", capture_output=True, text=True, timeout=600,
    )
    bench = None
    if test.returncode == 0:
        bench = subprocess.run(
            ["cargo", "run", "--release", "--manifest-path", manifest, "--bin", "compare", "--", "500000", "24"],
            cwd="/work", capture_output=True, text=True, timeout=600,
        )
    return {
        "platform": platform.platform(), "cpu_model": model,
        "cpu_flags_relevant": [f for f in flags.split() if f in {"pclmulqdq", "vpclmulqdq", "avx2", "avx512f", "bmi2", "adx"}],
        "affinity": sorted(os.sched_getaffinity(0)),
        "cargo_version": cargo_version, "rustc_version": rustc_version,
        "baseline_sha256": hashlib.sha256(pathlib.Path("/work/study/baseline_src/f2m.rs").read_bytes()).hexdigest(),
        "candidate_sha256": hashlib.sha256(pathlib.Path("/work/src/binary_ecc/f2m.rs").read_bytes()).hexdigest(),
        "test_exit": test.returncode, "test_stdout": test.stdout, "test_stderr": test.stderr,
        "bench_exit": None if bench is None else bench.returncode,
        "bench_stdout": None if bench is None else bench.stdout,
        "bench_stderr": None if bench is None else bench.stderr,
    }


@app.local_entrypoint()
def main():
    result = run_panel.remote()
    (STUDY / "modal_result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"test_exit={result['test_exit']} bench_exit={result['bench_exit']} cpu={result['cpu_model']}")

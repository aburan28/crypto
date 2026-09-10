#!/usr/bin/env python3
"""Measure clean Phase-B tool builds; never run a PDP instance."""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
import platform
import shutil
import stat
import subprocess

import run_koblitz_blind_pdp_phase_b as phase_b


SCHEMA = "koblitz_pdp_phase_b_tool_build_receipt.v1"
RUST_LOCK_PATH = "research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-20-rust-build/Cargo.lock"
RUST_PATHS = ["Cargo.toml", RUST_LOCK_PATH, "src", "examples/koblitz_pdp_export.rs", "examples/koblitz_pdp_backend.rs"]
CMS_COMMIT = "7ae1b4a74259cdce223a584281fb8f090bbd3eed"
CMS_DEPENDENCIES = {
    "cadical": "878a600ea92ce3a9173b069ab294429c851acf85",
    "cadiback": "acc8bb55b98e35afe0dc228662bf41d8341b8a9b",
}
CMS_UNUSED_SUBMODULES = {
    "tests/pbsugar": "3441462ed0af5f3d7d2b6a7ce55938ae305f34f4",
    "tests/simp-checks/simplifiy_testfiles": "345531cca666687830df6f57192e51d1f6ac52cc",
    "tests/xor_cnf_tests": "a0963d220a22e7f1433d7fb3c23d8f34a67a5611",
    "utils/OutputCheck": "3232d0937dcabf59df7ceb776f27b4d5cff6194d",
    "utils/cnf-utils": "b8f1daf34293180613bfc081f919f327a622a561",
    "utils/gtest": "1d17ea141d2c11b8917d2c7d029f1c4e2b9769b2",
    "utils/licensecheck": "74128e614361de8e57854d432fb0d22c28624737",
    "utils/lingeling-ala": "3c733d8d1c75055d4cc1a722ecaccd6130086b1f",
    "utils/minisat": "b58f67e5c1d4cc0be8ad588095a41182dfdfd259",
    "utils/minisat_only_elim_and_subsume": "a114747d6509a53c24c7c9a862d1b9111bc874da",
    "utils/sha1-sat": "e6347c118f183b1ec29d30d05e3f1dc39fccb949",
}
CMS_FLAGS = ["-DCMAKE_BUILD_TYPE=Release", "-DBUILD_SHARED_LIBS=OFF", "-DSTATIC_BINARY=ON", "-DENABLE_TESTING=OFF", "-DNOMPI=ON", "-DNOBREAKID=ON", "-DFETCHCONTENT_FULLY_DISCONNECTED=ON"]
BOUNDARY = {
    "cost_scope": "Source archival, extraction, tool-version probes, configure and compilation; charged once per build, separately from per-instance solver costs",
    "cpu_scope": "RUSAGE_CHILDREN includes waited compiler and linker descendants, including parallel jobs; total core-seconds are not single-core elapsed time",
    "memory_scope": "peak_rss_bytes is the largest child high-water RSS, not simultaneous aggregate RSS of a parallel process tree",
    "excluded": ["prior toolchain and system dependency installation", "prior source acquisition", "aggregate simultaneous build memory", "outer Python build-driver overhead"],
    "single_core_elapsed_seconds": None,
    "full_cost_gate_passed": False,
}


def git(source: Path, *args: str) -> str:
    return subprocess.run(["git", *args], cwd=source, capture_output=True, text=True, check=True).stdout.strip()


def clean_commit(source: Path, expected: str | None = None) -> str:
    commit = git(source, "rev-parse", "HEAD")
    phase_b.require_hex40(commit, "build source revision")
    if git(source, "status", "--porcelain=v1", "--untracked-files=all") or (expected and commit != expected):
        raise phase_b.PhaseBError(f"build source must be clean at {expected or commit}: {source}")
    links = {}
    for line in git(source, "ls-files", "--stage").splitlines():
        metadata, path = line.split("\t", 1)
        mode, object_id, _ = metadata.split()
        if mode == "160000":
            links[path] = object_id
    allowed = CMS_UNUSED_SUBMODULES if expected == CMS_COMMIT else {}
    if links != allowed:
        raise phase_b.PhaseBError("source submodules differ from the fixed unused-test exclusion list")
    return commit


def rust_source_objects(source: Path, revision: str = "HEAD") -> dict[str, str]:
    if revision != "HEAD":
        phase_b.require_hex40(revision, "Rust source revision")
    return {path: git(source, "rev-parse", f"{revision}:{path}") for path in RUST_PATHS}


def required_roles(tool: str) -> list[str]:
    if tool == "rust":
        return ["archive-source", "extract-source", "install-lock", "version-cargo", "version-rustc", "vendor", "build", "copy-exporter", "copy-backend"]
    return ["archive-source", "extract-source", "archive-cadical", "extract-cadical", "archive-cadiback", "extract-cadiback", "version-cmake", "version-cc", "version-cxx", "configure", "build", "copy-cryptominisat"]


def tool_identity(path: Path) -> dict:
    # System toolchains may use hard links or rustup's argv[0] symlinks.
    # Hash the resolved file while retaining the invocation name.
    resolved = path.resolve(strict=True)
    if not stat.S_ISREG(resolved.stat().st_mode) or not os.access(path, os.X_OK):
        raise phase_b.PhaseBError(f"build tool is not an executable regular file: {path}")
    data = resolved.read_bytes()
    return {"path": str(path.absolute()), "resolved_path": str(resolved), "bytes": len(data), "sha256": phase_b.sha256_bytes(data)}


def resource_summary(processes: list[dict]) -> dict:
    metrics = [item["metrics"] for item in processes]
    return {
        "total_core_seconds": round(math.fsum(item["total_core_seconds"] for item in metrics), 12),
        "summed_process_wall_seconds": round(math.fsum(item["wall_seconds"] for item in metrics), 12),
        "largest_child_peak_rss_bytes": max((item["peak_rss_bytes"] for item in metrics), default=0),
        "single_core_elapsed_seconds": None,
        "aggregate_process_tree_peak_rss_bytes": None,
    }


def build(tool: str, source: Path, output: Path, jobs: int, timeout: int) -> dict:
    source, output = source.resolve(), output.resolve()
    if tool not in ("rust", "cryptominisat") or type(jobs) is not int or not 1 <= jobs <= 64 or timeout <= 0:
        raise phase_b.PhaseBError("invalid tool, job count, or build timeout")
    if output.exists() or output.is_symlink():
        raise phase_b.PhaseBError(f"build output must be new: {output}")
    revision = clean_commit(source, CMS_COMMIT if tool == "cryptominisat" else None)
    implementation_state = phase_b.git_state()
    phase_b_helper_sha256 = phase_b.sha256_file(Path(phase_b.__file__), "Phase-B helper")
    objects = rust_source_objects(source) if tool == "rust" else None
    dependencies = {}
    if tool == "cryptominisat":
        for name, expected in CMS_DEPENDENCIES.items():
            dependency = source / f"build/_deps/{name}-src"
            dependencies[name] = {"path": str(dependency), "commit": clean_commit(dependency, expected)}
    names = ["git", "tar", "cp", "cargo", "rustc"] if tool == "rust" else ["git", "tar", "cp", "cmake", "cc", "c++", "make"]
    commands = {}
    for name in names:
        found = shutil.which(name)
        if found is None:
            raise phase_b.PhaseBError(f"required build tool missing: {name}")
        commands[name] = tool_identity(Path(found))
    output.mkdir(parents=False)
    evidence = output / "evidence"
    evidence.mkdir()
    source_root = output / "source"
    source_root.mkdir()
    build_root = output / "build"
    environment = phase_b.safe_child_environment()
    environment["PATH"] = os.pathsep.join(dict.fromkeys([str(Path(item["path"]).parent) for item in commands.values()] + environment["PATH"].split(os.pathsep)))
    # Vendoring may read the installed source cache. Compilation later uses a
    # fresh CARGO_HOME, a frozen vendor configuration and a new target directory.
    if tool == "rust":
        environment.update({"HOME": str(Path.home()), "CARGO_HOME": str(Path(os.environ.get("CARGO_HOME", Path.home() / ".cargo")).resolve()), "RUSTC": commands["rustc"]["path"], "CARGO_INCREMENTAL": "0"})
    processes = []

    def measured(role: str, command: list[str], cwd: Path) -> dict:
        record = phase_b.run_metered(role=role, command=command, cwd=cwd, task_root=evidence, input_paths=[], timeout=timeout, meter=phase_b.DEFAULT_METER, environment=environment, markers=set())
        if record["returncode"] != 0 or record["timed_out"]:
            raise phase_b.PhaseBError(f"{tool} {role} failed; retain {evidence} and use a new output directory")
        record["environment"] = dict(environment)
        processes.append(phase_b.compact_process(record))
        return record

    archives = {}

    def archive(name: str, checkout: Path, commit: str, destination: Path, paths: list[str]) -> None:
        tar_path = output / f"{name}.tar"
        measured(f"archive-{name}", [commands["git"]["path"], "archive", "--format=tar", f"--output={tar_path}", commit, *paths], checkout)
        archives[name] = phase_b.executable_identity(tar_path, "source archive", executable=False)
        measured(f"extract-{name}", [commands["tar"]["path"], "-xf", str(tar_path), "-C", str(destination)], output)

    archive("source", source, revision, source_root, RUST_PATHS if tool == "rust" else [])
    if tool == "rust":
        measured("install-lock", [commands["cp"]["path"], str(source_root / RUST_LOCK_PATH), str(source_root / "Cargo.lock")], output)
        shutil.copyfile(source_root / "Cargo.lock", evidence / "Cargo.lock")
        measured("version-cargo", [commands["cargo"]["path"], "--version", "--verbose"], source_root)
        measured("version-rustc", [commands["rustc"]["path"], "--version", "--verbose"], source_root)
        vendor = measured("vendor", [commands["cargo"]["path"], "vendor", "--locked", "--offline", "--versioned-dirs", str(output / "vendor")], source_root)
        (source_root / ".cargo").mkdir()
        (source_root / ".cargo/config.toml").write_text(vendor["stdout_text"])
        (evidence / "cargo-config.toml").write_text(vendor["stdout_text"])
        cargo_home = output / "cargo-home"
        cargo_home.mkdir()
        environment["CARGO_HOME"] = str(cargo_home)
        measured("build", [commands["cargo"]["path"], "build", "--release", "--locked", "--offline", "--jobs", str(jobs), "--target-dir", str(build_root), "--example", "koblitz_pdp_export", "--example", "koblitz_pdp_backend"], source_root)
        if (source_root / "Cargo.lock").read_bytes() != (evidence / "Cargo.lock").read_bytes():
            raise phase_b.PhaseBError("the frozen Rust lock snapshot changed during compilation")
        binaries = {name: build_root / "release/examples" / binary for name, binary in {"exporter": "koblitz_pdp_export", "backend": "koblitz_pdp_backend"}.items()}
    else:
        for name, dependency in dependencies.items():
            destination = output / name
            destination.mkdir()
            archive(name, Path(dependency["path"]), dependency["commit"], destination, [])
        for label, executable in (("cmake", "cmake"), ("cc", "cc"), ("cxx", "c++")):
            measured(f"version-{label}", [commands[executable]["path"], "--version"], source_root)
        measured("configure", [commands["cmake"]["path"], "-S", str(source_root), "-B", str(build_root), "-G", "Unix Makefiles", *CMS_FLAGS, f"-DCMAKE_C_COMPILER={commands['cc']['path']}", f"-DCMAKE_CXX_COMPILER={commands['c++']['path']}", f"-DFETCHCONTENT_SOURCE_DIR_CADICAL={output / 'cadical'}", f"-DFETCHCONTENT_SOURCE_DIR_CADIBACK={output / 'cadiback'}"], output)
        shutil.copyfile(build_root / "CMakeCache.txt", evidence / "CMakeCache.txt")
        measured("build", [commands["cmake"]["path"], "--build", str(build_root), "--target", "cryptominisat5-bin", "--parallel", str(jobs)], output)
        binaries = {"cryptominisat": build_root / "cryptominisat5"}
    packaged = output / "bin"
    packaged.mkdir()
    for name, binary in list(binaries.items()):
        destination = packaged / binary.name
        measured(f"copy-{name}", [commands["cp"]["path"], str(binary), str(destination)], output)
        os.chmod(destination, 0o755)
        binaries[name] = destination
    clean_commit(source, revision)
    for name, identity in commands.items():
        if tool_identity(Path(identity["path"])) != identity:
            raise phase_b.PhaseBError(f"build tool changed during compilation: {name}")
    for dependency in dependencies.values():
        clean_commit(Path(dependency["path"]), dependency["commit"])
    if phase_b.git_state() != implementation_state:
        raise phase_b.PhaseBError("Phase-B implementation state changed during the tool build")
    if phase_b.sha256_file(Path(phase_b.__file__), "Phase-B helper") != phase_b_helper_sha256:
        raise phase_b.PhaseBError("Phase-B helper changed during the tool build")
    receipt = {
        "schema": SCHEMA, "status": "completed", "tool": tool,
        "source_commit": revision, "source_clean": True, "rust_source_objects": objects,
        "dependency_commits": {name: item["commit"] for name, item in dependencies.items()},
        "excluded_test_submodules": CMS_UNUSED_SUBMODULES if tool == "cryptominisat" else {},
        "source_archives": archives, "requested_parallel_jobs": jobs,
        "toolchain": commands, "platform": {"system": platform.system(), "machine": platform.machine(), "release": platform.release()},
        "environment": environment, "build_processes": processes,
        "binaries": {name: phase_b.executable_identity(path, name) for name, path in binaries.items()},
        "resources": resource_summary(processes), "accounting_boundary": BOUNDARY,
        "evidence_inventory": phase_b.all_regular_inventory(evidence),
        "build_driver_sha256": phase_b.sha256_file(Path(__file__), "build driver"),
        "process_meter_sha256": phase_b.sha256_file(phase_b.DEFAULT_METER, "process meter"),
        "phase_b_helper_sha256": phase_b_helper_sha256,
        "implementation_state": implementation_state,
    }
    receipt["receipt_payload_sha256"] = phase_b.canonical_sha256(receipt)
    phase_b.write_json_new(output / "receipt.json", receipt)
    validate_receipt(
        output / "receipt.json", tool, receipt["binaries"], objects,
        implementation_state, require_clean_implementation=False,
    )
    return receipt


def validate_receipt(
    path: Path,
    tool: str,
    identities: dict,
    source_objects: dict | None = None,
    expected_implementation_state: dict | None = None,
    require_clean_implementation: bool = False,
) -> dict:
    value, raw = phase_b.read_json(path, "tool build receipt")
    payload = dict(value)
    self_hash = payload.pop("receipt_payload_sha256", None)
    if phase_b.canonical_sha256(payload) != self_hash:
        raise phase_b.PhaseBError("tool build receipt payload hash differs")
    trusted_implementation = {
        "build_driver_sha256": phase_b.sha256_file(Path(__file__), "build driver"),
        "process_meter_sha256": phase_b.sha256_file(phase_b.DEFAULT_METER, "process meter"),
        "phase_b_helper_sha256": phase_b.sha256_file(Path(phase_b.__file__), "Phase-B helper"),
    }
    for field, expected in trusted_implementation.items():
        if value.get(field) != expected:
            label = field.removesuffix("_sha256").replace("_", " ")
            raise phase_b.PhaseBError(
                f"tool build receipt {label} differs from the trusted current implementation"
            )
    implementation_state = value.get("implementation_state")
    if (
        not isinstance(implementation_state, dict)
        or set(implementation_state) != {"commit", "dirty", "porcelain"}
        or not isinstance(implementation_state.get("dirty"), bool)
        or not isinstance(implementation_state.get("porcelain"), list)
        or not all(isinstance(line, str) for line in implementation_state["porcelain"])
        or implementation_state["dirty"] != bool(implementation_state["porcelain"])
    ):
        raise phase_b.PhaseBError("tool build receipt has invalid implementation provenance")
    phase_b.require_hex40(implementation_state.get("commit"), "tool build implementation revision")
    if expected_implementation_state is not None and implementation_state != expected_implementation_state:
        raise phase_b.PhaseBError("tool build receipt differs from the expected implementation state")
    if require_clean_implementation and (
        implementation_state["dirty"] or implementation_state["porcelain"]
    ):
        raise phase_b.PhaseBError("production tool build requires a clean implementation state")
    if value.get("schema") != SCHEMA or value.get("status") != "completed" or value.get("tool") != tool or value.get("source_clean") is not True:
        raise phase_b.PhaseBError("tool build receipt is not a clean completed build of the requested tool")
    phase_b.require_hex40(value.get("source_commit"), "tool build source revision")
    if value.get("accounting_boundary") != BOUNDARY:
        raise phase_b.PhaseBError("tool build receipt changed its accounting boundary")
    if type(value.get("requested_parallel_jobs")) is not int or not 1 <= value["requested_parallel_jobs"] <= 64:
        raise phase_b.PhaseBError("tool build receipt has invalid parallel job count")
    expected_binaries = ("exporter", "backend") if tool == "rust" else ("cryptominisat",)
    if set(value.get("binaries", {})) != set(expected_binaries):
        raise phase_b.PhaseBError("tool build receipt has the wrong output binaries")
    for name in expected_binaries:
        if any(value["binaries"][name].get(key) != identities[name].get(key) for key in ("sha256", "bytes")):
            raise phase_b.PhaseBError("tool build receipt does not bind the requested executable")
    if tool == "cryptominisat" and (value["source_commit"] != CMS_COMMIT or value.get("dependency_commits") != CMS_DEPENDENCIES):
        raise phase_b.PhaseBError("CryptoMiniSat build source or dependency pins differ")
    if value.get("excluded_test_submodules") != (CMS_UNUSED_SUBMODULES if tool == "cryptominisat" else {}):
        raise phase_b.PhaseBError("tool build omitted submodules outside the fixed unused-test exclusion list")
    if tool == "rust" and (set(value.get("rust_source_objects") or {}) != set(RUST_PATHS) or (source_objects is not None and value["rust_source_objects"] != source_objects)):
        raise phase_b.PhaseBError("Rust build inputs differ from the solver implementation")
    evidence = path.parent / "evidence"
    if phase_b.all_regular_inventory(evidence) != value.get("evidence_inventory"):
        raise phase_b.PhaseBError("tool build evidence inventory differs")
    processes = value.get("build_processes")
    if not isinstance(processes, list) or [item.get("role") for item in processes] != required_roles(tool):
        raise phase_b.PhaseBError("tool build receipt lacks required separately metered stages")
    validate_recipe(value)
    expected_files = set()
    for process in processes:
        role = process["role"]
        for suffix in ("intent.json", "metrics.json", "stdout", "stderr"):
            expected_files.add(f"{role}.{suffix}")
        metrics, _ = phase_b.read_json(evidence / f"{role}.metrics.json", "build process metrics")
        intent, _ = phase_b.read_json(evidence / f"{role}.intent.json", "build process intent")
        if any(metrics.get(key) != process.get(key) for key in ("command", "returncode", "timed_out", "orphan_group_terminated", "metrics")) or intent.get("command") != process.get("command") or intent.get("environment") != process.get("environment"):
            raise phase_b.PhaseBError("tool build process differs from its raw receipt or intent")
        if metrics.get("returncode") != 0 or metrics.get("timed_out") is not False or metrics.get("orphan_group_terminated") is not False:
            raise phase_b.PhaseBError("tool build process did not complete cleanly")
        resources = metrics["metrics"]
        for key in ("wall_seconds", "user_seconds", "system_seconds", "total_core_seconds", "single_core_seconds", "peak_rss_bytes"):
            number = resources.get(key)
            if type(number) not in (int, float) or not math.isfinite(number) or number < 0:
                raise phase_b.PhaseBError("tool build resources contain invalid values")
        if not math.isclose(resources["total_core_seconds"], resources["user_seconds"] + resources["system_seconds"], abs_tol=1e-9, rel_tol=0) or resources["single_core_seconds"] != resources["total_core_seconds"] or resources.get("meter") != "fresh-process getrusage(RUSAGE_CHILDREN)":
            raise phase_b.PhaseBError("tool build resource accounting is inconsistent")
        for suffix in ("stdout", "stderr"):
            if phase_b.sha256_file(evidence / f"{role}.{suffix}", "build output") != process.get(f"{suffix}_sha256"):
                raise phase_b.PhaseBError("tool build output digest differs")
    if tool == "cryptominisat":
        expected_files.add("CMakeCache.txt")
    else:
        expected_files.add("cargo-config.toml")
        expected_files.add("Cargo.lock")
        if (evidence / "cargo-config.toml").read_bytes() != (evidence / "vendor.stdout").read_bytes():
            raise phase_b.PhaseBError("Cargo configuration differs from the metered vendor output")
    if {item["path"] for item in value["evidence_inventory"]} != expected_files:
        raise phase_b.PhaseBError("tool build evidence has unexpected or missing files")
    if value.get("resources") != resource_summary(processes):
        raise phase_b.PhaseBError("tool build aggregate resources are inconsistent")
    return {"receipt_sha256": phase_b.sha256_bytes(raw), "receipt": value}


def validate_recipe(value: dict) -> None:
    """Check measured commands implement this recipe, not just named stages."""
    tool = value["tool"]
    commands = value["toolchain"]
    names = {"git", "tar", "cp", "cargo", "rustc"} if tool == "rust" else {"git", "tar", "cp", "cmake", "cc", "c++", "make"}
    if set(commands) != names:
        raise phase_b.PhaseBError("tool build toolchain identities differ from the recipe")
    for identity in commands.values():
        phase_b.require_hex64(identity.get("sha256"), "build tool identity")
    archives = value["source_archives"]
    if set(archives) != ({"source"} if tool == "rust" else {"source", *CMS_DEPENDENCIES}):
        raise phase_b.PhaseBError("tool build source archives differ from the recipe")
    output = Path(archives["source"]["path"]).parent
    source_root, build_root = output / "source", output / "build"
    expected = {}
    for name, archive in archives.items():
        phase_b.require_hex64(archive.get("sha256"), "build source archive identity")
        commit = value["source_commit"] if name == "source" else CMS_DEPENDENCIES[name]
        destination = source_root if name == "source" else output / name
        paths = RUST_PATHS if tool == "rust" else []
        expected[f"archive-{name}"] = [commands["git"]["path"], "archive", "--format=tar", f"--output={output / (name + '.tar')}", commit, *paths]
        expected[f"extract-{name}"] = [commands["tar"]["path"], "-xf", str(output / (name + ".tar")), "-C", str(destination)]
    jobs = str(value["requested_parallel_jobs"])
    if tool == "rust":
        expected["install-lock"] = [commands["cp"]["path"], str(source_root / RUST_LOCK_PATH), str(source_root / "Cargo.lock")]
        for name in ("cargo", "rustc"):
            expected[f"version-{name}"] = [commands[name]["path"], "--version", "--verbose"]
        expected["vendor"] = [commands["cargo"]["path"], "vendor", "--locked", "--offline", "--versioned-dirs", str(output / "vendor")]
        expected["build"] = [commands["cargo"]["path"], "build", "--release", "--locked", "--offline", "--jobs", jobs, "--target-dir", str(build_root), "--example", "koblitz_pdp_export", "--example", "koblitz_pdp_backend"]
        binary_paths = {name: build_root / "release/examples" / binary for name, binary in {"exporter": "koblitz_pdp_export", "backend": "koblitz_pdp_backend"}.items()}
    else:
        for label, executable in (("cmake", "cmake"), ("cc", "cc"), ("cxx", "c++")):
            expected[f"version-{label}"] = [commands[executable]["path"], "--version"]
        expected["configure"] = [commands["cmake"]["path"], "-S", str(source_root), "-B", str(build_root), "-G", "Unix Makefiles", *CMS_FLAGS, f"-DCMAKE_C_COMPILER={commands['cc']['path']}", f"-DCMAKE_CXX_COMPILER={commands['c++']['path']}", f"-DFETCHCONTENT_SOURCE_DIR_CADICAL={output / 'cadical'}", f"-DFETCHCONTENT_SOURCE_DIR_CADIBACK={output / 'cadiback'}"]
        expected["build"] = [commands["cmake"]["path"], "--build", str(build_root), "--target", "cryptominisat5-bin", "--parallel", jobs]
        binary_paths = {"cryptominisat": build_root / "cryptominisat5"}
    for name, binary in binary_paths.items():
        expected[f"copy-{name}"] = [commands["cp"]["path"], str(binary), str(output / "bin" / binary.name)]
    for process in value["build_processes"]:
        if process["command"] != expected[process["role"]]:
            raise phase_b.PhaseBError("tool build command differs from the frozen build recipe")
        if tool == "rust" and process["role"] == "build" and process["environment"].get("CARGO_HOME") != str(output / "cargo-home"):
            raise phase_b.PhaseBError("Rust compilation did not use its fresh Cargo home")


def copy_capsule(path: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(path, destination / "receipt.json")
    shutil.copytree(path.parent / "evidence", destination / "evidence")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tool", choices=("rust", "cryptominisat"))
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=1800)
    args = parser.parse_args()
    try:
        result = build(args.tool, args.source.resolve(), args.output.resolve(), args.jobs, args.timeout)
        print(phase_b.pretty_bytes({"receipt": str(args.output / "receipt.json"), "resources": result["resources"], "binaries": result["binaries"]}).decode(), end="")
    except (OSError, subprocess.CalledProcessError, phase_b.PhaseBError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()

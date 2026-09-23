#!/usr/bin/env python3
"""Run the fixed, standalone Boolean schedule protocol once into a new directory."""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import tempfile
import time

HERE = Path(__file__).resolve().parent


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def child(command: list[str], out: Path, err: Path, timeout: float) -> dict:
    started = time.perf_counter_ns()
    with out.open("wb") as stdout, err.open("wb") as stderr:
        proc = subprocess.Popen(command, stdout=stdout, stderr=stderr)
        timed_out = False
        while True:
            pid, status, usage = os.wait4(proc.pid, os.WNOHANG)
            if pid:
                break
            if (time.perf_counter_ns() - started) / 1e9 >= timeout:
                timed_out = True
                proc.kill()
                _, status, usage = os.wait4(proc.pid, 0)
                break
            time.sleep(0.01)
        proc.returncode = os.waitstatus_to_exitcode(status)
    return {
        "command": command,
        "exit_code": proc.returncode,
        "timed_out": timed_out,
        "process_wall_ns": time.perf_counter_ns() - started,
        "user_seconds": usage.ru_utime,
        "system_seconds": usage.ru_stime,
        "peak_rss_bytes": usage.ru_maxrss * (1 if platform.system() == "Darwin" else 1024),
        "meter": "wait4 for this fresh worker, not cumulative RUSAGE_CHILDREN",
        "stdout_sha256": sha(out),
        "stderr_sha256": sha(err),
    }


def worker_command(executable, n, seed, batch, family, protocol):
    command = [str(executable), str(n), str(seed), str(batch), family, str(protocol["repetitions"])]
    if protocol.get("enable_census"):
        command.append("with-census")
    if protocol.get("balanced_order"):
        if not protocol.get("enable_census"):
            raise ValueError("balanced order requires the census arm set")
        command.append("balanced-order")
    return command


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, default=HERE / "protocol.json")
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    protocol = json.loads(args.protocol.read_text())
    for name in ("worker.rs", "protocol.json", "run.py", "analyze.py"):
        shutil.copyfile(args.protocol if name == "protocol.json" else HERE / name, out / name)
    source_hashes = {name: sha(out / name) for name in ("worker.rs", "protocol.json", "run.py", "analyze.py")}
    rustc = shutil.which("rustc")
    if not rustc:
        raise SystemExit("rustc is required; no dependencies are downloaded")
    metadata = {
        "schema_version": 1,
        "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "source_hashes": source_hashes,
        "rustc": subprocess.check_output([rustc, "--version", "--verbose"], text=True),
        "git_base": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=HERE, text=True).strip(),
        "git_status": subprocess.check_output(["git", "status", "--short"], cwd=HERE, text=True),
        "host": {"system": platform.system(), "machine": platform.machine(), "platform": platform.platform(), "cpus": os.cpu_count()},
        "scope": protocol["scope"],
        "complete": False,
    }
    dump(out / "metadata.json", metadata)
    receipts = []
    with tempfile.TemporaryDirectory(prefix="boolean-schedule-") as temporary:
        build = Path(temporary)
        executable, tests = build / "worker", build / "tests"
        for target, flags in [(executable, []), (tests, ["--test"])]:
            command = [rustc, "--edition", "2021", "-O", *flags, str(out / "worker.rs"), "-o", str(target)]
            result = subprocess.run(command, capture_output=True, timeout=60)
            (out / (target.name + "-compile.txt")).write_bytes(result.stdout + result.stderr)
            if result.returncode:
                raise SystemExit("compile failed; retained compiler output")
        metadata["worker_sha256"] = sha(executable)
        metadata["test_executable_sha256"] = sha(tests)
        test_stdout, test_stderr = build / "tests.stdout", build / "tests.stderr"
        test_receipt = child([str(tests)], test_stdout, test_stderr, 15)
        test_receipt["stdout"] = test_stdout.read_text()
        test_receipt["stderr"] = test_stderr.read_text()
        dump(out / "test_receipt.json", test_receipt)
        if test_receipt["exit_code"] != 0:
            raise SystemExit("correctness tests failed; campaign not launched")
        dump(out / "metadata.json", metadata)
        started = time.monotonic()
        for n in protocol["variables"]:
            with (out / f"raw-n{n}.jsonl").open("x") as raw:
                for split in ("discovery", "holdout"):
                    for seed in protocol[f"{split}_seeds"]:
                        for family in protocol["families"]:
                            for batch in protocol["batches"]:
                                if time.monotonic() - started > protocol["limits"]["campaign_seconds"]:
                                    raise SystemExit("campaign cap reached; incomplete evidence retained")
                                cell = f"n{n}-{split}-{seed}-{family}-b{batch}"
                                command = worker_command(executable, n, seed, batch, family, protocol)
                                stdout, stderr = build / "cell.stdout", build / "cell.stderr"
                                receipt = child(command, stdout, stderr, protocol["limits"]["worker_seconds"])
                                receipt.update(cell=cell, split=split, n=n, family=family, batch=batch, seed=seed)
                                receipts.append(receipt)
                                dump(out / "receipts.json", receipts)
                                if receipt["exit_code"] != 0:
                                    shutil.copyfile(stdout, out / "failed.stdout")
                                    shutil.copyfile(stderr, out / "failed.stderr")
                                    raise SystemExit(f"{cell} failed; censored evidence retained")
                                # Preserve raw worker bytes in the JSON records, allowing
                                # exact stdout-hash reconstruction by the analyzer.
                                for line in stdout.read_text().splitlines(keepends=True):
                                    row = json.loads(line)
                                    row.update(cell=cell, split=split, raw_line=line)
                                    raw.write(json.dumps(row, separators=(",", ":")) + "\n")
                                raw.flush()
            print(f"completed n={n}: {len(receipts)} fixed cells", flush=True)
        metadata["complete"] = True
        metadata["ended_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        metadata["campaign_seconds"] = time.monotonic() - started
        dump(out / "metadata.json", metadata)
    subprocess.run([shutil.which("python3"), str(out / "analyze.py"), str(out)], check=True)
    dump(out / "manifest.json", {
        "files": {str(path.relative_to(out)): sha(path) for path in sorted(out.rglob("*")) if path.is_file() and path.name != "manifest.json"},
        "scope": protocol["scope"],
    })


if __name__ == "__main__":
    main()

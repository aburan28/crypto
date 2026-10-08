#!/usr/bin/env python3
"""The `ic` tool programme's runner (IC_TOOL_PROGRAM.md §5), shared by
every round.

Every timed process is one `ic price --single-target` on one suite row,
run through `tools/isolated_bench.py run --wait --cpus 2` with
`RAYON_NUM_THREADS=1`, after PSI `some avg10` has fallen below 4.0, as
in ledger §23.  A refused start is logged and retried after 15 s.  A
contended or failed process is kept and run again, at most twice; the
first clean run is the row's figure, and every attempt stays on disk.
An existing output is never overwritten, so every command resumes where
the last one stopped.

Arms are named binaries.  `interleave` runs each row once per arm per
round, alternating the arms' order from round to round (A B, B A, ...).

The module is imported by each round's `run.py`; it has no command line
of its own.
"""
from __future__ import annotations

import functools
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TOOL = ROOT / "tools" / "isolated_bench.py"
SUITE = HERE.parent / "suite" / "v1"
CPUS = "2"
RETRIES = 2
REFUSAL_WAIT_S, REFUSAL_LIMIT = 15, 60
PSI_READY, PSI_WAIT_MAX_S = 4.0, 120
ENV = {**os.environ, "RAYON_NUM_THREADS": "1"}


def sh(cmd: list[str]) -> str:
    try:
        return subprocess.run(cmd, capture_output=True, text=True, check=False).stdout.strip()
    except FileNotFoundError:
        return ""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def suite_rows(tier: str = "S") -> list[dict]:
    rows = json.loads((SUITE / "SUITE.json").read_text())["rows"]
    return [r for r in rows if r["tier"] == tier]


@functools.cache
def curve_slug(a: int, n: int) -> str:
    """The ICV1 slug (AGENTS.md §11) of the suite's Koblitz curve with
    coefficient `a` over `GF(2^n)`, from `docs/curves/registry.json`."""
    sys.path.insert(0, str(ROOT / "scripts"))
    import curve_id  # noqa: E402

    curve = curve_id.resolve(f"K_{a} / GF(2^{n})")
    if curve is None:
        raise SystemExit(f"docs/curves/registry.json has no Koblitz curve a={a} over GF(2^{n})")
    return curve["slug"]


def slug_rows(rows: list[dict]) -> list[dict]:
    """Rows for a new run tree, keyed `<slug>/<recipe>-<target>`, as AGENTS.md
    §11 names run files.  The suite's frozen id stays as `suite_id`."""
    return [{**r, "suite_id": r["id"], "id": f"{curve_slug(r['a'], r['n'])}/{r['id'].split('-', 1)[1]}"}
            for r in rows]


def psi_some_avg10() -> float:
    worst = 0.0
    for kind in ("cpu", "memory"):
        try:
            for line in Path(f"/proc/pressure/{kind}").read_text().splitlines():
                if line.startswith("some"):
                    worst = max(worst, float(line.split("avg10=")[1].split()[0]))
        except OSError:
            pass
    return worst


def wait_for_quiet() -> float:
    start = time.monotonic()
    while psi_some_avg10() >= PSI_READY and time.monotonic() - start < PSI_WAIT_MAX_S:
        time.sleep(1)
    return time.monotonic() - start


def stem(out: Path) -> Path:
    name = out.name
    for suffix in (".price.json", ".json"):
        if name.endswith(suffix):
            return out.with_name(name[: -len(suffix)])
    return out


def record_path(out: Path) -> Path:
    return stem(out).with_name(stem(out).name + ".isolation.jsonl")


def run_record(out: Path) -> dict | None:
    rec = record_path(out)
    return json.loads(rec.read_text().splitlines()[-1])["run"] if rec.exists() else None


def clean(out: Path) -> bool:
    run = run_record(out)
    return run is not None and run["exit_status"] == 0 and not run["contended"]


def launch(cmd: list[str], out: Path, log_dir: Path, env: dict | None = None) -> None:
    """One timed process through the isolation tool, after PSI has fallen."""
    rec = record_path(out)
    label = str(out.relative_to(log_dir))
    err = stem(out).with_name(stem(out).name + ".stderr")
    for attempt in range(REFUSAL_LIMIT):
        waited = wait_for_quiet()
        tmp = err.with_suffix(".stderr.tmp")
        with open(tmp, "w") as e:
            subprocess.run([sys.executable, str(TOOL), "run", "--wait", "--cpus", CPUS, "--out", str(rec),
                            "--label", label, "--", *cmd], env=env or ENV, stdout=subprocess.DEVNULL,
                           stderr=e, check=False)
        if rec.exists():
            tmp.replace(err)
            with open(log_dir / "psi-waits.log", "a") as log:
                log.write(f"{label} {waited:.1f}\n")
            return
        with open(log_dir / "refusals.log", "a") as log:
            log.write(f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {label} attempt {attempt + 1} "
                      f"(after {waited:.1f} s of PSI wait): {tmp.read_text().strip()}\n")
        tmp.unlink()
        time.sleep(REFUSAL_WAIT_S)
    raise SystemExit(f"{label}: the machine stayed busy for {REFUSAL_LIMIT} tries; stopped (resumable)")


def load(out: Path) -> dict:
    try:
        return json.loads(out.read_text())
    except (OSError, json.JSONDecodeError):
        return {"status": "no report"}


def price_cmd(binary: Path, row: dict, out: Path, extra: tuple[str, ...] = ()) -> list[str]:
    return [str(binary), "price", "--params", str(SUITE / row["params"]), "--json", "--out", str(out),
            "--single-target", "--rho-seed", str(row["rho_seed"]), *extra]


def price(binary: Path, row: dict, out: Path, log_dir: Path, extra: tuple[str, ...] = (),
          env: dict | None = None) -> dict:
    """The first clean run of up to 1 + RETRIES; every attempt is kept."""
    attempt = out
    for k in range(RETRIES + 1):
        if k:
            attempt = out.with_name(stem(out).name + f"-retry{k}.price.json")
        if not attempt.exists():
            attempt.parent.mkdir(parents=True, exist_ok=True)
            launch(price_cmd(binary, row, attempt, extra), attempt, log_dir, env)
        rep = load(attempt)
        if clean(attempt) and rep.get("status") == "complete":
            return rep
    return rep


def figure_path(out: Path) -> Path:
    """The attempt that is the row's figure: the first clean, complete one."""
    for k in range(RETRIES + 1):
        attempt = out if k == 0 else out.with_name(stem(out).name + f"-retry{k}.price.json")
        if attempt.exists() and clean(attempt) and load(attempt).get("status") == "complete":
            return attempt
    return out


def interleave(arms: dict[str, Path], rows: list[dict], rounds: int, out_dir: Path,
               envs: dict[str, dict] | None = None, extra: tuple[str, ...] = ()) -> None:
    """Every row, every round, every arm; the arms' order alternates by round."""
    names = list(arms)
    for k in range(1, rounds + 1):
        order = names if k % 2 else names[::-1]
        for row in rows:
            for arm in order:
                out = out_dir / arm / row["id"] / f"r{k}.price.json"
                rep = price(arms[arm], row, out, out_dir, extra, (envs or {}).get(arm))
                status = rep.get("status")
                cold = rep.get("median", {})
                note = (f"setup {cold['s_setup']:.3f} online {cold['s_ic_online']:.4f}"
                        if status == "complete" else "")
                print(f"r{k} {row['id']} {arm}: {status} {note}", flush=True)


def untimed(binary: Path, row: dict, out: Path, extra: tuple[str, ...] = ()) -> dict:
    """Counts only, under taskset: for pins, never for time."""
    if not out.exists():
        out.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(["taskset", "-c", CPUS, *price_cmd(binary, row, out, extra)], env=ENV,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
    return load(out)


def outputs(rep: dict) -> dict:
    """Everything a pin compares: what the run decided, never how long it took."""
    return {
        "status": rep.get("status"),
        "counts": rep.get("counts"),
        "recovered_ic": rep.get("certificates", {}).get("ic", {}).get("scalar"),
        "recovered_rho": rep.get("certificates", {}).get("rho", {}).get("scalar"),
        "rho_counts": rep.get("rho_counts"),
        "all_verified": rep.get("all_verified"),
        "ic_and_rho_agree": rep.get("ic_and_rho_agree"),
    }


def host_manifest(out: Path, binaries: dict[str, dict]) -> dict:
    if out.exists():
        return json.loads(out.read_text())
    cpu = {}
    for line in Path("/proc/cpuinfo").read_text().splitlines():
        if ":" in line:
            k, v = (s.strip() for s in line.split(":", 1))
            if k in ("model name", "flags") and k not in cpu:
                cpu[k] = v
    flags = cpu.get("flags", "").split()
    wanted = [f for f in flags if f in ("popcnt", "avx2", "pclmulqdq", "bmi2", "gfni", "vpclmulqdq")
              or f.startswith("avx512")]
    mem = next((line for line in Path("/proc/meminfo").read_text().splitlines() if line.startswith("MemTotal")), "")
    thp = {}
    for name in ("enabled", "defrag"):
        try:
            thp[name] = Path(f"/sys/kernel/mm/transparent_hugepage/{name}").read_text().strip()
        except OSError:
            thp[name] = None
    lock = ROOT / "Cargo.lock"
    doc = {
        "repository_commit": sh(["git", "-C", str(ROOT), "rev-parse", "HEAD"]),
        "tree_status": sh(["git", "-C", str(ROOT), "status", "--porcelain"]),
        "rustc": sh(["rustc", "--version"]),
        "cargo_lock_sha256": sha256(lock) if lock.exists() else None,
        "cpu_model": cpu.get("model name"),
        "cpu_flags_relevant": wanted,
        "logical_cores": os.cpu_count(),
        "memory": mem,
        "transparent_hugepage": thp,
        "os": platform.platform(),
        "arch": platform.machine(),
        "glibc": " ".join(platform.libc_ver()),
        "binaries": binaries,
        "isolation_tool": {"path": "tools/isolated_bench.py", "sha256": sha256(TOOL)},
        "pinning": ("tools/isolated_bench.py run --wait --cpus 2 with RAYON_NUM_THREADS=1, after PSI some "
                    "avg10 falls below 4.0; the tool's defaults (settle 2 s, other processes at most 0.10 "
                    "CPUs, PSI at most 5); pins untimed under taskset -c 2"),
        "smt_siblings": {c: Path(f"/sys/devices/system/cpu/cpu{c}/topology/thread_siblings_list").read_text().strip()
                         for c in range(os.cpu_count() or 0)},
        "uptime": sh(["uptime"]),
        "hardware_class": "one x86-64 cloud container; no claim for Arm64, GPUs or other hosts",
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(doc, indent=1) + "\n")
    return doc

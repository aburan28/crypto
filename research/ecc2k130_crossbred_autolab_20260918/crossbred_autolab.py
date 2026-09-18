#!/usr/bin/env python3
"""Local AutoLab control plane for the ECC2K-130 Crossbred α thread."""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PROTOCOL_PATH = HERE / "protocol.json"
EVIDENCE_DIR = HERE / "evidence"
TASK_ID = "TASK-ECC2K130-CROSSBRED-AUTOLAB-20260918"
CARGO_BIN = Path.home() / ".cargo" / "bin"

# Positional columns of the oracle-cost table in examples/crossbred_bench.rs.
# Do not split the header: the |F| cell contains pipes.
CROSSBRED_COLS = (
    "n",
    "m",
    "ell",
    "v",
    "F",
    "Q_enum",
    "Q_word",
    "Q_over_Q_enum",
    "deg",
    "targets",
    "reference",
    "agree",
    "brute",
    "f4",
    "crossbred",
    "xb_brute",
    "xb_f4",
    "D",
    "k",
    "kernel",
    "filters",
    "xb_wall",
)


def runs_dir() -> Path:
    return Path(os.environ.get("CROSSBRED_AUTOLAB_RUNS", str(HERE / "runs")))


def current_path() -> Path:
    return runs_dir() / "current.json"


def lock_path() -> Path:
    return runs_dir() / "autolab.lock"


def load_protocol() -> dict[str, Any]:
    protocol = json.loads(PROTOCOL_PATH.read_text())
    if protocol.get("task_id") != TASK_ID:
        raise SystemExit(f"protocol task_id mismatch: {protocol.get('task_id')}")
    return protocol


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def cargo_env() -> dict[str, str]:
    env = os.environ.copy()
    cargo = str(CARGO_BIN)
    env["PATH"] = cargo + os.pathsep + env.get("PATH", "")
    env["CARGO_TERM_COLOR"] = "never"
    env["CARGO_HOME"] = str(Path.home() / ".cargo")
    return env


def markdown_cells(line: str) -> list[str]:
    return [p.strip() for p in line.strip().strip("|").split("|")]


def parse_dash(value: str) -> Any:
    if value in ("—", "-", "", "–"):
        return None
    lowered = value.lower()
    if lowered in ("yes", "true"):
        return True
    if lowered in ("no", "false"):
        return False
    try:
        if "." in value:
            return float(value)
        return int(value)
    except ValueError:
        return value


def parse_crossbred_output(text: str) -> list[dict[str, Any]]:
    """Parse the oracle-cost markdown table from crossbred_bench."""
    rows: list[dict[str, Any]] = []
    in_table = False
    for line in text.splitlines():
        if line.startswith("| n |"):
            in_table = True
            continue
        if not in_table:
            continue
        if line.startswith("|--"):
            continue
        if line.startswith("|"):
            parts = markdown_cells(line)
            if len(parts) < 8:
                continue
            rec = {
                CROSSBRED_COLS[i]: parse_dash(parts[i])
                for i in range(min(len(CROSSBRED_COLS), len(parts)))
            }
            rec["raw"] = line
            rec["reason"] = parts[-1] if len(parts) > len(CROSSBRED_COLS) - 1 else None
            rows.append(rec)
            continue
        if line.startswith("X1"):
            break
        if in_table and line.strip() == "":
            # Footer follows a blank line. Keep scanning until X1 in case
            # cargo mixed a T4 eprint between rows.
            continue
    return rows


def parse_ffd_table(text: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    in_table = False
    for line in text.splitlines():
        if line.startswith("| n |"):
            in_table = True
            continue
        if in_table and line.startswith("|--"):
            continue
        if in_table and line.startswith("|"):
            parts = markdown_cells(line)
            if len(parts) < 10:
                continue
            no_fall, _, trials = parts[8].partition("/")
            rows.append(
                {
                    "n": parse_dash(parts[0]),
                    "ell": parse_dash(parts[1]),
                    "m": parse_dash(parts[2]),
                    "vars": parse_dash(parts[3]),
                    "eqs": parse_dash(parts[4]),
                    "deg": parse_dash(parts[5]),
                    "ffd_min": parse_dash(parts[6]),
                    "ffd_max": parse_dash(parts[7]),
                    "no_fall": parse_dash(no_fall),
                    "trials": parse_dash(trials),
                    "mean_syz_d2": parse_dash(parts[9]),
                }
            )
            continue
        if in_table and not line.startswith("|"):
            break
    return rows


def usable_rungs(*summaries: dict[str, Any]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for summary in summaries:
        source = summary.get("schema", "")
        for row in summary.get("oracle_rows", []):
            if row.get("usable_for_fit") and row.get("Q_over_Q_enum"):
                tagged = dict(row)
                if "crossbred_x3" in source or row.get("a") is not None:
                    tagged["frame"] = "x3"
                else:
                    tagged["frame"] = "x1"
                out.append(tagged)
    return out


def fit_alpha(rungs: list[dict[str, Any]]) -> dict[str, Any]:
    if len(rungs) < 4:
        return {
            "fit": None,
            "reason": "fewer than four usable rungs",
            "n_rungs": len(rungs),
            "class": "measurement",
        }
    xs = [float(r["ell"]) for r in rungs]
    ys = [math.log2(float(r["Q_over_Q_enum"])) for r in rungs]
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    den = sum((x - mx) ** 2 for x in xs)
    if den == 0:
        return {"fit": None, "reason": "zero variance in ell", "n_rungs": n, "class": "measurement"}
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den
    alpha = slope + 2.0
    if slope <= -0.5:
        verdict = "success"
        klass = "advance"
    elif slope >= -0.1:
        verdict = "falsified"
        klass = "engineering"
    else:
        verdict = "inconclusive"
        klass = "measurement"
    return {
        "fit": True,
        "n_rungs": n,
        "slope": slope,
        "alpha": alpha,
        "verdict": verdict,
        "class": klass,
    }


def load_frozen(protocol: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    x1 = json.loads((REPO / protocol["frozen"]["x1"]).read_text())
    x3 = json.loads((REPO / protocol["frozen"]["x3"]).read_text())
    return x1, x3


def frame_fits(x1: dict[str, Any], x3: dict[str, Any]) -> dict[str, Any]:
    """OLS is per frame. Mixing chained X1 with symmetrised X3 is inadmissible."""
    x1_rungs = usable_rungs(x1)
    x3_rungs = usable_rungs(x3)
    frames = {"x1": fit_alpha(x1_rungs), "x3": fit_alpha(x3_rungs)}
    any_fit = any(f.get("fit") for f in frames.values())
    return {
        "fit": True if any_fit else None,
        "reason": None if any_fit else "no frame has four usable rungs",
        "frames": frames,
        "class": "measurement",
    }


def plan(protocol: dict[str, Any]) -> dict[str, Any]:
    x1, x3 = load_frozen(protocol)
    rungs = usable_rungs(x1, x3)
    report = {
        "task_id": TASK_ID,
        "incumbent": {
            "x1": protocol["frozen"]["x1_verdict"],
            "x3": protocol["frozen"]["x3_verdict"],
            "usable_rungs": [
                {
                    "frame": r.get("frame"),
                    "n": r.get("n"),
                    "ell": r.get("ell"),
                    "Q_over_Q_enum": r.get("Q_over_Q_enum"),
                    "a": r.get("a"),
                }
                for r in rungs
            ],
            "fit": frame_fits(x1, x3),
        },
        "beats": [
            {"beat_id": beat_id, "priority": spec.get("priority"), "label": spec.get("label")}
            for beat_id, spec in protocol["beats"].items()
        ],
        "next": None,
        "next_note": "X5 FFD beats are measured. X4 e2e is not an AutoLab beat.",
        "inadmissible": [
            "Changing T4 after seeing a cell",
            "Reporting a two-point sketch as α",
            "GPU search while filters=0",
            "Quoting wall-clock as S",
        ],
    }
    return report


def pid_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    return True


def clear_stale_lock() -> None:
    path = lock_path()
    if not path.exists():
        return
    try:
        pid = int(path.read_text().strip().split()[0])
    except (ValueError, IndexError, OSError):
        path.unlink(missing_ok=True)
        return
    if not pid_alive(pid):
        path.unlink(missing_ok=True)


class RunLock:
    def __enter__(self) -> "RunLock":
        d = runs_dir()
        d.mkdir(parents=True, exist_ok=True)
        clear_stale_lock()
        self.path = lock_path()
        try:
            self.fd = os.open(self.path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError as exc:
            raise SystemExit(f"autolab locked: {self.path}") from exc
        os.write(self.fd, f"{os.getpid()}\n".encode())
        return self

    def __exit__(self, *exc: Any) -> None:
        os.close(self.fd)
        try:
            self.path.unlink()
        except FileNotFoundError:
            pass


def preflight() -> dict[str, Any]:
    clear_stale_lock()
    env = cargo_env()
    rustc = subprocess.run(["rustc", "--version"], capture_output=True, text=True, env=env)
    cargo = subprocess.run(["cargo", "--version"], capture_output=True, text=True, env=env)
    bench = REPO / "examples/crossbred_bench.rs"
    ffd = REPO / "examples/ffd_probe.rs"
    x1 = REPO / "experiments/ecc2k130_crossbred_x1_20260918/summary.json"
    x3 = REPO / "experiments/ecc2k130_crossbred_x3_20260918/summary.json"
    ok = (
        rustc.returncode == 0
        and cargo.returncode == 0
        and bench.is_file()
        and ffd.is_file()
        and x1.is_file()
        and x3.is_file()
    )
    report = {
        "ok": ok,
        "rustc": rustc.stdout.strip() or rustc.stderr.strip(),
        "cargo": cargo.stdout.strip() or cargo.stderr.strip(),
        "host": platform.node(),
        "crossbred_bench": str(bench.relative_to(REPO)),
        "ffd_probe": str(ffd.relative_to(REPO)),
        "frozen_x1": str(x1.relative_to(REPO)),
        "frozen_x3": str(x3.relative_to(REPO)),
    }
    if not ok:
        raise SystemExit(json.dumps(report, indent=2))
    return report


def run_command(command: list[str], extra_env: dict[str, str] | None, log: Path) -> str:
    env = cargo_env()
    if extra_env:
        env.update(extra_env)
    with log.open("w") as fh:
        proc = subprocess.run(
            command,
            cwd=REPO,
            env=env,
            stdout=fh,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    text = log.read_text(errors="replace")
    if proc.returncode != 0:
        raise SystemExit(f"command failed ({proc.returncode}): {' '.join(command)}")
    return text


def claim_smoke(beat: dict[str, Any], text: str) -> dict[str, Any]:
    rows = parse_crossbred_output(text)
    expect = beat["expect"]
    match = next((r for r in rows if r.get("n") == expect["n"] and r.get("m") == expect["m"]), None)
    if match is None:
        return {"status": "FAIL", "reason": "expected row missing", "rows": rows}
    q = match.get("Q_over_Q_enum")
    ok_agree = match.get("agree") is True
    ok_q = isinstance(q, (int, float)) and abs(q - expect["Q_over_Q_enum"]) <= expect["rel_tol"] * expect[
        "Q_over_Q_enum"
    ]
    status = "PASS" if ok_agree and ok_q else "FAIL"
    return {
        "status": status,
        "row": {
            k: match.get(k)
            for k in ("n", "m", "ell", "v", "F", "Q_enum", "Q_word", "Q_over_Q_enum", "agree", "filters")
        },
        "expect": expect,
        "class": "accounting",
    }


def claim_fit(protocol: dict[str, Any]) -> dict[str, Any]:
    x1, x3 = load_frozen(protocol)
    rungs = usable_rungs(x1, x3)
    fit = frame_fits(x1, x3)
    expect_null = protocol["beats"]["fit.alpha"]["expect"]["fit"] is None
    status = "PASS" if (fit.get("fit") is None) == expect_null else "FAIL"
    return {
        "status": status,
        "rungs": [
            {
                "frame": r.get("frame"),
                "n": r.get("n"),
                "ell": r.get("ell"),
                "Q_over_Q_enum": r.get("Q_over_Q_enum"),
                "a": r.get("a"),
            }
            for r in rungs
        ],
        "fit": fit,
        "class": "measurement",
        "note": "Frozen X1+X3 only, OLS per frame. Mixing frames is inadmissible. A two-point sketch is not a fit.",
    }


def claim_ffd(text: str, beat: dict[str, Any]) -> dict[str, Any]:
    rows = parse_ffd_table(text)
    m4 = [r for r in rows if r.get("m") == 4]
    ell1 = [r for r in m4 if r.get("ell") == 1]
    comparable = [r for r in m4 if isinstance(r.get("ell"), int) and r["ell"] > 1]

    def maxima_and_growth(group: list[dict[str, Any]]) -> tuple[list[Any], bool]:
        growing = False
        maxima: list[Any] = []
        if len(group) >= 2:
            ordered = sorted(group, key=lambda r: r["n"] or 0)
            maxima = [r.get("ffd_max") for r in ordered]
            if all(isinstance(v, int) for v in maxima):
                growing = maxima[-1] > maxima[0]
        elif group:
            maxima = [group[0].get("ffd_max")]
        return maxima, growing

    maxima_all, growing_all = maxima_and_growth(m4)
    maxima_cmp, growing_cmp = maxima_and_growth(comparable)
    # ell = 1 has no FB bits, so S3 collapses (degree 1). That row is
    # not an m=4 Semaev system; H1 is read from ell > 1. Mixing it into
    # the growth flag is accounting, not a falsification of H1.
    return {
        "status": "PASS" if rows else "FAIL",
        "rows": rows,
        "m4_ffd_grows_with_n": growing_cmp,
        "m4_ffd_maxima": maxima_cmp,
        "m4_ffd_grows_with_n_including_ell_1": growing_all,
        "m4_ffd_maxima_including_ell_1": maxima_all,
        "m4_ell_eq_1": ell1,
        "falsifier": beat["expect"]["falsifier"],
        "class": "measurement",
        "note": beat.get("note"),
    }


def launch(beat_id: str) -> dict[str, Any]:
    protocol = load_protocol()
    if beat_id not in protocol["beats"]:
        raise SystemExit(f"unknown beat: {beat_id}")
    beat = protocol["beats"][beat_id]
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_id = f"{beat_id.replace('.', '-')}-{stamp}"
    run_dir = runs_dir() / run_id
    with RunLock():
        run_dir.mkdir(parents=True)
        (run_dir / "logs").mkdir()
        (run_dir / "artifacts").mkdir()
        write_json(run_dir / "inputs" / "protocol.json", protocol)
        started = dt.datetime.now(dt.timezone.utc).isoformat()
        if beat_id == "fit.alpha":
            claim = claim_fit(protocol)
            text = ""
        else:
            command = list(beat["command"])
            text = run_command(command, beat.get("env"), run_dir / "logs" / "bench.txt")
            if beat_id.startswith("x5."):
                claim = claim_ffd(text, beat)
            else:
                claim = claim_smoke(beat, text)
        claim.update(
            {
                "task_id": TASK_ID,
                "beat_id": beat_id,
                "run_id": run_id,
                "started_utc": started,
                "finished_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                "host": platform.node(),
                "log_sha256": sha256_file(run_dir / "logs" / "bench.txt") if text else None,
            }
        )
        write_json(run_dir / "artifacts" / "claim.json", claim)
        state = {"status": claim["status"], "beat_id": beat_id, "run_id": run_id, "dir": str(run_dir)}
        write_json(run_dir / "state.json", state)
        write_json(current_path(), state)
    return claim


def status() -> dict[str, Any]:
    path = current_path()
    if not path.exists():
        return {"status": "idle"}
    return json.loads(path.read_text())


def table_excerpt(text: str) -> str:
    lines = text.splitlines()
    start = next((i for i, line in enumerate(lines) if line.startswith("| n |")), None)
    if start is None:
        return text[-4000:]
    return "\n".join(lines[start:]) + "\n"


def verify() -> dict[str, Any]:
    st = status()
    if st.get("status") == "idle":
        raise SystemExit("no current run")
    run_dir = Path(st["dir"])
    claim_path = run_dir / "artifacts" / "claim.json"
    claim = json.loads(claim_path.read_text())
    log = run_dir / "logs" / "bench.txt"
    expected = claim.get("log_sha256")
    if expected:
        actual = sha256_file(log)
        if actual != expected:
            raise SystemExit(f"log hash mismatch: {actual} != {expected}")
    try:
        claim_rel = str(claim_path.relative_to(REPO))
    except ValueError:
        claim_rel = str(claim_path)
    return {"status": "PASS", "claim": claim_rel, "log_sha256": expected}


def _rel(path: Path) -> str:
    try:
        return str(path.relative_to(REPO))
    except ValueError:
        return str(path)


def promote_run(run_dir: Path) -> dict[str, Any]:
    state = json.loads((run_dir / "state.json").read_text())
    if state.get("status") != "PASS":
        raise SystemExit(f"refusing to promote status={state.get('status')}")
    dest = EVIDENCE_DIR / state["run_id"]
    dest.mkdir(parents=True, exist_ok=True)
    shutil.copy2(run_dir / "artifacts" / "claim.json", dest / "claim.json")
    shutil.copy2(run_dir / "state.json", dest / "state.json")
    log = run_dir / "logs" / "bench.txt"
    if log.exists():
        (dest / "logs").mkdir(exist_ok=True)
        (dest / "logs" / "table.txt").write_text(table_excerpt(log.read_text(errors="replace")))
        (dest / "logs" / "log_sha256.txt").write_text(sha256_file(log) + "\n")
    return {"status": "PASS", "evidence": _rel(dest), "run_id": state["run_id"]}


def promote(all_runs: bool = False) -> dict[str, Any]:
    if all_runs:
        promoted = []
        for state_path in sorted(runs_dir().glob("*/state.json")):
            state = json.loads(state_path.read_text())
            if state.get("status") == "PASS":
                promoted.append(promote_run(state_path.parent))
        if not promoted:
            raise SystemExit("no PASS runs to promote")
        return {"status": "PASS", "promoted": promoted}
    st = status()
    if st.get("status") == "idle":
        raise SystemExit("no current run")
    if st.get("status") != "PASS":
        raise SystemExit(f"refusing to promote status={st.get('status')}")
    return promote_run(Path(st["dir"]))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser("plan")
    sub.add_parser("preflight")
    launch_p = sub.add_parser("launch")
    launch_p.add_argument("--beat", required=True)
    sub.add_parser("status")
    sub.add_parser("verify")
    promote_p = sub.add_parser("promote")
    promote_p.add_argument("--all", action="store_true")
    args = parser.parse_args()
    if args.cmd == "plan":
        print(json.dumps(plan(load_protocol()), indent=2))
    elif args.cmd == "preflight":
        print(json.dumps(preflight(), indent=2))
    elif args.cmd == "launch":
        print(json.dumps(launch(args.beat), indent=2))
    elif args.cmd == "status":
        print(json.dumps(status(), indent=2))
    elif args.cmd == "verify":
        print(json.dumps(verify(), indent=2))
    elif args.cmd == "promote":
        print(json.dumps(promote(all_runs=args.all), indent=2))


if __name__ == "__main__":
    main()

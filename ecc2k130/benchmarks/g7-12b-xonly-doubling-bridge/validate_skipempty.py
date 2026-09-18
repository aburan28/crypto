#!/usr/bin/env python3
"""Resource, builtin, checkpoint and independent affine-x gates for skip-empty."""
from pathlib import Path
import json, re, subprocess
from experiment import ROOT, OUT, sha

CONTROL = ROOT / "build/g7-xonly-mod72-control"
CANDIDATE = ROOT / "build/g7-xonly-mod72-skipempty"
ORACLE_SRC = OUT / "check_bridge13_mod72_state.cpp"
ORACLE = ROOT / "build/check-bridge13-mod72-state"
WALK = re.compile(
    r"packed kernel: (\d+) registers/thread, (\d+) local bytes/thread, "
    r"(\d+) shared bytes/block"
)


def run(command, log, timeout=900, cwd=ROOT):
    log = Path(log)
    assert not log.exists(), log
    started = time.monotonic()
    with log.open("x") as stream:
        result = subprocess.run(command, cwd=cwd, stdout=stream, stderr=subprocess.STDOUT, timeout=timeout)
    return {
        "command": [str(part) for part in command],
        "returncode": result.returncode,
        "elapsed_seconds": time.monotonic() - started,
        "log": str(log.relative_to(OUT)),
        "log_sha256": sha(log),
    }


def resources(text):
    match = WALK.search(text)
    assert match, text[-2000:]
    registers, local, shared = map(int, match.groups())
    occupancy = "50.0% (3 blocks/SM, 256 threads/block" in text
    return {
        "registers": registers,
        "local_bytes": local,
        "shared_bytes": shared,
        "three_block_occupancy": occupancy,
        "passed": registers <= 80 and local == 0 and shared <= 32400 and occupancy,
    }


def main():
    rows = {}
    for label, binary in ("control", CONTROL), ("candidate", CANDIDATE):
        assert binary.exists(), binary
        test = run([str(binary), "--test"], OUT / f"skipempty-{label}-builtin-test.log")
        test_text = (OUT / test["log"]).read_text()
        test["passed"] = test["returncode"] == 0 and "FAILED" not in test_text
        ckpt = ROOT / "build" / f"g7-xonly-mod72-{label}-state.ckpt"
        assert not ckpt.exists(), ckpt
        state = run(
            [str(binary), "--packed", "--curve", "131", "--threads", "128",
             "--steps", "7", "--launches", "1", "--verify", "0",
             "--run-id", "62580", "--checkpoint", str(ckpt)],
            OUT / f"skipempty-{label}-state-client.log",
        )
        state_text = (OUT / state["log"]).read_text()
        state["resources"] = resources(state_text)
        state["checkpoint"] = str(ckpt.relative_to(ROOT))
        state["checkpoint_sha256"] = sha(ckpt)
        state["passed"] = (
            state["returncode"] == 0 and "MISMATCH" not in state_text
            and "14336 iterations" in state_text and "0 dropped" in state_text
            and state["resources"]["passed"]
        )
        rows[label] = {"binary_sha256": sha(binary), "test": test, "state": state}

    if not ORACLE.exists():
        compile_oracle = run(
            ["g++", "-O3", "-std=c++17", "-I", str(ROOT.parent), str(ORACLE_SRC), "-o", str(ORACLE)],
            OUT / "skipempty-oracle-compile.log",
            cwd=ROOT.parent,
        )
        assert compile_oracle["returncode"] == 0, compile_oracle
        rows["oracle_compile"] = compile_oracle

    for label in ("control", "candidate"):
        ckpt = ROOT / "build" / f"g7-xonly-mod72-{label}-state.ckpt"
        oracle = run([str(ORACLE), str(ckpt)], OUT / f"skipempty-{label}-state-oracle.log", cwd=ROOT.parent)
        text = (OUT / oracle["log"]).read_text().strip()
        oracle["output"] = text
        oracle["passed"] = oracle["returncode"] == 0 and text.endswith("mismatches=0")
        rows[label]["oracle"] = oracle

    rows["identical_checkpoints"] = (
        rows["control"]["state"]["checkpoint_sha256"] == rows["candidate"]["state"]["checkpoint_sha256"]
    )
    rows["passed"] = (
        rows["control"]["test"]["passed"] and rows["candidate"]["test"]["passed"]
        and rows["control"]["state"]["passed"] and rows["candidate"]["state"]["passed"]
        and rows["control"]["oracle"]["passed"] and rows["candidate"]["oracle"]["passed"]
        and rows["identical_checkpoints"]
    )
    (OUT / "skipempty-validation.json").write_text(json.dumps(rows, indent=2) + "\n")
    print(json.dumps(rows, indent=2), flush=True)
    raise SystemExit(0 if rows["passed"] else 1)


if __name__ == "__main__":
    main()

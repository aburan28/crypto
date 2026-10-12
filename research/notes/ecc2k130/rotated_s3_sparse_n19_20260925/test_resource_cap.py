#!/usr/bin/env python3
"""Harmless Linux controls for the external wall and process-group RSS gates."""
from __future__ import annotations

import json
import sys
import time
import tempfile
from pathlib import Path

import run


def expect_stop(root: Path, label: str, child_code: str, cap_seconds: float,
                cap_bytes: int, expected_stop: str) -> dict:
    receipt = {"attempts": []}
    expected = root / (label + ".never")
    try:
        run.child(receipt, root, label, [sys.executable, "-c", child_code],
                  expected, cap_seconds, cap_bytes)
    except RuntimeError:
        pass
    else:
        raise AssertionError(f"{label}: bounded child was admitted")
    attempt, = receipt["attempts"]
    assert attempt["phase"] == label
    assert attempt["resource_stop"] == expected_stop, attempt
    assert attempt["rss_cap_bytes"] == cap_bytes
    assert attempt["group_quiesced"] is True
    assert attempt["stdout_sha256"] == run.sha(root / f"{label}.stdout.txt")
    assert attempt["stderr_sha256"] == run.sha(root / f"{label}.stderr.txt")
    if expected_stop == "PROCESS_GROUP_RSS_CAP":
        assert attempt["sampled_group_peak_rss_bytes"] >= cap_bytes
    else:
        assert attempt["external_timeout"]
    return {"phase": label, "stop": attempt["resource_stop"],
            "observed_peak_bytes": attempt["sampled_group_peak_rss_bytes"]}


def main() -> None:
    assert sys.platform == "linux"
    with tempfile.TemporaryDirectory(prefix="n19-resource-gate-") as name:
        root = Path(name)
        wall = expect_stop(root, "wall", "import time; time.sleep(5)",
                           0.2, 128 * 1024 * 1024, "EXTERNAL_WALL_CAP")
        # Touch each page so the allocation is reflected in resident memory.
        allocation = ("import time; x=bytearray(64*1024*1024); "
                      "x[::4096]=b'X'*len(x[::4096]); time.sleep(5)")
        memory = expect_stop(root, "memory", allocation, 5,
                             32 * 1024 * 1024, "PROCESS_GROUP_RSS_CAP")
        descendant = ("import subprocess,sys,time; "
                      "subprocess.Popen([sys.executable,'-c',"
                      "'import time; x=bytearray(64*1024*1024); "
                      "x[::4096]=b\\'X\\'*len(x[::4096]); time.sleep(5)']); "
                      "time.sleep(5)")
        descendant_row = expect_stop(root, "descendant", descendant, 5,
                                     32 * 1024 * 1024, "PROCESS_GROUP_RSS_CAP")
        assert descendant_row["observed_peak_bytes"] >= 32 * 1024 * 1024
        lingering = ("import subprocess,sys; "
                     "subprocess.Popen([sys.executable,'-c',"
                     "'import sys,time; [(sys.stdout.write(\\'x\\'), "
                     "sys.stdout.flush(), time.sleep(.02)) for _ in range(250)]'])")
        receipt = {"attempts": []}
        try:
            run.child(receipt, root, "lingering", [sys.executable, "-c", lingering],
                      root / "lingering.never", 5, 128 * 1024 * 1024)
        except RuntimeError:
            pass
        else:
            raise AssertionError("lingering child unexpectedly admitted")
        lingering_attempt, = receipt["attempts"]
        assert lingering_attempt["resource_stop"] is None, lingering_attempt
        assert lingering_attempt["exit_code"] == 0
        assert lingering_attempt["group_quiesced"] is True
        recorded_stdout = lingering_attempt["stdout_sha256"]
        time.sleep(.2)
        assert run.sha(root / "lingering.stdout.txt") == recorded_stdout
    print(json.dumps({"decision": "PASS_RESOURCE_CAP_CONTROLS",
                      "controls": [wall, memory, descendant_row,
                                   {"phase": "lingering", "group_quiesced": True}]},
                     sort_keys=True))


if __name__ == "__main__":
    main()

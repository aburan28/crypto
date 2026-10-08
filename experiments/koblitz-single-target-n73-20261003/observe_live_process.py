#!/usr/bin/env python3
"""Sample the active IC child RSS with macOS proc_pidinfo; no solver inputs."""
import ctypes
import ctypes.util
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path


HERE = Path(__file__).resolve().parent
protocol = json.loads((HERE / "protocol.json").read_text())
run_dir = HERE / "runs" / protocol["run_id"]
run_record = json.loads((run_dir / "run.json").read_text())
pid = int(json.loads((run_dir / "ic_execution.json").read_text())["pid"])


class TaskInfo(ctypes.Structure):
    _fields_ = [
        ("virtual_size", ctypes.c_uint64),
        ("resident_size", ctypes.c_uint64),
        ("total_user", ctypes.c_uint64),
        ("total_system", ctypes.c_uint64),
        ("threads_user", ctypes.c_uint64),
        ("threads_system", ctypes.c_uint64),
        ("policy", ctypes.c_int32),
        ("faults", ctypes.c_int32),
        ("pageins", ctypes.c_int32),
        ("cow_faults", ctypes.c_int32),
        ("messages_sent", ctypes.c_int32),
        ("messages_received", ctypes.c_int32),
        ("syscalls_mach", ctypes.c_int32),
        ("syscalls_unix", ctypes.c_int32),
        ("csw", ctypes.c_int32),
        ("threadnum", ctypes.c_int32),
        ("numrunning", ctypes.c_int32),
        ("priority", ctypes.c_int32),
    ]


libproc = ctypes.CDLL(ctypes.util.find_library("proc") or "/usr/lib/libproc.dylib", use_errno=True)
proc_pidinfo = libproc.proc_pidinfo
proc_pidinfo.argtypes = [ctypes.c_int32, ctypes.c_int32, ctypes.c_uint64, ctypes.c_void_p, ctypes.c_int32]
proc_pidinfo.restype = ctypes.c_int32
output = run_dir / "live_resource_observations.jsonl"

with output.open("a") as stream:
    while True:
        task = TaskInfo()
        received = proc_pidinfo(pid, 4, 0, ctypes.byref(task), ctypes.sizeof(task))
        if received != ctypes.sizeof(task):
            break
        row = {
            "kind": "live_single_target_process_resource_observation",
            "observed_at_utc": datetime.now(timezone.utc).isoformat(),
            "pid": pid,
            "method": "macOS proc_pidinfo(PROC_PIDTASKINFO); exact-process point sample, not peak",
            "resident_rss_bytes": task.resident_size,
            "virtual_size_bytes": task.virtual_size,
            "running_threads": task.numrunning,
        }
        stream.write(json.dumps(row, sort_keys=True) + "\n")
        stream.flush()
        os.fsync(stream.fileno())
        time.sleep(30)

#!/usr/bin/env python3
"""Emit reproducible LOCAL evidence; never infer production certification.

The live GPU, cloud failover and campaign-scale gates are explicitly unrun.
--require-production deliberately exits 2 until a separate live review exists.
"""
import argparse
import datetime
import hashlib
from pathlib import Path
import platform
import subprocess
import sys
import time

from protocol import atomicJson, sha256File

ROOT = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", required=True)
    ap.add_argument("--require-production", action="store_true")
    args = ap.parse_args()
    output = Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    logs = output.parent / (output.stem + "-logs")
    logs.mkdir(exist_ok=True)
    files = [p for directory in ("include", "src", "generated", "aws")
             for p in (ROOT / directory).rglob("*")
             if p.is_file() and p.suffix in (".h", ".cuh", ".cu", ".cpp", ".py", ".sh")]
    files += [ROOT / "Makefile", ROOT / "aws/campaign.json"]
    hashes = {str(p.relative_to(ROOT)): sha256File(p) for p in sorted(files)}
    h = hashlib.sha256()
    for path, digest in sorted(hashes.items()):
        h.update((path + "\0" + digest + "\n").encode())
    report = {
        "schema": 1, "scope": "local-software-readiness", "classification": "engineering-correctness",
        "startedAt": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "gitCommit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "worktreeDiffersFromCommit": bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT)),
        "platform": platform.platform(), "python": platform.python_version(),
        "sourceDigest": h.hexdigest(), "files": hashes,
        "clientSha256": sha256File(ROOT / "ecc2k130-cpu"),
        "fixtureBinarySha256": sha256File(ROOT / "build/test-production"),
        "commands": [], "localPassed": False, "productionCertified": False,
        "externalGates": [
            {"gate": "exact-production-CUDA-build-and-device-differential-tests", "status": "NOT_RUN"},
            {"gate": "full-length-DP34-replay-and-multidevice-checkpoint-resume", "status": "NOT_RUN"},
            {"gate": "live-S3-lease-failover-network-partition-and-power-loss", "status": "NOT_RUN"},
            {"gate": "campaign-scale-external-merge-capacity-and-replay-latency", "status": "NOT_RUN"},
            {"gate": "RDS-seed-only-adapter-and-collision-resolver", "status": "NOT_CERTIFIED_NOT_IN_MAIN"},
        ],
        "notes": ["Raw DP payload remains 32 bytes; manifests supply version, campaign binding and SHA-256.",
                  "Local directory tests do not establish live S3 behavior or GPU correctness.",
                  "Passing software checks is not an ECC2K-130 discrete-log solution or a speedup claim."],
    }
    commands = [
        ("reference-arithmetic", ["./ecc2k130-cpu", "--test"], 180),
        ("persistence-faults", ["./build/test-production"], 90),
        ("storage-and-end-to-end", [sys.executable, "-m", "unittest", "discover", "-s", "aws", "-p", "test_certification.py", "-v"], 240),
        ("packed-host-and-timing", ["make", "test-packed-network", "test-timing"], 240),
    ]
    for name, cmd, limit in commands:
        started = time.monotonic()
        log = logs / (name + ".log")
        with log.open("wb") as fh:
            try:
                proc = subprocess.run(cmd, cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT, timeout=limit)
                rc = proc.returncode
            except subprocess.TimeoutExpired:
                fh.write(b"\nCERTIFICATION TIMEOUT\n")
                rc = 124
        report["commands"].append({"name": name, "argv": cmd, "returncode": rc,
                                   "seconds": time.monotonic() - started,
                                   "log": str(log.relative_to(output.parent)), "logSha256": sha256File(log)})
        atomicJson(output, report)
        print(name + ": " + ("PASS" if rc == 0 else "FAIL"), flush=True)
    report["localPassed"] = all(c["returncode"] == 0 for c in report["commands"])
    report["completedAt"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    atomicJson(output, report)
    print("Local software: %s; production: NOT CERTIFIED" % ("PASS" if report["localPassed"] else "FAIL"))
    return 1 if not report["localPassed"] else 2 if args.require_production else 0


if __name__ == "__main__":
    raise SystemExit(main())

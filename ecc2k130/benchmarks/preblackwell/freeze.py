#!/usr/bin/env python3
"""Compact a Modal bench JSON into the preblackwell receipt shape."""
from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

REF_6000_B = 15.115792


def compact(src: dict, family: str, sms: int, prior_b: float, prior_kind: str) -> dict:
    samples = []
    for s in src["samples"]:
        raw = s.get("raw", "")
        km = re.search(r"packed kernel: (\d+) registers/thread, (\d+) local bytes/thread", raw)
        dm = re.search(r"device: ([^\n]+)", raw)
        wm = re.search(r"backend cuda-packed131: (\d+) threads", raw)
        samples.append({
            "valid": s.get("valid"),
            "rateM": s.get("rate"),
            "packedClmad": s.get("packedClmad"),
            "expectedPackedClmad": s.get("expectedPackedClmad"),
            "registers": int(km.group(1)) if km else None,
            "localBytes": int(km.group(2)) if km else None,
            "deviceLine": dm.group(1) if dm else None,
            "threads": int(wm.group(1)) if wm else None,
        })
    ident = src.get("identity", {})
    gpu_state = ident.get("gpuState", "")
    um = re.search(r"GPU-[0-9a-f-]+", gpu_state)
    med = src["rate"]
    return {
        "kind": "preblackwell complete-scalar-iteration throughput",
        "family": family,
        "gpu": src["gpu"],
        "cc": src["cc"],
        "valid": src["valid"],
        "packedClmad": src["packedClmad"],
        "batch": src["batch"],
        "blockThreads": src["threads"],
        "minBlocks": src["minBlocks"],
        "workers": "automatic",
        "automaticThreads": samples[0].get("threads") if samples else None,
        "steps": src["steps"],
        "launches": src["launches"],
        "repeats": src["repeats"],
        "sms": sms,
        "sourceSha256": ident.get("sourceSha256"),
        "binarySha256": ident.get("binarySha256"),
        "activeBackend": ident.get("activeBackend"),
        "cudaImageVersion": ident.get("cudaImageVersion"),
        "gpuUuid": um.group(0) if um else None,
        "gpuState": gpu_state.strip(),
        "compiler": "nvcc 13.3.73",
        "ratesM": [s["rate"] for s in src["samples"]],
        "minRateM": src["minRate"],
        "maxRateM": src["maxRate"],
        "medianM": med,
        "medianB": med / 1000.0,
        "prior": {"kind": prior_kind, "B": prior_b},
        "ratioToPrior": (med / 1000.0) / prior_b if prior_b else None,
        "ratioTo6000": (med / 1000.0) / REF_6000_B,
        "perSmM": med / sms,
        "samples": samples,
    }


def from_log(text: str) -> dict | None:
    decoder = json.JSONDecoder()
    found = None
    i = 0
    while i < len(text):
        start = text.find("{", i)
        if start < 0:
            break
        try:
            obj, end = decoder.raw_decode(text, start)
        except json.JSONDecodeError:
            i = start + 1
            continue
        if isinstance(obj, dict) and obj.get("samples") and "rate" in obj:
            found = obj
        i = end
    return found


def main() -> None:
    if len(sys.argv) != 7:
        raise SystemExit(
            "usage: freeze.py SRCJSON_OR_LOG DESTJSON family sms priorB priorKind"
        )
    src_path, dest, family, sms, prior_b, prior_kind = (
        Path(sys.argv[1]),
        Path(sys.argv[2]),
        sys.argv[3],
        int(sys.argv[4]),
        float(sys.argv[5]),
        sys.argv[6],
    )
    raw = src_path.read_text(encoding="utf-8", errors="replace")
    try:
        src = json.loads(raw)
    except json.JSONDecodeError:
        src = from_log(raw)
        if src is None:
            raise SystemExit("no bench JSON in " + str(src_path))
    obj = compact(src, family, sms, prior_b, prior_kind)
    dest.write_text(json.dumps(obj, indent=2) + "\n")
    print(
        dest,
        "medianB=",
        obj["medianB"],
        "clmad=",
        obj["packedClmad"],
        "sha=",
        hashlib.sha256(dest.read_bytes()).hexdigest()[:12],
    )


if __name__ == "__main__":
    main()

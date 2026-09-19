#!/usr/bin/env python3
"""Turn a Modal ::bench log into the slim receipt B200.md cites."""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REF_B = 15.115792
PRIOR_SMS = 148
PRIOR_B = REF_B * PRIOR_SMS / 188.0


def last_json_object(text: str) -> dict:
    decoder = json.JSONDecoder()
    last = None
    i = 0
    while True:
        j = text.find("{", i)
        if j < 0:
            break
        try:
            obj, end = decoder.raw_decode(text, j)
        except json.JSONDecodeError:
            i = j + 1
            continue
        if isinstance(obj, dict):
            last = obj
        i = max(end, j + 1)
    if not isinstance(last, dict):
        raise SystemExit("no JSON object in log")
    return last


def first_match(pattern: str, text: str, flags=0) -> str | None:
    m = re.search(pattern, text, flags)
    return m.group(1) if m else None


def sample_view(sample: dict) -> dict:
    raw = sample.get("raw") or ""
    sms = first_match(r"(\d+) SMs", raw)
    threads = first_match(r"backend cuda-packed131: (\d+) threads", raw)
    if threads is None:
        threads = first_match(r", (\d+) packed threads resident", raw)
    return {
        "valid": bool(sample.get("valid")),
        "rateM": sample.get("rate"),
        "packedClmad": sample.get("packedClmad"),
        "expectedPackedClmad": sample.get("expectedPackedClmad"),
        "registers": int(first_match(r"packed kernel: (\d+) registers", raw) or 0) or None,
        "localBytes": int(first_match(r"(\d+) local bytes/thread", raw) or 0)
        if "local bytes/thread" in raw else None,
        "deviceLine": first_match(r"^(NVIDIA .+)$", raw, re.M),
        "threads": int(threads) if threads else None,
        "sms": int(sms) if sms else None,
        "error": sample.get("error"),
    }


def freeze(raw: dict, arm: str) -> dict:
    gpu = raw.get("gpu") or ""
    if "B200" not in gpu:
        raise SystemExit("refusing to freeze a non-B200 gpu name: %r" % gpu)
    cc = str(raw.get("cc") or "")
    if cc not in ("100", "10.0"):
        raise SystemExit("refusing to freeze cc=%r (want sm_100)" % cc)
    if not raw.get("valid"):
        raise SystemExit("benchmark invalid: %s" % raw.get("error"))
    samples = [sample_view(s) for s in (raw.get("samples") or [])]
    if len(samples) != 3 or not all(s["valid"] for s in samples):
        raise SystemExit("need three valid repeats, got %r" % samples)
    rates_m = [s["rateM"] for s in samples]
    median_m = raw["rate"]
    median_b = median_m / 1000.0
    sms = next((s["sms"] for s in samples if s["sms"]), None)
    identity = raw.get("identity") or {}
    packed = bool(raw.get("packedClmad"))
    return {
        "kind": "b200 complete-scalar-iteration throughput",
        "family": "b200",
        "gpu": gpu,
        "cc": "100",
        "valid": True,
        "packedClmad": packed,
        "arm": arm,
        "batch": raw.get("batch"),
        "blockThreads": raw.get("threads"),
        "minBlocks": raw.get("minBlocks"),
        "workers": "automatic",
        "requestedWorkers": raw.get("workers"),
        "automaticThreads": next((s["threads"] for s in samples if s["threads"]), None),
        "steps": raw.get("steps"),
        "launches": raw.get("launches"),
        "repeats": raw.get("repeats"),
        "sms": sms,
        "sourceSha256": identity.get("sourceSha256"),
        "binarySha256": identity.get("binarySha256"),
        "activeBackend": identity.get("activeBackend"),
        "cudaImageVersion": identity.get("cudaImageVersion"),
        "gpuState": identity.get("gpuState"),
        "compiler": identity.get("compiler"),
        "ratesM": rates_m,
        "minRateM": raw.get("minRate"),
        "maxRateM": raw.get("maxRate"),
        "medianM": median_m,
        "medianB": median_b,
        "prior": {"kind": "sm-count-scale", "sms": PRIOR_SMS, "B": PRIOR_B},
        "ratioToPrior": median_b / PRIOR_B,
        "ratioTo6000": median_b / REF_B,
        "perSmM": (median_m / sms) if sms else None,
        "samples": samples,
    }


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit("usage: freeze.py ARM LOG.json-or-log OUT.json")
    arm, src, dst = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
    text = src.read_text(errors="replace")
    try:
        raw = json.loads(text)
    except json.JSONDecodeError:
        raw = last_json_object(text)
    dst.write_text(json.dumps(freeze(raw, arm), indent=2) + "\n")
    print("wrote", dst)


if __name__ == "__main__":
    main()

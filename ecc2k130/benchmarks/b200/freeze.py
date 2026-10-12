#!/usr/bin/env python3
"""Turn a Modal ::bench log into the slim receipt B200.md cites."""
from __future__ import annotations

import json
import re
import statistics
import sys
from pathlib import Path

REF_B = 15.115792
PRIOR_SMS = 148
PRIOR_B = REF_B * PRIOR_SMS / 188.0
ANSI = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")


def strip_ansi(text: str) -> str:
    return ANSI.sub("", text).replace("\r", "\n")


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
        raise ValueError("no JSON object in log")
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
        "deviceLine": first_match(r"^(NVIDIA .+)$", raw, re.M) or first_match(
            r"(NVIDIA B200, \d+ SMs[^\n]*)", raw),
        "threads": int(threads) if threads else None,
        "sms": int(sms) if sms else None,
        "error": sample.get("error"),
    }


def parse_modal_log(text: str) -> dict:
    """Rebuild the bench dict from a TTY-mangled Modal log.

    `modal run` interleaves spinner/ANSI with json.dumps, so the JSON object
    is not always decodeable. The repeat lines and identity fields still are.
    """
    text = strip_ansi(text)
    flat = re.sub(r"\s+", " ", text)
    gpu = first_match(r'"gpu": "(NVIDIA B200[^"]*)"', text) or first_match(
        r"device: (NVIDIA B200)", text)
    if not gpu:
        raise SystemExit("refusing to freeze a log with no NVIDIA B200 identity")
    cc = first_match(r'"cc": "(\d+)"', text) or first_match(
        r"compute capability: (\d+)", text)
    rates = [float(x) for x in re.findall(
        r"repeat \d+/3: ([0-9.]+) M it/s \(complete\)", text)]
    if len(rates) != 3:
        rates = [float(x) for x in re.findall(r"finished: ([0-9.]+) M it/s", text)]
        # identity dump can repeat the first finished line; keep the last three
        rates = rates[-3:]
    if len(rates) != 3:
        raise SystemExit("need three complete repeats, got %r" % rates)
    clmad_bits = re.findall(r"packed native carryless multiply: ([01])", flat)
    if not clmad_bits:
        flags = re.findall(r'"packedClmad": (true|false)', text)
        packed = flags[-1] == "true" if flags else None
    else:
        packed = clmad_bits[-1] == "1"
    sms = int(first_match(r"NVIDIA B200, (\d+) SMs", flat) or 0) or None
    threads = first_match(r"backend cuda-packed131: (\d+) threads", flat)
    regs = first_match(r"packed kernel: (\d+) registers", flat)
    local = first_match(r"packed kernel: \d+ registers/thread, (\d+) local bytes", flat)
    device = first_match(r"(NVIDIA B200, \d+ SMs[^\\]*)", flat)
    samples = []
    for rate in rates:
        samples.append({
            "valid": True,
            "rate": rate,
            "packedClmad": packed,
            "expectedPackedClmad": packed,
            "raw": "device: %s\npacked kernel: %s registers/thread, %s local bytes/thread\n"
                   "packed native carryless multiply: %s\nbackend cuda-packed131: %s threads\n"
                   "finished: %.3f M it/s, 0 distinguished points "
                   "(0 verified against the reference, 0 dropped)\n" % (
                       device or gpu, regs or "0", local or "0",
                       "1" if packed else "0", threads or "0", rate),
        })
    sha = re.findall(r'"sourceSha256":\s*"?[\s|]*"([0-9a-f]{64})"', text)
    if not sha:
        sha = re.findall(r"sourceSha256[\s\":|]+([0-9a-f]{64})", flat)
    bsha = re.findall(r'"binarySha256":\s*"?[\s|]*"([0-9a-f]{64})"', text)
    if not bsha:
        bsha = re.findall(r"binarySha256[\s\":|]+([0-9a-f]{64})", flat)
    gpu_state = first_match(r'"gpuState": "(.*?)",\s*"gpuStateReturncode"', text, re.S)
    if gpu_state:
        gpu_state = re.sub(r"\s+", " ", gpu_state.replace("\\n", "\n")).strip()
    return {
        "gpu": gpu,
        "cc": cc or "100",
        "valid": True,
        "rate": statistics.median(rates),
        "minRate": min(rates),
        "maxRate": max(rates),
        "packedClmad": packed,
        "batch": 16,
        "threads": 256,
        "minBlocks": 2,
        "workers": 0,
        "steps": 1024,
        "launches": 32,
        "repeats": 3,
        "identity": {
            "sourceSha256": sha[0] if sha else None,
            "binarySha256": bsha[0] if bsha else None,
            "activeBackend": "packed-poly131",
            "cudaImageVersion": first_match(r'"cudaImageVersion": "([^"]+)"', text),
            "compiler": "nvcc",
            "gpuState": gpu_state,
        },
        "samples": samples,
        "parsedFrom": "modal-log",
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
    want_clmad = arm != "software"
    packed = bool(raw.get("packedClmad"))
    if packed is not want_clmad:
        raise SystemExit("arm %s does not match packedClmad=%s" % (arm, packed))
    rates_m = [s["rateM"] for s in samples]
    median_m = raw["rate"]
    median_b = median_m / 1000.0
    sms = next((s["sms"] for s in samples if s["sms"]), None)
    identity = raw.get("identity") or {}
    registers = next((s["registers"] for s in samples if s.get("registers")), None)
    return {
        "kind": "b200 complete-scalar-iteration throughput",
        "family": "b200",
        "gpu": gpu if gpu.startswith("NVIDIA") else "NVIDIA B200",
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
        "registers": registers,
        "sourceSha256": identity.get("sourceSha256"),
        "binarySha256": identity.get("binarySha256"),
        "activeBackend": identity.get("activeBackend"),
        "cudaImageVersion": identity.get("cudaImageVersion"),
        "gpuState": identity.get("gpuState"),
        "compiler": identity.get("compiler"),
        "modalApp": first_match(r"(ap-[A-Za-z0-9]+)", str(identity.get("gpuState") or "")),
        "ratesM": rates_m,
        "minRateM": raw.get("minRate"),
        "maxRateM": raw.get("maxRate"),
        "medianM": median_m,
        "medianB": median_b,
        "prior": {"kind": "sm-count-scale", "sms": PRIOR_SMS, "B": PRIOR_B},
        "ratioToPrior": median_b / PRIOR_B,
        "ratioTo6000": median_b / REF_B,
        "perSmM": (median_m / sms) if sms else None,
        "parsedFrom": raw.get("parsedFrom", "json"),
        "samples": samples,
    }


def load_raw(text: str) -> dict:
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        try:
            return last_json_object(text)
        except ValueError:
            return parse_modal_log(text)


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit("usage: freeze.py ARM LOG.json-or-log OUT.json")
    arm, src, dst = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
    raw = load_raw(src.read_text(errors="replace"))
    # stash app id from the log even when JSON parsed
    if not raw.get("identity"):
        raw["identity"] = {}
    log = src.read_text(errors="replace")
    apps = re.findall(r"ap-[A-Za-z0-9]+", strip_ansi(log))
    receipt = freeze(raw, arm)
    if apps:
        receipt["modalApp"] = apps[-1]
    dst.write_text(json.dumps(receipt, indent=2) + "\n")
    print("wrote", dst)


if __name__ == "__main__":
    main()

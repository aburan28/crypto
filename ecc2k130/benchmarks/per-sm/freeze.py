#!/usr/bin/env python3
"""Turn a Modal ::waves log into the slim receipt WAVES.md cites."""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ANSI = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")
PIPE = 85.0
REF_PER_SM = 15.115792 * 1000.0 / 188.0  # 80.403…


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
        if isinstance(obj, dict) and ("rows" in obj or "waves" in obj):
            last = obj
        i = max(end, j + 1)
    if not isinstance(last, dict):
        raise ValueError("no waves JSON object in log")
    return last


def first_match(pattern: str, text: str, flags=0) -> str | None:
    m = re.search(pattern, text, flags)
    return m.group(1) if m else None


def parse_modal_log(text: str) -> dict:
    text = strip_ansi(text)
    gpu = first_match(r'"gpu": "([^"]+)"', text)
    cc = first_match(r'"cc": "(\d+)"', text)
    sms = int(first_match(r'"sms": (\d+)', text)
              or first_match(r"(\d+) SMs", text) or 0)
    auto = int(first_match(r'"automaticThreads": (\d+)', text)
               or first_match(r"automatic occupancy: (\d+) threads", text) or 0)
    rows = []
    for m in re.finditer(
            r"=== wave (\d+): (\d+) workers ===\s*"
            r"(?:  repeat \d+/\d+: ([0-9.]+) M it/s \(complete\)\s*){3}",
            text):
        wave = int(m.group(1))
        workers = int(m.group(2))
        rates = [float(x) for x in re.findall(
            r"repeat \d+/\d+: ([0-9.]+) M it/s \(complete\)", m.group(0))]
        if len(rates) != 3:
            continue
        median = sorted(rates)[1]
        rows.append({
            "wave": wave,
            "workers": workers,
            "valid": True,
            "rate": median,
            "minRate": min(rates),
            "maxRate": max(rates),
            "ratesM": rates,
            "perSmM": (median / sms) if sms else None,
            "medianB": median / 1000.0,
        })
    if len(rows) < 2:
        raise SystemExit("need at least two complete waves, got %r" % rows)
    return {
        "gpu": gpu,
        "cc": cc,
        "sms": sms,
        "automaticThreads": auto,
        "valid": True,
        "rows": rows,
        "parsedFrom": "modal-log",
    }


def expected_name(modal: str) -> list[str]:
    return {
        "RTX-PRO-6000": ["RTX PRO 6000", "RTX 6000"],
        "L40S": ["L40S"],
        "B200": ["B200"],
        "H200": ["H200"],
        "H100!": ["H100"],
        "L4": ["L4"],
    }.get(modal, [modal.replace("-", " ")])


def freeze(raw: dict, modal: str) -> dict:
    gpu = raw.get("gpu") or ""
    tokens = expected_name(modal)
    if not any(t in gpu for t in tokens):
        raise SystemExit("refusing to freeze %s log with gpu %r" % (modal, gpu))
    if modal == "H100!" and "H200" in gpu:
        raise SystemExit("H100! log upgraded to H200")
    if not raw.get("valid"):
        raise SystemExit("sweep invalid: %s" % raw.get("error"))
    rows = []
    for src in raw.get("rows") or []:
        if not src.get("valid"):
            raise SystemExit("invalid wave %r" % src)
        rate = float(src["rate"])
        sms = raw.get("sms")
        per = src.get("perSmM") or ((rate / sms) if sms else None)
        rows.append({
            "wave": int(src["wave"]),
            "workers": int(src["workers"]),
            "medianM": rate,
            "medianB": rate / 1000.0,
            "minRateM": src.get("minRate"),
            "maxRateM": src.get("maxRate"),
            "ratesM": src.get("ratesM") or [
                s.get("rate") for s in (src.get("samples") or []) if s.get("rate")
            ],
            "perSmM": per,
            "ratioToPipe": (per / PIPE) if per else None,
        })
    if not rows:
        raise SystemExit("no valid wave rows")
    four = next((r for r in rows if r["wave"] == 4), None)
    for row in rows:
        row["ratioToFour"] = (
            (row["perSmM"] / four["perSmM"]) if four and row["perSmM"] else None)
        row["ratioToRefSm"] = (
            (row["perSmM"] / REF_PER_SM) if row["perSmM"] else None)
        row["class"] = "engineering"
    best = max(rows, key=lambda r: r["perSmM"] or 0)
    return {
        "kind": "per-SM occupancy wave sweep",
        "modal": modal,
        "gpu": gpu,
        "cc": str(raw.get("cc") or ""),
        "valid": True,
        "sms": raw.get("sms"),
        "automaticThreads": raw.get("automaticThreads"),
        "residentBlocks": raw.get("residentBlocks"),
        "blockThreads": raw.get("blockThreads") or 256,
        "pipeM": PIPE,
        "referencePerSmM": REF_PER_SM,
        "rows": rows,
        "bestWave": best["wave"],
        "bestPerSmM": best["perSmM"],
        "fourBeaten": bool(
            four and best["wave"] != 4 and best["perSmM"]
            and best["perSmM"] > four["perSmM"] * 1.01),
        "parsedFrom": raw.get("parsedFrom", "json"),
    }


def load_raw(text: str) -> dict:
    try:
        obj = json.loads(text)
        if isinstance(obj, dict) and obj.get("rows"):
            return obj
    except json.JSONDecodeError:
        pass
    try:
        return last_json_object(strip_ansi(text))
    except ValueError:
        return parse_modal_log(text)


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit("usage: freeze.py MODAL-TYPE LOG OUT.json")
    modal, src, dst = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
    raw = load_raw(src.read_text(errors="replace"))
    receipt = freeze(raw, modal)
    apps = re.findall(r"ap-[A-Za-z0-9]+", strip_ansi(src.read_text(errors="replace")))
    if apps:
        receipt["modalApp"] = apps[-1]
    dst.write_text(json.dumps(receipt, indent=2) + "\n")
    print("wrote", dst)


if __name__ == "__main__":
    main()

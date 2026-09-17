#!/usr/bin/env python3
"""Recover a bench receipt from a Modal log that may contain spinner junk."""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ANSI = re.compile(r"\x1b\[[0-9;?]*[ -/]*[@-~]|\x1b][^\x07]*\x07|\r")
FINISHED = re.compile(r"finished:\s+(\d+(?:\.\d+)?)\s+M it/s")
FLAG = re.compile(r"packed top clmad:\s+([01])")


def clean(text: str) -> str:
    return ANSI.sub("", text)


def from_json(text: str) -> dict | None:
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


def from_lines(text: str) -> dict | None:
    rates = [float(m.group(1)) for m in FINISHED.finditer(text)]
    flags = [int(m.group(1)) for m in FLAG.finditer(text)]
    if len(rates) < 3:
        return None
    rates = rates[-3:]
    top = flags[-1] if flags else None
    samples = [{"valid": True, "rate": r} for r in rates]
    return {
        "valid": True,
        "rate": sorted(rates)[1],
        "samples": samples,
        "packedTopClmad": bool(top) if top is not None else None,
        "expectedPackedTopClmad": bool(top) if top is not None else None,
        "recoveredFromLog": True,
        "sampleRatesM": rates,
    }


def extract(log_path: Path) -> dict:
    text = clean(log_path.read_text(encoding="utf-8", errors="replace"))
    obj = from_json(text) or from_lines(text)
    if obj is None:
        raise SystemExit("no recoverable receipt in " + str(log_path))
    return obj


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: extract.py LOG JSON")
    obj = extract(Path(sys.argv[1]))
    Path(sys.argv[2]).write_text(json.dumps(obj, indent=2) + "\n")
    print(
        "wrote",
        sys.argv[2],
        "valid=",
        obj.get("valid"),
        "rate=",
        obj.get("rate"),
        "topClmad=",
        obj.get("packedTopClmad"),
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Split a callgrind profile of `ic price` by measurement phase.

`ic_measurement` issues a valgrind `DUMP_STATS_AT` request at every phase
boundary, labelled with the phase that just ended. A run therefore leaves
one part file per segment: `<out>.1`, `<out>.2`, and so on. Each part's
trigger names its phase, and its `summary:` line holds its counts.

This sums the parts by label.
- The online phases' counts are then exact.
- The first `generic_ic_setup` part runs from process start to the first
  online boundary: start-up, the unit's calibration, the curve, and the
  index calculus's whole set-up.
- Later `setup` parts are the work between intervals (replays, rho's set-up).
- The `Program termination` part is what follows the last boundary.

The first set-up part is also annotated by function, exclusive costs.

    python3 callgrind_phases.py <dir> <stem> [--top 30] > phases.json

`<stem>` is the `--callgrind-out-file` path's basename, without the part
suffix.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path

EVENTS_RE = re.compile(r"^events:\s+(.*)$")
SUMMARY_RE = re.compile(r"^(?:summary|totals):\s+(.*)$")
TRIGGER_RE = re.compile(r"^desc: Trigger:\s+(.*)$")


def read_part(path: Path) -> tuple[str, dict]:
    events, totals, trigger = [], None, "unknown"
    with open(path, errors="replace") as f:
        for line in f:
            if m := TRIGGER_RE.match(line):
                trigger = m.group(1).strip()
                trigger = trigger.removeprefix("Client Request: ")
            elif m := EVENTS_RE.match(line):
                events = m.group(1).split()
            elif m := SUMMARY_RE.match(line):
                totals = [int(x) for x in m.group(1).split()]
    counts = {e: (totals[i] if totals and i < len(totals) else 0) for i, e in enumerate(events)}
    return trigger, counts


def part_number(path: Path) -> int:
    suffix = path.name.rsplit(".", 1)[-1]
    return int(suffix) if suffix.isdigit() else 0


def annotate(part: Path, top: int) -> list[dict]:
    text = subprocess.run(["callgrind_annotate", "--inclusive=no", "--threshold=99", str(part)],
                          capture_output=True, text=True, check=False).stdout
    rows, events, in_table = [], [], False
    for line in text.splitlines():
        if line.startswith("Events shown:"):
            events = line.split(":", 1)[1].split()
        if "file:function" in line:
            in_table = True
            continue
        if not in_table or not line.strip() or line.startswith("-"):
            continue
        # "1,028,390,508 (26.93%) 16,136 ( 8.64%) ... ???:fn [obj]": drop the
        # percentages, whose padding splits them, then read the counts.
        cells = re.sub(r"\(\s*[\d.]+%\)", " ", line).split()
        numbers, i = [], 0
        while i < len(cells) and len(numbers) < len(events):
            c = cells[i]
            if c == ".":
                numbers.append(0)
            elif re.fullmatch(r"[\d,]+", c):
                numbers.append(int(c.replace(",", "")))
            else:
                break
            i += 1
        if len(numbers) < len(events):
            continue
        name = " ".join(cells[i:])
        name = re.sub(r"\s*\[.*\]$", "", name)
        rows.append({"function": name.split(":", 1)[-1] if name.startswith("???:") else name,
                     **dict(zip(events, numbers))})
        if len(rows) >= top:
            break
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("dir", type=Path)
    ap.add_argument("stem")
    ap.add_argument("--top", type=int, default=30)
    args = ap.parse_args()
    parts = sorted((p for p in args.dir.glob(args.stem + ".*") if p.name.rsplit(".", 1)[-1].isdigit()),
                   key=part_number)
    by_label: dict[str, dict] = {}
    first_setup = None
    for p in parts:
        label, counts = read_part(p)
        slot = by_label.setdefault(label, {"parts": 0, "counts": {}})
        slot["parts"] += 1
        for e, v in counts.items():
            slot["counts"][e] = slot["counts"].get(e, 0) + v
        if first_setup is None and label == "generic_ic_setup":
            first_setup = p
    total = {}
    for slot in by_label.values():
        for e, v in slot["counts"].items():
            total[e] = total.get(e, 0) + v
    doc = {
        "stem": args.stem, "parts": len(parts), "total": total,
        "by_phase": dict(sorted(by_label.items(), key=lambda kv: -kv[1]["counts"].get("Ir", 0))),
        "first_setup_part": first_setup.name if first_setup else None,
        "first_setup_part_counts": read_part(first_setup)[1] if first_setup else None,
        "first_setup_part_top_functions": annotate(first_setup, args.top) if first_setup else [],
    }
    print(json.dumps(doc, indent=1))


if __name__ == "__main__":
    main()

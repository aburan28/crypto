#!/usr/bin/env python3
"""Validate and summarize the bounded split-forward GPU comparison."""

import argparse
import csv
import hashlib
import json
import math
import pathlib
import re
import struct


HEADER = struct.Struct("<8sII")
MAGIC = b"ECC2KDT3"
VERSION = 3
RECORD_BYTES = 32


def read_corpus(path):
    data = pathlib.Path(path).read_bytes()
    if len(data) < HEADER.size:
        raise ValueError(f"{path}: shorter than the {HEADER.size}-byte v3 header")
    magic, version, stride = HEADER.unpack_from(data)
    if magic != MAGIC or version != VERSION or stride != RECORD_BYTES:
        raise ValueError(
            f"{path}: expected {MAGIC!r}/v{VERSION}/{RECORD_BYTES}, "
            f"got {magic!r}/v{version}/{stride}"
        )
    payload = data[HEADER.size :]
    if len(payload) % stride:
        raise ValueError(f"{path}: {len(payload)} payload bytes is not a whole record count")
    records = [payload[i : i + stride] for i in range(0, len(payload), stride)]
    digest = hashlib.sha256(b"".join(sorted(records))).hexdigest()
    return {
        "magic": magic.decode("ascii"),
        "version": version,
        "recordBytes": stride,
        "records": len(records),
        "payloadBytes": len(payload),
        "sortedSha256": digest,
    }


def parse_build(path):
    text = pathlib.Path(path).read_text(errors="replace")
    marker = "Function properties for _ZN12eccPacked1314walkE10WalkParamsIjEPj"
    lines = text.splitlines()
    hits = [i for i, line in enumerate(lines) if marker in line]
    if len(hits) != 1:
        raise ValueError(f"{path}: packed walk resource marker count {len(hits)}")
    block = "\n".join(lines[hits[0] : hits[0] + 7])
    frame = re.search(
        r"(\d+) bytes stack frame, (\d+) bytes spill stores, (\d+) bytes spill loads",
        block,
    )
    registers = re.search(r"Used (\d+) registers", block)
    if not frame or not registers:
        raise ValueError(f"{path}: incomplete packed walk resource block")
    return {
        "registers": int(registers.group(1)),
        "stackFrameBytes": int(frame.group(1)),
        "spillStoreBytes": int(frame.group(2)),
        "spillLoadBytes": int(frame.group(3)),
        "raw": block,
    }


def parse_verify(path, expected_mode):
    text = pathlib.Path(path).read_text(errors="replace")
    finished = re.findall(
        r"finished: ([0-9.]+) M it/s, (\d+) distinguished points "
        r"\((\d+) verified against the reference, (\d+) dropped\)",
        text,
    )
    markers = re.findall(r"^packed table split forward: ([01])$", text, re.MULTILINE)
    errors = []
    if len(finished) != 1:
        errors.append(f"finished marker count {len(finished)}")
    if markers != [str(expected_mode)]:
        errors.append(f"identity markers {markers!r}, expected {[str(expected_mode)]!r}")
    if re.search(r"MISMATCH|OVERFLOW", text):
        errors.append("mismatch or overflow marker")
    row = finished[0] if len(finished) == 1 else (None, None, None, None)
    if row[2] != "300" or row[3] != "0":
        errors.append(f"expected 300 verified and 0 dropped, got {row[2]}/{row[3]}")
    return {
        "valid": not errors,
        "errors": errors,
        "rateMps": float(row[0]) if row[0] is not None else None,
        "distinguishedPoints": int(row[1]) if row[1] is not None else None,
        "verified": int(row[2]) if row[2] is not None else None,
        "dropped": int(row[3]) if row[3] is not None else None,
        "mode": int(markers[0]) if len(markers) == 1 else None,
    }


def read_samples(path):
    rows = []
    with pathlib.Path(path).open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            row["pair"] = int(row["pair"])
            row["order"] = int(row["order"])
            row["rateMps"] = float(row["rateMps"])
            if not math.isfinite(row["rateMps"]) or row["rateMps"] <= 0:
                raise ValueError(f"non-positive/non-finite rate in {row}")
            rows.append(row)
    return rows


def decide_samples(rows):
    errors = []
    warmup = [r for r in rows if r["phase"] == "warmup"]
    screen = [r for r in rows if r["phase"] == "screen"]
    confirm = [r for r in rows if r["phase"] == "confirm"]
    if [r["variant"] for r in warmup] != ["control", "candidate"]:
        errors.append("warmups must be exactly control,candidate")
    if [r["variant"] for r in screen] != ["control", "candidate", "control"]:
        errors.append("screen must be exactly control,candidate,control")
    qualified = False
    threshold = None
    screen_ratio = None
    if len(screen) == 3:
        controls = [screen[0]["rateMps"], screen[2]["rateMps"]]
        threshold = max(controls) * 1.005
        screen_ratio = screen[1]["rateMps"] / max(controls)
        qualified = screen[1]["rateMps"] >= threshold
    expected_confirmation = [
        (1, 1, "control"), (1, 2, "candidate"),
        (2, 1, "candidate"), (2, 2, "control"),
        (3, 1, "control"), (3, 2, "candidate"),
    ]
    got_confirmation = [(r["pair"], r["order"], r["variant"]) for r in confirm]
    if confirm and got_confirmation != expected_confirmation:
        errors.append(f"confirmation order is incomplete/different: {got_confirmation!r}")
    if not qualified and confirm:
        errors.append("unqualified candidate received confirmation")
    pairs = []
    if got_confirmation == expected_confirmation:
        for pair in (1, 2, 3):
            group = {r["variant"]: r["rateMps"] for r in confirm if r["pair"] == pair}
            pairs.append(group["candidate"] / group["control"])
    if errors:
        decision = "invalid"
    elif not qualified:
        decision = "unqualified_stop"
    elif not confirm:
        decision = "qualified_pending_confirmation"
    else:
        decision = "confirmation_complete"
    return {
        "valid": not errors,
        "errors": errors,
        "qualified": qualified,
        "qualificationThresholdMps": threshold,
        "screenRatioToFasterControl": screen_ratio,
        "decision": decision,
        "pairedConfirmationRatios": pairs,
        "samples": rows,
    }


def summarize(results, include_timing=True):
    results = pathlib.Path(results)
    errors = []
    output = {
        "schema": "ecc2k130_split_forward_screen.v1",
        "claimClass": "same-walk kernel engineering",
        "target": "candidate >= 1.005 * faster bracket control",
    }
    try:
        output["build"] = {
            name: parse_build(results / f"build-{name}.log")
            for name in ("control", "candidate")
        }
    except (OSError, ValueError) as exc:
        errors.append(str(exc))
    try:
        output["verify"] = {
            "control": parse_verify(results / "verify-control.log", 0),
            "candidate": parse_verify(results / "verify-candidate.log", 1),
        }
        for name, row in output["verify"].items():
            if not row["valid"]:
                errors.append(f"{name} replay: {row['errors']}")
    except OSError as exc:
        errors.append(str(exc))
    try:
        corpora = {
            name: read_corpus(results / f"dp-{name}.bin")
            for name in ("control", "candidate")
        }
        output["corpora"] = corpora
        output["corpusIdentity"] = (
            corpora["control"]["records"] == corpora["candidate"]["records"]
            and corpora["control"]["sortedSha256"] == corpora["candidate"]["sortedSha256"]
        )
        if not output["corpusIdentity"]:
            errors.append("header-aware corpus identities differ")
    except (OSError, ValueError) as exc:
        errors.append(str(exc))
    if include_timing:
        try:
            output["timing"] = decide_samples(read_samples(results / "samples.tsv"))
            if not output["timing"]["valid"]:
                errors.extend(output["timing"]["errors"])
        except (OSError, ValueError) as exc:
            errors.append(str(exc))
    for optional in ("host.txt", "source-files.sha256", "binary-sha256.txt"):
        path = results / optional
        if path.exists():
            output[optional] = path.read_text(errors="replace").strip()
    output["valid"] = not errors
    output["errors"] = errors
    if not output["valid"]:
        output["decision"] = "invalid"
    elif include_timing:
        output["decision"] = output["timing"]["decision"]
    else:
        output["decision"] = "preflight_pass"
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results")
    parser.add_argument("--out", required=True)
    parser.add_argument("--preflight", action="store_true", help="validate builds, replay and v3 corpora before timing")
    args = parser.parse_args()
    summary = summarize(args.results, include_timing=not args.preflight)
    pathlib.Path(args.out).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "valid": summary["valid"],
        "decision": summary["decision"],
        "errors": summary["errors"],
        "qualified": summary.get("timing", {}).get("qualified"),
        "screenRatio": summary.get("timing", {}).get("screenRatioToFasterControl"),
        "pairedRatios": summary.get("timing", {}).get("pairedConfirmationRatios"),
    }, sort_keys=True))
    return 0 if summary["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

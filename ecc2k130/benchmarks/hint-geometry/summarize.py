#!/usr/bin/env python3
"""Validate and summarize the B16/B32/B64 reconverged-hints sweep."""

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import pathlib
import re
import statistics


HERE = pathlib.Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("batch_hints_summary", HERE.parent / "batch-hints" / "summarize.py")
BASE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(BASE)

ARMS = ("b16", "b32", "b64")
GEOMETRY = {
    "b16": {"batch": 16, "blockThreads": 512, "runtimeThreads": 96256},
    "b32": {"batch": 32, "blockThreads": 256, "runtimeThreads": 48128},
    "b64": {"batch": 64, "blockThreads": 128, "runtimeThreads": 24064},
}
LIVE_SLOTS = 1_540_096
STEPS = 1024
LAUNCHES = {"warmup": 16, "screen": 32, "confirm": 64}


def parse_verify(path, arm):
    row = BASE.parse_verify(path, 1)
    text = pathlib.Path(path).read_text(errors="replace")
    expected = GEOMETRY[arm]
    backend = re.findall(
        r"^backend cuda-packed131: (\d+) threads x (\d+) slots x 1 lanes = (\d+) walks,",
        text,
        re.MULTILINE,
    )
    if len(backend) != 1:
        row["errors"].append(f"backend geometry marker count {len(backend)}")
    else:
        runtime, batch, walks = map(int, backend[0])
        if runtime != expected["runtimeThreads"] or batch != expected["batch"]:
            row["errors"].append(f"backend geometry {runtime}x{batch}, expected {expected['runtimeThreads']}x{expected['batch']}")
        if walks != expected["runtimeThreads"] * expected["batch"]:
            row["errors"].append(f"backend walk count {walks} is inconsistent")
        row.update(runtimeThreads=runtime, batch=batch, liveSlots=walks)
    row["valid"] = not row["errors"]
    return row


def read_samples(path):
    rows = []
    with pathlib.Path(path).open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            row["pair"] = int(row["pair"])
            row["order"] = int(row["order"])
            row["rateMps"] = float(row["rateMps"])
            if row["phase"] not in LAUNCHES or row["variant"] not in ARMS or not math.isfinite(row["rateMps"]) or row["rateMps"] <= 0:
                raise ValueError(f"invalid timing row {row}")
            if not re.fullmatch(r"[0-9a-f]{64}", row["logSha256"]):
                raise ValueError(f"invalid timing-log digest in {row}")
            rows.append(row)
    return rows


def validate_sample_logs(results, rows):
    errors = []
    for row in rows:
        arm = row["variant"]
        name = f"{row['phase']}-{row['pair']}-{row['order']}-{arm}.log"
        path = pathlib.Path(results) / name
        try:
            data = path.read_bytes()
        except OSError as exc:
            errors.append(str(exc))
            continue
        digest = hashlib.sha256(data).hexdigest()
        if digest != row["logSha256"]:
            errors.append(f"{name}: sha256 {digest}, TSV has {row['logSha256']}")
        text = data.decode(errors="replace")
        expected = GEOMETRY[arm]
        backend = re.findall(
            r"^backend cuda-packed131: (\d+) threads x (\d+) slots x 1 lanes = (\d+) walks,",
            text,
            re.MULTILINE,
        )
        wanted = [(str(expected["runtimeThreads"]), str(expected["batch"]), str(LIVE_SLOTS))]
        if backend != wanted:
            errors.append(f"{name}: backend geometry markers {backend!r}, expected {wanted!r}")
        finished = re.findall(
            r"finished: ([0-9.]+) M it/s, (\d+) distinguished points "
            r"\((\d+) verified against the reference, (\d+) dropped\)",
            text,
        )
        if len(finished) != 1:
            errors.append(f"{name}: finished marker count {len(finished)}")
        else:
            rate, _, _, dropped = finished[0]
            if float(rate) != row["rateMps"]:
                errors.append(f"{name}: final rate {rate}, TSV has {row['rateMps']}")
            if dropped != "0":
                errors.append(f"{name}: timed row dropped {dropped} reports")
        updates = [int(value) for value in re.findall(r"\s(\d+) iterations\s", text)]
        expected_work = LIVE_SLOTS * STEPS * LAUNCHES[row["phase"]]
        if not updates or updates[-1] != expected_work:
            errors.append(f"{name}: final work {updates[-1] if updates else None}, expected {expected_work}")
    return errors


def decide_samples(rows):
    errors = []
    warmup = [r for r in rows if r["phase"] == "warmup"]
    screen = [r for r in rows if r["phase"] == "screen"]
    confirm = [r for r in rows if r["phase"] == "confirm"]
    if [r["variant"] for r in warmup] != list(ARMS):
        errors.append("warmups must be exactly b16,b32,b64")
    if [(r["pair"], r["order"]) for r in warmup] != [(0, 1), (0, 2), (0, 3)]:
        errors.append("warmup pair/order keys differ")
    expected_screen = ["b16", "b32", "b64", "b64", "b32", "b16"]
    if [r["variant"] for r in screen] != expected_screen:
        errors.append(f"screen order differs: {[r['variant'] for r in screen]}")
    expected_screen_keys = [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3)]
    if [(r["pair"], r["order"]) for r in screen] != expected_screen_keys:
        errors.append("screen pass/order keys differ")
    medians = {}
    if len(screen) == 6:
        for arm in ARMS:
            medians[arm] = statistics.median(r["rateMps"] for r in screen if r["variant"] == arm)
    top = max(("b32", "b64"), key=lambda arm: medians.get(arm, 0.0)) if medians else None
    ratio = medians[top] / medians["b16"] if top and medians.get("b16") else None
    qualified = bool(ratio is not None and ratio >= 1.005)
    expected_confirmation = []
    if top:
        expected_confirmation = [
            (1, 1, "b16"), (1, 2, top),
            (2, 1, top), (2, 2, "b16"),
            (3, 1, "b16"), (3, 2, top),
        ]
    got = [(r["pair"], r["order"], r["variant"]) for r in confirm]
    if confirm and got != expected_confirmation:
        errors.append(f"confirmation order is incomplete/different: {got!r}")
    if confirm and not qualified:
        errors.append("unqualified arm received confirmation")
    paired = []
    if got and got == expected_confirmation:
        for pair in (1, 2, 3):
            group = {r["variant"]: r["rateMps"] for r in confirm if r["pair"] == pair}
            paired.append(group[top] / group["b16"])
    if errors:
        decision = "invalid"
    elif not qualified:
        decision = "reference_retained"
    elif not confirm:
        decision = "qualified_pending_confirmation"
    else:
        decision = "confirmation_complete"
    return {
        "valid": not errors,
        "errors": errors,
        "screenMediansMps": medians,
        "topCandidate": top,
        "screenRatioToB16": ratio,
        "qualified": qualified,
        "decision": decision,
        "pairedConfirmationRatios": paired,
        "samples": rows,
    }


def summarize(results, include_timing=True):
    root = pathlib.Path(results)
    errors = []
    out = {
        "schema": "ecc2k130_hint_geometry_sweep.v1",
        "claimClass": "same-walk kernel engineering",
        "geometry": GEOMETRY,
        "screenLaunches": 32,
        "screenSlotPopulations": 49283072,
        "qualificationRatio": 1.005,
    }
    try:
        out["build"] = {arm: BASE.parse_build(root / f"build-{arm}.log") for arm in ARMS}
    except (OSError, ValueError) as exc:
        errors.append(str(exc))
    try:
        out["verify"] = {arm: parse_verify(root / f"verify-{arm}.log", arm) for arm in ARMS}
        for arm, row in out["verify"].items():
            if not row["valid"]:
                errors.append(f"{arm} replay: {row['errors']}")
    except OSError as exc:
        errors.append(str(exc))
    try:
        corpora = {arm: BASE.read_corpus(root / f"dp-{arm}.bin") for arm in ARMS}
        out["corpora"] = corpora
        identities = {(row["records"], row["sortedSha256"]) for row in corpora.values()}
        out["corpusIdentity"] = len(identities) == 1
        if not out["corpusIdentity"]:
            errors.append("header-aware corpus identities differ")
    except (OSError, ValueError) as exc:
        errors.append(str(exc))
    if include_timing:
        try:
            samples = read_samples(root / "samples.tsv")
            out["timing"] = decide_samples(samples)
            out["timing"]["logErrors"] = validate_sample_logs(root, samples)
            if out["timing"]["logErrors"]:
                errors.extend(out["timing"]["logErrors"])
            if not out["timing"]["valid"]:
                errors.extend(out["timing"]["errors"])
        except (OSError, ValueError) as exc:
            errors.append(str(exc))
    for name in ("host.txt", "geometry.tsv", "source-files.sha256", "binary-sha256.txt"):
        path = root / name
        if path.exists():
            out[name] = path.read_text(errors="replace").strip()
        else:
            errors.append(f"missing evidence manifest: {path}")
    out["valid"] = not errors
    out["errors"] = errors
    out["decision"] = "invalid" if errors else (out["timing"]["decision"] if include_timing else "preflight_pass")
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results")
    parser.add_argument("--out", required=True)
    parser.add_argument("--preflight", action="store_true")
    args = parser.parse_args()
    result = summarize(args.results, include_timing=not args.preflight)
    pathlib.Path(args.out).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    timing = result.get("timing", {})
    print(json.dumps({"valid": result["valid"], "decision": result["decision"], "errors": result["errors"],
                      "topCandidate": timing.get("topCandidate"), "qualified": timing.get("qualified"),
                      "screenRatio": timing.get("screenRatioToB16"), "pairedRatios": timing.get("pairedConfirmationRatios")},
                     sort_keys=True))
    return 0 if result["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

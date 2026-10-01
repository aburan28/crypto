#!/usr/bin/env python3
"""Validate per-lane versus block-compacted exact v3 cycle hints."""

import argparse
import hashlib
import importlib.util
import json
import pathlib
import re


HERE = pathlib.Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "batch_hints_summary", HERE.parent / "batch-hints" / "summarize.py"
)
BASE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(BASE)


def parse_verify(path, expected_block):
    row = BASE.parse_verify(path, 1)
    text = pathlib.Path(path).read_text(errors="replace")
    fast2 = re.findall(r"^packed cycle fast2: ([01])$", text, re.MULTILINE)
    markers = re.findall(r"^packed table block hints: ([01]), queue (\d+)$", text, re.MULTILINE)
    if fast2 != ["1"]:
        row["errors"].append(f"fast2 markers {fast2!r}, expected ['1']")
    if markers != [(str(expected_block), "512")]:
        row["errors"].append(
            f"block-hint markers {markers!r}, expected {[(str(expected_block), '512')]!r}"
        )
    backend = re.findall(
        r"^backend cuda-packed131: (\d+) threads x (\d+) slots x 1 lanes = (\d+) walks,",
        text,
        re.MULTILINE,
    )
    if len(backend) != 1:
        row["errors"].append(f"backend geometry marker count {len(backend)}")
    else:
        threads, batch, walks = map(int, backend[0])
        if batch != 16 or walks != threads * batch:
            row["errors"].append(
                f"backend geometry {threads}x{batch}={walks}, expected B16 consistency"
            )
        row.update(runtimeThreads=threads, batch=batch, liveSlots=walks)
    row["fast2"] = int(fast2[0]) if len(fast2) == 1 else None
    row["blockHints"] = int(markers[0][0]) if len(markers) == 1 else None
    row["hintQueue"] = int(markers[0][1]) if len(markers) == 1 else None
    row["valid"] = not row["errors"]
    return row


def validate_block_timing_markers(results, rows):
    errors = []
    for row in rows:
        name = f"{row['phase']}-{row['pair']}-{row['order']}-{row['variant']}.log"
        try:
            text = (pathlib.Path(results) / name).read_text(errors="replace")
        except OSError as exc:
            errors.append(str(exc))
            continue
        fast2 = re.findall(r"^packed cycle fast2: ([01])$", text, re.MULTILINE)
        markers = re.findall(r"^packed table block hints: ([01]), queue (\d+)$", text, re.MULTILINE)
        expected = "1" if row["variant"] == "candidate" else "0"
        if fast2 != ["1"]:
            errors.append(f"{name}: fast2 markers {fast2!r}, expected ['1']")
        if markers != [(expected, "512")]:
            errors.append(f"{name}: block-hint markers {markers!r}, expected {[(expected, '512')]!r}")
    return errors


def summarize(results, include_timing=True):
    root = pathlib.Path(results)
    errors = []
    output = {
        "schema": "ecc2k130_block_hints_screen.v1",
        "claimClass": "same-walk kernel engineering",
        "common": {
            "batch": 16,
            "blockThreads": 512,
            "minBlocks": 1,
            "splitForward": 1,
            "batchHints": 1,
            "fast2": 1,
            "hintQueue": 512,
        },
        "target": "block hints >= 1.005 * faster bracketing control",
    }
    try:
        output["build"] = {
            name: BASE.parse_build(root / f"build-{name}.log")
            for name in ("control", "candidate")
        }
    except (OSError, ValueError) as exc:
        errors.append(str(exc))
    try:
        output["verify"] = {
            "control": parse_verify(root / "verify-control.log", 0),
            "candidate": parse_verify(root / "verify-candidate.log", 1),
        }
        for name, row in output["verify"].items():
            if not row["valid"]:
                errors.append(f"{name} reference replay: {row['errors']}")
    except OSError as exc:
        errors.append(str(exc))
    try:
        corpora = {
            name: BASE.read_corpus(root / f"dp-{name}.bin")
            for name in ("control", "candidate")
        }
        output["corpora"] = corpora
        output["corpusIdentity"] = (
            corpora["control"]["records"] == corpora["candidate"]["records"]
            and corpora["control"]["sortedSha256"]
            == corpora["candidate"]["sortedSha256"]
        )
        if not output["corpusIdentity"]:
            errors.append("header-aware table-v3 corpus identities differ")
    except (OSError, ValueError) as exc:
        errors.append(str(exc))
    if include_timing:
        try:
            samples = BASE.read_samples(root / "samples.tsv")
            output["timing"] = BASE.decide_samples(samples)
            output["timing"]["logErrors"] = (
                BASE.validate_sample_logs(root, samples)
                + validate_block_timing_markers(root, samples)
            )
            if output["timing"]["logErrors"]:
                errors.extend(output["timing"]["logErrors"])
            if not output["timing"]["valid"]:
                errors.extend(output["timing"]["errors"])
        except (OSError, ValueError) as exc:
            errors.append(str(exc))
    for name in (
        "host.txt",
        "source-files.sha256",
        "binary-sha256.txt",
        "dp-identity.txt",
    ):
        path = root / name
        if path.exists():
            output[name] = path.read_text(errors="replace").strip()
            output[name + "Sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
        else:
            errors.append(f"missing evidence manifest: {path}")
    output["valid"] = not errors
    output["errors"] = errors
    if errors:
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
    parser.add_argument("--preflight", action="store_true")
    args = parser.parse_args()
    result = summarize(args.results, include_timing=not args.preflight)
    pathlib.Path(args.out).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    timing = result.get("timing", {})
    print(
        json.dumps(
            {
                "valid": result["valid"],
                "decision": result["decision"],
                "errors": result["errors"],
                "qualified": timing.get("qualified"),
                "screenRatio": timing.get("screenRatioToFasterControl"),
                "pairedRatios": timing.get("pairedConfirmationRatios"),
            },
            sort_keys=True,
        )
    )
    return 0 if result["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

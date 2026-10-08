#!/usr/bin/env python3
"""Summarise the n=61 L=65536 compact-orbit vs KS v2 batched-rho panel.

Reads the run directory written by growing_n_n61_L65536_20260930.sh and
writes growing_n_summary.json beside it.  Every ratio is compact / rho on the
same block pair.  Walls come from /usr/bin/time -l on a shared macOS host and
are not isolated (tools/isolated_bench.py is Linux-only); retired instructions
are the whole-process hardware counter reported by the same tool.
"""

from __future__ import annotations

import hashlib
import json
import statistics
import sys
from pathlib import Path

L = 65536
N = 61


def parse_time(path: Path) -> dict:
    text = path.read_text().split("\n")
    first = text[0].split()
    record = {"real_s": float(first[0]), "user_s": float(first[2]), "sys_s": float(first[4])}
    keys = {
        "maximum resident set size": "max_rss_bytes",
        "instructions retired": "instructions_retired",
        "cycles elapsed": "cycles_elapsed",
        "involuntary context switches": "involuntary_context_switches",
        "peak memory footprint": "peak_memory_footprint_bytes",
        "page faults": "page_faults",
        "swaps": "swaps",
    }
    for line in text[1:]:
        for label, key in keys.items():
            if line.strip().endswith(label):
                record[key] = int(line.split()[0])
    return record


def read_jsonl(path: Path):
    with path.open() as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def ic_block(out: Path, b: int) -> dict:
    summary = json.loads((out / f"ic_n61_b{b}.summary.json").read_text())
    scalars, points, probes, verified = [], [], 0, True
    digest = hashlib.sha256()
    for record in read_jsonl(out / f"ic_n61_b{b}.jsonl"):
        if record.get("kind") != "compact_orbit_dlp_target":
            continue
        scalars.append(record["published_fixture_scalar"])
        points.append(tuple(record["target"]))
        digest.update(json.dumps({k: v for k, v in record.items() if "_ms" not in k}, sort_keys=True).encode())
        probes += record.get("probes") or 0
        verified &= bool(record.get("recovered_matches_published")) and bool(record.get("group_verified"))
    return {
        "time": parse_time(out / f"ic_n61_b{b}.time"),
        "load": (out / f"ic_n61_b{b}.load").read_text().strip(),
        "summary": summary,
        "target_probes_total": probes,
        "all_targets_verified": verified,
        "untimed_record_sha256": digest.hexdigest(),
        "scalars": scalars,
        "points": points,
    }


def ks_block(out: Path, b: int) -> dict:
    scalars, points, verified, summary = [], [], True, None
    digest = hashlib.sha256()
    for record in read_jsonl(out / f"ks_n61_b{b}.jsonl"):
        if record.get("kind") == "rho_ks_batch_fixture":
            scalars.append(record["published_fixture_scalar"])
            points.append(tuple(record["published_q"]))
            digest.update(json.dumps({k: v for k, v in record.items() if "_ms" not in k}, sort_keys=True).encode())
            verified &= record["recovered_fixture_scalar"] == record["published_fixture_scalar"] and record["verified"]
        elif record.get("kind") == "rho_ks_batch_summary":
            summary = record
    return {
        "time": parse_time(out / f"ks_n61_b{b}.time"),
        "load": (out / f"ks_n61_b{b}.load").read_text().strip(),
        "summary": summary,
        "all_targets_verified": verified,
        "untimed_record_sha256": digest.hexdigest(),
        "scalars": scalars,
        "points": points,
    }


def replay_summary(out: Path) -> dict:
    blocks = {}
    for path in sorted(out.glob("replay_b*/independent_replay.json")):
        field = json.loads(path.read_text())["fields"]["61"]
        blocks[path.parent.name] = {"records": field["records"], "pass": field["pass"], "fail": field["fail"]}
    return {
        "blocks": blocks,
        "records": sum(v["records"] for v in blocks.values()),
        "pass": sum(v["pass"] for v in blocks.values()),
        "fail": sum(v["fail"] for v in blocks.values()),
    }


def main() -> int:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).with_suffix("").parent / "growing_n_n61_L65536_20260930"
    K = int((out / "chosen_K_n61.txt").read_text())
    corpus = [int(x) for x in (out / f"scalars_n61-ks-growing-{L}-v1.txt").read_text().split()]
    tune = {}
    for path in sorted(out.glob("tune_n61_K*.done")):
        k = path.stem.split("_K")[1].split(".")[0]
        state = path.read_text().strip()
        entry = {"state": state}
        if state == "ok":
            entry["time"] = parse_time(out / f"tune_n61_K{k}.time")
            entry["summary"] = json.loads((out / f"tune_n61_K{k}.summary.json").read_text())
        tune[k] = entry

    blocks = []
    b = 0
    while (out / f"ic_n61_b{b}.done").exists() or (out / f"ks_n61_b{b}.done").exists():
        row = {"block": b, "order": "ic_then_rho" if b % 2 == 0 else "rho_then_ic"}
        ic = ic_block(out, b) if (out / f"ic_n61_b{b}.done").exists() else None
        ks = ks_block(out, b) if (out / f"ks_n61_b{b}.done").exists() else None
        if ic:
            row["ic"] = {k: v for k, v in ic.items() if k not in ("scalars", "points")}
            row["ic_matches_corpus"] = ic["scalars"] == corpus
        if ks:
            row["ks"] = {k: v for k, v in ks.items() if k not in ("scalars", "points")}
            row["ks_matches_corpus"] = ks["scalars"] == corpus
        if ic and ks:
            row["same_target_points"] = ic["points"] == ks["points"]
            it, kt = ic["time"], ks["time"]
            row["wall_ratio"] = it["real_s"] / kt["real_s"]
            row["user_ratio"] = it["user_s"] / kt["user_s"]
            if it.get("instructions_retired") and kt.get("instructions_retired"):
                row["instructions_ratio"] = it["instructions_retired"] / kt["instructions_retired"]
            row["paired_complete"] = True
        blocks.append(row)
        b += 1

    paired = [r for r in blocks if r.get("paired_complete")]

    def med(key):
        values = [r[key] for r in paired if key in r]
        return {"median": statistics.median(values), "range": [min(values), max(values)], "values": values} if values else None

    result = {
        "L": L, "n": N, "K": K,
        "comparator": "frozen Kuhn-Struik batched signed-Frobenius rho v2 (koblitz_rho_batch_ks_v2_n61, d=4); historical weaker comparator, superseded as the operative reference by the 2026-09-29 strong-rho ladder",
        "tune": tune,
        "blocks": blocks,
        "paired_blocks": len(paired),
        "wall_ratio": med("wall_ratio"),
        "user_ratio": med("user_ratio"),
        "instructions_ratio": med("instructions_ratio"),
        "all_paired_compact_faster_wall": bool(paired) and all(r["wall_ratio"] < 1 for r in paired),
        "all_targets_verified": all(
            r.get(arm, {}).get("all_targets_verified", True) for r in blocks for arm in ("ic", "ks")
        ),
        "all_scalars_match_corpus": all(
            r.get(key, True) for r in blocks for key in ("ic_matches_corpus", "ks_matches_corpus", "same_target_points")
        ),
        "untimed_records_identical_across_blocks": {
            arm: len({r[arm]["untimed_record_sha256"] for r in blocks if arm in r}) == 1 for arm in ("ic", "ks")
        },
        "replay": replay_summary(out),
        "timing_note": "walls from /usr/bin/time -l on a shared host (load per run recorded); tools/isolated_bench.py unsupported on macOS, so walls are not AGENTS.md section-10 evidence",
    }
    (out / "growing_n_summary.json").write_text(json.dumps(result, indent=1) + "\n")
    print(json.dumps({k: result[k] for k in ("K", "paired_blocks", "wall_ratio", "user_ratio", "instructions_ratio", "all_targets_verified", "all_scalars_match_corpus")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

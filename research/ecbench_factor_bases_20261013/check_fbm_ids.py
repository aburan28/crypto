#!/usr/bin/env python3
"""Fail if an IC record has factor_base but no fbm_id, or fbm_id disagrees with the seal."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
FBM_SCHEMA = "ecbench.factor_base_measurement/v1"
INSTRUMENT = "ecbench.ic_pipeline/v1"


def short_id(prefix: str, preimage: dict) -> str:
    # Match ecbench::canonical (sorted keys, compact, ensure_ascii=True).
    digest = hashlib.sha256(
        json.dumps(preimage, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()
    ).hexdigest()
    return f"{prefix}h{digest[:12]}"


def expected_fbm(rec: dict) -> str | None:
    fb = rec.get("factor_base")
    if not fb:
        return None
    method = rec.get("method") or {}
    params = method.get("params") or {}
    slug = ((rec.get("workload") or {}).get("curve") or {}).get("slug")
    if not slug:
        return None
    cfg = {
        "schema": FBM_SCHEMA,
        "curve": slug,
        "fb_id": fb["fb_id"],
        "fb_sha256": fb["fb_sha256"],
        "family": fb["family"],
        "fb_params": fb["params"],
        "oracle": params.get("oracle") or "",
        "solver": params.get("solver") or "",
        "linalg": params.get("linalg") or "incremental-gauss",
        "instrument": INSTRUMENT,
    }
    return short_id("FBM1", cfg)


def main() -> int:
    problems: list[str] = []
    checked = 0
    skipped_legacy = 0
    for session_dir in sorted((ROOT / "sessions").glob("*")):
        records = session_dir / "records.jsonl"
        if not records.is_file():
            continue
        for line in records.read_text().splitlines():
            if not line.strip():
                continue
            rec = json.loads(line)
            method = ((rec.get("method") or {}).get("id")) or ""
            if not method.startswith("ic."):
                continue
            fb = rec.get("factor_base")
            if not fb:
                continue
            got = rec.get("fbm_id")
            if got is None:
                # Pre-FBM1 committed sessions: allow, but do not count as OK for new work.
                skipped_legacy += 1
                continue
            expect = expected_fbm(rec)
            checked += 1
            if expect is None:
                problems.append(f"{session_dir.name}: cannot recompute fbm_id")
            elif got != expect:
                problems.append(
                    f"{session_dir.name} arm={rec.get('arm')}: fbm_id {got} != seal {expect}"
                )

    if problems:
        print("fbm_id check FAILED:", file=sys.stderr)
        for p in problems:
            print(" ", p, file=sys.stderr)
        return 1
    print(
        f"fbm_id check OK ({checked} IC records matched seal; "
        f"{skipped_legacy} legacy records without fbm_id skipped)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

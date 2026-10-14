#!/usr/bin/env python3
"""Fail if any IC record's fb_id disagrees with the catalog dump for that plug-in."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def main() -> int:
    catalog = json.loads((ROOT / "catalog.json").read_text())
    by_plugin_slug: dict[tuple[str, str], str] = {}
    for e in catalog["entries"]:
        by_plugin_slug[(e["plug_in"], e["slug"])] = e["fb_id"]

    problems: list[str] = []
    checked = 0
    for session_dir in sorted((ROOT / "sessions").glob("*")):
        records = session_dir / "records.jsonl"
        if not records.is_file():
            continue
        for line in records.read_text().splitlines():
            if not line.strip():
                continue
            rec = json.loads(line)
            fb = rec.get("factor_base") or {}
            method = ((rec.get("method") or {}).get("id")) or ""
            if not method.startswith("ic."):
                continue
            params = (rec.get("method") or {}).get("params") or {}
            plugin = params.get("factor_base")
            workload = rec.get("workload") or {}
            slug = ((workload.get("curve") or {}).get("slug")
                    or (rec.get("curve") or {}).get("slug"))
            got = fb.get("fb_id")
            if not plugin or not slug or not got:
                problems.append(f"{session_dir.name}: IC record missing plugin/slug/fb_id")
                continue
            expect = by_plugin_slug.get((plugin, slug))
            if expect is None:
                # arm may have run on a sibling curve in a multi-curve spec
                continue
            checked += 1
            if got != expect:
                problems.append(
                    f"{session_dir.name} arm={rec.get('arm')}: fb_id {got} != catalog {expect} "
                    f"for {plugin} on {slug}"
                )

    if problems:
        print("fb_id check FAILED:", file=sys.stderr)
        for p in problems:
            print(" ", p, file=sys.stderr)
        return 1
    print(f"fb_id check OK ({checked} IC records matched catalog)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

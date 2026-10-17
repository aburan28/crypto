#!/usr/bin/env python3
"""Join Sage factor_base_dor JSON to FB1 / FBM1; refuse orphan DoR cells.

A Sage cell must carry enough identity to resolve a catalog dump:

  - curve slug (ICV1) or a curve constructor matching the catalog
  - points_sha256 and/or fb_id
  - structure / constructor digest when present

Orphan DoR (no matching catalog dump) exits non-zero. On success, writes
FBM1-stamped measurement cells (schema ecbench.factor_base_measurement_cell/v1).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
INSTRUMENT = "sage.factor_base_dor/v1"
FBM_SCHEMA = "ecbench.factor_base_measurement/v1"
CELL_SCHEMA = "ecbench.factor_base_measurement_cell/v1"


def _canonical(obj: object) -> bytes:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()


def short_id(prefix: str, preimage: dict) -> tuple[str, str]:
    digest = hashlib.sha256(_canonical(preimage)).hexdigest()
    return f"{prefix}h{digest[:12]}", digest


def load_catalog() -> dict:
    return json.loads((ROOT / "catalog.json").read_text())


def load_dump(rel: str) -> dict:
    path = ROOT / rel
    if not path.is_file():
        raise FileNotFoundError(f"catalog git dump missing: {rel}")
    return json.loads(path.read_text())


def match_entry(catalog: dict, cell: dict) -> dict:
    slug = cell.get("curve_slug") or cell.get("slug") or (cell.get("curve") or {}).get("slug")
    fb_id = cell.get("fb_id")
    points_sha = cell.get("points_sha256") or (cell.get("factor_base") or {}).get("points_sha256")
    structure = cell.get("structure") or cell.get("constructor") or cell.get("family")
    constructor = cell.get("curve") if isinstance(cell.get("curve"), dict) else None

    candidates = []
    for e in catalog["entries"]:
        if slug and e.get("slug") != slug:
            continue
        if constructor and e.get("curve") != constructor:
            # allow n/a mismatch only when slug already matched
            if not slug:
                continue
        if fb_id and e.get("fb_id") != fb_id:
            continue
        if structure and e.get("plug_in", "").split(":")[0] not in (
            structure,
            structure.replace("_", "-"),
        ):
            # structure names from Sage (standard, random, subfield, …) map loosely;
            # require points_sha256 when structure alone is ambiguous.
            if not points_sha and not fb_id:
                continue
        dump = load_dump(e["git"]) if e.get("git") else None
        if points_sha and dump:
            got = dump["factor_base"]["points_sha256"]
            if got != points_sha:
                continue
        if dump and fb_id and dump["factor_base"]["fb_id"] != fb_id:
            continue
        candidates.append((e, dump))

    if not candidates:
        raise SystemExit(
            "orphan DoR refused: no catalog FB1 dump matches "
            f"slug={slug!r} fb_id={fb_id!r} points_sha256={points_sha!r} structure={structure!r}"
        )
    if len(candidates) > 1 and not (fb_id or points_sha):
        raise SystemExit(
            "orphan DoR refused: ambiguous catalog match; export curve slug + "
            "fb_id or points_sha256 with the Sage cell"
        )
    return candidates[0]


def seal_fbm(entry: dict, dump: dict, cell: dict) -> dict:
    fb = dump["factor_base"]
    oracle = cell.get("oracle") or entry.get("oracle") or "sage.factor_base_dor"
    solver = cell.get("solver") or "sage.groebner"
    linalg = cell.get("linalg") or "sage.macaulay"
    cfg = {
        "schema": FBM_SCHEMA,
        "curve": dump["curve"]["slug"],
        "fb_id": fb["fb_id"],
        "fb_sha256": fb["fb_sha256"],
        "family": fb["family"],
        "fb_params": fb["params"],
        "oracle": oracle,
        "solver": solver,
        "linalg": linalg,
        "instrument": INSTRUMENT,
    }
    fbm_id, fbm_sha = short_id("FBM1", cfg)
    return {"fbm_id": fbm_id, "fbm_sha256": fbm_sha, "config": cfg}


def cell_out(sealed: dict, sage_cell: dict) -> dict:
    extra = {
        k: sage_cell[k]
        for k in (
            "structure",
            "d_star",
            "d_star_mean",
            "d_star_max",
            "decomposable_fraction",
            "field",
            "dim",
            "targets",
            "dmax",
        )
        if k in sage_cell
    }
    return {
        "schema": CELL_SCHEMA,
        "fbm_id": sealed["fbm_id"],
        "fbm_sha256": sealed["fbm_sha256"],
        "config": sealed["config"],
        "solving_degree_mean": sage_cell.get("d_star_mean") or sage_cell.get("d_star"),
        "solving_degree_max": sage_cell.get("d_star_max") or sage_cell.get("d_star"),
        "status": sage_cell.get("status", "sage"),
        "extra": extra,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("sage_json", type=Path, help="Sage factor_base_dor output JSON")
    ap.add_argument("--out", type=Path, help="write joined FBM1 cells here")
    args = ap.parse_args()

    raw = json.loads(args.sage_json.read_text())
    cells = raw if isinstance(raw, list) else raw.get("cells") or raw.get("results") or [raw]
    catalog = load_catalog()
    joined = []
    for cell in cells:
        if not isinstance(cell, dict):
            continue
        entry, dump = match_entry(catalog, cell)
        if dump is None:
            raise SystemExit(f"catalog entry {entry.get('id')} has no git dump to join")
        sealed = seal_fbm(entry, dump, cell)
        joined.append(cell_out(sealed, cell))

    if not joined:
        print("no Sage cells to join", file=sys.stderr)
        return 1

    text = json.dumps({"schema": "ecbench.sage_dor_join/v1", "cells": joined}, indent=2) + "\n"
    if args.out:
        args.out.write_text(text)
        print(f"joined {len(joined)} Sage DoR cell(s) -> {args.out}", file=sys.stderr)
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

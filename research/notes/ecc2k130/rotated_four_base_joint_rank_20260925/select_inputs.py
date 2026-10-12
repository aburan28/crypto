#!/usr/bin/env python3
"""Freeze one unconditional same-Q SHA training/holdout order; no archive read."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import resource
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ARITHMETIC = HERE.parent / "rotated_subspace_support_20260925/gate.py"
DOMAIN = "ECC2K130-ROTATED-JOINT-RANK-20260925-v1"
Q = 130873
H = (385982, 301867)
TRAIN = 256
HOLDOUT = 64


def save(path: Path, value) -> None:
    path.write_text(json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rss() -> int:
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return raw if sys.platform == "darwin" else raw * 1024


def arithmetic():
    spec = importlib.util.spec_from_file_location("joint_input_arithmetic", ARITHMETIC)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def score(k: int) -> bytes:
    return hashlib.sha256(f"{DOMAIN}/{k}".encode("ascii")).digest()


def generate():
    start = time.monotonic()
    mod = arithmetic()
    field = mod.Field(19, [0, 1, 2, 5])
    field.rabin_prime_degree()
    curve = mod.Curve(field)
    assert field.poly == 0x80027 and mod.source_group_order(19) == 4 * Q
    assert curve.on_curve(H) and curve.scalar(H, Q) is None
    ordered = sorted(range(1, Q), key=lambda k: (score(k), k))
    assert len(ordered) == Q - 1 and len(set(ordered)) == Q - 1
    chosen = ordered[:TRAIN + HOLDOUT]
    assert len(chosen) == TRAIN + HOLDOUT and len(set(chosen)) == len(chosen)
    points = [list(curve.scalar(H, k)) for k in chosen]
    assert all(curve.on_curve(tuple(p)) for p in points)
    training = {"schema": "rotated_joint_training_v1", "domain": DOMAIN,
                "q": Q, "generator": list(H),
                "targets": [{"case_id": f"tr-{i:03d}", "k": k, "point": p}
                            for i, (k, p) in enumerate(zip(chosen[:TRAIN], points[:TRAIN]))]}
    point_only = {"schema": "rotated_joint_point_only_v1", "domain": DOMAIN,
                  "q": Q, "generator": list(H),
                  "targets": [{"case_id": f"ho-{i:03d}", "point": p}
                              for i, p in enumerate(points[TRAIN:])]}
    sealed = {"schema": "rotated_joint_sealed_v1", "domain": DOMAIN,
              "targets": [{"case_id": f"ho-{i:03d}", "k": k}
                          for i, k in enumerate(chosen[TRAIN:])]}
    receipt = {"status": "success", "domain": DOMAIN, "q": Q,
               "training_count": TRAIN, "holdout_count": HOLDOUT,
               "first_320_score_sha256": hashlib.sha256(b"".join(score(k) for k in chosen)).hexdigest(),
               "wall_seconds": time.monotonic() - start,
               "peak_rss_bytes": rss(),
               "field_operations": dict(field.operations),
               "curve_operations": dict(curve.operations)}
    assert receipt["wall_seconds"] <= 60 and receipt["peak_rss_bytes"] <= 512 * 1024 * 1024
    return training, point_only, sealed, receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists()
    args.out.mkdir(parents=True)
    training, point_only, sealed, receipt = generate()
    for name, value in (("training.json", training), ("point_only.json", point_only),
                        ("sealed_labels.json", sealed)):
        save(args.out / name, value)
    receipt["input_sha256"] = {name: sha(args.out / name)
                               for name in ("training.json", "point_only.json", "sealed_labels.json")}
    save(args.out / "construction_receipt.json", receipt)


if __name__ == "__main__":
    main()

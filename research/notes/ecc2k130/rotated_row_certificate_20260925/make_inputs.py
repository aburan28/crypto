#!/usr/bin/env python3
"""Freeze deterministic, public synthetic n131 planted tuples; no PDP search."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import resource
import signal
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / "rotated_subspace_support_20260925" / "gate.py"
DOMAIN = "ECC2K130-ROTATED-ROW-20260925-v1"
CELLS = ((5, 25), (6, 21))


def prior():
    spec = importlib.util.spec_from_file_location("frozen_rotated_gate", OLD)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def digest_label(m: int, d: int, t: int, i: int, counter: int) -> str:
    return f"{DOMAIN}/{m}/{d}/{t}/{i}/{counter}"


def source_point(mod, curve, basis, m: int, d: int, t: int, i: int,
                 forced_mask: int | None = None, forced_counter: int | None = None):
    f = curve.f
    if t == 1 and i == 0:
        return {"slot": i, "mask": 0, "counter": None, "rejected": [],
                "sign_bit": 0, "base_point": [0, 1], "source_point": [0, 1]}
    rejected = []
    for counter in range(128):
        if forced_mask is not None:
            assert forced_counter is not None
            counter = forced_counter
            mask = forced_mask
        else:
            seed = digest_label(m, d, t, i, counter).encode("ascii")
            mask = int.from_bytes(hashlib.sha256(seed).digest(), "big") % (1 << d)
        if mask == 0:
            rejected.append({"counter": counter, "reason": "zero"})
        else:
            x = mod.x_from_mask(basis, mask)
            assert x != 0
            rhs = x ^ f.square(f.inverse(x))
            if f.trace(rhs) != 0:
                rejected.append({"counter": counter, "reason": "no_rational_lift"})
            else:
                z = mod.half_trace(f, rhs)
                p0 = (x, f.mul(x, z))
                assert curve.on_curve(p0)
                sign_source = digest_label(m, d, t, i, counter) + "/sign"
                sign_bit = hashlib.sha256(sign_source.encode("ascii")).digest()[-1] & 1
                signed = curve.neg(p0) if sign_bit else p0
                p = signed
                for _ in range(i):
                    p = curve.tau(p)
                assert curve.on_curve(p)
                return {"slot": i, "mask": mask, "counter": counter,
                        "rejected": rejected, "sign_bit": sign_bit,
                        "base_point": list(signed), "source_point": list(p)}
        if forced_mask is not None:
            raise AssertionError("forced shared mask not rational")
    raise AssertionError((m, d, t, i, "no rational mask within 128 counters"))


def peak_rss_bytes() -> int:
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return raw if sys.platform == "darwin" else raw * 1024


def deadline(_signum, _frame):
    raise TimeoutError("n131 input cell exceeded hard 30-second cap")


def create_inputs() -> tuple[dict, list[dict]]:
    mod = prior()
    f = mod.Field(131, mod.MODELS[131]["low"])
    f.rabin_prime_degree()
    curve = mod.Curve(f)
    assert mod.source_group_order(131) == 4 * mod.Q131
    conjugates = mod.normal_conjugates(f, 3)
    torsion = mod.four_torsion(curve)
    cells = []
    costs = []
    inverse_four = pow(4, -1, mod.Q131)
    signal.signal(signal.SIGALRM, deadline)
    for m, d in CELLS:
        wall_start, cpu_start = time.monotonic(), time.process_time()
        before_f, before_c = dict(f.operations), dict(curve.operations)
        signal.setitimer(signal.ITIMER_REAL, 30)
        bases = mod.subspace_basis(conjugates, m, d)
        cell = {"m": m, "d": d, "tuples": []}
        for t in range(4):
            slots = []
            for i in range(m):
                if t == 0 and i > 0:
                    first = slots[0]
                    slot = source_point(mod, curve, bases[0], m, d, t, i,
                                        forced_mask=first["mask"],
                                        forced_counter=first["counter"])
                else:
                    slot = source_point(mod, curve, bases[0], m, d, t, i)
                assert mod.x_from_mask(bases[i], slot["mask"]) == slot["source_point"][0]
                slots.append(slot)
            total = None
            for slot in slots:
                total = curve.add(total, tuple(slot["source_point"]))
            projected = curve.scalar(total, 4)
            q = curve.scalar(projected, inverse_four)
            tors = curve.add(total, curve.neg(q))
            assert tors in torsion and curve.scalar(q, mod.Q131) is None
            assert curve.add(q, tors) == total
            cell["tuples"].append({"tuple": t, "slots": slots, "raw_sum": mod.pjson(total),
                                   "synthetic_Q": mod.pjson(q), "torsion_index": torsion.index(tors),
                                   "projected_Q": mod.pjson(curve.scalar(q, 4))})
        cells.append(cell)
        signal.setitimer(signal.ITIMER_REAL, 0)
        cost = {"m": m, "d": d, "wall_seconds": time.monotonic() - wall_start,
                "cpu_seconds": time.process_time() - cpu_start,
                "peak_rss_bytes": peak_rss_bytes(),
                "operations": {key: f.operations.get(key, 0) - before_f.get(key, 0) +
                               curve.operations.get(key, 0) - before_c.get(key, 0)
                               for key in sorted(set(f.operations) | set(curve.operations))}}
        assert cost["peak_rss_bytes"] <= 512 * 1024 * 1024
        costs.append(cost)
    return {"domain": DOMAIN, "kind": "public_synthetic_planted_row_inputs",
            "n": 131, "polynomial": f.poly, "q": mod.Q131,
            "lambda": mod.LAMBDA131, "beta": 3,
            "torsion": [mod.pjson(t) for t in torsion], "cells": cells}, costs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--receipt", required=True, type=Path)
    args = parser.parse_args()
    assert not args.out.exists() and not args.receipt.exists()
    wall_start, cpu_start = time.monotonic(), time.process_time()
    receipt = {"kind": "public_synthetic_input_construction", "status": "failed"}
    try:
        obj, costs = create_inputs()
        raw = json.dumps(obj, sort_keys=True, separators=(",", ":")) + "\n"
        args.out.write_text(raw)
        receipt.update({"status": "success", "cells": costs,
                        "input_sha256": hashlib.sha256(raw.encode()).hexdigest()})
    except Exception as exc:
        receipt.update({"error": repr(exc), "traceback": traceback.format_exc()})
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        receipt.update({"wall_seconds": time.monotonic() - wall_start,
                        "cpu_seconds": time.process_time() - cpu_start,
                        "peak_rss_bytes": peak_rss_bytes()})
        args.receipt.write_text(json.dumps(receipt, sort_keys=True, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    main()

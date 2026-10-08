#!/usr/bin/env python3
"""A check made after R04's run, not declared in PROTOCOL.md: is the
presence filter's admitted fraction what its sizing predicts?

    tar -xJf runs.tar.xz && python3 filter_check.py > filter_check.json

The folded pair table's presence filter is one bit per stored pair under
one hash, with `2^bits` bits, where `bits` is the bit length of four times
the stored pairs, clamped to 6..32 (`PairSumTable`'s builders in
`src/cryptanalysis/koblitz_index_calculus.rs`). A key that is not in the
table passes with the fraction of bits set, `1 - (1 - 2^-bits)^stored`.
True hits are a few in a hundred thousand summands, so the predicted
admitted fraction is that rate. The measured fraction is analyse.py's:
the probe arm's admitted keys over its scanned summands, the median over
processes. The admitted stage's nanoseconds per admitted key use the
same processes and the counter rate each one measured.
"""
from __future__ import annotations

import json
import math
import os
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "harness"))
import bench  # noqa: E402

RUNS = (HERE / os.environ.get("IC_RUNS", "runs")).resolve()


def filter_bits(stored: int) -> int:
    return min(max((max(stored, 1) * 4).bit_length(), 6), 32)


def size_row(d: Path, rs: list[dict]) -> dict:
    stored, measured, hits, per_key = set(), [], [], []
    for r in rs:
        for path in sorted((d / "probes" / r["id"]).glob("r*.price.json")):
            if "-retry" in path.name:
                continue
            p = bench.figure_path(path)
            rep = json.loads(p.read_text()) if p.exists() else None
            if rep is None or rep.get("status") != "complete" or not bench.clean(p) or "scan_probes" not in rep:
                continue
            probes = rep["scan_probes"]
            counts = probes["counts"]
            stored.add(rep["counts"]["pass"]["build"]["stored_pairs"])
            measured.append(counts["admitted"] / counts["summands"])
            hits.append(counts["pairs"] / counts["summands"])
            per_key.append(probes["stages"]["admitted"] * 1e9 / probes["tsc_hz"] / counts["admitted"])
    if len(stored) != 1:
        raise SystemExit(f"stored pairs differ within a size: {sorted(stored)}")
    n = stored.pop()
    bits = filter_bits(n)
    predicted = 1 - (1 - 2.0 ** -bits) ** n
    got = statistics.median(measured)
    return {"stored_pairs": n, "filter_bits": bits, "filter_bits_per_stored_pair": round(2 ** bits / n, 3),
            "predicted_pass_rate": round(predicted, 4), "measured_admitted_per_summand": round(got, 4),
            "measured_over_predicted": round(got / predicted, 4),
            "true_hits_per_summand": statistics.median(hits),
            "admitted_stage_ns_per_admitted_key": round(statistics.median(per_key), 1), "processes": len(measured)}


def main() -> None:
    rows = [r for r in bench.slug_rows(bench.suite_rows("S")) if r["recipe_seed"] == 201]
    by: dict[tuple[int, int], list[dict]] = {}
    for r in rows:
        by.setdefault((r["a"], r["n"]), []).append(r)
    sizes = []
    for (a, n), rs in sorted(by.items(), key=lambda kv: kv[1][0]["r"]):
        row = {"slug": bench.curve_slug(a, n), "log2_r": round(math.log2(rs[0]["r"]), 3)}
        row.update(size_row(RUNS / "compare", rs))
        sizes.append(row)
    ratios = [s["measured_over_predicted"] for s in sizes]
    doc = {
        "what_this_is": "a check after R04's run, not declared: the presence filter's predicted pass rate "
                        "against the admitted fraction R04 measured",
        "sizes": sizes,
        "measured_over_predicted_range": [min(ratios), max(ratios)],
    }
    json.dump(doc, sys.stdout, indent=1, sort_keys=False)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()

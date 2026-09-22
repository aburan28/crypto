#!/usr/bin/env python3
"""Flatten `isogeny experiment` JSON reports into one per-vertex CSV.

The committed `experiments/pooled_vertex_data.csv` was originally produced
by hand, which made it impossible to regenerate when the isogeny code
changed.  This script makes it reproducible:

    python3 scripts/pool_isogeny_vertices.py \
        experiments/06_experiment_14bit_10trials.json \
        experiments/06_experiment_16bit_10trials.json \
        experiments/07_experiment_18bit_10trials.json \
        experiments/13_22bit_dpoly.json \
        experiments/14b_30bit.json \
        experiments/15b_40bit.json \
        experiments/16c_50bit.json \
        > experiments/pooled_vertex_data.csv

One row per vertex of every curve's isogeny class, tagged with the bit
width and the index of the starting curve it came from.

Pass at most ONE artifact per bit width: several widths have runs at two
rho caps (40-bit at 2^18 and 2^23, 50-bit at 2^18 and 2^28) over the same
vertices, and pooling both counts each vertex twice.  The caps differ
across widths -- the cap has to scale with sqrt(n) to be reachable at all
-- so `rho_cap` is carried per row rather than left implicit; a success
rate is only comparable against another row at the same cap.
"""
import csv
import json
import sys

COLUMNS = [
    "bits", "rho_cap", "curve_idx", "class_size", "class_capped", "p", "a",
    "b", "trace", "n", "n_parity", "fundamental_disc", "conductor",
    "rho_iters", "rho_success", "mov_feasible", "smart_applies",
    "glv_speedup",
]


def rows(report):
    cfg = report["config"]
    bits = cfg["bits"]
    cap = cfg["rho_max_iters"]
    # A class that exactly reaches max_graph_nodes was truncated by the
    # walk's node budget, so its size is a censored lower bound, not a
    # measurement.  At 50 bits both classes hit it.
    node_cap = cfg["max_graph_nodes"]
    for idx, curve in enumerate(report["per_curve"]):
        for vertex in curve["class"]:
            p, a, b = vertex["curve"]
            cm = vertex["cm"]
            yield {
                "bits": bits,
                "rho_cap": cap,
                "curve_idx": idx,
                "class_size": curve["class_size"],
                "class_capped": int(curve["class_size"] >= node_cap),
                "p": p,
                "a": a,
                "b": b,
                "trace": cm["trace"],
                "n": cm["order"],
                "n_parity": cm["order"] % 2,
                "fundamental_disc": cm["fundamental_disc"],
                "conductor": cm["conductor"],
                "rho_iters": vertex["rho_iters"],
                "rho_success": int(vertex["rho_success"]),
                "mov_feasible": int(vertex["mov_feasible"]),
                "smart_applies": int(vertex["smart_applies"]),
                "glv_speedup": int(vertex["glv_speedup_available"]),
            }


def main(paths):
    if not paths:
        sys.exit(__doc__)
    writer = csv.DictWriter(sys.stdout, fieldnames=COLUMNS, lineterminator="\n")
    writer.writeheader()
    for path in paths:
        with open(path) as handle:
            writer.writerows(rows(json.load(handle)))


if __name__ == "__main__":
    main(sys.argv[1:])

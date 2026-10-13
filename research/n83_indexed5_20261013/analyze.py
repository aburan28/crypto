#!/usr/bin/env python3
"""Independent generic N83 group and archived-base replay for indexed5 runs."""

import json
import sys
from pathlib import Path

N = 83
MOD = (1 << N) | (1 << 7) | (1 << 4) | (1 << 2) | 1
R = 8569786107849059
LAMBDA = 2998004823950457
G = (6614244651398184860292015, 8868260410029249391539529)


def product(a, b):
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a & (1 << N):
            a ^= MOD
    return out


def inverse(a):
    if a == 0:
        raise ValueError("division by zero")
    u, v, g, h = a, MOD, 1, 0
    while u != 1:
        shift = u.bit_length() - v.bit_length()
        if shift < 0:
            u, v = v, u
            g, h = h, g
            shift = -shift
        u ^= v << shift
        g ^= h << shift
    while g.bit_length() > N:
        g ^= MOD << (g.bit_length() - N - 1)
    return g


def divide(a, b):
    return product(a, inverse(b))


def on_curve(point):
    if point is None:
        return True
    x, y = point
    return product(y, y) ^ product(x, y) == product(product(x, x), x) ^ product(x, x) ^ 1


def add(left, right):
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2:
        if (y1 ^ y2) == x1 or x1 == 0:
            return None
        slope = x1 ^ divide(y1, x1)
        x3 = product(slope, slope) ^ slope ^ 1
        return (x3, product(x1, x1) ^ product(slope ^ 1, x3))
    slope = divide(y1 ^ y2, x1 ^ x2)
    x3 = product(slope, slope) ^ slope ^ x1 ^ x2 ^ 1
    return (x3, product(slope, x1 ^ x3) ^ x3 ^ y1)


def multiply(point, scalar):
    acc = None
    while scalar:
        if scalar & 1:
            acc = add(acc, point)
        point = add(point, point)
        scalar >>= 1
    return acc


def point_from_json(value):
    if value is None:
        return None
    point = tuple(int(v) for v in value)
    if len(point) != 2 or any(v < 0 or v >= 1 << N for v in point):
        raise ValueError("invalid N83 point encoding")
    return point


def one_json(path):
    lines = Path(path).read_text().splitlines()
    if len(lines) != 1:
        raise ValueError(f"expected exactly one row: {path}")
    return json.loads(lines[0])


def replay_base(dump, summary):
    coords = [point_from_json(v) for v in dump["factor_base_point_coordinates"]]
    labels = dump["factor_base_point_labels"]
    reps = [point_from_json(v) for v in dump["factor_base_representatives"]]
    if (len(coords) != len(labels) or len(coords) != 166 * len(reps)
            or len(coords) != summary["factor_base_points"]
            or len(reps) != summary["orbit_columns"]
            or dump["base_hash"] != summary["base_hash"]):
        raise ValueError("archived base identity or dimensions mismatch")
    seen = set()
    for column, rep in enumerate(reps):
        if rep is None or not on_curve(rep) or multiply(rep, R) is not None:
            raise ValueError(f"representative {column} fails curve/subgroup check")
        current, coefficient = rep, 1
        for rotation in range(N):
            positive = 166 * column + 2 * rotation
            negative = positive + 1
            if coords[positive] != current or coords[negative] != (current[0], current[0] ^ current[1]):
                raise ValueError(f"orbit mismatch at {column}/{rotation}")
            if labels[positive] != [column, coefficient] or labels[negative] != [column, (R - coefficient) % R]:
                raise ValueError(f"label mismatch at {column}/{rotation}")
            for item in coords[positive:negative + 1]:
                if item in seen:
                    raise ValueError(f"duplicate base point at {column}/{rotation}")
                seen.add(item)
            current = (product(current[0], current[0]), product(current[1], current[1]))
            coefficient = coefficient * LAMBDA % R
        if current != rep or coefficient != 1:
            raise ValueError(f"orbit does not close at {column}")
    return coords


def verify_arm(study, phase, arm, answer):
    root = study / "runs" / phase / arm
    if not root.exists():
        return {"arm": arm, "status": "NOT_RUN"}
    if not (root / "exit_status.txt").exists():
        return {"arm": arm, "status": "INTERRUPTED_NO_EXIT_RECEIPT"}
    exit_status = int((root / "exit_status.txt").read_text().strip())
    usage = one_json(root / "usage.json")
    # A shell reports a signalled child's negative return code modulo 256.
    if usage["exit_status"] % 256 != exit_status or usage["peak_rss_unit"] != "bytes":
        raise ValueError(f"{arm}: direct-child usage mismatch")
    if exit_status:
        status = ("TIMEOUT" if usage["timed_out"] and exit_status == 124 else
                  "STOPPED_RESOURCE_GUARD" if usage["exit_status"] == -15 else "FAILED")
        return {"arm": arm, "status": status,
                "exit_status": exit_status, "usage": usage}
    summary = one_json(root / "stdout.jsonl")
    target = one_json(root / "targets_result.jsonl")
    if (summary["n"] != N or summary["a"] != 1 or summary["targets"] != 1
            or summary["targets_solved"] != 1 or summary["targets_failed"] != 0
            or summary["rank"] != summary["orbit_columns"]
            or summary["root_table_layout"] != "indexed83_5byte_exact"
            or summary["root_table_slot_bytes"] != 5
            or summary["root_table_indexed5_enabled"] is not True
            or summary["root_table_dense_enabled"] is not False
            or summary["root_table_load75_enabled"] is not False
            or summary["root_table_allocated_bytes"] > 10 * 1024 ** 3):
        raise ValueError(f"{arm}: rank, target or indexed-table boundary mismatch")
    if (summary["root_state_layout"] != "packed83_11byte_roots"
            or summary["root_state_slot_bytes"] != 28
            or summary["root_state_allocated_bytes"] < summary["regular_states"] * summary["root_state_slot_bytes"]):
        raise ValueError(f"{arm}: state layout receipt mismatch")
    q = point_from_json(answer["point"])
    scalar = int(target["recovered_scalar"])
    if (point_from_json(target["target"]) != q or target["fixture_index"] != 0
            or scalar != int(answer["scalar"]) or not target["group_verified"]
            or multiply(G, scalar) != q):
        raise ValueError(f"{arm}: independent scalar replay failed")
    if not on_curve(q) or multiply(q, R) is not None:
        raise ValueError(f"{arm}: public point fails curve/subgroup check")
    dump = one_json(root / "base_dump.jsonl")
    base = replay_base(dump, summary)
    indices = target["point_indices"]
    if len(indices) != 4:
        raise ValueError(f"{arm}: relation is not four-point")
    relation = None
    for index in indices:
        relation = add(relation, base[index])
    if relation != q:
        raise ValueError(f"{arm}: four-point relation replay failed")
    online = float(target["target_ms"])
    if abs(online - float(target["target_phase_sum_ms"])) > 1.0:
        raise ValueError(f"{arm}: online phase sum mismatch")
    return {"arm": arm, "status": "PASS_INDEPENDENT_REPLAY",
            "factor_base_points": len(base), "orbit_columns": len(dump["factor_base_representatives"]),
            "base_hash": summary["base_hash"], "regular_states": summary["regular_states"],
            "root_state_allocated_bytes": summary["root_state_allocated_bytes"],
            "root_table_allocated_bytes": summary["root_table_allocated_bytes"],
            "root_table_entries": summary["root_table_entries"],
            "rank": summary["rank"], "target_probes": target["probes"],
            "online_ms_exploratory": online, "setup_ms": summary["timing_ms"]["setup_total"],
            "rank_ms": summary["timing_ms"]["rank_stage"],
            "linear_algebra_ms": summary["timing_ms"]["linear_algebra"],
            "peak_rss_bytes": usage["peak_rss"], "point_indices": indices,
            "x_codes": target["x_codes"], "recovered_scalar": scalar}


def main():
    if len(sys.argv) != 3 or sys.argv[2] not in ("control", "primary"):
        raise SystemExit("usage: analyze.py STUDY_DIR control|primary")
    study, phase = Path(sys.argv[1]), sys.argv[2]
    answer = one_json(study / "answer.jsonl")
    arms = ("k512_indexed5",) if phase == "control" else ("k2500_indexed5",)
    rows = [verify_arm(study, phase, arm, answer) for arm in arms]
    parity = None
    if phase == "control" and rows[0]["status"] == "PASS_INDEPENDENT_REPLAY":
        prior = one_json(study / "../n83_state28_20261013/analysis_control.json")
        prior_packed = next(row for row in prior["arms"] if row["arm"] == "k512_packed")
        if prior_packed["status"] != "PASS_INDEPENDENT_REPLAY":
            raise ValueError("archived K512 packed control is unverified")
        fields = ("base_hash", "factor_base_points", "orbit_columns", "regular_states",
                  "point_indices", "x_codes", "target_probes", "recovered_scalar")
        parity = all(rows[0][field] == prior_packed[field] for field in fields)
        if not parity:
            raise ValueError("indexed5/archived-packed first-witness parity failed")
    result = {"schema": "n83.indexed5-analysis/v1", "phase": phase, "arms": rows,
              "first_witness_parity": parity}
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()

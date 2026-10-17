"""Reduce a bounded slice of the frozen new Type I quotient points.

The two previously checked constructors are loaded only through their
definition/preparation portions.  Their archived 141-class output loops are
not run and their receipts are not rewritten.  This wrapper freezes a new
slice, its cubic input divisors, and its fixed-origin representatives.
"""

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_prelude(name):
    path = HERE / name
    source = path.read_text()
    marker = "start = time.monotonic_ns()"
    assert source.count(marker) == 1
    namespace = {"__file__": str(path), "__name__": name.removesuffix(".py") + "_prelude"}
    exec(compile(source.split(marker)[0], str(path), "exec"), namespace)
    return namespace, source


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--count", type=int, default=64)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    input_path = HERE / "type_i_line_full_new_quotient_points.json"
    points = json.loads(input_path.read_text())
    assert len(points) == 9784
    assert 0 <= args.start < len(points) and 0 < args.count <= len(points) - args.start

    quotient, _ = load_prelude("line_quotient_cubic_probe.py")
    reduction, reduction_source = load_prelude("line_quotient_fixed_origin.py")
    density_path = HERE / "type_i_character_density.json"
    density = json.loads(density_path.read_text())
    coefficients = [int(value) for value in density["sextic_coefficients_decimal"]]

    def character(point):
        value = 0
        for coefficient in reversed(coefficients):
            value = (value * int(point[0].polynomial()[0]) + coefficient) % int(quotient["p"])
        assert value != 0
        return "+" if pow(value, (int(quotient["p"]) - 1) // 2, int(quotient["p"])) == 1 else "-"

    E, L, Integer = quotient["E"], quotient["L"], quotient["Integer"]
    r, O, inv2 = quotient["r"], quotient["O"], quotient["inv2"]
    rows = []
    start_ns = time.monotonic_ns()
    for index in range(args.start, args.start + args.count):
        xx, yy = points[index]
        Q = E(L(Integer(xx)), L(Integer(yy)))
        assert Q != O and r * Q == O
        H = inv2 * Q
        assert 2 * H == Q
        rows.append({
            "index": index,
            "cell": character(Q) + character(H),
            "quotient": [str(xx), str(yy)],
            "positive_input_residual": quotient["residual"](Q),
            "negative_input_residual": quotient["residual"](H),
            "status": "verified_prime_field_cubic_pair",
        })
        if len(rows) % 16 == 0:
            print("NEW_QUOTIENT_INPUT_PROGRESS", len(rows), args.count, flush=True)

    reduction["source"] = {"rows": rows}
    reduction["representatives"] = []
    loop = reduction_source[
        reduction_source.index('for item in source["rows"]:'):
        reduction_source.index("\nsuccesses =")
    ]
    exec(compile(loop, str(HERE / "line_quotient_fixed_origin.py"), "exec"), reduction)
    representatives = reduction["representatives"]
    successes = sum(row["status"] == "verified_fixed_origin_degree_three"
                    for row in representatives)
    output = {
        "status": "bounded_new_quotient_fixed_origin_reduction",
        "field_prime": str(quotient["p"]),
        "subgroup_order": str(r),
        "start_index": args.start,
        "input_classes": args.count,
        "reduced_classes": successes,
        "unreduced_classes": args.count - successes,
        "origin_base_indices": [0, 1, 2],
        "origin_x_polynomial": reduction["encode_poly"](reduction["origin_h"]),
        "origin_y_polynomial": reduction["encode_poly"](reduction["origin_y"]),
        "class_identity": "G_cut=B+T+U and F_cut=A+U+D imply [D-T]=[B-A]",
        "rows": rows,
        "representatives": representatives,
        "elapsed_ns_unisolated": time.monotonic_ns() - start_ns,
        "source_sha256": digest(Path(__file__)),
        "input_sha256": digest(input_path),
        "quotient_constructor_sha256": digest(HERE / "line_quotient_cubic_probe.py"),
        "reduction_constructor_sha256": digest(HERE / "line_quotient_fixed_origin.py"),
        "plane_receipt_sha256": digest(HERE / "type_i_plane_quartic.json"),
        "density_receipt_sha256": digest(density_path),
        "boundary": "Exact fixed-origin cubic cuts for the declared new quotient slice; independent native replay and larger slices are separate checks.",
    }
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print("NEW_QUOTIENT_REDUCTION_PASS", json.dumps({
        "start_index": args.start, "input_classes": args.count,
        "reduced_classes": successes, "unreduced_classes": args.count - successes,
        "output_sha256": digest(args.out),
    }, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()

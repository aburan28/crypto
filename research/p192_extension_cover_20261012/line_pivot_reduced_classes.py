"""Reduce a frozen slice of exact pivot-fan quotient classes on the quartic."""

import hashlib
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prelude(name):
    path = HERE / name
    source = path.read_text()
    marker = "start = time.monotonic_ns()"
    assert source.count(marker) == 1
    namespace = {"__file__": str(path), "__name__": name[:-3] + "_prelude"}
    exec(compile(source.split(marker)[0], str(path), "exec"), namespace)
    return namespace, source


def main():
    inputs = [HERE / f"type_i_residual_pivot{i:03d}_{'three' if i == 0 else 'cross'}.json"
              for i in (0, 1)]
    selected = set()
    samples = []
    for pivot, path in enumerate(inputs):
        fan = json.loads(path.read_text())
        assert fan.get("pivot_index", 0) == pivot
        count = 0
        for row_index, row in enumerate(fan["rows"]):
            for side, coordinates in enumerate(row["projected_residuals"]):
                key = tuple(map(int, coordinates))
                if key in selected:
                    continue
                selected.add(key)
                samples.append((pivot, row_index, side, key))
                count += 1
                if count == 32:
                    break
            if count == 32:
                break
        assert count == 32

    quotient, _ = prelude("line_quotient_cubic_probe.py")
    reduction, reduction_source = prelude("line_quotient_fixed_origin.py")
    density_path = HERE / "type_i_character_density.json"
    density = json.loads(density_path.read_text())
    coefficients = [int(value) for value in density["sextic_coefficients_decimal"]]
    E, L, Integer = quotient["E"], quotient["L"], quotient["Integer"]
    r, O, inv2 = quotient["r"], quotient["O"], quotient["inv2"]

    def character(point):
        x = int(point[0].polynomial()[0])
        value = 0
        for coefficient in reversed(coefficients):
            value = (value*x + coefficient) % int(quotient["p"])
        assert value != 0
        return "+" if pow(value, (int(quotient["p"]) - 1)//2,
                           int(quotient["p"])) == 1 else "-"

    rows = []
    start = time.monotonic_ns()
    for index, (pivot, row_index, side, (xx, yy)) in enumerate(samples):
        Q = E(L(Integer(xx)), L(Integer(yy)))
        assert Q != O and r * Q == O
        H = inv2 * Q
        assert 2 * H == Q
        rows.append({
            "index": index, "pivot": pivot, "fan_row_index": row_index,
            "residual_side": side, "quotient": [str(xx), str(yy)],
            "cell": character(H) + character(Q),
            "positive_input_residual": quotient["residual"](Q),
            "negative_input_residual": quotient["residual"](H),
            "status": "verified_prime_field_cubic_pair",
        })
        if (index + 1) % 8 == 0:
            print("PIVOT_CUBIC_PROGRESS", index + 1, len(samples), flush=True)

    reduction["source"] = {"rows": rows}
    reduction["representatives"] = []
    loop = reduction_source[
        reduction_source.index('for item in source["rows"]:'):
        reduction_source.index("\nsuccesses =")
    ]
    exec(compile(loop, str(HERE / "line_quotient_fixed_origin.py"), "exec"), reduction)
    representatives = reduction["representatives"]
    assert all(row["status"] == "verified_fixed_origin_degree_three"
               for row in representatives)
    output = {
        "status": "verified_pivot_fan_fixed_origin_slice",
        "field_prime": str(quotient["p"]), "subgroup_order": str(r),
        "classes_per_fan": 32, "input_classes": len(samples),
        "reduced_classes": len(representatives), "rows": rows,
        "representatives": representatives,
        "origin_x_polynomial": reduction["encode_poly"](reduction["origin_h"]),
        "origin_y_polynomial": reduction["encode_poly"](reduction["origin_y"]),
        "class_identity": "G_cut=B+T+U and F_cut=A+U+D imply [D-T]=[B-A]",
        "elapsed_ns_unisolated": time.monotonic_ns() - start,
        "source_sha256": sha(__file__),
        "fan_sha256": [sha(path) for path in inputs],
        "quotient_constructor_sha256": sha(HERE / "line_quotient_cubic_probe.py"),
        "reduction_constructor_sha256": sha(HERE / "line_quotient_fixed_origin.py"),
        "plane_receipt_sha256": sha(HERE / "type_i_plane_quartic.json"),
        "density_receipt_sha256": sha(density_path),
        "boundary": "First 32 distinct exact quotient points from each pivot fan. Cubic sections certify fixed-origin degree-three classes; the full fans and independent native replay remain separate gates.",
    }
    path = HERE / "type_i_two_pivot_reduced64.json"
    path.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print("PIVOT_REDUCED_PASS", len(representatives), sha(path), flush=True)


if __name__ == "__main__":
    main()

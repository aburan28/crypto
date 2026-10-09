"""Summarize fixed-input ARM field-multiply process pairs."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re
import statistics


LINE = re.compile(
    r"n=(\d+) fold=(true|false) repetition=(\d+) ns_per_mul=([0-9.]+) checksum=(\d+)")
NAMES = ("table1b", "fold1b", "fold2b", "table2b")


def main(directory: Path, output: Path) -> None:
    rows = {}
    for name in NAMES:
        found = [LINE.fullmatch(line).groups()
                 for line in (directory / f"{name}.txt").read_text().splitlines()]
        assert len(found) == 18
        expected_fold = name.startswith("fold")
        values = {}
        for n in (41, 53):
            selected = [entry for entry in found if int(entry[0]) == n]
            assert len(selected) == 9
            assert {int(entry[2]) for entry in selected} == set(range(9))
            assert all((entry[1] == "true") == expected_fold for entry in selected)
            assert len({entry[4] for entry in selected}) == 1
            values[str(n)] = {
                "median_ns_per_mul_after_warmup": statistics.median(
                    float(entry[3]) for entry in selected if int(entry[2]) > 0),
                "checksum": selected[0][4]}
        rows[name] = values
    curves = {}
    for n in (41, 53):
        key = str(n)
        assert len({rows[name][key]["checksum"] for name in NAMES}) == 1
        table = [rows[name][key]["median_ns_per_mul_after_warmup"]
                 for name in ("table1b", "table2b")]
        fold = [rows[name][key]["median_ns_per_mul_after_warmup"]
                for name in ("fold1b", "fold2b")]
        curves[key] = {
            "table_process_median_ns_per_mul": table,
            "fold_process_median_ns_per_mul": fold,
            "exploratory_table_over_fold_geometric_ratio": (
                math.sqrt(table[0] * table[1] / (fold[0] * fold[1]))),
            "table_aa_abs_log_ratio": abs(math.log(table[0] / table[1])),
            "fold_aa_abs_log_ratio": abs(math.log(fold[0] / fold[1])),
            "same_checksum": True}
    result = {"kind": "arm_pmull_fold_fixed_operand_kernel_diagnostic_v1",
              "architecture": "arm64", "operands_per_repetition": 1_048_576,
              "repetitions_per_process": 9,
              "process_order": list(NAMES),
              "controlled_speedup": None,
              "curves": curves}
    output.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    print(json.dumps(curves, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    main(args.directory, args.output)

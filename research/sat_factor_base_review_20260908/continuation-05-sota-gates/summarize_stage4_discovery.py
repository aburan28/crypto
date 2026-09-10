#!/usr/bin/env python3
"""Verify the charged public GGMP factor/curve census and freeze its winner."""

from __future__ import annotations

import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RUN = ROOT / "stage-4-ggmp-discovery-20260909"


def main() -> None:
    rows = []
    for curve_a in [0, 1]:
        for factor_index in range(6):
            name = f"a{curve_a}-f{factor_index}"
            metrics = json.loads((RUN / f"{name}.metrics.json").read_text())
            manifest_path = RUN / f"candidate-{name}" / "manifest.json"
            if manifest_path.exists():
                manifest = json.loads(manifest_path.read_text())
                count = manifest["factor_base_geometry"]["distinct_curve_points"]
                status = "admitted"
                factor = manifest["factor_base_predicate"]["factor_bitmask"]
                exponents = manifest["factor_base_predicate"]["linearised_exponents"]
            else:
                stderr = (RUN / f"{name}.stderr").read_text()
                match = re.search(r"degenerate factor base: only (\d+) curve point", stderr)
                assert match, f"missing geometry failure for {name}"
                count = int(match.group(1))
                status = "rejected_degenerate"
                factor = None
                exponents = None
            rows.append(
                {
                    "curve_a": curve_a,
                    "factor_index": factor_index,
                    "status": status,
                    "distinct_curve_points": count,
                    "factor_bitmask": factor,
                    "linearised_exponents": exponents,
                    "returncode": metrics["returncode"],
                    "wall_seconds": metrics["metrics"]["wall_seconds"],
                    "total_core_seconds": metrics["metrics"]["total_core_seconds"],
                    "peak_rss_bytes": metrics["metrics"]["peak_rss_bytes"],
                }
            )
    winner = min(rows, key=lambda row: (-row["distinct_curve_points"], row["curve_a"], row["factor_index"]))
    assert winner["curve_a"] == 0 and winner["factor_index"] == 0
    assert winner["distinct_curve_points"] == 63
    assert sum(row["status"] == "admitted" for row in rows) == 6
    assert sum(row["status"] == "rejected_degenerate" for row in rows) == 6
    output = {
        "schema": "koblitz_ggmp_public_factor_discovery.v1",
        "search_space_complete": True,
        "uses_target_scalar": False,
        "uses_discrete_log_labels": False,
        "enumerates_target_subgroup": False,
        "rows": rows,
        "selection_rule": "maximum distinct_curve_points, then minimum curve_a, then minimum factor_index",
        "selected": winner,
        "total_discovery_core_seconds": sum(row["total_core_seconds"] for row in rows),
        "peak_discovery_rss_bytes": max(row["peak_rss_bytes"] for row in rows),
    }
    (ROOT / "stage-4-discovery-summary.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Fail-closed verifier for GitHub-hosted Koblitz toy reproduction artifacts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    rows = []
    for secret in [53, 101]:
        ic = json.loads((args.directory / f"ic-{secret}.json").read_text())
        ic_meter = json.loads((args.directory / f"ic-{secret}.metrics.json").read_text())
        rho = json.loads((args.directory / f"rho-auto-{secret}.json").read_text())
        rho_meter = json.loads((args.directory / f"rho-auto-{secret}.metrics.json").read_text())
        assert ic_meter["returncode"] == 0 and not ic_meter["timed_out"]
        assert rho_meter["returncode"] == 0 and not rho_meter["timed_out"]
        assert ic["target"] == rho["target"]
        assert ic["report"]["verified_unknown_scalar_recovery"] is True
        assert ic["report"]["recovered_scalar"] == str(secret)
        assert ic["report"]["direct_relation"] is False
        assert ic["report"]["sat_invalid_models"] == 0
        assert ic["factor_base_predicate"]["factor_base_logs_constructed"] is False
        assert ic["factor_base_predicate"]["uses_discrete_log_labels"] is False
        assert ic["factor_base_predicate"]["enumerates_target_subgroup"] is False
        assert rho["automorphism_optimized"] is True
        assert rho["verified_unknown_scalar_recovery"] is True
        rows.append(
            {
                "secret": secret,
                "index_calculus": {
                    "relations": ic["report"]["relations"],
                    "trials": ic["report"]["trials"],
                    "conflicts": ic["report"]["sat_conflicts"],
                    "wall_seconds": ic_meter["metrics"]["wall_seconds"],
                    "core_seconds": ic_meter["metrics"]["total_core_seconds"],
                    "peak_rss_bytes": ic_meter["metrics"]["peak_rss_bytes"],
                },
                "automorphism_rho": {
                    "iterations": rho["iterations"],
                    "reported_group_additions": rho["reported_group_additions"],
                    "wall_seconds": rho_meter["metrics"]["wall_seconds"],
                    "core_seconds": rho_meter["metrics"]["total_core_seconds"],
                    "peak_rss_bytes": rho_meter["metrics"]["peak_rss_bytes"],
                },
            }
        )
    result = {
        "schema": "koblitz_external_host_reproduction.v1",
        "environment": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": platform.python_version(),
        },
        "rows": rows,
        "all_checks_pass": True,
        "scope": "independent GitHub-hosted environment reproduction of toy correctness and accounting; not an external novelty review",
    }
    (args.directory / "external-reproduction-summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

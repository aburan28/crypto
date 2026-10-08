#!/usr/bin/env python3
"""Create a compact, source-closed Stage-23 terminal-evidence bundle."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
from koblitz_stage23_terminal_evidence import EvidenceError, package_bundle


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--outer-root", type=Path, required=True)
    parser.add_argument("--outer-metrics", type=Path, required=True)
    parser.add_argument("--project-verification-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = package_bundle(
            run_root=args.run_root,
            outer_root=args.outer_root,
            outer_metrics=args.outer_metrics,
            project_verification_root=args.project_verification_root,
            output=args.output,
        )
    except (EvidenceError, OSError, ValueError, KeyError) as error:
        parser.exit(1, f"stage23-terminal-package: {error}\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

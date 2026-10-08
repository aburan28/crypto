#!/usr/bin/env python3
"""Verify a compact Stage-23 terminal-evidence bundle from archive-relative bytes."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
from koblitz_stage23_terminal_evidence import EvidenceError, verify_bundle


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    args = parser.parse_args()
    try:
        result = verify_bundle(args.bundle)
    except (EvidenceError, OSError, ValueError, KeyError) as error:
        parser.exit(1, f"stage23-terminal-verify: {error}\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

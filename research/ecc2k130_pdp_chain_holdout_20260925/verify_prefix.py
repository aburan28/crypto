#!/usr/bin/env python3
"""Check that the predeclared 2,048-attempt extension kept its 512 prefix."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("initial", type=Path)
    parser.add_argument("extended", type=Path)
    args = parser.parse_args()
    a, b = (json.loads(p.read_text()) for p in (args.initial, args.extended))
    assert a["parameters"]["attempts"] == 512
    assert b["parameters"]["attempts"] == 2048
    for key in ("target_coefficients", "masks", "source_targets",
                "transported_targets", "target_keys"):
        assert a[key] == b[key][:512], key
    assert a["geometry"]["bases"] == b["geometry"]["bases"]
    for name in a["arms"]:
        assert a["arms"][name]["cases"] == b["arms"][name]["cases"][:512], name
    print("PASS: 512 masked targets, masks, bases and eight arm outcomes "
          "match the frozen 2,048-attempt prefix")


if __name__ == "__main__":
    main()

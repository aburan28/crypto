#!/usr/bin/env python3
"""Ledger §23's one-target parameter files (PROTOCOL.md, "Recipe", "Targets").

Each file is §20's recipe for its size, built by §20's own rules
(`research/ic_exponent_20260926/make_params.py`), with one public
hash-to-curve target in place of §20's 32 known-answer targets:

- the four sizes §20 swept keep §20's chosen (columns, descent summands),
  and `make_params.params` rebuilds §20's frozen M1 recipe exactly;
- the two sizes §20 never ran take §20's model optimum, read in actual
  columns of 8*ceil(c/8), as §20 read its sweep grid;
- recipe seed 201 (the base and the collection, as §20's M1);
- target T<i>: `public_hash_seed` 23000 + i, whose logarithm nobody
  constructs (`ic workflow`'s public-target domain).

    python3 make_params.py <a> <n> <i> <out.json>
"""
from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
S20 = HERE.parents[1] / "research" / "ic_exponent_20260926" / "make_params.py"
# §20's rules, unchanged, under their own name (this file shares theirs).
_spec = importlib.util.spec_from_file_location("s20_make_params", S20)
s20 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s20)

s20.R.update({(0, 57): 275295876199, (1, 59): 25179555920633})

RECIPE_SEED = 201
TARGET_SEED_BASE = 23000
# (a, n): (columns, descent summands, source), in the order of r.
RECIPES = {
    (1, 47): (56, 3, "§20 sweep"),
    (0, 57): (8 * math.ceil(61 / 8), 2, "§20 model optimum (61 columns)"),
    (0, 41): (80, 2, "§20 sweep"),
    (0, 53): (264, 2, "§20 sweep"),
    (1, 59): (8 * math.ceil(251 / 8), 2, "§20 model optimum (251 columns)"),
    (0, 61): (320, 2, "§20 sweep"),
}
SIZES = list(RECIPES)


def target_seed(i: int) -> int:
    return TARGET_SEED_BASE + i


def params(a: int, n: int, i: int) -> dict:
    columns, m, _ = RECIPES[(a, n)]
    p = s20.params(a, n, columns, m, RECIPE_SEED)
    p["name"] = f"k{a}n{n}-c{columns}-m{m}-s{RECIPE_SEED}-T{i:02d}"
    p["targets"] = [{"public_hash_seed": target_seed(i)}]
    return p


def main() -> None:
    a, n, i = (int(x) for x in sys.argv[1:4])
    with open(sys.argv[4], "x") as f:
        json.dump(params(a, n, i), f, indent=1)
        f.write("\n")


if __name__ == "__main__":
    main()

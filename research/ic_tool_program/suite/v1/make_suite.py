#!/usr/bin/env python3
"""Suite v1 of the `ic` tool programme (IC_TOOL_PROGRAM.md §4).

Writes one `ic price --single-target` parameter file per row under
`params/`, and `SUITE.json` with every file's SHA-256 and each row's rho
seed.  Every recipe comes from §20's own rules
(`research/ic_exponent_20260926/make_params.py`), loaded through §23's
module so that the two rounds cannot drift apart.

- **The S suite.** The eleven Koblitz curves, in order of r. Each has
  §20's measurement seeds 201-204 (`M1`-`M4`) and two public
  hash-to-curve targets per seed: `T01`-`T08` of §23's law
  (`public_hash_seed` 23000 + i), two to a seed in order. The recipe is
  §20's frozen choice at its nine swept sizes and §20's model optimum
  at the two it never ran, as in §23.  `M1`'s rows at §23's six sizes
  are byte-for-byte §23's `T01`/`T02` files.
- **The smoke tier.** `E_0/GF(2^31)` by §20's rules at 8 columns,
  m = 2, seed 201, targets `T01`-`T02`. It is not evidence.

Each row's rho seed is 0x230000 + i, as in §23.

    python3 make_suite.py           # writes params/ and SUITE.json; refuses to overwrite
    python3 make_suite.py --check   # re-derives every file and compares bytes and hashes
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
S23_DIR = ROOT / "research" / "ic_single_target_20260930"
_spec = importlib.util.spec_from_file_location("s23_make_params", S23_DIR / "make_params.py")
s23 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s23)
s20 = s23.s20

# E_0/GF(2^31): #E = 4 * 373 * 1439393, the last factor prime.
s20.R.setdefault((0, 31), 1439393)

# (a, n): (columns, descent summands, source), in the order of r.
RECIPES = {
    (1, 19): (8, 2, "§20 sweep"),
    (1, 23): (8, 2, "§20 sweep"),
    (1, 45): (8, 2, "§20 sweep"),
    (0, 37): (8, 2, "§20 sweep"),
    (1, 43): (16, 2, "§20 sweep"),
    (1, 47): (56, 3, "§20 sweep"),
    (0, 57): (8 * math.ceil(61 / 8), 2, "§20 model optimum (61 columns)"),
    (0, 41): (80, 2, "§20 sweep"),
    (0, 53): (264, 2, "§20 sweep"),
    (1, 59): (8 * math.ceil(251 / 8), 2, "§20 model optimum (251 columns)"),
    (0, 61): (320, 2, "§20 sweep"),
}
SMOKE = {(0, 31): (8, 2, "§20's rules at 8 columns; smoke only")}
SEEDS = (201, 202, 203, 204)
TARGETS_PER_SEED = 2
TARGET_SEED_BASE = 23000
RHO_SEED_BASE = 0x230000


def params(a: int, n: int, columns: int, m: int, seed: int, i: int) -> dict:
    """§23's construction, with the recipe seed free."""
    p = s20.params(a, n, columns, m, seed)
    p["name"] = f"k{a}n{n}-c{columns}-m{m}-s{seed}-T{i:02d}"
    p["targets"] = [{"public_hash_seed": TARGET_SEED_BASE + i}]
    return p


def text(p: dict) -> str:
    return json.dumps(p, indent=1) + "\n"


def rows() -> list[dict]:
    out = []
    for tier, recipes, seeds in (("S", RECIPES, SEEDS), ("smoke", SMOKE, SEEDS[:1])):
        for (a, n), (columns, m, source) in recipes.items():
            r = s20.R[(a, n)]
            for j, seed in enumerate(seeds):
                for t in range(TARGETS_PER_SEED):
                    i = TARGETS_PER_SEED * j + t + 1
                    rel = f"params/{tier}/k{a}n{n}/M{seed - 200}-T{i:02d}.json"
                    out.append({
                        "id": f"k{a}n{n}-M{seed - 200}-T{i:02d}", "tier": tier, "a": a, "n": n, "r": r,
                        "log2_r": round(math.log2(r), 3), "columns": columns, "descent_summands": m,
                        "recipe_source": source, "recipe_seed": seed, "target": i,
                        "public_hash_seed": TARGET_SEED_BASE + i, "rho_seed": RHO_SEED_BASE + i,
                        "params": rel, "text": text(params(a, n, columns, m, seed, i)),
                    })
    return out


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def manifest(rs: list[dict]) -> dict:
    return {
        "suite": "ic tool programme suite v1",
        "plan": "research/notes/index-calculus/IC_TOOL_PROGRAM.md §4",
        "generator": {"path": "research/ic_tool_program/suite/v1/make_suite.py",
                      "sha256": sha256(Path(__file__).read_bytes())},
        "recipe_rules": {"path": "research/ic_exponent_20260926/make_params.py",
                         "sha256": sha256((ROOT / "research/ic_exponent_20260926/make_params.py").read_bytes())},
        "command": "ic price --params <params> --json --out <out> --single-target --rho-seed <rho_seed>",
        "rows": [{k: v for k, v in r.items() if k != "text"} | {"sha256": sha256(r["text"].encode())}
                 for r in rs],
    }


def main() -> None:
    rs = rows()
    doc = manifest(rs)
    if "--check" in sys.argv:
        frozen = json.loads((HERE / "SUITE.json").read_text())
        bad = [r["id"] for r in rs if (HERE / r["params"]).read_text() != r["text"]]
        bad += ["SUITE.json"] if frozen["rows"] != doc["rows"] else []
        for r in frozen["rows"]:
            if sha256((HERE / r["params"]).read_bytes()) != r["sha256"]:
                bad.append(r["params"])
        print(json.dumps({"rows": len(rs), "mismatches": bad}, indent=1))
        raise SystemExit(1 if bad else 0)
    if (HERE / "SUITE.json").exists():
        raise SystemExit("SUITE.json exists; suite v1 is frozen (use --check)")
    for r in rs:
        path = HERE / r["params"]
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "x") as f:
            f.write(r["text"])
    (HERE / "SUITE.json").write_text(json.dumps(doc, indent=1) + "\n")
    print(f"{len(rs)} rows written")


if __name__ == "__main__":
    main()

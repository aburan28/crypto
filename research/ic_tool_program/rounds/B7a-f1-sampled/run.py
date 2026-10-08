#!/usr/bin/env python3
"""B7a's measurements (PROTOCOL.md, with amendment 1).

    IC_BASE=<base> IC_BASE_COMMIT=<sha> IC_CAND=<B7a build> IC_CAND_COMMIT=<sha> \\
    IC_RUNS=<run tree> python3 run.py --steps <accepted>,B7a <step> [<step> ...]

Measurements 2-4 are Track B's, run by `harness/bround.py`: `conformance`,
`pin`, `translate` and `timing` (the base against B7a, F0, on M1's 22
rows, five rounds ABAB). B7a adds:

- `f1` (measurement 5): F1 on M1's 22 rows, three rounds, isolated. Each
  row's document is its v2 translation (`ic check --translate`, as
  `translate` writes it) with `method.fidelity` set to `F1`.
- `partial` (measurement 7): F1 forced to partial tables, `--f1-partial`
  at 1/2 and 1/4, on M1's rows at the six largest sizes, one round.
- `analyse`: everything, as JSON on stdout (`analyse.py`).

Every output stays on disk and every step resumes where it stopped. A
process that is contended or fails is kept and run again, at most twice,
as the programme's runner does.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "harness"))
import bench  # noqa: E402
import bround  # noqa: E402

F1_ROUNDS = 3
PARTIAL_FRACTIONS = (0.5, 0.25)
LARGEST = 6


def f1_document(row: dict) -> Path:
    """The row's v2 translation at `fidelity: F1`, written once."""
    path = bround.runs() / "f1" / "docs" / f"{row['id']}.json"
    if not path.exists():
        doc = bround.translate_row(row)
        doc["method"]["fidelity"] = "F1"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(doc, indent=1) + "\n")
    return path


def attempt_path(out: Path, k: int) -> Path:
    return out if k == 0 else out.with_name(bench.stem(out).name + f"-retry{k}.price.json")


def figure(out: Path) -> Path:
    """The first clean, extrapolated attempt: the process's figure."""
    for k in range(bench.RETRIES + 1):
        attempt = attempt_path(out, k)
        if attempt.exists() and bench.clean(attempt) and bench.load(attempt).get("status") == "extrapolated":
            return attempt
    return out


def f1_price(doc: Path, out: Path, log_dir: Path, extra: tuple[str, ...] = ()) -> dict:
    """The first clean run of up to 1 + RETRIES; every attempt is kept."""
    rep: dict = {}
    for k in range(bench.RETRIES + 1):
        attempt = attempt_path(out, k)
        if not attempt.exists():
            attempt.parent.mkdir(parents=True, exist_ok=True)
            cmd = [str(bround.arm("IC_CAND")), "price", "--params", str(doc), "--json", "--out",
                   str(attempt), *extra]
            bench.launch(cmd, attempt, log_dir)
        rep = bench.load(attempt)
        if bench.clean(attempt) and rep.get("status") == "extrapolated":
            return rep
    return rep


def note(rep: dict) -> str:
    cold = rep.get("extrapolated", {}).get("ic", {}).get("cold", {}).get("units", {})
    return f"median {cold.get('median', float('nan')):.4g} units" if cold else ""


def f1(_steps: str) -> dict:
    d = bround.runs() / "f1"
    for k in range(1, F1_ROUNDS + 1):
        for row in bround.m1():
            out = d / f"r{k}" / f"{row['id']}.price.json"
            rep = f1_price(f1_document(row), out, d)
            print(f"r{k} {row['id']} F1: {rep.get('status')} {note(rep)}", flush=True)
    return {"done": True}


def largest_rows() -> list[dict]:
    sizes = sorted({(r["a"], r["n"], r["r"]) for r in bround.m1()}, key=lambda s: s[2])[-LARGEST:]
    keep = {(a, n) for a, n, _ in sizes}
    return [r for r in bround.m1() if (r["a"], r["n"]) in keep]


def partial(_steps: str) -> dict:
    d = bround.runs() / "partial"
    for g in PARTIAL_FRACTIONS:
        for row in largest_rows():
            out = d / f"g{g}" / f"{row['id']}.price.json"
            rep = f1_price(f1_document(row), out, d, ("--f1-partial", str(g)))
            print(f"g={g} {row['id']}: {rep.get('status')} {note(rep)}", flush=True)
    return {"done": True}


def analyse(steps: str) -> dict:
    sys.path.insert(0, str(HERE))
    import analyse as b7a  # noqa: E402

    return b7a.analyse(steps)


STEPS = {"conformance": bround.conformance, "pin": bround.pin, "translate": bround.translate,
         "timing": bround.timing, "f1": f1, "partial": partial, "analyse": analyse,
         "manifest": bround.manifest}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--steps", required=True, help="the conformance steps, e.g. B0,B1,B2,B2b,B3,B7a")
    p.add_argument("step", nargs="+", choices=list(STEPS))
    args = p.parse_args()
    bround.runs().mkdir(parents=True, exist_ok=True)
    for step in args.step:
        result = STEPS[step](args.steps)
        print(json.dumps({step: result}, indent=1), flush=True)


if __name__ == "__main__":
    main()

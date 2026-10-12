#!/usr/bin/env python3
"""B2b's measurement 5, the subfield sweep: `kic` paired with `rho-negation`
on one curve over GF(2^k) for every k in 2..8 and odd e >= 3 with
16 <= n = k·e <= 40, two known-answer targets each (PROTOCOL.md).

    python3 sweep.py --ic <binary> --out <dir>

The documents are built here, in the conformance generator's own
arithmetic (`../../conformance/v2-b2b/make_cases.py`), which shares
nothing with the Rust tool.  Every answer the tool reports is replayed
here as [d]G = Q.  Untimed: run it under the benchmark lock
(`tools/isolated_bench.py busy --wait -- timeout ...`).

The curve for (k, e): y^2 + xy = x^3 + ax^2 + b under the repository's
modulus for n, with (a, b) the first over GF(2^k) — a in 0, 1, then the
rest ascending; b ascending — whose least subfield is GF(2^k) and whose
#E has a prime factor r > h with r^2 not dividing #E.  A pair with no
such curve is listed as such.

The recipe: ledger §20's rules with e in the place of n
(`make_cases.recipe`), at max(8, ceil(r^(1/3) / 2e)) columns, which
gives the suite's own |F| ≈ r^(1/3) at its sizes, and m = 2.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import subprocess
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location(
    "conformance_v2_b2b_cases", HERE.parents[1] / "conformance" / "v2-b2b" / "make_cases.py")
mc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mc)
v2, b2 = mc.v2, mc.b2

LABEL = "ic-b2b-sweep"
TARGETS = 2
TIMEOUT_S = 900
PANIC_STATUS = 101


def pairs() -> list[tuple[int, int]]:
    return [(k, e) for k in range(2, 9) for e in range(3, 41, 2) if 16 <= k * e <= 40]


def largest_prime_factor(m: int, label: str) -> int:
    return max(q for q, _ in mc.factorise(m, label))


def sweep_curve(k: int, e: int):
    n = k * e
    f = v2.curve_id.find_irreducible_sparse(n)
    C0 = v2.Curve(n, f, 0, 1)
    elements = mc.subfield(C0.k, k)
    a_order = [0, 1] + [x for x in elements if x > 1]
    for a in a_order:
        for b in elements:
            if b == 0 or mc.least_subfield(C0.k, a, b) != k:
                continue
            C = v2.Curve(n, f, a, b)
            N = mc.order_over_extension((1 << k) + 1 - mc.count_over_subfield(C, k, elements), 1 << k, e)
            r = largest_prime_factor(N, f"{LABEL}/{k}/{e}/{a}/{b}")
            if r > N // r and (N // r) % r:
                return mc.subfield_curve(n, k, a, b, f"{LABEL}/{k}/{e}")
    return None


def documents() -> tuple[list[dict], list[str]]:
    """The runs, and the pairs (k, e) on which no curve has r > h.  For a composite e, #E(GF(q^e)) is
    a multiple of #E(GF(q^d)) for every d dividing e, which can leave no such curve at all."""
    out, none = [], []
    for k, e in pairs():
        found = sweep_curve(k, e)
        if found is None:
            none.append(f"k{k}e{e}")
            continue
        C, N, r, h = found
        G = C.subgroup_point(f"{LABEL}/{k}/{e}/generator", h, r)
        columns = max(8, math.ceil(round(r ** (1 / 3), 6) / (2 * e)))
        for t in range(TARGETS):
            d = b2.known_log(f"{LABEL}/{k}/{e}/known_log/{t}", r)
            doc = {
                "schema_version": 2,
                "name": f"B2b sweep: k = {k}, e = {e}, n = {k * e}, target {t}",
                "field": {"kind": "binary", "degree": C.n, "modulus": v2.hx(C.f)},
                "curve": {"form": "binary_weierstrass", "a": v2.hx(C.a), "b": v2.hx(C.b)},
                "subgroup": {"order": str(r), "cofactor": str(h), "generator": v2.point_doc(G)},
                "target": {"known_log": str(d)},
                "method": {"solve": "paired",
                           "index_calculus": {"pipeline": "kic", "recipe": mc.recipe(e, r, columns)},
                           "rho": {"pipeline": "auto", "seed": mc.RHO_SEED}},
            }
            out.append({"k": k, "e": e, "n": k * e, "r": r, "h": h, "columns": columns, "target": t,
                        "known_log": d, "curve": C, "generator": G, "doc": doc})
    return out, none


def replay(row: dict, scalar: str | None) -> bool:
    if scalar is None:
        return False
    d = int(scalar)
    C, G = row["curve"], row["generator"]
    return 0 < d < row["r"] and C.mul(d, G) == C.mul(row["known_log"], G)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ic", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    (args.out / "params").mkdir()
    rows = []
    runs, no_curve = documents()
    for row in runs:
        stem = f"k{row['k']}e{row['e']}-t{row['target']}"
        params = args.out / "params" / f"{stem}.json"
        params.write_text(json.dumps(row["doc"], indent=1) + "\n")
        report_path = args.out / f"{stem}.report.json"
        started = time.monotonic()
        try:
            proc = subprocess.run([str(args.ic), "price", "--params", str(params), "--json", "--out",
                                   str(report_path), "--repeats", "1", "--repeats-fast", "1"],
                                  capture_output=True, text=True, timeout=TIMEOUT_S,
                                  env={"RAYON_NUM_THREADS": "1", "PATH": "/usr/bin:/bin"})
            status, stderr = proc.returncode, proc.stderr
        except subprocess.TimeoutExpired:
            status, stderr = "timeout", ""
        wall = time.monotonic() - started
        report = None
        if report_path.exists():
            try:
                report = json.loads(report_path.read_text())
            except json.JSONDecodeError:
                report = None
        get = (lambda *ks: _dig(report, ks)) if report else (lambda *ks: None)
        findings = []
        if status == PANIC_STATUS or "panicked" in stderr:
            findings.append("panic")
        if status not in (0, 1):
            findings.append(f"exit {status}")
        if report is None and status in (0, 1):
            findings.append("no report")
        complete = status == 0 and get("status") == "complete"
        if complete:
            if get("result", "verified") is not True:
                findings.append("complete but not verified")
            if get("result", "scalar") != str(row["known_log"]) or not replay(row, get("result", "scalar")):
                findings.append("wrong answer")
            if (get("route", "ic", "pipeline"), get("route", "rho", "pipeline")) != ("kic", "rho-negation"):
                findings.append("route")
        rows.append({k: row[k] for k in ("k", "e", "n", "r", "h", "columns", "target", "known_log")}
                    | {"params_sha256": hashlib.sha256(params.read_bytes()).hexdigest(), "exit": status,
                       "status": get("status"), "scalar": get("result", "scalar"),
                       "verified": get("result", "verified"), "wall_s": round(wall, 3),
                       "ic_online_ms": get("median", "ic_online_wall_ms"),
                       "rho_online_ms": get("median", "rho_online_wall_ms"), "findings": findings,
                       "stderr_tail": stderr[-400:] if findings else ""})
        print(f"k={row['k']} e={row['e']} n={row['n']} t={row['target']}: exit {status}, "
              f"{get('status')}, {wall:.1f} s{', FINDINGS ' + str(findings) if findings else ''}", flush=True)
    summary = {"runs": len(rows), "no_curve": no_curve,
               "complete": sum(r["status"] == "complete" for r in rows),
               "not_complete": [f"k{r['k']}e{r['e']}-t{r['target']}" for r in rows if r["status"] != "complete"],
               "findings": sum(bool(r["findings"]) for r in rows),
               "sweep_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               "generator_sha256": hashlib.sha256(Path(_spec.origin).read_bytes()).hexdigest()}
    (args.out / "results.json").write_text(json.dumps({"summary": summary, "rows": rows}, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    raise SystemExit(1 if summary["findings"] else 0)


def _dig(doc, keys):
    for k in keys:
        if not isinstance(doc, dict) or k not in doc:
            return None
        doc = doc[k]
    return doc


if __name__ == "__main__":
    main()

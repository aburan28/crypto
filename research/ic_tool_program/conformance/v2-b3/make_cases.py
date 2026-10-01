#!/usr/bin/env python3
"""Conformance cases for B3: rho on two-word binary fields, and rho alone
(`research/ic_tool_program/rounds/B3-two-word-rho/PROTOCOL.md`).

    python3 make_cases.py           # writes params/, cases.json and SHA256SUMS; refuses to overwrite
    python3 make_cases.py --check   # re-derives every file and compares the bytes

B1's cases and documents (`../v2/`) are frozen, so B3 adds its cases
here, beside them:
- most are B1's own documents with their method or target edited by the
  runner (`copy` with `set`), so their curves, generators and targets are
  the ones B1 checked;
- one new document, the translation of a suite row at `n = 61`, is
  built here with B1's generator code (`../v2/make_cases.py`), whose
  arithmetic shares nothing with the Rust tool.

`../run.py` runs every step's cases, these among them, with the `until`
rule applied.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
V2 = HERE.parent / "v2"
_spec = importlib.util.spec_from_file_location("conformance_v2_cases", V2 / "make_cases.py")
v2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(v2)

STEP = "B3"
RHO_SEED = 0x230000 + 1  # the suite's rho seed for target T01, and the gate document's
GATE = "{cases}/gate-m83-T001.json"
CHALLENGE = "{cases}/ecc2k130-challenge.json"
SMOKE = "{cases}/C009-translation-of-smoke-row.json"
CURVE_A = "{cases}/C010-curve-a-other-generator.json"
ROW61 = "params/S/k0n61/M1-T01.json"
COMPARED = ["rho_counts", "certificates.rho.scalar", "certificates.rho.target", "certificates.rho.curve"]


def rho_method(seed: int, cap: int | None = None) -> dict:
    rho = {"pipeline": "auto", "seed": seed}
    if cap is not None:
        rho["max_iterations"] = cap
    return {"solve": "rho", "fidelity": "F0", "rho": rho}


def case(cid: str, purpose: str, files: dict, argv: list[str], expect: dict, timeout: int = 120,
         until: str | None = None) -> dict:
    out = {"id": cid, "step": STEP, "checks": purpose, "files": files, "argv": argv,
           "env": {"RAYON_NUM_THREADS": "1"}, "expect": expect, "timeout_s": timeout}
    if until:
        out["until"] = until
    return out


def price(params: str = "{tmp}/params.json", report: str = "{tmp}/report.json") -> list[str]:
    return ["price", "--params", params, "--json", "--out", report]


def build() -> tuple[dict[str, dict], list[dict]]:
    files: dict[str, dict] = {}
    cases: list[dict] = []
    c010 = json.loads((V2 / "params" / "C010-curve-a-other-generator.json").read_text())
    k10 = c010["target"]["known_log"]
    rho_admitted = {"pipeline": "rho-koblitz", "admitted": True}
    kic_one_word = {"pipeline": "kic", "admitted": False, "gate": "field-wider-than-one-word"}

    cases.append(case(
        "C052-rho-alone-curve-a", "B3's solve: rho on curve A, with C010's generator and known answer",
        {"params.json": {"copy": CURVE_A, "set": {"method": rho_method(7)}}}, price(),
        {"exit": 0, "json_file": "{tmp}/report.json",
         "json_paths": {"status": "complete", "operation": "rho_single_target", "fidelity": "F0",
                        "result.scalar": k10, "result.verified": True, "result.known_answer": True,
                        "rho.pipeline": "rho-koblitz"},
         "json_contains": {"route.considered": [rho_admitted]}}))

    cases.append(case(
        "C053-gate-rho-capped", "B3's routing at n = 83: rho alone runs on the gate curve; a run stopped by "
        "its step cap reports no logarithm and exits 1",
        {"params.json": {"copy": GATE, "set": {"method": rho_method(RHO_SEED, 200000)}}}, price(),
        {"exit": 1, "json_file": "{tmp}/report.json",
         "json_paths": {"status": "not_recovered", "result.verified": False, "rho.step_cap": 200000,
                        "rho.counters.walk_operations": 200000, "rho.pipeline": "rho-koblitz"},
         "json_contains": {"route.considered": [kic_one_word, rho_admitted]}}))

    cases.append(case(
        "C054-gate-paired-has-no-ic-route", "C027's successor from B3 (the until rule): rho-koblitz admits the "
        "gate curve, kic does not, so paired is refused for want of an index-calculus route",
        {"params.json": {"copy": GATE}}, price() + ["--repeats", "1", "--repeats-fast", "1"],
        {"exit": 3, "json_file": "{tmp}/report.json",
         "json_paths": {"refusal.code": "no-ic-route", "refusal.class": "unsupported"},
         "json_contains": {"route.considered": [kic_one_word, rho_admitted]}},
        until="B3b"))

    cases.append(case(
        "C055-challenge-past-two-words", "B3's rho gate at n = 131: field-wider-than-two-words",
        {"params.json": {"copy": CHALLENGE}},
        ["check", "--params", "{tmp}/params.json", "--json", "--out", "{tmp}/report.json"],
        {"exit": 0, "json_file": "{tmp}/report.json", "json_paths": {"status": "checks_passed"},
         "json_contains": {"route.considered": [
             {"pipeline": "rho-koblitz", "admitted": False, "gate": "field-wider-than-two-words"}]}},
        until="B4"))

    cases.append(case(
        "C056-gate-rho-hashed-target", "rho alone past one word takes a point or a known logarithm; v1's "
        "hashed target there waits for B4",
        {"params.json": {"copy": GATE, "set": {"target": {"public_hash_seed": 1},
                                               "method": rho_method(RHO_SEED, 1000)}}}, price(),
        {"exit": 3, "json_file": "{tmp}/report.json",
         "json_paths": {"refusal.code": "not-yet-supported", "refusal.class": "unsupported"}},
        until="B4"))

    smoke_alone = case(
        "C057-smoke-rho-alone-is-the-paired-rho", "rho alone (the two-word walk) and the paired price's rho arm "
        "(ParallelRho) on the smoke row's translation: the same counts, operations, logarithm and certificate",
        {"params.json": {"copy": SMOKE, "set": {"method": rho_method(RHO_SEED)}},
         "paired.json": {"copy": SMOKE}}, price(),
        {"exit": 0, "json_file": "{tmp}/report.json",
         "json_paths": {"status": "complete", "result.verified": True},
         "same_outputs_as": {"argv": price("{tmp}/paired.json", "{tmp}/paired-report.json")
                             + ["--repeats", "1", "--repeats-fast", "1"],
                             "json_file": "{tmp}/paired-report.json", "paths": COMPARED}})
    cases.append(smoke_alone)

    # C058: the same at n = 61, against the v1 row itself.
    row = v2.v1_row(ROW61)
    n = row["curve"]["degree"]
    assert (n, row["curve"]["curve_a"]) == (61, 0)
    f = v2.curve_id.find_irreducible_sparse(n)
    assert v2.irreducible(f)
    curve = v2.Curve(n, f, 0, 1)
    order = v2.koblitz_order(0, n)
    r = v2.largest_prime_factor(order)
    h = order // r
    assert v2.prime_exact(r) and r > h and h % r
    assert v2.curve_id.binary_id(n, f, 0, 1, order, end="-7")["slug"] == "icv1-f2m61-t158598901-ab42b6c5"
    seed = row["targets"][0]["public_hash_seed"]
    doc = v2.document(f"C058: suite v1 row {ROW61}, translated to v2, rho alone", curve,
                      {"form": "koblitz", "a": 0}, r, h, {"rule": "koblitz_search_v1"},
                      {"public_hash_seed": seed}, v2.recipe_of(row), RHO_SEED)
    doc["method"] = rho_method(RHO_SEED)
    files["C058-n61-rho-alone.json"] = doc
    cases.append(case(
        "C058-n61-rho-alone-is-the-paired-rho", "C057 at n = 61, against the v1 row's single-target price",
        {"params.json": {"copy": "{here}/C058-n61-rho-alone.json"}, "v1.json": {"copy": "{suite}/" + ROW61}},
        price(),
        {"exit": 0, "json_file": "{tmp}/report.json",
         "json_paths": {"status": "complete", "result.verified": True},
         "same_outputs_as": {"argv": price("{tmp}/v1.json", "{tmp}/v1-report.json")
                             + ["--single-target", "--rho-seed", str(RHO_SEED), "--repeats", "1",
                                "--repeats-fast", "1"],
                             "json_file": "{tmp}/v1-report.json", "paths": COMPARED}},
        timeout=900))
    return files, cases


def texts() -> dict[str, str]:
    files, cases = build()
    out = {f"params/{name}": json.dumps(doc, indent=1) + "\n" for name, doc in files.items()}
    out["cases.json"] = json.dumps({
        "suite": "ic tool programme conformance suite v2, B3's cases",
        "design": "research/ic_tool_program/design/schema-v2.md §9-§10; "
                  "research/ic_tool_program/rounds/B3-two-word-rho/PROTOCOL.md",
        "includes": "../v1/cases.json (step B0) and ../v2/cases.json (B1), run first, with the until rule",
        "rules": [
            "The rules of ../v2/cases.json hold. {here} is this directory's params/; {cases} is still "
            "../v2/params/. ../run.py runs these after v1's and v2's.",
            "The until rule: a case whose `until` step is at or before the step run through is not run; "
            "its successor carries the new expectation. B1's files are frozen, so the dated note for "
            "C027 is here: from B3 (declared 2026-10-01) C027 is retired and C054 succeeds it.",
            "Steps, in order: B0, B1, B2, B3, B3b, B4, B5, B6, B7. B3b is the index calculus on "
            "two-word fields, declared separately.",
        ],
        "generator": {"path": "research/ic_tool_program/conformance/v2-b3/make_cases.py",
                      "sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        "cases": cases,
    }, indent=1, ensure_ascii=False) + "\n"
    out["SHA256SUMS"] = "".join(f"{hashlib.sha256(t.encode()).hexdigest()}  {rel}\n"
                                for rel, t in sorted(out.items()))
    return out


def main() -> None:
    out = texts()
    if "--check" in sys.argv:
        bad = [rel for rel, t in out.items() if not (HERE / rel).exists() or (HERE / rel).read_text() != t]
        print(json.dumps({"files": len(out), "mismatches": bad}, indent=1))
        raise SystemExit(1 if bad else 0)
    if (HERE / "cases.json").exists():
        raise SystemExit("cases.json exists; B3's cases are frozen (use --check)")
    for rel, t in out.items():
        path = HERE / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "x") as fh:
            fh.write(t)
    print(f"{len(out)} files written")


if __name__ == "__main__":
    main()

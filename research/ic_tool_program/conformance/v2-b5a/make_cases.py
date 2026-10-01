#!/usr/bin/env python3
"""Conformance cases for B5a: extension fields `GF(p^k)`
(`research/ic_tool_program/rounds/B5a-extension-fields/PROTOCOL.md`).

    python3 make_cases.py           # writes params/, cases.json and SHA256SUMS; refuses to overwrite
    python3 make_cases.py --check   # re-derives every file and compares the bytes

The arithmetic is B5a's instance generator's
(`../../rounds/B5a-extension-fields/instances.py`): `GF(p^k)` and its
curves in Python, sharing nothing with the Rust tool.  Each instance in
its `instances.json` is re-checked here before a document is written:
- the modulus is irreducible, by Rabin's test;
- the curve is non-singular;
- the generator is `[h]P` for a point hashed from a public label: not the
  identity, with `[r]G = O` and `r` prime by B1's exact test;
- `r > 4√q`, so `#E` is the only multiple of `r` in the Hasse interval,
  and that multiple is the recorded order `h·r` (design §4.5's method 3;
  the search found the order independently, by baby steps and giant
  steps);
- `r²` does not divide `#E`;
- the ICV1 slug is the reference implementation's (`scripts/curve_id.py`).

Each known-answer target is a known multiple of the generator, its
multiplier hashed from a public label; each public point `T001` is `[h]`
times a hashed point, and nobody knows its logarithm.

`../run.py` runs every step's cases, these among them, with the `until`
and `supersedes` rules applied.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
B2_DIR = HERE.parent / "v2-b2"
ROUND = HERE.parents[1] / "rounds" / "B5a-extension-fields"
_spec = importlib.util.spec_from_file_location("b5a_instances", ROUND / "instances.py")
inst = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(inst)
b2, curve_id = inst.b2, inst.curve_id

STEP = "B5a"
LABEL = "ic-conformance-v2-b5a"
RHO_SEED = 0x230000 + 1  # the suite's rho seed for target T01, as B1 to B4 use
TIMEOUT = 1800
C050_DOC = "C050-cubic-extension.json"  # B2's, copied byte for byte (C103, C107)
C050_KNOWN_LOG = "554837"


def arr(u) -> list[str]:
    return [str(c) for c in u]


def known_log(iid: str, r: int) -> int:
    digest = hashlib.sha256(f"{LABEL}/{iid}/known_log".encode()).digest()
    return int.from_bytes(digest, "big") % (r - 1) + 1


def subgroup_point(E, label: str, h: int, r: int):
    """[h]P for the first labelled P with [h]P not the identity; checked to have order r."""
    j = 0
    while True:
        G = E.mul(h, E.point(b2.labelled(label, j)))
        if G is not None:
            assert E.on_curve(G) and E.mul(r, G) is None
            return G
        j += 1


def checked(rec: dict):
    """An instance re-checked (module docstring): its field, curve and generator."""
    p, k = int(rec["p"]), rec["k"]
    mod = [int(c) for c in rec["modulus"]]
    order, r, h = int(rec["order"]), int(rec["r"]), int(rec["h"])
    assert inst.rabin_irreducible(p, mod), rec["id"]
    F = inst.Fpk(p, mod)
    E = inst.Curve(F, tuple(int(x) for x in rec["a"]), tuple(int(x) for x in rec["b"]))
    assert not E.singular(), rec["id"]
    assert order == h * r and b2.PRIME_TEST(r) and order % (r * r), rec["id"]
    G = subgroup_point(E, f"{LABEL}/{rec['id']}/generator", h, r)
    q = F.q
    s = math.isqrt(q)
    lo, hi = q + 1 - 2 * s - 2, q + 1 + 2 * s + 2
    assert r > 4 * s + 4 and [m for m in range(lo - lo % r, hi + 1, r) if m >= lo] == [order], rec["id"]
    ident = curve_id.extension_id(p, k, mod, list(E.a), list(E.b), order)
    assert ident["slug"] == rec["slug"] and ident["icv1"] == rec["icv1"], rec["id"]
    return F, E, G


def doc(name: str, rec: dict, G, target: dict, method: dict, curve: dict | None = None) -> dict:
    assert len(name) <= 120, name
    return {"schema_version": 2, "name": name,
            "field": {"kind": "prime_extension", "p": rec["p"], "degree": rec["k"], "modulus": rec["modulus"]},
            "curve": curve or {"form": "short_weierstrass", "a": rec["a"], "b": rec["b"]},
            "subgroup": {"order": rec["r"], "cofactor": rec["h"], "generator": {"x": arr(G[0]), "y": arr(G[1])}},
            "target": target, "method": method}


def case(cid: str, purpose: str, files: dict, argv: list[str], expect: dict, timeout: int = TIMEOUT,
         **extra) -> dict:
    return {"id": cid, "step": STEP, "checks": purpose, "files": files, "argv": argv,
            "env": {"RAYON_NUM_THREADS": "1"}, "expect": expect, "timeout_s": timeout, **extra}


def price(params: str = "{tmp}/params.json", report: str = "{tmp}/report.json") -> list[str]:
    return ["price", "--params", params, "--json", "--out", report, "--repeats", "1", "--repeats-fast", "1"]


def check(params: str = "{tmp}/params.json", report: str = "{tmp}/report.json") -> list[str]:
    return ["check", "--params", params, "--json", "--out", report]


def general_form(F, E, label: str):
    """`E` written as y² + a1xy + a3y = x³ + a2x² + a4x + a6 by x = X + s, y = Y + λX + μ,
    with s, λ, μ hashed from `label`, and the map of points to it."""
    s, lam, mu = (F.element(f"{label}/{v}") for v in ("s", "lambda", "mu"))
    two, three = F.scalar(2), F.scalar(3)
    a1, a3 = F.mul(two, lam), F.mul(two, mu)
    a2 = F.sub(F.mul(three, s), F.mul(lam, lam))
    a4 = F.sub(F.add(F.mul(three, F.mul(s, s)), E.a), F.mul(two, F.mul(lam, mu)))
    a6 = F.sub(F.add(F.add(F.mul(s, F.mul(s, s)), F.mul(E.a, s)), E.b), F.mul(mu, mu))

    def to_general(P):
        X = F.sub(P[0], s)
        Y = F.sub(F.sub(P[1], F.mul(lam, X)), mu)
        lhs = F.add(F.add(F.mul(Y, Y), F.mul(a1, F.mul(X, Y))), F.mul(a3, Y))
        rhs = F.add(F.add(F.add(F.mul(X, F.mul(X, X)), F.mul(a2, F.mul(X, X))), F.mul(a4, X)), a6)
        assert lhs == rhs, "the general form's equation"
        return X, Y

    curve = {"form": "general_weierstrass", "a1": arr(a1), "a2": arr(a2), "a3": arr(a3), "a4": arr(a4),
             "a6": arr(a6)}
    return curve, to_general


def build() -> tuple[dict[str, dict | str], list[dict]]:
    recs = {r["id"]: r for r in json.loads((ROUND / "instances.json").read_text())["instances"]}
    files: dict[str, dict | str] = {}
    cases: list[dict] = []
    paired, rho = b2.paired_auto(RHO_SEED), b2.rho_alone(RHO_SEED)
    known: dict[str, tuple] = {}

    # Every instance: its known-answer document, and measurement 5's public point.
    for iid, rec in recs.items():
        F, E, G = checked(rec)
        r, h = int(rec["r"]), int(rec["h"])
        k = known_log(iid, r)
        method = paired if iid.startswith("G") else rho
        files[f"{iid}-known.json"] = doc(f"B5a {iid}: {rec['slug']}, a known logarithm's multiple", rec, G,
                                         {"known_log": str(k)}, method)
        T = subgroup_point(E, f"{LABEL}/{iid}/T001", h, r)
        files[f"{iid}-T001.json"] = doc(f"B5a {iid}-T001: {rec['slug']}, a public point for measurement 5",
                                        rec, G, {"point": {"x": arr(T[0]), "y": arr(T[1])}}, method)
        known[iid] = (F, E, G, k)

    def considered(*rows):
        return [{"pipeline": p, "admitted": a, **({"gate": g} if g else {})} for p, a, g in rows]

    gaudry_rho = considered(("ic-gaudry-cubic", True, None), ("rho-negation", True, None))
    complete = {"status": "complete", "result.verified": True}

    # C103: C050's successor; C050 names B5 as its `until`, and B5 is now two steps.
    files[C050_DOC] = (B2_DIR / "params" / C050_DOC).read_text()
    cases.append(case(
        "C103-cubic-extension-rho", "C050's successor from B5a: its GF(1009^3) document, with a general "
        "modulus and h = 1524, recovers its known logarithm with rho-negation",
        {"params.json": {"copy": f"{{here}}/{C050_DOC}"}},
        ["price", "--params", "{tmp}/params.json", "--json", "--out", "{tmp}/report.json"],
        {"exit": 0, "json_file": "{tmp}/report.json",
         "json_paths": {**complete, "route.rho.pipeline": "rho-negation", "result.scalar": C050_KNOWN_LOG,
                        "result.known_answer": True}},
        timeout=300, supersedes="C050-cubic-extension"))

    # C104-C106: Gaudry's index calculus at F0, paired with rho-negation.
    for cid, iid in (("C104", "G1"), ("C105", "G2"), ("C106", "G3")):
        rec, k = recs[iid], known[iid][3]
        cases.append(case(
            f"{cid}-gaudry-cubic-{iid.lower()}", f"ic-gaudry-cubic paired with rho-negation at F0 on {iid}, "
            f"GF({rec['p']}^3) with t^3 - c and a prime group order of 2^{rec['log2_r']}",
            {"params.json": {"copy": f"{{here}}/{iid}-known.json"}}, price(),
            {"exit": 0, "json_file": "{tmp}/report.json",
             "json_paths": {**complete, "fidelity": "F0", "result.scalar": str(k), "result.known_answer": True,
                            "all_verified": True, "ic_and_rho_agree": True,
                            "route.ic.pipeline": "ic-gaudry-cubic", "route.rho.pipeline": "rho-negation",
                            "speedup_eligible": False, "curve_id.slug": rec["slug"]},
             "json_contains": {"route.considered": gaudry_rho, "disclosures": [{"code": "study-pipeline"}]}}))

    # C107-C109: the index calculus's gates, each with rho admitted.
    for cid, purpose, files_, gate in (
            ("C107-cubic-modulus-not-binomial", "C050's document paired: its general cubic modulus is refused by "
             "ic-gaudry-cubic as modulus-not-binomial, and rho-negation admits it",
             {"params.json": {"copy": f"{{here}}/{C050_DOC}", "set": {"method": b2.paired_auto(7)}}},
             "modulus-not-binomial"),
            ("C108-gaudry-cofactor-not-one", "H1 paired: a cofactor of 395 is refused by ic-gaudry-cubic as "
             "cofactor-not-one, and rho-negation admits it",
             {"params.json": {"copy": "{here}/H1-known.json", "set": {"method": paired}}}, "cofactor-not-one"),
            ("C109-quadratic-extension-paired", "E2 paired: GF(p^2) is refused by ic-gaudry-cubic as "
             "extension-degree-not-three, and rho-negation admits it",
             {"params.json": {"copy": "{here}/E2-known.json", "set": {"method": paired}}},
             "extension-degree-not-three")):
        cases.append(case(
            cid, purpose, files_, price(),
            {"exit": 3, "json_file": "{tmp}/report.json",
             "json_paths": {"refusal.code": "no-ic-route", "refusal.class": "unsupported"},
             "json_contains": {"route.considered": considered(("ic-gaudry-cubic", False, gate),
                                                              ("rho-negation", True, None))}},
            timeout=300))

    # C110-C113: rho alone on every width and degree.
    for cid, iid, purpose, pipe, extra in (
            ("C110", "E2", "the one-word arithmetic at its widest, q = (2^31 - 1)^2", "rho-negation", []),
            ("C111", "E5", "a degree-5 field, q ≈ 2^55", "rho-negation", []),
            ("C112", "E11", "a degree-11 field, past eight coefficients", "rho-negation", []),
            ("C113", "B2", "q ≈ 2^70, past one word: rho-negation is refused and rho-bignum runs", "rho-bignum",
             [("rho-negation", False, "field-wider-than-one-word"), ("rho-bignum", True, None)])):
        rec, k = recs[iid], known[iid][3]
        expect = {"exit": 0, "json_file": "{tmp}/report.json",
                  "json_paths": {**complete, "route.rho.pipeline": pipe, "result.scalar": str(k),
                                 "result.known_answer": True, "curve_id.slug": rec["slug"]}}
        if extra:
            expect["json_contains"] = {"route.considered": considered(*extra)}
        cases.append(case(
            f"{cid}-rho-{iid.lower()}", f"{iid} under solve: rho: {purpose}, the logarithm recovered and verified",
            {"params.json": {"copy": f"{{here}}/{iid}-known.json"}},
            ["price", "--params", "{tmp}/params.json", "--json", "--out", "{tmp}/report.json"], expect))

    # C114: G1's curve in general Weierstrass form, with both points mapped to it.
    F, E, G, k = known["G1"]
    curve, to_general = general_form(F, E, f"{LABEL}/G1/general")
    Gg, Qg = to_general(G), to_general(E.mul(k, G))
    files["C114-general-weierstrass-g1.json"] = doc(
        "C114: G1's curve in general Weierstrass form, its points mapped", recs["G1"], Gg,
        {"point": {"x": arr(Qg[0]), "y": arr(Qg[1])}}, rho, curve)
    cases.append(case(
        "C114-general-weierstrass-extension", "a general Weierstrass curve over GF(271^3), isomorphic to G1's: "
        "converted, the conversion recorded, and the logarithm of the mapped target recovered by rho",
        {"params.json": {"copy": "{here}/C114-general-weierstrass-g1.json"}},
        ["price", "--params", "{tmp}/params.json", "--json", "--out", "{tmp}/report.json"],
        {"exit": 0, "json_file": "{tmp}/report.json",
         "json_paths": {**complete, "result.scalar": str(k), "route.rho.pipeline": "rho-negation",
                        "conversion.from": "general_weierstrass", "conversion.to": "short_weierstrass"}}))

    # C115, C116: G1 under check, and at F2.
    cases.append(case(
        "C115-gaudry-check", "G1 under check: valid, with ic-gaudry-cubic and rho-negation admitted",
        {"params.json": {"copy": "{here}/G1-known.json"}}, check(),
        {"exit": 0, "json_file": "{tmp}/report.json", "json_paths": {"status": "checks_passed"},
         "json_contains": {"route.considered": gaudry_rho}}, timeout=120))
    cases.append(case(
        "C116-gaudry-estimates", "G1 at fidelity F2: estimated, with an estimate for each arm",
        {"params.json": {"copy": "{here}/G1-known.json", "set": {"method.fidelity": "F2"}}}, price(),
        {"exit": 0, "json_file": "{tmp}/report.json",
         "json_paths": {"status": "estimated", "estimate.ic.pipeline": "ic-gaudry-cubic",
                        "estimate.rho.pipeline": "rho-negation"}}, timeout=120))

    # C117: the pipeline named on a binary document.
    cases.append(case(
        "C117-gaudry-on-a-binary-field", "ic-gaudry-cubic named on the smoke row's binary document: refused "
        "with its gate, no-pipeline-for-field",
        {"params.json": {"copy": "{cases}/C009-translation-of-smoke-row.json",
                         "set": {"method.index_calculus": {"pipeline": "ic-gaudry-cubic", "recipe": "auto"}}}},
        price(),
        {"exit": 3, "json_file": "{tmp}/report.json",
         "json_paths": {"refusal.code": "no-pipeline-for-field", "refusal.class": "unsupported"}}, timeout=120))

    # C118: Gaudry's S4 solve by recipe.
    k = known["G1"][3]
    cases.append(case(
        "C118-gaudry-groebner-oracle", "G1 with recipe {oracle: groebner}: Gaudry's symmetrised S4 solve "
        "recovers the logarithm, verified with the matched rho's",
        {"params.json": {"copy": "{here}/G1-known.json",
                         "set": {"method.index_calculus": {"pipeline": "ic-gaudry-cubic",
                                                           "recipe": {"oracle": "groebner"}}}}}, price(),
        {"exit": 0, "json_file": "{tmp}/report.json",
         "json_paths": {**complete, "result.scalar": str(k), "all_verified": True,
                        "route.ic.pipeline": "ic-gaudry-cubic", "ic.recipe.oracle": "groebner"}}))
    cases.sort(key=lambda c: c["id"])
    return files, cases


def texts() -> dict[str, str]:
    curve_id.selftest()
    files, cases = build()
    out = {f"params/{name}": (d if isinstance(d, str) else json.dumps(d, indent=1) + "\n")
           for name, d in files.items()}
    out["cases.json"] = json.dumps({
        "suite": "ic tool programme conformance suite v2, B5a's cases",
        "design": "research/ic_tool_program/design/extension-fields.md; "
                  "research/ic_tool_program/rounds/B5a-extension-fields/PROTOCOL.md",
        "includes": "../v1/cases.json (step B0), ../v2/cases.json (B1) and the later steps' sets, run first, "
                    "with the until and supersedes rules",
        "rules": [
            "The rules of ../v2/cases.json hold. {here} is this directory's params/; {cases} is still "
            "../v2/params/. ../run.py runs these after the earlier steps' cases.",
            "The supersedes rule: C050 names B5 as its `until`, and B5 is now two steps, B5a and B5b, so C103 "
            "supersedes C050 from B5a on.",
            "params/C050-cubic-extension.json is B2's file, copied byte for byte for C103 and C107.",
            "params/*-T001.json are no case's: they are B5a's measurement 5's public points, frozen here with "
            "the instances' known-answer documents. Every instance comes from "
            "../../rounds/B5a-extension-fields/instances.json, re-checked by this generator.",
        ],
        "generator": {"path": "research/ic_tool_program/conformance/v2-b5a/make_cases.py",
                      "sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        "instances": {"path": "research/ic_tool_program/rounds/B5a-extension-fields/instances.json",
                      "sha256": hashlib.sha256((ROUND / "instances.json").read_bytes()).hexdigest()},
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
        raise SystemExit("cases.json exists; B5a's cases are frozen (use --check)")
    for rel, t in out.items():
        path = HERE / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "x") as fh:
            fh.write(t)


if __name__ == "__main__":
    main()

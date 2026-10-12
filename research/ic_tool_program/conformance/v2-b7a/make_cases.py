#!/usr/bin/env python3
"""Conformance cases for B7a: the F1 sampled level (C071-C077;
`research/ic_tool_program/rounds/B7a-f1-sampled/PROTOCOL.md`).

    python3 make_cases.py           # writes cases.json and SHA256SUMS; refuses to overwrite
    python3 make_cases.py --check   # re-derives every file and compares the bytes

Written at B7a's declaration, before any B7a code.  Every case copies a
document an earlier step froze and changes its method, so this set
needs no arithmetic of its own: the earlier generators asserted each
document's properties.  An F1 figure varies from run to run, so the
cases check the report's shape and its codes, not its numbers.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CONF = HERE.parent
STEP = "B7a"
LABEL = "an extrapolation from samples (F1); not a measurement; it does not discharge the m = 83 gate"

SOURCES = {
    "{cases}/C010-curve-a-other-generator.json": CONF / "v2" / "params" / "C010-curve-a-other-generator.json",
    "{cases}/gate-m83-T001.json": CONF / "v2" / "params" / "gate-m83-T001.json",
    "{here}/../../v2-b2/params/C047-secp256k1-over-budget.json":
        CONF / "v2-b2" / "params" / "C047-secp256k1-over-budget.json",
    "{here}/../../v2-b2/params/C038-generic-binary-n23.json":
        CONF / "v2-b2" / "params" / "C038-generic-binary-n23.json",
    "{here}/../../v2-b2b/params/C059-subfield-curve-paired.json":
        CONF / "v2-b2b" / "params" / "C059-subfield-curve-paired.json",
}


def build() -> list[dict]:
    for path in SOURCES.values():
        assert path.is_file(), f"{path} is an earlier step's frozen document"
    c010 = json.loads(SOURCES["{cases}/C010-curve-a-other-generator.json"].read_text())
    cases: list[dict] = []

    def case(cid: str, purpose: str, source: str, set_: dict, expect: dict, timeout: int = 600) -> None:
        cases.append({
            "id": cid, "step": STEP, "checks": purpose,
            "files": {"params.json": {"copy": source, "set": set_}},
            "argv": ["price", "--params", "{tmp}/params.json", "--json", "--out", "{tmp}/report.json"],
            "env": {"RAYON_NUM_THREADS": "1"},
            "expect": {"json_file": "{tmp}/report.json", **expect},
            "timeout_s": timeout,
        })

    def extrapolated(paths: dict) -> dict:
        return {"exit": 0, "json_paths": {"status": "extrapolated", "fidelity": "F1", "label": LABEL,
                                          "speedup_eligible": False, **paths}}

    c010_doc = "{cases}/C010-curve-a-other-generator.json"
    case("C071-curve-a-paired-f1",
         "B7a: kic and rho-koblitz extrapolated from samples on curve A; the yield counted, the costs sampled",
         c010_doc, {"method.fidelity": "F1"},
         extrapolated({"phases.ic.collect.how": "counted_and_sampled", "extrapolated.ic.complete": True,
                       "extrapolated.rho.complete": True, "route.ic.pipeline": "kic",
                       "route.rho.pipeline": "rho-koblitz"}))
    case("C072-curve-a-rho-f1", "B7a: rho alone at F1, its step cost measured and its expected steps counted",
         c010_doc, {"method": {"solve": "rho", "fidelity": "F1", "rho": {"pipeline": "auto", "seed": 7}}},
         extrapolated({"rho.pipeline": "rho-koblitz", "extrapolated.rho.complete": True}))
    case("C073-curve-a-kic-alone-f1", "B7a: kic alone at F1",
         c010_doc, {"method": {"solve": "index_calculus", "fidelity": "F1",
                               "index_calculus": c010["method"]["index_calculus"]}},
         extrapolated({"operation": "ic_single_target", "extrapolated.ic.complete": True}))
    case("C074-gate-rho-f1", "B7a: rho-koblitz at F1 on the m = 83 gate's field (an extrapolation, never the gate)",
         "{cases}/gate-m83-T001.json",
         {"method": {"solve": "rho", "fidelity": "F1", "rho": {"pipeline": "auto", "seed": 2293761}}},
         extrapolated({"rho.pipeline": "rho-koblitz", "extrapolated.rho.complete": True}))
    case("C075-secp256k1-rho-f1", "B7a: rho-bignum at F1 on secp256k1, where F0 is over any budget",
         "{here}/../../v2-b2/params/C047-secp256k1-over-budget.json", {"method.fidelity": "F1"},
         extrapolated({"rho.pipeline": "rho-bignum", "extrapolated.rho.complete": True}))
    case("C076-imported-ic-has-no-f1", "B7a: ic-binary-s4 has no F1 model",
         "{here}/../../v2-b2/params/C038-generic-binary-n23.json", {"method.fidelity": "F1"},
         {"exit": 3, "json_paths": {"refusal.code": "no-f1-model", "refusal.class": "unsupported"}})
    case("C077-subfield-curve-f1", "B7a: kic on a curve over GF(4), and rho-negation, at F1",
         "{here}/../../v2-b2b/params/C059-subfield-curve-paired.json", {"method.fidelity": "F1"},
         extrapolated({"route.ic.pipeline": "kic", "route.rho.pipeline": "rho-negation",
                       "extrapolated.ic.complete": True, "extrapolated.rho.complete": True}))
    return cases


def texts() -> dict[str, str]:
    cases = build()
    out = {"cases.json": json.dumps({
        "suite": "ic tool programme conformance suite v2, B7a's cases",
        "design": "research/ic_tool_program/design/f1-sampled.md (with amendment 1); "
                  "research/ic_tool_program/rounds/B7a-f1-sampled/PROTOCOL.md",
        "includes": "../v1 (B0), ../v2 (B1), ../v2-b2 (B2), ../v2-b2b (B2b), ../v2-b3 (B3); ../run.py runs every step's cases",
        "rules": [
            "The rules of ../v2/cases.json hold. {here} is this directory's params/, which is empty: every "
            "case copies an earlier step's frozen document. {cases} is ../v2/params/.",
            "An F1 figure varies from run to run, so the expectations are the report's codes and shape: "
            "status, fidelity, the label, the phases' kinds and the completeness flags.",
        ],
        "sources_sha256": {k: hashlib.sha256(v.read_bytes()).hexdigest() for k, v in SOURCES.items()},
        "generator": {"path": "research/ic_tool_program/conformance/v2-b7a/make_cases.py",
                      "sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        "cases": cases,
    }, indent=1, ensure_ascii=False) + "\n"}
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
        raise SystemExit("cases.json exists; B7a's cases are frozen (use --check)")
    for rel, t in out.items():
        with open(HERE / rel, "x") as fh:
            fh.write(t)
    print(f"{len(out)} files written")


if __name__ == "__main__":
    main()

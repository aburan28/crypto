#!/usr/bin/env python3
"""Conformance cases for B3b: the index calculus on two-word binary fields
(`research/ic_tool_program/rounds/B3b-two-word-kic/PROTOCOL.md`).

    python3 make_cases.py           # writes params/, cases.json and SHA256SUMS; refuses to overwrite
    python3 make_cases.py --check   # re-derives every file and compares the bytes

The documents are built with B1's generator code (`../v2/make_cases.py`),
whose arithmetic shares nothing with the Rust tool:
- each curve's order comes from the trace recurrence, and `r` and `h` are
  stated here and checked against it, with `r`'s primality exact;
- each generator is `[h]P` for a point hashed from a public label;
- each target is a known multiple of the generator, the multiplier hashed
  from a public label, so every case checks an exact answer.

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

STEP = "B3b"
RHO_SEED = 0x230000 + 1  # the suite's rho seed for target T01, as B1 and B3 use
GATE_MODULUS = (1 << 83) | (1 << 45) | (1 << 2) | (1 << 1) | 1  # AGENTS.md §8a
LABEL = "ic-conformance-v2-b3b"
GATE = "{cases}/gate-m83-T001.json"
SMOKE = "{cases}/C009-translation-of-smoke-row.json"
ROW61 = "params/S/k0n61/M1-T01.json"
TIMEOUT = 1800

# (case id, a, n, modulus or None for the repository's rule, r, h, purpose)
INSTANCES = [
    ("C078", 0, 67, None, 76589041, 1926828573412, "the smallest two-word instance"),
    ("C079", 1, 67, None, 68352708293, 2159006662, "a 36-bit subgroup at n = 67"),
    ("C080", 0, 79, None, 14377452373, 42042421294532, "n = 79, a 41-bit cofactor"),
    ("C081", 1, 71, None, 588353361747061, 4013206, "a 49-bit subgroup at n = 71"),
    ("C082", 1, 83, GATE_MODULUS, 8569786107849059, 1128547018,
     # Amendment 2: shortened, since the name it makes must stay within
     # schema v2's 120 characters.
     "the gate's field and modulus; K_1's 53-bit subgroup, not the gate's curve"),
]


def case(cid: str, purpose: str, files: dict, argv: list[str], expect: dict, timeout: int = TIMEOUT,
         until: str | None = None, supersedes: str | None = None) -> dict:
    out = {"id": cid, "step": STEP, "checks": purpose, "files": files, "argv": argv,
           "env": {"RAYON_NUM_THREADS": "1"}, "expect": expect, "timeout_s": timeout}
    if until:
        out["until"] = until
    if supersedes:
        out["supersedes"] = supersedes
    return out


def price(params: str = "{tmp}/params.json", report: str = "{tmp}/report.json") -> list[str]:
    return ["price", "--params", params, "--json", "--out", report, "--repeats", "1", "--repeats-fast", "1"]


def known_log(slug: str, r: int) -> int:
    digest = hashlib.sha256(f"{LABEL}/{slug}/known_log".encode()).digest()
    return int.from_bytes(digest, "big") % (r - 1) + 1


def build() -> tuple[dict[str, dict], list[dict]]:
    files: dict[str, dict] = {}
    cases: list[dict] = []
    for cid, a, n, modulus, r, h, purpose in INSTANCES:
        f = modulus if modulus is not None else v2.curve_id.find_irreducible_sparse(n)
        assert v2.irreducible(f)
        curve = v2.Curve(n, f, a, 1)
        order = v2.koblitz_order(a, n)
        assert order == r * h and v2.prime_exact(r) and h % r, (a, n)
        ids = v2.curve_id.binary_id(n, f, a, 1, order, end="-7")
        ident, slug = ids["icv1"], ids["slug"]
        gen = curve.subgroup_point(f"{LABEL}/{ident}/generator", h, r)
        k = known_log(ident, r)
        name = f"{cid}: {slug}, {purpose}"
        fname = f"{cid}-kic-two-word-n{n}.json"
        files[fname] = v2.document(name, curve, {"form": "koblitz", "a": a}, r, h, gen,
                                   {"known_log": str(k)}, "auto", RHO_SEED, routed=True)
        # Measurement 5's public target: a point hashed from a public
        # label, as the gate's T001 is, whose logarithm nobody knows.
        t001 = curve.subgroup_point(f"{LABEL}/{ident}/T001", h, r)
        files[f"{cid}-T001-n{n}.json"] = v2.document(
            f"{cid}-T001: {slug}, a public target for B3b's measurement 5", curve, {"form": "koblitz", "a": a},
            r, h, gen, {"point": v2.point_doc(t001)}, "auto", RHO_SEED, routed=True)
        cases.append(case(
            f"{cid}-kic-two-word-n{n}", f"B3b's kic at F0 on a two-word field: {purpose}",
            {"params.json": {"copy": f"{{here}}/{fname}"}}, price(),
            {"exit": 0, "json_file": "{tmp}/report.json",
             "json_paths": {"status": "complete", "fidelity": "F0", "result.scalar": str(k),
                            "result.verified": True, "result.known_answer": True, "all_verified": True,
                            "ic_and_rho_agree": True, "curve_id.icv1": ident, "ic.words": 2},
             "json_contains": {"route.considered": [{"pipeline": "kic", "admitted": True},
                                                    {"pipeline": "rho-koblitz", "admitted": True}]}}))
    # C083: C054's successor (B3 marked C054 `until: B3b`). With two-word
    # kernels kic admits the gate curve by width; at F0 its estimate at
    # r = 2^81 exceeds any budget, so a paired price is refused as over
    # budget, naming F1.
    cases.append(case(
        "C083-gate-paired-over-budget", "C054's successor from B3b: kic and rho-koblitz both admit the gate "
        "curve by width; F0's estimate exceeds the budget, so the paired price is refused as over budget",
        {"params.json": {"copy": GATE}}, price(),
        {"exit": 4, "json_file": "{tmp}/report.json",
         "json_paths": {"refusal.code": "over-budget", "refusal.class": "over_budget"},
         "json_contains": {"route.considered": [{"pipeline": "kic", "admitted": True},
                                                {"pipeline": "rho-koblitz", "admitted": True}]}},
        timeout=300))

    # C084, C085: sameness end to end. `--kic-wide` runs the two-word
    # pipeline where the one-word one would; the outputs must agree.
    compared = ["counts", "certificates.ic.scalar", "certificates.rho.scalar", "rho_counts", "all_verified",
                "ic_and_rho_agree"]
    cases.append(case(
        "C084-smoke-two-word-is-one-word", "the two-word pipeline on the smoke row's translation gives the "
        "one-word pipeline's counts, logarithm and certificates",
        {"params.json": {"copy": SMOKE}}, price() + ["--kic-wide"],
        {"exit": 0, "json_file": "{tmp}/report.json",
         "json_paths": {"status": "complete", "ic.words": 2, "all_verified": True},
         "same_outputs_as": {"argv": price("{tmp}/params.json", "{tmp}/one-word.json"),
                             "json_file": "{tmp}/one-word.json", "paths": compared}},
        timeout=300))
    # C085: the same at n = 61, the suite's largest size and a wide-tail
    # field, on the translation of its v1 row (built as B3's C058 is).
    row = v2.v1_row(ROW61)
    n = row["curve"]["degree"]
    assert (n, row["curve"]["curve_a"]) == (61, 0)
    f = v2.curve_id.find_irreducible_sparse(n)
    curve = v2.Curve(n, f, 0, 1)
    order = v2.koblitz_order(0, n)
    r = v2.largest_prime_factor(order)
    h = order // r
    assert v2.prime_exact(r) and r > h and h % r
    assert v2.curve_id.binary_id(n, f, 0, 1, order, end="-7")["slug"] == "icv1-f2m61-t158598901-ab42b6c5"
    files["C085-n61-paired.json"] = v2.document(
        f"C085: suite v1 row {ROW61}, translated to v2, paired", curve, {"form": "koblitz", "a": 0}, r, h,
        {"rule": "koblitz_search_v1"}, {"public_hash_seed": row["targets"][0]["public_hash_seed"]},
        v2.recipe_of(row), RHO_SEED)
    cases.append(case(
        "C085-n61-two-word-is-one-word", "C084 at n = 61, the suite's largest size, on its v1 row's translation",
        {"params.json": {"copy": "{here}/C085-n61-paired.json"}}, price() + ["--kic-wide"],
        {"exit": 0, "json_file": "{tmp}/report.json",
         "json_paths": {"status": "complete", "ic.words": 2, "all_verified": True},
         "same_outputs_as": {"argv": price("{tmp}/params.json", "{tmp}/one-word.json"),
                             "json_file": "{tmp}/one-word.json", "paths": compared}},
        timeout=900))

    # C086, C087 (amendment 2, 2026-10-01): two earlier cases read kic's
    # width gate past one word, which B3b moves. C031's `until` names B4
    # and C053 has none, so the until rule cannot retire them; each is
    # superseded from B3b instead, with the rest of its expectation kept.
    cases.append(case(
        "C086-challenge-kic-past-two-words", "C031's successor from B3b: the challenge file is valid, and kic, "
        "with two-word kernels, refuses n = 131 as wider than two words, as rho-koblitz has since B3",
        {"params.json": {"copy": "{cases}/ecc2k130-challenge.json"}},
        ["check", "--params", "{tmp}/params.json", "--json", "--out", "{tmp}/report.json"],
        {"exit": 0, "json_file": "{tmp}/report.json", "json_paths": {"status": "checks_passed"},
         "json_contains": {"checks": [{"code": "order-composite", "status": "pass", "exact": False}],
                           "disclosures": [{"code": "primality-screen"}],
                           "route.considered": [{"pipeline": "kic", "admitted": False,
                                                 "gate": "field-wider-than-two-words"}]}},
        timeout=120, until="B4", supersedes="C031-challenge-checks"))
    cases.append(case(
        "C087-gate-rho-capped-kic-admitted", "C053's successor from B3b: rho alone on the gate curve, stopped by "
        "its step cap, reports no logarithm and exits 1; kic now admits the gate's field by width",
        {"params.json": {"copy": GATE, "set": {"method": {
            "solve": "rho", "fidelity": "F0",
            "rho": {"pipeline": "auto", "seed": RHO_SEED, "max_iterations": 200000}}}}},
        ["price", "--params", "{tmp}/params.json", "--json", "--out", "{tmp}/report.json"],
        {"exit": 1, "json_file": "{tmp}/report.json",
         "json_paths": {"status": "not_recovered", "result.verified": False, "rho.step_cap": 200000,
                        "rho.counters.walk_operations": 200000, "rho.pipeline": "rho-koblitz"},
         "json_contains": {"route.considered": [{"pipeline": "kic", "admitted": True},
                                                {"pipeline": "rho-koblitz", "admitted": True}]}},
        timeout=120, supersedes="C053-gate-rho-capped"))
    return files, cases


def texts() -> dict[str, str]:
    files, cases = build()
    out = {f"params/{name}": json.dumps(doc, indent=1) + "\n" for name, doc in files.items()}
    out["cases.json"] = json.dumps({
        "suite": "ic tool programme conformance suite v2, B3b's cases",
        "design": "research/ic_tool_program/design/two-word-kic.md; "
                  "research/ic_tool_program/rounds/B3b-two-word-kic/PROTOCOL.md",
        "includes": "../v1/cases.json (step B0), ../v2/cases.json (B1) and the later steps' sets, run first, "
                    "with the until rule",
        "rules": [
            "The rules of ../v2/cases.json hold. {here} is this directory's params/; {cases} is still "
            "../v2/params/. ../run.py runs these after the earlier steps' cases.",
            "The until rule: B3 marked C054 `until: B3b`. From B3b (declared 2026-10-01) C054 is retired and "
            "C083 succeeds it. (Amendment 2 filled in this date, which the declaration left as a placeholder.)",
            "The supersedes rule (amendment 2, 2026-10-01): a case may name an earlier step's case that it "
            "supersedes; while the superseding case is run, the superseded one is not. B3b moves kic's width "
            "gate past one word, which C031 (until B4) and C053 (no until) read, so C086 supersedes C031 and "
            "C087 supersedes C053. B1's and B3's files stay frozen, so the note is here.",
            "`--kic-wide` is B3b's test switch: it runs the two-word pipeline where the one-word pipeline "
            "would run. C084 and C085 compare the two through `same_outputs_as`.",
            "params/C078-T001-n67.json to params/C082-T001-n83.json are no case's: they are B3b's "
            "measurement 5's public targets, frozen here with the cases' documents.",
        ],
        "generator": {"path": "research/ic_tool_program/conformance/v2-b3b/make_cases.py",
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
        raise SystemExit("cases.json exists; B3b's cases are frozen (use --check)")
    for rel, t in out.items():
        path = HERE / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "x") as fh:
            fh.write(t)


if __name__ == "__main__":
    main()

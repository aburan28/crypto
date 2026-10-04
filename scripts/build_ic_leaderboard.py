#!/usr/bin/env python3
"""Build the index-calculus leaderboard from frozen evidence.

    python3 scripts/build_ic_leaderboard.py                 # write the three outputs
    python3 scripts/build_ic_leaderboard.py --check         # fail if any is stale
    python3 scripts/build_ic_leaderboard.py --artifact P    # also a publishable copy

Outputs, all derived and never measured here:

- docs/ic/leaderboard.json    the rows, with the SHA-256 of every source file
- docs/ic/LEADERBOARD.md      the same tables in Markdown
- docs/ic-leaderboard.html    the page

The board has no single unit.  It is three tables (A, B, C), each priced in
its own unit against its own reference and floor, and the outputs state
those at the head of each table; see FAMILIES.

Every curve is named by its ICV1 slug (docs/curves/ICV1.md) through the
curve registry.  Every number is read from a frozen file named in
`SOURCES`; the only arithmetic here is averaging the repetitions a report
already holds, dividing a phase's operations by sqrt(r) to express it in S,
and taking a phase's share of its row's total.
"""
from __future__ import annotations

import argparse
import hashlib
import html
import json
import math
import re
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import curve_id as cid  # noqa: E402

REPO = cid.REPO
OUT_JSON = REPO / "docs/ic/leaderboard.json"
OUT_MD = REPO / "docs/ic/LEADERBOARD.md"
OUT_HTML = REPO / "docs/ic-leaderboard.html"

# The ledger section the page is current through (AGENTS.md §7b).  A round
# that adds a later section must update the page, or confirm that the round
# changes nothing on it, and then raise this number; `--check` fails until it
# does, because it cannot otherwise see evidence the builder does not read.
LEDGER = "research/notes/index-calculus/RESEARCH_IC_BOUNDARY_LEDGER.md"
LEDGER_COVERED_THROUGH = 23

SOURCES = {
    "round5": "docs/ic/runs/ic-boundary-ledger-round5-2026-09-22.json",
    "matched_rho": "research/ic_rho_reference_20260923/reprice/ic-boundary-ledger-round5-2026-09-22.json",
    "matched_rho_analysis": "research/ic_rho_reference_20260923/analysis.json",
    "koblitz_s20": "research/ic_exponent_20260926/analysis.json",
    "koblitz_s21": "research/ic_constructions_20260926/analysis.json",
    "koblitz_s22": "research/ic_descent_20260930/analysis-isolated.json",
    "koblitz_s23": "research/ic_single_target_20260930/analysis.json",
    "oracles": "docs/ic/runs/ic-oracle-pricing-lifted-2026-09-21.json",
    "n37_rank_columns": "research/ecbench_n37_rank_columns_20261004/RESULT.json",
    "n37_k8_k16": "research/ecbench_n37_k8_k16_20261004/DECISION.json",
    "n37_online_ir": "research/ecbench_n37_online_ir_20261004/DECISION.json",
    "n37_online_k8_claim": "research/ecbench_n37_online_ir_20261004/candidate_claims/ic-k8-0.json",
    "n37_online_k16_claim": "research/ecbench_n37_online_ir_20261004/candidate_claims/ic-k16-0.json",
    "registry": "docs/curves/registry.json",
}
S22_RUNS = "research/ic_descent_20260930/runs-isolated/main"
S23_RUNS = "research/ic_single_target_20260930/runs"

# Phases, in pipeline order, as the page groups them.
PHASES = [
    ("precomputation", "Reusable set-up", "Base, pair table, relations and logs, built once and reusable across targets (one-target Koblitz rows, §23)"),
    ("setup", "Curve and subgroup", "Point count, factoring #E, generator, cofactor clearing"),
    ("constructions", "Workflow constructions", "Orbit maps, point index, collectors and solvers built around the work (Koblitz rows)"),
    ("factor_base", "Factor base", "Choose and materialise the base points, or their orbit representatives"),
    ("table", "Pair table", "Every pair sum of the base, stored for meet-in-the-middle lookups"),
    ("relations", "Relations (PDP)", "Draw targets and decompose each over the base: the point decomposition problem"),
    ("linear_algebra", "Linear algebra", "Eliminate over Z/rZ until the target's column, or every base log, is pinned"),
    ("descent", "Descent", "Decompose each target over the logged base (Koblitz rows; the ladder puts the target in the matrix)"),
    ("verify", "Verification", "Re-add every relation and check [d]G = Q"),
]

VARIANTS = {
    "semaev_s3_roots_m2": "Semaev S3 roots, m = 2",
    "direct_subtraction_m2": "direct subtraction, m = 2",
    "semaev_s4_pairs_and_solve_m3": "Semaev S4 pairs-and-solve, m = 3",
    "mitm_m2": "pair table, m = 2",
    "mitm_m3": "pair table, m = 3",
    "mitm_m2_negfold": "negation-folded table, m = 2",
    "mitm_m3_negfold": "negation-folded table, m = 3",
    "mitm_m2_negfold_walk": "folded table, walked targets, m = 2",
    "mitm_m3_negfold_walk": "folded table, walked targets, m = 3",
    "mitm_m2_negfold_walk_balanced": "folded, walked, base at the family optimum, m = 2",
}
REGIME_NAME = {"prime": "Prime field", "char2": "Random binary", "koblitz": "Koblitz",
               "koblitz_batch": "Koblitz, 32-target batch"}

# The board mixes three unit families.  Their S, × reference and × floor are
# different quantities: a different unit (a count of group-addition
# equivalents, or a clock reading in batched additions), a different
# reference (one-target rho on the matched instance, one-target rho measured
# as built, or 32-target batch rho) and a different floor (a one-target
# bound in additions, a one-target bound in steps, or its 32-target batch
# share).  Every table, chart and summary below reads inside one family.
# The texts are the unit and reference definitions of the sources named in
# SOURCES: ledger §12 and §18 (table B), §20 and §22 (table C), §23 (table A);
# research/ic_exponent_20260926/PROTOCOL.md and
# research/ic_single_target_20260930/PROTOCOL.md.
FAMILY_OF = {"koblitz": "A", "prime": "B", "char2": "B", "koblitz_batch": "C"}
FAMILIES = {
    "A": {
        "title": "Koblitz, one target", "regimes": ["koblitz"], "section": "§23",
        "status": "primary comparison: one target, as AGENTS.md's one-target rule requires",
        "targets": "one unseen public point per process, 64 per size",
        "unit": "one batched affine addition (`add_many` over 1,024 subgroup points), per √r",
        "unit_how": "A clock reading, not a count: `ic price` times each phase exclusively, on one thread, "
                    "and divides by this unit, measured around each repetition (§20, §23).",
        "reference": "One-target rho on the same point, measured as built (its set-up plus its walk), "
                     "in the same unit and the same process; its answer replayed and verified (§23).",
        "floor": "√(π/4n) rho steps per √r: generic, one target, no precomputation, A = 2n by signed "
                 "Frobenius and negation, n the extension degree. Stated in steps and not rescaled to "
                 "units, so × floor here divides units by steps.",
        "unit_short": "batched-addition units", "ref_short": "one-target rho",
    },
    "B": {
        "title": "Prime and random-binary ladder", "regimes": ["prime", "char2"], "section": "§18",
        "status": "primary comparison: one target per instance",
        "targets": "one planted target per instance",
        "unit": "one group-addition equivalent (GAE), per √r",
        "unit_how": "A count, not a clock: exact native operation counts converted to additions at the "
                    "dimensionless ratios pinned per instance in `docs/ic/calibration.json` (ns per native "
                    "operation over ns per group addition; all 13 instances used the pinned ratios; "
                    "`ic boundary`, ledger §12). No conversion to table A's or C's unit is applied or "
                    "claimed.",
        "reference": "Matched Pollard rho on the same instance, one target: the negation-map r-adding "
                     "walk (A = 2), counted in GAE, verified [d]G = Q (§18).",
        "floor": "√(π/4) ≈ 0.886 GAE per √r: generic, one target, no precomputation, A = 2 usable "
                 "automorphisms (negation).",
        "unit_short": "GAE", "ref_short": "matched negation-map rho",
    },
    "C": {
        "title": "Koblitz, 32-target batch", "regimes": ["koblitz_batch"], "section": "§19–§22",
        "status": "historical diagnostic under AGENTS.md's one-target rule: 32 targets solved together",
        "targets": "32 targets solved together, priced per target",
        "unit": "one batched affine addition (`add_many` over 1,024 subgroup points), per √r, per target",
        "unit_how": "A clock reading, not a count: `ic price` times each phase exclusively, on one thread, "
                    "and divides by this unit (§20, §22 isolated re-run).",
        "reference": "Batch rho on the same 32 targets (signed Frobenius and negation), its counted "
                     "operations priced at the canonical step, one batched addition plus the table "
                     "canonicalisation, in the same process (§20). Collisions shared across the 32 "
                     "targets are in its cost per target.",
        "floor": "L(32)·√(π/4n) steps per √r per target, with L(32) = 0.19869: the one-target floor of "
                 "table A times Kuhn–Struik's batch share. Stated in steps and not rescaled to units, "
                 "and a smaller object than table A's floor, so × floor in C and in A do not compare.",
        "unit_short": "batched-addition units, per target", "ref_short": "32-target batch rho",
    },
}


def rng(lo: float, hi: float) -> str:
    return g3(lo) if lo == hi else f"{g3(lo)}–{g3(hi)}"


def step_line(f: dict) -> str | None:
    """What one rho step costs in the family's unit, where the unit prices a step."""
    parts = []
    if f.get("rho_step_units_as_built"):
        parts.append(f"rho's walk as built costs {rng(*f['rho_step_units_as_built'])} units a step")
    if f.get("canonical_step_units"):
        parts.append(f"the canonical step costs {rng(*f['canonical_step_units'])}")
    return "; ".join(parts) or None


def sha256(rel: str) -> str:
    return hashlib.sha256((REPO / rel).read_bytes()).hexdigest()


def load(rel: str):
    return json.loads((REPO / rel).read_text())


class Names:
    def __init__(self) -> None:
        self.reg = load(SOURCES["registry"])
        self.by = {}
        for c in self.reg["curves"]:
            for a in c["aliases"] + c["standard_names"] + [c["slug"]]:
                self.by[cid.normalise_alias(a)] = c

    def __call__(self, legacy: str) -> dict:
        c = self.by.get(cid.normalise_alias(legacy))
        if c is None:
            raise SystemExit(f"{legacy!r} is not in docs/curves/registry.json")
        return c


def ec1_of(c: dict) -> str | None:
    reps = c.get("representations") or []
    return reps[0]["ec1"] if reps else None


# ---------------------------------------------------------------------------
# The prime and random-binary ladder: Round 5, priced against the matched rho.
# ---------------------------------------------------------------------------


def ladder_rows(names: Names) -> list[dict]:
    r5 = load(SOURCES["round5"])
    mr = {i["instance"]: i for i in load(SOURCES["matched_rho"])["instances"]
          if i["regime"] in ("prime", "char2")}
    rows = []
    for inst in r5["ledger"]["instances"]:
        if inst["regime"] not in ("prime", "char2"):
            continue
        legacy = inst["curve"]["name"]
        c = names(legacy)
        ref = mr[legacy]
        rho_s = ref["matched_reference"]["mean_s"]
        sqrt_r = math.sqrt(inst["r"])
        by_variant: dict[str, list] = {}
        for v in inst["variants"]:
            by_variant.setdefault(v["name"], []).append(v)
        recipes = []
        for name, reps in by_variant.items():
            mean = lambda f: statistics.fmean(f(x) for x in reps)  # noqa: E731
            phases = {
                "factor_base": mean(lambda x: x["factor_base"]["gae"]) / sqrt_r,
                "relations": mean(lambda x: x["relations"]["gae"]) / sqrt_r,
                "linear_algebra": mean(lambda x: x["linear_algebra"]["gae"]) / sqrt_r,
                "verify": mean(lambda x: x["verify"]["gae"]) / sqrt_r,
            }
            s = mean(lambda x: x["s"])
            recipes.append({
                "variant": name, "label": VARIANTS.get(name, name),
                "m": reps[0]["summands"], "base_points": reps[0]["signed_points"],
                "columns": reps[0]["columns"], "table": reps[0]["table"],
                "targets": reps[0]["targets"], "trials": mean(lambda x: x["trials"]),
                "relations": mean(lambda x: x["relations_found"]),
                "yield_over_ceiling": mean(lambda x: x["yield_over_ceiling_exact"]),
                "s": s, "phases_s": phases, "ratio_rho": s / rho_s,
                "ratio_floor": s / math.sqrt(math.pi / 4),
                "verified": all(x["verified"] for x in reps), "repetitions": len(reps),
            })
        recipes.sort(key=lambda x: x["s"])
        rows.append({
            "regime": inst["regime"], "slug": c["slug"], "ec1": ec1_of(c),
            "legacy": legacy, "log2_r": inst["log2_r"], "cofactor": inst["cofactor"],
            "reference": "matched rho, negation map (A = 2), one target",
            "reference_s": rho_s, "floor_s": math.sqrt(math.pi / 4),
            "targets": 1, "recipes": recipes, "best": recipes[0],
            "source": SOURCES["round5"], "round": "Round 5 (2026-09-22), reference §18 (2026-09-23)",
            "class": "accounting",
        })
    return rows


# ---------------------------------------------------------------------------
# The Koblitz thread: §20's sweep at nine sizes, §21's constructions built once.
# ---------------------------------------------------------------------------

KOBLITZ_GROUPS = {
    "setup": ("setup",),
    "constructions": ("select_projection", "collect_setup", "logs_setup", "descent_setup"),
    "factor_base": ("select",),
    "table": ("build",),
    "relations": ("collect",),
    "linear_algebra": ("la",),
    "descent": ("descent",),
    "verify": ("verify", "verify_final", "other"),
}


def koblitz_phase_shares(a: int, n: int) -> tuple[dict, dict]:
    """Shares of each phase in §22's candidate reports (the isolated
    re-run): the median over five rounds within each set, then the mean
    over the four sets."""
    per_set = []
    counts = None
    for s in ("M1", "M2", "M3", "M4"):
        d = REPO / S22_RUNS / f"k{a}n{n}" / s
        reports = sorted(d.glob("r*-candidate.price.json"))
        if not reports:
            continue
        docs = [json.loads(p.read_text()) for p in reports]
        counts = counts or docs[0]["counts"]
        med = {}
        for g, keys in KOBLITZ_GROUPS.items():
            med[g] = statistics.median(sum(x["median"]["phases_units"].get(k, 0.0) for k in keys)
                                       for x in docs)
        tot = sum(med.values())
        per_set.append({g: v / tot for g, v in med.items()})
    shares = {g: statistics.fmean(p[g] for p in per_set) for g in KOBLITZ_GROUPS}
    return shares, counts or {}


def koblitz_rows(names: Names) -> list[dict]:
    """The Koblitz thread on current evidence (§22, isolated), with §21 and
    §22's own baseline as the before marks and §20's one-target cold ratio.

    These rows solve 32 targets together and are priced per target against
    batch rho on the same 32.  AGENTS.md ("Primary comparison uses one
    target") makes that a historical diagnostic, not a primary comparison,
    so each row carries its target count and the one-target figure §20
    measured beside it."""
    s20 = {x["curve"]: x for x in load(SOURCES["koblitz_s20"])["sizes"]}
    s21 = {x["curve"]: x for x in load(SOURCES["koblitz_s21"])["sizes"]}
    rows = []
    for x in load(SOURCES["koblitz_s22"])["sizes"]:
        c = names(x["curve"])
        a, n = x["a"], x["n"]
        before = s20[x["curve"]]
        shares, counts = koblitz_phase_shares(a, n)
        s = x["s_after"]
        floor_s = before["sets"][0]["floor"]  # L(32)·√(π/4n), per target
        logs = counts.get("logs", {})
        la = logs.get("linear_algebra", {})
        recipe = {
            "variant": "koblitz_collection_m3_aimed",
            "label": f"aimed m = 3 collection, m = {before['chosen']['descent_summands']} descent",
            "m": 3, "descent_summands": before["chosen"]["descent_summands"],
            "base_points": before["sets"][0]["points"], "columns": before["chosen"]["columns"],
            "table": before["sets"][0]["tier"], "stored_pairs": before["sets"][0]["stored_pairs"],
            "targets": "32 targets, batch", "relations": before["sets"][0]["relations"],
            "summands_scanned": before["sets"][0]["summands_scanned"],
            "descent_trials_per_target": before["sets"][0]["descent_trials"] / 32,
            "la_rows_in": la.get("filter", {}).get("rows_in"),
            "la_core_dimension": la.get("core_dimension"),
            "trials": counts.get("collect", {}).get("trials"),
            "s": s, "s_before": x["s_before"],
            "phases_s": {g: shares[g] * s for g in shares},
            "ratio_rho": x["ratio_after"],
            "ratio_rho_ci": [x["ratio_after_ci"]["lo"], x["ratio_after_ci"]["hi"]],
            "ratio_rho_before": x["ratio_before"],
            "ratio_rho_section21": s21[x["curve"]]["ratio_after"],
            "ratio_rho_section20": before["ratio"],
            "one_target_cold_ratio_section20": before["cold_ratio_M1"],
            "speedup": x["speedup"]["geomean"], "speedup_ci": [x["speedup"]["lo"], x["speedup"]["hi"]],
            "canonical_step_units": before["canonical_step_units"],
            "ratio_floor": s / floor_s, "verified": bool(x["pins"] and x["control1"]),
        }
        rows.append({
            "regime": "koblitz_batch", "slug": c["slug"], "ec1": ec1_of(c), "legacy": x["curve"],
            "log2_r": x["log2_r"], "cofactor": None,
            "reference": "batch rho, signed Frobenius and negation, the same 32 targets "
                         "(a 32-target batch: a historical diagnostic under AGENTS.md's one-target rule)",
            "reference_s": x["s_rho_s20"], "floor_s": floor_s, "targets": 32,
            "recipes": [recipe], "best": recipe, "source": SOURCES["koblitz_s22"],
            "round": "§22 isolated re-run (2026-09-30); before marks §21 and §22's baseline on main",
            "class": "engineering",
        })
    return rows


def koblitz_one_target_rows(names: Names) -> list[dict]:
    """Ledger §23: one unseen public point per process, 64 per size.  The
    primary comparison under AGENTS.md's one-target rule.  S here is cold,
    the reusable set-up plus the online interval, against one-target rho's
    set-up plus its walk; the online speedup, the canonical-step reading and
    the Bernstein–Lange precomputation model ride beside it, never alone."""
    rows = []
    for x in load(SOURCES["koblitz_s23"])["sizes"]:
        c = names(x["curve"])
        a, n = x["a"], x["n"]
        par = json.loads((REPO / S23_RUNS / f"k{a}n{n}" / "T01.params.json").read_text())
        columns = int(re.search(r"-c(\d+)-", par["name"]).group(1))
        d = x["diagnostics"]
        sh = d["online_phase_shares"]
        online, setup = x["s_ic_online_mean"], x["s_setup"]
        s = setup + online
        floor_s = math.sqrt(math.pi / (4 * n))
        artefact = d["s_curve_construction"] > x["s_rho_online_mean"]
        recipe = {
            "variant": "koblitz_one_target",
            "label": f"aimed m = {par['summands']} collection, m = {par['descent_summands']} descent; one target",
            "m": par["summands"], "descent_summands": par["descent_summands"],
            "base_points": par["factor_base"]["spec"]["points"], "columns": columns,
            "targets": "one public point per process, 64 per size",
            "s": s, "s_online": online, "s_setup": setup,
            "phases_s": {"precomputation": setup,
                         "descent": online * (sh["target_query"] + sh["target_pdp"]
                                              + sh["target_descent"] + sh["target_relation_check"]),
                         "verify": online * sh["recovery_check"]},
            "ratio_rho": x["cold_ratio_ic_over_rho"]["value"],
            "ratio_rho_ci": list(x["cold_ratio_ic_over_rho"]["ci95"]),
            "online_speedup": x["online_speedup"]["mean_ratio"],
            "online_speedup_ci": list(x["online_speedup"]["ci95"]),
            "online_speedup_canonical_step": x["online_speedup_rho_model"]["mean_ratio"],
            "ic_online_over_precomputation_model": x["precomputation_boundary_model"]["ic_online_over_bl"],
            "break_even_targets": x["break_even_targets"]["value"],
            "construction_artefact": artefact,
            "canonical_step_units": x["canonical_step_units"],
            "rho_step_units_as_built": x["rho_units_per_step_median"],
            "cold_ratio_without_curve_construction": d["cold_ratio_without_curve_construction"],
            "ratio_floor": s / floor_s,
            "verified": bool(x["all_rows_pass_the_check"]),
        }
        rows.append({
            "regime": "koblitz", "slug": c["slug"], "ec1": x.get("curve_id") or ec1_of(c),
            "legacy": x["curve"], "log2_r": x["log2_r"], "cofactor": None,
            "reference": "one-target rho on the same point, its set-up plus its walk (§23)",
            "reference_s": x["s_rho_online_mean"] + x["s_rho_setup"], "floor_s": floor_s,
            "targets": 1, "recipes": [recipe], "best": recipe, "source": SOURCES["koblitz_s23"],
            "round": "§23 one unseen point, online and cold (2026-09-30)", "class": "accounting",
        })
    return rows


# ---------------------------------------------------------------------------
# The decomposition oracles, per target.
# ---------------------------------------------------------------------------

ORACLE_LABEL = {
    "enumerate_m2": "enumerate", "enumerate_m3": "enumerate",
    "meet_in_the_middle_m2": "pair table", "meet_in_the_middle_m3": "pair table",
    "semaev_s4_pairs_and_solve": "S4 pairs-and-solve",
    "matrix_f4_splitting": "matrix F4", "cdcl_sat_native_xor": "CDCL SAT",
}


def oracle_rows(names: Names) -> list[dict]:
    rows = []
    for cell in load(SOURCES["oracles"])["oracle_pricing"]["cells"]:
        c = names(cell["curve"]["name"])
        rows.append({
            "slug": c["slug"], "n": cell["n"], "m": cell["m"], "dimension": cell["dimension"],
            "unknowns": cell["unknowns"], "equations": cell["equations"],
            "ffd": [cell["ffd_min"], cell["ffd_max"]], "hit_rate": cell["hit_rate"],
            "log2_r": cell["log2_r"], "disagreements": cell["disagreements"],
            "oracles": {ORACLE_LABEL.get(o["oracle"], o["oracle"]): {
                "gae_per_target": o["gae_per_target_mean"], "found": o["found"],
                "refuted": o.get("refuted", 0), "inconclusive": o.get("inconclusive", 0)}
                for o in cell["oracles"]},
        })
    return rows


# ---------------------------------------------------------------------------
# Exponents.
# ---------------------------------------------------------------------------


def exponents() -> dict:
    fits = load(SOURCES["round5"])["ledger"]["fits"]
    lead = {"prime": "mitm_m2_negfold_walk_balanced", "char2": "mitm_m2_negfold_walk_balanced"}
    out: dict = {}
    for regime, variant in lead.items():
        ph = {f["phase"]: {"alpha": f["alpha"], "r2": f["r_squared"], "points": f["points"]}
              for f in fits if f["regime"] == regime and f["variant"] == variant}
        out[regime] = {"variant": variant, "phases": ph,
                       "variants_total": {f["variant"]: {"alpha": f["alpha"], "r2": f["r_squared"]}
                                          for f in fits if f["regime"] == regime
                                          and f["phase"] == "total" and f["variant"] != "rho_reference"}}
    rho = load(SOURCES["matched_rho_analysis"])["tables"]["rho_exponent_round5_ladder"]
    for regime in ("prime", "char2"):
        out[regime]["rho_matched"] = {"alpha": rho[f"{regime}/rho_matched"]["alpha"],
                                      "r2": rho[f"{regime}/rho_matched"]["r_squared"]}
    s22 = load(SOURCES["koblitz_s22"])["exponent_refit"]
    out["koblitz"] = {"quantity": "ratio to batch rho times sqrt(n), against r, four largest sizes",
                      "beta": s22["after"]["beta"], "ci": [s22["after"]["lo"], s22["after"]["hi"]],
                      "points": s22["after"]["points"], "law": 1 / 6, "model": 0.141}
    return out


def roster(names: Names, measured: set[str]) -> list[dict]:
    out = []
    for c in names.reg["curves"]:
        p = c["params"]
        if c["family"] == "prime":
            field = f"GF(p), p = {p['p']}" if len(p["p"]) < 40 else f"GF(p), {int(p['p']).bit_length()}-bit p"
            coeff = f"a = {p['a']}, b = {p['b']}" if len(p["p"]) < 40 else "standard coefficients"
        else:
            m = p.get("n", p.get("m"))
            field = f"GF(2^{m})"
            a = p.get("a")
            coeff = (f"a = {a}, b = 1" if c["family"] == "koblitz"
                     else f"a = {a}, b = {p.get('b')}")
        out.append({
            "slug": c["slug"], "family": c["family"], "field": field, "coefficients": coeff,
            "order_bits": int(c["order"]).bit_length(), "standard": c["standard_names"],
            "ec1": [r["ec1"] for r in c.get("representations", [])],
            "ec1_unresolved": c.get("ec1_unresolved"), "legacy": c["aliases"][:4],
            "on_board": c["slug"] in measured,
        })
    return out


def families(board: list[dict]) -> dict:
    out = {}
    for code, f in FAMILIES.items():
        rows = [r for r in board if r["family"] == code]
        canon = [r["best"]["canonical_step_units"] for r in rows if "canonical_step_units" in r["best"]]
        built = [r["best"]["rho_step_units_as_built"] for r in rows if "rho_step_units_as_built" in r["best"]]
        out[code] = {**f, "rows": len(rows),
                     "canonical_step_units": [min(canon), max(canon)] if canon else None,
                     "rho_step_units_as_built": [min(built), max(built)] if built else None}
    return out


def build() -> dict:
    names = Names()
    rank_columns = load(SOURCES["n37_rank_columns"])
    if (rank_columns["status"] != "independently_replayed_l0_bounded_diagnostic"
            or rank_columns["selection"]["decision"] != "COUNTED_ENGINEERING_LEAD"
            or rank_columns["selection"]["selected_k"] != 16
            or rank_columns["online_speedup"] is not None
            or rank_columns["fully_priced_cold_speedup"] is not None):
        raise SystemExit("review the n37 bounded result before updating the leaderboard")
    selected = next(c for c in rank_columns["candidates"] if c["folded_columns"] == 16)
    bounded_diagnostic = {
        "curve": names("icv1-f2m37-tm534059-32aad96b")["slug"],
        "source": SOURCES["n37_rank_columns"],
        "class": "counted engineering lead over K42; no admitted IC/rho speedup",
        "candidate_id": selected["candidate_id"],
        "usable_points": selected["usable_points"],
        "folded_columns": selected["folded_columns"],
        "measured_runs": selected["measured_runs"],
        "mean_cold_s_lower_bound": selected["mean_cold_s_lower_bound"],
        "cold_counted_over_rho_diagnostic": selected["cold_counted_over_rho_diagnostic"],
        "cold_counted_over_k42_diagnostic": selected["cold_counted_over_k42_diagnostic"],
        "online_speedup": None,
        "fully_priced_cold_speedup": None,
    }
    k8_decision = load(SOURCES["n37_k8_k16"])
    if (k8_decision["schema"] != "ecbench.k8_k16_callgrind_decision/v1"
            or k8_decision["decision"] != "select_k8_for_larger_field_gate"
            or k8_decision["workloads"] != 16
            or k8_decision["profiles"] != 64
            or not k8_decision["all_archived_and_profiled_scalars_verified"]
            or k8_decision["comparisons"][0]["numerator"] != "ic-k16"
            or k8_decision["comparisons"][0]["denominator"] != "ic-k8"
            or k8_decision["comparisons"][1]["denominator"] != "rho-strong"):
        raise SystemExit("review the untouched K8/K16 instruction result before updating the leaderboard")
    k8_prior = next(c for c in rank_columns["candidates"] if c["folded_columns"] == 8)
    k8_confirmation = {
        "curve": names("icv1-f2m37-tm534059-32aad96b")["slug"],
        "source": SOURCES["n37_k8_k16"],
        "candidate_id": k8_prior["candidate_id"],
        "usable_points": k8_prior["usable_points"],
        "folded_columns": k8_prior["folded_columns"],
        "workloads": k8_decision["workloads"],
        "profiles": k8_decision["profiles"],
        "decision": k8_decision["decision"],
        "k16_over_k8_ir": k8_decision["comparisons"][0]["ratio_of_sums"],
        "k16_over_k8_ir_ci95": k8_decision["comparisons"][0]["bootstrap_95"],
        "k8_over_rho_ir": k8_decision["comparisons"][1]["ratio_of_sums"],
        "k8_over_rho_ir_ci95": k8_decision["comparisons"][1]["bootstrap_95"],
        "online_speedup": None,
        "cold_wall_speedup": None,
    }
    online_ir_decision = load(SOURCES["n37_online_ir"])
    online_k8_claim = load(SOURCES["n37_online_k8_claim"])
    online_k16_claim = load(SOURCES["n37_online_k16_claim"])
    if (online_ir_decision["schema"] != "ecbench.n37_online_ir_decision/v1"
            or online_ir_decision["decision"] != "prioritize_k16_for_isolated_n37_online_wall_gate"
            or online_ir_decision["workloads"] != 16
            or online_ir_decision["profiles"] != 64
            or not online_ir_decision["all_archived_and_profiled_scalars_verified"]
            or online_ir_decision["online"]["aa_max_absolute_relative_deviation"] != 0
            or online_ir_decision["online"]["comparisons"][0]["numerator"] != "ic-k8"
            or online_ir_decision["online"]["comparisons"][0]["denominator"] != "ic-k16"
            or online_ir_decision["online"]["comparisons"][2]["denominator"] != "rho-strong"
            or online_ir_decision["online"]["comparisons"][0]["bootstrap_95"][0] <= 1.10):
        raise SystemExit("review the target-only n37 instruction decision before updating the leaderboard")
    for claim, count, columns in [(online_k8_claim, 592, 8), (online_k16_claim, 1184, 16)]:
        fb = claim["candidate_manifest"]["factor_base"]["inventory"]
        if (not claim["independent_validation"] or not claim["ic_scalar_verified"]
                or not claim["rho_scalar_verified"]
                or claim["isolation_levels"]["ic"] != "L0"
                or claim["isolation_levels"]["rho"] != "L0"
                or fb["usable_point_count"] != count
                or fb["effective_columns"] != columns):
            raise SystemExit("review the exact n37 IC1 claims before updating the leaderboard")
    online_ir_diagnostic = {
        "curve": names("icv1-f2m37-tm534059-32aad96b")["slug"],
        "source": SOURCES["n37_online_ir"],
        "k8_candidate_id": online_k8_claim["candidate_id"],
        "k16_candidate_id": online_k16_claim["candidate_id"],
        "workloads": online_ir_decision["workloads"],
        "profiles": online_ir_decision["profiles"],
        "decision": online_ir_decision["decision"],
        "k8_over_k16_online_ir": online_ir_decision["online"]["comparisons"][0]["ratio_of_sums"],
        "k8_over_k16_online_ir_ci95": online_ir_decision["online"]["comparisons"][0]["bootstrap_95"],
        "k16_over_rho_online_ir": online_ir_decision["online"]["comparisons"][2]["ratio_of_sums"],
        "k16_over_rho_online_ir_ci95": online_ir_decision["online"]["comparisons"][2]["bootstrap_95"],
        "k8_over_k16_cold_ir": online_ir_decision["complete_solve"]["comparisons"][0]["ratio_of_sums"],
        "online_wall_speedup": None,
    }
    ladder, kob, kob1 = ladder_rows(names), koblitz_rows(names), koblitz_one_target_rows(names)
    board = ladder + kob1 + kob
    for r in board:
        r["family"] = FAMILY_OF[r["regime"]]
    measured = {r["slug"] for r in board}
    leaders = {}
    for regime in ("prime", "char2", "koblitz", "koblitz_batch"):
        rows = [r for r in board if r["regime"] == regime]
        ranked = [r for r in rows if not r["best"].get("construction_artefact")] or rows
        best = min(ranked, key=lambda r: r["best"]["ratio_rho"])
        top = max(rows, key=lambda r: r["log2_r"])
        leaders[regime] = {"slug": best["slug"], "log2_r": best["log2_r"],
                           "ratio_rho": best["best"]["ratio_rho"], "s": best["best"]["s"],
                           "recipe": best["best"]["label"],
                           "one_target": best["best"].get("one_target_cold_ratio_section20"),
                           "online_speedup": best["best"].get("online_speedup"),
                           "largest": {"slug": top["slug"], "log2_r": top["log2_r"],
                                       "ratio_rho": top["best"]["ratio_rho"]}}
    return {
        "schema_version": 1,
        "generated_by": "scripts/build_ic_leaderboard.py",
        "unit": "none: the board has three families, each with its own unit, reference and floor "
                "(families[*].unit, .reference, .floor); S, ratio to reference and ratio to floor "
                "are comparable within a family and nowhere else",
        "what_this_is": "A view of frozen whole-pipeline ECDLP measurements, curve by curve, "
                        "each against its own family's matched reference and generic floor. "
                        "No new measurement.",
        "what_this_is_not": ["a measurement", "a speedup", "a claim about any curve not on it"],
        "class": "accounting",
        "sources": {k: {"path": v, "sha256": sha256(v)} for k, v in SOURCES.items()},
        "families": families(board),
        "leaders": leaders, "board": board, "oracles": oracle_rows(names),
        "bounded_diagnostics": [bounded_diagnostic],
        "n37_k8_k16_confirmation": k8_confirmation,
        "n37_online_ir_diagnostic": online_ir_diagnostic,
        "exponents": exponents(), "roster": roster(names, measured),
        "phases": [{"id": i, "name": n, "what": w} for i, n, w in PHASES],
    }


# ---------------------------------------------------------------------------
# Rendering.
# ---------------------------------------------------------------------------


def g3(v: float | None) -> str:
    if v is None:
        return "—"
    if v == 0:
        return "0"
    if abs(v) >= 1000:
        return f"{v:,.0f}"
    return f"{v:.3g}"


def times(v: float) -> str:
    return f"{g3(v)}×"


def esc(s) -> str:
    return html.escape(str(s))


def rich(s) -> str:
    """Escaped text with `code` spans, for the family definitions shared with the Markdown."""
    return re.sub(r"`([^`]+)`", r"<code>\1</code>", esc(s))


def markdown(doc: dict) -> str:
    fams = doc["families"]
    L = ["# Index-calculus leaderboard", "",
         "Generated by `scripts/build_ic_leaderboard.py` from frozen evidence; the page is "
         "[`docs/ic-leaderboard.html`](../ic-leaderboard.html).  Curves are named by ICV1 slug "
         "([`docs/curves/ICV1.md`](../curves/ICV1.md)).  No row is below its reference.", "",
         "**There is no single unit on this page.**  The board is three tables, A, B and C.  Each is "
         "priced in its own unit, against its own reference and floor, stated at its head.  `S`, "
         "`× reference` and `× floor` are comparable within a table and nowhere else.  The same slug "
         "can appear in two tables with different `S`, because the unit, the reference and the number "
         "of targets differ.", ""]
    for code in ("A", "B", "C"):
        f = fams[code]
        rows = [r for r in doc["board"] if r["family"] == code]
        multi = len(f["regimes"]) > 1
        L += [f"## Table {code}: {f['title']} ({f['section']})", "",
              f"- **Unit:** {f['unit']}.  {f['unit_how']}",
              f"- **Reference:** {f['reference']}",
              f"- **Floor:** {f['floor']}",
              f"- **Targets:** {f['targets']}",
              f"- **Status:** {f['status']}"]
        steps = step_line(f)
        if steps:
            L.append(f"- **Rho step:** {steps}")
        head = (["regime"] if multi else []) + ["curve", "log₂ r", "recipe", "S, IC", "S, reference",
                                                "× reference", "S, floor", "× floor", "ok"]
        align = (([":--"] if multi else []) + [":--", "--:", ":--", "--:", "--:", "--:", "--:", "--:", ":--"])
        L += ["", "| " + " | ".join(head) + " |", "|" + "|".join(align) + "|"]
        for r in rows:
            b = r["best"]
            cells = ([REGIME_NAME[r["regime"]]] if multi else []) + [
                f"`{r['slug']}`", f"{r['log2_r']:.1f}", b["label"], g3(b["s"]), g3(r["reference_s"]),
                f"**{times(b['ratio_rho'])}**", g3(r["floor_s"]), times(b["ratio_floor"]),
                "✓" if b["verified"] else "✗"]
            L.append("| " + " | ".join(cells) + " |")
        L.append("")
    L += ["## Phases, as a share of S", "",
          "Shares are fractions of each row's own `S`, so they carry no unit; read each row in the "
          "table named in its first column.", "",
          "| table | curve | " + " | ".join(n for _, n, _ in PHASES) + " |",
          "|:--|:--|" + "--:|" * len(PHASES)]
    for code in ("A", "B", "C"):
        for r in (r for r in doc["board"] if r["family"] == code):
            ph, s = r["best"]["phases_s"], r["best"]["s"]
            L.append(f"| {code} | `{r['slug']}` | " + " | ".join(
                (f"{100 * ph[i] / s:.1f}%" if i in ph else "—") for i, _, _ in PHASES) + " |")
    d = doc["bounded_diagnostics"][0]
    L += ["", "## Bounded n37 diagnostic outside tables A–C", "",
          f"The separately calibrated `{d['curve']}` shared-rank K16 candidate "
          f"(`{d['candidate_id']}`) has {d['usable_points']:,} usable points and "
          f"{d['folded_columns']} folded columns. Across {d['measured_runs']} verified "
          f"one-target runs, its mean cold counted `S` lower bound is "
          f"{d['mean_cold_s_lower_bound']:.3f}; its counted IC/rho quotient is "
          f"{d['cold_counted_over_rho_diagnostic']:.3f} and K16/K42 is "
          f"{d['cold_counted_over_k42_diagnostic']:.3f}. Native work is unpriced "
          "for both arms, and L0 timing cannot establish an online speedup. "
          "This row is intentionally outside the three fully priced unit families; "
          f"read the [frozen decision](../../{d['source']}).", ""]
    c = doc["n37_k8_k16_confirmation"]
    L += [f"The untouched {c['workloads']}-target K8/K16 confirmation selected "
          f"`{c['candidate_id']}` for the cold implementation route: K16/K8 "
          f"whole-solve Callgrind Ir is {c['k16_over_k8_ir']:.3f} "
          f"[{c['k16_over_k8_ir_ci95'][0]:.3f}, {c['k16_over_k8_ir_ci95'][1]:.3f}], "
          f"and K8/rho is {c['k8_over_rho_ir']:.3f} "
          f"[{c['k8_over_rho_ir_ci95'][0]:.3f}, {c['k8_over_rho_ir_ci95'][1]:.3f}]. "
          "All 64 profiles and 320 independent measured replays verified. "
          "The unit is simulated whole-solve instructions, not isolated online "
          "wall time; K16 remains a target-only candidate. "
          f"Read the [raw instruction decision](../../{c['source']}).", ""]
    q = doc["n37_online_ir_diagnostic"]
    L += [f"The fresh {q['workloads']}-target panel reverses the base choice inside "
          "the **target-only simulated-instruction** interval: "
          f"`{q['k16_candidate_id']}` is prioritized for an isolated n37 online "
          f"wall test because K8/K16 target-only Callgrind Ir is "
          f"{q['k8_over_k16_online_ir']:.3f} "
          f"[{q['k8_over_k16_online_ir_ci95'][0]:.3f}, "
          f"{q['k8_over_k16_online_ir_ci95'][1]:.3f}]. "
          f"K16/rho is {q['k16_over_rho_online_ir']:.3f} "
          f"[{q['k16_over_rho_online_ir_ci95'][0]:.3f}, "
          f"{q['k16_over_rho_online_ir_ci95'][1]:.3f}] in the same instruction "
          f"unit, while K8/K16 complete-solve Ir is {q['k8_over_k16_cold_ir']:.3f}. "
          "All 64 profile scalars and 320 independent measured replays verified. "
          "Mac L0 timing leaves the primary online wall speedup unknown; both "
          "bases remain live at n41/n53. "
          f"Read the [raw target-only decision](../../{q['source']}).", ""]
    L += ["", "## Sources", ""]
    for k, v in doc["sources"].items():
        L.append(f"- `{v['path']}` — sha256 `{v['sha256'][:16]}…`")
    return "\n".join(L) + "\n"


STYLE = """
/* Layout: one reading column; the board and the step table are the page, everything else hangs off them. */
:root {
  --ground: #e7ebef; --surface: #fbfcfd; --sunk: #f0f3f6; --ink: #0e1922; --ink-2: #354855;
  --muted: #576976; --rule: #ccd5dc; --data: #1f6fd0; --data-soft: #d7e4f6; --bound: #b3372a;
  --bound-soft: #f3dcd8;
  --ph-precomputation: #b9c3cb; --ph-setup: #9aa7b1; --ph-constructions: #c4a35a; --ph-factor_base: #2f8f83; --ph-table: #6fb7a8;
  --ph-relations: #1f6fd0; --ph-linear_algebra: #8b5bc4; --ph-descent: #d0782f; --ph-verify: #4d5b66;
  --serif: "Newsreader", "Iowan Old Style", Georgia, serif;
  --sans: "IBM Plex Sans", system-ui, -apple-system, "Segoe UI", sans-serif;
  --mono: "IBM Plex Mono", ui-monospace, "SF Mono", Menlo, monospace;
}
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --ground: #0d1318; --surface: #151d24; --sunk: #101820; --ink: #e3eaef; --ink-2: #b6c5cf;
  --muted: #8799a6; --rule: #25323a; --data: #4e92de; --data-soft: #1c2e42; --bound: #d9604e;
  --bound-soft: #38221f;
  --ph-precomputation: #3d4a54; --ph-setup: #6c7a85; --ph-constructions: #b8954a; --ph-factor_base: #3aa698; --ph-table: #7cc4b5;
  --ph-relations: #4e92de; --ph-linear_algebra: #a57ad8; --ph-descent: #e08c45; --ph-verify: #93a3ae;
  color-scheme: dark; } }
:root[data-theme="dark"] {
  --ground: #0d1318; --surface: #151d24; --sunk: #101820; --ink: #e3eaef; --ink-2: #b6c5cf;
  --muted: #8799a6; --rule: #25323a; --data: #4e92de; --data-soft: #1c2e42; --bound: #d9604e;
  --bound-soft: #38221f;
  --ph-precomputation: #3d4a54; --ph-setup: #6c7a85; --ph-constructions: #b8954a; --ph-factor_base: #3aa698; --ph-table: #7cc4b5;
  --ph-relations: #4e92de; --ph-linear_algebra: #a57ad8; --ph-descent: #e08c45; --ph-verify: #93a3ae;
  color-scheme: dark; }
* { box-sizing: border-box; }
body { margin: 0; padding-inline: 16px; padding-block: 36px 64px; background: var(--ground); color: var(--ink);
  font-family: var(--sans); font-size: 15px; line-height: 1.55; }
.wrap { max-width: 1080px; margin: 0 auto; display: flex; flex-direction: column; gap: 28px; }
header { display: flex; flex-direction: column; gap: 12px; }
.eyebrow { font-family: var(--mono); font-size: 11px; letter-spacing: .13em; text-transform: uppercase; color: var(--muted); }
h1 { font-family: var(--serif); font-weight: 500; font-size: clamp(30px, 5.2vw, 44px); line-height: 1.08; margin: 0; text-wrap: balance; }
h2 { font-family: var(--serif); font-weight: 500; font-size: 22px; margin: 0; text-wrap: balance; }
.verdict { font-family: var(--serif); font-size: clamp(17px, 2.2vw, 20px); line-height: 1.45; color: var(--ink-2); max-width: 64ch; margin: 0; }
.facts { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 1px; background: var(--rule);
  border: 1px solid var(--rule); border-radius: 3px; overflow: hidden; }
.fact { background: var(--surface); padding: 14px 16px; display: flex; flex-direction: column; gap: 3px; min-width: 0; }
.fact dt { font-family: var(--mono); font-size: 10.5px; letter-spacing: .1em; text-transform: uppercase; color: var(--muted); }
.fact dd { margin: 0; font-family: var(--mono); font-size: 20px; font-weight: 500; font-variant-numeric: tabular-nums; }
.fact p { margin: 0; font-size: 12.5px; color: var(--muted); line-height: 1.4; overflow-wrap: anywhere; }
.panel { background: var(--surface); border: 1px solid var(--rule); border-radius: 4px; display: flex; flex-direction: column; gap: 12px; padding-block: 18px; }
.panel > * { padding-inline: 20px; }
.panel p { margin: 0; max-width: 72ch; color: var(--ink-2); font-size: 14px; }
.scroll { overflow-x: auto; padding-inline: 0; }
table { border-collapse: collapse; width: 100%; font-size: 13px; font-variant-numeric: tabular-nums; }
th, td { padding: 7px 10px; text-align: left; border-bottom: 1px solid var(--rule); vertical-align: top; }
th { font-family: var(--mono); font-size: 10.5px; font-weight: 500; letter-spacing: .06em; text-transform: uppercase; color: var(--muted); white-space: nowrap; }
td.n, th.n { text-align: right; white-space: nowrap; }
td.recipe { min-width: 210px; }
tr.lead td { background: var(--data-soft); }
code, .slug { font-family: var(--mono); font-size: 12px; }
.slug { white-space: nowrap; }
.chip { display: inline-block; font-family: var(--mono); font-size: 10px; letter-spacing: .08em; text-transform: uppercase;
  padding: 1px 6px; border: 1px solid var(--rule); border-radius: 2px; color: var(--ink-2); white-space: nowrap; }
.ratio { font-family: var(--mono); font-weight: 600; color: var(--bound); }
.bar { display: flex; height: 12px; min-width: 160px; border-radius: 2px; overflow: hidden; background: var(--sunk); }
.bar span { display: block; height: 100%; }
.legend { display: flex; flex-wrap: wrap; gap: 6px 16px; font-size: 12px; color: var(--ink-2); }
.legend i { display: inline-block; width: 10px; height: 10px; border-radius: 2px; margin-right: 6px; vertical-align: -1px; }
.muted { color: var(--muted); }
.defs { margin-block: 0; display: flex; flex-direction: column; gap: 6px; }
.defs div { display: grid; grid-template-columns: 84px 1fr; gap: 12px; }
.defs dt { font-family: var(--mono); font-size: 10.5px; letter-spacing: .08em; text-transform: uppercase; color: var(--muted); padding-top: 2px; }
.defs dd { margin: 0; font-size: 13px; color: var(--ink-2); max-width: 90ch; }
@media (max-width: 560px) { .defs div { grid-template-columns: 1fr; gap: 1px; } }
svg text { fill: var(--ink-2); font-family: var(--mono); font-size: 10px; }
svg .axis { stroke: var(--rule); }
svg .rho { stroke: var(--bound); stroke-width: 1.5; }
details summary { cursor: pointer; font-family: var(--mono); font-size: 12px; color: var(--data); }
details summary:focus-visible, a:focus-visible { outline: 2px solid var(--data); outline-offset: 2px; }
a { color: var(--data); }
footer { font-size: 12.5px; color: var(--muted); display: flex; flex-direction: column; gap: 6px; }
.site-nav { display: flex; flex-wrap: wrap; gap: 4px 18px; padding-block: 2px 10px; border-bottom: 1px solid var(--rule);
  font-family: var(--mono); font-size: 12px; letter-spacing: .02em; }
.site-nav a { color: var(--ink-2); text-decoration: none; padding-block: 4px; }
.site-nav a:hover { color: var(--data); }
.site-nav a[aria-current="page"] { color: var(--ink); border-bottom: 2px solid var(--data); }
footer p { margin: 0; max-width: 90ch; }
"""


def phase_bar(ph: dict, s: float) -> str:
    segs = []
    for pid, name, _ in PHASES:
        v = ph.get(pid)
        if not v:
            continue
        pct = 100 * v / s
        segs.append(f'<span style="width:{pct:.2f}%;background:var(--ph-{pid})" '
                    f'title="{esc(name)}: {pct:.1f}%"></span>')
    return f'<div class="bar" role="img" aria-label="phase shares">{"".join(segs)}</div>'


def ratio_chart(board: list[dict]) -> str:
    """Each curve's best recipe on a log axis of ratio to its reference."""
    W, H, left, right = 980, 60 + 26 * 3, 150, 20
    lo, hi = 0.0, 3.0   # log10 of 1x .. 1000x
    x = lambda r: left + (math.log10(r) - lo) / (hi - lo) * (W - left - right)  # noqa: E731
    out = [f'<svg viewBox="0 0 {W} {H}" width="100%" role="img" '
           f'aria-label="Ratio to reference per curve, log scale" style="min-width:640px">']
    for t in (1, 3, 10, 30, 100, 300, 1000):
        out.append(f'<line class="axis" x1="{x(t):.1f}" y1="18" x2="{x(t):.1f}" y2="{H - 18}"/>'
                   f'<text x="{x(t):.1f}" y="{H - 4}" text-anchor="middle">{t}×</text>')
    out.append(f'<line class="rho" x1="{x(1):.1f}" y1="14" x2="{x(1):.1f}" y2="{H - 18}"/>'
               f'<text x="{x(1) + 4:.1f}" y="12">reference = 1×</text>')
    for i, regime in enumerate(("koblitz", "prime", "char2")):
        y = 36 + i * 26
        out.append(f'<text x="0" y="{y + 4}">{esc(REGIME_NAME[regime])} ({FAMILY_OF[regime]})</text>')
        for r in (r for r in board if r["regime"] == regime):
            v = r["best"]["ratio_rho"]
            out.append(f'<circle cx="{x(v):.1f}" cy="{y}" r="5" fill="var(--data)" fill-opacity=".55" '
                       f'stroke="var(--data)"><title>{esc(r["slug"])}: {times(v)} at 2^{r["log2_r"]:.1f}</title></circle>')
    out.append("</svg>")
    return "".join(out)


def page(doc: dict, standalone: bool) -> str:
    board = doc["board"]
    L = doc["leaders"]
    by_regime = {g: sorted((r for r in board if r["regime"] == g), key=lambda r: r["best"]["ratio_rho"])
                 for g in ("koblitz", "prime", "char2", "koblitz_batch")}
    head = ('<title>Index Calculus Leaderboard</title>\n'
            '<link rel="preconnect" href="https://fonts.googleapis.com">\n'
            '<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>\n'
            '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,400;6..72,500'
            '&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500;600&display=swap">\n'
            f"<style>{STYLE}</style>\n")
    P = []
    K1 = L["koblitz"]
    kob_rows = by_regime["koblitz"]
    online_lo = min(r["best"]["online_speedup"] for r in kob_rows)
    online_hi = max(r["best"]["online_speedup"] for r in kob_rows)
    cold = [r["best"]["ratio_rho"] for r in kob_rows if not r["best"]["construction_artefact"]]
    model = [r["best"]["ic_online_over_precomputation_model"] for r in kob_rows]
    site = "https://aburan28.github.io/crypto/"
    nav = [("Overview", site), ("Scoreboard", site + "scoreboard/"),
           ("Leaderboard", None), ("Where things stand", site + "scoreboard/ic-current-state.html"),
           ("Algorithm lab", site + "scoreboard/algorithm-lab.html"), ("Campaign status", site + "status/"),
           ("Repository", "https://github.com/aburan28/crypto")]
    P.append('<div class="wrap"><nav class="site-nav" aria-label="Site">' + "".join(
        f'<a href="#board" aria-current="page">{esc(t)}</a>' if u is None else f'<a href="{esc(u)}">{esc(t)}</a>'
        for t, u in nav) + '</nav>')
    P.append('<header>'
             '<span class="eyebrow">ECDLP · index calculus · whole pipeline · accounting view, no new measurement</span>'
             '<h1>Index Calculus Leaderboard</h1>'
             '<p class="verdict">Every curve this repository has priced end to end, with the best index-calculus '
             'recipe measured on it, against Pollard rho on the same curve and the same point. <strong>The units '
             'differ between the three tables below, and so do the references and the floors: S and the ratios '
             'are comparable inside a table and nowhere else.</strong> <strong>Cold, set-up included, no row is '
             'below its reference.</strong> On one unseen point of a Koblitz curve (table A) '
             f'the index calculus\'s online interval is {g3(online_lo)}–{g3(online_hi)}× faster than rho once '
             'its reusable set-up exists, but the set-up is paid first: cold it is '
             f'{g3(min(cold))}–{g3(max(cold))}× slower, and a generic walk given the same precomputation would be '
             f'{g3(min(model))}–{g3(max(model))}× faster online again (a model). Prime and random-binary curves (table B), '
             f'one target each: {times(L["prime"]["ratio_rho"])} to {times(L["prime"]["largest"]["ratio_rho"])} and '
             f'{times(L["char2"]["ratio_rho"])} to {times(L["char2"]["largest"]["ratio_rho"])}.</p>'
             '<dl class="facts">')
    def unit_ref(code: str) -> str:
        f = doc["families"][code]
        return f'<br>× {esc(f["ref_short"])}, in {esc(f["unit_short"])} (table {code})'

    P.append(f'<div class="fact"><dt>Koblitz · one point, cold</dt><dd>{times(K1["ratio_rho"])}</dd>'
             f'<p>at log₂ r = {K1["log2_r"]:.1f}, <span class="slug">{esc(K1["slug"])}</span>; online alone '
             f'{times(K1["online_speedup"])} faster than rho; largest size {times(K1["largest"]["ratio_rho"])} cold'
             f'{unit_ref("A")}</p></div>')
    for g in ("prime", "char2"):
        ld = L[g]
        big = ld["largest"]
        P.append(f'<div class="fact"><dt>Closest · {esc(REGIME_NAME[g])}</dt><dd>{times(ld["ratio_rho"])}</dd>'
                 f'<p>at log₂ r = {ld["log2_r"]:.1f}, <span class="slug">{esc(ld["slug"])}</span>; '
                 f'at the largest size, log₂ r = {big["log2_r"]:.1f}: {times(big["ratio_rho"])}'
                 f'{unit_ref("B")}</p></div>')
    P.append(f'<div class="fact"><dt>Curves</dt><dd>{len({r["slug"] for r in board})} priced</dd>'
             f'<p>of {len(doc["roster"])} named in the repository; every one by its ICV1 slug</p></div></dl></header>')
    d = doc["bounded_diagnostics"][0]
    P.append(f'<section class="card" id="bounded-n37-diagnostic"><h2>New n37 counted diagnostic, outside the priced tables</h2>'
             f'<p>On <code>{esc(d["curve"])}</code>, shared-rank K16 has {d["usable_points"]:,} actual usable '
             f'points and {d["folded_columns"]} folded columns. In {d["measured_runs"]} verified '
             f'one-target runs, mean cold counted S is {d["mean_cold_s_lower_bound"]:.3f} as a lower '
             f'bound; its counted IC/rho quotient is {d["cold_counted_over_rho_diagnostic"]:.3f}, '
             f'and K16/K42 is {d["cold_counted_over_k42_diagnostic"]:.3f}. Both arms leave native '
             'work unpriced, and L0 timing gives no admitted online speedup. These figures are not '
             'comparable to tables A–C. '
             f'<a href="https://github.com/aburan28/crypto/blob/main/{esc(d["source"])}">Frozen decision</a>. '
             'A separate whole-solve Callgrind census finds K16/K42 = 0.367 and K16/rho = 4.461 '
             'in simulated instructions on eight paired points; native wall and online speed remain unknown. '
             '<a href="https://github.com/aburan28/crypto/blob/main/research/ecbench_callgrind_solve_20261004/RESULT.md">Raw instruction replay</a>.</p></section>')
    c = doc["n37_k8_k16_confirmation"]
    P.append(f'<section class="card" id="bounded-n37-k8-confirmation"><h2>Untouched n37 K8/K16 instruction confirmation</h2>'
             f'<p>On <code>{esc(c["curve"])}</code>, K8 has {c["usable_points"]:,} usable points '
             f'and {c["folded_columns"]} folded columns. Across {c["workloads"]} new same-point '
             f'public workloads, K16/K8 whole-solve Callgrind Ir is {c["k16_over_k8_ir"]:.3f} '
             f'[{c["k16_over_k8_ir_ci95"][0]:.3f}, {c["k16_over_k8_ir_ci95"][1]:.3f}], '
             f'selecting K8 for the cold implementation route. K8/rho is '
             f'{c["k8_over_rho_ir"]:.3f} [{c["k8_over_rho_ir_ci95"][0]:.3f}, '
             f'{c["k8_over_rho_ir_ci95"][1]:.3f}] in the same instruction unit. '
             'All 64 profiles and 320 independently replayed measured records verified. '
             'This is not isolated online wall speed; K16 remains a target-only candidate. '
             f'<a href="https://github.com/aburan28/crypto/blob/main/{esc(c["source"])}">Frozen decision</a>.</p></section>')
    q = doc["n37_online_ir_diagnostic"]
    P.append(f'<section class="card" id="bounded-n37-online-ir"><h2>Fresh n37 K8/K16 target-only instruction gate</h2>'
             f'<p>On <code>{esc(q["curve"])}</code>, K8/K16 target-only Callgrind Ir is '
             f'{q["k8_over_k16_online_ir"]:.3f} '
             f'[{q["k8_over_k16_online_ir_ci95"][0]:.3f}, '
             f'{q["k8_over_k16_online_ir_ci95"][1]:.3f}], so the frozen gate '
             'prioritizes K16 for an isolated n37 online wall test. K16/rho '
             f'target-only Ir is {q["k16_over_rho_online_ir"]:.3f} '
             f'[{q["k16_over_rho_online_ir_ci95"][0]:.3f}, '
             f'{q["k16_over_rho_online_ir_ci95"][1]:.3f}], but K8/K16 '
             f'complete-solve Ir is {q["k8_over_k16_cold_ir"]:.3f}. '
             'All 64 profiles and 320 independent native replays verified. '
             'Callgrind Ir is simulated instruction work, not isolated wall '
             'speed; n41/n53 and ECC2K-130 transfer stay open. '
             f'<a href="https://github.com/aburan28/crypto/blob/main/{esc(q["source"])}">Frozen decision and raw replay</a>.</p></section>')

    head_cells = ('<th class="n">#</th><th>curve</th><th class="n">log₂ r</th><th>recipe</th>'
                  '<th class="n">m</th><th class="n">|F|</th><th class="n">K</th><th class="n">S, IC</th>'
                  '<th class="n">S, ref</th><th class="n">× ref</th><th class="n">S, floor</th>'
                  '<th class="n">× floor</th><th>phases</th><th>ok</th>')
    ncols = 14

    def board_rows(g: str) -> None:
        for i, r in enumerate(by_regime[g], 1):
            b = r["best"]
            ci = (f'<br><span class="muted">[{g3(b["ratio_rho_ci"][0])}, {g3(b["ratio_rho_ci"][1])}]</span>'
                  if "ratio_rho_ci" in b else "")
            if "online_speedup" in b:
                ci += (f'<br><span class="muted">online {g3(b["online_speedup"])}× faster '
                       f'[{g3(b["online_speedup_ci"][0])}, {g3(b["online_speedup_ci"][1])}]; '
                       f'{g3(b["online_speedup_canonical_step"])}× at the canonical step</span>'
                       f'<br><span class="muted">model: {g3(b["ic_online_over_precomputation_model"])}× slower online '
                       f'than a generic walk with the same set-up; break-even {g3(b["break_even_targets"])} targets</span>')
                if b["construction_artefact"]:
                    ci += (f'<br><span class="muted">curve construction exceeds rho\'s walk here (§23.9); '
                           f'without it {times(b["cold_ratio_without_curve_construction"])}</span>')
            if "ratio_rho_section21" in b:
                ci += (f'<br><span class="muted">was {times(b["ratio_rho_section21"])} (§21)</span>'
                       f'<br><span class="muted">derived one-target cold estimate: {times(b["one_target_cold_ratio_section20"])} '
                     f'(§20; not table A\'s measurement)</span>')
            rank = "—" if b.get("construction_artefact") else i
            P.append(f'<tr class="{"lead" if r["slug"] == L[g]["slug"] else ""}"><td class="n">{rank}</td>'
                     f'<td><span class="slug">{esc(r["slug"])}</span></td><td class="n">{r["log2_r"]:.1f}</td>'
                     f'<td class="recipe">{esc(b["label"])}</td><td class="n">{b["m"]}</td><td class="n">{b["base_points"]:,}</td>'
                     f'<td class="n">{b["columns"]:,}</td><td class="n">{g3(b["s"])}</td>'
                     f'<td class="n">{g3(r["reference_s"])}</td><td class="n"><span class="ratio">{times(b["ratio_rho"])}</span>{ci}</td>'
                     f'<td class="n">{g3(r["floor_s"])}</td><td class="n">{times(b["ratio_floor"])}</td><td>{phase_bar(b["phases_s"], b["s"])}</td>'
                     f'<td>{"✓" if b["verified"] else "✗"}</td></tr>')


    def family_defs(code: str) -> str:
        f = doc["families"][code]
        defs = [("Unit", f"{f['unit']}. {f['unit_how']}"), ("Reference", f["reference"]),
                ("Floor", f["floor"]), ("Targets", f["targets"]), ("Status", f["status"])]
        steps = step_line(f)
        if steps:
            defs.append(("Rho step", steps))
        return ('<dl class="defs">' + "".join(
            f'<div><dt>{esc(k)}</dt><dd>{rich(v)}</dd></div>' for k, v in defs) + '</dl>'
            f'<p class="muted">Read S, × ref and × floor in table {code} only.</p>')

    def family_table(code: str) -> None:
        f = doc["families"][code]
        multi = len(f["regimes"]) > 1
        P.append(f'<div class="scroll"><table><thead><tr>{head_cells}</tr></thead><tbody>')
        for g in f["regimes"]:
            if multi:
                P.append(f'<tr><td colspan="{ncols}"><span class="chip">{esc(REGIME_NAME[g])}</span></td></tr>')
            board_rows(g)
        P.append('</tbody></table></div>')

    # The board: the index, the chart, and what the three tables are.
    P.append('<section class="panel" id="board"><h2>The board</h2>'
             '<p>One row per curve: its best recipe, ranked within its table by the cold ratio to that table\'s '
             'reference, the whole method from set-up to a verified logarithm. <strong>There is no single unit '
             'here.</strong> Tables A, B and C below are priced in different units against different references '
             'and floors, stated at the head of each; S, × ref and × floor are comparable inside a table and '
             'nowhere else. Tables A and B solve one target per row, as AGENTS.md\'s one-target rule requires: '
             'A on one unseen public point per process against one-target rho on the same point (§23, 64 points '
             'per size), B against the matched negation-map rho (ledger §18). Table C is a historical 32-target '
             'diagnostic. For the Koblitz rows of A the online interval, the reading with rho at the canonical '
             'step and the generic-precomputation model sit beside the cold ratio, as §23 requires. The bar '
             'splits S into phases. A walked-target row on the ladder is resolved to about a factor 1.8 either '
             'way at three repetitions (ledger §13.7).</p>'
             f'<div class="scroll">{ratio_chart(board)}</div>'
             '<p class="muted">Each dot is a ratio to its own table\'s reference, in that table\'s unit; the '
             'three rows of dots are not comparable with one another.</p>'
             '<div class="legend">' + "".join(
                 f'<span><i style="background:var(--ph-{pid})"></i>{esc(n)}</span>' for pid, n, _ in PHASES)
             + '</div></section>')
    for code, anchor, title in (("A", "board-A", "Table A · Koblitz, one target"),
                                ("B", "board-B", "Table B · prime and random-binary ladder")):
        f = doc["families"][code]
        P.append(f'<section class="panel" id="{anchor}"><h2>{esc(title)} <span class="chip">{esc(f["section"])}</span></h2>'
                 + family_defs(code))
        family_table(code)
        P.append('</section>')

    # Table C: the Koblitz batch, a historical diagnostic.
    fc = doc["families"]["C"]
    P.append(f'<section class="panel" id="batch"><h2>Table C · Koblitz, 32 targets in a batch '
             f'<span class="chip">{esc(fc["section"])}</span></h2>'
             '<p>The collection thread\'s rounds before §23 solved 32 targets together and priced them per target '
             'against batch rho on the same 32 (§19–§22). Under the one-target rule that is a historical '
             'diagnostic, not a primary comparison; it is kept here with its §21 figure as the before mark and '
             'the one-target cold ratio §20 derived from this run (shared phases plus one mean descent, against '
             'single-target rho), which table A\'s measured figure replaces. The phase split is the most detailed '
             'the thread has.</p>'
             + family_defs("C"))
    family_table("C")
    P.append('</section>')

    # Every recipe, every curve (the ladder).
    ladder = [r for r in board if r["regime"] in ("prime", "char2")]
    variants = []
    for r in ladder:
        for v in r["recipes"]:
            if v["variant"] not in variants:
                variants.append(v["variant"])
    P.append('<section class="panel" id="recipes"><h2>Every recipe on every ladder curve</h2>'
             '<p>Table B only: ratio to the matched rho for each variant the ladder ran, three repetitions averaged. '
             'The minimum per curve is the row on the board. The balanced recipe sizes the base at the '
             'family optimum F = (c·#E·t/k)<sup>1/3</sup>, which is 0.75·r<sup>1/6</sup> in S on a '
             'prime-order curve (ledger §11.2): a rising cost against a flat reference, so this family '
             'has no crossover to find at any size.</p><div class="scroll"><table><thead><tr><th>recipe</th>')
    for r in ladder:
        P.append(f'<th class="n" title="{esc(r["slug"])}">{"p" if r["regime"] == "prime" else "b"}·{r["log2_r"]:.1f}</th>')
    P.append('</tr></thead><tbody>')
    for vname in variants:
        P.append(f'<tr><td>{esc(VARIANTS.get(vname, vname))}</td>')
        for r in ladder:
            v = next((x for x in r["recipes"] if x["variant"] == vname), None)
            best = v is r["best"]
            cell = times(v["ratio_rho"]) if v else "—"
            P.append(f'<td class="n">{"<strong>" + cell + "</strong>" if best else cell}</td>')
        P.append('</tr>')
    P.append('</tbody></table></div><p class="muted">Column heads: p = prime field, b = random binary, '
             'then log₂ r. Hover a head for the curve\'s slug.</p></section>')

    # Step by step.
    ex = doc["exponents"]
    head_kob = next(r for r in board if r["regime"] == "koblitz_batch"
                    and r["slug"] == L["koblitz_batch"]["slug"])["best"]
    head_pr = next(r for r in board if r["slug"] == L["prime"]["slug"])["best"]
    rows_steps = [
        ("setup", "#E by Schoof / Koblitz recurrence, factor, clear the cofactor", "one-off; Koblitz: 25–171 units a base point at 2^27.8 and up", "—", "—"),
        ("factor_base", "Base of |F| points: small abscissae (prime), a subspace (binary), Frobenius orbits (Koblitz)", "|F|; the balanced optimum puts |F| ∝ r^{1/3}", "prime", "factor_base"),
        ("constructions", "Orbit map, point index and decomposition classes around the work (Koblitz); built once per base since §21", "per base point; 0.4–16% of S after §21 (ledger §21.3)", "—", "—"),
        ("table", "All |F|(|F|+1)/2 pair sums, folded by negation or by ⟨σ, −1⟩; on the ladder counted inside the factor base", "|F|² / (4t) stored pairs, t the fold", "—", "—"),
        ("relations", "Draw a target, ask the PDP oracle whether it is a sum of m base points, keep verified hits", "trials ≈ (K+1) / ceiling; per trial one table probe (m = 2) or |F| (m = 3)", "prime", "relations"),
        ("linear_algebra", "Relations as rows over Z/rZ; dense Gauss–Jordan (ladder) or filter + block Wiedemann (Koblitz)", "dense K³, sparse ~K·weight", "prime", "linear_algebra"),
        ("descent", "Each target decomposed over the logged base: an m = 2 or 3 whole-base scan", "≈ 2r / |F|² probes per target; does not amortise over targets", "—", "—"),
        ("verify", "Re-add every relation; check [d]G = Q for every target", "one scalar multiplication per target", "—", "—"),
    ]
    P.append('<section class="panel" id="steps"><h2>The ECDLP, step by step</h2>'
             '<p>What each phase computes, what it costs, and how fast it grows. The exponent α is the fit of that '
             'phase\'s operations to r<sup>α</sup> on the Round-5 ladder (table B, in GAE) for the balanced recipe; rho grows as '
             'r<sup>1/2</sup>, so a phase with α above ½ loses ground at every size. The share columns are the '
             'phase\'s part of S on each family\'s leading row.</p><div class="scroll"><table><thead><tr>'
             '<th>phase</th><th>what it does</th><th>cost model</th><th class="n">α prime</th><th class="n">α binary</th>'
             f'<th class="n">share · Koblitz batch 2^{L["koblitz_batch"]["log2_r"]:.0f} (§22)</th><th class="n">share · prime 2^{L["prime"]["log2_r"]:.1f}</th></tr></thead><tbody>')
    for pid, what, model, _, fitkey in rows_steps:
        name = next(n for i, n, _ in PHASES if i == pid)
        ap = ex["prime"]["phases"].get(fitkey, {}).get("alpha") if fitkey != "—" else None
        ab = ex["char2"]["phases"].get(fitkey, {}).get("alpha") if fitkey != "—" else None
        sk = head_kob["phases_s"].get(pid)
        sp = head_pr["phases_s"].get(pid)
        share_k = f"{100 * sk / head_kob['s']:.1f}%" if sk else "—"
        share_p = f"{100 * sp / head_pr['s']:.1f}%" if sp else "—"
        P.append(f'<tr><td><i style="display:inline-block;width:9px;height:9px;border-radius:2px;background:var(--ph-{pid});margin-right:6px"></i>{esc(name)}</td>'
                 f'<td>{esc(what)}</td><td>{esc(model)}</td>'
                 f'<td class="n">{g3(ap) if ap is not None else "—"}</td><td class="n">{g3(ab) if ab is not None else "—"}</td>'
                 f'<td class="n">{share_k}</td><td class="n">{share_p}</td></tr>')
    P.append(f'<tr><td><strong>Whole method</strong></td><td>all of the above, cold</td><td>family optimum S ∝ r<sup>1/6</sup></td>'
             f'<td class="n"><strong>{g3(ex["prime"]["phases"]["total"]["alpha"])}</strong></td>'
             f'<td class="n"><strong>{g3(ex["char2"]["phases"]["total"]["alpha"])}</strong></td><td class="n">100%</td><td class="n">100%</td></tr>'
             f'<tr><td>Matched rho</td><td>negation-map r-adding walk, distinguished points</td><td>√(πr/2A)</td>'
             f'<td class="n">{g3(ex["prime"]["rho_matched"]["alpha"])}</td><td class="n">{g3(ex["char2"]["rho_matched"]["alpha"])}</td><td class="n">—</td><td class="n">—</td></tr>')
    P.append('</tbody></table></div>'
             f'<p>At the Koblitz batch headline (§22) the pipeline stored {head_kob["stored_pairs"]:,} pairs, collected '
             f'{head_kob["relations"]:,} relations from {head_kob["summands_scanned"]:,} scanned summands, reduced '
             f'{head_kob["la_rows_in"]} rows to a Wiedemann core of dimension {head_kob["la_core_dimension"]}, and '
             f'spent {g3(head_kob["descent_trials_per_target"])} descent trials per target. Koblitz top-end exponent, '
             f'ratio·√n against r over the four largest sizes: β = {g3(ex["koblitz"]["beta"])} '
             f'[{g3(ex["koblitz"]["ci"][0])}, {g3(ex["koblitz"]["ci"][1])}], against the law\'s 1/6 and the model\'s '
             f'{ex["koblitz"]["model"]}; it cannot yet tell them apart. The ladder\'s rho fits sit under ½ because '
             'set-up still shows at these sizes (ledger §18.5).</p></section>')

    # PDP oracles.
    oracle_names = []
    for o in doc["oracles"]:
        for k in o["oracles"]:
            if k not in oracle_names:
                oracle_names.append(k)
    P.append('<section class="panel" id="pdp"><h2>The point decomposition problem, per target</h2>'
             '<p>Every oracle answers the same question on the same Semaev systems: is this target a sum of m base '
             'points? Cost is group-addition equivalents per target, found and refuted averaged; FFD is the first '
             'fall degree of the Weil-restricted system. The pair table answers in 0–25 units where F4 needs 10<sup>2</sup>–10<sup>6</sup> '
             'and SAT more, which is why the board\'s recipes use it. Zero disagreements on every cell.</p>'
             '<div class="scroll"><table><thead><tr><th>curve</th><th class="n">n</th><th class="n">m</th>'
             '<th class="n">ℓ</th><th class="n">unknowns</th><th class="n">FFD</th><th class="n">hit rate</th>'
             + "".join(f'<th class="n">{esc(k)}</th>' for k in oracle_names) + '</tr></thead><tbody>')
    for o in doc["oracles"]:
        ffd = "—" if o["ffd"][0] is None else (str(o["ffd"][0]) if o["ffd"][0] == o["ffd"][1] else f"{o['ffd'][0]}–{o['ffd'][1]}")
        P.append(f'<tr><td><span class="slug">{esc(o["slug"])}</span></td><td class="n">{o["n"]}</td><td class="n">{o["m"]}</td>'
                 f'<td class="n">{o["dimension"]}</td><td class="n">{o["unknowns"]}</td><td class="n">{ffd}</td>'
                 f'<td class="n">{o["hit_rate"]:.2f}</td>')
        for k in oracle_names:
            v = o["oracles"].get(k)
            if v is None:
                P.append('<td class="n">—</td>')
            elif v["inconclusive"] and not (v["found"] or v["refuted"]):
                P.append('<td class="n muted">no answer</td>')
            else:
                P.append(f'<td class="n">{g3(v["gae_per_target"])}</td>')
        P.append('</tr>')
    P.append(f'</tbody></table></div><p class="muted">Source: <code>{esc(SOURCES["oracles"])}</code>, a stage '
             'diagnostic: one phase priced alone is not a speed (AGENTS.md §2).</p></section>')

    # The curve roster.
    fam_order = ("koblitz", "subfield", "binary", "prime")
    P.append('<section class="panel" id="curves"><h2>Every curve the repository names</h2>'
             '<p>From the curve registry. The slug is the curve\'s name in text; the EC1 alias identifies the exact '
             'representation measured (subgroup and generator included) for comparisons across repositories. '
             'Legacy spellings are what older notes and frozen reports called the curve.</p>')
    for fam in fam_order:
        rows = [c for c in doc["roster"] if c["family"] == fam]
        if not rows:
            continue
        P.append(f'<details {"open" if fam == "koblitz" else ""}><summary>{esc(fam)} · {len(rows)} curves · '
                 f'{sum(c["on_board"] for c in rows)} on the board</summary><div class="scroll"><table><thead><tr>'
                 '<th>slug</th><th>field</th><th>coefficients</th><th class="n">#E bits</th><th>EC1</th>'
                 '<th>also called</th><th>board</th></tr></thead><tbody>')
        for c in rows:
            ec1 = ", ".join(c["ec1"]) if c["ec1"] else f'<span class="muted">{esc(c["ec1_unresolved"])}</span>'
            also = ", ".join(c["standard"] + c["legacy"])
            P.append(f'<tr><td><span class="slug">{esc(c["slug"])}</span></td><td>{esc(c["field"])}</td>'
                     f'<td>{esc(c["coefficients"])}</td><td class="n">{c["order_bits"]}</td>'
                     f'<td><code>{ec1 if c["ec1"] else ""}</code>{"" if c["ec1"] else ec1}</td>'
                     f'<td>{esc(also)}</td><td>{"✓" if c["on_board"] else ""}</td></tr>')
        P.append('</tbody></table></div></details>')
    P.append('</section>')

    # Footer.
    P.append('<footer><p><strong>What this is not.</strong> A view of frozen measurements, class accounting: no row '
             'here was measured for this page, and nothing is a speedup. Records kept in other units, wall-clock '
             'panels and batched producers awaiting independent replay, are in '
             '<a href="https://github.com/aburan28/crypto/blob/main/docs/ic/BOUNDARY_TARGETS.md">docs/ic/BOUNDARY_TARGETS.md</a>; the canonical page is the '
             '<a href="https://aburan28.github.io/crypto/scoreboard/">index-calculus scoreboard</a>. No IC pipeline has been priced '
             'end to end on ECC2K-130 or on the m = 83 confidence gate (AGENTS.md §8a).</p>'
             f'<p>Current through ledger §{LEDGER_COVERED_THROUGH}. Built by '
             '<code>scripts/build_ic_leaderboard.py</code> from: '
             + ", ".join(f'<code>{esc(v["path"])}</code> ({v["sha256"][:12]})' for v in doc["sources"].values())
             + '. Curves are named by ICV1 slug (<a href="https://github.com/aburan28/crypto/blob/main/docs/curves/ICV1.md">docs/curves/ICV1.md</a>).</p></footer></div>')
    body = "\n".join(P)
    if standalone:
        return ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
                '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
                '<meta name="description" content="Every curve priced end to end by index calculus in this '
                'repository, with the best recipe, the phase split and the ratio to matched Pollard rho.">\n'
                f"{head}</head>\n<body>\n{body}\n</body>\n</html>\n")
    return head + body + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--artifact", type=Path)
    args = ap.parse_args()
    doc = build()
    outs = {OUT_JSON: json.dumps(doc, indent=1, ensure_ascii=False) + "\n",
            OUT_MD: markdown(doc), OUT_HTML: page(doc, standalone=True)}
    if args.check:
        sections = [int(m) for m in re.findall(r"^## (\d+)\.", (REPO / LEDGER).read_text(), re.M)]
        latest = max(sections, default=0)
        if latest > LEDGER_COVERED_THROUGH:
            print(f"{LEDGER} has §{latest}; the leaderboard is current through §{LEDGER_COVERED_THROUGH}. "
                  "Point SOURCES at the round's frozen files and regenerate, or confirm the round "
                  "changes nothing on the page; then raise LEDGER_COVERED_THROUGH (AGENTS.md §7b).",
                  file=sys.stderr)
            return 1
        stale = [str(p.relative_to(REPO)) for p, t in outs.items()
                 if not p.exists() or p.read_text() != t]
        for s in stale:
            print(f"{s} is stale; run scripts/build_ic_leaderboard.py", file=sys.stderr)
        return 1 if stale else 0
    for p, t in outs.items():
        p.write_text(t)
        print(f"wrote {p.relative_to(REPO)}")
    if args.artifact:
        args.artifact.write_text(page(doc, standalone=False))
        print(f"wrote {args.artifact}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

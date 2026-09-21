#!/usr/bin/env python3
"""Update docs/ic/boundary_targets.json from a frozen boundary-ledger run.

usage: boundary_ledger_update.py <run.json> <boundary_targets.json> <run path for evidence> [--date=YYYY-MM-DD] [--dry]

Moves the binary and prime relation_yield / rank / end_to_end_dlp / vs_rho
records to the run's figures (the previous current block goes to history),
adds the Koblitz oracle-pricing and operation-counted vs_rho metrics, and
appends the ledger's agent priority.  Written for the 2026-09-21 run; a later
round should check the row texts it writes before running it again.
"""
import json, math, sys, copy
from collections import OrderedDict

run_path, ledger_path, evidence_path = sys.argv[1], sys.argv[2], sys.argv[3]
dry = "--dry" in sys.argv
run = json.load(open(run_path))
ledger = json.load(open(ledger_path, encoding="utf-8"))
L = run["ledger"]
DATE = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--date=")), "2026-09-21")
NOTE = "RESEARCH_IC_BOUNDARY_LEDGER.md"

def mean(xs):
    xs = [x for x in xs if x is not None and not (isinstance(x, float) and math.isnan(x))]
    return sum(xs) / len(xs) if xs else None

def r3(x):
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return None
    return float(f"{x:.4g}")

def variants_of(inst):
    seen = OrderedDict()
    for v in inst["variants"]:
        seen.setdefault(v["name"], []).append(v)
    return seen

def instance_rows(regime):
    rows = []
    for inst in L["instances"]:
        if inst["regime"] != regime:
            continue
        sqrt_r = math.sqrt(inst["r"])
        A = inst["automorphisms_generic"]
        ref = [x for x in inst["rho"] if (A <= 2 or x["automorphisms"] == A)]
        for name, runs in variants_of(inst).items():
            m = lambda k: mean([r[k] for r in runs])
            rows.append(OrderedDict([
                ("instance", inst["curve"]["name"]),
                ("log2_r", r3(inst["log2_r"])),
                ("r", inst["r"]),
                ("group_order", inst["group_order"]),
                ("cofactor", inst["cofactor"]),
                ("variant", name),
                ("m", runs[0]["summands"]),
                ("signed_points", runs[0]["signed_points"]),
                ("columns", runs[0]["columns"]),
                ("dimension", runs[0]["dimension"]),
                ("trials", r3(m("trials"))),
                ("relations", r3(m("relations_found"))),
                ("trials_per_relation", r3(m("trials") / max(m("relations_found"), 1))),
                ("decomposition_probability_ceiling", r3(runs[0]["decomposition_probability_ceiling"])),
                ("yield_over_ceiling", r3(m("yield_over_ceiling"))),
                ("matrix_rows", r3(m("relations_found"))),
                ("matrix_columns", runs[0]["columns"] + 1),
                ("rank", r3(m("rank"))),
                ("la_row_ops", r3(mean([r["linear_algebra"]["native"].get("row_ops", 0) for r in runs]))),
                ("la_gae", r3(mean([r["linear_algebra"]["gae"] for r in runs]))),
                ("factor_base_gae", r3(mean([r["factor_base"]["gae"] for r in runs]))),
                ("relations_gae", r3(mean([r["relations"]["gae"] for r in runs]))),
                ("total_gae", r3(m("total_gae"))),
                ("S", r3(m("s"))),
                ("S_wall", r3(m("s_wall"))),
                ("rho_S_mean", r3(inst["rho_s_mean"])),
                ("rho_walk_S_mean", r3(mean([x["s_walk"] for x in ref]))),
                ("S_over_rho", r3(m("s") / inst["rho_s_mean"])),
                ("floor_S", r3(inst["floor_s"])),
                ("S_over_floor", r3(m("s") / inst["floor_s"])),
                ("verified_all", all(r["verified"] for r in runs)),
                ("rho_verified_all", inst["rho_verified_all"]),
                ("repeats", len(runs)),
            ]))
    return rows

def fits_of(regime, phase=None):
    out = []
    for f in L["fits"]:
        if f["regime"] == regime and (phase is None or f["phase"] == phase):
            out.append(OrderedDict([("variant", f["variant"]), ("phase", f["phase"]), ("alpha_r", r3(f["alpha"])), ("r_squared", r3(f["r_squared"])), ("alpha_group_order", r3(f["alpha_group_order"])), ("points", f["points"])]))
    return out

def push_history(stage):
    cur = stage.get("current")
    if cur is not None:
        h = copy.deepcopy(cur)
        h["superseded_on"] = DATE
        stage.setdefault("history", []).insert(0, h)

def evidence(stage):
    ev = stage.setdefault("evidence", [])
    for e in (evidence_path, NOTE):
        if e not in ev:
            ev.append(e)

regimes = ledger["regimes"]

# ── binary (regime A) and prime (regime C): the same four stages ──
for key, regime_name in (("binary", "char2"), ("prime", "prime")):
    rec = regimes[key]["records"]
    rows = instance_rows(regime_name)
    if not rows:
        continue
    largest = max(rows, key=lambda r: r["r"])
    best = min(rows, key=lambda r: r["S_over_rho"])
    n_or_bits = "n" if key == "binary" else "bits"
    label = (lambda r: f"n={int(round(math.log2(r['group_order'])))}") if key == "binary" else (lambda r: f"{r['log2_r']:.1f} bits")

    st = rec["relation_yield"]
    push_history(st)
    st["status"] = "record"
    st["current"] = OrderedDict([
        ("summary", f"Operation-counted yield on the frozen {regime_name} ladder ({len(set(r['instance'] for r in rows))} instances, up to r=2^{largest['log2_r']:.1f}): relations per trial measured against the counting ceiling C(F+m-1,m)/#E for every variant; natural targets only"),
        ("date", DATE),
        ("metrics", OrderedDict([
            ("timing_class", "operation_counted"),
            ("target_mix", {"natural": "every trial", "planted_sat": 0, "proven_unsat": 0}),
            ("yield_over_ceiling_range", [r3(min(r["yield_over_ceiling"] for r in rows)), r3(max(r["yield_over_ceiling"] for r in rows))]),
            ("trials_per_relation_by_instance", [OrderedDict([("instance", r["instance"]), ("variant", r["variant"]), ("m", r["m"]), ("signed_points", r["signed_points"]), ("trials_per_relation", r["trials_per_relation"]), ("pr_decomposition_measured", r3(1.0 / r["trials_per_relation"]) if r["trials_per_relation"] else None), ("ceiling", r["decomposition_probability_ceiling"]), ("yield_over_ceiling", r["yield_over_ceiling"])]) for r in rows]),
            ("trials_exponent_fits", fits_of(regime_name, "trials")),
        ])),
    ])
    evidence(st)

    st = rec["rank"]
    push_history(st)
    st["status"] = "record"
    st["current"] = OrderedDict([
        ("summary", f"Relation-matrix LA priced per instance on the {regime_name} ladder: dense incremental Gauss-Jordan over Z/rZ, rank recomputed after every row, stop when the target column is pinned; multiply-subtracts counted and converted"),
        ("date", DATE),
        ("metrics", OrderedDict([
            ("sparse_or_dense", "dense incremental reduced echelon form (IncrementalGauss)"),
            ("rank_accumulation", "reduced after every relation; terminal when the logarithm's column is pinned"),
            ("by_instance", [OrderedDict([("instance", r["instance"]), ("variant", r["variant"]), ("matrix_dims", {"rows": r["matrix_rows"], "cols": r["matrix_columns"]}), ("rank", r["rank"]), ("row_ops", r["la_row_ops"]), ("la_charged_gae", r["la_gae"]), ("la_S", r3(r["la_gae"] / math.sqrt(r["r"])))]) for r in rows]),
            ("la_exponent_fits", fits_of(regime_name, "linear_algebra")),
        ])),
    ])
    evidence(st)

    st = rec["end_to_end_dlp"]
    push_history(st)
    st["status"] = "record"
    st["current"] = OrderedDict([
        ("summary", f"Known-answer DLP recovered and verified ([d]G = Q) by every variant on every instance of the {regime_name} ladder, largest {label(largest)} (r=2^{largest['log2_r']:.1f}); every phase counted in one unit"),
        ("date", DATE),
        ("metrics", OrderedDict([
            (f"{n_or_bits}_max", int(round(math.log2(largest["group_order"]))) if key == "binary" else r3(largest["log2_r"])),
            ("claim_boundary", "synthetic_known_answer"),
            ("stage_counts", "factor_base / relations / linear_algebra / verify per variant in the frozen run (native counts, wall_ns, gae)"),
            ("recovered_d_verified", all(r["verified_all"] for r in rows)),
            ("variants", sorted(set(r["variant"] for r in rows))),
        ])),
    ])
    evidence(st)

    st = rec["vs_rho"]
    push_history(st)
    st["status"] = "not_achieved"
    st["current"] = OrderedDict([
        ("summary", f"Whole-process operation-counted cost against a counted Pollard rho on the same instance: best variant {best['variant']} at {best['S_over_rho']}x rho ({label(best)}); no variant below rho at r >= 2^20"),
        ("date", DATE),
        ("metrics", OrderedDict([
            ("timing_class", "operation_counted_whole_process"),
            ("automorphism_discount", "reference walk uses A=1 (no negation map); floor stated at A=2"),
            ("by_instance", [OrderedDict([("instance", r["instance"]), ("log2_r", r["log2_r"]), ("variant", r["variant"]), ("S", r["S"]), ("rho_S", r["rho_S_mean"]), ("rho_walk_S", r["rho_walk_S_mean"]), ("S_over_rho", r["S_over_rho"]), ("S_over_floor", r["S_over_floor"]), ("verified", r["verified_all"]), ("rho_verified", r["rho_verified_all"])]) for r in rows]),
            ("total_exponent_fits", fits_of(regime_name, "total")),
            ("rho_health", "every rho run recovered and verified; steps against sqrt(pi r/2) recorded per run"),
        ])),
        ("verdict", False),
    ])
    evidence(st)

# ── koblitz: decomposition (FFD cells) and vs_rho (counted ratios) ──
rec = regimes["koblitz"]["records"]
rows = instance_rows("koblitz")
if rows:
    st = rec["vs_rho"]
    st["current"]["metrics"]["operation_counted_whole_process_2026_09_21"] = OrderedDict([
        ("timing_class", "operation_counted_whole_process"),
        ("automorphism_discount", "signed Frobenius walk, A=2n, measured (koblitz_signed_frobenius_rho_reference)"),
        ("by_instance", [OrderedDict([("instance", r["instance"]), ("log2_r", r["log2_r"]), ("cofactor", r["cofactor"]), ("variant", r["variant"]), ("signed_points", r["signed_points"]), ("columns", r["columns"]), ("S", r["S"]), ("rho_S", r["rho_S_mean"]), ("S_over_rho", r["S_over_rho"]), ("S_over_floor", r["S_over_floor"]), ("verified", r["verified_all"]), ("rho_verified", r["rho_verified_all"])]) for r in rows]),
        ("total_exponent_fits", fits_of("koblitz", "total")),
        ("reading", "single-target whole-process counts on a materialised base; the n=53 record above is a wall-clock charged class and is unchanged"),
    ])
    evidence(st)

op = run.get("oracle_pricing")
if op:
    st = rec["decomposition"]
    cells = []
    for c in op["cells"]:
        cells.append(OrderedDict([
            ("n", c["n"]), ("dimension", c["dimension"]), ("m", c["m"]), ("signed_points", c["signed_points"]), ("signed_orbits", c["signed_orbits"]),
            ("unknowns", c["unknowns"]), ("equations", c["equations"]), ("system_degree", c["system_degree"]), ("eq_var_ratio", r3(c["eq_var_ratio"])),
            ("ffd_min", c["ffd_min"]), ("ffd_max", c["ffd_max"]), ("ffd_no_fall", c["ffd_no_fall"]), ("macaulay", c["macaulay"]),
            ("targets", c["targets"]), ("hit_rate", r3(c["hit_rate"])), ("disagreements", c["disagreements"]),
            ("oracles", [OrderedDict([("oracle", o["oracle"]), ("native_unit", o["native_unit"]), ("found", o["found"]), ("refuted", o["refuted"]), ("inconclusive", o["inconclusive"]), ("native_median_refuted", r3(o["native_median_refuted"])), ("ms_median_refuted", r3(o["ms_median_refuted"])), ("gae_per_target_mean", r3(o["gae_per_target_mean"])), ("projected_relation_phase_S", r3(p["s_projected"])), ("projected_S_over_floor", r3(p["s_over_rho_floor"])), ("extrapolation", True)]) for o, p in zip(c["oracles"], c["projected"])]),
        ]))
    st["current"]["metrics"]["oracle_pricing_2026_09_21"] = OrderedDict([
        ("what", "every decomposition oracle priced per target on the same Semaev systems, with FFD over the target draws; F4 word XORs exact (F4_WORD_OPS_TOTAL), SAT conflicts, enumerate additions, pair-table probes, S4 pairs"),
        ("all_agree", op["all_agree"]),
        ("cells", cells),
    ])
    evidence(st)

ledger["updated"] = DATE
pri = ledger["agent_priorities"]
if not any("RESEARCH_IC_BOUNDARY_LEDGER" in p.get("beat", "") for p in pri):
    pri.append(OrderedDict([
        ("beat", "Iterate against RESEARCH_IC_BOUNDARY_LEDGER.md (ic boundary): in any regime, S/S_rho < 1 at r >= 2^20 with every phase counted, or a fitted total exponent below the reference's over >= 4 sizes, or yield/ceiling > 1.5 on a base outside every proper subgroup"),
        ("rank", len(pri) + 1),
        ("regime", "all"),
        ("stage", "vs_rho"),
    ]))

if dry:
    print(json.dumps({k: regimes[k]["records"]["vs_rho"]["current"]["summary"] if k != "koblitz" else "koblitz: metrics added" for k in regimes}, indent=1))
    print("prime relation_yield metrics keys:", list(regimes["prime"]["records"]["relation_yield"]["current"]["metrics"].keys()))
    print("priorities:", len(pri))
else:
    with open(ledger_path, "w", encoding="utf-8") as fh:
        json.dump(ledger, fh, indent=2, ensure_ascii=False)
        fh.write("\n")
    print("written", ledger_path)

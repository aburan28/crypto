#!/usr/bin/env python3
"""Bind already computed native results to presentation tables and figures.

This script selects fields and formats text/units. It does not run experiments,
calculate statistics, form ratios, estimate uncertainty, check mathematical
certificates, or recompute phase/model accounting. Its inputs are copied native
records, with the native independent audit and analysis gates already present.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BRANCH = "codex/n83-rho-speed-20261012"
STUDY = "research/n83_matched_20261014"


def number(value, precision=6):
    if value is None:
        return "not established"
    return f"{value:.{precision}f}" if isinstance(value, (float, int)) else str(value)


def interval(row):
    bounds = row.get("interval_95")
    return "not admitted" if not bounds else "[" + ", ".join(number(x) for x in bounds) + "]"


def triple(values):
    return ", ".join(number(v, 6) for v in values)


def link(name, path, directory=False):
    mode = "tree" if directory else "blob"
    url = f"https://github.com/aburan28/crypto/{mode}/{BRANCH}/{STUDY}/{path}"
    return f'<link href="{url}" color="#087F83">{name}</link>'


def source_note(*sources):
    return "Sources: " + ", ".join(link(name, path, directory) for name, path, directory in sources) + "."


def rewrite(value, replacements):
    if isinstance(value, str):
        return replacements.get(value, value)
    if isinstance(value, list):
        return [rewrite(v, replacements) for v in value]
    if isinstance(value, dict):
        return {k: rewrite(v, replacements) for k, v in value.items()}
    return value


def table(section, rows, headers=None, ordinal=0):
    target = [b for b in section["blocks"] if b["type"] == "table"][ordinal]
    target["rows"] = rows
    if headers:
        target["headers"] = headers


def chart(name, kind, title, ylabel, population, series, **extra):
    value = {"name": name, "kind": kind, "title": title, "ylabel": ylabel,
             "population": population,
             "source": "run-20261014T1940Z · records/postprocess-20261015T0018Z/statistics.json",
             "series": series}
    value.update(extra)
    return value


def bar_series(label, x, y, offset):
    return {"label": label, "x": x, "y": y,
            "display_offset": offset, "display_width": 0.35}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--statistics", type=Path, required=True)
    parser.add_argument("--analysis-input", type=Path, required=True)
    parser.add_argument("--replay", type=Path, required=True)
    parser.add_argument("--supplement", type=Path, required=True)
    parser.add_argument("--template", type=Path,
                        default=ROOT.parent/"source-preparation/report_content.template.json")
    parser.add_argument("--presentation-metadata", type=Path, required=True)
    args = parser.parse_args()
    stats = json.loads(args.statistics.read_text())
    inputs = json.loads(args.analysis_input.read_text())
    replay = json.loads(args.replay.read_text())
    supplement = json.loads(args.supplement.read_text())
    meta = json.loads(args.presentation_metadata.read_text())
    if stats.get("all_required_rows_valid") is not True or replay.get("status") != "PASS":
        raise SystemExit("Native independent-analysis gates are not complete; report binding withheld")
    # The gate above recognizes native receipts; it does not independently verify
    # their scientific content or their file bindings.
    cold = stats["primary_actual_cold_singleton_process_wall"]
    shared = stats["primary_prepared_batch_process_wall"]
    online = stats["secondary_prepared_online_target_means"]
    crossed = stats["secondary_crossed_online_sensitivity"]
    geometric = stats["secondary_geometric_mean_paired_target_online_ratio"]
    model = stats["modeled_break_even"]
    aa = stats["aa"]["wall_ratio_description"]
    rows = stats["per_target"]
    costs = stats["costs"]
    cold_pairs = stats["display_cold_pairs"]
    cohorts = stats["display_cohorts"]
    process_rows = inputs["process_rows"]
    display_status = "NATIVE_RECORDS_VERIFIED_LAYOUT_ONLY"

    replacements = {
        "PENDING:NATIVE_RECORDS_VERIFIED_LAYOUT_ONLY": display_status,
        "14 October 2026 — PENDING:ACTUAL_CAMPAIGN_ID": "14–15 October 2026 · run-20261014T1940Z",
        "PENDING:OBSERVED_PRIMARY_VALUES_AND_CONDITIONAL_SCOPE":
            f"Actual IC/rho mean complete-process ratio: {number(cold['ic_over_rho'])} for eight cold singleton pairs; "
            f"{number(shared['ic_over_rho'])} for three genuine 32-target cohort pairs. Values above 1 favor rho. "
            "The 32-target query comparison and crossover model remain secondary.",
        "PENDING:NATIVE_PRIMARY_BATCH_AND_COLD_CONCLUSION_WITH_COUNTS":
            f"All 32 solver processes completed and all 218 answers passed independent scalar replay. "
            f"Actual cold singleton means are IC {number(cold['arithmetic_mean_ic_s'])} s and rho {number(cold['arithmetic_mean_rho_s'])} s; "
            f"actual 32-target cohort means are IC {number(shared['arithmetic_mean_ic_s'])} s and rho {number(shared['arithmetic_mean_rho_s'])} s.",
        "PENDING:FAILURE_COMPLETENESS_AND_AA_ADMISSION_SCOPE": meta["completion_and_admission"],
        "PENDING:EC1_UID_POINTSET_HASH_SOURCE_BINDINGS": meta["identity_details"],
        "PENDING:ACTUAL_COMPLETED_FAILURES_GENERATION_RECEIPT_AND_SEEDS": meta["fixtures_and_seeds"],
        "PENDING:SETUP_PHASE_VALUES_RESIDUAL_AND_EXTERNAL_AUDIT_METERS": meta["phase_and_audit_boundary"],
        "PENDING:NATIVE_COLD_MEANS_MEDIANS_MIN_MAX_INTERVAL_AND_CORRECTNESS":
            f"Native paired eight-target percentile 95% interval for the mean-cost ratio: {interval(cold)}. "
            f"IC process median {number(costs['ic']['actual_cold_singleton_process']['wall_s']['median'])} s "
            f"(min {number(costs['ic']['actual_cold_singleton_process']['wall_s']['min'])}, "
            f"max {number(costs['ic']['actual_cold_singleton_process']['wall_s']['max'])}); "
            f"rho median {number(costs['rho']['actual_cold_singleton_process']['wall_s']['median'])} s "
            f"(min {number(costs['rho']['actual_cold_singleton_process']['wall_s']['min'])}, "
            f"max {number(costs['rho']['actual_cold_singleton_process']['wall_s']['max'])}). "
            "All eight original point pairs are included; given retained K600, without OS-cache flushing.",
        "PENDING:NATIVE_SHARED_RATIO_CI_COMPLETE_AND_RESOURCE_STATUS":
            f"Native mean actual process ratio IC/rho {number(shared['ic_over_rho'])}, "
            f"paired-cohort 95% interval {interval(shared)}. Each arm has exactly three observed whole-process costs. "
            "No process timeout, failure or target omission is converted into a zero-cost observation.",
        "PENDING:NATIVE_TARGETBLOCK_AND_CROSSED_SENSITIVITY_INTERVALS_WITH_FAILURES":
            f"Secondary ratio of arithmetic mean target costs IC/rho {number(online['ic_over_rho'])}, "
            f"32-target-block 95% interval {interval(online)}. Crossed target/cohort sensitivity interval {interval(crossed)}. "
            f"Descriptive geometric mean of paired target-mean ratios {number(geometric['ic_over_rho'])}. "
            "These clocks have unequal validation boundaries and do not establish the complete-process winner.",
        "PENDING:ACTUAL_RESOURCE_FAILURE_SUMMARY_WITH_RAW_LINKS": meta["resources_and_failures"],
        "PENDING:NATIVE_AA_VALUES_AND_CONDITIONAL_WALL_INTERPRETATION":
            f"A/A process A/B ratios: {triple(stats['aa']['A_over_B_wall_ratios'])}. "
            f"Native median {number(aa['median'])}, minimum {number(aa['min'])}, maximum {number(aa['max'])}. "
            "Five pairs use identical points, seeds, policies and exported semantic work, independently checked. "
            "No noise correction is subtracted from observed method costs.",
        "PENDING:VERIFICATION_SCOPE_FAILURE_PRESERVATION_AND_SOURCE_FREEZE": meta["verification_details"],
        "PENDING:EVIDENCE_BOUNDED_METHOD_CHOICE_FOR_SINGLETON_AND_32_TARGET_WORKLOADS": meta["method_choice"],
        "PENDING:SOURCES_INDEX_REPRODUCTION_AND_BENCHMARK_COPY_LINKS": meta["reproduction_details"],
        "PENDING:MODELED_THRESHOLD_INTERVAL_AND_ACTUAL_PREFIX_DISCLOSURE": meta["model_details"],
    }
    content = rewrite(json.loads(args.template.read_text()), replacements)
    sections = content["sections"]
    table(sections[0], [
        ["Eight fresh singleton pairs", number(cold["arithmetic_mean_ic_s"]), number(cold["arithmetic_mean_rho_s"]), number(cold["ic_over_rho"]) + " / " + interval(cold)],
        ["Three real 32-target cohort pairs", number(shared["arithmetic_mean_ic_s"]), number(shared["arithmetic_mean_rho_s"]), number(shared["ic_over_rho"]) + " / " + interval(shared)],
    ])
    table(sections[4], [[r["target_id"], number(r["ic"]["elapsed_process_wall_s"]),
                          number(r["rho"]["elapsed_process_wall_s"]),
                          number(r["process_wall_ic_over_rho"]) + " / verified"] for r in cold_pairs])
    by_cohort = {(r["cohort_id"], r["arm"]): r for r in cohorts}
    table(sections[5], [[str(cohort),
                         number(by_cohort[cohort, "ic"]["process"]["elapsed_process_wall_s"]) + " / " + number(by_cohort[cohort, "ic"]["preparation_s"]),
                         number(by_cohort[cohort, "rho"]["process"]["elapsed_process_wall_s"]) + " / " + number(by_cohort[cohort, "rho"]["preparation_s"]),
                         "Both verified; all 32 points"] for cohort in (1, 2, 3)],
          headers=["Cohort", "IC process / setup (s)", "Rho process / setup (s)", "Complete status"])
    table(sections[6], [[label, number(costs["ic"][scope][stat]), number(costs["rho"][scope][stat]), population]
                        for scope, population in (("32_target_three_repeat_online_means_s", "32 target means, three repeats each"),
                                                   ("pooled96_online_draws_s", "96 draws sharing three preparations"))
                        for stat, label in (("mean", "Mean query (s)"), ("median", "Median query (s)"),
                                            ("p90", "p90 query (s)"), ("p95", "p95 query (s)"),
                                            ("max", "Maximum query (s)"))])
    for section, part in ((sections[7], rows[:16]), (sections[8], rows[16:])):
        table(section, [[r["target_id"], triple(r["ic_online_s"]), triple(r["rho_online_s"]),
                         number(r["ic"]["mean"]) + " / " + number(r["rho"]["mean"]) + "\nratio " + number(r["paired_mean_ic_over_rho"])] for r in part],
              headers=["Target", "IC queries (s)", "Rho queries (s)", "Native means (s) / IC/rho"])
        next(block for block in section["blocks"] if block["type"] == "table")["widths"] = [0.65,1.45,1.45,1.3]
    table(sections[9], [[f"{arm.upper()} / {label}",
                         "CPU mean " + number(costs[arm][scope]["cpu_total_s"]["mean"]) + " s\nRSS min / max " + number(costs[arm][scope]["peak_rss_bytes"]["min"], 0) + " / " + number(costs[arm][scope]["peak_rss_bytes"]["max"], 0) + " bytes",
                         "All original processes verified; costs retained"]
                        for label, scope in (("8 cold", "actual_cold_singleton_process"),
                                             ("3 shared", "prepared_cohort_process")) for arm in ("ic", "rho")])
    for block in sections[9]["blocks"]:
        if block.get("type") == "paragraph" and block.get("text", "").startswith("Timeouts, OOMs"):
            block["text"] = "Each run retains exit/signal status, timeout flags, the raw memory-limit indicator, scalar-validation status and internal rho capped-walk counters. A null memory-limit indicator remains unavailable. Internal rho lane restarts are charged algorithmic work; zero whole-process timeouts does not imply zero internal restarts."
    main_rho_work = supplement["main104_rho_work_counter_totals"]
    all_rho_work = supplement["all114_rho_work_counter_totals"]
    aa_work = next(g["rho_work_counter_totals"] for g in supplement["group_summaries"]
                   if g["arm"] == "rho" and g["kind"] == "aa")
    sections[9]["blocks"].append({"type":"small", "text":
        "Native rho work counters: the 104 main answers retain " + str(main_rho_work["totals"]["walk_steps"]) +
        " walk steps and " + str(main_rho_work["totals"]["capped_walks"]) +
        " internal capped-walk restarts; the 10 A/A answers retain " + str(aa_work["totals"]["capped_walks"]) +
        " restarts. All 114 rho answers retain " + str(all_rho_work["totals"]["capped_walks"]) +
        " restarts. These algorithmic events and their work remain charged. They are distinct from whole-process timeout or global-cap failure."})
    table(sections[10], meta["model_table_rows"])
    table(sections[11], meta["host_table_rows"])
    table(sections[12], meta["requirement_rows"])
    table(sections[13], meta["measurement_changes"])
    if meta.get("verifier_revision_details"):
        sections[12]["blocks"].append({"type":"small", "text":meta["verifier_revision_details"]})

    # Citations are fixed report references to native records, never results.
    for index, section in enumerate(sections):
        section["source_note"] = meta["section_source_notes"][str(index)]
        for block in section["blocks"]:
            if block["type"] == "figure":
                block["caption"] = meta["figure_captions"][Path(block["path"]).stem]
    phase_labels = [
        ("target_query_ms", "Query-origin preparation"),
        ("target_pdp_ms", "Point decomposition excluding relation checks"),
        ("target_relation_check_ms", "Relation checks"),
        ("target_descent_ms", "Scalar descent / recovery"),
        ("target_recovery_check_ms", "Recovery group check"),
    ]
    phase_rows = [[label, number(supplement["ic104_phase_ms"][key]["mean"]),
                   number(supplement["ic104_phase_ms"][key]["median"]),
                   number(supplement["ic104_phase_ms"][key]["p95"]),
                   number(supplement["ic104_phase_ms"][key]["max"])]
                  for key, label in phase_labels]
    residual = supplement["ic104_target_minus_phase_sum_ms"]
    phase_rows.append(["Target clock minus five-phase sum", number(residual["mean"]),
                       number(residual["median"]), number(residual["p95"]), number(residual["max"])])
    phase_section = {
        "title":"Recorded phase accounting and residuals",
        "blocks":[
            {"type":"paragraph", "text":"All 104 main IC answers retain the five producer stage fields: 96 queries in shared cohorts and eight singleton queries. The native supplement computes the signed target-clock-minus-stage-sum residual. These fields describe the producer clock; they do not repair its omitted Q-validation/input/output boundary."},
            {"type":"table", "headers":["Native IC field (ms)","Mean","Median","p95","Maximum"],
             "rows":phase_rows,"widths":[2,1,1,1,1],"padding":4},
            {"type":"small", "text":meta["phase_details"]},
            {"type":"paragraph", "text":"Rho retains exact nanosecond phase and residual counters in its own raw target records, with A/A cases identified separately. Process-level supervisor-minus-producer residuals, CPU and RSS remain in the native supplement with source-file hashes and JSON/JSONL pointers. No phase sum or residual is recalculated by the PDF renderer."}
        ],
        "source_note":meta["phase_source_note"]
    }
    sections.insert(4, phase_section)
    (ROOT/"report_content.json").write_text(json.dumps(content, indent=2, ensure_ascii=False) + "\n")

    figures = []
    positions = list(range(1, 9))
    figures.append(chart("cold_pairs", "grouped_bars", "Actual one-target cold process costs", "Whole-process wall (s; log scale)",
                          "Eight actual paired public targets · retained K600 input", [
                              bar_series("IC", positions, [r["ic"]["elapsed_process_wall_s"] for r in cold_pairs], -0.18),
                              bar_series("Optimized rho", positions, [r["rho"]["elapsed_process_wall_s"] for r in cold_pairs], 0.18)],
                          xlabel="Predeclared cold target (01–08)", xticks=positions, yscale="log",
                          scope_note="Primary complete-process boundary; both arms independently verified"))
    figures.append(chart("shared_cohorts", "grouped_bars", "Actual genuine 32-target process costs", "Whole-process wall (s; log scale)",
                          "Three paired preparation cohorts · 32 points per arm/process", [
                              bar_series("IC", [1,2,3], [by_cohort[k,"ic"]["process"]["elapsed_process_wall_s"] for k in (1,2,3)], -0.18),
                              bar_series("Optimized rho", [1,2,3], [by_cohort[k,"rho"]["process"]["elapsed_process_wall_s"] for k in (1,2,3)], 0.18)],
                          xlabel="Preparation cohort", xticks=[1,2,3], yscale="log",
                          scope_note="Primary: actual shared infrastructure and every target solve charged"))
    target_positions = list(range(1, 33))
    figures.append(chart("query_targets", "points", "All 32 target-mean producer query costs", "Producer query wall (s; log scale)",
                          "32 paired target blocks · three shared preparations · unequal validation boundaries", [
                              {"label":"IC target mean", "x":target_positions, "y":[r["ic"]["mean"] for r in rows], "marker":"o"},
                              {"label":"Rho target mean", "x":target_positions, "y":[r["rho"]["mean"] for r in rows], "marker":"x"}],
                          xlabel="Frozen target index", xticks=[1,4,8,12,16,20,24,28,32], yscale="log",
                          figsize=[9.8,4.0],
                          scope_note="Secondary query comparison; not an end-to-end speedup"))
    ordered_resources = [("cold", str(i), arm, cold_pairs[i-1][arm]) for i in range(1, 9) for arm in ("ic","rho")]
    ordered_resources += [("cohort", str(i), arm, by_cohort[i,arm]["process"]) for i in (1,2,3) for arm in ("ic","rho")]
    resource_positions = list(range(1,12))
    figures.append(chart("process_memory", "grouped_bars", "Measured process peak resident memory", "wait4 peak RSS (bytes; log scale)",
                          "Eight cold pairs and three shared pairs · process peaks shown separately", [
                              bar_series("IC", resource_positions, [r[3]["peak_rss_bytes"] for r in ordered_resources if r[2]=="ic"], -0.18),
                              bar_series("Rho", resource_positions, [r[3]["peak_rss_bytes"] for r in ordered_resources if r[2]=="rho"], 0.18)],
                          xlabel="Cold targets C1–C8; genuine shared cohorts B1–B3", xticks=resource_positions, yscale="log",
                          xticklabels=["C1","C2","C3","C4","C5","C6","C7","C8","B1","B2","B3"],
                          figsize=[9.8,4.0],
                          scope_note="Observed RSS, not RLIMIT_AS allocation and not additive across sequential processes"))
    points = model["display_model_points"]
    figures.append(chart("break_even_model", "lines", "Additive setup + per-target model (not an observed crossover)", "Modeled cost (s)",
                          "Native model · target count above 32 is unexecuted extrapolation", [
                              {"label":"IC model", "x":[p["targets"] for p in points], "y":[p["ic_model_s"] for p in points]},
                              {"label":"Rho model", "x":[p["targets"] for p in points], "y":[p["rho_model_s"] for p in points]}],
                          xlabel="Modeled target count N (log scale)", xscale="log", yscale="log",
                          figsize=[9.8,4.5],
                          vertical_markers=[{"x":32,"label":"Actual shared workload N=32"}],
                          scope_note=meta["model_figure_scope"]))
    figures.append(chart("aa_noise", "points", "Identical optimized-rho A/A process pairs", "A/B whole-process wall ratio", "Five same-point, same-seed pairs · independent semantic gate passes", [
                              {"label":"Observed A/B", "x":[1,2,3,4,5], "y":stats["aa"]["A_over_B_wall_ratios"]}],
                          xlabel="Predeclared A/A pair", xticks=[1,2,3,4,5],
                          horizontal_markers=[{"y":1,"label":"Equal process time"}],
                          scope_note="Descriptive noise diagnostic; no confidence band and no correction"))
    (ROOT/"figure_data.json").write_text(json.dumps({"evidence_status":display_status,
         "research_results_recomputed":False, "native_statistics_source":str(args.statistics), "figures":figures}, indent=2)+"\n")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Aggregate the double-large-prime `ic bench` reports into results.json
and a markdown table.

    python3 summarise.py runs/<stamp> [runs/<stamp> ...]

One cell is (degree, columns, large_primes, |D|/|F|); its repeats (one
planted target each) are summed, so "trials per complete relation" is
total trials over total complete relations, not a mean of ratios.
Operation counts are the metric; wall time is a practicality note.
"""
import gzip, json, pathlib, re, sys

here = pathlib.Path(__file__).resolve().parent


def load(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        return json.load(f)


def native(phase, key):
    return phase.get("native", {}).get(key, 0)


def cell_key(row):
    spec = row["spec"]
    op = spec["oracle_params"]
    lp = int(op.get("large_primes", "0"))
    factor = float(op.get("lp_domain_factor", "0")) if lp else 0.0
    columns = int(spec["factor_base_params"]["columns"])
    degree = int(re.search(r"f2m(\d+)", row["instance"]).group(1)) if "f2m" in row["instance"] else None
    return (degree, columns, lp, factor)


def main(run_dirs):
    cells = {}
    sources = []
    for run_dir in run_dirs:
        run_dir = pathlib.Path(run_dir)
        for report in sorted(run_dir.glob("report_*.json*")):
            doc = load(report)
            sources.append(str(report.relative_to(here)) if report.is_relative_to(here) else str(report))
            degree = None
            m = re.search(r"report_dlp_n(\d+)_c(\d+)", report.name)
            if m:
                degree = int(m.group(1))
            for row in doc["rows"]:
                key = cell_key(row)
                if key[0] is None:
                    key = (degree,) + key[1:]
                cell = cells.setdefault(
                    key,
                    {
                        "degree": key[0],
                        "columns": key[1],
                        "large_primes": key[2],
                        "domain_factor": key[3],
                        "instance": row["instance"],
                        "log2_r": row["log2_r"],
                        "targets": 0,
                        "verified_targets": 0,
                        "signed_points": row["factor_base"]["signed_points"],
                        "trials": 0,
                        "complete_relations": 0,
                        "direct_three_sum": 0,
                        "direct_two_sum": 0,
                        "merged_relations": 0,
                        "partials_fed": 0,
                        "partials_one_large": 0,
                        "partials_mixed": 0,
                        "partials_two_large": 0,
                        "partials_stored_at_end": 0,
                        "pivot_reductions": 0,
                        "merge_depth_total": 0,
                        "merge_sources_total": 0,
                        "rejected": 0,
                        "capped": 0,
                        "verification_failures": 0,
                        "verification_scalar_mults": 0,
                        "relation_adds": 0,
                        "relation_doubles": 0,
                        "relation_scalar_mults": 0,
                        "relation_lookups": 0,
                        "relation_canonicalisations": 0,
                        "relation_gae": 0.0,
                        "setup_adds": 0,
                        "setup_gae": 0.0,
                        "classic_table_entries": 0,
                        "lp_table_entries": 0,
                        "lp_table_adds": 0,
                        "domain_orbits": 0,
                        "domain_points": 0,
                        "base_gae": 0.0,
                        "la_rows": 0,
                        "la_rank": 0,
                        "la_row_ops": 0,
                        "la_gae": 0.0,
                        "total_gae": 0.0,
                        "s_values": [],
                        "rho_s": row.get("rho_s"),
                        "filter_wall_seconds": 0.0,
                        "wall_seconds": 0.0,
                    },
                )
                rel = row["decomposition"]["cost"]
                setup = row["decomposition"]["setup"]
                la = row["linear_algebra"]
                cell["targets"] += 1
                cell["verified_targets"] += int(bool(row["verified"]))
                cell["trials"] += row["decomposition"]["targets_tried"]
                cell["complete_relations"] += row["decomposition"]["relations_found"]
                for name, key2 in [
                    ("direct_three_sum", "lp_direct_three_sum"),
                    ("direct_two_sum", "lp_direct_two_sum"),
                    ("merged_relations", "lp_merged_relations"),
                    ("partials_fed", "lp_partials_fed"),
                    ("partials_one_large", "lp_partials_one_large"),
                    ("partials_mixed", "lp_partials_mixed"),
                    ("partials_two_large", "lp_partials_two_large"),
                    ("partials_stored_at_end", "lp_partials_stored"),
                    ("pivot_reductions", "lp_pivot_reductions"),
                    ("merge_depth_total", "lp_merge_depth_total"),
                    ("merge_sources_total", "lp_merge_sources_total"),
                    ("rejected", "lp_rejected"),
                    ("capped", "lp_capped"),
                    ("verification_failures", "lp_verification_failures"),
                    ("verification_scalar_mults", "lp_verification_scalar_mults"),
                    ("relation_lookups", "lookups"),
                    ("relation_canonicalisations", "canonicalisations"),
                ]:
                    cell[name] += native(rel, key2)
                cell["filter_wall_seconds"] += native(rel, "lp_filter_wall_ns") / 1e9
                cell["relation_adds"] += rel["group_ops"]["adds"]
                cell["relation_doubles"] += rel["group_ops"]["doubles"]
                cell["relation_scalar_mults"] += rel["group_ops"]["scalar_mults"]
                cell["relation_gae"] += rel["gae"]
                cell["setup_adds"] += setup["group_ops"]["adds"]
                cell["setup_gae"] += setup["gae"]
                cell["classic_table_entries"] += native(setup, "pair_table_entries")
                cell["lp_table_entries"] += native(setup, "lp_table_entries")
                cell["lp_table_adds"] += native(setup, "lp_table_adds")
                cell["domain_orbits"] = native(setup, "lp_domain_orbits")
                cell["domain_points"] = native(setup, "lp_domain_points")
                cell["base_gae"] += row["factor_base"]["cost"]["gae"]
                cell["la_rows"] += la["rows"]
                cell["la_rank"] += la["rank"]
                cell["la_row_ops"] += la["work"]
                cell["la_gae"] += la["cost"]["gae"]
                cell["total_gae"] += row["total_gae"]
                cell["s_values"].append(row["s"])
                cell["wall_seconds"] += row["wall_seconds"]
    out = []
    for key in sorted(cells):
        c = cells[key]
        t = max(c["targets"], 1)
        rels = max(c["complete_relations"], 1)
        c["trials_per_complete_relation"] = c["trials"] / rels
        c["relation_adds_per_complete_relation"] = (c["relation_adds"] + c["relation_doubles"]) / rels
        c["relation_gae_per_complete_relation"] = c["relation_gae"] / rels
        c["setup_adds_per_target"] = c["setup_adds"] / t
        c["lp_table_adds_per_target"] = c["lp_table_adds"] / t
        c["lp_table_entries_per_target"] = c["lp_table_entries"] / t
        c["classic_table_entries_per_target"] = c["classic_table_entries"] / t
        c["partials_per_complete_relation"] = c["partials_fed"] / rels
        c["merge_yield"] = c["merged_relations"] / c["partials_fed"] if c["partials_fed"] else None
        c["merged_share"] = c["merged_relations"] / rels
        c["mean_merge_depth"] = (
            c["merge_depth_total"] / c["merged_relations"] if c["merged_relations"] else None
        )
        c["mean_merge_sources"] = (
            c["merge_sources_total"] / c["merged_relations"] if c["merged_relations"] else None
        )
        c["la_rows_per_target"] = c["la_rows"] / t
        c["la_row_ops_per_target"] = c["la_row_ops"] / t
        c["total_gae_per_target"] = c["total_gae"] / t
        c["s_mean"] = sum(c["s_values"]) / len(c["s_values"])
        c["s_over_rho"] = c["s_mean"] / c["rho_s"] if c.get("rho_s") else None
        c["domain_over_base"] = c["domain_points"] / c["signed_points"] if c["signed_points"] else None
        c["all_verified"] = c["verified_targets"] == c["targets"]
        out.append(c)
    (here / "results.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "what_this_is": "Double-large-prime collector on the compact-orbit Koblitz pipeline: per cell (degree, columns, large_primes, |D|/|F|), repeats summed. Operation counts; wall time is a practicality note on a shared host.",
                "sources": sources,
                "cells": out,
            },
            indent=1,
        )
        + "\n"
    )
    print(markdown(out))


def fmt(v, digits=3):
    if v is None:
        return "—"
    if isinstance(v, float):
        if v == 0:
            return "0"
        if abs(v) >= 1e5:
            return f"{v:.3g}"
        if abs(v) >= 100:
            return f"{v:.0f}"
        return f"{v:.{digits}g}"
    return str(v)


def markdown(cells):
    lines = [
        "| n | cols | \\|F\\| | LP | \\|D\\|/\\|F\\| | targets | trials | complete rel. | trials / rel. | partials fed | stored at end | merge yield | merged share | rel. adds / rel. | setup adds / target | LP table entries | LA rows / target | row ops / target | total GAE / target | S | S / rho | rejected / capped / replay fail | verified |",
        "|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|",
    ]
    for c in cells:
        lines.append(
            "| {n} | {cols} | {F} | {lp} | {df} | {t} | {trials} | {rels} | {tpr} | {pf} | {stored} | {my} | {ms} | {apr} | {sa} | {lte} | {lar} | {laops} | {tg} | {s} | {sr} | {bad} | {ok} |".format(
                bad=f"{c['rejected']} / {c['capped']} / {c['verification_failures']}",
                n=c["degree"],
                cols=c["columns"],
                F=c["signed_points"],
                lp=c["large_primes"],
                df=fmt(c["domain_over_base"], 3) if c["large_primes"] else "—",
                t=c["targets"],
                trials=c["trials"],
                rels=c["complete_relations"],
                tpr=fmt(c["trials_per_complete_relation"]),
                pf=c["partials_fed"],
                stored=c["partials_stored_at_end"],
                my=fmt(c["merge_yield"]),
                ms=fmt(c["merged_share"]),
                apr=fmt(c["relation_adds_per_complete_relation"]),
                sa=fmt(c["setup_adds_per_target"]),
                lte=fmt(c["lp_table_entries_per_target"]),
                lar=fmt(c["la_rows_per_target"]),
                laops=fmt(c["la_row_ops_per_target"]),
                tg=fmt(c["total_gae_per_target"]),
                s=fmt(c["s_mean"]),
                sr=fmt(c["s_over_rho"]),
                ok="yes" if c["all_verified"] else "NO",
            )
        )
    return "\n".join(lines)


if __name__ == "__main__":
    main(sys.argv[1:] or [p for p in sorted((here / "runs").glob("*")) if p.is_dir()])

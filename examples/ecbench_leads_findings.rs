//! Interpret sealed follow-up counts; never run or alter an experiment.
use serde_json::{json, Value};
use std::{collections::BTreeMap, fmt::Write as _, fs, path::Path};

pub fn main() -> Result<(), String> {
    let root = Path::new("research/ecbench_leads_20261015");
    let summary: Value =
        serde_json::from_slice(&fs::read(root.join("summary.json")).map_err(|e| e.to_string())?)
            .map_err(|e| e.to_string())?;
    let rows = summary["rows"].as_array().ok_or("missing rows")?;
    let mut prime = BTreeMap::<String, Vec<&Value>>::new();
    let mut orbits = BTreeMap::<String, BTreeMap<String, &Value>>::new();
    let mut wide = Vec::new();
    for row in rows {
        let curve = row["curve"].as_str().unwrap().to_string();
        let arm = row["arm"].as_str().unwrap();
        match row["session"].as_str().unwrap() {
            "prime-balance" if arm.contains("mitm2") => prime.entry(curve).or_default().push(row),
            "orbit-ladder" => {
                orbits.entry(curve).or_default().insert(arm.into(), row);
            }
            "wide-n83" => {
                let c = &row["counter_sums"];
                let attempts = row["records"].as_u64().unwrap();
                let b = c["factor_base_signed_points"]
                    .as_u64()
                    .ok_or("missing base support")?
                    / attempts;
                let r = 2_417_851_639_230_796_216_685_689u128;
                let tuples = u128::from(b) * (u128::from(b) + 1) * (u128::from(b) + 2) / 6;
                let ceiling = tuples as f64 / (r - 1) as f64;
                wide.push(json!({"arm":arm,"attempts":attempts,"signed_points_per_attempt":b,
                    "orbit_columns":c["factor_base_columns"].as_u64().unwrap()/attempts,
                    "total_queries":c["trials"],"point_checked_relations":c["point_checked_relations"].as_u64().unwrap_or(0),
                    "aggregate_rank":c["rank"],"triples_with_repetition":tuples.to_string(),
                    "uniform_query_counting_ceiling":ceiling,"union_bound_256_uniform_queries":256.*ceiling,
                    "model_scope":"existing multiset counting argument; a conditional model for uniform a; no success-probability claim for a fixed deterministic hash stream"}));
            }
            _ => {}
        }
    }
    let mut best = Vec::new();
    for (curve, candidates) in prime {
        let chosen = candidates
            .into_iter()
            .filter(|r| r["verified"] == r["records"])
            .min_by(|a, b| {
                a["mean_s"]
                    .as_f64()
                    .unwrap()
                    .total_cmp(&b["mean_s"].as_f64().unwrap())
            })
            .ok_or("no complete prime arm")?;
        best.push(json!({"curve":curve,"arm":chosen["arm"],"targets":chosen["records"],"mean_s":chosen["mean_s"],
            "counted_over_rho":chosen["paired_complete_counted_ratio_candidate_over_rho"],"selection_scope":"descriptive minimum on this panel, not a validated tuning rule"}));
    }
    let mut folded = Vec::new();
    for (curve, arms) in orbits {
        let a = arms["exact-c4"];
        let b = arms["orbit-c4"];
        if a["factor_base_ids"] != b["factor_base_ids"] {
            return Err("direct and folded bases differ".into());
        }
        folded.push(json!({"curve":curve,"direct_mean_s":a["mean_s"],"folded_mean_s":b["mean_s"],
            "direct_over_folded_count_quotient":a["mean_s"].as_f64().unwrap()/b["mean_s"].as_f64().unwrap(),
            "folded_over_rho_count_quotient":b["paired_complete_counted_ratio_candidate_over_rho"],
            "same_factor_base_ids":true,"fully_priced_ratio":null}));
    }
    let held: Vec<_> = rows
        .iter()
        .filter(|r| {
            r["session"] == "orbit-heldout" && !r["arm"].as_str().unwrap().starts_with("rho")
        })
        .collect();
    let c4 = held
        .iter()
        .find(|r| r["arm"] == "orbit-c4")
        .ok_or("heldout c4")?;
    let c8 = held
        .iter()
        .find(|r| r["arm"] == "orbit-c8")
        .ok_or("heldout c8")?;
    let heldout = json!({"curve":c4["curve"],"c4_mean_s":c4["mean_s"],"c8_mean_s":c8["mean_s"],
        "c4_over_c8_count_quotient":c4["mean_s"].as_f64().unwrap()/c8["mean_s"].as_f64().unwrap(),
        "c8_over_rho_count_quotient":c8["paired_complete_counted_ratio_candidate_over_rho"],"fully_priced_ratio":null});
    let attempts: u64 = rows.iter().map(|r| r["records"].as_u64().unwrap()).sum();
    let verified: u64 = rows.iter().map(|r| r["verified"].as_u64().unwrap()).sum();
    let exhausted: u64 = rows
        .iter()
        .map(|r| r["statuses"]["exhausted"].as_u64().unwrap_or(0))
        .sum();
    let errors: u64 = rows
        .iter()
        .map(|r| r["statuses"]["error"].as_u64().unwrap_or(0))
        .sum();
    let rank_rows: Vec<_> = rows
        .iter()
        .filter(|r| {
            matches!(
                r["session"].as_str(),
                Some("orbit-ladder" | "orbit-heldout")
            )
        })
        .collect();
    let rank_total: u64 = rank_rows
        .iter()
        .map(|r| r["records"].as_u64().unwrap())
        .sum();
    let rank_ic: u64 = rank_rows
        .iter()
        .filter(|r| !r["arm"].as_str().unwrap().starts_with("rho"))
        .map(|r| r["verified"].as_u64().unwrap())
        .sum();
    let paragraphs=vec![
        format!("All {attempts} attempts completed: {verified} verified recoveries, {exhausted} exhausted searches and {errors} retained legacy-width error. Of the {rank_total} attempts in the two rank panels, {rank_ic} are IC recoveries and {} are paired rho/A-A controls.",rank_total-rank_ic),
        format!("Both IC arms recovered all four held-out degree 37 public points, eight IC executions in total. Increasing the base from four to eight orbit columns changed mean counted S from {:.6} to {:.6}; the four/eight count quotient is {:.6}. The eight-column/rho count quotient remains {:.6}. These are incomplete-cost diagnostics, not speedups.",c4["mean_s"].as_f64().unwrap(),c8["mean_s"].as_f64().unwrap(),heldout["c4_over_c8_count_quotient"].as_f64().unwrap(),heldout["c8_over_rho_count_quotient"].as_f64().unwrap()),
        "The exact E0 degree 83 method now executes ordinary IC queries at both base sizes, but all four runs reached 256 queries without a point relation or rank gain. This discharges an implementation refusal; it does not pass the N83 target-recovery gate.".into(),
        "The zero-hit N83 result is weak evidence about mathematics: under the existing uniform-query multiset-counting model, these tiny bases cover a vanishing fraction of the subgroup. The model is not a calibrated probability for the fixed hash stream. A scalable algebraic base/witness improvement is still required.".into(),
        format!("Prime base balance changes across models: {}. These are descriptive minima on this panel, not a validated tuning rule or exponent.",best.iter().map(|r|format!("{}: {}",r["curve"].as_str().unwrap(),r["arm"].as_str().unwrap())).collect::<Vec<_>>().join("; ")),
        "Next discriminating work: price the new method's native operations and compare the column-count choice on further unseen points; retain all relative phases in any proposed algebraic correspondence; demand growing base dimension and successful exact E0 N83 recovery before any N131 inference. Recent main also has versioned full-price methods; their host-specific calibration is a follow-up input, not a price silently imported into this AMD run.".into(),
    ];
    let value = json!({"schema":"ecbench-leads-findings/v1","prime_descriptive_minima":best,"folding":folded,"heldout":heldout,"wide":wide,
        "paragraphs":paragraphs,"new_theorem":false,"full_cost_ratio":null,"wall_speedup":null});
    fs::write(
        root.join("findings.json"),
        format!("{}\n", serde_json::to_string_pretty(&value).unwrap()),
    )
    .map_err(|e| e.to_string())?;
    let mut report = String::from("# What the follow-up establishes\n\n");
    for p in value["paragraphs"].as_array().unwrap() {
        writeln!(report, "{}\n", p.as_str().unwrap()).unwrap();
    }
    report.push_str("## Prime balance on fresh targets\n\n| Exact curve | Descriptive minimum | Mean counted S | Counted / rho | Targets |\n|---|---|---:|---:|---:|\n");
    for b in value["prime_descriptive_minima"].as_array().unwrap() {
        writeln!(
            report,
            "|{}|{}|{:.6}|{:.6}|{}|",
            b["curve"].as_str().unwrap(),
            b["arm"].as_str().unwrap(),
            b["mean_s"].as_f64().unwrap(),
            b["counted_over_rho"].as_f64().unwrap(),
            b["targets"]
        )
        .unwrap();
    }
    report.push_str("\n## Direct versus folded tables, identical bases\n\n| Exact curve | Direct S | Folded S | Direct/folded counted quotient | Folded/rho counted quotient |\n|---|---:|---:|---:|---:|\n");
    for f in value["folding"].as_array().unwrap() {
        writeln!(
            report,
            "|{}|{:.6}|{:.6}|{:.6}|{:.6}|",
            f["curve"].as_str().unwrap(),
            f["direct_mean_s"].as_f64().unwrap(),
            f["folded_mean_s"].as_f64().unwrap(),
            f["direct_over_folded_count_quotient"].as_f64().unwrap(),
            f["folded_over_rho_count_quotient"].as_f64().unwrap()
        )
        .unwrap();
    }
    report.push_str("\n## N83 coverage model, separated from observations\n\nThe already-used multiset count C(B+2,3)/(r-1) bounds the number of reachable sums per uniformly randomized nonzero a, conditional on b. Distinct tuples can collide. The following union bounds are model calculations, not measured rates or calibrated probabilities for the frozen deterministic hash stream.\n\n| Columns | Signed points | Queries across two attempts | Checked relations | 256-query model union bound |\n|---|---:|---:|---:|---:|\n");
    for w in value["wide"].as_array().unwrap() {
        writeln!(
            report,
            "|{}|{}|{}|{}|{:.8e}|",
            w["orbit_columns"],
            w["signed_points_per_attempt"],
            w["total_queries"],
            w["point_checked_relations"],
            w["union_bound_256_uniform_queries"].as_f64().unwrap()
        )
        .unwrap();
    }
    report.push_str("\nSources: [sealed-record summary](summary.json), [per-target rows](per-target.json), [counting-ceiling implementation](../../src/cryptanalysis/ic_boundary.rs), [Galbraith et al., On Index Calculus Algorithms for Subfield Curves](https://sacworkshop.org/SAC20/files/preproceedings/18-IndexCalculus.pdf), [Caminata and Gorla, Solving degree, last fall degree, and related invariants](https://arxiv.org/abs/2112.05579), [recent main full-price implementation](https://github.com/aburan28/crypto/blob/fd8bcacf4cc097cc0dcb55ba068ad94a3d3202ae/src/cryptanalysis/ecbench/full_price.rs), and the [parent29-reference literature review](../ecbench_open_questions_20261014/literature/REPORT.md). No new mathematical theorem is claimed; the counting formula and its assumptions are inherited, and no fixed-stream probability bound is asserted.\n");
    fs::write(root.join("FINDINGS.md"), report).map_err(|e| e.to_string())?;
    let path = root.join("REPORT.md");
    let report = fs::read_to_string(&path).map_err(|e| e.to_string())?;
    if !report.contains("[Measured findings and next discriminating steps](FINDINGS.md)") {
        fs::write(path,report.replacen("# Pursuing the ecbench leads\n\n",
            "# Pursuing the ecbench leads\n\n[Measured findings and next discriminating steps](FINDINGS.md); [v1 reader correction and preserved failure](AUDIT-READER-AMENDMENT.md).\n\n",1)).map_err(|e|e.to_string())?;
    }
    Ok(())
}

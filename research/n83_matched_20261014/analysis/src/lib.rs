//! Prospective, native accounting for the 32-target / three-cohort matched study.
//! Ratios are IC/rho. All repeated draws remain attached to their target/cohort.
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};

pub const TARGETS: usize = 32;
pub const COHORTS: usize = 3;
pub const COLD_TARGETS: usize = 8;
pub const AA_PAIRS: usize = 5;
pub const DRAWS: usize = 20_000;
pub const SEED: u64 = 202610146799;

#[derive(Clone, Copy)]
struct SplitMix64(u64);
impl SplitMix64 {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^ (z >> 31)
    }
    fn index(&mut self, n: usize) -> usize {
        let modulus = n as u64;
        let threshold = modulus.wrapping_neg() % modulus;
        loop {
            let value = self.next();
            if value >= threshold {
                return (value % modulus) as usize;
            }
        }
    }
}

fn positive(v: &Value) -> Option<f64> {
    v.as_f64().filter(|n| n.is_finite() && *n > 0.)
}
fn nonnegative(v: &Value) -> Option<f64> {
    v.as_f64().filter(|n| n.is_finite() && *n >= 0.)
}
fn mean(v: &[f64]) -> Option<f64> {
    (!v.is_empty()).then(|| v.iter().sum::<f64>() / v.len() as f64)
}
fn ratio(a: f64, b: f64) -> Option<f64> {
    let r = a / b;
    (a >= 0. && b > 0. && r.is_finite()).then_some(r)
}

/// Linear interpolation between adjacent order statistics, including midpoint median.
pub fn quantile(sorted: &[f64], p: f64) -> Option<f64> {
    if sorted.is_empty() || !(0. ..=1.).contains(&p) {
        return None;
    }
    let index = p * (sorted.len() - 1) as f64;
    let lo = index.floor() as usize;
    let hi = index.ceil() as usize;
    Some(sorted[lo] + (sorted[hi] - sorted[lo]) * (index - lo as f64))
}

fn descriptive(v: &[f64]) -> Value {
    let mut sorted = v.to_vec();
    sorted.sort_by(f64::total_cmp);
    json!({"n":v.len(),"mean":mean(v),"median":quantile(&sorted,0.5),"p90":quantile(&sorted,0.9),"p95":quantile(&sorted,0.95),"min":sorted.first(),"max":sorted.last(),
        "quantile_rule":"linear interpolation on p*(n-1); median averages middle two for even n","tail_scope":"descriptive observed order statistics; no extreme-tail or population-tail guarantee"})
}

fn status_ok(row: &Value, batch: bool) -> bool {
    let verified = if batch {
        row["all_verified"].as_bool()
    } else {
        row["verified"].as_bool()
    };
    let acceptable_status = matches!(
        row["status"].as_str(),
        Some(
            "VERIFIED"
                | "VERIFIED_SOLVE"
                | "VERIFIED_SOLVED"
                | "COMPLETED_VERIFIED"
                | "COMPLETED"
                | "SUCCESS"
        )
    );
    verified == Some(true)
        && acceptable_status
        && row["wrong_answer"].as_bool() != Some(true)
        && row["timeout"].as_bool() != Some(true)
        && row.get("exit_code").map_or(true, |v| v.as_i64() == Some(0))
}

fn comparison(pairs: &[(f64, f64)], valid: bool, seed: u64, unit: &str, scope: &str) -> Value {
    if !valid || pairs.is_empty() {
        return json!({"valid":false,"ic_over_rho":null,"interval_95":null,"reason":"Required original matrix is incomplete, duplicated, unverified, failed, or invalid; no completed-subset primary estimate.","bootstrap_unit":unit,"scope":scope});
    }
    let estimate = ratio(
        pairs.iter().map(|x| x.0).sum(),
        pairs.iter().map(|x| x.1).sum(),
    );
    let mut rng = SplitMix64(seed);
    let mut draws = Vec::with_capacity(DRAWS);
    for _ in 0..DRAWS {
        let (mut a, mut b) = (0., 0.);
        for _ in pairs {
            let p = pairs[rng.index(pairs.len())];
            a += p.0;
            b += p.1;
        }
        if let Some(r) = ratio(a, b) {
            draws.push(r);
        }
    }
    draws.sort_by(f64::total_cmp);
    json!({"valid":estimate.is_some() && draws.len()==DRAWS,"ic_over_rho":estimate,"interval_95":[quantile(&draws,0.025),quantile(&draws,0.975)],
        "arithmetic_mean_ic_s":mean(&pairs.iter().map(|x|x.0).collect::<Vec<_>>()),"arithmetic_mean_rho_s":mean(&pairs.iter().map(|x|x.1).collect::<Vec<_>>()),
        "estimand":"ratio of arithmetic mean costs, IC/rho; below one favors IC","bootstrap_draws":DRAWS,"bootstrap_seed":seed,"bootstrap_unit":unit,"method":"paired percentile bootstrap, deterministic SplitMix64 with unbiased index draws","scope":scope})
}

fn strings(v: &Value) -> Option<Vec<String>> {
    v.as_array()?
        .iter()
        .map(|v| Some(v.as_str()?.to_owned()))
        .collect()
}
fn digest(v: &Value) -> bool {
    v.as_str()
        .is_some_and(|s| s.len() == 64 && s.chars().all(|c| c.is_ascii_hexdigit()))
}

fn process_summaries(rows: &[&Value]) -> Value {
    let numbers = |key: &str| {
        rows.iter()
            .filter_map(|r| nonnegative(&r[key]))
            .collect::<Vec<_>>()
    };
    json!({"rows":rows.len(),"verified_rows":rows.iter().filter(|r|status_ok(r,true)).count(),"wall_s":descriptive(&numbers("elapsed_process_wall_s")),
        "cpu_total_s":descriptive(&numbers("cpu_total_s")),"cpu_user_s":descriptive(&numbers("cpu_user_s")),"cpu_system_s":descriptive(&numbers("cpu_system_s")),
        "peak_rss_bytes":descriptive(&numbers("peak_rss_bytes")),"preparation_ms":descriptive(&numbers("producer_preparation_ms")),
        "scope":"all supplied observed rows, including consumed failure costs; no failed row removed; descriptive costs are not successful-completion mean estimates"})
}

fn process_pair_costs<'a>(
    rows: &[&'a Value],
    key: &str,
    keys: &[String],
    count: usize,
) -> (Vec<(f64, f64)>, Vec<String>) {
    let mut issues = Vec::new();
    let mut slots = BTreeMap::<(String, String), &Value>::new();
    for (index, row) in rows.iter().enumerate() {
        let id = if key == "cohort_id" {
            row[key].as_u64().map(|n| n.to_string())
        } else {
            row[key].as_str().map(str::to_owned)
        };
        let arm = row["arm"].as_str();
        let Some(id) = id else {
            issues.push(format!("process row{index}:missing {key}"));
            continue;
        };
        let Some(arm @ ("ic" | "rho")) = arm else {
            issues.push(format!("process row{index}:unknown arm"));
            continue;
        };
        if !keys.contains(&id) {
            issues.push(format!("process row{index}:unexpected {key}={id}"));
        }
        if slots.insert((id.clone(), arm.to_owned()), row).is_some() {
            issues.push(format!("duplicate process {id}/{arm}"));
        }
        if row["target_count"].as_u64() != Some(count as u64) {
            issues.push(format!("{id}/{arm}:incorrect target_count"));
        }
        if !status_ok(row, true) {
            issues.push(format!("{id}/{arm}:not verified uncensored completion"));
        }
        if positive(&row["elapsed_process_wall_s"]).is_none() {
            issues.push(format!("{id}/{arm}:invalid wall time"));
        }
    }
    if rows.len() != keys.len() * 2 {
        issues.push(format!(
            "expected {} process rows; observed {}",
            keys.len() * 2,
            rows.len()
        ));
    }
    let mut pairs = Vec::new();
    for id in keys {
        match (
            slots.get(&(id.clone(), "ic".into())),
            slots.get(&(id.clone(), "rho".into())),
        ) {
            (Some(a), Some(b)) => {
                if let (Some(a), Some(b)) = (
                    positive(&a["elapsed_process_wall_s"]),
                    positive(&b["elapsed_process_wall_s"]),
                ) {
                    pairs.push((a, b));
                }
            }
            _ => issues.push(format!("{id}:missing paired process arm")),
        }
    }
    (pairs, issues)
}

/// All rows are retained by the input receipt. No censoring imputation, fastest
/// seed selection, new tuning, or post-hoc target removal occurs here.
pub fn summarize(input: &Value) -> Value {
    let processes: Vec<&Value> = input["process_rows"]
        .as_array()
        .map(|a| a.iter().collect())
        .unwrap_or_default();
    let targets: Vec<&Value> = input["target_rows"]
        .as_array()
        .map(|a| a.iter().collect())
        .unwrap_or_default();
    let mut global = Vec::<String>::new();
    if input["schema"] != "n83_matched_statistics_input_v1" {
        global.push("input schema mismatch".into());
    }
    let audit = &input["independent_verification"];
    if audit["status"] != "PASS"
        || audit["checked_target_rows"].as_u64() != Some(218)
        || audit["invalid_rows"].as_u64() != Some(0)
        || audit["checked_aa_semantic_pairs"].as_u64() != Some(AA_PAIRS as u64)
        || audit["aa_semantic_gate"] != true
        || !digest(&audit["campaign_manifest_sha256"])
    {
        global.push(
            "require independently checked218-row and5AA-semantic-pair PASS receipt bound to campaign manifest".into(),
        );
    }
    if input["metadata"]["independently_verified"] != true
        || !digest(&input["metadata"]["protocol_sha256"])
        || !digest(&input["metadata"]["campaign_manifest_sha256"])
        || !digest(&input["metadata"]["independent_verification_receipt_sha256"])
        || audit["campaign_manifest_sha256"] != input["metadata"]["campaign_manifest_sha256"]
        || audit["protocol_sha256"] != input["metadata"]["protocol_sha256"]
    {
        global.push("independent gate must be true and bound to matching protocol/manifest plus actual receipt digest".into());
    }
    let ids = strings(&input["metadata"]["target_ids"]).unwrap_or_default();
    let cold_ids = strings(&input["metadata"]["cold_target_ids"]).unwrap_or_default();
    if ids.len() != TARGETS || ids.iter().collect::<BTreeSet<_>>().len() != TARGETS {
        global.push("require32 unique frozen target_ids".into());
    }
    if cold_ids.len() != COLD_TARGETS
        || cold_ids != ids.iter().take(COLD_TARGETS).cloned().collect::<Vec<_>>()
    {
        global.push("require first8 frozen cold_target_ids".into());
    }
    let batch_rows = processes
        .iter()
        .copied()
        .filter(|r| r["kind"] == "cohort")
        .collect::<Vec<_>>();
    let cold_rows = processes
        .iter()
        .copied()
        .filter(|r| r["kind"] == "cold")
        .collect::<Vec<_>>();
    let aa_rows = processes
        .iter()
        .copied()
        .filter(|r| r["kind"] == "aa")
        .collect::<Vec<_>>();
    if batch_rows.len() + cold_rows.len() + aa_rows.len() != processes.len() {
        global.push("unknown process kind".into());
    }
    let cohort_ids = (1..=COHORTS).map(|n| n.to_string()).collect::<Vec<_>>();
    let (batch_pairs, batch_issues) =
        process_pair_costs(&batch_rows, "cohort_id", &cohort_ids, TARGETS);
    let (cold_pairs, cold_issues) = process_pair_costs(&cold_rows, "target_id", &cold_ids, 1);
    let (mut aa, aa_valid) = summarize_aa(&aa_rows);
    aa["row_matrix_valid"] = json!(aa_valid);
    aa["independent_semantic_pairs_checked"] = audit["checked_aa_semantic_pairs"].clone();
    aa["valid"] = json!(aa_valid && global.is_empty());

    let mut query_issues = Vec::new();
    let mut slots = BTreeMap::<(String, usize, String), f64>::new();
    for (index, row) in targets.iter().enumerate() {
        let (Some(id), Some(cohort), Some(arm @ ("ic" | "rho"))) = (
            row["target_id"].as_str(),
            row["cohort_id"].as_u64(),
            row["arm"].as_str(),
        ) else {
            query_issues.push(format!("query row{index}:missing id/cohort/arm"));
            continue;
        };
        if !ids.iter().any(|x| x == id) || !(1..=COHORTS as u64).contains(&cohort) {
            query_issues.push(format!("query row{index}:outside frozen matrix"));
        }
        if !status_ok(row, false) {
            query_issues.push(format!(
                "{id}/{cohort}/{arm}:not verified uncensored completion"
            ));
        }
        let Some(time) = positive(&row["online_wall_s"]) else {
            query_issues.push(format!("{id}/{cohort}/{arm}:invalid time"));
            continue;
        };
        if slots
            .insert((id.into(), cohort as usize, arm.into()), time)
            .is_some()
        {
            query_issues.push(format!("duplicate query {id}/{cohort}/{arm}"));
        }
    }
    if targets.len() != TARGETS * COHORTS * 2 {
        query_issues.push(format!(
            "require{} query rows; observed{}",
            TARGETS * COHORTS * 2,
            targets.len()
        ));
    }
    let mut matrix = Vec::<[(f64, f64); COHORTS]>::new();
    let mut per_target = Vec::new();
    for id in &ids {
        let mut paired = [(0., 0.); COHORTS];
        let mut complete = true;
        for (c, slot) in paired.iter_mut().enumerate() {
            match (
                slots.get(&(id.clone(), c + 1, "ic".into())),
                slots.get(&(id.clone(), c + 1, "rho".into())),
            ) {
                (Some(a), Some(b)) => *slot = (*a, *b),
                _ => {
                    complete = false;
                    query_issues.push(format!("{id}/{}:missing paired target arm", c + 1));
                }
            }
        }
        let a = paired
            .iter()
            .map(|p| p.0)
            .filter(|v| *v > 0.)
            .collect::<Vec<_>>();
        let b = paired
            .iter()
            .map(|p| p.1)
            .filter(|v| *v > 0.)
            .collect::<Vec<_>>();
        per_target.push(json!({"target_id":id,"complete":complete,"ic_online_s":a,"rho_online_s":b,"ic":descriptive(&a),"rho":descriptive(&b),"paired_mean_ic_over_rho":if complete {ratio(mean(&a).unwrap(),mean(&b).unwrap())} else {None}}));
        matrix.push(paired);
    }
    let main_valid = global.is_empty()
        && batch_issues.is_empty()
        && query_issues.is_empty()
        && cold_issues.is_empty();
    let batch_valid = global.is_empty() && batch_issues.is_empty() && query_issues.is_empty();
    let cold_valid = global.is_empty() && cold_issues.is_empty();
    let query_valid = global.is_empty() && query_issues.is_empty() && batch_issues.is_empty();
    let target_pairs = matrix
        .iter()
        .map(|r| {
            (
                r.iter().map(|x| x.0).sum::<f64>() / COHORTS as f64,
                r.iter().map(|x| x.1).sum::<f64>() / COHORTS as f64,
            )
        })
        .collect::<Vec<_>>();
    let online = comparison(&target_pairs,query_valid,SEED+2,"32 paired target blocks, retaining each target's three cohort observations", "Conditional on these three preparation/cohort epochs, retained K600, frozen policies and this host; repeated IC queries are timing replication, not independent IC algorithms or rank seeds. Producer online boundaries differ: stock IC excludes input validation/serialization while rho online includes Q validation. Secondary diagnostic, not the end-to-end ranking.");
    let crossed = crossed_online(&matrix, query_valid);
    let geometric = if query_valid {
        Some((target_pairs.iter().map(|(a, b)| (a / b).ln()).sum::<f64>() / TARGETS as f64).exp())
    } else {
        None
    };
    let model = break_even(&batch_rows, &matrix, query_valid && batch_valid);
    let metrics = ["ic","rho"].iter().map(|arm| {
        let cohort= batch_rows.iter().copied().filter(|r|r["arm"]==*arm).collect::<Vec<_>>();
        let cold= cold_rows.iter().copied().filter(|r|r["arm"]==*arm).collect::<Vec<_>>();
        let mut all=cohort.clone(); all.extend_from_slice(&cold);
        let pooled=targets.iter().filter(|r|r["arm"]==*arm).filter_map(|r|positive(&r["online_wall_s"])).collect::<Vec<_>>();
        let means=target_pairs.iter().map(|p|if *arm=="ic" {p.0} else {p.1}).collect::<Vec<_>>();
        ((*arm).to_owned(),json!({"prepared_cohort_process":process_summaries(&cohort),"actual_cold_singleton_process":process_summaries(&cold),"all11_preparations":process_summaries(&all),"pooled96_online_draws_s":descriptive(&pooled),"32_target_three_repeat_online_means_s":if query_valid {descriptive(&means)} else {Value::Null}}))
    }).collect::<BTreeMap<_,_>>();
    json!({"schema":"n83_matched_statistics_v1","ratio_direction":"IC/rho; less than1 favors IC","all_required_rows_valid":main_valid && aa_valid,"overall_primary_winner_eligible":main_valid,
        "observed_counts":{"process_rows":processes.len(),"cohort_process_rows":batch_rows.len(),"cold_process_rows":cold_rows.len(),"aa_process_rows":aa_rows.len(),"cohort_target_rows":targets.len()},
        "requirements":{"targets":TARGETS,"cohorts":COHORTS,"cold_target_pairs":COLD_TARGETS,"aa_pairs":AA_PAIRS,"expected_process_rows":COHORTS*2+COLD_TARGETS*2+AA_PAIRS*2,"expected_query_rows":TARGETS*COHORTS*2},
        "issues":{"global":global,"batch":batch_issues,"cold":cold_issues,"query":query_issues},
        "primary_prepared_batch_process_wall":comparison(&batch_pairs,batch_valid,SEED,"3 actual paired whole-process cohorts, each32 targets","All startup, reusable setup, target work, built-in verification and output charged once per actual32-target process. Only three setup/timing clusters, fixed IC rank29; percentile interval has very limited preparation-variation resolution and no independent-host scope."),
        "primary_actual_cold_singleton_process_wall":comparison(&cold_pairs,cold_valid,SEED+1,"8 actual paired fresh-process target blocks","Fresh empty algorithm setup and target state for each one-target child; retained factor-base file remains an input. These are eight measured cold targets, not32 or96 cold observations. OS caches and original retained-base construction are outside this boundary."),
        "secondary_prepared_online_target_means":online,"secondary_crossed_online_sensitivity":crossed,"secondary_geometric_mean_paired_target_online_ratio":{"ic_over_rho":geometric,"scope":"descriptive geometric mean of32 ratios after within-target arithmetic averaging of three repeats; not the practical end-to-end headline"},
        "modeled_break_even":model,"costs":metrics,"per_target":per_target,"aa":aa,
        "process_rows":input["process_rows"],"target_rows":input["target_rows"],"auxiliary_target_rows":input["auxiliary_target_rows"],
        "display_cohorts":display_cohorts(&batch_rows,&targets),"display_cold_pairs":display_cold_pairs(&cold_rows,&cold_ids),
        "input_metadata":input["metadata"],"independent_verification":audit,"external_verification_accounting":input.get("external_verification_accounting").cloned().unwrap_or(json!({"scope":"independent audit outside solver intervals; require separately measured per-arm checker cost before any audit-inclusive end-to-end number","available":false})),
        "failure_policy":"All original attempts remain in input and descriptive consumed-cost summaries. No zero/missing imputation or post-hoc completed-subset primary. Overall winner requires all208planned main answers verified. Each narrower predeclared batch/cold comparison still requires its complete original matrix and the independent full218-row audit; a retained8-target cold metric does not imply32-target study completion. Failure at finite extraction exhaustion is not a mathematical no-decomposition proof.",
        "limitations":["One shared VM; physical-host neighbors and independent-host replication are not established.","Conditional target law excludes retained-base and prior exposed signed-Frobenius orbits; targets are sampled without replacement.","Three shared-setup repeats are not96 independent setups. IC rank29 is fixed, so setup seed-distribution inference is absent.","P90/p95/max at32 targets are observed descriptive order statistics, not precise extreme-tail estimates.","No common operation-equivalent cost, global-optimum claim or asymptotic conclusion."]})
}

fn display_cohorts(processes: &[&Value], targets: &[&Value]) -> Vec<Value> {
    let mut output = Vec::new();
    for cohort in 1..=COHORTS {
        for arm in ["ic", "rho"] {
            let matching = processes
                .iter()
                .filter(|r| r["cohort_id"].as_u64() == Some(cohort as u64) && r["arm"] == arm)
                .collect::<Vec<_>>();
            let process = if matching.len() == 1 {
                *matching[0]
            } else {
                &Value::Null
            };
            let prep = nonnegative(&process["producer_preparation_ms"]).map(|p| p / 1000.);
            let mut sum = 0.;
            let mut prefixes = Vec::new();
            for row in targets
                .iter()
                .filter(|r| r["cohort_id"].as_u64() == Some(cohort as u64) && r["arm"] == arm)
            {
                let time = positive(&row["online_wall_s"]);
                if let Some(t) = time {
                    sum += t;
                }
                prefixes.push(json!({"targets_so_far":prefixes.len()+1,"target_id":row["target_id"],"online_wall_s":time,"cumulative_online_s":sum,"producer_accounted_setup_plus_queries_s":prep.map(|p|p+sum),"verified":row["verified"]}));
            }
            output.push(json!({"cohort_id":cohort,"arm":arm,"process":process,"preparation_s":prep,"summed_online_s":sum,"process_minus_preparation_minus_online_s":positive(&process["elapsed_process_wall_s"]).zip(prep).map(|(wall,prep)|wall-prep-sum),"prefixes":prefixes,"prefix_scope":"Derived producer-accounted prefix using target rows in their supplied producer order; not actual observed prefix wall. Excludes validation/output/launch/residual costs outside online clocks. Whole-process cost is measured only at actual32-target cohort size."}));
        }
    }
    output
}

fn display_cold_pairs(rows: &[&Value], ids: &[String]) -> Vec<Value> {
    ids.iter().map(|id|{
        let find=|arm:&str|rows.iter().find(|r|r["target_id"]==*id && r["arm"]==arm).copied().cloned().unwrap_or(Value::Null);
        let ic=find("ic");
        let rho=find("rho");
        json!({"target_id":id,"ic":ic,"rho":rho,"process_wall_ic_over_rho":positive(&ic["elapsed_process_wall_s"]).zip(positive(&rho["elapsed_process_wall_s"])).and_then(|(a,b)|ratio(a,b)),"scope":"Actual paired singleton process costs; producer and phase details preserved inside raw process rows. Eight target pairs only."})
    }).collect()
}

fn crossed_online(matrix: &[[(f64, f64); COHORTS]], valid: bool) -> Value {
    if !valid {
        return json!({"valid":false,"ic_over_rho":null,"interval_95":null});
    }
    let mut rng = SplitMix64(SEED + 3);
    let mut draws = Vec::with_capacity(DRAWS);
    for _ in 0..DRAWS {
        let cohorts = (0..COHORTS).map(|_| rng.index(COHORTS)).collect::<Vec<_>>();
        let (mut a, mut b) = (0., 0.);
        for _ in 0..matrix.len() {
            let row = &matrix[rng.index(matrix.len())];
            for c in &cohorts {
                a += row[*c].0;
                b += row[*c].1;
            }
        }
        draws.push(a / b);
    }
    draws.sort_by(f64::total_cmp);
    json!({"valid":true,"ic_over_rho":matrix.iter().flat_map(|r|r.iter()).map(|p|p.0).sum::<f64>()/matrix.iter().flat_map(|r|r.iter()).map(|p|p.1).sum::<f64>(),"interval_95":[quantile(&draws,0.025),quantile(&draws,0.975)],"bootstrap_draws":DRAWS,"bootstrap_seed":SEED+3,"bootstrap_unit":"independently resample32 paired target blocks and3 paired cohort blocks, retaining the crossed cells","scope":"sensitivity to target and these three shared preparation/timing epochs; three cohort levels provide limited resolution; no host-general confidence interval"})
}

#[derive(Debug, PartialEq)]
struct ModelOutcome {
    first_n: Option<f64>,
    eventual_n: Option<f64>,
    last_n: Option<f64>,
    kind: &'static str,
}

fn threshold(prep_ic: f64, prep_rho: f64, query_ic: f64, query_rho: f64) -> ModelOutcome {
    let slope = query_rho - query_ic;
    let intercept = prep_ic - prep_rho;
    if slope > 0. {
        let n = (intercept / slope).ceil().max(1.);
        ModelOutcome {
            first_n: Some(n),
            eventual_n: Some(n),
            last_n: None,
            kind: "finite_eventual_ic_advantage",
        }
    } else if slope == 0. && intercept <= 0. {
        ModelOutcome {
            first_n: Some(1.),
            eventual_n: Some(1.),
            last_n: None,
            kind: "ic_at_most_rho_all_positive_n_equal_per_target_cost",
        }
    } else if slope < 0. && intercept - slope <= 0. {
        ModelOutcome {
            first_n: Some(1.),
            eventual_n: None,
            last_n: Some((intercept / slope).floor()),
            kind: "ic_initial_setup_advantage_lost_as_batch_grows",
        }
    } else {
        ModelOutcome {
            first_n: None,
            eventual_n: None,
            last_n: None,
            kind: "no_ic_advantage_for_positive_n_in_linear_model",
        }
    }
}

fn break_even(rows: &[&Value], matrix: &[[(f64, f64); COHORTS]], valid: bool) -> Value {
    let mut preps = [(0., 0.); COHORTS];
    let mut jumps = [0.; COHORTS];
    let mut prep_valid = valid;
    for cohort in 1..=COHORTS {
        for arm in ["ic", "rho"] {
            let list = rows
                .iter()
                .filter(|r| r["cohort_id"].as_u64() == Some(cohort as u64) && r["arm"] == arm)
                .collect::<Vec<_>>();
            let value = if list.len() == 1 {
                positive(&list[0]["producer_preparation_ms"]).map(|p| p / 1000.)
            } else {
                None
            };
            if let Some(p) = value {
                if arm == "ic" {
                    preps[cohort - 1].0 = p;
                } else {
                    let per_target =
                        positive(&list[0]["producer_per_target_preparation_ms"]).map(|j| j / 1000.);
                    if let Some(j) = per_target {
                        let shared = p - TARGETS as f64 * j;
                        if shared > 0. {
                            preps[cohort - 1].1 = shared;
                            jumps[cohort - 1] = j;
                        } else {
                            prep_valid = false;
                        }
                        if let Some(recorded) =
                            nonnegative(&list[0]["producer_shared_preparation_ms"])
                        {
                            if (recorded / 1000. - shared).abs() > 1e-6 {
                                prep_valid = false;
                            }
                        }
                    } else {
                        prep_valid = false;
                    }
                }
            } else {
                prep_valid = false;
            }
        }
    }
    if !prep_valid {
        return json!({"valid":false,"first_n_with_model_ic_at_most_rho":null,"reason":"Complete matching preparation/query clocks and rho shared versus32 seed-specific jump preparation required; no extrapolation from incomplete or misclassified costs."});
    }
    let a = preps.iter().map(|x| x.0).sum::<f64>() / COHORTS as f64;
    let b = preps.iter().map(|x| x.1).sum::<f64>() / COHORTS as f64;
    let c = matrix
        .iter()
        .flat_map(|r| r.iter())
        .map(|x| x.0)
        .sum::<f64>()
        / (matrix.len() * COHORTS) as f64;
    let online_d = matrix
        .iter()
        .flat_map(|r| r.iter())
        .map(|x| x.1)
        .sum::<f64>()
        / (matrix.len() * COHORTS) as f64;
    let jump_mean = jumps.iter().sum::<f64>() / COHORTS as f64;
    let d = online_d + jump_mean;
    let point = threshold(a, b, c, d);
    let mut rng = SplitMix64(SEED + 4);
    let mut finite = Vec::new();
    let mut nonfinite = 0;
    for _ in 0..DRAWS {
        let cohorts = (0..COHORTS).map(|_| rng.index(COHORTS)).collect::<Vec<_>>();
        let prep_a = cohorts.iter().map(|i| preps[*i].0).sum::<f64>() / COHORTS as f64;
        let prep_b = cohorts.iter().map(|i| preps[*i].1).sum::<f64>() / COHORTS as f64;
        let (mut qa, mut qb) = (0., 0.);
        for _ in 0..matrix.len() {
            let row = &matrix[rng.index(matrix.len())];
            for i in &cohorts {
                qa += row[*i].0;
                qb += row[*i].1 + jumps[*i];
            }
        }
        if let Some(n) = threshold(
            prep_a,
            prep_b,
            qa / (matrix.len() * COHORTS) as f64,
            qb / (matrix.len() * COHORTS) as f64,
        )
        .eventual_n
        {
            finite.push(n);
        } else {
            nonfinite += 1;
        }
    }
    finite.sort_by(f64::total_cmp);
    let extended_quantile = |p: f64| -> Value {
        let index = (p * (DRAWS - 1) as f64).ceil() as usize;
        finite
            .get(index)
            .copied()
            .map(|n| json!({"kind":"finite","targets":n}))
            .unwrap_or(json!({"kind":"no_eventual_ic_advantage_in_linear_model","targets":null}))
    };
    json!({"valid":true,"model":"P_IC + N*mean_online_IC versus P_rho_shared + N*(mean_online_rho + mean_seed_specific_G_jump_preparation_rho)","first_n_with_model_ic_at_most_rho":point.first_n,"first_n_with_eventual_model_ic_at_most_rho":point.eventual_n,"last_n_with_initial_model_ic_advantage":point.last_n,"point_kind":point.kind,
        "preparation_ic_s":a,"preparation_rho_shared_s":b,"mean_online_ic_s":c,"mean_online_rho_s":online_d,"rho_seed_specific_jump_per_target_s":jump_mean,"rho_per_target_online_plus_jump_s":d,"rho_per_target_minus_ic_s":d-c,
        "bootstrap_draws":DRAWS,"bootstrap_seed":SEED+4,"bootstrap_finite_draws":finite.len(),"bootstrap_nonfinite_draws":nonfinite,"bootstrap_nonfinite_fraction":nonfinite as f64/DRAWS as f64,"extended_interval_95":[extended_quantile(0.025),extended_quantile(0.975)],
        "finite_only_descriptive_quantiles":{"median":quantile(&finite,0.5),"p025":quantile(&finite,0.025),"p975":quantile(&finite,0.975),"scope":"conditional on finite draws; not an unconditional confidence interval; nonfinite mass retained above"},
        "measured_batch_size":TARGETS,"point_outside_measured_batch":point.first_n.map(|n|n>TARGETS as f64),
        "display_model_points":([1,2,4,8,16,32,64,128,256,512,1024,2048].iter().map(|n|json!({"targets":n,"ic_model_s":a+*n as f64*c,"rho_model_s":b+*n as f64*d,"ic_over_rho":(a+*n as f64*c)/(b+*n as f64*d),"measured_full_process_batch_size":*n==TARGETS,"extrapolated_beyond_measured_batch":*n>TARGETS})).collect::<Vec<_>>()),
        "bootstrap_nonfinite_meaning":"draw has no eventual IC advantage in this additive linear model, including an initial setup advantage that later disappears; retain its mass in unconditional extended interval",
        "scope":"MODEL ONLY: no measured execution above32 targets. Charges rho seed-specific G jump preparation per target while retaining curve/normal-basis preparation once. Uses producer preparation/query clocks with different online validation boundaries and excludes process/output residuals, independent audit and original retained-base construction. A nonpositive observed online slope does not prove universal absence of a crossover. A finite extrapolated threshold does not establish an observed batch win."})
}

fn summarize_aa(rows: &[&Value]) -> (Value, bool) {
    let mut slots = BTreeMap::<(u64, String), &Value>::new();
    let mut issues = Vec::<String>::new();
    for row in rows {
        let (Some(pair), Some(leg @ ("A" | "B"))) = (row["pair_id"].as_u64(), row["leg"].as_str())
        else {
            issues.push("AA row missing pair_id/leg".into());
            continue;
        };
        if !(1..=AA_PAIRS as u64).contains(&pair) || row["arm"] != "rho" || !status_ok(row, true) {
            issues.push(format!("AA{pair}/{leg}:invalid row"));
        }
        if slots.insert((pair, leg.into()), row).is_some() {
            issues.push(format!("AA{pair}/{leg}:duplicate"));
        }
    }
    if rows.len() != AA_PAIRS * 2 {
        issues.push("require10AA processes".into());
    }
    let mut ratios = Vec::new();
    for pair in 1..=AA_PAIRS as u64 {
        match (
            slots.get(&(pair, "A".into())),
            slots.get(&(pair, "B".into())),
        ) {
            (Some(a), Some(b)) => {
                if a["seed"].as_u64().is_none()
                    || a["target_id"].as_str().is_none()
                    || a["seed"] != b["seed"]
                    || a["target_id"] != b["target_id"]
                {
                    issues.push(format!("AA{pair}:missing or mismatched seed/target"));
                }
                if let (Some(a), Some(b)) = (
                    positive(&a["elapsed_process_wall_s"]),
                    positive(&b["elapsed_process_wall_s"]),
                ) {
                    ratios.push(a / b);
                } else {
                    issues.push(format!("AA{pair}:invalid time"));
                }
            }
            _ => issues.push(format!("AA{pair}:missing leg")),
        }
    }
    let valid = issues.is_empty();
    (
        json!({"valid":valid,"issues":issues,"A_over_B_wall_ratios":ratios,"wall_ratio_description":descriptive(&ratios),"scope":"Five same-target same-seed rho/rho process pairs; local timing diagnostic, not a confidence interval or noise ceiling. No correction is subtracted from IC/rho measurements."}),
        valid,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn panel() -> Value {
        let ids = (1..=TARGETS)
            .map(|i| format!("target_{i:02}"))
            .collect::<Vec<_>>();
        let mut process = Vec::new();
        let mut target = Vec::new();
        for cohort in 1..=COHORTS {
            for arm in ["ic", "rho"] {
                process.push(json!({"kind":"cohort","arm":arm,"cohort_id":cohort,"target_count":TARGETS,"status":"VERIFIED","all_verified":true,"elapsed_process_wall_s":if arm=="ic" {100.+cohort as f64} else {32.+cohort as f64},"producer_preparation_ms":if arm=="ic" {100000.} else {200.},"producer_per_target_preparation_ms":if arm=="ic" {0.}else{1.},"producer_shared_preparation_ms":if arm=="ic" {100000.} else {168.},"cpu_total_s":30.,"cpu_user_s":29.,"cpu_system_s":1.,"peak_rss_bytes":1000}));
                for id in &ids {
                    target.push(json!({"target_id":id,"cohort_id":cohort,"arm":arm,"status":"VERIFIED","verified":true,"online_wall_s":if arm=="ic" {0.2} else {1.}}));
                }
            }
        }
        for (i, id) in ids.iter().take(COLD_TARGETS).enumerate() {
            for arm in ["ic", "rho"] {
                process.push(json!({"kind":"cold","target_id":id,"target_count":1,"arm":arm,"status":"VERIFIED","all_verified":true,"elapsed_process_wall_s":if arm=="ic" {100.+i as f64} else {2.+i as f64},"producer_preparation_ms":if arm=="ic" {100000.} else {200.}}));
            }
        }
        for pair in 1..=AA_PAIRS {
            for leg in ["A", "B"] {
                process.push(json!({"kind":"aa","arm":"rho","pair_id":pair,"leg":leg,"target_id":ids[0],"seed":pair as u64,"status":"VERIFIED","all_verified":true,"target_count":1,"elapsed_process_wall_s":2.}));
            }
        }
        let hash = "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef";
        json!({"schema":"n83_matched_statistics_input_v1","metadata":{"target_ids":ids,"cold_target_ids":ids.iter().take(COLD_TARGETS).collect::<Vec<_>>(),"independently_verified":true,"protocol_sha256":hash,"campaign_manifest_sha256":hash,"independent_verification_receipt_sha256":hash},"process_rows":process,"target_rows":target,"independent_verification":{"status":"PASS","checked_target_rows":218,"checked_aa_semantic_pairs":5,"aa_semantic_gate":true,"invalid_rows":0,"campaign_manifest_sha256":hash,"protocol_sha256":hash}})
    }

    #[test]
    fn correct_whole_process_and_nested_estimands() {
        let r = summarize(&panel());
        assert_eq!(r["all_required_rows_valid"], true);
        assert!(
            (r["primary_prepared_batch_process_wall"]["ic_over_rho"]
                .as_f64()
                .unwrap()
                - 3.)
                .abs()
                < 1e-12
        );
        assert!(
            (r["primary_actual_cold_singleton_process_wall"]["ic_over_rho"]
                .as_f64()
                .unwrap()
                - 103.5 / 5.5)
                .abs()
                < 1e-12
        );
        assert!(
            (r["secondary_prepared_online_target_means"]["ic_over_rho"]
                .as_f64()
                .unwrap()
                - 0.2)
                .abs()
                < 1e-12
        );
        assert_eq!(r["requirements"]["expected_query_rows"], 192);
        assert_eq!(r["requirements"]["expected_process_rows"], 32);
        assert_eq!(
            r["modeled_break_even"]["first_n_with_model_ic_at_most_rho"],
            125.
        );
        assert_eq!(
            r["modeled_break_even"]["point_outside_measured_batch"],
            true
        );
    }

    #[test]
    fn conventional_even_median_and_tail_rule() {
        assert_eq!(quantile(&[1., 2., 3., 4.], 0.5), Some(2.5));
        assert_eq!(quantile(&[1., 2., 3., 4.], 0.9), Some(3.7));
        assert_eq!(quantile(&[], 0.5), None);
    }

    #[test]
    fn missing_seed_draw_cannot_be_imputed_or_completed_subset() {
        let mut p = panel();
        p["target_rows"].as_array_mut().unwrap().remove(0);
        let r = summarize(&p);
        assert_eq!(r["all_required_rows_valid"], false);
        assert!(r["primary_prepared_batch_process_wall"]["ic_over_rho"].is_null());
        assert!(r["secondary_prepared_online_target_means"]["ic_over_rho"].is_null());
        assert_eq!(
            r["primary_actual_cold_singleton_process_wall"]["valid"],
            true
        );
        assert_eq!(r["overall_primary_winner_eligible"], false);
    }

    #[test]
    fn duplicate_target_and_process_reject_matrices() {
        let mut p = panel();
        let duplicate = p["target_rows"][0].clone();
        p["target_rows"].as_array_mut().unwrap().push(duplicate);
        let r = summarize(&p);
        assert_eq!(r["all_required_rows_valid"], false);
        assert_eq!(r["primary_prepared_batch_process_wall"]["valid"], false);
        let mut p = panel();
        let duplicate = p["process_rows"][0].clone();
        p["process_rows"].as_array_mut().unwrap().push(duplicate);
        assert_eq!(
            summarize(&p)["primary_prepared_batch_process_wall"]["valid"],
            false
        );
    }

    #[test]
    fn failure_cost_retained_without_primary_success_claim() {
        let mut p = panel();
        p["process_rows"][0]["status"] = json!("TIMEOUT");
        p["process_rows"][0]["all_verified"] = json!(false);
        p["process_rows"][0]["elapsed_process_wall_s"] = json!(2100.);
        let r = summarize(&p);
        assert_eq!(r["primary_prepared_batch_process_wall"]["valid"], false);
        assert_eq!(
            r["costs"]["ic"]["prepared_cohort_process"]["wall_s"]["max"],
            2100.
        );
        assert_eq!(
            r["costs"]["ic"]["prepared_cohort_process"]["verified_rows"],
            2
        );
    }

    #[test]
    fn false_verified_flag_or_wrong_answer_rejects() {
        let mut p = panel();
        p["target_rows"][0]["verified"] = json!(false);
        assert_eq!(
            summarize(&p)["primary_prepared_batch_process_wall"]["valid"],
            false
        );
        let mut p = panel();
        p["process_rows"][7]["wrong_answer"] = json!(true);
        assert_eq!(
            summarize(&p)["primary_actual_cold_singleton_process_wall"]["valid"],
            false
        );
    }

    #[test]
    fn nonfinite_threshold_mass_is_not_discarded() {
        let mut p = panel();
        for row in p["target_rows"].as_array_mut().unwrap() {
            if row["arm"] == "ic" {
                row["online_wall_s"] = json!(2.);
            }
        }
        let r = summarize(&p);
        assert!(r["modeled_break_even"]["first_n_with_model_ic_at_most_rho"].is_null());
        assert_eq!(r["modeled_break_even"]["bootstrap_nonfinite_draws"], DRAWS);
        assert_eq!(r["modeled_break_even"]["bootstrap_nonfinite_fraction"], 1.);
        assert_eq!(
            r["modeled_break_even"]["extended_interval_95"][1]["kind"],
            "no_eventual_ic_advantage_in_linear_model"
        );
        assert!(r["modeled_break_even"]["finite_only_descriptive_quantiles"]["median"].is_null());
    }

    #[test]
    fn bootstrap_is_reproducible_and_keeps_pairs() {
        let pairs = [(10., 1.), (20., 2.), (30., 3.)];
        let a = comparison(&pairs, true, SEED, "paired", "test");
        let b = comparison(&pairs, true, SEED, "paired", "test");
        assert_eq!(a, b);
        assert_eq!(a["interval_95"], json!([10., 10.]));
    }

    #[test]
    fn aa_missing_seed_and_wrong_order_reject() {
        let mut p = panel();
        let rows = p["process_rows"].as_array_mut().unwrap();
        for row in rows {
            if row["kind"] == "aa" {
                row.as_object_mut().unwrap().remove("seed");
            }
        }
        assert_eq!(summarize(&p)["aa"]["valid"], false);
        let mut p = panel();
        p["metadata"]["cold_target_ids"]
            .as_array_mut()
            .unwrap()
            .swap(0, 1);
        assert_eq!(
            summarize(&p)["primary_actual_cold_singleton_process_wall"]["valid"],
            false
        );
    }

    #[test]
    fn missing_independent_gate_withholds_all_claims() {
        let mut p = panel();
        p.as_object_mut()
            .unwrap()
            .remove("independent_verification");
        let r = summarize(&p);
        assert_eq!(r["all_required_rows_valid"], false);
        assert_eq!(r["primary_prepared_batch_process_wall"]["valid"], false);
        assert_eq!(
            r["primary_actual_cold_singleton_process_wall"]["valid"],
            false
        );
        assert_eq!(r["modeled_break_even"]["valid"], false);
    }

    #[test]
    fn false_gate_and_mismatched_protocol_are_rejected() {
        let mut p = panel();
        p["metadata"]["independently_verified"] = json!(false);
        assert_eq!(summarize(&p)["all_required_rows_valid"], false);
        let mut p = panel();
        p["independent_verification"]["protocol_sha256"] = json!("a".repeat(64));
        assert_eq!(summarize(&p)["all_required_rows_valid"], false);
    }

    #[test]
    fn aa_semantic_receipt_is_mandatory_for_primary_gates() {
        let mut p = panel();
        p["independent_verification"]
            .as_object_mut()
            .unwrap()
            .remove("checked_aa_semantic_pairs");
        let r = summarize(&p);
        assert_eq!(r["all_required_rows_valid"], false);
        assert_eq!(r["primary_prepared_batch_process_wall"]["valid"], false);
        assert_eq!(r["aa"]["valid"], false);
        assert_eq!(
            r["primary_actual_cold_singleton_process_wall"]["valid"],
            false
        );
    }

    #[test]
    fn slower_query_with_lower_setup_is_initial_advantage_not_eventual_crossover() {
        assert_eq!(
            threshold(0., 10., 2., 1.),
            ModelOutcome {
                first_n: Some(1.),
                eventual_n: None,
                last_n: Some(10.),
                kind: "ic_initial_setup_advantage_lost_as_batch_grows"
            }
        );
        assert_eq!(threshold(0., 1., 2., 2.).eventual_n, Some(1.));
        assert_eq!(threshold(2., 1., 2., 2.).first_n, None);
        assert_eq!(threshold(2., 1., 2., 1.).first_n, None);
    }

    #[test]
    fn rho_jump_work_scales_per_target_not_fixed_batch_intercept() {
        let mut p = panel();
        for row in p["process_rows"].as_array_mut().unwrap() {
            if row["kind"] == "cohort" && row["arm"] == "rho" {
                row["producer_preparation_ms"] = json!(32100.);
                row["producer_shared_preparation_ms"] = json!(100.);
                row["producer_per_target_preparation_ms"] = json!(1000.);
            }
        }
        let r = summarize(&p);
        assert!(
            (r["modeled_break_even"]["preparation_rho_shared_s"]
                .as_f64()
                .unwrap()
                - 0.1)
                .abs()
                < 1e-12
        );
        assert_eq!(
            r["modeled_break_even"]["rho_seed_specific_jump_per_target_s"],
            1.
        );
        assert_eq!(
            r["modeled_break_even"]["rho_per_target_online_plus_jump_s"],
            2.
        );
        assert_eq!(
            r["modeled_break_even"]["first_n_with_model_ic_at_most_rho"],
            56.
        );
    }
}

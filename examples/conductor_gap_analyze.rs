//! Independent-target analysis for the frozen large-prime conductor-gap follow-up.
//! Reads sealed ecbench records; never launches a solver or prices wall time.
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::path::Path;

use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::ecbench::runner::read_records;
use crypto_lib::cryptanalysis::ecbench::stats::{mean, quantile, Resampler};
use serde_json::{json, Value};

const SAMPLES: usize = 20_000;
const SEED: u64 = 20261012;

fn cross_curve(a: &[f64], b: &[f64]) -> Value {
    let ratio = mean(b).unwrap() / mean(a).unwrap();
    let mut rng = Resampler::new(SEED);
    let mut boot = Vec::with_capacity(SAMPLES);
    for _ in 0..SAMPLES {
        let sa: f64 = (0..a.len()).map(|_| a[rng.below(a.len())]).sum();
        let sb: f64 = (0..b.len()).map(|_| b[rng.below(b.len())]).sum();
        boot.push((sb / b.len() as f64) / (sa / a.len() as f64));
    }
    let mut pool = a.iter().chain(b).copied().collect::<Vec<_>>();
    let observed = ratio.ln().abs();
    let mut extreme = 0;
    for _ in 0..SAMPLES {
        for i in (1..pool.len()).rev() {
            let j = rng.below(i + 1);
            pool.swap(i, j);
        }
        let ra = mean(&pool[..a.len()]).unwrap();
        let rb = mean(&pool[a.len()..]).unwrap();
        if (rb / ra).ln().abs() + 1e-12 >= observed {
            extreme += 1;
        }
    }
    json!({
        "ratio_floor_over_crater": ratio,
        "ci95": [quantile(&boot, 0.025), quantile(&boot, 0.975)],
        "ci99": [quantile(&boot, 0.005), quantile(&boot, 0.995)],
        "permutation_p_two_sided": (extreme + 1) as f64 / (SAMPLES + 1) as f64,
        "resamples": SAMPLES, "seed": SEED,
        "method": "independent workload bootstrap on each fixed curve; label permutation",
        "cross_curve_targets_transported": false,
    })
}

fn arg(args: &[String], name: &str) -> String {
    let i = args
        .iter()
        .position(|s| s == name)
        .expect("required argument");
    args.get(i + 1).expect("argument value").clone()
}

fn analyze(dir: &Path, ell: u64) -> Result<Value, String> {
    let spec: Value =
        serde_json::from_slice(&std::fs::read(dir.join("spec.json")).map_err(|e| e.to_string())?)
            .map_err(|e| e.to_string())?;
    let curves = spec["workloads"]["curves"].as_array().ok_or("curves")?;
    if curves.len() != 2 {
        return Err("the protocol requires exactly two fixed curves".into());
    }
    let records = read_records(dir)?;
    let measured = records.iter().filter(|r| !r.warmup).collect::<Vec<_>>();
    let failed = measured.iter().filter(|r| r.outcome.status != "verified").map(|r| {
        json!({"run_id":r.run_id,"arm":r.arm,"status":r.outcome.status,"error":r.outcome.error})
    }).collect::<Vec<_>>();
    let mut arms = BTreeMap::<String, Value>::new();
    let mut aa_pairs = 0;
    let mut aa_mismatch = Vec::new();
    let mut table = Vec::new();
    for arm in ["rho-neg", "rho-neg-aa", "ic-mitm"] {
        let mut role_rows = Vec::new();
        let mut samples = Vec::<Vec<f64>>::new();
        for (ci, role) in ["crater", "floor"].iter().enumerate() {
            let selected = measured
                .iter()
                .filter(|r| {
                    r.arm == arm
                        && serde_json::to_value(&r.workload.curve_spec).unwrap() == curves[ci]
                })
                .copied()
                .collect::<Vec<_>>();
            let mut workloads = BTreeMap::<String, Vec<f64>>::new();
            let mut statuses = BTreeMap::<String, usize>::new();
            let mut phases = BTreeMap::<String, f64>::new();
            let mut unpriced = BTreeSet::new();
            let mut bounded = false;
            let mut ids = Vec::new();
            let mut levels = BTreeMap::<String, usize>::new();
            for r in &selected {
                *statuses.entry(r.outcome.status.clone()).or_default() += 1;
                *levels
                    .entry(r.isolation.level.name().to_string())
                    .or_default() += 1;
                ids.push(r.run_id.clone());
                bounded |= r.cost.lower_bound;
                unpriced.extend(r.cost.unpriced.clone());
                if r.outcome.status == "verified" {
                    let s = r.cost.s.ok_or("verified run has no S")?;
                    if !s.is_finite() || s <= 0.0 {
                        return Err("invalid S".into());
                    }
                    workloads
                        .entry(r.workload.workload_id.clone())
                        .or_default()
                        .push(s);
                    for p in &r.phases {
                        *phases.entry(p.name.clone()).or_default() += p.gae;
                    }
                }
                if arm == "rho-neg" {
                    let controls = selected_control(&measured, r);
                    aa_pairs += 1;
                    if controls.len() != 1
                        || controls[0].cost != r.cost
                        || controls[0].counters != r.counters
                        || controls[0].outcome.recovered != r.outcome.recovered
                    {
                        aa_mismatch.push(r.run_id.clone());
                    }
                }
            }
            let valid = selected.iter().all(|r| r.outcome.status == "verified")
                && selected.len() == 8
                && workloads.len() == 8;
            let values = workloads
                .values()
                .map(|v| mean(v).unwrap())
                .collect::<Vec<_>>();
            samples.push(values.clone());
            let mean_s = valid.then(|| mean(&values)).flatten();
            let n = selected.len() as f64;
            let phase_means = phases
                .into_iter()
                .map(|(k, v)| (k, v / n.max(1.0)))
                .collect::<BTreeMap<_, _>>();
            let slug = selected.first().map(|r| r.workload.curve.slug.clone());
            let row = json!({
                "role": role, "slug": slug, "measured":selected.len(),
                "workloads":workloads.len(),"statuses":statuses,"valid":valid,
                "mean_s":mean_s,"min_s":values.iter().copied().reduce(f64::min),
                "max_s":values.iter().copied().reduce(f64::max),
                "lower_bound":bounded,"unpriced":unpriced,"isolation_levels":levels,
                "phase_mean_gae":phase_means,"run_ids":ids,
                "factor_base":selected.first().and_then(|r|r.factor_base.clone()),
                "floor_s":selected.first().map(|r|r.boundaries.floor_s),
            });
            table.push(json!({"arm":arm,"row":row.clone()}));
            role_rows.push(row);
        }
        let admissible = role_rows.iter().all(|r| r["valid"] == true);
        arms.insert(arm.to_string(),json!({"curves":role_rows,
            "cross_curve":if admissible {cross_curve(&samples[0],&samples[1])} else {json!({"status":"inadmissible: measured failures or missing workloads"})}
        }));
    }
    let slugs = table
        .iter()
        .filter_map(|r| r["row"]["slug"].as_str())
        .collect::<BTreeSet<_>>();
    if slugs.len() != 2 {
        return Err("curve roles did not resolve uniquely".into());
    }
    let mut by_curve = BTreeMap::new();
    for role in ["crater", "floor"] {
        let rho = arms["rho-neg"]["curves"]
            .as_array()
            .unwrap()
            .iter()
            .find(|r| r["role"] == role)
            .unwrap();
        let ic = arms["ic-mitm"]["curves"]
            .as_array()
            .unwrap()
            .iter()
            .find(|r| r["role"] == role)
            .unwrap();
        by_curve.insert(role, json!({"ic_over_rho_count_diagnostic":ic["mean_s"].as_f64().zip(rho["mean_s"].as_f64()).map(|(i,r)|i/r),"ratio_of_lower_bounds":rho["lower_bound"]==true || ic["lower_bound"]==true}));
    }
    let source = std::fs::read(dir.join("records.jsonl")).map_err(|e| e.to_string())?;
    Ok(
        json!({"ell":ell,"session":dir,"session_id":records.first().map(|r|&r.session_id),
        "records":records.len(),"measured":measured.len(),"failed":failed,
        "records_sha256":sha256_hex(&source),"arms":arms,"ic_vs_rho":by_curve,
        "aa":{"pairs":aa_pairs,"mismatches":aa_mismatch},
        "status":"observations on two fixed curves; no population hardness or speedup claim"}),
    )
}

fn selected_control<'a>(
    measured: &[&'a crypto_lib::cryptanalysis::ecbench::record::Record],
    r: &crypto_lib::cryptanalysis::ecbench::record::Record,
) -> Vec<&'a crypto_lib::cryptanalysis::ecbench::record::Record> {
    measured
        .iter()
        .filter(|x| {
            x.arm == "rho-neg-aa"
                && x.workload.workload_id == r.workload.workload_id
                && x.round == r.round
        })
        .copied()
        .collect()
}

fn main() -> Result<(), String> {
    let args = std::env::args().collect::<Vec<_>>();
    let root = arg(&args, "--sessions");
    let out = arg(&args, "--out");
    let mut sessions = Vec::new();
    for (name, ell) in [("ell8209-q29", 8209), ("ell131101-q37", 131101)] {
        sessions.push(analyze(&Path::new(&root).join(name), ell)?);
    }
    let mut markdown = String::from("# Frozen conductor-gap count diagnostic\n\nIndependent target samples on each fixed curve. Ratios of lower bounds are diagnostics; no wall-time or intrinsic-hardness claim.\n\n| ell | arm | floor/crater | 95% interval | 99% interval | permutation p |\n|---:|---|---:|---|---|---:|\n");
    for s in &sessions {
        for arm in ["rho-neg", "rho-neg-aa", "ic-mitm"] {
            let x = &s["arms"][arm]["cross_curve"];
            writeln!(
                &mut markdown,
                "| {} | {} | {} | {} | {} | {} |",
                s["ell"],
                arm,
                x["ratio_floor_over_crater"],
                x["ci95"],
                x["ci99"],
                x["permutation_p_two_sided"]
            )
            .unwrap();
        }
    }
    markdown.push_str("\n| ell | curve | arm | measured / verified | mean counted S | S range | lower bound |\n|---:|---|---|---|---:|---|---|\n");
    for s in &sessions {
        for arm in ["rho-neg", "rho-neg-aa", "ic-mitm"] {
            for r in s["arms"][arm]["curves"].as_array().unwrap() {
                writeln!(
                    &mut markdown,
                    "| {} | {} ({}) | {} | {} / {} | {} | {}–{} | {} |",
                    s["ell"],
                    r["slug"].as_str().unwrap(),
                    r["role"].as_str().unwrap(),
                    arm,
                    r["measured"],
                    r["statuses"]["verified"],
                    r["mean_s"],
                    r["min_s"],
                    r["max_s"],
                    r["lower_bound"]
                )
                .unwrap();
            }
        }
    }
    let doc = json!({"schema":"conductor-gap-analysis/v1","protocol":"research/ecbench_conductor_gap_20261012/PROTOCOL.md","resamples":SAMPLES,"seed":SEED,"sessions":sessions});
    std::fs::write(
        format!("{out}.json"),
        serde_json::to_vec_pretty(&doc).unwrap(),
    )
    .map_err(|e| e.to_string())?;
    std::fs::write(format!("{out}.md"), markdown).map_err(|e| e.to_string())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn constant_equal_samples_have_null_interval_and_p_one() {
        let x = cross_curve(&[2.0; 8], &[2.0; 8]);
        assert_eq!(x["ratio_floor_over_crater"], 1.0);
        assert_eq!(x["ci99"], json!([1.0, 1.0]));
        assert_eq!(x["permutation_p_two_sided"], 1.0);
    }
    #[test]
    fn separated_samples_preserve_direction_and_bootstrap_scale() {
        let x = cross_curve(&[2.0; 8], &[4.0; 8]);
        assert_eq!(x["ratio_floor_over_crater"], 2.0);
        assert_eq!(x["ci99"], json!([2.0, 2.0]));
        assert!(x["permutation_p_two_sided"].as_f64().unwrap() < 0.01);
    }
}

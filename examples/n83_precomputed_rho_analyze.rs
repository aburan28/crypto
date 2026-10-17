//! Verifier and analyzer for one job of the N83 precomputed-rho control.
//!
//! Usage: n83_precomputed_rho_analyze <job-dir> <ic-label> <analysis.json>
//!
//! `<job-dir>` holds `targets.jsonl`, `ic/records.jsonl`, `ic/stdout.log`
//! (the IC setup summary) and `ic/time.txt` (GNU time -v), `rho_plain/w<NN>.json` (one strong-rho receipt per target),
//! and `bl_a/records.jsonl`, `bl_b/records.jsonl` (precomputed-rho records).
//! Every recovered scalar of every arm is replayed here with the generic
//! library curve (`[d]G == Q`), independently of the producers, and the arms
//! must agree target by target. Online-interval statistics: mean and median
//! with 95% percentile bootstrap intervals over targets (10,000 resamples,
//! fixed seed), and paired IC/arm ratios resampled by target. Operation
//! counts are priced at 0, 0.1 and 1 group-addition equivalents per
//! unpriced native operation (AGENTS.md 5a).

use std::fs;

use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use rand::{rngs::StdRng, Rng, SeedableRng};
use serde_json::{json, Value};

const RESAMPLES: usize = 10_000;
const WEIGHTS: [f64; 3] = [0.0, 0.1, 1.0];

fn read_jsonl(path: &str) -> Result<Vec<Value>, String> {
    let text = fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| serde_json::from_str(l).map_err(|e| format!("{path}: {e}")))
        .collect()
}

fn as_u128(v: &Value) -> Option<u128> {
    match v {
        Value::String(s) => s.parse().ok(),
        Value::Number(n) => n.as_u64().map(u128::from).or_else(|| n.to_string().parse().ok()),
        _ => None,
    }
}

fn parse_targets(path: &str) -> Result<Vec<[u128; 2]>, String> {
    let text = fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    let mut out = Vec::new();
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let inner = line.trim().trim_start_matches('[').trim_end_matches(']');
        let parts: Vec<u128> = inner
            .split(',')
            .map(|p| p.trim().trim_matches('"').parse::<u128>())
            .collect::<Result<_, _>>()
            .map_err(|_| format!("{path}: bad target line"))?;
        if parts.len() != 2 {
            return Err(format!("{path}: bad target line"));
        }
        out.push([parts[0], parts[1]]);
    }
    Ok(out)
}

/// GNU `time -v` maximum resident set size, in bytes.
fn time_max_rss(path: &str) -> Option<u64> {
    let text = fs::read_to_string(path).ok()?;
    text.lines()
        .find_map(|l| l.trim().strip_prefix("Maximum resident set size (kbytes):"))
        .and_then(|v| v.trim().parse::<u64>().ok())
        .map(|k| k * 1024)
}

fn time_field(path: &str, label: &str) -> Option<String> {
    let text = fs::read_to_string(path).ok()?;
    text.lines()
        .find_map(|l| l.trim().strip_prefix(label).map(|v| v.trim().to_string()))
}

struct Curve {
    curve: KoblitzCurve,
}

impl Curve {
    fn verify(&self, d: u64, q: [u128; 2]) -> bool {
        let point = BinaryPoint::Affine {
            x: F2mElement::from_biguint(&BigUint::from(q[0]), 83),
            y: F2mElement::from_biguint(&BigUint::from(q[1]), 83),
        };
        self.curve.mul(self.curve.generator(), &BigUint::from(d)) == point
    }
}

#[derive(Clone)]
struct Arm {
    name: String,
    online_ms: Vec<f64>,
    scalars: Vec<Option<u64>>,
    statuses: Vec<String>,
    /// (group additions, unpriced native operations) per target, if emitted.
    counts: Vec<Option<(f64, f64)>>,
    count_note: String,
}

fn mean(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / v.len() as f64
}

fn median(v: &[f64]) -> f64 {
    let mut s = v.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = s.len();
    if n % 2 == 1 {
        s[n / 2]
    } else {
        0.5 * (s[n / 2 - 1] + s[n / 2])
    }
}

fn percentile_ci(mut v: Vec<f64>) -> [f64; 2] {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let at = |p: f64| v[((p * (v.len() - 1) as f64).round()) as usize];
    [at(0.025), at(0.975)]
}

/// Bootstrap over target indices; `stat` maps a resampled index set to a value.
fn bootstrap(n: usize, seed: u64, stat: impl Fn(&[usize]) -> f64) -> [f64; 2] {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut values = Vec::with_capacity(RESAMPLES);
    let mut idx = vec![0usize; n];
    for _ in 0..RESAMPLES {
        for slot in idx.iter_mut() {
            *slot = rng.gen_range(0..n);
        }
        values.push(stat(&idx));
    }
    percentile_ci(values)
}

fn pick(v: &[f64], idx: &[usize]) -> Vec<f64> {
    idx.iter().map(|&i| v[i]).collect()
}

fn arm_stats(arm: &Arm, seed: u64) -> Value {
    let v = &arm.online_ms;
    json!({
        "arm": arm.name,
        "n": v.len(),
        "online_ms": v,
        "mean_ms": mean(v),
        "mean_ms_ci95": bootstrap(v.len(), seed, |i| mean(&pick(v, i))),
        "median_ms": median(v),
        "median_ms_ci95": bootstrap(v.len(), seed ^ 1, |i| median(&pick(v, i))),
        "min_ms": v.iter().cloned().fold(f64::INFINITY, f64::min),
        "max_ms": v.iter().cloned().fold(0.0, f64::max),
        "statuses": arm.statuses,
        "count_note": arm.count_note,
    })
}

/// Paired comparison: `num` over `den` (for example IC over precomputed rho).
fn paired(num: &Arm, den: &Arm, seed: u64) -> Value {
    let (a, b) = (&num.online_ms, &den.online_ms);
    let ratios: Vec<f64> = a.iter().zip(b).map(|(x, y)| x / y).collect();
    let n = ratios.len();
    let band: Vec<Value> = WEIGHTS
        .iter()
        .map(|&w| {
            let cost = |arm: &Arm| -> Option<Vec<f64>> {
                arm.counts
                    .iter()
                    .map(|c| c.map(|(g, u)| g + w * u))
                    .collect()
            };
            match (cost(num), cost(den)) {
                (Some(x), Some(y)) => json!({
                    "weight_gae_per_native_op": w,
                    "mean_ratio": mean(&x) / mean(&y),
                    "median_ratio": median(&x) / median(&y),
                    "numerator_mean_gae": mean(&x),
                    "denominator_mean_gae": mean(&y),
                }),
                _ => json!({"weight_gae_per_native_op": w, "status": "COUNT_NOT_EMITTED"}),
            }
        })
        .collect();
    json!({
        "numerator": num.name,
        "denominator": den.name,
        "per_target_ratio": ratios,
        "median_of_paired_ratios": median(&ratios),
        "median_of_paired_ratios_ci95": bootstrap(n, seed, |i| median(&pick(&ratios, i))),
        "ratio_of_means": mean(a) / mean(b),
        "ratio_of_means_ci95": bootstrap(n, seed ^ 2, |i| mean(&pick(a, i)) / mean(&pick(b, i))),
        "ratio_of_medians": median(a) / median(b),
        "ratio_of_medians_ci95": bootstrap(n, seed ^ 3, |i| median(&pick(a, i)) / median(&pick(b, i))),
        "targets_where_numerator_slower": ratios.iter().filter(|&&r| r > 1.0).count(),
        "operation_count_band": band,
    })
}

fn load_ic(dir: &str, targets: &[[u128; 2]]) -> Result<(Arm, Value), String> {
    let records = read_jsonl(&format!("{dir}/records.jsonl"))?;
    let rows: Vec<&Value> = records
        .iter()
        .filter(|r| r["kind"] == "compact_orbit_dlp_target")
        .collect();
    // Both frozen producers print their setup summary on stdout.
    let stdout = read_jsonl(&format!("{dir}/stdout.log")).unwrap_or_default();
    let summary = records
        .iter()
        .chain(stdout.iter())
        .find(|r| r["kind"] == "compact_orbit_dlp_summary")
        .ok_or("IC summary missing")?;
    if rows.len() != targets.len() {
        return Err(format!("IC rows {} != targets {}", rows.len(), targets.len()));
    }
    let mut arm = Arm {
        name: "ic".into(),
        online_ms: vec![],
        scalars: vec![],
        statuses: vec![],
        counts: vec![],
        count_note: "group additions = 4 per relation check + 1 recovery check (lower bound); unpriced = probes + root_calls where emitted".into(),
    };
    for (row, t) in rows.iter().zip(targets) {
        let q = &row["published_q"];
        if as_u128(&q[0]) != Some(t[0]) || as_u128(&q[1]) != Some(t[1]) {
            return Err("IC target order differs from targets.jsonl".into());
        }
        let ms = row["online_ms"]
            .as_f64()
            .or_else(|| row["target_ms"].as_f64())
            .ok_or("IC online time missing")?;
        arm.online_ms.push(ms);
        arm.scalars.push(row["recovered_scalar"].as_u64());
        arm.statuses.push(if row["group_verified"] == true { "PRODUCER_VERIFIED".into() } else { "PRODUCER_UNVERIFIED".into() });
        let probes = row["probes"].as_f64();
        let roots = row["root_calls"].as_f64().unwrap_or(0.0);
        let checks = row["relation_checks"].as_f64().unwrap_or(1.0);
        arm.counts.push(probes.map(|p| (4.0 * checks + 1.0, p + roots)));
    }
    let timing = &summary["timing_ms"];
    let process = timing["process_total"].as_f64();
    let targets_total = timing["targets_total"].as_f64();
    let info = json!({
        "candidate_columns": summary["orbit_columns"],
        "factor_base_points": summary["factor_base_points"],
        "regular_states": summary["regular_states"],
        "root_table_entries": summary["root_table_entries"],
        "root_state_allocated_bytes": summary["root_state_allocated_bytes"],
        "root_table_allocated_bytes": summary["root_table_allocated_bytes"],
        "summary_peak_rss_bytes": summary["peak_rss_bytes"],
        "time_v_max_rss_bytes": time_max_rss(&format!("{dir}/time.txt")),
        "time_v_user_seconds": time_field(&format!("{dir}/time.txt"), "User time (seconds):"),
        "time_v_wall": time_field(&format!("{dir}/time.txt"), "Elapsed (wall clock) time (h:mm:ss or m:ss):"),
        "timing_ms": timing,
        "target_independent_ms": process.zip(targets_total).map(|(p, t)| p - t),
        "rank": summary["rank"],
        "targets_solved": summary["targets_solved"],
        "threads": summary["threads"],
        "target_threads": summary["target_threads"],
    });
    Ok((arm, info))
}

fn load_rho_plain(dir: &str, targets: &[[u128; 2]]) -> Result<Arm, String> {
    let mut arm = Arm {
        name: "rho_plain_fold".into(),
        online_ms: vec![],
        scalars: vec![],
        statuses: vec![],
        counts: vec![],
        count_note: "group additions = charges.group_additions; unpriced = canonicalizations + partition_hashes + table_queries + table_inserts".into(),
    };
    for (i, t) in targets.iter().enumerate() {
        let path = format!("{dir}/w{i:02}.json");
        let text = fs::read_to_string(&path).map_err(|e| format!("{path}: {e}"))?;
        let r: Value = serde_json::from_str(text.trim()).map_err(|e| format!("{path}: {e}"))?;
        let q = &r["public_target_q"];
        if as_u128(&q[0]) != Some(t[0]) || as_u128(&q[1]) != Some(t[1]) {
            return Err(format!("{path}: target differs"));
        }
        arm.online_ms.push(r["online_ms"].as_f64().ok_or("rho online_ms missing")?);
        arm.scalars.push(r["recovered_scalar"].as_str().and_then(|s| s.parse().ok()));
        arm.statuses.push(r["status"].as_str().unwrap_or("MISSING").into());
        let c = &r["outcome"]["charges_setup_and_solve"];
        let g = c["group_additions"].as_f64();
        let u = ["canonicalizations", "partition_hashes", "table_queries", "table_inserts"]
            .iter()
            .map(|k| c[*k].as_f64())
            .sum::<Option<f64>>();
        arm.counts.push(g.zip(u));
    }
    Ok(arm)
}

fn load_bl(dir: &str, name: &str, targets: &[[u128; 2]]) -> Result<(Arm, Value), String> {
    let records = read_jsonl(&format!("{dir}/records.jsonl"))?;
    let setup = records
        .iter()
        .find(|r| r["record"] == "setup")
        .ok_or(format!("{name}: setup record missing"))?
        .clone();
    let end = records.iter().find(|r| r["record"] == "end").cloned();
    let rows: Vec<&Value> = records.iter().filter(|r| r["record"] == "target").collect();
    if rows.len() != targets.len() {
        return Err(format!("{name}: rows {} != targets {}", rows.len(), targets.len()));
    }
    let mut arm = Arm {
        name: name.into(),
        online_ms: vec![],
        scalars: vec![],
        statuses: vec![],
        counts: vec![],
        count_note: "group additions = steps + escapes (doublings) + walks (start additions); unpriced = canonicalizations (steps + escapes + walks) + partition hashes (steps) + table lookups (distinguished points); scalar multiplications for walk setup and the recovery check are not counted".into(),
    };
    for (row, t) in rows.iter().zip(targets) {
        let q = &row["public_target_q"];
        if as_u128(&q[0]) != Some(t[0]) || as_u128(&q[1]) != Some(t[1]) {
            return Err(format!("{name}: target order differs"));
        }
        arm.online_ms.push(row["online_ms"].as_f64().ok_or("online_ms missing")?);
        arm.scalars.push(row["recovered_scalar"].as_str().and_then(|s| s.parse().ok()));
        arm.statuses.push(row["status"].as_str().unwrap_or("MISSING").into());
        let w = &row["online_walks"];
        let (steps, walks, escapes, dps) = (
            w["steps"].as_f64().unwrap_or(0.0),
            w["walks"].as_f64().unwrap_or(0.0),
            w["escapes"].as_f64().unwrap_or(0.0),
            w["distinguished"].as_f64().unwrap_or(0.0),
        );
        arm.counts.push(Some((steps + escapes + walks, (steps + escapes + walks) + steps + dps)));
    }
    let info = json!({
        "setup_ms": setup["setup_ns"].as_f64().map(|n| n / 1e6),
        "setup_cpu_seconds": setup["setup_cpu_seconds"],
        "setup_phase_ns": setup["setup_phase_ns"],
        "table_walks": setup["table_walks"],
        "table_walk_rate_steps_per_second": setup["table_walk_rate_steps_per_second"],
        "table_distinct_entries": setup["table_distinct_entries"],
        "table_raw_entries": setup["table_raw_entries"],
        "table_bytes": setup["table_bytes"],
        "table_stopped_by": setup["table_stopped_by"],
        "predicted_online_steps": setup["predicted_online_steps"],
        "config": setup["config"],
        "executable_sha256": setup["executable_sha256"],
        "worker_source_sha256": setup["worker_source_sha256"],
        "peak_rss_bytes": end.as_ref().map(|e| e["peak_rss_bytes"].clone()),
        "process_cpu_seconds": end.as_ref().map(|e| e["process_cpu_seconds"].clone()),
        "online_steps": rows.iter().map(|r| r["online_walks"]["steps"].clone()).collect::<Vec<_>>(),
        "recovery_source": rows.iter().map(|r| r["recovery_source"].clone()).collect::<Vec<_>>(),
    });
    Ok((arm, info))
}

fn run(job: &str, ic_label: &str, out: &str) -> Result<bool, String> {
    let curve = KoblitzCurve::new(1, 83).ok_or("curve")?;
    if curve.subgroup_order.to_u64() != Some(8_569_786_107_849_059) {
        return Err("unexpected curve".into());
    }
    let curve = Curve { curve };
    let targets = parse_targets(&format!("{job}/targets.jsonl"))?;
    let (ic, ic_info) = load_ic(&format!("{job}/ic"), &targets)?;
    let rho = load_rho_plain(&format!("{job}/rho_plain"), &targets)?;
    let (bl_a, a_info) = load_bl(&format!("{job}/bl_a"), "bl_a_time_matched", &targets)?;
    let (bl_b, b_info) = load_bl(&format!("{job}/bl_b"), "bl_b_memory_matched", &targets)?;
    // Amendment 1 (run 2 onward): time-matched to the measured IC setup wall.
    let a2 = if fs::metadata(format!("{job}/bl_a2/records.jsonl")).is_ok() {
        Some(load_bl(&format!("{job}/bl_a2"), "bl_a2_wall_equal", &targets)?)
    } else {
        None
    };
    let mut arms = vec![&ic, &rho, &bl_a, &bl_b];
    if let Some((a2_arm, _)) = &a2 {
        arms.push(a2_arm);
    }
    // Optional sidecar: the scalars the hosted single-target run recovered
    // for these public points (read only after every arm has finished).
    let reference: Option<Vec<u64>> = fs::read_to_string(format!("{job}/reference_scalars.txt"))
        .ok()
        .map(|t| t.lines().filter(|l| !l.trim().is_empty()).map(|l| l.trim().parse().unwrap()).collect());
    // Independent replay and cross-arm agreement.
    let mut verification = Vec::new();
    let mut all_ok = true;
    for (i, t) in targets.iter().enumerate() {
        let mut row = json!({"target_index": i});
        let mut agreed: Option<u64> = None;
        let mut agree = true;
        for arm in &arms {
            let d = arm.scalars[i];
            let ok = d.is_some_and(|d| curve.verify(d, *t));
            all_ok &= ok;
            if let Some(d) = d {
                agree &= agreed.is_none_or(|a| a == d);
                agreed = Some(d);
            }
            row[&arm.name] = json!({"recovered_scalar": d.map(|d| d.to_string()), "generic_replay": ok});
        }
        all_ok &= agree;
        row["arms_agree"] = json!(agree);
        if let Some(reference) = &reference {
            let same = agreed == reference.get(i).copied();
            all_ok &= same;
            row["matches_hosted_reference"] = json!(same);
        }
        verification.push(row);
    }
    let mut analysis = json!({
        "schema": "n83-precomputed-rho-analysis/v1",
        "job": job,
        "ic_label": ic_label,
        "targets": targets.len(),
        "all_scalars_verified_and_agreeing": all_ok,
        "hosted_reference_checked": reference.is_some(),
        "verification": verification,
        "bootstrap": {"resamples": RESAMPLES, "interval": "95% percentile, targets resampled with replacement", "seeds": "fixed per statistic"},
        "arms": {
            "ic": arm_stats(&ic, 11),
            "rho_plain_fold": arm_stats(&rho, 12),
            "bl_a_time_matched": arm_stats(&bl_a, 13),
            "bl_b_memory_matched": arm_stats(&bl_b, 14),
        },
        "setup": {"ic": ic_info, "bl_a_time_matched": a_info, "bl_b_memory_matched": b_info},
        "paired": {
            "rho_plain_over_ic": paired(&rho, &ic, 21),
            "ic_over_bl_a": paired(&ic, &bl_a, 22),
            "ic_over_bl_b": paired(&ic, &bl_b, 23),
            "rho_plain_over_bl_a": paired(&rho, &bl_a, 24),
        },
    });
    if let Some((a2_arm, a2_info)) = &a2 {
        analysis["arms"]["bl_a2_wall_equal"] = arm_stats(a2_arm, 15);
        analysis["setup"]["bl_a2_wall_equal"] = a2_info.clone();
        analysis["paired"]["ic_over_bl_a2"] = paired(&ic, a2_arm, 25);
    }
    fs::write(out, format!("{}\n", serde_json::to_string_pretty(&analysis).unwrap()))
        .map_err(|e| e.to_string())?;
    Ok(all_ok)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 {
        eprintln!("usage: n83_precomputed_rho_analyze <job-dir> <ic-label> <analysis.json>");
        std::process::exit(2);
    }
    match run(&args[1], &args[2], &args[3]) {
        Ok(true) => println!("PASS: every recovered scalar verified and all arms agree"),
        Ok(false) => {
            println!("FAIL: a scalar failed replay or arms disagree");
            std::process::exit(1);
        }
        Err(e) => {
            eprintln!("error: {e}");
            std::process::exit(2);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn median_and_percentiles() {
        assert_eq!(median(&[3.0, 1.0, 2.0]), 2.0);
        assert_eq!(median(&[4.0, 1.0, 2.0, 3.0]), 2.5);
        let ci = percentile_ci((0..=1000).map(|i| i as f64).collect());
        assert_eq!(ci, [25.0, 975.0]);
    }

    #[test]
    fn bootstrap_of_constant_is_constant() {
        let v = vec![5.0; 16];
        assert_eq!(bootstrap(16, 1, |i| mean(&pick(&v, i))), [5.0, 5.0]);
    }

    #[test]
    fn generic_replay_accepts_the_answer_and_rejects_others() {
        let curve = Curve { curve: KoblitzCurve::new(1, 83).unwrap() };
        // Hosted workload 0 of gha-38597929540-1 (public target, both arms' answer).
        let q = [6274346871052084705058109u128, 7777137791885914859810330];
        assert!(curve.verify(4565490491554819, q));
        assert!(!curve.verify(4565490491554818, q));
    }

    #[test]
    fn gnu_time_rss_parses() {
        let dir = std::env::temp_dir().join(format!("n83bl-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let p = dir.join("time.txt");
        fs::write(&p, "\tMaximum resident set size (kbytes): 1024\n").unwrap();
        assert_eq!(time_max_rss(p.to_str().unwrap()), Some(1024 * 1024));
    }
}

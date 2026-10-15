//! Read one frozen native campaign. Never pool processors or accept censored timings.
use isogeny_algos::{json::Json, sha256::sha256_hex};
use std::collections::BTreeMap;
use std::io::Write;

fn median(xs: &[f64]) -> Result<f64, String> {
    if xs.len() != 5 || xs.iter().any(|x| !x.is_finite() || *x <= 0.0) {
        return Err("expected five finite positive round ratios".into());
    }
    let mut values = xs.to_vec();
    values.sort_by(f64::total_cmp);
    Ok(values[2])
}
fn summarize(j: &Json) -> Result<Json, String> {
    if j.get("schema").and_then(Json::as_str) != Some("isogeny-perf.campaign/v1") {
        return Err("unknown campaign schema".into());
    }
    if j.get("status").and_then(Json::as_str) != Some("PASS")
        || j.get("exact_outputs_all_runs") != Some(&Json::Bool(true))
        || j.get("repeatable_counts_each_binary") != Some(&Json::Bool(true))
        || j.get("wall_speedup_eligible") != Some(&Json::Bool(true))
    {
        return Ok(Json::obj(vec![
            ("schema", Json::str("isogeny-perf.cpu-analysis/v1")),
            ("status", Json::str("CENSORED")),
            ("qualified_wall_results", Json::Arr(vec![])),
            ("reason", Json::str("campaign lacks complete uncontended timing, exact output or repeated-count receipts")),
        ]));
    }
    let Some(Json::Arr(pairs)) = j.get("qualified_pairs") else {
        return Err("missing qualified pairs".into());
    };
    let mut groups: BTreeMap<(String, String, String), (BTreeMap<i64, f64>, BTreeMap<i64, f64>)> =
        BTreeMap::new();
    for pair in pairs {
        let text = |key| {
            pair.get(key)
                .and_then(Json::as_str)
                .ok_or_else(|| format!("missing {key}"))
        };
        if !["192", "224", "256", "511"].contains(&text("bits")?)
            || !["7", "31"].contains(&text("ell")?)
            || !isogeny_algos::find::bmss::Method::ALL
                .iter()
                .any(|m| Some(m.name()) == pair.get("algorithm").and_then(Json::as_str))
        {
            return Err("unexpected matrix cell".into());
        }
        let key = (
            text("bits")?.into(),
            text("ell")?.into(),
            text("algorithm")?.into(),
        );
        let round = pair
            .get("round")
            .and_then(Json::as_f64)
            .ok_or("missing round")?;
        if !round.is_finite() || round.fract() != 0.0 || !(0.0..5.0).contains(&round) {
            return Err("invalid round".into());
        }
        let ratio = pair
            .get("a_over_b")
            .and_then(Json::as_f64)
            .ok_or("missing ratio")?;
        let a = pair
            .get("a_median_ns")
            .and_then(Json::as_f64)
            .ok_or("missing baseline median")?;
        let b = pair
            .get("b_median_ns")
            .and_then(Json::as_f64)
            .ok_or("missing candidate median")?;
        if !a.is_finite()
            || !b.is_finite()
            || a <= 0.0
            || b <= 0.0
            || (ratio - a / b).abs() > 1e-12 * (a / b)
        {
            return Err("inconsistent or nonpositive sample medians".into());
        }
        let arm = groups.entry(key).or_default();
        let rows = match text("phase")? {
            "aa" => &mut arm.0,
            "ab" => &mut arm.1,
            _ => return Err("unknown phase".into()),
        };
        if rows.insert(round as i64, ratio).is_some() {
            return Err("duplicate round".into());
        }
    }
    if groups.len() != 64 {
        return Err("missing requested matrix cells".into());
    }
    let mut results = vec![];
    for ((bits, ell, algorithm), (aa, ab)) in groups {
        let aa: Vec<_> = aa.values().copied().collect();
        let ab: Vec<_> = ab.values().copied().collect();
        let aa_median = median(&aa)?;
        let ab_median = median(&ab)?;
        let minimum = |xs: &[f64]| xs.iter().copied().fold(f64::INFINITY, f64::min);
        let maximum = |xs: &[f64]| xs.iter().copied().fold(0.0, f64::max);
        // A/A can move in either direction. Keep that observed envelope;
        // this threshold is descriptive, not a statistical confidence interval.
        let noise = maximum(&aa).max(1.0 / minimum(&aa));
        results.push(Json::obj(vec![
            ("bits", Json::str(bits)),
            ("ell", Json::str(ell)),
            ("algorithm", Json::str(algorithm)),
            (
                "aa_round_ratios",
                Json::Arr(aa.iter().copied().map(Json::Float).collect()),
            ),
            (
                "ab_round_ratios",
                Json::Arr(ab.iter().copied().map(Json::Float).collect()),
            ),
            ("aa_median", Json::Float(aa_median)),
            ("aa_min", Json::Float(minimum(&aa))),
            ("aa_max", Json::Float(maximum(&aa))),
            ("ab_median", Json::Float(ab_median)),
            ("ab_min", Json::Float(minimum(&ab))),
            ("ab_max", Json::Float(maximum(&ab))),
            ("observed_noise_envelope", Json::Float(noise)),
            ("median_exceeds_aa_envelope", Json::Bool(ab_median > noise)),
        ]));
    }
    Ok(Json::obj(vec![
        ("schema", Json::str("isogeny-perf.cpu-analysis/v1")), ("status", Json::str("PASS")),
        ("unit", Json::str("per-case baseline/candidate ratio of constructor sample medians, five paired rounds")),
        ("qualification", Json::str("native reserved-core receipt; guest-level contention excluded; VM host tenants and frequency remain unobserved")),
        ("qualified_wall_results", Json::Arr(results)),
    ]))
}
fn main() -> Result<(), String> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 3 {
        return Err("summarize-cpu SUMMARY.json FRESH_ANALYSIS.json".into());
    }
    let bytes = std::fs::read(&args[1]).map_err(|e| e.to_string())?;
    let data = Json::parse(std::str::from_utf8(&bytes).map_err(|e| e.to_string())?)?;
    let mut result = summarize(&data)?;
    let Json::Obj(fields) = &mut result else {
        unreachable!()
    };
    fields.push((
        "source_summary_sha256".into(),
        Json::str(sha256_hex(&bytes)),
    ));
    let mut out = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&args[2])
        .map_err(|e| e.to_string())?;
    writeln!(out, "{}", result.dump()).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn shared_host_is_censored_and_missing_rounds_fail() {
        let j = Json::obj(vec![
            ("schema", Json::str("isogeny-perf.campaign/v1")),
            ("status", Json::str("PASS")),
        ]);
        assert_eq!(
            summarize(&j).unwrap().get("status").and_then(Json::as_str),
            Some("CENSORED")
        );
        assert!(median(&[1.0, 2.0]).is_err());
        assert_eq!(median(&[1.0, 1.1, 0.9, 1.05, 1.02]).unwrap(), 1.02);
        assert!(median(&[1.0, 1.0, 1.0, 1.0, f64::NAN]).is_err());
    }
}

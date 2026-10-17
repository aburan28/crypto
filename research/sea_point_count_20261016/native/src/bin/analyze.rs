//! Validate paired outputs before computing timing ratios. Never admit contention.
use isogeny_algos::{
    api::{self, PrimeCurve},
    bigint::Big,
    int::Int,
    json::Json,
};
use std::{collections::BTreeMap, fs, path::Path};

fn key<'a>(j: &'a Json, k: &str) -> &'a Json {
    j.get(k).unwrap_or_else(|| panic!("missing {k}"))
}
fn number(j: &Json) -> f64 {
    j.as_f64()
        .or_else(|| j.as_str().and_then(|s| s.parse().ok()))
        .expect("number")
}
fn median(v: &[f64]) -> f64 {
    let mut s = v.to_vec();
    s.sort_by(f64::total_cmp);
    if s.len() % 2 == 1 {
        s[s.len() / 2]
    } else {
        (s[s.len() / 2 - 1] + s[s.len() / 2]) / 2.
    }
}
fn valid_pair(a: &Json, b: &Json) -> Result<(), String> {
    for row in [a, b] {
        if key(row, "exit_code").as_f64() != Some(0.) {
            return Err("nonzero exit".into());
        }
        if key(key(row, "result"), "verified") != &Json::Bool(true) {
            return Err("unverified result".into());
        }
        for field in ["calls", "elapsed_ns"] {
            let value = number(key(key(row, "result"), field));
            if !value.is_finite() || value <= 0. {
                return Err(format!("invalid positive measurement: {field}"));
            }
        }
    }
    for name in ["curve", "order", "calls", "seed_schedule"] {
        if key(key(a, "result"), name) != key(key(b, "result"), name) {
            return Err(format!("mismatched {name}"));
        }
    }
    Ok(())
}
fn per_call(row: &Json) -> f64 {
    number(key(key(row, "result"), "elapsed_ns")) / number(key(key(row, "result"), "calls")) / 1000.
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 3 {
        return Err("usage: analyze RUN_DIRECTORY REGISTRY_JSON".into());
    }
    let dir = Path::new(&args[1]);
    let registry = Json::parse(&fs::read_to_string(&args[2])?)?;
    let Json::Arr(curves) = key(&registry, "curves") else {
        return Err("registry curves".into());
    };
    let mut groups: BTreeMap<(String, String, i64), Vec<Json>> = BTreeMap::new();
    for line in fs::read_to_string(dir.join("runs.jsonl"))?.lines() {
        let row = Json::parse(line)?;
        let slug = key(key(&row, "result"), "curve");
        if !curves.iter().any(|c| c.get("slug") == Some(slug)) {
            return Err(format!("unregistered model: {}", slug.dump()).into());
        }
        let group = (
            key(&row, "kind").as_str().unwrap().into(),
            key(&row, "case").as_str().unwrap().into(),
            number(key(&row, "round")) as i64,
        );
        groups.entry(group).or_default().push(row);
    }
    let mut summaries = vec![];
    let mut certificates = vec![];
    let mut md=String::from("| Workload | A median us/call | B median us/call | A/B median paired ratio | Paired range | A/A max deviation | Admitted A/B pairs |\n| --- | ---: | ---: | ---: | --- | ---: | ---: |\n");
    for case in [
        "bsgs24",
        "bsgs32",
        "bsgs40",
        "sea40_cached",
        "sea192_cached",
        "count192_cold",
        "cm256_control",
    ] {
        let mut aa = vec![];
        let mut av = vec![];
        let mut bv = vec![];
        let mut ratios = vec![];
        let mut excluded = 0;
        let mut assetups = vec![];
        let mut bsetups = vec![];
        let mut curve = Json::Null;
        let mut signature: Option<Vec<Json>> = None;
        for kind in ["aa", "ab"] {
            for round in 0..5 {
                let rows = groups
                    .get_mut(&(kind.into(), case.into(), round))
                    .ok_or("missing pair")?;
                rows.sort_by_key(|r| number(key(r, "slot")) as i64);
                if rows.len() != 2 {
                    return Err("pair must have two rows".into());
                }
                valid_pair(&rows[0], &rows[1])?;
                for row in rows.iter() {
                    let result = key(row, "result");
                    if key(result, "case").as_str() != Some(case) {
                        return Err("result workload differs from requested workload".into());
                    }
                    let current: Vec<Json> = ["curve", "order", "calls", "seed_schedule"]
                        .into_iter()
                        .map(|name| key(result, name).clone())
                        .collect();
                    if signature.as_ref().is_some_and(|s| s != &current) {
                        return Err("workload changed between rounds".into());
                    }
                    signature = Some(current);
                }
                curve = key(key(&rows[0], "result"), "curve").clone();
                if rows
                    .iter()
                    .any(|r| key(r, "contended") != &Json::Bool(false))
                {
                    excluded += 1;
                    continue;
                }
                let (a, b) = (per_call(&rows[0]), per_call(&rows[1]));
                if kind == "aa" {
                    aa.push(a / b)
                } else {
                    av.push(a);
                    bv.push(b);
                    ratios.push(a / b);
                    assetups.push(number(key(key(&rows[0], "result"), "setup_ns")) / 1e6);
                    bsetups.push(number(key(key(&rows[1], "result"), "setup_ns")) / 1e6);
                }
            }
        }
        if aa.is_empty() || ratios.is_empty() {
            return Err(format!("no admitted pairs for {case}; retain raw evidence").into());
        }
        let noise = aa.iter().map(|r| (r - 1.).abs()).fold(0., f64::max);
        let low = ratios.iter().copied().fold(f64::INFINITY, f64::min);
        let high = ratios.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let ma = median(&av);
        let mb = median(&bv);
        let ratio = median(&ratios);
        // Independent from SEA congruence recombination and the BSGS search:
        // prove the registered reference via point-order divisors and Hasse.
        // This runs after timing, so validation cost is not silently charged
        // to one arm or used to disturb a timed sample.
        let entry = curves
            .iter()
            .find(|c| c.get("slug") == Some(&curve))
            .unwrap();
        if key(entry, "order") != &signature.as_ref().unwrap()[1] {
            return Err("measured answer differs from registered reference order".into());
        }
        let model = Json::parse(key(entry, "model_json").as_str().unwrap())?;
        let integer = |v: &Json| Int::from_big(&Big::from_dec(v.as_str().unwrap()));
        let c = PrimeCurve::new(
            integer(key(&model, "p")),
            integer(key(&model, "a")),
            integer(key(&model, "b")),
        )?;
        let order = integer(key(entry, "order"));
        let checked = api::check_order(&c, &order, 1);
        if !checked.certificate.certified || !checked.twist_certificate.certified {
            return Err(format!("reference order could not be certified for {case}").into());
        }
        certificates.push(Json::obj(vec![
            ("case", Json::str(case)),
            ("curve", curve.clone()),
            ("order", Json::str(&order)),
            ("curve_certified", Json::Bool(checked.certificate.certified)),
            (
                "twist_certified",
                Json::Bool(checked.twist_certificate.certified),
            ),
            ("curve_divisor", Json::str(checked.certificate.divisor)),
            (
                "twist_divisor",
                Json::str(checked.twist_certificate.divisor),
            ),
            ("via_twist", Json::Bool(checked.certificate.via_twist)),
            (
                "assumption",
                Json::str("factor primality is checked by probable-prime tests"),
            ),
        ]));
        md.push_str(&format!(
            "| {case} | {ma:.3} | {mb:.3} | {ratio:.5} | {low:.5}..{high:.5} | {:.3}% | {} / 5 |\n",
            100. * noise,
            ratios.len()
        ));
        summaries.push(Json::obj(vec![
            ("case", Json::str(case)),
            ("curve", curve),
            ("a_median_us", Json::Float(ma)),
            ("b_median_us", Json::Float(mb)),
            (
                "a_min_us",
                Json::Float(av.iter().copied().fold(f64::INFINITY, f64::min)),
            ),
            (
                "b_min_us",
                Json::Float(bv.iter().copied().fold(f64::INFINITY, f64::min)),
            ),
            ("ratio_median", Json::Float(ratio)),
            ("ratio_min", Json::Float(low)),
            ("ratio_max", Json::Float(high)),
            ("aa_max_deviation", Json::Float(noise)),
            ("aa_pairs", Json::Num(aa.len() as i64)),
            ("ab_pairs", Json::Num(ratios.len() as i64)),
            ("excluded_pairs", Json::Num(excluded)),
            ("a_setup_median_ms", Json::Float(median(&assetups))),
            ("b_setup_median_ms", Json::Float(median(&bsetups))),
            (
                "outside_aa_envelope",
                Json::Bool((ratio - 1.).abs() > noise),
            ),
        ]));
    }
    fs::write(dir.join("summary.json"), Json::Arr(summaries).dump())?;
    fs::write(
        dir.join("reference-certificates.json"),
        Json::Arr(certificates).dump(),
    )?;
    fs::write(dir.join("summary.md"), &md)?;
    print!("{md}");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn row(order: &str, calls: i64, verified: bool, exit: i64) -> Json {
        Json::obj(vec![
            ("exit_code", Json::Num(exit)),
            (
                "result",
                Json::obj(vec![
                    ("curve", Json::str("fixture")),
                    ("order", Json::str(order)),
                    ("calls", Json::Num(calls)),
                    ("seed_schedule", Json::str("1 + i % 4")),
                    ("elapsed_ns", Json::str("1000")),
                    ("verified", Json::Bool(verified)),
                ]),
            ),
        ])
    }
    #[test]
    fn rejects_wrong_answers_different_work_and_failed_runs() {
        let a = row("101", 5, true, 0);
        assert!(valid_pair(&a, &a).is_ok());
        for b in [
            row("102", 5, true, 0),
            row("101", 6, true, 0),
            row("101", 5, false, 0),
            row("101", 5, true, 1),
        ] {
            assert!(valid_pair(&a, &b).is_err());
        }
        let zero_calls = row("101", 0, true, 0);
        assert!(valid_pair(&zero_calls, &zero_calls).is_err());
    }
}

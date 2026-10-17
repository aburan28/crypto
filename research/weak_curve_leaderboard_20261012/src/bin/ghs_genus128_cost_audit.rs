//! Exact GF(4) monic-polynomial counts and conditional genus-128 cost replay.
//! This prices a model; it does not construct divisor places on a GHS cover.

use num_bigint::{BigInt, BigUint, Sign};
use num_traits::{One, ToPrimitive, Zero};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{env, error::Error, fs};

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const G: usize = 128;
const Q: u32 = 4;
const TRANSFER: &str = "../../docs/curves/jacobian-transfers.json";
const COST_SOURCE: &str = "../../src/cryptanalysis/jacobian_transfer/cost.rs";
const OUTPUT: &str = "ghs_genus128_cost_audit.json";

fn hash(path: &str) -> Result<String> {
    Ok(format!("{:x}", Sha256::digest(fs::read(path)?)))
}

fn mobius(mut n: u32) -> i32 {
    let mut result = 1;
    let mut prime = 2;
    while prime * prime <= n {
        if n % prime == 0 {
            n /= prime;
            if n % prime == 0 {
                return 0;
            }
            result = -result;
        }
        prime += 1;
    }
    if n > 1 {
        -result
    } else {
        result
    }
}

fn irreducible_count(degree: u32) -> BigUint {
    let mut numerator = BigInt::zero();
    for divisor in 1..=degree {
        if degree % divisor != 0 {
            continue;
        }
        let mu = mobius(divisor);
        if mu == 0 {
            continue;
        }
        let power = BigInt::from_biguint(Sign::Plus, BigUint::from(Q).pow(degree / divisor));
        numerator += BigInt::from(mu) * power;
    }
    assert_eq!((&numerator % degree).to_i32(), Some(0));
    (numerator / degree).to_biguint().unwrap()
}

fn log2_big(value: &BigUint) -> f64 {
    value.to_f64().expect("counts fit in finite f64").log2()
}

fn rounded(value: f64) -> f64 {
    (value * 1_000_000.0).round() / 1_000_000.0
}

fn logadd(left: f64, right: f64) -> f64 {
    let hi = left.max(right);
    let lo = left.min(right);
    hi + (1.0 + (lo - hi).exp2()).log2()
}

fn combinations_with_repetition(items: &BigUint, k_max: usize) -> Vec<BigUint> {
    let mut out = Vec::with_capacity(k_max + 1);
    out.push(BigUint::one());
    for k in 1..=k_max {
        let next = (&out[k - 1] * (items + BigUint::from(k - 1))) / BigUint::from(k);
        out.push(next);
    }
    out
}

fn advance_counts(previous: &[BigUint], degree: usize, irreducibles: &BigUint) -> Vec<BigUint> {
    let combos = combinations_with_repetition(irreducibles, G / degree);
    (0..=G)
        .map(|target| {
            (0..=target / degree).fold(BigUint::zero(), |acc, k| {
                acc + &previous[target - k * degree] * &combos[k]
            })
        })
        .collect()
}

fn sweep(model: &Value, subgroup_bits: f64) -> Value {
    let mut counts = vec![BigUint::zero(); G + 1];
    counts[0] = BigUint::one();
    let mut total_irreducibles = BigUint::zero();
    let exponent = model["karatsuba_exponent"].as_f64().unwrap();
    let group_op = model["jacobian_group_op_mults_per_genus_squared"]
        .as_f64()
        .unwrap()
        .log2()
        + 2.0 * (G as f64).log2()
        + exponent * (Q as f64).log2().log2();
    let linear_op = exponent * subgroup_bits.log2();
    let mut rows = Vec::new();
    let mut selected = BigUint::zero();
    let mut selected_base = BigUint::zero();
    for degree in 1..=G {
        let irreducibles = irreducible_count(degree as u32);
        total_irreducibles += &irreducibles;
        counts = advance_counts(&counts, degree, &irreducibles);
        if degree == 15 {
            selected = counts[G].clone();
            selected_base = total_irreducibles.clone();
        }
        if !counts[G].is_zero() {
            let log2_base = log2_big(&total_irreducibles) - 1.0;
            let log2_smooth = log2_big(&counts[G]) - (G as f64) * (Q as f64).log2();
            let log2_relations = log2_base - log2_smooth;
            let log2_linear = 2.0 * log2_base + (G as f64).log2();
            let total = logadd(log2_relations + group_op, log2_linear + linear_op);
            rows.push(json!({
                "degree_bound":degree,
                "log2_nominal_base":rounded(log2_base),
                "log2_smooth_fraction":rounded(log2_smooth),
                "log2_relation_queries":rounded(log2_relations),
                "log2_linear_work":rounded(log2_linear),
                "log2_total":rounded(total)
            }));
        }
        if degree == G || log2_big(&total_irreducibles) - 1.0 >= 160.0 {
            break;
        }
    }
    let best = rows
        .iter()
        .min_by(|a, b| {
            a["log2_total"]
                .as_f64()
                .unwrap()
                .total_cmp(&b["log2_total"].as_f64().unwrap())
        })
        .unwrap()
        .clone();
    assert_eq!(best["degree_bound"], 15);
    assert_eq!(selected_base, BigUint::from(97_842_104u64));
    assert_eq!(
        selected.to_string(),
        "1374288545819365998587940548442445318355469588992205247403709040067002"
    );
    let mut best = best;
    best["nominal_half_base_places"] = json!(48_921_052u64);
    best["smooth_polynomial_count"] = json!(selected.to_string());
    best["all_monic_polynomials"] = json!(BigUint::from(Q).pow(G as u32).to_string());
    json!({"tested_degree_bounds":rows.len(), "best":best, "sweep":rows})
}

fn build() -> Result<Value> {
    let source: Value = serde_json::from_slice(&fs::read(TRANSFER)?)?;
    let model = &source["cost_model"];
    let mut curves = Vec::new();
    for name in ["c2pnb176w1", "c2pnb272w1"] {
        let curve = source["curves"]
            .as_array()
            .unwrap()
            .iter()
            .find(|curve| {
                curve["standard_names"]
                    .as_array()
                    .is_some_and(|names| names.iter().any(|v| v == name))
            })
            .expect("named curve");
        let route = curve["routes"]
            .as_array()
            .unwrap()
            .iter()
            .find(|route| {
                route["kind"] == "ghs_weil_descent" && route["invariants"]["base_degree"] == 2
            })
            .expect("GF(4) GHS route");
        let inv = &route["invariants"];
        let dest = &route["destination"];
        let cost = &route["cost"];
        assert_eq!(route["evidence"], "exact_invariant");
        assert_eq!(route["class"], "formula_advantage");
        assert_eq!(route["transport"], "size_admits");
        assert_eq!(inv["magic_number"], 8);
        assert_eq!(inv["type_i"], false);
        assert_eq!(inv["mmt_condition_8"], true);
        assert_eq!(inv["relative_trace_b_zero"], false);
        assert_eq!(dest["genus"], "128");
        assert_eq!(dest["field_bits"], 2.0);
        assert_eq!(cost["destination"]["solver"], "enge_gaudry_smoothness");
        assert_eq!(cost["destination"]["enge_gaudry_hypothesis_holds"], true);
        let audit = sweep(model, curve["subgroup"]["bits"].as_f64().unwrap());
        let delta = (audit["best"]["log2_total"].as_f64().unwrap()
            - cost["destination"]["log2_cost"].as_f64().unwrap())
        .abs();
        assert!(delta < 0.001, "source formula disagreement: {delta}");
        curves.push(json!({
            "name":name, "curve_id":curve["icv1"],
            "subgroup_order":curve["subgroup"]["order"],
            "subgroup_bits":curve["subgroup"]["bits"],
            "relative_degree":inv["relative_degree"],
            "reference_rho_log2":curve["reference_log2"],
            "source_formula_log2":cost["destination"]["log2_cost"],
            "audit":audit,
            "cover_model":Value::Null, "selected_subgroup_image":Value::Null,
            "actual_usable_base":Value::Null, "verified_relation_yield":Value::Null,
            "useful_rank":Value::Null, "target_descent_cost":Value::Null
        }));
    }
    Ok(json!({
        "schema":"weakcurves.ghs_genus128_cost_audit.v1",
        "status":"conditional_cost_model_replayed",
        "producer_language":"rust",
        "field_size":Q, "genus":G,
        "model_assumptions":[
            "half of base irreducibles yield usable places",
            "monic-polynomial smoothness represents divisor smoothness",
            "reported operation prices and sparse linear algebra apply"
        ],
        "source_sha256":{
            "docs/curves/jacobian-transfers.json":hash(TRANSFER)?,
            "src/cryptanalysis/jacobian_transfer/cost.rs":hash(COST_SOURCE)?,
            "research/weak_curve_leaderboard_20261012/src/bin/ghs_genus128_cost_audit.rs":
                hash("src/bin/ghs_genus128_cost_audit.rs")?
        },
        "curves":curves
    }))
}

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();
    let value = build()?;
    match args.get(1).map(String::as_str) {
        Some("--write") if args.len() == 2 => {
            fs::write(
                OUTPUT,
                format!("{}\n", serde_json::to_string_pretty(&value)?),
            )?;
        }
        Some("--check") if args.len() == 2 => {
            let saved: Value = serde_json::from_slice(&fs::read(OUTPUT)?)?;
            assert_eq!(value, saved, "GHS genus-128 model receipt drifted");
        }
        _ => return Err("usage: ghs_genus128_cost_audit --write|--check".into()),
    }
    println!("GHS genus-128 cost model: two routes; b=15; 48,921,052 nominal half-base places");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn mobius_and_irreducibles() {
        assert_eq!(mobius(1), 1);
        assert_eq!(mobius(4), 0);
        assert_eq!(mobius(6), 1);
        assert_eq!(irreducible_count(1), BigUint::from(4u32));
        assert_eq!(irreducible_count(2), BigUint::from(6u32));
    }
    #[test]
    fn genus128_model_replays_registry() {
        let value = build().unwrap();
        assert_eq!(value["curves"][0]["audit"]["best"]["degree_bound"], 15);
        assert_eq!(value["curves"][1]["audit"]["best"]["degree_bound"], 15);
        assert_eq!(value["curves"][0]["audit"]["tested_degree_bounds"], 84);
    }
}

//! Bounded six-summand SAT access control on two S3-derived planted witnesses.
//! Open-stage models are field models, not measured natural relations.

use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::binary_semaev_chain_sat::ChainedS3Encoding;
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use crypto_lib::cryptanalysis::sat::SolveResult;
use num_bigint::BigUint;
use serde_json::{json, Value};
use std::error::Error;
use std::time::Instant;

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const FIXTURE: &[u8] = include_bytes!("../research/koblitz_n83_factor_base_sweep_20261008/verification/hosted-m6-witness-ladder-20261017/fixture.json");
const FIXTURE_BLAKE3: &str = "f07b515f10b4cfd0b0d82267a22c1d46018d278fc922de7cd7b73c77e61ed43a";
const CURVE: &str = "icv1-f2m83-tm6151469093347-debefd74";

fn decimal(value: &Value) -> Result<BigUint> {
    let raw = value.as_str().ok_or("decimal string required")?;
    if raw.is_empty() || !raw.bytes().all(|byte| byte.is_ascii_digit()) {
        return Err("invalid decimal string".into());
    }
    let parsed: BigUint = raw.parse()?;
    if parsed.to_string() != raw || parsed.bits() > 83 {
        return Err("noncanonical field integer".into());
    }
    Ok(parsed)
}

fn decimals(value: &Value, count: usize) -> Result<Vec<BigUint>> {
    let array = value.as_array().ok_or("array required")?;
    if array.len() != count {
        return Err("array length differs".into());
    }
    array.iter().map(decimal).collect()
}

fn point(value: &Value) -> Result<BinaryPoint> {
    let parts = value.as_array().ok_or("affine point required")?;
    if parts.len() != 2 {
        return Err("affine point arity differs".into());
    }
    Ok(BinaryPoint::Affine {
        x: F2mElement::from_biguint(&decimal(&parts[0])?, 83),
        y: F2mElement::from_biguint(&decimal(&parts[1])?, 83),
    })
}

fn x_of(point: &BinaryPoint) -> Result<BigUint> {
    let BinaryPoint::Affine { x, .. } = point else {
        return Err("identity in planted witness".into());
    };
    Ok(x.to_biguint())
}

fn pin_block(enc: &mut ChainedS3Encoding, block: usize, value: &BigUint) {
    for bit in 0..83 {
        let wire = (block * 83 + bit + 1) as i32;
        enc.solver
            .add_clause(vec![if value.bit(bit as u64) { wire } else { -wire }]);
    }
}

fn run(role: &str, stage: &str, conflict_cap: u64) -> Result<Value> {
    if !["gray_pool", "hash_pool"].contains(&role)
        || !["all_x", "summand_x", "open"].contains(&stage)
        || !(1..=20_000).contains(&conflict_cap)
        || blake3::hash(FIXTURE).to_hex().as_str() != FIXTURE_BLAKE3
    {
        return Err("frozen control argument or fixture differs".into());
    }
    let fixture: Value = serde_json::from_slice(FIXTURE)?;
    if fixture["schema"] != "n83.m6-planted-witness-ladder/v1"
        || fixture["status"] != "PASS_S3_BOUND_GROUP_WITNESSES"
        || fixture["curve_slug"] != CURVE
        || fixture["summands"] != 6
        || fixture["base_summands"] != 4
        || fixture["residual_summands"] != 2
        || !fixture["natural_relations"].is_null()
        || !fixture["selected_best_total_runtime"].is_null()
    {
        return Err("fixture contract differs".into());
    }
    let arms = fixture["arms"].as_array().ok_or("arms array required")?;
    let selected: Vec<&Value> = arms.iter().filter(|arm| arm["role"] == role).collect();
    if arms.len() != 2 || selected.len() != 1 {
        return Err("arm identity differs".into());
    }
    let arm = selected[0];
    let base_count = 600;
    let residual_count = if role == "gray_pool" { 97 } else { 600 };
    if arm["identity_mask"] != 0
        || arm["residual_signed_orbits"].as_u64() != Some(residual_count as u64)
    {
        return Err("arm domain shape differs".into());
    }
    let base = decimals(&arm["base_reps"], base_count)?;
    let residual = decimals(&arm["residual_reps"], residual_count)?;
    let xs = decimals(&arm["ordered_xs"], 6)?;
    let partial_xs = decimals(&arm["partial_xs"], 4)?;
    let target_x = decimal(&arm["target_x"])?;
    if xs[0] >= xs[1]
        || xs[..4].iter().any(|x| !base.contains(x))
        || xs[4..].iter().any(|x| !residual.contains(x))
    {
        return Err("selected x coordinates differ from domains".into());
    }
    let selected_points = arm["selected_points"].as_array().ok_or("selected points")?;
    if selected_points.len() != 6 {
        return Err("selected point count differs".into());
    }
    let curve = KoblitzCurve::known_n83_k0().ok_or("pinned N83 curve")?;
    if curve.label() != CURVE {
        return Err("curve identity differs".into());
    }
    let mut sum = BinaryPoint::Infinity;
    for (index, raw) in selected_points.iter().enumerate() {
        let term = point(raw)?;
        if !curve.curve.is_on_curve(&term)
            || curve.mul(&term, &curve.subgroup_order) != BinaryPoint::Infinity
            || x_of(&term)? != xs[index]
        {
            return Err("selected point failed original-group replay".into());
        }
        sum = curve.add(&sum, &term);
        let current_x = x_of(&sum)?;
        if (1..5).contains(&index) && current_x != partial_xs[index - 1] {
            return Err("partial group sum differs".into());
        }
    }
    if sum != point(&arm["target"])? || x_of(&sum)? != target_x {
        return Err("planted target differs from group sum".into());
    }
    let mut domains = vec![base.clone(), base.clone(), base.clone(), base];
    domains.push(residual.clone());
    domains.push(residual);
    let started = Instant::now();
    eprintln!(
        "{}",
        json!({"phase":"build_started","elapsed_wall_ms_L0":0.0})
    );
    let target = F2mElement::from_biguint(&target_x, 83);
    let mut enc = ChainedS3Encoding::build(
        6,
        &curve.curve.irreducible,
        &curve.curve.b,
        Some(&target),
        0,
        150_000,
    )
    .map_err(|error| format!("circuit build: {error:?}"))?;
    eprintln!(
        "{}",
        json!({"phase":"circuit_built",
        "elapsed_wall_ms_L0":started.elapsed().as_secs_f64()*1000.0,
        "variables":enc.solver.n_vars(), "clauses":enc.solver.n_clauses(),
        "xors":enc.solver.n_xors()})
    );
    let domain_size = enc
        .constrain_frobenius_orbit_domains(&domains, 150_000, 12_000_000)
        .map_err(|error| format!("domain build: {error:?}"))?;
    let expected_domain_clauses = if role == "gray_pool" {
        237_280
    } else {
        323_796
    };
    if domain_size.clauses != expected_domain_clauses {
        return Err("domain count differs from frozen natural-search arm".into());
    }
    enc.order_first_pair(150_000)
        .map_err(|error| format!("first-pair order: {error:?}"))?;
    eprintln!(
        "{}",
        json!({"phase":"domains_built",
        "elapsed_wall_ms_L0":started.elapsed().as_secs_f64()*1000.0,
        "variables":enc.solver.n_vars(), "clauses":enc.solver.n_clauses(),
        "xors":enc.solver.n_xors(), "domain_clauses":domain_size.clauses})
    );
    if stage != "open" {
        for (index, x) in xs.iter().enumerate() {
            pin_block(&mut enc, index, x);
        }
    }
    if stage == "all_x" {
        for (index, x) in partial_xs.iter().enumerate() {
            pin_block(&mut enc, 6 + index, x);
        }
    }
    let fixed_x_blocks = match stage {
        "all_x" => 10,
        "summand_x" => 6,
        _ => 0,
    };
    eprintln!(
        "{}",
        json!({"phase":"fixed_units_installed",
        "elapsed_wall_ms_L0":started.elapsed().as_secs_f64()*1000.0,
        "fixed_x_blocks":fixed_x_blocks, "variables":enc.solver.n_vars(),
        "clauses":enc.solver.n_clauses(), "xors":enc.solver.n_xors()})
    );
    enc.solver.conflict_budget = conflict_cap;
    let outcome = enc.solver.solve();
    let model_verified = if outcome == SolveResult::Sat {
        enc.verify_model()
            && (stage == "open"
                || enc
                    .decode_summands()
                    .iter()
                    .zip(&xs)
                    .all(|(observed, expected)| observed.to_biguint() == *expected))
    } else {
        false
    };
    Ok(json!({
        "schema":"n83.m6-planted-witness-control/v1", "scope":"known_answer_solver_access_only",
        "role":role, "stage":stage,
        "status":match outcome {
            SolveResult::Sat if model_verified && stage == "open" => "SAT_verified_field_model",
            SolveResult::Sat if model_verified => "SAT_verified_planted_x",
            SolveResult::Sat => "PRODUCER_FAILURE_invalid_model",
            SolveResult::Unsat => "PRODUCER_FAILURE_unsat_planted_encoding",
            SolveResult::Unknown => "UNKNOWN_conflict_cap",
        },
        "fixture_blake3":FIXTURE_BLAKE3, "base_object_blake3":arm["base_object_blake3"],
        "pool_object_blake3":arm["pool_object_blake3"],
        "target":arm["target"], "identity_mask":0, "fixed_x_blocks":fixed_x_blocks,
        "domain_clauses":domain_size.clauses,
        "variables":enc.solver.n_vars(), "clauses":enc.solver.n_clauses(),
        "xors":enc.solver.n_xors(), "conflict_cap":conflict_cap,
        "conflicts":enc.solver.conflicts(), "model_verified":model_verified,
        "elapsed_wall_ms_L0":started.elapsed().as_secs_f64()*1000.0,
        "natural_relations":Value::Null, "total_index_calculus_runtime_ms":Value::Null,
        "selected_best_total_runtime":Value::Null,
    }))
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 {
        return Err("usage: koblitz_n83_m6_witness_ladder ROLE STAGE CONFLICT_CAP".into());
    }
    println!(
        "{}",
        serde_json::to_string(&run(&args[1], &args[2], args[3].parse()?)?)?
    );
    Ok(())
}

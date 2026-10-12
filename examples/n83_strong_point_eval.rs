//! Exact public-point strong rho control for the legacy n=83, a=1, K600 IC lead.
//!
//! Usage: n83_strong_point_eval <one-point.jsonl> <walk-seed> [receipt.json]
//! The input must contain exactly one JSON array of two decimal coordinate
//! strings. No known-answer scalar or validation sidecar is accepted or opened.
//!
//! Curve/normal-basis/G-only-jump preparation is target independent. Online time
//! includes target membership/subgroup checks, Q-dependent walk starts, the
//! distinguished-point solve and a separate generic arithmetic scalar replay.
//! All preparation, validation and output boundaries are reported explicitly.
//! This worker does not certify hardware isolation or infer a speedup by itself.
#![recursion_limit = "256"]

use std::fs;
use std::io::Write;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use crypto_lib::cryptanalysis::koblitz_strong_rho::{
    RawPointG, StrongRhoCharges, StrongRhoG, StrongRhoParams,
};
use crypto_lib::cryptanalysis::wide_gf2m::WideGf2;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use rand::{rngs::StdRng, SeedableRng};
use serde_json::{json, Value};

type ExactStrong = StrongRhoG<WideGf2, u64>;
const BASE_SNAPSHOT: &str = "a8194d616ed8c04489ec46177d369f1b1670360a";
const LEGACY_R: u64 = 8_569_786_107_849_059;
const LEGACY_G: [u128; 2] = [
    6_614_244_651_398_184_860_292_015,
    8_868_260_410_029_249_391_539_529,
];

fn parse_public_point(text: &str, n: u32) -> Result<[u128; 2], String> {
    let mut rows = text.lines().filter(|line| !line.trim().is_empty());
    let row = rows.next().ok_or("public point input is empty")?;
    if rows.next().is_some() {
        return Err("exactly one public point is required".into());
    }
    let strings: [String; 2] = serde_json::from_str(row)
        .map_err(|_| "input must be a JSON array of two decimal strings".to_owned())?;
    let x = strings[0].parse::<u128>().map_err(|_| "bad x coordinate")?;
    let y = strings[1].parse::<u128>().map_err(|_| "bad y coordinate")?;
    let bound = 1u128 << n;
    if x >= bound || y >= bound {
        return Err("coordinate is not a canonical field word".into());
    }
    Ok([x, y])
}

fn point_strings(point: &BinaryPoint) -> Result<[String; 2], String> {
    match point {
        BinaryPoint::Infinity => Err("affine point required".into()),
        BinaryPoint::Affine { x, y } => {
            Ok([x.to_biguint().to_string(), y.to_biguint().to_string()])
        }
    }
}

fn make_strong(curve: &KoblitzCurve) -> Result<ExactStrong, String> {
    let r = curve
        .subgroup_order
        .to_u64()
        .ok_or("subgroup scalar does not fit u64")?;
    let lambda = curve.lambda.to_u64().ok_or("lambda does not fit u64")?;
    Ok(ExactStrong::from_parts(
        WideGf2::new(&curve.curve.irreducible),
        u128::from(curve.a),
        r,
        lambda,
        RawPointG::<u128>::from_binary_wide(curve.generator()),
    ))
}

fn assert_legacy_curve(curve: &KoblitzCurve) -> Result<(), String> {
    if curve.n != 83
        || curve.a != 1
        || curve.curve.irreducible.low_terms != [0, 2, 4, 7]
        || curve.subgroup_order.to_u64() != Some(LEGACY_R)
        || curve.cofactor.to_u64() != Some(1_128_547_018)
        || point_strings(curve.generator())? != [LEGACY_G[0].to_string(), LEGACY_G[1].to_string()]
    {
        return Err("constructor does not match the exact legacy K600 curve and generator".into());
    }
    Ok(())
}

fn charges_json(c: StrongRhoCharges) -> Value {
    json!({
        "group_additions": c.group_additions,
        "scalar_multiplications": c.scalar_multiplications,
        "canonicalizations": c.canonicalizations,
        "partition_hashes": c.partition_hashes,
        "table_queries": c.table_queries,
        "table_inserts": c.table_inserts,
        "failed_collisions": c.failed_collisions,
    })
}

fn peak_rss_bytes() -> Option<u64> {
    // SAFETY: getrusage fills the initialized plain-data structure for this process.
    let mut usage: libc::rusage = unsafe { std::mem::zeroed() };
    if unsafe { libc::getrusage(libc::RUSAGE_SELF, &mut usage) } != 0 {
        return None;
    }
    let raw = u64::try_from(usage.ru_maxrss).ok()?;
    #[cfg(target_os = "macos")]
    {
        Some(raw)
    }
    #[cfg(not(target_os = "macos"))]
    {
        Some(raw.saturating_mul(1024))
    }
}

fn evaluate(input: &str, seed: u64) -> Result<(Value, bool), String> {
    let process_started = Instant::now();
    let started_unix_ms = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|e| e.to_string())?
        .as_millis();
    let input_bytes = fs::read(input).map_err(|e| format!("read public point: {e}"))?;
    let input_text = std::str::from_utf8(&input_bytes).map_err(|_| "input is not UTF-8")?;
    let coordinates = parse_public_point(input_text, 83)?;
    let input_blake3 = blake3::hash(&input_bytes).to_hex().to_string();

    let prepare_started = Instant::now();
    // This constructor is deliberately different from known_n83_k1(): that
    // helper uses modulus [0,1,2,45] and a different generator.
    let curve = KoblitzCurve::new(1, 83).ok_or("legacy curve construction failed")?;
    assert_legacy_curve(&curve)?;
    let exact_id = curve.curve_id();
    let gf = WideGf2::new(&curve.curve.irreducible);
    let kernel = gf.kernel_name();
    let rho = make_strong(&curve)?;
    let mut setup_charges = StrongRhoCharges::default();
    // G-only jumps can be prepared before Q arrives. Walk-start RNG and jump
    // RNG domains are separate, both derived only from the declared seed.
    let mut jump_material = b"n83-exact-public-strong-jumps-v1\0".to_vec();
    jump_material.extend_from_slice(&seed.to_le_bytes());
    let digest = blake3::hash(&jump_material);
    let jump_seed = u64::from_le_bytes(digest.as_bytes()[..8].try_into().unwrap());
    let jumps = rho.jumps(jump_seed, &mut setup_charges);
    let prepare_ns = prepare_started.elapsed().as_nanos();
    let params = StrongRhoParams::default();
    let mut rng = StdRng::seed_from_u64(seed);
    let predicted_step_cap = (rho.ideal_steps().ceil() as u64)
        .saturating_mul(params.step_cap_factor)
        .max(1_000_000);

    let online_started = Instant::now();
    let query_started = Instant::now();
    let q = BinaryPoint::Affine {
        x: F2mElement::from_biguint(&BigUint::from(coordinates[0]), 83),
        y: F2mElement::from_biguint(&BigUint::from(coordinates[1]), 83),
    };
    let on_curve = curve.curve.is_on_curve(&q);
    let in_subgroup = on_curve && curve.mul(&q, &curve.subgroup_order) == BinaryPoint::Infinity;
    let raw_q = RawPointG::<u128>::from_binary_wide(&q);
    let query_ns = query_started.elapsed().as_nanos();

    let walk_started = Instant::now();
    let outcome = if in_subgroup {
        rho.solve(raw_q, &jumps, &mut rng, &params, setup_charges)
    } else {
        None
    };
    let walk_ns = walk_started.elapsed().as_nanos();
    let validation_started = Instant::now();
    let recovered = outcome.as_ref().map(|value| value.scalar);
    let replay_verified =
        recovered.is_some_and(|d| curve.mul(curve.generator(), &BigUint::from(d)) == q);
    let validation_ns = validation_started.elapsed().as_nanos();
    let online_ns = online_started.elapsed().as_nanos();
    let phase_sum_ns = query_ns + walk_ns + validation_ns;
    let outcome_json = outcome.as_ref().map(|out| {
        json!({
            "walk_steps": out.walk_steps,
            "walks": out.walks,
            "fruitless_cycles": out.fruitless,
            "fruitless_escapes": out.escapes,
            "capped_walks": out.capped,
            "wasted_merges": out.wasted_merges,
            "table_entries": out.table_entries,
            "table_payload_lower_bound_bytes": out.table_entries * (1 + 2 * 16 + 2 * 8),
            "ideal_steps": out.ideal_steps,
            "automorphism_size": out.automorphisms,
            "charges_setup_and_solve": charges_json(out.charges),
        })
    });
    let status = if !in_subgroup {
        "REJECTED_PUBLIC_POINT"
    } else if outcome.is_none() {
        "UNKNOWN_STEP_CAP"
    } else if !replay_verified {
        "FAILED_SCALAR_REPLAY"
    } else {
        "PASS_VERIFIED_PUBLIC_TARGET"
    };
    let executable = std::env::current_exe().ok().and_then(|p| fs::read(p).ok());
    let source = include_bytes!("n83_strong_point_eval.rs");
    let receipt = json!({
        "schema": "n83-exact-public-strong-evaluation/v1",
        "status": status,
        "run_id": format!("N83A1StrongW{seed}Q{}", &input_blake3[..12]),
        "canonical_candidate_id": null,
        "canonical_workload_id": null,
        "started_unix_ms": started_unix_ms,
        "record_class": "verified_answer_exploratory_wall",
        "host_isolation_receipt": null,
        "controlled_speedup": null,
        "worker_count": 1,
        "target_count": 1,
        "known_answer_input": false,
        "validation_sidecar_opened": false,
        "base_snapshot": BASE_SNAPSHOT,
        "compiled_source_revision": option_env!("N83_SOURCE_COMMIT"),
        "worker_source_blake3": blake3::hash(source).to_hex().to_string(),
        "executable_blake3": executable.as_deref().map(|b| blake3::hash(b).to_hex().to_string()),
        "public_input_blake3": input_blake3,
        "reference_grade": "strong_wide_field_u64_scalar",
        "curve_id": exact_id.icv1,
        "curve_slug": exact_id.slug,
        "n": 83,
        "a": 1,
        "b": 1,
        "field_modulus_low_terms": [0, 2, 4, 7],
        "field_encoding": "polynomial words as canonical decimal strings",
        "subgroup_order": LEGACY_R.to_string(),
        "cofactor": curve.cofactor.to_string(),
        "lambda": rho.lambda().to_string(),
        "generator": point_strings(curve.generator())?,
        "public_target_q": [coordinates[0].to_string(), coordinates[1].to_string()],
        "published_q": [coordinates[0].to_string(), coordinates[1].to_string()],
        "target_on_curve": on_curve,
        "target_in_subgroup": in_subgroup,
        "recovered_scalar": recovered.map(|d| d.to_string()),
        "generic_scalar_replay_verified": replay_verified,
        "group_verified": replay_verified,
        "reference_group_validation": replay_verified,
        "arithmetic_kernel": kernel,
        "walk_seed": seed,
        "jump_seed": jump_seed,
        "rho_policy": {
            "quotient": "signed_frobenius",
            "lanes": params.lanes,
            "distinguished_bits": params.dp_bits,
            "step_cap_factor": params.step_cap_factor,
            "step_cap": predicted_step_cap,
            "escape_fruitless": params.escape_fruitless,
            "jump_count": jumps.len(),
            "G_only_jumps_target_independent": true,
        },
        "preparation_ns_excluded": prepare_ns,
        "preparation_charges": charges_json(setup_charges),
        "online_interval": {
            "start_event": "before_public_target_membership_and_subgroup_validation",
            "stop_event": "after_generic_scalar_replay",
            "included_stages": ["target_validation", "walk_starts", "walk_and_collision", "generic_scalar_replay"],
            "excluded_stages": ["input_loading", "curve_normal_basis_G_only_jumps", "receipt_serialization", "receipt_io"],
        },
        "online_phase_ns": {
            "target_query_and_validation": query_ns,
            "walk_and_collision": walk_ns,
            "generic_scalar_replay": validation_ns,
        },
        "online_phase_sum_ns": phase_sum_ns,
        "online_interval_ns": online_ns,
        "online_timer_boundary_residual_ns": online_ns.saturating_sub(phase_sum_ns),
        "online_ms": online_ns as f64 / 1_000_000.0,
        "online_total_ms": online_ns as f64 / 1_000_000.0,
        "target_ms": online_ns as f64 / 1_000_000.0,
        "outcome": outcome_json,
        "failure_counter_note": if outcome.is_none() {
            Some("library solve returns no outcome on cap; actual steps/charges are unavailable, not zero")
        } else { None },
        "peak_rss_bytes": peak_rss_bytes(),
        "peak_rss_source": "getrusage RUSAGE_SELF ru_maxrss; platform-scaled bytes",
        "process_ns_before_receipt_io": process_started.elapsed().as_nanos(),
        "claim_boundary": "public synthetic exact legacy n83 a1 target; answer correctness and exploratory local timing only; no admitted speedup",
    });
    Ok((receipt, replay_verified))
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if !(3..=4).contains(&args.len()) {
        eprintln!("usage: n83_strong_point_eval <one-point.jsonl> <walk-seed> [receipt.json]");
        std::process::exit(2);
    }
    let result = args[2]
        .parse::<u64>()
        .map_err(|_| "walk seed must be u64".to_owned())
        .and_then(|seed| evaluate(&args[1], seed));
    let (receipt, success) = match result {
        Ok(value) => value,
        Err(message) => (
            json!({"schema":"n83-exact-public-strong-evaluation/v1", "status":"INPUT_OR_SETUP_FAILURE", "error":message}),
            false,
        ),
    };
    let bytes = format!("{}\n", serde_json::to_string(&receipt).unwrap());
    if let Some(path) = args.get(3) {
        let write_result = fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(path)
            .and_then(|mut file| file.write_all(bytes.as_bytes()));
        if let Err(error) = write_result {
            eprintln!("receipt write failed: {error}");
            std::process::exit(2);
        }
    }
    print!("{bytes}");
    if !success {
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crypto_lib::cryptanalysis::koblitz_strong_rho::StrongRho;

    #[test]
    fn point_input_refuses_scalar_sidecars_noncanonical_and_multiple_targets() {
        assert!(parse_public_point("3141592653589793", 83).is_err());
        assert!(parse_public_point(r#"{"fixture_scalar":42}"#, 83).is_err());
        assert!(parse_public_point("[1,2]", 83).is_err());
        assert!(parse_public_point("[\"1\",\"2\"]\n[\"3\",\"4\"]", 83).is_err());
        let high = format!("[\"{}\",\"2\"]", 1u128 << 83);
        assert!(parse_public_point(&high, 83).is_err());
        assert_eq!(parse_public_point("[\"1\",\"2\"]\n", 83).unwrap(), [1, 2]);
    }

    #[test]
    fn wide_field_control_reproduces_narrow_reference_on_small_curves() {
        for n in [7, 17] {
            let curve = KoblitzCurve::new(1, n).expect("small control curve");
            let narrow = StrongRho::new(&curve);
            let wide = make_strong(&curve).unwrap();
            assert_eq!(wide.modulus(), narrow.modulus());
            assert_eq!(wide.lambda(), narrow.lambda());
            assert_eq!(wide.automorphisms(), narrow.automorphisms());
            for seed in [1, 29] {
                let scalar = 37 % narrow.modulus();
                let q_narrow = narrow.scalar_mul(narrow.generator(), scalar);
                let q_wide = wide.scalar_mul(wide.generator(), scalar);
                let mut charges_narrow = StrongRhoCharges::default();
                let mut charges_wide = StrongRhoCharges::default();
                let jumps_narrow = narrow.jumps(seed ^ 0x5eed, &mut charges_narrow);
                let jumps_wide = wide.jumps(seed ^ 0x5eed, &mut charges_wide);
                // The n=7 signed quotient has only about five nonzero
                // classes; dp_bits=4 may select none of them. Store every
                // point in this tiny equivalence control. Production n=83
                // parameters remain the unchanged strong defaults.
                let params = StrongRhoParams {
                    lanes: 4,
                    dp_bits: 0,
                    ..StrongRhoParams::default()
                };
                let narrow_result = narrow
                    .solve(
                        q_narrow,
                        &jumps_narrow,
                        &mut StdRng::seed_from_u64(seed),
                        &params,
                        charges_narrow,
                    )
                    .expect("narrow known-answer control");
                let wide_result = wide
                    .solve(
                        q_wide,
                        &jumps_wide,
                        &mut StdRng::seed_from_u64(seed),
                        &params,
                        charges_wide,
                    )
                    .expect("wide known-answer control");
                assert_eq!(wide_result.scalar, scalar);
                assert_eq!(narrow_result.scalar, wide_result.scalar);
                assert_eq!(narrow_result.walk_steps, wide_result.walk_steps);
                assert_eq!(narrow_result.walks, wide_result.walks);
                assert_eq!(narrow_result.fruitless, wide_result.fruitless);
                assert_eq!(narrow_result.capped, wide_result.capped);
                assert_eq!(narrow_result.wasted_merges, wide_result.wasted_merges);
                assert_eq!(narrow_result.table_entries, wide_result.table_entries);
                assert_eq!(narrow_result.charges, wide_result.charges);
                let q_general = curve.mul(curve.generator(), &BigUint::from(scalar));
                assert_eq!(
                    curve.mul(curve.generator(), &BigUint::from(wide_result.scalar)),
                    q_general,
                );
            }
        }
    }
}

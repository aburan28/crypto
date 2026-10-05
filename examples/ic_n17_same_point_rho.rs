//! One disclosed public-point strong-rho diagnostic paired to the F5 target.
//! This is not a frozen worker or a controlled CPU speedup measurement.
use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use crypto_lib::cryptanalysis::koblitz_strong_rho::{
    RawPoint, StrongRho, StrongRhoCharges, StrongRhoParams,
};
use num_bigint::BigUint;
use rand::{rngs::StdRng, SeedableRng};
use serde_json::json;
use std::time::Instant;

const TARGET: [u64; 2] = [61889, 74818];
const GENERATOR: [u64; 2] = [43693, 23339];
const ORDER: u64 = 65587;
const JUMP_SEED: u64 = 2026100502;
const WALK_SEED: u64 = 2026100503;

fn ns(start: Instant, end: Instant) -> u64 {
    end.duration_since(start)
        .as_nanos()
        .try_into()
        .expect("clock overflow")
}
fn encoded(point: &BinaryPoint) -> Option<[u64; 2]> {
    match point {
        BinaryPoint::Infinity => None,
        BinaryPoint::Affine { x, y } => Some([
            x.raw_bits().first().copied().unwrap_or(0),
            y.raw_bits().first().copied().unwrap_or(0),
        ]),
    }
}
fn charges(value: StrongRhoCharges) -> serde_json::Value {
    json!({"group_additions":value.group_additions,
        "scalar_multiplications":value.scalar_multiplications,
        "canonicalizations":value.canonicalizations,
        "partition_hashes":value.partition_hashes,
        "table_queries":value.table_queries,
        "table_inserts":value.table_inserts,
        "failed_collisions":value.failed_collisions})
}

fn main() {
    let setup_started = Instant::now();
    let curve = KoblitzCurve::new(1, 17).expect("exact n17 curve");
    assert_eq!(curve.subgroup_order, BigUint::from(ORDER));
    assert_eq!(encoded(curve.generator()), Some(GENERATOR));
    let rho = StrongRho::new(&curve);
    assert_eq!(rho.modulus(), ORDER);
    let params = StrongRhoParams::default();
    assert_eq!(
        (params.lanes, params.dp_bits, params.step_cap_factor),
        (32, 4, 2000)
    );
    let mut prep_charges = StrongRhoCharges::default();
    let jumps = rho.jumps(JUMP_SEED, &mut prep_charges);
    let setup_end = Instant::now();

    // Q is supplied as public coordinates. Validation and conversion are
    // target dependent and charged before the first walk operation.
    let online_start = Instant::now();
    let target = BinaryPoint::Affine {
        x: F2mElement::from_biguint(&BigUint::from(TARGET[0]), 17),
        y: F2mElement::from_biguint(&BigUint::from(TARGET[1]), 17),
    };
    assert!(curve.curve.is_on_curve(&target));
    assert_eq!(
        curve.mul(&target, &BigUint::from(ORDER)),
        BinaryPoint::Infinity
    );
    let q = RawPoint::from_binary(&target);
    let walk_start = Instant::now();
    let mut rng = StdRng::seed_from_u64(WALK_SEED);
    let outcome = rho.solve(q, &jumps, &mut rng, &params, prep_charges);
    let walk_end = Instant::now();

    let common = json!({"schema_version":1,"question":"n17-disclosed-same-point-strong-rho-v1",
        "target":TARGET,"generator":GENERATOR,"subgroup_order":ORDER,
        "target_count":1,"fixture_scalar_generated":false,"worker_count":1,
        "rho_policy":"signed-frobenius-distinguished-point-normal-basis-batch-inversion",
        "lanes":params.lanes,"dp_bits":params.dp_bits,"step_cap_factor":params.step_cap_factor,
        "jump_seed":JUMP_SEED,"walk_seed":WALK_SEED,
        "setup_ns":ns(setup_started,setup_end),
        "target_validation_and_conversion_ns":ns(online_start,walk_start),
        "walk_ns":ns(walk_start,walk_end),
        "jump_setup_charges":charges(prep_charges),
        "peak_memory_bytes":null,"candidate_id":null,"workload_id":null,"run_id":null,
        "source_bound_execution_admitted":false,"fresh_paired_qualification":false,
        "headline_eligible":false,"promotion_eligible":false,
        "online_speedup":null,"wall_timing_class":"ordinary-host-exploratory; isolation unverified"});
    match outcome {
        Some(found) => {
            // Independent polynomial-basis group arithmetic, charged inside
            // the one-target online interval after normal-basis rho stops.
            assert_eq!(
                curve.mul(curve.generator(), &BigUint::from(found.scalar)),
                target
            );
            let online_end = Instant::now();
            assert_eq!(
                ns(online_start, walk_start) + ns(walk_start, walk_end) + ns(walk_end, online_end),
                ns(online_start, online_end)
            );
            let mut row = common;
            row["status"] = json!("VERIFIED_SAME_POINT_STRONG_RHO_DIAGNOSTIC");
            row["verified_recovery"] = json!(true);
            row["scalar"] = json!(found.scalar);
            row["recovery_check_ns"] = json!(ns(walk_end, online_end));
            row["online_ns"] = json!(ns(online_start, online_end));
            row["walk_steps"] = json!(found.walk_steps);
            row["walks"] = json!(found.walks);
            row["distinguished_point_entries_at_stop"] = json!(found.table_entries);
            row["automorphisms"] = json!(found.automorphisms);
            row["whole_solve_charges_including_jump_setup"] = charges(found.charges);
            row["recovery_check_scalar_multiplications"] = json!(1);
            println!("{row}");
        }
        None => {
            let mut row = common;
            row["status"] = json!("STEP_CAP_EXHAUSTED");
            row["verified_recovery"] = json!(false);
            row["scalar"] = serde_json::Value::Null;
            row["recovery_check_ns"] = serde_json::Value::Null;
            row["online_ns"] = json!(ns(online_start, walk_end));
            println!("{row}");
            std::process::exit(1);
        }
    }
}

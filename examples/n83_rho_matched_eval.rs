//! Matched one-target rho workloads with genuinely reusable preprocessing.
//!
//! Usage: n83_rho_matched_eval POINTS.jsonl SEEDS.json CONFIG.json TARGETS.jsonl
//! Every solve has a fresh RNG, walk state and distinguished-point table. Only
//! curve/normal-basis data and G-only jump tables are prepared before Q arrives.
//! Point files contain decimal coordinate strings; seeds never contain scalars.
use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use crypto_lib::cryptanalysis::koblitz_strong_rho::{
    RawPointG, RhoField, StrongRhoCharges, StrongRhoG, StrongRhoParams,
};
use crypto_lib::cryptanalysis::semaev_decomp::Gf2_128;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use rand::{rngs::StdRng, SeedableRng};
use serde::Deserialize;
use serde_json::{json, Value};
use std::fs;
use std::io::Write;
use std::time::Instant;

const LEGACY_R: u64 = 8_569_786_107_849_059;
const LEGACY_G: [u128; 2] = [
    6_614_244_651_398_184_860_292_015,
    8_868_260_410_029_249_391_539_529,
];

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Config {
    backend: String,
    scalar_kernel: String,
    lanes: usize,
    dp_bits: u32,
    escape_fruitless: bool,
    step_cap_factor: u64,
}
impl Config {
    fn parse(bytes: &[u8]) -> Result<Self, String> {
        let c: Self = serde_json::from_slice(bytes).map_err(|e| e.to_string())?;
        if c.backend != "fold"
            || c.scalar_kernel != "reference"
            || c.lanes != 256
            || c.dp_bits != 6
            || c.escape_fruitless
            || c.step_cap_factor != 2000
        {
            return Err(
                "configuration must be the frozen Fold/reference/256/DP6/escape-off/2000 policy"
                    .into(),
            );
        }
        Ok(c)
    }
    fn params(&self) -> StrongRhoParams {
        StrongRhoParams {
            lanes: self.lanes,
            dp_bits: self.dp_bits,
            escape_fruitless: self.escape_fruitless,
            step_cap_factor: self.step_cap_factor,
        }
    }
    fn json(&self) -> Value {
        json!({"backend":self.backend,"scalar_kernel":self.scalar_kernel,"lanes":self.lanes,"dp_bits":self.dp_bits,"escape_fruitless":self.escape_fruitless,"step_cap_factor":self.step_cap_factor})
    }
}

struct FoldField(Gf2_128);
impl RhoField for FoldField {
    type E = u128;
    fn degree(&self) -> u32 {
        self.0.n
    }
    #[inline(always)]
    fn mul(&self, a: u128, b: u128) -> u128 {
        self.0.mul(a, b)
    }
    #[inline(always)]
    fn sqr(&self, a: u128) -> u128 {
        self.0.sqr(a)
    }
    fn batch_inv(&self, xs: &mut [u128], scratch: &mut Vec<u128>) {
        self.0.batch_inv(xs, scratch);
    }
}
fn sha256(bytes: &[u8]) -> String {
    hex::encode(crypto_lib::hash::sha256(bytes))
}
fn charges(c: StrongRhoCharges) -> Value {
    json!({"group_additions":c.group_additions,"scalar_multiplications":c.scalar_multiplications,"canonicalizations":c.canonicalizations,"partition_hashes":c.partition_hashes,"table_queries":c.table_queries,"table_inserts":c.table_inserts,"failed_collisions":c.failed_collisions})
}
fn parse_points(text: &str) -> Result<Vec<[u128; 2]>, String> {
    let points = text
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|line| {
            let s: [String; 2] = serde_json::from_str(line)
                .map_err(|_| "point needs two decimal coordinate strings")?;
            let x = s[0].parse::<u128>().map_err(|_| "invalid x")?;
            let y = s[1].parse::<u128>().map_err(|_| "invalid y")?;
            if x >= 1u128 << 83 || y >= 1u128 << 83 {
                return Err("noncanonical field coordinate".to_string());
            }
            Ok([x, y])
        })
        .collect::<Result<Vec<_>, String>>()?;
    if points.len() != 1 && points.len() != 32 {
        return Err("a singleton or exactly 32 independent public points is required".into());
    }
    let unique: std::collections::BTreeSet<_> = points.iter().collect();
    if unique.len() != points.len() {
        return Err("duplicate public target".into());
    }
    Ok(points)
}
fn parse_seeds(bytes: &[u8], count: usize) -> Result<Vec<u64>, String> {
    let seeds: Vec<u64> = serde_json::from_slice(bytes).map_err(|e| format!("seed list: {e}"))?;
    if seeds.len() != count {
        return Err("one u64 walk seed per public point is required".into());
    }
    Ok(seeds)
}
fn verify_curve(curve: &KoblitzCurve) -> Result<(), String> {
    let g = RawPointG::<u128>::from_binary_wide(curve.generator());
    if curve.n != 83
        || curve.a != 1
        || curve.curve.irreducible.low_terms != [0, 2, 4, 7]
        || curve.subgroup_order.to_u64() != Some(LEGACY_R)
        || curve.cofactor.to_u64() != Some(1_128_547_018)
        || curve.lambda.to_u64() != Some(2_998_004_823_950_457)
        || curve.curve.b.to_biguint() != BigUint::from(1u8)
        || g != (RawPointG::Affine {
            x: LEGACY_G[0],
            y: LEGACY_G[1],
        })
    {
        return Err("curve/generator differs from the exact legacy fixture".into());
    }
    Ok(())
}
fn run(args: &[String]) -> Result<bool, String> {
    if args.len() != 4 {
        return Err(
            "usage: n83_rho_matched_eval POINTS.jsonl SEEDS.json CONFIG.json TARGETS.jsonl".into(),
        );
    }
    let process_started = Instant::now();
    let point_bytes = fs::read(&args[0]).map_err(|e| e.to_string())?;
    let points = parse_points(std::str::from_utf8(&point_bytes).map_err(|e| e.to_string())?)?;
    let seed_bytes = fs::read(&args[1]).map_err(|e| e.to_string())?;
    let seeds = parse_seeds(&seed_bytes, points.len())?;
    let config_bytes = fs::read(&args[2]).map_err(|e| e.to_string())?;
    let config = Config::parse(&config_bytes)?;
    let mut out = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&args[3])
        .map_err(|e| e.to_string())?;
    let preparation_started = Instant::now();
    let curve = KoblitzCurve::new(1, 83).ok_or("legacy curve construction failed")?;
    verify_curve(&curve)?;
    let rho = StrongRhoG::<FoldField, u64>::from_parts(
        FoldField(Gf2_128::new(&curve.curve.irreducible)),
        1,
        LEGACY_R,
        curve.lambda.to_u64().ok_or("lambda does not fit u64")?,
        RawPointG::<u128>::from_binary_wide(curve.generator()),
    );
    let shared_curve_normal_basis_ns = preparation_started.elapsed().as_nanos();
    let jump_started = Instant::now();
    let jump_tables: Vec<_> = seeds
        .iter()
        .map(|seed| {
            let mut material = b"n83-exact-public-strong-jumps-v1\0".to_vec();
            material.extend_from_slice(&seed.to_le_bytes());
            let hash = blake3::hash(&material);
            let jump_seed = u64::from_le_bytes(hash.as_bytes()[..8].try_into().unwrap());
            let mut c = StrongRhoCharges::default();
            let table = rho.jumps(jump_seed, &mut c);
            (table, c, jump_seed)
        })
        .collect();
    let all_g_only_jumps_ns = jump_started.elapsed().as_nanos();
    let preparation_ns = preparation_started.elapsed().as_nanos();
    let source = include_bytes!("n83_rho_matched_eval.rs");
    let provenance = json!({
        "compiled_source_revision":option_env!("N83_SOURCE_COMMIT"),
        "worker_source_sha256":sha256(source),
        "compiled_source_sha256":{
            "examples/n83_rho_matched_eval.rs":sha256(source),
            "src/cryptanalysis/koblitz_strong_rho.rs":sha256(include_bytes!("../src/cryptanalysis/koblitz_strong_rho.rs")),
            "src/cryptanalysis/semaev_decomp.rs":sha256(include_bytes!("../src/cryptanalysis/semaev_decomp.rs")),
        },
        "public_input_sha256":sha256(&point_bytes),"seed_input_sha256":sha256(&seed_bytes),
        "config_input_sha256":sha256(&config_bytes),
        "executable_sha256":std::env::current_exe().ok().and_then(|p|fs::read(p).ok()).map(|b|sha256(&b)),
    });
    let params = config.params();
    let targets_started = Instant::now();
    let mut solved = 0usize;
    let mut online_sum_ns = 0u128;
    for (index, coordinates) in points.iter().enumerate() {
        let (jumps, jump_charges, jump_seed) = &jump_tables[index];
        let mut rng = StdRng::seed_from_u64(seeds[index]);
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
        // solve owns its DP table, lane and trail vectors locally. No state
        // carrying Q or collisions survives this call into another target.
        let outcome = if in_subgroup {
            rho.solve(raw_q, jumps, &mut rng, &params, *jump_charges)
        } else {
            None
        };
        let walk_ns = walk_started.elapsed().as_nanos();
        let replay_started = Instant::now();
        let recovered = outcome.as_ref().map(|o| o.scalar);
        let verified = recovered.is_some_and(|d| {
            d > 0 && d < LEGACY_R && curve.mul(curve.generator(), &BigUint::from(d)) == q
        });
        let replay_ns = replay_started.elapsed().as_nanos();
        let online_ns = online_started.elapsed().as_nanos();
        online_sum_ns += online_ns;
        if verified {
            solved += 1;
        }
        let phase_sum_ns = query_ns + walk_ns + replay_ns;
        let outcome_json = outcome.as_ref().map(|o| json!({
            "walk_steps":o.walk_steps,"walks":o.walks,"fruitless_cycles":o.fruitless,
            "fruitless_escapes":o.escapes,"capped_walks":o.capped,"wasted_merges":o.wasted_merges,
            "table_entries":o.table_entries,"ideal_steps":o.ideal_steps,"automorphism_size":o.automorphisms,
            "charges_jump_preparation_and_solve":charges(o.charges),
        }));
        let record = json!({
            "schema":"n83-rho-matched-target/v1","kind":"n83_rho_matched_target",
            "status":if verified {"PASS_VERIFIED_PUBLIC_TARGET"} else if !in_subgroup {"REJECTED_PUBLIC_POINT"} else if outcome.is_none() {"UNKNOWN_STEP_CAP"} else {"FAILED_SCALAR_REPLAY"},
            "fixture_index":index,"target_count":1,"cohort_target_count":points.len(),
            "published_q":[coordinates[0].to_string(),coordinates[1].to_string()],
            "target":[coordinates[0].to_string(),coordinates[1].to_string()],
            "curve_slug":curve.curve_id().slug,"field_modulus_low_terms":[0,2,4,7],
            "subgroup_order":LEGACY_R.to_string(),"generator":[LEGACY_G[0].to_string(),LEGACY_G[1].to_string()],
            "walk_seed":seeds[index],"jump_seed":jump_seed,"optimization_config":config.json(),
            "known_answer_input":false,"validation_sidecar_opened":false,
            "group_verified":verified,"verified":verified,"exit_code":if verified {0}else{1},
            "target_on_curve":on_curve,"target_in_subgroup":in_subgroup,
            "recovered_scalar":recovered.map(|d|d.to_string()),
            "online_total_ms":online_ns as f64 / 1e6,"target_ms":online_ns as f64 / 1e6,
            "online_interval_ns":online_ns,"online_phase_sum_ns":phase_sum_ns,
            "online_timer_boundary_residual_ns":online_ns.saturating_sub(phase_sum_ns),
            "online_phase_ns":{"target_query_and_validation":query_ns,"walk_and_collision":walk_ns,"generic_scalar_replay":replay_ns},
            "online_interval":{"start_event":"before_public_target_membership_and_subgroup_validation","stop_event":"after_generic_scalar_replay","excluded":["public_input_loading","shared_curve_normal_basis","G_only_jumps","receipt_hashing_serialization_and_IO"]},
            "fresh_target_dependent_DP_walk_and_rng_state":true,"cross_target_collision_sharing":false,
            "preparation_charges_for_this_seed_only":charges(*jump_charges),
            "outcome":outcome_json,"provenance":provenance,
        });
        writeln!(out, "{record}").map_err(|e| e.to_string())?;
        out.flush().map_err(|e| e.to_string())?;
    }
    let summary = json!({
        "kind":"n83_rho_matched_summary","schema":"n83-rho-matched-summary/v1",
        "targets":points.len(),"targets_solved":solved,"targets_failed":points.len()-solved,
        "all_verified":solved==points.len(),"optimization_config":config.json(),
        "curve_slug":curve.curve_id().slug,"target_dependent_state_shared":false,
        "target_independent_curve_normal_basis_constructions":1,"G_only_jump_tables":jump_tables.len(),
        "timing_ms":{"setup_total":preparation_ns as f64/1e6,"shared_curve_normal_basis":shared_curve_normal_basis_ns as f64/1e6,"all_G_only_jumps":all_g_only_jumps_ns as f64/1e6,"targets_online_sum":online_sum_ns as f64/1e6,"targets_total":targets_started.elapsed().as_secs_f64()*1000.,"process_total":process_started.elapsed().as_secs_f64()*1000.},
        "producer_preparation_ms":preparation_ns as f64/1e6,"setup_ms":preparation_ns as f64/1e6,
        "provenance":provenance,
    });
    println!("{summary}");
    out.sync_all().map_err(|e| e.to_string())?;
    Ok(solved == points.len())
}
fn main() {
    let result = run(&std::env::args().skip(1).collect::<Vec<_>>());
    match result {
        Ok(true) => {}
        Ok(false) => std::process::exit(1),
        Err(e) => {
            eprintln!("n83 matched rho: {e}");
            std::process::exit(2);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn config() -> Value {
        json!({"backend":"fold","scalar_kernel":"reference","lanes":256,"dp_bits":6,"escape_fruitless":false,"step_cap_factor":2000})
    }
    #[test]
    fn policy_drift_and_unknown_fields_rejected() {
        let mut c = config();
        c["lanes"] = json!(64);
        assert!(Config::parse(c.to_string().as_bytes()).is_err());
        let mut c = config();
        c["known_scalar"] = json!(1);
        assert!(Config::parse(c.to_string().as_bytes()).is_err());
    }
    #[test]
    fn numeric_coordinates_and_scalars_rejected() {
        assert!(parse_points("[1,2]").is_err());
        assert!(parse_points("123456").is_err());
        assert!(parse_points("[\"9671406556917033397649408\",\"0\"]").is_err());
    }
    #[test]
    fn seed_point_count_must_match() {
        assert!(parse_seeds(b"[1,2]", 1).is_err());
        assert!(parse_seeds(b"[1]", 1).is_ok());
    }
    #[test]
    fn exact_legacy_constructor_is_pinned() {
        verify_curve(&KoblitzCurve::new(1, 83).unwrap()).unwrap();
    }
}

//! Configured scalar-kernel strong rho control for the exact legacy n=83/a=1 curve.
//!
//! Usage: n83_rho_speed_eval <point.jsonl> <walk-seed> <config.json> [receipt.json]
//! The input must contain exactly one JSON array of two decimal coordinate
//! strings. No known-answer scalar or validation sidecar is accepted or opened.
//!
//! Curve/normal-basis/G-only-jump preparation is target independent. Online time
//! includes target membership/subgroup checks, Q-dependent walk starts, the
//! distinguished-point solve and a separate generic arithmetic scalar replay.
//! A strict config selects reference or exact Barrett scalars and a bounded policy.
//! The production StrongRhoG algorithm is unchanged. All boundaries are explicit.
//! This worker does not certify hardware isolation or infer a speedup by itself.
#![recursion_limit = "256"]

use std::fs;
use std::io::Write;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use crypto_lib::cryptanalysis::koblitz_strong_rho::{
    RawPointG, RhoField, RhoScalar, StrongRhoCharges, StrongRhoG, StrongRhoParams,
};
use crypto_lib::cryptanalysis::semaev_decomp::Gf2_128;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use rand::distributions::uniform::{SampleBorrow, SampleUniform, UniformInt, UniformSampler};
use rand::{rngs::StdRng, Rng, SeedableRng};
use serde::Deserialize;
use serde_json::{json, Value};

#[derive(Clone, Copy, Debug, Deserialize)]
#[serde(rename_all = "lowercase")]
enum Backend {
    Fold,
}

impl Backend {
    fn name(self) -> &'static str {
        match self {
            Self::Fold => "fold",
        }
    }
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[serde(rename_all = "lowercase")]
enum ScalarKernel {
    Reference,
    Barrett,
}

impl ScalarKernel {
    fn name(self) -> &'static str {
        match self {
            Self::Reference => "reference",
            Self::Barrett => "barrett",
        }
    }
}

#[derive(Clone, Copy, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Config {
    backend: Backend,
    scalar_kernel: ScalarKernel,
    lanes: usize,
    dp_bits: u32,
    escape_fruitless: bool,
    step_cap_factor: u64,
}

impl Config {
    fn parse(bytes: &[u8]) -> Result<Self, String> {
        let config: Self = serde_json::from_slice(bytes)
            .map_err(|error| format!("invalid strict config: {error}"))?;
        if ![64, 128, 256].contains(&config.lanes)
            || ![6, 8, 10].contains(&config.dp_bits)
            || config.escape_fruitless
            || config.step_cap_factor != 2_000
        {
            return Err(
                "config requires fold backend, reference|barrett scalars, lanes 64|128|256, dp_bits 6|8|10, escape_fruitless false and step_cap_factor 2000".into(),
            );
        }
        Ok(config)
    }

    fn params(self) -> StrongRhoParams {
        StrongRhoParams {
            lanes: self.lanes,
            dp_bits: self.dp_bits,
            escape_fruitless: self.escape_fruitless,
            step_cap_factor: self.step_cap_factor,
        }
    }

    fn as_json(self) -> Value {
        json!({
            "backend": self.backend.name(),
            "scalar_kernel": self.scalar_kernel.name(),
            "lanes": self.lanes,
            "dp_bits": self.dp_bits,
            "escape_fruitless": self.escape_fruitless,
            "step_cap_factor": self.step_cap_factor,
        })
    }
}

/// The IC producer's existing wide-field implementation, behind the same
/// RhoField contract. In particular, zeros remain zero in Montgomery inversion
/// and no replacement walk, canonicalization or scalar arithmetic is introduced.
struct FoldField(Gf2_128);

/// Canonical residue arithmetic used only by this worker. Every hot-path input
/// is in [0,m); the production walk preserves that invariant. Subtraction,
/// inversion and infrequent multiplication retain the production exact code.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
#[repr(transparent)]
struct ExactScalar(u64);

impl std::fmt::Display for ExactScalar {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

/// Delegation includes both one-shot methods: Rng::gen_range uses them rather
/// than the constructed sampler. This keeps rejection draws and RNG positions
/// identical to u64 for every supported range.
struct ExactUniform(UniformInt<u64>);

impl SampleUniform for ExactScalar {
    type Sampler = ExactUniform;
}

impl UniformSampler for ExactUniform {
    type X = ExactScalar;
    fn new<B1: SampleBorrow<Self::X>, B2: SampleBorrow<Self::X>>(low: B1, high: B2) -> Self {
        Self(UniformInt::new(low.borrow().0, high.borrow().0))
    }
    fn new_inclusive<B1: SampleBorrow<Self::X>, B2: SampleBorrow<Self::X>>(
        low: B1,
        high: B2,
    ) -> Self {
        Self(UniformInt::new_inclusive(low.borrow().0, high.borrow().0))
    }
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Self::X {
        ExactScalar(self.0.sample(rng))
    }
    fn sample_single<R: Rng + ?Sized, B1: SampleBorrow<Self::X>, B2: SampleBorrow<Self::X>>(
        low: B1,
        high: B2,
        rng: &mut R,
    ) -> Self::X {
        ExactScalar(UniformInt::<u64>::sample_single(
            low.borrow().0,
            high.borrow().0,
            rng,
        ))
    }
    fn sample_single_inclusive<
        R: Rng + ?Sized,
        B1: SampleBorrow<Self::X>,
        B2: SampleBorrow<Self::X>,
    >(
        low: B1,
        high: B2,
        rng: &mut R,
    ) -> Self::X {
        ExactScalar(UniformInt::<u64>::sample_single_inclusive(
            low.borrow().0,
            high.borrow().0,
            rng,
        ))
    }
}

#[derive(Clone, Copy, Debug)]
struct ExactMulCtx {
    modulus: u64,
    reciprocal: u128,
}

/// The upper 128 bits of the exact 256-bit product. Four 64x64 products
/// reconstruct its four limbs; the middle carry is at most two.
#[inline(always)]
fn mul_hi_128(a: u128, b: u128) -> u128 {
    let a0 = u128::from(a as u64);
    let a1 = a >> 64;
    let b0 = u128::from(b as u64);
    let b1 = b >> 64;
    let p00 = a0 * b0;
    let p01 = a0 * b1;
    let p10 = a1 * b0;
    let p11 = a1 * b1;
    let carry = (p00 >> 64) + u128::from(p01 as u64) + u128::from(p10 as u64);
    p11 + (p01 >> 64) + (p10 >> 64) + (carry >> 64)
}

impl ExactMulCtx {
    fn new(modulus: u64) -> Self {
        assert!(modulus > 1, "modulus must exceed 1");
        let m = u128::from(modulus);
        // floor(2^128/m), without representing 2^128 in u128. For m>=2
        // the result is at most 2^127; the addition therefore cannot overflow.
        let reciprocal = u128::MAX / m + u128::from(u128::MAX % m == m - 1);
        Self {
            modulus,
            reciprocal,
        }
    }

    #[inline(always)]
    fn mul(self, a: u64, b: u64) -> u64 {
        debug_assert!(a < self.modulus && b < self.modulus);
        let m = u128::from(self.modulus);
        let x = u128::from(a) * u128::from(b);
        let q = mul_hi_128(x, self.reciprocal);
        // q <= floor(x/m), so this subtraction is exact. Because x<2^128,
        // floor(x/m)-1 <= q; hence the residual is below 2m. One correction
        // suffices for every 2<=m<=u64::MAX, not only this 53-bit subgroup.
        let residual = x - q * m;
        let reduced = if residual >= m {
            residual - m
        } else {
            residual
        };
        debug_assert!(reduced < m);
        reduced as u64
    }
}

impl RhoScalar for ExactScalar {
    type MulCtx = ExactMulCtx;
    fn zero() -> Self {
        Self(0)
    }
    fn one() -> Self {
        Self(1)
    }
    fn from_biguint(v: &BigUint) -> Option<Self> {
        <u64 as RhoScalar>::from_biguint(v).map(Self)
    }
    fn to_f64(self) -> f64 {
        self.0 as f64
    }
    fn to_u128(self) -> u128 {
        u128::from(self.0)
    }
    fn bit_len(self) -> u32 {
        64 - self.0.leading_zeros()
    }
    fn bit(self, i: u32) -> bool {
        (self.0 >> i) & 1 == 1
    }
    fn mul_ctx(m: Self) -> Self::MulCtx {
        ExactMulCtx::new(m.0)
    }
    #[inline(always)]
    fn mul_in(ctx: &Self::MulCtx, a: Self, b: Self) -> Self {
        Self(ctx.mul(a.0, b.0))
    }
    fn mul_mod(a: Self, b: Self, m: Self) -> Self {
        Self(<u64 as RhoScalar>::mul_mod(a.0, b.0, m.0))
    }
    #[inline(always)]
    fn add_mod(a: Self, b: Self, m: Self) -> Self {
        debug_assert!(a.0 < m.0 && b.0 < m.0);
        let (sum, carry) = a.0.overflowing_add(b.0);
        Self(if carry || sum >= m.0 {
            sum.wrapping_sub(m.0)
        } else {
            sum
        })
    }
    fn sub_mod(a: Self, b: Self, m: Self) -> Self {
        Self(<u64 as RhoScalar>::sub_mod(a.0, b.0, m.0))
    }
    fn inverse_mod(v: Self, m: Self) -> Option<Self> {
        <u64 as RhoScalar>::inverse_mod(v.0, m.0).map(Self)
    }
}

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

fn arithmetic_kernel(scalar: ScalarKernel) -> Value {
    #[cfg(target_arch = "x86_64")]
    let hardware = std::arch::is_x86_feature_detected!("pclmulqdq").then_some("pclmulqdq");
    #[cfg(target_arch = "aarch64")]
    let hardware = std::arch::is_aarch64_feature_detected!("aes").then_some("pmull");
    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    let hardware: Option<&str> = None;
    json!({
        "backend": "fold",
        "carryless_multiply": hardware.unwrap_or("portable"),
        "reduction": if hardware.is_some() { "two_clmul_tail_folds" } else { "byte_table_reduction" },
        "squaring": "bit_spread_then_reduce",
        "scalar_kernel": scalar.name(),
        "scalar_multiply": match scalar { ScalarKernel::Reference => "production_u64_u128_remainder", ScalarKernel::Barrett => "exact_reciprocal_high128_one_correction" },
        "scalar_add": match scalar { ScalarKernel::Reference => "production_u128_remainder", ScalarKernel::Barrett => "canonical_overflow_aware_subtract" },
    })
}
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

fn make_strong<S: RhoScalar>(
    curve: &KoblitzCurve,
    field: FoldField,
) -> Result<StrongRhoG<FoldField, S>, String> {
    let r = curve
        .subgroup_order
        .to_u64()
        .ok_or("subgroup scalar does not fit u64")?;
    let lambda = curve.lambda.to_u64().ok_or("lambda does not fit u64")?;
    Ok(StrongRhoG::from_parts(
        field,
        u128::from(curve.a),
        S::from_biguint(&BigUint::from(r)).ok_or("scalar conversion failed")?,
        S::from_biguint(&BigUint::from(lambda)).ok_or("lambda conversion failed")?,
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

fn evaluate_scalar<S: RhoScalar>(
    input: &str,
    seed: u64,
    config: Config,
    config_bytes: &[u8],
) -> Result<(Value, bool), String> {
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
    let field = FoldField(Gf2_128::new(&curve.curve.irreducible));
    let kernel = arithmetic_kernel(config.scalar_kernel);
    let rho = make_strong::<S>(&curve, field)?;
    let mut setup_charges = StrongRhoCharges::default();
    // G-only jumps can be prepared before Q arrives. Walk-start RNG and jump
    // RNG domains are separate, both derived only from the declared seed.
    let mut jump_material = b"n83-exact-public-strong-jumps-v1\0".to_vec();
    jump_material.extend_from_slice(&seed.to_le_bytes());
    let digest = blake3::hash(&jump_material);
    let jump_seed = u64::from_le_bytes(digest.as_bytes()[..8].try_into().unwrap());
    let jumps = rho.jumps(jump_seed, &mut setup_charges);
    let prepare_ns = prepare_started.elapsed().as_nanos();
    let params = config.params();
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
        recovered.is_some_and(|d| curve.mul(curve.generator(), &BigUint::from(d.to_u128())) == q);
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
    let source = include_bytes!("n83_rho_speed_eval.rs");
    let receipt = json!({
        "schema": "n83-rho-speed-evaluation/v1",
        "status": status,
        "run_id": format!("N83A1RhoSpeed{}L{}DP{}E{}W{seed}Q{}", config.scalar_kernel.name(), config.lanes, config.dp_bits, u8::from(config.escape_fruitless), &input_blake3[..12]),
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
        "worker_source_sha256": sha256(source),
        "compiled_source_sha256": {
            "examples/n83_rho_speed_eval.rs": sha256(source),
            "src/cryptanalysis/koblitz_strong_rho.rs": sha256(include_bytes!("../src/cryptanalysis/koblitz_strong_rho.rs")),
            "src/cryptanalysis/wide_gf2m.rs": sha256(include_bytes!("../src/cryptanalysis/wide_gf2m.rs")),
            "src/cryptanalysis/semaev_decomp.rs": sha256(include_bytes!("../src/cryptanalysis/semaev_decomp.rs")),
        },
        "optimization_backend": config.backend.name(),
        "scalar_kernel": config.scalar_kernel.name(),
        "optimization_config": config.as_json(),
        "config_input_sha256": sha256(config_bytes),
        "config_input_blake3": blake3::hash(config_bytes).to_hex().to_string(),
        "executable_blake3": executable.as_deref().map(|b| blake3::hash(b).to_hex().to_string()),
        "executable_sha256": executable.as_deref().map(sha256),
        "public_input_blake3": input_blake3,
        "public_input_sha256": sha256(&input_bytes),
        "reference_grade": "strong_configured_fold_exact_u64_scalar",
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
            "excluded_stages": ["input_and_config_loading", "curve_normal_basis_G_only_jumps", "receipt_hashing", "receipt_serialization", "receipt_io"],
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

fn evaluate(input: &str, seed: u64, config_path: &str) -> Result<(Value, bool), String> {
    let config_bytes = fs::read(config_path).map_err(|error| format!("read config: {error}"))?;
    let config = Config::parse(&config_bytes)?;
    match config.scalar_kernel {
        ScalarKernel::Reference => evaluate_scalar::<u64>(input, seed, config, &config_bytes),
        ScalarKernel::Barrett => evaluate_scalar::<ExactScalar>(input, seed, config, &config_bytes),
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if !(4..=5).contains(&args.len()) {
        eprintln!(
            "usage: n83_rho_speed_eval <one-point.jsonl> <walk-seed> <config.json> [receipt.json]"
        );
        std::process::exit(2);
    }
    let result = args[2]
        .parse::<u64>()
        .map_err(|_| "walk seed must be u64".to_owned())
        .and_then(|seed| evaluate(&args[1], seed, &args[3]));
    let (receipt, success) = match result {
        Ok(value) => value,
        Err(message) => (
            json!({"schema":"n83-rho-speed-evaluation/v1", "status":"INPUT_OR_SETUP_FAILURE", "error":message}),
            false,
        ),
    };
    let bytes = format!("{}\n", serde_json::to_string(&receipt).unwrap());
    if let Some(path) = args.get(4) {
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
    use rand::distributions::{Distribution, Uniform};
    use rand::RngCore;

    #[test]
    fn strict_configs_accept_exact_grid_and_reject_other_input() {
        for kernel in ["reference", "barrett"] {
            for lanes in [64, 128, 256] {
                for dp in [6, 8, 10] {
                    let config = format!(
                        r#"{{"backend":"fold","scalar_kernel":"{kernel}","lanes":{lanes},"dp_bits":{dp},"escape_fruitless":false,"step_cap_factor":2000}}"#
                    );
                    assert!(Config::parse(config.as_bytes()).is_ok());
                }
            }
        }
        let valid = r#"{"backend":"fold","scalar_kernel":"barrett","lanes":64,"dp_bits":6,"escape_fruitless":false,"step_cap_factor":2000}"#;
        for invalid in [
            "{}".to_owned(),
            valid.replace("fold", "wide"),
            valid.replace("barrett", "fast"),
            valid.replace("64", "32"),
            valid.replace("\"dp_bits\":6", "\"dp_bits\":4"),
            valid.replace("false", "true"),
            valid.replace("2000", "1000"),
            valid.replace("\"scalar_kernel\":\"barrett\",", ""),
            valid.replace(
                "\"scalar_kernel\":\"barrett\",",
                "\"scalar_kernel\":\"barrett\",\"scalar_kernel\":\"reference\",",
            ),
            valid.replace("}", ",\"scalar\":12345}"),
            format!("{valid}{{}}"),
        ] {
            assert!(Config::parse(invalid.as_bytes()).is_err(), "{invalid}");
        }
    }

    #[test]
    fn public_input_refuses_scalar_sidecars_multiple_points_and_noncanonical_words() {
        for bad in [
            "12345",
            r#"{"fixture_scalar":12345}"#,
            "[1,2]",
            "[\"1\",\"2\"]\n[\"3\",\"4\"]",
        ] {
            assert!(parse_public_point(bad, 83).is_err());
        }
        assert!(parse_public_point(&format!("[\"{}\",\"2\"]", 1u128 << 83), 83).is_err());
        assert_eq!(parse_public_point("[\"1\",\"2\"]", 83).unwrap(), [1, 2]);
    }

    #[test]
    fn upper_product_matches_biguint_on_limb_edges_and_random_words() {
        let edge = [
            0,
            1,
            u128::from(u64::MAX),
            1u128 << 64,
            (1u128 << 64) + 1,
            1u128 << 127,
            u128::MAX,
        ];
        for a in edge {
            for b in edge {
                let expected = ((BigUint::from(a) * BigUint::from(b)) >> 128usize)
                    .to_u128()
                    .unwrap();
                assert_eq!(mul_hi_128(a, b), expected);
            }
        }
        let mut rng = StdRng::seed_from_u64(0x6e38336869676831);
        for _ in 0..16_384 {
            let a = rng.gen::<u128>();
            let b = rng.gen::<u128>();
            let expected = ((BigUint::from(a) * BigUint::from(b)) >> 128usize)
                .to_u128()
                .unwrap();
            assert_eq!(mul_hi_128(a, b), expected);
        }
    }

    #[test]
    fn reciprocal_reduction_and_addition_are_exact_with_one_correction_at_all_widths() {
        let mut moduli = vec![
            2,
            3,
            4,
            5,
            7,
            255,
            256,
            257,
            (1u64 << 32) - 1,
            1u64 << 32,
            (1u64 << 32) + 1,
            (1u64 << 50) - 1,
            1u64 << 50,
            (1u64 << 50) + 1,
            LEGACY_R,
            (1u64 << 63) - 1,
            1u64 << 63,
            (1u64 << 63) + 1,
            u64::MAX - 1,
            u64::MAX,
        ];
        let mut rng = StdRng::seed_from_u64(0x6e38337265647563);
        for _ in 0..256 {
            moduli.push(rng.gen_range(2..=u64::MAX));
        }
        let bound = BigUint::from(1u8) << 128usize;
        let mut saw_correction = false;
        let mut saw_carry = false;
        for m in moduli {
            let ctx = ExactMulCtx::new(m);
            assert_eq!(BigUint::from(ctx.reciprocal), &bound / BigUint::from(m));
            let edges = [0, 1, m / 2, m - 2, m - 1];
            let mut pairs: Vec<(u64, u64)> = edges
                .iter()
                .flat_map(|&a| edges.iter().map(move |&b| (a, b)))
                .collect();
            for _ in 0..512 {
                pairs.push((rng.gen_range(0..m), rng.gen_range(0..m)));
            }
            for (a, b) in pairs {
                let x = u128::from(a) * u128::from(b);
                let q = mul_hi_128(x, ctx.reciprocal);
                let exact_q = x / u128::from(m);
                assert!(q <= exact_q && exact_q - q <= 1);
                let residual = x - q * u128::from(m);
                assert!(residual < 2 * u128::from(m));
                saw_correction |= residual >= u128::from(m);
                assert_eq!(ctx.mul(a, b), (x % u128::from(m)) as u64);
                let expected_sum = ((u128::from(a) + u128::from(b)) % u128::from(m)) as u64;
                assert_eq!(
                    ExactScalar::add_mod(ExactScalar(a), ExactScalar(b), ExactScalar(m)).0,
                    expected_sum
                );
                saw_carry |= a.overflowing_add(b).1;
                assert_eq!(
                    ExactScalar::sub_mod(ExactScalar(a), ExactScalar(b), ExactScalar(m)).0,
                    <u64 as RhoScalar>::sub_mod(a, b, m)
                );
                assert_eq!(
                    ExactScalar::mul_mod(ExactScalar(a), ExactScalar(b), ExactScalar(m)).0,
                    <u64 as RhoScalar>::mul_mod(a, b, m)
                );
            }
        }
        assert!(saw_correction, "tests exercise the correction branch");
        assert!(saw_carry, "tests exercise overflow-aware addition");
    }

    #[test]
    fn delegated_inverse_and_scalar_conversions_match_reference() {
        for m in [2, 3, 8, 17, LEGACY_R, u64::MAX] {
            for a in [0, 1, 2 % m, m / 2, m - 1] {
                assert_eq!(
                    ExactScalar::inverse_mod(ExactScalar(a), ExactScalar(m)).map(|v| v.0),
                    <u64 as RhoScalar>::inverse_mod(a, m)
                );
            }
        }
        for a in [0, 1, LEGACY_R, u64::MAX] {
            let s = ExactScalar::from_biguint(&BigUint::from(a)).unwrap();
            assert_eq!(s.0, a);
            assert_eq!(s.to_u128(), u128::from(a));
            assert_eq!(s.to_string(), a.to_string());
            assert_eq!(s.bit_len(), <u64 as RhoScalar>::bit_len(a));
            for bit in 0..64 {
                assert_eq!(s.bit(bit), <u64 as RhoScalar>::bit(a, bit));
            }
        }
        assert!(ExactScalar::from_biguint(&(BigUint::from(1u8) << 64usize)).is_none());
    }

    #[test]
    fn uniform_samplers_and_one_shot_ranges_consume_identical_u64_rng_streams() {
        for (low, high) in [(0, 2), (1, LEGACY_R), (0, (1u64 << 63) + 1), (17, u64::MAX)] {
            let mut left = StdRng::seed_from_u64(0x6e383372616e6431);
            let mut right = left.clone();
            let u = Uniform::new(low, high);
            let e = Uniform::new(ExactScalar(low), ExactScalar(high));
            let ui = Uniform::new_inclusive(low, high);
            let ei = Uniform::new_inclusive(ExactScalar(low), ExactScalar(high));
            for _ in 0..4096 {
                assert_eq!(u.sample(&mut left), e.sample(&mut right).0);
                assert_eq!(ui.sample(&mut left), ei.sample(&mut right).0);
                assert_eq!(
                    left.gen_range(low..high),
                    right.gen_range(ExactScalar(low)..ExactScalar(high)).0
                );
                assert_eq!(
                    left.gen_range(low..=high),
                    right.gen_range(ExactScalar(low)..=ExactScalar(high)).0
                );
                assert_eq!(left.next_u64(), right.next_u64());
            }
        }
    }

    #[test]
    fn exact_n83_full_solve_jumps_trajectory_counters_and_charges_match_reference() {
        let curve = KoblitzCurve::new(1, 83).unwrap();
        assert_legacy_curve(&curve).unwrap();
        let reference =
            make_strong::<u64>(&curve, FoldField(Gf2_128::new(&curve.curve.irreducible))).unwrap();
        let exact =
            make_strong::<ExactScalar>(&curve, FoldField(Gf2_128::new(&curve.curve.irreducible)))
                .unwrap();
        let mut cr = StrongRhoCharges::default();
        let mut ce = StrongRhoCharges::default();
        // Correctness-only inputs; these seeds and Q are outside timed corpora.
        let jr = reference.jumps(0x6e38337363616c61, &mut cr);
        let je = exact.jumps(0x6e38337363616c61, &mut ce);
        assert_eq!(cr, ce);
        for (a, b) in jr.iter().zip(&je) {
            assert_eq!(a.point, b.point);
            assert_eq!(a.a, b.a.0);
        }
        for scalar in [1, 2, 12345, LEGACY_R - 1] {
            let point = reference.scalar_mul(reference.generator(), scalar);
            assert_eq!(
                point,
                exact.scalar_mul(exact.generator(), ExactScalar(scalar))
            );
            let sr = crypto_lib::cryptanalysis::koblitz_strong_rho::WalkStateG {
                point,
                a: scalar,
                b: 0,
            };
            let se = crypto_lib::cryptanalysis::koblitz_strong_rho::WalkStateG {
                point,
                a: ExactScalar(scalar),
                b: ExactScalar(0),
            };
            let a = reference.canonicalize(sr);
            let b = exact.canonicalize(se);
            assert_eq!(a.point, b.point);
            assert_eq!(a.a, b.a.0);
            assert_eq!(a.b, b.b.0);
        }
        let q = reference.scalar_mul(reference.generator(), 12345);
        let params = StrongRhoParams {
            lanes: 64,
            dp_bits: 6,
            escape_fruitless: false,
            step_cap_factor: 2000,
        };
        let mut rr = StdRng::seed_from_u64(0x6e3833736f6c7632);
        let mut re = rr.clone();
        let a = reference
            .solve(q, &jr, &mut rr, &params, cr)
            .expect("reference solve");
        let b = exact
            .solve(q, &je, &mut re, &params, ce)
            .expect("exact solve");
        assert_eq!(a.scalar, 12345);
        assert_eq!(a.scalar, b.scalar.0);
        assert_eq!(a.walk_steps, b.walk_steps);
        assert_eq!(a.walks, b.walks);
        assert_eq!(a.fruitless, b.fruitless);
        assert_eq!(a.escapes, b.escapes);
        assert_eq!(a.capped, b.capped);
        assert_eq!(a.wasted_merges, b.wasted_merges);
        assert_eq!(a.table_entries, b.table_entries);
        assert_eq!(a.ideal_steps, b.ideal_steps);
        assert_eq!(a.automorphisms, b.automorphisms);
        assert_eq!(a.charges, b.charges);
        assert_eq!(rr.next_u64(), re.next_u64());
        assert_eq!(
            RawPointG::<u128>::from_binary_wide(
                &curve.mul(curve.generator(), &BigUint::from(a.scalar))
            ),
            q
        );
    }
}

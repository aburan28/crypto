//! `cryptopro-b-glv-bench`: GLV-2 with an isogeny-chain endomorphism against
//! width-`w` NAF on GOST R 34.10-2001 CryptoPro-B, for
//! `research/cryptopro_b_glv_chain_20261005/`.
//!
//! Every arm runs on the same variable-time Jacobian `a = -3` arithmetic of
//! [`crypto_lib::ecc::cryptopro_b_point`]:
//!
//! * **baseline** — width-`w` NAF (`scalar_mul_wnaf`);
//! * **GLV-2 total**, one arm per chain in `CHAINS` (`scalar_mul_glv`):
//!   decomposition + `φ(P)` + interleaved width-`w` NAF, everything
//!   input-dependent inside the timed call;
//! * **stage diagnostics**, labelled as such and never a speed: the
//!   decomposition alone and `φ(P)` alone, per chain.
//!
//! Protocol.  A seeded RNG draws `M` (scalar, point) pairs.  Before any
//! timing the binary checks, on every pair, that each GLV arm equals the
//! baseline, that the baseline equals the crate's textbook `BigUint`
//! implementation, that `φ(P) = λ·P` for each chain, and that each
//! decomposition satisfies `k1 + k2·λ ≡ k (mod n)` within its Babai bound.
//! Then `R` rounds time every arm over all `M` pairs, the arm order
//! alternating between rounds, and the table reports the median and the
//! minimum over rounds of the per-operation time in nanoseconds: one
//! table, one unit, a correctness column.  The same record is written as
//! JSON with `--json`.  Nothing here is a constant-time implementation.

use clap::Parser;
use crypto_lib::ecc::cryptopro_b_chain_consts::{Chain, CHAINS};
use crypto_lib::ecc::cryptopro_b_point::{
    scalar_mul_glv, scalar_mul_wnaf, CryptoProBAffine, CryptoProBJacobian, GlvContext,
};
use crypto_lib::ecc::CurveParams;
use num_bigint::{BigUint, RandBigInt};
use num_traits::Zero;
use rand::rngs::StdRng;
use rand::SeedableRng;
use serde_json::{json, Value};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

#[derive(Parser, Debug)]
#[command(
    name = "cryptopro-b-glv-bench",
    about = "GLV-2 with an isogeny-chain endomorphism vs width-w NAF on GOST CryptoPro-B"
)]
struct Args {
    /// Number M of random (scalar, point) pairs.
    #[arg(long, default_value_t = 256)]
    pairs: usize,
    /// Number R of timing rounds.
    #[arg(long, default_value_t = 20)]
    rounds: usize,
    /// wNAF width, used by every arm.
    #[arg(long, default_value_t = 5)]
    w: u32,
    /// RNG seed for the pairs.
    #[arg(long, default_value_t = 20261005)]
    seed: u64,
    /// Write the full record here as JSON.
    #[arg(long)]
    json: Option<PathBuf>,
}

/// What an arm times, per pair index.
#[derive(Clone, Copy)]
enum Kind {
    Baseline,
    GlvTotal(usize),
    Decompose(usize),
    Phi(usize),
    /// The baseline again, as a separate arm: the A/A noise floor.
    BaselineAgain,
}

struct Arm {
    label: String,
    class: &'static str,
    kind: Kind,
    chain: Option<usize>,
    correctness: String,
    correct: bool,
}

fn nonzero_below(rng: &mut StdRng, n: &BigUint) -> BigUint {
    loop {
        let k = rng.gen_biguint_below(n);
        if !k.is_zero() {
            return k;
        }
    }
}

/// `5·5·7` for a chain's steps.
fn chain_shape(chain: &Chain) -> String {
    chain
        .steps
        .iter()
        .map(|s| s.ell.to_string())
        .collect::<Vec<_>>()
        .join("·")
}

/// The parenthetical of the chain's frozen name, e.g.
/// `element 4 + 1 w, degree 175`.
fn chain_element(chain: &Chain) -> String {
    match (chain.name.find('('), chain.name.rfind(')')) {
        (Some(a), Some(b)) if b > a => chain.name[a + 1..b].to_string(),
        _ => format!("degree {}", chain.degree),
    }
}

fn median(values: &[f64]) -> f64 {
    let mut v = values.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).expect("finite timings"));
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

fn min(values: &[f64]) -> f64 {
    values.iter().copied().fold(f64::INFINITY, f64::min)
}

fn shell(cmd: &str, args: &[&str]) -> Option<String> {
    let out = Command::new(cmd).args(args).output().ok()?;
    if !out.status.success() {
        return None;
    }
    Some(String::from_utf8_lossy(&out.stdout).trim().to_string())
}

fn cpu_info() -> (String, Vec<String>) {
    let text = std::fs::read_to_string("/proc/cpuinfo").unwrap_or_default();
    let mut model = String::from("unknown");
    let mut flags = Vec::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("model name") {
            if let Some((_, v)) = rest.split_once(':') {
                model = v.trim().to_string();
            }
        }
        if let Some(rest) = line.strip_prefix("flags") {
            if let Some((_, v)) = rest.split_once(':') {
                let of_interest = [
                    "sse4_2",
                    "popcnt",
                    "avx",
                    "avx2",
                    "avx512f",
                    "bmi1",
                    "bmi2",
                    "adx",
                    "pclmulqdq",
                    "aes",
                ];
                flags = v
                    .split_whitespace()
                    .filter(|f| of_interest.contains(f))
                    .map(str::to_string)
                    .collect();
            }
        }
        if model != "unknown" && !flags.is_empty() {
            break;
        }
    }
    (model, flags)
}

fn yes_no(b: bool) -> &'static str {
    if b {
        "yes"
    } else {
        "no"
    }
}

fn main() {
    let args = Args::parse();
    let (m, rounds, w) = (args.pairs, args.rounds, args.w);
    assert!(m > 0 && rounds > 0, "need at least one pair and one round");

    // ── curve-constant setup (untimed, like Montgomery-form curve constants) ──
    let curve = CurveParams::gost_cryptopro_b();
    let g = CryptoProBAffine::from_textbook(&curve.generator()).expect("generator");
    let ctxs: Vec<GlvContext> = CHAINS.iter().map(|c| GlvContext::new(c)).collect();

    // ── inputs ──
    let mut rng = StdRng::seed_from_u64(args.seed);
    let mut points: Vec<CryptoProBAffine> = Vec::with_capacity(m);
    let mut scalars: Vec<BigUint> = Vec::with_capacity(m);
    for _ in 0..m {
        let r = nonzero_below(&mut rng, &curve.n);
        points.push(
            scalar_mul_wnaf(&g, &r, w)
                .to_affine()
                .expect("r·G with 0 < r < n is not the identity"),
        );
        scalars.push(nonzero_below(&mut rng, &curve.n));
    }

    // ── correctness, on every pair, before any timing ──
    let a_fe = curve.a_fe();
    let baseline: Vec<CryptoProBJacobian> = (0..m)
        .map(|i| scalar_mul_wnaf(&points[i], &scalars[i], w))
        .collect();
    let mut baseline_vs_textbook = true;
    for i in 0..m {
        let want = points[i]
            .to_textbook()
            .scalar_mul_vartime(&scalars[i], &a_fe);
        baseline_vs_textbook &= baseline[i].to_textbook() == want;
    }
    struct ChainCheck {
        glv_ok: bool,
        phi_ok: bool,
        decomp_ok: bool,
    }
    let checks: Vec<ChainCheck> = ctxs
        .iter()
        .map(|ctx| {
            let mut c = ChainCheck {
                glv_ok: true,
                phi_ok: true,
                decomp_ok: true,
            };
            for i in 0..m {
                let glv = scalar_mul_glv(ctx, &points[i], &scalars[i], w);
                c.glv_ok &= glv.eq_point(&baseline[i]);
                let phi = ctx.chain.apply(&points[i]);
                c.phi_ok &= phi.eq_point(&scalar_mul_wnaf(&points[i], &ctx.lambda, w));
                let (k1, k2) = ctx.decompose(&scalars[i]);
                c.decomp_ok &= ctx.decomposition_is_valid(&scalars[i], &k1, &k2);
            }
            c
        })
        .collect();

    // ── arms ──
    let mut arms: Vec<Arm> = vec![Arm {
        label: format!("baseline: width-{w} NAF (affine odd-multiple table, mixed additions)"),
        class: "reference",
        kind: Kind::Baseline,
        chain: None,
        correctness: format!(
            "baseline == textbook BigUint reference on all {m} pairs: {}",
            yes_no(baseline_vs_textbook)
        ),
        correct: baseline_vs_textbook,
    }];
    for (c, chain) in CHAINS.iter().enumerate() {
        let shape = chain_shape(chain);
        let element = chain_element(chain);
        let bound = chain.babai_bound_bits;
        arms.push(Arm {
            label: format!(
                "GLV-2 total via chain {shape} ({element}): decomposition + phi(P) + interleaved width-{w} NAF"
            ),
            class: "engineering",
            kind: Kind::GlvTotal(c),
            chain: Some(c),
            correctness: format!(
                "GLV == baseline on all {m} pairs: {}",
                yes_no(checks[c].glv_ok)
            ),
            correct: checks[c].glv_ok,
        });
        arms.push(Arm {
            label: format!(
                "stage diagnostic: decomposition alone, basis of chain {shape} (num-bigint Babai rounding)"
            ),
            class: "stage diagnostic",
            kind: Kind::Decompose(c),
            chain: Some(c),
            correctness: format!(
                "k1 + k2·λ ≡ k (mod n) and |k1|,|k2| ≤ {bound} bits on all {m} scalars: {}",
                yes_no(checks[c].decomp_ok)
            ),
            correct: checks[c].decomp_ok,
        });
        arms.push(Arm {
            label: format!(
                "stage diagnostic: phi(P) alone, chain {shape} ({element}; projective Horner, no inversion)"
            ),
            class: "stage diagnostic",
            kind: Kind::Phi(c),
            chain: Some(c),
            correctness: format!(
                "phi(P) == λ·P on all {m} points: {}",
                yes_no(checks[c].phi_ok)
            ),
            correct: checks[c].phi_ok,
        });
    }
    arms.push(Arm {
        label: "A/A control: the baseline arm timed again as a separate arm (identical code and inputs; its ratio to the baseline is the noise floor)".to_string(),
        class: "A/A control",
        kind: Kind::BaselineAgain,
        chain: None,
        correctness: format!(
            "same routine as the baseline, so baseline == textbook BigUint reference on all {m} pairs: {}",
            yes_no(baseline_vs_textbook)
        ),
        correct: baseline_vs_textbook,
    });

    // ── timing ──
    let run = |kind: Kind, i: usize| match kind {
        Kind::Baseline | Kind::BaselineAgain => {
            black_box(scalar_mul_wnaf(&points[i], &scalars[i], w));
        }
        Kind::GlvTotal(c) => {
            black_box(scalar_mul_glv(&ctxs[c], &points[i], &scalars[i], w));
        }
        Kind::Decompose(c) => {
            black_box(ctxs[c].decompose(&scalars[i]));
        }
        Kind::Phi(c) => {
            black_box(ctxs[c].chain.apply(&points[i]));
        }
    };
    let mut per_round: Vec<Vec<f64>> = vec![Vec::with_capacity(rounds); arms.len()];
    let started = Instant::now();
    for round in 0..rounds {
        let order: Vec<usize> = if round % 2 == 0 {
            (0..arms.len()).collect()
        } else {
            (0..arms.len()).rev().collect()
        };
        for a in order {
            let t0 = Instant::now();
            for i in 0..m {
                run(arms[a].kind, i);
            }
            per_round[a].push(t0.elapsed().as_nanos() as f64 / m as f64);
        }
    }
    let timing_wall_s = started.elapsed().as_secs_f64();

    let medians: Vec<f64> = per_round.iter().map(|r| median(r)).collect();
    let mins: Vec<f64> = per_round.iter().map(|r| min(r)).collect();
    let (base_med, base_min) = (medians[0], mins[0]);

    // ── host and build record ──
    let (cpu_model, cpu_flags) = cpu_info();
    let rustc = shell("rustc", &["-V"]).unwrap_or_else(|| "unavailable".into());
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    let git_commit = shell("git", &["-C", manifest_dir, "rev-parse", "HEAD"]);
    let git_dirty =
        shell("git", &["-C", manifest_dir, "status", "--porcelain"]).map(|s| !s.is_empty());
    let profile = if cfg!(debug_assertions) {
        "debug"
    } else {
        "release"
    };
    let kernel = std::fs::read_to_string("/proc/sys/kernel/osrelease")
        .map(|s| s.trim().to_string())
        .unwrap_or_else(|_| "unknown".into());
    let logical_cpus = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(0);
    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);

    // ── the table ──
    println!("# cryptopro-b-glv-bench");
    println!();
    println!(
        "curve: GOST R 34.10-2001 CryptoPro-B (RFC 4357 section 11.4.4), p = 2^255 + 3225, prime order n, cofactor 1"
    );
    println!("arithmetic: variable-time Jacobian a = -3 (dbl-2001-b, add-2007-bl, madd-2007-bl), Montgomery-form field, shared by every arm");
    println!(
        "protocol: M = {m} seeded (scalar, point) pairs, R = {rounds} rounds, w = {w}, seed = {}, arm order alternates between rounds; timed output is a Jacobian point (no final inversion in any arm)",
        args.seed
    );
    println!(
        "host: {cpu_model}; flags: {}; logical cpus: {logical_cpus}; kernel: {kernel}; {} {}",
        cpu_flags.join(" "),
        std::env::consts::OS,
        std::env::consts::ARCH
    );
    println!(
        "build: {rustc}; profile {profile}; commit {}{}",
        git_commit.clone().unwrap_or_else(|| "unknown".into()),
        match git_dirty {
            Some(true) => " (dirty)",
            _ => "",
        }
    );
    println!("timing wall time: {timing_wall_s:.1} s");
    println!();
    println!("| arm | class | ns per scalar multiplication, median of {rounds} rounds | min of {rounds} rounds | ratio baseline/arm (median) | ratio baseline/arm (min) | correctness |");
    println!("|:--|:--|--:|--:|--:|--:|:--|");
    for (a, arm) in arms.iter().enumerate() {
        let (r_med, r_min) = match arm.kind {
            Kind::Baseline => ("1.000".to_string(), "1.000".to_string()),
            Kind::GlvTotal(_) | Kind::BaselineAgain => (
                format!("{:.3}", base_med / medians[a]),
                format!("{:.3}", base_min / mins[a]),
            ),
            Kind::Decompose(_) | Kind::Phi(_) => ("— (stage)".to_string(), "— (stage)".to_string()),
        };
        println!(
            "| {} | {} | {:.0} | {:.0} | {} | {} | {} |",
            arm.label, arm.class, medians[a], mins[a], r_med, r_min, arm.correctness
        );
    }
    println!();
    let all_correct = arms.iter().all(|a| a.correct);
    println!("all correctness checks passed: {}", yes_no(all_correct));

    // ── the JSON record ──
    let arms_json: Vec<Value> = arms
        .iter()
        .enumerate()
        .map(|(a, arm)| {
            let kind = match arm.kind {
                Kind::Baseline => "baseline",
                Kind::GlvTotal(_) => "glv_total",
                Kind::Decompose(_) => "stage_decomposition",
                Kind::Phi(_) => "stage_phi",
                Kind::BaselineAgain => "aa_control",
            };
            let ratios = match arm.kind {
                Kind::Baseline | Kind::GlvTotal(_) | Kind::BaselineAgain => json!({
                    "median": base_med / medians[a],
                    "min": base_min / mins[a],
                }),
                _ => Value::Null,
            };
            json!({
                "label": arm.label,
                "class": arm.class,
                "kind": kind,
                "chain": arm.chain.map(|c| CHAINS[c].name),
                "ns_per_op_per_round": per_round[a],
                "median_ns_per_op": medians[a],
                "min_ns_per_op": mins[a],
                "ratio_baseline_over_arm": ratios,
                "correctness": arm.correctness,
                "correct": arm.correct,
            })
        })
        .collect();
    let chains_json: Vec<Value> = CHAINS
        .iter()
        .map(|c| {
            json!({
                "name": c.name,
                "degree": c.degree,
                "steps": c.steps.iter().map(|s| s.ell).collect::<Vec<u32>>(),
                "babai_bound_bits": c.babai_bound_bits,
                "test_vectors": c.vectors.len(),
            })
        })
        .collect();
    let record = json!({
        "bench": "cryptopro-b-glv-bench",
        "curve": "GOST R 34.10-2001 CryptoPro-B (RFC 4357 section 11.4.4)",
        "unit": "nanoseconds per scalar multiplication (per operation for stage diagnostics)",
        "protocol": {
            "pairs_m": m,
            "rounds_r": rounds,
            "w": w,
            "seed": args.seed,
            "arm_order": "forward on even rounds, reversed on odd rounds",
            "statistic": "median and minimum over rounds of (round wall time / M)",
            "timed_output": "Jacobian point; the final affine conversion is outside every arm",
            "untimed_setup": "GlvContext per chain (Montgomery-form chain constants, lambda, reduced basis) built once",
            "timing_wall_seconds": timing_wall_s,
        },
        "host": {
            "cpu_model": cpu_model,
            "cpu_flags_of_interest": cpu_flags,
            "logical_cpus": logical_cpus,
            "os": std::env::consts::OS,
            "arch": std::env::consts::ARCH,
            "kernel": kernel,
        },
        "build": {
            "rustc": rustc,
            "profile": profile,
            "crate": env!("CARGO_PKG_NAME"),
            "crate_version": env!("CARGO_PKG_VERSION"),
            "git_commit": git_commit,
            "git_dirty": git_dirty,
        },
        "chains": chains_json,
        "correctness": {
            "baseline_equals_textbook_on_all_pairs": baseline_vs_textbook,
            "per_chain": checks.iter().zip(CHAINS.iter()).map(|(c, chain)| json!({
                "chain": chain.name,
                "glv_equals_baseline_on_all_pairs": c.glv_ok,
                "phi_equals_lambda_p_on_all_points": c.phi_ok,
                "decomposition_valid_on_all_scalars": c.decomp_ok,
            })).collect::<Vec<Value>>(),
            "all": all_correct,
        },
        "arms": arms_json,
        "timestamp_unix": timestamp,
    });
    if let Some(path) = &args.json {
        let text = serde_json::to_string_pretty(&record).expect("serialisable record");
        std::fs::write(path, text + "\n")
            .unwrap_or_else(|e| panic!("write {}: {e}", path.display()));
        println!("json record written to {}", path.display());
    }
    if !all_correct {
        eprintln!("correctness check failed; the timings above are not evidence");
        std::process::exit(1);
    }
}

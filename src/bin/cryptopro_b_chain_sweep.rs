//! `cryptopro-b-chain-sweep`: every isogeny-chain realisation of the cheap
//! endomorphisms of GOST R 34.10-2001 CryptoPro-B, with the
//! scalar-multiplication parameters around them, against width-`w` NAF,
//! for `research/cryptopro_b_chain_sweep_20261006/`.
//!
//! The chains come from a frozen JSON chain set (format
//! `endosweep-chainsweep/1`, loaded by
//! [`crypto_lib::ecc::cryptopro_b_chain_set`]).  Every arm runs on the same
//! variable-time Jacobian `a = -3` arithmetic of
//! [`crypto_lib::ecc::cryptopro_b_point`]:
//!
//! * **baseline** — width-`w` NAF for every `w` in `--widths`, with the
//!   odd-multiple table affine (one batched inversion, mixed additions) or
//!   Jacobian (no inversion, full additions);
//! * **GLV grid** — on `--grid-chain`: decomposition + `φ(P)` + interleaved
//!   two-scalar width-`w` NAF, for both chain evaluators (`generic`, the
//!   evaluator of `research/cryptopro_b_glv_chain_20261005`; `optimised`,
//!   Jacobian steps with an affine first step), both table modes and every
//!   width;
//! * **GLV per element** — the cheapest ordering (by the chain set's model)
//!   of every element in the set, at `--element-config`;
//! * **continuity** — `--continuity-chain` at `--continuity-config`, the
//!   configuration of the earlier measurement;
//! * **stage diagnostics**, labelled as such and never a speed: `φ(P)` alone
//!   for every chain in the set and both evaluators, and the decomposition
//!   alone for every element;
//! * **A/A control** — the reference baseline timed again as its own arm.
//!
//! Before any timing the binary checks, on every pair, that every baseline
//! arm equals the crate's textbook `BigUint` scalar multiplication, that
//! every GLV arm equals it too, that every `φ(P)` arm equals `λ·P` for its
//! chain's `λ`, that every decomposition satisfies `k1 + k2·λ ≡ k (mod n)`
//! within its Babai bound, and that every chain reproduces the test vectors
//! frozen with it.  The exit status is non-zero if any check fails.
//!
//! With `--count-only` (built with `--features cryptopro-b-opcount`) the
//! binary counts field operations per arm instead of timing, and writes them
//! as JSON; `--counts` merges such a file into a timing run's table, after
//! checking that it was made from the same chain set, pairs and seed.
//! Counted operations are deterministic and host-independent, so they are
//! the primary columns; wall time is the practicality columns.
//!
//! Nothing here is a constant-time implementation.

use clap::Parser;
use crypto_lib::ecc::cryptopro_b_chain_set::{ChainSet, ChainSpec};
use crypto_lib::ecc::cryptopro_b_point::{
    scalar_mul_glv_with, scalar_mul_wnaf, scalar_mul_wnaf_with, ChainEval, CryptoProBAffine,
    CryptoProBJacobian, GlvContext, TableMode,
};
use crypto_lib::ecc::CurveParams;
use crypto_lib::hash::sha256::sha256;
use num_bigint::{BigUint, RandBigInt};
use num_traits::Zero;
use rand::rngs::StdRng;
use rand::SeedableRng;
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

#[derive(Parser, Debug)]
#[command(
    name = "cryptopro-b-chain-sweep",
    about = "Every chain realisation and scalar-multiplication parameter of the CryptoPro-B GLV endomorphisms"
)]
struct Args {
    /// The frozen chain set (JSON, format endosweep-chainsweep/1).
    #[arg(long)]
    chains: PathBuf,
    /// Refuse to run unless the chain set's sha256 is this (hex).
    #[arg(long)]
    expect_sha256: Option<String>,
    /// Number M of random (scalar, point) pairs.
    #[arg(long, default_value_t = 256)]
    pairs: usize,
    /// Number R of timing rounds.
    #[arg(long, default_value_t = 20)]
    rounds: usize,
    /// RNG seed for the pairs.
    #[arg(long, default_value_t = 20261006)]
    seed: u64,
    /// wNAF widths of the baseline and the GLV grid.
    #[arg(long, value_delimiter = ',', default_value = "3,4,5,6,7")]
    widths: Vec<u32>,
    /// Chain whose full evaluator × table × width grid is timed (default:
    /// the least optimised model count in the set).
    #[arg(long)]
    grid_chain: Option<String>,
    /// "EVAL/TABLE/W" at which every element's cheapest ordering is timed.
    #[arg(long, default_value = "optimised/jacobian/5")]
    element_config: String,
    /// "TABLE/W" of the baseline that is the reference for every ratio and
    /// is timed a second time as the A/A control.
    #[arg(long, default_value = "jacobian/5")]
    reference: String,
    /// Chain id of the continuity arm.
    #[arg(long, default_value = "4+1w/5.5.7")]
    continuity_chain: String,
    /// "EVAL/TABLE/W" of the continuity arm.
    #[arg(long, default_value = "generic/affine/5")]
    continuity_config: String,
    /// Count field operations per arm instead of timing (needs the
    /// `cryptopro-b-opcount` feature).
    #[arg(long)]
    count_only: bool,
    /// Counts JSON from a `--count-only` run on the same inputs.
    #[arg(long)]
    counts: Option<PathBuf>,
    /// Write the full record here as JSON.
    #[arg(long)]
    json: Option<PathBuf>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Kind {
    Baseline {
        table: TableMode,
        w: u32,
    },
    Glv {
        chain: usize,
        eval: ChainEval,
        table: TableMode,
        w: u32,
    },
    Phi {
        chain: usize,
        eval: ChainEval,
    },
    Decompose {
        chain: usize,
    },
    AaControl {
        table: TableMode,
        w: u32,
    },
}

struct Arm {
    key: String,
    label: String,
    class: &'static str,
    kind: Kind,
    correctness: String,
    correct: bool,
}

fn eval_name(e: ChainEval) -> &'static str {
    match e {
        ChainEval::Generic => "generic",
        ChainEval::Optimised => "optimised",
    }
}

fn table_name(t: TableMode) -> &'static str {
    match t {
        TableMode::Affine => "affine",
        TableMode::Jacobian => "jacobian",
    }
}

fn parse_eval(s: &str) -> ChainEval {
    match s {
        "generic" => ChainEval::Generic,
        "optimised" => ChainEval::Optimised,
        _ => panic!("evaluator must be generic or optimised, not {s:?}"),
    }
}

fn parse_table(s: &str) -> TableMode {
    match s {
        "affine" => TableMode::Affine,
        "jacobian" => TableMode::Jacobian,
        _ => panic!("table mode must be affine or jacobian, not {s:?}"),
    }
}

fn parse_width(s: &str) -> u32 {
    let w: u32 = s
        .parse()
        .unwrap_or_else(|_| panic!("width {s:?} is not a number"));
    assert!((2..=7).contains(&w), "width must be in 2..=7");
    w
}

/// "EVAL/TABLE/W".
fn parse_config(s: &str) -> (ChainEval, TableMode, u32) {
    let parts: Vec<&str> = s.split('/').collect();
    assert_eq!(parts.len(), 3, "configuration {s:?} is not EVAL/TABLE/W");
    (
        parse_eval(parts[0]),
        parse_table(parts[1]),
        parse_width(parts[2]),
    )
}

/// "TABLE/W".
fn parse_baseline(s: &str) -> (TableMode, u32) {
    let parts: Vec<&str> = s.split('/').collect();
    assert_eq!(
        parts.len(),
        2,
        "baseline configuration {s:?} is not TABLE/W"
    );
    (parse_table(parts[0]), parse_width(parts[1]))
}

fn model_m_eq(c: &ChainSpec, eval: ChainEval) -> u64 {
    let m = match eval {
        ChainEval::Generic => c.model_generic,
        ChainEval::Optimised => c.model_optimised,
    };
    m.mul + m.sqr
}

fn nonzero_below(rng: &mut StdRng, n: &BigUint) -> BigUint {
    loop {
        let k = rng.gen_biguint_below(n);
        if !k.is_zero() {
            return k;
        }
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

/// Counted field operations of one arm, per scalar multiplication (or per
/// stage call), over all pairs.
#[derive(Clone, Copy, Debug, Default)]
struct Counted {
    mul: f64,
    sqr: f64,
    inv: f64,
    add: f64,
    m_eq_min: f64,
    m_eq_max: f64,
}

impl Counted {
    fn m_eq(&self) -> f64 {
        self.mul + self.sqr
    }
}

fn main() {
    let args = Args::parse();
    let (m, rounds) = (args.pairs, args.rounds);
    assert!(m > 0 && rounds > 0, "need at least one pair and one round");
    if args.count_only && !cfg!(feature = "cryptopro-b-opcount") {
        eprintln!("--count-only needs a build with --features cryptopro-b-opcount");
        std::process::exit(2);
    }
    if !args.count_only && cfg!(feature = "cryptopro-b-opcount") {
        eprintln!(
            "this build counts field operations (feature cryptopro-b-opcount); time with a default build"
        );
        std::process::exit(2);
    }

    // ── the frozen chain set ──
    let bytes = std::fs::read(&args.chains)
        .unwrap_or_else(|e| panic!("read {}: {e}", args.chains.display()));
    let chains_sha256 = hex::encode(sha256(&bytes));
    if let Some(want) = &args.expect_sha256 {
        if !want.eq_ignore_ascii_case(&chains_sha256) {
            eprintln!(
                "chain set sha256 {chains_sha256} differs from the expected {want}; refusing to run"
            );
            std::process::exit(2);
        }
    }
    let set = ChainSet::parse(std::str::from_utf8(&bytes).expect("UTF-8 chain set"))
        .unwrap_or_else(|e| panic!("chain set: {e}"));
    let ctxs: Vec<GlvContext> = set.chains.iter().map(ChainSpec::glv_context).collect();
    let index = |id: &str| -> usize {
        set.chains
            .iter()
            .position(|c| c.id == id)
            .unwrap_or_else(|| panic!("chain {id:?} is not in the set"))
    };
    let grid_chain = match &args.grid_chain {
        Some(id) => index(id),
        None => (0..set.chains.len())
            .min_by_key(|&c| {
                (
                    model_m_eq(&set.chains[c], ChainEval::Optimised),
                    set.chains[c].id.clone(),
                )
            })
            .expect("non-empty set"),
    };
    // each element's cheapest ordering under the optimised model
    let mut by_element: BTreeMap<(i64, i64), usize> = BTreeMap::new();
    for (c, ch) in set.chains.iter().enumerate() {
        let e = by_element.entry(ch.element).or_insert(c);
        let cur = &set.chains[*e];
        if (model_m_eq(ch, ChainEval::Optimised), &ch.id)
            < (model_m_eq(cur, ChainEval::Optimised), &cur.id)
        {
            *e = c;
        }
    }
    let mut element_chains: Vec<usize> = by_element.values().copied().collect();
    element_chains.sort_by_key(|&c| {
        (
            model_m_eq(&set.chains[c], ChainEval::Optimised),
            set.chains[c].id.clone(),
        )
    });
    let (el_eval, el_table, el_w) = parse_config(&args.element_config);
    let (ref_table, ref_w) = parse_baseline(&args.reference);
    let cont_chain = index(&args.continuity_chain);
    let (cont_eval, cont_table, cont_w) = parse_config(&args.continuity_config);
    for &w in &args.widths {
        parse_width(&w.to_string());
    }

    // ── inputs ──
    let curve = CurveParams::gost_cryptopro_b();
    let g = CryptoProBAffine::from_textbook(&curve.generator()).expect("generator");
    let mut rng = StdRng::seed_from_u64(args.seed);
    let mut points: Vec<CryptoProBAffine> = Vec::with_capacity(m);
    let mut scalars: Vec<BigUint> = Vec::with_capacity(m);
    for _ in 0..m {
        let r = nonzero_below(&mut rng, &curve.n);
        points.push(
            scalar_mul_wnaf(&g, &r, 5)
                .to_affine()
                .expect("r·G with 0 < r < n is not the identity"),
        );
        scalars.push(nonzero_below(&mut rng, &curve.n));
    }

    // ── arms ──
    let mut arms: Vec<Arm> = Vec::new();
    let tables = [TableMode::Affine, TableMode::Jacobian];
    let evals = [ChainEval::Generic, ChainEval::Optimised];
    for &table in &tables {
        for &w in &args.widths {
            let is_ref = table == ref_table && w == ref_w;
            arms.push(Arm {
                key: format!("baseline/{}/w{w}", table_name(table)),
                label: format!(
                    "baseline: width-{w} NAF, {} table{}",
                    match table {
                        TableMode::Affine => "affine (one batched inversion, mixed additions)",
                        TableMode::Jacobian => "Jacobian (no inversion, full additions)",
                    },
                    if is_ref { " — the reference" } else { "" }
                ),
                class: if is_ref {
                    "reference"
                } else {
                    "baseline variant"
                },
                kind: Kind::Baseline { table, w },
                correctness: String::new(),
                correct: false,
            });
        }
    }
    let glv_label = |c: usize, eval: ChainEval, table: TableMode, w: u32, what: &str| -> String {
        let ch = &set.chains[c];
        format!(
            "GLV-2 total{what} via chain {} ({}, norm {}, {} evaluator): decomposition + phi(P) + interleaved width-{w} NAF, {} tables",
            ch.shape(),
            ch.element_label(),
            ch.norm,
            eval_name(eval),
            table_name(table)
        )
    };
    let mut glv_keys: Vec<(usize, ChainEval, TableMode, u32)> = Vec::new();
    for &eval in &evals {
        for &table in &tables {
            for &w in &args.widths {
                glv_keys.push((grid_chain, eval, table, w));
                arms.push(Arm {
                    key: format!(
                        "glv/{}/{}/{}/w{w}",
                        set.chains[grid_chain].id,
                        eval_name(eval),
                        table_name(table)
                    ),
                    label: glv_label(grid_chain, eval, table, w, " (grid)"),
                    class: "engineering",
                    kind: Kind::Glv {
                        chain: grid_chain,
                        eval,
                        table,
                        w,
                    },
                    correctness: String::new(),
                    correct: false,
                });
            }
        }
    }
    let mut extra =
        |c: usize, eval: ChainEval, table: TableMode, w: u32, what: &str, arms: &mut Vec<Arm>| {
            if glv_keys.contains(&(c, eval, table, w)) {
                return;
            }
            glv_keys.push((c, eval, table, w));
            arms.push(Arm {
                key: format!(
                    "glv/{}/{}/{}/w{w}",
                    set.chains[c].id,
                    eval_name(eval),
                    table_name(table)
                ),
                label: glv_label(c, eval, table, w, what),
                class: "engineering",
                kind: Kind::Glv {
                    chain: c,
                    eval,
                    table,
                    w,
                },
                correctness: String::new(),
                correct: false,
            });
        };
    for &c in &element_chains {
        extra(c, el_eval, el_table, el_w, " (per element)", &mut arms);
    }
    extra(
        cont_chain,
        cont_eval,
        cont_table,
        cont_w,
        " (continuity: the earlier measurement's configuration)",
        &mut arms,
    );
    for c in 0..set.chains.len() {
        for &eval in &evals {
            let ch = &set.chains[c];
            arms.push(Arm {
                key: format!("phi/{}/{}", ch.id, eval_name(eval)),
                label: format!(
                    "stage diagnostic: phi(P) alone, chain {} ({}, norm {}), {} evaluator (model {} M_eq)",
                    ch.shape(),
                    ch.element_label(),
                    ch.norm,
                    eval_name(eval),
                    model_m_eq(ch, eval)
                ),
                class: "stage diagnostic",
                kind: Kind::Phi { chain: c, eval },
                correctness: String::new(),
                correct: false,
            });
        }
    }
    for &c in &element_chains {
        let ch = &set.chains[c];
        arms.push(Arm {
            key: format!("decompose/{}", ch.id),
            label: format!(
                "stage diagnostic: decomposition alone, basis of {} (num-bigint Babai rounding)",
                ch.element_label()
            ),
            class: "stage diagnostic",
            kind: Kind::Decompose { chain: c },
            correctness: String::new(),
            correct: false,
        });
    }
    arms.push(Arm {
        key: format!("aa/{}/w{ref_w}", table_name(ref_table)),
        label: "A/A control: the reference baseline timed again as a separate arm (identical code and inputs; its ratio to the reference is the noise floor)".to_string(),
        class: "A/A control",
        kind: Kind::AaControl {
            table: ref_table,
            w: ref_w,
        },
        correctness: String::new(),
        correct: false,
    });

    // ── correctness, on every pair, before any timing ──
    let a_fe = curve.a_fe();
    let truth: Vec<_> = (0..m)
        .map(|i| {
            points[i]
                .to_textbook()
                .scalar_mul_vartime(&scalars[i], &a_fe)
        })
        .collect();
    let reference: Vec<CryptoProBJacobian> = (0..m)
        .map(|i| scalar_mul_wnaf(&points[i], &scalars[i], 5))
        .collect();
    let reference_ok = (0..m).all(|i| reference[i].to_textbook() == truth[i]);
    // λ·P per chain, from the verified wNAF routine
    let lambda_p: Vec<Vec<CryptoProBJacobian>> = ctxs
        .iter()
        .map(|ctx| {
            (0..m)
                .map(|i| scalar_mul_wnaf(&points[i], &ctx.lambda, 5))
                .collect()
        })
        .collect();
    // frozen test vectors, per chain, both evaluators
    let vectors_ok: Vec<bool> = set
        .chains
        .iter()
        .zip(ctxs.iter())
        .map(|(ch, ctx)| {
            ch.vectors.iter().all(|v| {
                let (p, phi_p, k_p) = ChainSpec::vector_points(v);
                let phi_want = CryptoProBJacobian::from_affine(&phi_p);
                let kp_want = CryptoProBJacobian::from_affine(&k_p);
                let (k1, k2) = ctx.decompose(&v.k);
                evals
                    .iter()
                    .all(|&e| ctx.chain.apply_with(&p, e).eq_point(&phi_want))
                    && k1 == v.k1
                    && k2 == v.k2
                    && scalar_mul_glv_with(
                        ctx,
                        &p,
                        &v.k,
                        5,
                        TableMode::Affine,
                        ChainEval::Optimised,
                    )
                    .eq_point(&kp_want)
            })
        })
        .collect();
    for arm in arms.iter_mut() {
        let (ok, text) = match arm.kind {
            Kind::Baseline { table, w } | Kind::AaControl { table, w } => {
                let ok = reference_ok
                    && (0..m).all(|i| {
                        scalar_mul_wnaf_with(&points[i], &scalars[i], w, table).to_textbook()
                            == truth[i]
                    });
                (
                    ok,
                    format!("== textbook BigUint k·P on all {m} pairs: {}", yes_no(ok)),
                )
            }
            Kind::Glv {
                chain,
                eval,
                table,
                w,
            } => {
                let ok = reference_ok
                    && (0..m).all(|i| {
                        scalar_mul_glv_with(&ctxs[chain], &points[i], &scalars[i], w, table, eval)
                            .eq_point(&reference[i])
                    });
                (
                    ok,
                    format!("GLV == textbook k·P on all {m} pairs: {}", yes_no(ok)),
                )
            }
            Kind::Phi { chain, eval } => {
                let ok = vectors_ok[chain]
                    && (0..m).all(|i| {
                        ctxs[chain]
                            .chain
                            .apply_with(&points[i], eval)
                            .eq_point(&lambda_p[chain][i])
                    });
                (
                    ok,
                    format!(
                        "phi(P) == λ·P on all {m} points and the {} frozen test vectors: {}",
                        set.chains[chain].vectors.len(),
                        yes_no(ok)
                    ),
                )
            }
            Kind::Decompose { chain } => {
                let ctx = &ctxs[chain];
                let ok = (0..m).all(|i| {
                    let (k1, k2) = ctx.decompose(&scalars[i]);
                    ctx.decomposition_is_valid(&scalars[i], &k1, &k2)
                });
                (
                    ok,
                    format!(
                        "k1 + k2·λ ≡ k (mod n) and |k1|,|k2| ≤ {} bits on all {m} scalars: {}",
                        ctx.babai_bound_bits,
                        yes_no(ok)
                    ),
                )
            }
        };
        arm.correct = ok;
        arm.correctness = text;
    }
    let all_correct = arms.iter().all(|a| a.correct) && vectors_ok.iter().all(|&v| v);

    let run = |kind: Kind, i: usize| match kind {
        Kind::Baseline { table, w } | Kind::AaControl { table, w } => {
            black_box(scalar_mul_wnaf_with(&points[i], &scalars[i], w, table));
        }
        Kind::Glv {
            chain,
            eval,
            table,
            w,
        } => {
            black_box(scalar_mul_glv_with(
                &ctxs[chain],
                &points[i],
                &scalars[i],
                w,
                table,
                eval,
            ));
        }
        Kind::Phi { chain, eval } => {
            black_box(ctxs[chain].chain.apply_with(&points[i], eval));
        }
        Kind::Decompose { chain } => {
            black_box(ctxs[chain].decompose(&scalars[i]));
        }
    };

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
    let host = json!({
        "cpu_model": cpu_model,
        "cpu_flags_of_interest": cpu_flags,
        "logical_cpus": logical_cpus,
        "os": std::env::consts::OS,
        "arch": std::env::consts::ARCH,
        "kernel": kernel,
    });
    let build = json!({
        "rustc": rustc,
        "profile": profile,
        "crate": env!("CARGO_PKG_NAME"),
        "crate_version": env!("CARGO_PKG_VERSION"),
        "git_commit": git_commit,
        "git_dirty": git_dirty,
        "opcount_feature": cfg!(feature = "cryptopro-b-opcount"),
    });
    let inputs = json!({
        "chains_path": args.chains.display().to_string(),
        "chains_sha256": chains_sha256,
        "chains_in_set": set.chains.len(),
        "pairs_m": m,
        "seed": args.seed,
        "widths": args.widths,
        "grid_chain": set.chains[grid_chain].id,
        "element_config": args.element_config,
        "element_chains": element_chains.iter().map(|&c| set.chains[c].id.clone()).collect::<Vec<_>>(),
        "reference": args.reference,
        "continuity_chain": args.continuity_chain,
        "continuity_config": args.continuity_config,
    });

    // ── counting mode ──
    if args.count_only {
        let counted = count_arms(&arms, m, &run);
        let rows: Vec<Value> = arms
            .iter()
            .zip(counted.iter())
            .map(|(arm, c)| {
                json!({
                    "key": arm.key,
                    "mul": c.mul, "sqr": c.sqr, "inv": c.inv, "add": c.add,
                    "m_eq": c.m_eq(), "m_eq_min": c.m_eq_min, "m_eq_max": c.m_eq_max,
                    "correct": arm.correct,
                })
            })
            .collect();
        let record = json!({
            "bench": "cryptopro-b-chain-sweep",
            "mode": "count",
            "unit": "field operations per scalar multiplication (per call for stage diagnostics), mean over all pairs",
            "inputs": inputs,
            "host": host,
            "build": build,
            "all_correct": all_correct,
            "arms": rows,
            "timestamp_unix": timestamp,
        });
        println!("# cryptopro-b-chain-sweep (counting run)");
        println!();
        println!("| arm | M | S | inversions | M_eq = M + S (mean) | min | max | correct |");
        println!("|:--|--:|--:|--:|--:|--:|--:|:--|");
        for (arm, c) in arms.iter().zip(counted.iter()) {
            println!(
                "| {} | {:.2} | {:.2} | {:.3} | {:.2} | {:.0} | {:.0} | {} |",
                arm.key,
                c.mul,
                c.sqr,
                c.inv,
                c.m_eq(),
                c.m_eq_min,
                c.m_eq_max,
                yes_no(arm.correct)
            );
        }
        println!();
        println!("all correctness checks passed: {}", yes_no(all_correct));
        write_json(&args.json, &record);
        if !all_correct {
            eprintln!("correctness check failed; the counts above are not evidence");
            std::process::exit(1);
        }
        return;
    }

    // ── counts from a counting run, if given ──
    let counts: Option<BTreeMap<String, Counted>> = args.counts.as_ref().map(|path| {
        let text = std::fs::read_to_string(path)
            .unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
        let v: Value = serde_json::from_str(&text).expect("counts JSON");
        let vin = &v["inputs"];
        assert_eq!(
            vin["chains_sha256"],
            json!(chains_sha256),
            "counts were made from another chain set"
        );
        assert_eq!(
            vin["pairs_m"],
            json!(m),
            "counts were made from another number of pairs"
        );
        assert_eq!(
            vin["seed"],
            json!(args.seed),
            "counts were made from another seed"
        );
        assert_eq!(
            v["all_correct"],
            json!(true),
            "the counting run failed a correctness check"
        );
        v["arms"]
            .as_array()
            .expect("arms")
            .iter()
            .map(|a| {
                let f = |k: &str| a[k].as_f64().expect("number");
                (
                    a["key"].as_str().expect("key").to_string(),
                    Counted {
                        mul: f("mul"),
                        sqr: f("sqr"),
                        inv: f("inv"),
                        add: f("add"),
                        m_eq_min: f("m_eq_min"),
                        m_eq_max: f("m_eq_max"),
                    },
                )
            })
            .collect()
    });
    if let Some(c) = &counts {
        for arm in &arms {
            assert!(
                c.contains_key(&arm.key),
                "the counting run has no arm {}",
                arm.key
            );
        }
    }

    // ── timing ──
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
    let ref_idx = arms
        .iter()
        .position(|a| a.class == "reference")
        .expect("reference arm");
    let (ref_med, ref_min) = (medians[ref_idx], mins[ref_idx]);
    let ref_count = counts.as_ref().map(|c| c[&arms[ref_idx].key].m_eq());
    let is_total = |k: Kind| !matches!(k, Kind::Phi { .. } | Kind::Decompose { .. });

    // ── the table ──
    println!("# cryptopro-b-chain-sweep");
    println!();
    println!(
        "curve: GOST R 34.10-2001 CryptoPro-B (RFC 4357 section 11.4.4), p = 2^255 + 3225, prime order n, cofactor 1; CM by the maximal order of discriminant {} (class number {})",
        set.d_k,
        set.class_number.map(|h| h.to_string()).unwrap_or_else(|| "?".into())
    );
    println!("arithmetic: variable-time Jacobian a = -3 (dbl-2001-b, add-2007-bl, madd-2007-bl), Montgomery-form field, shared by every arm");
    println!(
        "chain set: {} ({} chains, sha256 {chains_sha256}); grid chain {}; per-element configuration {}; reference baseline {}",
        args.chains.display(),
        set.chains.len(),
        set.chains[grid_chain].id,
        args.element_config,
        args.reference
    );
    println!(
        "protocol: M = {m} seeded (scalar, point) pairs, R = {rounds} rounds, seed = {}, arm order alternates between rounds; timed output is a Jacobian point (no final inversion in any arm)",
        args.seed
    );
    println!(
        "host: {}; flags: {}; logical cpus: {logical_cpus}; kernel: {}; {} {}",
        host["cpu_model"].as_str().unwrap_or("?"),
        cpu_flags.join(" "),
        host["kernel"].as_str().unwrap_or("?"),
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
    println!(
        "counted operations: {}",
        match &args.counts {
            Some(p) => format!("from {} (a counting build on the same inputs)", p.display()),
            None => "not supplied (columns marked —)".into(),
        }
    );
    println!("timing wall time: {timing_wall_s:.1} s");
    println!();
    println!("| arm | class | counted M_eq per op (mean; M + S) | M | S | inversions | ratio reference/arm (counted) | ns per op, median of {rounds} | min of {rounds} | ratio reference/arm (median) | ratio reference/arm (min) | correctness |");
    println!("|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|");
    for (a, arm) in arms.iter().enumerate() {
        let c = counts.as_ref().map(|c| c[&arm.key]);
        let cnum = |f: fn(&Counted) -> f64, digits: usize| match c {
            Some(c) => format!("{:.*}", digits, f(&c)),
            None => "—".into(),
        };
        let total = is_total(arm.kind);
        let rc = match (c, ref_count, total) {
            (Some(c), Some(r), true) => format!("{:.3}", r / c.m_eq()),
            (_, _, false) => "— (stage)".into(),
            _ => "—".into(),
        };
        let (rm, rn) = if total {
            (
                format!("{:.3}", ref_med / medians[a]),
                format!("{:.3}", ref_min / mins[a]),
            )
        } else {
            ("— (stage)".into(), "— (stage)".into())
        };
        println!(
            "| {} | {} | {} | {} | {} | {} | {} | {:.0} | {:.0} | {} | {} | {} |",
            arm.label,
            arm.class,
            cnum(Counted::m_eq_fn, 1),
            cnum(|c| c.mul, 1),
            cnum(|c| c.sqr, 1),
            cnum(|c| c.inv, 0),
            rc,
            medians[a],
            mins[a],
            rm,
            rn,
            arm.correctness
        );
    }
    println!();
    println!("all correctness checks passed: {}", yes_no(all_correct));

    let arms_json: Vec<Value> = arms
        .iter()
        .enumerate()
        .map(|(a, arm)| {
            let c = counts.as_ref().map(|c| c[&arm.key]);
            let total = is_total(arm.kind);
            json!({
                "key": arm.key,
                "label": arm.label,
                "class": arm.class,
                "ns_per_op_per_round": per_round[a],
                "median_ns_per_op": medians[a],
                "min_ns_per_op": mins[a],
                "ratio_reference_over_arm": if total { json!({
                    "median": ref_med / medians[a],
                    "min": ref_min / mins[a],
                    "counted": match (c, ref_count) { (Some(c), Some(r)) => json!(r / c.m_eq()), _ => Value::Null },
                }) } else { Value::Null },
                "counted": c.map(|c| json!({
                    "mul": c.mul, "sqr": c.sqr, "inv": c.inv, "add": c.add,
                    "m_eq": c.m_eq(), "m_eq_min": c.m_eq_min, "m_eq_max": c.m_eq_max,
                })),
                "correctness": arm.correctness,
                "correct": arm.correct,
            })
        })
        .collect();
    let record = json!({
        "bench": "cryptopro-b-chain-sweep",
        "mode": "timing",
        "curve": "GOST R 34.10-2001 CryptoPro-B (RFC 4357 section 11.4.4)",
        "unit": "nanoseconds per scalar multiplication (per call for stage diagnostics); counted field operations per op from the counting run",
        "protocol": {
            "rounds_r": rounds,
            "arm_order": "forward on even rounds, reversed on odd rounds",
            "statistic": "median and minimum over rounds of (round wall time / M)",
            "timed_output": "Jacobian point; the final affine conversion is outside every arm",
            "untimed_setup": "GlvContext per chain (Montgomery-form chain constants, u^-1, lambda, reduced basis) built once",
            "timing_wall_seconds": timing_wall_s,
        },
        "inputs": inputs,
        "host": host,
        "build": build,
        "chains": set.chains.iter().map(|c| json!({
            "id": c.id,
            "element": [c.element.0, c.element.1],
            "norm": c.norm,
            "order": c.order,
            "babai_bound_bits": c.babai_bound_bits,
            "model_m_eq": {"generic": model_m_eq(c, ChainEval::Generic), "optimised": model_m_eq(c, ChainEval::Optimised)},
        })).collect::<Vec<Value>>(),
        "correctness": {
            "reference_wnaf_equals_textbook_on_all_pairs": reference_ok,
            "test_vectors_per_chain": set.chains.iter().zip(vectors_ok.iter()).map(|(c, ok)| json!({"chain": c.id, "ok": ok})).collect::<Vec<Value>>(),
            "all": all_correct,
        },
        "arms": arms_json,
        "timestamp_unix": timestamp,
    });
    write_json(&args.json, &record);
    if !all_correct {
        eprintln!("correctness check failed; the timings above are not evidence");
        std::process::exit(1);
    }
}

impl Counted {
    fn m_eq_fn(c: &Counted) -> f64 {
        c.m_eq()
    }
}

fn write_json(path: &Option<PathBuf>, record: &Value) {
    if let Some(path) = path {
        let text = serde_json::to_string_pretty(record).expect("serialisable record");
        std::fs::write(path, text + "\n")
            .unwrap_or_else(|e| panic!("write {}: {e}", path.display()));
        println!("json record written to {}", path.display());
    }
}

/// Per arm: run every pair once with the counters reset before each pair.
#[cfg(feature = "cryptopro-b-opcount")]
fn count_arms(arms: &[Arm], m: usize, run: &dyn Fn(Kind, usize)) -> Vec<Counted> {
    use crypto_lib::ecc::cryptopro_b_field::opcount;
    arms.iter()
        .map(|arm| {
            let mut c = Counted {
                m_eq_min: f64::INFINITY,
                m_eq_max: 0.0,
                ..Counted::default()
            };
            for i in 0..m {
                opcount::reset();
                run(arm.kind, i);
                let o = opcount::read();
                c.mul += o.mul as f64;
                c.sqr += o.sqr as f64;
                c.inv += o.inv as f64;
                c.add += o.add as f64;
                let e = (o.mul + o.sqr) as f64;
                c.m_eq_min = c.m_eq_min.min(e);
                c.m_eq_max = c.m_eq_max.max(e);
            }
            let mf = m as f64;
            c.mul /= mf;
            c.sqr /= mf;
            c.inv /= mf;
            c.add /= mf;
            c
        })
        .collect()
}

#[cfg(not(feature = "cryptopro-b-opcount"))]
fn count_arms(_arms: &[Arm], _m: usize, _run: &dyn Fn(Kind, usize)) -> Vec<Counted> {
    unreachable!("--count-only is refused without the cryptopro-b-opcount feature")
}

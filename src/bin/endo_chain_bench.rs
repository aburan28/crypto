//! `endo-chain-bench`: GLV-2 via a frozen isogeny-chain endomorphism against
//! width-`w` NAF, on any short Weierstrass curve over a prime field, for
//! `research/endo_native_bench_20261016/`.
//!
//! The curve-generic form of `cryptopro-b-chain-sweep` (whose code and
//! results it leaves untouched): the chain set (format
//! `endosweep-chainsweep/1`, loaded by
//! [`crypto_lib::ecc::endo_generic_chain_set`]) carries the curve's own `p`,
//! `a`, `b`, `n`, cofactor, generator, chains, `λ`, reduced bases and test
//! vectors, and every arm runs on the same variable-time Jacobian arithmetic
//! of [`crypto_lib::ecc::endo_generic_point`] over a run-time Montgomery
//! field of 4, 6 or 13 limbs, with the curve's own doubling formula:
//!
//! * **baseline** — width-`w` NAF for every `w` in `--widths`, table affine
//!   (one batched inversion, mixed additions) or Jacobian (full additions);
//! * **GLV grid** — on `--grid-chain` (default: the least `optimised` model
//!   count): both chain evaluators × both table modes × every width;
//! * **GLV per element** — every element's cheapest ordering at
//!   `--element-config`;
//! * **continuity** (optional) — `--continuity-chain` at
//!   `--continuity-config`;
//! * **stage diagnostics**, never a speed: `φ(P)` alone for every chain and
//!   both evaluators, and the decomposition alone for every element;
//! * **A/A control** — the reference baseline timed again as its own arm.
//!
//! Arm keys, pair generation (seeded `StdRng`, `r ← [1, n)`, `P = r·G`,
//! `k ← [1, n)`), the counting protocol and the timing loop are those of
//! `cryptopro-b-chain-sweep`, so that on the CryptoPro-B chain set the two
//! harnesses can be compared arm by arm.
//!
//! Before any timing the binary checks on every pair that every baseline arm
//! equals the crate's textbook `BigUint` `k·P`, every GLV arm too, every
//! `φ(P)` equals `λ·P`, every decomposition is valid within its Babai bound,
//! and every chain reproduces its frozen test vectors; it exits non-zero if
//! any check fails.  `--count-only` (built with `--features
//! endo-chain-opcount`) counts field operations per arm instead of timing;
//! `--counts` merges such a file into a timing run's table.
//!
//! Nothing here is a constant-time implementation.

use clap::Parser;
use crypto_lib::ecc::endo_generic_chain_set::{ChainSet, ChainSpec, PreparedSet};
use crypto_lib::ecc::endo_generic_point::{
    jacobian_to_biguint, scalar_mul_glv, scalar_mul_wnaf, textbook_mul, Affine, ChainEval,
    Jacobian, TableMode,
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
    name = "endo-chain-bench",
    about = "GLV-2 via a frozen isogeny-chain endomorphism against width-w NAF, on any short Weierstrass curve"
)]
struct Args {
    /// The frozen chain set (JSON, format endosweep-chainsweep/1).
    #[arg(long)]
    chains: PathBuf,
    /// Refuse to run unless the chain set's sha256 is this (hex).
    #[arg(long)]
    expect_sha256: Option<String>,
    /// The chain sweep's scalar-multiplication model (JSON with
    /// `baseline` and `glv` expected counts), set beside the counts.
    #[arg(long)]
    model: Option<PathBuf>,
    /// Number M of random (scalar, point) pairs.
    #[arg(long, default_value_t = 256)]
    pairs: usize,
    /// Number R of timing rounds.
    #[arg(long, default_value_t = 20)]
    rounds: usize,
    /// RNG seed for the pairs.
    #[arg(long, default_value_t = 20261016)]
    seed: u64,
    /// wNAF widths of the baseline and the GLV grid.
    #[arg(long, value_delimiter = ',', default_value = "3,4,5,6,7")]
    widths: Vec<u32>,
    /// Chain whose evaluator × table × width grid is run (default: the
    /// least optimised model count in the set, ties by id).
    #[arg(long)]
    grid_chain: Option<String>,
    /// "EVAL/TABLE/W" at which every element's cheapest ordering is run.
    #[arg(long, default_value = "optimised/jacobian/5")]
    element_config: String,
    /// "TABLE/W" of the reference baseline (every ratio's numerator; timed
    /// again as the A/A control).
    #[arg(long, default_value = "jacobian/5")]
    reference: String,
    /// Chain id of an optional continuity arm.
    #[arg(long)]
    continuity_chain: Option<String>,
    /// "EVAL/TABLE/W" of the continuity arm.
    #[arg(long, default_value = "generic/affine/5")]
    continuity_config: String,
    /// Take the generator from the crate's curve zoo instead of the file
    /// (`gost-cryptopro-b`; for a chain set without `curve.G`).
    #[arg(long)]
    zoo_generator: Option<String>,
    /// Count field operations per arm instead of timing (needs the
    /// `endo-chain-opcount` feature).
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
    model_m_eq: Option<f64>,
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

fn parse_config(s: &str) -> (ChainEval, TableMode, u32) {
    let parts: Vec<&str> = s.split('/').collect();
    assert_eq!(parts.len(), 3, "configuration {s:?} is not EVAL/TABLE/W");
    (
        parse_eval(parts[0]),
        parse_table(parts[1]),
        parse_width(parts[2]),
    )
}

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

fn write_json(path: &Option<PathBuf>, record: &Value) {
    if let Some(path) = path {
        let text = serde_json::to_string_pretty(record).expect("serialisable record");
        std::fs::write(path, text + "\n")
            .unwrap_or_else(|e| panic!("write {}: {e}", path.display()));
        println!("json record written to {}", path.display());
    }
}

fn main() {
    let args = Args::parse();
    assert!(
        args.pairs > 0 && args.rounds > 0,
        "need at least one pair and one round"
    );
    if args.count_only && !cfg!(feature = "endo-chain-opcount") {
        eprintln!("--count-only needs a build with --features endo-chain-opcount");
        std::process::exit(2);
    }
    if !args.count_only && cfg!(feature = "endo-chain-opcount") {
        eprintln!("this build counts field operations (feature endo-chain-opcount); time with a default build");
        std::process::exit(2);
    }
    let bytes = std::fs::read(&args.chains)
        .unwrap_or_else(|e| panic!("read {}: {e}", args.chains.display()));
    let sha = hex::encode(sha256(&bytes));
    if let Some(want) = &args.expect_sha256 {
        if !want.eq_ignore_ascii_case(&sha) {
            eprintln!("chain set sha256 {sha} differs from the expected {want}; refusing to run");
            std::process::exit(2);
        }
    }
    let set = ChainSet::parse(std::str::from_utf8(&bytes).expect("UTF-8 chain set"))
        .unwrap_or_else(|e| panic!("chain set: {e}"));
    let fallback = match args.zoo_generator.as_deref() {
        None => None,
        Some("gost-cryptopro-b") => {
            let c = CurveParams::gost_cryptopro_b();
            assert!(
                c.p == set.p && c.a == set.a && c.b == set.b && c.n == set.n,
                "--zoo-generator gost-cryptopro-b on a chain set of another curve"
            );
            Some((c.gx, c.gy))
        }
        Some(other) => panic!("unknown --zoo-generator {other:?}"),
    };
    match set.limbs() {
        4 => run::<4>(&args, &set, &sha, fallback),
        6 => run::<6>(&args, &set, &sha, fallback),
        13 => run::<13>(&args, &set, &sha, fallback),
        l => {
            eprintln!("p needs {l} limbs; this binary is built for 4, 6 and 13");
            std::process::exit(2);
        }
    }
}

fn model_arms(path: &Option<PathBuf>) -> Option<Value> {
    path.as_ref().map(|p| {
        let text =
            std::fs::read_to_string(p).unwrap_or_else(|e| panic!("read {}: {e}", p.display()));
        serde_json::from_str(&text).expect("model JSON")
    })
}

fn model_lookup(
    model: &Option<Value>,
    section: &str,
    chain: Option<&str>,
    cfg: &str,
) -> Option<f64> {
    let m = model.as_ref()?;
    let s = &m[section];
    let v = match chain {
        Some(c) => &s[c][cfg],
        None => &s[cfg],
    };
    v.get("M_eq").and_then(Value::as_f64).or_else(|| v.as_f64())
}

fn run<const L: usize>(
    args: &Args,
    set: &ChainSet,
    chains_sha256: &str,
    fallback: Option<(BigUint, BigUint)>,
) {
    let (m, rounds) = (args.pairs, args.rounds);
    let PreparedSet {
        curve,
        generator: g,
        contexts: ctxs,
    } = set
        .prepare::<L>(fallback)
        .unwrap_or_else(|e| panic!("chain set: {e}"));
    let c = &curve;
    let model = model_arms(&args.model);
    let index = |id: &str| -> usize {
        set.chains
            .iter()
            .position(|c| c.id == id)
            .unwrap_or_else(|| panic!("chain {id:?} is not in the set"))
    };
    let grid_chain = match &args.grid_chain {
        Some(id) => index(id),
        None => (0..set.chains.len())
            .min_by_key(|&i| {
                (
                    model_m_eq(&set.chains[i], ChainEval::Optimised),
                    set.chains[i].id.clone(),
                )
            })
            .expect("non-empty set"),
    };
    let mut by_element: BTreeMap<(i64, i64), usize> = BTreeMap::new();
    for (i, ch) in set.chains.iter().enumerate() {
        let e = by_element.entry(ch.element).or_insert(i);
        let cur = &set.chains[*e];
        if (model_m_eq(ch, ChainEval::Optimised), &ch.id)
            < (model_m_eq(cur, ChainEval::Optimised), &cur.id)
        {
            *e = i;
        }
    }
    let mut element_chains: Vec<usize> = by_element.values().copied().collect();
    element_chains.sort_by_key(|&i| {
        (
            model_m_eq(&set.chains[i], ChainEval::Optimised),
            set.chains[i].id.clone(),
        )
    });
    let (el_eval, el_table, el_w) = parse_config(&args.element_config);
    let (ref_table, ref_w) = parse_baseline(&args.reference);
    let cont = args.continuity_chain.as_deref().map(|id| {
        let (e, t, w) = parse_config(&args.continuity_config);
        (index(id), e, t, w)
    });
    for &w in &args.widths {
        parse_width(&w.to_string());
    }

    // ── inputs: the pair generation of cryptopro-b-chain-sweep ──
    let mut rng = StdRng::seed_from_u64(args.seed);
    let mut points: Vec<Affine<L>> = Vec::with_capacity(m);
    let mut scalars: Vec<BigUint> = Vec::with_capacity(m);
    for _ in 0..m {
        let r = nonzero_below(&mut rng, &set.n);
        points.push(
            c.to_affine(&scalar_mul_wnaf(c, &g, &r, 5, TableMode::Affine))
                .expect("r·G with 0 < r < n is not the identity"),
        );
        scalars.push(nonzero_below(&mut rng, &set.n));
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
                model_m_eq: model_lookup(
                    &model,
                    "baseline",
                    None,
                    &format!("{}/w{w}", table_name(table)),
                ),
                correctness: String::new(),
                correct: false,
            });
        }
    }
    let glv_label = |i: usize, eval: ChainEval, table: TableMode, w: u32, what: &str| -> String {
        let ch = &set.chains[i];
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
    let mut push_glv =
        |i: usize, eval: ChainEval, table: TableMode, w: u32, what: &str, arms: &mut Vec<Arm>| {
            if glv_keys.contains(&(i, eval, table, w)) {
                return;
            }
            glv_keys.push((i, eval, table, w));
            let id = &set.chains[i].id;
            let cfg = format!("{}/{}/w{w}", eval_name(eval), table_name(table));
            arms.push(Arm {
                key: format!("glv/{id}/{cfg}"),
                label: glv_label(i, eval, table, w, what),
                class: "engineering",
                kind: Kind::Glv {
                    chain: i,
                    eval,
                    table,
                    w,
                },
                model_m_eq: model_lookup(&model, "glv", Some(id), &cfg),
                correctness: String::new(),
                correct: false,
            });
        };
    for &eval in &evals {
        for &table in &tables {
            for &w in &args.widths {
                push_glv(grid_chain, eval, table, w, " (grid)", &mut arms);
            }
        }
    }
    for &i in &element_chains {
        push_glv(i, el_eval, el_table, el_w, " (per element)", &mut arms);
    }
    if let Some((i, e, t, w)) = cont {
        push_glv(i, e, t, w, " (continuity)", &mut arms);
    }
    for (i, ch) in set.chains.iter().enumerate() {
        for &eval in &evals {
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
                kind: Kind::Phi { chain: i, eval },
                model_m_eq: Some(model_m_eq(ch, eval) as f64),
                correctness: String::new(),
                correct: false,
            });
        }
    }
    for &i in &element_chains {
        let ch = &set.chains[i];
        arms.push(Arm {
            key: format!("decompose/{}", ch.id),
            label: format!(
                "stage diagnostic: decomposition alone, basis of {} (num-bigint Babai rounding)",
                ch.element_label()
            ),
            class: "stage diagnostic",
            kind: Kind::Decompose { chain: i },
            model_m_eq: Some(0.0),
            correctness: String::new(),
            correct: false,
        });
    }
    arms.push(Arm {
        key: format!("aa/{}/w{ref_w}", table_name(ref_table)),
        label: "A/A control: the reference baseline timed again as a separate arm (identical code and inputs; its ratio to the reference is the noise floor)".to_string(),
        class: "A/A control",
        kind: Kind::AaControl { table: ref_table, w: ref_w },
        model_m_eq: model_lookup(&model, "baseline", None, &format!("{}/w{ref_w}", table_name(ref_table))),
        correctness: String::new(),
        correct: false,
    });

    // ── correctness, on every pair, before any timing ──
    let truth: Vec<Option<(BigUint, BigUint)>> = (0..m)
        .map(|i| textbook_mul(c, &points[i], &scalars[i]))
        .collect();
    let reference: Vec<Jacobian<L>> = (0..m)
        .map(|i| scalar_mul_wnaf(c, &points[i], &scalars[i], 5, TableMode::Affine))
        .collect();
    let reference_ok = (0..m).all(|i| jacobian_to_biguint(c, &reference[i]) == truth[i]);
    let lambda_p: Vec<Vec<Jacobian<L>>> = ctxs
        .iter()
        .map(|ctx| {
            (0..m)
                .map(|i| scalar_mul_wnaf(c, &points[i], &ctx.lambda, 5, TableMode::Affine))
                .collect()
        })
        .collect();
    let vectors_ok: Vec<bool> = set
        .chains
        .iter()
        .zip(ctxs.iter())
        .map(|(ch, ctx)| {
            ch.vectors.iter().all(|v| {
                let p = c.affine(&v.p.0, &v.p.1).expect("canonical");
                let (k1, k2) = ctx.decompose(&v.k);
                c.is_on_curve(&p)
                    && evals.iter().all(|&e| {
                        jacobian_to_biguint(c, &ctx.chain.apply_with(c, &p, e))
                            == Some(v.phi_p.clone())
                    })
                    && k1 == v.k1
                    && k2 == v.k2
                    && jacobian_to_biguint(
                        c,
                        &scalar_mul_glv(
                            c,
                            ctx,
                            &p,
                            &v.k,
                            5,
                            TableMode::Affine,
                            ChainEval::Optimised,
                        ),
                    ) == Some(v.k_p.clone())
                    && textbook_mul(c, &p, &v.k) == Some(v.k_p.clone())
            })
        })
        .collect();
    for arm in arms.iter_mut() {
        let (ok, text) = match arm.kind {
            Kind::Baseline { table, w } | Kind::AaControl { table, w } => {
                let ok = reference_ok
                    && (0..m).all(|i| {
                        jacobian_to_biguint(
                            c,
                            &scalar_mul_wnaf(c, &points[i], &scalars[i], w, table),
                        ) == truth[i]
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
                        c.eq_point(
                            &scalar_mul_glv(
                                c,
                                &ctxs[chain],
                                &points[i],
                                &scalars[i],
                                w,
                                table,
                                eval,
                            ),
                            &reference[i],
                        )
                    });
                (
                    ok,
                    format!("GLV == textbook k·P on all {m} pairs: {}", yes_no(ok)),
                )
            }
            Kind::Phi { chain, eval } => {
                let ok = vectors_ok[chain]
                    && (0..m).all(|i| {
                        c.eq_point(
                            &ctxs[chain].chain.apply_with(c, &points[i], eval),
                            &lambda_p[chain][i],
                        )
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
            black_box(scalar_mul_wnaf(c, &points[i], &scalars[i], w, table));
        }
        Kind::Glv {
            chain,
            eval,
            table,
            w,
        } => {
            black_box(scalar_mul_glv(
                c,
                &ctxs[chain],
                &points[i],
                &scalars[i],
                w,
                table,
                eval,
            ));
        }
        Kind::Phi { chain, eval } => {
            black_box(ctxs[chain].chain.apply_with(c, &points[i], eval));
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
        "opcount_feature": cfg!(feature = "endo-chain-opcount"),
        "limbs": L,
    });
    let curve_json = json!({
        "target": set.target,
        "p": format!("0x{:x}", set.p),
        "a": format!("0x{:x}", set.a),
        "b": format!("0x{:x}", set.b),
        "n": format!("0x{:x}", set.n),
        "p_bits": set.p.bits(),
        "n_bits": set.n.bits(),
        "cofactor": set.cofactor.to_string(),
        "D_K": set.d_k,
        "class_number": set.class_number,
        "doubling": c.dbl.describe(),
        "inversion_M_S": c.f.inv_cost(),
        "generator": [format!("0x{:x}", c.f.to_biguint(&g.x)), format!("0x{:x}", c.f.to_biguint(&g.y))],
        "generator_source": if set.generator.is_some() { set.generator_source.clone() } else if args.zoo_generator.is_some() { format!("crate curve zoo ({})", args.zoo_generator.clone().unwrap_or_default()) } else { "first test vector's P".into() },
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
        "element_chains": element_chains.iter().map(|&i| set.chains[i].id.clone()).collect::<Vec<_>>(),
        "reference": args.reference,
        "continuity_chain": args.continuity_chain,
        "continuity_config": if args.continuity_chain.is_some() { json!(args.continuity_config) } else { Value::Null },
        "model_path": args.model.as_ref().map(|p| p.display().to_string()),
    });

    // ── counting mode ──
    if args.count_only {
        let counted = count_arms(&arms, m, &run);
        let rows: Vec<Value> = arms
            .iter()
            .zip(counted.iter())
            .map(|(arm, k)| {
                json!({
                    "key": arm.key,
                    "mul": k.mul, "sqr": k.sqr, "inv": k.inv, "add": k.add,
                    "m_eq": k.m_eq(), "m_eq_min": k.m_eq_min, "m_eq_max": k.m_eq_max,
                    "model_m_eq": arm.model_m_eq,
                    "correct": arm.correct,
                })
            })
            .collect();
        let record = json!({
            "bench": "endo-chain-bench",
            "mode": "count",
            "unit": "field operations per scalar multiplication (per call for stage diagnostics), mean over all pairs",
            "curve": curve_json,
            "inputs": inputs,
            "host": host,
            "build": build,
            "all_correct": all_correct,
            "arms": rows,
            "timestamp_unix": timestamp,
        });
        println!("# endo-chain-bench (counting run): {}", set.target);
        println!();
        println!(
            "curve: {} ({} -bit p, {}-bit n, cofactor {}, D_K = {}); doubling {}; inversion {} M + {} S",
            set.target,
            set.p.bits(),
            set.n.bits(),
            set.cofactor,
            set.d_k,
            c.dbl.describe(),
            c.f.inv_cost().0,
            c.f.inv_cost().1
        );
        println!(
            "chain set: {} (sha256 {chains_sha256}); M = {m} pairs, seed {}",
            args.chains.display(),
            args.seed
        );
        println!();
        println!("| arm | M | S | inversions | M_eq = M + S (mean) | min | max | model M_eq | counted / model − 1 | correct |");
        println!("|:--|--:|--:|--:|--:|--:|--:|--:|--:|:--|");
        for (arm, k) in arms.iter().zip(counted.iter()) {
            let (mm, dev) = match arm.model_m_eq {
                Some(x) if x > 0.0 => (
                    format!("{x:.2}"),
                    format!("{:+.3} %", 100.0 * (k.m_eq() / x - 1.0)),
                ),
                Some(x) => (format!("{x:.0}"), "—".into()),
                None => ("—".into(), "—".into()),
            };
            println!(
                "| {} | {:.2} | {:.2} | {:.3} | {:.2} | {:.0} | {:.0} | {mm} | {dev} | {} |",
                arm.key,
                k.mul,
                k.sqr,
                k.inv,
                k.m_eq(),
                k.m_eq_min,
                k.m_eq_max,
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
    if let Some(k) = &counts {
        for arm in &arms {
            assert!(
                k.contains_key(&arm.key),
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
    let ref_count = counts.as_ref().map(|k| k[&arms[ref_idx].key].m_eq());
    let is_total = |k: Kind| !matches!(k, Kind::Phi { .. } | Kind::Decompose { .. });

    println!("# endo-chain-bench: {}", set.target);
    println!();
    println!(
        "curve: {} ({}-bit p, {}-bit n, cofactor {}, D_K = {}, class number {}); arithmetic: variable-time Jacobian, {}, add-2007-bl, madd-2007-bl, Fermat inversion {} M + {} S, run-time Montgomery field of {L} limbs, shared by every arm",
        set.target,
        set.p.bits(),
        set.n.bits(),
        set.cofactor,
        set.d_k,
        set.class_number.map(|h| h.to_string()).unwrap_or_else(|| "?".into()),
        c.dbl.describe(),
        c.f.inv_cost().0,
        c.f.inv_cost().1
    );
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
        let k = counts.as_ref().map(|k| k[&arm.key]);
        let cnum = |f: fn(&Counted) -> f64, digits: usize| match k {
            Some(k) => format!("{:.*}", digits, f(&k)),
            None => "—".into(),
        };
        let total = is_total(arm.kind);
        let rc = match (k, ref_count, total) {
            (Some(k), Some(r), true) => format!("{:.3}", r / k.m_eq()),
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
            cnum(|k| k.m_eq(), 1),
            cnum(|k| k.mul, 1),
            cnum(|k| k.sqr, 1),
            cnum(|k| k.inv, 0),
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
            let k = counts.as_ref().map(|k| k[&arm.key]);
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
                    "counted": match (k, ref_count) { (Some(k), Some(r)) => json!(r / k.m_eq()), _ => Value::Null },
                }) } else { Value::Null },
                "counted": k.map(|k| json!({
                    "mul": k.mul, "sqr": k.sqr, "inv": k.inv, "add": k.add,
                    "m_eq": k.m_eq(), "m_eq_min": k.m_eq_min, "m_eq_max": k.m_eq_max,
                })),
                "model_m_eq": arm.model_m_eq,
                "correctness": arm.correctness,
                "correct": arm.correct,
            })
        })
        .collect();
    let record = json!({
        "bench": "endo-chain-bench",
        "mode": "timing",
        "unit": "nanoseconds per scalar multiplication (per call for stage diagnostics); counted field operations per op from the counting run",
        "curve": curve_json,
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
        "chains": set.chains.iter().map(|ch| json!({
            "id": ch.id,
            "element": [ch.element.0, ch.element.1],
            "norm": ch.norm,
            "order": ch.order,
            "babai_bound_bits": ch.babai_bound_bits,
            "model_m_eq": {"generic": model_m_eq(ch, ChainEval::Generic), "optimised": model_m_eq(ch, ChainEval::Optimised)},
        })).collect::<Vec<Value>>(),
        "correctness": {
            "reference_wnaf_equals_textbook_on_all_pairs": reference_ok,
            "test_vectors_per_chain": set.chains.iter().zip(vectors_ok.iter()).map(|(ch, ok)| json!({"chain": ch.id, "ok": ok})).collect::<Vec<Value>>(),
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

/// Per arm: run every pair once with the counters reset before each pair.
#[cfg(feature = "endo-chain-opcount")]
fn count_arms(arms: &[Arm], m: usize, run: &dyn Fn(Kind, usize)) -> Vec<Counted> {
    use crypto_lib::ecc::endo_generic_field::opcount;
    arms.iter()
        .map(|arm| {
            let mut k = Counted {
                m_eq_min: f64::INFINITY,
                m_eq_max: 0.0,
                ..Counted::default()
            };
            for i in 0..m {
                opcount::reset();
                run(arm.kind, i);
                let o = opcount::read();
                k.mul += o.mul as f64;
                k.sqr += o.sqr as f64;
                k.inv += o.inv as f64;
                k.add += o.add as f64;
                let e = (o.mul + o.sqr) as f64;
                k.m_eq_min = k.m_eq_min.min(e);
                k.m_eq_max = k.m_eq_max.max(e);
            }
            let mf = m as f64;
            k.mul /= mf;
            k.sqr /= mf;
            k.inv /= mf;
            k.add /= mf;
            k
        })
        .collect()
}

#[cfg(not(feature = "endo-chain-opcount"))]
fn count_arms(_arms: &[Arm], _m: usize, _run: &dyn Fn(Kind, usize)) -> Vec<Counted> {
    unreachable!("--count-only is refused without the endo-chain-opcount feature")
}

//! `ecdlp-nist` — the ECDLP across every NIST curve.
//!
//! ```text
//! ecdlp-nist list                         # the fifteen curves with their audit
//! ecdlp-nist audit K-233                  # one curve's structural checks
//! ecdlp-nist solve P-256 --secret 123456789 --interval-bits 32
//! ecdlp-nist solve B-571 --random-bits 24 --method kangaroo --threads 4
//! ecdlp-nist solve toy-k23a1 --random-bits 21          # whole-group rho, Frobenius folded
//! ecdlp-nist solve toy-k23a1 --random-bits 21 --no-fold
//! ecdlp-nist smart --bits 256 --seed 7                 # construct an anomalous curve and break it
//! ecdlp-nist bench P-521                               # additions per second on this host
//! ```
//!
//! Every subcommand takes `--json` for a machine-readable report.  See
//! `docs/ECDLP_NIST.md`.

use std::process::ExitCode;
use std::time::Instant;

use clap::{Args, Parser, Subcommand};
use num_bigint::{BigUint, RandBigInt};
use rand::rngs::StdRng;
use rand::SeedableRng;
use serde_json::json;

use crypto_lib::cryptanalysis::ecdlp_nist::{
    self as nist, anomalous, audit, curve_by_name, curve_names, nist_curves, solve, Audit, Curve,
    EcdlpGroup, Group, Method, SolveOptions,
};

#[derive(Parser)]
#[command(
    name = "ecdlp-nist",
    version,
    about = "ECDLP solvers and structural audit across the fifteen NIST curves",
    long_about = "Solve Q = [k]G on P-192…P-521, K-163…K-571 and B-163…B-571 by the \
                  branch the curve calls for (Smart on anomalous prime curves, \
                  Pohlig–Hellman on composite orders, BSGS/kangaroo on an interval, \
                  folded parallel rho otherwise), and audit each curve for the \
                  structural conditions those attacks need.  Bounded instances run \
                  on the real curves; whole-group searches are refused as infeasible \
                  unless forced, and no result here is a statement about deployed keys."
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
    /// Emit a machine-readable JSON report.
    #[arg(long, global = true)]
    json: bool,
}

#[derive(Subcommand)]
enum Command {
    /// List the NIST curves with their audit summary.
    List,
    /// Structural audit of one curve.
    Audit(CurveArg),
    /// Plant (or accept) a target point and solve for its scalar.
    Solve(SolveArgs),
    /// Construct an anomalous curve by CM and break it with Smart's attack.
    Smart(SmartArgs),
    /// Measure group additions per second on a curve.
    Bench(BenchArgs),
}

#[derive(Args)]
struct CurveArg {
    /// Curve name: P-256, K-163, B-571, secp256r1, toy-k23a1, anomalous-128-3, …
    curve: String,
}

#[derive(Args)]
struct SolveArgs {
    /// Curve name (see `list`).
    curve: String,
    /// Plant Q = [k]G with this scalar (decimal or 0x-hex).
    #[arg(long, conflicts_with_all = ["random_bits", "point"])]
    secret: Option<String>,
    /// Plant Q = [k]G with a random k of this many bits (seeded).
    #[arg(long, conflicts_with = "point")]
    random_bits: Option<u32>,
    /// Solve an externally supplied point `x,y` (hex).  Needs an interval
    /// unless the curve is a toy.
    #[arg(long, value_name = "X,Y")]
    point: Option<String>,
    /// Interval lower bound `lo` (decimal or 0x-hex); default 0.
    #[arg(long)]
    interval_lo: Option<String>,
    /// Interval width as a bit size: k ∈ [lo, lo + 2^bits).
    #[arg(long)]
    interval_bits: Option<u32>,
    /// Interval width as a number (decimal or 0x-hex).
    #[arg(long, conflicts_with = "interval_bits")]
    interval_width: Option<String>,
    /// auto, smart, pohlig-hellman, bsgs, kangaroo, rho.
    #[arg(long, default_value = "auto")]
    method: String,
    #[arg(long, default_value_t = 1)]
    threads: usize,
    #[arg(long, default_value_t = 0x5EED)]
    seed: u64,
    /// Iteration budget for rho / kangaroo.
    #[arg(long)]
    max_iterations: Option<u64>,
    /// Distinguished-point bits (default: derived from the expected cost).
    #[arg(long)]
    dp_bits: Option<u32>,
    /// Disable negation / Frobenius folding in the rho (control arm).
    #[arg(long)]
    no_fold: bool,
    /// Start a whole-group rho even when the audit calls it infeasible.
    #[arg(long)]
    force: bool,
}

#[derive(Args)]
struct SmartArgs {
    /// Size of the anomalous prime.
    #[arg(long, default_value_t = 128)]
    bits: u32,
    /// Selects the CM discriminant and the prime.
    #[arg(long, default_value_t = 0)]
    seed: u64,
    /// The planted scalar (default: random, seeded).
    #[arg(long)]
    secret: Option<String>,
}

#[derive(Args)]
struct BenchArgs {
    curve: String,
    /// Additions to time.
    #[arg(long, default_value_t = 2000)]
    ops: u64,
}

fn parse_uint(s: &str) -> Result<BigUint, String> {
    let t = s.trim();
    if let Some(h) = t.strip_prefix("0x").or_else(|| t.strip_prefix("0X")) {
        BigUint::parse_bytes(h.as_bytes(), 16)
    } else {
        BigUint::parse_bytes(t.as_bytes(), 10)
    }
    .ok_or_else(|| format!("not a number: {s:?}"))
}

fn audit_row(a: &Audit) -> String {
    format!(
        "{:<6} {:<8} {:>5} {:>5} {:>3}  {:<5} {:<5} {:<5}  {:>6}  {:>7.1}  {:>7.1}  {}",
        a.name,
        a.family.tag(),
        a.field_bits,
        a.order_bits,
        a.cofactor,
        if a.anomalous { "YES" } else { "no" },
        if a.order_is_prime { "yes" } else { "NO" },
        a.embedding_degree
            .map(|k| k.to_string())
            .unwrap_or_else(|| format!(">{}", a.embedding_degree_bound)),
        a.rho_class_size,
        a.rho_log2_expected_unfolded,
        a.rho_log2_expected,
        a.recommended
    )
}

fn audit_header() -> String {
    format!(
        "{:<6} {:<8} {:>5} {:>5} {:>3}  {:<5} {:<5} {:<5}  {:>6}  {:>7}  {:>7}  {}",
        "curve",
        "family",
        "q",
        "n",
        "h",
        "anom",
        "prime",
        "embed",
        "class",
        "log2ρ",
        "folded",
        "method"
    )
}

fn cmd_list(json: bool) -> Result<(), String> {
    let curves = nist_curves();
    let audits: Vec<Audit> = curves.iter().map(audit).collect();
    if json {
        let rows: Vec<serde_json::Value> = curves
            .iter()
            .zip(&audits)
            .map(|(c, a)| json!({"curve": c.describe(), "audit": a}))
            .collect();
        println!(
            "{}",
            serde_json::to_string_pretty(&json!({
                "operation": "list",
                "count": rows.len(),
                "curves": rows,
                "other_names": curve_names(),
            }))
            .unwrap()
        );
        return Ok(());
    }
    println!("{}", audit_header());
    for a in &audits {
        println!("{}", audit_row(a));
    }
    println!();
    println!("q, n: bits of the field and of the generator's order; h: cofactor; anom: #E(F_q) = q (Smart's attack);");
    println!("prime: generator order prime (else Pohlig–Hellman); embed: embedding degree if ≤ {} (MOV);", nist::EMBEDDING_DEGREE_BOUND);
    println!(
        "class: automorphisms folded by the rho (2 = negation, 2m = negation × Frobenius on K-m);"
    );
    println!(
        "log2ρ: log2 of expected rho additions unfolded √(πn/2); folded: with the class folded."
    );
    println!(
        "Also accepted: {}",
        curve_names()
            .into_iter()
            .filter(|n| n.starts_with("toy") || n.starts_with("anom"))
            .collect::<Vec<_>>()
            .join(", ")
    );
    Ok(())
}

fn cmd_audit(name: &str, json: bool) -> Result<(), String> {
    let c = curve_by_name(name)?;
    let a = audit(&c);
    if json {
        println!(
            "{}",
            serde_json::to_string_pretty(
                &json!({"operation": "audit", "curve": c.describe(), "audit": a})
            )
            .unwrap()
        );
        return Ok(());
    }
    println!("{}", audit_header());
    println!("{}", audit_row(&a));
    println!();
    println!("standard:      {}", a.standard);
    println!("trace:         {}", a.trace);
    println!("group order:   {}", c.group_order);
    println!("subgroup n:    {}", c.order());
    if let Some(l) = &a.frobenius_eigenvalue {
        println!("Frobenius λ:   {l}");
    }
    if let Some(p) = a.extension_degree_prime {
        println!("m prime:       {p} (composite m would open Weil descent)");
    }
    if let Some(k) = a.koblitz_order_consistent {
        println!(
            "Lucas order:   {}",
            if k { "matches h·n" } else { "MISMATCH" }
        );
    }
    if !a.order_factors.is_empty() {
        println!("order factors: {:?}", a.order_factors);
    }
    for n in &a.notes {
        println!("note:          {n}");
    }
    Ok(())
}

fn cmd_solve(args: &SolveArgs, json: bool) -> Result<(), String> {
    let curve = curve_by_name(&args.curve)?;
    let method: Method = args.method.parse()?;
    let mut rng = StdRng::seed_from_u64(args.seed ^ 0xC0FFEE);
    let lo = match &args.interval_lo {
        Some(s) => parse_uint(s)?,
        None => BigUint::from(0u32),
    };
    let width = match (&args.interval_bits, &args.interval_width) {
        (Some(b), _) => Some(BigUint::from(1u32) << *b),
        (None, Some(w)) => Some(parse_uint(w)?),
        (None, None) => None,
    };
    // Target: planted or supplied.
    let (target, planted): (nist::AnyPoint, Option<BigUint>) = if let Some(pt) = &args.point {
        let (x, y) = pt
            .split_once(',')
            .ok_or_else(|| "--point wants X,Y in hex".to_string())?;
        (curve.parse_point(x, y)?, None)
    } else {
        let k = match (&args.secret, args.random_bits) {
            (Some(s), _) => parse_uint(s)? % curve.order(),
            (None, Some(bits)) => {
                let w = width.clone().unwrap_or_else(|| BigUint::from(1u32) << bits);
                // Random k in [lo, lo + min(w, 2^bits)) ∩ [0, n).
                let span = w.min(BigUint::from(1u32) << bits);
                (&lo + rng.gen_biguint_below(&span)) % curve.order()
            }
            (None, None) => {
                return Err("give --secret K, --random-bits B, or --point X,Y".into());
            }
        };
        (curve.plant(&k), Some(k))
    };
    let interval = width.map(|w| (lo.clone(), w));
    let opts = SolveOptions {
        method,
        interval,
        threads: args.threads.max(1),
        seed: args.seed,
        max_iterations: args.max_iterations,
        dp_bits: args.dp_bits,
        fold_automorphisms: !args.no_fold,
        force: args.force,
        ..SolveOptions::default()
    };
    let report = solve(&curve, &target, &opts);
    let (tx, ty) = curve.point_hex(&target).unwrap_or_default();
    let planted_ok = planted
        .as_ref()
        .map(|k| report.scalar_value.as_ref() == Some(k));
    if json {
        println!(
            "{}",
            serde_json::to_string_pretty(&json!({
                "operation": "solve",
                "target": {"x": tx, "y": ty},
                "planted_scalar": planted.as_ref().map(|k| k.to_string()),
                "planted_recovered": planted_ok,
                "report": report,
            }))
            .unwrap()
        );
    } else {
        println!("curve:     {} ({})", report.curve, report.family.tag());
        println!("target:    ({tx}, {ty})");
        println!("method:    {}", report.method);
        match &report.scalar {
            Some(k) => println!(
                "scalar:    {k}  (0x{})  verified={}",
                report.scalar_hex.as_deref().unwrap_or(""),
                report.verified
            ),
            None => println!("scalar:    not found"),
        }
        if let Some(k) = &planted {
            println!("planted:   {k}  recovered={}", planted_ok.unwrap_or(false));
        }
        println!(
            "cost:      {} group additions in {} ms on {} thread(s); expected ≈ {:.3e}",
            report.group_ops, report.elapsed_ms, report.threads, report.expected_ops
        );
        if report.table_size > 0 {
            println!("table:     {} entries", report.table_size);
        }
        for n in &report.notes {
            println!("note:      {n}");
        }
        if let Some(f) = &report.failure {
            println!("failure:   {f}");
        }
    }
    if report.failure.is_some() && report.scalar.is_none() {
        return Err(String::new());
    }
    Ok(())
}

fn cmd_smart(args: &SmartArgs, json: bool) -> Result<(), String> {
    let t0 = Instant::now();
    let cp = anomalous::generate(args.bits, args.seed)?;
    let construct_ms = t0.elapsed().as_millis();
    let curve = curve_by_name(&format!("anomalous-{}-{}", args.bits, args.seed))?;
    let mut rng = StdRng::seed_from_u64(args.seed ^ 0x5A4D);
    let k = match &args.secret {
        Some(s) => parse_uint(s)? % curve.order(),
        None => rng.gen_biguint_below(curve.order()),
    };
    let target = curve.plant(&k);
    let report = solve(&curve, &target, &SolveOptions::default());
    let ok = report.scalar_value.as_ref() == Some(&k);
    if json {
        println!(
            "{}",
            serde_json::to_string_pretty(&json!({
                "operation": "smart",
                "curve": {
                    "name": cp.name, "p": cp.p.to_string(), "a": cp.a.to_string(), "b": cp.b.to_string(),
                    "gx": cp.gx.to_string(), "gy": cp.gy.to_string(), "order": cp.n.to_string(),
                    "bits": cp.p.bits(), "construct_ms": construct_ms,
                },
                "planted_scalar": k.to_string(),
                "planted_recovered": ok,
                "report": report,
            }))
            .unwrap()
        );
    } else {
        println!(
            "anomalous curve ({} bits, {}): y² = x³ + {}x + {} over F_p",
            cp.p.bits(),
            cp.name,
            cp.a,
            cp.b
        );
        println!("p = #E(F_p) = {}", cp.p);
        println!("constructed in {construct_ms} ms; planted k = {k}");
        println!("method:  {}", report.method);
        match &report.scalar {
            Some(s) => println!(
                "recovered k = {s}  (match={ok}) in {} ms",
                report.elapsed_ms
            ),
            None => println!(
                "not recovered: {}",
                report.failure.as_deref().unwrap_or("?")
            ),
        }
        for n in &report.notes {
            println!("note:    {n}");
        }
    }
    if ok {
        Ok(())
    } else {
        Err(String::new())
    }
}

fn bench_group<G: EcdlpGroup>(g: &G, ops: u64) -> (f64, u64) {
    let gen = g.generator();
    let mut p = g.double(&gen);
    let t0 = Instant::now();
    for _ in 0..ops {
        p = g.add(&p, &gen);
    }
    let dt = t0.elapsed().as_secs_f64();
    (ops as f64 / dt, g.key(&p))
}

fn cmd_bench(args: &BenchArgs, json: bool) -> Result<(), String> {
    let c: Curve = curve_by_name(&args.curve)?;
    let (rate, _) = match &c.group {
        Group::Prime(g) => bench_group(g.as_ref(), args.ops),
        Group::Binary(g) => bench_group(g.as_ref(), args.ops),
    };
    let a = audit(&c);
    let folded_secs = 2f64.powf(a.rho_log2_expected) / rate;
    if json {
        println!(
            "{}",
            serde_json::to_string_pretty(&json!({
                "operation": "bench",
                "curve": c.name,
                "additions_per_second": rate,
                "rho_log2_expected_folded": a.rho_log2_expected,
                "single_thread_rho_seconds_log2": folded_secs.log2(),
            }))
            .unwrap()
        );
    } else {
        println!(
            "{}: {:.0} affine additions/s (this host, one thread, BigUint/F2m arithmetic)",
            c.name, rate
        );
        println!(
            "whole-group rho would need 2^{:.1} additions ≈ 2^{:.1} s single-threaded here (≈ 2^{:.1} years)",
            a.rho_log2_expected,
            folded_secs.log2(),
            (folded_secs / 3.15e7).log2()
        );
    }
    Ok(())
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let r = match &cli.command {
        Command::List => cmd_list(cli.json),
        Command::Audit(a) => cmd_audit(&a.curve, cli.json),
        Command::Solve(a) => cmd_solve(a, cli.json),
        Command::Smart(a) => cmd_smart(a, cli.json),
        Command::Bench(a) => cmd_bench(a, cli.json),
    };
    match r {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            if !e.is_empty() {
                eprintln!("error: {e}");
            }
            ExitCode::FAILURE
        }
    }
}

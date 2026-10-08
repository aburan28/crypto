//! Native checks behind the ranked P-224 hidden-weakness hypotheses.
//!
//! Writes one JSON file per check into `--out`.  Protocol and report:
//! `research/p224_hypothesis_ranking_20261008/`.

use clap::Parser;
use crypto_lib::cryptanalysis::b_seed_profile::{profile_on_primes, profile_primes};
use crypto_lib::cryptanalysis::p224_hypotheses as h;
use crypto_lib::ecc::curve::CurveParams;
use num_bigint::{BigInt, BigUint};
use rayon::prelude::*;
use serde::Serialize;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode};
use std::time::Instant;

#[derive(Parser, Serialize)]
#[command(
    about = "Seed, p±1, CM-discriminant, twist, embedding-degree and b-residue checks for P-224"
)]
struct Cli {
    /// Output directory.
    #[arg(long)]
    out: PathBuf,
    /// Trial-division bound before Pollard–Brent.
    #[arg(long, default_value_t = 1 << 20)]
    trial_bound: u32,
    /// Pollard–Brent iterations allowed per factored integer.
    #[arg(long, default_value_t = 1 << 26)]
    rho_budget: u64,
    /// Prime bounds Q for the b-residue profile of P-224 and P-256.
    #[arg(long, value_delimiter = ',', default_values_t = [100_000u64, 1_000_000, 10_000_000, 100_000_000])]
    bounds: Vec<u64>,
    /// Prime bound for the null calibration panel.
    #[arg(long, default_value_t = 1_000_000)]
    null_bound: u64,
    /// Number of null coefficients.
    #[arg(long, default_value_t = 1000)]
    null_draws: u32,
    /// Domain label hashed into every null coefficient.
    #[arg(long, default_value = "p224-b-null/v1")]
    null_label: String,
}

#[derive(Serialize)]
struct ResidueRow {
    curve: String,
    prime_bound: u64,
    n_primes: usize,
    max_prime: u64,
    ks_statistic: f64,
    ks_pvalue: f64,
    decile_counts: [u64; 10],
    /// `(1 + #{null D ≥ D}) / (draws + 1)`, at the null bound only.
    empirical_pvalue: Option<f64>,
}

#[derive(Serialize)]
struct NullPanel {
    label: String,
    field: String,
    prime_bound: u64,
    draws: u32,
    ks_statistics: Vec<f64>,
    ks_pvalues: Vec<f64>,
    below_0_001: usize,
    below_0_01: usize,
    below_0_05: usize,
    mean_pvalue: f64,
}

#[derive(Serialize)]
struct PublishedCheck {
    curve: String,
    source: String,
    quantity: String,
    factors: Vec<String>,
    /// `null` if the product or a primality test failed.
    weakest_certainty: Option<h::Certainty>,
}

#[derive(Serialize)]
struct RunRecord {
    args: Cli,
    git_revision: String,
    git_dirty: bool,
    rustc: String,
    threads: usize,
    seconds: Vec<(String, f64)>,
}

fn write_json<T: Serialize>(dir: &Path, name: &str, value: &T) -> Result<(), String> {
    let path = dir.join(name);
    let mut bytes = serde_json::to_vec_pretty(value).map_err(|e| e.to_string())?;
    bytes.push(b'\n');
    fs::write(&path, bytes).map_err(|e| format!("write {}: {e}", path.display()))
}

fn command_output(program: &str, args: &[&str]) -> String {
    Command::new(program)
        .args(args)
        .output()
        .ok()
        .filter(|o| o.status.success())
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_else(|| "unknown".to_string())
}

/// SafeCurves' factorisations ("CM field discriminants" and "twist
/// security" pages), checked here independently of how they were found.
fn published_checks() -> Vec<PublishedCheck> {
    let p224 = CurveParams::p224();
    let p256 = CurveParams::p256();
    let abs_delta = |c: &CurveParams| {
        let t = h::trace(c);
        (BigInt::from(&c.p * 4u8) - &t * &t).magnitude().clone()
    };
    let twist = |c: &CurveParams| &c.p + 1u8 + &c.p + 1u8 - &c.n;
    let rows: [(&str, &str, BigUint, Vec<&str>); 3] = [
        (
            "P-224",
            "4p - t^2 (= 9 |D_K|)",
            abs_delta(&p224),
            vec![
                "3",
                "3",
                "3",
                "29",
                "79",
                "7523",
                "40927",
                "11549194661",
                "388425074903852603481408727235725774844022299",
            ],
        ),
        (
            "P-256",
            "4p - t^2 (= |D_K|)",
            abs_delta(&p256),
            vec![
                "3",
                "5",
                "456597257999",
                "1428624589419343516204097",
                "46523541035814968339936406074986559003387",
            ],
        ),
        (
            "P-224",
            "twist order p + 1 + t",
            twist(&p224),
            vec![
                "3",
                "3",
                "11",
                "47",
                "3015283",
                "40375823",
                "267983539294927",
                "177594041488131583478651368420021457",
            ],
        ),
    ];
    rows.into_iter()
        .map(|(curve, quantity, n, factors)| PublishedCheck {
            curve: curve.to_string(),
            source: "https://safecurves.cr.yp.to/ (disc.html, twist.html)".to_string(),
            quantity: quantity.to_string(),
            weakest_certainty: h::verify_factorization(&n, &factors),
            factors: factors.into_iter().map(str::to_string).collect(),
        })
        .collect()
}

fn run(cli: Cli) -> Result<(), String> {
    fs::create_dir_all(&cli.out).map_err(|e| format!("create {}: {e}", cli.out.display()))?;
    let out = cli.out.clone();
    let mut seconds = Vec::new();
    let all = Instant::now();

    let started = Instant::now();
    let seeds: Vec<h::SeedCheck> = h::nist_seeded_curves().iter().map(h::check_seed).collect();
    write_json(&out, "seeds.json", &seeds)?;
    seconds.push(("seeds".to_string(), started.elapsed().as_secs_f64()));

    let started = Instant::now();
    let mut fields: Vec<CurveParams> = h::nist_seeded_curves()
        .into_iter()
        .map(|c| c.params)
        .collect();
    fields.push(CurveParams::brainpool_p224r1());
    let sides: Vec<h::PrimeSide> = fields
        .par_iter()
        .flat_map(|c| h::prime_sides(c.name, &c.p, cli.trial_bound, cli.rho_budget))
        .collect();
    write_json(&out, "prime_sides.json", &sides)?;
    seconds.push(("prime_sides".to_string(), started.elapsed().as_secs_f64()));

    let curves = [CurveParams::p224(), CurveParams::p256()];
    let started = Instant::now();
    let cm: Vec<h::CmDiscriminant> = curves
        .par_iter()
        .map(|c| h::cm_discriminant(c, cli.trial_bound, cli.rho_budget))
        .collect();
    write_json(&out, "cm_discriminant.json", &cm)?;
    seconds.push((
        "cm_discriminant".to_string(),
        started.elapsed().as_secs_f64(),
    ));

    let started = Instant::now();
    let twists: Vec<h::TwistOrder> = curves
        .par_iter()
        .map(|c| h::twist_order(c, cli.trial_bound, cli.rho_budget))
        .collect();
    write_json(&out, "twist.json", &twists)?;
    seconds.push(("twist".to_string(), started.elapsed().as_secs_f64()));

    let started = Instant::now();
    let embedding: Vec<h::EmbeddingDegree> = curves
        .par_iter()
        .map(|c| h::embedding_degree(c, cli.trial_bound, cli.rho_budget))
        .collect();
    write_json(&out, "embedding_degree.json", &embedding)?;
    seconds.push((
        "embedding_degree".to_string(),
        started.elapsed().as_secs_f64(),
    ));

    let started = Instant::now();
    write_json(&out, "published_factorizations.json", &published_checks())?;
    seconds.push((
        "published_factorizations".to_string(),
        started.elapsed().as_secs_f64(),
    ));

    let started = Instant::now();
    let null_primes = profile_primes(cli.null_bound);
    let null_b = h::null_coefficients(&cli.null_label, &CurveParams::p224().p, cli.null_draws);
    let null: Vec<(f64, f64)> = null_b
        .par_iter()
        .map(|b| {
            let r = profile_on_primes(b, &null_primes);
            (r.ks_statistic, r.ks_pvalue)
        })
        .collect();
    let pvalues: Vec<f64> = null.iter().map(|x| x.1).collect();
    let panel = NullPanel {
        label: cli.null_label.clone(),
        field: "P-224 prime".to_string(),
        prime_bound: cli.null_bound,
        draws: cli.null_draws,
        ks_statistics: null.iter().map(|x| x.0).collect(),
        below_0_001: pvalues.iter().filter(|&&p| p < 0.001).count(),
        below_0_01: pvalues.iter().filter(|&&p| p < 0.01).count(),
        below_0_05: pvalues.iter().filter(|&&p| p < 0.05).count(),
        mean_pvalue: pvalues.iter().sum::<f64>() / pvalues.len().max(1) as f64,
        ks_pvalues: pvalues,
    };
    seconds.push(("b_null_panel".to_string(), started.elapsed().as_secs_f64()));

    let started = Instant::now();
    let mut rows = Vec::new();
    for c in &curves {
        for &bound in &cli.bounds {
            let primes = profile_primes(bound);
            let r = profile_on_primes(&c.b, &primes);
            let empirical = (bound == cli.null_bound).then(|| {
                let above = panel
                    .ks_statistics
                    .iter()
                    .filter(|&&d| d >= r.ks_statistic)
                    .count();
                (1 + above) as f64 / (panel.draws + 1) as f64
            });
            rows.push(ResidueRow {
                curve: c.name.to_string(),
                prime_bound: bound,
                n_primes: r.n_primes,
                max_prime: r.max_prime,
                ks_statistic: r.ks_statistic,
                ks_pvalue: r.ks_pvalue,
                decile_counts: r.decile_counts,
                empirical_pvalue: empirical,
            });
        }
    }
    write_json(&out, "b_residues.json", &rows)?;
    write_json(&out, "b_null_panel.json", &panel)?;
    seconds.push(("b_residues".to_string(), started.elapsed().as_secs_f64()));
    seconds.push(("total".to_string(), all.elapsed().as_secs_f64()));

    let record = RunRecord {
        git_revision: command_output("git", &["rev-parse", "HEAD"]),
        git_dirty: !command_output("git", &["status", "--porcelain", "--untracked-files=no"])
            .is_empty(),
        rustc: command_output("rustc", &["-V"]),
        threads: rayon::current_num_threads(),
        seconds,
        args: cli,
    };
    write_json(&out, "run.json", &record)?;
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("p224_hypotheses: {e}");
            ExitCode::FAILURE
        }
    }
}

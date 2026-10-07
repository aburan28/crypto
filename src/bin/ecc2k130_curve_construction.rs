//! Run the curve-construction study and freeze its report.
//!
//! ```text
//! cargo run --release --bin ecc2k130_curve_construction -- \
//!     --out research/notes/ecc2k130/curve_construction_20261006/results/curve_construction.json \
//!     --legacy experiments/ecc2k130_hyperelliptic_cover_boundary.json
//! ```

use clap::Parser;
use crypto_lib::cryptanalysis::curve_construction::{run, Report, RunConfig};
use serde::Serialize;
use serde_json::Value;
use std::{fs, path::PathBuf, time::Instant};

#[derive(Parser)]
#[command(
    about = "Can a curve over F_2 carrying the ECC2K-130 subgroup be built? Four routes, one report."
)]
struct Args {
    /// Write the JSON report here.
    #[arg(long)]
    out: PathBuf,
    /// Legacy frozen JSON to replay against (the Python-era boundary file).
    #[arg(long)]
    legacy: Option<PathBuf>,
    #[arg(long, default_value_t = 400)]
    samples: usize,
    #[arg(long, default_value_t = 20261006)]
    seed: u64,
    #[arg(long, default_value_t = 500)]
    max_prime: u64,
    /// Skip the genus-5 hyperelliptic enumeration (2^20 models).
    #[arg(long)]
    skip_genus5: bool,
}

#[derive(Serialize)]
struct LegacyCheck {
    file: String,
    compared: usize,
    max_abs_diff: f64,
    tolerance: f64,
    exact_values_match: bool,
    agrees: bool,
}

#[derive(Serialize)]
struct Frozen<'a> {
    #[serde(flatten)]
    report: &'a Report,
    legacy_check: Option<LegacyCheck>,
    runtime_seconds: f64,
    command: String,
}

fn num(v: &Value, path: &[&str]) -> f64 {
    let mut cur = v;
    for p in path {
        cur = &cur[*p];
    }
    cur.as_f64()
        .unwrap_or_else(|| panic!("legacy field {path:?} missing"))
}

fn compare(report: &Report, legacy: &Value, file: String) -> LegacyCheck {
    let tol = 0.006;
    let mut diffs = Vec::new();
    let t = &report.legacy_replay.target;
    diffs.push((t.log2_rho_reference - num(legacy, &["target", "log2_rho_reference"])).abs());
    diffs.push((t.log2_rho_plain - num(legacy, &["target", "log2_rho_plain"])).abs());
    let ex = &report.legacy_replay.window.exact_genus_130_jac_a;
    let lex = &legacy["boundary_E_window"]["exact_cell_if_A_were_a_jacobian"];
    for (mine, key) in [
        (ex.log2_factor_base, "log2_factor_base"),
        (ex.log2_smooth_probability, "log2_smooth_probability"),
        (ex.log2_relations, "log2_relations"),
        (ex.log2_linear_algebra, "log2_linear_algebra"),
        (ex.log2_total, "log2_total"),
    ] {
        diffs.push((mine - lex[key].as_f64().expect("number")).abs());
    }
    let mut exact = lex["smoothness_bound_b"].as_u64() == Some(ex.smoothness_bound_b as u64)
        && legacy["target"]["r"].as_str() == Some(t.r.as_str())
        && legacy["boundary_D_genus_floor"]["points_on_A_equal_r"].as_bool()
            == Some(report.legacy_replay.trace_zero.points_on_a_equal_r);
    let cells = legacy["boundary_E_window"]["cells"]
        .as_array()
        .expect("cells");
    for c in &report.legacy_replay.window.cells {
        if let Some(lc) = cells
            .iter()
            .find(|x| x["genus"].as_u64() == Some(c.genus as u64))
        {
            diffs.push((c.log2_total - lc["log2_total"].as_f64().expect("number")).abs());
            exact &= lc["smoothness_bound_b"].as_u64() == Some(c.smoothness_bound_b as u64);
        }
    }
    let lw = &legacy["boundary_E_window"]["window_high_between"];
    exact &= lw[0].as_u64() == Some(report.legacy_replay.window.crossover_between.0 as u64)
        && lw[1].as_u64() == Some(report.legacy_replay.window.crossover_between.1 as u64);
    let max = diffs.iter().cloned().fold(0.0, f64::max);
    LegacyCheck {
        file,
        compared: diffs.len(),
        max_abs_diff: max,
        tolerance: tol,
        exact_values_match: exact,
        agrees: exact && max <= tol,
    }
}

fn main() {
    let args = Args::parse();
    let cfg = RunConfig {
        ghs_samples: args.samples,
        seed: args.seed,
        max_prime: args.max_prime,
        genus5_hyperelliptic: !args.skip_genus5,
    };
    let start = Instant::now();
    let report = run(&cfg);
    let legacy_check = args.legacy.as_ref().map(|p| {
        let v: Value = serde_json::from_str(&fs::read_to_string(p).expect("read legacy file"))
            .expect("parse legacy JSON");
        compare(&report, &v, p.display().to_string())
    });
    let frozen = Frozen {
        report: &report,
        legacy_check,
        runtime_seconds: start.elapsed().as_secs_f64(),
        command: std::env::args().collect::<Vec<_>>().join(" "),
    };
    if let Some(dir) = args.out.parent() {
        fs::create_dir_all(dir).expect("create output directory");
    }
    fs::write(
        &args.out,
        serde_json::to_string_pretty(&frozen).expect("serialise") + "\n",
    )
    .expect("write report");
    println!("{}", report.verdict);
    if let Some(c) = &frozen.legacy_check {
        println!(
            "legacy replay: {} values, max |diff| {:.4}, exact fields match: {}, agrees: {}",
            c.compared, c.max_abs_diff, c.exact_values_match, c.agrees
        );
    }
    println!(
        "wrote {} in {:.1}s",
        args.out.display(),
        frozen.runtime_seconds
    );
}

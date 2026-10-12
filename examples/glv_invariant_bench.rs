//! **Endomorphism-invariant factor bases across curve families —
//! measurement.**
//!
//! Companion to `crypto_lib::cryptanalysis::glv_invariant_base`,
//! `gls_fp2`, and `RESEARCH_GLV_INVARIANT_FACTOR_BASES.md`.  For every
//! family, size and seed it generates an instance with a certified
//! order, runs the framework pipeline over the folded base and over the
//! control (the same points folded by negation alone), and prices both
//! against a counted negation-folded rho on the same instance.  For the
//! type-C families (`d7`, `d8`) it also measures the degree-2
//! endomorphism: its eigenvalue order and how much of a base it keeps in
//! the base.
//!
//! ```text
//! cargo run --release --example glv_invariant_bench -- \
//!     --families j0,j1728,generic,d7,d8 --bits 16,20,24 --seeds 2 \
//!     --oracles subtract,mitm --json experiments/23_glv_invariant_pilot.json
//! cargo run --release --example glv_invariant_bench -- \
//!     --families gls --bits 8,10,12 --seeds 2 --oracles subtract \
//!     --json experiments/23_glv_invariant_gls_pilot.json
//! ```
//!
//! `--bits` is the subgroup size for the prime families and the size
//! of `p` for `gls` (whose group has about `2·bits` bits).  Every row
//! carries `S` over the whole pipeline, the fold measured as points per
//! column, and whether the planted logarithm came back.

use std::env;
use std::fs;

use crypto_lib::cryptanalysis::gls_fp2::{generate_gls_instance, Fp2Curve, GlsInstance};
use crypto_lib::cryptanalysis::glv_invariant_base::{
    endomorphism_overlap, generate_cm_instance, glv_orbit_base, velu_degree2_endomorphisms,
    verify_endomorphism, AutomorphismGroup, CmFamily, OverlapReport,
};
use crypto_lib::cryptanalysis::ic_boundary::{
    calibrate_group, calibrate_row_ops, rho_reference_negation, Calibration, CountedGroup,
    GroupOps, PrimeInstance,
};
use crypto_lib::cryptanalysis::ic_framework::plugins::{
    GlsLineBase, GlvOrbitBase, MitmOracle, SubtractOracle,
};
use crypto_lib::cryptanalysis::ic_framework::stages::{DecompositionOracle, InstanceCtx, Targets};
use crypto_lib::cryptanalysis::ic_framework::{run_pipeline, PipelineSpec, RunReport};
use serde::Serialize;

#[derive(Serialize)]
struct Row {
    family: String,
    bits: u32,
    seed: u64,
    instance: String,
    log2_r: f64,
    r: u64,
    cofactor: u64,
    oracle: String,
    /// `size` handed to the base builder (seed abscissae; GLS: none).
    seed_abscissae: u64,
    planted: u64,
    rho_runs: usize,
    rho_s_mean: f64,
    rho_all_verified: bool,
    folded: RunReport,
    control: RunReport,
    /// `control.s / folded.s`.
    s_ratio: f64,
    /// `control.columns / folded.columns`: the fold.
    column_ratio: f64,
    /// `control trials / folded trials`.
    trial_ratio: f64,
}

#[derive(Serialize)]
struct TypeCRow {
    family: String,
    bits: u32,
    seed: u64,
    instance: String,
    r: u64,
    endomorphisms_found: usize,
    overlap: Option<OverlapReport>,
    verified: bool,
}

#[derive(Serialize)]
struct Output {
    what_this_is: &'static str,
    what_this_is_not: [&'static str; 3],
    command: String,
    host: serde_json::Value,
    rows: Vec<Row>,
    type_c: Vec<TypeCRow>,
}

fn base_size(r: u64) -> u64 {
    ((r as f64).cbrt() as u64).max(8)
}

fn rho_mean<G: CountedGroup>(
    g: &G,
    gen: G::Elt,
    target: G::Elt,
    r: u64,
    seed: u64,
    runs: usize,
) -> (f64, bool) {
    let mut sum = 0.0;
    let mut ok = true;
    for k in 0..runs {
        let res = rho_reference_negation(g, gen, target, r, seed ^ (k as u64 * 0x9E37), 1 << 34);
        sum += res.s;
        ok &= res.verified;
    }
    (sum / runs.max(1) as f64, ok)
}

fn calibration<G: CountedGroup>(g: &G, points: &[G::Elt], r: u64) -> Calibration {
    let mut calib = Calibration::default();
    calibrate_group(g, points, &mut calib);
    calibrate_row_ops(r, &mut calib);
    calib
}

fn spec_for(base: &str, oracle: &str, size: Option<u64>, fold: bool, seed: u64) -> PipelineSpec {
    let mut spec = PipelineSpec {
        factor_base: base.into(),
        oracle: oracle.into(),
        targets: Targets::Walk,
        max_trials: 50_000_000,
        seed,
        ..Default::default()
    };
    if let Some(s) = size {
        spec.factor_base_params.set("size", s.to_string());
    }
    if !fold {
        spec.factor_base_params.set("no_fold", "1");
    }
    if oracle == "mitm" {
        spec.oracle_params.set("negation_folded", "1");
    }
    spec
}

#[allow(clippy::too_many_arguments)]
fn run_prime(
    family: CmFamily,
    bits: u32,
    seed: u64,
    oracle_name: &str,
    rho_runs: usize,
    rows: &mut Vec<Row>,
    type_c: &mut Vec<TypeCRow>,
    inst: &PrimeInstance,
) {
    let g = inst.generator_point();
    let mut ops = GroupOps::default();
    let planted = 1 + (seed.wrapping_mul(0x9E37_79B9)) % (inst.r - 1);
    let target = inst.curve.mul(&mut ops, g, planted);
    let ctx = InstanceCtx {
        group: &inst.curve,
        generator: g,
        target,
        r: inst.r,
        cofactor: inst.cofactor,
        group_order: inst.group_order,
        name: inst.name.clone(),
        field_degree: None,
    };
    let sample: Vec<_> = (1..=8u64).map(|k| inst.curve.mul(&mut ops, g, k)).collect();
    let calib = calibration(&inst.curve, &sample, inst.r);
    let (rho_s, rho_ok) = rho_mean(&inst.curve, g, target, inst.r, seed, rho_runs);
    let size = base_size(inst.r);
    let base = GlvOrbitBase { instance: inst };
    let mut reports = Vec::new();
    for fold in [true, false] {
        let spec = spec_for("glv-orbit", oracle_name, Some(size), fold, seed);
        let mut subtract = SubtractOracle;
        let mut mitm = MitmOracle::new(2);
        let oracle: &mut dyn DecompositionOracle<_> = if oracle_name == "mitm" {
            &mut mitm
        } else {
            &mut subtract
        };
        let report = run_pipeline(&ctx, &spec, &base, oracle, planted, &calib, Some(rho_s))
            .expect("the configuration runs");
        eprintln!(
            "  {} {} fold={fold}: base {} cols {} pts/col {:.1} trials {} S {:.2} correct {}",
            inst.name,
            oracle_name,
            report.factor_base.signed_points,
            report.factor_base.columns,
            report.factor_base.points_per_column,
            report.decomposition.targets_tried,
            report.s,
            report.verified
        );
        reports.push(report);
    }
    let control = reports.pop().unwrap();
    let folded = reports.pop().unwrap();
    rows.push(Row {
        family: family.name().into(),
        bits,
        seed,
        instance: inst.name.clone(),
        log2_r: (inst.r as f64).log2(),
        r: inst.r,
        cofactor: inst.cofactor,
        oracle: oracle_name.into(),
        seed_abscissae: size,
        planted,
        rho_runs,
        rho_s_mean: rho_s,
        rho_all_verified: rho_ok,
        s_ratio: control.s / folded.s,
        column_ratio: control.factor_base.columns as f64 / folded.factor_base.columns as f64,
        trial_ratio: control.decomposition.targets_tried as f64
            / folded.decomposition.targets_tried.max(1) as f64,
        folded,
        control,
    });
    if matches!(family, CmFamily::D7 | CmFamily::D8) && oracle_name == "subtract" {
        let endos = velu_degree2_endomorphisms(inst);
        let (fb, _) =
            glv_orbit_base(inst, size as usize, AutomorphismGroup::Negation, true).unwrap();
        let (overlap, verified) = match endos.first() {
            Some(phi) => (
                Some(endomorphism_overlap(
                    &inst.curve,
                    &fb,
                    phi,
                    inst.r,
                    inst.group_order,
                )),
                verify_endomorphism(&inst.curve, g, inst.r, phi, 20, seed).is_ok(),
            ),
            None => (None, false),
        };
        if let Some(o) = &overlap {
            eprintln!(
                "  {} type C: {} found, ord_r(λ) = {}, {}/{} images in base (chance {:.4})",
                inst.name,
                endos.len(),
                o.eigenvalue_order,
                o.images_in_base,
                o.base_points,
                o.chance_fraction
            );
        }
        type_c.push(TypeCRow {
            family: family.name().into(),
            bits,
            seed,
            instance: inst.name.clone(),
            r: inst.r,
            endomorphisms_found: endos.len(),
            overlap,
            verified,
        });
    }
}

fn run_gls(
    bits: u32,
    seed: u64,
    oracle_name: &str,
    rho_runs: usize,
    rows: &mut Vec<Row>,
    inst: &GlsInstance,
) {
    let g = inst.generator;
    let mut ops = GroupOps::default();
    let planted = 1 + (seed.wrapping_mul(0x9E37_79B9)) % (inst.r - 1);
    let target = inst.curve.mul(&mut ops, g, planted);
    let ctx: InstanceCtx<Fp2Curve> = InstanceCtx {
        group: &inst.curve,
        generator: g,
        target,
        r: inst.r,
        cofactor: inst.cofactor,
        group_order: inst.group_order,
        name: inst.name.clone(),
        field_degree: Some(2),
    };
    let sample: Vec<_> = (1..=8u64).map(|k| inst.curve.mul(&mut ops, g, k)).collect();
    let calib = calibration(&inst.curve, &sample, inst.r);
    let (rho_s, rho_ok) = rho_mean(&inst.curve, g, target, inst.r, seed, rho_runs);
    let base = GlsLineBase { instance: inst };
    let mut reports = Vec::new();
    for fold in [true, false] {
        let spec = spec_for("gls-line", oracle_name, None, fold, seed);
        let mut subtract = SubtractOracle;
        let mut mitm = MitmOracle::new(2);
        let oracle: &mut dyn DecompositionOracle<_> = if oracle_name == "mitm" {
            &mut mitm
        } else {
            &mut subtract
        };
        let report = run_pipeline(&ctx, &spec, &base, oracle, planted, &calib, Some(rho_s))
            .expect("the configuration runs");
        eprintln!(
            "  {} {} fold={fold}: base {} cols {} pts/col {:.1} trials {} S {:.2} correct {}",
            inst.name,
            oracle_name,
            report.factor_base.signed_points,
            report.factor_base.columns,
            report.factor_base.points_per_column,
            report.decomposition.targets_tried,
            report.s,
            report.verified
        );
        reports.push(report);
    }
    let control = reports.pop().unwrap();
    let folded = reports.pop().unwrap();
    rows.push(Row {
        family: "gls".into(),
        bits,
        seed,
        instance: inst.name.clone(),
        log2_r: (inst.r as f64).log2(),
        r: inst.r,
        cofactor: inst.cofactor,
        oracle: oracle_name.into(),
        seed_abscissae: 0,
        planted,
        rho_runs,
        rho_s_mean: rho_s,
        rho_all_verified: rho_ok,
        s_ratio: control.s / folded.s,
        column_ratio: control.factor_base.columns as f64 / folded.factor_base.columns as f64,
        trial_ratio: control.decomposition.targets_tried as f64
            / folded.decomposition.targets_tried.max(1) as f64,
        folded,
        control,
    });
}

fn markdown(rows: &[Row]) -> String {
    let mut out = String::from(
        "| family | log2 r | h | oracle | base | cols fold | cols control | pts/col | trials fold | trials control | S fold | S control | S control / S fold | rho S | S fold / rho | correct |\n|:--|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|\n",
    );
    for r in rows {
        out.push_str(&format!(
            "| {} | {:.1} | {} | {} | {} | {} | {} | {:.1} | {} | {} | {:.2} | {:.2} | {:.2} | {:.2} | {:.1}× | {} |\n",
            r.family,
            r.log2_r,
            r.cofactor,
            r.oracle,
            r.folded.factor_base.signed_points,
            r.folded.factor_base.columns,
            r.control.factor_base.columns,
            r.folded.factor_base.points_per_column,
            r.folded.decomposition.targets_tried,
            r.control.decomposition.targets_tried,
            r.folded.s,
            r.control.s,
            r.s_ratio,
            r.rho_s_mean,
            r.folded.s / r.rho_s_mean,
            if r.folded.verified && r.control.verified {
                "yes"
            } else {
                "NO"
            },
        ));
    }
    out
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let mut families = vec!["j0".to_string(), "j1728".into(), "generic".into()];
    let mut bits: Vec<u32> = vec![16, 20];
    let mut seeds = 1u64;
    let mut oracles = vec!["subtract".to_string()];
    let mut rho_runs = 8usize;
    let mut max_cofactor = 8u64;
    let mut json: Option<String> = None;
    let mut i = 0;
    while i < args.len() {
        let next = |i: &mut usize| -> String {
            *i += 1;
            args[*i].clone()
        };
        match args[i].as_str() {
            "--families" => families = next(&mut i).split(',').map(str::to_string).collect(),
            "--bits" => {
                bits = next(&mut i)
                    .split(',')
                    .map(|s| s.parse().unwrap())
                    .collect()
            }
            "--seeds" => seeds = next(&mut i).parse().unwrap(),
            "--oracles" => oracles = next(&mut i).split(',').map(str::to_string).collect(),
            "--rho-runs" => rho_runs = next(&mut i).parse().unwrap(),
            "--max-cofactor" => max_cofactor = next(&mut i).parse().unwrap(),
            "--json" => json = Some(next(&mut i)),
            other => panic!("unknown argument {other}"),
        }
        i += 1;
    }
    let mut rows = Vec::new();
    let mut type_c = Vec::new();
    for family in &families {
        for &b in &bits {
            for seed in 1..=seeds {
                if family == "gls" {
                    let inst = match generate_gls_instance(b, seed, max_cofactor.max(16)) {
                        Ok(i) => i,
                        Err(e) => {
                            eprintln!("gls {b}: {e}");
                            continue;
                        }
                    };
                    eprintln!(
                        "gls p_bits={b} seed={seed}: {} r=2^{:.1} h={}",
                        inst.name,
                        (inst.r as f64).log2(),
                        inst.cofactor
                    );
                    for o in &oracles {
                        run_gls(b, seed, o, rho_runs, &mut rows, &inst);
                    }
                    continue;
                }
                let fam = CmFamily::parse(family).unwrap();
                let inst = match generate_cm_instance(fam, b, seed, max_cofactor) {
                    Ok(i) => i,
                    Err(e) => {
                        eprintln!("{family} {b}: {e}");
                        continue;
                    }
                };
                eprintln!(
                    "{family} bits={b} seed={seed}: {} r=2^{:.1} h={}",
                    inst.name,
                    (inst.r as f64).log2(),
                    inst.cofactor
                );
                for o in &oracles {
                    run_prime(fam, b, seed, o, rho_runs, &mut rows, &mut type_c, &inst);
                }
            }
        }
    }
    println!("{}", markdown(&rows));
    if !type_c.is_empty() {
        println!("\n| family | log2 r | degree-2 maps | ord_r(λ) | images in base | base points | chance fraction | verified |\n|:--|--:|--:|--:|--:|--:|--:|:--|");
        for t in &type_c {
            let (ord, inside, pts, chance) = match &t.overlap {
                Some(o) => (
                    o.eigenvalue_order.to_string(),
                    o.images_in_base.to_string(),
                    o.base_points.to_string(),
                    format!("{:.5}", o.chance_fraction),
                ),
                None => ("—".into(), "—".into(), "—".into(), "—".into()),
            };
            println!(
                "| {} | {:.1} | {} | {} | {} | {} | {} | {} |",
                t.family,
                (t.r as f64).log2(),
                t.endomorphisms_found,
                ord,
                inside,
                pts,
                chance,
                t.verified
            );
        }
    }
    if let Some(path) = json {
        let out = Output {
            what_this_is: "Endomorphism-invariant factor bases (glv-orbit, gls-line) against the negation-only control on the same points, every phase priced in group-addition equivalents, S = total / sqrt(r), with a counted negation-folded rho on the same instance and target.",
            what_this_is_not: [
                "not a speed claim: the fold divides the columns by the automorphism group's order over 2 and the relation count moves with its floor (AGENTS.md section 3: engineering)",
                "not a claim about any deployed curve: toy instances with certified orders",
                "square roots and Legendre symbols in the base build are counted but unpriced on these generated curves, as in `ic bench` on a curve outside the pinned table",
            ],
            command: format!("glv_invariant_bench {}", args.join(" ")),
            host: serde_json::json!({
                "os": std::env::consts::OS,
                "arch": std::env::consts::ARCH,
                "threads": std::thread::available_parallelism().map(|n| n.get()).unwrap_or(0),
            }),
            rows,
            type_c,
        };
        fs::write(&path, serde_json::to_string_pretty(&out).unwrap()).unwrap();
        eprintln!("wrote {path}");
    }
}

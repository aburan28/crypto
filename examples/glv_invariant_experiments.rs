//! **Experiments E1–E15 of `RESEARCH_GLV_INVARIANT_FACTOR_BASES.md` —
//! measurement.**
//!
//! ```text
//! cargo run --release --example glv_invariant_experiments -- --exp e1 --bits 16,20,24,28,32 --seeds 6 --json experiments/23_glv_invariant_e1.json
//! cargo run --release --example glv_invariant_experiments -- --exp e2 --bits 8,10,12,14,16 --seeds 6 --json experiments/23_glv_invariant_e2.json
//! cargo run --release --example glv_invariant_experiments -- --exp e3 --bits 8,10,12 --seeds 4 --json experiments/23_glv_invariant_e3.json
//! cargo run --release --example glv_invariant_experiments -- --exp e4 --bits 16,20,24,28 --seeds 4 --json experiments/23_glv_invariant_e4.json
//! cargo run --release --example glv_invariant_experiments -- --exp e5 --bits 7,8,9,10,11 --seeds 4 --json experiments/23_glv_invariant_e5.json
//! cargo run --release --example glv_invariant_experiments -- --exp e6 --bits 16,20,24,28 --seeds 6 --json experiments/23_glv_invariant_e6.json
//! cargo run --release --example glv_invariant_experiments -- --exp e7 --bits 18,22,26 --seeds 4 --json experiments/23_glv_invariant_e7.json
//! ```
//!
//! Every experiment's rows are `serde_json` objects; the tables are
//! printed from the frozen files by `scripts/glv_invariant_tables.py`.
//! What each experiment measures, its arms and its boundary are in the
//! note's §4; the E1/E2/E5/E7 rows come from
//! `glv_invariant_experiments::full_rank_stream`, one target stream
//! feeding the folded and the control matrix at once.

use std::env;
use std::fs;
use std::time::Instant;

use crypto_lib::cryptanalysis::ext_curve::{take_field_counters, Fp2, Fp3};
use crypto_lib::cryptanalysis::fghr_line::{
    fghr_line_base, fghr_polynomials, generate_fghr_instance, FghrFold, FghrOracle, YLineS4Oracle,
};
use crypto_lib::cryptanalysis::gls_fp2::{
    generate_gls_instance_of, gls_line_base_by, GlsFamily, GlsFold, GlsInstance,
};
use crypto_lib::cryptanalysis::glv_invariant_base::{
    automorphism_generators, eigenvalue_order, endomorphism_overlap, generate_cm_instance,
    glv_orbit_base, rho_reference_folded, velu_degree2_endomorphisms, velu_degree3_endomorphisms,
    verify_endomorphism, AutomorphismGroup, CmFamily, Endomorphism, EndomorphismClasses, Negation,
};
use crypto_lib::cryptanalysis::glv_invariant_experiments::{
    classes_for, e1_stream, full_rank_stream, full_rank_stream_until, StopRule, StreamReport,
};
use crypto_lib::cryptanalysis::ic_boundary::{
    koblitz_factor_base, koblitz_instance, rho_reference_negation, BinaryGroup, ColumnFold,
    CountedGroup, GroupOps, OracleCounters, PrimeInstance,
};
use crypto_lib::cryptanalysis::ic_framework::plugins::{MitmOracle, SubtractOracle};
use crypto_lib::cryptanalysis::ic_framework::stages::{DecompositionOracle, InstanceCtx, Params};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    all_factors_of_x_n_minus_1, build_frobenius_factor_base_from_divisor,
};
use crypto_lib::cryptanalysis::line_oracle::LineOracle;
use crypto_lib::cryptanalysis::line_s4_oracle::LineS4Oracle;
use crypto_lib::cryptanalysis::orbit_pair_table::{
    OrbitMitmOracle, PrimePowerKey, SubfieldOrbitKey,
};
use crypto_lib::cryptanalysis::subfield_fp3::{
    generate_subfield_instance, subfield_line_base, SubfieldFold,
};
use serde_json::{json, Value};

struct Opts {
    exp: String,
    bits: Vec<u32>,
    seeds: u64,
    rho_runs: usize,
    /// Steps a single rho reference may spend before it is reported as
    /// capped (`verified: false`); a cap is never evidence.
    rho_max_steps: u64,
    max_trials: u64,
    json: Option<String>,
}

fn base_size(r: u64) -> usize {
    ((r as f64).cbrt() as usize).max(8)
}

fn planted_for(seed: u64, r: u64) -> u64 {
    1 + seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) % (r - 1)
}

fn stream_json(rep: &StreamReport) -> Value {
    serde_json::to_value(rep).unwrap()
}

fn prime_ctx(
    inst: &PrimeInstance,
    planted: u64,
) -> (
    InstanceCtx<'_, crypto_lib::cryptanalysis::ic_boundary::PrimeCurve>,
    crypto_lib::cryptanalysis::ic_boundary::PrimePoint,
) {
    let g = inst.generator_point();
    let mut ops = GroupOps::default();
    let target = inst.curve.mul(&mut ops, g, planted);
    (
        InstanceCtx {
            group: &inst.curve,
            generator: g,
            target,
            r: inst.r,
            cofactor: inst.cofactor,
            group_order: inst.group_order,
            name: inst.name.clone(),
            field_degree: None,
        },
        target,
    )
}

// ── E1: automorphism fold at full rank ─────────────────────────────

fn e1(o: &Opts) -> Vec<Value> {
    let mut rows = Vec::new();
    for family in [CmFamily::J0, CmFamily::J1728, CmFamily::Generic] {
        for &bits in &o.bits {
            if family == CmFamily::Generic && bits > 26 {
                eprintln!("e1 generic {bits}: point counting is O(p); skipped above 26 bits");
                continue;
            }
            for seed in 1..=o.seeds {
                let inst = match generate_cm_instance(family, bits, seed, 8) {
                    Ok(i) => i,
                    Err(e) => {
                        eprintln!("e1 {} {bits} seed {seed}: {e}", family.name());
                        continue;
                    }
                };
                let size = base_size(inst.r);
                let (folded, frep) =
                    glv_orbit_base(&inst, size, AutomorphismGroup::Auto, true).unwrap();
                let (control, _) =
                    glv_orbit_base(&inst, size, AutomorphismGroup::Auto, false).unwrap();
                let planted = planted_for(seed, inst.r);
                let (ctx, target) = prime_ctx(&inst, planted);
                let mut oracles: Vec<(&str, Box<dyn DecompositionOracle<_>>)> =
                    vec![("mitm", Box::new(MitmOracle::new(2)))];
                if bits <= 22 {
                    oracles.push(("subtract", Box::new(SubtractOracle)));
                }
                for (name, mut oracle) in oracles {
                    let started = Instant::now();
                    let mut params = Params::default();
                    params.set("negation_folded", "1");
                    let mut prep_ops = GroupOps::default();
                    oracle
                        .prepare(&ctx, &folded, &params, &mut prep_ops)
                        .unwrap();
                    let rep = e1_stream(
                        &inst.curve,
                        ctx.generator,
                        target,
                        inst.r,
                        inst.cofactor,
                        planted,
                        &folded,
                        &control,
                        seed,
                        o.max_trials,
                        |ops, ctr, pt| oracle.decompose(&ctx, &folded, ops, ctr, pt),
                    )
                    .unwrap();
                    eprintln!(
                        "e1 {} 2^{:.1} {name}: cols {}/{} deficiency {}/{} square rel {:?}/{:?} (ratio {:?}) first-pin {:?}/{:?} trials {} ok {}/{} [{:.1}s]",
                        family.name(),
                        (inst.r as f64).log2(),
                        rep.folded.columns,
                        rep.control.columns,
                        rep.folded.deficiency_total,
                        rep.control.deficiency_total,
                        rep.folded.square_relations,
                        rep.control.square_relations,
                        rep.square_ratio,
                        rep.folded.first_pin_relations,
                        rep.control.first_pin_relations,
                        rep.trials,
                        rep.folded.verified,
                        rep.control.verified,
                        started.elapsed().as_secs_f64()
                    );
                    rows.push(json!({
                        "experiment": "e1",
                        "family": family.name(),
                        "bits": bits,
                        "seed": seed,
                        "instance": inst.name,
                        "log2_r": (inst.r as f64).log2(),
                        "r": inst.r,
                        "cofactor": inst.cofactor,
                        "seed_abscissae": size,
                        "oracle": name,
                        "prep_group_ops": prep_ops,
                        "fold_generators": frep.generators,
                        "planted": planted,
                        "stream": stream_json(&rep),
                        "wall_seconds": started.elapsed().as_secs_f64(),
                    }));
                }
            }
        }
    }
    rows
}

// ── E2: GLS line against Koblitz orbit, same unit ──────────────────

fn gls_stream(
    inst: &GlsInstance,
    folded: &crypto_lib::cryptanalysis::ic_boundary::FactorBase<
        crypto_lib::cryptanalysis::gls_fp2::Fp2Point,
    >,
    control: &crypto_lib::cryptanalysis::ic_boundary::FactorBase<
        crypto_lib::cryptanalysis::gls_fp2::Fp2Point,
    >,
    seed: u64,
    max_trials: u64,
    use_subtract: bool,
) -> (StreamReport, u64, u64) {
    let planted = planted_for(seed, inst.r);
    let mut ops = GroupOps::default();
    let target = inst.curve.mul(&mut ops, inst.generator, planted);
    let ctx = InstanceCtx {
        group: &inst.curve,
        generator: inst.generator,
        target,
        r: inst.r,
        cofactor: inst.cofactor,
        group_order: inst.group_order,
        name: inst.name.clone(),
        field_degree: Some(2),
    };
    let mut line: LineOracle<Fp2> = LineOracle::new(inst.curve.f, inst.line);
    let mut subtract = SubtractOracle;
    line.prepare(&ctx, folded, &Params::default(), &mut GroupOps::default())
        .unwrap();
    let rep = full_rank_stream(
        &inst.curve,
        inst.generator,
        target,
        inst.r,
        inst.cofactor,
        planted,
        folded,
        control,
        seed,
        max_trials,
        2,
        None,
        |ops, ctr, pt| {
            if use_subtract {
                subtract.decompose(&ctx, folded, ops, ctr, pt)
            } else {
                line.decompose(&ctx, folded, ops, ctr, pt)
            }
        },
    )
    .unwrap();
    let totals = line.solver_totals().unwrap();
    (rep, totals.ops, totals.calls)
}

/// The line oracle against `subtract`, target by target (AGENTS.md §6).
fn gls_agreement(
    inst: &GlsInstance,
    fb: &crypto_lib::cryptanalysis::ic_boundary::FactorBase<
        crypto_lib::cryptanalysis::gls_fp2::Fp2Point,
    >,
    targets: u64,
) -> (u64, u64, u64) {
    let mut ops = GroupOps::default();
    let target = inst.curve.mul(&mut ops, inst.generator, 3);
    let ctx = InstanceCtx {
        group: &inst.curve,
        generator: inst.generator,
        target,
        r: inst.r,
        cofactor: inst.cofactor,
        group_order: inst.group_order,
        name: inst.name.clone(),
        field_degree: Some(2),
    };
    let mut line: LineOracle<Fp2> = LineOracle::new(inst.curve.f, inst.line);
    line.prepare(&ctx, fb, &Params::default(), &mut ops)
        .unwrap();
    let mut subtract = SubtractOracle;
    let (mut ca, mut cb) = (OracleCounters::default(), OracleCounters::default());
    let (mut agree, mut disagree, mut hits) = (0u64, 0u64, 0u64);
    for k in 1..=targets {
        let pt = inst.curve.mul(&mut ops, inst.generator, k);
        let a = line.decompose(&ctx, fb, &mut ops, &mut ca, pt);
        let b = subtract.decompose(&ctx, fb, &mut ops, &mut cb, pt);
        if a.is_some() == b.is_some() {
            agree += 1;
        } else {
            disagree += 1;
        }
        if a.is_some() {
            hits += 1;
        }
    }
    (agree, disagree, hits)
}

fn koblitz_divisor_for(n: u32, target_dim: u32) -> Option<Vec<usize>> {
    let factors = all_factors_of_x_n_minus_1(n);
    let degs: Vec<u32> = factors.iter().map(|&f| 63 - f.leading_zeros()).collect();
    let count = factors.len().min(16);
    let mut candidates: Vec<(u32, Vec<usize>)> = Vec::new();
    for mask in 1u32..(1u32 << count) {
        let idx: Vec<usize> = (0..count).filter(|&i| (mask >> i) & 1 == 1).collect();
        let total: u32 = idx.iter().map(|&i| degs[i]).sum();
        if total >= 4 && total < n && (!n.is_multiple_of(total) || total == 1) {
            candidates.push((total.abs_diff(target_dim), idx));
        }
    }
    candidates.sort();
    candidates.into_iter().map(|(_, idx)| idx).next()
}

fn e2(o: &Opts) -> Vec<Value> {
    let mut rows = Vec::new();
    // GLS arm.
    for &p_bits in &o.bits {
        for seed in 1..=o.seeds {
            let inst = match generate_gls_instance_of(GlsFamily::Generic, p_bits, seed, 64) {
                Ok(i) => i,
                Err(e) => {
                    eprintln!("e2 gls {p_bits} seed {seed}: {e}");
                    continue;
                }
            };
            let started = Instant::now();
            let (folded, frep) = gls_line_base_by(&inst, GlsFold::Psi).unwrap();
            let (control, _) = gls_line_base_by(&inst, GlsFold::Negation).unwrap();
            let (rep, oracle_muls, oracle_calls) =
                gls_stream(&inst, &folded, &control, seed, o.max_trials, false);
            let agreement = if p_bits <= 10 {
                let (a, d, h) = gls_agreement(&inst, &folded, 300);
                Some(json!({"targets": 300, "agree": a, "disagree": d, "hits": h}))
            } else {
                None
            };
            eprintln!(
                "e2 gls p=2^{p_bits} r=2^{:.1}: cols {}/{} square {:?}/{:?} ratio {:?} trials {} muls/call {:.0} ok {}/{} [{:.1}s]",
                (inst.r as f64).log2(),
                rep.folded.columns,
                rep.control.columns,
                rep.folded.square_relations,
                rep.control.square_relations,
                rep.square_ratio,
                rep.trials,
                oracle_muls as f64 / oracle_calls.max(1) as f64,
                rep.folded.verified,
                rep.control.verified,
                started.elapsed().as_secs_f64()
            );
            rows.push(json!({
                "experiment": "e2",
                "family": "gls",
                "p_bits": p_bits,
                "seed": seed,
                "instance": inst.name,
                "log2_r": (inst.r as f64).log2(),
                "r": inst.r,
                "cofactor": inst.cofactor,
                "eigenvalue_order": 4,
                "oracle": "line-resultant",
                "oracle_fp_muls": oracle_muls,
                "oracle_calls": oracle_calls,
                "fold_generators": frep.generators,
                "agreement_with_subtract": agreement,
                "stream": stream_json(&rep),
                "wall_seconds": started.elapsed().as_secs_f64(),
            }));
        }
    }
    // Koblitz arm: the ledger's fold, on the same driver.
    // Degrees whose x^n − 1 has factors of intermediate degree, so an
    // invariant subspace of about r^{1/3} abscissae exists; prime n with
    // ord_n(2) = n − 1 (17, 19, 29) have only the (n − 1)-dimensional one.
    for n in 13u32..=33 {
        let Some(inst) = koblitz_instance(1, n).or_else(|| koblitz_instance(0, n)) else {
            eprintln!("e2 koblitz n={n}: no instance");
            continue;
        };
        let Some(kc) = inst.koblitz.as_ref() else {
            continue;
        };
        let target_dim = ((inst.r as f64).log2() / 3.0).ceil() as u32 + 2;
        let Some(idx) = koblitz_divisor_for(n, target_dim) else {
            eprintln!("e2 koblitz n={n}: no divisor");
            continue;
        };
        let Some(frob) = build_frobenius_factor_base_from_divisor(kc, &idx) else {
            eprintln!("e2 koblitz n={n}: divisor {idx:?} gives no base");
            continue;
        };
        if frob.points.len() > 60_000 {
            eprintln!("e2 koblitz n={n}: the invariant subspace of dimension {} carries {} points, too many for a pair table", frob.ell, frob.points.len());
            continue;
        }
        let Some(folded) = koblitz_factor_base(
            &inst,
            &frob,
            ColumnFold::SignedFrobeniusOrbit,
            "orbit".into(),
        ) else {
            continue;
        };
        let Some(control) =
            koblitz_factor_base(&inst, &frob, ColumnFold::Abscissa, "abscissa".into())
        else {
            continue;
        };
        let group = BinaryGroup(&inst.fast);
        for seed in 1..=o.seeds {
            let started = Instant::now();
            let planted = planted_for(seed, inst.r);
            let mut ops = GroupOps::default();
            let target = group.mul(&mut ops, inst.generator, planted);
            let ctx = InstanceCtx {
                group: &group,
                generator: inst.generator,
                target,
                r: inst.r,
                cofactor: inst.cofactor,
                group_order: inst.group_order,
                name: inst.name.clone(),
                field_degree: Some(n),
            };
            let mut oracle = MitmOracle::new(2);
            let mut params = Params::default();
            params.set("negation_folded", "1");
            oracle.prepare(&ctx, &folded, &params, &mut ops).unwrap();
            let rep = full_rank_stream(
                &group,
                inst.generator,
                target,
                inst.r,
                inst.cofactor,
                planted,
                &folded,
                &control,
                seed,
                o.max_trials,
                2,
                None,
                |ops, ctr, pt| oracle.decompose(&ctx, &folded, ops, ctr, pt),
            );
            let rep = match rep {
                Ok(r) => r,
                Err(e) => {
                    eprintln!("e2 koblitz n={n}: {e}");
                    break;
                }
            };
            eprintln!(
                "e2 koblitz n={n} r=2^{:.1} dim {}: cols {}/{} pts/col {:.1} square {:?}/{:?} ratio {:?} ok {}/{} [{:.1}s]",
                (inst.r as f64).log2(),
                frob.ell,
                rep.folded.columns,
                rep.control.columns,
                rep.folded.points_per_column,
                rep.folded.square_relations,
                rep.control.square_relations,
                rep.square_ratio,
                rep.folded.verified,
                rep.control.verified,
                started.elapsed().as_secs_f64()
            );
            rows.push(json!({
                "experiment": "e2",
                "family": "koblitz",
                "n": n,
                "seed": seed,
                "instance": inst.name,
                "log2_r": (inst.r as f64).log2(),
                "r": inst.r,
                "cofactor": inst.cofactor,
                "eigenvalue_order": n,
                "divisor": idx,
                "subspace_dimension": frob.ell,
                "oracle": "mitm",
                "stream": stream_json(&rep),
                "wall_seconds": started.elapsed().as_secs_f64(),
            }));
        }
    }
    rows
}

// ── E3: composite groups on twisted CM curves ──────────────────────

fn e3(o: &Opts) -> Vec<Value> {
    let mut rows = Vec::new();
    for family in [GlsFamily::J0, GlsFamily::J1728] {
        for &p_bits in &o.bits {
            for seed in 1..=o.seeds {
                let max_cofactor = if family == GlsFamily::J1728 {
                    16u64 << p_bits
                } else {
                    64
                };
                let inst = match generate_gls_instance_of(family, p_bits, seed, max_cofactor) {
                    Ok(i) => i,
                    Err(e) => {
                        eprintln!("e3 {} {p_bits} seed {seed}: {e}", family.name());
                        continue;
                    }
                };
                let started = Instant::now();
                let aut = inst.automorphism.as_ref().unwrap();
                let mut folds = serde_json::Map::new();
                for fold in [
                    GlsFold::Negation,
                    GlsFold::Psi,
                    GlsFold::Aut,
                    GlsFold::PsiAut,
                ] {
                    let (fb, rep) = gls_line_base_by(&inst, fold).unwrap();
                    folds.insert(
                        fold.name().into(),
                        json!({"columns": fb.columns, "points": fb.points.len(), "points_per_column": rep.points_per_orbit, "points_added_by_closure": rep.points_added}),
                    );
                }
                let (both, _) = gls_line_base_by(&inst, GlsFold::PsiAut).unwrap();
                let (psi_only, _) = gls_line_base_by(&inst, GlsFold::Psi).unwrap();
                let (control, _) = gls_line_base_by(&inst, GlsFold::Negation).unwrap();
                // A stream needs more distinct targets than relations; a
                // j = 1728 twist has r ≈ p, fewer than its control's columns
                // times the coupon-collector factor, so only the folds are
                // measured there.
                let enough_targets = inst.r >= 16 * control.columns as u64;
                let same_points = both.points.len() == control.points.len() && enough_targets;
                let (stream_both, muls_b, calls_b) = if same_points {
                    let (r, m, c) = gls_stream(&inst, &both, &control, seed, o.max_trials, false);
                    (Some(stream_json(&r)), m, c)
                } else {
                    (None, 0, 0)
                };
                let stream_psi = enough_targets
                    .then(|| gls_stream(&inst, &psi_only, &control, seed, o.max_trials, false).0);
                eprintln!(
                    "e3 {} p=2^{p_bits} r=2^{:.1} h={}: aut ord {} psi ord {} coincide {} | folds {}",
                    family.name(),
                    (inst.r as f64).log2(),
                    inst.cofactor,
                    eigenvalue_order(aut.eigenvalue, inst.r),
                    eigenvalue_order(inst.psi.eigenvalue, inst.r),
                    aut.eigenvalue == inst.psi.eigenvalue || aut.eigenvalue == inst.r - inst.psi.eigenvalue,
                    serde_json::to_string(&folds).unwrap()
                );
                rows.push(json!({
                    "experiment": "e3",
                    "family": family.name(),
                    "p_bits": p_bits,
                    "seed": seed,
                    "instance": inst.name,
                    "log2_r": (inst.r as f64).log2(),
                    "r": inst.r,
                    "cofactor": inst.cofactor,
                    "psi_eigenvalue_order": eigenvalue_order(inst.psi.eigenvalue, inst.r),
                    "aut_eigenvalue_order": eigenvalue_order(aut.eigenvalue, inst.r),
                    "aut_is_plus_minus_psi_on_subgroup": aut.eigenvalue == inst.psi.eigenvalue || aut.eigenvalue == inst.r - inst.psi.eigenvalue,
                    "folds": folds,
                    "stream_psi_aut_vs_negation": stream_both,
                    "oracle_fp_muls": muls_b,
                    "oracle_calls": calls_b,
                    "stream_psi_vs_negation": stream_psi.as_ref().map(stream_json),
                    "enough_targets_for_a_stream": enough_targets,
                    "wall_seconds": started.elapsed().as_secs_f64(),
                }));
            }
        }
    }
    rows
}

// ── E4: type C, degree 2 and 3 ─────────────────────────────────────

fn e4(o: &Opts) -> Vec<Value> {
    let mut rows = Vec::new();
    for (family, degree) in [
        (CmFamily::D7, 2u64),
        (CmFamily::D8, 2),
        (CmFamily::J1728, 2),
        (CmFamily::D11, 3),
        (CmFamily::J0, 3),
    ] {
        for &bits in &o.bits {
            for seed in 1..=o.seeds {
                let inst = match generate_cm_instance(family, bits, seed, 16) {
                    Ok(i) => i,
                    Err(e) => {
                        eprintln!("e4 {} {bits} seed {seed}: {e}", family.name());
                        continue;
                    }
                };
                let size = base_size(inst.r);
                let (fb, _) =
                    glv_orbit_base(&inst, size, AutomorphismGroup::Negation, true).unwrap();
                let g = inst.generator_point();
                let endos: Vec<Box<dyn Endomorphism<_>>> = if degree == 2 {
                    velu_degree2_endomorphisms(&inst)
                        .into_iter()
                        .map(|e| Box::new(e) as Box<dyn Endomorphism<_>>)
                        .collect()
                } else {
                    velu_degree3_endomorphisms(&inst)
                        .into_iter()
                        .map(|e| Box::new(e) as Box<dyn Endomorphism<_>>)
                        .collect()
                };
                let mut maps = Vec::new();
                for e in &endos {
                    let check = verify_endomorphism(&inst.curve, g, inst.r, e.as_ref(), 20, seed);
                    let overlap = endomorphism_overlap(
                        &inst.curve,
                        &fb,
                        e.as_ref(),
                        inst.r,
                        inst.group_order,
                    );
                    maps.push(json!({
                        "name": e.name(),
                        "degree": e.degree(),
                        "verified": check.is_ok(),
                        "overlap": overlap,
                    }));
                }
                eprintln!(
                    "e4 {} deg {degree} 2^{:.1}: {} maps, ord {:?}, images in base {:?} of {}",
                    family.name(),
                    (inst.r as f64).log2(),
                    endos.len(),
                    endos
                        .iter()
                        .map(|e| eigenvalue_order(e.eigenvalue(), inst.r))
                        .collect::<Vec<_>>(),
                    maps.iter()
                        .map(|m| m["overlap"]["images_in_base"].as_u64().unwrap_or(0))
                        .collect::<Vec<_>>(),
                    fb.points.len()
                );
                rows.push(json!({
                    "experiment": "e4",
                    "family": family.name(),
                    "degree": degree,
                    "bits": bits,
                    "seed": seed,
                    "instance": inst.name,
                    "log2_r": (inst.r as f64).log2(),
                    "r": inst.r,
                    "cofactor": inst.cofactor,
                    "base_points": fb.points.len(),
                    "maps": maps,
                }));
            }
        }
    }
    rows
}

// ── E5: subfield curves over F_{p³} ────────────────────────────────

fn e5(o: &Opts) -> Vec<Value> {
    let mut rows = Vec::new();
    for j0 in [false, true] {
        for &p_bits in &o.bits {
            for seed in 1..=o.seeds {
                let ratio = if j0 { 16u64 << p_bits } else { 8 };
                let inst = match generate_subfield_instance(p_bits, seed, j0, ratio) {
                    Ok(i) => i,
                    Err(e) => {
                        eprintln!(
                            "e5 {} {p_bits} seed {seed}: {e}",
                            if j0 { "j0" } else { "generic" }
                        );
                        continue;
                    }
                };
                let started = Instant::now();
                let mut folds = serde_json::Map::new();
                let fold_list: Vec<SubfieldFold> = if j0 {
                    vec![
                        SubfieldFold::Negation,
                        SubfieldFold::Frobenius,
                        SubfieldFold::Zeta,
                        SubfieldFold::FrobeniusZeta,
                    ]
                } else {
                    vec![SubfieldFold::Negation, SubfieldFold::Frobenius]
                };
                for fold in &fold_list {
                    let (fb, rep) = subfield_line_base(&inst, *fold).unwrap();
                    folds.insert(fold.name().into(), json!({"columns": fb.columns, "points": fb.points.len(), "points_per_column": rep.points_per_orbit}));
                }
                // On j = 0 the π-eigenline carrying the target subgroup is the
                // F_p-points of a cubic twist: the base is the whole subgroup,
                // every target decomposes in about p/2 ways and x_R lies on the
                // line, so the descent degenerates.  Folds are recorded; the
                // relation stream is not an index calculus there.
                let stream = if j0 {
                    None
                } else {
                    let (folded, _) = subfield_line_base(&inst, SubfieldFold::Frobenius).unwrap();
                    let (control, _) = subfield_line_base(&inst, SubfieldFold::Negation).unwrap();
                    let planted = planted_for(seed, inst.r);
                    let mut ops = GroupOps::default();
                    let target = inst.curve.mul(&mut ops, inst.generator, planted);
                    let ctx = InstanceCtx {
                        group: &inst.curve,
                        generator: inst.generator,
                        target,
                        r: inst.r,
                        cofactor: inst.cofactor,
                        group_order: inst.group_order,
                        name: inst.name.clone(),
                        field_degree: Some(3),
                    };
                    let mut line: LineOracle<Fp3> = LineOracle::new(inst.curve.f, inst.line);
                    line.prepare(&ctx, &folded, &Params::default(), &mut ops)
                        .unwrap();
                    let rep = full_rank_stream(
                        &inst.curve,
                        inst.generator,
                        target,
                        inst.r,
                        inst.cofactor,
                        planted,
                        &folded,
                        &control,
                        seed,
                        o.max_trials,
                        2,
                        None,
                        |ops, ctr, pt| line.decompose(&ctx, &folded, ops, ctr, pt),
                    )
                    .unwrap();
                    let totals = line.solver_totals().unwrap();
                    eprintln!(
                        "e5 generic p=2^{p_bits} r=2^{:.1} h/N1={}: square {:?}/{:?} ratio {:?} rank@cols {:?}/{:?} trials {} ok {}/{} degenerate {} [{:.1}s]",
                        (inst.r as f64).log2(),
                        inst.cofactor / inst.base_order,
                        rep.folded.square_relations,
                        rep.control.square_relations,
                        rep.square_ratio,
                        rep.folded.rank_fraction_at_columns,
                        rep.control.rank_fraction_at_columns,
                        rep.trials,
                        rep.folded.verified,
                        rep.control.verified,
                        line.degenerate,
                        started.elapsed().as_secs_f64()
                    );
                    Some(json!({
                        "stream": stream_json(&rep),
                        "oracle": "line-resultant",
                        "oracle_fp_muls": totals.ops,
                        "oracle_calls": totals.calls,
                        "oracle_degenerate": line.degenerate,
                    }))
                };
                eprintln!(
                    "e5 {} p=2^{p_bits} r=2^{:.1} h/N1={}: folds {}",
                    if j0 { "j0" } else { "generic" },
                    (inst.r as f64).log2(),
                    inst.cofactor / inst.base_order,
                    serde_json::to_string(&folds).unwrap()
                );
                rows.push(json!({
                    "experiment": "e5",
                    "family": if j0 { "subfield-j0" } else { "subfield" },
                    "p_bits": p_bits,
                    "seed": seed,
                    "instance": inst.name,
                    "log2_r": (inst.r as f64).log2(),
                    "r": inst.r,
                    "base_order": inst.base_order,
                    "group_order": inst.group_order,
                    "cofactor": inst.cofactor,
                    "frobenius_eigenvalue_order": eigenvalue_order(inst.frobenius.eigenvalue, inst.r),
                    "zeta_eigenvalue_order": inst.zeta.as_ref().map(|z| eigenvalue_order(z.eigenvalue, inst.r)),
                    "zeta_is_frobenius_power": inst.zeta.as_ref().map(|z| {
                        let l = inst.frobenius.eigenvalue;
                        z.eigenvalue == l || z.eigenvalue == ((l as u128 * l as u128) % inst.r as u128) as u64
                    }),
                    "folds": folds,
                    "relation_stream": stream,
                    "wall_seconds": started.elapsed().as_secs_f64(),
                }));
            }
        }
    }
    rows
}

// ── E6: the matched folded rho ─────────────────────────────────────

fn e6(o: &Opts) -> Vec<Value> {
    let mut rows = Vec::new();
    for family in [CmFamily::J0, CmFamily::J1728] {
        for &bits in &o.bits {
            for seed in 1..=o.seeds {
                let inst = match generate_cm_instance(family, bits, seed, 8) {
                    Ok(i) => i,
                    Err(e) => {
                        eprintln!("e6 {} {bits} seed {seed}: {e}", family.name());
                        continue;
                    }
                };
                let gens = automorphism_generators(&inst, AutomorphismGroup::Auto).unwrap();
                let refs: Vec<&dyn Endomorphism<_>> = gens.iter().map(|b| b.as_ref()).collect();
                let g = inst.generator_point();
                let planted = planted_for(seed, inst.r);
                let mut ops = GroupOps::default();
                let target = inst.curve.mul(&mut ops, g, planted);
                let started = Instant::now();
                let mut walks = Vec::new();
                let (mut s_neg, mut s_fold, mut st_neg, mut st_fold, mut all_ok) =
                    (0.0, 0.0, 0u64, 0u64, true);
                for k in 0..o.rho_runs as u64 {
                    let ws = seed ^ (k * 0x9E37);
                    let n = rho_reference_negation(&inst.curve, g, target, inst.r, ws, 1 << 40);
                    let f =
                        rho_reference_folded(&inst.curve, g, target, inst.r, ws, 1 << 40, &refs)
                            .unwrap();
                    all_ok &= n.verified && f.verified;
                    s_neg += n.s;
                    s_fold += f.s;
                    st_neg += n.steps;
                    st_fold += f.steps;
                    walks.push(json!({"walk_seed": ws, "negation": n, "folded": f}));
                }
                let runs = o.rho_runs.max(1) as f64;
                eprintln!(
                    "e6 {} 2^{:.1}: A={} S neg {:.3} folded {:.3} ratio {:.2} (steps ratio {:.2}, expected sqrt(A/2) {:.2}) ok {} [{:.1}s]",
                    family.name(),
                    (inst.r as f64).log2(),
                    refs.len().max(1),
                    s_neg / runs,
                    s_fold / runs,
                    s_neg / s_fold,
                    st_neg as f64 / st_fold.max(1) as f64,
                    (if family == CmFamily::J0 { 3.0f64 } else { 2.0 }).sqrt(),
                    all_ok,
                    started.elapsed().as_secs_f64()
                );
                rows.push(json!({
                    "experiment": "e6",
                    "family": family.name(),
                    "bits": bits,
                    "seed": seed,
                    "instance": inst.name,
                    "log2_r": (inst.r as f64).log2(),
                    "r": inst.r,
                    "cofactor": inst.cofactor,
                    "planted": planted,
                    "rho_runs": o.rho_runs,
                    "negation_s_mean": s_neg / runs,
                    "folded_s_mean": s_fold / runs,
                    "s_ratio": s_neg / s_fold,
                    "steps_ratio": st_neg as f64 / st_fold.max(1) as f64,
                    "expected_ratio": (if family == CmFamily::J0 { 3.0f64 } else { 2.0 }).sqrt(),
                    "all_verified": all_ok,
                    "walks": walks,
                    "wall_seconds": started.elapsed().as_secs_f64(),
                }));
            }
        }
    }
    rows
}

// ── E7: three summands and orbit duplicates ────────────────────────

fn e7(o: &Opts) -> Vec<Value> {
    let mut rows = Vec::new();
    for &bits in &o.bits {
        for seed in 1..=o.seeds {
            let inst = match generate_cm_instance(CmFamily::J0, bits, seed, 8) {
                Ok(i) => i,
                Err(e) => {
                    eprintln!("e7 j0 {bits} seed {seed}: {e}");
                    continue;
                }
            };
            // A smaller base for three summands: about r^{1/4} seeds.
            let size = ((inst.r as f64).powf(0.25) as usize).max(6);
            let (folded, _) = glv_orbit_base(&inst, size, AutomorphismGroup::Auto, true).unwrap();
            let (control, _) = glv_orbit_base(&inst, size, AutomorphismGroup::Auto, false).unwrap();
            let gens = automorphism_generators(&inst, AutomorphismGroup::Auto).unwrap();
            let refs: Vec<&dyn Endomorphism<_>> = gens.iter().map(|b| b.as_ref()).collect();
            let g = inst.generator_point();
            let classes = classes_for(&inst.curve, g, inst.r, &refs).unwrap();
            let planted = planted_for(seed, inst.r);
            let (ctx, target) = prime_ctx(&inst, planted);
            for m in [2u32, 3] {
                let started = Instant::now();
                let mut oracle = MitmOracle::new(m);
                let mut params = Params::default();
                params.set("negation_folded", "1");
                let mut prep_ops = GroupOps::default();
                oracle
                    .prepare(&ctx, &folded, &params, &mut prep_ops)
                    .unwrap();
                let rep = full_rank_stream(
                    &inst.curve,
                    g,
                    target,
                    inst.r,
                    inst.cofactor,
                    planted,
                    &folded,
                    &control,
                    seed,
                    o.max_trials,
                    m,
                    Some(&classes),
                    |ops, ctr, pt| oracle.decompose(&ctx, &folded, ops, ctr, pt),
                )
                .unwrap();
                eprintln!(
                    "e7 j0 2^{:.1} m={m}: cols {}/{} square {:?}/{:?} ratio {:?} zero-support {}/{} single {}/{} orbit dups {} of {} ok {}/{} [{:.1}s]",
                    (inst.r as f64).log2(),
                    rep.folded.columns,
                    rep.control.columns,
                    rep.folded.square_relations,
                    rep.control.square_relations,
                    rep.square_ratio,
                    rep.folded.zero_support_rows,
                    rep.control.zero_support_rows,
                    rep.folded.single_column_rows,
                    rep.control.single_column_rows,
                    rep.orbit_duplicates,
                    rep.trials,
                    rep.folded.verified,
                    rep.control.verified,
                    started.elapsed().as_secs_f64()
                );
                rows.push(json!({
                    "experiment": "e7",
                    "family": "j0",
                    "bits": bits,
                    "seed": seed,
                    "instance": inst.name,
                    "log2_r": (inst.r as f64).log2(),
                    "r": inst.r,
                    "cofactor": inst.cofactor,
                    "seed_abscissae": size,
                    "summands": m,
                    "group_order": classes.group_order,
                    "prep_group_ops": prep_ops,
                    "stream": stream_json(&rep),
                    "wall_seconds": started.elapsed().as_secs_f64(),
                }));
            }
        }
    }
    rows
}

// ── E8: three summands on the Frobenius line, every phase priced ───

/// The measured conversion that prices an `F_p` inversion in `F_p`
/// multiplications on this host (AGENTS.md §2: foreign units convert
/// by a measured factor, recorded): the median of three wall-time
/// ratios of `inv_mod` to `mulm` over `2·10⁵` operands each.
fn inversion_in_multiplications(p: u64) -> f64 {
    use crypto_lib::cryptanalysis::glv_invariant_base::mulm;
    use crypto_lib::cryptanalysis::residual_walk::inv_mod;
    let n = 200_000u64;
    let mut ratios = Vec::new();
    for _ in 0..3 {
        let mut acc = 1u64;
        let t0 = Instant::now();
        for k in 1..=n {
            acc = mulm(acc, k % (p - 1) + 1, p);
        }
        let t_mul = t0.elapsed().as_secs_f64();
        std::hint::black_box(acc);
        let mut acc = 0u64;
        let t1 = Instant::now();
        for k in 1..=n {
            acc ^= inv_mod(k % (p - 1) + 1, p);
        }
        let t_inv = t1.elapsed().as_secs_f64();
        std::hint::black_box(acc);
        ratios.push(t_inv / t_mul.max(1e-12));
    }
    ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
    ratios[1]
}

/// E8: the three-summand pair-table oracle on the Frobenius line of a
/// generic subfield curve, both arms on one stream, with every phase
/// counted in `F_p` multiplications and inversions (the base build,
/// the pair table, the stream, the linear algebra's row operations)
/// and the matched rho — negation walk and Frobenius-folded walk —
/// counted the same way on the same instance.
fn e8(o: &Opts) -> Vec<Value> {
    let mut rows = Vec::new();
    for &p_bits in &o.bits {
        for seed in 1..=o.seeds {
            let inst = match generate_subfield_instance(p_bits, seed, false, 8) {
                Ok(i) => i,
                Err(e) => {
                    eprintln!("e8 {p_bits} seed {seed}: {e}");
                    continue;
                }
            };
            let started = Instant::now();
            let k_inv = inversion_in_multiplications(inst.p);
            take_field_counters();
            let (folded, frep) = subfield_line_base(&inst, SubfieldFold::Frobenius).unwrap();
            let (build_muls_f, build_invs_f) = take_field_counters();
            let (control, _) = subfield_line_base(&inst, SubfieldFold::Negation).unwrap();
            let (build_muls_c, build_invs_c) = take_field_counters();
            let neg = Negation { r: inst.r };
            let refs: Vec<&dyn Endomorphism<_>> = vec![&neg, &inst.frobenius];
            let classes = classes_for(&inst.curve, inst.generator, inst.r, &refs).unwrap();
            let planted = planted_for(seed, inst.r);
            let mut ops = GroupOps::default();
            let target = inst.curve.mul(&mut ops, inst.generator, planted);
            let ctx = InstanceCtx {
                group: &inst.curve,
                generator: inst.generator,
                target,
                r: inst.r,
                cofactor: inst.cofactor,
                group_order: inst.group_order,
                name: inst.name.clone(),
                field_degree: Some(3),
            };
            take_field_counters();
            // One pair table over the shared points serves both arms;
            // a pipeline would build it once per arm at the same cost,
            // so it is charged to each arm in full.
            let mut oracle = MitmOracle::new(3);
            let mut params = Params::default();
            params.set("negation_folded", "1");
            let mut table_ops = GroupOps::default();
            oracle
                .prepare(&ctx, &folded, &params, &mut table_ops)
                .unwrap();
            let (table_muls, table_invs) = take_field_counters();
            let rep = full_rank_stream(
                &inst.curve,
                inst.generator,
                target,
                inst.r,
                inst.cofactor,
                planted,
                &folded,
                &control,
                seed,
                o.max_trials,
                3,
                Some(&classes),
                |ops, ctr, pt| oracle.decompose(&ctx, &folded, ops, ctr, pt),
            )
            .unwrap();
            let (stream_muls, stream_invs) = take_field_counters();
            let mut walks = Vec::new();
            let (mut st_neg, mut st_fold, mut all_ok) = (0u64, 0u64, true);
            for k in 0..o.rho_runs as u64 {
                let ws = seed ^ (k * 0x9E37);
                take_field_counters();
                let n = rho_reference_negation(
                    &inst.curve,
                    inst.generator,
                    target,
                    inst.r,
                    ws,
                    o.rho_max_steps,
                );
                let (nm, ni) = take_field_counters();
                let f = rho_reference_folded(
                    &inst.curve,
                    inst.generator,
                    target,
                    inst.r,
                    ws,
                    o.rho_max_steps,
                    &refs,
                )
                .unwrap();
                let (fm, fi) = take_field_counters();
                all_ok &= n.verified && f.verified;
                st_neg += n.steps;
                st_fold += f.steps;
                walks.push(json!({
                    "walk_seed": ws,
                    "negation": n,
                    "negation_fp_muls": nm,
                    "negation_fp_invs": ni,
                    "folded": f,
                    "folded_fp_muls": fm,
                    "folded_fp_invs": fi,
                }));
            }
            eprintln!(
                "e8 p=2^{p_bits} r=2^{:.1} h/N1={}: cols {}/{} square {:?}/{:?} ratio {:?} rank@cols {:?}/{:?} single {}/{} dups {} of {} stream muls {} table muls {} rho steps ratio {:.2} ok {}/{}/{} [{:.1}s]",
                (inst.r as f64).log2(),
                inst.cofactor / inst.base_order,
                rep.folded.columns,
                rep.control.columns,
                rep.folded.square_relations,
                rep.control.square_relations,
                rep.square_ratio,
                rep.folded.rank_fraction_at_columns,
                rep.control.rank_fraction_at_columns,
                rep.folded.single_column_rows,
                rep.control.single_column_rows,
                rep.orbit_duplicates,
                rep.trials,
                stream_muls,
                table_muls,
                st_neg as f64 / st_fold.max(1) as f64,
                rep.folded.verified,
                rep.control.verified,
                all_ok,
                started.elapsed().as_secs_f64()
            );
            rows.push(json!({
                "experiment": "e8",
                "family": "subfield",
                "p_bits": p_bits,
                "seed": seed,
                "instance": inst.name,
                "log2_r": (inst.r as f64).log2(),
                "r": inst.r,
                "base_order": inst.base_order,
                "group_order": inst.group_order,
                "cofactor": inst.cofactor,
                "summands": 3,
                "oracle": "mitm",
                "group_order_of_fold": classes.group_order,
                "points_per_column_folded": frep.points_per_orbit,
                "unit": "F_p multiplications; an inversion is priced at inversion_in_multiplications, measured on this host; a row operation of the linear algebra is one Z/rZ multiplication, priced as one",
                "inversion_in_multiplications": k_inv,
                "base_build": {
                    "folded": {"fp_muls": build_muls_f, "fp_invs": build_invs_f},
                    "control": {"fp_muls": build_muls_c, "fp_invs": build_invs_c},
                },
                "pair_table": {"group_ops": table_ops, "fp_muls": table_muls, "fp_invs": table_invs},
                "stream_cost": {"fp_muls": stream_muls, "fp_invs": stream_invs},
                "stream": stream_json(&rep),
                "rho_runs": o.rho_runs,
                "rho_steps_ratio": st_neg as f64 / st_fold.max(1) as f64,
                "rho_expected_ratio": ((classes.group_order as f64) / 2.0).sqrt(),
                "rho_all_verified": all_ok,
                "walks": walks,
                "wall_seconds": started.elapsed().as_secs_f64(),
            }));
        }
    }
    rows
}

// ── E9: the matched folded rho on the F_{p²} and F_{p³} groups ─────

fn e9(o: &Opts) -> Vec<Value> {
    let mut rows = Vec::new();
    // Runs the two walks on one instance and returns the row's shared part.
    fn walks_on<G: CountedGroup>(
        g: &G,
        generator: G::Elt,
        target: G::Elt,
        r: u64,
        seed: u64,
        runs: usize,
        max_steps: u64,
        refs: &[&dyn Endomorphism<G>],
    ) -> (Vec<Value>, f64, f64, f64, bool) {
        let mut walks = Vec::new();
        let (mut s_neg, mut s_fold, mut st_neg, mut st_fold, mut ok) = (0.0, 0.0, 0u64, 0u64, true);
        for k in 0..runs as u64 {
            let ws = seed ^ (k * 0x9E37);
            let n = rho_reference_negation(g, generator, target, r, ws, max_steps);
            let f = rho_reference_folded(g, generator, target, r, ws, max_steps, refs).unwrap();
            ok &= n.verified && f.verified;
            s_neg += n.s;
            s_fold += f.s;
            st_neg += n.steps;
            st_fold += f.steps;
            walks.push(json!({"walk_seed": ws, "negation": n, "folded": f}));
        }
        let runs = runs.max(1) as f64;
        (
            walks,
            s_neg / runs,
            s_fold / runs,
            st_neg as f64 / st_fold.max(1) as f64,
            ok,
        )
    }
    for family in [GlsFamily::Generic, GlsFamily::J0, GlsFamily::J1728] {
        for &p_bits in &o.bits {
            for seed in 1..=o.seeds {
                let max_cofactor = if family == GlsFamily::J1728 {
                    16u64 << p_bits
                } else {
                    64
                };
                let inst = match generate_gls_instance_of(family, p_bits, seed, max_cofactor) {
                    Ok(i) => i,
                    Err(e) => {
                        eprintln!("e9 gls {family:?} {p_bits} seed {seed}: {e}");
                        continue;
                    }
                };
                let started = Instant::now();
                let neg = Negation { r: inst.r };
                let mut refs: Vec<&dyn Endomorphism<_>> = vec![&neg, &inst.psi];
                if let Some(a) = &inst.automorphism {
                    refs.push(a);
                }
                let classes = classes_for(&inst.curve, inst.generator, inst.r, &refs).unwrap();
                let planted = planted_for(seed, inst.r);
                let mut ops = GroupOps::default();
                let target = inst.curve.mul(&mut ops, inst.generator, planted);
                let (walks, s_neg, s_fold, steps_ratio, ok) = walks_on(
                    &inst.curve,
                    inst.generator,
                    target,
                    inst.r,
                    seed,
                    o.rho_runs,
                    o.rho_max_steps,
                    &refs,
                );
                let expected = (classes.group_order as f64 / 2.0).sqrt();
                eprintln!(
                    "e9 gls {} p=2^{p_bits} r=2^{:.1}: A={} S neg {:.3} folded {:.3} ratio {:.2} (steps ratio {:.2}, expected {:.2}) ok {} [{:.1}s]",
                    inst.family,
                    (inst.r as f64).log2(),
                    classes.group_order,
                    s_neg,
                    s_fold,
                    s_neg / s_fold,
                    steps_ratio,
                    expected,
                    ok,
                    started.elapsed().as_secs_f64()
                );
                rows.push(json!({
                    "experiment": "e9",
                    "family": format!("gls-{}", inst.family),
                    "p_bits": p_bits,
                    "seed": seed,
                    "instance": inst.name,
                    "log2_r": (inst.r as f64).log2(),
                    "r": inst.r,
                    "cofactor": inst.cofactor,
                    "generators": refs.iter().map(|e| e.name()).collect::<Vec<_>>(),
                    "group_order_of_fold": classes.group_order,
                    "rho_runs": o.rho_runs,
                    "negation_s_mean": s_neg,
                    "folded_s_mean": s_fold,
                    "s_ratio": s_neg / s_fold,
                    "steps_ratio": steps_ratio,
                    "expected_ratio": expected,
                    "all_verified": ok,
                    "walks": walks,
                    "wall_seconds": started.elapsed().as_secs_f64(),
                }));
            }
        }
    }
    for j0 in [false, true] {
        for &p_bits in &o.bits {
            for seed in 1..=o.seeds {
                let ratio = if j0 { 16u64 << p_bits } else { 8 };
                let inst = match generate_subfield_instance(p_bits, seed, j0, ratio) {
                    Ok(i) => i,
                    Err(e) => {
                        eprintln!("e9 subfield j0={j0} {p_bits} seed {seed}: {e}");
                        continue;
                    }
                };
                let started = Instant::now();
                let neg = Negation { r: inst.r };
                let mut refs: Vec<&dyn Endomorphism<_>> = vec![&neg, &inst.frobenius];
                if let Some(z) = &inst.zeta {
                    refs.push(z);
                }
                let classes = classes_for(&inst.curve, inst.generator, inst.r, &refs).unwrap();
                let planted = planted_for(seed, inst.r);
                let mut ops = GroupOps::default();
                let target = inst.curve.mul(&mut ops, inst.generator, planted);
                let (walks, s_neg, s_fold, steps_ratio, ok) = walks_on(
                    &inst.curve,
                    inst.generator,
                    target,
                    inst.r,
                    seed,
                    o.rho_runs,
                    o.rho_max_steps,
                    &refs,
                );
                let expected = (classes.group_order as f64 / 2.0).sqrt();
                eprintln!(
                    "e9 subfield{} p=2^{p_bits} r=2^{:.1}: A={} S neg {:.3} folded {:.3} ratio {:.2} (steps ratio {:.2}, expected {:.2}) ok {} [{:.1}s]",
                    if j0 { "-j0" } else { "" },
                    (inst.r as f64).log2(),
                    classes.group_order,
                    s_neg,
                    s_fold,
                    s_neg / s_fold,
                    steps_ratio,
                    expected,
                    ok,
                    started.elapsed().as_secs_f64()
                );
                rows.push(json!({
                    "experiment": "e9",
                    "family": if j0 { "subfield-j0" } else { "subfield" },
                    "p_bits": p_bits,
                    "seed": seed,
                    "instance": inst.name,
                    "log2_r": (inst.r as f64).log2(),
                    "r": inst.r,
                    "cofactor": inst.cofactor,
                    "generators": refs.iter().map(|e| e.name()).collect::<Vec<_>>(),
                    "group_order_of_fold": classes.group_order,
                    "rho_runs": o.rho_runs,
                    "negation_s_mean": s_neg,
                    "folded_s_mean": s_fold,
                    "s_ratio": s_neg / s_fold,
                    "steps_ratio": steps_ratio,
                    "expected_ratio": expected,
                    "all_verified": ok,
                    "walks": walks,
                    "wall_seconds": started.elapsed().as_secs_f64(),
                }));
            }
        }
    }
    rows
}

// ── E11: the algebraic S₄ oracle on the Frobenius line ─────────────

/// The S₄ line oracle against the pair table, target by target
/// (AGENTS.md §6): both asked whether `[k]G` decomposes over the base.
fn s4_agreement(
    inst: &crypto_lib::cryptanalysis::subfield_fp3::SubfieldInstance,
    fb: &crypto_lib::cryptanalysis::ic_boundary::FactorBase<
        crypto_lib::cryptanalysis::subfield_fp3::Fp3Point,
    >,
    targets: u64,
    seed: u64,
) -> Value {
    let mut ops = GroupOps::default();
    let target = inst.curve.mul(&mut ops, inst.generator, 3);
    let ctx = InstanceCtx {
        group: &inst.curve,
        generator: inst.generator,
        target,
        r: inst.r,
        cofactor: inst.cofactor,
        group_order: inst.group_order,
        name: inst.name.clone(),
        field_degree: Some(3),
    };
    let mut s4 = LineS4Oracle::new(inst.curve.f, inst.line, seed);
    s4.prepare(&ctx, fb, &Params::default(), &mut ops).unwrap();
    let mut mitm = MitmOracle::new(3);
    let mut params = Params::default();
    params.set("negation_folded", "1");
    mitm.prepare(&ctx, fb, &params, &mut ops).unwrap();
    let (mut ca, mut cb) = (OracleCounters::default(), OracleCounters::default());
    let (mut agree, mut s4_only, mut mitm_only, mut hits) = (0u64, 0u64, 0u64, 0u64);
    for k in 1..=targets {
        let pt = inst.curve.mul(&mut ops, inst.generator, k);
        let a = s4.decompose(&ctx, fb, &mut ops, &mut ca, pt);
        let b = mitm.decompose(&ctx, fb, &mut ops, &mut cb, pt);
        if let Some(idx) = &a {
            let sum = idx.iter().fold(inst.curve.identity(), |acc, &i| {
                inst.curve.add(&mut ops, acc, fb.points[i])
            });
            assert_eq!(sum, pt, "an S₄ decomposition did not sum to its target");
            hits += 1;
        }
        match (a.is_some(), b.is_some()) {
            (true, true) | (false, false) => agree += 1,
            (true, false) => s4_only += 1,
            (false, true) => mitm_only += 1,
        }
    }
    json!({
        "targets": targets,
        "agree": agree,
        "disagree": s4_only + mitm_only,
        "s4_only": s4_only,
        "mitm_only": mitm_only,
        "hits": hits,
        "unsolved": s4.stats.unsolved,
        "border_unreachable": s4.stats.border_unreachable,
    })
}

/// E11: E8's stream with the algebraic oracle in place of the pair table.
fn e11(o: &Opts) -> Vec<Value> {
    let mut rows = Vec::new();
    for &p_bits in &o.bits {
        for seed in 1..=o.seeds {
            let inst = match generate_subfield_instance(p_bits, seed, false, 8) {
                Ok(i) => i,
                Err(e) => {
                    eprintln!("e11 {p_bits} seed {seed}: {e}");
                    continue;
                }
            };
            let started = Instant::now();
            let k_inv = inversion_in_multiplications(inst.p);
            take_field_counters();
            let (folded, frep) = subfield_line_base(&inst, SubfieldFold::Frobenius).unwrap();
            let (build_muls_f, build_invs_f) = take_field_counters();
            let (control, _) = subfield_line_base(&inst, SubfieldFold::Negation).unwrap();
            let (build_muls_c, build_invs_c) = take_field_counters();
            let neg = Negation { r: inst.r };
            let refs: Vec<&dyn Endomorphism<_>> = vec![&neg, &inst.frobenius];
            let classes = classes_for(&inst.curve, inst.generator, inst.r, &refs).unwrap();
            let planted = planted_for(seed, inst.r);
            let mut ops = GroupOps::default();
            let target = inst.curve.mul(&mut ops, inst.generator, planted);
            let ctx = InstanceCtx {
                group: &inst.curve,
                generator: inst.generator,
                target,
                r: inst.r,
                cofactor: inst.cofactor,
                group_order: inst.group_order,
                name: inst.name.clone(),
                field_degree: Some(3),
            };
            let prep_started = Instant::now();
            let mut oracle = LineS4Oracle::new(inst.curve.f, inst.line, seed);
            let mut prep_ops = GroupOps::default();
            oracle
                .prepare(&ctx, &folded, &Params::default(), &mut prep_ops)
                .unwrap();
            let prep_wall = prep_started.elapsed().as_secs_f64();
            take_field_counters();
            let rep = full_rank_stream(
                &inst.curve,
                inst.generator,
                target,
                inst.r,
                inst.cofactor,
                planted,
                &folded,
                &control,
                seed,
                o.max_trials,
                3,
                Some(&classes),
                |ops, ctr, pt| oracle.decompose(&ctx, &folded, ops, ctr, pt),
            )
            .unwrap();
            let (stream_muls, stream_invs) = take_field_counters();
            let st = &oracle.stats;
            let totals = oracle.solver_totals().unwrap();
            let agreement = if p_bits <= 9 {
                Some(s4_agreement(&inst, &folded, 300, seed))
            } else {
                None
            };
            let mut walks = Vec::new();
            let (mut st_neg, mut st_fold, mut all_ok) = (0u64, 0u64, true);
            for k in 0..o.rho_runs as u64 {
                let ws = seed ^ (k * 0x9E37);
                take_field_counters();
                let n = rho_reference_negation(
                    &inst.curve,
                    inst.generator,
                    target,
                    inst.r,
                    ws,
                    o.rho_max_steps,
                );
                let (nm, ni) = take_field_counters();
                let f = rho_reference_folded(
                    &inst.curve,
                    inst.generator,
                    target,
                    inst.r,
                    ws,
                    o.rho_max_steps,
                    &refs,
                )
                .unwrap();
                let (fm, fi) = take_field_counters();
                all_ok &= n.verified && f.verified;
                st_neg += n.steps;
                st_fold += f.steps;
                walks.push(json!({
                    "walk_seed": ws,
                    "negation": n,
                    "negation_fp_muls": nm,
                    "negation_fp_invs": ni,
                    "folded": f,
                    "folded_fp_muls": fm,
                    "folded_fp_invs": fi,
                }));
            }
            eprintln!(
                "e11 p=2^{p_bits} r=2^{:.1} h/N1={}: cols {}/{} square {:?}/{:?} ratio {:?} rank@cols {:?}/{:?} trials {} solves {} unsolved {} muls/solve {:.0} stream muls {} agreement {} ok {}/{}/{} [{:.1}s]",
                (inst.r as f64).log2(),
                inst.cofactor / inst.base_order,
                rep.folded.columns,
                rep.control.columns,
                rep.folded.square_relations,
                rep.control.square_relations,
                rep.square_ratio,
                rep.folded.rank_fraction_at_columns,
                rep.control.rank_fraction_at_columns,
                rep.trials,
                st.solves,
                st.unsolved,
                st.fp_muls as f64 / st.solves.max(1) as f64,
                stream_muls,
                agreement
                    .as_ref()
                    .map(|a| format!("{}/{}", a["agree"], a["targets"]))
                    .unwrap_or_else(|| "—".into()),
                rep.folded.verified,
                rep.control.verified,
                all_ok,
                started.elapsed().as_secs_f64()
            );
            rows.push(json!({
                "experiment": "e11",
                "family": "subfield",
                "p_bits": p_bits,
                "seed": seed,
                "instance": inst.name,
                "log2_r": (inst.r as f64).log2(),
                "r": inst.r,
                "base_order": inst.base_order,
                "group_order": inst.group_order,
                "cofactor": inst.cofactor,
                "summands": 3,
                "oracle": "line-s4-macaulay",
                "group_order_of_fold": classes.group_order,
                "points_per_column_folded": frep.points_per_orbit,
                "unit": "F_p multiplications; an inversion is priced at inversion_in_multiplications, measured on this host; the solver counts its own multiplications; a row operation of the linear algebra is one Z/rZ multiplication, priced as one",
                "inversion_in_multiplications": k_inv,
                "base_build": {
                    "folded": {"fp_muls": build_muls_f, "fp_invs": build_invs_f},
                    "control": {"fp_muls": build_muls_c, "fp_invs": build_invs_c},
                },
                "oracle_setup": {"wall_seconds": prep_wall},
                "solver": {
                    "calls": st.solves,
                    "fp_muls": st.fp_muls,
                    "macaulay_muls": st.macaulay_muls,
                    "unsolved": st.unsolved,
                    "border_unreachable": st.border_unreachable,
                    "retried_at_degree_11": st.retried_at_degree_11,
                    "retried_at_degree_12": st.retried_at_degree_12,
                    "retried_at_degree_13": st.retried_at_degree_13,
                    "e_solutions": st.e_solutions,
                    "split_cubics": st.split_cubics,
                    "root_triples": oracle.solutions,
                    "unliftable_systems": oracle.unliftable,
                    "wall_ns": totals.wall_ns,
                },
                "stream_cost": {"fp_muls": stream_muls, "fp_invs": stream_invs},
                "stream": stream_json(&rep),
                "agreement_with_mitm": agreement,
                "rho_runs": o.rho_runs,
                "rho_steps_ratio": st_neg as f64 / st_fold.max(1) as f64,
                "rho_expected_ratio": ((classes.group_order as f64) / 2.0).sqrt(),
                "rho_all_verified": all_ok,
                "walks": walks,
                "wall_seconds": started.elapsed().as_secs_f64(),
            }));
        }
    }
    rows
}

// ── E12: the pair table over orbit representatives ─────────────────

/// One stream on a subfield line with a given three-summand oracle,
/// field multiplications counted around the table build and the stream.
#[allow(clippy::too_many_arguments)]
fn e12_subfield_arm(
    inst: &crypto_lib::cryptanalysis::subfield_fp3::SubfieldInstance,
    ctx: &InstanceCtx<'_, crypto_lib::cryptanalysis::subfield_fp3::Fp3Curve>,
    folded: &crypto_lib::cryptanalysis::ic_boundary::FactorBase<
        crypto_lib::cryptanalysis::subfield_fp3::Fp3Point,
    >,
    control: &crypto_lib::cryptanalysis::ic_boundary::FactorBase<
        crypto_lib::cryptanalysis::subfield_fp3::Fp3Point,
    >,
    classes: &EndomorphismClasses<'_, crypto_lib::cryptanalysis::subfield_fp3::Fp3Curve>,
    planted: u64,
    seed: u64,
    max_trials: u64,
    oracle: &mut dyn DecompositionOracle<crypto_lib::cryptanalysis::subfield_fp3::Fp3Curve>,
) -> (Value, Value, StreamReport) {
    take_field_counters();
    let mut table_ops = GroupOps::default();
    let mut params = Params::default();
    params.set("negation_folded", "1");
    let started = Instant::now();
    oracle
        .prepare(ctx, folded, &params, &mut table_ops)
        .unwrap();
    let table_wall = started.elapsed().as_secs_f64();
    let (tm, ti) = take_field_counters();
    let rep = full_rank_stream(
        &inst.curve,
        inst.generator,
        ctx.target,
        inst.r,
        inst.cofactor,
        planted,
        folded,
        control,
        seed,
        max_trials,
        3,
        Some(classes),
        |ops, ctr, pt| oracle.decompose(ctx, folded, ops, ctr, pt),
    )
    .unwrap();
    let (sm, si) = take_field_counters();
    (
        json!({"group_ops": table_ops, "fp_muls": tm, "fp_invs": ti, "wall_seconds": table_wall}),
        json!({"fp_muls": sm, "fp_invs": si}),
        rep,
    )
}

/// The negation table and the orbit table on one target set, target by
/// target (AGENTS.md §6).
fn e12_agreement<G: CountedGroup>(
    ctx: &InstanceCtx<'_, G>,
    fb: &crypto_lib::cryptanalysis::ic_boundary::FactorBase<G::Elt>,
    a: &mut dyn DecompositionOracle<G>,
    b: &mut dyn DecompositionOracle<G>,
    targets: u64,
) -> Value {
    let mut ops = GroupOps::default();
    let (mut ca, mut cb) = (OracleCounters::default(), OracleCounters::default());
    let (mut agree, mut a_only, mut b_only, mut hits) = (0u64, 0u64, 0u64, 0u64);
    for k in 2..targets + 2 {
        let pt = ctx.group.mul(&mut ops, ctx.generator, k);
        let x = a.decompose(ctx, fb, &mut ops, &mut ca, pt);
        let y = b.decompose(ctx, fb, &mut ops, &mut cb, pt);
        if let Some(idx) = &x {
            let sum = idx.iter().fold(ctx.group.identity(), |acc, &i| {
                ctx.group.add(&mut ops, acc, fb.points[i])
            });
            assert_eq!(
                sum, pt,
                "an orbit-table decomposition did not sum to its target"
            );
            hits += 1;
        }
        match (x.is_some(), y.is_some()) {
            (true, true) | (false, false) => agree += 1,
            (true, false) => a_only += 1,
            (false, true) => b_only += 1,
        }
    }
    json!({"targets": targets, "agree": agree, "orbit_only": a_only, "negation_only": b_only, "hits": hits})
}

/// E12 on the subfield line, three summands: E8's stream once over the
/// negation table and once over the orbit table, every phase counted.
fn e12(o: &Opts) -> Vec<Value> {
    let mut rows = Vec::new();
    for &p_bits in &o.bits {
        for seed in 1..=o.seeds {
            let inst = match generate_subfield_instance(p_bits, seed, false, 8) {
                Ok(i) => i,
                Err(e) => {
                    eprintln!("e12 {p_bits} seed {seed}: {e}");
                    continue;
                }
            };
            let started = Instant::now();
            let k_inv = inversion_in_multiplications(inst.p);
            take_field_counters();
            let (folded, frep) = subfield_line_base(&inst, SubfieldFold::Frobenius).unwrap();
            let (build_muls_f, build_invs_f) = take_field_counters();
            let (control, _) = subfield_line_base(&inst, SubfieldFold::Negation).unwrap();
            let (build_muls_c, build_invs_c) = take_field_counters();
            let neg = Negation { r: inst.r };
            let refs: Vec<&dyn Endomorphism<_>> = vec![&neg, &inst.frobenius];
            let classes = classes_for(&inst.curve, inst.generator, inst.r, &refs).unwrap();
            let planted = planted_for(seed, inst.r);
            let mut ops = GroupOps::default();
            let target = inst.curve.mul(&mut ops, inst.generator, planted);
            let ctx = InstanceCtx {
                group: &inst.curve,
                generator: inst.generator,
                target,
                r: inst.r,
                cofactor: inst.cofactor,
                group_order: inst.group_order,
                name: inst.name.clone(),
                field_degree: Some(3),
            };
            let mut neg_oracle = MitmOracle::new(3);
            let (neg_table, neg_stream, neg_rep) = e12_subfield_arm(
                &inst,
                &ctx,
                &folded,
                &control,
                &classes,
                planted,
                seed,
                o.max_trials,
                &mut neg_oracle,
            );
            let orbit_classes =
                EndomorphismClasses::new(&inst.curve, inst.generator, inst.r, refs.clone())
                    .unwrap();
            let mut orbit_oracle = OrbitMitmOracle::new(3, SubfieldOrbitKey, orbit_classes);
            let (orbit_table, orbit_stream, orbit_rep) = e12_subfield_arm(
                &inst,
                &ctx,
                &folded,
                &control,
                &classes,
                planted,
                seed,
                o.max_trials,
                &mut orbit_oracle,
            );
            let entries = orbit_oracle
                .table()
                .map(|t| (t.entries, t.representatives))
                .unwrap();
            let agreement = if p_bits <= 9 {
                let ac =
                    EndomorphismClasses::new(&inst.curve, inst.generator, inst.r, refs.clone())
                        .unwrap();
                let mut a = OrbitMitmOracle::new(3, SubfieldOrbitKey, ac);
                let mut b = MitmOracle::new(3);
                let mut params = Params::default();
                params.set("negation_folded", "1");
                let mut tmp = GroupOps::default();
                a.prepare(&ctx, &folded, &params, &mut tmp).unwrap();
                b.prepare(&ctx, &folded, &params, &mut tmp).unwrap();
                Some(e12_agreement(&ctx, &folded, &mut a, &mut b, 200))
            } else {
                None
            };
            let mut walks = Vec::new();
            let (mut st_neg, mut st_fold, mut all_ok) = (0u64, 0u64, true);
            for k in 0..o.rho_runs as u64 {
                let ws = seed ^ (k * 0x9E37);
                take_field_counters();
                let n = rho_reference_negation(
                    &inst.curve,
                    inst.generator,
                    target,
                    inst.r,
                    ws,
                    o.rho_max_steps,
                );
                let (nm, ni) = take_field_counters();
                let f = rho_reference_folded(
                    &inst.curve,
                    inst.generator,
                    target,
                    inst.r,
                    ws,
                    o.rho_max_steps,
                    &refs,
                )
                .unwrap();
                let (fm, fi) = take_field_counters();
                all_ok &= n.verified && f.verified;
                st_neg += n.steps;
                st_fold += f.steps;
                walks.push(json!({
                    "walk_seed": ws,
                    "negation": n,
                    "negation_fp_muls": nm,
                    "negation_fp_invs": ni,
                    "folded": f,
                    "folded_fp_muls": fm,
                    "folded_fp_invs": fi,
                }));
            }
            eprintln!(
                "e12 p=2^{p_bits} r=2^{:.1}: entries orbit {} | table muls orbit {} neg {} | stream muls orbit {} neg {} | square {:?}/{:?} vs {:?}/{:?} | collisions {} mismatches {} | ok {}/{}/{}/{}/{} [{:.1}s]",
                (inst.r as f64).log2(),
                entries.0,
                orbit_table["fp_muls"],
                neg_table["fp_muls"],
                orbit_stream["fp_muls"],
                neg_stream["fp_muls"],
                orbit_rep.folded.square_relations,
                orbit_rep.control.square_relations,
                neg_rep.folded.square_relations,
                neg_rep.control.square_relations,
                orbit_oracle.stats.key_collisions,
                orbit_oracle.stats.image_mismatches,
                orbit_rep.folded.verified,
                orbit_rep.control.verified,
                neg_rep.folded.verified,
                neg_rep.control.verified,
                all_ok,
                started.elapsed().as_secs_f64()
            );
            rows.push(json!({
                "experiment": "e12",
                "family": "subfield",
                "p_bits": p_bits,
                "seed": seed,
                "instance": inst.name,
                "log2_r": (inst.r as f64).log2(),
                "r": inst.r,
                "base_order": inst.base_order,
                "group_order": inst.group_order,
                "cofactor": inst.cofactor,
                "summands": 3,
                "group_order_of_fold": classes.group_order,
                "points_per_column_folded": frep.points_per_orbit,
                "points": folded.points.len(),
                "unit": "F_p multiplications; an inversion is priced at inversion_in_multiplications, measured on this host; the orbit key's two uncounted multiplications per call are added; a row operation of the linear algebra is one Z/rZ multiplication, priced as one",
                "inversion_in_multiplications": k_inv,
                "base_build": {
                    "folded": {"fp_muls": build_muls_f, "fp_invs": build_invs_f},
                    "control": {"fp_muls": build_muls_c, "fp_invs": build_invs_c},
                },
                "negation_table": {"table": neg_table, "stream_cost": neg_stream, "stream": stream_json(&neg_rep)},
                "orbit_table": {
                    "table": orbit_table,
                    "entries": entries.0,
                    "representatives": entries.1,
                    "build_uncounted_muls": orbit_oracle.build_uncounted_muls(),
                    "stream_uncounted_muls": orbit_oracle.uncounted_muls() - orbit_oracle.build_uncounted_muls(),
                    "probe_stats": orbit_oracle.stats,
                    "stream_cost": orbit_stream,
                    "stream": stream_json(&orbit_rep),
                },
                "agreement": agreement,
                "rho_runs": o.rho_runs,
                "rho_steps_ratio": st_neg as f64 / st_fold.max(1) as f64,
                "rho_expected_ratio": ((classes.group_order as f64) / 2.0).sqrt(),
                "rho_all_verified": all_ok,
                "walks": walks,
                "wall_seconds": started.elapsed().as_secs_f64(),
            }));
        }
    }
    rows
}

/// Measured conversions on a prime-field curve: an `F_p` multiplication
/// and a square root, each in affine additions (median of three timed
/// batches each).  AGENTS.md §2: a foreign unit converts at a measured
/// factor, and the factor is recorded with the row.
fn prime_conversions(inst: &PrimeInstance) -> (f64, f64) {
    use crypto_lib::cryptanalysis::glv_invariant_base::mulm;
    let curve = &inst.curve;
    let g = inst.generator_point();
    let p = curve.p;
    let med = |f: &mut dyn FnMut() -> f64| {
        let mut v: Vec<f64> = (0..3).map(|_| f()).collect();
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v[1]
    };
    let n = 100_000u64;
    let add_ns = med(&mut || {
        let mut ops = GroupOps::default();
        let mut acc = g;
        let t = Instant::now();
        for _ in 0..n {
            acc = curve.add(&mut ops, acc, g);
        }
        std::hint::black_box(acc);
        t.elapsed().as_nanos() as f64 / n as f64
    });
    let mul_ns = med(&mut || {
        let mut acc = 3u64;
        let t = Instant::now();
        for k in 1..=n {
            acc = mulm(acc, k % (p - 1) + 1, p);
        }
        std::hint::black_box(acc);
        t.elapsed().as_nanos() as f64 / n as f64
    });
    let sqrt_ns = med(&mut || {
        let mut acc = 0u64;
        let t = Instant::now();
        for k in 1..=(n / 10) {
            acc ^= curve.sqrt(k % (p - 1) + 1).unwrap_or(0);
        }
        std::hint::black_box(acc);
        t.elapsed().as_nanos() as f64 / (n / 10) as f64
    });
    (mul_ns / add_ns, sqrt_ns / add_ns)
}

/// E12 on F_p, two summands: the regime where the pair table is the
/// cost (§5.2).  Unit: group additions; key multiplications, row
/// operations and square roots convert at measured factors.
fn e12p(o: &Opts) -> Vec<Value> {
    let mut rows = Vec::new();
    for (family, exponent) in [(CmFamily::J0, 3u32), (CmFamily::J1728, 2)] {
        for &bits in &o.bits {
            for seed in 1..=o.seeds {
                let inst = match generate_cm_instance(family, bits, seed, 8) {
                    Ok(i) => i,
                    Err(e) => {
                        eprintln!("e12p {} {bits} seed {seed}: {e}", family.name());
                        continue;
                    }
                };
                let started = Instant::now();
                let (mul_in_adds, sqrt_in_adds) = prime_conversions(&inst);
                let size = base_size(inst.r);
                let (folded, frep) =
                    glv_orbit_base(&inst, size, AutomorphismGroup::Auto, true).unwrap();
                let (control, _) =
                    glv_orbit_base(&inst, size, AutomorphismGroup::Auto, false).unwrap();
                let gens = automorphism_generators(&inst, AutomorphismGroup::Auto).unwrap();
                let refs: Vec<&dyn Endomorphism<_>> = gens.iter().map(|b| b.as_ref()).collect();
                let planted = planted_for(seed, inst.r);
                let (ctx, target) = prime_ctx(&inst, planted);
                let g = ctx.generator;
                let run = |oracle: &mut dyn DecompositionOracle<_>| {
                    let mut params = Params::default();
                    params.set("negation_folded", "1");
                    let mut prep_ops = GroupOps::default();
                    oracle
                        .prepare(&ctx, &folded, &params, &mut prep_ops)
                        .unwrap();
                    let rep = full_rank_stream(
                        &inst.curve,
                        g,
                        target,
                        inst.r,
                        inst.cofactor,
                        planted,
                        &folded,
                        &control,
                        seed,
                        o.max_trials,
                        2,
                        None,
                        |ops, ctr, pt| oracle.decompose(&ctx, &folded, ops, ctr, pt),
                    )
                    .unwrap();
                    (prep_ops, rep)
                };
                let mut neg_oracle = MitmOracle::new(2);
                let (neg_prep, neg_rep) = run(&mut neg_oracle);
                let classes =
                    EndomorphismClasses::new(&inst.curve, g, inst.r, refs.clone()).unwrap();
                let w = classes.group_order;
                let mut orbit_oracle = OrbitMitmOracle::new(2, PrimePowerKey { exponent }, classes);
                let (orbit_prep, orbit_rep) = run(&mut orbit_oracle);
                let entries = orbit_oracle
                    .table()
                    .map(|t| (t.entries, t.representatives, t.build_keys))
                    .unwrap();
                let agreement = if bits <= 20 {
                    let ac =
                        EndomorphismClasses::new(&inst.curve, g, inst.r, refs.clone()).unwrap();
                    let mut a = OrbitMitmOracle::new(2, PrimePowerKey { exponent }, ac);
                    let mut b = MitmOracle::new(2);
                    let mut params = Params::default();
                    params.set("negation_folded", "1");
                    let mut tmp = GroupOps::default();
                    a.prepare(&ctx, &folded, &params, &mut tmp).unwrap();
                    b.prepare(&ctx, &folded, &params, &mut tmp).unwrap();
                    Some(e12_agreement(&ctx, &folded, &mut a, &mut b, 3000))
                } else {
                    None
                };
                let mut walks = Vec::new();
                let (mut st_neg, mut st_fold, mut all_ok) = (0u64, 0u64, true);
                for k in 0..o.rho_runs as u64 {
                    let ws = seed ^ (k * 0x9E37);
                    let n =
                        rho_reference_negation(&inst.curve, g, target, inst.r, ws, o.rho_max_steps);
                    let f = rho_reference_folded(
                        &inst.curve,
                        g,
                        target,
                        inst.r,
                        ws,
                        o.rho_max_steps,
                        &refs,
                    )
                    .unwrap();
                    all_ok &= n.verified && f.verified;
                    st_neg += n.steps;
                    st_fold += f.steps;
                    walks.push(json!({"walk_seed": ws, "negation": n, "folded": f}));
                }
                eprintln!(
                    "e12p {} 2^{:.1}: points {} w {} | table adds orbit {} neg {} | entries orbit {} | stream adds orbit {} neg {} | square {:?}/{:?} vs {:?}/{:?} | collisions {} mismatches {} | mul/add {:.3} sqrt/add {:.2} | ok {}/{}/{}/{}/{} [{:.1}s]",
                    family.name(),
                    (inst.r as f64).log2(),
                    folded.points.len(),
                    w,
                    orbit_prep.adds,
                    neg_prep.adds,
                    entries.0,
                    orbit_rep.group_ops.adds,
                    neg_rep.group_ops.adds,
                    orbit_rep.folded.square_relations,
                    orbit_rep.control.square_relations,
                    neg_rep.folded.square_relations,
                    neg_rep.control.square_relations,
                    orbit_oracle.stats.key_collisions,
                    orbit_oracle.stats.image_mismatches,
                    mul_in_adds,
                    sqrt_in_adds,
                    orbit_rep.folded.verified,
                    orbit_rep.control.verified,
                    neg_rep.folded.verified,
                    neg_rep.control.verified,
                    all_ok,
                    started.elapsed().as_secs_f64()
                );
                rows.push(json!({
                    "experiment": "e12p",
                    "family": family.name(),
                    "bits": bits,
                    "seed": seed,
                    "instance": inst.name,
                    "log2_r": (inst.r as f64).log2(),
                    "r": inst.r,
                    "cofactor": inst.cofactor,
                    "summands": 2,
                    "group_order_of_fold": w,
                    "seed_abscissae": size,
                    "points": folded.points.len(),
                    "points_per_column_folded": frep.points_per_orbit,
                    "unit": "group additions and doublings; F_p multiplications (orbit keys, generator maps on hits, linear-algebra row operations) at mul_in_adds and square roots (base build) at sqrt_in_adds, both measured on this host",
                    "mul_in_adds": mul_in_adds,
                    "sqrt_in_adds": sqrt_in_adds,
                    "base_build": {
                        "folded": {"group_ops": folded.cost.group_ops, "sqrt_solves": folded.cost.get("sqrt_solves"), "endomorphism_maps": folded.cost.get("endomorphism_maps")},
                        "control": {"group_ops": control.cost.group_ops, "sqrt_solves": control.cost.get("sqrt_solves"), "endomorphism_maps": control.cost.get("endomorphism_maps")},
                    },
                    "negation_table": {"table_group_ops": neg_prep, "stream": stream_json(&neg_rep)},
                    "orbit_table": {
                        "table_group_ops": orbit_prep,
                        "entries": entries.0,
                        "representatives": entries.1,
                        "build_keys": entries.2,
                        "build_uncounted_muls": orbit_oracle.build_uncounted_muls(),
                        "stream_uncounted_muls": orbit_oracle.uncounted_muls() - orbit_oracle.build_uncounted_muls(),
                        "probe_stats": orbit_oracle.stats,
                        "stream": stream_json(&orbit_rep),
                    },
                    "agreement": agreement,
                    "rho_runs": o.rho_runs,
                    "rho_steps_ratio": st_neg as f64 / st_fold.max(1) as f64,
                    "rho_expected_ratio": ((w as f64) / 2.0).sqrt(),
                    "rho_all_verified": all_ok,
                    "walks": walks,
                    "wall_seconds": started.elapsed().as_secs_f64(),
                }));
            }
        }
    }
    rows
}

// ── E13: the 2-torsion symmetry of the system with the fold ────────

/// E13: on the `Y`-line of a subfield curve with a rational 2-torsion
/// point, the `D₃`-symmetrised oracle (16 solutions a system) on three
/// arms — fold `⟨−1, π, τ_T⟩` (12 a column), fold `⟨−1, π⟩` (6), negation
/// (2) — and, at `p ≤ 2^10`, the `S₃`-symmetrised Macaulay oracle (64) on
/// the two folded bases and the same targets, every phase counted, beside
/// the matched rho.
fn e13(o: &Opts) -> Vec<Value> {
    let mut rows = Vec::new();
    for &p_bits in &o.bits {
        for seed in 1..=o.seeds {
            let inst = match generate_fghr_instance(p_bits, seed, 8) {
                Ok(i) => i,
                Err(e) => {
                    eprintln!("e13 {p_bits} seed {seed}: {e}");
                    continue;
                }
            };
            let sub = &inst.sub;
            let started = Instant::now();
            let k_inv = inversion_in_multiplications(sub.p);
            take_field_counters();
            let (fold12, r12) = fghr_line_base(&inst, FghrFold::FrobeniusTranslation).unwrap();
            let (b12m, b12i) = take_field_counters();
            let (fold6, r6) = fghr_line_base(&inst, FghrFold::Frobenius).unwrap();
            let (b6m, b6i) = take_field_counters();
            let (control, _) = fghr_line_base(&inst, FghrFold::Negation).unwrap();
            let (b2m, b2i) = take_field_counters();
            let polys = fghr_polynomials(&inst).unwrap();
            let neg = Negation { r: sub.r };
            let refs: Vec<&dyn Endomorphism<_>> = vec![&neg, &sub.frobenius];
            let classes = classes_for(&sub.curve, sub.generator, sub.r, &refs).unwrap();
            let planted = planted_for(seed, sub.r);
            let mut ops = GroupOps::default();
            let target = sub.curve.mul(&mut ops, sub.generator, planted);
            let ctx = InstanceCtx {
                group: &sub.curve,
                generator: sub.generator,
                target,
                r: sub.r,
                cofactor: sub.cofactor,
                group_order: sub.group_order,
                name: sub.name.clone(),
                field_degree: Some(3),
            };
            let stream = |folded: &crypto_lib::cryptanalysis::ic_boundary::FactorBase<
                crypto_lib::cryptanalysis::subfield_fp3::Fp3Point,
            >,
                          oracle: &mut dyn DecompositionOracle<
                crypto_lib::cryptanalysis::subfield_fp3::Fp3Curve,
            >| {
                take_field_counters();
                let rep = full_rank_stream_until(
                    &sub.curve,
                    sub.generator,
                    target,
                    sub.r,
                    sub.cofactor,
                    planted,
                    folded,
                    &control,
                    seed,
                    o.max_trials,
                    3,
                    Some(&classes),
                    StopRule::FoldedSquareBothPinned,
                    |ops, ctr, pt| oracle.decompose(&ctx, folded, ops, ctr, pt),
                )
                .unwrap();
                let (m, i) = take_field_counters();
                (rep, m, i)
            };
            let mut d3a = FghrOracle::new(&inst, &polys, seed);
            let (rep_a, ma, ia) = stream(&fold12, &mut d3a);
            let mut d3b = FghrOracle::new(&inst, &polys, seed);
            let (rep_b, mb, ib) = stream(&fold6, &mut d3b);
            // Both bases under the S₃ oracle too, so that the combined
            // ratio (fold 6 with S₃ against fold 12 with D₃) is measured,
            // not composed from the two separate ones.
            let s3_run = |base| {
                (p_bits <= 10).then(|| {
                    let mut s3 = YLineS4Oracle::new(&inst, &polys, seed);
                    let (rep, m, i) = stream(base, &mut s3);
                    json!({
                        "stream": stream_json(&rep),
                        "stream_cost": {"fp_muls": m, "fp_invs": i},
                        "solver": {"calls": s3.stats.solves, "fp_muls": s3.stats.fp_muls, "unsolved": s3.stats.unsolved, "quotient_dim_total": s3.stats.quotient_dim_total, "unliftable": s3.unliftable},
                    })
                })
            };
            let s3_stream = s3_run(&fold12);
            let s3_fold6 = s3_run(&fold6);
            let agreement = if p_bits <= 9 {
                let mut d3 = FghrOracle::new(&inst, &polys, seed ^ 5);
                let mut s3 = YLineS4Oracle::new(&inst, &polys, seed ^ 5);
                let mut mitm = MitmOracle::new(3);
                let mut params = Params::default();
                params.set("negation_folded", "1");
                let mut tmp = GroupOps::default();
                mitm.prepare(&ctx, &fold12, &params, &mut tmp).unwrap();
                let mut c = [
                    OracleCounters::default(),
                    OracleCounters::default(),
                    OracleCounters::default(),
                ];
                let (mut hits, mut d3_dis, mut s3_dis) = (0u64, 0u64, 0u64);
                let n = 200u64;
                for k in 2..n + 2 {
                    let pt = sub.curve.mul(&mut tmp, sub.generator, k);
                    let a = d3.decompose(&ctx, &fold12, &mut tmp, &mut c[0], pt);
                    let b = s3.decompose(&ctx, &fold12, &mut tmp, &mut c[1], pt);
                    let m = mitm.decompose(&ctx, &fold12, &mut tmp, &mut c[2], pt);
                    for idx in [&a, &b].into_iter().flatten() {
                        let sum = idx.iter().fold(sub.curve.identity(), |acc, &i| {
                            sub.curve.add(&mut tmp, acc, fold12.points[i])
                        });
                        assert_eq!(
                            sum, pt,
                            "an algebraic decomposition did not sum to its target"
                        );
                    }
                    hits += m.is_some() as u64;
                    d3_dis += (a.is_some() != m.is_some()) as u64;
                    s3_dis += (b.is_some() != m.is_some()) as u64;
                }
                Some(
                    json!({"targets": n, "pair_table_hits": hits, "d3_disagreements": d3_dis, "s3_disagreements": s3_dis, "d3_muls_per_call": d3.stats.fp_muls as f64 / n as f64, "s3_muls_per_call": s3.stats.fp_muls as f64 / n as f64, "s3_quotient_per_call": s3.stats.quotient_dim_total as f64 / n as f64}),
                )
            } else {
                None
            };
            let mut walks = Vec::new();
            let (mut st_neg, mut st_fold, mut all_ok) = (0u64, 0u64, true);
            for k in 0..o.rho_runs as u64 {
                let ws = seed ^ (k * 0x9E37);
                take_field_counters();
                let n = rho_reference_negation(
                    &sub.curve,
                    sub.generator,
                    target,
                    sub.r,
                    ws,
                    o.rho_max_steps,
                );
                let (nm, ni) = take_field_counters();
                let f = rho_reference_folded(
                    &sub.curve,
                    sub.generator,
                    target,
                    sub.r,
                    ws,
                    o.rho_max_steps,
                    &refs,
                )
                .unwrap();
                let (fm, fi) = take_field_counters();
                all_ok &= n.verified && f.verified;
                st_neg += n.steps;
                st_fold += f.steps;
                walks.push(json!({
                    "walk_seed": ws,
                    "negation": n,
                    "negation_fp_muls": nm,
                    "negation_fp_invs": ni,
                    "folded": f,
                    "folded_fp_muls": fm,
                    "folded_fp_invs": fi,
                }));
            }
            eprintln!(
                "e13 p=2^{p_bits} r=2^{:.1} h={}: cols {}/{}/{} | square fold12 {:?}/{:?} fold6 {:?}/{:?} | D3 muls/call {:.0} deg max {} unsolved {} | S3 {} | agreement {} | ok {}/{}/{} [{:.1}s]",
                (sub.r as f64).log2(),
                sub.cofactor,
                fold12.columns,
                fold6.columns,
                control.columns,
                rep_a.folded.square_relations,
                rep_a.control.square_relations,
                rep_b.folded.square_relations,
                rep_b.control.square_relations,
                d3a.stats.fp_muls as f64 / d3a.stats.calls.max(1) as f64,
                d3a.stats.resultant_degree_max,
                d3a.stats.unsolved,
                s3_stream
                    .as_ref()
                    .map(|v| format!("muls/call {:.0}", v["solver"]["fp_muls"].as_f64().unwrap() / v["solver"]["calls"].as_f64().unwrap().max(1.0)))
                    .unwrap_or_else(|| "—".into()),
                agreement
                    .as_ref()
                    .map(|a| format!("d3 {}/{} s3 {}/{}", a["d3_disagreements"], a["targets"], a["s3_disagreements"], a["targets"]))
                    .unwrap_or_else(|| "—".into()),
                rep_a.folded.verified && rep_a.control.verified,
                rep_b.folded.verified,
                all_ok,
                started.elapsed().as_secs_f64()
            );
            rows.push(json!({
                "experiment": "e13",
                "family": "subfield-2torsion",
                "p_bits": p_bits,
                "seed": seed,
                "instance": sub.name,
                "log2_r": (sub.r as f64).log2(),
                "r": sub.r,
                "base_order": sub.base_order,
                "group_order": sub.group_order,
                "cofactor": sub.cofactor,
                "x0": inst.x0,
                "c": inst.c,
                "summands": 3,
                "unit": "F_p multiplications; an inversion is priced at inversion_in_multiplications, measured on this host; each solver counts its own multiplications; a row operation of the linear algebra is one Z/rZ multiplication, priced as one; the polynomials' once-per-curve set-up is charged to every arm",
                "inversion_in_multiplications": k_inv,
                "points": fold12.points.len(),
                "points_per_column": {"fold12": r12.points_per_orbit, "fold6": r6.points_per_orbit, "control": 2.0},
                "columns": {"fold12": fold12.columns, "fold6": fold6.columns, "control": control.columns},
                "base_build": {
                    "fold12": {"fp_muls": b12m, "fp_invs": b12i},
                    "fold6": {"fp_muls": b6m, "fp_invs": b6i},
                    "control": {"fp_muls": b2m, "fp_invs": b2i},
                },
                "polynomials": {"y_terms": polys.y_terms, "s3_terms": polys.s3_terms, "d3_terms": polys.d3_terms, "setup_fp_muls": polys.setup_muls},
                "d3_fold12": {"stream": stream_json(&rep_a), "stream_cost": {"fp_muls": ma, "fp_invs": ia}, "solver": d3a.stats},
                "d3_fold6": {"stream": stream_json(&rep_b), "stream_cost": {"fp_muls": mb, "fp_invs": ib}, "solver": d3b.stats},
                "s3_fold12": s3_stream,
                "s3_fold6": s3_fold6,
                "agreement": agreement,
                "rho_runs": o.rho_runs,
                "rho_steps_ratio": st_neg as f64 / st_fold.max(1) as f64,
                "rho_expected_ratio": ((classes.group_order as f64) / 2.0).sqrt(),
                "rho_all_verified": all_ok,
                "walks": walks,
                "wall_seconds": started.elapsed().as_secs_f64(),
            }));
        }
    }
    rows
}

// ── E15: does the §6.5 degeneracy transfer to the ECC2K-130 family? ──

/// The prime factorisation of a small cofactor, by trial division.
fn small_factors(mut h: u64) -> Vec<u64> {
    let mut out = Vec::new();
    let mut d = 2u64;
    while d * d <= h {
        while h.is_multiple_of(d) {
            out.push(d);
            h /= d;
        }
        d += 1;
    }
    if h > 1 {
        out.push(h);
    }
    out
}

/// The multiplicative order of 2 modulo an odd `n`.
fn ord2(n: u32) -> u32 {
    let (mut x, mut k) = (2 % n, 1);
    while x != 1 {
        x = x * 2 % n;
        k += 1;
    }
    k
}

/// E15: two-summand streams on `E_0: y² + xy = x³ + 1` (the ECC2K-130
/// family) over `GF(2^n)`, the signed-Frobenius-orbit base against the
/// abscissa control, as E2 ran them on `E_1`.  The question is the §6.5
/// one: `E_0(F_2) ≅ Z/4` is cofactor and Frobenius-fixed, as `E(F_p)` is
/// on a subfield curve — does the block degeneracy follow?  Each row
/// records the cofactor's factorisation, whether the cofactor is exactly
/// `E_0(F_2)` (the challenge's shape), `ord_n(2)` and the intermediate
/// subfields (AGENTS.md §8b), the structural deficiency `D` and the
/// single-column rows of both arms.
fn e15(o: &Opts) -> Vec<Value> {
    let mut rows = Vec::new();
    let degrees: Vec<u32> = if o.bits.is_empty() {
        (13..=33).collect()
    } else {
        o.bits.clone()
    };
    for n in degrees {
        let Some(inst) = koblitz_instance(0, n) else {
            eprintln!("e15 n={n}: no instance");
            continue;
        };
        let Some(kc) = inst.koblitz.as_ref() else {
            continue;
        };
        let subfields: Vec<u32> = (2..n).filter(|d| n.is_multiple_of(*d)).collect();
        let cofactor_factors = small_factors(inst.cofactor);
        let target_dim = ((inst.r as f64).log2() / 3.0).ceil() as u32 + 2;
        let Some(idx) = koblitz_divisor_for(n, target_dim) else {
            eprintln!(
                "e15 n={n}: no proper invariant subspace (ord_n(2) = {})",
                ord2(n)
            );
            continue;
        };
        let Some(frob) = build_frobenius_factor_base_from_divisor(kc, &idx) else {
            eprintln!("e15 n={n}: divisor {idx:?} gives no base");
            continue;
        };
        if frob.points.len() > 60_000 {
            eprintln!(
                "e15 n={n}: the invariant subspace of dimension {} carries {} points, too many for a pair table",
                frob.ell,
                frob.points.len()
            );
            continue;
        }
        let Some(folded) = koblitz_factor_base(
            &inst,
            &frob,
            ColumnFold::SignedFrobeniusOrbit,
            "orbit".into(),
        ) else {
            continue;
        };
        let Some(control) =
            koblitz_factor_base(&inst, &frob, ColumnFold::Abscissa, "abscissa".into())
        else {
            continue;
        };
        let group = BinaryGroup(&inst.fast);
        for seed in 1..=o.seeds {
            let started = Instant::now();
            let planted = planted_for(seed, inst.r);
            let mut ops = GroupOps::default();
            let target = group.mul(&mut ops, inst.generator, planted);
            let ctx = InstanceCtx {
                group: &group,
                generator: inst.generator,
                target,
                r: inst.r,
                cofactor: inst.cofactor,
                group_order: inst.group_order,
                name: inst.name.clone(),
                field_degree: Some(n),
            };
            let mut oracle = MitmOracle::new(2);
            let mut params = Params::default();
            params.set("negation_folded", "1");
            oracle.prepare(&ctx, &folded, &params, &mut ops).unwrap();
            let rep = match full_rank_stream(
                &group,
                inst.generator,
                target,
                inst.r,
                inst.cofactor,
                planted,
                &folded,
                &control,
                seed,
                o.max_trials,
                2,
                None,
                |ops, ctr, pt| oracle.decompose(&ctx, &folded, ops, ctr, pt),
            ) {
                Ok(r) => r,
                Err(e) => {
                    eprintln!("e15 n={n}: {e}");
                    break;
                }
            };
            eprintln!(
                "e15 E_0 n={n} r=2^{:.1} h={} {:?} dim {}: cols {}/{} D {}/{} single-column rows {}/{} square {:?}/{:?} ok {}/{} [{:.1}s]",
                (inst.r as f64).log2(),
                inst.cofactor,
                cofactor_factors,
                frob.ell,
                rep.folded.columns,
                rep.control.columns,
                rep.folded.deficiency_total,
                rep.control.deficiency_total,
                rep.folded.single_column_rows,
                rep.control.single_column_rows,
                rep.folded.square_relations,
                rep.control.square_relations,
                rep.folded.verified,
                rep.control.verified,
                started.elapsed().as_secs_f64()
            );
            rows.push(json!({
                "experiment": "e15",
                "family": "koblitz-e0",
                "curve": "E_0: y^2 + xy = x^3 + 1",
                "n": n,
                "seed": seed,
                "instance": inst.name,
                "irreducible": format!("{:?}", inst.irreducible),
                "log2_r": (inst.r as f64).log2(),
                "r": inst.r,
                "group_order": inst.group_order,
                "cofactor": inst.cofactor,
                "cofactor_factors": cofactor_factors,
                "cofactor_is_e0_f2": inst.cofactor == 4,
                "ord_n_2": ord2(n),
                "intermediate_subfield_degrees": subfields,
                "eigenvalue_order": n,
                "divisor": idx,
                "subspace_dimension": frob.ell,
                "oracle": "mitm",
                "stream": stream_json(&rep),
                "wall_seconds": started.elapsed().as_secs_f64(),
            }));
        }
    }
    rows
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let mut o = Opts {
        exp: "e1".into(),
        bits: vec![16, 20],
        seeds: 2,
        rho_runs: 8,
        rho_max_steps: 1 << 28,
        max_trials: 50_000_000,
        json: None,
    };
    let mut i = 0;
    while i < args.len() {
        let next = |i: &mut usize| -> String {
            *i += 1;
            args[*i].clone()
        };
        match args[i].as_str() {
            "--exp" => o.exp = next(&mut i),
            "--bits" => {
                o.bits = next(&mut i)
                    .split(',')
                    .map(|s| s.parse().unwrap())
                    .collect()
            }
            "--seeds" => o.seeds = next(&mut i).parse().unwrap(),
            "--rho-runs" => o.rho_runs = next(&mut i).parse().unwrap(),
            "--rho-max-steps" => o.rho_max_steps = next(&mut i).parse().unwrap(),
            "--max-trials" => o.max_trials = next(&mut i).parse().unwrap(),
            "--json" => o.json = Some(next(&mut i)),
            other => panic!("unknown argument {other}"),
        }
        i += 1;
    }
    let started = Instant::now();
    let rows = match o.exp.as_str() {
        "e1" => e1(&o),
        "e2" => e2(&o),
        "e3" => e3(&o),
        "e4" => e4(&o),
        "e5" => e5(&o),
        "e6" => e6(&o),
        "e7" => e7(&o),
        "e8" => e8(&o),
        "e9" => e9(&o),
        "e11" => e11(&o),
        "e12" => e12(&o),
        "e12p" => e12p(&o),
        "e13" => e13(&o),
        "e15" => e15(&o),
        other => panic!("unknown experiment {other}; try e1..e9, e11, e12, e12p, e13"),
    };
    let out = json!({
        "what_this_is": format!("Experiment {} of research/notes/index-calculus/RESEARCH_GLV_INVARIANT_FACTOR_BASES.md: rows as measured, every logarithm checked against the planted one.", o.exp),
        "what_this_is_not": [
            "not a speed claim: relation counts and S are reported against their floors and a counted rho; the fold is an engineering lever unless a ratio to the floor falls",
            "not a claim about any deployed curve: toy instances with certified orders",
        ],
        "command": format!("glv_invariant_experiments {}", args.join(" ")),
        "host": {"os": std::env::consts::OS, "arch": std::env::consts::ARCH, "threads": std::thread::available_parallelism().map(|n| n.get()).unwrap_or(0)},
        "wall_seconds": started.elapsed().as_secs_f64(),
        "rows": rows,
    });
    println!("{}", serde_json::to_string(&out["rows"]).unwrap().len());
    if let Some(path) = o.json {
        fs::write(&path, serde_json::to_string_pretty(&out).unwrap()).unwrap();
        eprintln!(
            "wrote {path} ({} rows, {:.0}s)",
            out["rows"].as_array().unwrap().len(),
            started.elapsed().as_secs_f64()
        );
    }
}

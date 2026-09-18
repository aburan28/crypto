//! Crossbred against matrix-F4 and exhaustive search, on the Boolean
//! systems of Koblitz-curve index calculus.
//!
//! Three engines answer the same question — *does the target decompose
//! over the factor base?* — and this prices them in one unit.
//!
//! ```bash
//! cargo run --release --example crossbred_bench            # chained default ladder
//! cargo run --release --example crossbred_bench -- 9 11 13 # chosen n, m = 2
//! cargo run --release --example crossbred_bench -- --sym --no-sweep 9:3 15:3
//! ```
//!
//! `--sym` prices the Artin–Schreier `u`-frame (`build_symmetrised_system`
//! on `F_u`). `--a 0|1` selects `K_a`. `--no-sweep` skips the `(D, k)`
//! grid print.
//!
//! ## The unit, and the boundary
//!
//! Everything is reported in **bit operations**, with one stated
//! conversion: a 64-bit word XOR is 64 bit operations.  That is a
//! definition, not a measurement, and it is the only conversion used.
//!
//! The boundary for the engine comparison is exhaustive search over the
//! Boolean system: `2^v` points, each evaluating `T` monomials.  The
//! boundary for the index-calculus exponent (X1 / X3) is
//! `Q_enum = C(|F|, m−1)` word-ops at one word-op per pair.  That
//! undercounts enumeration, so `Q / Q_enum` is conservative against an
//! advance.
//!
//! Note what this does **not** claim.  Per `RESEARCH_SEMAEV_DECOMPOSITION.md`,
//! a faster decomposition oracle does not move the exponent of the
//! attack — relation collection stays `Θ(2^n)` with any oracle that is
//! polynomial in the factor base.  This table prices the oracle, which
//! is one phase, and says so.

use crypto_lib::binary_ecc::BinaryPoint;
use crypto_lib::cryptanalysis::crossbred::{
    extract_crossbred, format_sweep, solve_crossbred, sweep, CrossbredParams, SearchOptions,
};
use crypto_lib::cryptanalysis::koblitz_groebner::{
    build_decomposition_system, matrix_f4_f2_counted, solve_boolean_system, FieldStructure,
    SolveOptions, SolverEngine,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    build_frobenius_factor_base, KoblitzCurve,
};
use crypto_lib::cryptanalysis::koblitz_symmetrised::{
    build_symmetrised_factor_base, build_symmetrised_system, divisor_for_dimension,
};
use crypto_lib::cryptanalysis::pq_groebner_f2::F2BoolPoly;
use num_bigint::BigUint;
use std::time::Instant;

/// Bit operations per 64-bit word XOR.  The one conversion in this file.
const BITS_PER_WORD_OP: u64 = 64;

/// Exhaustive search over the Boolean system, as the reference oracle.
/// Returns the roots and the bit-operation count.
fn brute_force(system: &[F2BoolPoly], n_vars: usize) -> (Vec<u64>, u64) {
    let terms: u64 = system.iter().map(|p| p.terms.len() as u64).sum();
    let roots = (0..(1u64 << n_vars))
        .filter(|&p| system.iter().all(|q| q.eval(p) == 0))
        .collect();
    (roots, (1u64 << n_vars) * terms)
}

/// T4 of `RESEARCH_ECC2K130_ROUTE_TARGETS.md`, applied to the probe.
/// Smallest `D` then smallest `k` whose extraction determines the
/// remaining variables *and* whose search does not exhaust.
fn pick_determining_params(
    equations: &[F2BoolPoly],
    n_vars: usize,
    sys_deg: u32,
) -> Option<(usize, u32)> {
    let search_opts = SearchOptions::default();
    let k_hi = n_vars.min(search_opts.max_enumerated_bits as usize);
    if k_hi <= 2 {
        return None;
    }
    for d in sys_deg..=(sys_deg + 2) {
        for k in 2..k_hi {
            let solved = n_vars - k;
            if solved == 0 || solved > 64 {
                continue;
            }
            let params = CrossbredParams {
                macaulay_degree: d,
                enumerated: k,
                target_degree: 1,
                max_rows: 20_000,
            };
            let Some(xb) = extract_crossbred(equations, n_vars, &params) else {
                continue;
            };
            if xb.stats.kernel_dim < solved {
                continue;
            }
            let (_, stats) = solve_crossbred(equations, &xb, &search_opts);
            if !stats.exhausted {
                return Some((k, d));
            }
        }
    }
    None
}

fn q_enum_word_ops(m: usize, factor_base_points: u64) -> u64 {
    // One word-op per enumerated pair (m = 3) or point (m = 2). This
    // undercounts real enumeration, so Q / Q_enum is biased against
    // calling Crossbred an advance.
    if m <= 2 {
        factor_base_points
    } else {
        factor_base_points.saturating_mul(factor_base_points.saturating_sub(1)) / 2
    }
}

fn system_degree(equations: &[F2BoolPoly]) -> u32 {
    equations
        .iter()
        .flat_map(|e| e.terms.iter())
        .map(|t| t.mask.count_ones())
        .max()
        .unwrap_or(0)
}

fn f4_verdict(system: &[F2BoolPoly], n_vars: usize, max_degree: u32) -> (bool, u64, u32) {
    let base = system
        .iter()
        .flat_map(|p| p.terms.iter())
        .map(|t| t.mask.count_ones())
        .max()
        .unwrap_or(2)
        .max(2);
    let mut total = 0u64;
    for d in base..=max_degree.max(base) {
        match matrix_f4_f2_counted(system, n_vars, d) {
            Some((rows, ops)) => {
                total += ops;
                let refuted = rows
                    .iter()
                    .any(|p| p.terms.len() == 1 && p.terms[0].mask == 0);
                if refuted {
                    return (true, total, d);
                }
                if d == max_degree.max(base) {
                    return (false, total, d);
                }
            }
            None => return (false, total, d),
        }
    }
    (false, total, base)
}

struct Cli {
    symmetrised: bool,
    no_sweep: bool,
    probe_only: bool,
    curve_a: u8,
    cases: Vec<(u32, usize)>,
}

fn parse_cli() -> Cli {
    let raw: Vec<String> = std::env::args().skip(1).collect();
    let mut symmetrised = false;
    let mut no_sweep = false;
    let mut probe_only = false;
    let mut curve_a: u8 = 0;
    let mut rest: Vec<String> = Vec::new();
    let mut i = 0;
    while i < raw.len() {
        match raw[i].as_str() {
            "--sym" => symmetrised = true,
            "--no-sweep" => no_sweep = true,
            "--probe-only" => probe_only = true,
            "--a" => {
                i += 1;
                curve_a = raw
                    .get(i)
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(0);
            }
            other => rest.push(other.to_string()),
        }
        i += 1;
    }
    let cases = if rest.is_empty() {
        if symmetrised {
            vec![(9, 2), (9, 3), (15, 2), (15, 3)]
        } else {
            vec![(9, 2), (13, 2), (7, 3), (9, 3)]
        }
    } else {
        rest.iter()
            .filter_map(|a| {
                if let Some((n, m)) = a.split_once(':') {
                    Some((n.parse().ok()?, m.parse().ok()?))
                } else {
                    Some((a.parse().ok()?, 2))
                }
            })
            .collect()
    };
    Cli {
        symmetrised,
        no_sweep,
        probe_only,
        curve_a,
        cases,
    }
}

fn print_table_header() {
    println!(
        "| n | m | ℓ | v | |F| | Q_enum | Q_word | Q/Q_enum | deg | targets | reference | agree | brute (bit ops) | F4 (bit ops) | \
         crossbred (bit ops) | xb/brute | xb/F4 | D | k | kernel | filters | xb wall |"
    );
    println!(
        "|--:|--:|--:|--:|----:|-------:|-------:|---------:|----:|--------:|:----------|:-----:|----------------:|-------------:|\
         --------------------:|---------:|------:|--:|--:|-------:|--------:|--------:|"
    );
}

fn print_skip_row(n: u32, m: usize, ell: usize, v: usize, sys_deg: u32, reason: &str) {
    println!(
        "| {n} | {m} | {ell} | {v} | — | — | — | — | {sys_deg} | — | — | — | — | — | — \
         | — | — | — | — | — | — | {reason} |"
    );
}

fn print_footer() {
    println!();
    println!("X1/X3 unit: Q_word is 64-bit word XORs per target (extraction + search).");
    println!("Q_enum is C(|F|, m-1) word-ops at one word-op per pair (m=3) or point (m=2);");
    println!("that undercounts enumeration, so Q/Q_enum is conservative against an advance.");
    println!("`--sym` uses |F_u| and ℓ = dim V from divisor_for_dimension(n, (n+1)/m).");
    println!();
    println!("The `agree` column is the correctness gate, and it has two strengths.");
    println!("  `exhaustive, ==` : crossbred's root set equals exhaustive search's, exactly.");
    println!("  `F4 capped, ⊆`   : past v = 24 exhaustive search does not fit, and F4+splitting");
    println!("                     runs under its default max_solutions cap, so its list is a");
    println!("                     prefix of the root set.  The check is then two one-sided");
    println!("                     ones: every crossbred root verifies against the original");
    println!("                     equations, and every root F4 found is one crossbred found.");
    println!("A row that does not say `yes` is not a result.");
    println!();
    println!("What this prices is the decomposition oracle only.  Relation collection stays");
    println!("Θ(2^n) with any oracle polynomial in the factor base — see");
    println!("RESEARCH_SEMAEV_DECOMPOSITION.md, `What this does not buy`.");
}

fn price_built_system(
    n: u32,
    m: usize,
    ell: usize,
    f_pts: u64,
    v: usize,
    sys_deg: u32,
    k: usize,
    deg: u32,
    build: &dyn Fn(&BinaryPoint) -> Option<(Vec<F2BoolPoly>, usize)>,
    points: &[BinaryPoint],
) {
    let params = CrossbredParams {
        macaulay_degree: deg,
        enumerated: k,
        target_degree: 1,
        max_rows: 20_000,
    };
    let brute_ok = v <= 24;
    let q_enum = q_enum_word_ops(m, f_pts);
    let (mut brute_ops, mut f4_ops, mut xb_ops, mut xb_word) = (0u64, 0u64, 0u64, 0u64);
    let (mut kernel, mut filters) = (0usize, 0usize);
    let mut agree = true;
    let mut targets = 0u32;
    let mut wall = std::time::Duration::ZERO;

    for point in points {
        let Some((equations, nv)) = build(point) else {
            continue;
        };
        if nv != v {
            agree = false;
            continue;
        }
        targets += 1;

        let (roots, exact) = if brute_ok {
            let (roots, bops) = brute_force(&equations, v);
            brute_ops += bops;
            (roots, true)
        } else {
            let (mut roots, _) = solve_boolean_system(
                &equations,
                v,
                &SolveOptions {
                    engine: SolverEngine::MatrixF4 { max_degree: 4 },
                    ..Default::default()
                },
            );
            roots.sort_unstable();
            (roots, false)
        };

        let (refuted, fops, _) = f4_verdict(&equations, v, 4);
        f4_ops += fops * BITS_PER_WORD_OP;
        if refuted && !roots.is_empty() {
            agree = false;
        }

        let t0 = Instant::now();
        let Some(xb) = extract_crossbred(&equations, v, &params) else {
            agree = false;
            continue;
        };
        let (mut got, stats) = solve_crossbred(&equations, &xb, &SearchOptions::default());
        wall += t0.elapsed();
        kernel = xb.stats.kernel_dim;
        filters = xb.stats.filters;
        let words = xb.stats.word_ops
            + stats.transform_word_ops
            + stats.filter_word_ops
            + stats.solve_row_ops;
        xb_word += words;
        xb_ops += words * BITS_PER_WORD_OP;

        if stats.exhausted || xb.stats.kernel_dim < v.saturating_sub(k) {
            agree = false;
        }
        let mut expect = roots;
        got.sort_unstable();
        expect.sort_unstable();
        if !got
            .iter()
            .all(|&r| equations.iter().all(|p| p.eval(r) == 0))
        {
            agree = false;
        }
        if exact {
            if got != expect {
                agree = false;
            }
        } else if !expect.iter().all(|r| got.contains(r)) {
            agree = false;
        }
    }

    if targets == 0 {
        print_skip_row(n, m, ell, v, sys_deg, "no usable targets");
        return;
    }
    let per = |x: u64| x / targets as u64;
    let (b, f, x) = (per(brute_ops), per(f4_ops), per(xb_ops));
    let show = |val: u64| {
        if val == 0 {
            "—".to_string()
        } else {
            val.to_string()
        }
    };
    let ratio = |num: u64, den: u64| {
        if den > 0 {
            format!("{:.3}", num as f64 / den as f64)
        } else {
            "—".to_string()
        }
    };
    let qw = per(xb_word);
    println!(
        "| {n} | {m} | {ell} | {v} | {f_pts} | {q_enum} | {qw} | {} | {sys_deg} | {targets} | {} | {} | {} | {} | {x} | {} | {} \
         | {deg} | {k} | {kernel} | {filters} | {:.1} ms |",
        ratio(qw, q_enum),
        if brute_ok {
            "exhaustive, =="
        } else {
            "F4 capped, ⊆"
        },
        if agree { "yes" } else { "NO" },
        show(b),
        show(f),
        ratio(x, b),
        ratio(x, f),
        wall.as_secs_f64() * 1000.0 / targets as f64,
    );
}

fn candidate_points(kc: &KoblitzCurve, max_tries: u32) -> Vec<BinaryPoint> {
    let g = kc.generator().clone();
    (1..=max_tries)
        .filter_map(|t| {
            let point = kc.mul(&g, &BigUint::from(t * 7 + 1));
            match point {
                BinaryPoint::Infinity => None,
                p => Some(p),
            }
        })
        .collect()
}

fn first_probe(
    points: &[BinaryPoint],
    build: &dyn Fn(&BinaryPoint) -> Option<(Vec<F2BoolPoly>, usize)>,
) -> Option<(Vec<F2BoolPoly>, usize)> {
    points.iter().find_map(build)
}

fn main() {
    let cli = parse_cli();
    let frame = if cli.symmetrised {
        "symmetrised u / F_u"
    } else {
        "chained x / F_x"
    };

    if !cli.no_sweep {
        println!();
        if cli.symmetrised {
            println!("=== Crossbred space by (D, k), K_{}/F_2^9, m = 2, {frame} ===", cli.curve_a);
            println!();
            if let Some(kc) = KoblitzCurve::new(cli.curve_a, 9) {
                let target_dim = (9u32 + 1).div_ceil(2);
                if let Some(div) = divisor_for_dimension(9, target_dim) {
                    if let Some(fb) = build_symmetrised_factor_base(&kc, &div) {
                        let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
                        let g = kc.generator().clone();
                        let target = kc.mul(&g, &BigUint::from(53u32));
                        if let Some(sys) = build_symmetrised_system(&kc, &fb, &target, 2, &st) {
                            println!(
                                "unknowns v = {}, equations = {}, ℓ = {}, |F_u| = {}",
                                sys.n_vars,
                                sys.equations.len(),
                                fb.ell,
                                fb.points.len()
                            );
                            println!();
                            let cells = sweep(
                                &sys.equations,
                                sys.n_vars,
                                &[2, 3, 4],
                                &(2..sys.n_vars.min(8)).collect::<Vec<_>>(),
                                20_000,
                            );
                            print!("{}", format_sweep(&cells));
                        } else {
                            eprintln!("skip sweep: no symmetrised system at n=9 m=2");
                        }
                    }
                }
            } else {
                eprintln!("skip sweep: no K_{}/F_2^9", cli.curve_a);
            }
        } else {
            println!("=== Crossbred space by (D, k), K_0 / F_2^9, m = 2 ===");
            println!();
            let kc = KoblitzCurve::new(0, 9).expect("K_0 / F_2^9");
            let fb = build_frobenius_factor_base(&kc, 0).expect("factor base");
            let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
            let g = kc.generator().clone();
            let target = kc.mul(&g, &BigUint::from(53u32));
            let x_r = match &target {
                BinaryPoint::Affine { x, .. } => x.clone(),
                BinaryPoint::Infinity => unreachable!(),
            };
            let sys =
                build_decomposition_system(&fb.subspace_basis, &x_r, &kc.curve.b, 2, &st).unwrap();
            println!(
                "unknowns v = {}, equations = {}, ℓ = {}",
                sys.n_vars,
                sys.equations.len(),
                fb.ell
            );
            println!();
            let cells = sweep(
                &sys.equations,
                sys.n_vars,
                &[2, 3, 4],
                &(2..sys.n_vars.min(8)).collect::<Vec<_>>(),
                20_000,
            );
            print!("{}", format_sweep(&cells));
        }
    }

    println!();
    println!("=== Oracle cost per target, three engines, same verdicts ({frame}, K_{}) ===", cli.curve_a);
    println!();
    print_table_header();

    for &(n, m) in &cli.cases {
        let kc = match KoblitzCurve::new(cli.curve_a, n) {
            Some(c) => c,
            None => {
                eprintln!("skip n={n} m={m}: no KoblitzCurve K_{}", cli.curve_a);
                continue;
            }
        };
        let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
        let points = candidate_points(&kc, 40);

        if cli.symmetrised {
            let target_dim = (n + 1).div_ceil(m as u32);
            let Some(div) = divisor_for_dimension(n, target_dim) else {
                eprintln!("skip n={n} m={m}: no divisor for dim {target_dim}");
                continue;
            };
            let Some(fb) = build_symmetrised_factor_base(&kc, &div) else {
                eprintln!("skip n={n} m={m}: no F_u");
                continue;
            };
            let ell = fb.ell;
            let f_pts = fb.points.len() as u64;
            let build = |p: &BinaryPoint| {
                build_symmetrised_system(&kc, &fb, p, m, &st).map(|s| (s.equations, s.n_vars))
            };
            let probe7 = kc.mul(kc.generator(), &BigUint::from(7u32));
            let Some((probe_eqs, v)) = build(&probe7).or_else(|| first_probe(&points, &build))
            else {
                eprintln!("skip n={n} m={m}: no symmetrised probe system");
                continue;
            };
            let sys_deg = system_degree(&probe_eqs);
            eprintln!(
                "probe n={n} m={m} v={v} deg={sys_deg} |F_u|={f_pts} ell={ell} eqs={}",
                probe_eqs.len()
            );
            if cli.probe_only {
                println!(
                    "| {n} | {m} | {ell} | {v} | {f_pts} | {} | — | — | {sys_deg} | — | probe | — | — | — | — | — | — | — | — | — | — | probe only |",
                    q_enum_word_ops(m, f_pts)
                );
                continue;
            }
            let (k, deg) = match pick_determining_params(&probe_eqs, v, sys_deg) {
                Some(c) => {
                    eprintln!(
                        "T4 pick n={n} m={m} v={v} D={} k={} |F_u|={f_pts} ell={ell}",
                        c.1, c.0
                    );
                    c
                }
                None => {
                    print_skip_row(n, m, ell, v, sys_deg, "no determining space");
                    continue;
                }
            };
            price_built_system(n, m, ell, f_pts, v, sys_deg, k, deg, &build, &points);
        } else {
            let fb = match build_frobenius_factor_base(&kc, 0) {
                Some(f) => f,
                None => {
                    eprintln!("skip n={n} m={m}: no factor base");
                    continue;
                }
            };
            let ell = fb.ell as usize;
            let f_pts = fb.points.len() as u64;
            let b = kc.curve.b.clone();
            let build = |p: &BinaryPoint| {
                let x_r = match p {
                    BinaryPoint::Affine { x, .. } => x.clone(),
                    BinaryPoint::Infinity => return None,
                };
                build_decomposition_system(&fb.subspace_basis, &x_r, &b, m, &st)
                    .map(|s| (s.equations, s.n_vars))
            };
            // X1 freeze used generator×7 as the T4 probe. Keep it.
            let probe7 = kc.mul(kc.generator(), &BigUint::from(7u32));
            let Some((probe_eqs, v)) = build(&probe7).or_else(|| first_probe(&points, &build))
            else {
                eprintln!("skip n={n} m={m}: no decomposition system");
                continue;
            };
            let sys_deg = system_degree(&probe_eqs);
            let (k, deg) = match pick_determining_params(&probe_eqs, v, sys_deg) {
                Some(c) => {
                    eprintln!("T4 pick n={n} m={m} v={v} D={} k={}", c.1, c.0);
                    c
                }
                None => {
                    print_skip_row(n, m, ell, v, sys_deg, "no determining space");
                    continue;
                }
            };
            // Chained X1 used scalars t = 1..=12. Keep that set so the
            // frozen X1 rows stay comparable; extra candidates exist only
            // so a probe still exists if the first few are infinity.
            let chained_points: Vec<BinaryPoint> = points.iter().take(12).cloned().collect();
            price_built_system(
                n,
                m,
                ell,
                f_pts,
                v,
                sys_deg,
                k,
                deg,
                &build,
                &chained_points,
            );
        }
    }

    print_footer();
}

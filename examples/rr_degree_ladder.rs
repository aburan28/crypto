//! Refutation-degree ladder for the Riemann–Roch (Nagao) norm form at `m = 3`,
//! beside Semaev's direct `S₄` and a random control, on the same draws.
//!
//! One process per cell `(K_a, n, ℓ)`.  Each draw picks a random `ℓ`-dimensional
//! subspace `V` and a uniform target `R = [k]G` in the prime-order subgroup
//! (redrawn while `x(R) ∈ V`).  Three systems in the summand bits `x_i ∈ V`:
//!
//! * `rr`: the norm form of `f = x² + αx + γ + βy ∈ L(4O)`, `f(−R) = 0`, with
//!   `N(X) = (X + x_R)(X + x₁)(X + x₂)(X + x₃)` matched coefficient by coefficient.
//!   `β` is eliminated through the Artin–Schreier equation of the `X³` coefficient
//!   (`β = H(e₁) + ε`, one free bit `ε`, plus `Tr(e₁) = 0`), `α` through the constant
//!   coefficient.  Left: `3ℓ + 1` unknowns, `2n + 1` equations of Boolean degree 4.
//! * `x4`: `S₄(x₁, x₂, x₃, x_R)` Weil-descended: `3ℓ` unknowns, `n` equations, degree 6.
//! * `ctrl`: a random system of `rr`'s shape (unknowns, equations, degree, density).
//!
//! Roots are counted by group arithmetic, independently of every Macaulay code path:
//! `rr`'s roots are the ordered triples of factor-base points summing to `R`; `x4`'s are
//! the `x`-triples at which `S₄` vanishes in the field.  Only a system with no root is
//! measured, by `solving_degree` up to `--d-max` (`--d-max-x4` for the cheaper `x4`).  Output: one JSON line per draw and
//! arm, in the order `x4`, `rr`, `ctrl`, written as soon as it exists, plus a provisional
//! `at_least` line after every unresolved degree; the last line per (draw, arm) counts.
//!
//!     rr_degree_ladder --a 1 --n 17 --ell 4 --unsat 4 --max-draws 256 --d-max 8 \
//!         --seed 20260930 --out cell.jsonl
//!
//! `--check` brute-forces the `rr` root count over its cube and the `x4` system at every
//! point of `V³`, and compares both with the group-arithmetic counts (tiny cells only).
use crypto_lib::binary_ecc::{BinaryPoint, F2mElement, IrreduciblePoly};
use crypto_lib::cryptanalysis::koblitz_bench::{random_control_system, random_subspace_basis};
use crypto_lib::cryptanalysis::koblitz_groebner::{
    solving_profile_sparse, system_degree, FieldStructure, SymElement,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use crypto_lib::cryptanalysis::koblitz_symmetrised::{build_direct_x_system, plain_terms};
use crypto_lib::cryptanalysis::pq_groebner_f2::{F2BoolMono, F2BoolPoly};
use num_bigint::BigUint;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::collections::HashSet;
use std::io::Write;
use std::time::Instant;

fn word(e: &F2mElement) -> u64 {
    e.raw_bits().first().copied().unwrap_or(0)
}

fn trace(e: &F2mElement, n: u32, irr: &IrreduciblePoly) -> u64 {
    let (mut acc, mut t) = (e.clone(), e.clone());
    for _ in 1..n {
        t = t.square(irr);
        acc = acc.add(&t);
    }
    word(&acc)
}

/// Half trace, `n` odd: `H(c)² + H(c) = c` when `Tr(c) = 0`.
fn half_trace(e: &F2mElement, n: u32, irr: &IrreduciblePoly) -> F2mElement {
    let (mut acc, mut t) = (e.clone(), e.clone());
    for _ in 0..(n - 1) / 2 {
        t = t.square_k_times(2, irr);
        acc = acc.add(&t);
    }
    acc
}

/// The affine points of `K_a` with abscissa `x`.
fn points_over(kc: &KoblitzCurve, x: &F2mElement) -> Vec<BinaryPoint> {
    let (n, irr) = (kc.n, &kc.curve.irreducible);
    if x.is_zero() {
        return vec![BinaryPoint::Affine {
            x: x.clone(),
            y: F2mElement::one(n),
        }];
    }
    // y = x z, z² + z = x + a + 1/x²
    let c = x
        .add(&kc.curve.a)
        .add(&x.square(irr).flt_inverse(irr).expect("x ≠ 0"));
    if trace(&c, n, irr) != 0 {
        return vec![];
    }
    let z = half_trace(&c, n, irr);
    let y = x.mul(&z, irr);
    let pts = vec![
        BinaryPoint::Affine {
            x: x.clone(),
            y: y.clone(),
        },
        BinaryPoint::Affine {
            x: x.clone(),
            y: y.add(x),
        },
    ];
    assert!(pts.iter().all(|p| kc.curve.is_on_curve(p)), "point over x");
    pts
}

fn neg(p: &BinaryPoint) -> BinaryPoint {
    match p {
        BinaryPoint::Infinity => BinaryPoint::Infinity,
        BinaryPoint::Affine { x, y } => BinaryPoint::Affine {
            x: x.clone(),
            y: x.add(y),
        },
    }
}

fn span(basis: &[F2mElement], n: u32) -> Vec<F2mElement> {
    let mut out = vec![F2mElement::zero(n)];
    for b in basis {
        let more: Vec<F2mElement> = out.iter().map(|v| v.add(b)).collect();
        out.extend(more);
    }
    out
}

/// Ordered triples of factor-base points (every point with `x ∈ V`, the 2-torsion
/// point included) whose sum is `R`.
fn rr_root_count(
    kc: &KoblitzCurve,
    base: &[BinaryPoint],
    in_v: &HashSet<u64>,
    r: &BinaryPoint,
) -> u64 {
    let mut count = 0u64;
    for p1 in base {
        let r1 = kc.add(r, &neg(p1));
        for p2 in base {
            match kc.add(&r1, &neg(p2)) {
                BinaryPoint::Affine { x, .. } if in_v.contains(&word(&x)) => count += 1,
                _ => {}
            }
        }
    }
    count
}

fn pow(e: &F2mElement, k: u32, irr: &IrreduciblePoly) -> F2mElement {
    let mut acc = F2mElement::one(e.m_value());
    for _ in 0..k {
        acc = acc.mul(e, irr);
    }
    acc
}

/// `S₄(x₁, x₂, x₃, x_R)` evaluated in the field from the same term list the symbolic
/// system uses.
fn s4_field(xs: [&F2mElement; 4], terms: &[Vec<u32>], irr: &IrreduciblePoly) -> F2mElement {
    let n = xs[0].m_value();
    let mut acc = F2mElement::zero(n);
    for t in terms {
        let mut prod = F2mElement::one(n);
        for (i, &e) in t.iter().enumerate() {
            if e > 0 {
                prod = prod.mul(&pow(xs[i], e, irr), irr);
            }
        }
        acc = acc.add(&prod);
    }
    acc
}

/// `x`-triples in `V³` at which `S₄` vanishes in the field.
fn x4_root_count(
    v: &[F2mElement],
    x_r: &F2mElement,
    terms: &[Vec<u32>],
    irr: &IrreduciblePoly,
) -> u64 {
    let mut count = 0u64;
    for x1 in v {
        for x2 in v {
            for x3 in v {
                if s4_field([x1, x2, x3, x_r], terms, irr).is_zero() {
                    count += 1;
                }
            }
        }
    }
    count
}

fn sym_sqrt(e: &SymElement, st: &FieldStructure) -> SymElement {
    let mut t = e.clone();
    for _ in 1..st.n {
        t = t.square(st);
    }
    t
}

fn sym_half_trace(e: &SymElement, st: &FieldStructure) -> SymElement {
    let (mut acc, mut t) = (e.clone(), e.clone());
    for _ in 0..(st.n - 1) / 2 {
        t = t.square(st).square(st);
        acc = acc.add(&t);
    }
    acc
}

/// The trace as a Boolean polynomial; every other coordinate of `Σ e^{2^j}` must vanish
/// identically, which is checked.
fn sym_trace(e: &SymElement, st: &FieldStructure) -> F2BoolPoly {
    let (mut acc, mut t) = (e.clone(), e.clone());
    for _ in 1..st.n {
        t = t.square(st);
        acc = acc.add(&t);
    }
    assert!(
        acc.coords[1..].iter().all(|c| c.is_zero()),
        "trace lies in F_2"
    );
    acc.coords[0].clone()
}

/// The eliminated norm form: `2n + 1` equations in `3ℓ + 1` unknowns.
fn build_rr(
    kc: &KoblitzCurve,
    basis: &[F2mElement],
    r_pt: &BinaryPoint,
    st: &FieldStructure,
) -> Option<(Vec<F2BoolPoly>, usize)> {
    let (n, irr, ell) = (kc.n, &kc.curve.irreducible, basis.len());
    let n_vars = 3 * ell + 1;
    if n_vars > 64 || n % 2 == 0 {
        return None;
    }
    let BinaryPoint::Affine { x: r, y: s } = r_pt else {
        return None;
    };
    let c = |v: &F2mElement| SymElement::constant(v, n, n_vars);
    let one = F2mElement::one(n);
    let x: Vec<SymElement> = (0..3)
        .map(|i| SymElement::from_subspace_vars(basis, i * ell, n, n_vars))
        .collect();
    let rc = c(r);
    // elementary symmetric functions of {x_R, x₁, x₂, x₃}
    let s1 = x[0].add(&x[1]).add(&x[2]);
    let x12 = x[0].mul(&x[1], st);
    let s2 = x12.add(&x[0].mul(&x[2], st)).add(&x[1].mul(&x[2], st));
    let s3 = x12.mul(&x[2], st);
    let e1 = rc.add(&s1);
    let e2 = rc.mul(&s1, st).add(&s2);
    let e3 = rc.mul(&s2, st).add(&s3);
    let e4 = rc.mul(&s3, st);
    // β² + β = e₁  ⇒  β = H(e₁) + ε, with Tr(e₁) = 0 imposed
    let mut eps = SymElement::zero(n, n_vars);
    eps.coords[0] = F2BoolPoly::from_monos(vec![F2BoolMono::var((3 * ell) as u32)], n_vars);
    let beta = sym_half_trace(&e1, st).add(&eps);
    // γ² + β² = e₄, γ = x_R² + x_R α + (x_R + y_R) β
    //   ⇒  α = x_R⁻¹ √e₄ + x_R + x_R⁻¹ (x_R + y_R + 1) β
    let r_inv = r.flt_inverse(irr)?;
    let k1 = r_inv.mul(&r.add(s).add(&one), irr);
    let alpha = c(&r_inv)
        .mul(&sym_sqrt(&e4, st), st)
        .add(&rc)
        .add(&c(&k1).mul(&beta, st));
    let gamma = c(&r.square(irr))
        .add(&rc.mul(&alpha, st))
        .add(&c(&r.add(s)).mul(&beta, st));
    // identities that the eliminations must make exact
    let res4 = gamma.square(st).add(&beta.square(st)).add(&e4);
    assert!(
        res4.coords.iter().all(|p| p.is_zero()),
        "constant coefficient identity"
    );
    let res1 = beta.square(st).add(&beta).add(&e1);
    let tr = sym_trace(&e1, st);
    assert!(
        res1.coords[1..].iter().all(|p| p.is_zero()) && res1.coords[0] == tr,
        "X³ coefficient identity up to Tr(e₁)"
    );
    // X²:  α² + αβ + a β² = e₂        X¹:  βγ = e₃
    let ab = alpha.mul(&beta, st);
    let eq2 = alpha
        .square(st)
        .add(&ab)
        .add(&c(&kc.curve.a).mul(&beta.square(st), st))
        .add(&e2);
    let eq3 = beta.mul(&gamma, st).add(&e3);
    let mut eqs: Vec<F2BoolPoly> = eq2.coords.into_iter().chain(eq3.coords).collect();
    eqs.push(tr);
    eqs.retain(|p| !p.is_zero());
    Some((eqs, n_vars))
}

fn eval(p: &F2BoolPoly, a: u64) -> bool {
    p.terms
        .iter()
        .fold(false, |acc, t| acc ^ (t.mask & a == t.mask))
}

fn brute_roots(polys: &[F2BoolPoly], n_vars: usize) -> u64 {
    (0u64..1 << n_vars)
        .filter(|&a| polys.iter().all(|p| !eval(p, a)))
        .count() as u64
}

fn terms_per_eq(polys: &[F2BoolPoly]) -> usize {
    (polys.iter().map(|p| p.terms.len()).sum::<usize>() as f64 / polys.len() as f64).round()
        as usize
}

struct Measured {
    outcome: String,
    secs: f64,
}

/// `solving_degree`'s scan, one degree at a time, calling `partial(d, secs)` after every
/// unresolved degree below `d_max` so that a lower bound survives a kill.
fn measure(
    polys: &[F2BoolPoly],
    n_vars: usize,
    d_max: u32,
    mut partial: impl FnMut(u32, f64),
) -> Measured {
    let started = Instant::now();
    let (mut solved, mut built, mut refuted) = (None, None, false);
    for d in system_degree(polys).max(1)..=d_max {
        let Some(prof) = solving_profile_sparse(polys, n_vars, d) else {
            break;
        };
        built = Some(d);
        refuted = prof.refuted;
        if prof.resolves() {
            solved = Some(d);
            break;
        }
        if d < d_max {
            partial(d, started.elapsed().as_secs_f64());
        }
    }
    let outcome = match solved {
        Some(degree) => format!(r#"{{"kind":"resolved","degree":{degree},"refuted":{refuted}}}"#),
        None if built == Some(d_max) => format!(r#"{{"kind":"at_least","degree":{}}}"#, d_max + 1),
        None => format!(
            r#"{{"kind":"caps","built":{}}}"#,
            built.map_or("null".to_string(), |b| b.to_string())
        ),
    };
    Measured {
        outcome,
        secs: started.elapsed().as_secs_f64(),
    }
}

fn put(out: &mut Option<std::fs::File>, line: &str) {
    match out.as_mut() {
        Some(f) => {
            writeln!(f, "{line}").expect("write");
            f.flush().expect("flush");
        }
        None => println!("{line}"),
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let get = |name: &str| {
        args.iter()
            .position(|a| a == name)
            .and_then(|i| args.get(i + 1))
            .cloned()
    };
    let num = |name: &str, default: u64| -> u64 {
        get(name).map_or(default, |v| {
            v.parse().unwrap_or_else(|_| panic!("bad {name}"))
        })
    };
    let a = num("--a", 1) as u8;
    let n = num("--n", 17) as u32;
    let ell = num("--ell", 3) as usize;
    let want_unsat = num("--unsat", 4) as u32;
    let max_draws = num("--max-draws", 256) as u32;
    let d_max = num("--d-max", 7) as u32;
    let d_max_x4 = num("--d-max-x4", 9) as u32;
    let seed = num("--seed", 20260930);
    let check = args.iter().any(|s| s == "--check");
    assert!((2..=8).contains(&ell), "need 2 ≤ ℓ ≤ 8");
    let mut out = get("--out").map(|p| {
        assert!(
            !std::path::Path::new(&p).exists(),
            "--out exists; never overwritten"
        );
        std::fs::File::create(p).expect("create --out")
    });

    let kc = KoblitzCurve::new(a, n).unwrap_or_else(|| panic!("no curve K_{a}/2^{n}"));
    let irr = &kc.curve.irreducible;
    let r64: u64 = kc.subgroup_order.iter_u64_digits().next().unwrap_or(0);
    assert!(
        kc.subgroup_order.bits() <= 63 && r64 > 2,
        "subgroup order out of range"
    );
    let st = FieldStructure::new(n, irr);
    let terms = plain_terms(3).expect("S₄ terms");
    let cell_seed = seed ^ (u64::from(a) << 48) ^ (u64::from(n) << 40) ^ ((ell as u64) << 32);
    let mut rng = StdRng::seed_from_u64(cell_seed);
    let cell = format!("K{a}n{n}l{ell}");
    let g = kc.generator().clone();
    let (mut unsat, mut sat, mut skipped) = (0u32, 0u32, 0u32);

    for draw in 0..max_draws {
        if unsat >= want_unsat {
            break;
        }
        let basis = random_subspace_basis(n, ell, &mut rng);
        let k = 1 + rng.gen::<u64>() % (r64 - 1);
        let target = kc.mul(&g, &BigUint::from(k));
        let BinaryPoint::Affine { x: x_r, .. } = &target else {
            unreachable!("k in [1, r-1]")
        };
        let v = span(&basis, n);
        let in_v: HashSet<u64> = v.iter().map(word).collect();
        if in_v.contains(&word(x_r)) {
            skipped += 1;
            continue;
        }
        let base: Vec<BinaryPoint> = v.iter().flat_map(|x| points_over(&kc, x)).collect();
        let started = Instant::now();
        let (rr, rr_vars) = build_rr(&kc, &basis, &target, &st).expect("rr builds");
        let (x4, x4_vars) = build_direct_x_system(&basis, x_r, 3, &st).expect("x4 builds");
        let rr_roots = rr_root_count(&kc, &base, &in_v, &target);
        let x4_roots = x4_root_count(&v, x_r, &terms, irr);
        if check {
            let rr_brute = brute_roots(&rr, rr_vars);
            let x4_brute = brute_roots(&x4, x4_vars);
            eprintln!(
                "{cell} draw {draw}: rr roots group {rr_roots} brute {rr_brute}; x4 roots field {x4_roots} brute {x4_brute}"
            );
            assert_eq!(rr_roots, rr_brute, "rr root count");
            assert_eq!(x4_roots, x4_brute, "x4 root count");
        }
        let head = format!(
            r#""cell":"{cell}","a":{a},"n":{n},"ell":{ell},"m":3,"seed":{seed},"draw":{draw},"v_basis":{:?},"k":{k},"x_r":{},"base_points":{}"#,
            basis.iter().map(word).collect::<Vec<_>>(),
            word(x_r),
            base.len(),
        );
        // One line per arm, written as soon as it exists, in the order x4, rr, ctrl, and
        // after every unresolved degree a provisional `at_least` line (`"partial":true`):
        // the last line per (draw, arm) is the result, so a cell killed by its limit keeps
        // every arm it finished and a lower bound for the one it was on.
        let mut emit = |arm: &str,
                        polys: &[F2BoolPoly],
                        n_vars: usize,
                        roots: Option<u64>,
                        d_max: u32| {
            let shape = format!(
                r#""arm":"{arm}","n_vars":{n_vars},"n_eqs":{},"degree":{},"terms_per_eq":{}"#,
                polys.len(),
                system_degree(polys),
                terms_per_eq(polys)
            );
            let roots_field = roots.map_or(String::new(), |r| format!(r#","roots":{r}"#));
            let body = if roots == Some(0) || roots.is_none() {
                let m = measure(polys, n_vars, d_max, |d, secs| {
                    put(
                        &mut out,
                        &format!(
                            r#"{{{head},{shape}{roots_field},"outcome":{{"kind":"at_least","degree":{},"partial":true}},"secs":{secs:.3}}}"#,
                            d + 1
                        ),
                    );
                });
                format!(r#","outcome":{},"secs":{:.3}"#, m.outcome, m.secs)
            } else {
                r#","outcome":{"kind":"satisfiable"}"#.to_string()
            };
            put(&mut out, &format!("{{{head},{shape}{roots_field}{body}}}"));
        };
        emit("x4", &x4, x4_vars, Some(x4_roots), d_max_x4);
        emit("rr", &rr, rr_vars, Some(rr_roots), d_max);
        if rr_roots == 0 {
            let ctrl = random_control_system(
                rr_vars,
                rr.len(),
                system_degree(&rr),
                terms_per_eq(&rr),
                cell_seed ^ 0x0C01_7201_u64.wrapping_mul(u64::from(draw) + 1),
            );
            emit("ctrl", &ctrl, rr_vars, None, d_max);
            unsat += 1;
        } else {
            sat += 1;
        }
        let _ = started;
    }
    eprintln!("{cell}: {unsat} unsat, {sat} sat, {skipped} skipped (x_R ∈ V)");
}

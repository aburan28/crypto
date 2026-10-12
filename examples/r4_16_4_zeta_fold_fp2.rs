//! Round-4 lead 16.4: the free automorphism fold of the Gaudry base over
//! `F_{p²}`, against a 3-torsion-saturated base, on `j = 0` curves.
//!
//! Curve family: `E : y² = x³ + b` over `F_{p²}`, `p ≡ 1 (mod 3)`,
//! `b ∈ F_{p²} ∖ F_p` a square, `#E = 3r` with `r` prime.  `T = (0, √b)` is
//! a rational point of order 3 (the flex on `x = 0`), and
//! `φ : E → E' : y² = x³ − 27b`, `(x, y) ↦ ((x³ + 4b)/x², y(x³ − 8b)/x³)` is
//! the 3-isogeny with kernel `⟨T⟩`.  The group order is certified against
//! the six sextic-twist candidates `p² + 1 − Tr(ζ₆ᵏ π²)` (Cornacchia,
//! `4p = a² + 3v²`) by `[N]P = O` on random points; `φ` is checked to land
//! on `E'`, to be additive, and to kill `T`.
//!
//! Arms (all three have about `p` points, so they are compared at equal
//! base size, as Prop. 15.2 states the law):
//!
//! * **(a)** `F = {P ∈ E : x(P) ∈ F_p^*}`, columns = orbits under `⟨−1⟩`.
//! * **(b)** the same points, columns = orbits under `⟨−1, ζ⟩`,
//!   `ζ(x, y) = (ωx, y)`.  The base is checked to be `ζ`-closed with
//!   `Closure::Strict` (refused if not).
//! * **(c)** `F_c = φ⁻¹({Q ∈ E'(F_{p²}) : x(Q) ∈ F_p^*})`, the base that
//!   is closed under translation by `⟨T⟩` and is a line in the quotient
//!   coordinate `x' = x ∘ φ`.  On `E` it folds by `⟨−1, ζ, τ_T⟩` (18 points
//!   a column).  Through `φ` (an isomorphism `⟨G⟩ → ⟨G'⟩`) the same arm is
//!   run on `E'` with the base `F'' = {Q : x(Q) ∈ F_p^*, Q ∈ φ(E) = ⟨G'⟩}`
//!   folded by `⟨−1, ζ'⟩`; the column counts of the two views are checked
//!   to agree.
//!
//! Experiments (`--exp`):
//!
//! * `census` — exact counts on `E`: every 2-sum and every 3-sum of the base
//!   that equals each target (pair table + one sweep), raw and as distinct
//!   rows in each arm's columns.
//! * `pipeline` — the repository's two-summand Weil-descent oracle
//!   (`line_oracle::LineOracle`, line `c = 1`) on arms (a)/(b) over `E` and on
//!   arm (c) over `E'`; hits, `F_p` multiplications per call, agreement with
//!   an exhaustive subtraction check.
//! * `stream` — `glv_invariant_experiments::full_rank_stream` with that
//!   oracle: (b) against (a) on one stream, and (c) folded against (c)
//!   negation-only on one stream; relations and trials to full rank, and the
//!   recovered logarithm.
//! * `sweep` — a three-summand oracle written here (the pipeline has none
//!   over `F_{p²}` on a line: two descended equations in three `F_p`
//!   unknowns are positive-dimensional): sweep the first summand over the
//!   base, decompose `R − P₁` with the line oracle.  On (c) the sweep runs on
//!   `E'` over `F''`, i.e. over `F_c / ⟨T⟩`.
//!
//! Counted unit: `F_p` multiplications (the extension-field arithmetic's
//! thread-local counter plus the line solver's own count).  Inversions are
//! reported, not priced.
//!
//!     cargo run --release --example r4_16_4_zeta_fold_fp2 -- \
//!         --exp census --bits 9,10,11 --seeds 3 --targets 300 --json out.json

use std::collections::{BTreeMap, HashSet};
use std::env;
use std::fs;
use std::time::Instant;

use crypto_lib::cryptanalysis::ext_curve::{
    random_point, take_field_counters, ExtCurve, ExtDiagonalAutomorphism, ExtField, ExtPoint, Fp2,
    Fp2El,
};
use crypto_lib::cryptanalysis::glv_invariant_base::{
    cornacchia_4p, factor_u64, fold_by_endomorphisms, verify_endomorphism, Closure, Endomorphism,
    Negation,
};
use crypto_lib::cryptanalysis::glv_invariant_experiments::full_rank_stream;
use crypto_lib::cryptanalysis::ic_boundary::{CountedGroup, FactorBase, GroupOps, OracleCounters};
use crypto_lib::cryptanalysis::ic_framework::stages::{DecompositionOracle, InstanceCtx, Params};
use crypto_lib::cryptanalysis::line_oracle::LineOracle;
use crypto_lib::cryptanalysis::residual_walk::{inv_mod, is_prime_u64, pow_mod, sqrt_mod};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};

type Curve = ExtCurve<Fp2>;
type Pt = ExtPoint<Fp2El>;
type Fb = FactorBase<Pt>;

// ── Instance ────────────────────────────────────────────────────────

struct Inst {
    p: u64,
    f: Fp2,
    b: Fp2El,
    e: Curve,
    ep: Curve,
    n: u64,
    r: u64,
    h: u64,
    torsion: bool,
    split: Option<(u64, u64)>,
    trace: i64,
    g: Pt,
    gp: Pt,
    t: Pt,
    omega: u64,
    zeta: ExtDiagonalAutomorphism<Fp2El>,
    zetap: ExtDiagonalAutomorphism<Fp2El>,
    k3: u64,
}

/// The six candidate orders `p² + 1 − Tr(ζ₆ᵏ π₁²)` of `j = 0` curves over
/// `F_{p²}`, with `π₁ = (a + v√−3)/2`, `4p = a² + 3v²`.
fn j0_fp2_candidates(p: u64) -> Result<Vec<(u64, i64)>, String> {
    let (a, v) = cornacchia_4p(p, 3).ok_or("no Cornacchia solution for 4p = a² + 3v²")?;
    let (a, v) = (a as i128, v as i128);
    // (A + B√−3)/2 = π₁².
    let mut big_a = (a * a - 3 * v * v) / 2;
    let mut big_b = a * v;
    let p2 = (p as i128) * (p as i128);
    let mut out = Vec::new();
    for _ in 0..6 {
        if big_a * big_a + 3 * big_b * big_b != 4 * p2 {
            return Err("norm of the Frobenius candidate is not p²".into());
        }
        out.push(((p2 + 1 - big_a) as u64, big_a as i64));
        let (na, nb) = ((big_a - 3 * big_b) / 2, (big_a + big_b) / 2);
        big_a = na;
        big_b = nb;
    }
    Ok(out)
}

fn inf(f: &Fp2) -> Pt {
    ExtPoint::infinity(f.zero())
}

/// `φ_b : y² = x³ + b → y² = x³ − 27b`, kernel `x = 0`.
fn phi(f: &Fp2, b: Fp2El, pt: Pt) -> Pt {
    if pt.infinity || f.is_zero(pt.x) {
        return inf(f);
    }
    let x2 = f.sqr(pt.x);
    let x3 = f.mul(x2, pt.x);
    let inv3 = f.inv(x3).expect("x ≠ 0");
    let xx = f.mul(f.mul(f.add(x3, f.scale(b, 4)), pt.x), inv3);
    let yy = f.mul(f.mul(pt.y, f.sub(x3, f.scale(b, 8))), inv3);
    ExtPoint::affine(xx, yy)
}

/// `E' → E`: `φ_{b'}` to `y² = x³ + 729b`, then `(x, y) ↦ (x/9, y/27)`.
fn psi(f: &Fp2, bp: Fp2El, pt: Pt) -> Pt {
    let q = phi(f, bp, pt);
    if q.infinity {
        return q;
    }
    let p = f.p;
    ExtPoint::affine(f.scale(q.x, inv_mod(9, p)), f.scale(q.y, inv_mod(27, p)))
}

fn zeta_on(
    c: &Curve,
    gen: Pt,
    r: u64,
    omega: u64,
    seed: u64,
) -> Result<ExtDiagonalAutomorphism<Fp2El>, String> {
    let f = &c.f;
    let s = sqrt_mod(r - 3, r).ok_or("−3 is not a square mod r")?;
    let half = inv_mod(2, r);
    let cands = [
        ((r - 1 + s) % r * half) % r,
        (((r - 1) + (r - s)) % r * half) % r,
    ];
    let img = ExtPoint::affine(f.scale(gen.x, omega), gen.y);
    let mut ops = GroupOps::default();
    let lambda = cands
        .into_iter()
        .find(|&l| c.mul(&mut ops, gen, l) == img)
        .ok_or("neither root of λ² + λ + 1 is ζ's eigenvalue")?;
    let z = ExtDiagonalAutomorphism {
        cx: f.from_fp(omega),
        cy: f.one(),
        eigenvalue: lambda,
        order: 3,
        label: "zeta3",
    };
    verify_endomorphism(c, gen, r, &z, 20, seed)?;
    Ok(z)
}

/// `#E = Norm(α − 1)·Norm(α + 1)` for some `α = ζ₆ʲ π₁`: the order of the
/// base change of an `F_p` curve (and its quadratic twist).
fn base_change_split(p: u64, n: u64) -> Option<(u64, u64)> {
    let (a, v) = cornacchia_4p(p, 3)?;
    let (mut big_a, mut big_b) = (a as i128, v as i128);
    for _ in 0..6 {
        let n1 = ((big_a - 2) * (big_a - 2) + 3 * big_b * big_b) / 4;
        let n2 = ((big_a + 2) * (big_a + 2) + 3 * big_b * big_b) / 4;
        if n1 * n2 == n as i128 {
            return Some((n1 as u64, n2 as u64));
        }
        let (na, nb) = ((big_a - 3 * big_b) / 2, (big_a + big_b) / 2);
        big_a = na;
        big_b = nb;
    }
    None
}

/// `torsion = true`: `b` a square (rational `T = (0, √b)` of order 3),
/// `9 ∤ #E`, `r` the largest prime factor.  `torsion = false`: `b` a
/// non-square (no rational 3-torsion), `#E = h·r` with `h ≤ 16`.
fn build(p_bits: u32, seed: u64, torsion: bool) -> Result<Inst, String> {
    let mut rng =
        StdRng::seed_from_u64(seed ^ 0x1640_F2F2 ^ ((p_bits as u64) << 40) ^ u64::from(torsion));
    let mut ops = GroupOps::default();
    for _ in 0..2000 {
        let p = loop {
            let c = rng.gen_range((1u64 << (p_bits - 1))..(1u64 << p_bits)) | 1;
            if c % 3 == 1 && is_prime_u64(c) {
                break c;
            }
        };
        let f = Fp2::new(p)?;
        let cands = j0_fp2_candidates(p)?;
        for _ in 0..12 {
            let b: Fp2El = [rng.gen_range(0..p), rng.gen_range(1..p)];
            if f.is_square(b) != torsion {
                continue;
            }
            let e = Curve { f, a: f.zero(), b };
            let pts: Vec<Pt> = (0..4).map(|_| random_point(&e, &mut rng)).collect();
            let mut matches: Vec<(u64, i64)> = cands
                .iter()
                .copied()
                .filter(|&(n, _)| pts.iter().all(|&q| e.mul(&mut ops, q, n).infinity))
                .collect();
            matches.dedup();
            if matches.is_empty() {
                return Err(format!("no sextic-twist candidate order fits at p = {p}"));
            }
            if matches.len() > 1 {
                continue;
            }
            let (n, trace) = matches[0];
            let factors = factor_u64(n);
            let &(r, er) = factors.last().ok_or("n = 1")?;
            let h = n / r;
            if torsion {
                if n % 3 != 0 {
                    return Err(format!("b square but 3 ∤ #E at p = {p}"));
                }
                // r ≤ (√p + 1)² here (the order is a base-change product);
                // keep the curves whose r is a large part of one factor.
                if n % 9 == 0 || r <= 3 || er != 1 || r < p / 4 {
                    continue;
                }
            } else if h > 16 || er != 1 || r <= 3 {
                continue;
            }
            let g = loop {
                let q = e.mul(&mut ops, random_point(&e, &mut rng), h);
                if !q.infinity {
                    break q;
                }
            };
            if !e.mul(&mut ops, g, r).infinity {
                return Err("[r]G ≠ O".into());
            }
            let omega = (2..p)
                .map(|x| pow_mod(x, (p - 1) / 3, p))
                .find(|&w| w != 1)
                .ok_or("no cube root of unity")?;
            let zeta = zeta_on(&e, g, r, omega, seed)?;
            let split = base_change_split(p, n);
            if !torsion {
                return Ok(Inst {
                    p,
                    f,
                    b,
                    e,
                    ep: e,
                    n,
                    r,
                    h,
                    torsion,
                    split,
                    trace,
                    g,
                    gp: g,
                    t: inf(&f),
                    omega,
                    zeta,
                    zetap: zeta,
                    k3: 0,
                });
            }
            let t = ExtPoint::affine(f.zero(), f.sqrt(b).ok_or("b not a square")?);
            if !e.is_on_curve(t) || !e.mul(&mut ops, t, 3).infinity {
                return Err("T is not a rational 3-torsion point".into());
            }
            let bp = f.neg(f.scale(b, 27));
            let ep = Curve {
                f,
                a: f.zero(),
                b: bp,
            };
            if !phi(&f, b, t).infinity {
                return Err("φ(T) ≠ O".into());
            }
            for _ in 0..10 {
                let x = random_point(&e, &mut rng);
                let y = random_point(&e, &mut rng);
                let (px, py) = (phi(&f, b, x), phi(&f, b, y));
                if !ep.is_on_curve(px) {
                    return Err("φ(P) is not on E'".into());
                }
                if phi(&f, b, e.add(&mut ops, x, y)) != ep.add(&mut ops, px, py) {
                    return Err("φ is not additive".into());
                }
            }
            let gp = phi(&f, b, g);
            if gp.infinity || !ep.mul(&mut ops, gp, r).infinity {
                return Err("φ(G) is not of order r".into());
            }
            let zetap = zeta_on(&ep, gp, r, omega, seed ^ 1)?;
            return Ok(Inst {
                p,
                f,
                b,
                e,
                ep,
                n,
                r,
                h,
                torsion,
                split,
                trace,
                g,
                gp,
                t,
                omega,
                zeta,
                zetap,
                // φ(E) is the prime-to-3 part of E' (9 ∤ #E), where [3] is invertible.
                k3: inv_mod(3, n / 3),
            });
        }
    }
    Err(format!(
        "no instance at {p_bits} bits (torsion = {torsion})"
    ))
}

// ── Bases ───────────────────────────────────────────────────────────

/// Translation by `T`: on `[h]`-projected logarithms it acts as `1`.
struct TranslateT {
    t: Pt,
}

impl Endomorphism<Curve> for TranslateT {
    fn name(&self) -> String {
        "translate-T".into()
    }
    fn degree(&self) -> u64 {
        1
    }
    fn eigenvalue(&self) -> u64 {
        1
    }
    fn apply(&self, g: &Curve, p: Pt) -> Pt {
        g.add(&mut GroupOps::default(), p, self.t)
    }
}

fn fold(
    c: &Curve,
    r: u64,
    h: u64,
    seed: Vec<Pt>,
    gens: &[&dyn Endomorphism<Curve>],
    d: &str,
) -> Fb {
    let (fb, _) = fold_by_endomorphisms(
        c,
        r,
        h,
        seed,
        gens,
        Closure::Strict,
        32,
        |q| c.key(q),
        |q| c.f.pack(q.x),
        d.to_string(),
    )
    .unwrap_or_else(|e| panic!("fold {d}: {e}"));
    fb
}

fn line_points(c: &Curve, p: u64) -> Vec<Pt> {
    (1..p).flat_map(|t| c.lift_x(c.f.from_fp(t))).collect()
}

struct Bases {
    fa: Fb,          // E, x ∈ F_p^*, ⟨−1⟩
    fb: Fb,          // same points, ⟨−1, ζ⟩
    cz: Fb,          // E', F'' ⊂ ⟨G'⟩, ⟨−1, ζ'⟩
    cn: Fb,          // E', F'', ⟨−1⟩
    c18: Option<Fb>, // E, F_c = φ⁻¹(F''), ⟨−1, ζ, τ_T⟩
    c6: Option<Fb>,  // E, F_c, ⟨−1, τ_T⟩
    info: Value,
}

fn build_bases(inst: &Inst, with_fc_on_e: bool) -> Bases {
    let started = Instant::now();
    let (e, ep, f, r) = (&inst.e, &inst.ep, &inst.f, inst.r);
    let neg = Negation { r };
    let seed_a = line_points(e, inst.p);
    let mut ops = GroupOps::default();
    let in_g = seed_a
        .iter()
        .filter(|&&q| e.mul(&mut ops, q, r).infinity)
        .count();
    let fa = fold(e, r, inst.h, seed_a.clone(), &[&neg], "x ∈ F_p^*, ⟨−1⟩");
    let fb = fold(
        e,
        r,
        inst.h,
        seed_a,
        &[&neg, &inst.zeta],
        "x ∈ F_p^*, ⟨−1, ζ⟩",
    );
    // F'' on E': x ∈ F_p^*, filtered to ⟨G'⟩ = φ(E); the cube test is
    // checked against [r]Q = O on every point.
    let beta = if inst.torsion {
        f.sqrt(ep.b).expect("b' square")
    } else {
        f.zero()
    };
    let cube_exp = (inst.p * inst.p - 1) / 3;
    let (mut agree, mut disagree, mut line_pts) = (0u64, 0u64, 0u64);
    let mut seed_c = Vec::new();
    for q in if inst.torsion {
        line_points(ep, inst.p)
    } else {
        Vec::new()
    } {
        line_pts += 1;
        let truth = ep.mul(&mut ops, q, inst.n / 3).infinity;
        let z = f.sub(q.y, beta);
        let cube = f.pow(z, cube_exp) == f.one();
        if cube == truth {
            agree += 1;
        } else {
            disagree += 1;
        }
        if truth {
            seed_c.push(q);
        }
    }
    // Without rational 3-torsion there is no arm (c); the slots repeat (b)
    // and (a) so the drivers stay uniform, and every (c) field is marked
    // invalid in the output.
    if !inst.torsion {
        seed_c = line_points(e, inst.p);
    }
    let cz = fold(
        ep,
        r,
        inst.h,
        seed_c.clone(),
        &[&neg, &inst.zetap],
        "E': F'', ⟨−1, ζ'⟩",
    );
    let cn = fold(ep, r, inst.h, seed_c, &[&neg], "E': F'', ⟨−1⟩");
    // F_c on E: preimages [3⁻¹]ψ(Q) + ⟨T⟩.
    let tt = TranslateT { t: inst.t };
    let mut fc_pts = Vec::new();
    let mut fc_ok = true;
    if with_fc_on_e {
        let mut keys = HashSet::new();
        let mt = e.neg(inst.t);
        for &q in &cz.points {
            let p1 = e.mul(&mut ops, psi(f, ep.b, q), inst.k3);
            let img = phi(f, inst.b, p1);
            if cz.index_of_key(ep.key(&img)).is_none() {
                fc_ok = false;
            }
            for s in [p1, e.add(&mut ops, p1, inst.t), e.add(&mut ops, p1, mt)] {
                if keys.insert(e.key(&s)) {
                    fc_pts.push(s);
                }
            }
        }
        if fc_pts.len() != 3 * cz.points.len() {
            fc_ok = false;
        }
    }
    let (c18, c6) = if with_fc_on_e {
        (
            Some(fold(
                e,
                r,
                inst.h,
                fc_pts.clone(),
                &[&neg, &inst.zeta, &tt],
                "F_c, ⟨−1, ζ, τ_T⟩",
            )),
            Some(fold(e, r, inst.h, fc_pts, &[&neg, &tt], "F_c, ⟨−1, τ_T⟩")),
        )
    } else {
        (None, None)
    };
    let info = json!({
        "points_a": fa.points.len(),
        "fraction_of_F_in_G": in_g as f64 / fa.points.len() as f64,
        "U_a": fa.columns, "U_b": fb.columns,
        "points_per_col_a": fa.points.len() as f64 / fa.columns as f64,
        "points_per_col_b": fb.points.len() as f64 / fb.columns as f64,
        "zeta_closed_strict": true,
        "c_arm_valid": inst.torsion,
        "E_prime_line_points": line_pts,
        "F2_points_on_E_prime": cz.points.len(),
        "cube_test_agree": agree, "cube_test_disagree": disagree,
        "U_c_on_E_prime": cz.columns, "U_c_neg_on_E_prime": cn.columns,
        "Fc_on_E_built": with_fc_on_e,
        "Fc_points_on_E": c18.as_ref().map(|b| b.points.len()),
        "Fc_preimage_checks_ok": fc_ok,
        "U_c18_on_E": c18.as_ref().map(|b| b.columns), "U_c6_on_E": c6.as_ref().map(|b| b.columns),
        "build_wall_s_contended": started.elapsed().as_secs_f64(),
    });
    Bases {
        fa,
        fb,
        cz,
        cn,
        c18,
        c6,
        info,
    }
}

fn inst_json(inst: &Inst) -> Value {
    json!({
        "family": if inst.torsion {
            "y^2 = x^3 + b over F_{p^2}, b not in F_p, b a square (T = (0, sqrt b) rational), 9 does not divide #E"
        } else {
            "y^2 = x^3 + b over F_{p^2}, b a non-square (no rational 3-torsion), #E = h r, h <= 16"
        },
        "p": inst.p, "b": inst.b, "nu": inst.f.nu,
        "group_order": inst.n, "r": inst.r, "log2_r": (inst.r as f64).log2(),
        "cofactor": inst.h, "rational_3_torsion": inst.torsion,
        "order_is_base_change_product_N1_N2": inst.split, "trace_fp2": inst.trace, "omega": inst.omega,
        "lambda_zeta": inst.zeta.eigenvalue, "lambda_zeta_prime": inst.zetap.eigenvalue,
    })
}

fn random_targets(inst: &Inst, n: usize, seed: u64) -> Vec<Pt> {
    let mut rng = StdRng::seed_from_u64(seed ^ 0x7A12_6E75);
    let mut ops = GroupOps::default();
    (0..n)
        .map(|_| inst.e.mul(&mut ops, inst.g, rng.gen_range(1..inst.r)))
        .collect()
}

// ── Projection of relations ─────────────────────────────────────────

fn proj_row(idx: &[usize], fb: &Fb, r: u64) -> Vec<(usize, u64)> {
    let mut m: BTreeMap<usize, u64> = BTreeMap::new();
    for &i in idx {
        let e = m.entry(fb.col_of[i]).or_insert(0);
        *e = (*e + fb.coef_of[i]) % r;
    }
    m.into_iter().filter(|&(_, c)| c != 0).collect()
}

// ── Census: exact 2- and 3-sums on E ────────────────────────────────

struct PairTable {
    v: Vec<(u64, u32, u32)>,
}

impl PairTable {
    fn new(c: &Curve, pts: &[Pt]) -> Self {
        let mut ops = GroupOps::default();
        let mut v = Vec::with_capacity(pts.len() * (pts.len() + 1) / 2);
        for i in 0..pts.len() {
            for j in i..pts.len() {
                let s = c.add(&mut ops, pts[i], pts[j]);
                v.push((c.key(&s), i as u32, j as u32));
            }
        }
        v.sort_unstable();
        Self { v }
    }
    fn get(&self, key: u64) -> &[(u64, u32, u32)] {
        let lo = self.v.partition_point(|x| x.0 < key);
        let hi = self.v.partition_point(|x| x.0 <= key);
        &self.v[lo..hi]
    }
}

/// Per target: raw counts and distinct rows under each projection.
fn census_base(
    c: &Curve,
    pts: &[Pt],
    projs: &[&Fb],
    r: u64,
    targets: &[Pt],
) -> Vec<(usize, Vec<usize>, usize, Vec<usize>)> {
    let table = PairTable::new(c, pts);
    let mut ops = GroupOps::default();
    let mut out = Vec::new();
    for &rt in targets {
        let pairs = table.get(c.key(&rt));
        let n2_raw = pairs.len();
        let n2p: Vec<usize> = projs
            .iter()
            .map(|fb| {
                pairs
                    .iter()
                    .map(|&(_, i, j)| proj_row(&[i as usize, j as usize], fb, r))
                    .collect::<HashSet<_>>()
                    .len()
            })
            .collect();
        let mut triples: HashSet<[u32; 3]> = HashSet::new();
        for (i0, &p0) in pts.iter().enumerate() {
            let s = c.add(&mut ops, rt, c.neg(p0));
            for &(_, i, j) in table.get(c.key(&s)) {
                let mut t = [i0 as u32, i, j];
                t.sort_unstable();
                triples.insert(t);
            }
        }
        let n3p: Vec<usize> = projs
            .iter()
            .map(|fb| {
                triples
                    .iter()
                    .map(|t| proj_row(&[t[0] as usize, t[1] as usize, t[2] as usize], fb, r))
                    .collect::<HashSet<_>>()
                    .len()
            })
            .collect();
        out.push((n2_raw, n2p, triples.len(), n3p));
    }
    out
}

fn mean(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / v.len().max(1) as f64
}

/// Paired bootstrap of `Σnum / Σden` over targets: (point, lo, hi).
fn boot_ratio(num: &[f64], den: &[f64], seed: u64) -> (f64, f64, f64) {
    let n = num.len();
    let point = num.iter().sum::<f64>() / den.iter().sum::<f64>();
    let mut rng = StdRng::seed_from_u64(seed ^ 0xB007);
    let mut reps: Vec<f64> = (0..2000)
        .map(|_| {
            let (mut a, mut b) = (0.0, 0.0);
            for _ in 0..n {
                let k = rng.gen_range(0..n);
                a += num[k];
                b += den[k];
            }
            if b > 0.0 {
                a / b
            } else {
                f64::INFINITY
            }
        })
        .collect();
    reps.sort_by(|x, y| x.partial_cmp(y).unwrap());
    (point, reps[50], reps[1949])
}

fn census(bits: &[u32], seeds: u64, ntargets: usize, torsion: bool) -> Vec<Value> {
    let mut rows = Vec::new();
    for &b in bits {
        for s in 1..=seeds {
            let started = Instant::now();
            let inst = match build(b, s, torsion) {
                Ok(i) => i,
                Err(e) => {
                    eprintln!("census {b} seed {s}: {e}");
                    continue;
                }
            };
            let bases = build_bases(&inst, true);
            let targets = random_targets(&inst, ntargets, s);
            let r = inst.r;
            let ab = census_base(
                &inst.e,
                &bases.fa.points,
                &[&bases.fa, &bases.fb],
                r,
                &targets,
            );
            let (c18, c6) = (bases.c18.as_ref().unwrap(), bases.c6.as_ref().unwrap());
            let cc = census_base(&inst.e, &c18.points, &[c18, c6], r, &targets);
            type Cen = [(usize, Vec<usize>, usize, Vec<usize>)];
            let n2 = |v: &Cen, k: usize| -> Vec<f64> { v.iter().map(|x| x.1[k] as f64).collect() };
            let n3 = |v: &Cen, k: usize| -> Vec<f64> { v.iter().map(|x| x.3[k] as f64).collect() };
            let hit = |v: &[f64]| v.iter().filter(|&&x| x > 0.0).count() as f64 / v.len() as f64;
            let a2r: Vec<f64> = ab.iter().map(|x| x.0 as f64).collect();
            let a3r: Vec<f64> = ab.iter().map(|x| x.2 as f64).collect();
            let c2r: Vec<f64> = cc.iter().map(|x| x.0 as f64).collect();
            let c3r: Vec<f64> = cc.iter().map(|x| x.2 as f64).collect();
            // Projections: ab = [a, b], cc = [c18, c6].
            let (a2p, b2p, a3p, b3p) = (n2(&ab, 0), n2(&ab, 1), n3(&ab, 0), n3(&ab, 1));
            let (c2p, c62p, c3p, c63p) = (n2(&cc, 0), n2(&cc, 1), n3(&cc, 0), n3(&cc, 1));
            let (ua, ub, uc) = (
                bases.fa.columns as f64,
                bases.fb.columns as f64,
                c18.columns as f64,
            );
            // Relations per unknown per target, (c) over (b): the law's
            // class factor and column fold together.
            let rpu2 = boot_ratio(
                &c2p.iter().map(|x| x / uc).collect::<Vec<_>>(),
                &b2p.iter().map(|x| x / ub).collect::<Vec<_>>(),
                s,
            );
            let rpu3 = boot_ratio(
                &c3p.iter().map(|x| x / uc).collect::<Vec<_>>(),
                &b3p.iter().map(|x| x / ub).collect::<Vec<_>>(),
                s + 7,
            );
            let d2 = boot_ratio(&c2p, &b2p, s + 11);
            let ratio_f = c18.points.len() as f64 / bases.fa.points.len() as f64;
            let d3 = boot_ratio(&c3p, &b3p, s + 13);
            let row = json!({
                "experiment": "census", "bits": b, "seed": s,
                "instance": inst_json(&inst), "bases": bases.info,
                "targets": ntargets,
                "mean_n2_raw_F": mean(&a2r), "mean_n2_proj_a": mean(&a2p),
                "mean_n2_proj_b": mean(&b2p),
                "mean_n3_raw_F": mean(&a3r), "mean_n3_proj_a": mean(&a3p),
                "mean_n3_proj_b": mean(&b3p),
                "mean_n2_raw_Fc": mean(&c2r), "mean_n2_proj_c18": mean(&c2p),
                "mean_n3_raw_Fc": mean(&c3r), "mean_n3_proj_c18": mean(&c3p),
                "mean_n2_proj_c6": mean(&c62p),
                "mean_n3_proj_c6": mean(&c63p),
                "hit2_b": hit(&b2p), "hit2_c": hit(&c2p), "hit3_b": hit(&b3p), "hit3_c": hit(&c3p),
                "ratio_distinct_rows_c_over_b_m2": d2, "ratio_distinct_rows_c_over_b_m3": d3,
                "law_distinct_rows_c_over_b": [1.0/3.0, 1.0/9.0],
                "ratio_rows_per_unknown_c_over_b_m2": rpu2,
                "ratio_rows_per_unknown_c_over_b_m3": rpu3,
                "law_rows_per_unknown_c_over_b": [1.0, 1.0/3.0],
                "U_b_over_U_a": ub / ua, "U_c_over_U_b": uc / ub,
                // The law is stated at equal |F|; |F_c| and |F| differ by a
                // few percent per curve, so the size-matched ratios divide
                // the m-sum counts by (|F_c|/|F|)^m and U by |F_c|/|F|.
                "base_size_ratio_Fc_over_F": ratio_f,
                "size_matched_distinct_rows_c_over_b": [d2.0 / ratio_f.powi(2), d3.0 / ratio_f.powi(3)],
                "size_matched_rows_per_unknown_c_over_b": [rpu2.0 / ratio_f, rpu3.0 / ratio_f.powi(2)],
                "wall_s_contended": started.elapsed().as_secs_f64(),
            });
            eprintln!(
                "census {b} seed {s} p={} r=2^{:.1}: |F|={} U a/b/c={}/{}/{} n2 b {:.3} c {:.3} (ratio {:.3} [{:.3},{:.3}]) n3 b {:.2} c {:.2} (ratio {:.3} [{:.3},{:.3}]) rows/unknown c/b m2 {:.2} m3 {:.2} [{:.1}s]",
                inst.p, (inst.r as f64).log2(), bases.fa.points.len(),
                ua, ub, uc, mean(&b2p), mean(&c2p), d2.0, d2.1, d2.2,
                mean(&b3p), mean(&c3p), d3.0, d3.1, d3.2, rpu2.0, rpu3.0,
                started.elapsed().as_secs_f64()
            );
            rows.push(row);
        }
    }
    rows
}

// ── Pipeline: the line oracle, m = 2 ────────────────────────────────

fn ctx_for<'a>(c: &'a Curve, gen: Pt, target: Pt, inst: &Inst) -> InstanceCtx<'a, Curve> {
    InstanceCtx {
        group: c,
        generator: gen,
        target,
        r: inst.r,
        cofactor: inst.h,
        group_order: inst.n,
        name: format!("j0-fp2-p{}", inst.p),
        field_degree: Some(2),
    }
}

struct CallStats {
    hit: Vec<f64>,
    muls: Vec<f64>,
    invs: Vec<f64>,
    exhaustive_agree: u64,
    exhaustive_disagree: u64,
    unliftable: u64,
    degenerate: u64,
}

/// Decompose every target with the line oracle; `map` sends a target on
/// `E` to the curve the base lives on (identity or `φ`), and its cost is
/// charged to the call.
fn line_calls(
    c: &Curve,
    gen: Pt,
    inst: &Inst,
    fb: &Fb,
    targets: &[Pt],
    map: &dyn Fn(Pt) -> Pt,
    check: usize,
) -> CallStats {
    let ctx = ctx_for(c, gen, gen, inst);
    let mut line: LineOracle<Fp2> = LineOracle::new(c.f, c.f.one());
    let mut ops = GroupOps::default();
    line.prepare(&ctx, fb, &Params::default(), &mut ops)
        .expect("line base");
    let mut ctr = OracleCounters::default();
    let mut st = CallStats {
        hit: vec![],
        muls: vec![],
        invs: vec![],
        exhaustive_agree: 0,
        exhaustive_disagree: 0,
        unliftable: 0,
        degenerate: 0,
    };
    for (k, &rt) in targets.iter().enumerate() {
        take_field_counters();
        let s0 = line.solver_totals().unwrap().ops;
        let pt = map(rt);
        let res = line.decompose(&ctx, fb, &mut ops, &mut ctr, pt);
        let (m, i) = take_field_counters();
        let s1 = line.solver_totals().unwrap().ops;
        st.hit.push(if res.is_some() { 1.0 } else { 0.0 });
        st.muls.push((m + s1 - s0) as f64);
        st.invs.push(i as f64);
        if let Some(v) = &res {
            let sum = c.add(&mut ops, fb.points[v[0]], fb.points[v[1]]);
            assert_eq!(sum, pt, "a decomposition that does not sum to its target");
        }
        if k < check {
            let exists = fb.points.iter().any(|&q| {
                let s = c.add(&mut ops, pt, c.neg(q));
                fb.index_of_key(c.key(&s)).is_some()
            });
            if exists == res.is_some() {
                st.exhaustive_agree += 1;
            } else {
                st.exhaustive_disagree += 1;
            }
        }
    }
    st.unliftable = line.unliftable;
    st.degenerate = line.degenerate;
    st
}

fn pipeline(bits: &[u32], seeds: u64, ntargets: usize, torsion: bool) -> Vec<Value> {
    let mut rows = Vec::new();
    for &b in bits {
        for s in 1..=seeds {
            let started = Instant::now();
            let inst = match build(b, s, torsion) {
                Ok(i) => i,
                Err(e) => {
                    eprintln!("pipeline {b} seed {s}: {e}");
                    continue;
                }
            };
            let bases = build_bases(&inst, false);
            let targets = random_targets(&inst, ntargets, s);
            let check = ntargets.min(if b <= 12 { ntargets } else { 200 });
            let id = |x: Pt| x;
            let (f, bb) = (inst.f, inst.b);
            let tors = inst.torsion;
            let ph = move |x: Pt| if tors { phi(&f, bb, x) } else { x };
            let ab = line_calls(&inst.e, inst.g, &inst, &bases.fb, &targets, &id, check);
            let cc = line_calls(&inst.ep, inst.gp, &inst, &bases.cz, &targets, &ph, check);
            let (ua, ub, uc) = (
                bases.fa.columns as f64,
                bases.fb.columns as f64,
                bases.cz.columns as f64,
            );
            // cost per relation = Σmuls / Σhits; total = U × that.
            let tot = |u: f64, st: &CallStats| {
                u * st.muls.iter().sum::<f64>() / st.hit.iter().sum::<f64>()
            };
            let (ta, tb, tc) = (tot(ua, &ab), tot(ub, &ab), tot(uc, &cc));
            // Paired bootstrap of (c)/(b) total: (U_c Σm_c/Σh_c)/(U_b Σm_b/Σh_b).
            let n = targets.len();
            let mut rng = StdRng::seed_from_u64(s ^ 0xC0B0);
            let mut reps: Vec<f64> = (0..2000)
                .map(|_| {
                    let (mut mb, mut hb, mut mc, mut hc) = (0.0, 0.0, 0.0, 0.0);
                    for _ in 0..n {
                        let k = rng.gen_range(0..n);
                        mb += ab.muls[k];
                        hb += ab.hit[k];
                        mc += cc.muls[k];
                        hc += cc.hit[k];
                    }
                    (uc * mc / hc) / (ub * mb / hb)
                })
                .collect();
            reps.sort_by(|x, y| x.partial_cmp(y).unwrap());
            let hb = mean(&ab.hit);
            let hc = mean(&cc.hit);
            let row = json!({
                "experiment": "pipeline", "bits": b, "seed": s,
                "instance": inst_json(&inst), "bases": bases.info, "targets": ntargets,
                "oracle": "line_oracle::LineOracle<Fp2>, line c = 1 (x in F_p), two summands",
                "U_a": ua, "U_b": ub, "U_c": uc,
                "hit_rate_ab": hb, "hit_rate_c": hc,
                "muls_per_call_ab": mean(&ab.muls), "muls_per_call_c": mean(&cc.muls),
                "invs_per_call_ab": mean(&ab.invs), "invs_per_call_c": mean(&cc.invs),
                "cost_per_relation_ab": mean(&ab.muls) / hb,
                "cost_per_relation_c": mean(&cc.muls) / hc,
                "total_a": ta, "total_b": tb, "total_c": tc,
                "total_b_over_a": tb / ta, "total_c_over_b": tc / tb,
                "total_c_over_b_ci95": [reps[50], reps[1949]],
                "base_size_ratio_Fc_over_F": 3.0 * bases.cz.points.len() as f64 / bases.fa.points.len() as f64,
                "size_matched_total_c_over_b": (tc / tb) * 3.0 * bases.cz.points.len() as f64 / bases.fa.points.len() as f64,
                "law_total_b_over_a": 1.0/3.0, "law_total_c_over_b": 1.0,
                "c_arm_valid": inst.torsion,
                "agreement_ab": [ab.exhaustive_agree, ab.exhaustive_disagree],
                "agreement_c": [cc.exhaustive_agree, cc.exhaustive_disagree],
                "unliftable_ab": ab.unliftable, "unliftable_c": cc.unliftable,
                "degenerate_ab": ab.degenerate, "degenerate_c": cc.degenerate,
                "wall_s_contended": started.elapsed().as_secs_f64(),
            });
            eprintln!(
                "pipeline {b} seed {s} p={} r=2^{:.1}: U {}/{}/{} hit b {:.3} c {:.3} muls/call b {:.0} c {:.0} total b/a {:.3} c/b {:.3} [{:.3},{:.3}] agree {}:{} {}:{} [{:.1}s]",
                inst.p, (inst.r as f64).log2(), ua, ub, uc, hb, hc,
                mean(&ab.muls), mean(&cc.muls), tb / ta, tc / tb, reps[50], reps[1949],
                ab.exhaustive_agree, ab.exhaustive_disagree, cc.exhaustive_agree, cc.exhaustive_disagree,
                started.elapsed().as_secs_f64()
            );
            rows.push(row);
        }
    }
    rows
}

// ── Stream: full rank, logarithm recovered ──────────────────────────

fn stream(bits: &[u32], seeds: u64, max_trials: u64, torsion: bool) -> Vec<Value> {
    let mut rows = Vec::new();
    for &b in bits {
        for s in 1..=seeds {
            let started = Instant::now();
            let inst = match build(b, s, torsion) {
                Ok(i) => i,
                Err(e) => {
                    eprintln!("stream {b} seed {s}: {e}");
                    continue;
                }
            };
            let bases = build_bases(&inst, false);
            let mut rng = StdRng::seed_from_u64(s ^ 0x5EED);
            let planted = rng.gen_range(2..inst.r);
            let mut ops = GroupOps::default();
            let q = inst.e.mul(&mut ops, inst.g, planted);
            let qp = phi(&inst.f, inst.b, q);
            let run = |c: &Curve, gen: Pt, tgt: Pt, folded: &Fb, control: &Fb| {
                let ctx = ctx_for(c, gen, tgt, &inst);
                let mut line: LineOracle<Fp2> = LineOracle::new(c.f, c.f.one());
                line.prepare(&ctx, folded, &Params::default(), &mut GroupOps::default())
                    .unwrap();
                take_field_counters();
                let rep = full_rank_stream(
                    c,
                    gen,
                    tgt,
                    inst.r,
                    inst.h,
                    planted,
                    folded,
                    control,
                    s,
                    max_trials,
                    2,
                    None,
                    |o, ctr, pt| line.decompose(&ctx, folded, o, ctr, pt),
                )
                .unwrap();
                let (m, i) = take_field_counters();
                let solver = line.solver_totals().unwrap().ops;
                (rep, m + solver, i)
            };
            let (r_ab, m_ab, i_ab) = run(&inst.e, inst.g, q, &bases.fb, &bases.fa);
            let (r_c, m_c, i_c) = if inst.torsion {
                run(&inst.ep, inst.gp, qp, &bases.cz, &bases.cn)
            } else {
                run(&inst.e, inst.g, q, &bases.cz, &bases.cn)
            };
            let arm = |a: &crypto_lib::cryptanalysis::glv_invariant_experiments::ArmReport| {
                json!({
                    "columns": a.columns, "points": a.points,
                    "square_relations": a.square_relations, "square_trials": a.square_trials,
                    "first_pin_relations": a.first_pin_relations, "first_pin_trials": a.first_pin_trials,
                    "rank_final": a.rank_final, "deficiency_total": a.deficiency_total,
                    "recovered": a.recovered, "verified": a.verified,
                })
            };
            let sq = |a: &crypto_lib::cryptanalysis::glv_invariant_experiments::ArmReport| {
                a.square_trials.map(|t| t as f64)
            };
            let tb = sq(&r_ab.folded);
            let tc = sq(&r_c.folded);
            let row = json!({
                "experiment": "stream", "bits": b, "seed": s,
                "instance": inst_json(&inst), "planted": planted,
                "ab_stream": {
                    "trials": r_ab.trials, "relations": r_ab.relations, "hit_rate": r_ab.hit_rate,
                    "field_muls_total": m_ab, "field_invs_total": i_ab,
                    "b": arm(&r_ab.folded), "a": arm(&r_ab.control),
                    "square_ratio_a_over_b": r_ab.square_ratio, "column_ratio": r_ab.column_ratio,
                },
                "c_stream": {
                    "trials": r_c.trials, "relations": r_c.relations, "hit_rate": r_c.hit_rate,
                    "field_muls_total": m_c, "field_invs_total": i_c,
                    "c_zeta": arm(&r_c.folded), "c_neg": arm(&r_c.control),
                    "square_ratio_neg_over_zeta": r_c.square_ratio, "column_ratio": r_c.column_ratio,
                },
                "c_arm_valid": inst.torsion,
                "square_trials_c_over_b": match (tc, tb) { (Some(c), Some(b)) => Some(c / b), _ => None },
                "law_square_trials_c_over_b": 1.0,
                "wall_s_contended": started.elapsed().as_secs_f64(),
            });
            eprintln!(
                "stream {b} seed {s} p={}: ab trials {} rel {} | b sq rel {:?} trials {:?} ok {} | a sq rel {:?} ok {} | c trials {} sq rel {:?} trials {:?} ok {} | cneg sq rel {:?} ok {} [{:.1}s]",
                inst.p, r_ab.trials, r_ab.relations,
                r_ab.folded.square_relations, r_ab.folded.square_trials, r_ab.folded.verified,
                r_ab.control.square_relations, r_ab.control.verified,
                r_c.trials, r_c.folded.square_relations, r_c.folded.square_trials, r_c.folded.verified,
                r_c.control.square_relations, r_c.control.verified,
                started.elapsed().as_secs_f64()
            );
            rows.push(row);
        }
    }
    rows
}

// ── Sweep: a three-summand oracle from the two-summand line oracle ──

struct SweepStats {
    calls: Vec<f64>,
    muls: Vec<f64>,
    rows: Vec<f64>,
    raw: Vec<f64>,
}

/// For each target `R` (on the base's curve): decompose `R − P₁` for every
/// base point `P₁`; collect distinct triples and distinct projected rows.
fn sweep_m3(c: &Curve, gen: Pt, inst: &Inst, fb: &Fb, targets: &[Pt]) -> SweepStats {
    let ctx = ctx_for(c, gen, gen, inst);
    let mut line: LineOracle<Fp2> = LineOracle::new(c.f, c.f.one());
    let mut ops = GroupOps::default();
    line.prepare(&ctx, fb, &Params::default(), &mut ops)
        .unwrap();
    let mut ctr = OracleCounters::default();
    let mut st = SweepStats {
        calls: vec![],
        muls: vec![],
        rows: vec![],
        raw: vec![],
    };
    for &rt in targets {
        take_field_counters();
        let s0 = line.solver_totals().unwrap().ops;
        let mut triples: HashSet<[usize; 3]> = HashSet::new();
        let mut calls = 0u64;
        for (i, &p1) in fb.points.iter().enumerate() {
            let s = c.add(&mut ops, rt, c.neg(p1));
            if s.infinity {
                continue;
            }
            calls += 1;
            if let Some(v) = line.decompose(&ctx, fb, &mut ops, &mut ctr, s) {
                let mut t = [i, v[0], v[1]];
                t.sort_unstable();
                triples.insert(t);
            }
        }
        let (m, _) = take_field_counters();
        let s1 = line.solver_totals().unwrap().ops;
        let rows: HashSet<Vec<(usize, u64)>> =
            triples.iter().map(|t| proj_row(t, fb, inst.r)).collect();
        st.calls.push(calls as f64);
        st.muls.push((m + s1 - s0) as f64);
        st.rows.push(rows.len() as f64);
        st.raw.push(triples.len() as f64);
    }
    st
}

fn sweep(bits: &[u32], seeds: u64, ntargets: usize, torsion: bool) -> Vec<Value> {
    let mut rows = Vec::new();
    for &b in bits {
        for s in 1..=seeds {
            let started = Instant::now();
            let inst = match build(b, s, torsion) {
                Ok(i) => i,
                Err(e) => {
                    eprintln!("sweep {b} seed {s}: {e}");
                    continue;
                }
            };
            let bases = build_bases(&inst, false);
            let targets = random_targets(&inst, ntargets, s ^ 0x33);
            let tp: Vec<Pt> = targets.iter().map(|&x| phi(&inst.f, inst.b, x)).collect();
            let sb = sweep_m3(&inst.e, inst.g, &inst, &bases.fb, &targets);
            let sc = sweep_m3(&inst.ep, inst.gp, &inst, &bases.cz, &tp);
            let (ua, ub, uc) = (
                bases.fa.columns as f64,
                bases.fb.columns as f64,
                bases.cz.columns as f64,
            );
            let sum = |v: &[f64]| v.iter().sum::<f64>();
            let cpr_b = sum(&sb.muls) / sum(&sb.rows);
            let cpr_c = sum(&sc.muls) / sum(&sc.rows);
            // Paired bootstrap of the (c)/(b) total, quotient sweep.
            let n = targets.len();
            let mut rng = StdRng::seed_from_u64(s ^ 0x5133);
            let mut reps: Vec<f64> = (0..2000)
                .map(|_| {
                    let (mut mb, mut rb, mut mc, mut rc) = (0.0, 0.0, 0.0, 0.0);
                    for _ in 0..n {
                        let k = rng.gen_range(0..n);
                        mb += sb.muls[k];
                        rb += sb.rows[k];
                        mc += sc.muls[k];
                        rc += sc.rows[k];
                    }
                    (uc * mc / rc) / (ub * mb / rb)
                })
                .collect();
            reps.sort_by(|x, y| x.partial_cmp(y).unwrap());
            let q = (uc * cpr_c) / (ub * cpr_b);
            let row = json!({
                "experiment": "sweep_m3", "bits": b, "seed": s,
                "instance": inst_json(&inst), "bases": bases.info, "targets": ntargets,
                "oracle": "sweep P1 over the base, LineOracle on R - P1 (one pair per call)",
                "U_a": ua, "U_b": ub, "U_c": uc,
                "calls_per_target_b": mean(&sb.calls), "calls_per_target_c": mean(&sc.calls),
                "muls_per_target_b": mean(&sb.muls), "muls_per_target_c": mean(&sc.muls),
                "distinct_rows_per_target_b": mean(&sb.rows), "distinct_rows_per_target_c": mean(&sc.rows),
                "distinct_triples_per_target_b": mean(&sb.raw), "distinct_triples_per_target_c": mean(&sc.raw),
                "cost_per_relation_b": cpr_b, "cost_per_relation_c_quotient": cpr_c,
                "total_a": ua * cpr_b, "total_b": ub * cpr_b, "total_c_quotient": uc * cpr_c,
                "total_b_over_a": ub / ua,
                "total_c_quotient_over_b": q, "total_c_quotient_over_b_ci95": [reps[50], reps[1949]],
                "total_c_naive_over_b_derived": 3.0 * q,
                "base_size_ratio_Fc_over_F": 3.0 * bases.cz.points.len() as f64 / bases.fa.points.len() as f64,
                "size_matched_total_c_quotient_over_b": q * 3.0 * bases.cz.points.len() as f64 / bases.fa.points.len() as f64,
                "law_total_c_over_b": 3.0,
                "wall_s_contended": started.elapsed().as_secs_f64(),
            });
            eprintln!(
                "sweep {b} seed {s} p={}: U {}/{}/{} calls b {:.0} c {:.0} rows b {:.1} c {:.1} muls b {:.3e} c {:.3e} total c/b quotient {:.3} [{:.3},{:.3}] naive {:.3} [{:.1}s]",
                inst.p, ua, ub, uc, mean(&sb.calls), mean(&sc.calls), mean(&sb.rows), mean(&sc.rows),
                mean(&sb.muls), mean(&sc.muls), q, reps[50], reps[1949], 3.0 * q,
                started.elapsed().as_secs_f64()
            );
            rows.push(row);
        }
    }
    rows
}

// ── Driver ──────────────────────────────────────────────────────────

fn main() {
    let args: Vec<String> = env::args().collect();
    let get = |k: &str| {
        args.iter()
            .position(|a| a == k)
            .and_then(|i| args.get(i + 1).cloned())
    };
    let exp = get("--exp").unwrap_or_else(|| "census".into());
    let bits: Vec<u32> = get("--bits")
        .unwrap_or_else(|| "10".into())
        .split(',')
        .map(|x| x.parse().expect("bits"))
        .collect();
    let seeds: u64 = get("--seeds").map_or(3, |x| x.parse().unwrap());
    let targets: usize = get("--targets").map_or(300, |x| x.parse().unwrap());
    let max_trials: u64 = get("--max-trials").map_or(2_000_000, |x| x.parse().unwrap());
    let family = get("--family").unwrap_or_else(|| "torsion".into());
    let torsion = match family.as_str() {
        "torsion" => true,
        "prime" => false,
        other => panic!("unknown --family {other}; torsion or prime"),
    };
    let started = Instant::now();
    let rows = match exp.as_str() {
        "census" => census(&bits, seeds, targets, torsion),
        "pipeline" => pipeline(&bits, seeds, targets, torsion),
        "stream" => stream(&bits, seeds, max_trials, torsion),
        "sweep" => sweep(&bits, seeds, targets, torsion),
        other => panic!("unknown --exp {other}"),
    };
    let out = json!({
        "lead": "16.4", "experiment": exp, "family": family, "bits": bits, "seeds": seeds,
        "targets": targets, "rows": rows,
        "unit": "F_p multiplications (ext_curve thread-local counter + LineOracle solver count); inversions reported, not priced; wall times are contended",
        "wall_s_contended": started.elapsed().as_secs_f64(),
    });
    let text = serde_json::to_string_pretty(&out).unwrap();
    match get("--json") {
        Some(path) => fs::write(&path, text).expect("write json"),
        None => println!("{text}"),
    }
}

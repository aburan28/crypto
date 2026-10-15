//! Factor-base shape at `m ≥ 3`: goals G1, G2, G3 and G5 of
//! `research/notes/index-calculus/RESEARCH_FACTOR_BASE_SHAPE_SEARCH.md` §7.
//!
//! The instrument is the chained-`S₃` ladder of `RESEARCH_DREG_MEASUREMENT.md`
//! (`koblitz_bench::ladder_measure_m`: the system of
//! `koblitz_groebner::build_decomposition_system`, an exact solution count
//! first, then `solving_degree_from` on unsatisfiable draws only), with the
//! subspace drawn by the shape builders of `factor_base_closure_sweep.rs`
//! instead of uniformly, every target tagged with its trace class
//! `t = Tr(√b/x_R)` (`b = 1`, so `t = Tr(1/x_R)`), and only class-0 targets
//! measured.  The system is built here, summand by summand, so that each
//! summand can carry its own basis (G2) and the field modulus can vary (G3);
//! for one basis and the sparse modulus it is term for term the ladder's
//! system (`selftest` checks the solution count against an exhaustive
//! evaluation of the built equations).
//!
//! Subcommands (one JSON line per measured target on stdout):
//!
//! ```text
//! g1 --cell n:l:dmax --shapes random,random2,gp,fp,coord,subfield --instances K --targets T --seed S
//! g3 --cell n:l:dmax --targets T --seed S            (coord window over the modulus panel)
//! g2 --cell n:l:dmax --twist s --shapes random,fp --instances K --targets T --seed S
//! g5 --cell n:l:dmax --betas R --instances K --targets T --seed S   (m = 2, Möbius pole)
//! selftest
//! score FILE.jsonl …                                (PROTOCOL.md §4 rules)
//! ```
//!
//! Stage diagnostic: refutation degrees on toy fields.  No `S`, no speedup.

use crypto_lib::binary_ecc::{F2mElement, IrreduciblePoly};
use crypto_lib::cryptanalysis::descent_expansion::enumerate_irreducibles;
use crypto_lib::cryptanalysis::koblitz_bench::random_subspace_basis;
use crypto_lib::cryptanalysis::koblitz_groebner::{
    solving_degree_from, sym_semaev_s3, FieldStructure, SymElement,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::find_irreducible_sparse;
use crypto_lib::cryptanalysis::pq_groebner_f2::F2BoolPoly;
use crypto_lib::cryptanalysis::semaev_decomp::Gf2;
use num_bigint::BigUint;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::time::Instant;

// ── field helpers on u64 words ───────────────────────────────────────────

struct Field {
    n: u32,
    irr: IrreduciblePoly,
    gf: Gf2,
    st: FieldStructure,
}

impl Field {
    fn new(irr: IrreduciblePoly) -> Self {
        let n = irr.degree;
        Field {
            n,
            gf: Gf2::new(&irr),
            st: FieldStructure::new(n, &irr),
            irr,
        }
    }
    fn elem(&self, w: u64) -> F2mElement {
        F2mElement::from_biguint(&BigUint::from(w), self.n)
    }
    fn trace(&self, w: u64) -> u64 {
        let mut acc = 0u64;
        let mut t = w;
        for _ in 0..self.n {
            acc ^= t;
            t = self.gf.sqr(t);
        }
        debug_assert!(acc <= 1);
        acc
    }
    fn mask(&self) -> u64 {
        (1u64 << self.n) - 1
    }
    fn modulus(&self) -> String {
        let mut terms: Vec<String> = vec![format!("z^{}", self.irr.degree)];
        let mut low = self.irr.low_terms.clone();
        low.sort_unstable_by(|a, b| b.cmp(a));
        for k in low {
            terms.push(match k {
                0 => "1".into(),
                1 => "z".into(),
                k => format!("z^{k}"),
            });
        }
        terms.join("+")
    }
}

fn rank(words: &[u64]) -> usize {
    let mut echelon: Vec<u64> = Vec::new();
    for &v in words {
        let r = echelon.iter().fold(v, |r, &e| r.min(r ^ e));
        if r != 0 {
            echelon.push(r);
            echelon.sort_unstable_by(|a, b| b.cmp(a));
        }
    }
    echelon.len()
}

fn span(words: &[u64]) -> Vec<u64> {
    (0..1u64 << words.len())
        .map(|c| {
            (0..words.len())
                .filter(|t| (c >> t) & 1 == 1)
                .fold(0, |acc, t| acc ^ words[t])
        })
        .collect()
}

fn rand_not_01(rng: &mut StdRng, f: &Field) -> u64 {
    loop {
        let w = rng.gen::<u64>() & f.mask();
        if w > 1 {
            return w;
        }
    }
}

// ── shapes (builders of factor_base_closure_sweep.rs, on words) ─────────

/// One draw of a shape: `l` independent words, or `None` if the shape does
/// not exist at `(n, l)`.
fn draw_shape(name: &str, f: &Field, l: usize, rng: &mut StdRng) -> Option<(Vec<u64>, String)> {
    let n = f.n;
    for _ in 0..64 {
        let (cols, detail): (Vec<u64>, String) = match name {
            // The ladder's own random ℓ (`random_subspace_basis`); random2 is
            // a second, independent arm of the same builder (negative control).
            "random" | "random2" => {
                let b = random_subspace_basis(n, l, rng);
                (
                    b.iter()
                        .map(|e| e.raw_bits().first().copied().unwrap_or(0))
                        .collect(),
                    "uniform".into(),
                )
            }
            "coord" => ((0..l).map(|i| 1u64 << i).collect(), "g=z".into()),
            "gp" => {
                let g = rand_not_01(rng, f);
                let mut acc = 1u64;
                let mut cols = Vec::with_capacity(l);
                for _ in 0..l {
                    cols.push(acc);
                    acc = f.gf.mul(acc, g);
                }
                (cols, format!("g={g:#x}"))
            }
            "fp" => {
                let a = rand_not_01(rng, f);
                (
                    (0..l).map(|i| f.gf.sqr_k(a, i as u32)).collect(),
                    format!("alpha={a:#x} Tr={}", f.trace(a)),
                )
            }
            "subfield" => {
                if !(n as usize).is_multiple_of(l) || l as u32 == n {
                    return None;
                }
                // F_{2^l} = ker(x ↦ x^{2^l} + x); greedy basis from the
                // enumerated fixed points (n ≤ 16 here).
                if n > 20 {
                    return None;
                }
                let fixed: Vec<u64> = (1..=f.mask())
                    .filter(|&w| f.gf.sqr_k(w, l as u32) == w)
                    .collect();
                let mut cols = Vec::new();
                for w in fixed {
                    let mut cand = cols.clone();
                    cand.push(w);
                    if rank(&cand) == cand.len() {
                        cols = cand;
                    }
                    if cols.len() == l {
                        break;
                    }
                }
                (cols, format!("F_2^{l}"))
            }
            _ => panic!("unknown shape {name}"),
        };
        if cols.len() == l && rank(&cols) == l {
            return Some((cols, detail));
        }
        if matches!(name, "coord" | "subfield") {
            return None;
        }
    }
    None
}

/// `dim V·V`, `dim V ∩ σV`, and whether `V ⊂ ker Tr`.
fn profile(f: &Field, basis: &[u64]) -> (usize, usize, bool) {
    let mut prods = Vec::new();
    for &a in basis {
        for &b in basis {
            prods.push(f.gf.mul(a, b));
        }
    }
    let vv = rank(&prods);
    let sigma: Vec<u64> = basis.iter().map(|&a| f.gf.sqr(a)).collect();
    let mut both = basis.to_vec();
    both.extend(&sigma);
    let cap = 2 * basis.len() - rank(&both);
    let tz = basis.iter().all(|&a| f.trace(a) == 0);
    (vv, cap, tz)
}

// ── systems and exact counts ─────────────────────────────────────────────

/// The `m = 3` chained system `S₃(x₁, x₂, u)`, `S₃(u, x₃, x_R)` with summand
/// `i` in the span of `bases[i]`.  Layout as `build_decomposition_system`:
/// summand coordinates first, then `u`'s `n` bits.
fn system_m3(f: &Field, bases: [&[u64]; 3], x_r: u64, b: u64) -> (Vec<F2BoolPoly>, usize) {
    let l: usize = bases[0].len();
    assert!(bases.iter().all(|v| v.len() == l));
    let n = f.n;
    let n_vars = 3 * l + n as usize;
    let summands: Vec<SymElement> = (0..3)
        .map(|i| {
            let elems: Vec<F2mElement> = bases[i].iter().map(|&w| f.elem(w)).collect();
            SymElement::from_subspace_vars(&elems, i * l, n, n_vars)
        })
        .collect();
    let u = SymElement::from_free_vars(3 * l, n, n_vars);
    let target = SymElement::constant(&f.elem(x_r), n, n_vars);
    let be = f.elem(b);
    let mut eqs = sym_semaev_s3(&summands[0], &summands[1], &u, &be, &f.st);
    eqs.extend(sym_semaev_s3(&u, &summands[2], &target, &be, &f.st));
    (eqs, n_vars)
}

/// Exact solutions of [`system_m3`]: `(x₁, x₂, x₃, u)` with `x_i` in the
/// span of `bases[i]`, `u ∈ F_{2^n}`.
fn count_m3(f: &Field, bases: [&[u64]; 3], x_r: u64, b: u64) -> u64 {
    let gf = &f.gf;
    let s3 = |a: u64, c: u64, d: u64| {
        let ac = gf.mul(a, c);
        gf.mul(gf.sqr(a ^ c), gf.sqr(d)) ^ gf.mul(ac, d) ^ gf.sqr(ac) ^ b
    };
    let (s1, s2, s3span) = (span(bases[0]), span(bases[1]), span(bases[2]));
    let squares: Vec<u64> = (0..=f.mask()).map(|u| gf.sqr(u)).collect();
    let mut count = 0u64;
    for &x1 in &s1 {
        for &x2 in &s2 {
            let (lead, mid) = (gf.sqr(x1 ^ x2), gf.mul(x1, x2));
            let tail = gf.sqr(mid) ^ b;
            for (u, &uu) in squares.iter().enumerate() {
                let u = u as u64;
                if gf.mul(lead, uu) ^ gf.mul(mid, u) ^ tail != 0 {
                    continue;
                }
                count += s3span.iter().filter(|&&x3| s3(u, x3, x_r) == 0).count() as u64;
            }
        }
    }
    count
}

/// The `m = 2` system in the Möbius coordinate `u = 1/(x + β)` (G5):
/// `F = {P : u(P) ∈ V}`.  With `X_i = β + 1/u_i`,
/// `(u₁u₂)²·S₃(X₁, X₂, x₃) = (u₁² + u₂²)x₃² + x₃·u₁u₂·N + N² + b(u₁u₂)²`,
/// `N = β²u₁u₂ + β(u₁ + u₂) + 1`, written so every term is Boolean degree
/// ≤ 2: `u₁u₂N = β²(u₁u₂)² + β(u₁²u₂ + u₁u₂²) + u₁u₂`,
/// `N² = β⁴(u₁u₂)² + β²(u₁² + u₂²) + 1`.  `beta = None` is the plain
/// subspace in `x` (`S₃(x₁, x₂, x₃)`), the reference arm.
fn system_m2(
    f: &Field,
    basis: &[u64],
    beta: Option<u64>,
    x3: u64,
    b: u64,
) -> (Vec<F2BoolPoly>, usize) {
    let l = basis.len();
    let n = f.n;
    let n_vars = 2 * l;
    let elems: Vec<F2mElement> = basis.iter().map(|&w| f.elem(w)).collect();
    let u1 = SymElement::from_subspace_vars(&elems, 0, n, n_vars);
    let u2 = SymElement::from_subspace_vars(&elems, l, n, n_vars);
    let be = f.elem(b);
    let Some(beta) = beta else {
        let t = SymElement::constant(&f.elem(x3), n, n_vars);
        return (sym_semaev_s3(&u1, &u2, &t, &be, &f.st), n_vars);
    };
    let st = &f.st;
    let c = |w: u64| f.elem(w);
    let u1s = u1.square(st);
    let u2s = u2.square(st);
    let p = u1.mul(&u2, st); // u1 u2
    let ps = p.square(st); // (u1 u2)^2
    let mixed = u1s.mul(&u2, st).add(&u1.mul(&u2s, st)); // u1²u2 + u1u2²
    let sum_sq = u1s.add(&u2s);
    let gf = &f.gf;
    let b2 = gf.sqr(beta);
    let b4 = gf.sqr(b2);
    let x3s = gf.sqr(x3);
    // x3 * [β²(u1u2)² + β(u1²u2+u1u2²) + u1u2]
    let un = ps
        .mul_constant(&c(b2), st)
        .add(&mixed.mul_constant(&c(beta), st))
        .add(&p);
    let mut acc = sum_sq.mul_constant(&c(x3s), st);
    acc = acc.add(&un.mul_constant(&c(x3), st));
    // N² = β⁴(u1u2)² + β²(u1²+u2²) + 1
    acc = acc.add(&ps.mul_constant(&c(b4), st));
    acc = acc.add(&sum_sq.mul_constant(&c(b2), st));
    acc = acc.add(&SymElement::constant(&c(1), n, n_vars));
    acc = acc.add(&ps.mul_constant(&c(b), st));
    (acc.coords, n_vars)
}

/// Exact solutions of [`system_m2`] by enumeration of `V²` (including the
/// `u = 0` points, which the polynomial system admits).
fn count_m2(f: &Field, basis: &[u64], beta: Option<u64>, x3: u64, b: u64) -> u64 {
    let gf = &f.gf;
    let s = span(basis);
    let x3s = gf.sqr(x3);
    let mut count = 0u64;
    for &u1 in &s {
        for &u2 in &s {
            let v = match beta {
                None => {
                    let ac = gf.mul(u1, u2);
                    gf.mul(gf.sqr(u1 ^ u2), x3s) ^ gf.mul(ac, x3) ^ gf.sqr(ac) ^ b
                }
                Some(be) => {
                    let p = gf.mul(u1, u2);
                    let nn = gf.mul(gf.sqr(be), p) ^ gf.mul(be, u1 ^ u2) ^ 1;
                    gf.mul(gf.sqr(u1 ^ u2), x3s)
                        ^ gf.mul(x3, gf.mul(p, nn))
                        ^ gf.sqr(nn)
                        ^ gf.mul(b, gf.sqr(p))
                }
            };
            if v == 0 {
                count += 1;
            }
        }
    }
    count
}

// ── measurement ──────────────────────────────────────────────────────────

/// `(json, degree, censored)`: resolved at `degree`, or not resolved on a
/// fully built degree-`d_max` matrix (`censored`, degree is the lower bound
/// `d_max + 1`), or the caps stopped it (`degree = None`).
fn measure(eqs: &[F2BoolPoly], n_vars: usize, d_max: u32) -> (String, Option<u32>, bool) {
    let (d, profs) = solving_degree_from(eqs, n_vars, 1, d_max);
    let built = profs.last().map(|p| p.degree);
    match d {
        Some(degree) => {
            let refuted = profs.last().map(|p| p.refuted).unwrap_or(false);
            (
                format!(r#"{{"kind":"resolved","degree":{degree},"refuted":{refuted}}}"#),
                Some(degree),
                false,
            )
        }
        None if built == Some(d_max) => (
            format!(r#"{{"kind":"at_least","degree":{}}}"#, d_max + 1),
            Some(d_max + 1),
            true,
        ),
        None => (
            format!(
                r#"{{"kind":"caps_hit","built":{}}}"#,
                built.map_or("null".into(), |b| b.to_string())
            ),
            None,
            false,
        ),
    }
}

fn seed_of(seed: u64, parts: &[u64]) -> u64 {
    // SplitMix64 over the parts: stable, order-sensitive.
    let mut h = seed ^ 0x9E37_79B9_7F4A_7C15;
    for &p in parts {
        h = h.wrapping_add(p).wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = h;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        h = z ^ (z >> 31);
    }
    h
}

fn shape_id(name: &str) -> u64 {
    name.bytes()
        .fold(0u64, |h, c| h.wrapping_mul(131).wrapping_add(c as u64))
}

struct Args(Vec<String>);
impl Args {
    fn get(&self, k: &str) -> Option<String> {
        self.0
            .iter()
            .position(|a| a == k)
            .and_then(|i| self.0.get(i + 1))
            .cloned()
    }
    fn num(&self, k: &str, d: u64) -> u64 {
        self.get(k).and_then(|v| v.parse().ok()).unwrap_or(d)
    }
    fn cell(&self) -> (u32, usize, u32) {
        let c = self.get("--cell").expect("--cell n:l:dmax");
        let p: Vec<u64> = c.split(':').map(|x| x.parse().unwrap()).collect();
        (p[0] as u32, p[1] as usize, p[2] as u32)
    }
    fn list(&self, k: &str, d: &str) -> Vec<String> {
        self.get(k)
            .unwrap_or_else(|| d.into())
            .split(',')
            .map(String::from)
            .collect()
    }
}

/// Draw targets for one basis triple until `targets` class-0 unsatisfiable
/// ones are measured; prints one line per measured target and a summary
/// line with the draw accounting.  `tag` is the JSON prefix (cell, arm …).
#[allow(clippy::too_many_arguments)]
fn run_targets_m3(
    f: &Field,
    bases: [&[u64]; 3],
    d_max: u32,
    targets: usize,
    rng: &mut StdRng,
    tag: &str,
    max_draws: usize,
) {
    let b = 1u64;
    let (mut drawn, mut class1, mut sat0, mut measured) = (0usize, 0usize, 0usize, 0usize);
    // Targets are drawn without replacement: at n ≤ 13 a with-replacement
    // stream repeats targets, and a repeat is not an independent draw.
    let mut seen = std::collections::HashSet::new();
    while measured < targets && drawn < max_draws && (seen.len() as u64) < f.mask() {
        let x_r = loop {
            let w = rng.gen::<u64>() & f.mask();
            if w != 0 && seen.insert(w) {
                break w;
            }
        };
        drawn += 1;
        // t = Tr(√b / x_R) = Tr(b / x_R²) = Tr(1/x_R) at b = 1.
        let class = f.trace(f.gf.inv(x_r));
        if class == 1 {
            class1 += 1;
            continue;
        }
        let sols = count_m3(f, bases, x_r, b);
        if sols > 0 {
            sat0 += 1;
            continue;
        }
        let started = Instant::now();
        let (eqs, n_vars) = system_m3(f, bases, x_r, b);
        let (oj, _, _) = measure(&eqs, n_vars, d_max);
        println!(
            r#"{{{tag},"x_r":{x_r},"class":0,"solutions":0,"n_vars":{n_vars},"n_eqs":{},"outcome":{oj},"secs":{:.3}}}"#,
            eqs.len(),
            started.elapsed().as_secs_f64()
        );
        measured += 1;
    }
    println!(
        r#"{{{tag},"summary":true,"drawn":{drawn},"class1_skipped":{class1},"class0_satisfiable":{sat0},"measured":{measured}}}"#
    );
}

fn g1(a: &Args) {
    let (n, l, d_max) = a.cell();
    let seed = a.num("--seed", 20261015);
    let instances = a.num("--instances", 2) as usize;
    let targets = a.num("--targets", 2) as usize;
    let max_draws = a.num("--max-draws", 4096) as usize;
    let f = Field::new(find_irreducible_sparse(n).expect("field"));
    for shape in a.list("--shapes", "random,random2,gp,fp,coord") {
        let inst = if matches!(shape.as_str(), "coord" | "subfield") {
            1
        } else {
            instances
        };
        let per_targets = if inst == 1 {
            targets * instances
        } else {
            targets
        };
        for k in 0..inst {
            let mut rng = StdRng::seed_from_u64(seed_of(
                seed,
                &[n as u64, l as u64, shape_id(&shape), k as u64],
            ));
            let Some((basis, detail)) = draw_shape(&shape, &f, l, &mut rng) else {
                println!(
                    r#"{{"goal":"G1","cell":"n{n}l{l}","shape":"{shape}","instance":{k},"absent":true}}"#
                );
                continue;
            };
            let (vv, cap, tz) = profile(&f, &basis);
            let tag = format!(
                r#""goal":"G1","cell":"n{n}l{l}","n":{n},"l":{l},"d_max":{d_max},"modulus":"{}","shape":"{shape}","instance":{k},"detail":"{detail}","basis":{basis:?},"dim_VV":{vv},"dim_V_cap_sigmaV":{cap},"trace_zero":{tz}"#,
                f.modulus()
            );
            run_targets_m3(
                &f,
                [&basis, &basis, &basis],
                d_max,
                per_targets,
                &mut rng,
                &tag,
                max_draws,
            );
        }
    }
}

/// G3's modulus panel: every trinomial, the lightest pentanomial and the two
/// heaviest moduli in `enumerate_irreducibles` order (the closure sweep's
/// panel), or every irreducible when there are at most `all_below` of them.
fn panel(n: u32, all_below: usize) -> Vec<IrreduciblePoly> {
    let all = enumerate_irreducibles(n, 1 << 20);
    if all.len() <= all_below {
        return all;
    }
    let weight = |p: &IrreduciblePoly| p.low_terms.len() + 1;
    let mut out: Vec<IrreduciblePoly> = all.iter().filter(|p| weight(p) == 3).cloned().collect();
    if let Some(p) = all.iter().find(|p| weight(p) == 5) {
        out.push(p.clone());
    }
    let wmax = all.iter().map(weight).max().unwrap();
    out.extend(all.iter().filter(|p| weight(p) == wmax).take(2).cloned());
    out
}

fn g3(a: &Args) {
    let (n, l, d_max) = a.cell();
    let seed = a.num("--seed", 20261015);
    let targets = a.num("--targets", 4) as usize;
    let max_draws = a.num("--max-draws", 4096) as usize;
    let all_below = a.num("--all-below", 8) as usize;
    for (i, irr) in panel(n, all_below).into_iter().enumerate() {
        let f = Field::new(irr);
        let basis: Vec<u64> = (0..l).map(|i| 1u64 << i).collect();
        let (vv, cap, tz) = profile(&f, &basis);
        let mut rng = StdRng::seed_from_u64(seed_of(seed, &[n as u64, l as u64, 3, i as u64]));
        let tag = format!(
            r#""goal":"G3","cell":"n{n}l{l}","n":{n},"l":{l},"d_max":{d_max},"modulus":"{}","weight":{},"shape":"coord","instance":{i},"basis":{basis:?},"dim_VV":{vv},"dim_V_cap_sigmaV":{cap},"trace_zero":{tz}"#,
            f.modulus(),
            f.irr.low_terms.len() + 1
        );
        run_targets_m3(
            &f,
            [&basis, &basis, &basis],
            d_max,
            targets,
            &mut rng,
            &tag,
            max_draws,
        );
    }
}

fn g2(a: &Args) {
    let (n, l, d_max) = a.cell();
    let seed = a.num("--seed", 20261015);
    let instances = a.num("--instances", 2) as usize;
    let targets = a.num("--targets", 2) as usize;
    let max_draws = a.num("--max-draws", 4096) as usize;
    let twists: Vec<u32> = a
        .list("--twists", "1,2")
        .iter()
        .map(|s| {
            if s == "l" {
                l as u32
            } else {
                s.parse().unwrap()
            }
        })
        .collect();
    let f = Field::new(find_irreducible_sparse(n).expect("field"));
    for shape in a.list("--shapes", "random,fp") {
        for k in 0..instances {
            let mut rng = StdRng::seed_from_u64(seed_of(
                seed,
                &[n as u64, l as u64, 2, shape_id(&shape), k as u64],
            ));
            let (basis, detail) = draw_shape(&shape, &f, l, &mut rng).expect("shape");
            let target_seed: u64 = rng.gen();
            let tw = |s: u32| -> Vec<u64> { basis.iter().map(|&w| f.gf.sqr_k(w, s)).collect() };
            // Arms: plain chain, then for each s the literal instrument
            // ("second": summand 2 squared s times) and the ordered tuple
            // ("ordered": summand j squared j·s times).  Every arm sees the
            // same target stream.
            let mut arms: Vec<(String, [Vec<u64>; 3])> = vec![(
                "plain".into(),
                [basis.clone(), basis.clone(), basis.clone()],
            )];
            for &s in &twists {
                arms.push((
                    format!("second_s{s}"),
                    [basis.clone(), tw(s), basis.clone()],
                ));
                arms.push((format!("ordered_s{s}"), [basis.clone(), tw(s), tw(2 * s)]));
            }
            for (arm, bs) in &arms {
                let mut trng = StdRng::seed_from_u64(target_seed);
                let tag = format!(
                    r#""goal":"G2","cell":"n{n}l{l}","n":{n},"l":{l},"d_max":{d_max},"shape":"{shape}","instance":{k},"detail":"{detail}","arm":"{arm}","bases":[{:?},{:?},{:?}]"#,
                    bs[0], bs[1], bs[2]
                );
                run_targets_m3(
                    &f,
                    [&bs[0], &bs[1], &bs[2]],
                    d_max,
                    targets,
                    &mut trng,
                    &tag,
                    max_draws,
                );
            }
        }
    }
}

fn g5(a: &Args) {
    let (n, l, d_max) = a.cell();
    let seed = a.num("--seed", 20261015);
    let instances = a.num("--instances", 2) as usize;
    let targets = a.num("--targets", 4) as usize;
    let n_betas = a.num("--betas", 3) as usize;
    let max_draws = a.num("--max-draws", 4096) as usize;
    let f = Field::new(find_irreducible_sparse(n).expect("field"));
    let b = 1u64;
    for k in 0..instances {
        let mut rng = StdRng::seed_from_u64(seed_of(seed, &[n as u64, l as u64, 5, k as u64]));
        let (basis, _) = draw_shape("random", &f, l, &mut rng).expect("basis");
        let mut arms: Vec<(String, Option<u64>)> = vec![
            ("x_subspace".into(), None),
            ("beta0".into(), Some(0)),
            ("beta1".into(), Some(1)),
        ];
        for r in 0..n_betas {
            let be = rand_not_01(&mut rng, &f);
            arms.push((format!("beta_rand{r}"), Some(be)));
        }
        let target_seed: u64 = rng.gen();
        for (arm, beta) in &arms {
            let mut trng = StdRng::seed_from_u64(target_seed);
            let (mut drawn, mut class1, mut sat0, mut measured) = (0usize, 0usize, 0usize, 0usize);
            let tag = format!(
                r#""goal":"G5","cell":"n{n}l{l}","n":{n},"l":{l},"d_max":{d_max},"instance":{k},"basis":{basis:?},"arm":"{arm}","beta":{}"#,
                beta.map_or("null".into(), |x| x.to_string())
            );
            let mut seen = std::collections::HashSet::new();
            while measured < targets && drawn < max_draws && (seen.len() as u64) < f.mask() {
                let x3 = loop {
                    let w = trng.gen::<u64>() & f.mask();
                    if w != 0 && seen.insert(w) {
                        break w;
                    }
                };
                drawn += 1;
                if f.trace(f.gf.inv(x3)) == 1 {
                    class1 += 1;
                    continue;
                }
                if count_m2(&f, &basis, *beta, x3, b) > 0 {
                    sat0 += 1;
                    continue;
                }
                let started = Instant::now();
                let (eqs, n_vars) = system_m2(&f, &basis, *beta, x3, b);
                let (oj, _, _) = measure(&eqs, n_vars, d_max);
                println!(
                    r#"{{{tag},"x3":{x3},"class":0,"solutions":0,"n_vars":{n_vars},"outcome":{oj},"secs":{:.3}}}"#,
                    started.elapsed().as_secs_f64()
                );
                measured += 1;
            }
            println!(
                r#"{{{tag},"summary":true,"drawn":{drawn},"class1_skipped":{class1},"class0_satisfiable":{sat0},"measured":{measured}}}"#
            );
        }
    }
}

/// Exact counts against exhaustive evaluation of the built equations, and
/// the system against `build_decomposition_system` for one basis.
fn selftest() {
    use crypto_lib::cryptanalysis::koblitz_groebner::build_decomposition_system;
    let mut rng = StdRng::seed_from_u64(7);
    let mut checked = 0;
    for &(n, l) in &[(5u32, 2usize), (7, 2)] {
        let f = Field::new(find_irreducible_sparse(n).unwrap());
        for _ in 0..6 {
            let v = draw_shape("random", &f, l, &mut rng).unwrap().0;
            let w = draw_shape("fp", &f, l, &mut rng).unwrap().0;
            let x_r = rand_not_01(&mut rng, &f);
            for bases in [[&v, &v, &v], [&v, &w, &v]] {
                let (eqs, nv) = system_m3(&f, [bases[0], bases[1], bases[2]], x_r, 1);
                let ex = (0..1u64 << nv)
                    .filter(|&p| eqs.iter().all(|e| e.eval(p) == 0))
                    .count() as u64;
                let c = count_m3(&f, [bases[0], bases[1], bases[2]], x_r, 1);
                assert_eq!(ex, c, "m=3 count n={n} l={l}");
                checked += 1;
            }
            // Same system as the ladder's for one basis.
            let elems: Vec<F2mElement> = v.iter().map(|&w| f.elem(w)).collect();
            let sys =
                build_decomposition_system(&elems, &f.elem(x_r), &f.elem(1), 3, &f.st).unwrap();
            let (eqs, _) = system_m3(&f, [&v, &v, &v], x_r, 1);
            assert_eq!(sys.equations, eqs, "ladder system n={n}");
        }
    }
    for &(n, l) in &[(7u32, 3usize), (9, 4)] {
        let f = Field::new(find_irreducible_sparse(n).unwrap());
        for _ in 0..8 {
            let v = draw_shape("random", &f, l, &mut rng).unwrap().0;
            let x3 = rand_not_01(&mut rng, &f);
            for beta in [None, Some(0), Some(1), Some(rand_not_01(&mut rng, &f))] {
                let (eqs, nv) = system_m2(&f, &v, beta, x3, 1);
                assert!(eqs
                    .iter()
                    .all(|e| e.terms.iter().all(|m| m.mask.count_ones() <= 2)));
                let ex = (0..1u64 << nv)
                    .filter(|&p| eqs.iter().all(|e| e.eval(p) == 0))
                    .count() as u64;
                assert_eq!(
                    ex,
                    count_m2(&f, &v, beta, x3, 1),
                    "m=2 count n={n} beta={beta:?}"
                );
                checked += 1;
            }
            // The pole arm decomposes x3 exactly when some (β+1/u1, β+1/u2)
            // with u_i ∈ V \ {0} satisfies S₃, or a u = 0 point does.
            let be = rand_not_01(&mut rng, &f);
            let s = span(&v);
            let gf = &f.gf;
            let x3s = gf.sqr(x3);
            let direct = s
                .iter()
                .flat_map(|&a| s.iter().map(move |&c| (a, c)))
                .filter(|&(a, c)| a != 0 && c != 0)
                .filter(|&(a, c)| {
                    let (x1, x2) = (be ^ gf.inv(a), be ^ gf.inv(c));
                    let ac = gf.mul(x1, x2);
                    gf.mul(gf.sqr(x1 ^ x2), x3s) ^ gf.mul(ac, x3) ^ gf.sqr(ac) ^ 1 == 0
                })
                .count() as u64;
            let zero_pts = s
                .iter()
                .flat_map(|&a| s.iter().map(move |&c| (a, c)))
                .filter(|&(a, c)| a == 0 || c == 0)
                .filter(|&(a, c)| {
                    let p = gf.mul(a, c);
                    let nn = gf.mul(gf.sqr(be), p) ^ gf.mul(be, a ^ c) ^ 1;
                    gf.mul(gf.sqr(a ^ c), x3s)
                        ^ gf.mul(x3, gf.mul(p, nn))
                        ^ gf.sqr(nn)
                        ^ gf.mul(1, gf.sqr(p))
                        == 0
                })
                .count() as u64;
            assert_eq!(direct + zero_pts, count_m2(&f, &v, Some(be), x3, 1));
            checked += 1;
        }
    }
    println!(r#"{{"selftest":"ok","checks":{checked}}}"#);
}

// ── scoring (the decision rules of PROTOCOL.md §4) ───────────────────────

#[derive(Default, Clone)]
struct Group {
    degrees: Vec<u32>,
    censored: usize,
    caps: usize,
    sat0: usize,
    unsat_seen: usize,
    inst: std::collections::BTreeMap<u64, Vec<u32>>,
}

impl Group {
    fn mean(&self) -> Option<f64> {
        (!self.degrees.is_empty())
            .then(|| self.degrees.iter().sum::<u32>() as f64 / self.degrees.len() as f64)
    }
    fn inst_spread(&self) -> Option<f64> {
        let ms: Vec<f64> = self
            .inst
            .values()
            .filter(|v| !v.is_empty())
            .map(|v| v.iter().sum::<u32>() as f64 / v.len() as f64)
            .collect();
        (ms.len() >= 2).then(|| {
            ms.iter().cloned().fold(f64::MIN, f64::max)
                - ms.iter().cloned().fold(f64::MAX, f64::min)
        })
    }
}

fn score(files: &[String]) {
    use serde_json::Value;
    use std::collections::BTreeMap;
    let mut groups: BTreeMap<(String, String, String), Group> = BTreeMap::new();
    for path in files {
        let text = std::fs::read_to_string(path).expect("read");
        for line in text.lines() {
            let Ok(v) = serde_json::from_str::<Value>(line) else {
                continue;
            };
            let Some(goal) = v["goal"].as_str() else {
                continue;
            };
            if v["absent"].as_bool() == Some(true) {
                continue;
            }
            let cell = v["cell"].as_str().unwrap_or("").to_string();
            let arms: Vec<String> = match goal {
                "G1" => vec![v["shape"].as_str().unwrap().to_string()],
                "G3" => vec![v["modulus"].as_str().unwrap().to_string()],
                "G2" => vec![format!(
                    "{}/{}",
                    v["shape"].as_str().unwrap(),
                    v["arm"].as_str().unwrap()
                )],
                "G5" => {
                    let a = v["arm"].as_str().unwrap().to_string();
                    if a.starts_with("beta_rand") {
                        vec![a, "beta_rand".into()]
                    } else {
                        vec![a]
                    }
                }
                _ => continue,
            };
            let inst = v["instance"].as_u64().unwrap_or(0);
            for arm in arms {
                let g = groups
                    .entry((goal.to_string(), cell.clone(), arm))
                    .or_default();
                if v["summary"].as_bool() == Some(true) {
                    g.sat0 += v["class0_satisfiable"].as_u64().unwrap_or(0) as usize;
                    g.unsat_seen += v["measured"].as_u64().unwrap_or(0) as usize;
                    continue;
                }
                let o = &v["outcome"];
                match o["kind"].as_str() {
                    Some("resolved") | Some("at_least") => {
                        let d = o["degree"].as_u64().unwrap() as u32;
                        if o["kind"] == "at_least" {
                            g.censored += 1;
                        }
                        g.degrees.push(d);
                        g.inst.entry(inst).or_default().push(d);
                    }
                    _ => g.caps += 1,
                }
            }
        }
    }
    println!("| goal | cell | arm | draws | mean D (censored at bound) | censored | caps | degrees | class-0 sat. fraction | instance-mean spread |");
    println!("|:--|:--|:--|--:|--:|--:|--:|:--|--:|--:|");
    for ((goal, cell, arm), g) in &groups {
        let frac = if g.sat0 + g.unsat_seen > 0 {
            format!("{:.3}", g.sat0 as f64 / (g.sat0 + g.unsat_seen) as f64)
        } else {
            "-".into()
        };
        println!(
            "| {goal} | {cell} | {arm} | {} | {} | {} | {} | {} | {frac} | {} |",
            g.degrees.len(),
            g.mean().map_or("-".into(), |m| format!("{m:.3}")),
            g.censored,
            g.caps,
            g.degrees
                .iter()
                .map(|d| d.to_string())
                .collect::<Vec<_>>()
                .join(" "),
            g.inst_spread().map_or("-".into(), |s| format!("{s:.3}")),
        );
    }
    let get =
        |goal: &str, cell: &str, arm: &str| groups.get(&(goal.into(), cell.into(), arm.into()));
    let cells = |goal: &str| -> Vec<String> {
        let mut c: Vec<String> = groups
            .keys()
            .filter(|k| k.0 == goal)
            .map(|k| k.1.clone())
            .collect();
        c.dedup();
        c
    };
    println!();
    // G1
    let mut established: BTreeMap<String, usize> = BTreeMap::new();
    for cell in cells("G1") {
        let Some(r) = get("G1", &cell, "random") else {
            continue;
        };
        let null = get("G1", &cell, "random2").and_then(|r2| Some(r.mean()? - r2.mean()?));
        let noisy = null.is_some_and(|x| x.abs() >= 0.5);
        for shape in ["gp", "fp", "coord", "subfield"] {
            let Some(x) = get("G1", &cell, shape) else {
                continue;
            };
            let (Some(mr), Some(mx)) = (r.mean(), x.mean()) else {
                continue;
            };
            let gap = mr - mx;
            let status = if x.censored == 0 && x.caps == 0 && gap >= 0.5 && !noisy {
                "established"
            } else if x.censored == 0 && r.censored == 0 && gap < 0.5 {
                "excluded"
            } else if x.censored > 0 && r.censored == 0 && gap < 0.5 {
                "excluded by censoring"
            } else if gap >= 0.5 && noisy {
                "not distinguishable from the negative control"
            } else {
                "undetermined"
            };
            if status == "established" && (shape == "gp" || shape == "fp") {
                *established.entry(shape.into()).or_default() += 1;
            }
            println!(
                "G1 {cell} {shape}: gap = mean(random) - mean({shape}) = {gap:.3} -> {status}; negative control random-random2 = {}",
                null.map_or("-".into(), |x| format!("{x:.3}"))
            );
        }
    }
    let best = established.values().copied().max().unwrap_or(0);
    println!(
        "G1 verdict (scaled cells): {}",
        if best >= 2 {
            "survives"
        } else if best == 0 {
            "killed (no cell with the gap)"
        } else {
            "not killed, not supported"
        }
    );
    // G3
    let mut any_survive = false;
    let mut any_over = false;
    for cell in cells("G3") {
        let ms: Vec<f64> = groups
            .iter()
            .filter(|(k, _)| k.0 == "G3" && k.1 == cell)
            .filter_map(|(_, g)| g.mean())
            .collect();
        if ms.len() < 2 {
            continue;
        }
        let spread = ms.iter().cloned().fold(f64::MIN, f64::max)
            - ms.iter().cloned().fold(f64::MAX, f64::min);
        let noise = get("G1", &cell, "random").and_then(|g| g.inst_spread());
        let over = spread > 0.3 + 1e-9;
        any_over |= over;
        let beyond = over && noise.is_some_and(|n| spread > n);
        any_survive |= beyond;
        println!(
            "G3 {cell}: panel spread {spread:.3} over {} moduli; G1 random instance-mean spread {}",
            ms.len(),
            noise.map_or("-".into(), |x| format!("{x:.3}"))
        );
    }
    println!(
        "G3 verdict: {}",
        if !any_over {
            "killed (flat in the panel to 0.3)"
        } else if any_survive {
            "survives"
        } else {
            "not killed, within noise"
        }
    );
    // G2
    let mut killed = false;
    for ((goal, cell, arm), g) in &groups {
        if goal != "G2" || arm.ends_with("/plain") {
            continue;
        }
        let shape = arm.split('/').next().unwrap();
        let Some(p) = get("G2", cell, &format!("{shape}/plain")) else {
            continue;
        };
        match (g.mean(), p.mean()) {
            (Some(a), Some(b)) => {
                let rise = a - b;
                let censored = g.censored + p.censored > 0;
                if rise >= 1.0 && !censored {
                    killed = true;
                }
                println!(
                    "G2 {cell} {arm}: rise {rise:.3} (censored draws: {})",
                    g.censored + p.censored
                );
            }
            _ => println!("G2 {cell} {arm}: no unsatisfiable class-0 draw to compare"),
        }
    }
    println!(
        "G2 verdict: {}",
        if killed {
            "killed (rise >= 1)"
        } else {
            "survives at these sizes (no rise >= 1)"
        }
    );
    // G5
    let mut over = false;
    for cell in cells("G5") {
        let ms: Vec<f64> = ["beta0", "beta1", "beta_rand"]
            .iter()
            .filter_map(|a| get("G5", &cell, a).and_then(|g| g.mean()))
            .collect();
        let rs: Vec<f64> = ["beta_rand0", "beta_rand1", "beta_rand2"]
            .iter()
            .filter_map(|a| get("G5", &cell, a).and_then(|g| g.mean()))
            .collect();
        let sp = |v: &[f64]| {
            v.iter().cloned().fold(f64::MIN, f64::max) - v.iter().cloned().fold(f64::MAX, f64::min)
        };
        if ms.len() < 2 {
            continue;
        }
        let s = sp(&ms);
        over |= s > 0.3 + 1e-9;
        println!(
            "G5 {cell}: beta-spread {s:.3}; random-beta arms' own spread {}",
            if rs.len() >= 2 {
                format!("{:.3}", sp(&rs))
            } else {
                "-".into()
            }
        );
    }
    println!(
        "G5 verdict: {}",
        if over {
            "not killed (spread > 0.3 somewhere)"
        } else {
            "killed (flat in beta to 0.3)"
        }
    );
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let a = Args(args.clone());
    match args.get(1).map(String::as_str) {
        Some("g1") => g1(&a),
        Some("g2") => g2(&a),
        Some("g3") => g3(&a),
        Some("g5") => g5(&a),
        Some("selftest") => selftest(),
        Some("score") => score(&args[2..]),
        _ => eprintln!("usage: factor_base_shape_ladder g1|g2|g3|g5|selftest …"),
    }
}

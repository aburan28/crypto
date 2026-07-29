//! Thread 22 — exhaustive genus-2 Jacobian census over small proxy primes.
//!
//! Background (see `RESEARCH_AUTOLAB_LOG.md`, entries 2026-07-08, 2026-07-09,
//! 2026-07-21, 2026-07-26 and `RESEARCH_MESTRE_HOWE.md`):
//!
//!   * For p ≡ 1 (mod 6) the j = 0 curves y² = x³ + b fall into exactly six
//!     sextic-twist classes E_0, …, E_5 with traces t_0, …, t_5.
//!   * `howe_sextic_twists_all15.gp` checks Howe's gluing conditions on each of
//!     the C(6,2) = 15 pairs:
//!         (H1) #E_i ≠ #E_j
//!         (H2) x³+b_i and x³+b_j have the same F_p-factorisation pattern
//!              (same Galois module structure on the 2-torsion)
//!         (H3) gcd(#E_i, #E_j) = 1
//!     For secp256k1, 5 of the 15 pairs satisfy (H1)+(H2)+(H3).
//!   * The 2026-07-26 run closed the Rosenhain route for secp256k1 (cubic
//!     residue obstruction) and proposed *Thread 22*: instead of using the
//!     Rosenhain/CHLRS formula, search directly for a genus-2 curve C/F_p with
//!         charpoly(Frob | Jac C) = P_{E_i}(T) · P_{E_j}(T).
//!
//! This program answers Thread 22 **exhaustively** for small proxy primes.
//!
//! Method
//! ------
//! Every genus-2 curve over F_p (p > 5) has a model y² = f(x) with f squarefree
//! of degree 5 or 6.  Up to the substitutions x ↦ ux+v, y ↦ wy the leading
//! coefficient only matters modulo squares and the sub-leading coefficient can
//! be translated away.  So the two families
//!
//!     deg 6:  lead · (x⁶ + a₄x⁴ + a₃x³ + a₂x² + a₁x + a₀),  lead ∈ {1, ν}
//!     deg 5:  lead · (x⁵ + a₃x³ + a₂x² + a₁x + a₀),          lead ∈ {1, ν}
//!
//! (ν a fixed non-residue) cover *every* genus-2 curve over F_p.  For each we
//! compute N₁ = #C(F_p) and N₂ = #C(F_p²) and hence
//!
//!     s₁ = p + 1 − N₁,        s₂ = (N₂ − p² − 1 + s₁²)/2,
//!     charpoly = T⁴ − s₁T³ + s₂T² − p·s₁T + p².
//!
//! The target for the pair (i, j) is
//!
//!     (T²−t_iT+p)(T²−t_jT+p)  ⟹  s₁ = t_i + t_j,  s₂ = t_i·t_j + 2p.
//!
//! Since s₁ is cheap (one pass over F_p) and s₂ is expensive (one pass over
//! F_p²), we prefilter on s₁ and only then compute s₂.  This preserves
//! exhaustiveness: a curve whose s₁ misses every target cannot match.
//!
//! Run:  cargo run --release --example howe_genus2_census -- 7 13 19 31

use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

// ---------------------------------------------------------------------------
// F_p and F_p² = F_p[w]/(w² − ν) tables
// ---------------------------------------------------------------------------

struct Ctx {
    p: u64,
    nu: u64,
    /// Quadratic character on F_p, indexed by value.
    chi: Vec<i64>,
    /// Quadratic character on F_p², indexed by a*p + b for the element a + b·w.
    chi2: Vec<i64>,
    /// xpow[k][x] = x^k in F_p, k = 0..=6.
    xpow: Vec<Vec<u64>>,
    /// x2pow[k][idx] = (idx)^k in F_p², stored as (a, b), k = 0..=6.
    x2pow: Vec<Vec<(u64, u64)>>,
}

fn pow_mod(mut b: u64, mut e: u64, m: u64) -> u64 {
    let mut r = 1u64;
    b %= m;
    while e > 0 {
        if e & 1 == 1 {
            r = r * b % m;
        }
        b = b * b % m;
        e >>= 1;
    }
    r
}

impl Ctx {
    fn new(p: u64) -> Ctx {
        // Quadratic character on F_p.
        let mut chi = vec![-1i64; p as usize];
        chi[0] = 0;
        for x in 1..p {
            chi[(x * x % p) as usize] = 1;
        }
        // Smallest non-residue.
        let nu = (2..p).find(|&v| chi[v as usize] == -1).expect("p > 2");

        // F_p²: element a + b·w with w² = ν, indexed a*p + b.
        let mul2 = |x: (u64, u64), y: (u64, u64)| -> (u64, u64) {
            (
                (x.0 * y.0 + x.1 % p * y.1 % p * nu) % p,
                (x.0 * y.1 + x.1 * y.0) % p,
            )
        };
        let n2 = (p * p) as usize;
        let mut is_sq2 = vec![false; n2];
        for idx in 1..n2 as u64 {
            let e = (idx / p, idx % p);
            let s = mul2(e, e);
            is_sq2[(s.0 * p + s.1) as usize] = true;
        }
        let mut chi2 = vec![-1i64; n2];
        chi2[0] = 0;
        for idx in 1..n2 {
            if is_sq2[idx] {
                chi2[idx] = 1;
            }
        }

        // Power tables.
        let mut xpow = vec![vec![0u64; p as usize]; 7];
        for x in 0..p {
            let mut acc = 1u64;
            for k in 0..7 {
                xpow[k][x as usize] = acc;
                acc = acc * x % p;
            }
        }
        let mut x2pow = vec![vec![(0u64, 0u64); n2]; 7];
        for idx in 0..n2 as u64 {
            let e = (idx / p, idx % p);
            let mut acc = (1u64, 0u64);
            for k in 0..7 {
                x2pow[k][idx as usize] = acc;
                acc = mul2(acc, e);
            }
        }

        Ctx { p, nu, chi, chi2, xpow, x2pow }
    }
}

// ---------------------------------------------------------------------------
// Polynomial helpers over F_p (coefficients low-to-high)
// ---------------------------------------------------------------------------

fn poly_trim(a: &mut Vec<u64>) {
    while a.len() > 1 && *a.last().unwrap() == 0 {
        a.pop();
    }
}

fn poly_rem(a: &[u64], b: &[u64], p: u64) -> Vec<u64> {
    let mut r = a.to_vec();
    poly_trim(&mut r);
    let mut d = b.to_vec();
    poly_trim(&mut d);
    if d.len() == 1 && d[0] == 0 {
        return r;
    }
    let inv_lead = pow_mod(*d.last().unwrap(), p - 2, p);
    while r.len() >= d.len() && !(r.len() == 1 && r[0] == 0) {
        let shift = r.len() - d.len();
        let c = *r.last().unwrap() * inv_lead % p;
        if c != 0 {
            for i in 0..d.len() {
                r[shift + i] = (r[shift + i] + p - d[i] * c % p) % p;
            }
        }
        poly_trim(&mut r);
        if r.len() < d.len() {
            break;
        }
        // Guard against a stuck loop (should not happen with exact arithmetic).
        if c == 0 {
            r.pop();
            poly_trim(&mut r);
        }
    }
    r
}

/// Is `f` squarefree over F_p?  (deg f ≥ 1, p > deg f so f' ≠ 0 for our models)
fn is_squarefree(f: &[u64], p: u64) -> bool {
    let mut df: Vec<u64> = (1..f.len()).map(|i| f[i] * (i as u64) % p).collect();
    if df.is_empty() {
        return false;
    }
    poly_trim(&mut df);
    if df.len() == 1 && df[0] == 0 {
        return false;
    }
    let mut a = f.to_vec();
    let mut b = df;
    poly_trim(&mut a);
    loop {
        let r = poly_rem(&a, &b, p);
        let mut r = r;
        poly_trim(&mut r);
        if r.len() == 1 && r[0] == 0 {
            break;
        }
        a = b;
        b = r;
    }
    b.len() == 1 && b[0] != 0
}

// ---------------------------------------------------------------------------
// Sextic twists of the j = 0 family
// ---------------------------------------------------------------------------

/// Returns (b_k, t_k, roots_of_x3_plus_b, #E_k) for the six sextic twists.
fn sextic_twists(ctx: &Ctx) -> Vec<(u64, i64, usize, i64)> {
    let p = ctx.p;
    // A generator of F_p*.
    let mut fac = vec![];
    let mut m = p - 1;
    let mut d = 2;
    while d * d <= m {
        if m % d == 0 {
            fac.push(d);
            while m % d == 0 {
                m /= d;
            }
        }
        d += 1;
    }
    if m > 1 {
        fac.push(m);
    }
    let g = (2..p)
        .find(|&g| fac.iter().all(|&q| pow_mod(g, (p - 1) / q, p) != 1))
        .expect("generator exists");

    let mut out = vec![];
    for k in 0..6u32 {
        let b = pow_mod(g, k as u64, p);
        // t = -sum_x chi(x^3 + b)
        let mut s = 0i64;
        let mut roots = 0usize;
        for x in 0..p {
            let v = (ctx.xpow[3][x as usize] + b) % p;
            s += ctx.chi[v as usize];
            if v == 0 {
                roots += 1;
            }
        }
        let t = -s;
        let card = p as i64 + 1 - t;
        out.push((b, t, roots, card));
    }
    out
}

fn gcd(a: i64, b: i64) -> i64 {
    if b == 0 { a.abs() } else { gcd(b, a % b) }
}

// ---------------------------------------------------------------------------
// Census
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
struct Hit {
    deg: usize,
    lead: u64,
    coeffs: Vec<u64>, // low-to-high, monic part
}

fn census(ctx: &Ctx, targets: &HashMap<(i64, i64), Vec<(usize, usize)>>) -> HashMap<(i64, i64), Hit> {
    let p = ctx.p;
    let pu = p as usize;
    let s1_targets: HashSet<i64> = targets.keys().map(|&(s1, _)| s1).collect();

    // ---- degree 6 family: lead * (x^6 + a4 x^4 + a3 x^3 + a2 x^2 + a1 x + a0)
    let deg6: Vec<HashMap<(i64, i64), Hit>> = (0..p)
        .into_par_iter()
        .map(|a4| {
            let mut found: HashMap<(i64, i64), Hit> = HashMap::new();
            let mut cnt = vec![0u32; pu];
            let mut cnt2 = vec![0u32; pu * pu];
            for a3 in 0..p {
                for a2 in 0..p {
                    for a1 in 0..p {
                        // g(x) = x^6 + a4 x^4 + a3 x^3 + a2 x^2 + a1 x
                        for c in cnt.iter_mut() {
                            *c = 0;
                        }
                        for x in 0..pu {
                            let v = (ctx.xpow[6][x]
                                + a4 * ctx.xpow[4][x]
                                + a3 * ctx.xpow[3][x]
                                + a2 * ctx.xpow[2][x]
                                + a1 * ctx.xpow[1][x])
                                % p;
                            cnt[v as usize] += 1;
                        }
                        // S(a0) = sum_v cnt[v] * chi[v + a0]
                        let mut survivors: Vec<(u64, i64)> = vec![];
                        for a0 in 0..pu {
                            let mut s = 0i64;
                            for v in 0..pu {
                                let c = cnt[v];
                                if c != 0 {
                                    let w = if v + a0 >= pu { v + a0 - pu } else { v + a0 };
                                    s += c as i64 * ctx.chi[w];
                                }
                            }
                            // lead=1 → s1 = -(S+1);  lead=nu → s1 = S+1
                            let s1a = -(s + 1);
                            let s1b = s + 1;
                            if s1_targets.contains(&s1a) || s1_targets.contains(&s1b) {
                                let f = vec![a0 as u64, a1, a2, a3, a4, 0, 1];
                                if is_squarefree(&f, p) {
                                    survivors.push((a0 as u64, s));
                                }
                            }
                        }
                        if survivors.is_empty() {
                            continue;
                        }
                        // Build the F_p² count table once per (a4,a3,a2,a1).
                        for c in cnt2.iter_mut() {
                            *c = 0;
                        }
                        for idx in 0..pu * pu {
                            let e6 = ctx.x2pow[6][idx];
                            let e4 = ctx.x2pow[4][idx];
                            let e3 = ctx.x2pow[3][idx];
                            let e2 = ctx.x2pow[2][idx];
                            let e1 = ctx.x2pow[1][idx];
                            let ra = (e6.0 + a4 * e4.0 + a3 * e3.0 + a2 * e2.0 + a1 * e1.0) % p;
                            let rb = (e6.1 + a4 * e4.1 + a3 * e3.1 + a2 * e2.1 + a1 * e1.1) % p;
                            cnt2[(ra * p + rb) as usize] += 1;
                        }
                        for (a0, s) in survivors {
                            let mut s2sum = 0i64;
                            for a in 0..pu {
                                let a_sh = {
                                    let t = a + a0 as usize;
                                    if t >= pu { t - pu } else { t }
                                };
                                for b in 0..pu {
                                    let c = cnt2[a * pu + b];
                                    if c != 0 {
                                        s2sum += c as i64 * ctx.chi2[a_sh * pu + b];
                                    }
                                }
                            }
                            let n2 = (p * p) as i64 + s2sum + 2;
                            for &(lead, s1) in
                                [(1u64, -(s + 1)), (ctx.nu, s + 1)].iter()
                            {
                                let s2 = (n2 - (p * p) as i64 - 1 + s1 * s1) / 2;
                                if targets.contains_key(&(s1, s2)) {
                                    found.entry((s1, s2)).or_insert(Hit {
                                        deg: 6,
                                        lead,
                                        coeffs: vec![a0, a1, a2, a3, a4, 0, 1],
                                    });
                                }
                            }
                        }
                    }
                }
            }
            found
        })
        .collect();

    // ---- degree 5 family: lead * (x^5 + a3 x^3 + a2 x^2 + a1 x + a0)
    let deg5: Vec<HashMap<(i64, i64), Hit>> = (0..p)
        .into_par_iter()
        .map(|a3| {
            let mut found: HashMap<(i64, i64), Hit> = HashMap::new();
            let mut cnt = vec![0u32; pu];
            let mut cnt2 = vec![0u32; pu * pu];
            for a2 in 0..p {
                for a1 in 0..p {
                    for c in cnt.iter_mut() {
                        *c = 0;
                    }
                    for x in 0..pu {
                        let v = (ctx.xpow[5][x]
                            + a3 * ctx.xpow[3][x]
                            + a2 * ctx.xpow[2][x]
                            + a1 * ctx.xpow[1][x])
                            % p;
                        cnt[v as usize] += 1;
                    }
                    let mut survivors: Vec<(u64, i64)> = vec![];
                    for a0 in 0..pu {
                        let mut s = 0i64;
                        for v in 0..pu {
                            let c = cnt[v];
                            if c != 0 {
                                let w = if v + a0 >= pu { v + a0 - pu } else { v + a0 };
                                s += c as i64 * ctx.chi[w];
                            }
                        }
                        // lead=1 → s1 = -S;  lead=nu → s1 = S
                        if s1_targets.contains(&(-s)) || s1_targets.contains(&s) {
                            let f = vec![a0 as u64, a1, a2, a3, 0, 1];
                            if is_squarefree(&f, p) {
                                survivors.push((a0 as u64, s));
                            }
                        }
                    }
                    if survivors.is_empty() {
                        continue;
                    }
                    for c in cnt2.iter_mut() {
                        *c = 0;
                    }
                    for idx in 0..pu * pu {
                        let e5 = ctx.x2pow[5][idx];
                        let e3 = ctx.x2pow[3][idx];
                        let e2 = ctx.x2pow[2][idx];
                        let e1 = ctx.x2pow[1][idx];
                        let ra = (e5.0 + a3 * e3.0 + a2 * e2.0 + a1 * e1.0) % p;
                        let rb = (e5.1 + a3 * e3.1 + a2 * e2.1 + a1 * e1.1) % p;
                        cnt2[(ra * p + rb) as usize] += 1;
                    }
                    for (a0, s) in survivors {
                        let mut s2sum = 0i64;
                        for a in 0..pu {
                            let a_sh = {
                                let t = a + a0 as usize;
                                if t >= pu { t - pu } else { t }
                            };
                            for b in 0..pu {
                                let c = cnt2[a * pu + b];
                                if c != 0 {
                                    s2sum += c as i64 * ctx.chi2[a_sh * pu + b];
                                }
                            }
                        }
                        let n2 = (p * p) as i64 + s2sum + 1;
                        for &(lead, s1) in [(1u64, -s), (ctx.nu, s)].iter() {
                            let s2 = (n2 - (p * p) as i64 - 1 + s1 * s1) / 2;
                            if targets.contains_key(&(s1, s2)) {
                                found.entry((s1, s2)).or_insert(Hit {
                                    deg: 5,
                                    lead,
                                    coeffs: vec![a0, a1, a2, a3, 0, 1],
                                });
                            }
                        }
                    }
                }
            }
            found
        })
        .collect();

    let mut all: HashMap<(i64, i64), Hit> = HashMap::new();
    for m in deg6.into_iter().chain(deg5.into_iter()) {
        for (k, v) in m {
            all.entry(k).or_insert(v);
        }
    }
    all
}

/// Full census: every realisable genus-2 Frobenius char poly (s1, s2) over F_p.
/// No s1 prefilter, so this costs a full F_p²-pass per model — small p only.
fn census_all(ctx: &Ctx) -> HashSet<(i64, i64)> {
    let p = ctx.p;
    let pu = p as usize;

    let parts: Vec<HashSet<(i64, i64)>> = (0..p)
        .into_par_iter()
        .map(|a4| {
            let mut out: HashSet<(i64, i64)> = HashSet::new();
            let mut cnt = vec![0u32; pu];
            let mut cnt2 = vec![0u32; pu * pu];
            for a3 in 0..p {
                for a2 in 0..p {
                    for a1 in 0..p {
                        for c in cnt.iter_mut() {
                            *c = 0;
                        }
                        for x in 0..pu {
                            let v = (ctx.xpow[6][x]
                                + a4 * ctx.xpow[4][x]
                                + a3 * ctx.xpow[3][x]
                                + a2 * ctx.xpow[2][x]
                                + a1 * ctx.xpow[1][x])
                                % p;
                            cnt[v as usize] += 1;
                        }
                        for c in cnt2.iter_mut() {
                            *c = 0;
                        }
                        for idx in 0..pu * pu {
                            let e6 = ctx.x2pow[6][idx];
                            let e4 = ctx.x2pow[4][idx];
                            let e3 = ctx.x2pow[3][idx];
                            let e2 = ctx.x2pow[2][idx];
                            let e1 = ctx.x2pow[1][idx];
                            let ra = (e6.0 + a4 * e4.0 + a3 * e3.0 + a2 * e2.0 + a1 * e1.0) % p;
                            let rb = (e6.1 + a4 * e4.1 + a3 * e3.1 + a2 * e2.1 + a1 * e1.1) % p;
                            cnt2[(ra * p + rb) as usize] += 1;
                        }
                        for a0 in 0..pu {
                            let f = vec![a0 as u64, a1, a2, a3, a4, 0, 1];
                            if !is_squarefree(&f, p) {
                                continue;
                            }
                            let mut s = 0i64;
                            for v in 0..pu {
                                let c = cnt[v];
                                if c != 0 {
                                    let w = if v + a0 >= pu { v + a0 - pu } else { v + a0 };
                                    s += c as i64 * ctx.chi[w];
                                }
                            }
                            let mut s2sum = 0i64;
                            for a in 0..pu {
                                let a_sh = {
                                    let t = a + a0;
                                    if t >= pu { t - pu } else { t }
                                };
                                for b in 0..pu {
                                    let c = cnt2[a * pu + b];
                                    if c != 0 {
                                        s2sum += c as i64 * ctx.chi2[a_sh * pu + b];
                                    }
                                }
                            }
                            let n2 = (p * p) as i64 + s2sum + 2;
                            for &s1 in [-(s + 1), s + 1].iter() {
                                let s2 = (n2 - (p * p) as i64 - 1 + s1 * s1) / 2;
                                out.insert((s1, s2));
                            }
                        }
                    }
                }
            }
            out
        })
        .collect();

    // degree-5 models
    let parts5: Vec<HashSet<(i64, i64)>> = (0..p)
        .into_par_iter()
        .map(|a3| {
            let mut out: HashSet<(i64, i64)> = HashSet::new();
            let mut cnt = vec![0u32; pu];
            let mut cnt2 = vec![0u32; pu * pu];
            for a2 in 0..p {
                for a1 in 0..p {
                    for c in cnt.iter_mut() {
                        *c = 0;
                    }
                    for x in 0..pu {
                        let v = (ctx.xpow[5][x]
                            + a3 * ctx.xpow[3][x]
                            + a2 * ctx.xpow[2][x]
                            + a1 * ctx.xpow[1][x])
                            % p;
                        cnt[v as usize] += 1;
                    }
                    for c in cnt2.iter_mut() {
                        *c = 0;
                    }
                    for idx in 0..pu * pu {
                        let e5 = ctx.x2pow[5][idx];
                        let e3 = ctx.x2pow[3][idx];
                        let e2 = ctx.x2pow[2][idx];
                        let e1 = ctx.x2pow[1][idx];
                        let ra = (e5.0 + a3 * e3.0 + a2 * e2.0 + a1 * e1.0) % p;
                        let rb = (e5.1 + a3 * e3.1 + a2 * e2.1 + a1 * e1.1) % p;
                        cnt2[(ra * p + rb) as usize] += 1;
                    }
                    for a0 in 0..pu {
                        let f = vec![a0 as u64, a1, a2, a3, 0, 1];
                        if !is_squarefree(&f, p) {
                            continue;
                        }
                        let mut s = 0i64;
                        for v in 0..pu {
                            let c = cnt[v];
                            if c != 0 {
                                let w = if v + a0 >= pu { v + a0 - pu } else { v + a0 };
                                s += c as i64 * ctx.chi[w];
                            }
                        }
                        let mut s2sum = 0i64;
                        for a in 0..pu {
                            let a_sh = {
                                let t = a + a0;
                                if t >= pu { t - pu } else { t }
                            };
                            for b in 0..pu {
                                let c = cnt2[a * pu + b];
                                if c != 0 {
                                    s2sum += c as i64 * ctx.chi2[a_sh * pu + b];
                                }
                            }
                        }
                        let n2 = (p * p) as i64 + s2sum + 1;
                        for &s1 in [-s, s].iter() {
                            let s2 = (n2 - (p * p) as i64 - 1 + s1 * s1) / 2;
                            out.insert((s1, s2));
                        }
                    }
                }
            }
            out
        })
        .collect();

    let mut all = HashSet::new();
    for s in parts.into_iter().chain(parts5.into_iter()) {
        all.extend(s);
    }
    all
}

/// Conjecture test over ALL split isogeny classes E1 × E2 over F_p (not just
/// the j = 0 sextic twists).  For a prime field every integer t with
/// |t| ≤ ⌊2√p⌋ is the trace of some elliptic curve (Deuring), so we range over
/// all unordered pairs t1 ≠ t2.
///
/// Conjecture (T22): the class contains a Jacobian  ⟺  N₁ = p+1−(t1+t2) ≥ 0
/// and |t1 − t2| ≥ 2.
///
/// Rationale for the second clause: a ppav isogenous to E1 × E2 that is not a
/// product is a gluing of E1 and E2 along an n-torsion subgroup for some
/// n ≥ 2, which needs E1[n] ≅ E2[n] as Galois modules, forcing n | (t1 − t2).
/// |t1 − t2| = 1 leaves no admissible n.
fn conjecture_test(ctx: &Ctx) {
    let p = ctx.p as i64;
    let realisable = census_all(ctx);
    let bound = (2.0 * (p as f64).sqrt()).floor() as i64;

    let (mut ok, mut viol_a, mut viol_b, mut n_pairs) = (0, 0, 0, 0);
    let mut viol_examples: Vec<String> = vec![];
    for t1 in -bound..=bound {
        for t2 in (t1 + 1)..=bound {
            n_pairs += 1;
            let s1 = t1 + t2;
            let s2 = t1 * t2 + 2 * p;
            let n1_ok = p + 1 - s1 >= 0;
            let pred = n1_ok && (t2 - t1) >= 2;
            let got = realisable.contains(&(s1, s2));
            if pred == got {
                ok += 1;
            } else if pred && !got {
                viol_a += 1;
                if viol_examples.len() < 6 {
                    viol_examples.push(format!("predicted-yes/absent  (t1,t2)=({},{})", t1, t2));
                }
            } else {
                viol_b += 1;
                if viol_examples.len() < 6 {
                    viol_examples.push(format!("predicted-no/present  (t1,t2)=({},{})", t1, t2));
                }
            }
        }
    }
    println!("\n  --- conjecture T22 over ALL trace pairs |t| ≤ {} ---", bound);
    println!("  distinct realisable genus-2 char polys over F_{}: {}", p, realisable.len());
    println!("  pairs tested: {}   agree: {}   pred-yes-but-absent: {}   pred-no-but-present: {}",
             n_pairs, ok, viol_a, viol_b);
    for e in &viol_examples {
        println!("    counterexample: {}", e);
    }
}

fn poly_str(lead: u64, coeffs: &[u64], p: u64) -> String {
    let mut terms: Vec<String> = vec![];
    for (i, &c) in coeffs.iter().enumerate().rev() {
        let c = c * lead % p;
        if c == 0 {
            continue;
        }
        terms.push(match i {
            0 => format!("{}", c),
            1 => if c == 1 { "x".into() } else { format!("{}*x", c) },
            _ => if c == 1 { format!("x^{}", i) } else { format!("{}*x^{}", c, i) },
        });
    }
    if terms.is_empty() {
        "0".into()
    } else {
        terms.join(" + ")
    }
}

fn run(p: u64) {
    assert!(p % 6 == 1, "need p ≡ 1 (mod 6) for six sextic twists");
    let ctx = Ctx::new(p);
    let tw = sextic_twists(&ctx);

    println!("\n================================================================");
    println!("p = {}   (p mod 6 = {}), non-residue ν = {}", p, p % 6, ctx.nu);
    println!("================================================================");
    println!("  k |    b_k |  t_k |  #E_k | #roots(x³+b_k)");
    for (k, &(b, t, roots, card)) in tw.iter().enumerate() {
        println!("  {} | {:6} | {:4} | {:5} | {}", k, b, t, card, roots);
    }
    let tsum: i64 = tw.iter().map(|x| x.1).sum();
    println!("  Σ t_k = {}  (must be 0 by twist symmetry)", tsum);

    // Howe conditions and target char polys.
    let mut targets: HashMap<(i64, i64), Vec<(usize, usize)>> = HashMap::new();
    let mut pairs = vec![];
    for i in 0..6 {
        for j in (i + 1)..6 {
            let (_, ti, ri, ci) = tw[i];
            let (_, tj, rj, cj) = tw[j];
            let h1 = ci != cj;
            let h2 = ri == rj;
            let h3 = gcd(ci, cj) == 1;
            let s1 = ti + tj;
            let s2 = ti * tj + 2 * p as i64;
            targets.entry((s1, s2)).or_default().push((i, j));
            pairs.push((i, j, h1, h2, h3, s1, s2));
        }
    }

    let t0 = std::time::Instant::now();
    let found = census(&ctx, &targets);
    let elapsed = t0.elapsed();

    println!("\n  exhaustive census over 2·p⁵ + 2·p⁴ = {} models, {:.1}s",
             2 * p.pow(5) + 2 * p.pow(4), elapsed.as_secs_f64());
    println!("\n  pair | H1 H2 H3 | Howe |   s1 |    s2 | genus-2 curve realising P_i·P_j");
    println!("  -----+----------+------+------+-------+---------------------------------");
    let mut agree = 0;
    let mut howe_yes_real_no = 0;
    let mut howe_no_real_yes = 0;
    for &(i, j, h1, h2, h3, s1, s2) in &pairs {
        let howe = h1 && h2 && h3;
        let hit = found.get(&(s1, s2));
        let desc = match hit {
            Some(h) => format!("y² = {}   [deg {}]", poly_str(h.lead, &h.coeffs, p), h.deg),
            None => "— none exists —".to_string(),
        };
        println!(
            "  ({},{}) |  {}  {}  {} | {:>4} | {:4} | {:5} | {}",
            i, j,
            if h1 { "Y" } else { "n" },
            if h2 { "Y" } else { "n" },
            if h3 { "Y" } else { "n" },
            if howe { "YES" } else { "no" },
            s1, s2, desc
        );
        match (howe, hit.is_some()) {
            (true, true) | (false, false) => agree += 1,
            (true, false) => howe_yes_real_no += 1,
            (false, true) => howe_no_real_yes += 1,
        }
    }
    let n_howe = pairs.iter().filter(|x| x.2 && x.3 && x.4).count();
    let n_real = pairs
        .iter()
        .filter(|&&(_, _, _, _, _, s1, s2)| found.contains_key(&(s1, s2)))
        .count();
    println!("\n  Howe-glueable pairs : {}/15", n_howe);
    println!("  realised as Jac     : {}/15", n_real);
    println!("  agreement           : {}/15  (Howe⇒realised misses: {}, realised-but-not-Howe: {})",
             agree, howe_yes_real_no, howe_no_real_yes);
}

fn main() {
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    // `--all-pairs` additionally runs the conjecture test over every split
    // isogeny class (expensive: full F_p² pass per model, so p ≤ 19).
    let all_pairs = args.iter().any(|a| a == "--all-pairs");
    args.retain(|a| !a.starts_with("--"));
    let primes: Vec<u64> = if args.is_empty() {
        vec![7, 13, 19, 31]
    } else {
        args.iter().map(|s| s.parse().expect("prime")).collect()
    };
    println!("Thread 22 — exhaustive genus-2 Jacobian census for j=0 sextic twists");
    for p in primes {
        run(p);
        if all_pairs {
            let ctx = Ctx::new(p);
            let t0 = std::time::Instant::now();
            conjecture_test(&ctx);
            println!("  ({:.1}s)", t0.elapsed().as_secs_f64());
        }
    }
}

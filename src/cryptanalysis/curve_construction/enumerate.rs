//! Exhaustive enumeration of curves of small genus over `F_2`.
//!
//! Every smooth curve of genus 2, 3 or 4 over `F_2` appears below with at
//! least one model:
//!
//! * hyperelliptic: `y² + h(x)·y = f(x)`, `deg h ≤ g+1`, `deg f ≤ 2g+2`, smooth in
//!   weighted projective space (exact criterion, both charts);
//! * genus 3 non-hyperelliptic: smooth plane quartics, all `2^15` forms;
//! * genus 4 non-hyperelliptic: canonical curves `Q = K = 0` in `P³`, with `Q`
//!   one of the three quadric types over `F_2` (hyperbolic, elliptic, cone) and
//!   `K` running over cubics modulo `Q·(linear)`.
//!
//! Genus 5 is enumerated for hyperelliptic models only.  For each model the
//! characteristic polynomial of Frobenius is read off its point counts, and
//! every target polynomial is tested for divisibility.  A miss needs no
//! smoothness check — completeness of the enumeration is what makes a miss a
//! proof — but every hit is checked for smoothness before it is reported.

use super::gf2k::SmallField;
use super::lpoly::{char_poly_from_counts, divide_exact};
use rayon::prelude::*;
use serde::Serialize;
use std::collections::BTreeMap;

/// A target abelian variety over `F_2`, by its monic characteristic polynomial.
#[derive(Clone, Debug)]
pub struct Target {
    pub name: String,
    pub poly: Vec<i128>,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct TargetResult {
    pub target: String,
    /// Smooth models whose characteristic polynomial the target divides.
    pub hit_models: u64,
    /// Cofactor `P_C / P_target` → number of hit models with that cofactor.
    pub cofactors: BTreeMap<String, u64>,
    /// Up to eight hit models, family-specific encoding.
    pub examples: Vec<String>,
}

#[derive(Clone, Debug, Serialize)]
pub struct FamilyResult {
    pub family: String,
    pub genus: usize,
    pub models_enumerated: u64,
    /// Models with a valid (Newton-integral) characteristic polynomial; for the
    /// hyperelliptic family, the models that pass the exact smoothness test.
    pub models_admitted: u64,
    pub complete_for_genus: bool,
    pub targets: Vec<TargetResult>,
}

fn poly_key(p: &[i128]) -> String {
    p.iter()
        .map(|c| c.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn record(
    results: &mut [TargetResult],
    targets: &[Target],
    p: &[i128],
    model: &str,
    smooth: impl Fn() -> bool,
) {
    for (t, res) in targets.iter().zip(results.iter_mut()) {
        if t.poly.len() > p.len() {
            continue;
        }
        if let Some(cof) = divide_exact(p, &t.poly) {
            if !smooth() {
                continue;
            }
            res.hit_models += 1;
            *res.cofactors.entry(poly_key(&cof)).or_insert(0) += 1;
            if res.examples.len() < 8 {
                res.examples.push(model.to_string());
            }
        }
    }
}

fn merge(into: &mut [TargetResult], from: Vec<TargetResult>) {
    for (a, b) in into.iter_mut().zip(from) {
        a.hit_models += b.hit_models;
        for (k, v) in b.cofactors {
            *a.cofactors.entry(k).or_insert(0) += v;
        }
        for e in b.examples {
            if a.examples.len() < 8 {
                a.examples.push(e);
            }
        }
    }
}

fn empty_results(targets: &[Target]) -> Vec<TargetResult> {
    targets
        .iter()
        .map(|t| TargetResult {
            target: t.name.clone(),
            ..Default::default()
        })
        .collect()
}

// ---------------------------------------------------------------------------
// F_2[x] as u64 bit vectors
// ---------------------------------------------------------------------------

fn pdeg(a: u64) -> i32 {
    63 - a.leading_zeros() as i32
}

fn pmul(mut a: u64, mut b: u64) -> u64 {
    let mut r = 0u64;
    while b != 0 {
        if b & 1 == 1 {
            r ^= a;
        }
        b >>= 1;
        a <<= 1;
    }
    r
}

fn pmod(mut a: u64, m: u64) -> u64 {
    let dm = pdeg(m);
    while a != 0 && pdeg(a) >= dm {
        a ^= m << (pdeg(a) - dm);
    }
    a
}

fn pgcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let r = pmod(a, b);
        a = b;
        b = r;
    }
    a
}

/// Formal derivative in characteristic 2: odd-degree terms shift down.
fn pderiv(a: u64) -> u64 {
    (a >> 1) & 0x5555_5555_5555_5555
}

/// Exact smoothness of `y² + h y = f` (`deg h ≤ g+1`, `deg f ≤ 2g+2`) on both
/// charts of the weighted projective model.  Singular affine points sit over
/// common roots of `h` and `h'²f + f'²`; the point(s) at infinity are singular
/// iff `h_{g+1} = 0` and `h_g f_{2g+2} + f_{2g+1} = 0`.
pub fn hyperelliptic_smooth(h: u64, f: u64, g: usize) -> bool {
    if h == 0 {
        return false;
    }
    let hp = pderiv(h);
    let fp = pderiv(f);
    if pgcd(h, pmul(pmul(hp, hp), f) ^ pmul(fp, fp)) != 1 {
        return false;
    }
    let htop = (h >> (g + 1)) & 1;
    if htop == 0 {
        let hg = (h >> g) & 1;
        let ftop = (f >> (2 * g + 2)) & 1;
        let fsub = (f >> (2 * g + 1)) & 1;
        if (hg & ftop) ^ fsub == 0 {
            return false;
        }
    }
    true
}

struct HyperField {
    q: usize,
    nh: usize,
    nf: usize,
    hval: Vec<u16>,
    fval: Vec<u16>,
    contrib: Vec<u8>,
    k_even: bool,
}

fn eval_table(field: &SmallField, nbits: usize) -> Vec<u16> {
    let q = field.q;
    let mut out = vec![0u16; q << nbits];
    for x in 0..q {
        let mut pw = vec![1u16; nbits];
        for i in 1..nbits {
            pw[i] = field.mul(pw[i - 1], x as u16);
        }
        let row = &mut out[x << nbits..(x + 1) << nbits];
        for poly in 1..(1usize << nbits) {
            row[poly] = row[poly & (poly - 1)] ^ pw[poly.trailing_zeros() as usize];
        }
    }
    out
}

impl HyperField {
    fn new(k: u32, g: usize) -> Self {
        let field = SmallField::new(k);
        let q = field.q;
        let (bh, bf) = (g + 2, 2 * g + 3);
        let mut contrib = vec![0u8; q * q];
        for hv in 0..q {
            for fv in 0..q {
                contrib[hv * q + fv] = if hv == 0 {
                    1
                } else {
                    let ih = field.inv(hv as u16);
                    let c = field.mul(fv as u16, field.mul(ih, ih));
                    2 * (1 - field.trace(c))
                };
            }
        }
        HyperField {
            q,
            nh: 1 << bh,
            nf: 1 << bf,
            hval: eval_table(&field, bh),
            fval: eval_table(&field, bf),
            contrib,
            k_even: k.is_multiple_of(2),
        }
    }

    #[inline]
    fn count(&self, h: usize, f: usize, g: usize) -> i64 {
        let mut c = 0i64;
        for x in 0..self.q {
            let hv = self.hval[x * self.nh + h] as usize;
            let fv = self.fval[x * self.nf + f] as usize;
            c += self.contrib[hv * self.q + fv] as i64;
        }
        if (h >> (g + 1)) & 1 == 1 {
            if (f >> (2 * g + 2)) & 1 == 0 || self.k_even {
                c += 2;
            }
        } else {
            c += 1;
        }
        c
    }
}

/// Every hyperelliptic model of genus `g` over `F_2`.
pub fn hyperelliptic(g: usize, targets: &[Target]) -> FamilyResult {
    let fields: Vec<HyperField> = (1..=g as u32).map(|k| HyperField::new(k, g)).collect();
    let nh = 1usize << (g + 2);
    let nf = 1usize << (2 * g + 3);
    let partial: Vec<(u64, Vec<TargetResult>)> = (1..nh)
        .into_par_iter()
        .map(|h| {
            let mut res = empty_results(targets);
            let mut admitted = 0u64;
            for f in 0..nf {
                if !hyperelliptic_smooth(h as u64, f as u64, g) {
                    continue;
                }
                admitted += 1;
                let counts: Vec<i64> = fields.iter().map(|fl| fl.count(h, f, g)).collect();
                let Some(p) = char_poly_from_counts(&counts, 2) else {
                    continue;
                };
                record(&mut res, targets, &p, &format!("h={h:#x},f={f:#x}"), || {
                    true
                });
            }
            (admitted, res)
        })
        .collect();
    let mut results = empty_results(targets);
    let mut admitted = 0;
    for (a, r) in partial {
        admitted += a;
        merge(&mut results, r);
    }
    FamilyResult {
        family: "hyperelliptic".into(),
        genus: g,
        models_enumerated: ((nh - 1) * nf) as u64,
        models_admitted: admitted,
        complete_for_genus: g == 2,
        targets: results,
    }
}

// ---------------------------------------------------------------------------
// projective points and forms
// ---------------------------------------------------------------------------

/// Normalised points of `P^{dim}(F_{2^k})` (first non-zero coordinate is 1).
pub fn projective_points(field: &SmallField, dim: usize) -> Vec<Vec<u16>> {
    let q = field.q;
    let mut pts = Vec::new();
    for lead in 0..=dim {
        let free = dim - lead;
        let total = q.pow(free as u32);
        for idx in 0..total {
            let mut p = vec![0u16; dim + 1];
            p[lead] = 1;
            let mut r = idx;
            for slot in p.iter_mut().skip(lead + 1) {
                *slot = (r % q) as u16;
                r /= q;
            }
            pts.push(p);
        }
    }
    pts
}

/// Monomials of degree `deg` in `nvars` variables, as exponent vectors.
pub fn monomials(nvars: usize, deg: usize) -> Vec<Vec<usize>> {
    fn rec(nvars: usize, left: usize, cur: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if cur.len() == nvars - 1 {
            cur.push(left);
            out.push(cur.clone());
            cur.pop();
            return;
        }
        for e in (0..=left).rev() {
            cur.push(e);
            rec(nvars, left - e, cur, out);
            cur.pop();
        }
    }
    let mut out = Vec::new();
    rec(nvars, deg, &mut Vec::new(), &mut out);
    out
}

fn mono_value(field: &SmallField, pt: &[u16], e: &[usize]) -> u16 {
    let mut v = 1u16;
    for (x, &k) in pt.iter().zip(e) {
        for _ in 0..k {
            v = field.mul(v, *x);
        }
    }
    v
}

/// Value of a form (bitmask over `monos`) at a point.
pub fn form_value(field: &SmallField, monos: &[Vec<usize>], form: u64, pt: &[u16]) -> u16 {
    let mut v = 0u16;
    for (i, e) in monos.iter().enumerate() {
        if (form >> i) & 1 == 1 {
            v ^= mono_value(field, pt, e);
        }
    }
    v
}

/// Partial derivative `∂/∂x_var` of a form, as a list of (exponent vector) terms.
fn partial_terms(monos: &[Vec<usize>], form: u64, var: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    for (i, e) in monos.iter().enumerate() {
        if (form >> i) & 1 == 1 && e[var] % 2 == 1 {
            let mut d = e.clone();
            d[var] -= 1;
            out.push(d);
        }
    }
    out
}

fn terms_value(field: &SmallField, terms: &[Vec<usize>], pt: &[u16]) -> u16 {
    terms
        .iter()
        .fold(0u16, |acc, e| acc ^ mono_value(field, pt, e))
}

/// Nibble tables: 16-entry XOR tables per 4-bit chunk of a form's bitmask.
fn nibble_tables(vals: &[u16]) -> Vec<[u16; 16]> {
    vals.chunks(4)
        .map(|chunk| {
            let mut t = [0u16; 16];
            for (m, slot) in t.iter_mut().enumerate() {
                for (b, v) in chunk.iter().enumerate() {
                    if (m >> b) & 1 == 1 {
                        *slot ^= v;
                    }
                }
            }
            t
        })
        .collect()
}

#[inline]
fn nibble_eval(tabs: &[[u16; 16]], form: usize) -> u16 {
    let mut v = 0u16;
    for (i, t) in tabs.iter().enumerate() {
        v ^= t[(form >> (4 * i)) & 15];
    }
    v
}

// ---------------------------------------------------------------------------
// plane quartics
// ---------------------------------------------------------------------------

pub fn quartic_monomials() -> Vec<Vec<usize>> {
    monomials(3, 4)
}

/// `#C(F_{2^k})` for a plane curve given by a form.
pub fn plane_count(field: &SmallField, monos: &[Vec<usize>], form: u64) -> i64 {
    projective_points(field, 2)
        .iter()
        .filter(|p| form_value(field, monos, form, p) == 0)
        .count() as i64
}

/// No singular point over `F_{2^k}`, `k ≤ kmax` (for a plane quartic every
/// singular point is defined over an extension of degree ≤ 4).
pub fn plane_smooth(monos: &[Vec<usize>], form: u64, kmax: u32) -> bool {
    let partials: Vec<Vec<Vec<usize>>> = (0..3).map(|v| partial_terms(monos, form, v)).collect();
    for k in 1..=kmax {
        let field = SmallField::new(k);
        for pt in projective_points(&field, 2) {
            if form_value(&field, monos, form, &pt) == 0
                && partials.iter().all(|t| terms_value(&field, t, &pt) == 0)
            {
                return false;
            }
        }
    }
    true
}

/// Every plane quartic over `F_2`.  Returns the family result and, per target,
/// every smooth hit form (for orbit analysis).
pub fn plane_quartics(targets: &[Target]) -> (FamilyResult, Vec<Vec<u64>>) {
    let monos = quartic_monomials();
    let nforms = 1usize << monos.len();
    let mut counts = vec![[0i64; 3]; nforms];
    for k in 1..=3u32 {
        let field = SmallField::new(k);
        for pt in projective_points(&field, 2) {
            let vals: Vec<u16> = monos.iter().map(|e| mono_value(&field, &pt, e)).collect();
            let tabs = nibble_tables(&vals);
            for (form, c) in counts.iter_mut().enumerate().skip(1) {
                if nibble_eval(&tabs, form) == 0 {
                    c[k as usize - 1] += 1;
                }
            }
        }
    }
    let mut results = empty_results(targets);
    let mut hit_forms: Vec<Vec<u64>> = vec![Vec::new(); targets.len()];
    let mut admitted = 0u64;
    for form in 1..nforms {
        let Some(p) = char_poly_from_counts(&counts[form], 2) else {
            continue;
        };
        admitted += 1;
        for (ti, t) in targets.iter().enumerate() {
            if t.poly.len() <= p.len()
                && divide_exact(&p, &t.poly).is_some()
                && plane_smooth(&monos, form as u64, 6)
            {
                hit_forms[ti].push(form as u64);
            }
        }
        record(
            &mut results,
            targets,
            &p,
            &format!("form={form:#06x}"),
            || plane_smooth(&monos, form as u64, 6),
        );
    }
    (
        FamilyResult {
            family: "plane quartic".into(),
            genus: 3,
            models_enumerated: (nforms - 1) as u64,
            models_admitted: admitted,
            complete_for_genus: true,
            targets: results,
        },
        hit_forms,
    )
}

/// Apply `x ↦ A·x` (A a 3×3 matrix over F_2, rows as bitmasks) to a quartic form.
pub fn substitute_quartic(monos: &[Vec<usize>], form: u64, a: [[u8; 3]; 3]) -> u64 {
    // dense polynomials in 3 variables, degree ≤ 4: index e0*25 + e1*5 + e2
    let idx = |e: &[usize]| e[0] * 25 + e[1] * 5 + e[2];
    let mul = |p: &[u8; 125], q: &[u8; 125]| {
        let mut r = [0u8; 125];
        for i in 0..125 {
            if p[i] == 0 {
                continue;
            }
            let (a0, a1, a2) = (i / 25, (i / 5) % 5, i % 5);
            for j in 0..125 {
                if q[j] == 0 {
                    continue;
                }
                let (b0, b1, b2) = (j / 25, (j / 5) % 5, j % 5);
                if a0 + b0 < 5 && a1 + b1 < 5 && a2 + b2 < 5 {
                    r[(a0 + b0) * 25 + (a1 + b1) * 5 + a2 + b2] ^= 1;
                }
            }
        }
        r
    };
    let lin: Vec<[u8; 125]> = (0..3)
        .map(|row| {
            let mut p = [0u8; 125];
            for (col, unit) in [[1usize, 0, 0], [0, 1, 0], [0, 0, 1]].iter().enumerate() {
                if a[row][col] == 1 {
                    p[idx(unit)] ^= 1;
                }
            }
            p
        })
        .collect();
    let mut total = [0u8; 125];
    for (i, e) in monos.iter().enumerate() {
        if (form >> i) & 1 == 0 {
            continue;
        }
        let mut term = [0u8; 125];
        term[0] = 1;
        for (var, &k) in e.iter().enumerate() {
            for _ in 0..k {
                term = mul(&term, &lin[var]);
            }
        }
        for j in 0..125 {
            total[j] ^= term[j];
        }
    }
    let mut out = 0u64;
    for (i, e) in monos.iter().enumerate() {
        if total[idx(e)] == 1 {
            out |= 1 << i;
        }
    }
    out
}

/// `GL_3(F_2)`, all 168 matrices.
pub fn gl3_f2() -> Vec<[[u8; 3]; 3]> {
    let mut out = Vec::new();
    for bits in 0u32..512 {
        let a: [[u8; 3]; 3] =
            std::array::from_fn(|i| std::array::from_fn(|j| ((bits >> (3 * i + j)) & 1) as u8));
        let det = (a[0][0] & (a[1][1] & a[2][2] ^ a[1][2] & a[2][1]))
            ^ (a[0][1] & (a[1][0] & a[2][2] ^ a[1][2] & a[2][0]))
            ^ (a[0][2] & (a[1][0] & a[2][1] ^ a[1][1] & a[2][0]));
        if det == 1 {
            out.push(a);
        }
    }
    out
}

// ---------------------------------------------------------------------------
// canonical genus-4 curves: Q ∩ K in P³
// ---------------------------------------------------------------------------

/// The three quadric types over F_2 carrying a canonical genus-4 curve.
pub fn quadrics() -> Vec<(&'static str, Vec<(Vec<usize>, u8)>)> {
    let m = |e: [usize; 4]| e.to_vec();
    vec![
        (
            "hyperbolic x0x1+x2x3",
            vec![(m([1, 1, 0, 0]), 1), (m([0, 0, 1, 1]), 1)],
        ),
        (
            "elliptic x0x1+x2^2+x2x3+x3^2",
            vec![
                (m([1, 1, 0, 0]), 1),
                (m([0, 0, 2, 0]), 1),
                (m([0, 0, 1, 1]), 1),
                (m([0, 0, 0, 2]), 1),
            ],
        ),
        (
            "cone x0x1+x2^2",
            vec![(m([1, 1, 0, 0]), 1), (m([0, 0, 2, 0]), 1)],
        ),
    ]
}

/// Cubic monomials in x0..x3 other than the pivots x0²x1, x0x1², x0x1x2,
/// x0x1x3 — a complement to Q·(linear) for all three quadrics above.
pub fn free_cubic_monomials() -> Vec<Vec<usize>> {
    let pivots = [[2, 1, 0, 0], [1, 2, 0, 0], [1, 1, 1, 0], [1, 1, 0, 1]];
    monomials(4, 3)
        .into_iter()
        .filter(|e| !pivots.iter().any(|p| p.as_slice() == e.as_slice()))
        .collect()
}

fn quadric_value(field: &SmallField, q: &[(Vec<usize>, u8)], pt: &[u16]) -> u16 {
    q.iter()
        .fold(0u16, |acc, (e, _)| acc ^ mono_value(field, pt, e))
}

fn quadric_points(field: &SmallField, q: &[(Vec<usize>, u8)]) -> Vec<Vec<u16>> {
    projective_points(field, 3)
        .into_iter()
        .filter(|p| quadric_value(field, q, p) == 0)
        .collect()
}

/// No point of `Q ∩ K` over `F_{2^k}`, `k ≤ kmax`, where the Jacobian
/// `[∇Q; ∇K]` has rank < 2.
pub fn ci_smooth(q: &[(Vec<usize>, u8)], monos: &[Vec<usize>], cubic: u64, kmax: u32) -> bool {
    let qforms: Vec<Vec<usize>> = q.iter().map(|(e, _)| e.clone()).collect();
    let qmask = (1u64 << qforms.len()) - 1;
    let gq: Vec<Vec<Vec<usize>>> = (0..4).map(|v| partial_terms(&qforms, qmask, v)).collect();
    let gk: Vec<Vec<Vec<usize>>> = (0..4).map(|v| partial_terms(monos, cubic, v)).collect();
    for k in 1..=kmax {
        let field = SmallField::new(k);
        for pt in quadric_points(&field, q) {
            if form_value(&field, monos, cubic, &pt) != 0 {
                continue;
            }
            let a: Vec<u16> = gq.iter().map(|t| terms_value(&field, t, &pt)).collect();
            let b: Vec<u16> = gk.iter().map(|t| terms_value(&field, t, &pt)).collect();
            let mut rank2 = false;
            for i in 0..4 {
                for j in i + 1..4 {
                    if field.mul(a[i], b[j]) ^ field.mul(a[j], b[i]) != 0 {
                        rank2 = true;
                    }
                }
            }
            if !rank2 {
                return false;
            }
        }
    }
    true
}

/// Every canonical genus-4 model `Q ∩ K` over `F_2`.
pub fn canonical_genus4(targets: &[Target]) -> FamilyResult {
    let monos = free_cubic_monomials();
    assert_eq!(monos.len(), 16);
    let ncub = 1usize << 16;
    let mut results = empty_results(targets);
    let mut admitted = 0u64;
    let mut enumerated = 0u64;
    for (qname, q) in quadrics() {
        let mut counts = vec![[0i64; 4]; ncub];
        for k in 1..=4u32 {
            let field = SmallField::new(k);
            for pt in quadric_points(&field, &q) {
                let vals: Vec<u16> = monos.iter().map(|e| mono_value(&field, &pt, e)).collect();
                let tabs = nibble_tables(&vals);
                for (cubic, c) in counts.iter_mut().enumerate().skip(1) {
                    if nibble_eval(&tabs, cubic) == 0 {
                        c[k as usize - 1] += 1;
                    }
                }
            }
        }
        for cubic in 1..ncub {
            enumerated += 1;
            let Some(p) = char_poly_from_counts(&counts[cubic], 2) else {
                continue;
            };
            admitted += 1;
            record(
                &mut results,
                targets,
                &p,
                &format!("{qname}; K={cubic:#06x}"),
                || ci_smooth(&q, &monos, cubic as u64, 6),
            );
        }
    }
    FamilyResult {
        family: "canonical (2,3) complete intersection".into(),
        genus: 4,
        models_enumerated: enumerated,
        models_admitted: admitted,
        complete_for_genus: true,
        targets: results,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn smoothness_criterion_on_known_curves() {
        // y² + y = x⁵ is a smooth genus-2 curve over F_2.
        assert!(hyperelliptic_smooth(0b1, 1 << 5, 2));
        assert!(
            !hyperelliptic_smooth(0, 1 << 5, 2),
            "h = 0 is purely inseparable"
        );
    }

    #[test]
    fn klein_quartic_is_smooth_with_known_counts() {
        let monos = quartic_monomials();
        let klein = [[3usize, 1, 0], [0, 3, 1], [1, 0, 3]]
            .iter()
            .map(|e| 1u64 << monos.iter().position(|m| m.as_slice() == e).unwrap())
            .fold(0, |a, b| a | b);
        let counts: Vec<i64> = (1..=3)
            .map(|k| plane_count(&SmallField::new(k), &monos, klein))
            .collect();
        assert_eq!(counts, vec![3, 5, 24]);
        assert!(plane_smooth(&monos, klein, 4));
        assert_eq!(gl3_f2().len(), 168);
    }

    #[test]
    fn genus_two_counts_match_direct_enumeration() {
        // y² + y = x⁵ + x³: count points over F_2..F_4 directly and via tables.
        let (h, f, g) = (1usize, (1usize << 5) | (1 << 3), 2usize);
        for k in 1..=2u32 {
            let field = SmallField::new(k);
            let hf = HyperField::new(k, g);
            let mut direct = 0i64;
            for x in 0..field.q as u16 {
                let fx = field.mul(field.mul(field.mul(x, x), x), field.mul(x, x))
                    ^ field.mul(field.mul(x, x), x);
                for y in 0..field.q as u16 {
                    if field.mul(y, y) ^ y == fx {
                        direct += 1;
                    }
                }
            }
            direct += 1; // deg h = 0 < g + 1: one point at infinity
            assert_eq!(hf.count(h, f, g), direct);
        }
    }
}

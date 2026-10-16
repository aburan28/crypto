//! CM-guided volcano interpolation (Bröker–Lauter–Sutherland, Algorithm 2.1).
//!
//! This bounded reference implementation uses odd fundamental discriminants with
//! split ell, v=1 or 2, and certified traces at primes below 2^40. Class actions use
//! Frobenius eigenkernels, not a precomputed modular-polynomial table. The floor
//! class group is enumerated with its projection to the surface class group.
//! These deliberately simple choices do not implement the paper's asymptotic
//! parameter selection or its optimized polycyclic presentations.
use super::crt::ExplicitCrt;
use crate::bigint::Big;
use crate::cm::{hilbert_class_poly, reduced_forms};
use crate::curve::{self, Curve, Pt};
use crate::field::{is_prime, Field, Rng, Zp};
use crate::find::{divpoly::division_poly, sea::ring_mul};
use crate::int::Int;
use crate::kernel::kohel::kohel;
use crate::path::relation::Form;
use crate::poly;
use std::collections::{HashMap, HashSet};

/// Arithmetic and trial factorization are bounded by this prime limit.
pub const MAX_CM_PRIME: u64 = (1 << 40) - 1;
const MAX_DISCRIMINANT: i64 = 16_000_000;
const MAX_CLASSES: usize = 4096;
const MAX_AUX_PRIME: u64 = 43;

#[derive(Clone, Debug)]
struct Generator {
    q: u64,
    b: i128,
    top: Form,
    floor: Form,
}
#[derive(Clone, Debug)]
struct ClassWalk {
    classes: Vec<Form>,
    index: HashMap<Form, usize>,
    steps: Vec<(usize, usize)>,
    projection: Vec<Form>,
}

fn class_walk(d: i128, top_d: i128, gs: &[Generator], floor: bool) -> Result<ClassWalk, String> {
    let id = Form::identity(d);
    let mut w = ClassWalk {
        classes: vec![id],
        index: HashMap::from([(id, 0)]),
        steps: vec![],
        projection: vec![Form::identity(top_d)],
    };
    let mut i = 0;
    while i < w.classes.len() {
        for (gidx, g) in gs.iter().enumerate() {
            let c = w.classes[i].compose(if floor { &g.floor } else { &g.top });
            let projection = w.projection[i].compose(&g.top);
            if let Some(&k) = w.index.get(&c) {
                if w.projection[k] != projection {
                    return Err("inconsistent class-group projection".into());
                }
            } else {
                if w.classes.len() == MAX_CLASSES {
                    return Err("class group exceeds resource limit".into());
                }
                w.index.insert(c, w.classes.len());
                w.classes.push(c);
                w.projection.push(projection);
                w.steps.push((i, gidx));
            }
        }
        i += 1;
    }
    Ok(w)
}

fn prime_form(d: i128, q: u64) -> Option<Form> {
    let q = q as i128;
    (1..2 * q)
        .step_by(2)
        .find(|&b| (b * b - d).rem_euclid(4 * q) == 0)
        .map(|b| Form {
            a: q,
            b,
            c: (b * b - d) / (4 * q),
        })
}
fn fundamental(d: i64) -> bool {
    if d >= -4 || d.rem_euclid(4) != 1 {
        return false;
    }
    let n = -d;
    let mut q = 3;
    while q * q <= n {
        if n % (q * q) == 0 {
            return false;
        }
        q += 2;
    }
    true
}

/// A reusable CM setup. Coefficient matrices use c[i][j] for X^i Y^j.
#[derive(Clone, Debug)]
pub struct CmPlan {
    ell: usize,
    d: i64,
    v: u64,
    h: Vec<Int>,
    generators: Vec<Generator>,
    top: ClassWalk,
    floor: ClassWalk,
    horizontal: Form,
}

impl CmPlan {
    /// Select the first suitable fundamental discriminant with |D| <= 20000.
    /// Supports odd prime degrees 3..=31; failure is an explicit resource limit.
    pub fn new(ell: usize) -> Result<Self, String> {
        validate_degree(ell)?;
        for a in (7..=20_000i64).step_by(4) {
            let d = -a;
            if fundamental(d)
                && a * (ell * ell) as i64 <= MAX_DISCRIMINANT
                && d.rem_euclid(ell as i64) != 0
                && prime_form(d as i128, ell as u64).is_some()
                && reduced_forms(d).len() >= ell + 2
            {
                return Self::with_discriminant(ell, d);
            }
        }
        Err("no discriminant found within the CM search limits".into())
    }

    /// Use a specified odd fundamental D, with ell split and h(D) >= ell+2.
    /// The existing numerical H_D routine supplies candidates; each finite-field
    /// run certifies its seed's trace and maximal endomorphism ring algebraically.
    pub fn with_discriminant(ell: usize, d: i64) -> Result<Self, String> {
        validate_degree(ell)?;
        if d >= -4 || d < -MAX_DISCRIMINANT / (ell * ell) as i64 || !fundamental(d) {
            return Err("D must be odd fundamental, D < -4, with |ell^2 D| <= 16000000".into());
        }
        if d.rem_euclid(ell as i64) == 0 {
            return Err("ell must split in D".into());
        }
        let horizontal = prime_form(d as i128, ell as u64)
            .ok_or("ell must split in D")?
            .reduce();
        let r = d * (ell * ell) as i64;
        let h_top = reduced_forms(d).len();
        let h_floor = reduced_forms(r).len();
        if h_top < ell + 2 || h_floor != h_top * (ell - 1) || h_floor > MAX_CLASSES {
            return Err("class numbers do not satisfy the CM interpolation limits".into());
        }
        let mut generators = vec![];
        let mut top = class_walk(d as i128, d as i128, &generators, false)?;
        let mut floor = class_walk(r as i128, d as i128, &generators, true)?;
        for q in (3..=MAX_AUX_PRIME).step_by(2) {
            if q == ell as u64 || !is_prime(q) || d.rem_euclid(q as i64) == 0 {
                continue;
            }
            if let Some(g) = prime_form(d as i128, q) {
                let b = (ell as i128 * g.b).rem_euclid(2 * q as i128);
                let bottom = Form {
                    a: q as i128,
                    b,
                    c: (b * b - r as i128) / (4 * q as i128),
                }
                .reduce();
                if top.index.contains_key(&g.reduce()) && floor.index.contains_key(&bottom) {
                    continue;
                }
                generators.push(Generator {
                    q,
                    b: g.b,
                    top: g.reduce(),
                    floor: bottom,
                });
                top = class_walk(d as i128, d as i128, &generators, false)?;
                floor = class_walk(r as i128, d as i128, &generators, true)?;
                if top.classes.len() == h_top && floor.classes.len() == h_floor {
                    break;
                }
            }
        }
        if top.classes.len() != h_top || floor.classes.len() != h_floor {
            return Err("auxiliary primes <= 43 do not generate both class groups".into());
        }
        Ok(Self {
            ell,
            d,
            v: if d.rem_euclid(8) == 1 { 2 } else { 1 },
            h: hilbert_class_poly(d)?.coeffs,
            generators,
            top,
            floor,
            horizontal,
        })
    }

    pub fn discriminant(&self) -> i64 {
        self.d
    }
    pub fn degree(&self) -> usize {
        self.ell
    }
    pub fn class_numbers(&self) -> (usize, usize) {
        (self.top.classes.len(), self.floor.classes.len())
    }

    /// Validates 4p=t^2-v^2 ell^2 D, t=2 mod ell, and ell^3 not dividing #E.
    fn check_prime(&self, p: u64, t: u64) -> Result<(), String> {
        let ell = self.ell as u64;
        let norm = (t as u128) * (t as u128)
            + (self.v * ell) as u128 * (self.v * ell) as u128 * (-self.d) as u128;
        if p <= MAX_AUX_PRIME
            || p > MAX_CM_PRIME
            || !is_prime(p)
            || norm != 4 * p as u128
            || t % ell != 2
            || t.is_multiple_of(p)
            || ((-self.d) as u64).is_multiple_of(p)
            || t > p
        {
            return Err("prime/trace does not meet the bounded CM norm conditions".into());
        }
        let n = p + 1 - t;
        if !n.is_multiple_of(ell * ell) || n.is_multiple_of(ell * ell * ell) {
            return Err("CM descent requires exactly ell^2 dividing the point count".into());
        }
        Ok(())
    }

    /// Deterministic norm-equation search sufficient for a certified height bound.
    pub fn crt_primes(&self) -> Result<Vec<(u64, u64)>, String> {
        let bound = Big::from_u64(1).shl(coefficient_bound_bits(self.ell)? + 2);
        let mut product = Big::from_u64(1);
        let mut primes = vec![];
        let mut t = 2u64;
        loop {
            let norm = t * t + self.v * self.v * (self.ell * self.ell) as u64 * (-self.d) as u64;
            if norm / 4 > MAX_CM_PRIME {
                break;
            }
            if norm.is_multiple_of(4) {
                let p = norm / 4;
                if self.check_prime(p, t).is_ok() {
                    primes.push((p, t));
                    product = product.mul_small(p);
                    if product > bound {
                        return Ok(primes);
                    }
                }
            }
            t += self.ell as u64;
        }
        Err(format!(
            "insufficient CM primes below 2^40 for ell={}, D={} ({} primes found)",
            self.ell,
            self.d,
            primes.len()
        ))
    }

    /// Compute Phi_ell mod p by CM surface/floor walks and interpolation.
    pub fn phi_mod_p(&self, p: u64, t: u64, rng: &mut Rng) -> Result<Vec<Vec<u64>>, String> {
        self.check_prime(p, t)?;
        let f = Zp::new(p);
        let h: Vec<_> = self.h.iter().map(|c| c.mod_u64(p)).collect();
        let roots = poly::roots(&f, &h, rng);
        if roots.len() != self.top.classes.len() {
            return Err("H_D does not split with distinct roots".into());
        }
        let j = *roots
            .iter()
            .find(|&&j| j != 0 && j != 1728 % p)
            .ok_or("no nonspecial CM seed")?;
        let mut seed = curve::from_j(&f, j);
        let n = p + 1 - t;
        if !full_torsion(&f, &seed, self.ell) {
            let z = (2..p)
                .find(|&x| f.legendre(x) == -1)
                .ok_or("no twist parameter")?;
            seed = Curve::new(f.mul(seed.a, f.sq(z)), f.mul(seed.b, f.mul(f.sq(z), z)));
        }
        // Fundamental D and conductor(pi)=v*ell: full rational ell-torsion
        // forces conductor(End(E)) | v. If v=2, full E[2] forces it odd.
        if !full_torsion(&f, &seed, self.ell) || (self.v == 2 && !full_torsion(&f, &seed, 2)) {
            return Err("CM seed failed maximal-order torsion certification".into());
        }
        certify_trace(&f, &seed, n, self.ell as u64, rng)?;
        let top = self.walk_curves(&f, t, seed, &self.top)?;
        let js: Vec<_> = top.iter().map(|e| curve::jinv(&f, e)).collect();
        let surface: HashSet<_> = js.iter().copied().collect();
        if surface.len() != js.len() || js.iter().any(|&j| poly::eval(&f, &h, j) != 0) {
            return Err("surface class action failed H_D/distinctness checks".into());
        }
        // One descending isogeny suffices: the kernel of Cl(ell^2 D)->Cl(D)
        // groups the floor orbit into the ell-1 descendants of each surface j.
        let mut bottom = None;
        for _ in 0..128 {
            let pt = seed.random_point(&f, rng);
            let pt = curve::pmul(&f, &seed, &pt, (n / (self.ell * self.ell) as u64) as u128);
            if pt == Pt::Inf || curve::pmul(&f, &seed, &pt, self.ell as u128) != Pt::Inf {
                continue;
            }
            let (kernel, _) = curve::kernel_poly_from_point(&f, &seed, &pt, self.ell as u64);
            let cod = kohel(&f, &seed, &kernel, self.ell as u64).cod;
            if !surface.contains(&curve::jinv(&f, &cod)) {
                bottom = Some(cod);
                break;
            }
        }
        let floor = self.walk_curves(
            &f,
            t,
            bottom.ok_or("descent sampling budget exhausted")?,
            &self.floor,
        )?;
        let mut seen = HashSet::new();
        let mut children = vec![vec![]; top.len()];
        for (i, e) in floor.iter().enumerate() {
            let j = curve::jinv(&f, e);
            if surface.contains(&j) || !seen.insert(j) {
                return Err("floor orbit is not distinct from the surface".into());
            }
            let k = *self
                .top
                .index
                .get(&self.floor.projection[i])
                .ok_or("missing projected class")?;
            children[k].push(j);
        }
        let mut values = vec![];
        for (i, c) in self.top.classes.iter().enumerate() {
            if children[i].len() != self.ell - 1 {
                return Err("incorrect number of descending neighbors".into());
            }
            let mut neighbors = children[i].clone();
            for g in [self.horizontal, self.horizontal.inverse()] {
                let k = self.top.index[&c.compose(&g)];
                neighbors.push(js[k]); // multiplicities are intentional
            }
            values.push(poly::from_roots(&f, &neighbors));
        }
        let size = self.ell + 2;
        let mut c = Vec::with_capacity(size);
        for i in 0..size {
            c.push(super::lagrange(
                &f,
                &js[..size],
                &values[..size].iter().map(|v| v[i]).collect::<Vec<_>>(),
            ));
        }
        for i in 0..size {
            for j in 0..size {
                if c[i][j] != c[j][i] {
                    return Err("CM interpolation failed symmetry".into());
                }
            }
            for (j, val) in js.iter().zip(&values) {
                if poly::eval(&f, &c[i], *j) != val[i] {
                    return Err("CM interpolation failed validation point".into());
                }
            }
        }
        if c[self.ell + 1][0] != 1 || c[self.ell][self.ell] != p - 1 {
            return Err("CM interpolation failed leading-coefficient checks".into());
        }
        Ok(c)
    }

    fn walk_curves(
        &self,
        f: &Zp,
        t: u64,
        seed: Curve<u64>,
        walk: &ClassWalk,
    ) -> Result<Vec<Curve<u64>>, String> {
        let mut curves = vec![seed];
        for &(parent, gidx) in &walk.steps {
            let g = &self.generators[gidx];
            // Both orders give the same eigenvalue under the paired ideals:
            // b_floor=ell*b_top and pi=(t+v*ell*sqrt(D))/2.
            let q = g.q as i128;
            let lambda = ((t as i128 - self.v as i128 * self.ell as i128 * g.b) * (q + 1) / 2)
                .rem_euclid(q) as u64;
            curves.push(eigenstep(f, &curves[parent], g.q, lambda)?);
        }
        Ok(curves)
    }

    /// Stream one finite-field matrix at a time into explicit CRT modulo m>0.
    /// The modulus may be composite or share factors with the CRT primes.
    pub fn phi_mod_m(&self, modulus: &Int, rng: &mut Rng) -> Result<Vec<Vec<Int>>, String> {
        let primes = self.crt_primes()?;
        let size = self.ell + 2;
        let mut crt = ExplicitCrt::new(
            &primes.iter().map(|&(p, _)| p).collect::<Vec<_>>(),
            modulus,
            size * size,
            coefficient_bound_bits(self.ell)?,
        )?;
        for (p, t) in primes {
            let c = self
                .phi_mod_p(p, t, rng)
                .map_err(|e| format!("CM prime {p}, trace {t}: {e}"))?;
            crt.push(p, &c.into_iter().flatten().collect::<Vec<_>>())?;
        }
        Ok(crt.finish()?.chunks(size).map(|r| r.to_vec()).collect())
    }

    /// Recover signed integer coefficients using a target larger than twice B.
    pub fn phi_integer(&self, rng: &mut Rng) -> Result<Vec<Vec<Int>>, String> {
        let half = Int::one().shl(coefficient_bound_bits(self.ell)? + 1);
        let modulus = half.shl(1);
        Ok(self
            .phi_mod_m(&modulus, rng)?
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|c| if c > half { c - &modulus } else { c })
                    .collect()
            })
            .collect())
    }
}

// The exact order of each sampled point divides #E. Full rational E[ell]
// implies ell divides the smaller invariant factor, so ell*lcm(point orders)
// also divides #E. A divisor >4*sqrt(p) makes N unique in the Hasse interval.
// Sampling affects termination only; it cannot cause a false certificate.
fn certify_trace(f: &Zp, e: &Curve<u64>, n: u64, ell: u64, rng: &mut Rng) -> Result<(), String> {
    if f.p < 1 << 22 {
        return if curve::order(f, e, rng) == n {
            Ok(())
        } else {
            Err("CM seed has the wrong trace".into())
        };
    }
    let mut factors = vec![];
    let mut rest = n;
    let mut q = 2;
    while q * q <= rest {
        if rest.is_multiple_of(q) {
            factors.push(q);
            while rest.is_multiple_of(q) {
                rest /= q;
            }
        }
        q += if q == 2 { 1 } else { 2 };
    }
    if rest > 1 {
        factors.push(rest);
    }
    let mut multiple = ell;
    for _ in 0..64 {
        let pt = e.random_point(f, rng);
        if curve::pmul(f, e, &pt, n as u128) != Pt::Inf {
            return Err("CM seed has the wrong trace".into());
        }
        let mut order = n;
        for &q in &factors {
            while order.is_multiple_of(q) && curve::pmul(f, e, &pt, (order / q) as u128) == Pt::Inf
            {
                order /= q;
            }
        }
        let (mut a, mut b) = (multiple, order);
        while b != 0 {
            (a, b) = (b, a % b);
        }
        multiple = multiple / a * order;
        let divisor = ell as u128 * multiple as u128;
        if divisor * divisor > 16 * f.p as u128 {
            return Ok(());
        }
    }
    Err("trace certificate sampling budget exhausted".into())
}

fn validate_degree(ell: usize) -> Result<(), String> {
    if !(3..=31).contains(&ell) || !is_prime(ell as u64) {
        return Err("CM implementation supports odd prime degrees 3..=31".into());
    }
    Ok(())
}

/// Integer upper bound in bits from log |Phi_ell| <= 6 ell log ell + 18 ell.
/// ceil(log2 ell) and 27 > 18/log(2) avoid floating-point underestimation.
pub fn coefficient_bound_bits(ell: usize) -> Result<usize, String> {
    validate_degree(ell)?;
    Ok(6 * ell * (usize::BITS - (ell - 1).leading_zeros()) as usize + 27 * ell)
}

fn full_torsion(f: &Zp, e: &Curve<u64>, ell: usize) -> bool {
    let psi = if ell == 2 {
        vec![e.b, e.a, 0, 1]
    } else {
        poly::monic(f, &division_poly(f, e, ell))
    };
    let x = poly::x_poly(f);
    if !poly::sub(f, &poly::powmod(f, &x, f.p as u128, &psi), &x).is_empty() {
        return false;
    }
    ell == 2 || poly::powmod(f, &vec![e.b, e.a, 0, 1], ((f.p - 1) / 2) as u128, &psi) == vec![1]
}

fn eigenstep(f: &Zp, e: &Curve<u64>, q: u64, lambda: u64) -> Result<Curve<u64>, String> {
    if lambda == 0 {
        return Err("zero Frobenius eigenvalue".into());
    }
    let psi = poly::monic(f, &division_poly(f, e, q as usize));
    let (a, b) = ring_mul(f, e, &psi, lambda).ok_or("Frobenius scalar multiplication failed")?;
    let xp = poly::powmod(f, &poly::x_poly(f), f.p as u128, &psi);
    let yp = poly::powmod(f, &vec![e.b, e.a, 0, 1], ((f.p - 1) / 2) as u128, &psi);
    let h = poly::gcd(
        f,
        &poly::gcd(f, &psi, &poly::sub(f, &xp, &a)),
        &poly::sub(f, &yp, &b),
    );
    if poly::deg(f, &h) != (q as isize - 1) / 2 {
        return Err(format!(
            "incorrect eigenkernel degree at auxiliary prime {q}"
        ));
    }
    let cod = kohel(f, e, &h, q).cod;
    if !curve::is_smooth(f, &cod) {
        return Err("singular class-action codomain".into());
    }
    Ok(cod)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::find::modpoly::{integer_coeffs, Phi};

    #[test]
    fn talk_fixture_and_noncyclic_surface_match_hecke() {
        let plan = CmPlan::with_discriminant(5, -151).unwrap();
        assert_eq!(plan.class_numbers(), (7, 28));
        // The talk's example: 4*4451 = 52^2 + 2^2*5^2*151.
        for seed in [1, 17, 0x5107] {
            assert_eq!(
                plan.phi_mod_p(4451, 52, &mut Rng::new(seed)).unwrap(),
                Phi::compute(&Zp::new(4451), 5).c
            );
        }
        // h(-231)=12 with a noncyclic class group, and extra interpolation points.
        let plan = CmPlan::with_discriminant(5, -231).unwrap();
        assert_eq!(plan.class_numbers(), (12, 48));
        assert!(plan
            .top
            .classes
            .iter()
            .all(|c| c.pow(6) == Form::identity(-231)));
        let (p, t) = plan.crt_primes().unwrap()[0];
        assert_eq!(
            plan.phi_mod_p(p, t, &mut Rng::new(3)).unwrap(),
            Phi::compute(&Zp::new(p), 5).c
        );
    }

    #[test]
    fn selected_degrees_and_v_one_match_hecke() {
        for ell in [3, 7, 11, 13, 17, 19, 23, 29, 31] {
            let plan = CmPlan::new(ell).unwrap();
            for (p, t) in plan
                .crt_primes()
                .unwrap()
                .into_iter()
                .take(if ell <= 11 { 2 } else { 1 })
            {
                assert_eq!(
                    plan.phi_mod_p(p, t, &mut Rng::new(p)).unwrap(),
                    Phi::compute(&Zp::new(p), ell).c,
                    "ell={ell} D={} p={p}",
                    plan.d
                );
            }
        }
        let plan = CmPlan::with_discriminant(3, -131).unwrap();
        assert_eq!(plan.v, 1);
        let (p, t) = plan.crt_primes().unwrap()[0];
        assert_eq!(
            plan.phi_mod_p(p, t, &mut Rng::new(9)).unwrap(),
            Phi::compute(&Zp::new(p), 3).c
        );
    }

    #[test]
    fn large_prime_uses_exact_trace_certificate() {
        let plan = CmPlan::with_discriminant(5, -151).unwrap();
        let (p, t) = (4102..20000u64)
            .step_by(5)
            .filter_map(|t| {
                let norm = t * t + 4 * 25 * 151;
                let p = norm / 4;
                (norm.is_multiple_of(4) && plan.check_prime(p, t).is_ok()).then_some((p, t))
            })
            .next()
            .unwrap();
        assert!(p >= 1 << 22);
        assert_eq!(
            plan.phi_mod_p(p, t, &mut Rng::new(123)).unwrap(),
            Phi::compute(&Zp::new(p), 5).c
        );
        let f = Zp::new(p);
        let j = poly::roots(
            &f,
            &plan.h.iter().map(|c| c.mod_u64(p)).collect(),
            &mut Rng::new(1),
        )[0];
        let e = curve::from_j(&f, j);
        assert!(certify_trace(&f, &e, 1, 5, &mut Rng::new(1)).is_err());
    }

    #[test]
    fn explicit_cm_crt_matches_integer_q_expansion() {
        for ell in [3, 5] {
            let plan = CmPlan::new(ell).unwrap();
            let p = plan.crt_primes().unwrap()[0].0;
            let expected = integer_coeffs(ell);
            // A composite target sharing a factor with the CRT product.
            let modulus = Int::from(p * 12);
            let actual = plan.phi_mod_m(&modulus, &mut Rng::new(7)).unwrap();
            let expected: Vec<Vec<_>> = expected
                .iter()
                .map(|row| row.iter().map(|c| c.modulo(&modulus)).collect())
                .collect();
            assert_eq!(actual, expected, "ell={ell}");
        }
    }

    #[test]
    fn rejects_invalid_parameters_and_wrong_norm() {
        for ell in [0, 1, 2, 4, 33, usize::MAX] {
            assert!(CmPlan::new(ell).is_err());
        }
        for d in [i64::MIN, -3, 0, -15, -23, -99, -150, -1_000_000] {
            assert!(CmPlan::with_discriminant(5, d).is_err(), "D={d}");
        }
        let plan = CmPlan::with_discriminant(5, -151).unwrap();
        for (p, t) in [
            (4451, 51),
            (4453, 52),
            (MAX_CM_PRIME + 2, 52),
            (0, u64::MAX),
        ] {
            assert!(plan.phi_mod_p(p, t, &mut Rng::new(1)).is_err());
        }
        assert!(plan.phi_mod_m(&Int::zero(), &mut Rng::new(1)).is_err());
    }
}

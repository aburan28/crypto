//! sidh-toy: full SIDH key exchange in Rust.
//!
//!   p = 431 = 2^4 · 3^3 − 1   (textbook SIDH micro example)
//!   F_{p²} = F_p[i] / (i² + 1)
//!   E_0 : y² = x³ + 6x² + x       (supersingular over F_{p²})
//!   Alice torsion: 2^4 = 16       secret key m_A ∈ [0, 16)
//!   Bob   torsion: 3^3 = 27       secret key m_B ∈ [0, 27)
//!
//! Protocol:
//!   keygen_A: pick m_A; R_A = P_A + m_A · Q_A (order 16);
//!             chain 4 deg-2 isogenies; output (A_A, φ_A(P_B), φ_A(Q_B), φ_A(P_B − Q_B))
//!   keygen_B: symmetric with 3-isogenies on 3^3 torsion
//!   shared_A: R_AB = φ_B(P_A) + m_A · φ_B(Q_A); chain 4 deg-2 isogenies; output j(E_AB)
//!   shared_B: R_BA = φ_A(P_B) + m_B · φ_A(Q_B); chain 3 deg-3 isogenies; output j(E_BA)
//!   key agreement ⇔ j(E_AB) == j(E_BA)
//!
//! Caveat: SIDH was broken by Castryck–Decru in 2022. See cd_attack_sketch.rs
//! for the algorithm document.

use num_bigint::{BigInt, RandBigInt, Sign};
use num_integer::Integer;
use num_traits::{One, Zero};
use rand::SeedableRng;
use rand::Rng;
use rand_chacha::ChaCha8Rng;

// ============================================================================
// F_p
// ============================================================================

#[derive(Clone)]
struct Fp { p: BigInt }

impl Fp {
    fn r(&self, x: &BigInt) -> BigInt {
        let r = x % &self.p;
        if r.sign() == Sign::Minus { r + &self.p } else { r }
    }
    fn add(&self, a: &BigInt, b: &BigInt) -> BigInt { self.r(&(a + b)) }
    fn sub(&self, a: &BigInt, b: &BigInt) -> BigInt { self.r(&(a - b)) }
    fn mul(&self, a: &BigInt, b: &BigInt) -> BigInt { self.r(&(a * b)) }
    fn inv(&self, a: &BigInt) -> BigInt {
        let eg = self.r(a).extended_gcd(&self.p);
        assert!(eg.gcd.is_one()); self.r(&eg.x)
    }
    fn is_square(&self, a: &BigInt) -> bool {
        if a.is_zero() { return true; }
        a.modpow(&((&self.p - 1) / 2), &self.p).is_one()
    }
    fn sqrt(&self, a: &BigInt) -> Option<BigInt> {
        if a.is_zero() { return Some(BigInt::zero()); }
        if !self.is_square(a) { return None; }
        Some(a.modpow(&((&self.p + 1) / 4), &self.p))
    }
}

// ============================================================================
// F_{p²} = F_p[i] / (i² + 1)
// ============================================================================

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct F2 { a: BigInt, b: BigInt }

#[derive(Clone)]
struct Fp2 { fp: Fp }

impl Fp2 {
    fn zero(&self) -> F2 { F2 { a: BigInt::zero(), b: BigInt::zero() } }
    fn one(&self) -> F2 { F2 { a: BigInt::one(), b: BigInt::zero() } }
    fn from_int(&self, n: i64) -> F2 { F2 { a: self.fp.r(&BigInt::from(n)), b: BigInt::zero() } }
    fn is_zero(&self, x: &F2) -> bool { x.a.is_zero() && x.b.is_zero() }

    fn add(&self, x: &F2, y: &F2) -> F2 {
        F2 { a: self.fp.add(&x.a, &y.a), b: self.fp.add(&x.b, &y.b) }
    }
    fn sub(&self, x: &F2, y: &F2) -> F2 {
        F2 { a: self.fp.sub(&x.a, &y.a), b: self.fp.sub(&x.b, &y.b) }
    }
    fn neg(&self, x: &F2) -> F2 {
        F2 { a: self.fp.sub(&BigInt::zero(), &x.a), b: self.fp.sub(&BigInt::zero(), &x.b) }
    }
    fn mul(&self, x: &F2, y: &F2) -> F2 {
        let ac = self.fp.mul(&x.a, &y.a);
        let bd = self.fp.mul(&x.b, &y.b);
        let ad = self.fp.mul(&x.a, &y.b);
        let bc = self.fp.mul(&x.b, &y.a);
        F2 { a: self.fp.sub(&ac, &bd), b: self.fp.add(&ad, &bc) }
    }
    fn sq(&self, x: &F2) -> F2 { self.mul(x, x) }
    fn inv(&self, x: &F2) -> F2 {
        let n = self.fp.add(&self.fp.mul(&x.a, &x.a), &self.fp.mul(&x.b, &x.b));
        let ni = self.fp.inv(&n);
        F2 {
            a: self.fp.mul(&x.a, &ni),
            b: self.fp.sub(&BigInt::zero(), &self.fp.mul(&x.b, &ni)),
        }
    }
    fn pow(&self, x: &F2, e: &BigInt) -> F2 {
        let mut r = self.one(); let mut b = x.clone(); let mut k = e.clone();
        while k.sign() == Sign::Plus {
            if k.is_odd() { r = self.mul(&r, &b); }
            b = self.sq(&b);
            k >>= 1;
        }
        r
    }
    fn is_square(&self, x: &F2) -> bool {
        if self.is_zero(x) { return true; }
        let exp = (&self.fp.p * &self.fp.p - BigInt::one()) / 2;
        let r = self.pow(x, &exp);
        r == self.one()
    }
    /// Square root in F_{p²}, p ≡ 3 (mod 4). Returns None if non-square.
    fn sqrt(&self, x: &F2) -> Option<F2> {
        if self.is_zero(x) { return Some(self.zero()); }
        if x.b.is_zero() {
            if let Some(rt) = self.fp.sqrt(&x.a) {
                return Some(F2 { a: rt, b: BigInt::zero() });
            }
            let neg = self.fp.sub(&BigInt::zero(), &x.a);
            let rt = self.fp.sqrt(&neg)?;
            return Some(F2 { a: BigInt::zero(), b: rt });
        }
        let two_inv = self.fp.inv(&BigInt::from(2));
        let norm = self.fp.add(&self.fp.mul(&x.a, &x.a), &self.fp.mul(&x.b, &x.b));
        let s = self.fp.sqrt(&norm)?;
        for sign in [s.clone(), self.fp.sub(&BigInt::zero(), &s)] {
            let half = self.fp.mul(&self.fp.add(&x.a, &sign), &two_inv);
            if let Some(a) = self.fp.sqrt(&half) {
                if a.is_zero() { continue; }
                let inv = self.fp.inv(&self.fp.mul(&BigInt::from(2), &a));
                let b = self.fp.mul(&x.b, &inv);
                let cand = F2 { a, b };
                if self.sq(&cand) == *x { return Some(cand); }
            }
        }
        None
    }
}

// ============================================================================
// Affine Montgomery point (with sign), for basis construction & scalar mult.
// ============================================================================

#[derive(Clone, Debug, PartialEq, Eq)]
enum Aff { Inf, P(F2, F2) }   // P(x, y), y² = x³ + Ax² + x

#[derive(Clone)]
struct AffMont { a: F2 }

impl AffMont {
    fn lift(&self, x: &F2, fp2: &Fp2) -> Option<Aff> {
        let x2 = fp2.sq(x);
        let x3 = fp2.mul(&x2, x);
        let rhs = fp2.add(&fp2.add(&x3, &fp2.mul(&self.a, &x2)), x);
        fp2.sqrt(&rhs).map(|y| Aff::P(x.clone(), y))
    }
    fn neg(&self, p: &Aff, fp2: &Fp2) -> Aff {
        match p { Aff::Inf => Aff::Inf, Aff::P(x, y) => Aff::P(x.clone(), fp2.neg(y)) }
    }
    fn add(&self, p: &Aff, q: &Aff, fp2: &Fp2) -> Aff {
        let (xp, yp, xq, yq) = match (p, q) {
            (Aff::Inf, _) => return q.clone(),
            (_, Aff::Inf) => return p.clone(),
            (Aff::P(xp, yp), Aff::P(xq, yq)) => (xp, yp, xq, yq),
        };
        if xp == xq {
            if fp2.add(yp, yq) == fp2.zero() { return Aff::Inf; }
            // Doubling: λ = (3x² + 2Ax + 1)/(2y)
            let num = fp2.add(
                &fp2.add(&fp2.mul(&fp2.from_int(3), &fp2.sq(xp)),
                         &fp2.mul(&fp2.from_int(2), &fp2.mul(&self.a, xp))),
                &fp2.one());
            let den = fp2.mul(&fp2.from_int(2), yp);
            let lambda = fp2.mul(&num, &fp2.inv(&den));
            let xr = fp2.sub(
                &fp2.sub(&fp2.sq(&lambda), &self.a),
                &fp2.mul(&fp2.from_int(2), xp));
            let yr = fp2.sub(&fp2.mul(&lambda, &fp2.sub(xp, &xr)), yp);
            Aff::P(xr, yr)
        } else {
            let lambda = fp2.mul(&fp2.sub(yq, yp), &fp2.inv(&fp2.sub(xq, xp)));
            let xr = fp2.sub(
                &fp2.sub(&fp2.sub(&fp2.sq(&lambda), &self.a), xp), xq);
            let yr = fp2.sub(&fp2.mul(&lambda, &fp2.sub(xp, &xr)), yp);
            Aff::P(xr, yr)
        }
    }
    fn mul(&self, k: &BigInt, p: &Aff, fp2: &Fp2) -> Aff {
        let mut r = Aff::Inf; let mut q = p.clone(); let mut k = k.clone();
        while k.sign() == Sign::Plus {
            if k.is_odd() { r = self.add(&r, &q, fp2); }
            q = self.add(&q, &q, fp2);
            k >>= 1;
        }
        r
    }
}

// ============================================================================
// X-only projective Montgomery — used for the isogeny chain.
// ============================================================================

#[derive(Clone, Debug, PartialEq, Eq)]
struct Proj { x: F2, z: F2 }

impl Proj {
    fn inf(fp2: &Fp2) -> Self { Proj { x: fp2.one(), z: fp2.zero() } }
    fn aff(x: F2, fp2: &Fp2) -> Self { Proj { x, z: fp2.one() } }
    fn is_inf(&self) -> bool { self.z.a.is_zero() && self.z.b.is_zero() }
    fn norm(&self, fp2: &Fp2) -> F2 { fp2.mul(&self.x, &fp2.inv(&self.z)) }
}

#[derive(Clone)]
struct Mont { a: F2 }

impl Mont {
    fn xdbl(&self, p: &Proj, fp2: &Fp2) -> Proj {
        if p.is_inf() { return Proj::inf(fp2); }
        let xpz = fp2.add(&p.x, &p.z);
        let xmz = fp2.sub(&p.x, &p.z);
        let u = fp2.sq(&xpz);
        let v = fp2.sq(&xmz);
        let d = fp2.sub(&u, &v);
        let x3 = fp2.mul(&u, &v);
        let a24 = fp2.mul(
            &fp2.add(&self.a, &fp2.from_int(2)),
            &fp2.inv(&fp2.from_int(4)));
        let z3 = fp2.mul(&d, &fp2.add(&v, &fp2.mul(&a24, &d)));
        Proj { x: x3, z: z3 }
    }
    fn xadd(&self, p: &Proj, q: &Proj, pmq: &Proj, fp2: &Fp2) -> Proj {
        if p.is_inf() { return q.clone(); }
        if q.is_inf() { return p.clone(); }
        let u = fp2.mul(&fp2.sub(&p.x, &p.z), &fp2.add(&q.x, &q.z));
        let v = fp2.mul(&fp2.add(&p.x, &p.z), &fp2.sub(&q.x, &q.z));
        Proj {
            x: fp2.mul(&pmq.z, &fp2.sq(&fp2.add(&u, &v))),
            z: fp2.mul(&pmq.x, &fp2.sq(&fp2.sub(&u, &v))),
        }
    }
    fn ladder(&self, k: &BigInt, p: &Proj, fp2: &Fp2) -> Proj {
        if k.is_zero() || p.is_inf() { return Proj::inf(fp2); }
        let one = BigInt::one();
        let mut r0 = Proj::inf(fp2); let mut r1 = p.clone();
        for i in (0..k.bits()).rev() {
            let bit = ((k >> i) & &one).is_one();
            if bit {
                r0 = self.xadd(&r0, &r1, p, fp2);
                r1 = self.xdbl(&r1, fp2);
            } else {
                r1 = self.xadd(&r0, &r1, p, fp2);
                r0 = self.xdbl(&r0, fp2);
            }
        }
        r0
    }
}

// ============================================================================
// Isogenies
// ============================================================================

/// Degree-2 isogeny.  Kernel must be (X_K : Z_K) of exact order 2 with
/// X_K ≠ 0 (i.e., not the special (0:1) 2-torsion).
/// Codomain (SIKE / CLN 2016): A' = 2 − 4·x_K².
/// Push (Costello-Longa-Naehrig form):
///   φ(P) = (X_P·(X_P·X_K − Z_P·Z_K) : Z_P·(X_P·Z_K − Z_P·X_K))
fn isog2(_curve: &Mont, kernel: &Proj, push: &mut [Proj], fp2: &Fp2) -> F2 {
    let x_k = kernel.norm(fp2);
    let xk_sq = fp2.sq(&x_k);
    let a_new = fp2.sub(&fp2.from_int(2), &fp2.mul(&fp2.from_int(4), &xk_sq));
    for p in push.iter_mut() {
        if p.is_inf() { continue; }
        let xx = fp2.mul(&p.x, &kernel.x);
        let zz = fp2.mul(&p.z, &kernel.z);
        let xz = fp2.mul(&p.x, &kernel.z);
        let zx = fp2.mul(&p.z, &kernel.x);
        let new_x = fp2.mul(&p.x, &fp2.sub(&xx, &zz));
        let new_z = fp2.mul(&p.z, &fp2.sub(&xz, &zx));
        *p = Proj { x: new_x, z: new_z };
    }
    a_new
}

/// Odd-ℓ isogeny (works for any odd ℓ ≥ 3). Same algorithm as csidh_toy:
/// Mont↔Edwards conversion, Moody–Shumow codomain, Costello-Hisil point eval.
fn isog_odd(
    curve: &Mont, kernel: &Proj, ell: u64, push: &mut [Proj], fp2: &Fp2,
) -> F2 {
    assert!(ell >= 3 && ell % 2 == 1);
    let two = fp2.from_int(2);
    let d_x = fp2.sub(&curve.a, &two);
    let d_z = fp2.add(&curve.a, &two);
    let mut prod_x = fp2.sub(&kernel.x, &kernel.z);
    let mut prod_z = fp2.add(&kernel.x, &kernel.z);

    // accumulators for each pushed point
    let mut acc: Vec<Option<(F2, F2)>> = push.iter().map(|p| {
        if p.is_inf() { None } else {
            Some((
                fp2.sub(&fp2.mul(&p.x, &kernel.x), &fp2.mul(&p.z, &kernel.z)),
                fp2.sub(&fp2.mul(&p.x, &kernel.z), &fp2.mul(&p.z, &kernel.x)),
            ))
        }
    }).collect();

    let mut m = vec![kernel.clone(), curve.xdbl(kernel, fp2), Proj::inf(fp2)];
    for i in 1..(ell / 2) {
        let i_mod = (i % 3) as usize;
        let im1 = ((i - 1) % 3) as usize;
        let im2 = ((i + 1) % 3) as usize;
        if i >= 2 { m[i_mod] = curve.xadd(&m[im1], kernel, &m[im2], fp2); }
        let m_i = m[i_mod].clone();
        prod_x = fp2.mul(&prod_x, &fp2.sub(&m_i.x, &m_i.z));
        prod_z = fp2.mul(&prod_z, &fp2.add(&m_i.x, &m_i.z));
        for (a, p) in acc.iter_mut().zip(push.iter()) {
            if let Some((ax, az)) = a {
                let t0 = fp2.mul(&fp2.sub(&p.x, &p.z), &fp2.add(&m_i.x, &m_i.z));
                let t1 = fp2.mul(&fp2.sub(&m_i.x, &m_i.z), &fp2.add(&p.x, &p.z));
                *ax = fp2.mul(ax, &fp2.add(&t0, &t1));
                *az = fp2.mul(az, &fp2.sub(&t0, &t1));
            }
        }
    }
    let ell_big = BigInt::from(ell);
    let eight = BigInt::from(8);
    let new_d_x = fp2.mul(&fp2.pow(&d_x, &ell_big), &fp2.pow(&prod_x, &eight));
    let new_d_z = fp2.mul(&fp2.pow(&d_z, &ell_big), &fp2.pow(&prod_z, &eight));
    let a_new = fp2.mul(
        &fp2.mul(&two, &fp2.add(&new_d_x, &new_d_z)),
        &fp2.inv(&fp2.sub(&new_d_z, &new_d_x)));
    for (a, p) in acc.iter().zip(push.iter_mut()) {
        if let Some((ax, az)) = a {
            *p = Proj {
                x: fp2.mul(&p.x, &fp2.sq(ax)),
                z: fp2.mul(&p.z, &fp2.sq(az)),
            };
        }
    }
    a_new
}

// ============================================================================
// Torsion basis & SIDH protocol
// ============================================================================

fn j_inv(a: &F2, fp2: &Fp2) -> F2 {
    let a2 = fp2.sq(a);
    let inner = fp2.sub(&a2, &fp2.from_int(3));
    let cube = fp2.mul(&inner, &fp2.sq(&inner));
    let num = fp2.mul(&fp2.from_int(256), &cube);
    let den = fp2.sub(&a2, &fp2.from_int(4));
    fp2.mul(&num, &fp2.inv(&den))
}

/// Build a basis {P, Q} of E_0[N] with N = ℓ^e | p+1, plus P − Q for the
/// 3-point ladder analog. Returns (P, Q, P-Q) as affine points.
fn make_basis(
    curve: &AffMont, big_cofactor: &BigInt, n: &BigInt, fp2: &Fp2, rng: &mut ChaCha8Rng,
) -> (Aff, Aff, Aff) {
    let find_order_n = |rng: &mut ChaCha8Rng| -> Aff {
        loop {
            let xa = rng.gen_bigint_range(&BigInt::zero(), &fp2.fp.p);
            let xb = rng.gen_bigint_range(&BigInt::zero(), &fp2.fp.p);
            let x = F2 { a: xa, b: xb };
            if let Some(pt) = curve.lift(&x, fp2) {
                let scaled = curve.mul(big_cofactor, &pt, fp2);
                if matches!(scaled, Aff::Inf) { continue; }
                // Confirm order divides N but isn't smaller — check [N/ℓ]·scaled ≠ ∞ for the prime ℓ.
                let full = curve.mul(n, &scaled, fp2);
                if !matches!(full, Aff::Inf) { continue; }
                return scaled;
            }
        }
    };

    let p_pt = find_order_n(rng);
    // Find Q linearly independent from P: brute-force check no scalar k∈[1,N) has [k]·P = Q.
    loop {
        let q_pt = find_order_n(rng);
        let mut indep = true;
        let mut acc = Aff::Inf;
        let mut k = BigInt::zero();
        while &k < n {
            if acc == q_pt { indep = false; break; }
            acc = curve.add(&acc, &p_pt, fp2);
            k += 1;
        }
        if !indep { continue; }
        let pmq = curve.add(&p_pt, &curve.neg(&q_pt, fp2), fp2);
        return (p_pt, q_pt, pmq);
    }
}

/// Alice's keygen: m_A ∈ [0, 2^a). Returns (A_A, x_φ(P_B), x_φ(Q_B), x_φ(P_B − Q_B)).
fn keygen_alice(
    e0_aff: &AffMont, e0_mont: &Mont,
    p_a: &Aff, q_a: &Aff, _pmq_a: &Aff,
    p_b: &Aff, q_b: &Aff, pmq_b: &Aff,
    m_a: &BigInt, a_steps: u32, fp2: &Fp2,
) -> (F2, Proj, Proj, Proj) {
    // R_A = P_A + m_A · Q_A  (in affine, then go x-only)
    let mq = e0_aff.mul(m_a, q_a, fp2);
    let r_a = e0_aff.add(p_a, &mq, fp2);
    let r_a_x = match &r_a {
        Aff::Inf => panic!("R_A is infinity"),
        Aff::P(x, _) => Proj::aff(x.clone(), fp2),
    };
    // Auxiliary points pushed through the chain (x-only).
    let mut push_pts: Vec<Proj> = vec![
        r_a_x,
        match p_b { Aff::P(x,_) => Proj::aff(x.clone(), fp2), _ => Proj::inf(fp2) },
        match q_b { Aff::P(x,_) => Proj::aff(x.clone(), fp2), _ => Proj::inf(fp2) },
        match pmq_b { Aff::P(x,_) => Proj::aff(x.clone(), fp2), _ => Proj::inf(fp2) },
    ];
    let mut cur = e0_mont.clone();
    for step in 0..a_steps {
        // Extract a kernel of order exactly 2 from the first push point R_A.
        let pow = 1u64 << (a_steps - step - 1);
        let small_k = cur.ladder(&BigInt::from(pow), &push_pts[0], fp2);
        let a_new = isog2(&cur, &small_k, &mut push_pts, fp2);
        cur = Mont { a: a_new };
    }
    // After chain, push_pts[0] (the R_A kernel) is annihilated. The other three
    // are the images of P_B, Q_B, P_B − Q_B on Alice's curve.
    (cur.a, push_pts[1].clone(), push_pts[2].clone(), push_pts[3].clone())
}

/// Bob's keygen: m_B ∈ [0, 3^b). Returns (A_B, x_φ(P_A), x_φ(Q_A), x_φ(P_A − Q_A)).
fn keygen_bob(
    e0_aff: &AffMont, e0_mont: &Mont,
    p_a: &Aff, q_a: &Aff, pmq_a: &Aff,
    p_b: &Aff, q_b: &Aff, _pmq_b: &Aff,
    m_b: &BigInt, b_steps: u32, fp2: &Fp2,
) -> (F2, Proj, Proj, Proj) {
    let mq = e0_aff.mul(m_b, q_b, fp2);
    let r_b = e0_aff.add(p_b, &mq, fp2);
    let r_b_x = match &r_b {
        Aff::Inf => panic!("R_B is infinity"),
        Aff::P(x, _) => Proj::aff(x.clone(), fp2),
    };
    let mut push_pts: Vec<Proj> = vec![
        r_b_x,
        match p_a { Aff::P(x,_) => Proj::aff(x.clone(), fp2), _ => Proj::inf(fp2) },
        match q_a { Aff::P(x,_) => Proj::aff(x.clone(), fp2), _ => Proj::inf(fp2) },
        match pmq_a { Aff::P(x,_) => Proj::aff(x.clone(), fp2), _ => Proj::inf(fp2) },
    ];
    let mut cur = e0_mont.clone();
    for step in 0..b_steps {
        let pow = 3u64.pow(b_steps - step - 1);
        let small_k = cur.ladder(&BigInt::from(pow), &push_pts[0], fp2);
        let a_new = isog_odd(&cur, &small_k, 3, &mut push_pts, fp2);
        cur = Mont { a: a_new };
    }
    (cur.a, push_pts[1].clone(), push_pts[2].clone(), push_pts[3].clone())
}

/// Alice's shared-secret pass: reconstruct kernel on Bob's curve, chain a 2-isogenies.
/// Requires knowing the y-coordinate of (φ_B(P_A) + m_A · φ_B(Q_A)), which from x-only
/// we recover by lifting on E_B.
fn shared_alice(
    a_b: &F2,
    phi_pa: &Proj, phi_qa: &Proj, phi_pmqa: &Proj,
    m_a: &BigInt, a_steps: u32, fp2: &Fp2,
) -> F2 {
    let cur_aff = AffMont { a: a_b.clone() };
    let xa_p = phi_pa.norm(fp2);
    let xa_q = phi_qa.norm(fp2);
    let xa_pmq = phi_pmqa.norm(fp2);
    let pa_lift = cur_aff.lift(&xa_p, fp2).expect("φ(P_A) didn't lift on E_B");
    let qa_lift = cur_aff.lift(&xa_q, fp2).expect("φ(Q_A) didn't lift");
    let pmq_lift = cur_aff.lift(&xa_pmq, fp2).expect("φ(P−Q) didn't lift");
    // Need (P + m·Q) on E_B with y-sign matching the published P−Q relation.
    // Resolve sign of Q_A by checking that (P_A - Q_A) matches.
    let try_q = |q: Aff| -> Option<Aff> {
        let diff = cur_aff.add(&pa_lift, &cur_aff.neg(&q, fp2), fp2);
        if let (Aff::P(x1, _), Aff::P(x2, _)) = (&diff, &pmq_lift) {
            if x1 == x2 { return Some(q); }
        }
        None
    };
    let qa_correct = try_q(qa_lift.clone()).unwrap_or_else(|| {
        try_q(cur_aff.neg(&qa_lift, fp2)).expect("Q_A sign disambiguation failed")
    });
    let r = cur_aff.add(&pa_lift, &cur_aff.mul(m_a, &qa_correct, fp2), fp2);
    let r_x = match r { Aff::P(x,_) => Proj::aff(x, fp2), _ => panic!("kernel ∞") };

    let mut push_pts = vec![r_x];
    let mut cur = Mont { a: a_b.clone() };
    for step in 0..a_steps {
        let pow = 1u64 << (a_steps - step - 1);
        let small_k = cur.ladder(&BigInt::from(pow), &push_pts[0], fp2);
        let a_new = isog2(&cur, &small_k, &mut push_pts, fp2);
        cur = Mont { a: a_new };
    }
    j_inv(&cur.a, fp2)
}

fn shared_bob(
    a_a: &F2,
    phi_pb: &Proj, phi_qb: &Proj, phi_pmqb: &Proj,
    m_b: &BigInt, b_steps: u32, fp2: &Fp2,
) -> F2 {
    let cur_aff = AffMont { a: a_a.clone() };
    let xb_p = phi_pb.norm(fp2);
    let xb_q = phi_qb.norm(fp2);
    let xb_pmq = phi_pmqb.norm(fp2);
    let pb_lift = cur_aff.lift(&xb_p, fp2).expect("φ(P_B) didn't lift");
    let qb_lift = cur_aff.lift(&xb_q, fp2).expect("φ(Q_B) didn't lift");
    let pmq_lift = cur_aff.lift(&xb_pmq, fp2).expect("φ(P_B−Q_B) didn't lift");
    let try_q = |q: Aff| -> Option<Aff> {
        let diff = cur_aff.add(&pb_lift, &cur_aff.neg(&q, fp2), fp2);
        if let (Aff::P(x1,_), Aff::P(x2,_)) = (&diff, &pmq_lift) {
            if x1 == x2 { return Some(q); }
        }
        None
    };
    let qb_correct = try_q(qb_lift.clone()).unwrap_or_else(|| {
        try_q(cur_aff.neg(&qb_lift, fp2)).expect("Q_B sign disambiguation failed")
    });
    let r = cur_aff.add(&pb_lift, &cur_aff.mul(m_b, &qb_correct, fp2), fp2);
    let r_x = match r { Aff::P(x,_) => Proj::aff(x, fp2), _ => panic!("kernel ∞") };

    let mut push_pts = vec![r_x];
    let mut cur = Mont { a: a_a.clone() };
    for step in 0..b_steps {
        let pow = 3u64.pow(b_steps - step - 1);
        let small_k = cur.ladder(&BigInt::from(pow), &push_pts[0], fp2);
        let a_new = isog_odd(&cur, &small_k, 3, &mut push_pts, fp2);
        cur = Mont { a: a_new };
    }
    j_inv(&cur.a, fp2)
}

// ============================================================================
// Main
// ============================================================================

fn main() {
    let p = BigInt::from(431u64);
    let fp = Fp { p: p.clone() };
    let fp2 = Fp2 { fp };
    let a_exp = 4u32; let b_exp = 3u32;
    let n_a = BigInt::from(1u64 << a_exp);          // 16
    let n_b = BigInt::from(3u64.pow(b_exp));        // 27
    let cof_a = BigInt::from(3u64.pow(b_exp));      // = 27, cofactor for 2^a torsion
    let cof_b = BigInt::from(1u64 << a_exp);        // = 16, cofactor for 3^b torsion
    println!("[*] SIDH-toy");
    println!("    p = {}  (= 2^{a_exp} · 3^{b_exp} − 1),   p+1 = {}", p, &p + 1);
    println!("    Alice torsion: 2^{a_exp} = {}", n_a);
    println!("    Bob   torsion: 3^{b_exp} = {}", n_b);
    println!();

    // E_0: y² = x³ + 6x² + x
    let e0_aff = AffMont { a: fp2.from_int(6) };
    let e0_mont = Mont { a: fp2.from_int(6) };
    println!("[*] Starting curve E_0:  A = 6   j(E_0) = {:?}", j_inv(&e0_mont.a, &fp2));

    let mut rng = ChaCha8Rng::seed_from_u64(0xD15EA5E);
    let (p_a, q_a, pmq_a) = make_basis(&e0_aff, &cof_a, &n_a, &fp2, &mut rng);
    let (p_b, q_b, pmq_b) = make_basis(&e0_aff, &cof_b, &n_b, &fp2, &mut rng);
    println!("[*] Built torsion bases P_A, Q_A (order {}) and P_B, Q_B (order {}).", n_a, n_b);
    println!();

    // Secret keys
    let mut rng_a = ChaCha8Rng::seed_from_u64(0xA11CE);
    let mut rng_b = ChaCha8Rng::seed_from_u64(0xB0B);
    let m_a = BigInt::from(rng_a.gen_range(1u64..(1u64 << a_exp)));
    let m_b = BigInt::from(rng_b.gen_range(1u64..3u64.pow(b_exp)));
    println!("[*] Alice secret m_A = {}", m_a);
    println!("[*] Bob   secret m_B = {}", m_b);

    let (a_a, phi_pb, phi_qb, phi_pmqb) = keygen_alice(
        &e0_aff, &e0_mont, &p_a, &q_a, &pmq_a, &p_b, &q_b, &pmq_b,
        &m_a, a_exp, &fp2);
    let (a_b, phi_pa, phi_qa, phi_pmqa) = keygen_bob(
        &e0_aff, &e0_mont, &p_a, &q_a, &pmq_a, &p_b, &q_b, &pmq_b,
        &m_b, b_exp, &fp2);
    println!();
    println!("[*] Alice public key:");
    println!("    E_A: A = {:?}", a_a);
    println!("    j(E_A) = {:?}", j_inv(&a_a, &fp2));
    println!("[*] Bob   public key:");
    println!("    E_B: A = {:?}", a_b);
    println!("    j(E_B) = {:?}", j_inv(&a_b, &fp2));
    println!();

    let j_alice = shared_alice(&a_b, &phi_pa, &phi_qa, &phi_pmqa, &m_a, a_exp, &fp2);
    let j_bob = shared_bob(&a_a, &phi_pb, &phi_qb, &phi_pmqb, &m_b, b_exp, &fp2);
    println!("[*] Alice's shared j = {:?}", j_alice);
    println!("[*] Bob's   shared j = {:?}", j_bob);
    println!("[*] AGREE: {}", j_alice == j_bob);
}

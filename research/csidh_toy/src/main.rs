//! csidh-toy: a small CSIDH key exchange implementation.
//!
//! Math:
//!   - p = 4·ℓ_1·ℓ_2·...·ℓ_n − 1, with all ℓ_i small distinct odd primes
//!   - p ≡ 3 (mod 4)
//!   - E_0: y² = x³ + x is supersingular over F_p, in Montgomery form A=0
//!   - All curves in the action graph are supersingular over F_p with #E = p+1
//!   - Secret key: vector (e_1,...,e_n) ∈ [-B, B]^n
//!   - Action of (e_i): apply |e_i| ℓ_i-isogenies in the direction sign(e_i)
//!   - Public key: the A-coefficient of the resulting Montgomery curve
//!   - Shared secret: j-invariant of the doubly-acted curve (commutative action)
//!
//! All EC arithmetic is Montgomery x-only on (X : Z) projective coords.
//! Isogenies use the Costello–Hisil 2017 explicit formula for odd-degree
//! Montgomery isogenies. Affine normalization at each step keeps code simple.

use num_bigint::{BigInt, RandBigInt, Sign};
use num_integer::Integer;
use num_traits::{One, Zero};
use rand::SeedableRng;
use rand::Rng;
use rand_chacha::ChaCha8Rng;

// ============================================================================
// Parameters
// ============================================================================

/// Small odd primes used in CSIDH. p+1 = 4·∏LIST.
/// 12 primes gives p ≈ 2^34 — small enough for fast BigInt math, big enough
/// that the action graph has plenty of curves.
const LIST: &[u64] = &[3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];

/// Secret key coordinate bound: e_i ∈ [-BOUND, BOUND].
const BOUND: i32 = 5;

// ============================================================================
// F_p arithmetic
// ============================================================================

#[derive(Clone)]
struct Fp {
    p: BigInt,
}

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
        assert!(eg.gcd.is_one(), "no inverse");
        self.r(&eg.x)
    }
    /// p ≡ 3 (mod 4): legendre symbol via a^((p−1)/2).
    fn is_square(&self, a: &BigInt) -> bool {
        if a.is_zero() { return true; }
        let leg = a.modpow(&((&self.p - 1) / 2), &self.p);
        leg.is_one()
    }
}

// ============================================================================
// Primality test (Miller–Rabin) — used to verify p is prime.
// ============================================================================

fn miller_rabin(n: &BigInt) -> bool {
    let one = BigInt::one();
    if n <= &one { return false; }
    if n == &BigInt::from(2) || n == &BigInt::from(3) { return true; }
    if n.is_even() { return false; }
    let n_minus_1: BigInt = n - 1;
    let mut d = n_minus_1.clone();
    let mut r: u32 = 0;
    while d.is_even() { d /= 2; r += 1; }
    for &w in &[2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let wb = BigInt::from(w);
        if &wb >= n { continue; }
        let mut x = wb.modpow(&d, n);
        if x.is_one() || x == n_minus_1 { continue; }
        let mut composite = true;
        for _ in 0..r.saturating_sub(1) {
            x = (&x * &x) % n;
            if x == n_minus_1 { composite = false; break; }
        }
        if composite { return false; }
    }
    true
}

// ============================================================================
// Montgomery curve x-only arithmetic
// ============================================================================

#[derive(Clone, Debug, PartialEq, Eq)]
struct Proj {
    x: BigInt,
    z: BigInt,
}

impl Proj {
    fn infinity() -> Self { Proj { x: BigInt::one(), z: BigInt::zero() } }
    fn affine(x: BigInt) -> Self { Proj { x, z: BigInt::one() } }
    fn is_inf(&self) -> bool { self.z.is_zero() }
}

/// Montgomery curve By² = x³ + Ax² + x, with B=1. We track only A.
#[derive(Clone)]
struct Mont {
    a: BigInt,
}

impl Mont {
    /// x-only doubling on M_A.
    fn xdbl(&self, p: &Proj, fp: &Fp) -> Proj {
        if p.is_inf() { return Proj::infinity(); }
        let xpz = fp.add(&p.x, &p.z);
        let xmz = fp.sub(&p.x, &p.z);
        let u = fp.mul(&xpz, &xpz);
        let v = fp.mul(&xmz, &xmz);
        let d = fp.sub(&u, &v);
        let x3 = fp.mul(&u, &v);
        let a24 = fp.mul(
            &fp.add(&self.a, &BigInt::from(2)),
            &fp.inv(&BigInt::from(4)),
        );
        let z3 = fp.mul(&d, &fp.add(&v, &fp.mul(&a24, &d)));
        Proj { x: x3, z: z3 }
    }

    /// Differential addition: P+Q given P, Q, P−Q (x-only).
    fn xadd(&self, p: &Proj, q: &Proj, pmq: &Proj, fp: &Fp) -> Proj {
        if p.is_inf() { return q.clone(); }
        if q.is_inf() { return p.clone(); }
        let u = fp.mul(&fp.sub(&p.x, &p.z), &fp.add(&q.x, &q.z));
        let v = fp.mul(&fp.add(&p.x, &p.z), &fp.sub(&q.x, &q.z));
        let upv = fp.add(&u, &v);
        let umv = fp.sub(&u, &v);
        let x = fp.mul(&pmq.z, &fp.mul(&upv, &upv));
        let z = fp.mul(&pmq.x, &fp.mul(&umv, &umv));
        Proj { x, z }
    }

    /// Montgomery ladder: compute [k]P (x-only) for k ≥ 0.
    fn ladder(&self, k: &BigInt, p: &Proj, fp: &Fp) -> Proj {
        if k.is_zero() || p.is_inf() { return Proj::infinity(); }
        let one = BigInt::one();
        let mut r0 = Proj::infinity();
        let mut r1 = p.clone();
        let bits = k.bits();
        for i in (0..bits).rev() {
            let bit = ((k >> i) & &one).is_one();
            if bit {
                r0 = self.xadd(&r0, &r1, p, fp);
                r1 = self.xdbl(&r1, fp);
            } else {
                r1 = self.xadd(&r0, &r1, p, fp);
                r0 = self.xdbl(&r0, fp);
            }
        }
        r0
    }

    /// Is `x` the x-coord of a point on M_A (rather than its twist)?
    fn is_on_curve(&self, x: &BigInt, fp: &Fp) -> bool {
        // rhs = x·(x² + A·x + 1)
        let x2 = fp.mul(x, x);
        let rhs = fp.mul(x, &fp.add(&fp.add(&x2, &fp.mul(&self.a, x)), &BigInt::one()));
        fp.is_square(&rhs)
    }
}

// ============================================================================
// ℓ-isogeny (odd ℓ), x-only.
// Algorithm: Montgomery↔Edwards conversion, Moody–Shumow Edwards-form
// codomain, Costello–Hisil Montgomery-form point evaluation.
// Matches the CSIDH reference C implementation (yx7.cc/code/csidh).
// ============================================================================

fn xisog(
    curve: &Mont,
    kernel: &Proj,
    ell: u64,
    push: Option<&Proj>,
    fp: &Fp,
) -> (BigInt, Option<Proj>) {
    assert!(ell >= 3 && ell % 2 == 1);
    let two = BigInt::from(2);

    // Mont → Edwards:  d = (A−2)/(A+2),  stored projectively as (d_x : d_z).
    let d_x = fp.sub(&curve.a, &two);
    let d_z = fp.add(&curve.a, &two);

    // Edwards-form codomain accumulator: starts with i=1 contribution from K.
    //   prod = ((X_K − Z_K) : (X_K + Z_K))
    let mut prod_x = fp.sub(&kernel.x, &kernel.z);
    let mut prod_z = fp.add(&kernel.x, &kernel.z);

    // Montgomery-form point evaluation: init Q from input P and K.
    //   Q = ( P.x·K.x − P.z·K.z  :  P.x·K.z − P.z·K.x )
    let mut q_acc = push.map(|p| Proj {
        x: fp.sub(&fp.mul(&p.x, &kernel.x), &fp.mul(&p.z, &kernel.z)),
        z: fp.sub(&fp.mul(&p.x, &kernel.z), &fp.mul(&p.z, &kernel.x)),
    });

    // M[0] = K, M[1] = 2K, M[2] holds (i+1)K computed during loop.
    let mut m: Vec<Proj> = vec![
        kernel.clone(),
        curve.xdbl(kernel, fp),
        Proj::infinity(),
    ];

    // Loop: for i = 1 .. (ℓ−1)/2 (exclusive upper, since C uses `i < ℓ/2`).
    // At iteration i, M[i%3] represents (i+1)·K.
    for i in 1..(ell / 2) {
        let i_mod = (i % 3) as usize;
        let im1_mod = ((i - 1) % 3) as usize;
        let im2_mod = ((i + 1) % 3) as usize; // (i-2) ≡ (i+1) (mod 3)

        if i >= 2 {
            m[i_mod] = curve.xadd(&m[im1_mod], kernel, &m[im2_mod], fp);
        }
        let m_i = m[i_mod].clone();

        // Consume into Edwards prod
        prod_x = fp.mul(&prod_x, &fp.sub(&m_i.x, &m_i.z));
        prod_z = fp.mul(&prod_z, &fp.add(&m_i.x, &m_i.z));

        // Consume into Montgomery point accumulator
        if let (Some(q), Some(pp)) = (q_acc.as_mut(), push) {
            let t0 = fp.mul(&fp.sub(&pp.x, &pp.z), &fp.add(&m_i.x, &m_i.z));
            let t1 = fp.mul(&fp.sub(&m_i.x, &m_i.z), &fp.add(&pp.x, &pp.z));
            q.x = fp.mul(&q.x, &fp.add(&t0, &t1));
            q.z = fp.mul(&q.z, &fp.sub(&t0, &t1));
        }
    }

    // Edwards codomain finish:  d_new = d^ℓ · prod^8
    let ell_big = BigInt::from(ell);
    let eight = BigInt::from(8);
    let new_d_x = fp.mul(
        &d_x.modpow(&ell_big, &fp.p),
        &prod_x.modpow(&eight, &fp.p),
    );
    let new_d_z = fp.mul(
        &d_z.modpow(&ell_big, &fp.p),
        &prod_z.modpow(&eight, &fp.p),
    );

    // Edwards → Mont:  A = 2(d.x + d.z)/(d.z − d.x)
    let a_new = fp.mul(
        &fp.mul(&two, &fp.add(&new_d_x, &new_d_z)),
        &fp.inv(&fp.sub(&new_d_z, &new_d_x)),
    );

    // Montgomery point finish:  P' = (P.x · Q.x² : P.z · Q.z²)
    let pushed = if let (Some(q), Some(pp)) = (q_acc.as_ref(), push) {
        let qx2 = fp.mul(&q.x, &q.x);
        let qz2 = fp.mul(&q.z, &q.z);
        Some(Proj {
            x: fp.mul(&pp.x, &qx2),
            z: fp.mul(&pp.z, &qz2),
        })
    } else {
        None
    };

    (a_new, pushed)
}

// ============================================================================
// CSIDH group action
// ============================================================================

#[derive(Clone)]
struct Csidh {
    fp: Fp,
    p_plus_1: BigInt,
    primes: Vec<u64>,
}

impl Csidh {
    fn new() -> Self {
        let mut p_plus_1 = BigInt::from(4);
        for &l in LIST { p_plus_1 *= l; }
        let p: BigInt = &p_plus_1 - 1;
        assert!(miller_rabin(&p), "p = 4·∏LIST − 1 is not prime; adjust LIST");
        Csidh { fp: Fp { p }, p_plus_1, primes: LIST.to_vec() }
    }

    fn keygen(&self, rng: &mut ChaCha8Rng) -> Vec<i32> {
        (0..self.primes.len()).map(|_| rng.gen_range(-BOUND..=BOUND)).collect()
    }

    /// Apply the action of secret vector `e` to curve M_A, return new A.
    fn action(&self, a: &BigInt, e: &[i32]) -> BigInt {
        let mut e = e.to_vec();
        let mut a = a.clone();

        // Deterministic sampling seeded from inputs, so the test is reproducible.
        let mut seed: u64 = 1469598103934665603;
        for &ei in &e { seed ^= ei as u64; seed = seed.wrapping_mul(1099511628211); }
        for byte in a.to_bytes_le().1 {
            seed ^= byte as u64; seed = seed.wrapping_mul(1099511628211);
        }
        let mut rng = ChaCha8Rng::seed_from_u64(seed);

        loop {
            if e.iter().all(|&v| v == 0) { return a; }
            let curve = Mont { a: a.clone() };

            let x = rng.gen_bigint_range(&BigInt::from(1), &self.fp.p);
            let s: i32 = if curve.is_on_curve(&x, &self.fp) { 1 } else { -1 };
            let active: Vec<usize> = (0..self.primes.len())
                .filter(|&i| e[i] != 0 && e[i].signum() == s)
                .collect();
            if active.is_empty() { continue; }
            let k_prod: BigInt = active.iter()
                .map(|&i| BigInt::from(self.primes[i])).product();

            let p_pt = Proj::affine(x);
            let mut q_pt = curve.ladder(&(&self.p_plus_1 / &k_prod), &p_pt, &self.fp);
            if q_pt.is_inf() { continue; }

            let mut k_remaining = k_prod.clone();
            let mut cur_curve = curve;
            for &i in &active {
                let ell = self.primes[i];
                let cofactor = &k_remaining / BigInt::from(ell);
                let kernel = cur_curve.ladder(&cofactor, &q_pt, &self.fp);
                if kernel.is_inf() {
                    k_remaining /= ell;
                    continue;
                }
                let (a_new, q_new) = xisog(&cur_curve, &kernel, ell, Some(&q_pt), &self.fp);
                cur_curve = Mont { a: a_new };
                q_pt = q_new.unwrap();
                e[i] -= s;
                k_remaining /= ell;
                if q_pt.is_inf() { break; }
            }
            a = cur_curve.a;
        }
    }

    /// j-invariant of M_A: j = 256·(A² − 3)³ / (A² − 4)
    fn j_invariant(&self, a: &BigInt) -> BigInt {
        let a2 = self.fp.mul(a, a);
        let three = BigInt::from(3);
        let four = BigInt::from(4);
        let inner = self.fp.sub(&a2, &three);
        let cube = self.fp.mul(&inner, &self.fp.mul(&inner, &inner));
        let num = self.fp.mul(&BigInt::from(256), &cube);
        let den = self.fp.sub(&a2, &four);
        self.fp.mul(&num, &self.fp.inv(&den))
    }
}

fn fmt_vec(v: &[i32]) -> String {
    let parts: Vec<String> = v.iter().map(|e| format!("{:>2}", e)).collect();
    format!("[{}]", parts.join(", "))
}

fn main() {
    let csidh = Csidh::new();
    println!("[*] CSIDH-toy parameters");
    println!("    primes: {:?}", csidh.primes);
    println!("    p = {}  ({} bits)", csidh.fp.p, csidh.fp.p.bits());
    println!("    p+1 = {}", csidh.p_plus_1);
    println!("    starting curve E_0: A = 0  (y² = x³ + x)");
    println!();

    let mut rng_alice = ChaCha8Rng::seed_from_u64(0xA11CE);
    let mut rng_bob = ChaCha8Rng::seed_from_u64(0xB0B);

    let sk_alice = csidh.keygen(&mut rng_alice);
    let sk_bob = csidh.keygen(&mut rng_bob);
    println!("[*] Alice secret key e_A = {}", fmt_vec(&sk_alice));
    println!("[*] Bob   secret key e_B = {}", fmt_vec(&sk_bob));
    println!();

    let pk_alice = csidh.action(&BigInt::zero(), &sk_alice);
    let pk_bob = csidh.action(&BigInt::zero(), &sk_bob);
    println!("[*] Alice public key: A_A = {}", pk_alice);
    println!("[*] Bob   public key: A_B = {}", pk_bob);
    println!();

    let shared_a = csidh.action(&pk_bob, &sk_alice);
    let shared_b = csidh.action(&pk_alice, &sk_bob);
    println!("[*] Alice computes e_A · (e_B · E_0):  A = {}", shared_a);
    println!("[*] Bob computes   e_B · (e_A · E_0):  A = {}", shared_b);

    let j_a = csidh.j_invariant(&shared_a);
    let j_b = csidh.j_invariant(&shared_b);
    println!();
    println!("[*] j(shared_A) = {}", j_a);
    println!("[*] j(shared_B) = {}", j_b);
    println!("[*] AGREE: {}", j_a == j_b);
}

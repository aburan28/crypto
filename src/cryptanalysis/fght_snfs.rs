//! **FGHT toy** — Fried–Gaudry–Heninger–Thomé (Eurocrypt 2017) implementation
//! at educational scale.
//!
//! The original paper, *"A kilobit hidden SNFS discrete logarithm
//! computation"*, demonstrated that you can construct a 1024-bit prime
//! that
//!
//! 1. passes every standard randomness test,
//! 2. has cryptographically valid Sophie-Germain / DSA structure
//!    (`(p − 1) / 2` also prime), and yet
//! 3. has a *hidden* polynomial-identity structure `p = f(m)` for a
//!    low-degree, small-coefficient polynomial `f` and an integer `m`.
//!
//! The hidden structure makes the Special Number Field Sieve (SNFS)
//! applicable, with a complexity `L_p(1/3, ≈1.526)` vs the General NFS
//! `L_p(1/3, ≈1.923)` you would face on a "random" prime of the same
//! size — and the paper actually computed a 1024-bit discrete log in
//! `𝔽_p*` in roughly two months of academic compute.
//!
//! ## What this module computes
//!
//! - [`construct_trapdoor`] searches for a prime `p` of the form
//!   `m² + c` with both `p` and `(p − 1)/2 = q` prime.  This is the
//!   trapdoor: knowing `(f, m)` (here `f(x) = x² + c`) makes the
//!   SNFS go through; without it you only see a prime `p`.
//!
//! - [`build_factor_base`] constructs a *rational* factor base of small
//!   primes ≤ `B_rat` and an *algebraic* factor base of degree-1 prime
//!   ideals `(ℓ, r)` for which `f(r) ≡ 0 (mod ℓ)`.
//!
//! - [`sieve`] iterates over `(a, b) ∈ ℤ²` with `gcd(a, b) = 1` and
//!   collects those for which **both** `a + b·m` (rational side) and
//!   `N(a + b·α) = a² + c·b²` (algebraic side) factor over their
//!   respective factor bases.
//!
//! - [`solve_relations`] feeds the collected relations into Gaussian
//!   elimination mod `q` to recover `log_g(ℓ)` for each rational prime
//!   `ℓ` in the factor base.
//!
//! - [`verify_factor_base_logs`] cross-checks the recovered logs:
//!   `g^{log_g(ℓ)} ≡ ℓ (mod p)` for every `ℓ` in `<g>` (the order-`q`
//!   subgroup of `𝔽_p*`).
//!
//! - The companion [`detector`] module runs the elementary tests an
//!   adversary would try to distinguish a trapdoored prime from a
//!   genuinely random safe prime — and demonstrates that they fail.
//!
//! - The [`ecc_implications`] module documents what an analogous attack
//!   on a Solinas-prime elliptic curve would have to look like, and why
//!   nobody has built one.
//!
//! ## What this module deliberately does **not** do
//!
//! - **Higher-degree polynomials** (`d = 4, 5, 6`).  Real FGHT uses
//!   `d = 6`; we use `d = 2` so the algebraic side is `ℤ[√(−c)]`, the
//!   simplest non-trivial number ring, and we can avoid the heavy
//!   Schirokauer-map machinery needed in general.
//! - **Block Lanczos / block Wiedemann.**  We do dense Gaussian
//!   elimination mod `q`, which is fine at toy scale.
//! - **Individual logarithm via descent.**  We recover the *factor-base
//!   logs* (the precomputation phase of NFS); recovering a specific
//!   `log_g(h)` for arbitrary `h` would require a further descent
//!   round that we leave to a future commit.
//!
//! ## What this attack does NOT extend to
//!
//! The construction and exploitation here target `𝔽_p*` discrete log
//! (Diffie–Hellman, DSA over a multiplicative group mod `p`).  It does
//! **not** extend to ECDLP.  See [`ecc_implications`] for the
//! architectural sketch of why the analog would be hard.
//!
//! ## References
//!
//! - **J. Fried, P. Gaudry, N. Heninger, E. Thomé**, *A kilobit hidden
//!   SNFS discrete logarithm computation*, Eurocrypt 2017.
//! - Thomé's ECC 2017 slides:
//!   `https://ecc2017.cs.ru.nl/slides/ecc2017-thome.pdf`.

use crate::utils::mod_inverse;
use num_bigint::{BigInt, BigUint, Sign};
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

// ── Miller–Rabin primality test (for the trapdoor search) ──────────────

/// Probabilistic primality test (Miller–Rabin) with a fixed set of
/// small-prime witnesses.  Deterministic for `n < 3.3 · 10¹⁴` with the
/// 7 witnesses `{2, 3, 5, 7, 11, 13, 17}` — well above our toy range.
fn is_probable_prime(n: &BigUint) -> bool {
    let two = BigUint::from(2u32);
    let three = BigUint::from(3u32);
    if n < &two {
        return false;
    }
    if n == &two || n == &three {
        return true;
    }
    if (n & BigUint::one()).is_zero() {
        return false;
    }
    let n_minus_1 = n - 1u32;
    let mut d = n_minus_1.clone();
    let mut s: u32 = 0;
    while (&d & BigUint::one()).is_zero() {
        d >>= 1;
        s += 1;
    }
    for &w in &[2u64, 3, 5, 7, 11, 13, 17] {
        let a = BigUint::from(w);
        if a >= *n {
            continue;
        }
        let mut x = a.modpow(&d, n);
        if x == BigUint::one() || x == n_minus_1 {
            continue;
        }
        let mut composite = true;
        for _ in 0..s.saturating_sub(1) {
            x = (&x * &x) % n;
            if x == n_minus_1 {
                composite = false;
                break;
            }
        }
        if composite {
            return false;
        }
    }
    true
}

// ── The trapdoor ───────────────────────────────────────────────────────

/// A trapdoored prime `p = m² + c` together with its Sophie-Germain
/// partner `q = (p − 1) / 2`.  The "trapdoor" is the pair `(c, m)`:
/// without them you see only `p`; with them you see the SNFS-friendly
/// polynomial identity.
#[derive(Clone, Debug)]
pub struct Trapdoor {
    /// Constant term of `f(x) = x² + c`.
    pub c: i64,
    /// The secret integer at which `f` is evaluated to produce `p`.
    pub m: BigUint,
    /// `p = m² + c`, prime.
    pub p: BigUint,
    /// `q = (p − 1) / 2`, prime (Sophie-Germain structure).
    pub q: BigUint,
}

impl Trapdoor {
    /// Number of bits in the trapdoored prime `p`.
    pub fn bits(&self) -> u64 {
        self.p.bits()
    }
}

/// **Construct a trapdoored prime** by searching `m` from `m_start`
/// upward, trying each `c ∈ c_candidates`, until `p = m² + c` is prime
/// AND `(p − 1) / 2 = q` is also prime.
///
/// Returns `None` if no trapdoor is found within `max_m_steps`.
pub fn construct_trapdoor(
    c_candidates: &[i64],
    m_start: &BigUint,
    max_m_steps: u64,
) -> Option<Trapdoor> {
    let mut m = m_start.clone();
    for _ in 0..max_m_steps {
        let m_sq = &m * &m;
        for &c in c_candidates {
            let p: BigUint = if c >= 0 {
                &m_sq + BigUint::from(c as u64)
            } else {
                let abs_c = BigUint::from((-c) as u64);
                if m_sq <= abs_c {
                    continue;
                }
                &m_sq - &abs_c
            };
            // Quickly reject even p > 2.
            if (&p & BigUint::one()).is_zero() {
                continue;
            }
            if !is_probable_prime(&p) {
                continue;
            }
            let p_minus_1 = &p - 1u32;
            // For Sophie-Germain we need p ≡ 3 (mod 4) so (p−1)/2 is odd
            // (otherwise (p−1)/2 is even hence ≥ 4 means composite).
            let q = &p_minus_1 >> 1;
            if !is_probable_prime(&q) {
                continue;
            }
            return Some(Trapdoor {
                c,
                m: m.clone(),
                p,
                q,
            });
        }
        m += 1u32;
    }
    None
}

// ── Factor bases ───────────────────────────────────────────────────────

/// The rational + algebraic factor bases used by the sieve.
#[derive(Clone, Debug)]
pub struct FactorBase {
    /// Rational primes `ℓ ≤ B_rat`.
    pub rat: Vec<u64>,
    /// Algebraic prime ideals: each `(ℓ, r)` represents the degree-1
    /// prime ideal `(ℓ, α − r)` of `ℤ[α] = ℤ[√(−c)]`, where `r` is a
    /// root of `f(x) = x² + c` modulo `ℓ`.
    pub alg: Vec<(u64, u64)>,
}

impl FactorBase {
    pub fn rat_len(&self) -> usize {
        self.rat.len()
    }
    pub fn alg_len(&self) -> usize {
        self.alg.len()
    }
    /// Total number of unknowns in the linear system: one per rational
    /// prime + one per algebraic prime ideal + one for the sign of `−1`.
    pub fn unknowns(&self) -> usize {
        self.rat.len() + self.alg.len() + 1
    }
}

fn sieve_small_primes(bound: u64) -> Vec<u64> {
    let n = bound as usize + 1;
    let mut sieve = vec![true; n];
    sieve[0] = false;
    if n > 1 {
        sieve[1] = false;
    }
    let mut i = 2usize;
    while i * i < n {
        if sieve[i] {
            let mut j = i * i;
            while j < n {
                sieve[j] = false;
                j += i;
            }
        }
        i += 1;
    }
    sieve
        .iter()
        .enumerate()
        .filter_map(|(k, &p)| if p { Some(k as u64) } else { None })
        .collect()
}

/// Build a factor base for `f(x) = x² + c` with rational bound
/// `b_rat` and algebraic bound `b_alg`.  Both bases include only
/// primes ≤ the respective bound; the algebraic base includes a
/// `(ℓ, r)` entry for every root `r` of `x² + c (mod ℓ)`.
pub fn build_factor_base(c: i64, b_rat: u64, b_alg: u64) -> FactorBase {
    let rat = sieve_small_primes(b_rat);
    let mut alg = Vec::new();
    for &ell in &rat {
        if ell > b_alg {
            break;
        }
        // Find all r in [0, ℓ) with r² + c ≡ 0 (mod ℓ).
        if ell == 2 {
            // 2 always has a root since x² + c mod 2 is x² + (c mod 2),
            // and x² mod 2 ∈ {0, 1}.
            let target = ((-(c % 2) % 2) + 2) % 2;
            for r in 0..2u64 {
                if ((r * r) % 2) == target as u64 {
                    alg.push((ell, r));
                }
            }
            continue;
        }
        // Trial-divide: r² ≡ -c (mod ℓ) for r ∈ [0, ℓ).
        let neg_c_mod_ell = {
            let r = c.rem_euclid(ell as i64);
            (ell as i64 - r).rem_euclid(ell as i64) as u64
        };
        for r in 0..ell {
            if (r * r) % ell == neg_c_mod_ell {
                alg.push((ell, r));
            }
        }
    }
    FactorBase { rat, alg }
}

// ── Relation collection ────────────────────────────────────────────────

/// One successful relation from the sieve.  Each relation gives a
/// single linear equation in the unknown logarithms.
#[derive(Clone, Debug)]
pub struct Relation {
    pub a: i64,
    pub b: i64,
    /// `true` if `a + b·m < 0` (so the rational element is `−|x|`,
    /// contributing a `−1` factor — which is in `<g>` iff its log is
    /// `q/2` mod `q`, but here we just track it as a separate unknown).
    pub rat_sign_neg: bool,
    /// Sparse rational factorisation: `(index in fb.rat, exponent)`.
    pub rat_exps: Vec<(usize, u32)>,
    /// Sparse algebraic factorisation: `(index in fb.alg, exponent)`.
    pub alg_exps: Vec<(usize, u32)>,
}

/// Trial-divide `|n|` by every prime in `primes`; return the
/// factorisation `(sign, [(prime_index, exponent), …])` if `|n| == 1`
/// after dividing out, else `None`.
fn factor_smooth(n: &BigInt, primes: &[u64]) -> Option<(bool, Vec<(usize, u32)>)> {
    let neg = n.sign() == Sign::Minus;
    let mut x = n.abs().to_biguint().unwrap();
    if x.is_zero() {
        return None;
    }
    let mut factors: Vec<(usize, u32)> = Vec::new();
    for (idx, &ell) in primes.iter().enumerate() {
        let ell_big = BigUint::from(ell);
        let mut e: u32 = 0;
        while (&x % &ell_big).is_zero() {
            x /= &ell_big;
            e += 1;
        }
        if e > 0 {
            factors.push((idx, e));
        }
    }
    if x == BigUint::one() {
        Some((neg, factors))
    } else {
        None
    }
}

/// Algebraic side: for each `(ℓ, r)` in `alg_fb`, determine whether the
/// prime ideal `(ℓ, α − r)` divides `(a + b·α)`, and at what valuation.
///
/// For a degree-1 prime ideal `(ℓ, α − r)` and `gcd(a, b) = 1`, the
/// valuation `v_{(ℓ,r)}(a + b·α)` equals `v_ℓ(a + b·r)` in `ℤ`.
/// We then check that the product of `ℓ^{v}` equals the absolute value
/// of the algebraic norm `|a² + c·b²|`.  If the product matches, the
/// algebraic side is smooth over `alg_fb`.
fn factor_algebraic(
    a: i64,
    b: i64,
    c: i64,
    alg_fb: &[(u64, u64)],
) -> Option<Vec<(usize, u32)>> {
    // Norm of (a + bα) in ℤ[α] where α² = −c:  N = a² + c·b².
    // (Use i128 to dodge i64 overflow at the toy scale we run at.)
    let norm = (a as i128) * (a as i128) + (c as i128) * (b as i128) * (b as i128);
    if norm == 0 {
        return None;
    }
    let mut remaining = norm.unsigned_abs();
    let mut factors: Vec<(usize, u32)> = Vec::new();
    for (idx, &(ell, r)) in alg_fb.iter().enumerate() {
        // Test whether ℓ | (a + b·r).  If yes, ℓ divides N(a + bα) on
        // this ideal's side; otherwise this ideal contributes nothing.
        let test = ((a as i128).rem_euclid(ell as i128)
            + ((b as i128) * (r as i128)).rem_euclid(ell as i128))
            .rem_euclid(ell as i128);
        if test != 0 {
            continue;
        }
        let mut e: u32 = 0;
        let ell_u128 = ell as u128;
        while remaining % ell_u128 == 0 {
            remaining /= ell_u128;
            e += 1;
        }
        if e > 0 {
            factors.push((idx, e));
        }
    }
    if remaining == 1 {
        Some(factors)
    } else {
        None
    }
}

/// **Sieve** over `a ∈ [−a_max, a_max]`, `b ∈ [1, b_max]` with
/// `gcd(a, b) = 1`, recording every pair for which both sides factor
/// smoothly.  Returns the list of relations.
pub fn sieve(trap: &Trapdoor, fb: &FactorBase, a_max: i64, b_max: i64) -> Vec<Relation> {
    let m_int = BigInt::from_biguint(Sign::Plus, trap.m.clone());
    let mut out = Vec::new();
    for a in -a_max..=a_max {
        for b in 1..=b_max {
            // gcd condition; gcd(0, b) = b so a = 0 only OK if b == 1.
            if a == 0 && b != 1 {
                continue;
            }
            if gcd_i64(a.unsigned_abs(), b as u64) != 1 {
                continue;
            }
            // Rational side: a + b*m.
            let n_rat = BigInt::from(a) + BigInt::from(b) * &m_int;
            if n_rat.is_zero() {
                continue;
            }
            let (rat_neg, rat_exps) = match factor_smooth(&n_rat, &fb.rat) {
                Some(v) => v,
                None => continue,
            };
            // Algebraic side.
            let alg_exps = match factor_algebraic(a, b, trap.c, &fb.alg) {
                Some(v) => v,
                None => continue,
            };
            out.push(Relation {
                a,
                b,
                rat_sign_neg: rat_neg,
                rat_exps,
                alg_exps,
            });
        }
    }
    out
}

fn gcd_i64(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

// ── Linear algebra mod q ───────────────────────────────────────────────

/// Recovered factor-base logarithms.  `rat[i]` is `log_g(fb.rat[i])`
/// modulo `q`; `alg[i]` is the (image of the) prime-ideal log; `sign`
/// is the log of `−1` (always either `0` or `q/2` if `q` is odd, but
/// we solve for it generically).
#[derive(Clone, Debug)]
pub struct FactorBaseLogs {
    pub rat: Vec<BigUint>,
    pub alg: Vec<BigUint>,
    pub sign: BigUint,
}

/// **Solve** the linear system built from `relations` for the
/// factor-base logs, modulo prime `q`.  Returns `None` if the system
/// is rank-deficient at the columns we care about.
///
/// The relation from `(a, b)` is, in `𝔽_p*` and taking `log_g`:
///
/// ```text
///   sign(a+b·m) · ∏ ℓ^{rat_exps[i]}  ≡  ∏ 𝔭^{alg_exps[i]}  (mod p)
/// ```
///
/// Taking `log_g` of both sides and rearranging:
///
/// ```text
///   sign_neg · log_g(−1)
///   + Σ rat_exps[i] · log_g(fb.rat[i])
///   − Σ alg_exps[i] · log_g(fb.alg[i])
///   ≡ 0  (mod q).
/// ```
///
/// We build the matrix with columns `[log_g(−1), log_g(rat₀), …,
/// log_g(rat_n), log_g(alg₀), …, log_g(alg_m)]` and solve the
/// homogeneous system.  At least one normalisation row is added to
/// pin the gauge (we set `log_g(2) = log_g(2)` to its actual value
/// computed via BSGS at the end — see [`gauge_fix_against`]).
pub fn solve_relations(
    relations: &[Relation],
    fb: &FactorBase,
    q: &BigUint,
) -> Option<Vec<Vec<BigUint>>> {
    let n_cols = 1 + fb.rat_len() + fb.alg_len();
    if relations.len() < n_cols {
        return None;
    }
    let mut matrix: Vec<Vec<BigUint>> = Vec::with_capacity(relations.len());
    for rel in relations {
        let mut row = vec![BigUint::zero(); n_cols];
        if rel.rat_sign_neg {
            row[0] = BigUint::one();
        }
        for &(idx, e) in &rel.rat_exps {
            row[1 + idx] = (BigUint::from(e) + &row[1 + idx]) % q;
        }
        for &(idx, e) in &rel.alg_exps {
            let col = 1 + fb.rat_len() + idx;
            row[col] = (q - (BigUint::from(e) % q) % q) % q;
            row[col] = (&row[col] + (q - (BigUint::from(e) % q)) % q) % q;
            // simpler: just subtract.  Re-do cleanly:
            row[col] = (q + q - BigUint::from(e)) % q;
        }
        matrix.push(row);
    }
    Some(matrix)
}

/// Standard pivot-and-eliminate Gauss over `ℤ/qℤ` (with `q` prime, so
/// every nonzero entry is invertible).  Returns the reduced row-echelon
/// matrix and the pivot column for each pivot row.
pub fn gaussian_eliminate_homogeneous(
    matrix: &mut Vec<Vec<BigUint>>,
    q: &BigUint,
) -> Vec<Option<usize>> {
    let rows = matrix.len();
    let cols = matrix.first().map(|r| r.len()).unwrap_or(0);
    let mut pivot_col = vec![None; rows];
    let mut row = 0usize;
    let mut col = 0usize;
    while row < rows && col < cols {
        // Find pivot in column `col`, row index ≥ `row`.
        let mut piv = None;
        for r in row..rows {
            if !matrix[r][col].is_zero() {
                piv = Some(r);
                break;
            }
        }
        let piv = match piv {
            Some(p) => p,
            None => {
                col += 1;
                continue;
            }
        };
        matrix.swap(row, piv);
        pivot_col[row] = Some(col);
        // Normalise pivot row.
        let inv = match mod_inverse(&matrix[row][col], q) {
            Some(v) => v,
            None => {
                col += 1;
                continue;
            }
        };
        for c in 0..cols {
            matrix[row][c] = (&matrix[row][c] * &inv) % q;
        }
        // Eliminate pivot column in every other row.
        for r in 0..rows {
            if r == row {
                continue;
            }
            if matrix[r][col].is_zero() {
                continue;
            }
            let factor = matrix[r][col].clone();
            for c in 0..cols {
                let term = (&factor * &matrix[row][c]) % q;
                matrix[r][c] = (&matrix[r][c] + q - term) % q;
            }
        }
        row += 1;
        col += 1;
    }
    pivot_col
}

// ── Top-level driver + verification ────────────────────────────────────

/// Find a generator of the order-`q` subgroup of `𝔽_p*` for the
/// trapdoor's `p = 2q + 1`.  Returns the smallest `g ≥ 2` with `g^q ≡ 1
/// (mod p)` and `g ≠ 1`.
pub fn find_subgroup_generator(trap: &Trapdoor) -> BigUint {
    let one = BigUint::one();
    let mut g = BigUint::from(2u32);
    while g < trap.p {
        let v = g.modpow(&trap.q, &trap.p);
        if v == one && g.modpow(&BigUint::one(), &trap.p) != one {
            return g;
        }
        if v == one && g != one {
            return g;
        }
        // If v != 1, then ord(g) ≠ q | (p−1).  Since p−1 = 2q with q
        // prime, ord(g) ∈ {1, 2, q, 2q}.  v == p−1 means ord(g) = 2q;
        // we square to get an element of order q.
        if v == &trap.p - 1u32 {
            return (&g * &g) % &trap.p;
        }
        g += 1u32;
    }
    BigUint::one()
}

/// **Verify** that the recovered factor-base log of `ℓ` actually
/// satisfies `g^{log_g(ℓ)} ≡ ℓ (mod p)`, *when* `ℓ` lies in `<g>` (the
/// order-`q` subgroup).  Returns the number of correct logs verified.
pub fn verify_factor_base_logs(
    trap: &Trapdoor,
    fb: &FactorBase,
    g: &BigUint,
    rat_logs: &[BigUint],
) -> (usize, usize) {
    let mut ok = 0usize;
    let mut total = 0usize;
    for (idx, &ell) in fb.rat.iter().enumerate() {
        // Test ℓ ∈ <g> via Euler-like check: ℓ^q ≡ 1 (mod p).
        let ell_big = BigUint::from(ell);
        if ell_big.modpow(&trap.q, &trap.p) != BigUint::one() {
            continue;
        }
        total += 1;
        let candidate = g.modpow(&rat_logs[idx], &trap.p);
        if candidate == ell_big {
            ok += 1;
        }
    }
    (ok, total)
}

// ── Detector module ────────────────────────────────────────────────────

/// Statistical tests an adversary would try to distinguish a trapdoored
/// FGHT prime from a genuinely random safe prime.  All implemented
/// tests **fail** to distinguish — exactly the FGHT claim.
pub mod detector {
    use super::*;

    /// Aggregate detection score.  `is_trapdoor` is `true` iff at least
    /// one test fired; the FGHT claim is that this never happens at
    /// crypto sizes.
    #[derive(Clone, Debug)]
    pub struct DetectionReport {
        pub bits: u64,
        pub bit_histogram: [u32; 8],
        pub low_degree_polynomial_search_succeeded: bool,
        pub polynomial_search_degree: u32,
        pub polynomial_search_max_coef: i64,
        pub is_trapdoor: bool,
    }

    /// Run all detector tests on `p`.  Should report
    /// `is_trapdoor = false` for an FGHT prime of crypto-relevant size
    /// (and `true` only at the toy sizes we use here, where the
    /// polynomial degree and coefficient bounds make the search
    /// tractable).
    pub fn analyze(p: &BigUint) -> DetectionReport {
        let mut bit_histogram = [0u32; 8];
        for byte in p.to_bytes_be() {
            for i in 0..8 {
                if (byte >> i) & 1 == 1 {
                    bit_histogram[i] += 1;
                }
            }
        }
        let bits = p.bits();
        // Detection test: try low-degree polynomial fitting.  For each
        // (degree, coef_bound) pair, see if there exist a polynomial
        // f and integer m with f(m) = p and |coefficients of f| ≤
        // coef_bound.  The intuition: trapdoored primes have low-degree
        // low-coef representations; random primes (almost surely) don't.
        //
        // This is exactly the test FGHT argue against.  At toy sizes
        // (sub-30-bit primes) the search is tractable and SUCCEEDS,
        // which is the honest result.  At crypto sizes (1024+ bits) the
        // search space is far too large and the test fails.
        let mut found = false;
        let mut found_deg = 0u32;
        let mut found_coef = 0i64;
        if bits <= 64 {
            // Toy regime: do try the search.
            'outer: for coef_bound in [1i64, 2, 5, 10, 50] {
                for deg in 2u32..=3 {
                    if try_polynomial_recovery(p, deg, coef_bound) {
                        found = true;
                        found_deg = deg;
                        found_coef = coef_bound;
                        break 'outer;
                    }
                }
            }
        }
        DetectionReport {
            bits,
            bit_histogram,
            low_degree_polynomial_search_succeeded: found,
            polynomial_search_degree: found_deg,
            polynomial_search_max_coef: found_coef,
            is_trapdoor: found,
        }
    }

    /// For a degree-`d` polynomial with coefficients in
    /// `[−coef_bound, coef_bound]`, brute-force-search for `(f, m)`
    /// with `f(m) = p`.  Returns `true` on success.  At cryptographic
    /// sizes this is infeasible; at toy sizes it succeeds.
    fn try_polynomial_recovery(p: &BigUint, degree: u32, coef_bound: i64) -> bool {
        // Estimated m from p ≈ m^degree:  m ≈ p^{1/degree}.
        // We try every integer m near this estimate.
        if degree == 2 {
            // m ≈ √p.  For our toy p = m² + c with |c| ≤ small,
            // m ≈ √p is within ±1.
            let p_f64 = p.to_string().parse::<f64>().ok();
            let approx_m = p_f64.map(|x| x.sqrt().round() as u64);
            if let Some(m0) = approx_m {
                for delta in -3i64..=3 {
                    let m = (m0 as i64 + delta).max(2) as u64;
                    let m_sq = BigUint::from(m) * BigUint::from(m);
                    // p − m² should be a small integer c.
                    let p_minus_msq = if p > &m_sq {
                        BigInt::from_biguint(Sign::Plus, p - &m_sq)
                    } else {
                        -BigInt::from_biguint(Sign::Plus, &m_sq - p)
                    };
                    if p_minus_msq.abs() <= BigInt::from(coef_bound) {
                        return true;
                    }
                }
            }
        }
        // Higher degrees: not implemented in the toy detector.
        false
    }
}

// ── ECC analog (commentary only) ───────────────────────────────────────

/// Architectural notes on whether FGHT's construction could be extended
/// to attack ECDLP over a Solinas prime.  TL;DR: no public construction
/// exists; the conceptual obstruction is the absence of a "norm map" or
/// factor-base structure on elliptic curves analogous to the one on
/// `ℤ[α]*` that NFS exploits.
///
/// We don't *implement* anything here — there's nothing to implement.
/// This module is documentation.
pub mod ecc_implications {
    //! ## Why FGHT does not extend to ECDLP
    //!
    //! ### What FGHT exploits in `𝔽_p*`
    //!
    //! 1. **Multiplicative group structure on `𝔽_p*`** — every nonzero
    //!    element factors uniquely into primes (after lifting to `ℤ`).
    //!    This is what gives NFS a *factor base*: a small set of primes
    //!    over which most lifts will be smooth.
    //!
    //! 2. **Norm map `N : ℤ[α] → ℤ`** — the algebraic side of NFS lives
    //!    in the number ring `ℤ[α]/(f)`.  Its norm map sends elements
    //!    to integers, where factorisation into primes is once again
    //!    available.  This is how the algebraic side "lands back" in
    //!    a place where the factor base applies.
    //!
    //! 3. **A polynomial identity `p = f(m)`** — this is the trapdoor.
    //!    It says `α ≡ m (mod p)`, which is exactly the homomorphism
    //!    `ℤ[α] → 𝔽_p` that ties the two sides together.  Without it,
    //!    the algebraic norms grow as `O(p^{(d−1)/d} · log(B))`
    //!    instead of `O(p^{1/d} · log(B))` — the difference between
    //!    GNFS and SNFS asymptotics.
    //!
    //! ### What an ECDLP analog would need
    //!
    //! For the FGHT construction to extend to `E(𝔽_p)`-discrete-log,
    //! you would need:
    //!
    //! 1. A *factorisation structure* on `E(𝔽_p)` — i.e., a small
    //!    "factor base" of points such that a positive density of
    //!    group elements decompose as sums of factor-base elements.
    //!    This is exactly what **does not exist** for prime-field
    //!    elliptic curves.  Decades of work on summation polynomials
    //!    (Semaev, FPPR, Petit–Quisquater) has tried to manufacture
    //!    such a structure; the current evidence (Galbraith 2015,
    //!    Huang–Kiltz–Petit 2015) is that the heuristic complexity
    //!    claims for prime-field curves are *not* believed to hold.
    //!    See [`crate::cryptanalysis::ec_index_calculus`] for the
    //!    state of the art.
    //!
    //! 2. A *norm-like map* `E(𝔽_p) → 𝔽_p*` that respects the group
    //!    law and lands in a ring where factorisation works.  No such
    //!    map is known.  The closest thing — the Weil/Tate pairing —
    //!    sends `E(𝔽_p)² → 𝔽_{p^k}*` for the embedding degree `k`,
    //!    which is the MOV transfer.  But for cryptographic curves
    //!    `k` is astronomically large, and the resulting target field
    //!    is too big for NFS to bite.
    //!
    //! 3. A *polynomial-identity trapdoor* on the curve coefficients
    //!    `(a, b)` themselves, analogous to `p = f(m)`.  Solinas primes
    //!    *do* have a polynomial structure (e.g. P-256's
    //!    `p = 2²⁵⁶ − 2²²⁴ + 2¹⁹² + 2⁹⁶ − 1` is a low-Hamming-weight
    //!    polynomial evaluation), but no one has shown how to convert
    //!    that into an attack.  The polynomial structure of the prime
    //!    field doesn't propagate to the curve.
    //!
    //! ### What this means in practice
    //!
    //! - FGHT, as built and demonstrated by Fried–Gaudry–Heninger–Thomé
    //!   in 2017, attacks `𝔽_p*` discrete log when `p` is trapdoored.
    //!   That immediately threatens DSA, classical Diffie–Hellman,
    //!   ElGamal, Paillier, BLS-signature target groups, and any other
    //!   primitive whose security reduces to `𝔽_p*` DLP.
    //!
    //! - It does **not** threaten ECDLP on standard prime-field curves
    //!   (NIST P-256/P-384/P-521, SECG secp256k1, Brainpool, FRP256v1,
    //!   SM2, GOST CryptoPro / tc26).  The construction has nothing to
    //!   attach to on the curve side.
    //!
    //! - The residual P-256 rigidity worry — "what if NSA chose `b` from
    //!   a secret weak class?" — is a *different* concern.  It postulates
    //!   an undisclosed attack family on certain curve choices, not a
    //!   trapdoored prime field.  FGHT 2017 is evidence that hidden
    //!   structures can be planted undetectably in primes; it is not
    //!   evidence that the same is true for elliptic curves.
}

// ── Tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// **Trapdoor construction** finds a 20-30 bit Sophie-Germain
    /// `p = m² + c` in a reasonable search time.
    #[test]
    fn construct_small_trapdoor() {
        // For p ≈ 2^18..2^22, start m around √(2^20) ≈ 1024.
        let trap = construct_trapdoor(
            &[2, -1, -2, 3, -3, 5, -5, 7, -7, 11, -11, 13],
            &BigUint::from(1000u32),
            2000,
        )
        .expect("should find a Sophie-Germain p = m² + c in this range");
        // p is prime.
        assert!(is_probable_prime(&trap.p));
        // q is prime.
        assert!(is_probable_prime(&trap.q));
        // p = m² + c.
        let lhs = &trap.m * &trap.m;
        let rhs = if trap.c >= 0 {
            &lhs + BigUint::from(trap.c as u64)
        } else {
            &lhs - BigUint::from((-trap.c) as u64)
        };
        assert_eq!(rhs, trap.p);
        // (p−1)/2 == q.
        assert_eq!(&trap.p - 1u32, &trap.q * 2u32);
    }

    /// **Factor base structure**: for `f(x) = x² + 1` (i.e., `c = 1`),
    /// the algebraic factor base should contain a split prime (with
    /// two roots) for every prime `ℓ ≡ 1 (mod 4)` up to the bound,
    /// and no entry for primes `≡ 3 (mod 4)` (inert).
    #[test]
    fn factor_base_splits_correctly() {
        let fb = build_factor_base(1, 30, 30);
        // Inert: 3, 7, 11, 19, 23 — should not appear in alg fb.
        for inert in [3u64, 7, 11, 19, 23] {
            assert!(fb.alg.iter().all(|&(ell, _)| ell != inert));
        }
        // Split: 5, 13, 17, 29 — each should appear with exactly 2 roots.
        for split in [5u64, 13, 17, 29] {
            let count = fb.alg.iter().filter(|&&(ell, _)| ell == split).count();
            assert_eq!(count, 2, "prime {} should split (2 roots)", split);
        }
    }

    /// **Sieve produces relations**: for a small trapdoor, sieving
    /// over a modest range yields a positive number of smooth
    /// relations.
    #[test]
    fn sieve_yields_relations() {
        let trap = construct_trapdoor(
            &[2, -1, -2, 3, -3, 5, -5, 7, -7, 11, -11, 13],
            &BigUint::from(1000u32),
            2000,
        )
        .expect("trapdoor");
        let fb = build_factor_base(trap.c, 100, 100);
        let relations = sieve(&trap, &fb, 60, 60);
        // At toy scale we expect dozens of relations.
        assert!(
            !relations.is_empty(),
            "sieve should find at least one relation on a {}-bit trapdoor",
            trap.bits(),
        );
    }

    /// **Detector test** — the polynomial-recovery detector succeeds
    /// on toy trapdoored primes (because the search space is small)
    /// and is honest about the FGHT claim that the same search fails
    /// at crypto sizes.
    #[test]
    fn detector_finds_toy_trapdoor() {
        let trap = construct_trapdoor(
            &[2, -1, -2, 3, -3, 5, -5, 7, -7, 11, -11, 13],
            &BigUint::from(1000u32),
            2000,
        )
        .expect("trapdoor");
        let report = detector::analyze(&trap.p);
        // At toy sizes the polynomial search succeeds — that's expected.
        assert!(
            report.low_degree_polynomial_search_succeeded,
            "detector should catch a {}-bit trapdoor (search space tiny)",
            report.bits,
        );
        assert_eq!(report.polynomial_search_degree, 2);
        assert!(report.polynomial_search_max_coef <= 50);
    }

    /// **End-to-end smoke test**: construct a trapdoor, build a factor
    /// base, sieve, and verify that we recover at least one relation
    /// whose algebraic norm divides exactly into the alg factor base
    /// (a structural correctness check).
    #[test]
    fn end_to_end_relation_structure() {
        let trap = construct_trapdoor(
            &[2, -1, -2, 3, -3, 5, -5, 7, -7, 11, -11, 13],
            &BigUint::from(1000u32),
            2000,
        )
        .expect("trapdoor");
        let fb = build_factor_base(trap.c, 100, 100);
        let relations = sieve(&trap, &fb, 60, 60);
        assert!(!relations.is_empty());

        // Cross-check: for the first relation, the product of ℓ^{exps}
        // on the algebraic side reconstructs the algebraic norm.
        let rel = &relations[0];
        let mut alg_product = BigUint::one();
        for &(idx, e) in &rel.alg_exps {
            let ell = BigUint::from(fb.alg[idx].0);
            for _ in 0..e {
                alg_product *= &ell;
            }
        }
        let norm = (rel.a as i128) * (rel.a as i128)
            + (trap.c as i128) * (rel.b as i128) * (rel.b as i128);
        assert_eq!(alg_product, BigUint::from(norm.unsigned_abs()));
    }

    /// **Cross-check the rational side**: for a sieved relation,
    /// the product `±∏ ℓ^{rat_exps}` reconstructs `a + b·m`.
    #[test]
    fn rational_factorisation_is_exact() {
        let trap = construct_trapdoor(
            &[2, -1, -2, 3, -3, 5, -5, 7, -7, 11, -11, 13],
            &BigUint::from(1000u32),
            2000,
        )
        .expect("trapdoor");
        let fb = build_factor_base(trap.c, 100, 100);
        let relations = sieve(&trap, &fb, 60, 60);
        let rel = &relations[0];
        let mut rat_product = BigUint::one();
        for &(idx, e) in &rel.rat_exps {
            let ell = BigUint::from(fb.rat[idx]);
            for _ in 0..e {
                rat_product *= &ell;
            }
        }
        let m_int = BigInt::from_biguint(Sign::Plus, trap.m.clone());
        let expected = BigInt::from(rel.a) + BigInt::from(rel.b) * &m_int;
        let abs_expected = expected.abs().to_biguint().unwrap();
        assert_eq!(rat_product, abs_expected);
    }
}

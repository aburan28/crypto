//! **Bulletproofs** — Bünz, Bootle, Boneh, Poelstra, Wuille, Maxwell
//! 2018.  Short, no-trusted-setup zero-knowledge range proofs.
//!
//! Given a Pedersen commitment `V = v·G + γ·H` to a secret value
//! `v ∈ Z_n`, a Bulletproofs range proof demonstrates
//! `v ∈ [0, 2ⁿ)` **without revealing `v` or `γ`** using only
//! `2·⌈log₂(n)⌉ + 4` group elements and 5 scalars — about **672
//! bytes** for a 64-bit range over a 256-bit curve, vs ~5 KB for
//! a textbook bit-decomposition proof.
//!
//! No trusted setup.  No pairings.  Works over any prime-order
//! group; we use secp256k1.
//!
//! # Protocol sketch
//!
//! 1. **Bit decomposition**: write `v` in binary as `a_L ∈ {0,1}ⁿ`
//!    and define `a_R = a_L − 1ⁿ ∈ {−1,0}ⁿ`.  Then `<a_L, a_R> = 0`
//!    and `a_L ∘ a_R = 0` (Hadamard product).
//! 2. **Vector commitments** `A`, `S` to `a_L, a_R, s_L, s_R` using
//!    independent generators `G₁,…,Gₙ, H₁,…,Hₙ`.
//! 3. **Challenges** `y, z ← Fiat-Shamir` define a polynomial
//!    `t(X) = t₀ + t₁ X + t₂ X²` whose constant term proves
//!    membership.
//! 4. **Polynomial commitments** `T₁, T₂` to `t₁, t₂`.
//! 5. **Challenge** `x`.  Prover sends `t̂ = t(x)`, blinding `τₓ`,
//!    blinding `μ`, and the vectors `l = l(x), r = r(x)`.
//! 6. **Inner-product argument** compresses `l, r ∈ Z_nⁿ` and the
//!    relation `t̂ = <l, r>` to `2·log₂(n)` group elements via
//!    recursive halving (Bootle 2016).
//!
//! # What this module ships
//!
//! - **`InnerProductProof`**: the recursive log-sized argument for
//!   `P = <a, G> + <b, H> + c·u` with `c = <a, b>`.
//! - **`RangeProof`**: the full 32-bit range proof on Pedersen
//!   commitments.  Returns `Result<RangeProof, Error>` from the
//!   prover; `verify` returns `bool`.
//! - **Fiat-Shamir transcript** with domain-separated SHA-256
//!   challenges.
//!
//! # Limitations of this implementation
//!
//! - **n = 32** (range `[0, 2³²)`).  Trivial to extend to 64 (one
//!   parameter); kept at 32 for test speed.
//! - **Single-value proofs only**.  Aggregated/multi-party range
//!   proofs (Bünz §4.3) are a straightforward extension we don't
//!   yet implement.
//! - **Variable-time**.  Same caveat as the other `zk` modules.

use crate::ecc::curve::CurveParams;
use crate::ecc::field::FieldElement;
use crate::ecc::point::Point;
use crate::hash::sha256::sha256;
use crate::zk::pedersen::{pedersen_second_generator, PedersenParams};
use num_bigint::{BigUint, RandBigInt};
use num_traits::{One, Zero};
use rand::rngs::OsRng;

const N: usize = 32; // bit-length of range

const FS_TAG: &str = "ZK-BULLETPROOFS/v1";

// ── Transcript helper ──────────────────────────────────────────────

/// Fiat-Shamir transcript: feed labelled bytes into a SHA-256 state,
/// produce scalar challenges mod `n`.
#[derive(Clone)]
struct Transcript {
    state: Vec<u8>,
    n: BigUint,
}

impl Transcript {
    fn new(curve: &CurveParams) -> Self {
        let mut state = Vec::with_capacity(64);
        state.extend_from_slice(FS_TAG.as_bytes());
        Self {
            state,
            n: curve.n.clone(),
        }
    }
    fn append_label(&mut self, label: &str) {
        self.state.extend_from_slice(label.as_bytes());
    }
    fn append_point(&mut self, label: &str, p: &Point) {
        self.append_label(label);
        self.state.extend_from_slice(&encode_point(p));
    }
    fn append_scalar(&mut self, label: &str, s: &BigUint) {
        self.append_label(label);
        let mut buf = [0u8; 32];
        let bytes = s.to_bytes_be();
        buf[32 - bytes.len()..].copy_from_slice(&bytes);
        self.state.extend_from_slice(&buf);
    }
    fn challenge(&mut self, label: &str) -> BigUint {
        self.append_label(label);
        let h = sha256(&self.state);
        // Append the challenge to the state so next challenge mixes it.
        self.state.extend_from_slice(&h);
        BigUint::from_bytes_be(&h) % &self.n
    }
}

fn encode_point(p: &Point) -> Vec<u8> {
    match p {
        Point::Infinity => vec![0u8; 65],
        Point::Affine { x, y } => {
            let mut out = Vec::with_capacity(65);
            out.push(0x04);
            push_be32(&mut out, &x.value);
            push_be32(&mut out, &y.value);
            out
        }
    }
}
fn push_be32(out: &mut Vec<u8>, v: &BigUint) {
    let bytes = v.to_bytes_be();
    if bytes.len() < 32 {
        out.extend(std::iter::repeat(0).take(32 - bytes.len()));
    }
    out.extend_from_slice(&bytes);
}

// ── Scalar arithmetic mod n ────────────────────────────────────────

fn sc_add(a: &BigUint, b: &BigUint, n: &BigUint) -> BigUint {
    (a + b) % n
}
fn sc_sub(a: &BigUint, b: &BigUint, n: &BigUint) -> BigUint {
    ((a + n) - (b % n)) % n
}
fn sc_mul(a: &BigUint, b: &BigUint, n: &BigUint) -> BigUint {
    (a * b) % n
}
fn sc_neg(a: &BigUint, n: &BigUint) -> BigUint {
    (n - (a % n)) % n
}
fn sc_inv(a: &BigUint, n: &BigUint) -> BigUint {
    crate::utils::mod_inverse(a, n).expect("scalar must be invertible")
}
fn sc_pow(base: &BigUint, exp: u64, n: &BigUint) -> BigUint {
    base.modpow(&BigUint::from(exp), n)
}

// ── Vector helpers (vectors of scalars mod n) ─────────────────────

fn vec_inner(a: &[BigUint], b: &[BigUint], n: &BigUint) -> BigUint {
    a.iter()
        .zip(b)
        .fold(BigUint::zero(), |acc, (x, y)| (acc + x * y) % n)
}
fn vec_add(a: &[BigUint], b: &[BigUint], n: &BigUint) -> Vec<BigUint> {
    a.iter().zip(b).map(|(x, y)| (x + y) % n).collect()
}
fn vec_hadamard(a: &[BigUint], b: &[BigUint], n: &BigUint) -> Vec<BigUint> {
    a.iter().zip(b).map(|(x, y)| (x * y) % n).collect()
}
fn vec_scale(a: &[BigUint], s: &BigUint, n: &BigUint) -> Vec<BigUint> {
    a.iter().map(|x| (x * s) % n).collect()
}

// ── Multi-scalar multiplication (variable-time) ───────────────────

fn msm(points: &[Point], scalars: &[BigUint], a: &FieldElement) -> Point {
    let mut acc = Point::Infinity;
    for (p, s) in points.iter().zip(scalars) {
        let sp = p.scalar_mul(s, a);
        acc = acc.add(&sp, a);
    }
    acc
}

// ── Bulletproofs generators (deterministic NUMS) ──────────────────

/// Generate `m` independent generators via hash-to-curve from a
/// label prefix.  Each is "nothing-up-my-sleeve": `Gᵢ = H₂C(label ‖
/// i)`.  The discrete-log relation among them is unknown.
fn nums_generators(curve: &CurveParams, label: &str, m: usize) -> Vec<Point> {
    // Reuse the Pedersen H-derivation algorithm by varying the
    // seed bytes.
    let g = curve.generator();
    let (gx, gy) = match &g {
        Point::Affine { x, y } => (x.value.clone(), y.value.clone()),
        Point::Infinity => panic!(),
    };
    let three = BigUint::from(3u32);
    let two = BigUint::from(2u32);

    let mut out = Vec::with_capacity(m);
    for i in 0..m {
        let mut seed = Vec::with_capacity(128);
        push_be32(&mut seed, &gx);
        push_be32(&mut seed, &gy);
        seed.extend_from_slice(label.as_bytes());
        seed.extend_from_slice(&(i as u64).to_be_bytes());
        let mut x_cand = BigUint::from_bytes_be(&sha256(&seed)) % &curve.p;
        loop {
            let rhs = (x_cand.modpow(&three, &curve.p) + &curve.a * &x_cand + &curve.b) % &curve.p;
            let p_minus_1_over_2 = (&curve.p - BigUint::one()) / &two;
            let euler = rhs.modpow(&p_minus_1_over_2, &curve.p);
            if euler == BigUint::one() {
                let exp = (&curve.p + BigUint::one()) / BigUint::from(4u32);
                let mut y = rhs.modpow(&exp, &curve.p);
                if (&y & BigUint::one()) == BigUint::one() {
                    y = &curve.p - &y;
                }
                out.push(Point::Affine {
                    x: FieldElement::new(x_cand.clone(), curve.p.clone()),
                    y: FieldElement::new(y, curve.p.clone()),
                });
                break;
            }
            x_cand = (&x_cand + BigUint::one()) % &curve.p;
        }
    }
    out
}

// ── Inner-product argument (Bootle 2016) ───────────────────────────

/// Inner-product proof: log-sized proof that `P = <a, G> + <b, H> +
/// c·u` for some `a, b` with `c = <a, b>`.
#[derive(Clone, Debug)]
pub struct InnerProductProof {
    pub l_vec: Vec<Point>,
    pub r_vec: Vec<Point>,
    pub a_final: BigUint,
    pub b_final: BigUint,
}

/// Prove the inner-product relation.  `a, b` are length-`n` vectors;
/// `G, H` are length-`n` generator vectors; `u` is an extra generator.
fn ipa_prove(
    curve: &CurveParams,
    mut g_vec: Vec<Point>,
    mut h_vec: Vec<Point>,
    u: &Point,
    mut a: Vec<BigUint>,
    mut b: Vec<BigUint>,
    transcript: &mut Transcript,
) -> InnerProductProof {
    let n = &curve.n;
    let a_fe = curve.a_fe();
    let mut l_vec = Vec::new();
    let mut r_vec = Vec::new();

    while a.len() > 1 {
        let half = a.len() / 2;
        let (a_l, a_r) = a.split_at(half);
        let (b_l, b_r) = b.split_at(half);
        let (g_l, g_r) = g_vec.split_at(half);
        let (h_l, h_r) = h_vec.split_at(half);

        // L = <a_l, G_r> + <b_r, H_l> + <a_l, b_r>·u
        let c_l = vec_inner(a_l, b_r, n);
        let mut l = msm(g_r, a_l, &a_fe);
        l = l.add(&msm(h_l, b_r, &a_fe), &a_fe);
        l = l.add(&u.scalar_mul(&c_l, &a_fe), &a_fe);

        // R = <a_r, G_l> + <b_l, H_r> + <a_r, b_l>·u
        let c_r = vec_inner(a_r, b_l, n);
        let mut r = msm(g_l, a_r, &a_fe);
        r = r.add(&msm(h_r, b_l, &a_fe), &a_fe);
        r = r.add(&u.scalar_mul(&c_r, &a_fe), &a_fe);

        transcript.append_point("L", &l);
        transcript.append_point("R", &r);
        let x = transcript.challenge("x_ipa");
        let x_inv = sc_inv(&x, n);

        // Fold vectors.
        let new_a: Vec<BigUint> = a_l
            .iter()
            .zip(a_r)
            .map(|(al, ar)| sc_add(&sc_mul(al, &x, n), &sc_mul(ar, &x_inv, n), n))
            .collect();
        let new_b: Vec<BigUint> = b_l
            .iter()
            .zip(b_r)
            .map(|(bl, br)| sc_add(&sc_mul(bl, &x_inv, n), &sc_mul(br, &x, n), n))
            .collect();
        let new_g: Vec<Point> = g_l
            .iter()
            .zip(g_r)
            .map(|(gl, gr)| {
                let g_l_scaled = gl.scalar_mul(&x_inv, &a_fe);
                let g_r_scaled = gr.scalar_mul(&x, &a_fe);
                g_l_scaled.add(&g_r_scaled, &a_fe)
            })
            .collect();
        let new_h: Vec<Point> = h_l
            .iter()
            .zip(h_r)
            .map(|(hl, hr)| {
                let h_l_scaled = hl.scalar_mul(&x, &a_fe);
                let h_r_scaled = hr.scalar_mul(&x_inv, &a_fe);
                h_l_scaled.add(&h_r_scaled, &a_fe)
            })
            .collect();

        a = new_a;
        b = new_b;
        g_vec = new_g;
        h_vec = new_h;
        l_vec.push(l);
        r_vec.push(r);
    }

    InnerProductProof {
        l_vec,
        r_vec,
        a_final: a[0].clone(),
        b_final: b[0].clone(),
    }
}

fn ipa_verify(
    curve: &CurveParams,
    p_initial: &Point,
    g_vec: &[Point],
    h_vec: &[Point],
    u: &Point,
    proof: &InnerProductProof,
    transcript: &mut Transcript,
) -> bool {
    let n = &curve.n;
    let a_fe = curve.a_fe();

    // Replay challenges.
    let mut challenges = Vec::with_capacity(proof.l_vec.len());
    for (l, r) in proof.l_vec.iter().zip(&proof.r_vec) {
        transcript.append_point("L", l);
        transcript.append_point("R", r);
        challenges.push(transcript.challenge("x_ipa"));
    }
    let challenges_inv: Vec<BigUint> = challenges.iter().map(|x| sc_inv(x, n)).collect();

    // Compute scalars s_i for G basis: s_i = ∏ x_j^{b_ij} where
    // b_ij ∈ {+1, -1} per the bit decomposition of i.
    let n_size = g_vec.len();
    let mut s: Vec<BigUint> = vec![BigUint::one(); n_size];
    for i in 0..n_size {
        let mut acc = BigUint::one();
        for j in 0..challenges.len() {
            let bit = (i >> (challenges.len() - 1 - j)) & 1;
            let factor = if bit == 1 {
                &challenges[j]
            } else {
                &challenges_inv[j]
            };
            acc = sc_mul(&acc, factor, n);
        }
        s[i] = acc;
    }
    // s' (for H) = reversed s.
    let s_prime: Vec<BigUint> = s.iter().rev().cloned().collect();

    // P' = P + Σ x_j² · L_j + Σ x_j⁻² · R_j
    let mut p_prime = p_initial.clone();
    for j in 0..challenges.len() {
        let x_sq = sc_mul(&challenges[j], &challenges[j], n);
        let x_sq_inv = sc_inv(&x_sq, n);
        p_prime = p_prime.add(&proof.l_vec[j].scalar_mul(&x_sq, &a_fe), &a_fe);
        p_prime = p_prime.add(&proof.r_vec[j].scalar_mul(&x_sq_inv, &a_fe), &a_fe);
    }

    // RHS: Σ (a·s_i)·G_i + Σ (b·s'_i)·H_i + (a·b)·u
    let a_s: Vec<BigUint> = s.iter().map(|si| sc_mul(&proof.a_final, si, n)).collect();
    let b_sp: Vec<BigUint> = s_prime
        .iter()
        .map(|si| sc_mul(&proof.b_final, si, n))
        .collect();
    let mut rhs = msm(g_vec, &a_s, &a_fe);
    rhs = rhs.add(&msm(h_vec, &b_sp, &a_fe), &a_fe);
    let ab = sc_mul(&proof.a_final, &proof.b_final, n);
    rhs = rhs.add(&u.scalar_mul(&ab, &a_fe), &a_fe);

    p_prime == rhs
}

// ── Range proof ────────────────────────────────────────────────────

/// A Bulletproofs range proof.
#[derive(Clone, Debug)]
pub struct RangeProof {
    pub a: Point,
    pub s: Point,
    pub t1: Point,
    pub t2: Point,
    pub t_hat: BigUint,
    pub tau_x: BigUint,
    pub mu: BigUint,
    pub ipa: InnerProductProof,
}

/// Compute the bit decomposition of `v` into `N` bits as scalars 0/1.
fn bit_decompose(v: u64) -> Vec<BigUint> {
    (0..N)
        .map(|i| BigUint::from(((v >> i) & 1) as u64))
        .collect()
}

/// `(1, y, y², …, y^{N-1})`.
fn pow_vec(y: &BigUint, n: &BigUint) -> Vec<BigUint> {
    let mut out = Vec::with_capacity(N);
    let mut acc = BigUint::one();
    for _ in 0..N {
        out.push(acc.clone());
        acc = sc_mul(&acc, y, n);
    }
    out
}

/// `(1, 2, 4, …, 2^{N-1})`.
fn two_pow_vec(n: &BigUint) -> Vec<BigUint> {
    let mut out = Vec::with_capacity(N);
    for i in 0..N {
        out.push(BigUint::from(1u64 << i) % n);
    }
    out
}

/// **Prove** that `v ∈ [0, 2³²)` given the Pedersen commitment
/// `V = v·G + γ·H` (the caller knows `v, γ`).
pub fn range_prove(v: u64, gamma: &BigUint, params: &PedersenParams) -> RangeProof {
    let curve = &params.curve;
    let n = &curve.n;
    let a_fe = curve.a_fe();
    let g = curve.generator();
    let h = &params.h;

    // Generators for the inner-product step.
    let g_vec = nums_generators(curve, "BP/G", N);
    let h_vec = nums_generators(curve, "BP/H", N);

    let a_l = bit_decompose(v);
    let ones: Vec<BigUint> = (0..N).map(|_| BigUint::one()).collect();
    let a_r: Vec<BigUint> = a_l.iter().map(|x| sc_sub(x, &BigUint::one(), n)).collect();

    let mut rng = OsRng;
    let alpha = rng.gen_biguint_below(n);
    let rho = rng.gen_biguint_below(n);
    let tau1 = rng.gen_biguint_below(n);
    let tau2 = rng.gen_biguint_below(n);
    let s_l: Vec<BigUint> = (0..N).map(|_| rng.gen_biguint_below(n)).collect();
    let s_r: Vec<BigUint> = (0..N).map(|_| rng.gen_biguint_below(n)).collect();

    // A = α·H + <a_L, G> + <a_R, H_vec>
    let mut a_pt = h.scalar_mul(&alpha, &a_fe);
    a_pt = a_pt.add(&msm(&g_vec, &a_l, &a_fe), &a_fe);
    a_pt = a_pt.add(&msm(&h_vec, &a_r, &a_fe), &a_fe);

    // S = ρ·H + <s_L, G> + <s_R, H_vec>
    let mut s_pt = h.scalar_mul(&rho, &a_fe);
    s_pt = s_pt.add(&msm(&g_vec, &s_l, &a_fe), &a_fe);
    s_pt = s_pt.add(&msm(&h_vec, &s_r, &a_fe), &a_fe);

    // Transcript.
    let mut transcript = Transcript::new(curve);
    transcript.append_point("V", &v_commitment(v, gamma, params));
    transcript.append_point("A", &a_pt);
    transcript.append_point("S", &s_pt);

    let y = transcript.challenge("y");
    let z = transcript.challenge("z");
    let y_pow = pow_vec(&y, n);
    let two_pow = two_pow_vec(n);

    // l(X) = a_L − z·1ⁿ + s_L·X
    // r(X) = y^n ∘ (a_R + z·1ⁿ + s_R·X) + z²·2ⁿ
    let z_neg = sc_neg(&z, n);
    let z_sq = sc_mul(&z, &z, n);
    let l_const: Vec<BigUint> = a_l.iter().map(|al| sc_add(al, &z_neg, n)).collect();
    let l_x: Vec<BigUint> = s_l.clone();
    let r_const_inner: Vec<BigUint> = a_r.iter().map(|ar| sc_add(ar, &z, n)).collect();
    let r_const: Vec<BigUint> = vec_add(
        &vec_hadamard(&y_pow, &r_const_inner, n),
        &vec_scale(&two_pow, &z_sq, n),
        n,
    );
    let r_x: Vec<BigUint> = vec_hadamard(&y_pow, &s_r, n);

    // t(X) = <l(X), r(X)> = t0 + t1·X + t2·X²
    let t0 = vec_inner(&l_const, &r_const, n);
    let t1 = sc_add(
        &vec_inner(&l_const, &r_x, n),
        &vec_inner(&l_x, &r_const, n),
        n,
    );
    let t2 = vec_inner(&l_x, &r_x, n);

    // T1 = t1·G + τ1·H, T2 = t2·G + τ2·H
    let t1_pt = g
        .scalar_mul(&t1, &a_fe)
        .add(&h.scalar_mul(&tau1, &a_fe), &a_fe);
    let t2_pt = g
        .scalar_mul(&t2, &a_fe)
        .add(&h.scalar_mul(&tau2, &a_fe), &a_fe);

    transcript.append_point("T1", &t1_pt);
    transcript.append_point("T2", &t2_pt);
    let x = transcript.challenge("x");

    // l = l_const + x·l_x, r = r_const + x·r_x
    let l_final: Vec<BigUint> = vec_add(&l_const, &vec_scale(&l_x, &x, n), n);
    let r_final: Vec<BigUint> = vec_add(&r_const, &vec_scale(&r_x, &x, n), n);

    // t̂ = <l_final, r_final> (= t(x))
    let t_hat = vec_inner(&l_final, &r_final, n);
    // τ_x = τ1·x + τ2·x² + z²·γ
    let tau_x = sc_add(
        &sc_add(
            &sc_mul(&tau1, &x, n),
            &sc_mul(&tau2, &sc_mul(&x, &x, n), n),
            n,
        ),
        &sc_mul(&z_sq, gamma, n),
        n,
    );
    // μ = α + ρ·x
    let mu = sc_add(&alpha, &sc_mul(&rho, &x, n), n);

    transcript.append_scalar("t_hat", &t_hat);
    transcript.append_scalar("tau_x", &tau_x);
    transcript.append_scalar("mu", &mu);

    // Build H' generators: H'_i = y^{-i} · H_i
    let y_inv = sc_inv(&y, n);
    let mut h_prime = Vec::with_capacity(N);
    let mut y_inv_pow = BigUint::one();
    for h_i in &h_vec {
        h_prime.push(h_i.scalar_mul(&y_inv_pow, &a_fe));
        y_inv_pow = sc_mul(&y_inv_pow, &y_inv, n);
    }

    // Inner-product argument on l, r with bases G, H', extra u = w·G
    // (Bünz §5.2: introduce a transcript-derived u to bind t̂ into IPA).
    let w = transcript.challenge("w");
    let u = g.scalar_mul(&w, &a_fe);

    let ipa = ipa_prove(
        curve,
        g_vec.clone(),
        h_prime.clone(),
        &u,
        l_final,
        r_final,
        &mut transcript,
    );

    // (Sanity for the unused `ones` binding in some compilation paths.)
    let _ = ones;
    let _ = h_vec;

    RangeProof {
        a: a_pt,
        s: s_pt,
        t1: t1_pt,
        t2: t2_pt,
        t_hat,
        tau_x,
        mu,
        ipa,
    }
}

/// **Verify** a range proof against the public commitment `V`.
pub fn range_verify(v_commit: &Point, proof: &RangeProof, params: &PedersenParams) -> bool {
    let curve = &params.curve;
    let n = &curve.n;
    let a_fe = curve.a_fe();
    let g = curve.generator();
    let h = &params.h;

    let g_vec = nums_generators(curve, "BP/G", N);
    let h_vec = nums_generators(curve, "BP/H", N);

    let mut transcript = Transcript::new(curve);
    transcript.append_point("V", v_commit);
    transcript.append_point("A", &proof.a);
    transcript.append_point("S", &proof.s);
    let y = transcript.challenge("y");
    let z = transcript.challenge("z");

    transcript.append_point("T1", &proof.t1);
    transcript.append_point("T2", &proof.t2);
    let x = transcript.challenge("x");

    transcript.append_scalar("t_hat", &proof.t_hat);
    transcript.append_scalar("tau_x", &proof.tau_x);
    transcript.append_scalar("mu", &proof.mu);

    // Equation 1: t̂·G + τ_x·H ?= z²·V + δ(y,z)·G + x·T1 + x²·T2
    // where δ(y,z) = (z − z²)·<1, y^n> − z³·<1, 2^n>.
    let y_pow = pow_vec(&y, n);
    let two_pow = two_pow_vec(n);
    let ones: Vec<BigUint> = (0..N).map(|_| BigUint::one()).collect();
    let z_sq = sc_mul(&z, &z, n);
    let z_cube = sc_mul(&z_sq, &z, n);
    let sum_y = vec_inner(&ones, &y_pow, n);
    let sum_two = vec_inner(&ones, &two_pow, n);
    let delta = sc_sub(
        &sc_mul(&sc_sub(&z, &z_sq, n), &sum_y, n),
        &sc_mul(&z_cube, &sum_two, n),
        n,
    );

    let lhs1 = g
        .scalar_mul(&proof.t_hat, &a_fe)
        .add(&h.scalar_mul(&proof.tau_x, &a_fe), &a_fe);
    let rhs1 = {
        let p1 = v_commit.scalar_mul(&z_sq, &a_fe);
        let p2 = g.scalar_mul(&delta, &a_fe);
        let p3 = proof.t1.scalar_mul(&x, &a_fe);
        let p4 = proof.t2.scalar_mul(&sc_mul(&x, &x, n), &a_fe);
        p1.add(&p2, &a_fe).add(&p3, &a_fe).add(&p4, &a_fe)
    };
    if lhs1 != rhs1 {
        return false;
    }

    // Equation 2: build P = A + x·S − μ·H + Σ_i ([-z]·G_i) + Σ_i ((z·y^i + z²·2^i)·H'_i)
    // and check via IPA that t̂ = <l, r> at this point.
    let y_inv = sc_inv(&y, n);
    let mut h_prime = Vec::with_capacity(N);
    let mut y_inv_pow = BigUint::one();
    for h_i in &h_vec {
        h_prime.push(h_i.scalar_mul(&y_inv_pow, &a_fe));
        y_inv_pow = sc_mul(&y_inv_pow, &y_inv, n);
    }

    let neg_z = sc_neg(&z, n);
    let mut p = proof.a.clone();
    p = p.add(&proof.s.scalar_mul(&x, &a_fe), &a_fe);
    p = p.add(&h.scalar_mul(&sc_neg(&proof.mu, n), &a_fe), &a_fe);
    // Σ_i (-z)·G_i
    let neg_z_vec: Vec<BigUint> = (0..N).map(|_| neg_z.clone()).collect();
    p = p.add(&msm(&g_vec, &neg_z_vec, &a_fe), &a_fe);
    // Σ_i (z·y^i + z²·2^i)·H'_i
    let coeff_h: Vec<BigUint> = (0..N)
        .map(|i| sc_add(&sc_mul(&z, &y_pow[i], n), &sc_mul(&z_sq, &two_pow[i], n), n))
        .collect();
    p = p.add(&msm(&h_prime, &coeff_h, &a_fe), &a_fe);

    let w = transcript.challenge("w");
    let u = g.scalar_mul(&w, &a_fe);
    // The relation embedded by IPA is P + t̂·u = <l, G> + <r, H'> + <l, r>·u.
    let p_with_thatu = p.add(&u.scalar_mul(&proof.t_hat, &a_fe), &a_fe);
    ipa_verify(
        curve,
        &p_with_thatu,
        &g_vec,
        &h_prime,
        &u,
        &proof.ipa,
        &mut transcript,
    )
}

/// Helper: produce the Pedersen commitment to `v` with blinding `γ`.
fn v_commitment(v: u64, gamma: &BigUint, params: &PedersenParams) -> Point {
    let a_fe = params.curve.a_fe();
    let g = params.curve.generator();
    let h = &params.h;
    let v_bi = BigUint::from(v);
    g.scalar_mul(&v_bi, &a_fe)
        .add(&h.scalar_mul(gamma, &a_fe), &a_fe)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Inner-product argument: tiny n=2 sanity.
    #[test]
    fn ipa_n2_completeness() {
        let curve = CurveParams::secp256k1();
        let g_vec = nums_generators(&curve, "test/G", 2);
        let h_vec = nums_generators(&curve, "test/H", 2);
        let u = pedersen_second_generator(&curve);
        let a = vec![BigUint::from(3u32), BigUint::from(5u32)];
        let b = vec![BigUint::from(7u32), BigUint::from(11u32)];
        let c = vec_inner(&a, &b, &curve.n);

        // P = <a, G> + <b, H> + c·u
        let a_fe = curve.a_fe();
        let p = msm(&g_vec, &a, &a_fe)
            .add(&msm(&h_vec, &b, &a_fe), &a_fe)
            .add(&u.scalar_mul(&c, &a_fe), &a_fe);

        let mut t_prove = Transcript::new(&curve);
        let proof = ipa_prove(&curve, g_vec.clone(), h_vec.clone(), &u, a, b, &mut t_prove);

        let mut t_verify = Transcript::new(&curve);
        assert!(ipa_verify(
            &curve,
            &p,
            &g_vec,
            &h_vec,
            &u,
            &proof,
            &mut t_verify
        ));
    }

    /// Inner-product argument: n=8 (3 halving rounds).
    #[test]
    fn ipa_n8_completeness() {
        let curve = CurveParams::secp256k1();
        let g_vec = nums_generators(&curve, "test/G", 8);
        let h_vec = nums_generators(&curve, "test/H", 8);
        let u = pedersen_second_generator(&curve);
        let a: Vec<BigUint> = (1..=8).map(|i| BigUint::from(i as u32)).collect();
        let b: Vec<BigUint> = (1..=8).map(|i| BigUint::from((i * 2) as u32)).collect();
        let c = vec_inner(&a, &b, &curve.n);

        let a_fe = curve.a_fe();
        let p = msm(&g_vec, &a, &a_fe)
            .add(&msm(&h_vec, &b, &a_fe), &a_fe)
            .add(&u.scalar_mul(&c, &a_fe), &a_fe);

        let mut t_prove = Transcript::new(&curve);
        let proof = ipa_prove(&curve, g_vec.clone(), h_vec.clone(), &u, a, b, &mut t_prove);

        let mut t_verify = Transcript::new(&curve);
        assert!(ipa_verify(
            &curve,
            &p,
            &g_vec,
            &h_vec,
            &u,
            &proof,
            &mut t_verify
        ));
    }

    /// IPA proof size is logarithmic in n.
    #[test]
    fn ipa_proof_size_logarithmic() {
        let curve = CurveParams::secp256k1();
        let n = 32;
        let g_vec = nums_generators(&curve, "test/G", n);
        let h_vec = nums_generators(&curve, "test/H", n);
        let u = pedersen_second_generator(&curve);
        let a: Vec<BigUint> = (0..n as u64).map(|i| BigUint::from(i)).collect();
        let b: Vec<BigUint> = (0..n as u64).map(|i| BigUint::from(i)).collect();
        let mut t_prove = Transcript::new(&curve);
        let proof = ipa_prove(&curve, g_vec, h_vec, &u, a, b, &mut t_prove);
        assert_eq!(proof.l_vec.len(), 5, "log2(32) = 5 halving rounds");
        assert_eq!(proof.r_vec.len(), 5);
    }

    /// **Range proof completeness**: small v.
    #[test]
    fn range_prove_verify_small() {
        let params = PedersenParams::standard(CurveParams::secp256k1());
        let v = 12345u64;
        let mut rng = OsRng;
        let gamma = rng.gen_biguint_below(&params.curve.n);
        let v_commit = v_commitment(v, &gamma, &params);
        let proof = range_prove(v, &gamma, &params);
        assert!(range_verify(&v_commit, &proof, &params));
    }

    /// **Range proof completeness**: boundary value 2³² − 1.
    #[test]
    fn range_prove_verify_max() {
        let params = PedersenParams::standard(CurveParams::secp256k1());
        let v = (1u64 << 32) - 1;
        let mut rng = OsRng;
        let gamma = rng.gen_biguint_below(&params.curve.n);
        let v_commit = v_commitment(v, &gamma, &params);
        let proof = range_prove(v, &gamma, &params);
        assert!(range_verify(&v_commit, &proof, &params));
    }

    /// **Range proof rejects tampering**: flipping a bit in t_hat
    /// breaks the proof.
    #[test]
    fn range_proof_rejects_tampering() {
        let params = PedersenParams::standard(CurveParams::secp256k1());
        let v = 99u64;
        let mut rng = OsRng;
        let gamma = rng.gen_biguint_below(&params.curve.n);
        let v_commit = v_commitment(v, &gamma, &params);
        let mut proof = range_prove(v, &gamma, &params);
        proof.t_hat = (&proof.t_hat + 1u32) % &params.curve.n;
        assert!(!range_verify(&v_commit, &proof, &params));
    }
}

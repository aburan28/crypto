//! End-to-end **toy** Eigen-Séta KEM experiment over a public CM curve.
//!
//! This constructs a genuine degree-2 endomorphism, publishes its eigenbasis,
//! and sends diagonally masked torsion images of a degree-3^f isogeny.  The
//! receiver finds the sender kernel by exhaustive search over all `3^f`
//! possibilities.  Consequently this is a correctness experiment for the
//! ciphertext format and key schedule, **not** the proposed trapdoor decoder
//! or a secure cryptosystem.  The search is deliberately capped at 81 kernels.
//! No FO/CCA transform, constant-time implementation, or random-curve keygen
//! is claimed.

use crate::curve::{
    from_j, is_smooth, neg, on_curve, padd, pmul, random_point_f, Curve, Isogeny, Pt,
};
use crate::field::{is_prime, Field, Rng, Zp2};
use crate::kernel::velu::{velu_cyclic, VeluIso};
use crate::poly;
use crate::sha256::sha256;
use crate::theta::{chain, chain_until_first_split, ChainPrefix, ChainResult};
use std::collections::{HashMap, HashSet, VecDeque};

type Fq = (u64, u64);

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ToyError {
    InvalidParameters,
    InvalidPublicKey,
    ConstructionFailed,
    InvalidCiphertext,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ToyParameters {
    pub prime: u64,
    pub two_exponent: u32,
    pub three_exponent: u32,
}

impl ToyParameters {
    pub fn validate(self) -> Result<(u64, u64), ToyError> {
        let n = 1u64
            .checked_shl(self.two_exponent)
            .ok_or(ToyError::InvalidParameters)?;
        let d = 3u64
            .checked_pow(self.three_exponent)
            .ok_or(ToyError::InvalidParameters)?;
        if self.two_exponent < 3
            || self.three_exponent == 0
            || n > 1024
            || d > 81
            || self.prime > 100_000
            || !is_prime(self.prime)
            || self.prime % 4 != 3
            || !matches!(self.prime % 7, 3 | 5 | 6)
            || !(self.prime + 1).is_multiple_of(n * d)
            || n * n <= d * d * 7
        {
            return Err(ToyError::InvalidParameters);
        }
        Ok((n, d))
    }
}

/// Public toy level structure.  The 3-primary basis is published to avoid a
/// separate canonical basis algorithm in this bounded experiment.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ToyPublicKey {
    pub parameters: ToyParameters,
    pub curve: Curve<Fq>,
    pub eigen_p: Pt<Fq>,
    pub eigen_q: Pt<Fq>,
    pub sender_p: Pt<Fq>,
    pub sender_q: Pt<Fq>,
}

impl ToyPublicKey {
    pub fn validate(&self) -> Result<(), ToyError> {
        let (n, d) = self.parameters.validate()?;
        let f = Zp2::new(self.parameters.prime);
        let e = &self.curve;
        if !is_smooth(&f, e)
            || [self.eigen_p, self.eigen_q, self.sender_p, self.sender_q]
                .iter()
                .any(|p| !on_curve(&f, e, p))
            || !exact_order(&f, e, &self.eigen_p, n, 2)
            || !exact_order(&f, e, &self.eigen_q, n, 2)
            || pmul(&f, e, &self.eigen_p, (n / 2) as u128)
                == pmul(&f, e, &self.eigen_q, (n / 2) as u128)
            || !exact_order(&f, e, &self.sender_p, d, 3)
            || !exact_order(&f, e, &self.sender_q, d, 3)
        {
            return Err(ToyError::InvalidPublicKey);
        }
        let p3 = pmul(&f, e, &self.sender_p, (d / 3) as u128);
        let q3 = pmul(&f, e, &self.sender_q, (d / 3) as u128);
        if q3 == p3 || q3 == neg(&f, &p3) {
            return Err(ToyError::InvalidPublicKey);
        }
        Ok(())
    }
}

#[derive(Debug)]
pub struct ToySecretKey {
    pub public: ToyPublicKey,
    /// 2-torsion kernel of a CM endomorphism of norm 2.
    pub theta_kernel: Pt<Fq>,
    /// Scaling factor for the isomorphism from the 2-isogeny codomain.
    pub theta_scale: Fq,
    pub trace: i8,
    pub eigenvalues: (u64, u64),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ToyCiphertext {
    pub curve: Curve<Fq>,
    pub s: Pt<Fq>,
    pub t: Pt<Fq>,
}

pub type SharedKey = [u8; 32];

#[derive(Debug)]
pub struct ToyKaniTranscript {
    /// Number of explicit product-to-product quotients.
    pub product_steps: u32,
    /// The remaining first-half theta chain.  Its final quotient is generally
    /// a Jacobian; the unimplemented second half would split it.
    pub first_half: ChainResult<Fq>,
}

/// Evaluated prefix in the special-key regime `N > d²Δ'`.  An intermediate
/// product still needs explicit factor maps before the chain can continue.
#[derive(Debug)]
pub struct ToySingleChain {
    pub product_steps: u32,
    pub result: ChainResult<Fq>,
    pub remaining_steps: u32,
    pub c: u64,
    pub cp: u64,
}

fn exact_order(f: &Zp2, e: &Curve<Fq>, p: &Pt<Fq>, order: u64, prime: u64) -> bool {
    pmul(f, e, p, order as u128) == Pt::Inf && pmul(f, e, p, (order / prime) as u128) != Pt::Inf
}

fn torsion_point(f: &Zp2, e: &Curve<Fq>, order: u64, prime: u64, rng: &mut Rng) -> Option<Pt<Fq>> {
    for _ in 0..4096 {
        let r = random_point_f(f, e, rng);
        let p = pmul(f, e, &r, ((f.char() + 1) / order) as u128);
        if exact_order(f, e, &p, order, prime) {
            return Some(p);
        }
    }
    None
}

fn theta_eval(f: &Zp2, iso: &VeluIso<Zp2>, scale: Fq, point: &Pt<Fq>) -> Pt<Fq> {
    match iso.eval(f, point) {
        Pt::Inf => Pt::Inf,
        Pt::Aff(x, y) => {
            let u2 = f.sq(scale);
            Pt::Aff(f.div(x, u2), f.div(y, f.mul(u2, scale)))
        }
    }
}

/// Scale for an isomorphism E1 -> E2, (x,y) -> (x/u²,y/u³).
/// The toy CM walk uses only generic j, where both coefficients are nonzero.
fn toy_curve_iso_scale(f: &Zp2, e1: &Curve<Fq>, e2: &Curve<Fq>) -> Option<Fq> {
    if e1 == e2 {
        return Some(f.one());
    }
    if [e1.a, e1.b, e2.a, e2.b].iter().any(|&x| f.is_zero(x)) {
        return None;
    }
    let u2 = f.div(f.mul(e1.b, e2.a), f.mul(e2.b, e1.a));
    let u = f.sqrt(u2)?;
    (e1.a == f.mul(e2.a, f.sq(u2)) && e1.b == f.mul(e2.b, f.mul(f.sq(u2), u2))).then_some(u)
}

fn toy_iso_inverse(f: &Zp2, u: Fq, p: Pt<Fq>) -> Pt<Fq> {
    match p {
        Pt::Inf => Pt::Inf,
        Pt::Aff(x, y) => {
            let u2 = f.sq(u);
            Pt::Aff(f.mul(u2, x), f.mul(f.mul(u2, u), y))
        }
    }
}

fn toy_iso_forward(f: &Zp2, u: Fq, p: Pt<Fq>) -> Pt<Fq> {
    match p {
        Pt::Inf => Pt::Inf,
        Pt::Aff(x, y) => {
            let u2 = f.sq(u);
            Pt::Aff(f.div(x, u2), f.div(y, f.mul(u2, u)))
        }
    }
}

fn find_theta(f: &Zp2, e: &Curve<Fq>, rng: &mut Rng) -> Option<(Pt<Fq>, Fq, i8)> {
    let roots = poly::roots(f, &vec![e.b, e.a, f.zero(), f.one()], rng);
    for x in roots {
        let kernel = Pt::Aff(x, f.zero());
        let iso = velu_cyclic(f, e, &kernel, 2);
        if iso.cod.a == f.zero() || iso.cod.b == f.zero() {
            continue;
        }
        // For E(a,b) -> E(a',b'), x -> x/u^2, y -> y/u^3.
        let u2 = f.div(f.mul(iso.cod.b, e.a), f.mul(e.b, iso.cod.a));
        let Some(scale) = f.sqrt(u2) else { continue };
        if iso.cod.a != f.mul(f.sq(u2), e.a) || iso.cod.b != f.mul(f.mul(f.sq(u2), u2), e.b) {
            continue;
        }
        let r = torsion_point(f, e, 8, 2, rng)?;
        let tr = theta_eval(f, &iso, scale, &r);
        let tr2 = theta_eval(f, &iso, scale, &tr);
        let two_r = pmul(f, e, &r, 2);
        let plus = padd(f, e, &padd(f, e, &tr2, &neg(f, &tr)), &two_r);
        let minus = padd(f, e, &padd(f, e, &tr2, &tr), &two_r);
        if plus == Pt::Inf {
            return Some((kernel, scale, 1));
        }
        if minus == Pt::Inf {
            return Some((kernel, scale, -1));
        }
    }
    None
}

fn eigen_roots(n: u64, trace: i8) -> Option<(u64, u64)> {
    let roots: Vec<u64> = (0..n)
        .filter(|&r| (r * r + 2 + if trace == 1 { n - r } else { r }) % n == 0)
        .collect();
    let even = roots.iter().copied().find(|r| r % 2 == 0)?;
    let odd = roots.iter().copied().find(|r| r % 2 == 1)?;
    Some((even, odd))
}

fn eigen_point(
    f: &Zp2,
    e: &Curve<Fq>,
    iso: &VeluIso<Zp2>,
    scale: Fq,
    n: u64,
    own: u64,
    other: u64,
    rng: &mut Rng,
) -> Option<Pt<Fq>> {
    for _ in 0..4096 {
        let r = torsion_point(f, e, n, 2, rng)?;
        let p = padd(
            f,
            e,
            &theta_eval(f, iso, scale, &r),
            &neg(f, &pmul(f, e, &r, other as u128)),
        );
        if exact_order(f, e, &p, n, 2)
            && theta_eval(f, iso, scale, &p) == pmul(f, e, &p, own as u128)
        {
            return Some(p);
        }
    }
    None
}

/// Constructs the public CM toy key.  The degree-2 endomorphism is genuine;
/// the CM ring of this curve is public, so this key has no secrecy claim.
pub fn toy_keygen(parameters: ToyParameters, rng: &mut Rng) -> Result<ToySecretKey, ToyError> {
    let (n, d) = parameters.validate()?;
    let f = Zp2::new(parameters.prime);
    let j = f.neg(f.from_u64(3375));
    let e = from_j(&f, j);
    if !is_smooth(&f, &e) {
        return Err(ToyError::ConstructionFailed);
    }
    let (kernel, scale, trace) = find_theta(&f, &e, rng).ok_or(ToyError::ConstructionFailed)?;
    let iso = velu_cyclic(&f, &e, &kernel, 2);
    let (r0, r1) = eigen_roots(n, trace).ok_or(ToyError::ConstructionFailed)?;
    let p = eigen_point(&f, &e, &iso, scale, n, r0, r1, rng).ok_or(ToyError::ConstructionFailed)?;
    let q = eigen_point(&f, &e, &iso, scale, n, r1, r0, rng).ok_or(ToyError::ConstructionFailed)?;
    if pmul(&f, &e, &p, (n / 2) as u128) == pmul(&f, &e, &q, (n / 2) as u128) {
        return Err(ToyError::ConstructionFailed);
    }
    let pd = torsion_point(&f, &e, d, 3, rng).ok_or(ToyError::ConstructionFailed)?;
    let pd3 = pmul(&f, &e, &pd, (d / 3) as u128);
    let qd = (0..4096)
        .find_map(|_| {
            let candidate = torsion_point(&f, &e, d, 3, rng)?;
            let c3 = pmul(&f, &e, &candidate, (d / 3) as u128);
            if c3 != pd3 && c3 != neg(&f, &pd3) {
                Some(candidate)
            } else {
                None
            }
        })
        .ok_or(ToyError::ConstructionFailed)?;
    let sk = ToySecretKey {
        public: ToyPublicKey {
            parameters,
            curve: e,
            eigen_p: p,
            eigen_q: q,
            sender_p: pd,
            sender_q: qd,
        },
        theta_kernel: kernel,
        theta_scale: scale,
        trace,
        eigenvalues: (r0, r1),
    };
    sk.public.validate()?;
    Ok(sk)
}

fn odd_mask(seed: u64, order: u64, label: &[u8]) -> u64 {
    let mut input = b"Eigen-Seta toy mask v1".to_vec();
    input.extend_from_slice(label);
    input.extend_from_slice(&seed.to_le_bytes());
    let h = sha256(&input);
    (u64::from_le_bytes(h[..8].try_into().unwrap()) % (order / 2)) * 2 + 1
}

fn serialize_point(out: &mut Vec<u8>, p: Pt<Fq>) {
    match p {
        Pt::Inf => out.push(0),
        Pt::Aff((x0, x1), (y0, y1)) => {
            out.push(1);
            for value in [x0, x1, y0, y1] {
                out.extend_from_slice(&value.to_le_bytes());
            }
        }
    }
}

fn derive_key(seed: u64, ct: &ToyCiphertext) -> SharedKey {
    let mut input = b"Eigen-Seta toy shared key v1".to_vec();
    input.extend_from_slice(&seed.to_le_bytes());
    for value in [ct.curve.a.0, ct.curve.a.1, ct.curve.b.0, ct.curve.b.1] {
        input.extend_from_slice(&value.to_le_bytes());
    }
    serialize_point(&mut input, ct.s);
    serialize_point(&mut input, ct.t);
    sha256(&input)
}

fn encaps_for_seed(pk: &ToyPublicKey, seed: u64) -> Result<(ToyCiphertext, SharedKey), ToyError> {
    let (n, d) = pk.parameters.validate()?;
    if seed >= d {
        return Err(ToyError::InvalidParameters);
    }
    let f = Zp2::new(pk.parameters.prime);
    let k = padd(
        &f,
        &pk.curve,
        &pk.sender_p,
        &pmul(&f, &pk.curve, &pk.sender_q, seed as u128),
    );
    if !exact_order(&f, &pk.curve, &k, d, 3) {
        return Err(ToyError::ConstructionFailed);
    }
    let phi = velu_cyclic(&f, &pk.curve, &k, d);
    let a = odd_mask(seed, n, b"a");
    let b = odd_mask(seed, n, b"b");
    let s = pmul(&f, &phi.cod, &phi.eval(&f, &pk.eigen_p), a as u128);
    let t = pmul(&f, &phi.cod, &phi.eval(&f, &pk.eigen_q), b as u128);
    let ct = ToyCiphertext {
        curve: phi.cod,
        s,
        t,
    };
    Ok((ct, derive_key(seed, &ct)))
}

pub fn toy_encaps(
    pk: &ToyPublicKey,
    rng: &mut Rng,
) -> Result<(ToyCiphertext, SharedKey), ToyError> {
    pk.validate()?;
    let (_, d) = pk.parameters.validate()?;
    encaps_for_seed(pk, rng.below(d))
}

/// Public exhaustive decoder for the bounded toy.  This is also a direct
/// demonstration that the toy ciphertext and key are recoverable by anyone.
pub fn toy_public_recover(pk: &ToyPublicKey, ct: &ToyCiphertext) -> Result<SharedKey, ToyError> {
    pk.validate()?;
    let (_, d) = pk.parameters.validate()?;
    let f = Zp2::new(pk.parameters.prime);
    if !is_smooth(&f, &ct.curve)
        || !on_curve(&f, &ct.curve, &ct.s)
        || !on_curve(&f, &ct.curve, &ct.t)
        || ct.s == Pt::Inf
        || ct.t == Pt::Inf
    {
        return Err(ToyError::InvalidCiphertext);
    }
    for seed in 0..d {
        let (candidate, key) = encaps_for_seed(pk, seed)?;
        if candidate == *ct {
            return Ok(key);
        }
    }
    Err(ToyError::InvalidCiphertext)
}

/// Same bounded public recovery path, exposed with the usual KEM API shape.
pub fn toy_decaps(sk: &ToySecretKey, ct: &ToyCiphertext) -> Result<SharedKey, ToyError> {
    toy_public_recover(&sk.public, ct)
}

/// Evaluate the first half of the proposed Kani kernel obtained from the
/// masked ciphertext basis and the secret eigenvalues.  At toy parameters,
/// a 4-division table supplies the lifts.  The second half and the component
/// maps needed to recover the sender kernel remain unimplemented.
pub fn toy_kani_first_half(
    sk: &ToySecretKey,
    ct: &ToyCiphertext,
) -> Result<ToyKaniTranscript, ToyError> {
    let (product_steps, prefix, _, _) = toy_kani_eval(sk, ct, true, false, &[])?;
    let first_half = prefix.result;
    Ok(ToyKaniTranscript {
        product_steps,
        first_half,
    })
}

/// First half for the Rosati-dual matrix `F† = W F W`, whose N-torsion
/// kernel is the factor-swap of the original first-half kernel.  Evaluating
/// this quotient is a prerequisite, not the missing middle-Jacobian map.
pub fn toy_kani_dual_first_half(
    sk: &ToySecretKey,
    ct: &ToyCiphertext,
) -> Result<ToyKaniTranscript, ToyError> {
    let (product_steps, prefix, _, _) = toy_kani_eval(sk, ct, true, true, &[])?;
    Ok(ToyKaniTranscript {
        product_steps,
        first_half: prefix.result,
    })
}

/// Evaluate the *single* Kani chain when `N > d²Δ'`, including supplied
/// product points.  Returns at the first intermediate split because the
/// current theta evaluator cannot yet continue from a product.  This is a
/// diagnostic for the special-key regime, not trapdoor decapsulation.
pub fn toy_kani_single_chain(
    sk: &ToySecretKey,
    ct: &ToyCiphertext,
    extra: &[(Pt<Fq>, Pt<Fq>)],
) -> Result<ToySingleChain, ToyError> {
    let (product_steps, prefix, c, cp) = toy_kani_eval(sk, ct, false, false, extra)?;
    Ok(ToySingleChain {
        product_steps,
        result: prefix.result,
        remaining_steps: prefix.remaining_steps,
        c,
        cp,
    })
}

fn toy_kani_eval(
    sk: &ToySecretKey,
    ct: &ToyCiphertext,
    square_norm: bool,
    dual_graph: bool,
    extra: &[(Pt<Fq>, Pt<Fq>)],
) -> Result<(u32, ChainPrefix<Fq>, u64, u64), ToyError> {
    let pk = &sk.public;
    pk.validate()?;
    let (n, d) = pk.parameters.validate()?;
    if !(pk.parameters.prime + 1).is_multiple_of(4 * n) {
        return Err(ToyError::InvalidParameters);
    }
    let f = Zp2::new(pk.parameters.prime);
    if !is_smooth(&f, &ct.curve)
        || !on_curve(&f, &ct.curve, &ct.s)
        || !on_curve(&f, &ct.curve, &ct.t)
        || !exact_order(&f, &ct.curve, &ct.s, n, 2)
        || !exact_order(&f, &ct.curve, &ct.t, n, 2)
        || pmul(&f, &ct.curve, &ct.s, (n / 2) as u128)
            == pmul(&f, &ct.curve, &ct.t, (n / 2) as u128)
    {
        return Err(ToyError::InvalidCiphertext);
    }
    if extra
        .iter()
        .any(|(x, y)| !on_curve(&f, &ct.curve, x) || !on_curve(&f, &ct.curve, y))
    {
        return Err(ToyError::InvalidCiphertext);
    }
    let norm = if square_norm { n * n } else { n };
    let target = norm
        .checked_sub(d * d * 7)
        .ok_or(ToyError::InvalidParameters)?;
    let (c, cp) = (0..(n / 2))
        .rev()
        .map(|i| 2 * i)
        .find_map(|c| {
            (1..n)
                .step_by(2)
                .find(|&cp| c * c + cp * cp == target)
                .map(|cp| (c, cp))
        })
        .ok_or(ToyError::ConstructionFailed)?;
    let scalar = |r: u64| -> u64 {
        ((d as i128 * (2 * r as i128 - sk.trace as i128) + c as i128).rem_euclid(n as i128)) as u64
    };
    let mut graph = [
        (
            pmul(&f, &ct.curve, &ct.s, cp as u128),
            neg(
                &f,
                &pmul(&f, &ct.curve, &ct.s, scalar(sk.eigenvalues.0) as u128),
            ),
        ),
        (
            pmul(&f, &ct.curve, &ct.t, cp as u128),
            neg(
                &f,
                &pmul(&f, &ct.curve, &ct.t, scalar(sk.eigenvalues.1) as u128),
            ),
        ),
    ];
    if dual_graph {
        graph = graph.map(|(x, y)| (y, x));
    }
    let mut rng = Rng::new(0xe156_5e7a);
    let order = 4 * n;
    let basis_p =
        torsion_point(&f, &ct.curve, order, 2, &mut rng).ok_or(ToyError::ConstructionFailed)?;
    let p2 = pmul(&f, &ct.curve, &basis_p, (order / 2) as u128);
    let basis_q = (0..4096)
        .find_map(|_| {
            let q = torsion_point(&f, &ct.curve, order, 2, &mut rng)?;
            (pmul(&f, &ct.curve, &q, (order / 2) as u128) != p2).then_some(q)
        })
        .ok_or(ToyError::ConstructionFailed)?;
    // Enumerate the 4-division map on E[4N].  The same table lifts the four
    // components of the Kani graph; the cost is O(N^2) and is toy-only.
    let mut lifts = HashMap::new();
    let mut row = Pt::Inf;
    for _ in 0..n {
        let mut point = row;
        for _ in 0..n {
            lifts.insert(pmul(&f, &ct.curve, &point, 4), point);
            point = padd(&f, &ct.curve, &point, &basis_q);
        }
        row = padd(&f, &ct.curve, &row, &basis_p);
    }
    let lifted = graph.map(|(left, right)| {
        Ok::<_, ToyError>((
            *lifts.get(&left).ok_or(ToyError::ConstructionFailed)?,
            *lifts.get(&right).ok_or(ToyError::ConstructionFailed)?,
        ))
    });
    let [a, b] = lifted;
    let mut k = [a?, b?];
    let mut curves = [ct.curve, ct.curve];
    let mut pushed = extra.to_vec();
    let mut remaining = pk.parameters.two_exponent;
    let mut product_steps = 0;
    // When the 2-torsion kernel is diagonal, (x,y) -> (x+y,x-y) is a
    // degree-4 quotient with precisely that kernel.  It remains on E x E.
    // The generic theta chain cannot represent this first split quotient.
    while remaining > 1 {
        let g0 = (
            pmul(&f, &curves[0], &k[0].0, 1u128 << (remaining + 1)),
            pmul(&f, &curves[1], &k[0].1, 1u128 << (remaining + 1)),
        );
        let g1 = (
            pmul(&f, &curves[0], &k[1].0, 1u128 << (remaining + 1)),
            pmul(&f, &curves[1], &k[1].1, 1u128 << (remaining + 1)),
        );
        let g2 = (
            padd(&f, &curves[0], &g0.0, &g1.0),
            padd(&f, &curves[1], &g0.1, &g1.1),
        );
        let graph_scale = toy_curve_iso_scale(&f, &curves[0], &curves[1]).filter(|&u| {
            g0.0 != Pt::Inf
                && g1.0 != Pt::Inf
                && g0.0 != g1.0
                && toy_iso_inverse(&f, u, g0.1) == g0.0
                && toy_iso_inverse(&f, u, g1.1) == g1.0
        });
        if let Some(u) = graph_scale {
            k = k.map(|(x, y)| {
                let y = toy_iso_inverse(&f, u, y);
                (
                    padd(&f, &curves[0], &x, &y),
                    padd(&f, &curves[0], &x, &neg(&f, &y)),
                )
            });
            for (x, y) in &mut pushed {
                let (left, right) = (*x, toy_iso_inverse(&f, u, *y));
                *x = padd(&f, &curves[0], &left, &right);
                *y = padd(&f, &curves[0], &left, &neg(&f, &right));
            }
            curves[1] = curves[0];
        } else {
            let subgroup = [g0, g1, g2];
            let left = subgroup
                .iter()
                .find_map(|(x, y)| (*x != Pt::Inf && *y == Pt::Inf).then_some(*x));
            let right = subgroup
                .iter()
                .find_map(|(x, y)| (*x == Pt::Inf && *y != Pt::Inf).then_some(*y));
            let (Some(left), Some(right)) = (left, right) else {
                break;
            };
            let iso_left = velu_cyclic(&f, &curves[0], &left, 2);
            let iso_right = velu_cyclic(&f, &curves[1], &right, 2);
            k = k.map(|(x, y)| (iso_left.eval(&f, &x), iso_right.eval(&f, &y)));
            for (x, y) in &mut pushed {
                *x = iso_left.eval(&f, x);
                *y = iso_right.eval(&f, y);
            }
            curves = [iso_left.cod, iso_right.cod];
        }
        remaining -= 1;
        product_steps += 1;
    }
    let res = if square_norm {
        ChainPrefix {
            result: chain(&f, &curves[0], &curves[1], k, remaining, &pushed)
                .ok_or(ToyError::ConstructionFailed)?,
            remaining_steps: 0,
        }
    } else {
        chain_until_first_split(&f, &curves[0], &curves[1], k, remaining, &pushed)
            .ok_or(ToyError::ConstructionFailed)?
    };
    Ok((product_steps, res, c, cp))
}

/// Recover coordinates in a 3-power torsion basis by digit lifting.  This is
/// a small native analogue of the smooth-order DLP in the Sage reference;
/// it never enumerates sender ciphertexts.
fn bidlp3(
    f: &Zp2,
    e: &Curve<Fq>,
    r: Pt<Fq>,
    p: Pt<Fq>,
    q: Pt<Fq>,
    exponent: u32,
) -> Option<(u64, u64)> {
    let d = 3u64.checked_pow(exponent)?;
    let p3 = pmul(f, e, &p, (d / 3) as u128);
    let q3 = pmul(f, e, &q, (d / 3) as u128);
    let (mut x, mut y, mut place) = (0u64, 0u64, 1u64);
    for k in 0..exponent {
        let residual = padd(
            f,
            e,
            &r,
            &neg(
                f,
                &padd(f, e, &pmul(f, e, &p, x as u128), &pmul(f, e, &q, y as u128)),
            ),
        );
        let target = pmul(f, e, &residual, 3u64.pow(exponent - 1 - k) as u128);
        let digits = (0..3u64).find_map(|dx| {
            (0..3u64)
                .find(|&dy| {
                    target
                        == padd(
                            f,
                            e,
                            &pmul(f, e, &p3, dx as u128),
                            &pmul(f, e, &q3, dy as u128),
                        )
                })
                .map(|dy| (dx, dy))
        })?;
        x += digits.0 * place;
        y += digits.1 * place;
        place *= 3;
    }
    (padd(f, e, &pmul(f, e, &p, x as u128), &pmul(f, e, &q, y as u128)) == r).then_some((x, y))
}

fn toy_single_chain_action_x(
    sk: &ToySecretKey,
    ct: &ToyCiphertext,
    x1: Pt<Fq>,
    x2: Pt<Fq>,
) -> Result<(u64, [Fq; 2]), ToyError> {
    let f = Zp2::new(sk.public.parameters.prime);
    let payload = [(x1, Pt::Inf), (Pt::Inf, x1), (x2, Pt::Inf), (Pt::Inf, x2)];
    let quotient = toy_kani_single_chain(sk, ct, &payload)?;
    if quotient.remaining_steps != 0 || quotient.result.images.len() != 4 {
        return Err(ToyError::ConstructionFailed);
    }
    let null = quotient
        .result
        .nulls
        .last()
        .ok_or(ToyError::ConstructionFailed)?;
    let chart = crate::theta::find_split_chart(&f, null, &quotient.result.images, 50_000)
        .ok_or(ToyError::ConstructionFailed)?;
    let factor_curves = [
        crate::theta::curve_from_theta1(&f, &chart.factors.0)
            .ok_or(ToyError::ConstructionFailed)?,
        crate::theta::curve_from_theta1(&f, &chart.factors.1)
            .ok_or(ToyError::ConstructionFailed)?,
    ];
    let scales = [
        toy_curve_iso_scale(&f, &ct.curve, &factor_curves[0])
            .ok_or(ToyError::ConstructionFailed)?,
        toy_curve_iso_scale(&f, &ct.curve, &factor_curves[1])
            .ok_or(ToyError::ConstructionFailed)?,
    ];
    let mut coordinates = [[f.zero(); 2]; 4];
    for (i, image) in quotient.result.images.iter().enumerate() {
        let transformed = chart.apply(&f, image).ok_or(ToyError::ConstructionFailed)?;
        let factors =
            crate::theta::factor_product(&f, &transformed).ok_or(ToyError::ConstructionFailed)?;
        for j in 0..2 {
            let theta_null = if j == 0 {
                &chart.factors.0
            } else {
                &chart.factors.1
            };
            let theta_pt = if j == 0 { &factors.0 } else { &factors.1 };
            let x = crate::theta::x_on_theta1_curve(&f, theta_null, theta_pt)
                .ok_or(ToyError::ConstructionFailed)?;
            coordinates[i][j] = f.mul(f.sq(scales[j]), x);
        }
    }
    let expected = [x1, x2].map(|x| match pmul(&f, &ct.curve, &x, quotient.cp as u128) {
        Pt::Aff(x, _) => x,
        Pt::Inf => f.zero(),
    });
    let first_factor = (0..2)
        .find(|&first| {
            let second = 1 - first;
            coordinates[0][second] == expected[0]
                && coordinates[2][second] == expected[1]
                && coordinates[1][first] == expected[0]
                && coordinates[3][first] == expected[1]
        })
        .ok_or(ToyError::ConstructionFailed)?;
    Ok((
        quotient.c,
        [coordinates[0][first_factor], coordinates[2][first_factor]],
    ))
}

fn normalized_theta(f: &Zp2, v: crate::theta::Th<Fq>) -> Option<crate::theta::Th<Fq>> {
    let z = v.iter().position(|&x| !f.is_zero(x))?;
    let inv = f.inv(v[z]);
    Some(v.map(|x| f.mul(x, inv)))
}

fn point_x(p: Pt<Fq>) -> Option<Fq> {
    match p {
        Pt::Aff(x, _) => Some(x),
        Pt::Inf => None,
    }
}

/// Exponential toy oracle for the missing second Kani half.  On E[d]^2 it
/// inverts the first quotient by a complete image table, then enumerates the
/// 11,520 projective two-qubit Clifford changes to find middle-Jacobian
/// automorphisms with the prescribed off-diagonal `[c']` action.  This is
/// bounded to d <= 9 and is *not* a scalable or secure second-half evaluator.
fn toy_two_halves_action_x(
    sk: &ToySecretKey,
    ct: &ToyCiphertext,
    p: Pt<Fq>,
    q: Pt<Fq>,
) -> Result<Vec<(u64, [Option<Fq>; 2])>, ToyError> {
    let (n, d) = sk.public.parameters.validate()?;
    if d > 9 || n > d * d * 7 || n * n <= d * d * 7 {
        return Err(ToyError::InvalidParameters);
    }
    let f = Zp2::new(sk.public.parameters.prime);
    let mut points = Vec::with_capacity((d * d) as usize);
    for i in 0..d {
        for j in 0..d {
            points.push(padd(
                &f,
                &ct.curve,
                &pmul(&f, &ct.curve, &p, i as u128),
                &pmul(&f, &ct.curve, &q, j as u128),
            ));
        }
    }
    let mut extras = Vec::with_capacity(points.len() * points.len());
    for &a in &points {
        for &b in &points {
            extras.push((a, b));
        }
    }
    let (_, first, c, cp) = toy_kani_eval(sk, ct, true, false, &extras)?;
    if first.remaining_steps != 0 || first.result.images.len() != extras.len() {
        return Err(ToyError::ConstructionFailed);
    }
    let null = *first
        .result
        .nulls
        .last()
        .ok_or(ToyError::ConstructionFailed)?;
    if first.result.split.is_some() {
        return Err(ToyError::ConstructionFailed);
    }
    let mut table = HashMap::new();
    for (&v, &source) in first.result.images.iter().zip(&extras) {
        let key = normalized_theta(&f, v).ok_or(ToyError::ConstructionFailed)?;
        if let Some(previous) = table.insert(key, source) {
            let negative = (neg(&f, &source.0), neg(&f, &source.1));
            if previous != source && previous != negative {
                return Err(ToyError::ConstructionFailed);
            }
        }
    }
    if table.len() != extras.len().div_ceil(2) {
        return Err(ToyError::ConstructionFailed);
    }
    let imaginary = f.sqrt(f.neg(f.one())).ok_or(ToyError::ConstructionFailed)?;
    let apply = |path: &[u8], v: crate::theta::Th<Fq>| {
        path.iter().fold(v, |w, &g| {
            crate::theta::split_chart_step(&f, imaginary, g, &w)
        })
    };
    let identity =
        [0, 1, 2, 3].map(|j| [0, 1, 2, 3].map(|i| if i == j { f.one() } else { f.zero() }));
    let normalize_matrix = |m: [crate::theta::Th<Fq>; 4]| {
        let z = m.iter().flatten().copied().find(|&x| !f.is_zero(x))?;
        let inv = f.inv(z);
        Some(m.map(|col| col.map(|x| f.mul(x, inv))))
    };
    let len = points.len();
    let sample_indices = [d as usize * len, d as usize, len, 1usize];
    let expected_cp = [p, q].map(|x| point_x(pmul(&f, &ct.curve, &x, cp as u128)));
    let mut seen = HashSet::new();
    let mut todo = VecDeque::from([(identity, Vec::<u8>::new())]);
    let mut actions = HashSet::new();
    while let Some((matrix, path)) = todo.pop_front() {
        let matrix = normalize_matrix(matrix).ok_or(ToyError::ConstructionFailed)?;
        if !seen.insert(matrix) {
            continue;
        }
        if crate::theta::proj_eq(&f, &apply(&path, null), &null) {
            let mapped: Option<Vec<_>> = sample_indices
                .iter()
                .map(|&idx| {
                    table
                        .get(&normalized_theta(
                            &f,
                            apply(&path, first.result.images[idx]),
                        )?)
                        .copied()
                })
                .collect();
            if let Some(mapped) = mapped {
                let outputs: Vec<_> = mapped
                    .iter()
                    .map(|&(a, b)| {
                        (
                            pmul(&f, &ct.curve, &b, n as u128),
                            pmul(&f, &ct.curve, &a, n as u128),
                        )
                    })
                    .collect();
                if point_x(outputs[0].1) == expected_cp[0]
                    && point_x(outputs[1].0) == expected_cp[0]
                    && point_x(outputs[2].1) == expected_cp[1]
                    && point_x(outputs[3].0) == expected_cp[1]
                {
                    actions.insert((c, [point_x(outputs[0].0), point_x(outputs[2].0)]));
                }
            }
        }
        if seen.len() >= 50_000 {
            return Err(ToyError::ConstructionFailed);
        }
        for g in 0..6 {
            let mut next = path.clone();
            next.push(g);
            todo.push_back((
                matrix.map(|v| crate::theta::split_chart_step(&f, imaginary, g, &v)),
                next,
            ));
        }
    }
    if actions.is_empty() {
        return Err(ToyError::ConstructionFailed);
    }
    Ok(actions.into_iter().collect())
}

/// Experimental *secret-dependent* decoder for the CM fixtures.  It uses a
/// direct single Kani chain when `N>d²Δ'` and a full odd-torsion image table
/// for the two-halves toy regime.  Both are bounded research paths; the
/// latter is exponential and is not the scalable random-curve decoder.
pub fn toy_trapdoor_recover_seed(sk: &ToySecretKey, ct: &ToyCiphertext) -> Result<u64, ToyError> {
    let pk = &sk.public;
    pk.validate()?;
    let (n, d) = pk.parameters.validate()?;
    let f = Zp2::new(pk.parameters.prime);
    if !is_smooth(&f, &ct.curve)
        || !on_curve(&f, &ct.curve, &ct.s)
        || !on_curve(&f, &ct.curve, &ct.t)
    {
        return Err(ToyError::InvalidCiphertext);
    }
    let mut rng = Rng::new(0xe156_5e7b);
    let x1 = torsion_point(&f, &ct.curve, d, 3, &mut rng).ok_or(ToyError::InvalidCiphertext)?;
    let x3 = pmul(&f, &ct.curve, &x1, (d / 3) as u128);
    let x2 = (0..4096)
        .find_map(|_| {
            let q = torsion_point(&f, &ct.curve, d, 3, &mut rng)?;
            let q3 = pmul(&f, &ct.curve, &q, (d / 3) as u128);
            (q3 != x3 && q3 != neg(&f, &x3)).then_some(q)
        })
        .ok_or(ToyError::InvalidCiphertext)?;
    let actions = if n > d * d * 7 {
        let (c, alpha_x) = toy_single_chain_action_x(sk, ct, x1, x2)?;
        vec![(c, alpha_x.map(Some))]
    } else {
        toy_two_halves_action_x(sk, ct, x1, x2)?
    };
    for (c, alpha_x) in actions {
        let alpha_candidates = alpha_x.map(|x| match x {
            Some(x) => {
                let y = f
                    .sqrt(crate::curve::rhs(&f, &ct.curve, x))
                    .ok_or(ToyError::ConstructionFailed)?;
                Ok::<_, ToyError>(Pt::Aff(x, y))
            }
            None => Ok(Pt::Inf),
        });
        let [a1, a2] = alpha_candidates;
        let (a1, a2) = (a1?, a2?);
        for sign1 in [false, true] {
            for sign2 in [false, true] {
                let alpha1 = if sign1 { neg(&f, &a1) } else { a1 };
                let alpha2 = if sign2 { neg(&f, &a2) } else { a2 };
                let psi1 = padd(
                    &f,
                    &ct.curve,
                    &alpha1,
                    &neg(&f, &pmul(&f, &ct.curve, &x1, c as u128)),
                );
                let psi2 = padd(
                    &f,
                    &ct.curve,
                    &alpha2,
                    &neg(&f, &pmul(&f, &ct.curve, &x2, c as u128)),
                );
                for i in 0..d {
                    for j in 0..d {
                        let khat = padd(
                            &f,
                            &ct.curve,
                            &pmul(&f, &ct.curve, &x1, i as u128),
                            &pmul(&f, &ct.curve, &x2, j as u128),
                        );
                        if !exact_order(&f, &ct.curve, &khat, d, 3)
                            || padd(
                                &f,
                                &ct.curve,
                                &pmul(&f, &ct.curve, &psi1, i as u128),
                                &pmul(&f, &ct.curve, &psi2, j as u128),
                            ) != Pt::Inf
                        {
                            continue;
                        }
                        let dual = velu_cyclic(&f, &ct.curve, &khat, d);
                        let Some(u) = toy_curve_iso_scale(&f, &dual.cod, &pk.curve) else {
                            continue;
                        };
                        for x in [x1, x2, padd(&f, &ct.curve, &x1, &x2)] {
                            let image = toy_iso_forward(&f, u, dual.eval(&f, &x));
                            if !exact_order(&f, &pk.curve, &image, d, 3) {
                                continue;
                            }
                            let Some((px, qx)) = bidlp3(
                                &f,
                                &pk.curve,
                                image,
                                pk.sender_p,
                                pk.sender_q,
                                pk.parameters.three_exponent,
                            ) else {
                                continue;
                            };
                            if px % 3 == 0 {
                                continue;
                            }
                            let Some(inv) = (1..d).find(|&v| px * v % d == 1) else {
                                continue;
                            };
                            let seed = qx * inv % d;
                            if let Ok((candidate, _)) = encaps_for_seed(pk, seed) {
                                if candidate == *ct {
                                    return Ok(seed);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    Err(ToyError::InvalidCiphertext)
}

/// Secret-dependent toy decapsulation in the special single-chain and bounded
/// two-halves-oracle regimes.  No security claim: the CM endomorphism is
/// public, the sender space is tiny, and the two-halves path enumerates E[d]².
pub fn toy_trapdoor_decaps(sk: &ToySecretKey, ct: &ToyCiphertext) -> Result<SharedKey, ToyError> {
    let seed = toy_trapdoor_recover_seed(sk, ct)?;
    let (candidate, key) = encaps_for_seed(&sk.public, seed)?;
    (candidate == *ct)
        .then_some(key)
        .ok_or(ToyError::InvalidCiphertext)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn direct_action_x(
        sk: &ToySecretKey,
        ct: &ToyCiphertext,
        seed: u64,
        c: u64,
        points: [Pt<Fq>; 2],
    ) -> [Option<Fq>; 2] {
        let f = Zp2::new(sk.public.parameters.prime);
        let d = 3u64.pow(sk.public.parameters.three_exponent);
        let kernel = padd(
            &f,
            &sk.public.curve,
            &sk.public.sender_p,
            &pmul(&f, &sk.public.curve, &sk.public.sender_q, seed as u128),
        );
        let phi = velu_cyclic(&f, &sk.public.curve, &kernel, d);
        let dual_kernel = phi.eval(&f, &sk.public.sender_q);
        let dual = velu_cyclic(&f, &ct.curve, &dual_kernel, d);
        let u = toy_curve_iso_scale(&f, &dual.cod, &sk.public.curve).unwrap();
        let check = toy_iso_forward(&f, u, dual.eval(&f, &phi.eval(&f, &sk.public.eigen_p)));
        let want = pmul(&f, &sk.public.curve, &sk.public.eigen_p, d as u128);
        let dual_sign = if check == want {
            false
        } else if check == neg(&f, &want) {
            true
        } else {
            panic!("dual normalization");
        };
        let theta_iso = velu_cyclic(&f, &sk.public.curve, &sk.theta_kernel, 2);
        points.map(|point| {
            let mut pre = toy_iso_forward(&f, u, dual.eval(&f, &point));
            if dual_sign {
                pre = neg(&f, &pre);
            }
            let theta_pre = theta_eval(&f, &theta_iso, sk.theta_scale, &pre);
            let theta_prime = padd(
                &f,
                &sk.public.curve,
                &pmul(&f, &sk.public.curve, &theta_pre, 2),
                &if sk.trace == 1 { neg(&f, &pre) } else { pre },
            );
            let psi = phi.eval(&f, &theta_prime);
            point_x(padd(
                &f,
                &ct.curve,
                &psi,
                &pmul(&f, &ct.curve, &point, c as u128),
            ))
        })
    }

    #[test]
    fn cm_key_and_all_sender_indices_round_trip() {
        let parameters = ToyParameters {
            prime: 16127,
            two_exponent: 8,
            three_exponent: 2,
        };
        let mut rng = Rng::new(0x51e7a);
        let sk = toy_keygen(parameters, &mut rng).unwrap();
        let f = Zp2::new(parameters.prime);
        let iso = velu_cyclic(&f, &sk.public.curve, &sk.theta_kernel, 2);
        let n = 1u64 << parameters.two_exponent;
        for (p, lambda) in [
            (sk.public.eigen_p, sk.eigenvalues.0),
            (sk.public.eigen_q, sk.eigenvalues.1),
        ] {
            assert_eq!(
                theta_eval(&f, &iso, sk.theta_scale, &p),
                pmul(&f, &sk.public.curve, &p, lambda as u128)
            );
            assert!(exact_order(&f, &sk.public.curve, &p, n, 2));
        }
        for seed in 0..9 {
            let (ct, key) = encaps_for_seed(&sk.public, seed).unwrap();
            assert_eq!(toy_decaps(&sk, &ct), Ok(key));
            assert_eq!(toy_public_recover(&sk.public, &ct), Ok(key));
            let mut tampered = ct;
            tampered.s = neg(&f, &ct.s);
            assert_eq!(toy_decaps(&sk, &tampered), Err(ToyError::InvalidCiphertext));
        }
        let (ct, key) = toy_encaps(&sk.public, &mut rng).unwrap();
        assert_eq!(toy_decaps(&sk, &ct), Ok(key));
    }

    #[test]
    fn invalid_parameters_rejected() {
        let p = ToyParameters {
            prime: 16127,
            two_exponent: 8,
            three_exponent: 2,
        };
        assert!(p.validate().is_ok());
        assert_eq!(
            ToyParameters { prime: 15, ..p }.validate(),
            Err(ToyError::InvalidParameters)
        );
        assert_eq!(
            ToyParameters {
                three_exponent: 5,
                ..p
            }
            .validate(),
            Err(ToyError::InvalidParameters)
        );
    }

    #[test]
    fn ciphertext_kani_first_half_reaches_jacobian() {
        let parameters = ToyParameters {
            prime: 16127,
            two_exponent: 6,
            three_exponent: 2,
        };
        let mut rng = Rng::new(0x44c0);
        let sk = toy_keygen(parameters, &mut rng).unwrap();
        let f = Zp2::new(parameters.prime);
        for seed in 0..9 {
            let (ct, key) = encaps_for_seed(&sk.public, seed).unwrap();
            let transcript = toy_kani_first_half(&sk, &ct).unwrap();
            let dual = toy_kani_dual_first_half(&sk, &ct).unwrap();
            assert!(transcript.product_steps > 0);
            assert_eq!(
                transcript.product_steps as usize + transcript.first_half.nulls.len(),
                parameters.two_exponent as usize
            );
            assert_eq!(transcript.first_half.split, None);
            assert_eq!(dual.first_half.split, None);
            let jacobian_key = |null: &crate::theta::Th<Fq>| {
                let rh = crate::theta::rosenhain(&f, null).unwrap();
                let t = (2..20u64)
                    .map(|x| f.from_u64(x))
                    .find(|&x| ![f.zero(), f.one(), rh[0], rh[1], rh[2]].contains(&x))
                    .unwrap();
                let roots = [f.zero(), f.one(), rh[0], rh[1], rh[2]].map(|x| f.inv(f.sub(x, t)));
                let roots = [roots[0], roots[1], roots[2], roots[3], roots[4], f.zero()];
                crate::genus2::ic_key(&f, &crate::genus2::igusa_clebsch(&f, f.one(), &roots))
            };
            assert_eq!(
                jacobian_key(transcript.first_half.nulls.last().unwrap()),
                jacobian_key(dual.first_half.nulls.last().unwrap())
            );
            let left = transcript.first_half.nulls.last().unwrap();
            let right = dual.first_half.nulls.last().unwrap();
            assert!(crate::theta::proj_eq(&f, left, right));
            let x = torsion_point(&f, &ct.curve, 9, 3, &mut rng).unwrap();
            let x3 = pmul(&f, &ct.curve, &x, 3);
            let y = (0..4096)
                .find_map(|_| {
                    let q = torsion_point(&f, &ct.curve, 9, 3, &mut rng)?;
                    let q3 = pmul(&f, &ct.curve, &q, 3);
                    (q3 != x3 && q3 != neg(&f, &x3)).then_some(q)
                })
                .unwrap();
            let extras = [(x, Pt::Inf), (Pt::Inf, y), (x, y)];
            let swapped = extras.map(|(a, b)| (b, a));
            let (_, forward, _, _) = toy_kani_eval(&sk, &ct, true, false, &extras).unwrap();
            let (_, backward, _, _) = toy_kani_eval(&sk, &ct, true, true, &swapped).unwrap();
            for (a, b) in forward.result.images.iter().zip(&backward.result.images) {
                assert!(crate::theta::proj_eq(&f, a, b));
            }
            let actions = toy_two_halves_action_x(&sk, &ct, x, y).unwrap();
            let direct = direct_action_x(&sk, &ct, seed, actions[0].0, [x, y]);
            assert!(actions
                .iter()
                .any(|&(c, alpha)| c == actions[0].0 && alpha == direct));
            assert_eq!(toy_trapdoor_recover_seed(&sk, &ct), Ok(seed));
            assert_eq!(toy_trapdoor_decaps(&sk, &ct), Ok(key));
            assert_eq!(toy_decaps(&sk, &ct), Ok(key));
        }
    }

    #[test]
    fn special_key_single_chain_splits_and_payload_factors() {
        let parameters = ToyParameters {
            prime: 21503,
            two_exponent: 7,
            three_exponent: 1,
        };
        let mut rng = Rng::new(0x44c0);
        let sk = toy_keygen(parameters, &mut rng).unwrap();
        let f = Zp2::new(parameters.prime);
        for seed in 0..3 {
            let (ct, key) = encaps_for_seed(&sk.public, seed).unwrap();
            let x = torsion_point(&f, &ct.curve, 3, 3, &mut rng).unwrap();
            let transcript = toy_kani_single_chain(&sk, &ct, &[(x, Pt::Inf)]).unwrap();
            assert_eq!((transcript.c, transcript.cp), (8, 1));
            assert_eq!(transcript.result.images.len(), 1);
            assert!(transcript.result.split.is_some(), "seed {seed}");
            assert_eq!(transcript.remaining_steps, 0, "seed {seed}");
            let null = transcript.result.nulls.last().unwrap();
            let actual = transcript.result.split.unwrap();
            assert_eq!(
                actual,
                (
                    crate::curve::jinv(&f, &ct.curve),
                    crate::curve::jinv(&f, &ct.curve)
                )
            );
            let chart = crate::theta::find_split_chart(&f, null, &transcript.result.images, 50_000)
                .expect("factorization chart");
            let image = chart.apply(&f, &transcript.result.images[0]).unwrap();
            let (image_a, image_b) = crate::theta::factor_product(&f, &image).unwrap();
            for (factor_null, factor_image) in
                [(chart.factors.0, image_a), (chart.factors.1, image_b)]
            {
                let curve = crate::theta::curve_from_theta1(&f, &factor_null).unwrap();
                assert_eq!(
                    crate::curve::jinv(&f, &curve),
                    crate::curve::jinv(&f, &ct.curve)
                );
                let x_coord =
                    crate::theta::x_on_theta1_curve(&f, &factor_null, &factor_image).unwrap();
                let y_coord = f.sqrt(crate::curve::rhs(&f, &curve, x_coord)).unwrap();
                let point = Pt::Aff(x_coord, y_coord);
                assert!(on_curve(&f, &curve, &point));
                assert_eq!(pmul(&f, &curve, &point, 3), Pt::Inf);
            }
            assert_eq!(
                transcript.product_steps as usize
                    + transcript.result.nulls.len()
                    + transcript.remaining_steps as usize,
                parameters.two_exponent as usize
            );
            assert!(!crate::theta::proj_eq(
                &f,
                &transcript.result.images[0],
                transcript.result.nulls.last().unwrap()
            ));
            let x1 = torsion_point(&f, &ct.curve, 3, 3, &mut rng).unwrap();
            let x2 = (0..4096)
                .find_map(|_| {
                    let q = torsion_point(&f, &ct.curve, 3, 3, &mut rng)?;
                    (q != x1 && q != neg(&f, &x1)).then_some(q)
                })
                .unwrap();
            let (c, alpha_x) = toy_single_chain_action_x(&sk, &ct, x1, x2).unwrap();
            let kernel = padd(
                &f,
                &sk.public.curve,
                &sk.public.sender_p,
                &pmul(&f, &sk.public.curve, &sk.public.sender_q, seed as u128),
            );
            let phi = velu_cyclic(&f, &sk.public.curve, &kernel, 3);
            let dual_kernel = phi.eval(&f, &sk.public.sender_q);
            let dual = velu_cyclic(&f, &ct.curve, &dual_kernel, 3);
            let u = toy_curve_iso_scale(&f, &dual.cod, &sk.public.curve).unwrap();
            let check = toy_iso_forward(&f, u, dual.eval(&f, &phi.eval(&f, &sk.public.eigen_p)));
            let want = pmul(&f, &sk.public.curve, &sk.public.eigen_p, 3);
            let dual_sign = if check == want {
                false
            } else if check == neg(&f, &want) {
                true
            } else {
                panic!("dual normalization")
            };
            let theta_iso = velu_cyclic(&f, &sk.public.curve, &sk.theta_kernel, 2);
            for (idx, point) in [x1, x2].into_iter().enumerate() {
                let mut pre = toy_iso_forward(&f, u, dual.eval(&f, &point));
                if dual_sign {
                    pre = neg(&f, &pre);
                }
                let theta_pre = theta_eval(&f, &theta_iso, sk.theta_scale, &pre);
                let theta_prime = padd(
                    &f,
                    &sk.public.curve,
                    &pmul(&f, &sk.public.curve, &theta_pre, 2),
                    &if sk.trace == 1 { neg(&f, &pre) } else { pre },
                );
                let psi = phi.eval(&f, &theta_prime);
                let alpha = padd(&f, &ct.curve, &psi, &pmul(&f, &ct.curve, &point, c as u128));
                assert_eq!(
                    match alpha {
                        Pt::Aff(x, _) => x,
                        Pt::Inf => panic!("zero alpha"),
                    },
                    alpha_x[idx]
                );
            }
            assert_eq!(toy_trapdoor_recover_seed(&sk, &ct), Ok(seed));
            assert_eq!(toy_trapdoor_decaps(&sk, &ct), Ok(key));
            let mut tampered = ct;
            tampered.s = neg(&f, &ct.s);
            assert!(toy_trapdoor_decaps(&sk, &tampered).is_err());
        }
    }

    #[test]
    fn single_chain_decoder_multiple_keys() {
        let parameters = ToyParameters {
            prime: 21503,
            two_exponent: 7,
            three_exponent: 1,
        };
        for seed in [1, 17, 0x12345] {
            let mut rng = Rng::new(seed);
            let sk = toy_keygen(parameters, &mut rng).unwrap();
            for sender in 0..3 {
                let (ct, key) = encaps_for_seed(&sk.public, sender).unwrap();
                assert_eq!(
                    toy_trapdoor_decaps(&sk, &ct),
                    Ok(key),
                    "key seed {seed}, sender {sender}"
                );
            }
        }
    }

    #[test]
    fn two_halves_oracle_decoder_multiple_keys() {
        let parameters = ToyParameters {
            prime: 16127,
            two_exponent: 6,
            three_exponent: 2,
        };
        for key_seed in [1, 17] {
            let mut rng = Rng::new(key_seed);
            let sk = toy_keygen(parameters, &mut rng).unwrap();
            for sender in 0..9 {
                let (ct, key) = encaps_for_seed(&sk.public, sender).unwrap();
                assert_eq!(
                    toy_trapdoor_decaps(&sk, &ct),
                    Ok(key),
                    "key seed {key_seed}, sender {sender}"
                );
                if sender == 0 {
                    let f = Zp2::new(parameters.prime);
                    let mut tampered = ct;
                    tampered.s = neg(&f, &ct.s);
                    assert!(toy_trapdoor_decaps(&sk, &tampered).is_err());
                }
            }
        }
    }

    #[test]
    #[ignore = "bounded exploratory two-halves odd-torsion oracle; not a scalable decoder"]
    fn two_halves_middle_automorphism_oracle() {
        use std::collections::{HashSet, VecDeque};
        let parameters = ToyParameters {
            prime: 16127,
            two_exponent: 6,
            three_exponent: 2,
        };
        let mut rng = Rng::new(0x44c0);
        let sk = toy_keygen(parameters, &mut rng).unwrap();
        let (ct, _) = encaps_for_seed(&sk.public, 0).unwrap();
        let f = Zp2::new(parameters.prime);
        let p = torsion_point(&f, &ct.curve, 9, 3, &mut rng).unwrap();
        let p3 = pmul(&f, &ct.curve, &p, 3);
        let q = (0..4096)
            .find_map(|_| {
                let x = torsion_point(&f, &ct.curve, 9, 3, &mut rng)?;
                let x3 = pmul(&f, &ct.curve, &x, 3);
                (x3 != p3 && x3 != neg(&f, &p3)).then_some(x)
            })
            .unwrap();
        let points: Vec<_> = (0..9u64)
            .flat_map(|i| (0..9u64).map(move |j| (i, j)))
            .map(|(i, j)| {
                padd(
                    &f,
                    &ct.curve,
                    &pmul(&f, &ct.curve, &p, i as u128),
                    &pmul(&f, &ct.curve, &q, j as u128),
                )
            })
            .collect();
        let extras: Vec<_> = points
            .iter()
            .flat_map(|&x| points.iter().map(move |&y| (x, y)))
            .collect();
        let (_, first, c, cp) = toy_kani_eval(&sk, &ct, true, false, &extras).unwrap();
        assert_eq!(first.result.images.len(), 81 * 81);
        let null = *first.result.nulls.last().unwrap();
        let normalize = |v: crate::theta::Th<Fq>| {
            let k = v.iter().position(|&x| !f.is_zero(x)).unwrap();
            let inv = f.inv(v[k]);
            v.map(|x| f.mul(x, inv))
        };
        let mut table = HashMap::new();
        for (&v, &source) in first.result.images.iter().zip(&extras) {
            table.insert(normalize(v), source);
        }
        let imaginary = f.sqrt(f.neg(f.one())).unwrap();
        let apply = |path: &[u8], v: crate::theta::Th<Fq>| {
            path.iter().fold(v, |w, &g| {
                crate::theta::split_chart_step(&f, imaginary, g, &w)
            })
        };
        let identity =
            [0, 1, 2, 3].map(|j| [0, 1, 2, 3].map(|i| if i == j { f.one() } else { f.zero() }));
        let normalize_matrix = |m: [crate::theta::Th<Fq>; 4]| {
            let z = m
                .iter()
                .flatten()
                .copied()
                .find(|&x| !f.is_zero(x))
                .unwrap();
            let inv = f.inv(z);
            m.map(|col| col.map(|x| f.mul(x, inv)))
        };
        let sample = [points[9], points[1]];
        let sample_indices = [(9usize * 81), 9usize, 81usize, 1usize];
        let expected_cp = sample.map(|x| pmul(&f, &ct.curve, &x, cp as u128));
        let xcoord = |x: Pt<Fq>| match x {
            Pt::Aff(x, _) => x,
            Pt::Inf => f.zero(),
        };
        let mut seen = HashSet::new();
        let mut todo = VecDeque::from([(identity, Vec::<u8>::new())]);
        let mut candidates = Vec::new();
        while let Some((matrix, path)) = todo.pop_front() {
            let matrix = normalize_matrix(matrix);
            if !seen.insert(matrix) {
                continue;
            }
            if crate::theta::proj_eq(&f, &apply(&path, null), &null) {
                let mapped: Option<Vec<_>> = sample_indices
                    .iter()
                    .map(|&idx| {
                        table
                            .get(&normalize(apply(&path, first.result.images[idx])))
                            .copied()
                    })
                    .collect();
                if let Some(mapped) = mapped {
                    let outputs = mapped
                        .iter()
                        .map(|&(a, b)| (pmul(&f, &ct.curve, &b, 64), pmul(&f, &ct.curve, &a, 64)))
                        .collect::<Vec<_>>();
                    if xcoord(outputs[0].1) == xcoord(expected_cp[0])
                        && xcoord(outputs[1].0) == xcoord(expected_cp[0])
                        && xcoord(outputs[2].1) == xcoord(expected_cp[1])
                        && xcoord(outputs[3].0) == xcoord(expected_cp[1])
                    {
                        candidates.push((path.clone(), outputs));
                    }
                }
            }
            if seen.len() >= 50_000 {
                break;
            }
            for g in 0..6 {
                let mut next = path.clone();
                next.push(g);
                todo.push_back((
                    matrix.map(|v| crate::theta::split_chart_step(&f, imaginary, g, &v)),
                    next,
                ));
            }
        }
        eprintln!(
            "middle oracle: states={}, table={}, candidates={}, c={}, cp={}",
            seen.len(),
            table.len(),
            candidates.len(),
            c,
            cp
        );
        for (path, outputs) in &candidates {
            eprintln!(
                "candidate path length {} alpha x {:?} {:?}",
                path.len(),
                xcoord(outputs[0].0),
                xcoord(outputs[2].0)
            );
        }
        assert!(!candidates.is_empty());
        let sender_kernel = sk.public.sender_p;
        let phi = velu_cyclic(&f, &sk.public.curve, &sender_kernel, 9);
        let dual_kernel = phi.eval(&f, &sk.public.sender_q);
        let dual = velu_cyclic(&f, &ct.curve, &dual_kernel, 9);
        let u = toy_curve_iso_scale(&f, &dual.cod, &sk.public.curve).unwrap();
        let check = toy_iso_forward(&f, u, dual.eval(&f, &phi.eval(&f, &sk.public.eigen_p)));
        let want = pmul(&f, &sk.public.curve, &sk.public.eigen_p, 9);
        let dual_sign = if check == want {
            false
        } else if check == neg(&f, &want) {
            true
        } else {
            panic!("dual normalization")
        };
        let theta_iso = velu_cyclic(&f, &sk.public.curve, &sk.theta_kernel, 2);
        let direct: Vec<_> = sample
            .iter()
            .map(|&point| {
                let mut pre = toy_iso_forward(&f, u, dual.eval(&f, &point));
                if dual_sign {
                    pre = neg(&f, &pre);
                }
                let theta_pre = theta_eval(&f, &theta_iso, sk.theta_scale, &pre);
                let theta_prime = padd(
                    &f,
                    &sk.public.curve,
                    &pmul(&f, &sk.public.curve, &theta_pre, 2),
                    &if sk.trace == 1 { neg(&f, &pre) } else { pre },
                );
                padd(
                    &f,
                    &ct.curve,
                    &phi.eval(&f, &theta_prime),
                    &pmul(&f, &ct.curve, &point, c as u128),
                )
            })
            .collect();
        let matching: Vec<_> = candidates
            .iter()
            .filter(|(_, outputs)| {
                xcoord(outputs[0].0) == xcoord(direct[0])
                    && xcoord(outputs[2].0) == xcoord(direct[1])
            })
            .collect();
        eprintln!(
            "direct alpha x {:?} {:?}; matching candidates={}",
            xcoord(direct[0]),
            xcoord(direct[1]),
            matching.len()
        );
        assert!(!matching.is_empty());
    }
}

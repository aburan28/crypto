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
use crate::theta::{chain, ChainResult};
use std::collections::HashMap;

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
    let target = (n * n)
        .checked_sub(d * d * 7)
        .ok_or(ToyError::InvalidParameters)?;
    let (c, cp) = (2..n)
        .step_by(2)
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
    let graph = [
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
        if curves[0] == curves[1]
            && g0.0 != Pt::Inf
            && g1.0 != Pt::Inf
            && g0.0 != g1.0
            && g0.0 == g0.1
            && g1.0 == g1.1
        {
            k = k.map(|(x, y)| {
                (
                    padd(&f, &curves[0], &x, &y),
                    padd(&f, &curves[0], &x, &neg(&f, &y)),
                )
            });
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
            curves = [iso_left.cod, iso_right.cod];
        }
        remaining -= 1;
        product_steps += 1;
    }
    let res =
        chain(&f, &curves[0], &curves[1], k, remaining, &[]).ok_or(ToyError::ConstructionFailed)?;
    Ok(ToyKaniTranscript {
        product_steps,
        first_half: res,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

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
        for seed in 0..9 {
            let (ct, key) = encaps_for_seed(&sk.public, seed).unwrap();
            let transcript = toy_kani_first_half(&sk, &ct).unwrap();
            assert!(transcript.product_steps > 0);
            assert_eq!(
                transcript.product_steps as usize + transcript.first_half.nulls.len(),
                parameters.two_exponent as usize
            );
            assert_eq!(transcript.first_half.split, None);
            assert_eq!(toy_decaps(&sk, &ct), Ok(key));
        }
    }
}

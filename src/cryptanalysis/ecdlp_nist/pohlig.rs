//! Pohlig–Hellman on an [`EcdlpGroup`] whose declared order is composite.
//!
//! The order `N` of the base point is factored (trial division to `2^16`,
//! then Brent's rho — `curve_traits::arith::factor`); for each `q^e` the
//! problem is projected onto the `q^e`-subgroup and solved digit by digit
//! in the `q`-subgroup, with [`interval::bsgs`] for `q ≤ 2^40` and the
//! generic rho above that; the digits are combined by CRT.
//!
//! On every NIST curve the generator's order `n` is prime, so this path is
//! the identity: the audit records that, and the dispatcher only takes this
//! branch for toy full-group fixtures or user-supplied composite orders.

use std::time::Instant;

use num_bigint::BigUint;
use num_traits::{One, Zero};

use super::group::EcdlpGroup;
use super::interval;
use super::rho::{pollard_rho, RhoOptions};
use crate::cryptanalysis::curve_traits::arith::factor;

/// What a Pohlig–Hellman run did.
#[derive(Clone, Debug)]
pub struct PohligHellmanReport {
    pub scalar: Option<BigUint>,
    /// `(q, e)` for `N = ∏ q^e`.
    pub factors: Vec<(BigUint, u32)>,
    /// `(q^e, k mod q^e)` recovered.
    pub residues: Vec<(BigUint, BigUint)>,
    /// Composite cofactors the factoring budget could not split.
    pub unfactored: Vec<BigUint>,
    pub group_ops: u64,
    pub elapsed_ms: u128,
    pub failure: Option<String>,
}

/// A view of `g` whose generator is `base` of order `order`, with
/// automorphism folding switched off (the eigenvalue is only known modulo
/// the full subgroup order).
pub struct SubgroupView<'a, G: EcdlpGroup> {
    pub inner: &'a G,
    pub base: G::Elt,
    pub order: BigUint,
}

impl<G: EcdlpGroup> EcdlpGroup for SubgroupView<'_, G> {
    type Elt = G::Elt;
    fn name(&self) -> &str {
        self.inner.name()
    }
    fn order(&self) -> &BigUint {
        &self.order
    }
    fn generator(&self) -> G::Elt {
        self.base.clone()
    }
    fn identity(&self) -> G::Elt {
        self.inner.identity()
    }
    fn is_identity(&self, p: &G::Elt) -> bool {
        self.inner.is_identity(p)
    }
    fn add(&self, p: &G::Elt, q: &G::Elt) -> G::Elt {
        self.inner.add(p, q)
    }
    fn double(&self, p: &G::Elt) -> G::Elt {
        self.inner.double(p)
    }
    fn neg(&self, p: &G::Elt) -> G::Elt {
        self.inner.neg(p)
    }
    fn mul(&self, p: &G::Elt, k: &BigUint) -> G::Elt {
        self.inner.mul(p, k)
    }
    fn key(&self, p: &G::Elt) -> u64 {
        self.inner.key(p)
    }
    fn canonical_sign(&self, p: &G::Elt) -> (G::Elt, bool) {
        self.inner.canonical_sign(p)
    }
    fn coords_hex(&self, p: &G::Elt) -> Option<(String, String)> {
        self.inner.coords_hex(p)
    }
}

/// Solve a DLP in a subgroup of prime order `q` by BSGS or rho.
fn prime_order_dlp<G: EcdlpGroup>(
    view: &SubgroupView<'_, G>,
    target: &G::Elt,
    rho_opts: &RhoOptions,
    ops: &mut u64,
) -> Option<BigUint> {
    if view.is_identity(target) {
        return Some(BigUint::zero());
    }
    if view.order.bits() <= 40 {
        let rep = interval::bsgs(view, target, &BigUint::zero(), &view.order);
        *ops += rep.group_ops;
        rep.scalar
    } else {
        let rep = pollard_rho(view, target, rho_opts);
        *ops += rep.iterations;
        rep.scalar
    }
}

/// CRT for pairwise coprime moduli.
pub fn crt(residues: &[(BigUint, BigUint)]) -> Option<BigUint> {
    let mut x = BigUint::zero();
    let mut m = BigUint::one();
    for (modulus, r) in residues {
        // x' ≡ x (mod m), x' ≡ r (mod modulus).
        let m_inv = m.modpow(&euler_phi_minus_one(modulus)?, modulus);
        let diff = (r + modulus - (&x % modulus)) % modulus;
        let t = diff * m_inv % modulus;
        x += &m * t;
        m *= modulus;
    }
    Some(x)
}

/// Exponent for inversion mod a prime power `q^e`: `φ(q^e) − 1`.
fn euler_phi_minus_one(modulus: &BigUint) -> Option<BigUint> {
    // Caller passes q^e with q prime; recover q by factoring (tiny input).
    let f = factor(modulus, 10_000);
    if !f.complete() || f.primes.len() != 1 {
        return None;
    }
    let (q, e) = &f.primes[0];
    let phi = q.pow(e - 1) * (q - BigUint::one());
    Some(phi - BigUint::one())
}

/// Solve `Q = [k]G` where `G = g.generator()` has (possibly composite)
/// order `g.order()`.
pub fn pohlig_hellman<G: EcdlpGroup>(
    g: &G,
    target: &G::Elt,
    rho_opts: &RhoOptions,
) -> PohligHellmanReport {
    let t0 = Instant::now();
    let n = g.order().clone();
    let gen = g.generator();
    let f = factor(&n, 1 << 22);
    let mut report = PohligHellmanReport {
        scalar: None,
        factors: f.primes.clone(),
        residues: Vec::new(),
        unfactored: f.composites.iter().map(|(c, _)| c.clone()).collect(),
        group_ops: 0,
        elapsed_ms: 0,
        failure: None,
    };
    if !f.complete() {
        report.failure = Some("order has an unfactored composite part".into());
        report.elapsed_ms = t0.elapsed().as_millis();
        return report;
    }
    for (q, e) in &f.primes {
        let qe = q.pow(*e);
        let cof = &n / &qe;
        // Project onto the q^e-subgroup.
        let g_i = g.mul(&gen, &cof);
        let q_i = g.mul(target, &cof);
        report.group_ops += 4 * cof.bits();
        // Base of order q.
        let g_q = g.mul(&g_i, &q.pow(e - 1));
        let view = SubgroupView {
            inner: g,
            base: g_q,
            order: q.clone(),
        };
        let mut k_i = BigUint::zero();
        let mut q_pow = BigUint::one();
        for j in 0..*e {
            // h_j = [q^(e−1−j)] (Q_i − [k_i] G_i) lies in the q-subgroup.
            let partial = g.mul(&g_i, &k_i);
            let diff = g.add(&q_i, &g.neg(&partial));
            let h_j = g.mul(&diff, &q.pow(e - 1 - j));
            report.group_ops += 4 * q.bits() * (*e as u64);
            let Some(d_j) = prime_order_dlp(&view, &h_j, rho_opts, &mut report.group_ops) else {
                report.failure = Some(format!("digit {j} of the {q}-adic expansion not found"));
                report.elapsed_ms = t0.elapsed().as_millis();
                return report;
            };
            k_i += &q_pow * d_j;
            q_pow *= q;
        }
        report.residues.push((qe, k_i));
    }
    report.scalar = crt(&report.residues).filter(|k| g.mul(&gen, k) == *target);
    if report.scalar.is_none() {
        report.failure = Some("CRT combination failed verification".into());
    }
    report.elapsed_ms = t0.elapsed().as_millis();
    report
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::ecdlp_nist::toy;

    #[test]
    fn pohlig_hellman_on_toy_full_groups() {
        for (m, a) in [(17u32, 0u8), (19, 0)] {
            let g = toy::koblitz_full_group(m, a).expect("toy");
            let n = g.order().clone();
            for k in [1u32, 2, 3, 12_345, 100_000] {
                let kb = BigUint::from(k) % &n;
                let q = g.mul(&g.generator(), &kb);
                let rep = pohlig_hellman(&g, &q, &RhoOptions::default());
                assert_eq!(rep.scalar, Some(kb), "m={m} k={k}: {rep:?}");
                assert!(rep.factors.len() >= 2);
            }
        }
    }

    #[test]
    fn crt_recombines() {
        let r = vec![
            (BigUint::from(4u32), BigUint::from(3u32)),
            (BigUint::from(9u32), BigUint::from(5u32)),
            (BigUint::from(25u32), BigUint::from(7u32)),
        ];
        let x = crt(&r).unwrap();
        assert_eq!(&x % 4u32, BigUint::from(3u32));
        assert_eq!(&x % 9u32, BigUint::from(5u32));
        assert_eq!(&x % 25u32, BigUint::from(7u32));
    }
}

//! # Challenge 59 — Invalid-Curve Attack on ECDH
//!
//! Standard ECDH addition formulas use only the `a` coefficient of
//! the short-Weierstrass curve `y² = x³ + a x + b` — `b` is implicit
//! in any valid point `(x, y)`.  An attacker who skips that
//! validation can feed Alice a point from a *different* curve
//! (same `a`, different `b`) and ride its small subgroups to recover
//! Alice's secret key.
//!
//! The cryptopals curve has near-prime order — only a `2³` cofactor —
//! so on-curve subgroup confinement leaks at most 3 bits.  But three
//! attacker-chosen sister curves with rich small-factor structure
//! cover enough residues that the CRT reassembly lands on the full
//! key.
//!
//! ## The recipe
//!
//! 1. Pick a sister curve `E_b' : y² = x³ + a x + b'` whose order
//!    has a small factor `r` not present in the cryptopals curve.
//! 2. Generate a point `P` of order `r` on `E_b'`:
//!    - Pick random `x` in GF(p).
//!    - Compute `y = sqrt(x³ + a x + b') mod p` via Tonelli–Shanks.
//!    - `P := (order_of_E_b' / r) · (x, y)` (using addition formulas
//!      that ignore `b`).  Restart if `P == O`.
//! 3. Send `P` to Alice.  She'll compute `K = d · P` *on her own
//!    formulas*, but mathematically `K` lives in the order-`r`
//!    subgroup of `E_b'`.
//! 4. Brute-force `b ∈ [0, r)` by checking
//!    `HMAC-SHA256(b · P, msg) == tag` to recover `d mod r`.
//! 5. Repeat with more `r` values.  CRT-combine until the product
//!    exceeds the order of Alice's curve — and you have `d`.

use crate::cryptopals::set8_util::{
    biguint_to_bytes_be, crt_combine, hmac_sha256, parse_big,
};
use crate::cryptopals::Report;
use crate::cryptanalysis::ec_index_calculus::sqrt_mod_p;
use num_bigint::{BigInt, BigUint, ToBigInt};
use num_traits::{One, Zero};
use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;

/// Convert BigUint mod p (operates on signed BigInt internally to
/// keep negatives sane).
fn sub_mod(a: &BigUint, b: &BigUint, p: &BigUint) -> BigUint {
    let a_i = a.to_bigint().unwrap();
    let b_i = b.to_bigint().unwrap();
    let p_i = p.to_bigint().unwrap();
    let r = ((&a_i - &b_i) % &p_i + &p_i) % &p_i;
    r.to_biguint().unwrap()
}

fn mod_inv(a: &BigUint, p: &BigUint) -> BigUint {
    use num_integer::Integer;
    let a_i = a.to_bigint().unwrap();
    let p_i = p.to_bigint().unwrap();
    let g = a_i.extended_gcd(&p_i);
    let r = ((g.x % &p_i) + &p_i) % &p_i;
    r.to_biguint().unwrap()
}

/// Affine `(x, y)` or the point at infinity.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Pt {
    Inf,
    Aff(BigUint, BigUint),
}

/// Short-Weierstrass curve `y² = x³ + a x + b` over GF(p).
///
/// Note: addition is computed *without* reference to `b`, which is
/// the very property the invalid-curve attack exploits.
#[derive(Clone, Debug)]
pub struct Curve {
    pub p: BigUint,
    pub a: BigUint,
    pub b: BigUint,
}

impl Curve {
    pub fn add(&self, p1: &Pt, p2: &Pt) -> Pt {
        match (p1, p2) {
            (Pt::Inf, _) => p2.clone(),
            (_, Pt::Inf) => p1.clone(),
            (Pt::Aff(x1, y1), Pt::Aff(x2, y2)) => {
                if x1 == x2 && y1 != y2 {
                    return Pt::Inf;
                }
                let p = &self.p;
                let m = if x1 == x2 && y1 == y2 {
                    // doubling: m = (3 x² + a) / 2y
                    let num = ((BigUint::from(3u32) * x1 * x1) + &self.a) % p;
                    let den = mod_inv(&((y1 + y1) % p), p);
                    (num * den) % p
                } else {
                    // chord: m = (y2 - y1) / (x2 - x1)
                    let num = sub_mod(y2, y1, p);
                    let den = mod_inv(&sub_mod(x2, x1, p), p);
                    (num * den) % p
                };
                let x3 = sub_mod(&sub_mod(&(&m * &m % p), x1, p), x2, p);
                let y3 = sub_mod(&((&m * sub_mod(x1, &x3, p)) % p), y1, p);
                Pt::Aff(x3, y3)
            }
        }
    }

    pub fn scalar_mul(&self, k: &BigUint, q: &Pt) -> Pt {
        let mut result = Pt::Inf;
        let mut base = q.clone();
        let mut k = k.clone();
        while !k.is_zero() {
            if (&k & BigUint::one()).is_one() {
                result = self.add(&result, &base);
            }
            base = self.add(&base, &base);
            k >>= 1;
        }
        result
    }
}

/// Cryptopals' main curve and its sisters.
pub fn cryptopals_curve() -> (Curve, BigUint, Pt) {
    let p = parse_big("233970423115425145524320034830162017933");
    let a = sub_mod(&BigUint::zero(), &BigUint::from(95051u32), &p);
    let b = BigUint::from(11279326u32);
    let curve = Curve { p, a, b };
    let order = parse_big("29246302889428143187362802287225875743");
    let gx = BigUint::from(182u32);
    let gy = parse_big("85518893674295321206118380980485522083");
    (curve, order, Pt::Aff(gx, gy))
}

/// Three sister curves cryptopals supplies with their orders.
pub fn sister_curves() -> Vec<(BigUint, BigUint)> {
    // (b, order)
    vec![
        (
            BigUint::from(210u32),
            parse_big("233970423115425145550826547352470124412"),
        ),
        (
            BigUint::from(504u32),
            parse_big("233970423115425145544350131142039591210"),
        ),
        (
            BigUint::from(727u32),
            parse_big("233970423115425145545378039958152057148"),
        ),
    ]
}

/// Alice plays victim.  Holds secret `d`, accepts any peer point,
/// computes `K = d · P` *with no curve check*, and returns
/// `HMAC-SHA256(K_x_bytes, msg)`.
pub struct Alice {
    pub curve: Curve,
    pub d: BigUint,
}

impl Alice {
    pub fn oracle(&self, p: &Pt) -> ([u8; 32], &'static [u8]) {
        let k = self.curve.scalar_mul(&self.d, p);
        let msg: &[u8] = b"crazy flamboyant for the rap enjoyment";
        // Use BOTH x and y in the MAC key.  In a pure x-coordinate
        // shared-secret scheme (Montgomery ladder) we'd lose the sign
        // and end up with `±(d mod r)` ambiguity per residue, but
        // for short-Weierstrass ECDH the y is observable and
        // resolves the ambiguity.
        let key_bytes = match k {
            Pt::Inf => vec![0u8; 32],
            Pt::Aff(x, y) => {
                let mut b = biguint_to_bytes_be(&x, 16);
                b.extend_from_slice(&biguint_to_bytes_be(&y, 16));
                b
            }
        };
        (hmac_sha256(&key_bytes, msg), msg)
    }
}

fn rhs(curve: &Curve, x: &BigUint) -> BigUint {
    let p = &curve.p;
    let x2 = (x * x) % p;
    let x3 = (&x2 * x) % p;
    let ax = (&curve.a * x) % p;
    ((x3 + ax) + &curve.b) % p
}

/// Find a point of order *exactly* `r` on `curve` (which has order
/// `curve_order`).  We pick random `x`, compute `y` via
/// Tonelli-Shanks, then multiply by `curve_order / r`.  Restart if
/// the result is the identity.
pub fn random_point_of_order(
    curve: &Curve,
    curve_order: &BigUint,
    r: &BigUint,
    seed: u64,
) -> Pt {
    let mut rng = StdRng::seed_from_u64(seed);
    let cofactor = curve_order / r;
    let p = &curve.p;
    loop {
        let x = rng.gen::<u128>() % 1_000_000;
        let x_b = BigUint::from(x) % p;
        let rhs_val = rhs(curve, &x_b);
        let y = match sqrt_mod_p(&rhs_val, p) {
            Some(y) => y,
            None => continue,
        };
        let pt = Pt::Aff(x_b, y);
        let scaled = curve.scalar_mul(&cofactor, &pt);
        if scaled != Pt::Inf {
            return scaled;
        }
    }
}

/// Brute-force `b ∈ [0, r)` to recover `d mod r` from one MAC.
pub fn recover_residue_ec(
    curve: &Curve,
    pt: &Pt,
    r: &BigUint,
    msg: &[u8],
    tag: &[u8; 32],
) -> Option<BigUint> {
    let mut b = BigUint::zero();
    let mut cur = Pt::Inf;
    while &b < r {
        let key_bytes = match &cur {
            Pt::Inf => vec![0u8; 32],
            Pt::Aff(x, y) => {
                let mut bb = biguint_to_bytes_be(x, 16);
                bb.extend_from_slice(&biguint_to_bytes_be(y, 16));
                bb
            }
        };
        let got = hmac_sha256(&key_bytes, msg);
        if &got == tag {
            return Some(b);
        }
        cur = curve.add(&cur, pt);
        b += 1u32;
    }
    None
}

/// Collect residues across all three sister curves until the
/// product covers the main curve's order, then CRT-combine.
pub fn attack(alice: &Alice, main_curve: &Curve, main_order: &BigUint) -> Option<BigUint> {
    use crate::cryptopals::set8_util::small_factors;
    let mut residues: Vec<(BigUint, BigUint)> = Vec::new();
    let mut product = BigUint::one();
    let mut seed = 31_u64;
    'outer: for (b, ord) in sister_curves() {
        if &product >= main_order {
            break;
        }
        let curve = Curve {
            p: main_curve.p.clone(),
            a: main_curve.a.clone(),
            b,
        };
        // Small factors of this sister's order.
        let factors = small_factors(&ord, 1 << 16);
        for r in factors {
            if &product >= main_order {
                break 'outer;
            }
            if r <= BigUint::from(2u32) {
                continue;
            }
            if r > BigUint::from(1u64 << 17) {
                continue;
            }
            // Avoid duplicates: skip if a previous residue already used `r`.
            if residues.iter().any(|(_, rr)| rr == &r) {
                continue;
            }
            let pt = random_point_of_order(&curve, &ord, &r, seed);
            seed = seed.wrapping_add(1);
            let (tag, msg) = alice.oracle(&pt);
            let resi = match recover_residue_ec(&curve, &pt, &r, msg, &tag) {
                Some(b) => b,
                None => continue,
            };
            residues.push((resi, r.clone()));
            product *= r;
        }
    }
    if &product < main_order {
        return None;
    }
    let (a, _m) = crt_combine(&residues);
    let _ = BigInt::from(0);
    Some(a % main_order)
}

pub fn run() -> Report {
    let mut r = Report::new(59, "Invalid-Curve Attack on ECDH");
    let (curve, order, _g) = cryptopals_curve();
    let secret_d = parse_big("12345678901234567890123456789012345");
    let alice = Alice {
        curve: curve.clone(),
        d: secret_d.clone(),
    };
    r.line(format!("Alice's secret d : {}", secret_d));
    r.line(format!("curve order n    : {}", order));
    match attack(&alice, &curve, &order) {
        Some(d) => {
            r.line(format!("Recovered     d  : {}", d));
            r.line(format!("Match            : {}", d == secret_d));
            if d == secret_d {
                return r.succeed();
            }
        }
        None => r.line("Attack did not converge (not enough small factors)."),
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn generator_lands_on_infinity_at_order() {
        let (curve, order, g) = cryptopals_curve();
        let r = curve.scalar_mul(&order, &g);
        assert_eq!(r, Pt::Inf);
    }

    #[test]
    fn invalid_curve_attack_recovers_key() {
        let (curve, order, _g) = cryptopals_curve();
        let secret = parse_big("987654321987654321987654321");
        let alice = Alice {
            curve: curve.clone(),
            d: secret.clone(),
        };
        let d = attack(&alice, &curve, &order).expect("attack must succeed");
        assert_eq!(d, secret);
    }
}

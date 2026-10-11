//! Weil pairing of order n = 3^k on y^2 = x^3 + a2 x^2 + a4 x + a6 and a
//! constant-time discrete logarithm in mu_{3^k}.  The Miller loop runs over
//! the public bits of n; line functions are kept with their exact scale
//! factors (numerator/denominator accumulation, one inversion at the end),
//! which the Weil-pairing ratio needs because the embedding degree is 1.

use crate::ct::ConditionallySelectable;
use crate::field::Field;
use crate::scalar::Scalar;
use crate::weier::{WCurve, WPoint};
use num_bigint::BigUint;

/// f_{n,P}(Q) as (num, den), for odd n with P of exact order n.
fn miller<F: Field>(e: &WCurve<F>, p: &WPoint<F>, q: &WPoint<F>, n: &BigUint) -> (F, F) {
    let bits: Vec<bool> = (0..n.bits()).rev().map(|i| n.bit(i)).collect();
    let (xp, yp) = (p.x, p.y);
    let (xq, yq) = (q.x, q.y);
    let (mut x, mut y, mut z) = (xp, yp, F::one());
    let mut num = F::one();
    let mut den = F::one();
    let last = bits.len() - 1;
    for (idx, &bit) in bits.iter().enumerate().skip(1) {
        // doubling step
        let z2 = z.square();
        let y2 = y.square();
        let m = x.square().mul_small(3) + (e.a2 * x * z2).double() + e.a4 * z2.square();
        let l = (y * z).double() * (yq * z2 * z - y) - m * z * (xq * z2 - x);
        let sl = (y * z2.square()).double();
        let z3 = (y * z).double();
        let x3 = m.square() - e.a2 * z3.square() - (x * y2).mul_small(8);
        let y3 = m * ((x * y2).mul_small(4) - x3) - y2.square().mul_small(8);
        let z3s = z3.square();
        let v = xq * z3s - x3;
        num = num.square() * l * z3s;
        den = den.square() * sl * v;
        x = x3;
        y = y3;
        z = z3;
        if bit {
            if idx == last {
                // [n-1]P + P = O: vertical line through P, no vertical denominator
                num *= xq - xp;
            } else {
                let z2 = z.square();
                let h = xp * z2 - x;
                let r = yp * z2 * z - y;
                let l = (yq - yp) * h * z - r * (xq - xp);
                let sl = h * z;
                let z3 = h * z;
                let h2 = h.square();
                let x3 = r.square() - e.a2 * z3.square() - (x + xp * z2) * h2;
                let y3 = r * (x * h2 - x3) - y * h2 * h;
                let z3s = z3.square();
                let v = xq * z3s - x3;
                num = num * l * z3s;
                den = den * sl * v;
                x = x3;
                y = y3;
                z = z3;
            }
        }
    }
    (num, den)
}

/// e_n(P, Q)^2 (the square removes the (-1)^n sign and is a bijection on
/// mu_n for odd n).  Constant time in P and Q for fixed public n.
pub fn weil_sq<F: Field>(e: &WCurve<F>, p: &WPoint<F>, q: &WPoint<F>, n: &BigUint) -> F {
    let (n1, d1) = miller(e, p, q, n);
    let (n2, d2) = miller(e, q, p, n);
    let w = (n1 * d2) * (d1 * n2).invert();
    w.square()
}

/// x^(3^j).
pub fn cube_iter<F: Field>(x: &F, j: usize) -> F {
    let mut r = *x;
    for _ in 0..j {
        r = r.square() * r;
    }
    r
}

/// Is zeta (in mu_{3^k}) of full order 3^k?  Constant time.
pub fn full_order_3k<F: Field>(zeta: &F, k: usize) -> crate::ct::Choice {
    !cube_iter(zeta, k - 1).ct_eq(&F::one())
}

/// s with base^s = target in mu_{3^k} (base of full order).  Constant time:
/// fixed number of cubings, digits chosen by constant-time comparison.
/// Elements of mu_{3^k} have norm 1, so x^-1 = x^p = frob(x).
pub fn dlog_3k<F: Field>(base: &F, target: &F, k: usize) -> Scalar {
    let gamma = cube_iter(base, k - 1);
    let gamma2 = gamma.square();
    let mut pw = Vec::with_capacity(k);
    let mut cur = *base;
    for _ in 0..k {
        pw.push(cur);
        cur = cur.square() * cur;
    }
    let mut acc = F::one();
    let mut digits = vec![0u64; k];
    for i in 0..k {
        let h = cube_iter(&(*target * acc), k - 1 - i);
        let d1 = h.ct_eq(&gamma);
        let d2 = h.ct_eq(&gamma2);
        let inv = pw[i].frob();
        let t1 = F::conditional_select(&F::one(), &inv, d1);
        let t2 = F::conditional_select(&F::one(), &inv.square(), d2);
        acc = acc * t1 * t2;
        digits[i] = u64::conditional_select(&0, &1, d1) | u64::conditional_select(&0, &2, d2);
    }
    // s = sum d_i 3^i by Horner (fixed width, constant time)
    let limbs = (k * 16 / 10) / 64 + 2; // 3^k < 2^(1.585 k)
    let mut s = vec![0u64; limbs];
    for i in (0..k).rev() {
        let mut carry = digits[i] as u128;
        for w in s.iter_mut() {
            let t = (*w as u128) * 3 + carry;
            *w = t as u64;
            carry = t >> 64;
        }
    }
    Scalar { l: s }
}

//! Characteristic-2 division polynomials and a floor test for l-volcanoes.
//!
//! For y^2 + xy = x^3 + a2 x^2 + a6 (IEEE P1363 A.12.?): f_0 = 0, f_1 = 1,
//! f_2 = x, f_3 = x^4 + x^3 + a6, f_4 = x^6 + a6 x^2,
//! f_{2n+1} = f_{n+2} f_n^3 + f_{n-1} f_{n+1}^3,
//! f_{2n}   = (f_{n+2} f_{n-1}^2 + f_{n-2} f_{n+1}^2) f_n / x.
//!
//! Floor test. Let l | f be odd, a = t/2 (mod l) and k the order of +-a in
//! (Z/l)^*/{+-1}. Above the floor pi acts on E[l] as the scalar a, so all
//! (l^2 - 1)/2 x-coordinates of E[l] lie in F_{q^k}; on the floor pi - a is a
//! non-zero nilpotent on E[l], so only the eigenline's (l - 1)/2 do. Hence
//! deg gcd(f_l, x^(q^k) - x) is (l^2 - 1)/2 above the floor and (l - 1)/2 on it.

use crate::f2m::F2m;
use crate::iso3::{add, deg, divrem, gcd, mul, trim, Poly};
use std::collections::HashMap;

fn shift_div_x(p: Poly) -> Poly {
    assert_eq!(p[0], 0, "not divisible by x");
    trim(p[1..].to_vec())
}

fn fpoly(k: &F2m, a6: u64, n: usize, memo: &mut HashMap<usize, Poly>) -> Poly {
    if let Some(p) = memo.get(&n) {
        return p.clone();
    }
    let p = match n {
        0 => vec![0],
        1 => vec![1],
        2 => vec![0, 1],
        3 => vec![a6, 0, 0, 1, 1],
        4 => vec![0, 0, a6, 0, 0, 0, 1],
        _ if n % 2 == 1 => {
            let j = (n - 1) / 2;
            let fj = fpoly(k, a6, j, memo);
            let fj1 = fpoly(k, a6, j + 1, memo);
            let fj2 = fpoly(k, a6, j + 2, memo);
            let fjm = fpoly(k, a6, j - 1, memo);
            let a = mul(k, &fj2, &mul(k, &fj, &mul(k, &fj, &fj)));
            let b = mul(k, &fjm, &mul(k, &fj1, &mul(k, &fj1, &fj1)));
            add(&a, &b)
        }
        _ => {
            let j = n / 2;
            let fj = fpoly(k, a6, j, memo);
            let fj1 = fpoly(k, a6, j + 1, memo);
            let fj2 = fpoly(k, a6, j + 2, memo);
            let fjm = fpoly(k, a6, j - 1, memo);
            let fjm2 = fpoly(k, a6, j - 2, memo);
            let a = mul(k, &fj2, &mul(k, &fjm, &fjm));
            let b = mul(k, &fjm2, &mul(k, &fj1, &fj1));
            shift_div_x(mul(k, &add(&a, &b), &fj))
        }
    };
    memo.insert(n, p.clone());
    p
}

pub fn division_poly(k: &F2m, a6: u64, n: usize) -> Poly {
    fpoly(k, a6, n, &mut HashMap::new())
}

fn eval(k: &F2m, p: &Poly, x: u64) -> u64 {
    let mut r = 0u64;
    for &c in p.iter().rev() {
        r = k.mul(r, x) ^ c;
    }
    r
}

/// deg gcd(f_l, x^(q^kdeg) - x).
pub fn fixed_x_count(k: &F2m, a6: u64, ell: usize, kdeg: u32) -> usize {
    let psi = division_poly(k, a6, ell);
    assert_eq!(deg(&psi) as usize, (ell * ell - 1) / 2);
    let mut xp: Poly = vec![0, 1];
    for _ in 0..(k.m * kdeg) {
        // squaring is coefficient-wise in characteristic 2
        let mut sq = vec![0u64; 2 * xp.len() - 1];
        for (i, &c) in xp.iter().enumerate() {
            sq[2 * i] = k.sqr(c);
        }
        xp = divrem(k, &trim(sq), &psi).1;
    }
    let g = gcd(k, &psi, &add(&xp, &vec![0, 1]));
    deg(&g) as usize
}

#[derive(Debug)]
pub struct FloorTest {
    pub ell: usize,
    pub kernel_x_field_degree: u32,
    pub fixed_x: usize,
    pub on_floor: bool,
}

pub fn floor_test(k: &F2m, a6: u64, ell: usize, kdeg: u32) -> FloorTest {
    let n = fixed_x_count(k, a6, ell, kdeg);
    let all = (ell * ell - 1) / 2;
    let eig = (ell - 1) / 2;
    assert!(n == all || n == eig, "unexpected fixed-x count {n} for l={ell}");
    FloorTest { ell, kernel_x_field_degree: kdeg, fixed_x: n, on_floor: n == eig }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::f2m::{standard_modulus, BPt, BinCurve};

    #[test]
    fn division_polynomials_vanish_on_torsion() {
        let m = 13;
        let k = F2m::new(standard_modulus(m));
        for b in [5u64, 77, 1234, 4321] {
            let c = BinCurve { k: k.clone(), a2: 1, a6: b };
            for n in [3usize, 5, 7, 9, 11] {
                let psi = division_poly(&k, b, n);
                assert_eq!(deg(&psi) as usize, (n * n - 1) / 2);
                for x in 1..(1u64 << m) {
                    if let Some(p) = c.lift_x(x) {
                        let tors = c.mul(p, n as u128) == BPt::Inf;
                        assert_eq!(eval(&k, &psi, x) == 0, tors, "b={b} n={n} x={x}");
                    }
                }
            }
        }
    }
}

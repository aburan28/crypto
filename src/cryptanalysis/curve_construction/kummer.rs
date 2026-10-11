//! Route 2 at toy size: Kummer covers `y^ℓ = f` of `E_0` over `F_2`.
//!
//! `f ∈ F_2(E_0)` has divisor `Σ k_i (p_i)` supported on `E_0(F_2) ≅ Z/4` with
//! every `k_i ≢ 0 mod ℓ`, so the cover is geometrically cyclic of degree `ℓ`,
//! totally ramified at the `r` points `p_i`, genus `1 + (ℓ − 1)r/2`, and its
//! automorphisms are defined over `F_2(μ_ℓ)`: Frobenius acts on the group with
//! order `d = ord_ℓ(2)`.  `f` is built from Miller functions on `⟨P⟩`.  Point
//! counts give `P_C`; the run checks that it is a Weil polynomial, that `E_0`
//! splits off, that the Prym's polynomial lies in `Z[T^d]` (the `μ_d`-stability
//! the structural bound relies on), and whether `P_{A_ℓ(E_0)}` divides it.

use super::gf2k::SmallField;
use super::lpoly::{char_poly_from_counts, divide_exact, is_weil_polynomial};
use super::weil::{curve_poly, mult_order, to_i128, trace_zero_poly};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Pt {
    O,
    A(u8, u8),
}

/// A line `cx·x + cy·y + c0` with coefficients in `F_2`.
type Line = (u8, u8, u8);

const A_COEFF: u8 = 0; // E_0

fn add(p: Pt, q: Pt) -> Pt {
    match (p, q) {
        (Pt::O, _) => q,
        (_, Pt::O) => p,
        (Pt::A(x1, y1), Pt::A(x2, y2)) => {
            if x1 == x2 {
                if y2 == x1 ^ y1 || x1 == 0 {
                    return Pt::O;
                }
                let l = x1 ^ y1;
                let x3 = A_COEFF;
                Pt::A(x3, x1 ^ ((l ^ 1) & x3))
            } else {
                let l = y1 ^ y2;
                let x3 = 1 ^ A_COEFF;
                Pt::A(x3, (l & (x1 ^ x3)) ^ x3 ^ y1)
            }
        }
    }
}

fn line(p: Pt, q: Pt) -> Option<Line> {
    match (p, q) {
        (Pt::O, Pt::O) => None,
        (Pt::O, Pt::A(x, _)) | (Pt::A(x, _), Pt::O) => Some((1, 0, x)),
        (Pt::A(x1, y1), Pt::A(x2, y2)) => {
            if x1 == x2 && (y2 == x1 ^ y1) {
                return Some((1, 0, x1));
            }
            let l = if x1 == x2 { x1 ^ y1 } else { y1 ^ y2 };
            Some((l, 1, y1 ^ (l & x1)))
        }
    }
}

fn vertical(p: Pt) -> Option<Line> {
    match p {
        Pt::O => None,
        Pt::A(x, _) => Some((1, 0, x)),
    }
}

type Function = BTreeMap<Line, i64>;

fn mul_into(f: &mut Function, g: &Function, e: i64) {
    for (l, k) in g {
        *f.entry(*l).or_insert(0) += k * e;
    }
    f.retain(|_, k| *k != 0);
}

/// Miller functions `f_j`, `div f_j = j(P) − (jP) − (j − 1)(O)`, `j = 0..=4`.
fn miller(p: Pt) -> Vec<Function> {
    let mut out = vec![Function::new(), Function::new()];
    let mut jp = p;
    for _ in 1..4 {
        let mut next = out.last().expect("non-empty").clone();
        if let Some(l) = line(jp, p) {
            *next.entry(l).or_insert(0) += 1;
        }
        let jp1 = add(jp, p);
        if let Some(v) = vertical(jp1) {
            *next.entry(v).or_insert(0) -= 1;
        }
        next.retain(|_, k| *k != 0);
        out.push(next);
        jp = jp1;
    }
    out
}

fn eval(field: &SmallField, f: &Function, x: u16, y: u16) -> u16 {
    let (mut num, mut den) = (1u16, 1u16);
    for (&(cx, cy, c0), &e) in f {
        let v = (if cx == 1 { x } else { 0 }) ^ (if cy == 1 { y } else { 0 }) ^ c0 as u16;
        assert!(v != 0, "line vanishes off E_0(F_2)");
        for _ in 0..e.abs() {
            if e > 0 {
                num = field.mul(num, v);
            } else {
                den = field.mul(den, v);
            }
        }
    }
    field.mul(num, field.inv(den))
}

#[derive(Clone, Debug, Serialize)]
pub struct Cover {
    pub support_multiples_of_p: Vec<u8>,
    pub weights: Vec<i64>,
    pub genus: usize,
    pub counts: Vec<i64>,
    pub char_poly: Vec<i128>,
    pub weil: bool,
    pub e0_splits_off: bool,
    pub prym_in_z_t_d: bool,
    pub contains_target: bool,
}

#[derive(Clone, Debug, Serialize)]
pub struct KummerFamily {
    pub ell: u64,
    pub d: u64,
    pub branch_points: usize,
    pub genus: usize,
    pub covers: usize,
    pub all_weil: bool,
    pub all_e0_split: bool,
    pub all_prym_mu_d_stable: bool,
    pub covers_containing_target: usize,
    pub structural_bound_allows: bool,
    pub examples: Vec<Cover>,
}

/// Every Kummer cover `y^ℓ = f` of `E_0` branched at `r` points of `E_0(F_2)`.
pub fn family(ell: u64, r: usize) -> KummerFamily {
    let rational: Vec<Pt> = [Pt::O, Pt::A(0, 1), Pt::A(1, 0), Pt::A(1, 1)]
        .into_iter()
        .collect();
    // generator of E_0(F_2) ≅ Z/4
    let gen = rational
        .iter()
        .copied()
        .find(|&p| add(add(p, p), add(p, p)) == Pt::O && add(p, p) != Pt::O)
        .expect("E_0(F_2) is cyclic of order 4");
    let mut mult = vec![Pt::O];
    for _ in 1..4 {
        let last = *mult.last().expect("non-empty");
        mult.push(add(last, gen));
    }
    let fj = miller(gen);
    let d = mult_order(2, ell);
    let target = to_i128(&trace_zero_poly(-1, ell as u32));
    let pe = to_i128(&curve_poly(-1));
    let genus = 1 + (ell as usize - 1) * r / 2;
    let fields: Vec<SmallField> = (1..=genus as u32).map(SmallField::new).collect();

    let w = 2 * ell as i64 + 2;
    let mut seen: BTreeSet<(Vec<u8>, Vec<i64>)> = BTreeSet::new();
    let mut covers = Vec::new();
    for mask in 0u8..16 {
        if mask.count_ones() as usize != r {
            continue;
        }
        let support: Vec<u8> = (0..4u8).filter(|j| mask >> j & 1 == 1).collect();
        // weight vectors in [−w, w]^r, ordered by Σ|k|, first representative per class
        let mut cands: Vec<Vec<i64>> = Vec::new();
        let mut cur = vec![-w; r];
        loop {
            let ok = cur.iter().sum::<i64>() == 0
                && cur.iter().all(|k| k.rem_euclid(ell as i64) != 0)
                && support
                    .iter()
                    .zip(&cur)
                    .map(|(&j, &k)| k * j as i64)
                    .sum::<i64>()
                    .rem_euclid(4)
                    == 0;
            if ok {
                cands.push(cur.clone());
            }
            let mut i = 0;
            while i < r && cur[i] == w {
                cur[i] = -w;
                i += 1;
            }
            if i == r {
                break;
            }
            cur[i] += 1;
        }
        cands.sort_by_key(|v| (v.iter().map(|k| k.abs()).sum::<i64>(), v.clone()));
        for k in cands {
            let k0 = k[0].rem_euclid(ell as i64);
            let inv = (1..ell as i64)
                .find(|s| s * k0 % ell as i64 == 1)
                .expect("unit");
            let class: Vec<i64> = k.iter().map(|x| (x * inv).rem_euclid(ell as i64)).collect();
            if !seen.insert((support.clone(), class)) {
                continue;
            }
            // f = f_4^m · ∏ f_{j_i}^{−k_i},  m = Σ k_i j_i / 4
            let m: i64 = support
                .iter()
                .zip(&k)
                .map(|(&j, &kk)| kk * j as i64)
                .sum::<i64>()
                / 4;
            let mut f = Function::new();
            mul_into(&mut f, &fj[4], m);
            for (&j, &kk) in support.iter().zip(&k) {
                mul_into(&mut f, &fj[j as usize], -kk);
            }
            let branch: BTreeSet<u8> = support.iter().copied().collect();
            let counts: Vec<i64> = fields
                .iter()
                .map(|field| {
                    let order = (field.q - 1) as u64;
                    let nk = |v: u16| -> i64 {
                        if !order.is_multiple_of(ell) {
                            1
                        } else if field.pow(v, order / ell) == 1 {
                            ell as i64
                        } else {
                            0
                        }
                    };
                    let rational_term = |j: u8| if branch.contains(&j) { 1 } else { nk(1) };
                    let mut c = rational_term(0); // O
                    for x in 0..field.q as u16 {
                        let ys: Vec<u16> = if x == 0 {
                            vec![1]
                        } else {
                            let ix = field.inv(x);
                            let cc = x ^ (A_COEFF as u16) ^ field.mul(ix, ix);
                            match field.solve_as(cc) {
                                Some(z) => vec![field.mul(x, z), field.mul(x, z ^ 1)],
                                None => vec![],
                            }
                        };
                        for y in ys {
                            if x <= 1 && y <= 1 {
                                let j = mult
                                    .iter()
                                    .position(|&p| p == Pt::A(x as u8, y as u8))
                                    .expect("rational point")
                                    as u8;
                                c += rational_term(j);
                            } else {
                                c += nk(eval(field, &f, x, y));
                            }
                        }
                    }
                    c
                })
                .collect();
            let Some(p) = char_poly_from_counts(&counts, 2) else {
                covers.push(Cover {
                    support_multiples_of_p: support.clone(),
                    weights: k.clone(),
                    genus,
                    counts,
                    char_poly: vec![],
                    weil: false,
                    e0_splits_off: false,
                    prym_in_z_t_d: false,
                    contains_target: false,
                });
                continue;
            };
            let weil = is_weil_polynomial(&p, 2);
            let prym = divide_exact(&p, &pe);
            let stable = prym.as_ref().is_some_and(|q| {
                q.iter()
                    .enumerate()
                    .all(|(i, c)| *c == 0 || (i as u64).is_multiple_of(d))
            });
            let contains = divide_exact(&p, &target).is_some();
            covers.push(Cover {
                support_multiples_of_p: support.clone(),
                weights: k.clone(),
                genus,
                counts,
                char_poly: p,
                weil,
                e0_splits_off: prym.is_some(),
                prym_in_z_t_d: stable,
                contains_target: contains,
            });
        }
    }
    KummerFamily {
        ell,
        d,
        branch_points: r,
        genus,
        covers: covers.len(),
        all_weil: covers.iter().all(|c| c.weil),
        all_e0_split: covers.iter().all(|c| c.e0_splits_off),
        all_prym_mu_d_stable: covers.iter().all(|c| c.prym_in_z_t_d),
        covers_containing_target: covers.iter().filter(|c| c.contains_target).count(),
        structural_bound_allows: r as u64 >= 2 * d,
        examples: covers.into_iter().take(4).collect(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn miller_functions_have_the_right_shape() {
        let p = Pt::A(1, 0);
        assert_eq!(add(add(p, p), add(p, p)), Pt::O);
        let f = miller(p);
        // f_4 = (y + x + 1)² / x
        let mut expect = Function::new();
        expect.insert((1, 1, 1), 2);
        expect.insert((1, 0, 0), -1);
        assert_eq!(f[4], expect);
    }

    #[test]
    fn two_branch_point_covers_are_weil_and_split() {
        let fam = family(3, 2);
        assert!(fam.covers > 0);
        assert!(fam.all_weil && fam.all_e0_split && fam.all_prym_mu_d_stable);
    }
}

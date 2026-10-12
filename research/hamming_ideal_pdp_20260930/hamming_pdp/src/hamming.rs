//! The Hamming ideals of La Scala–Marchesin–Tiwari for the constraint
//! `wt(x) ≤ w` on a block of `n` Boolean variables, in the paper's three
//! lifted presentations, plus a monomial-ideal control.

use crate::boolpoly::{Mono, Poly};
use std::collections::HashMap;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Encoding {
    Sub,
    Mono,
    C,
    Fc,
    Qfc,
}

impl Encoding {
    pub fn name(&self) -> &'static str {
        match self {
            Encoding::Sub => "SUB",
            Encoding::Mono => "MONO",
            Encoding::C => "C",
            Encoding::Fc => "FC",
            Encoding::Qfc => "QFC",
        }
    }
    pub fn parse(s: &str) -> Option<Encoding> {
        match s {
            "SUB" => Some(Encoding::Sub),
            "MONO" => Some(Encoding::Mono),
            "C" => Some(Encoding::C),
            "FC" => Some(Encoding::Fc),
            "QFC" => Some(Encoding::Qfc),
            _ => None,
        }
    }
}

pub struct Alloc {
    pub next: usize,
}
impl Alloc {
    pub fn fresh(&mut self) -> usize {
        let v = self.next;
        self.next += 1;
        v
    }
}

fn floor_log2(n: usize) -> u32 {
    63 - (n as u64).leading_zeros()
}

/// `⌊log₂ n⌋`: the digit index needed to determine any weight in `[0, n]`.
pub fn digit_bound(n: usize) -> u32 {
    floor_log2(n)
}

/// Root digit constraints for `wt ≤ w`: `y_{2^k} = 0` for `2^k > w`, and the
/// ANF of `[Σ_{2^k ≤ w} t_k 2^k > w]` over the remaining digits.
fn root_constraints(root_digits: &[Option<Poly>], w: u32, l: u32) -> Vec<Poly> {
    let mut out = Vec::new();
    let mut low: Vec<(u32, Poly)> = Vec::new();
    for k in 0..=l {
        let y = match &root_digits[k as usize] {
            Some(y) => y.clone(),
            None => continue,
        };
        if (1u32 << k) > w {
            out.push(y);
        } else {
            low.push((k, y));
        }
    }
    // ANF over the low digits.
    let m = low.len();
    let mut truth = vec![false; 1 << m];
    for a in 0..(1usize << m) {
        let mut val = 0u32;
        for (i, (k, _)) in low.iter().enumerate() {
            if (a >> i) & 1 == 1 {
                val += 1 << k;
            }
        }
        truth[a] = val > w;
    }
    // Möbius transform.
    let mut anf = truth.clone();
    for i in 0..m {
        for a in 0..(1usize << m) {
            if (a >> i) & 1 == 1 {
                anf[a] ^= anf[a ^ (1 << i)];
            }
        }
    }
    let mut p = Poly::zero();
    for a in 0..(1usize << m) {
        if anf[a] {
            let mut term = Poly::one();
            for (i, (_, y)) in low.iter().enumerate() {
                if (a >> i) & 1 == 1 {
                    term = term.mul(y);
                }
            }
            p.add_assign(&term);
        }
    }
    if !p.is_zero() {
        out.push(p);
    }
    out
}

/// Node of the balanced decomposition tree: ESFs as polynomials.
/// `esf[d]` for `d ≥ 1` (index 0 unused = 1); `None` where the ESF is not carried.
struct Node {
    size: usize,
    /// C-Hamming: e_d for 1 ≤ d ≤ min(size, 2^L). FC/QFC: only powers of two, at index 2^k.
    esf: Vec<Option<Poly>>,
    /// QFC product variables c_S, keyed by the subset mask over k.
    prods: HashMap<u32, Poly>,
}

pub struct Built {
    pub eqs: Vec<Poly>,
    pub aux_vars: usize,
}

pub fn hamming_ideal(enc: Encoding, coord_vars: &[usize], w: u32, alloc: &mut Alloc) -> Built {
    let n = coord_vars.len();
    let l = digit_bound(n);
    let start = alloc.next;
    let mut eqs = Vec::new();
    let root = match enc {
        Encoding::Mono => {
            // all products of w+1 distinct variables
            let mut idx: Vec<usize> = (0..=w as usize).collect();
            loop {
                let mut m = Mono::ONE;
                for &i in &idx {
                    m = m.mul(&Mono::var(coord_vars[i]));
                }
                eqs.push(Poly { terms: vec![m] });
                // next combination
                let k = w as usize + 1;
                let mut i = k;
                loop {
                    if i == 0 {
                        return Built { eqs, aux_vars: 0 };
                    }
                    i -= 1;
                    if idx[i] < n - k + i {
                        idx[i] += 1;
                        for j in i + 1..k {
                            idx[j] = idx[j - 1] + 1;
                        }
                        break;
                    }
                }
            }
        }
        Encoding::Sub => return Built { eqs, aux_vars: 0 },
        Encoding::C => build_c(coord_vars, 0, n - 1, l, alloc, &mut eqs),
        Encoding::Fc => build_fc(coord_vars, 0, n - 1, l, false, alloc, &mut eqs),
        Encoding::Qfc => build_fc(coord_vars, 0, n - 1, l, true, alloc, &mut eqs),
    };
    let mut digits: Vec<Option<Poly>> = Vec::new();
    for k in 0..=l {
        let d = 1usize << k;
        digits.push(if d < root.esf.len() { root.esf[d].clone() } else { None });
    }
    eqs.extend(root_constraints(&digits, w, l));
    Built { eqs, aux_vars: alloc.next - start }
}

fn leaf(var: usize, l: u32) -> Node {
    let mut esf = vec![None; 2];
    esf[1] = Some(Poly::var(var));
    let _ = l;
    Node { size: 1, esf, prods: HashMap::new() }
}

/// Convolution Hamming ideal (Theorem 3.5): `y_d = Σ_k y_k^L y_{d-k}^R`.
fn build_c(cv: &[usize], a: usize, b: usize, l: u32, alloc: &mut Alloc, eqs: &mut Vec<Poly>) -> Node {
    if a == b {
        return leaf(cv[a], l);
    }
    let m = (a + b) / 2;
    let left = build_c(cv, a, m, l, alloc, eqs);
    let right = build_c(cv, m + 1, b, l, alloc, eqs);
    let size = b - a + 1;
    let dmax = size.min(1usize << l);
    let mut esf = vec![None; dmax + 1];
    for d in 1..=dmax {
        let y = Poly::var(alloc.fresh());
        let mut rhs = Poly::zero();
        for k in 0..=d {
            let lk = get_esf(&left, k);
            let rk = get_esf(&right, d - k);
            match (lk, rk) {
                (Some(x), Some(y)) => rhs.add_assign(&x.mul(&y)),
                _ => {}
            }
        }
        eqs.push(y.add(&rhs));
        esf[d] = Some(y);
    }
    Node { size, esf, prods: HashMap::new() }
}

fn get_esf(node: &Node, d: usize) -> Option<Poly> {
    if d == 0 {
        return Some(Poly::one());
    }
    if d > node.size || d >= node.esf.len() {
        return None;
    }
    node.esf[d].clone()
}

/// Factorized convolution (Theorem 4.5) and its quadratic form (Section 5).
fn build_fc(cv: &[usize], a: usize, b: usize, l: u32, quad: bool, alloc: &mut Alloc, eqs: &mut Vec<Poly>) -> Node {
    if a == b {
        return leaf(cv[a], l);
    }
    let m = (a + b) / 2;
    let mut left = build_fc(cv, a, m, l, quad, alloc, eqs);
    let mut right = build_fc(cv, m + 1, b, l, quad, alloc, eqs);
    let size = b - a + 1;
    let rho = floor_log2(size).min(l) + 1;
    let mut esf = vec![None; (1usize << (rho - 1)) + 1];
    for k in 0..rho {
        let d = 1usize << k;
        let y = Poly::var(alloc.fresh());
        let mut rhs = Poly::zero();
        for j in 0..=d {
            let lt = lucas_product(&mut left, j as u32, quad, alloc, eqs);
            let rt = lucas_product(&mut right, (d - j) as u32, quad, alloc, eqs);
            if let (Some(x), Some(z)) = (lt, rt) {
                rhs.add_assign(&x.mul(&z));
            }
        }
        eqs.push(y.add(&rhs));
        esf[d] = Some(y);
    }
    Node { size, esf, prods: HashMap::new() }
}

/// `e_j` of a node as the Lucas product `Π_{h ∈ bits(j)} e_{2^h}`; `None`
/// (zero) if some required power-of-two ESF is not carried. In the
/// quadratic form the product is a single product variable `c_S`.
fn lucas_product(node: &mut Node, j: u32, quad: bool, alloc: &mut Alloc, eqs: &mut Vec<Poly>) -> Option<Poly> {
    if j == 0 {
        return Some(Poly::one());
    }
    if j as usize > node.size {
        return None;
    }
    // every bit must be carried
    let mut bits = Vec::new();
    for h in 0..32 {
        if (j >> h) & 1 == 1 {
            let d = 1usize << h;
            if d >= node.esf.len() || node.esf[d].is_none() {
                return None;
            }
            bits.push(h);
        }
    }
    if bits.len() == 1 {
        return node.esf[1usize << bits[0]].clone();
    }
    if !quad {
        let mut p = Poly::one();
        for h in &bits {
            p = p.mul(node.esf[1usize << h].as_ref().unwrap());
        }
        return Some(p);
    }
    // quadratic: c_S = c_{S \ {max}} · y_{2^max}, introduced on demand
    if let Some(c) = node.prods.get(&j) {
        return Some(c.clone());
    }
    let top = *bits.last().unwrap();
    let rest = j & !(1 << top);
    let c_rest = lucas_product(node, rest, quad, alloc, eqs).unwrap();
    let y = node.esf[1usize << top].clone().unwrap();
    let c = Poly::var(alloc.fresh());
    eqs.push(c.add(&c_rest.mul(&y)));
    node.prods.insert(j, c.clone());
    Some(c)
}

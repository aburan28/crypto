//! Symbolic elements of `F_{2^n}` with Boolean-polynomial coordinates, and
//! the Weil descent of Semaev's `S₃`.

use crate::boolpoly::Poly;
use crate::gf2n::Field;

#[derive(Clone)]
pub struct SymElem {
    /// Polynomial-basis coordinates, `n` of them.
    pub coords: Vec<Poly>,
}

impl SymElem {
    pub fn constant(f: &Field, c: u64) -> SymElem {
        SymElem {
            coords: (0..f.n).map(|t| Poly::constant((c >> t) & 1 == 1)).collect(),
        }
    }
    /// `x = Σ_j u_j · b_j` for the given basis elements and variables.
    pub fn from_vars(f: &Field, basis: &[u64], vars: &[usize]) -> SymElem {
        assert_eq!(basis.len(), vars.len());
        let mut coords = vec![Poly::zero(); f.n as usize];
        for (b, v) in basis.iter().zip(vars) {
            for t in 0..f.n as usize {
                if (b >> t) & 1 == 1 {
                    coords[t].add_assign(&Poly::var(*v));
                }
            }
        }
        SymElem { coords }
    }
    pub fn add(&self, o: &SymElem) -> SymElem {
        SymElem {
            coords: self.coords.iter().zip(&o.coords).map(|(a, b)| a.add(b)).collect(),
        }
    }
    pub fn mul(&self, o: &SymElem, f: &Field) -> SymElem {
        let n = f.n as usize;
        let mut coords = vec![Poly::zero(); n];
        // z^{i+k} mod irr
        let mut zpow = vec![1u64; 2 * n];
        for i in 1..2 * n {
            zpow[i] = f.mul(zpow[i - 1], 2);
        }
        for i in 0..n {
            if self.coords[i].is_zero() {
                continue;
            }
            for k in 0..n {
                if o.coords[k].is_zero() {
                    continue;
                }
                let p = self.coords[i].mul(&o.coords[k]);
                let z = zpow[i + k];
                for t in 0..n {
                    if (z >> t) & 1 == 1 {
                        coords[t].add_assign(&p);
                    }
                }
            }
        }
        SymElem { coords }
    }
    pub fn sqr(&self, f: &Field) -> SymElem {
        let n = f.n as usize;
        let mut coords = vec![Poly::zero(); n];
        for i in 0..n {
            if self.coords[i].is_zero() {
                continue;
            }
            let z = f.pow(2, 2 * i as u64);
            for t in 0..n {
                if (z >> t) & 1 == 1 {
                    coords[t].add_assign(&self.coords[i]);
                }
            }
        }
        SymElem { coords }
    }
}

/// `S₃(x₁, x₂, x₃) = (x₁+x₂)² x₃² + x₁x₂x₃ + (x₁x₂)² + b` with `x₃` known,
/// `b = 1`: the `n` coordinate equations.
pub fn s3_descent(f: &Field, x1: &SymElem, x2: &SymElem, x3: u64) -> Vec<Poly> {
    let s = x1.add(x2).sqr(f);
    let t1 = s.mul(&SymElem::constant(f, f.sqr(x3)), f);
    let p = x1.mul(x2, f);
    let t2 = p.mul(&SymElem::constant(f, x3), f);
    let t3 = p.sqr(f);
    let one = SymElem::constant(f, 1);
    let total = t1.add(&t2).add(&t3).add(&one);
    total.coords.into_iter().filter(|c| !c.is_zero()).collect()
}

/// Numeric `S₃` for verification.
pub fn s3_value(f: &Field, x1: u64, x2: u64, x3: u64) -> u64 {
    let s = f.sqr(x1 ^ x2);
    let p = f.mul(x1, x2);
    f.mul(s, f.sqr(x3)) ^ f.mul(p, x3) ^ f.sqr(p) ^ 1
}

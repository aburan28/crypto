//! Boolean polynomials over `F_2` in the quotient by the field equations:
//! square-free monomials as bit masks, degrevlex with variable `0` largest.

use std::cmp::Ordering;

pub const W: usize = 8;
pub const MAX_VARS: usize = 64 * W;

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct Mono(pub [u64; W]);

impl Mono {
    pub const ONE: Mono = Mono([0; W]);
    #[inline]
    pub fn var(i: usize) -> Mono {
        let mut m = [0u64; W];
        m[i / 64] |= 1 << (i % 64);
        Mono(m)
    }
    #[inline]
    pub fn degree(&self) -> u32 {
        self.0.iter().map(|w| w.count_ones()).sum()
    }
    #[inline]
    pub fn mul(&self, o: &Mono) -> Mono {
        let mut m = [0u64; W];
        for k in 0..W {
            m[k] = self.0[k] | o.0[k];
        }
        Mono(m)
    }
    #[inline]
    pub fn divides(&self, o: &Mono) -> bool {
        (0..W).all(|k| self.0[k] & !o.0[k] == 0)
    }
    #[inline]
    pub fn div(&self, o: &Mono) -> Mono {
        let mut m = [0u64; W];
        for k in 0..W {
            m[k] = self.0[k] & !o.0[k];
        }
        Mono(m)
    }
    #[inline]
    pub fn has(&self, i: usize) -> bool {
        (self.0[i / 64] >> (i % 64)) & 1 == 1
    }
    #[inline]
    pub fn without(&self, i: usize) -> Mono {
        let mut m = self.0;
        m[i / 64] &= !(1 << (i % 64));
        Mono(m)
    }
    pub fn vars(&self) -> Vec<usize> {
        let mut v = Vec::new();
        for k in 0..W {
            let mut w = self.0[k];
            while w != 0 {
                let b = w.trailing_zeros() as usize;
                v.push(k * 64 + b);
                w &= w - 1;
            }
        }
        v
    }
}

/// Degrevlex, variable `0` largest: higher degree wins; at equal degree the
/// monomial *without* the smallest (highest-index) differing variable wins.
#[inline]
pub fn cmp_mono(a: &Mono, b: &Mono) -> Ordering {
    let da = a.degree();
    let db = b.degree();
    if da != db {
        return da.cmp(&db);
    }
    for k in (0..W).rev() {
        let d = a.0[k] ^ b.0[k];
        if d != 0 {
            let h = 63 - d.leading_zeros();
            return if (a.0[k] >> h) & 1 == 1 {
                Ordering::Less
            } else {
                Ordering::Greater
            };
        }
    }
    Ordering::Equal
}

/// Terms sorted descending; the first term is the leading monomial.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Poly {
    pub terms: Vec<Mono>,
}

impl Poly {
    pub fn zero() -> Poly {
        Poly { terms: vec![] }
    }
    pub fn one() -> Poly {
        Poly { terms: vec![Mono::ONE] }
    }
    pub fn var(i: usize) -> Poly {
        Poly { terms: vec![Mono::var(i)] }
    }
    pub fn constant(b: bool) -> Poly {
        if b {
            Poly::one()
        } else {
            Poly::zero()
        }
    }
    /// Sort and cancel duplicate monomials.
    pub fn from_terms(mut t: Vec<Mono>) -> Poly {
        t.sort_by(|a, b| cmp_mono(b, a));
        let mut out: Vec<Mono> = Vec::with_capacity(t.len());
        for m in t {
            if let Some(last) = out.last() {
                if *last == m {
                    out.pop();
                    continue;
                }
            }
            out.push(m);
        }
        Poly { terms: out }
    }
    #[inline]
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }
    pub fn is_one(&self) -> bool {
        self.terms.len() == 1 && self.terms[0] == Mono::ONE
    }
    pub fn lm(&self) -> Option<&Mono> {
        self.terms.first()
    }
    pub fn degree(&self) -> u32 {
        self.terms.first().map(|m| m.degree()).unwrap_or(0)
    }
    pub fn add(&self, o: &Poly) -> Poly {
        let mut out = Vec::with_capacity(self.terms.len() + o.terms.len());
        let (mut i, mut j) = (0, 0);
        while i < self.terms.len() && j < o.terms.len() {
            match cmp_mono(&self.terms[i], &o.terms[j]) {
                Ordering::Greater => {
                    out.push(self.terms[i]);
                    i += 1;
                }
                Ordering::Less => {
                    out.push(o.terms[j]);
                    j += 1;
                }
                Ordering::Equal => {
                    i += 1;
                    j += 1;
                }
            }
        }
        out.extend_from_slice(&self.terms[i..]);
        out.extend_from_slice(&o.terms[j..]);
        Poly { terms: out }
    }
    pub fn add_assign(&mut self, o: &Poly) {
        *self = self.add(o);
    }
    pub fn mul_mono(&self, m: &Mono) -> Poly {
        Poly::from_terms(self.terms.iter().map(|t| t.mul(m)).collect())
    }
    pub fn mul(&self, o: &Poly) -> Poly {
        if self.is_zero() || o.is_zero() {
            return Poly::zero();
        }
        let mut t = Vec::with_capacity(self.terms.len() * o.terms.len());
        for a in &self.terms {
            for b in &o.terms {
                t.push(a.mul(b));
            }
        }
        Poly::from_terms(t)
    }
    /// Substitute `x_i := b`.
    pub fn substitute(&self, i: usize, b: bool) -> Poly {
        let mut t = Vec::with_capacity(self.terms.len());
        for m in &self.terms {
            if m.has(i) {
                if b {
                    t.push(m.without(i));
                }
            } else {
                t.push(*m);
            }
        }
        Poly::from_terms(t)
    }
    /// Evaluate at a point given as a bitset over variables.
    pub fn eval(&self, point: &[u64; W]) -> bool {
        let mut acc = false;
        for m in &self.terms {
            if (0..W).all(|k| m.0[k] & !point[k] == 0) {
                acc ^= true;
            }
        }
        acc
    }
    pub fn max_var_index(&self) -> Option<usize> {
        self.terms.iter().flat_map(|m| m.vars()).max()
    }
}

impl std::fmt::Display for Poly {
    fn fmt(&self, fm: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.terms.is_empty() {
            return write!(fm, "0");
        }
        let parts: Vec<String> = self
            .terms
            .iter()
            .map(|m| {
                let v = m.vars();
                if v.is_empty() {
                    "1".to_string()
                } else {
                    v.iter().map(|i| format!("x{i}")).collect::<Vec<_>>().join("*")
                }
            })
            .collect();
        write!(fm, "{}", parts.join(" + "))
    }
}

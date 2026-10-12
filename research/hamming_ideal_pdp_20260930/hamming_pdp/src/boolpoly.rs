//! Boolean polynomials over `F_2` in the quotient by the field equations:
//! square-free monomials as bit masks, degrevlex with variable `0` largest.

use std::cell::Cell;
use std::cmp::Ordering;

/// 1024 Boolean variables. Frozen cells use fewer than 512. C-Hamming at
/// `n = 53` needs two 53-bit coordinate blocks plus the convolution
/// auxiliaries, which does not fit in 512, so the width is 16 limbs.
/// High limbs of a smaller system stay zero.
pub const W: usize = 16;
pub const MAX_VARS: usize = 64 * W;

thread_local! {
    /// How many `u64` limbs of a monomial can be nonzero on this thread.
    /// High limbs of every monomial built inside a `LimbsGuard` are zero, so
    /// degree, order and divisibility may ignore them. The guard is restored
    /// on drop: a later system with a wider variable set must not inherit it.
    static ACTIVE_LIMBS: Cell<usize> = const { Cell::new(W) };
}

#[inline]
pub fn active_limbs() -> usize {
    ACTIVE_LIMBS.with(|c| c.get())
}

/// Sets the limb width for the current thread and restores the previous
/// width when dropped.
pub struct LimbsGuard {
    prev: usize,
}

impl LimbsGuard {
    pub fn set(n: usize) -> Self {
        let n = n.clamp(1, W);
        let prev = ACTIVE_LIMBS.with(|c| c.replace(n));
        Self { prev }
    }
}

impl Drop for LimbsGuard {
    fn drop(&mut self) {
        ACTIVE_LIMBS.with(|c| c.set(self.prev));
    }
}

/// Smallest limb count that covers every nonzero word of `polys`.
pub fn limbs_covering(polys: &[Poly]) -> usize {
    let mut n = 1usize;
    for p in polys {
        for m in &p.terms {
            for k in (0..W).rev() {
                if m.0[k] != 0 {
                    n = n.max(k + 1);
                    break;
                }
            }
        }
    }
    n
}

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
        let n = active_limbs();
        let mut s = 0u32;
        for k in 0..n {
            s += self.0[k].count_ones();
        }
        s
    }
    #[inline]
    pub fn mul(&self, o: &Mono) -> Mono {
        let n = active_limbs();
        let mut m = [0u64; W];
        for k in 0..n {
            m[k] = self.0[k] | o.0[k];
        }
        Mono(m)
    }
    #[inline]
    pub fn divides(&self, o: &Mono) -> bool {
        let n = active_limbs();
        for k in 0..n {
            if self.0[k] & !o.0[k] != 0 {
                return false;
            }
        }
        true
    }
    #[inline]
    pub fn div(&self, o: &Mono) -> Mono {
        let n = active_limbs();
        let mut m = [0u64; W];
        for k in 0..n {
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
        let n = active_limbs();
        for k in 0..n {
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

/// Ascending key whose order is descending degrevlex (`cmp_mono(b, a)`).
/// Degree is reversed; at equal degree the high limb compares first, and a
/// set high bit makes the monomial smaller, so it sorts later.
#[inline]
fn mono_desc_key(m: &Mono) -> (u32, [u64; W]) {
    let n = active_limbs();
    let mut deg = 0u32;
    let mut limbs = [0u64; W];
    for k in 0..n {
        deg += m.0[k].count_ones();
        limbs[W - 1 - k] = m.0[k];
    }
    (!deg, limbs)
}

/// Degrevlex, variable `0` largest: higher degree wins; at equal degree the
/// monomial *without* the smallest (highest-index) differing variable wins.
#[inline]
pub fn cmp_mono(a: &Mono, b: &Mono) -> Ordering {
    let n = active_limbs();
    let mut da = 0u32;
    let mut db = 0u32;
    for k in 0..n {
        da += a.0[k].count_ones();
        db += b.0[k].count_ones();
    }
    if da != db {
        return da.cmp(&db);
    }
    for k in (0..n).rev() {
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
        Poly {
            terms: vec![Mono::ONE],
        }
    }
    pub fn var(i: usize) -> Poly {
        Poly {
            terms: vec![Mono::var(i)],
        }
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
        match t.len() {
            0 | 1 => return Poly { terms: t },
            n if n <= 24 => t.sort_unstable_by(|a, b| cmp_mono(b, a)),
            _ => t.sort_by_cached_key(mono_desc_key),
        }
        // In-place GF(2) stack cancel. Writes stay at indices `< r`, so
        // `t[r]` is still the sorted value.
        let mut w = 0usize;
        for r in 0..t.len() {
            if w > 0 && t[w - 1] == t[r] {
                w -= 1;
            } else {
                if w != r {
                    t[w] = t[r];
                }
                w += 1;
            }
        }
        t.truncate(w);
        Poly { terms: t }
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
        let n = active_limbs();
        for m in &self.terms {
            let mut hits = true;
            for k in 0..n {
                if m.0[k] & !point[k] != 0 {
                    hits = false;
                    break;
                }
            }
            if hits {
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
                    v.iter()
                        .map(|i| format!("x{i}"))
                        .collect::<Vec<_>>()
                        .join("*")
                }
            })
            .collect();
        write!(fm, "{}", parts.join(" + "))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn reference_from_terms(mut t: Vec<Mono>) -> Poly {
        t.sort_by(|a, b| cmp_mono(b, a));
        let mut out: Vec<Mono> = Vec::new();
        for m in t {
            if out.last() == Some(&m) {
                out.pop();
            } else {
                out.push(m);
            }
        }
        Poly { terms: out }
    }

    #[test]
    fn from_terms_matches_degrevlex_cancel() {
        let _g = LimbsGuard::set(2);
        let mut monos = vec![Mono::ONE];
        for i in 0..80 {
            monos.push(Mono::var(i));
            if i % 3 == 0 {
                monos.push(Mono::var(i));
                monos.push(Mono::var(i));
            }
            for j in (i + 1..80).step_by(7) {
                monos.push(Mono::var(i).mul(&Mono::var(j)));
            }
        }
        let got = Poly::from_terms(monos.clone());
        let want = reference_from_terms(monos);
        assert_eq!(got, want);
        let triple = Poly::from_terms(vec![Mono::var(3), Mono::var(3), Mono::var(3)]);
        assert_eq!(triple.terms, vec![Mono::var(3)]);
    }
}

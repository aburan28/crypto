//! Sparse multivariate polynomials over Z, resultants, and exact integer root finding.
//!
//! [`MPoly`] stores a polynomial in a fixed number of variables `x_0, .., x_{n-1}` as a map
//! from exponent vectors to nonzero [`Int`] coefficients, ordered lexicographically with
//! `x_0` most significant.  Everything here is exact.
//!
//! * Resultants ([`MPoly::resultant`]) are Sylvester determinants computed by fraction-free
//!   Bareiss elimination over the polynomial ring Z[other variables]; every Bareiss division
//!   is an exact polynomial division ([`MPoly::div_exact`]).
//! * Univariate helpers live in [`upoly`] (coefficient vectors, constant term first),
//!   including [`upoly::integer_roots`]: all integer roots of an integer polynomial,
//!   optionally restricted to `|r| <= bound`.  Method: make the polynomial squarefree
//!   (rigorously: a prime p not dividing the leading coefficient with `f mod p` squarefree
//!   proves `f` squarefree over Z; otherwise divide by an exact Z[x] gcd with `f'`), find all
//!   roots of `f mod p` (distinct-degree part `gcd(f, x^p - x)` and Cantor–Zassenhaus
//!   splitting), Hensel-lift each (simple) root to `p^e > 2 * bound` and test the symmetric
//!   representative exactly.  Since `f mod p` is squarefree every integer root in range is the
//!   lift of exactly one root mod p, so the output is complete and every returned value is
//!   checked by exact evaluation.  The only randomness (choice of splitting elements) affects
//!   running time, not the result.
use crate::int::Int;
use std::collections::BTreeMap;
use std::fmt;
use std::ops::{Add, Mul, Neg, Sub};

fn int(v: i64) -> Int {
    Int::from(v)
}

/// Sparse polynomial in `nvars` variables with integer coefficients.
#[derive(Clone, PartialEq, Eq)]
pub struct MPoly {
    nvars: usize,
    terms: BTreeMap<Vec<u32>, Int>,
}

impl fmt::Debug for MPoly {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.terms.is_empty() {
            return write!(f, "0");
        }
        let mut first = true;
        for (e, c) in self.terms.iter().rev() {
            if !first {
                write!(f, " + ")?;
            }
            first = false;
            write!(f, "{}", c)?;
            for (i, &k) in e.iter().enumerate() {
                if k == 1 {
                    write!(f, "*x{}", i)?;
                } else if k > 1 {
                    write!(f, "*x{}^{}", i, k)?;
                }
            }
        }
        Ok(())
    }
}

impl MPoly {
    pub fn zero(nvars: usize) -> MPoly {
        MPoly {
            nvars,
            terms: BTreeMap::new(),
        }
    }
    pub fn constant(nvars: usize, c: &Int) -> MPoly {
        MPoly::monomial(nvars, &vec![0; nvars], c)
    }
    /// The variable `x_i`.
    pub fn var(nvars: usize, i: usize) -> MPoly {
        let mut e = vec![0; nvars];
        e[i] = 1;
        MPoly::monomial(nvars, &e, &Int::one())
    }
    /// `c * x^e`.
    pub fn monomial(nvars: usize, e: &[u32], c: &Int) -> MPoly {
        assert_eq!(e.len(), nvars);
        let mut p = MPoly::zero(nvars);
        if !c.is_zero() {
            p.terms.insert(e.to_vec(), c.clone());
        }
        p
    }
    /// Sum of `c * x^e` over the given terms (repeated exponents are added).
    pub fn from_terms(nvars: usize, terms: &[(Vec<u32>, Int)]) -> MPoly {
        let mut p = MPoly::zero(nvars);
        for (e, c) in terms {
            p.add_term(e, c);
        }
        p
    }
    pub fn nvars(&self) -> usize {
        self.nvars
    }
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }
    pub fn num_terms(&self) -> usize {
        self.terms.len()
    }
    /// Terms in increasing lexicographic order of exponents.
    pub fn terms(&self) -> impl Iterator<Item = (&Vec<u32>, &Int)> {
        self.terms.iter()
    }
    /// Coefficient of `x^e` (zero if absent).
    pub fn coeff(&self, e: &[u32]) -> Int {
        self.terms.get(e).cloned().unwrap_or_else(Int::zero)
    }
    /// Constant term.
    pub fn constant_term(&self) -> Int {
        self.coeff(&vec![0; self.nvars])
    }
    /// Lex-leading term (exponent, coefficient).
    pub fn leading_term(&self) -> Option<(&Vec<u32>, &Int)> {
        self.terms.iter().next_back()
    }
    pub fn add_term(&mut self, e: &[u32], c: &Int) {
        if c.is_zero() {
            return;
        }
        match self.terms.get_mut(e) {
            Some(v) => {
                *v = &*v + c;
                if v.is_zero() {
                    self.terms.remove(e);
                }
            }
            None => {
                self.terms.insert(e.to_vec(), c.clone());
            }
        }
    }
    /// Degree in `x_i` (`-1` for the zero polynomial is represented as `None`).
    pub fn deg_in(&self, i: usize) -> Option<u32> {
        self.terms.keys().map(|e| e[i]).max()
    }
    pub fn total_degree(&self) -> Option<u32> {
        self.terms.keys().map(|e| e.iter().sum()).max()
    }
    /// Whether only variable `i` (or none) occurs.
    pub fn only_var(&self, i: usize) -> bool {
        self.terms
            .keys()
            .all(|e| e.iter().enumerate().all(|(j, &k)| j == i || k == 0))
    }
    pub fn scale(&self, c: &Int) -> MPoly {
        if c.is_zero() {
            return MPoly::zero(self.nvars);
        }
        MPoly {
            nvars: self.nvars,
            terms: self.terms.iter().map(|(e, v)| (e.clone(), v * c)).collect(),
        }
    }
    /// Exact division by an integer, `None` unless every coefficient is divisible.
    pub fn div_exact_int(&self, c: &Int) -> Option<MPoly> {
        let mut out = MPoly::zero(self.nvars);
        for (e, v) in &self.terms {
            let (q, r) = v.divrem_trunc(c);
            if !r.is_zero() {
                return None;
            }
            out.terms.insert(e.clone(), q);
        }
        Some(out)
    }
    /// Every coefficient reduced into `[0, m)`.
    pub fn modulo(&self, m: &Int) -> MPoly {
        let mut out = MPoly::zero(self.nvars);
        for (e, v) in &self.terms {
            out.add_term(e, &v.modulo(m));
        }
        out
    }
    /// gcd of the coefficients (non-negative; 0 for the zero polynomial).
    pub fn content(&self) -> Int {
        let mut g = Int::zero();
        for v in self.terms.values() {
            g = g.gcd(v);
            if g == Int::one() {
                break;
            }
        }
        g
    }
    /// `self / content`, with positive leading coefficient.
    pub fn primitive_part(&self) -> MPoly {
        if self.is_zero() {
            return self.clone();
        }
        let mut c = self.content();
        if self.leading_term().unwrap().1.is_neg() {
            c = -c;
        }
        self.div_exact_int(&c).unwrap()
    }
    /// Max-norm of the coefficients.
    pub fn norm_inf(&self) -> Int {
        self.terms
            .values()
            .map(|v| v.abs())
            .max()
            .unwrap_or_else(Int::zero)
    }
    /// Squared Euclidean norm of the coefficient vector.
    pub fn norm2_sq(&self) -> Int {
        let mut s = Int::zero();
        for v in self.terms.values() {
            s = &s + &(v * v);
        }
        s
    }
    /// Exact evaluation at a point.
    pub fn eval(&self, x: &[Int]) -> Int {
        assert_eq!(x.len(), self.nvars);
        // cache powers
        let mut pw: Vec<Vec<Int>> = vec![vec![Int::one()]; self.nvars];
        let mut s = Int::zero();
        for (e, c) in &self.terms {
            let mut t = c.clone();
            for (i, &k) in e.iter().enumerate() {
                if k == 0 {
                    continue;
                }
                while pw[i].len() <= k as usize {
                    let nx = pw[i].last().unwrap() * &x[i];
                    pw[i].push(nx);
                }
                t = &t * &pw[i][k as usize];
            }
            s = &s + &t;
        }
        s
    }
    /// Substitute `x_i = v` (the result keeps `nvars` variables, `x_i` no longer occurs).
    pub fn subs_value(&self, i: usize, v: &Int) -> MPoly {
        let mut pw = vec![Int::one()];
        let mut out = MPoly::zero(self.nvars);
        for (e, c) in &self.terms {
            let k = e[i] as usize;
            while pw.len() <= k {
                let nx = pw.last().unwrap() * v;
                pw.push(nx);
            }
            let mut e2 = e.clone();
            e2[i] = 0;
            out.add_term(&e2, &(c * &pw[k]));
        }
        out
    }
    /// `x_i -> c_i x_i` for all i (e.g. `p(xX, yY)`).
    pub fn scale_vars(&self, c: &[Int]) -> MPoly {
        assert_eq!(c.len(), self.nvars);
        let mut pw: Vec<Vec<Int>> = vec![vec![Int::one()]; self.nvars];
        let mut out = MPoly::zero(self.nvars);
        for (e, v) in &self.terms {
            let mut t = v.clone();
            for (i, &k) in e.iter().enumerate() {
                while pw[i].len() <= k as usize {
                    let nx = pw[i].last().unwrap() * &c[i];
                    pw[i].push(nx);
                }
                t = &t * &pw[i][k as usize];
            }
            out.add_term(e, &t);
        }
        out
    }
    /// Multiply by the monomial `x^e`.
    pub fn mul_monomial(&self, e: &[u32]) -> MPoly {
        MPoly {
            nvars: self.nvars,
            terms: self
                .terms
                .iter()
                .map(|(k, v)| (k.iter().zip(e).map(|(a, b)| a + b).collect(), v.clone()))
                .collect(),
        }
    }
    /// `x_i -> x_i + c` (Taylor shift in one variable).
    pub fn shift_var(&self, i: usize, c: &Int) -> MPoly {
        let co = self.coeffs_in(i);
        // Horner in x_i + c
        let xi = &MPoly::var(self.nvars, i) + &MPoly::constant(self.nvars, c);
        let mut acc = MPoly::zero(self.nvars);
        for q in co.iter().rev() {
            acc = &(&acc * &xi) + q;
        }
        acc
    }
    /// View as a polynomial in `x_i`: `self = sum_k co[k] x_i^k` (the `co[k]` do not involve
    /// `x_i`).  Empty for the zero polynomial.
    pub fn coeffs_in(&self, i: usize) -> Vec<MPoly> {
        let d = match self.deg_in(i) {
            Some(d) => d as usize,
            None => return vec![],
        };
        let mut co = vec![MPoly::zero(self.nvars); d + 1];
        for (e, v) in &self.terms {
            let mut e2 = e.clone();
            let k = e2[i] as usize;
            e2[i] = 0;
            co[k].terms.insert(e2, v.clone());
        }
        co
    }
    /// Inverse of [`MPoly::coeffs_in`].
    pub fn from_coeffs_in(nvars: usize, i: usize, co: &[MPoly]) -> MPoly {
        let mut out = MPoly::zero(nvars);
        for (k, q) in co.iter().enumerate() {
            for (e, v) in &q.terms {
                let mut e2 = e.clone();
                e2[i] += k as u32;
                out.add_term(&e2, v);
            }
        }
        out
    }
    /// Coefficient vector in `x_i` if no other variable occurs.
    pub fn to_univariate(&self, i: usize) -> Option<Vec<Int>> {
        if !self.only_var(i) {
            return None;
        }
        let d = match self.deg_in(i) {
            Some(d) => d as usize,
            None => return Some(vec![]),
        };
        let mut v = vec![Int::zero(); d + 1];
        for (e, c) in &self.terms {
            v[e[i] as usize] = c.clone();
        }
        Some(v)
    }
    pub fn from_univariate(nvars: usize, i: usize, c: &[Int]) -> MPoly {
        let mut out = MPoly::zero(nvars);
        for (k, v) in c.iter().enumerate() {
            let mut e = vec![0; nvars];
            e[i] = k as u32;
            out.add_term(&e, v);
        }
        out
    }
    /// Exact division: `Some(q)` with `self = q * d`, or `None` if `d` does not divide `self`
    /// in Z[x_0..x_{n-1}].  (Lex leading-term division; correct because Z[x] is a domain and
    /// lex is a monomial order, so `LT(q d) = LT(q) LT(d)`.)
    pub fn div_exact(&self, d: &MPoly) -> Option<MPoly> {
        assert!(!d.is_zero(), "division by zero polynomial");
        let (de, dc) = d
            .leading_term()
            .map(|(e, c)| (e.clone(), c.clone()))
            .unwrap();
        if d.terms.len() == 1 {
            let mut out = MPoly::zero(self.nvars);
            for (e, v) in &self.terms {
                if e.iter().zip(&de).any(|(a, b)| a < b) {
                    return None;
                }
                let (q, r) = v.divrem_trunc(&dc);
                if !r.is_zero() {
                    return None;
                }
                out.terms
                    .insert(e.iter().zip(&de).map(|(a, b)| a - b).collect(), q);
            }
            return Some(out);
        }
        let mut r = self.clone();
        let mut q = MPoly::zero(self.nvars);
        while let Some((re, rc)) = r.leading_term().map(|(e, c)| (e.clone(), c.clone())) {
            if re.iter().zip(&de).any(|(a, b)| a < b) {
                return None;
            }
            let (qc, rem) = rc.divrem_trunc(&dc);
            if !rem.is_zero() {
                return None;
            }
            let qe: Vec<u32> = re.iter().zip(&de).map(|(a, b)| a - b).collect();
            for (e, v) in &d.terms {
                let e2: Vec<u32> = e.iter().zip(&qe).map(|(a, b)| a + b).collect();
                r.add_term(&e2, &-(v * &qc));
            }
            debug_assert!(!r.terms.contains_key(&re));
            q.terms.insert(qe, qc);
        }
        Some(q)
    }
    /// Resultant with respect to `x_i` (a polynomial in the other variables), as the
    /// determinant of the Sylvester matrix, computed by fraction-free Bareiss elimination over
    /// Z[other variables].  Convention: `Res(f, g) = 0` if either is zero; if `deg f = 0` then
    /// `Res = f^{deg g}` (and symmetrically).
    pub fn resultant(&self, other: &MPoly, i: usize) -> MPoly {
        let n = self.nvars;
        if self.is_zero() || other.is_zero() {
            return MPoly::zero(n);
        }
        let a = self.coeffs_in(i);
        let b = other.coeffs_in(i);
        let (m, k) = (a.len() - 1, b.len() - 1);
        if m == 0 {
            return a[0].pow(k as u32);
        }
        if k == 0 {
            return b[0].pow(m as u32);
        }
        let s = m + k;
        let mut mat = vec![vec![MPoly::zero(n); s]; s];
        // rows 0..k: shifts of a (highest coefficient first); rows k..s: shifts of b
        for r in 0..k {
            for (t, c) in a.iter().rev().enumerate() {
                mat[r][r + t] = c.clone();
            }
        }
        for r in 0..m {
            for (t, c) in b.iter().rev().enumerate() {
                mat[k + r][r + t] = c.clone();
            }
        }
        det_bareiss(mat, n)
    }
    /// Resultant w.r.t. `x_i` of two polynomials in which only `x_i` and `x_keep` occur, as a
    /// dense univariate polynomial in `x_keep` (same Sylvester convention as
    /// [`MPoly::resultant`]); Bareiss over Z[x_keep] with dense coefficient vectors, which is
    /// much faster than the sparse general version.  `None` if another variable occurs.
    pub fn resultant_univariate(&self, other: &MPoly, i: usize, keep: usize) -> Option<Vec<Int>> {
        let ok = |p: &MPoly| {
            p.terms.keys().all(|e| {
                e.iter()
                    .enumerate()
                    .all(|(j, &k)| j == i || j == keep || k == 0)
            })
        };
        if !ok(self) || !ok(other) {
            return None;
        }
        if self.is_zero() || other.is_zero() {
            return Some(vec![]);
        }
        let to_u = |p: &MPoly| -> Vec<Vec<Int>> {
            p.coeffs_in(i)
                .iter()
                .map(|c| c.to_univariate(keep).unwrap())
                .collect()
        };
        let (a, b) = (to_u(self), to_u(other));
        let (m, k) = (a.len() - 1, b.len() - 1);
        let upow = |c: &Vec<Int>, e: usize| {
            let mut r = vec![Int::one()];
            for _ in 0..e {
                r = upoly::mul(&r, c);
            }
            r
        };
        if m == 0 {
            return Some(upow(&a[0], k));
        }
        if k == 0 {
            return Some(upow(&b[0], m));
        }
        let s = m + k;
        let mut mat = vec![vec![Vec::<Int>::new(); s]; s];
        for r in 0..k {
            for (t, c) in a.iter().rev().enumerate() {
                mat[r][r + t] = upoly::trim(c.clone());
            }
        }
        for r in 0..m {
            for (t, c) in b.iter().rev().enumerate() {
                mat[k + r][r + t] = upoly::trim(c.clone());
            }
        }
        Some(upoly::det_bareiss(mat))
    }
    pub fn pow(&self, e: u32) -> MPoly {
        let mut r = MPoly::constant(self.nvars, &Int::one());
        let mut b = self.clone();
        let mut e = e;
        while e > 0 {
            if e & 1 == 1 {
                r = &r * &b;
            }
            e >>= 1;
            if e > 0 {
                b = &b * &b;
            }
        }
        r
    }
    /// Formal partial derivative in `x_i`.
    pub fn derivative(&self, i: usize) -> MPoly {
        let mut out = MPoly::zero(self.nvars);
        for (e, v) in &self.terms {
            if e[i] > 0 {
                let mut e2 = e.clone();
                e2[i] -= 1;
                out.add_term(&e2, &(v * &Int::from(e[i] as u64)));
            }
        }
        out
    }
}

/// Determinant of a square matrix over Z[x_0..] by fraction-free Bareiss elimination (exact).
pub fn det_bareiss(mut m: Vec<Vec<MPoly>>, nvars: usize) -> MPoly {
    let s = m.len();
    if s == 0 {
        return MPoly::constant(nvars, &Int::one());
    }
    let mut sign = false;
    let mut prev = MPoly::constant(nvars, &Int::one());
    for k in 0..s - 1 {
        if m[k][k].is_zero() {
            // pick the smallest nonzero pivot below (fewest terms) to limit growth
            let piv = (k + 1..s)
                .filter(|&r| !m[r][k].is_zero())
                .min_by_key(|&r| m[r][k].num_terms());
            match piv {
                Some(r) => {
                    m.swap(k, r);
                    sign = !sign;
                }
                None => return MPoly::zero(nvars),
            }
        }
        for i in k + 1..s {
            for j in k + 1..s {
                let t = &(&m[i][j] * &m[k][k]) - &(&m[i][k] * &m[k][j]);
                m[i][j] = if k == 0 {
                    t
                } else {
                    t.div_exact(&prev).expect("Bareiss division must be exact")
                };
            }
            m[i][k] = MPoly::zero(nvars);
        }
        prev = m[k][k].clone();
    }
    let d = m[s - 1][s - 1].clone();
    if sign {
        -d
    } else {
        d
    }
}

impl Add<&MPoly> for &MPoly {
    type Output = MPoly;
    fn add(self, o: &MPoly) -> MPoly {
        assert_eq!(self.nvars, o.nvars);
        let (big, small) = if self.terms.len() >= o.terms.len() {
            (self, o)
        } else {
            (o, self)
        };
        let mut out = big.clone();
        for (e, v) in &small.terms {
            out.add_term(e, v);
        }
        out
    }
}
impl Sub<&MPoly> for &MPoly {
    type Output = MPoly;
    fn sub(self, o: &MPoly) -> MPoly {
        assert_eq!(self.nvars, o.nvars);
        let mut out = self.clone();
        for (e, v) in &o.terms {
            out.add_term(e, &-v);
        }
        out
    }
}
impl Neg for &MPoly {
    type Output = MPoly;
    fn neg(self) -> MPoly {
        MPoly {
            nvars: self.nvars,
            terms: self.terms.iter().map(|(e, v)| (e.clone(), -v)).collect(),
        }
    }
}
impl Neg for MPoly {
    type Output = MPoly;
    fn neg(self) -> MPoly {
        -&self
    }
}
impl Mul<&MPoly> for &MPoly {
    type Output = MPoly;
    fn mul(self, o: &MPoly) -> MPoly {
        assert_eq!(self.nvars, o.nvars);
        let mut out = MPoly::zero(self.nvars);
        let mut e = vec![0u32; self.nvars];
        for (e1, v1) in &self.terms {
            for (e2, v2) in &o.terms {
                for t in 0..self.nvars {
                    e[t] = e1[t] + e2[t];
                }
                out.add_term(&e, &(v1 * v2));
            }
        }
        out
    }
}
macro_rules! owned_ops {
    ($tr:ident, $m:ident) => {
        impl $tr<MPoly> for MPoly {
            type Output = MPoly;
            fn $m(self, o: MPoly) -> MPoly {
                (&self).$m(&o)
            }
        }
        impl $tr<&MPoly> for MPoly {
            type Output = MPoly;
            fn $m(self, o: &MPoly) -> MPoly {
                (&self).$m(o)
            }
        }
    };
}
owned_ops!(Add, add);
owned_ops!(Sub, sub);
owned_ops!(Mul, mul);

/// Univariate integer polynomials as coefficient vectors (constant term first, no trailing
/// zeros after [`upoly::trim`]).
pub mod upoly {
    use super::int;
    use crate::int::Int;

    pub fn trim(mut f: Vec<Int>) -> Vec<Int> {
        while f.last().is_some_and(|c| c.is_zero()) {
            f.pop();
        }
        f
    }
    /// Degree (`None` for zero).
    pub fn degree(f: &[Int]) -> Option<usize> {
        f.iter().rposition(|c| !c.is_zero())
    }
    pub fn eval(f: &[Int], x: &Int) -> Int {
        let mut s = Int::zero();
        for c in f.iter().rev() {
            s = &(&s * x) + c;
        }
        s
    }
    /// `f(x) mod m` in `[0, m)`.
    pub fn eval_mod(f: &[Int], x: &Int, m: &Int) -> Int {
        let mut s = Int::zero();
        for c in f.iter().rev() {
            s = (&(&s * x) + c).modulo(m);
        }
        s
    }
    pub fn derivative(f: &[Int]) -> Vec<Int> {
        (1..f.len()).map(|i| &f[i] * &Int::from(i as u64)).collect()
    }
    pub fn content(f: &[Int]) -> Int {
        f.iter().fold(Int::zero(), |g, c| g.gcd(c))
    }
    pub fn primitive_part(f: &[Int]) -> Vec<Int> {
        let f = trim(f.to_vec());
        if f.is_empty() {
            return f;
        }
        let mut c = content(&f);
        if f.last().unwrap().is_neg() {
            c = -c;
        }
        f.iter().map(|v| v.divrem_trunc(&c).0).collect()
    }
    pub fn mul(a: &[Int], b: &[Int]) -> Vec<Int> {
        if a.is_empty() || b.is_empty() {
            return vec![];
        }
        let mut out = vec![Int::zero(); a.len() + b.len() - 1];
        for (i, x) in a.iter().enumerate() {
            if x.is_zero() {
                continue;
            }
            for (j, y) in b.iter().enumerate() {
                out[i + j] = &out[i + j] + &(x * y);
            }
        }
        trim(out)
    }
    /// Exact division `a / b` in Z[x], `None` if not exact.
    pub fn div_exact(a: &[Int], b: &[Int]) -> Option<Vec<Int>> {
        let b = trim(b.to_vec());
        assert!(!b.is_empty(), "division by zero polynomial");
        let mut r = trim(a.to_vec());
        if r.is_empty() {
            return Some(vec![]);
        }
        if r.len() < b.len() {
            return None;
        }
        let db = b.len() - 1;
        let lb = b[db].clone();
        let mut q = vec![Int::zero(); r.len() - db];
        for i in (0..q.len()).rev() {
            let (c, rem) = r[i + db].divrem_trunc(&lb);
            if !rem.is_zero() {
                return None;
            }
            if !c.is_zero() {
                for j in 0..=db {
                    r[i + j] = &r[i + j] - &(&c * &b[j]);
                }
            }
            q[i] = c;
        }
        if r.iter().any(|c| !c.is_zero()) {
            return None;
        }
        Some(trim(q))
    }
    pub fn sub(a: &[Int], b: &[Int]) -> Vec<Int> {
        let n = a.len().max(b.len());
        let z = Int::zero();
        trim(
            (0..n)
                .map(|i| a.get(i).unwrap_or(&z) - b.get(i).unwrap_or(&z))
                .collect(),
        )
    }
    /// Determinant of a square matrix over Z[x] (dense entries) by fraction-free Bareiss.
    pub fn det_bareiss(mut m: Vec<Vec<Vec<Int>>>) -> Vec<Int> {
        let s = m.len();
        if s == 0 {
            return vec![Int::one()];
        }
        let mut sign = false;
        let mut prev = vec![Int::one()];
        for k in 0..s - 1 {
            if m[k][k].is_empty() {
                let piv = (k + 1..s)
                    .filter(|&r| !m[r][k].is_empty())
                    .min_by_key(|&r| m[r][k].len());
                match piv {
                    Some(r) => {
                        m.swap(k, r);
                        sign = !sign;
                    }
                    None => return vec![],
                }
            }
            for i in k + 1..s {
                for j in k + 1..s {
                    let t = sub(&mul(&m[i][j], &m[k][k]), &mul(&m[i][k], &m[k][j]));
                    m[i][j] = if k == 0 {
                        t
                    } else {
                        div_exact(&t, &prev).expect("Bareiss division must be exact")
                    };
                }
                m[i][k] = vec![];
            }
            prev = m[k][k].clone();
        }
        let d = m[s - 1][s - 1].clone();
        if sign {
            d.iter().map(|c| -c).collect()
        } else {
            d
        }
    }
    /// Pseudo-remainder `prem(a, b)` (`lc(b)^{da-db+1} a = q b + r`).
    pub fn prem(a: &[Int], b: &[Int]) -> Vec<Int> {
        let mut r = trim(a.to_vec());
        let b = trim(b.to_vec());
        let db = b.len() - 1;
        let lb = b[db].clone();
        while !r.is_empty() && r.len() > db {
            let dr = r.len() - 1;
            let c = r[dr].clone();
            for v in r.iter_mut() {
                *v = &*v * &lb;
            }
            for j in 0..=db {
                r[dr - db + j] = &r[dr - db + j] - &(&c * &b[j]);
            }
            r = trim(r);
        }
        r
    }
    /// gcd in Z[x] (primitive, positive leading coefficient) by the primitive PRS.
    pub fn gcd(a: &[Int], b: &[Int]) -> Vec<Int> {
        let mut a = primitive_part(a);
        let mut b = primitive_part(b);
        if a.is_empty() {
            return b;
        }
        if b.is_empty() {
            return a;
        }
        let cg = content(&[content(a.as_slice()), content(b.as_slice())]);
        if a.len() < b.len() {
            std::mem::swap(&mut a, &mut b);
        }
        while !b.is_empty() {
            let r = prem(&a, &b);
            a = b;
            b = primitive_part(&r);
        }
        let mut g = primitive_part(&a);
        if !cg.is_zero() && cg != Int::one() {
            g = g.iter().map(|v| v * &cg).collect();
        }
        g
    }

    // ---- arithmetic in F_p[x] (p < 2^62), coefficient vectors as u64 ----
    fn mm(a: u64, b: u64, p: u64) -> u64 {
        ((a as u128 * b as u128) % p as u128) as u64
    }
    fn pw(mut a: u64, mut e: u64, p: u64) -> u64 {
        let mut r = 1;
        while e > 0 {
            if e & 1 == 1 {
                r = mm(r, a, p);
            }
            a = mm(a, a, p);
            e >>= 1;
        }
        r
    }
    fn inv(a: u64, p: u64) -> u64 {
        pw(a, p - 2, p)
    }
    fn tr(mut f: Vec<u64>) -> Vec<u64> {
        while f.last() == Some(&0) {
            f.pop();
        }
        f
    }
    fn pmulmod(a: &[u64], b: &[u64], m: &[u64], p: u64) -> Vec<u64> {
        if a.is_empty() || b.is_empty() {
            return vec![];
        }
        let mut out = vec![0u128; a.len() + b.len() - 1];
        let pp = p as u128;
        for (i, &x) in a.iter().enumerate() {
            if x == 0 {
                continue;
            }
            for (j, &y) in b.iter().enumerate() {
                out[i + j] = (out[i + j] + x as u128 * y as u128) % pp;
            }
        }
        prem_p(&out.into_iter().map(|v| v as u64).collect::<Vec<_>>(), m, p)
    }
    fn prem_p(a: &[u64], m: &[u64], p: u64) -> Vec<u64> {
        let mut r = tr(a.to_vec());
        let dm = m.len() - 1;
        let li = inv(m[dm], p);
        while r.len() > dm {
            let dr = r.len() - 1;
            let c = mm(r[dr], li, p);
            if c != 0 {
                for j in 0..=dm {
                    let t = mm(c, m[j], p);
                    r[dr - dm + j] = (r[dr - dm + j] + p - t) % p;
                }
            }
            r.pop();
            r = tr(r);
        }
        r
    }
    fn pgcd(a: &[u64], b: &[u64], p: u64) -> Vec<u64> {
        let mut a = tr(a.to_vec());
        let mut b = tr(b.to_vec());
        while !b.is_empty() {
            let r = prem_p(&a, &b, p);
            a = b;
            b = r;
        }
        if let Some(&l) = a.last() {
            let li = inv(l, p);
            for v in a.iter_mut() {
                *v = mm(*v, li, p);
            }
        }
        a
    }
    fn ppowmod(base: &[u64], mut e: u64, m: &[u64], p: u64) -> Vec<u64> {
        let mut r = vec![1u64];
        let mut b = prem_p(base, m, p);
        while e > 0 {
            if e & 1 == 1 {
                r = pmulmod(&r, &b, m, p);
            }
            e >>= 1;
            if e > 0 {
                b = pmulmod(&b, &b, m, p);
            }
        }
        r
    }
    fn pdiv(a: &[u64], b: &[u64], p: u64) -> Vec<u64> {
        let mut r = tr(a.to_vec());
        let db = b.len() - 1;
        let li = inv(b[db], p);
        let mut q = vec![0u64; r.len().saturating_sub(db)];
        while r.len() > db {
            let dr = r.len() - 1;
            let c = mm(r[dr], li, p);
            q[dr - db] = c;
            for j in 0..=db {
                let t = mm(c, b[j], p);
                r[dr - db + j] = (r[dr - db + j] + p - t) % p;
            }
            r.pop();
        }
        q
    }
    /// All roots in F_p of a squarefree monic split polynomial g (product of distinct linears).
    fn split_roots(g: &[u64], p: u64, seed: &mut u64, out: &mut Vec<u64>) {
        let d = g.len() - 1;
        if d == 0 {
            return;
        }
        if d == 1 {
            out.push((p - g[0]) % p);
            return;
        }
        if p == 2 {
            for x in 0..2u64 {
                let v = g.iter().rev().fold(0u64, |s, &c| (s * x + c) % 2);
                if v == 0 {
                    out.push(x);
                }
            }
            return;
        }
        loop {
            *seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let a = (*seed >> 11) % p;
            let h = ppowmod(&[a, 1], (p - 1) / 2, g, p);
            let mut h1 = h.clone();
            if h1.is_empty() {
                h1.push(0);
            }
            h1[0] = (h1[0] + p - 1) % p;
            let f1 = pgcd(g, &h1, p);
            let d1 = f1.len().saturating_sub(1);
            if d1 > 0 && d1 < d {
                let f2 = pdiv(g, &f1, p);
                split_roots(&f1, p, seed, out);
                split_roots(&f2, p, seed, out);
                return;
            }
        }
    }
    /// Distinct roots of `f mod p` in F_p.
    /// (`p` an odd prime below 2^63, or 2; `f` given by residues in `[0, p)`.)
    pub fn roots_mod_p(f: &[u64], p: u64) -> Vec<u64> {
        let f = tr(f.to_vec());
        if f.len() <= 1 {
            return vec![];
        }
        let xp = ppowmod(&[0, 1], p, &f, p);
        let mut xpx = xp.clone();
        while xpx.len() < 2 {
            xpx.push(0);
        }
        xpx[1] = (xpx[1] + p - 1) % p;
        let g = pgcd(&f, &tr(xpx), p);
        let mut out = vec![];
        let mut seed = p ^ 0x5DEECE66D;
        split_roots(&g, p, &mut seed, &mut out);
        out.sort();
        out
    }
    /// Deterministic Miller–Rabin for `n < 2^64` (first 12 prime bases).
    pub fn is_prime_u64(n: u64) -> bool {
        if n < 2 {
            return false;
        }
        for q in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
            if n.is_multiple_of(q) {
                return n == q;
            }
        }
        let (mut d, mut s) = (n - 1, 0);
        while d % 2 == 0 {
            d /= 2;
            s += 1;
        }
        'w: for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
            let mut x = pw(a, d, n);
            if x == 1 || x == n - 1 {
                continue;
            }
            for _ in 1..s {
                x = mm(x, x, n);
                if x == n - 1 {
                    continue 'w;
                }
            }
            return false;
        }
        true
    }
    fn reduce_p(f: &[Int], p: u64) -> Vec<u64> {
        f.iter().map(|c| c.mod_u64(p)).collect()
    }

    /// Whether `f` is squarefree modulo the prime `p` (and `p` does not divide `lc(f)`).
    fn good_prime(f: &[Int], p: u64) -> bool {
        let fp = reduce_p(f, p);
        if fp.last() == Some(&0) {
            return false;
        }
        let dp: Vec<u64> = (1..fp.len()).map(|i| mm(fp[i], i as u64 % p, p)).collect();
        pgcd(&fp, &dp, p).len() == 1
    }

    /// All integer roots of `f` (sorted, without multiplicity), restricted to `|r| <= bound`
    /// when a bound is given (otherwise the Cauchy bound is used).  Panics on the zero
    /// polynomial (every integer is a root).  Exact and complete; see the module docs.
    pub fn integer_roots(f: &[Int], bound: Option<&Int>) -> Vec<Int> {
        let f = trim(f.to_vec());
        assert!(!f.is_empty(), "integer_roots of the zero polynomial");
        let mut out = vec![];
        // strip x^v
        let v = f.iter().position(|c| !c.is_zero()).unwrap();
        if v > 0 {
            out.push(Int::zero());
        }
        let f = primitive_part(&f[v..]);
        if f.len() <= 1 {
            return out;
        }
        // Cauchy bound: 1 + max |a_i / a_n|
        let lc = f.last().unwrap().abs();
        let mut cb = Int::zero();
        for c in &f[..f.len() - 1] {
            let q = c.abs().div_floor(&lc);
            if q > cb {
                cb = q;
            }
        }
        cb = &cb + &int(2);
        let b = match bound {
            Some(b) if b.abs() < cb => b.abs(),
            _ => cb,
        };
        if f.len() == 2 {
            // linear: -f0 / f1
            let (q, r) = (-&f[0]).divrem_trunc(&f[1]);
            if r.is_zero() && q.abs() <= b {
                out.push(q);
            }
            out.sort();
            return out;
        }
        // find a prime with f mod p squarefree
        let mut sq = f.clone();
        let mut p = (1u64 << 31) - 1;
        let mut tries = 0;
        let mut checked_sqfree = false;
        loop {
            while !is_prime_u64(p) {
                p -= 2;
            }
            if good_prime(&sq, p) {
                break;
            }
            tries += 1;
            p -= 2;
            if tries == 3 && !checked_sqfree {
                // possibly not squarefree over Z: divide by the exact gcd with f'
                let g = gcd(&sq, &derivative(&sq));
                if g.len() > 1 {
                    sq = primitive_part(&div_exact(&sq, &g).expect("gcd divides"));
                    if sq.len() <= 2 {
                        if sq.len() == 2 {
                            let (q, r) = (-&sq[0]).divrem_trunc(&sq[1]);
                            if r.is_zero() && q.abs() <= b {
                                out.push(q);
                            }
                        }
                        out.sort();
                        out.dedup();
                        return out;
                    }
                }
                checked_sqfree = true;
            }
        }
        let rp = roots_mod_p(&reduce_p(&sq, p), p);
        // Hensel lift to p^e > 2b
        let pi = Int::from(p);
        let two_b = &b * &int(2);
        let df = derivative(&sq);
        for r0 in rp {
            let mut m = pi.clone();
            let mut r = Int::from(r0);
            while m <= two_b {
                let m2 = &m * &m;
                let fr = eval_mod(&sq, &r, &m2);
                let dfr = eval_mod(&df, &r, &m2);
                let di = dfr.inv_mod(&m2).expect("simple root mod p");
                r = (&r - &(&fr * &di)).modulo(&m2);
                m = m2;
            }
            // symmetric representative
            let half = m.div_floor(&int(2));
            let c = if r > half { &r - &m } else { r };
            if c.abs() <= b && eval(&sq, &c).is_zero() {
                out.push(c);
            }
        }
        out.sort();
        out.dedup();
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::Rng;

    fn ip(v: &[i64]) -> Vec<Int> {
        v.iter().map(|&x| int(x)).collect()
    }
    fn x(n: usize, i: usize) -> MPoly {
        MPoly::var(n, i)
    }
    fn c(n: usize, v: i64) -> MPoly {
        MPoly::constant(n, &int(v))
    }

    #[test]
    fn ring_ops_and_eval() {
        let (a, b) = (x(2, 0), x(2, 1));
        let p = &(&(&a * &a) + &(&b * &c(2, 3))) - &c(2, 5); // x^2 + 3y - 5
        let q = &a - &b;
        let pq = &p * &q;
        assert_eq!(
            pq.eval(&[int(2), int(-1)]),
            &p.eval(&[int(2), int(-1)]) * &int(3)
        );
        assert_eq!(pq.div_exact(&q).unwrap(), p);
        assert_eq!(pq.div_exact(&p).unwrap(), q);
        assert!(p.div_exact(&q).is_none());
        assert_eq!(
            p.subs_value(1, &int(2)).to_univariate(0).unwrap(),
            ip(&[1, 0, 1])
        );
        let s = p.scale_vars(&[int(10), int(100)]);
        assert_eq!(s.coeff(&[2, 0]), int(100));
        assert_eq!(s.coeff(&[0, 1]), int(300));
        let sh = p.shift_var(0, &int(1)); // (x+1)^2 + 3y - 5
        assert_eq!(sh.eval(&[int(3), int(7)]), p.eval(&[int(4), int(7)]));
        let co = pq.coeffs_in(1);
        assert_eq!(MPoly::from_coeffs_in(2, 1, &co), pq);
        assert_eq!(pq.scale(&int(6)).content(), int(6));
        assert_eq!(
            pq.scale(&int(-6))
                .primitive_part()
                .leading_term()
                .unwrap()
                .1
                .signum(),
            1
        );
        assert_eq!(p.derivative(0), a.scale(&int(2)));
        assert_eq!(p.pow(3), &(&p * &p) * &p);
    }

    #[test]
    fn resultants() {
        // Res_y(y^2 - x, y - 2) = 4 - x  (Sylvester of y^2+0y-x and y-2)
        let (xx, yy) = (x(2, 0), x(2, 1));
        let f = &(&yy * &yy) - &xx;
        let g = &yy - &c(2, 2);
        let r = f.resultant(&g, 1);
        assert_eq!(r.to_univariate(0).unwrap(), ip(&[4, -1]));
        // Res_x(x^2 - 2, x^2 - 3) = (√2-√3)(√2+√3)(-√2-√3)(-√2+√3)... = 1
        let f1 = MPoly::from_univariate(1, 0, &ip(&[-2, 0, 1]));
        let f2 = MPoly::from_univariate(1, 0, &ip(&[-3, 0, 1]));
        assert_eq!(f1.resultant(&f2, 0).constant_term(), int(1));
        // common root -> 0
        let g1 = MPoly::from_univariate(1, 0, &ip(&[-6, 1, 1])); // (x+3)(x-2)
        let g2 = MPoly::from_univariate(1, 0, &ip(&[-2, 1])); // x-2
        assert!(g1.resultant(&g2, 0).is_zero());
        // trivariate: eliminating z then y from planted common root
        let mut rng = Rng::new(5);
        let n = 3;
        let root = [int(3), int(-2), int(5)];
        let mk = |rng: &mut Rng| {
            let mut p = MPoly::zero(n);
            for e in [
                [2, 0, 0],
                [0, 2, 0],
                [0, 0, 2],
                [1, 1, 0],
                [0, 1, 1],
                [1, 0, 1],
                [1, 0, 0],
                [0, 0, 1],
            ] {
                p.add_term(&e, &int((rng.next() % 41) as i64 - 20));
            }
            let v = p.eval(&root);
            &p - &MPoly::constant(n, &v)
        };
        let (p1, p2, p3) = (mk(&mut rng), mk(&mut rng), mk(&mut rng));
        let r1 = p1.resultant(&p2, 2);
        let r2 = p1.resultant(&p3, 2);
        let r = r1.resultant(&r2, 1);
        let u = r.to_univariate(0).unwrap();
        assert!(!u.is_empty());
        assert!(upoly::eval(&u, &root[0]).is_zero());
        assert!(upoly::integer_roots(&u, None).contains(&root[0]));
    }

    #[test]
    fn integer_roots_basic() {
        // (x-3)(x+5)(2x-1)(x^2+1)
        let mut f = ip(&[1]);
        for g in [ip(&[-3, 1]), ip(&[5, 1]), ip(&[-1, 2]), ip(&[1, 0, 1])] {
            f = upoly::mul(&f, &g);
        }
        assert_eq!(upoly::integer_roots(&f, None), ip(&[-5, 3]));
        assert_eq!(upoly::integer_roots(&f, Some(&int(4))), ip(&[3]));
        // multiple roots and x factor: x^2 (x-7)^3 (x+2)
        let mut g = ip(&[0, 0, 1]);
        for _ in 0..3 {
            g = upoly::mul(&g, &ip(&[-7, 1]));
        }
        g = upoly::mul(&g, &ip(&[2, 1]));
        assert_eq!(upoly::integer_roots(&g, None), ip(&[-2, 0, 7]));
        // no integer roots
        assert!(upoly::integer_roots(&ip(&[-2, 0, 1]), None).is_empty());
        // (x^2+1)^2 (x-1): non-squarefree non-linear factor
        let h = upoly::mul(&upoly::mul(&ip(&[1, 0, 1]), &ip(&[1, 0, 1])), &ip(&[-1, 1]));
        assert_eq!(upoly::integer_roots(&h, None), ip(&[1]));
    }

    #[test]
    fn integer_roots_large() {
        // degree 40, roots of ~300 bits, random cofactor with ~2000-bit coefficients
        let mut rng = Rng::new(17);
        let big = |rng: &mut Rng, limbs: usize| {
            let mut v = Int::zero();
            for _ in 0..limbs {
                v = &v.shl(64) + &Int::from(rng.next());
            }
            if rng.next() & 1 == 1 {
                -v
            } else {
                v
            }
        };
        let roots: Vec<Int> = (0..6).map(|_| big(&mut rng, 5)).collect();
        let mut f = ip(&[1]);
        for r in &roots {
            f = upoly::mul(&f, &[-r, int(1)]);
        }
        let cof: Vec<Int> = (0..35).map(|_| big(&mut rng, 30)).collect();
        f = upoly::mul(&f, &cof);
        let t = std::time::Instant::now();
        let got = upoly::integer_roots(&f, None);
        let mut want = roots.clone();
        want.sort();
        // the random cofactor could in principle have integer roots too: check containment
        for r in &want {
            assert!(got.contains(r));
        }
        for r in &got {
            assert!(upoly::eval(&f, r).is_zero());
        }
        eprintln!("integer_roots deg {}: {:?}", f.len() - 1, t.elapsed());
    }
}

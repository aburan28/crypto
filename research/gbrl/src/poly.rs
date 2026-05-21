use crate::field::Fp;
use crate::monomial::Monomial;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;

#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub struct Term {
    pub coef: Fp,
    pub mono: Monomial,
}

// Sparse multivariate polynomial. Terms stored in DESCENDING monomial order
// (so terms[0] is the leading term) with no zero coefficients and no duplicates.
#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub struct Poly {
    pub terms: Vec<Term>,
    pub nvars: usize,
}

impl Poly {
    pub fn zero(nvars: usize) -> Self { Poly { terms: vec![], nvars } }
    pub fn is_zero(&self) -> bool { self.terms.is_empty() }
    pub fn lt(&self) -> Option<&Term> { self.terms.first() }
    pub fn lm(&self) -> Option<&Monomial> { self.terms.first().map(|t| &t.mono) }
    pub fn lc(&self) -> Option<Fp> { self.terms.first().map(|t| t.coef) }
    pub fn nterms(&self) -> usize { self.terms.len() }

    pub fn total_degree(&self) -> u32 {
        self.terms.iter().map(|t| t.mono.degree()).max().unwrap_or(0)
    }

    // Build from an unsorted term list. Sorts, drops zero coefs, combines duplicates.
    pub fn from_terms(mut terms: Vec<Term>, nvars: usize) -> Self {
        terms.retain(|t| !t.coef.is_zero());
        terms.sort_by(|a, b| b.mono.cmp(&a.mono));
        let mut out: Vec<Term> = Vec::with_capacity(terms.len());
        for t in terms {
            if let Some(last) = out.last_mut() {
                if last.mono == t.mono {
                    last.coef = last.coef + t.coef;
                    if last.coef.is_zero() { out.pop(); }
                    continue;
                }
            }
            out.push(t);
        }
        Poly { terms: out, nvars }
    }

    // Merge-add two polynomials whose term lists are already sorted descending.
    pub fn add(&self, other: &Poly) -> Poly {
        debug_assert_eq!(self.nvars, other.nvars);
        let mut i = 0;
        let mut j = 0;
        let mut out: Vec<Term> = Vec::with_capacity(self.terms.len() + other.terms.len());
        while i < self.terms.len() && j < other.terms.len() {
            match self.terms[i].mono.cmp(&other.terms[j].mono) {
                Ordering::Greater => { out.push(self.terms[i].clone()); i += 1; }
                Ordering::Less    => { out.push(other.terms[j].clone()); j += 1; }
                Ordering::Equal   => {
                    let c = self.terms[i].coef + other.terms[j].coef;
                    if !c.is_zero() {
                        out.push(Term { coef: c, mono: self.terms[i].mono.clone() });
                    }
                    i += 1; j += 1;
                }
            }
        }
        while i < self.terms.len() { out.push(self.terms[i].clone()); i += 1; }
        while j < other.terms.len() { out.push(other.terms[j].clone()); j += 1; }
        Poly { terms: out, nvars: self.nvars }
    }

    pub fn neg(&self) -> Poly {
        Poly {
            terms: self.terms.iter().map(|t| Term { coef: -t.coef, mono: t.mono.clone() }).collect(),
            nvars: self.nvars,
        }
    }

    pub fn sub(&self, other: &Poly) -> Poly { self.add(&other.neg()) }

    // Multiply every term by (c * m). Preserves sort order.
    pub fn mul_term(&self, c: Fp, m: &Monomial) -> Poly {
        if c.is_zero() { return Poly::zero(self.nvars); }
        Poly {
            terms: self.terms.iter().map(|t| Term {
                coef: t.coef * c,
                mono: t.mono.mul(m),
            }).collect(),
            nvars: self.nvars,
        }
    }

    // Multiply by a scalar.
    pub fn scale(&self, c: Fp) -> Poly {
        self.mul_term(c, &Monomial::one(self.nvars))
    }

    // Full polynomial multiplication. Naïve O(|self| * |other|); fine for
    // small/medium Semaev systems and the toy benchmark instances. Use FFT-
    // backed schoolbook only if it becomes a measurable bottleneck.
    pub fn mul(&self, other: &Poly) -> Poly {
        if self.is_zero() || other.is_zero() { return Poly::zero(self.nvars); }
        let mut terms = Vec::with_capacity(self.terms.len() * other.terms.len());
        for a in &self.terms {
            for b in &other.terms {
                terms.push(Term {
                    coef: a.coef * b.coef,
                    mono: a.mono.mul(&b.mono),
                });
            }
        }
        Poly::from_terms(terms, self.nvars)
    }

    // Evaluate this polynomial at a concrete point. `point[i]` is the value
    // of variable x_i.
    pub fn eval(&self, point: &[Fp]) -> Fp {
        assert_eq!(point.len(), self.nvars);
        let mut sum = Fp::zero();
        for t in &self.terms {
            let mut val = t.coef;
            for (xi, &e) in point.iter().zip(&t.mono.exp) {
                for _ in 0..e {
                    val = val * *xi;
                }
            }
            sum = sum + val;
        }
        sum
    }

    // Substitute a constant for one variable AND drop that variable from the
    // representation. Output has nvars = self.nvars - 1. Used to specialize
    // Semaev S_{m+1}(x_R, x_1, ..., x_m) to a fixed target x_R.
    pub fn eliminate(&self, var_idx: usize, value: Fp) -> Poly {
        assert!(var_idx < self.nvars);
        let new_nvars = self.nvars - 1;
        let mut terms = Vec::with_capacity(self.terms.len());
        for t in &self.terms {
            let e = t.mono.exp[var_idx];
            let mut val_pow = Fp::one();
            for _ in 0..e {
                val_pow = val_pow * value;
            }
            let new_exp: Vec<u32> = t.mono.exp.iter().enumerate()
                .filter(|(i, _)| *i != var_idx)
                .map(|(_, e)| *e)
                .collect();
            terms.push(Term {
                coef: t.coef * val_pow,
                mono: Monomial { exp: new_exp },
            });
        }
        Poly::from_terms(terms, new_nvars)
    }

    // Embed a polynomial with `self.nvars` variables into a `target_nvars`-var
    // ambient by mapping var i → var_map[i]. Used to lift a univariate factor
    // base polynomial into the multi-variable decomposition system.
    pub fn embed(&self, target_nvars: usize, var_map: &[usize]) -> Poly {
        assert_eq!(var_map.len(), self.nvars);
        assert!(var_map.iter().all(|&i| i < target_nvars));
        let mut terms = Vec::with_capacity(self.terms.len());
        for t in &self.terms {
            let mut new_exp = vec![0u32; target_nvars];
            for (src, &dst) in var_map.iter().enumerate() {
                new_exp[dst] += t.mono.exp[src];
            }
            terms.push(Term { coef: t.coef, mono: Monomial { exp: new_exp } });
        }
        Poly::from_terms(terms, target_nvars)
    }
}

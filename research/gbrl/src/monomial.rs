use serde::{Deserialize, Serialize};
use std::cmp::Ordering;

#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
pub struct Monomial {
    pub exp: Vec<u32>,
}

impl Monomial {
    pub fn one(n: usize) -> Self { Monomial { exp: vec![0; n] } }

    pub fn var(i: usize, n: usize) -> Self {
        let mut e = vec![0; n];
        e[i] = 1;
        Monomial { exp: e }
    }

    pub fn degree(&self) -> u32 { self.exp.iter().sum() }
    pub fn nvars(&self) -> usize { self.exp.len() }

    pub fn mul(&self, other: &Monomial) -> Monomial {
        debug_assert_eq!(self.exp.len(), other.exp.len());
        Monomial { exp: self.exp.iter().zip(&other.exp).map(|(a, b)| a + b).collect() }
    }

    pub fn divides(&self, other: &Monomial) -> bool {
        debug_assert_eq!(self.exp.len(), other.exp.len());
        self.exp.iter().zip(&other.exp).all(|(a, b)| a <= b)
    }

    pub fn div(&self, other: &Monomial) -> Option<Monomial> {
        debug_assert_eq!(self.exp.len(), other.exp.len());
        if !other.divides(self) { return None; }
        Some(Monomial { exp: self.exp.iter().zip(&other.exp).map(|(a, b)| a - b).collect() })
    }

    pub fn lcm(&self, other: &Monomial) -> Monomial {
        debug_assert_eq!(self.exp.len(), other.exp.len());
        Monomial { exp: self.exp.iter().zip(&other.exp).map(|(a, b)| (*a).max(*b)).collect() }
    }

    pub fn gcd_is_one(&self, other: &Monomial) -> bool {
        debug_assert_eq!(self.exp.len(), other.exp.len());
        self.exp.iter().zip(&other.exp).all(|(a, b)| *a == 0 || *b == 0)
    }
}

// Graded reverse lex.
// First compare total degree; tie-break by reverse lex (rightmost nonzero of
// (self - other) being negative => self is greater).
impl Ord for Monomial {
    fn cmp(&self, other: &Self) -> Ordering {
        let d1 = self.degree();
        let d2 = other.degree();
        if d1 != d2 { return d1.cmp(&d2); }
        for i in (0..self.exp.len()).rev() {
            match self.exp[i].cmp(&other.exp[i]) {
                Ordering::Equal => continue,
                ord => return ord.reverse(),
            }
        }
        Ordering::Equal
    }
}

impl PartialOrd for Monomial {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> { Some(self.cmp(other)) }
}

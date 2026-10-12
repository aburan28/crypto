//! The Hamming-weight base `{x : wt(x) ≤ w}` in the ONB, restricted to the
//! abscissae that lift, and its Frobenius orbits.
//!
//! Squaring is a coordinate permutation, so `x ↦ x²` maps the base to itself
//! (weight is invariant) and the liftable abscissae to liftable ones
//! (`Tr(c²) = Tr(c)` for `c = x + 1/x²`).  An orbit of the abscissae is
//! `{x, x², x⁴, …}`; its size divides `n`, and is `n` except for `x = 0` and
//! the words with a period dividing `n`.

use crate::curve::lift;
use crate::field::Onb;
use crate::targets::enumerate_weight;

pub struct Base {
    /// Liftable abscissae of weight `≤ w`, ascending.
    pub xs: Vec<u64>,
    /// A rational ordinate for each.
    pub ys: Vec<u64>,
    /// Words of weight `≤ w` in all (liftable or not).
    pub words: usize,
}

pub fn build_base(f: &Onb, w: usize) -> Base {
    let mut words: Vec<u64> = enumerate_weight(f.n, w)
        .into_iter()
        .map(|v| v as u64)
        .collect();
    let n_words = words.len();
    words.sort_unstable();
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for x in words {
        if let Some(y) = lift(f, x) {
            xs.push(x);
            ys.push(y);
        }
    }
    Base {
        xs,
        ys,
        words: n_words,
    }
}

pub fn index_of(b: &Base, x: u64) -> Option<usize> {
    b.xs.binary_search(&x).ok()
}

pub struct Orbits {
    /// Base indices of the representatives (smallest word of each orbit).
    pub reps: Vec<usize>,
    /// Orbit size of each representative.
    pub sizes: Vec<usize>,
}

pub fn orbits(f: &Onb, b: &Base) -> Orbits {
    let m = b.xs.len();
    let mut canon = vec![0u64; m];
    for (i, &x) in b.xs.iter().enumerate() {
        let mut best = x;
        let mut v = x;
        for _ in 1..f.n {
            v = f.sqr_raw(v);
            best = best.min(v);
        }
        canon[i] = best;
    }
    let mut reps = Vec::new();
    let mut sizes = Vec::new();
    for i in 0..m {
        if canon[i] == b.xs[i] {
            // Orbit size: first k >= 1 with x^(2^k) = x.
            let mut k = 1usize;
            let mut v = f.sqr_raw(b.xs[i]);
            while v != b.xs[i] {
                v = f.sqr_raw(v);
                k += 1;
            }
            reps.push(i);
            sizes.push(k);
        }
    }
    Orbits { reps, sizes }
}

impl Orbits {
    /// `Σ` of the orbit sizes of the representatives: the number of first
    /// summands the representatives stand for.
    pub fn covered(&self) -> usize {
        self.sizes.iter().sum()
    }

    /// The orbits smaller than `n`, as `(representative word, size)`: the
    /// boundary orbits whose size the covering count must use.
    pub fn short_orbits(&self, b: &Base, n: usize) -> Vec<(u64, usize)> {
        self.reps
            .iter()
            .zip(&self.sizes)
            .filter(|(_, &s)| s < n)
            .map(|(&i, &s)| (b.xs[i], s))
            .collect()
    }
}

/// Whether every member of every orbit of a base element is in the base
/// (the base is Frobenius-stable).
pub fn stable(f: &Onb, b: &Base) -> bool {
    b.xs.iter().all(|&x| index_of(b, f.sqr_raw(x)).is_some())
}

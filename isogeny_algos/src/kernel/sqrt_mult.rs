//! Algorithm 1 of Bernstein–De Feo–Leroux–Smith 2020 ("Faster computation of isogenies of large
//! prime degree", §3): evaluate h_S(alpha) = prod_{s in S} (alpha - zeta^s) for S = (I + J) ∪ K in
//! O~(max(#I, #J, #K)) field operations instead of #S.
//!
//! Lemma 3.4: if I and J have no common differences (i1 - i2 != j1 - j2 for i1 != i2, j1 != j2),
//! then (i, j) -> i + j is a bijection I x J -> I + J and
//!     h_{I+J}(alpha) = Res_Z(h_I(Z), H_J(alpha, Z)),  h_I(Z) = prod_I (Z - zeta^i),
//!                                                    H_J(alpha, Z) = prod_J (alpha - zeta^j Z).
//! h_I is monic, so the resultant is prod_{i in I} H_J(alpha, zeta^i): H_J is built by a product
//! tree (O~(#J)) and evaluated at the fixed points zeta^i by a remainder tree whose node
//! reciprocals are precomputed once per parameter set (O~(#I + #J) per alpha). h_K(alpha) is
//! the direct product (Algorithm 1, line 4).
//!
//! [`Alg1::eval_pollard`] is the §3.2 form of the same resultant: the fixed polynomial h_I is
//! evaluated at the alpha-dependent points alpha / zeta^j (Pollard's arrangement, with a
//! remainder tree in place of Bluestein's chirp-z transform), then multiplied by zeta^{j #I}.
//!
//! Index systems: [`Alg1::progression`] is Example 3.7 (S an arithmetic progression),
//! [`Alg1::for_set`] a greedy decomposition of an arbitrary finite S into translates of one
//! progression I = {0, r, ..., (b-1) r} (used for the cyclotomic example of §1:
//! Phi_n(X) = h_S(X) for S = (Z/nZ)^*, zeta a primitive n-th root of unity).
use crate::field::Field;
use crate::poly::{self, Poly};
use crate::strassen::RemTree;
use std::collections::BTreeSet;

/// zeta^e for a signed exponent (zeta must be invertible when e < 0).
pub fn zpow<F: Field>(f: &F, zeta: F::E, e: i64) -> F::E {
    if e >= 0 {
        f.pow(zeta, e as u128)
    } else {
        f.pow(f.inv(zeta), e.unsigned_abs() as u128)
    }
}

/// The direct product prod_{s in S} (alpha - zeta^s): #S multiplications (reference).
pub fn hs_direct<F: Field>(f: &F, zeta: F::E, s: &[i64], alpha: F::E) -> F::E {
    let mut acc = f.one();
    for &e in s {
        acc = f.mul(acc, f.sub(alpha, zpow(f, zeta, e)));
    }
    acc
}

/// Whether I and J have no common differences (Lemma 3.4's hypothesis): i1 - i2 != j1 - j2 for
/// all i1 != i2 in I and j1 != j2 in J. O(#I^2 + #J^2).
pub fn no_common_differences(i: &[i64], j: &[i64]) -> bool {
    let di: BTreeSet<i64> = i
        .iter()
        .flat_map(|&x| i.iter().filter(move |&&y| y != x).map(move |&y| x - y))
        .collect();
    j.iter()
        .all(|&x| j.iter().all(|&y| x == y || !di.contains(&(x - y))))
}

/// Parameters (zeta, I, J, K) of Algorithm 1 with the precomputed tree over {zeta^i : i in I}.
pub struct Alg1<F: Field> {
    pub zeta: F::E,
    pub i: Vec<i64>,
    pub j: Vec<i64>,
    pub k: Vec<i64>,
    /// zeta^j for j in J
    zj: Vec<F::E>,
    /// zeta^k for k in K
    zk: Vec<F::E>,
    /// remainder tree over the points zeta^i (its root is h_I)
    tree_i: Option<RemTree<F>>,
}

impl<F: Field> Alg1<F> {
    /// Parameters for S = (I + J) ∪ K. Panics unless I, J have no common differences and
    /// K ∩ (I + J) is empty (the hypotheses of Proposition 3.5).
    pub fn new(f: &F, zeta: F::E, i: Vec<i64>, j: Vec<i64>, k: Vec<i64>) -> Self {
        assert!(
            i.is_empty() || j.is_empty() || Self::check_sets(&i, &j, &k),
            "Alg1: I, J must have no common differences and K must avoid I + J"
        );
        let zj = j.iter().map(|&e| zpow(f, zeta, e)).collect();
        let zk = k.iter().map(|&e| zpow(f, zeta, e)).collect();
        let tree_i = if i.is_empty() || j.is_empty() {
            None
        } else {
            let pts: Vec<F::E> = i.iter().map(|&e| zpow(f, zeta, e)).collect();
            Some(RemTree::from_roots(f, &pts))
        };
        Alg1 {
            zeta,
            i,
            j,
            k,
            zj,
            zk,
            tree_i,
        }
    }

    fn check_sets(i: &[i64], j: &[i64], k: &[i64]) -> bool {
        // distinct elements, no common differences, and K disjoint from I + J
        let si: BTreeSet<i64> = i.iter().copied().collect();
        let sj: BTreeSet<i64> = j.iter().copied().collect();
        if si.len() != i.len() || sj.len() != j.len() || !no_common_differences(i, j) {
            return false;
        }
        let sum: BTreeSet<i64> = i
            .iter()
            .flat_map(|&x| j.iter().map(move |&y| x + y))
            .collect();
        k.iter().all(|e| !sum.contains(e))
    }

    /// Example 3.7: S = {m + k r : 0 <= k < n}, b = floor(sqrt n), I = {i r : i < b},
    /// J = {m + j b r : j < b}, K = {m + k r : b^2 <= k < n}.
    pub fn progression(f: &F, zeta: F::E, m: i64, r: i64, n: u64) -> Self {
        assert!(r != 0);
        let b = crate::strassen::isqrt(n) as i64;
        let i: Vec<i64> = (0..b).map(|t| t * r).collect();
        let j: Vec<i64> = (0..b).map(|t| m + t * b * r).collect();
        let k: Vec<i64> = (b * b..n as i64).map(|t| m + t * r).collect();
        Alg1::new(f, zeta, i, j, k)
    }

    /// A decomposition of an arbitrary finite set S (distinct integers): for each step r in
    /// `steps` (and b = floor(sqrt #S) shrinking by halves), I = {0, r, ..., (b-1) r} and J the
    /// starts of disjoint translates j + I contained in S, chosen greedily left to right with
    /// consecutive starts at least b r apart (so I, J have no common differences); K = the rest.
    /// The candidate with the smallest max(#I, #J, #K) wins.
    pub fn for_set(f: &F, zeta: F::E, s: &[i64], steps: &[i64]) -> Self {
        let set: BTreeSet<i64> = s.iter().copied().collect();
        assert_eq!(set.len(), s.len(), "Alg1::for_set: S has repeated elements");
        let mut best: Option<(usize, Vec<i64>, Vec<i64>, Vec<i64>)> = None;
        let b0 = crate::strassen::isqrt(s.len() as u64) as i64;
        for &r in steps {
            assert!(r > 0);
            let mut b = b0.max(1);
            while b >= 1 {
                let mut covered = BTreeSet::new();
                let mut j = vec![];
                let mut last: Option<i64> = None;
                for &x in &set {
                    if covered.contains(&x) {
                        continue;
                    }
                    if let Some(l) = last {
                        // |j1 - j2| >= b r rules out every difference of I, |i1 - i2| <= (b-1) r
                        if x - l < b * r {
                            continue;
                        }
                    }
                    if (0..b).all(|t| set.contains(&(x + t * r)) && !covered.contains(&(x + t * r)))
                    {
                        for t in 0..b {
                            covered.insert(x + t * r);
                        }
                        j.push(x);
                        last = Some(x);
                    }
                }
                let i: Vec<i64> = (0..b).map(|t| t * r).collect();
                let k: Vec<i64> = set
                    .iter()
                    .copied()
                    .filter(|x| !covered.contains(x))
                    .collect();
                let cost = i.len().max(j.len()).max(k.len());
                if best.as_ref().is_none_or(|bb| cost < bb.0) {
                    best = Some((cost, i, j, k));
                }
                if b == 1 {
                    break;
                }
                b /= 2;
            }
        }
        let (_, i, j, k) = best.expect("Alg1::for_set: no step given");
        if j.is_empty() {
            // nothing decomposed: everything is in K
            return Alg1::new(f, zeta, vec![], vec![], s.to_vec());
        }
        Alg1::new(f, zeta, i, j, k)
    }

    /// The exponent set S = (I + J) ∪ K, sorted.
    pub fn set(&self) -> Vec<i64> {
        let mut s: Vec<i64> = self
            .i
            .iter()
            .flat_map(|&x| self.j.iter().map(move |&y| x + y))
            .chain(self.k.iter().copied())
            .collect();
        s.sort_unstable();
        s
    }

    /// H_J(alpha, Z) = prod_{j in J} (alpha - zeta^j Z) by a product tree.
    fn h_j(&self, f: &F, alpha: F::E) -> Poly<F> {
        let layer: Vec<Poly<F>> = self.zj.iter().map(|&z| vec![alpha, f.neg(z)]).collect();
        product(f, layer)
    }

    /// h_{I+J}(alpha) = Res_Z(h_I, H_J(alpha, Z)) = prod_{i in I} H_J(alpha, zeta^i).
    pub fn h_ij(&self, f: &F, alpha: F::E) -> F::E {
        match &self.tree_i {
            None => f.one(),
            Some(t) => t
                .eval(f, &self.h_j(f, alpha))
                .into_iter()
                .fold(f.one(), |a, v| f.mul(a, v)),
        }
    }

    /// h_K(alpha) = prod_{k in K} (alpha - zeta^k).
    pub fn h_k(&self, f: &F, alpha: F::E) -> F::E {
        self.zk
            .iter()
            .fold(f.one(), |a, &z| f.mul(a, f.sub(alpha, z)))
    }

    /// h_S(alpha) by Algorithm 1.
    pub fn eval(&self, f: &F, alpha: F::E) -> F::E {
        f.mul(self.h_ij(f, alpha), self.h_k(f, alpha))
    }

    /// h_S(alpha) by the §3.2 arrangement: prod_j zeta^{j #I} h_I(alpha / zeta^j), with h_I the
    /// root of the precomputed tree and a fresh remainder tree over the points alpha / zeta^j.
    /// Needs zeta invertible.
    pub fn eval_pollard(&self, f: &F, alpha: F::E) -> F::E {
        let hij = match &self.tree_i {
            None => f.one(),
            Some(t) => {
                let hi = t.root();
                let zji = f.batch_inv(&self.zj);
                let pts: Vec<F::E> = zji.iter().map(|&z| f.mul(alpha, z)).collect();
                let vals = RemTree::from_roots(f, &pts).eval(f, hi);
                let ni = self.i.len() as u128;
                vals.into_iter()
                    .zip(&self.zj)
                    .fold(f.one(), |a, (v, &z)| f.mul(a, f.mul(v, f.pow(z, ni))))
            }
        };
        f.mul(hij, self.h_k(f, alpha))
    }
}

/// Product of the polynomials by a balanced tree.
pub fn product<F: Field>(f: &F, mut layer: Vec<Poly<F>>) -> Poly<F> {
    while layer.len() > 1 {
        layer = layer
            .chunks(2)
            .map(|c| {
                if c.len() == 2 {
                    poly::mul_raw(f, &c[0], &c[1])
                } else {
                    c[0].clone()
                }
            })
            .collect();
    }
    layer.pop().unwrap_or_else(|| vec![f.one()])
}

/// Integer coefficients of the cyclotomic polynomial Phi_n (low to high), by exact division of
/// X^n - 1 by Phi_d for the proper divisors d of n (independent of Algorithm 1; n small).
pub fn cyclotomic_int(n: u64) -> Vec<i128> {
    let n = n as usize;
    let mut phis: Vec<Option<Vec<i128>>> = vec![None; n + 1];
    for d in 1..=n {
        if !n.is_multiple_of(d) {
            continue;
        }
        // X^d - 1 divided by Phi_e for e | d, e < d
        let mut num = vec![0i128; d + 1];
        num[0] = -1;
        num[d] = 1;
        for e in 1..d {
            if d % e == 0 {
                let den = phis[e].as_ref().expect("divisor computed first");
                num = exact_div(&num, den);
            }
        }
        phis[d] = Some(num);
    }
    phis[n].take().unwrap()
}

fn exact_div(a: &[i128], b: &[i128]) -> Vec<i128> {
    // b monic
    let (da, db) = (a.len() - 1, b.len() - 1);
    let mut r = a.to_vec();
    let mut q = vec![0i128; da - db + 1];
    for t in (0..=da - db).rev() {
        let c = r[t + db];
        q[t] = c;
        for s in 0..=db {
            r[t + s] -= c * b[s];
        }
    }
    debug_assert!(r.iter().all(|&x| x == 0), "not an exact division");
    q
}

//! 2^e-isogeny chains on Montgomery curves in x-only projective form, as in SIDH/SIKE (Costello–
//! Longa–Naehrig 2016; SIKE specification formulas): curves as (A24+ : C24) = (A + 2C : 4C),
//! 4-isogenies (and 2-isogenies) without inversions, and the kernel points found by an optimal
//! strategy (dynamic programming over "multiply by 4" and "push a point" costs).
use crate::field::Field;
use crate::kernel::montgomery::{xdbl_p, Proj24, XZ};

/// A 4-isogeny with kernel <P4> (P4 of order 4 with 2 P4 != (0,0)): codomain and the constants
/// (K1, K2, K3) of the evaluation.
#[derive(Clone, Copy, Debug)]
pub struct Four<E> {
    pub cod: Proj24<E>,
    pub k: [E; 3],
}

pub fn four_isogeny<F: Field>(f: &F, p4: XZ<F::E>) -> Option<Four<F::E>> {
    let (x, z) = p4;
    let k2 = f.sub(x, z);
    let k3 = f.add(x, z);
    if f.is_zero(k2) || f.is_zero(k3) {
        return None; // 2 P4 = (0, 0): use two 2-isogenies instead
    }
    let mut k1 = f.sq(z);
    k1 = f.add(k1, k1);
    let c24 = f.sq(k1);
    k1 = f.add(k1, k1);
    let mut a24 = f.sq(x);
    a24 = f.add(a24, a24);
    a24 = f.sq(a24);
    Some(Four { cod: (a24, c24), k: [k1, k2, k3] })
}

impl<E: Copy> Four<E> {
    pub fn eval<F: Field<E = E>>(&self, f: &F, p: XZ<E>) -> XZ<E> {
        let [k1, k2, k3] = self.k;
        let mut t0 = f.add(p.0, p.1);
        let mut t1 = f.sub(p.0, p.1);
        let x = f.mul(t0, k2);
        let z = f.mul(t1, k3);
        t0 = f.mul(f.mul(t0, t1), k1);
        t1 = f.sq(f.add(x, z));
        let zz = f.sq(f.sub(x, z));
        let xs = f.add(t0, t1);
        let t0b = f.sub(zz, t0);
        (f.mul(xs, t1), f.mul(zz, t0b))
    }
}

/// A 2-isogeny with kernel <P2> (P2 of order 2, P2 != (0,0)).
#[derive(Clone, Copy, Debug)]
pub struct Two<E> {
    pub cod: Proj24<E>,
    pub k: XZ<E>,
}

pub fn two_isogeny<F: Field>(f: &F, p2: XZ<F::E>) -> Option<Two<F::E>> {
    if f.is_zero(p2.0) {
        return None;
    }
    let a24 = f.sq(p2.0);
    let c24 = f.sq(p2.1);
    Some(Two { cod: (f.sub(c24, a24), c24), k: p2 })
}

impl<E: Copy> Two<E> {
    pub fn eval<F: Field<E = E>>(&self, f: &F, q: XZ<E>) -> XZ<E> {
        let (x2, z2) = self.k;
        let t0 = f.mul(f.add(x2, z2), f.sub(q.0, q.1));
        let t1 = f.mul(f.sub(x2, z2), f.add(q.0, q.1));
        (f.mul(q.0, f.add(t0, t1)), f.mul(q.1, f.sub(t0, t1)))
    }
}

/// [2^e] P.
pub fn xdbl_e<F: Field>(f: &F, k: Proj24<F::E>, mut p: XZ<F::E>, e: usize) -> XZ<F::E> {
    for _ in 0..e {
        p = xdbl_p(f, k, p);
    }
    p
}

/// Optimal strategy for n steps: split[n] = number i of steps done on the kernel [m^(n-i)] R
/// before continuing with the image of R, minimising S(n) = S(i) + S(n-i) + (n-i) mul + i eval.
pub fn optimal_splits(n: usize, mul: f64, eval: f64) -> Vec<usize> {
    let mut s = vec![0.0f64; n + 1];
    let mut split = vec![0usize; n + 1];
    for m in 2..=n {
        let mut best = f64::INFINITY;
        for i in 1..m {
            let c = s[i] + s[m - i] + (m - i) as f64 * mul + i as f64 * eval;
            if c < best {
                best = c;
                split[m] = i;
            }
        }
        s[m] = best;
    }
    split
}

#[derive(Clone, Copy, Debug, Default)]
pub struct TwoPowerStats {
    pub doublings: usize,
    pub evals: usize,
    pub isogenies: usize,
}

/// The 2^e-isogeny with kernel <R> (R of order exactly 2^e, e even, and no step's 4-torsion
/// point doubling to (0,0)) as e/2 four-isogenies with the given splits; `pts` are pushed through.
/// Returns the codomain, or None on a degenerate step.
pub fn four_chain<F: Field>(
    f: &F,
    curve: Proj24<F::E>,
    r: XZ<F::E>,
    e: usize,
    pts: &mut Vec<XZ<F::E>>,
    splits: &[usize],
    st: &mut TwoPowerStats,
) -> Option<Proj24<F::E>> {
    assert!(e % 2 == 0);
    fn rec<F: Field>(
        f: &F,
        curve: &mut Proj24<F::E>,
        r: XZ<F::E>,
        n: usize,
        stack: &mut Vec<XZ<F::E>>,
        splits: &[usize],
        st: &mut TwoPowerStats,
    ) -> Option<()> {
        if n == 1 {
            let iso = four_isogeny(f, r)?;
            for p in stack.iter_mut() {
                *p = iso.eval(f, *p);
            }
            st.evals += stack.len();
            st.isogenies += 1;
            *curve = iso.cod;
            return Some(());
        }
        let i = splits[n];
        let s = xdbl_e(f, *curve, r, 2 * (n - i));
        st.doublings += 2 * (n - i);
        stack.push(r);
        rec(f, curve, s, i, stack, splits, st)?;
        let r2 = stack.pop().unwrap();
        rec(f, curve, r2, n - i, stack, splits, st)
    }
    let mut c = curve;
    rec(f, &mut c, r, e / 2, pts, splits, st)?;
    Some(c)
}

/// The same 2^e-isogeny as e two-isogenies (strategy over "double" / "push").
pub fn two_chain<F: Field>(
    f: &F,
    curve: Proj24<F::E>,
    r: XZ<F::E>,
    e: usize,
    pts: &mut Vec<XZ<F::E>>,
    splits: &[usize],
    st: &mut TwoPowerStats,
) -> Option<Proj24<F::E>> {
    fn rec<F: Field>(
        f: &F,
        curve: &mut Proj24<F::E>,
        r: XZ<F::E>,
        n: usize,
        stack: &mut Vec<XZ<F::E>>,
        splits: &[usize],
        st: &mut TwoPowerStats,
    ) -> Option<()> {
        if n == 1 {
            let iso = two_isogeny(f, r)?;
            for p in stack.iter_mut() {
                *p = iso.eval(f, *p);
            }
            st.evals += stack.len();
            st.isogenies += 1;
            *curve = iso.cod;
            return Some(());
        }
        let i = splits[n];
        let s = xdbl_e(f, *curve, r, n - i);
        st.doublings += n - i;
        stack.push(r);
        rec(f, curve, s, i, stack, splits, st)?;
        let r2 = stack.pop().unwrap();
        rec(f, curve, r2, n - i, stack, splits, st)
    }
    let mut c = curve;
    rec(f, &mut c, r, e, pts, splits, st)?;
    Some(c)
}

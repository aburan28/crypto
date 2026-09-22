//! Quadratic Semaev decomposition oracle (`m = 2`) on packed Gray FES.
//!
//! [`crate::cryptanalysis::mq_fes::mq_fes_decompose_reference`] rebuilds the
//! Weil-restricted `S₃` system symbolically for every target and walks all
//! `2^{2ℓ}` assignments before lifting.  This oracle keeps the answer set and
//! changes three costs:
//!
//! 1. **Packed template.**  `S₃(x₁, x₂, x_R)` is affine in the bits of `x_R`
//!    (squaring is `F₂`-linear), so the system is `C + Σ_k r_k A_k` for
//!    target-independent quadratic systems `C, A_k`.  Each is stored once per
//!    factor base as one `u64` per monomial (bit `e` = equation `e`); a target
//!    costs at most `n` table XORs instead of a symbolic rebuild.
//! 2. **Swap symmetry.**  `S₃` depends on `x₁ + x₂` and `x₁x₂` only, so with
//!    `x₁ = Σ aᵢvᵢ`, `x₂ = Σ bᵢvᵢ` every equation is invariant under `a ↔ b`.
//!    Substitute `s = a ⊕ b`.  Unordered pairs with `s ≠ 0` are represented
//!    exactly once by `s_{<j} = 0, s_j = 1, a_j = 0` for `j` the lowest set bit
//!    of `s`; `s = 0` is the diagonal.  The case split walks
//!    `2^ℓ + Σ_j 2^{2ℓ-2-j} ≈ 2^{2ℓ-1}` points instead of `2^{2ℓ}`.
//! 3. **Lift on the fly.**  Each root is checked and lifted when the walk
//!    reaches it, so a decomposable target stops at its first liftable root
//!    and there is no root cap.
//!
//! The walk itself is the register-blocked libfes `L = 4` chunk.  None of
//! this moves the free-oracle floor: the question asked of the oracle is
//! unchanged.

use super::koblitz_groebner::{DecompositionSystem, FieldStructure};
use super::mq_fes::{idxq, Ffs};
use super::polynomial_reuse::DecompositionTemplate;
use crate::binary_ecc::{BinaryPoint, F2mElement};
use num_bigint::BigUint;
use std::collections::HashMap;
use std::sync::{Arc, Mutex, OnceLock};

/// A system of up to 64 quadratic Boolean equations, bit-sliced by equation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PackedQuad {
    pub n: usize,
    /// `lin[0]` is the constant, `lin[1 + i]` the coefficient of `xᵢ`.
    pub lin: Vec<u64>,
    /// `quad[idxq(j, i)]`, `j < i`, the coefficient of `xⱼxᵢ`.
    pub quad: Vec<u64>,
}

/// Image of an old variable under an affine substitution: `c ⊕ Σ_{p ∈ vars} y_p`.
#[derive(Clone, Copy, Debug)]
struct Affine {
    c: bool,
    vars: u64,
}

impl PackedQuad {
    pub fn zero(n: usize) -> Self {
        Self {
            n,
            lin: vec![0; n + 1],
            quad: vec![0; n * n.saturating_sub(1) / 2],
        }
    }

    /// Pack polynomials of degree ≤ 2, one per equation; `None` past degree 2.
    pub fn from_polys(polys: &[super::pq_groebner_f2::F2BoolPoly], n: usize) -> Option<Self> {
        if polys.len() > 64 || n > 64 {
            return None;
        }
        let mut out = Self::zero(n);
        for (eq, p) in polys.iter().enumerate() {
            let bit = 1u64 << eq;
            for t in &p.terms {
                let mask = t.mask;
                match mask.count_ones() {
                    0 => out.lin[0] ^= bit,
                    1 => {
                        let i = mask.trailing_zeros() as usize;
                        if i >= n {
                            return None;
                        }
                        out.lin[1 + i] ^= bit;
                    }
                    2 => {
                        let j = mask.trailing_zeros() as usize;
                        let i = 63 - mask.leading_zeros() as usize;
                        if i >= n {
                            return None;
                        }
                        out.quad[idxq(j, i)] ^= bit;
                    }
                    _ => return None,
                }
            }
        }
        Some(out)
    }

    pub fn xor_assign(&mut self, other: &Self) {
        debug_assert_eq!(self.n, other.n);
        for (a, b) in self.lin.iter_mut().zip(&other.lin) {
            *a ^= b;
        }
        for (a, b) in self.quad.iter_mut().zip(&other.quad) {
            *a ^= b;
        }
    }

    /// Equation values at `x`; zero iff `x` is a common root.
    pub fn eval(&self, x: u64) -> u64 {
        let mut v = self.lin[0];
        for i in 0..self.n {
            if (x >> i) & 1 == 0 {
                continue;
            }
            v ^= self.lin[1 + i];
            let row = idxq(0, i.max(1));
            for j in 0..i {
                if (x >> j) & 1 == 1 {
                    v ^= self.quad[row + j];
                }
            }
        }
        v
    }

    /// The system in the new variables `y` given `x_old[i] = map[i]`.
    fn substitute(&self, map: &[Affine], n_new: usize) -> Self {
        debug_assert_eq!(map.len(), self.n);
        let mut out = Self::zero(n_new);
        out.lin[0] = self.lin[0];
        let add_var = |out: &mut Self, vars: u64, w: u64| {
            let mut v = vars;
            while v != 0 {
                let p = v.trailing_zeros() as usize;
                out.lin[1 + p] ^= w;
                v &= v - 1;
            }
        };
        for i in 0..self.n {
            let w = self.lin[1 + i];
            if w == 0 {
                continue;
            }
            if map[i].c {
                out.lin[0] ^= w;
            }
            add_var(&mut out, map[i].vars, w);
        }
        for i in 1..self.n {
            let row = idxq(0, i);
            for j in 0..i {
                let w = self.quad[row + j];
                if w == 0 {
                    continue;
                }
                let (mi, mj) = (map[i], map[j]);
                if mi.c && mj.c {
                    out.lin[0] ^= w;
                }
                if mi.c {
                    add_var(&mut out, mj.vars, w);
                }
                if mj.c {
                    add_var(&mut out, mi.vars, w);
                }
                let mut vp = mi.vars;
                while vp != 0 {
                    let p = vp.trailing_zeros() as usize;
                    vp &= vp - 1;
                    let mut vq = mj.vars;
                    while vq != 0 {
                        let q = vq.trailing_zeros() as usize;
                        vq &= vq - 1;
                        if p == q {
                            out.lin[1 + p] ^= w;
                        } else {
                            out.quad[idxq(p.min(q), p.max(q))] ^= w;
                        }
                    }
                }
            }
        }
        out
    }

    /// libfes `Fq`/`Fl` tables with the two fictive variables padded in.
    fn to_tables(&self, fq: &mut [u64; 561], fl: &mut [u64; 34]) {
        let n = self.n;
        fq.fill(0);
        fl.fill(0);
        fl[..=n].copy_from_slice(&self.lin);
        fq[..self.quad.len()].copy_from_slice(&self.quad);
        for i in 0..n {
            fq[idxq(i, n)] = 0;
        }
        fq[idxq(0, n + 1)] = 0;
        for i in 1..n {
            fq[idxq(i, n + 1)] = fq[idxq(i - 1, i)];
        }
        fq[idxq(n, n + 1)] = 0;
    }
}

/// `C` and `A_k` of the affine-in-`x_R` Semaev system, packed.
#[derive(Debug)]
pub struct PackedSemaevTemplate {
    pub ell: usize,
    pub n_vars: usize,
    pub constant: PackedQuad,
    pub coefficients: Vec<PackedQuad>,
}

impl PackedSemaevTemplate {
    pub fn build(basis: &[F2mElement], b: &F2mElement, st: &FieldStructure) -> Option<Self> {
        let t = DecompositionTemplate::build(basis, b, 2, st)?;
        if !t.prefix.is_empty() {
            return None;
        }
        let constant = PackedQuad::from_polys(&t.constant, t.n_vars)?;
        let coefficients = t
            .coefficients
            .iter()
            .map(|c| PackedQuad::from_polys(c, t.n_vars))
            .collect::<Option<Vec<_>>>()?;
        // The case split is only complete for swap-invariant systems.
        let ell = t.ell;
        let swap: Vec<Affine> = (0..2 * ell)
            .map(|i| Affine {
                c: false,
                vars: 1u64 << ((i + ell) % (2 * ell)),
            })
            .collect();
        if t.n_vars != 2 * ell
            || std::iter::once(&constant)
                .chain(&coefficients)
                .any(|q| q.substitute(&swap, 2 * ell) != *q)
        {
            return None;
        }
        Some(Self {
            ell: t.ell,
            n_vars: t.n_vars,
            constant,
            coefficients,
        })
    }

    pub fn instantiate(&self, x_r_bits: u64) -> PackedQuad {
        let mut sys = self.constant.clone();
        let mut bits = x_r_bits;
        while bits != 0 {
            let k = bits.trailing_zeros() as usize;
            bits &= bits - 1;
            if let Some(a) = self.coefficients.get(k) {
                sys.xor_assign(a);
            }
        }
        sys
    }
}

fn template_cache() -> &'static Mutex<HashMap<Vec<u64>, Arc<PackedSemaevTemplate>>> {
    static CACHE: OnceLock<Mutex<HashMap<Vec<u64>, Arc<PackedSemaevTemplate>>>> = OnceLock::new();
    CACHE.get_or_init(|| Mutex::new(HashMap::new()))
}

fn bits(e: &F2mElement) -> u64 {
    e.raw_bits().first().copied().unwrap_or(0)
}

/// The packed template for this factor base, built on first use.
///
/// The key fixes the field (degree and `z^n mod f`), `b` and the ordered
/// subspace basis, which is everything [`DecompositionTemplate`] reads.
pub fn packed_template(
    basis: &[F2mElement],
    b: &F2mElement,
    st: &FieldStructure,
) -> Option<Arc<PackedSemaevTemplate>> {
    let n = st.n as usize;
    let mut key = Vec::with_capacity(basis.len() + 3 + n);
    key.push(st.n as u64);
    key.push(bits(b));
    key.extend(basis.iter().map(bits));
    key.extend(st.squares.iter().copied());
    if let Some(row) = st.reduced.last() {
        key.extend(row.iter().copied());
    }
    if let Some(t) = template_cache().lock().ok()?.get(&key) {
        return Some(t.clone());
    }
    let t = Arc::new(PackedSemaevTemplate::build(basis, b, st)?);
    template_cache().lock().ok()?.insert(key, t.clone());
    Some(t)
}

/// Visit every common root of `sys` in Gray order until `visit` returns
/// `true`.  Returns whether it stopped, and the assignments walked.
fn walk_visit(sys: &PackedQuad, visit: &mut dyn FnMut(u64) -> bool) -> (bool, u64) {
    let n = sys.n;
    if n < 4 {
        for x in 0..(1u64 << n) {
            if sys.eval(x) == 0 && visit(x) {
                return (true, x + 1);
            }
        }
        return (false, 1u64 << n);
    }
    let mut fq = [0u64; 561];
    let mut fl = [0u64; 34];
    sys.to_tables(&mut fq, &mut fl);
    walk_l4_visit(&fq, &mut fl, n, visit)
}

#[cold]
#[inline(never)]
fn hit(visit: &mut dyn FnMut(u64) -> bool, index: u64) -> bool {
    visit(index ^ (index >> 1))
}

/// Register-blocked `L = 4` Gray walk (as `gray_ffs_l4_find_one`) with a visitor.
fn walk_l4_visit(
    fq: &[u64; 561],
    fl: &mut [u64; 34],
    n: usize,
    visit: &mut dyn FnMut(u64) -> bool,
) -> (bool, u64) {
    const L: usize = 4;
    let mut ffs = Ffs::reset(n - L);
    let mut k1 = ffs.k1 + L as i32;
    let mut k2;
    let iterations = 1u64 << (n - L);
    let mut f0 = fl[0];
    let mut f1 = fl[1];
    let mut f2 = fl[2];
    let mut f3 = fl[3];
    let mut f4 = fl[4];
    for j in 0..iterations {
        let alpha = idxq(0, k1 as usize);
        ffs.step();
        k1 = ffs.k1 + L as i32;
        k2 = ffs.k2 + L as i32;
        let beta = (1 + k1) as usize;
        let gamma = idxq(k1 as usize, k2 as usize);
        let base = j << L;
        // Safety: as in `gray_ffs_l4_find_one` — alpha + 3 < 561, beta < 34.
        unsafe {
            let qa0 = *fq.get_unchecked(alpha);
            let qa1 = *fq.get_unchecked(alpha + 1);
            let qa2 = *fq.get_unchecked(alpha + 2);
            let qa3 = *fq.get_unchecked(alpha + 3);
            let q0 = *fq.get_unchecked(0);
            let q1 = *fq.get_unchecked(1);
            let q2 = *fq.get_unchecked(2);
            let q3 = *fq.get_unchecked(3);
            let q4 = *fq.get_unchecked(4);
            let q5 = *fq.get_unchecked(5);
            macro_rules! step {
                ($fa:ident, $qb:expr, $idx:expr) => {{
                    if f0 == 0 && hit(visit, $idx) {
                        return (true, $idx + 1);
                    }
                    $fa ^= $qb;
                    f0 ^= $fa;
                }};
            }
            step!(f1, qa0, base);
            step!(f2, qa1, base + 1);
            step!(f1, q0, base + 2);
            step!(f3, qa2, base + 3);
            step!(f1, q1, base + 4);
            step!(f2, q2, base + 5);
            step!(f1, q0, base + 6);
            step!(f4, qa3, base + 7);
            step!(f1, q3, base + 8);
            step!(f2, q4, base + 9);
            step!(f1, q0, base + 10);
            step!(f3, q5, base + 11);
            step!(f1, q1, base + 12);
            step!(f2, q2, base + 13);
            step!(f1, q0, base + 14);
            if f0 == 0 && hit(visit, base + 15) {
                return (true, base + 16);
            }
            let fb = fl.get_unchecked_mut(beta);
            *fb ^= *fq.get_unchecked(gamma);
            f0 ^= *fb;
        }
    }
    (false, iterations << L)
}

/// One sub-walk of the swap-symmetric case split.
struct SymCase {
    sys: PackedQuad,
    /// For each new variable, its role: `(is_a, coordinate)`.
    roles: Vec<(bool, usize)>,
    /// Fixed part of `s` (the single set bit `j`), zero on the diagonal.
    s_fixed: u64,
}

impl SymCase {
    /// Old assignment `a | b << ℓ` for new-variable assignment `y`.
    fn old_point(&self, y: u64, ell: usize) -> u64 {
        let mut a = 0u64;
        let mut s = self.s_fixed;
        for (p, &(is_a, i)) in self.roles.iter().enumerate() {
            if (y >> p) & 1 == 1 {
                if is_a {
                    a |= 1 << i;
                } else {
                    s |= 1 << i;
                }
            }
        }
        a | ((a ^ s) << ell)
    }
}

/// The `ℓ + 1` sub-systems covering every unordered pair `{a, b}` once.
fn symmetric_cases(sys: &PackedQuad, ell: usize) -> Vec<SymCase> {
    debug_assert_eq!(sys.n, 2 * ell);
    let mut cases = Vec::with_capacity(ell + 1);
    // Diagonal: b = a.
    {
        let mut map = vec![Affine { c: false, vars: 0 }; 2 * ell];
        for i in 0..ell {
            map[i] = Affine { c: false, vars: 1 << i };
            map[ell + i] = Affine { c: false, vars: 1 << i };
        }
        cases.push(SymCase {
            sys: sys.substitute(&map, ell),
            roles: (0..ell).map(|i| (true, i)).collect(),
            s_fixed: 0,
        });
    }
    for j in 0..ell {
        let mut roles = Vec::with_capacity(2 * ell - 2 - j);
        let mut a_var = vec![None; ell];
        for i in (0..ell).filter(|&i| i != j) {
            a_var[i] = Some(roles.len());
            roles.push((true, i));
        }
        let mut s_var = vec![None; ell];
        for i in j + 1..ell {
            s_var[i] = Some(roles.len());
            roles.push((false, i));
        }
        let mut map = vec![Affine { c: false, vars: 0 }; 2 * ell];
        for i in 0..ell {
            let a = a_var[i].map_or(0, |p| 1u64 << p);
            map[i] = Affine { c: false, vars: a };
            // b_i = a_i ⊕ s_i with s_{<j} = 0, s_j = 1, a_j = 0.
            map[ell + i] = match i.cmp(&j) {
                std::cmp::Ordering::Less => Affine { c: false, vars: a },
                std::cmp::Ordering::Equal => Affine { c: true, vars: 0 },
                std::cmp::Ordering::Greater => Affine {
                    c: false,
                    vars: a | (1u64 << s_var[i].unwrap()),
                },
            };
        }
        let n_new = roles.len();
        cases.push(SymCase {
            sys: sys.substitute(&map, n_new),
            roles,
            s_fixed: 1 << j,
        });
    }
    cases
}

/// Visit each unordered common root `{a, b}` of a swap-symmetric `m = 2`
/// system once, as the old assignment `a | b << ℓ`.  Returns whether the
/// visitor stopped the search and the assignments walked.
pub fn symmetric_pair_visit(
    sys: &PackedQuad,
    ell: usize,
    visit: &mut dyn FnMut(u64) -> bool,
) -> (bool, u64) {
    let mut walked = 0u64;
    for case in symmetric_cases(sys, ell) {
        let mut inner = |y: u64| visit(case.old_point(y, ell));
        let (stopped, w) = walk_visit(&case.sys, &mut inner);
        walked += w;
        if stopped {
            return (true, walked);
        }
    }
    (false, walked)
}

/// The quadratic Semaev decomposition oracle: same contract and answer set
/// as [`crate::cryptanalysis::mq_fes::mq_fes_decompose_reference`].
pub fn mq_fes_decompose(
    kc: &crate::cryptanalysis::koblitz_index_calculus::KoblitzCurve,
    fb: &crate::cryptanalysis::koblitz_index_calculus::FrobeniusFactorBase,
    index_of: &HashMap<(BigUint, BigUint), usize>,
    st: &FieldStructure,
    target: &BinaryPoint,
    m: usize,
) -> (
    Option<Vec<usize>>,
    crate::cryptanalysis::koblitz_index_calculus::SatDecompositionStats,
) {
    use super::mq_fes::profile;
    use crate::cryptanalysis::koblitz_index_calculus::{lift_candidate, SatDecompositionStats};

    let mut stats = SatDecompositionStats::default();
    if m != 2 {
        stats.exhausted = true;
        return (None, stats);
    }
    let x_r = match target {
        BinaryPoint::Affine { x, .. } => x.clone(),
        BinaryPoint::Infinity => return (None, stats),
    };
    let ell = fb.subspace_basis.len();
    // Sub-walks have at most 2ℓ − 1 variables; the Gray tables stop at 32.
    if 2 * ell > 33 || st.n > 64 {
        return super::mq_fes::mq_fes_decompose_reference(kc, fb, index_of, st, target, m);
    }
    profile::add_calls(1);
    let build_start = std::time::Instant::now();
    let Some(template) = packed_template(&fb.subspace_basis, &kc.curve.b, st) else {
        return super::mq_fes::mq_fes_decompose_reference(kc, fb, index_of, st, target, m);
    };
    let sys = template.instantiate(bits(&x_r));
    profile::add_build_ns(build_start.elapsed().as_nanos() as u64);
    stats.solver_calls = 1;

    let shape = DecompositionSystem {
        equations: Vec::new(),
        n_vars: 2 * ell,
        ell,
        m: 2,
    };
    let mut found = None;
    let mut roots = 0u64;
    let mut lift_ns = 0u64;
    let walk_start = std::time::Instant::now();
    let mut visit = |root: u64| {
        let t = std::time::Instant::now();
        roots += 1;
        stats.models += 1;
        let done = if sys.eval(root) != 0 {
            stats.spurious += 1;
            false
        } else {
            let xs: Vec<_> = (0..2)
                .map(|i| shape.summand_x(&fb.subspace_basis, root, i, kc.n))
                .collect();
            found = lift_candidate(kc, fb, index_of, &xs, target);
            found.is_some()
        };
        lift_ns += t.elapsed().as_nanos() as u64;
        done
    };
    let (_, walked) = symmetric_pair_visit(&sys, ell, &mut visit);
    let total = walk_start.elapsed().as_nanos() as u64;
    profile::add_walk(walked, total.saturating_sub(lift_ns), roots);
    profile::add_lift_ns(lift_ns);

    if found.is_some() {
        return (found, stats);
    }
    if stats.models == 0 {
        stats.refuted = true;
    } else {
        stats.exhausted = true;
    }
    (None, stats)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::mq_fes::{fill_fq_fl, gray_incremental_find_all, QuadraticForm};
    use crate::cryptanalysis::pq_groebner_f2::{F2BoolMono, F2BoolPoly};

    fn rng(seed: &mut u64) -> u64 {
        *seed ^= *seed << 13;
        *seed ^= *seed >> 7;
        *seed ^= *seed << 17;
        *seed
    }

    /// A random system invariant under swapping the halves `a ↔ b`.
    fn random_symmetric(ell: usize, m: usize, seed: &mut u64) -> Vec<F2BoolPoly> {
        let n = 2 * ell;
        let swap = |mask: u64| {
            let lo = mask & ((1 << ell) - 1);
            let hi = mask >> ell;
            hi | (lo << ell)
        };
        (0..m)
            .map(|_| {
                let mut monos = std::collections::BTreeSet::new();
                let mut toggle = |mask: u64| {
                    if !monos.insert(mask) {
                        monos.remove(&mask);
                    }
                };
                for mask in 0..(1u64 << n) {
                    if mask.count_ones() > 2 || mask > swap(mask) {
                        continue;
                    }
                    if rng(seed) % 3 == 0 {
                        toggle(mask);
                        if swap(mask) != mask {
                            toggle(swap(mask));
                        }
                    }
                }
                F2BoolPoly::from_monos(
                    monos.into_iter().map(F2BoolMono::from_mask).collect(),
                    n,
                )
            })
            .collect()
    }

    fn canonical(x: u64, ell: usize) -> (u64, u64) {
        let a = x & ((1 << ell) - 1);
        let b = x >> ell;
        (a.min(b), a.max(b))
    }

    #[test]
    fn packed_matches_quadratic_form_tables() {
        let mut seed = 0x1234_5678_9abc_def1u64;
        let polys = random_symmetric(4, 10, &mut seed);
        let packed = PackedQuad::from_polys(&polys, 8).unwrap();
        let forms: Vec<_> = polys.iter().map(|p| QuadraticForm::from_poly(p).unwrap()).collect();
        let (mut fq_a, mut fl_a) = ([0u64; 561], [0u64; 34]);
        let (mut fq_b, mut fl_b) = ([0u64; 561], [0u64; 34]);
        packed.to_tables(&mut fq_a, &mut fl_a);
        fill_fq_fl(&forms, 8, &mut fq_b, &mut fl_b);
        assert_eq!(fq_a, fq_b);
        assert_eq!(fl_a, fl_b);
        for x in 0..256u64 {
            let want = forms
                .iter()
                .enumerate()
                .fold(0u64, |acc, (e, f)| acc | ((f.eval(x) as u64) << e));
            assert_eq!(packed.eval(x), want);
        }
    }

    #[test]
    fn symmetric_split_covers_each_unordered_root_once() {
        let mut seed = 0xdead_beef_cafe_f00du64;
        for ell in 1..=7 {
            for _ in 0..6 {
                // Few equations so there are many roots to cover.
                let polys = random_symmetric(ell, 3, &mut seed);
                let packed = PackedQuad::from_polys(&polys, 2 * ell).unwrap();
                let forms: Vec<_> =
                    polys.iter().map(|p| QuadraticForm::from_poly(p).unwrap()).collect();
                let full = gray_incremental_find_all(&forms, usize::MAX).unwrap();
                let mut want: Vec<_> = full.iter().map(|&x| canonical(x, ell)).collect();
                want.sort_unstable();
                want.dedup();

                let mut got = Vec::new();
                let (stopped, walked) = symmetric_pair_visit(&packed, ell, &mut |x| {
                    assert_eq!(packed.eval(x), 0, "visited a non-root");
                    got.push(canonical(x, ell));
                    false
                });
                assert!(!stopped);
                let expect_walk = (1u64 << ell)
                    + (0..ell).map(|j| 1u64 << (2 * ell - 2 - j)).sum::<u64>();
                assert_eq!(walked, expect_walk);
                let unsorted = got.len();
                got.sort_unstable();
                got.dedup();
                assert_eq!(got.len(), unsorted, "an unordered root was visited twice");
                assert_eq!(got, want, "ell={ell}");
            }
        }
    }

    #[test]
    fn template_instantiation_matches_symbolic_rebuild() {
        use crate::cryptanalysis::koblitz_groebner::build_decomposition_system;
        use crate::cryptanalysis::koblitz_index_calculus::{build_frobenius_factor_base, KoblitzCurve};
        for (a, n) in [(1u8, 7u32), (0, 13), (1, 17)] {
            let kc = KoblitzCurve::new(a, n).expect("curve");
            let fb = build_frobenius_factor_base(&kc, 0).expect("factor base");
            let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
            let t = packed_template(&fb.subspace_basis, &kc.curve.b, &st).expect("template");
            let mut seed = 0x5eed_0000 ^ n as u64;
            for _ in 0..40 {
                let r = rng(&mut seed) & ((1u64 << n) - 1);
                let x_r = F2mElement::from_bit_positions(
                    &(0..n).filter(|k| (r >> k) & 1 == 1).collect::<Vec<_>>(),
                    n,
                );
                let sys = build_decomposition_system(&fb.subspace_basis, &x_r, &kc.curve.b, 2, &st)
                    .expect("system");
                let want = PackedQuad::from_polys(&sys.equations, sys.n_vars).unwrap();
                assert_eq!(t.instantiate(r), want, "n={n} r={r:#x}");
            }
        }
    }

    #[test]
    fn oracle_agrees_with_reference_on_every_target() {
        use crate::cryptanalysis::koblitz_index_calculus::{
            build_frobenius_factor_base, point_key, KoblitzCurve,
        };
        use crate::cryptanalysis::mq_fes::mq_fes_decompose_reference;
        for (a, n) in [(1u8, 7u32), (0, 13), (1, 17)] {
            let kc = KoblitzCurve::new(a, n).expect("curve");
            let fb = build_frobenius_factor_base(&kc, 0).expect("factor base");
            let st = FieldStructure::new(kc.n, &kc.curve.irreducible);
            let index_of: HashMap<_, _> = fb
                .points
                .iter()
                .enumerate()
                .map(|(i, p)| (point_key(p), i))
                .collect();
            let mut decomposed = 0;
            for t in 1..=160u64 {
                let target = kc.mul(kc.generator(), &BigUint::from(t * 7 + 3));
                let (got, _) = mq_fes_decompose(&kc, &fb, &index_of, &st, &target, 2);
                let (want, _) = mq_fes_decompose_reference(&kc, &fb, &index_of, &st, &target, 2);
                assert_eq!(got.is_some(), want.is_some(), "n={n} t={t}");
                if let Some(idxs) = got {
                    decomposed += 1;
                    let sum = idxs
                        .iter()
                        .fold(BinaryPoint::Infinity, |acc, &i| kc.add(&acc, &fb.points[i]));
                    assert_eq!(sum, target, "n={n} t={t}: relation does not verify");
                }
            }
            assert!(decomposed > 0, "n={n}: no target decomposed");
        }
    }
}

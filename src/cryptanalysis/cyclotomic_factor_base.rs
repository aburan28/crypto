//! # Cyclotomic-sparse Frobenius-stable factor bases (lead A2).
//!
//! When `p = 2n + 1` is prime and `ord_p(2) = n` (so `n ∈ {11, 23, 83,
//! 131}`), the field `F_{2^n}` is `F_2(γ)` for a primitive `p`-th root of
//! unity `γ`, and `μ_p ⊂ F_{2^n}^×`.  The **γ-weight-`w`** set
//!
//! ```text
//!     C_w = { Σ_{k ∈ S} γ^k : S ⊂ Z/p, |S| = w }
//! ```
//!
//! is closed under the field Frobenius (`x ↦ x²` acts as `S ↦ 2S`), has at
//! most `C(p, w)` elements, and — unlike the Hamming-weight set in the
//! type-II optimal normal basis `β_i = γ^{2^i} + γ^{−2^i}`, whose
//! membership is a cardinality constraint — admits a **polynomial**
//! membership description: `x ∈ C_1 ⟺ x^p = 1`, which in Boolean
//! coordinates is the product of the `popcount(p)` Frobenius images of
//! `x` (`x·x²·x⁴·x¹⁶ = 1` at `p = 23`), and `x ∈ C_w` is `x = Σ z_j` with
//! `z_j^p = 1`.
//!
//! This module builds both families as [`FrobeniusFactorBase`]s, reports
//! their census (abscissae, points, subgroup-usable points, signed orbits,
//! projected columns), and writes the Weil-descended decomposition
//! systems with either membership encoding so the falsifier in
//! `docs/ic/RESEARCH_IC_NOVEL_DIRECTIONS_20261007.md` §A2 — solving degree
//! and oracle wall of the cyclotomic system against a matched-size
//! Hamming-weight ONB system, plus natural 2-sum / 3-sum yields — can be
//! measured rather than argued.  `examples/a2_cyclotomic_base.rs` runs
//! it; `research/ic_techniques_20261012/a2_cyclotomic_and_delta/` holds
//! the evidence.
//!
//! Everything here is public synthetic / known-answer: the curves are the
//! Koblitz models `K_a` at toy degrees, targets are seeded multiples of
//! the public generator, and no logarithm is computed.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::time::Instant;

use num_bigint::BigUint;
use num_traits::ToPrimitive;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::binary_ecc::{BinaryPoint, F2mElement, IrreduciblePoly};

use super::koblitz_groebner::{
    f4_word_ops_thread, solve_boolean_system, solving_degree, sym_semaev_s3, system_degree,
    FieldStructure, SolveOptions, SolveStats, SolverEngine, SplitRule, SymElement, MAX_VARS,
};
use super::koblitz_index_calculus::{
    build_explicit_frobenius_orbit_factor_base, cofactor_project_factor_base, lift_candidate,
    order_of_2_mod_n, point_key, restrict_factor_base_to_orbits, FrobeniusFactorBase, KoblitzCurve,
};
use super::pq_groebner_f2::{F2BoolMono, F2BoolPoly};

// ── The cyclotomic structure of F_{2^n} ─────────────────────────────

/// Why a degree does not carry the construction.
fn check_degree(n: u32) -> Result<u32, String> {
    if !(3..=62).contains(&n) {
        return Err(format!("degree {n} outside the single-word range 3..=62"));
    }
    let p = 2 * n + 1;
    if !is_prime_u32(p) {
        return Err(format!("2n+1 = {p} is not prime at n = {n}"));
    }
    match order_of_2_mod_n(p) {
        Some(order) if order == n => Ok(p),
        Some(order) => Err(format!(
            "ord_{p}(2) = {order} ≠ {n}: γ ∉ F_(2^{n}), so the cyclotomic set is not inside the field"
        )),
        None => Err(format!("could not compute the order of 2 mod {p}")),
    }
}

fn is_prime_u32(v: u32) -> bool {
    if v < 2 {
        return false;
    }
    let mut d = 2u32;
    while d * d <= v {
        if v.is_multiple_of(d) {
            return false;
        }
        d += 1;
    }
    true
}

/// `C(n, k)` as `u64`, saturating.
pub fn binomial(n: u64, k: u64) -> u64 {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut acc: u128 = 1;
    for i in 0..k {
        acc = acc * (n - i) as u128 / (i + 1) as u128;
        if acc > u64::MAX as u128 {
            return u64::MAX;
        }
    }
    acc as u64
}

/// `F_{2^n} = F_2(γ)` with `γ` a primitive `(2n+1)`-th root of unity, in
/// the polynomial basis of `irr`, together with the type-II optimal
/// normal basis `β_i = γ^{2^i} + γ^{−2^i}` that the Hamming-weight control
/// base is written in.
#[derive(Clone, Debug)]
pub struct CyclotomicField {
    pub n: u32,
    /// `p = 2n + 1`.
    pub p: u32,
    pub irr: IrreduciblePoly,
    /// `γ`, the least (as an integer word) primitive `p`-th root of unity.
    pub gamma: F2mElement,
    /// `powers[k] = γ^k` for `k ∈ 0..p`.
    pub powers: Vec<F2mElement>,
    /// `onb[i] = γ^{2^i} + γ^{−2^i}`, `i ∈ 0..n` — a normal basis, so the
    /// Frobenius is the cyclic shift `onb[i] ↦ onb[i+1]`.
    pub onb: Vec<F2mElement>,
    /// Rows of the inverse of the ONB matrix: `onb_coordinates` reads the
    /// coordinate `i` of `x` as the parity of `x & onb_inverse_rows[i]`.
    onb_inverse_rows: Vec<u64>,
}

fn word(x: &F2mElement) -> u64 {
    x.raw_bits().first().copied().unwrap_or(0)
}

fn pow(x: &F2mElement, mut e: u64, irr: &IrreduciblePoly) -> F2mElement {
    let n = x.m_value();
    let mut base = x.clone();
    let mut acc = F2mElement::one(n);
    while e > 0 {
        if e & 1 == 1 {
            acc = acc.mul(&base, irr);
        }
        base = base.square(irr);
        e >>= 1;
    }
    acc
}

/// Rows of the inverse of the binary matrix whose columns are `columns`.
fn invert_f2(columns: &[u64], n: u32) -> Option<Vec<u64>> {
    let mut matrix: Vec<u64> = (0..n)
        .map(|row| {
            (0..n).fold(0u64, |bits, col| {
                bits | (((columns[col as usize] >> row) & 1) << col)
            })
        })
        .collect();
    let mut inverse: Vec<u64> = (0..n).map(|row| 1u64 << row).collect();
    for col in 0..n as usize {
        let pivot = (col..n as usize).find(|&row| (matrix[row] >> col) & 1 == 1)?;
        matrix.swap(col, pivot);
        inverse.swap(col, pivot);
        for row in 0..n as usize {
            if row != col && (matrix[row] >> col) & 1 == 1 {
                matrix[row] ^= matrix[col];
                inverse[row] ^= inverse[col];
            }
        }
    }
    Some(inverse)
}

impl CyclotomicField {
    /// Locate `γ` in the field `F_2[z]/(irr)` of degree `n`.
    pub fn new(irr: &IrreduciblePoly, n: u32) -> Result<Self, String> {
        let p = check_degree(n)?;
        if irr.degree != n {
            return Err("irreducible polynomial degree does not match n".into());
        }
        let group_order = (1u64 << n) - 1;
        if !group_order.is_multiple_of(p as u64) {
            return Err(format!("{p} does not divide 2^{n} − 1"));
        }
        let exponent = group_order / p as u64;
        // Every γ^k (k ≠ 0) is primitive since p is prime; take the least
        // word so the construction is canonical for the field.
        let mut gamma = None;
        for g in 2..(1u64 << n) {
            let candidate = pow(
                &F2mElement::from_biguint(&BigUint::from(g), n),
                exponent,
                irr,
            );
            if candidate != F2mElement::one(n) && !candidate.is_zero() {
                let mut least = candidate.clone();
                let mut current = candidate.clone();
                for _ in 1..p {
                    current = current.mul(&candidate, irr);
                    if word(&current) < word(&least) && current != F2mElement::one(n) {
                        least = current.clone();
                    }
                }
                gamma = Some(least);
                break;
            }
        }
        let gamma = gamma.ok_or("no primitive p-th root of unity found")?;
        let mut powers = Vec::with_capacity(p as usize);
        let mut current = F2mElement::one(n);
        for _ in 0..p {
            powers.push(current.clone());
            current = current.mul(&gamma, irr);
        }
        if current != F2mElement::one(n) {
            return Err("γ^p ≠ 1".into());
        }
        let onb: Vec<F2mElement> = (0..n)
            .map(|i| {
                let e = pow2_mod(i, p) as usize;
                powers[e].add(&powers[(p as usize - e) % p as usize])
            })
            .collect();
        let columns: Vec<u64> = onb.iter().map(word).collect();
        let onb_inverse_rows =
            invert_f2(&columns, n).ok_or("β_i = γ^(2^i) + γ^(−2^i) is not a basis")?;
        Ok(Self {
            n,
            p,
            irr: irr.clone(),
            gamma,
            powers,
            onb,
            onb_inverse_rows,
        })
    }

    /// The structure for a Koblitz curve's field.
    pub fn for_curve(kc: &KoblitzCurve) -> Result<Self, String> {
        if kc.k != 1 {
            return Err("cyclotomic bases are defined for the Koblitz models K_a (k = 1)".into());
        }
        Self::new(&kc.curve.irreducible, kc.n)
    }

    /// Hamming weight of `x` in the ONB coordinates.
    pub fn onb_coordinates(&self, x: &F2mElement) -> u64 {
        let w = word(x);
        self.onb_inverse_rows
            .iter()
            .enumerate()
            .fold(0u64, |acc, (i, row)| {
                acc | (((w & row).count_ones() as u64 & 1) << i)
            })
    }

    /// Every element `Σ_{k ∈ S} γ^k`, `|S| = w`, deduplicated and sorted.
    /// Capped at `C(p, w) ≤ 4·10⁶` subsets.
    pub fn gamma_weight_set(&self, w: usize) -> Result<Vec<F2mElement>, String> {
        self.weight_set(&self.powers, w)
    }

    /// Every element of ONB Hamming weight exactly `w`, sorted.
    pub fn onb_weight_set(&self, w: usize) -> Result<Vec<F2mElement>, String> {
        self.weight_set(&self.onb, w)
    }

    fn weight_set(&self, atoms: &[F2mElement], w: usize) -> Result<Vec<F2mElement>, String> {
        let count = binomial(atoms.len() as u64, w as u64);
        if w == 0 || w > atoms.len() {
            return Err(format!("weight {w} outside 1..={}", atoms.len()));
        }
        if count > 4_000_000 {
            return Err(format!(
                "C({}, {w}) = {count} subsets is above the 4·10⁶ cap",
                atoms.len()
            ));
        }
        let mut out: BTreeMap<u64, F2mElement> = BTreeMap::new();
        for_each_subset(atoms.len(), w, |subset| {
            let mut acc = F2mElement::zero(self.n);
            for &k in subset {
                acc = acc.add(&atoms[k]);
            }
            out.entry(word(&acc)).or_insert(acc);
        });
        Ok(out.into_values().collect())
    }

    /// Whether the set is closed under `x ↦ x²`.
    pub fn is_frobenius_stable(&self, set: &[F2mElement]) -> bool {
        let members: HashSet<u64> = set.iter().map(word).collect();
        set.iter()
            .all(|x| members.contains(&word(&x.square(&self.irr))))
    }

    /// The `n` Boolean coordinate equations of `z^p + 1 = 0` for a symbolic
    /// `z`: the product of the Frobenius images `z^{2^j}` over the set bits
    /// `j` of `p`, so the Boolean degree is `popcount(p)` times the degree
    /// of `z`'s coordinates.
    pub fn unity_membership(&self, z: &SymElement, st: &FieldStructure) -> Vec<F2BoolPoly> {
        let n_vars = z.coords[0].n_vars;
        let mut images = Vec::new();
        let mut current = z.clone();
        for j in 0..32 {
            if (self.p >> j) == 0 {
                break;
            }
            if (self.p >> j) & 1 == 1 {
                images.push(current.clone());
            }
            current = current.square(st);
        }
        let mut acc = images[0].clone();
        for image in &images[1..] {
            acc = acc.mul(image, st);
        }
        acc = acc.add(&SymElement::constant(
            &F2mElement::one(self.n),
            self.n,
            n_vars,
        ));
        acc.coords
    }
}

fn pow2_mod(i: u32, p: u32) -> u32 {
    let mut acc = 1u64;
    for _ in 0..i {
        acc = (acc * 2) % p as u64;
    }
    acc as u32
}

/// Visit every `k`-subset of `0..n` in lexicographic order.
pub fn for_each_subset(n: usize, k: usize, mut visit: impl FnMut(&[usize])) {
    if k == 0 || k > n {
        return;
    }
    let mut idx: Vec<usize> = (0..k).collect();
    loop {
        visit(&idx);
        // Advance.
        let mut i = k;
        loop {
            if i == 0 {
                return;
            }
            i -= 1;
            if idx[i] != i + n - k {
                break;
            }
            if i == 0 {
                return;
            }
        }
        idx[i] += 1;
        for j in i + 1..k {
            idx[j] = idx[j - 1] + 1;
        }
    }
}

// ── Factor bases and their census ───────────────────────────────────

/// Which family a base belongs to.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "family", rename_all = "snake_case")]
pub enum BaseFamily {
    /// `x = Σ_{k∈S} γ^k`, `|S| = w`.
    Cyclotomic { w: usize },
    /// Hamming weight exactly `w` in the type-II ONB.
    OnbWeight { w: usize },
}

/// What a base is, counted.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BaseCensus {
    pub family: BaseFamily,
    pub n: u32,
    pub a: u8,
    pub p: u32,
    /// `C(p, w)` (cyclotomic) or `C(n, w)` (ONB): the size before collisions.
    pub predicted_abscissae: u64,
    /// Distinct abscissae in the set.
    pub abscissae: usize,
    /// Closed under `x ↦ x²`.
    pub frobenius_stable: bool,
    /// Abscissae carrying a rational point.
    pub abscissae_with_points: usize,
    /// `|F|`: rational points above the set.
    pub points: usize,
    /// Points with `[r]P = O` (already in the prime-order subgroup).
    pub points_in_subgroup: usize,
    /// Frobenius orbits and signed Frobenius orbits (relation columns).
    pub orbits: usize,
    pub signed_orbits: usize,
    /// Distinct non-identity cofactor images `[h]P` and their signed
    /// Frobenius orbits — the subgroup-usable count after projection.
    pub projected_points: usize,
    pub projected_columns: usize,
    /// Orbits retained when the base was truncated to match another.
    pub retained_orbits: Option<usize>,
    pub build_ms: f64,
}

/// A built base with its census.
#[derive(Clone, Debug)]
pub struct CensusedBase {
    pub base: FrobeniusFactorBase,
    pub census: BaseCensus,
}

fn census(
    kc: &KoblitzCurve,
    cf: &CyclotomicField,
    family: BaseFamily,
    predicted: u64,
    set: &[F2mElement],
    base: &FrobeniusFactorBase,
    retained: Option<usize>,
    started: Instant,
) -> BaseCensus {
    let xs_with_points: HashSet<BigUint> = base
        .points
        .iter()
        .filter_map(|p| match p {
            BinaryPoint::Affine { x, .. } => Some(x.to_biguint()),
            BinaryPoint::Infinity => None,
        })
        .collect();
    let r = &kc.subgroup_order;
    let points_in_subgroup = base
        .points
        .par_iter()
        .filter(|p| kc.mul(p, r) == BinaryPoint::Infinity)
        .count();
    let (projected_points, projected_columns) = match cofactor_project_factor_base(kc, base) {
        Ok(projected) => (projected.points.len(), projected.signed_orbits.len()),
        Err(_) => (0, 0),
    };
    BaseCensus {
        family,
        n: kc.n,
        a: kc.a,
        p: cf.p,
        predicted_abscissae: predicted,
        abscissae: set.len(),
        frobenius_stable: cf.is_frobenius_stable(set),
        abscissae_with_points: xs_with_points.len(),
        points: base.points.len(),
        points_in_subgroup,
        orbits: base.orbits.len(),
        signed_orbits: base.signed_orbits.len(),
        projected_points,
        projected_columns,
        retained_orbits: retained,
        build_ms: started.elapsed().as_secs_f64() * 1e3,
    }
}

fn build_from_set(
    kc: &KoblitzCurve,
    cf: &CyclotomicField,
    family: BaseFamily,
    predicted: u64,
    set: Vec<F2mElement>,
    started: Instant,
) -> Result<CensusedBase, String> {
    if !cf.is_frobenius_stable(&set) {
        return Err("abscissa set is not Frobenius-stable".into());
    }
    let base = build_explicit_frobenius_orbit_factor_base(kc, &set)
        .ok_or("explicit-orbit constructor rejected the cyclotomic set")?;
    // The constructor closes the set under Frobenius; a stable set gains
    // nothing, and the census checks that too.
    if base.subspace.len() != set.len() {
        return Err(format!(
            "Frobenius closure grew the set from {} to {} abscissae",
            set.len(),
            base.subspace.len()
        ));
    }
    let c = census(kc, cf, family, predicted, &set, &base, None, started);
    Ok(CensusedBase { base, census: c })
}

/// The γ-weight-`w` base on `kc`.
pub fn build_cyclotomic_factor_base(
    kc: &KoblitzCurve,
    cf: &CyclotomicField,
    w: usize,
) -> Result<CensusedBase, String> {
    let started = Instant::now();
    let set = cf.gamma_weight_set(w)?;
    build_from_set(
        kc,
        cf,
        BaseFamily::Cyclotomic { w },
        binomial(cf.p as u64, w as u64),
        set,
        started,
    )
}

/// The ONB-Hamming-weight-`w` base on `kc`.
pub fn build_onb_weight_factor_base(
    kc: &KoblitzCurve,
    cf: &CyclotomicField,
    w: usize,
) -> Result<CensusedBase, String> {
    let started = Instant::now();
    let set = cf.onb_weight_set(w)?;
    build_from_set(
        kc,
        cf,
        BaseFamily::OnbWeight { w },
        binomial(cf.n as u64, w as u64),
        set,
        started,
    )
}

/// An ONB-weight base with (close to) `columns` signed Frobenius orbits:
/// the least weight `w` whose base has at least that many columns,
/// truncated to the first `columns` orbits by canonical abscissa.  The
/// control the §A2 falsifier compares against.
pub fn build_matched_onb_weight_factor_base(
    kc: &KoblitzCurve,
    cf: &CyclotomicField,
    columns: usize,
) -> Result<CensusedBase, String> {
    if columns == 0 {
        return Err("matched base needs at least one column".into());
    }
    let started = Instant::now();
    for w in 1..=(cf.n as usize / 2) {
        let built = build_onb_weight_factor_base(kc, cf, w)?;
        if built.base.signed_orbits.len() < columns {
            continue;
        }
        if built.base.signed_orbits.len() == columns {
            return Ok(built);
        }
        // Orbits ordered by their least abscissa word, deterministic.
        let mut keyed: Vec<(u64, usize)> = built
            .base
            .signed_orbits
            .iter()
            .enumerate()
            .map(|(o, members)| {
                let least = members
                    .iter()
                    .filter_map(|&i| match &built.base.points[i] {
                        BinaryPoint::Affine { x, .. } => Some(word(x)),
                        BinaryPoint::Infinity => None,
                    })
                    .min()
                    .unwrap_or(u64::MAX);
                (least, o)
            })
            .collect();
        keyed.sort_unstable();
        let keep: Vec<usize> = keyed.iter().take(columns).map(|&(_, o)| o).collect();
        let restricted = restrict_factor_base_to_orbits(kc, &built.base, &keep)
            .ok_or("orbit restriction produced no base")?;
        let set: Vec<F2mElement> = restricted.subspace.clone();
        let c = census(
            kc,
            cf,
            BaseFamily::OnbWeight { w },
            binomial(cf.n as u64, w as u64),
            &set,
            &restricted,
            Some(columns),
            started,
        );
        return Ok(CensusedBase {
            base: restricted,
            census: c,
        });
    }
    Err(format!(
        "no ONB weight up to n/2 reaches {columns} signed orbits"
    ))
}

// ── Weil-descended decomposition systems with membership encodings ──

/// How the summand abscissae are constrained to the base.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "encoding", rename_all = "snake_case")]
pub enum MembershipEncoding {
    /// `x_i = Σ_{j<w} z_{i,j}`, each `z_{i,j}` a full field unknown with
    /// `z^p = 1` (degree `popcount(p)`); `w·n` variables per summand.
    CyclotomicUnity { w: usize },
    /// `x_i = Σ_t v_{i,t} β_t` in the ONB, with every `(w+1)`-product of
    /// the `v_{i,t}` zero and `e_w(v_i) = 1` (exact weight); degree `w+1`;
    /// `n` variables per summand.
    OnbWeight { w: usize },
}

/// A decomposition system with a membership encoding.
#[derive(Clone, Debug)]
pub struct MembershipSystem {
    pub encoding: MembershipEncoding,
    pub equations: Vec<F2BoolPoly>,
    pub n_vars: usize,
    pub m: usize,
    /// Boolean unknowns per summand (`w·n` or `n`).
    pub vars_per_summand: usize,
    /// `x_i = Σ_t bit_{i,t} · summand_basis[t]`.
    pub summand_basis: Vec<F2mElement>,
    pub semaev_equations: usize,
    pub membership_equations: usize,
    /// Total Boolean degree of the system.
    pub degree: u32,
}

impl MembershipSystem {
    /// Abscissa of summand `i` at Boolean assignment `root`.
    pub fn summand_x(&self, root: u64, i: usize, n: u32) -> F2mElement {
        let mut acc = F2mElement::zero(n);
        for t in 0..self.vars_per_summand {
            if (root >> (i * self.vars_per_summand + t)) & 1 == 1 {
                acc = acc.add(&self.summand_basis[t]);
            }
        }
        acc
    }

    /// One-based ids of every variable, summands first and then the
    /// `(m − 2)·n` intermediate-point unknowns, for a WDSat branch list.
    /// WDSat assigns only the variables it is told to branch on (plus
    /// what XOR-Gauss propagation forces), and on the chained systems the
    /// intermediate unknowns are not forced, so a summand-only list
    /// leaves them unassigned and the model unparsable.
    pub fn branch_vars(&self) -> Vec<u32> {
        (1..=self.n_vars as u32).collect()
    }
}

/// Build the system for `R = P_1 + … + P_m` whose summands are constrained
/// by `encoding`.  `None` when the layout exceeds [`MAX_VARS`].
pub fn build_membership_system(
    cf: &CyclotomicField,
    st: &FieldStructure,
    b: &F2mElement,
    x_r: &F2mElement,
    encoding: MembershipEncoding,
    m: usize,
) -> Option<MembershipSystem> {
    if m < 2 {
        return None;
    }
    let n = cf.n;
    let nn = n as usize;
    let vars_per_summand = match encoding {
        MembershipEncoding::CyclotomicUnity { w } => w.checked_mul(nn)?,
        MembershipEncoding::OnbWeight { .. } => nn,
    };
    let n_vars = m
        .checked_mul(vars_per_summand)?
        .checked_add((m - 2).checked_mul(nn)?)?;
    if n_vars > MAX_VARS || vars_per_summand == 0 {
        return None;
    }
    let standard: Vec<F2mElement> = (0..n)
        .map(|i| F2mElement::from_bit_positions(&[i], n))
        .collect();
    let mut equations = Vec::new();
    let mut membership_equations = 0usize;
    let mut summands = Vec::with_capacity(m);
    let summand_basis: Vec<F2mElement>;
    match encoding {
        MembershipEncoding::CyclotomicUnity { w } => {
            summand_basis = (0..w).flat_map(|_| standard.iter().cloned()).collect();
            for i in 0..m {
                let mut x = SymElement::zero(n, n_vars);
                for j in 0..w {
                    let z = SymElement::from_subspace_vars(
                        &standard,
                        i * vars_per_summand + j * nn,
                        n,
                        n_vars,
                    );
                    let eqs = cf.unity_membership(&z, st);
                    membership_equations += eqs.len();
                    equations.extend(eqs);
                    x = x.add(&z);
                }
                summands.push(x);
            }
        }
        MembershipEncoding::OnbWeight { w } => {
            if w == 0 || w > nn {
                return None;
            }
            summand_basis = cf.onb.clone();
            for i in 0..m {
                let offset = i * nn;
                summands.push(SymElement::from_subspace_vars(&cf.onb, offset, n, n_vars));
                // Every (w+1)-product vanishes: HW ≤ w.
                if w < nn {
                    for_each_subset(nn, w + 1, |subset| {
                        let mask = subset
                            .iter()
                            .fold(0u64, |acc, &t| acc | (1u64 << (offset + t)));
                        equations.push(F2BoolPoly::from_monos(vec![F2BoolMono { mask }], n_vars));
                        membership_equations += 1;
                    });
                }
                // e_w(v) = 1: with HW ≤ w this is HW = w exactly.
                let mut monos = vec![F2BoolMono::one()];
                for_each_subset(nn, w, |subset| {
                    let mask = subset
                        .iter()
                        .fold(0u64, |acc, &t| acc | (1u64 << (offset + t)));
                    monos.push(F2BoolMono { mask });
                });
                equations.push(F2BoolPoly::from_monos(monos, n_vars));
                membership_equations += 1;
            }
        }
    }
    let target = SymElement::constant(x_r, n, n_vars);
    let before = equations.len();
    if m == 2 {
        equations.extend(sym_semaev_s3(&summands[0], &summands[1], &target, b, st));
    } else {
        let inter: Vec<SymElement> = (0..m - 2)
            .map(|i| SymElement::from_free_vars(m * vars_per_summand + i * nn, n, n_vars))
            .collect();
        equations.extend(sym_semaev_s3(&summands[0], &summands[1], &inter[0], b, st));
        for i in 0..m - 3 {
            equations.extend(sym_semaev_s3(
                &inter[i],
                &summands[i + 2],
                &inter[i + 1],
                b,
                st,
            ));
        }
        equations.extend(sym_semaev_s3(
            &inter[m - 3],
            &summands[m - 1],
            &target,
            b,
            st,
        ));
    }
    let semaev_equations = equations.len() - before;
    equations.retain(|e| !e.is_zero());
    let degree = system_degree(&equations);
    Some(MembershipSystem {
        encoding,
        equations,
        n_vars,
        m,
        vars_per_summand,
        summand_basis,
        semaev_equations,
        membership_equations,
        degree,
    })
}

// ── Measuring one system ───────────────────────────────────────────

/// The algebraic measurement of one target's system.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SystemMeasurement {
    pub target_index: usize,
    pub n_vars: usize,
    pub equations: usize,
    pub system_degree: u32,
    /// Lowest Macaulay degree at which the sparse profile resolves the
    /// system (refutes it or pins every variable), if reached by `d_max`.
    pub macaulay_solving_degree: Option<u32>,
    /// Highest Macaulay degree the profile actually built.
    pub macaulay_degree_reached: u32,
    pub macaulay_ms: f64,
    /// Matrix-F4 splitting solve.
    pub f4_ms: f64,
    pub f4_word_xors: u64,
    pub f4_max_degree_built: u32,
    pub f4_splits: usize,
    pub f4_reductions: usize,
    pub f4_exhausted: bool,
    pub boolean_roots: usize,
    /// Roots that lift to a genuine point decomposition over the base.
    pub lifted_decomposition: bool,
}

/// Macaulay solving-degree probe and matrix-F4 solve of one system,
/// lifting roots against `base`.
pub fn measure_system(
    kc: &KoblitzCurve,
    base: &FrobeniusFactorBase,
    index_of: &HashMap<(BigUint, BigUint), usize>,
    sys: &MembershipSystem,
    target: &BinaryPoint,
    target_index: usize,
    d_max: u32,
    engine: SolverEngine,
    node_budget: usize,
) -> SystemMeasurement {
    let started = Instant::now();
    let (solving, profiles) = if d_max >= sys.degree {
        solving_degree(&sys.equations, sys.n_vars, d_max)
    } else {
        (None, Vec::new())
    };
    let macaulay_ms = started.elapsed().as_secs_f64() * 1e3;
    let reached = profiles.iter().map(|p| p.degree).max().unwrap_or(0);

    let opts = SolveOptions {
        engine,
        max_solutions: 64,
        node_budget,
        split_rule: SplitRule::Auto,
    };
    let xors_before = f4_word_ops_thread();
    let started = Instant::now();
    let (roots, stats): (Vec<u64>, SolveStats) =
        solve_boolean_system(&sys.equations, sys.n_vars, &opts);
    let f4_ms = started.elapsed().as_secs_f64() * 1e3;
    let f4_word_xors = f4_word_ops_thread().saturating_sub(xors_before);
    let lifted = roots.iter().any(|&root| {
        let xs: Vec<F2mElement> = (0..sys.m).map(|i| sys.summand_x(root, i, kc.n)).collect();
        lift_candidate(kc, base, index_of, &xs, target).is_some()
    });
    SystemMeasurement {
        target_index,
        n_vars: sys.n_vars,
        equations: sys.equations.len(),
        system_degree: sys.degree,
        macaulay_solving_degree: solving,
        macaulay_degree_reached: reached,
        macaulay_ms,
        f4_ms,
        f4_word_xors,
        f4_max_degree_built: stats.max_degree_built,
        f4_splits: stats.splits,
        f4_reductions: stats.reductions,
        f4_exhausted: stats.exhausted,
        boolean_roots: roots.len(),
        lifted_decomposition: lifted,
    }
}

/// Seeded public targets `[k]G`, `k` drawn from a splitmix64 stream.
pub fn seeded_targets(kc: &KoblitzCurve, count: usize, seed: u64) -> Vec<(u64, BinaryPoint)> {
    let r = kc.subgroup_order.to_u64().unwrap_or(u64::MAX);
    let mut state = seed ^ 0xa2c1_c107_0d1c_0000;
    let mut out = Vec::with_capacity(count);
    let g = kc.generator();
    while out.len() < count {
        state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^= z >> 31;
        let k = 1 + z % (r - 1);
        out.push((k, kc.mul(g, &BigUint::from(k))));
    }
    out
}

/// Index map of a base, keyed as [`lift_candidate`] expects.
pub fn index_map(base: &FrobeniusFactorBase) -> HashMap<(BigUint, BigUint), usize> {
    base.points
        .iter()
        .enumerate()
        .map(|(i, p)| (point_key(p), i))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `K_1 / F_2^11` (`#E = 2·991`).  `K_0` at `n = 11` has `#E = 4·23²`,
    /// so the constructor rejects it: the `K_0` toy rung is `n = 23`.
    fn field_11() -> (KoblitzCurve, CyclotomicField) {
        let kc = KoblitzCurve::new(1, 11).expect("K_1 / F_2^11");
        let cf = CyclotomicField::for_curve(&kc).expect("23 is prime with ord_23(2) = 11");
        (kc, cf)
    }

    #[test]
    fn degrees_with_the_structure_include_the_ladder() {
        // The falsifier's toy rungs, and the exclusions that matter: a
        // prime 2n+1 with 2 a non-residue (ord = 2n, so γ ∉ F_{2^n}), a
        // composite 2n+1, and the landed n = 83 which needs two words.
        assert_eq!(check_degree(11), Ok(23));
        assert_eq!(check_degree(23), Ok(47));
        assert!(check_degree(5).unwrap_err().contains("ord_11(2) = 10"));
        assert!(check_degree(31).unwrap_err().contains("not prime"));
        assert!(check_degree(83).unwrap_err().contains("3..=62"));
        let have: Vec<u32> = (3..=62).filter(|&n| check_degree(n).is_ok()).collect();
        assert!(have.contains(&11) && have.contains(&23));
        for &n in &have {
            let p = 2 * n + 1;
            assert!(
                ((1u64 << n) - 1).is_multiple_of(p as u64),
                "μ_{p} ⊂ F_2^{n}"
            );
        }
    }

    #[test]
    fn gamma_weight_sets_are_frobenius_stable_with_predicted_size() {
        let (_, cf) = field_11();
        assert_eq!(cf.p, 23);
        assert_eq!(pow(&cf.gamma, 23, &cf.irr), F2mElement::one(11));
        assert_ne!(cf.gamma, F2mElement::one(11));
        // The Golay code [23,12,7] has minimum weight 7, so sums over
        // subsets of size ≤ 3 are pairwise distinct: |C_w| = C(23, w).
        for w in 1..=3 {
            let set = cf.gamma_weight_set(w).expect("set");
            assert_eq!(set.len() as u64, binomial(23, w as u64), "w = {w}");
            assert!(cf.is_frobenius_stable(&set), "w = {w}");
        }
        // ONB weight sets: a basis, so C(n, w) each, and Frobenius is a
        // coordinate shift.
        for w in 1..=3 {
            let set = cf.onb_weight_set(w).expect("set");
            assert_eq!(set.len() as u64, binomial(11, w as u64));
            assert!(cf.is_frobenius_stable(&set));
            for x in &set {
                assert_eq!(cf.onb_coordinates(x).count_ones() as usize, w);
                let shifted = cf.onb_coordinates(&x.square(&cf.irr));
                let coords = cf.onb_coordinates(x);
                assert_eq!(shifted, ((coords << 1) | (coords >> 10)) & ((1 << 11) - 1));
            }
        }
        // The two families are not nested the obvious way round: every
        // ONB basis element β_i = γ^{2^i} + γ^{−2^i} is a γ-weight-2
        // element, while no γ^k is an ONB basis element (and
        // γ^0 = 1 = Σ β_i has full ONB weight).
        let one = cf.gamma_weight_set(1).expect("μ_23");
        let two: HashSet<u64> = cf
            .gamma_weight_set(2)
            .expect("C_2")
            .iter()
            .map(word)
            .collect();
        assert_eq!(cf.onb_coordinates(&F2mElement::one(11)).count_ones(), 11);
        assert_eq!(
            one.iter()
                .filter(|x| cf.onb_coordinates(x).count_ones() == 1)
                .count(),
            0
        );
        for beta in &cf.onb {
            assert!(two.contains(&word(beta)), "β_i ∈ C_2");
        }
    }

    #[test]
    fn lifted_points_lie_on_the_curve_and_project_into_the_subgroup() {
        let (kc, cf) = field_11();
        for w in 1..=2 {
            let built = build_cyclotomic_factor_base(&kc, &cf, w).expect("base");
            let c = &built.census;
            assert!(c.frobenius_stable);
            assert_eq!(c.abscissae as u64, c.predicted_abscissae);
            for p in &built.base.points {
                assert!(kc.curve.is_on_curve(p));
                let projected = kc.mul(p, &kc.cofactor);
                assert!(
                    projected == BinaryPoint::Infinity
                        || kc.mul(&projected, &kc.subgroup_order) == BinaryPoint::Infinity
                );
            }
            assert_eq!(c.points, built.base.points.len());
            assert!(c.projected_points <= c.points);
            // Every signed orbit is closed: all its points share the orbit id.
            for (o, members) in built.base.signed_orbits.iter().enumerate() {
                for &i in members {
                    assert_eq!(built.base.signed_orbit_of[i].0, o);
                }
            }
        }
        // On K_1 every γ^k ≠ 1 lifts (Tr(γ^k) + Tr(1) + Tr(γ^{-k}) = 0 for
        // k ≠ 0 since Tr(γ^k) + Tr(γ^{-k}) = 1 when −1 is a non-residue)
        // and x = 1 does not: 22 abscissae, 44 points, the two Frobenius
        // orbits of μ_23 \ {1} as the two signed orbits.
        let w1 = build_cyclotomic_factor_base(&kc, &cf, 1).expect("w = 1");
        assert_eq!(w1.census.abscissae, 23);
        assert_eq!(w1.census.abscissae_with_points, 22);
        assert_eq!(w1.census.points, 44);
        assert_eq!(w1.census.signed_orbits, 2);
    }

    /// On `K_0` the same trace identity runs the other way: no `γ^k ≠ 1`
    /// lifts, the weight-1 base is the two 4-torsion points above `x = 1`
    /// and nothing is subgroup-usable.  `n = 23` is the `K_0` toy rung.
    #[test]
    fn k0_weight_one_base_is_the_four_torsion_above_x_equals_one() {
        let kc = KoblitzCurve::new(0, 23).expect("K_0 / F_2^23");
        let cf = CyclotomicField::for_curve(&kc).expect("47 is prime with ord_47(2) = 23");
        assert_eq!(cf.p, 47);
        let w1 = build_cyclotomic_factor_base(&kc, &cf, 1).expect("w = 1");
        assert_eq!(w1.census.abscissae, 47);
        assert_eq!(w1.census.abscissae_with_points, 1);
        assert_eq!(w1.census.points, 2);
        assert_eq!(w1.census.points_in_subgroup, 0);
        assert_eq!(w1.census.projected_points, 0);
        for p in &w1.base.points {
            assert_eq!(kc.mul(p, &BigUint::from(4u32)), BinaryPoint::Infinity);
        }
    }

    #[test]
    fn matched_onb_base_has_the_requested_columns() {
        let (kc, cf) = field_11();
        let cyc = build_cyclotomic_factor_base(&kc, &cf, 2).expect("w = 2");
        let columns = cyc.base.signed_orbits.len();
        assert!(columns > 0);
        let onb = build_matched_onb_weight_factor_base(&kc, &cf, columns).expect("matched");
        assert_eq!(onb.base.signed_orbits.len(), columns);
        assert!(onb.census.frobenius_stable);
        let BaseFamily::OnbWeight { w } = onb.census.family else {
            panic!("matched base is an ONB base");
        };
        for x in &onb.base.subspace {
            assert_eq!(cf.onb_coordinates(x).count_ones() as usize, w);
        }
    }

    #[test]
    fn unity_membership_is_satisfied_exactly_by_mu_p_at_n_11() {
        let (kc, cf) = field_11();
        let st = FieldStructure::new(11, &kc.curve.irreducible);
        let z = SymElement::from_subspace_vars(
            &(0..11)
                .map(|i| F2mElement::from_bit_positions(&[i], 11))
                .collect::<Vec<_>>(),
            0,
            11,
            11,
        );
        let eqs = cf.unity_membership(&z, &st);
        assert_eq!(eqs.len(), 11);
        assert_eq!(system_degree(&eqs), 4, "23 = 16 + 4 + 2 + 1");
        let members: HashSet<u64> = cf
            .gamma_weight_set(1)
            .expect("μ_23")
            .iter()
            .map(word)
            .collect();
        assert_eq!(members.len(), 23);
        let mut satisfied = 0;
        for x in 0..(1u64 << 11) {
            let holds = eqs.iter().all(|e| e.eval(x) == 0);
            assert_eq!(holds, members.contains(&x), "x = {x}");
            satisfied += holds as usize;
        }
        assert_eq!(satisfied, 23);
    }

    #[test]
    fn membership_systems_find_planted_decompositions() {
        let k1 = KoblitzCurve::new(1, 11).expect("K_1 / F_2^11");
        let cf = CyclotomicField::for_curve(&k1).expect("field");
        let st = FieldStructure::new(11, &k1.curve.irreducible);
        for (encoding, w) in [
            (MembershipEncoding::CyclotomicUnity { w: 1 }, 1usize),
            (MembershipEncoding::OnbWeight { w: 2 }, 2),
        ] {
            let built = match encoding {
                MembershipEncoding::CyclotomicUnity { .. } => {
                    build_cyclotomic_factor_base(&k1, &cf, w).expect("base")
                }
                MembershipEncoding::OnbWeight { .. } => {
                    build_onb_weight_factor_base(&k1, &cf, w).expect("base")
                }
            };
            let index_of = index_map(&built.base);
            let pts = &built.base.points;
            assert!(pts.len() >= 4);
            let target = k1.add(&pts[0], &pts[3]);
            let BinaryPoint::Affine { x: x_r, .. } = &target else {
                panic!("planted target is infinity");
            };
            let sys = build_membership_system(&cf, &st, &k1.curve.b, x_r, encoding, 2)
                .expect("system fits");
            assert_eq!(sys.n_vars, 22);
            let measured = measure_system(
                &k1,
                &built.base,
                &index_of,
                &sys,
                &target,
                0,
                2,
                SolverEngine::MatrixF4 { max_degree: 3 },
                4096,
            );
            assert!(measured.boolean_roots > 0, "{encoding:?}");
            assert!(measured.lifted_decomposition, "{encoding:?}: {measured:?}");
            // Every root respects the membership encoding.
            let opts = SolveOptions {
                engine: SolverEngine::MatrixF4 { max_degree: 3 },
                max_solutions: 64,
                node_budget: 4096,
                split_rule: SplitRule::Auto,
            };
            let (roots, _) = solve_boolean_system(&sys.equations, sys.n_vars, &opts);
            let allowed: HashSet<u64> = built.base.subspace.iter().map(word).collect();
            for root in roots {
                for i in 0..2 {
                    let x = sys.summand_x(root, i, 11);
                    assert!(
                        allowed.contains(&word(&x)),
                        "{encoding:?} root leaves the set"
                    );
                }
            }
        }
    }

    #[test]
    fn subset_enumeration_counts_binomials() {
        for (n, k) in [(5, 2), (7, 3), (11, 1), (6, 6)] {
            let mut count = 0u64;
            let mut last: Option<Vec<usize>> = None;
            for_each_subset(n, k, |s| {
                if let Some(prev) = &last {
                    assert!(prev.as_slice() < s, "lexicographic order");
                }
                last = Some(s.to_vec());
                count += 1;
            });
            assert_eq!(count, binomial(n as u64, k as u64));
        }
    }
}

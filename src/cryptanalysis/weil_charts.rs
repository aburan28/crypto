//! Exact two-summand S3 solving on a Frobenius union of small F2 spaces.
//!
//! Preserve the complete union: every unordered component pair is visited.
//! Normalize (i,j) by inverse Frobenius^i, so only j-i product tensors are
//! stored. Optional row-space projection exposes all linear consequences of
//! the coordinate equations before calling the existing F4 splitting solver.
//! Product ranks and XOR counters are diagnostics, not total attack costs.
use super::{
    koblitz_fast::{FastCurve, FastPoint},
    koblitz_groebner::{solve_boolean_system_filtered, SolveOptions, SolveStats},
    koblitz_index_calculus::{
        enumerate_decompose, lift_candidate, FrobeniusFactorBase, KoblitzCurve,
    },
    pq_groebner_f2::{F2BoolMono, F2BoolPoly},
    semaev_decomp::Gf2,
};
use crate::binary_ecc::{BinaryPoint, F2mElement};
use num_bigint::BigUint;
use std::collections::{BTreeSet, HashMap, HashSet};

fn bits(x: &F2mElement) -> u64 {
    x.raw_bits().first().copied().unwrap_or(0)
}

/// Canonical F2 basis, with highest-bit pivots, independent of input ordering.
pub fn canonical_basis(values: impl IntoIterator<Item = u64>) -> Vec<u64> {
    let mut pivots = [0u64; 64];
    for mut v in values {
        while v != 0 {
            let p = 63 - v.leading_zeros() as usize;
            if pivots[p] != 0 {
                v ^= pivots[p];
            } else {
                pivots[p] = v;
                break;
            }
        }
    }
    for p in 0..64 {
        if pivots[p] != 0 {
            for q in p + 1..64 {
                if pivots[q] >> p & 1 != 0 {
                    pivots[q] ^= pivots[p];
                }
            }
        }
    }
    pivots.into_iter().filter(|&v| v != 0).collect()
}

fn frob(f: &Gf2, mut x: u64, k: u32) -> u64 {
    for _ in 0..k {
        x = f.sqr(x);
    }
    x
}

fn combination(basis: &[u64], assignment: u64) -> u64 {
    basis.iter().enumerate().fold(
        0,
        |a, (i, b)| if assignment >> i & 1 != 0 { a ^ b } else { a },
    )
}

#[derive(Clone, Debug)]
struct PairTemplate {
    squares: Vec<u64>,
    products: Vec<u64>,
    product_squares: Vec<u64>,
    product_rank: usize,
}

/// 63 field equations plus at most two sets of seven support constraints.
struct CoefficientRows {
    data: [u64; 77],
    len: usize,
}

impl CoefficientRows {
    fn new(len: usize) -> Self {
        assert!(len <= 77);
        Self { data: [0; 77], len }
    }
    fn extend(&mut self, rows: impl IntoIterator<Item = u64>) {
        for row in rows {
            self.data[self.len] = row;
            self.len += 1;
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct WeilSolveStats {
    pub component_pairs: usize,
    pub duplicate_pairs: usize,
    pub algebraic_pairs: usize,
    pub projection_refutations: usize,
    pub linear_constraints: usize,
    pub projection_word_xors: u64,
    pub coefficient_field_muls: u64,
    pub coefficient_field_squares: u64,
    /// Exact assignments checked in small projected affine spaces. Charged
    /// alongside reductions to the query budget, but reported separately.
    pub residual_assignments: usize,
    pub solver: SolveStats,
}

/// Immutable target-independent preprocessing. Only complete, validated covers
/// are accepted. At most seven Boolean coordinates per component are supported
/// by this experimental single-word coefficient-row implementation.
#[derive(Clone, Debug)]
pub struct WeilChartPlan {
    field: Gf2,
    k: u32,
    a: u64,
    b: u64,
    charts: Vec<Vec<u64>>,
    templates: Vec<PairTemplate>,
    domain: Vec<u64>,
    project_linear: bool,
    points: Vec<BinaryPoint>,
    point_lifts: Option<(FastCurve, HashMap<u64, Vec<(FastPoint, usize)>>)>,
    point_constraints: Vec<Vec<u64>>,
    point_chart: Option<Box<WeilChartPlan>>,
}

impl WeilChartPlan {
    pub fn new(
        kc: &KoblitzCurve,
        fb: &FrobeniusFactorBase,
        seed: &[F2mElement],
        project_linear: bool,
    ) -> Result<Self, String> {
        if kc.n > 63 || seed.is_empty() || seed.len() > 7 {
            return Err("chart solver requires n<=63 and 1..=7 seed coordinates".into());
        }
        if seed
            .iter()
            .any(|x| x.raw_bits().len() > 1 || bits(x) >> kc.n != 0)
        {
            return Err("seed exceeds field width".into());
        }
        let field = Gf2::new(&kc.curve.irreducible);
        let first: Vec<_> = seed.iter().map(bits).collect();
        let canonical = canonical_basis(first.iter().copied());
        if canonical.len() != seed.len() {
            return Err("dependent seed basis".into());
        }
        let mut charts = Vec::new();
        let mut current = first;
        for _ in 0..kc.extension_degree() {
            charts.push(current.clone());
            current = current.iter().map(|&x| frob(&field, x, kc.k)).collect();
            if canonical_basis(current.iter().copied()) == canonical {
                break;
            }
        }
        if canonical_basis(current.iter().copied()) != canonical {
            return Err("Frobenius cover did not close".into());
        }
        let mut domain = BTreeSet::new();
        for basis in &charts {
            for mask in 0..1u64 << basis.len() {
                domain.insert(combination(basis, mask));
            }
        }
        let expected: BTreeSet<_> = fb.subspace.iter().map(bits).collect();
        if domain != expected {
            return Err("chart cover differs from factor-base domain".into());
        }
        let templates = charts
            .iter()
            .map(|basis| {
                let products: Vec<_> = charts[0]
                    .iter()
                    .flat_map(|&u| basis.iter().map(move |&v| (u, v)))
                    .map(|(u, v)| field.mul(u, v))
                    .collect();
                PairTemplate {
                    squares: charts[0]
                        .iter()
                        .chain(basis)
                        .map(|&x| field.sqr(x))
                        .collect(),
                    product_squares: products.iter().map(|&x| field.sqr(x)).collect(),
                    product_rank: canonical_basis(products.iter().copied()).len(),
                    products,
                }
            })
            .collect();
        let point_lifts = if let Some(fc) = FastCurve::new(&kc.curve) {
            let mut lifts: HashMap<u64, Vec<(FastPoint, usize)>> = HashMap::new();
            for (i, p) in fb.points.iter().enumerate() {
                let p = fc.lift(p);
                if p.infinity || !domain.contains(&p.x) || !fc.is_on_curve(p) {
                    return Err("invalid factor-base point".into());
                }
                lifts.entry(p.x).or_insert_with(|| Vec::with_capacity(2)).push((p, i));
            }
            Some((fc, lifts))
        } else {
            None
        };
        // A valid point's local coordinates lie in the span of the rational
        // abscissae present in that component. Keep these necessary equations
        // separate: the unrestricted algebraic visitor must retain S3 roots
        // that do not lift to rational factor-base points.
        let rational_xs: HashSet<_> = fb.points.iter().filter_map(|p| match p {
            BinaryPoint::Affine { x, .. } => Some(bits(x)),
            BinaryPoint::Infinity => None,
        }).collect();
        let mut constraint_cache: HashMap<Vec<u64>, Vec<u64>> = HashMap::new();
        let point_constraints: Vec<Vec<u64>> = charts.iter().map(|basis| {
            let support = canonical_basis((0..1u64 << basis.len())
                .filter(|&m| rational_xs.contains(&combination(basis, m))));
            constraint_cache.entry(support.clone()).or_insert_with(|| {
                canonical_basis((1..1u64 << basis.len())
                    .filter(|&mask| support.iter().all(|s| (mask & s).count_ones() % 2 == 0)))
            }).clone()
        }).collect();
        let support_basis = canonical_basis((0..1u64 << charts[0].len())
            .filter(|&m| rational_xs.contains(&combination(&charts[0], m))));
        let point_chart = if !support_basis.is_empty() && support_basis.len() < seed.len() {
            let compact_seed: Vec<_> = support_basis.iter()
                .map(|&m| field.to_element(combination(&charts[0], m))).collect();
            let mut compact_domain = BTreeSet::new();
            for mask in 0..1u64 << compact_seed.len() {
                let mut x = compact_seed.iter().enumerate()
                    .fold(0, |x, (i, b)| if mask >> i & 1 != 0 { x ^ bits(b) } else { x });
                for _ in 0..kc.extension_degree() {
                    compact_domain.insert(x);
                    x = frob(&field, x, kc.k);
                }
            }
            // A caller can supply a non-Frobenius-closed subset of points.
            // Compress only when every original point remains covered.
            if rational_xs.iter().all(|x| compact_domain.contains(x)) {
                let mut compact_fb = fb.clone();
                compact_fb.subspace = compact_domain.into_iter().map(|x| field.to_element(x)).collect();
                Some(Box::new(Self::new(kc, &compact_fb, &compact_seed, project_linear)?))
            } else {
                None
            }
        } else {
            None
        };
        Ok(Self {
            field,
            k: kc.k,
            a: bits(&kc.curve.a),
            b: bits(&kc.curve.b),
            charts,
            templates,
            domain: domain.into_iter().collect(),
            project_linear,
            points: fb.points.clone(),
            point_lifts,
            point_constraints,
            point_chart,
        })
    }

    pub fn component_count(&self) -> usize {
        self.charts.len()
    }
    pub fn seed_dimension(&self) -> usize {
        self.charts[0].len()
    }
    pub fn relative_product_ranks(&self) -> Vec<usize> {
        self.templates.iter().map(|t| t.product_rank).collect()
    }
    pub fn pair_weighted_product_rank(&self) -> f64 {
        let c = self.component_count();
        let sum: usize = self
            .templates
            .iter()
            .enumerate()
            .map(|(d, t)| (c - d) * t.product_rank)
            .sum();
        sum as f64 / (c * (c + 1) / 2) as f64
    }
    /// Tensor payload only; excludes allocator metadata and field tables.
    pub fn tensor_bytes(&self) -> usize {
        self.templates
            .iter()
            .map(|t| 8 * (t.squares.len() + t.products.len() + t.product_squares.len()))
            .sum()
    }
    pub fn component_elements(&self) -> Vec<Vec<u64>> {
        self.charts
            .iter()
            .map(|b| (0..1u64 << b.len()).map(|m| combination(b, m)).collect())
            .collect()
    }
    /// Dimensions of the spans of local coordinates with factor-base lifts.
    pub fn rational_support_dimensions(&self) -> Vec<usize> {
        self.point_constraints.iter().map(|cs| self.seed_dimension() - cs.len()).collect()
    }
    pub fn matches(&self, kc: &KoblitzCurve, fb: &FrobeniusFactorBase) -> bool {
        self.points == fb.points
            && self.field.n == kc.n
            && self.k == kc.k
            && self.a == bits(&kc.curve.a)
            && self.b == bits(&kc.curve.b)
            && self.field.irr
                == kc
                    .curve
                    .irreducible
                    .low_terms
                    .iter()
                    .fold(1u64 << kc.n, |x, &b| x | (1 << b))
            && self.domain
                == fb
                    .subspace
                    .iter()
                    .map(bits)
                    .collect::<BTreeSet<_>>()
                    .into_iter()
                    .collect::<Vec<_>>()
    }

    fn rows(&self, delta: usize, r: u64, stats: &mut WeilSolveStats) -> CoefficientRows {
        let t = &self.templates[delta];
        let ell = self.seed_dimension();
        let quadratic = ell * ell;
        let r2 = self.field.sqr(r);
        stats.coefficient_field_squares += 1;
        let mut rows = CoefficientRows::new(self.field.n as usize);
        let mut put = |mut coefficient: u64, column: usize| {
            while coefficient != 0 {
                let bit = coefficient.trailing_zeros() as usize;
                rows.data[bit] ^= 1 << column;
                coefficient &= coefficient - 1;
            }
        };
        for (j, (&p, &p2)) in t.products.iter().zip(&t.product_squares).enumerate() {
            put(p2 ^ self.field.mul(r, p), j);
        }
        for (j, &x2) in t.squares.iter().enumerate() {
            put(self.field.mul(r2, x2), quadratic + j);
        }
        put(self.b, quadratic + 2 * ell);
        stats.coefficient_field_muls += (quadratic + 2 * ell) as u64;
        rows
    }

    /// Visit each unordered algebraic pair once. Returning true from `accept`
    /// stops on a witness. Exhaustion is explicit and budgets cover the whole
    /// cover, rather than resetting for every component combination.
    pub fn visit_pairs(
        &self,
        x_target: u64,
        options: &SolveOptions,
        accept: impl FnMut(u64, u64) -> bool,
    ) -> WeilSolveStats {
        self.visit_pairs_inner(x_target, options, false, accept)
    }

    fn visit_pairs_inner(
        &self,
        x_target: u64,
        options: &SolveOptions,
        point_search: bool,
        mut accept: impl FnMut(u64, u64) -> bool,
    ) -> WeilSolveStats {
        let mut stats = WeilSolveStats::default();
        if x_target >> self.field.n != 0 {
            stats.solver.exhausted = true;
            return stats;
        }
        let e = self.field.n / self.k;
        let mut rotations = vec![x_target];
        for i in 1..e as usize {
            rotations.push(frob(&self.field, rotations[i - 1], self.k));
        }
        let mut seen = HashSet::new();
        let ell = self.seed_dimension();
        for i in 0..self.charts.len() {
            for j in i..self.charts.len() {
                stats.component_pairs += 1;
                let r = rotations[(e as usize - i) % e as usize];
                let mut rows = self.rows(j - i, r, &mut stats);
                if point_search {
                    rows.extend(self.point_constraints[i].iter().map(|row| row << (ell * ell)));
                    rows.extend(self.point_constraints[j].iter().map(|row| row << (ell * ell + ell)));
                }
                let prepared = prepare(rows, ell, self.project_linear, &mut stats);
                let Some(prepared) = prepared else {
                    stats.projection_refutations += 1;
                    continue;
                };
                let remaining = options.node_budget.saturating_sub(
                    stats.solver.reductions.saturating_add(stats.residual_assignments));
                if remaining == 0 {
                    stats.solver.exhausted = true;
                    return stats;
                }
                let opts = SolveOptions {
                    node_budget: remaining,
                    max_solutions: usize::MAX,
                    ..*options
                };
                let mut stop = false;
                let mut accept_assignment = |assignment| {
                        let x = combination(&self.charts[i], assignment);
                        let y = combination(&self.charts[j], assignment >> ell);
                        let pair = if x <= y { (x, y) } else { (y, x) };
                        if !seen.insert(pair) {
                            stats.duplicate_pairs += 1;
                            return false;
                        }
                        stats.algebraic_pairs += 1;
                        stop = accept(pair.0, pair.1);
                        stop
                };
                if let Some(rows) = &prepared.residual_rows {
                    for root in 0..1u64 << prepared.n_vars {
                        if root as usize >= remaining {
                            stats.solver.exhausted = true;
                            break;
                        }
                        stats.residual_assignments += 1;
                        let assignment = prepared.lift(root);
                        let mut monomials = (assignment << (ell * ell))
                            | (1u64 << (ell * ell + 2 * ell));
                        let right = assignment >> ell;
                        for bit in 0..ell {
                            if assignment >> bit & 1 != 0 {
                                monomials |= right << (bit * ell);
                            }
                        }
                        if rows.data[..rows.len].iter().all(|r| (r & monomials).count_ones() % 2 == 0)
                            && accept_assignment(assignment)
                        {
                            break;
                        }
                    }
                } else {
                    let (_, got) = solve_boolean_system_filtered(
                        &prepared.equations,
                        prepared.n_vars,
                        &opts,
                        |root| accept_assignment(prepared.lift(root)),
                    );
                    stats.solver.reductions += got.reductions;
                    stats.solver.infeasible_branches += got.infeasible_branches;
                    stats.solver.propagations += got.propagations;
                    stats.solver.splits += got.splits;
                    stats.solver.max_degree_built =
                        stats.solver.max_degree_built.max(got.max_degree_built);
                    stats.solver.oversize += got.oversize;
                    stats.solver.exhausted |= got.exhausted;
                }
                if stop || stats.solver.exhausted {
                    return stats;
                }
            }
        }
        stats
    }

    /// Search for a verified point decomposition; None rejects mismatched
    /// field/curve/domain metadata, rather than silently using an incomplete cover.
    pub fn decompose(
        &self,
        kc: &KoblitzCurve,
        fb: &FrobeniusFactorBase,
        index: &HashMap<(BigUint, BigUint), usize>,
        target: &BinaryPoint,
        options: &SolveOptions,
    ) -> Option<(Option<Vec<usize>>, WeilSolveStats)> {
        if !self.matches(kc, fb) {
            return None;
        }
        Some(self.decompose_validated(kc, fb, index, target, options))
    }

    /// The immutable curve, domain and point order must already have passed
    /// `matches` at pipeline entry. This avoids rebuilding a set every trial.
    pub(crate) fn decompose_validated(
        &self,
        kc: &KoblitzCurve,
        fb: &FrobeniusFactorBase,
        index: &HashMap<(BigUint, BigUint), usize>,
        target: &BinaryPoint,
        options: &SolveOptions,
    ) -> (Option<Vec<usize>>, WeilSolveStats) {
        if let Some(compact) = &self.point_chart {
            // Point indices and coverage were checked when this immutable
            // plan was constructed, and again by the parent pipeline entry.
            return compact.decompose_validated(kc, fb, index, target, options);
        }
        let BinaryPoint::Affine { x, .. } = target else {
            // Infinity has no S3 abscissa. Handle the group identity directly.
            return (
                enumerate_decompose(kc, fb, index, target, 2),
                WeilSolveStats::default(),
            );
        };
        let mut answer = None;
        let stats = self.visit_pairs_inner(bits(x), options, true, |a, b| {
            if let Some((fc, lifts)) = &self.point_lifts {
                let Some(left) = lifts.get(&a) else { return false; };
                let Some(right) = lifts.get(&b) else { return false; };
                for &(p, i) in left {
                    for &(q, j) in right {
                        if fc.add(p, q) == fc.lift(target) {
                            answer = Some(vec![i.min(j), i.max(j)]);
                            return true;
                        }
                    }
                }
                return false;
            }
            let xs = [
                F2mElement::from_biguint(&BigUint::from(a), kc.n),
                F2mElement::from_biguint(&BigUint::from(b), kc.n),
            ];
            answer = lift_candidate(kc, fb, index, &xs, target);
            answer.is_some()
        });
        (answer, stats)
    }
}

struct Prepared {
    equations: Vec<F2BoolPoly>,
    n_vars: usize,
    expressions: [u64; 14],
    n_original_vars: usize,
    residual_rows: Option<CoefficientRows>,
}
impl Prepared {
    fn lift(&self, assignment: u64) -> u64 {
        let with_constant = assignment | (1u64 << self.n_vars);
        self.expressions[..self.n_original_vars].iter().enumerate().fold(0, |a, (i, &e)| {
            a | (((e & with_constant).count_ones() as u64 & 1) << i)
        })
    }
}

fn rref(rows: &mut [u64], columns: usize, xors: &mut u64) -> u64 {
    let mut pivots = 0;
    let mut rank = 0;
    for c in 0..columns {
        let p = rank;
        let Some(found) = (p..rows.len()).find(|&i| rows[i] >> c & 1 != 0) else {
            continue;
        };
        rows.swap(p, found);
        for i in 0..rows.len() {
            if i != p && rows[i] >> c & 1 != 0 {
                rows[i] ^= rows[p];
                *xors += 1;
            }
        }
        pivots |= 1 << c;
        rank += 1;
        if rank == rows.len() {
            break;
        }
    }
    pivots
}

fn monomial_masks(ell: usize) -> Vec<u64> {
    let mut masks: Vec<_> = (0..ell)
        .flat_map(|i| (0..ell).map(move |j| (1 << i) | (1 << (ell + j))))
        .collect();
    masks.extend((0..2 * ell).map(|i| 1 << i));
    masks.push(0);
    masks
}

fn prepare(
    mut rows: CoefficientRows,
    ell: usize,
    project: bool,
    stats: &mut WeilSolveStats,
) -> Option<Prepared> {
    let nv = 2 * ell;
    let quadratic = ell * ell;
    if !project {
        let masks = monomial_masks(ell);
        let equations = rows.data[..rows.len]
            .iter().copied()
            .filter(|&r| r != 0)
            .map(|r| {
                F2BoolPoly::from_monos(
                    masks
                        .iter()
                        .enumerate()
                        .filter(|(i, _)| r >> i & 1 != 0)
                        .map(|(_, &m)| F2BoolMono::from_mask(m))
                        .collect(),
                    nv,
                )
            })
            .collect();
        let mut expressions = [0; 14];
        for (i, e) in expressions[..nv].iter_mut().enumerate() { *e = 1 << i; }
        return Some(Prepared {
            equations,
            n_vars: nv,
            expressions,
            n_original_vars: nv,
            residual_rows: None,
        });
    }
    let rank = rref(&mut rows.data[..rows.len], quadratic, &mut stats.projection_word_xors).count_ones() as usize;
    let linear = &mut rows.data[rank..rows.len];
    for row in linear.iter_mut() { *row >>= quadratic; }
    let pivots = rref(linear, nv, &mut stats.projection_word_xors);
    let linear_rank = pivots.count_ones() as usize;
    if linear[linear_rank..].contains(&(1 << nv)) {
        return None;
    }
    stats.linear_constraints += linear_rank;
    let free = ((1u64 << nv) - 1) ^ pivots;
    let nf = nv - linear_rank;
    let mut expressions = [0; 14];
    for (i, var) in (0..nv).filter(|i| free >> i & 1 != 0).enumerate() {
        expressions[var] = 1 << i;
    }
    for (row, pivot) in linear.iter().zip((0..nv).filter(|i| pivots >> i & 1 != 0)) {
        expressions[pivot] = ((row >> nv) & 1) << nf;
        for var in (0..nv).filter(|i| free >> i & 1 != 0) {
            if row >> var & 1 != 0 {
                expressions[pivot] ^= expressions[var];
            }
        }
    }
    if nf <= 12 {
        rows.len = rank;
        return Some(Prepared {
            equations: Vec::new(),
            n_vars: nf,
            expressions,
            n_original_vars: nv,
            residual_rows: Some(rows),
        });
    }
    let expand = |expression: u64| -> Vec<u64> {
        (0..=nf)
            .filter(|i| expression >> i & 1 != 0)
            .map(|i| if i == nf { 0 } else { 1 << i })
            .collect()
    };
    let masks = monomial_masks(ell);
    let mut equations = Vec::new();
    for row in &rows.data[..rank] {
        let mut terms = BTreeSet::new();
        for (column, &mask) in masks.iter().enumerate() {
            if row >> column & 1 == 0 {
                continue;
            }
            let mut product = vec![0];
            for var in (0..nv).filter(|i| mask >> i & 1 != 0) {
                let affine = expand(expressions[var]);
                product = product
                    .iter()
                    .flat_map(|&a| affine.iter().map(move |&b| a | b))
                    .collect();
            }
            for term in product {
                if !terms.insert(term) {
                    terms.remove(&term);
                }
            }
        }
        if terms.len() == 1 && terms.contains(&0) {
            return None;
        }
        if !terms.is_empty() {
            equations.push(F2BoolPoly::from_monos(
                terms.into_iter().map(F2BoolMono::from_mask).collect(),
                nf,
            ));
        }
    }
    Some(Prepared {
        equations,
        n_vars: nf,
        expressions,
        n_original_vars: nv,
        residual_rows: None,
    })
}

#[cfg(test)]
mod tests {
    use super::super::koblitz_index_calculus::{
        build_frobenius_union_factor_base, linearised_kernel_basis,
    };
    use super::*;

    #[test]
    fn packed_projection_handles_the_highest_coefficient_bit_and_all_rows() {
        let mut rows = CoefficientRows::new(63);
        rows.extend((0..14).map(|i| 1u64 << (49 + i)));
        assert_eq!(rows.len, 77);
        let prepared = prepare(rows, 7, true, &mut WeilSolveStats::default()).unwrap();
        assert_eq!(prepared.n_vars, 0);
        assert_eq!(prepared.lift(0), 0);
        let mut inconsistent = CoefficientRows::new(63);
        inconsistent.data[0] = 1u64 << 63;
        inconsistent.extend((0..14).map(|i| 1u64 << (49 + i)));
        assert!(prepare(inconsistent, 7, true, &mut WeilSolveStats::default()).is_none());
    }

    #[test]
    fn cached_lifts_match_general_witnesses_and_reject_reordered_bases() {
        let kc = KoblitzCurve::new(1, 7).unwrap();
        let seed: Vec<_> = [1u64, 2].into_iter()
            .map(|x| F2mElement::from_biguint(&BigUint::from(x), 7)).collect();
        let fb = build_frobenius_union_factor_base(&kc, &seed).unwrap();
        let plan = WeilChartPlan::new(&kc, &fb, &seed, true).unwrap();
        let index = fb.index_map();
        for k in 1..71u64 {
            let target = kc.mul(kc.generator(), &BigUint::from(k));
            let mut general = plan.clone();
            general.point_lifts = None;
            let expected = general.decompose(&kc, &fb, &index, &target, &SolveOptions::default()).unwrap();
            let actual = plan.decompose(&kc, &fb, &index, &target, &SolveOptions::default()).unwrap();
            assert_eq!(actual.0, expected.0, "k={k}");
            assert_eq!(actual.1.solver.reductions, expected.1.solver.reductions);
        }
        let mut reordered = fb.clone();
        reordered.points.swap(0, 1);
        assert!(!plan.matches(&kc, &reordered));
        assert!(plan.decompose(&kc, &reordered, &reordered.index_map(), kc.generator(), &SolveOptions::default()).is_none());
        let mut changed = fb.clone();
        changed.points.pop();
        assert!(!plan.matches(&kc, &changed));
    }
    fn truth(plan: &WeilChartPlan, r: u64) -> BTreeSet<(u64, u64)> {
        let f = &plan.field;
        let mut out = BTreeSet::new();
        for (i, &x) in plan.domain.iter().enumerate() {
            for &y in &plan.domain[i..] {
                let p = f.mul(x, y);
                if f.sqr(p) ^ f.mul(r, p) ^ f.mul(f.sqr(r), f.sqr(x ^ y)) ^ plan.b == 0 {
                    out.insert((x, y));
                }
            }
        }
        out
    }
    #[test]
    fn every_chart_pair_matches_exhaustive_s3_with_and_without_projection() {
        for (n, k, a, b) in [(5, 1, 0, 1), (7, 1, 1, 1), (15, 3, 1, 2)] {
            let kc = KoblitzCurve::subfield(k, n, a, b).unwrap();
            let seed: Vec<_> = [1u64, 2]
                .into_iter()
                .map(|x| F2mElement::from_biguint(&BigUint::from(x), n))
                .collect();
            let fb = build_frobenius_union_factor_base(&kc, &seed).unwrap();
            for project in [false, true] {
                let plan = WeilChartPlan::new(&kc, &fb, &seed, project).unwrap();
                for r in 0..1u64 << n.min(5) {
                    let mut actual = BTreeSet::new();
                    let stats = plan.visit_pairs(
                        r,
                        &SolveOptions {
                            node_budget: 100000,
                            ..Default::default()
                        },
                        |x, y| {
                            actual.insert((x, y));
                            false
                        },
                    );
                    assert!(!stats.solver.exhausted);
                    assert_eq!(
                        actual,
                        truth(&plan, r),
                        "n={n} k={k} r={r} projected={project}"
                    );
                    assert_eq!(
                        stats.component_pairs,
                        plan.component_count() * (plan.component_count() + 1) / 2
                    );
                }
            }
        }
    }
    #[test]
    fn scaled_subfield_products_stay_small_and_cover_is_exact() {
        let kc = KoblitzCurve::new(0, 9).unwrap();
        let f = Gf2::new(&kc.curve.irreducible);
        let sub = linearised_kernel_basis(&[0, 3], 9, &kc.curve.irreducible);
        let seed: Vec<_> = sub
            .iter()
            .map(|x| F2mElement::from_biguint(&BigUint::from(f.mul(2, bits(x))), 9))
            .collect();
        let fb = build_frobenius_union_factor_base(&kc, &seed).unwrap();
        let plan = WeilChartPlan::new(&kc, &fb, &seed, true).unwrap();
        assert!(plan.relative_product_ranks().iter().all(|&r| r == 3));
        assert_eq!(
            plan.component_elements()
                .into_iter()
                .flatten()
                .collect::<BTreeSet<_>>(),
            fb.subspace.iter().map(bits).collect()
        );
        for r in [0, 1, 7, 93, 271, 511] {
            let mut actual = BTreeSet::new();
            let stats = plan.visit_pairs(
                r,
                &SolveOptions {
                    node_budget: 100000,
                    ..Default::default()
                },
                |x, y| {
                    actual.insert((x, y));
                    false
                },
            );
            assert!(!stats.solver.exhausted);
            assert_eq!(actual, truth(&plan, r));
        }
        let bad = vec![seed[0].clone(), seed[0].clone()];
        assert!(WeilChartPlan::new(&kc, &fb, &bad, true).is_err());
        assert!(WeilChartPlan::new(&kc, &fb, &sub, true).is_err());
        let other = KoblitzCurve::new(1, 9).unwrap();
        assert!(!plan.matches(&other, &fb));
    }
    #[test]
    fn zero_budget_never_claims_complete_algebraic_search() {
        let kc = KoblitzCurve::new(1, 7).unwrap();
        let seed: Vec<_> = [1u64, 2]
            .into_iter()
            .map(|x| F2mElement::from_biguint(&BigUint::from(x), 7))
            .collect();
        let fb = build_frobenius_union_factor_base(&kc, &seed).unwrap();
        let plan = WeilChartPlan::new(&kc, &fb, &seed, false).unwrap();
        let stats = plan.visit_pairs(
            3,
            &SolveOptions {
                node_budget: 0,
                ..Default::default()
            },
            |_, _| false,
        );
        assert!(stats.solver.exhausted);
        assert_eq!(stats.solver.reductions, 0);
    }

    #[test]
    fn residual_assignment_budget_covers_all_components() {
        let kc = KoblitzCurve::new(1, 7).unwrap();
        let seed: Vec<_> = [1u64, 2, 4].into_iter()
            .map(|x| F2mElement::from_biguint(&BigUint::from(x), 7)).collect();
        let fb = build_frobenius_union_factor_base(&kc, &seed).unwrap();
        let plan = WeilChartPlan::new(&kc, &fb, &seed, true).unwrap();
        let mut exhausted = 0;
        for r in 0..16 {
            let mut roots = BTreeSet::new();
            let stats = plan.visit_pairs(r, &SolveOptions { node_budget: 1, ..Default::default() }, |x, y| {
                roots.insert((x, y));
                false
            });
            assert!(stats.residual_assignments + stats.solver.reductions <= 1);
            if stats.solver.exhausted {
                exhausted += 1;
            } else {
                assert_eq!(roots, truth(&plan, r));
            }
        }
        assert!(exhausted > 0);
    }

    #[test]
    fn rational_support_compression_preserves_every_subgroup_decomposition() {
        let kc = KoblitzCurve::new(1, 11).unwrap();
        let seed: Vec<_> = [1u64, 2, 4].into_iter()
            .map(|x| F2mElement::from_biguint(&BigUint::from(x), 11)).collect();
        let fb = build_frobenius_union_factor_base(&kc, &seed).unwrap();
        let plan = WeilChartPlan::new(&kc, &fb, &seed, true).unwrap();
        assert!(plan.rational_support_dimensions().iter().all(|&d| d == 1));
        assert_eq!(plan.point_chart.as_ref().unwrap().seed_dimension(), 1);
        let fc = FastCurve::new(&kc.curve).unwrap();
        let points: Vec<_> = fb.points.iter().map(|p| fc.lift(p)).collect();
        let expected: HashSet<_> = points.iter().enumerate()
            .flat_map(|(i, &p)| points[i..].iter().map(move |&q| (p,q)))
            .map(|(p,q)| fc.add(p,q)).collect();
        let index = fb.index_map();
        for k in 1..991u64 {
            let target = kc.mul(kc.generator(), &BigUint::from(k));
            let (w, stats) = plan.decompose(&kc, &fb, &index, &target, &SolveOptions::default()).unwrap();
            assert!(!stats.solver.exhausted);
            assert_eq!(w.is_some(), expected.contains(&fc.lift(&target)), "k={k}");
            if let Some(w) = w {
                assert_eq!(kc.add(&fb.points[w[0]], &fb.points[w[1]]), target);
            }
        }
        // The algebraic API still includes nonlifting roots of the full union.
        let mut actual = BTreeSet::new();
        let stats = plan.visit_pairs(13, &SolveOptions::default(), |x,y| { actual.insert((x,y)); false });
        assert!(!stats.solver.exhausted);
        assert_eq!(actual, truth(&plan, 13));
    }

    #[test]
    fn pipeline_rejects_wrong_arity_and_identity_lifts() {
        use super::super::koblitz_index_calculus::{
            koblitz_index_calculus_dlp_with_factor_base, KoblitzIcOptions,
        };
        let kc = KoblitzCurve::new(0, 9).unwrap();
        let seed: Vec<_> = [1u64, 2, 4]
            .into_iter()
            .map(|x| F2mElement::from_biguint(&BigUint::from(x), 9))
            .collect();
        let fb = build_frobenius_union_factor_base(&kc, &seed).unwrap();
        let plan = std::sync::Arc::new(WeilChartPlan::new(&kc, &fb, &seed, true).unwrap());
        let (ids, _) = plan
            .decompose(
                &kc,
                &fb,
                &fb.index_map(),
                &BinaryPoint::Infinity,
                &SolveOptions::default(),
            )
            .unwrap();
        let ids = ids.unwrap();
        assert_eq!(
            kc.add(&fb.points[ids[0]], &fb.points[ids[1]]),
            BinaryPoint::Infinity
        );
        let opts = KoblitzIcOptions {
            m: 3,
            weil_charts: Some(plan),
            ..Default::default()
        };
        assert!(
            koblitz_index_calculus_dlp_with_factor_base(&kc, kc.generator(), &fb, &opts).is_none()
        );
    }
}

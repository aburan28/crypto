//! Indexed S3-image filtering for residual-depth-1 Dickson components.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::path::PathBuf;
use std::process::ExitCode;
use std::time::Duration;

use clap::Parser;
use crypto_lib::cryptanalysis::f4_fp::{self, F4Options, Ordering, Poly};
use crypto_lib::cryptanalysis::p256_dickson_factor_base::CURVE_SLUG;
use crypto_lib::hash::sha256;
use serde::{Deserialize, Serialize};

const P: u64 = 7681;
const A: u64 = P - 3;
const B: u64 = 3506;
const DEPTHS: [u32; 4] = [5, 6, 7, 8];
const TERMINALS: [u64; 2] = [1, 272];
const MAX_DEGREE: u32 = 12;
const BUDGET_SECS: u64 = 5;
const ROUND6_SHA256: &str = "71d63031111ba48ff831e79430bc87e6c7bac34626d7d0de4acfb4eaa65d40f4";

#[derive(Parser)]
#[command(about = "Test an indexed quadratic S3 image filter on frozen Dickson cells")]
struct Cli {
    /// Round-6 canonical result JSON.
    #[arg(long)]
    round6: PathBuf,
    /// Compact deterministic JSON output; stdout when omitted.
    #[arg(long)]
    out: Option<PathBuf>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize)]
struct Affine {
    x: u64,
    y: u64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Point {
    Infinity,
    Affine(Affine),
}

#[derive(Clone, Debug)]
struct Pol {
    n: usize,
    terms: BTreeMap<Vec<u32>, u64>,
}

impl Pol {
    fn zero(n: usize) -> Self {
        Self {
            n,
            terms: BTreeMap::new(),
        }
    }

    fn constant(n: usize, value: u64) -> Self {
        let mut out = Self::zero(n);
        let value = value % P;
        if value != 0 {
            out.terms.insert(vec![0; n], value);
        }
        out
    }

    fn var(n: usize, index: usize) -> Self {
        let mut exponent = vec![0; n];
        exponent[index] = 1;
        let mut out = Self::zero(n);
        out.terms.insert(exponent, 1);
        out
    }

    fn add(&self, other: &Self) -> Self {
        let mut out = self.clone();
        for (exponent, coefficient) in &other.terms {
            let value = out.terms.entry(exponent.clone()).or_insert(0);
            *value = (*value + coefficient) % P;
        }
        out.terms.retain(|_, coefficient| *coefficient != 0);
        out
    }

    fn scale(&self, scalar: u64) -> Self {
        let mut out = Self::zero(self.n);
        for (exponent, coefficient) in &self.terms {
            let value = mulm(*coefficient, scalar);
            if value != 0 {
                out.terms.insert(exponent.clone(), value);
            }
        }
        out
    }

    fn sub(&self, other: &Self) -> Self {
        self.add(&other.scale(P - 1))
    }

    fn mul(&self, other: &Self) -> Self {
        let mut out = Self::zero(self.n);
        for (left_exp, left_coefficient) in &self.terms {
            for (right_exp, right_coefficient) in &other.terms {
                let exponent: Vec<u32> = left_exp
                    .iter()
                    .zip(right_exp)
                    .map(|(left, right)| left + right)
                    .collect();
                let value = out.terms.entry(exponent).or_insert(0);
                *value = (*value + mulm(*left_coefficient, *right_coefficient)) % P;
            }
        }
        out.terms.retain(|_, coefficient| *coefficient != 0);
        out
    }

    fn square(&self) -> Self {
        self.mul(self)
    }

    fn to_f4(&self) -> Poly {
        let terms: Vec<(Vec<u32>, u64)> = self
            .terms
            .iter()
            .map(|(exponent, coefficient)| (exponent.clone(), *coefficient))
            .collect();
        f4_fp::normalise(&terms, P, Ordering::Grevlex)
    }
}

#[derive(Clone, Debug, Deserialize)]
struct Round6Cell {
    terminal: u64,
    depth: u32,
    residual_depth: u32,
    target_kind: String,
    target: [u64; 2],
    components: usize,
    complete_components: usize,
    correct_components: usize,
    max_solving_degree: u32,
    total_field_ops: u64,
}

#[derive(Clone, Debug, Deserialize)]
struct Round6Result {
    schema: String,
    curve: String,
    prime: u64,
    a: u64,
    b: u64,
    cells: Vec<Round6Cell>,
}

#[derive(Clone, Debug, Serialize)]
struct RetainedRun {
    boundary1: u64,
    boundary2: u64,
    inconsistent: bool,
    complete: bool,
    timed_out: bool,
    pairs_above_bound: usize,
    solving_degree_max: u32,
    max_cols_to_solution: usize,
    field_ops: u64,
}

#[derive(Clone, Debug, Serialize)]
struct CellResult {
    terminal: u64,
    depth: u32,
    target_kind: String,
    target: [u64; 2],
    baseline_pairs: usize,
    left_boundaries: usize,
    liftable_factor_base_x: usize,
    quadratic_solves: usize,
    quadratic_roots_returned: usize,
    index_lookups: usize,
    linear_degeneracies: usize,
    universal_degeneracies: usize,
    retained_pairs: usize,
    exact_positive_pairs: usize,
    false_negatives: usize,
    false_positives: usize,
    exact_reference_additions: u64,
    index_build_multiplications: u64,
    filter_multiplications: u64,
    baseline_f4_field_ops: u64,
    retained_f4_field_ops: u64,
    candidate_field_ops: u64,
    candidate_over_baseline: f64,
    retained_max_solving_degree: u32,
    retained_max_columns: usize,
    retained_complete: usize,
    retained_consistent: usize,
    retained_correct: bool,
    retained_pair_sha256: String,
    exact_pair_sha256: String,
    retained_runs: Vec<RetainedRun>,
}

#[derive(Clone, Debug, Serialize)]
struct ExperimentResult {
    schema: String,
    curve: String,
    prime: u64,
    a: u64,
    b: u64,
    round6_sha256: String,
    tonelli_shanks_nonresidue: u64,
    cells: Vec<CellResult>,
}

#[derive(Clone, Copy, Debug)]
struct QuadraticRoots {
    roots: [u64; 2],
    len: usize,
    linear: bool,
    universal: bool,
}

impl QuadraticRoots {
    fn empty() -> Self {
        Self {
            roots: [0; 2],
            len: 0,
            linear: false,
            universal: false,
        }
    }
}

#[derive(Default)]
struct CountedField {
    multiplications: u64,
}

impl CountedField {
    fn mul(&mut self, left: u64, right: u64) -> u64 {
        self.multiplications += 1;
        mulm(left, right)
    }

    fn square(&mut self, value: u64) -> u64 {
        self.mul(value, value)
    }

    fn pow(&mut self, mut base: u64, mut exponent: u64) -> u64 {
        let mut result = 1;
        while exponent != 0 {
            if exponent & 1 == 1 {
                result = self.mul(result, base);
            }
            base = self.square(base);
            exponent >>= 1;
        }
        result
    }

    fn inverse(&mut self, value: u64) -> u64 {
        self.pow(value, P - 2)
    }

    fn sqrt_tonelli_shanks(&mut self, value: u64, nonresidue: u64) -> Option<u64> {
        if value == 0 {
            return Some(0);
        }
        if self.pow(value, (P - 1) / 2) != 1 {
            return None;
        }
        let mut q = P - 1;
        let mut s = 0u32;
        while q & 1 == 0 {
            q >>= 1;
            s += 1;
        }
        let mut c = self.pow(nonresidue, q);
        let mut x = self.pow(value, q.div_ceil(2));
        let mut t = self.pow(value, q);
        let mut m = s;
        while t != 1 {
            let mut i = 1u32;
            let mut power = self.square(t);
            while power != 1 && i < m {
                power = self.square(power);
                i += 1;
            }
            if i == m {
                return None;
            }
            let b = self.pow(c, 1u64 << (m - i - 1));
            x = self.mul(x, b);
            let b2 = self.square(b);
            t = self.mul(t, b2);
            c = b2;
            m = i;
        }
        Some(x)
    }
}

fn addm(left: u64, right: u64) -> u64 {
    (left + right) % P
}

fn negm(value: u64) -> u64 {
    if value == 0 {
        0
    } else {
        P - value
    }
}

fn subm(left: u64, right: u64) -> u64 {
    (left + P - right % P) % P
}

fn mulm(left: u64, right: u64) -> u64 {
    ((left as u128 * right as u128) % P as u128) as u64
}

fn powm(mut base: u64, mut exponent: u64) -> u64 {
    let mut result = 1;
    while exponent != 0 {
        if exponent & 1 == 1 {
            result = mulm(result, base);
        }
        base = mulm(base, base);
        exponent >>= 1;
    }
    result
}

fn inverse(value: u64) -> u64 {
    powm(value, P - 2)
}

fn find_nonresidue() -> u64 {
    (2..P)
        .find(|&value| powm(value, (P - 1) / 2) == P - 1)
        .expect("prime field has a quadratic nonresidue")
}

fn square_root_table() -> Vec<Vec<u64>> {
    let mut table = vec![Vec::new(); P as usize];
    for y in 0..P {
        table[mulm(y, y) as usize].push(y);
    }
    table
}

fn curve_rhs(x: u64) -> u64 {
    addm(addm(mulm(mulm(x, x), x), mulm(A, x)), B)
}

fn curve_points(square_roots: &[Vec<u64>]) -> Vec<Affine> {
    let mut points = Vec::new();
    for x in 0..P {
        for &y in &square_roots[curve_rhs(x) as usize] {
            points.push(Affine { x, y });
        }
    }
    points
}

fn add_points(left: Affine, right: Affine) -> Point {
    let slope = if left.x == right.x {
        if addm(left.y, right.y) == 0 || left.y == 0 {
            return Point::Infinity;
        }
        mulm(
            addm(mulm(3, mulm(left.x, left.x)), A),
            inverse(mulm(2, left.y)),
        )
    } else {
        mulm(subm(right.y, left.y), inverse(subm(right.x, left.x)))
    };
    let x = subm(subm(mulm(slope, slope), left.x), right.x);
    let y = subm(mulm(slope, subm(left.x, x)), left.y);
    Point::Affine(Affine { x, y })
}

fn iterate(mut value: u64, depth: u32) -> u64 {
    for _ in 0..depth {
        value = subm(mulm(value, value), 2);
    }
    value
}

fn roots(depth: u32, terminal: u64) -> Vec<u64> {
    (0..P)
        .filter(|&value| iterate(value, depth) == terminal)
        .collect()
}

fn signed_points(xs: &[u64], square_roots: &[Vec<u64>]) -> Vec<Affine> {
    xs.iter()
        .flat_map(|&x| {
            square_roots[curve_rhs(x) as usize]
                .iter()
                .copied()
                .map(move |y| Affine { x, y })
        })
        .collect()
}

fn component_contains_target(
    left: &[Affine],
    right: &[Affine],
    target: Affine,
    additions: &mut u64,
) -> bool {
    for &p1 in left {
        for &p2 in right {
            *additions += 1;
            if add_points(p1, p2) == Point::Affine(target) {
                return true;
            }
        }
    }
    false
}

#[cfg(test)]
fn s3_value(u: u64, v: u64, t: u64) -> u64 {
    let q = mulm(u, v);
    let sum = addm(u, v);
    let difference = subm(u, v);
    let first = mulm(mulm(t, t), mulm(difference, difference));
    let inner = addm(mulm(sum, addm(q, A)), mulm(2, B));
    let second = negm(mulm(mulm(2, t), inner));
    let third = subm(mulm(subm(q, A), subm(q, A)), mulm(mulm(4, B), sum));
    addm(addm(first, second), third)
}

fn coefficients(field: &mut CountedField, u: u64, t: u64) -> (u64, u64, u64) {
    let t2 = field.square(t);
    let u2 = field.square(u);
    let a2 = field.square(A);
    let au = field.mul(A, u);
    let two_b = addm(B, B);
    let av2 = field.square(subm(t, u));

    let mut inner_b = field.mul(t2, u);
    inner_b = addm(inner_b, field.mul(t, u2));
    inner_b = addm(inner_b, field.mul(t, A));
    inner_b = addm(inner_b, au);
    inner_b = addm(inner_b, two_b);
    let bv = negm(addm(inner_b, inner_b));

    let mut cv = field.mul(t2, u2);
    let t_inner = field.mul(t, addm(au, two_b));
    cv = subm(cv, addm(t_inner, t_inner));
    cv = addm(cv, a2);
    let bu = field.mul(B, u);
    cv = subm(cv, addm(addm(bu, bu), addm(bu, bu)));
    (av2, bv, cv)
}

fn solve_quadratic(
    field: &mut CountedField,
    a: u64,
    b: u64,
    c: u64,
    nonresidue: u64,
) -> QuadraticRoots {
    if a == 0 {
        if b == 0 {
            return QuadraticRoots {
                universal: c == 0,
                ..QuadraticRoots::empty()
            };
        }
        let inverse_b = field.inverse(b);
        let root = field.mul(negm(c), inverse_b);
        return QuadraticRoots {
            roots: [root, 0],
            len: 1,
            linear: true,
            universal: false,
        };
    }
    let b2 = field.square(b);
    let ac = field.mul(a, c);
    let four_ac = addm(addm(ac, ac), addm(ac, ac));
    let discriminant = subm(b2, four_ac);
    let Some(root_discriminant) = field.sqrt_tonelli_shanks(discriminant, nonresidue) else {
        return QuadraticRoots::empty();
    };
    let denominator = addm(a, a);
    let denominator_inverse = field.inverse(denominator);
    let first = field.mul(addm(negm(b), root_discriminant), denominator_inverse);
    let second = field.mul(subm(negm(b), root_discriminant), denominator_inverse);
    if first == second {
        QuadraticRoots {
            roots: [first, 0],
            len: 1,
            linear: false,
            universal: false,
        }
    } else {
        QuadraticRoots {
            roots: [first, second],
            len: 2,
            linear: false,
            universal: false,
        }
    }
}

fn add_chain(equations: &mut Vec<Pol>, n: usize, start: usize, terminal: u64) -> Pol {
    let current = Pol::var(n, start);
    equations.push(
        current
            .square()
            .sub(&Pol::constant(n, 2))
            .sub(&Pol::constant(n, terminal)),
    );
    current
}

fn system(boundary1: u64, boundary2: u64, target_x: u64) -> (Vec<Pol>, usize) {
    let block = 3usize;
    let n = 2 * block + 1;
    let mut equations = Vec::new();
    let x1 = add_chain(&mut equations, n, 0, boundary1);
    let x2 = add_chain(&mut equations, n, block, boundary2);
    let y1 = Pol::var(n, 1);
    let u1 = Pol::var(n, 2);
    let y2 = Pol::var(n, 4);
    let u2 = Pol::var(n, 5);
    for (x, y, u) in [(&x1, &y1, &u1), (&x2, &y2, &u2)] {
        equations.push(u.sub(&x.square()));
        equations.push(
            y.square()
                .sub(&u.mul(x))
                .sub(&x.scale(A))
                .sub(&Pol::constant(n, B)),
        );
    }
    let q = Pol::var(n, 6);
    equations.push(q.sub(&x1.mul(&x2)));
    let sum = x1.add(&x2);
    let difference = x1.sub(&x2);
    let first = difference.square().scale(mulm(target_x, target_x));
    let inner = sum
        .mul(&q.add(&Pol::constant(n, A)))
        .add(&Pol::constant(n, mulm(2, B)));
    let second = inner.scale(negm(mulm(2, target_x)));
    let third = q
        .sub(&Pol::constant(n, A))
        .square()
        .sub(&sum.scale(mulm(4, B)));
    equations.push(first.add(&second).add(&third));
    (equations, n)
}

fn run_retained(boundary1: u64, boundary2: u64, target_x: u64) -> RetainedRun {
    let (equations, n) = system(boundary1, boundary2, target_x);
    let input: Vec<Poly> = equations.iter().map(Pol::to_f4).collect();
    let report = f4_fp::f4(
        &input,
        n,
        P,
        &F4Options::new(Ordering::Grevlex, MAX_DEGREE)
            .with_budget(Duration::from_secs(BUDGET_SECS)),
    );
    let complete =
        !report.timed_out && report.pairs_above_bound == 0 && report.staircase_at_stop.is_none();
    RetainedRun {
        boundary1,
        boundary2,
        inconsistent: report.inconsistent,
        complete,
        timed_out: report.timed_out,
        pairs_above_bound: report.pairs_above_bound,
        solving_degree_max: report.solving_degree_max,
        max_cols_to_solution: report.max_cols_to_solution,
        field_ops: report.field_ops,
    }
}

fn pair_digest(pairs: &BTreeSet<(u64, u64)>) -> String {
    let mut text = String::new();
    for (left, right) in pairs {
        writeln!(text, "{left},{right}").expect("writing to String cannot fail");
    }
    hex::encode(sha256(text.as_bytes()))
}

fn frozen_target(terminal: u64, depth: u32, kind: &str) -> Option<[u64; 2]> {
    let pair = match (terminal, depth) {
        (1, 5) => ([9, 500], [1, 1057]),
        (1, 6) => ([15, 1783], [1, 1057]),
        (1, 7) => ([2, 1636], [1, 1057]),
        (1, 8) => ([1, 1057], [10, 323]),
        (272, 5) => ([1, 1057], [8, 1810]),
        (272, 6) => ([19, 2530], [1, 1057]),
        (272, 7) => ([2, 1636], [1, 1057]),
        (272, 8) => ([1, 1057], [53, 1838]),
        _ => return None,
    };
    match kind {
        "positive" => Some(pair.0),
        "negative" => Some(pair.1),
        _ => None,
    }
}

fn baseline_cell<'a>(
    round6: &'a Round6Result,
    terminal: u64,
    depth: u32,
    kind: &str,
) -> Result<&'a Round6Cell, String> {
    let matches: Vec<&Round6Cell> = round6
        .cells
        .iter()
        .filter(|cell| {
            cell.terminal == terminal
                && cell.depth == depth
                && cell.residual_depth == 1
                && cell.target_kind == kind
        })
        .collect();
    if matches.len() != 1 {
        return Err(format!(
            "expected one round-6 cell for terminal {terminal}, depth {depth}, {kind}; got {}",
            matches.len()
        ));
    }
    Ok(matches[0])
}

fn run_cell(
    baseline: &Round6Cell,
    square_roots: &[Vec<u64>],
    nonresidue: u64,
) -> Result<CellResult, String> {
    let target = Affine {
        x: baseline.target[0],
        y: baseline.target[1],
    };
    let branch_depth = baseline.depth - 1;
    let boundaries = roots(branch_depth, baseline.terminal);
    let expected_boundaries = 1usize << branch_depth;
    if boundaries.len() != expected_boundaries {
        return Err(format!(
            "terminal {}, depth {}: expected {expected_boundaries} boundaries, got {}",
            baseline.terminal,
            baseline.depth,
            boundaries.len()
        ));
    }

    let full_x: Vec<u64> = roots(baseline.depth, baseline.terminal)
        .into_iter()
        .filter(|&x| !square_roots[curve_rhs(x) as usize].is_empty())
        .collect();
    let mut index_build_multiplications = 0u64;
    let mut right_index = BTreeMap::new();
    let mut x_by_boundary: BTreeMap<u64, Vec<u64>> = boundaries
        .iter()
        .map(|&boundary| (boundary, Vec::new()))
        .collect();
    for x in full_x.iter().copied() {
        index_build_multiplications += 1;
        let boundary = subm(mulm(x, x), 2);
        if !x_by_boundary.contains_key(&boundary) {
            return Err(format!(
                "x={x} maps outside the frozen boundary set at terminal {}, depth {}",
                baseline.terminal, baseline.depth
            ));
        }
        right_index.insert(x, boundary);
        x_by_boundary.get_mut(&boundary).unwrap().push(x);
    }

    let point_sets: BTreeMap<u64, Vec<Affine>> = x_by_boundary
        .iter()
        .map(|(&boundary, xs)| (boundary, signed_points(xs, square_roots)))
        .collect();
    let mut exact_pairs = BTreeSet::new();
    let mut exact_reference_additions = 0u64;
    for &left in &boundaries {
        for &right in &boundaries {
            if component_contains_target(
                &point_sets[&left],
                &point_sets[&right],
                target,
                &mut exact_reference_additions,
            ) {
                exact_pairs.insert((left, right));
            }
        }
    }

    let mut field = CountedField::default();
    let mut retained_pairs = BTreeSet::new();
    let mut quadratic_solves = 0usize;
    let mut quadratic_roots_returned = 0usize;
    let mut index_lookups = 0usize;
    let mut linear_degeneracies = 0usize;
    let mut universal_degeneracies = 0usize;
    for &left in &boundaries {
        for &u in &x_by_boundary[&left] {
            quadratic_solves += 1;
            let (qa, qb, qc) = coefficients(&mut field, u, target.x);
            let roots = solve_quadratic(&mut field, qa, qb, qc, nonresidue);
            linear_degeneracies += usize::from(roots.linear);
            universal_degeneracies += usize::from(roots.universal);
            if roots.universal {
                for &right in &boundaries {
                    retained_pairs.insert((left, right));
                }
                continue;
            }
            for &v in &roots.roots[..roots.len] {
                quadratic_roots_returned += 1;
                index_lookups += 1;
                if let Some(&right) = right_index.get(&v) {
                    retained_pairs.insert((left, right));
                }
            }
        }
    }

    let false_negatives = exact_pairs.difference(&retained_pairs).count();
    let false_positives = retained_pairs.difference(&exact_pairs).count();
    let retained_runs: Vec<RetainedRun> = retained_pairs
        .iter()
        .map(|&(left, right)| run_retained(left, right, target.x))
        .collect();
    let retained_f4_field_ops = retained_runs.iter().map(|run| run.field_ops).sum();
    let candidate_field_ops =
        index_build_multiplications + field.multiplications + retained_f4_field_ops;
    let retained_complete = retained_runs.iter().filter(|run| run.complete).count();
    let retained_consistent = retained_runs.iter().filter(|run| !run.inconsistent).count();
    let retained_max_solving_degree = retained_runs
        .iter()
        .map(|run| run.solving_degree_max)
        .max()
        .unwrap_or(0);
    let retained_max_columns = retained_runs
        .iter()
        .map(|run| run.max_cols_to_solution)
        .max()
        .unwrap_or(0);
    let retained_correct = false_negatives == 0
        && false_positives == 0
        && universal_degeneracies == 0
        && retained_complete == retained_runs.len()
        && retained_consistent == retained_runs.len()
        && retained_max_solving_degree <= 3;

    Ok(CellResult {
        terminal: baseline.terminal,
        depth: baseline.depth,
        target_kind: baseline.target_kind.clone(),
        target: baseline.target,
        baseline_pairs: baseline.components,
        left_boundaries: boundaries.len(),
        liftable_factor_base_x: full_x.len(),
        quadratic_solves,
        quadratic_roots_returned,
        index_lookups,
        linear_degeneracies,
        universal_degeneracies,
        retained_pairs: retained_pairs.len(),
        exact_positive_pairs: exact_pairs.len(),
        false_negatives,
        false_positives,
        exact_reference_additions,
        index_build_multiplications,
        filter_multiplications: field.multiplications,
        baseline_f4_field_ops: baseline.total_field_ops,
        retained_f4_field_ops,
        candidate_field_ops,
        candidate_over_baseline: candidate_field_ops as f64 / baseline.total_field_ops as f64,
        retained_max_solving_degree,
        retained_max_columns,
        retained_complete,
        retained_consistent,
        retained_correct,
        retained_pair_sha256: pair_digest(&retained_pairs),
        exact_pair_sha256: pair_digest(&exact_pairs),
        retained_runs,
    })
}

fn validate_round6(round6: &Round6Result) -> Result<(), String> {
    if round6.schema != "p256.dickson_residual_scaling/v1"
        || round6.curve != CURVE_SLUG
        || round6.prime != P
        || round6.a != A
        || round6.b != B
    {
        return Err("round-6 header does not match the frozen corpus".into());
    }
    for terminal in TERMINALS {
        for depth in DEPTHS {
            for kind in ["positive", "negative"] {
                let cell = baseline_cell(round6, terminal, depth, kind)?;
                if Some(cell.target) != frozen_target(terminal, depth, kind)
                    || cell.components != 1usize << (2 * (depth - 1))
                    || cell.complete_components != cell.components
                    || cell.correct_components != cell.components
                    || cell.max_solving_degree != 3
                    || cell.total_field_ops == 0
                {
                    return Err(format!(
                        "round-6 cell failed frozen admission: terminal {terminal}, depth {depth}, {kind}"
                    ));
                }
            }
        }
    }
    Ok(())
}

fn run(cli: Cli) -> Result<(), String> {
    let round6_bytes = std::fs::read(&cli.round6).map_err(|error| error.to_string())?;
    let round6_sha256 = hex::encode(sha256(&round6_bytes));
    if round6_sha256 != ROUND6_SHA256 {
        return Err(format!(
            "round-6 SHA-256 mismatch: expected {ROUND6_SHA256}, got {round6_sha256}"
        ));
    }
    let round6: Round6Result =
        serde_json::from_slice(&round6_bytes).map_err(|error| error.to_string())?;
    validate_round6(&round6)?;
    let square_roots = square_root_table();
    let all_points = curve_points(&square_roots);
    if all_points.len() != 7520 {
        return Err(format!(
            "toy curve point count changed: expected 7520, got {}",
            all_points.len()
        ));
    }
    let nonresidue = find_nonresidue();

    let mut cells = Vec::new();
    for terminal in TERMINALS {
        for depth in DEPTHS {
            for kind in ["positive", "negative"] {
                cells.push(run_cell(
                    baseline_cell(&round6, terminal, depth, kind)?,
                    &square_roots,
                    nonresidue,
                )?);
            }
        }
    }
    let result = ExperimentResult {
        schema: "p256.dickson_s3_image_filter/v1".into(),
        curve: CURVE_SLUG.into(),
        prime: P,
        a: A,
        b: B,
        round6_sha256,
        tonelli_shanks_nonresidue: nonresidue,
        cells,
    };
    let text = serde_json::to_string_pretty(&result).map_err(|error| error.to_string())? + "\n";
    match cli.out {
        Some(path) => std::fs::write(path, text).map_err(|error| error.to_string())?,
        None => print!("{text}"),
    }
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quadratic_coefficients_equal_direct_s3() {
        let nonresidue = find_nonresidue();
        assert_eq!(powm(nonresidue, (P - 1) / 2), P - 1);
        for (u, t) in [(1, 2), (17, 301), (500, 500), (7679, 42)] {
            let mut field = CountedField::default();
            let (qa, qb, qc) = coefficients(&mut field, u, t);
            for v in [0, 1, 2, 99, 1024, 7680] {
                let polynomial = addm(addm(mulm(qa, mulm(v, v)), mulm(qb, v)), qc);
                assert_eq!(polynomial, s3_value(u, v, t));
            }
        }
    }

    #[test]
    fn counted_tonelli_shanks_recovers_all_squares() {
        let nonresidue = find_nonresidue();
        for value in 0..P {
            let mut field = CountedField::default();
            let square = mulm(value, value);
            let root = field.sqrt_tonelli_shanks(square, nonresidue).unwrap();
            assert_eq!(mulm(root, root), square);
        }
    }

    #[test]
    fn quadratic_solver_returns_exact_roots() {
        let nonresidue = find_nonresidue();
        for (u, t) in [(1, 2), (17, 301), (500, 500), (7679, 42)] {
            let mut field = CountedField::default();
            let (qa, qb, qc) = coefficients(&mut field, u, t);
            let roots = solve_quadratic(&mut field, qa, qb, qc, nonresidue);
            if roots.universal {
                for v in 0..P {
                    assert_eq!(s3_value(u, v, t), 0);
                }
            } else {
                let expected: BTreeSet<u64> = (0..P).filter(|&v| s3_value(u, v, t) == 0).collect();
                let actual: BTreeSet<u64> = roots.roots[..roots.len].iter().copied().collect();
                assert_eq!(actual, expected);
            }
        }
    }
}

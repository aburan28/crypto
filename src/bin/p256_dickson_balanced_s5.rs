//! Balanced four-summand degree screen on branched Dickson fibres.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::path::PathBuf;
use std::process::ExitCode;
use std::time::Duration;

use clap::Parser;
use crypto_lib::cryptanalysis::f4_fp::{self, F4Options, Ordering, Poly as F4Poly};
use crypto_lib::cryptanalysis::p256_dickson_factor_base::CURVE_SLUG;
use crypto_lib::ecc::CurveParams;
use crypto_lib::hash::sha256;
use num_traits::ToPrimitive;
use serde::{Deserialize, Serialize};

const P: u64 = 1151;
const A: u64 = P - 3;
const DEPTH: u32 = 5;
const TERMINALS: [u64; 2] = [0, 369];
const CELLS_PER_CLASS: usize = 4;
const MAX_DEGREE: u32 = 8;
const BUDGET_SECS: u64 = 20;
const ROUND9_SHA256: &str = "44f8412ff1c4dec61d0fd05bc2d66f8788e14da4063239e7431e6410da67565a";
const ROUND10_SHA256: &str = "1e14e88f624bb1dadee3ea398f33d4fe0c51b9ea0956d8a4377bbccd6810bc93";
const ROUND11_SHA256: &str = "7266c05093c815734db0a0319b946f76d003e836cce0911541be5ab4d50f2e60";
const IMAGE_MAX_DEGREE: u32 = 4;
const IMAGE_BUDGET_SECS: u64 = 5;

#[derive(Parser)]
#[command(about = "Compare balanced S5 F4 with exact liftability specialisation")]
struct Cli {
    /// Compact deterministic JSON output; stdout when omitted.
    #[arg(long)]
    out: Option<PathBuf>,
    /// Run the round-10 atomized replay from this hash-pinned round-9 result.
    #[arg(long)]
    atomize_round9: Option<PathBuf>,
    /// Run the round-11 pair-image atomization from this round-10 result.
    #[arg(long)]
    image_atomize_round10: Option<PathBuf>,
    /// Run the round-12 eight-summand image-growth screen after hash-checking round 11.
    #[arg(long)]
    scale_eight_round11: Option<PathBuf>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize)]
struct Affine {
    x: u64,
    y: u64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
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
            *value = addm(*value, *coefficient);
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
                *value = addm(*value, mulm(*left_coefficient, *right_coefficient));
            }
        }
        out.terms.retain(|_, coefficient| *coefficient != 0);
        out
    }

    fn square(&self) -> Self {
        self.mul(self)
    }

    fn degree(&self) -> u32 {
        self.terms
            .keys()
            .map(|exponent| exponent.iter().sum())
            .max()
            .unwrap_or(0)
    }

    #[cfg(test)]
    fn evaluate(&self, values: &[u64]) -> u64 {
        self.terms.iter().fold(0, |sum, (exponents, coefficient)| {
            let term = exponents
                .iter()
                .zip(values)
                .fold(*coefficient, |value, (&exponent, &variable)| {
                    mulm(value, powm(variable, u64::from(exponent)))
                });
            addm(sum, term)
        })
    }

    fn to_f4(&self) -> F4Poly {
        let terms: Vec<(Vec<u32>, u64)> = self
            .terms
            .iter()
            .map(|(exponent, coefficient)| (exponent.clone(), *coefficient))
            .collect();
        f4_fp::normalise(&terms, P, Ordering::Grevlex)
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct F4Run {
    arm: String,
    variables: usize,
    equations: usize,
    input_max_degree: u32,
    inconsistent: bool,
    complete: bool,
    correct: bool,
    timed_out: bool,
    pairs_above_bound: usize,
    solving_degree_max: u32,
    max_cols_to_solution: usize,
    field_ops: u64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct StagedRun {
    positive: bool,
    correct: bool,
    pair_image_solves: u64,
    join_solves: u64,
    roots_returned: u64,
    index_lookups: u64,
    field_multiplications: u64,
    maximum_local_degree: u32,
    left_image_size: usize,
    right_image_size: usize,
    witness: Option<[u64; 2]>,
    image_sha256: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct CellResult {
    boundaries: [u64; 4],
    expected_positive: bool,
    exhaustive_signed_tuples: usize,
    baseline: F4Run,
    specialised: F4Run,
    staged: StagedRun,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct BoundaryRecord {
    boundary: u64,
    liftable_x: Vec<u64>,
    signed_points: usize,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct TerminalResult {
    terminal: u64,
    boundaries: Vec<BoundaryRecord>,
    planted_source_boundaries: [u64; 4],
    target: [u64; 2],
    selected_positive_cells: usize,
    selected_negative_cells: usize,
    cells: Vec<CellResult>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct ExperimentResult {
    schema: String,
    curve: String,
    prime: u64,
    a: u64,
    b: u64,
    depth: u32,
    terminals: Vec<u64>,
    cells_per_class: usize,
    max_degree: u32,
    budget_seconds_per_arm_cell: u64,
    terminal_results: Vec<TerminalResult>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct AtomReference {
    affine: bool,
    identity: bool,
    full: bool,
    signed_tuples: usize,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct AtomRun {
    x: [u64; 4],
    reference: AtomReference,
    identity_chart_positive: bool,
    identity_s3_checks: u64,
    identity_chart_correct: bool,
    affine_f4: F4Run,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct AtomizedStagedRun {
    positive: bool,
    correct: bool,
    pair_image_solves: u64,
    join_solves: u64,
    roots_returned: u64,
    index_lookups: u64,
    field_multiplications: u64,
    left_identity: bool,
    right_identity: bool,
    identity_index_checks: u64,
    witness_chart: Option<String>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct AtomizedCellResult {
    boundaries: [u64; 4],
    expected_positive: bool,
    atoms: usize,
    signed_tuples: usize,
    affine_positive: bool,
    identity_positive: bool,
    combined_positive: bool,
    false_negative: bool,
    false_positive: bool,
    atom_runs: Vec<AtomRun>,
    staged: AtomizedStagedRun,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct AtomizedTerminalResult {
    terminal: u64,
    target: [u64; 2],
    parent_cells: usize,
    atomized_subcomponents: usize,
    signed_tuples: usize,
    cells: Vec<AtomizedCellResult>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct AtomizedExperimentResult {
    schema: String,
    curve: String,
    prime: u64,
    a: u64,
    b: u64,
    round9_sha256: String,
    max_degree: u32,
    budget_seconds_per_atom: u64,
    terminal_results: Vec<AtomizedTerminalResult>,
}

#[derive(Clone, Debug, Serialize)]
struct ImageAtomRun {
    left_image_x: u64,
    expected_positive: bool,
    coefficient_multiplications: u64,
    f4: F4Run,
}

#[derive(Clone, Debug, Serialize)]
struct ImageLeafRun {
    leaf_x: [u64; 4],
    left_affine_image: Vec<u64>,
    right_affine_image: Vec<u64>,
    image_build_solves: u64,
    image_build_roots: u64,
    image_build_multiplications: u64,
    identity_chart_positive: bool,
    identity_s3_checks: u64,
    image_atoms: Vec<ImageAtomRun>,
}

#[derive(Clone, Debug, Serialize)]
struct ImageCellResult {
    boundaries: [u64; 4],
    expected_positive: bool,
    leaf_atoms: usize,
    image_atoms: usize,
    signed_tuples: usize,
    affine_positive: bool,
    identity_positive: bool,
    combined_positive: bool,
    false_negative: bool,
    false_positive: bool,
    leaf_runs: Vec<ImageLeafRun>,
}

#[derive(Clone, Debug, Serialize)]
struct ImageTerminalResult {
    terminal: u64,
    target: [u64; 2],
    parent_cells: usize,
    leaf_atoms: usize,
    image_atoms: usize,
    signed_tuples: usize,
    cells: Vec<ImageCellResult>,
}

#[derive(Clone, Debug, Serialize)]
struct ImageExperimentResult {
    schema: String,
    curve: String,
    prime: u64,
    a: u64,
    b: u64,
    round10_sha256: String,
    max_degree: u32,
    budget_seconds_per_image_atom: u64,
    terminal_results: Vec<ImageTerminalResult>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ImageState {
    affine: BTreeSet<u64>,
    identity: bool,
}

#[derive(Clone, Debug, Serialize)]
struct S9CellResult {
    kind: String,
    boundaries: [u64; 8],
    expected_positive: bool,
    candidate_positive: bool,
    x_atoms: usize,
    signed_tuples: usize,
    quadratic_solves: u64,
    roots_returned: u64,
    final_index_lookups: u64,
    field_multiplications: u64,
    maximum_two_leaf_affine_image: usize,
    maximum_four_leaf_affine_image: usize,
    four_leaf_identity_images: u64,
    false_negative: bool,
    false_positive: bool,
    image_stream_sha256: String,
}

#[derive(Clone, Debug, Serialize)]
struct S9TerminalResult {
    terminal: u64,
    target: [u64; 2],
    nonempty_boundaries: usize,
    selected_boundaries: [u64; 8],
    cells: Vec<S9CellResult>,
}

#[derive(Clone, Debug, Serialize)]
struct S9ExperimentResult {
    schema: String,
    curve: String,
    prime: u64,
    a: u64,
    b: u64,
    round11_sha256: String,
    local_maximum_degree: u32,
    terminal_results: Vec<S9TerminalResult>,
}

#[derive(Clone, Copy)]
struct QuadraticRoots {
    roots: [u64; 2],
    len: usize,
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
}

fn addm(left: u64, right: u64) -> u64 {
    (left + right) % P
}

fn subm(left: u64, right: u64) -> u64 {
    (left + P - right % P) % P
}

fn mulm(left: u64, right: u64) -> u64 {
    ((left as u128 * right as u128) % P as u128) as u64
}

fn negm(value: u64) -> u64 {
    if value == 0 {
        0
    } else {
        P - value
    }
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

fn square_root_table() -> Vec<Vec<u64>> {
    let mut table = vec![Vec::new(); P as usize];
    for value in 0..P {
        table[mulm(value, value) as usize].push(value);
    }
    table
}

fn curve_rhs(x: u64, b: u64) -> u64 {
    addm(addm(mulm(mulm(x, x), x), mulm(A, x)), b)
}

fn add_points(left: Point, right: Point) -> Point {
    let (left, right) = match (left, right) {
        (Point::Infinity, value) | (value, Point::Infinity) => return value,
        (Point::Affine(left), Point::Affine(right)) => (left, right),
    };
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

fn negate(point: Affine) -> Affine {
    Affine {
        x: point.x,
        y: negm(point.y),
    }
}

fn add_four(points: [Affine; 4]) -> Point {
    points.into_iter().fold(Point::Infinity, |sum, point| {
        add_points(sum, Point::Affine(point))
    })
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

fn liftable_x(boundary: u64, b: u64, square_roots: &[Vec<u64>]) -> Vec<u64> {
    roots(1, boundary)
        .into_iter()
        .filter(|&x| !square_roots[curve_rhs(x, b) as usize].is_empty())
        .collect()
}

fn signed_points(xs: &[u64], b: u64, square_roots: &[Vec<u64>]) -> Vec<Affine> {
    xs.iter()
        .flat_map(|&x| {
            square_roots[curve_rhs(x, b) as usize]
                .iter()
                .copied()
                .map(move |y| Affine { x, y })
        })
        .collect()
}

fn classify(points: [&[Affine]; 4], target: Affine) -> (bool, usize) {
    let mut positive = false;
    let mut tuples = 0;
    for &p1 in points[0] {
        for &p2 in points[1] {
            for &p3 in points[2] {
                for &p4 in points[3] {
                    tuples += 1;
                    let sum = add_four([p1, p2, p3, p4]);
                    positive |=
                        sum == Point::Affine(target) || sum == Point::Affine(negate(target));
                }
            }
        }
    }
    (positive, tuples)
}

fn planted_target(nonempty: &[(u64, Vec<Affine>)]) -> Result<([u64; 4], Affine), String> {
    if nonempty.len() < 4 {
        return Err(format!(
            "need four nonempty boundaries, found {}",
            nonempty.len()
        ));
    }
    let chosen = [&nonempty[0], &nonempty[1], &nonempty[2], &nonempty[3]];
    for &p1 in &chosen[0].1 {
        for &p2 in &chosen[1].1 {
            for &p3 in &chosen[2].1 {
                for &p4 in &chosen[3].1 {
                    if let Point::Affine(target) = add_four([p1, p2, p3, p4]) {
                        return Ok(([chosen[0].0, chosen[1].0, chosen[2].0, chosen[3].0], target));
                    }
                }
            }
        }
    }
    Err("all planted signed-row tuples summed to infinity".into())
}

fn variable_s3(
    equations: &mut Vec<Pol>,
    x: &Pol,
    y: &Pol,
    z: &Pol,
    s: &Pol,
    q: &Pol,
    h: &Pol,
    k: &Pol,
    b: u64,
) {
    let n = x.n;
    equations.push(s.sub(&x.add(y)));
    equations.push(q.sub(&x.mul(y)));
    equations.push(h.sub(&z.square()));
    equations.push(k.sub(&z.mul(s)));
    let equation = k
        .square()
        .sub(&h.mul(q).scale(4))
        .sub(&k.mul(&q.add(&Pol::constant(n, A))).scale(2))
        .sub(&z.scale(mulm(4, b)))
        .add(&q.sub(&Pol::constant(n, A)).square())
        .sub(&s.scale(mulm(4, b)));
    equations.push(equation);
}

fn final_s3(equations: &mut Vec<Pol>, left: &Pol, right: &Pol, q: &Pol, target_x: u64, b: u64) {
    let n = left.n;
    equations.push(q.sub(&left.mul(right)));
    let sum = left.add(right);
    let difference = left.sub(right);
    let first = difference.square().scale(mulm(target_x, target_x));
    let inner = sum
        .mul(&q.add(&Pol::constant(n, A)))
        .add(&Pol::constant(n, mulm(2, b)));
    let second = inner.scale(negm(mulm(2, target_x)));
    let third = q
        .sub(&Pol::constant(n, A))
        .square()
        .sub(&sum.scale(mulm(4, b)));
    equations.push(first.add(&second).add(&third));
}

fn baseline_system(boundaries: [u64; 4], target_x: u64, b: u64) -> (Vec<Pol>, usize) {
    let n = 23;
    let mut equations = Vec::new();
    let mut xs = Vec::new();
    for (leaf, boundary) in boundaries.into_iter().enumerate() {
        let x = Pol::var(n, 3 * leaf);
        let y = Pol::var(n, 3 * leaf + 1);
        let u = Pol::var(n, 3 * leaf + 2);
        equations.push(
            x.square()
                .sub(&Pol::constant(n, 2))
                .sub(&Pol::constant(n, boundary)),
        );
        equations.push(u.sub(&x.square()));
        equations.push(
            y.square()
                .sub(&u.mul(&x))
                .sub(&x.scale(A))
                .sub(&Pol::constant(n, b)),
        );
        xs.push(x);
    }
    let zl = Pol::var(n, 12);
    let sl = Pol::var(n, 13);
    let ql = Pol::var(n, 14);
    let hl = Pol::var(n, 15);
    let kl = Pol::var(n, 16);
    variable_s3(&mut equations, &xs[0], &xs[1], &zl, &sl, &ql, &hl, &kl, b);
    let zr = Pol::var(n, 17);
    let sr = Pol::var(n, 18);
    let qr = Pol::var(n, 19);
    let hr = Pol::var(n, 20);
    let kr = Pol::var(n, 21);
    variable_s3(&mut equations, &xs[2], &xs[3], &zr, &sr, &qr, &hr, &kr, b);
    final_s3(&mut equations, &zl, &zr, &Pol::var(n, 22), target_x, b);
    (equations, n)
}

fn leaf_polynomial(n: usize, variable: usize, values: &[u64]) -> Result<Pol, String> {
    let x = Pol::var(n, variable);
    match values {
        [root] => Ok(x.sub(&Pol::constant(n, *root))),
        [first, second] => Ok(x
            .sub(&Pol::constant(n, *first))
            .mul(&x.sub(&Pol::constant(n, *second)))),
        _ => Err(format!(
            "residual-depth-1 leaf has {} liftable roots",
            values.len()
        )),
    }
}

fn specialised_system(
    leaves: [&[u64]; 4],
    target_x: u64,
    b: u64,
) -> Result<(Vec<Pol>, usize), String> {
    let n = 15;
    let mut equations = Vec::new();
    let xs: Vec<Pol> = (0..4).map(|index| Pol::var(n, index)).collect();
    for (index, values) in leaves.into_iter().enumerate() {
        equations.push(leaf_polynomial(n, index, values)?);
    }
    let zl = Pol::var(n, 4);
    let sl = Pol::var(n, 5);
    let ql = Pol::var(n, 6);
    let hl = Pol::var(n, 7);
    let kl = Pol::var(n, 8);
    variable_s3(&mut equations, &xs[0], &xs[1], &zl, &sl, &ql, &hl, &kl, b);
    let zr = Pol::var(n, 9);
    let sr = Pol::var(n, 10);
    let qr = Pol::var(n, 11);
    let hr = Pol::var(n, 12);
    let kr = Pol::var(n, 13);
    variable_s3(&mut equations, &xs[2], &xs[3], &zr, &sr, &qr, &hr, &kr, b);
    final_s3(&mut equations, &zl, &zr, &Pol::var(n, 14), target_x, b);
    Ok((equations, n))
}

fn run_f4(arm: &str, equations: &[Pol], n: usize, expected_positive: bool) -> F4Run {
    run_f4_bounded(
        arm,
        equations,
        n,
        expected_positive,
        MAX_DEGREE,
        BUDGET_SECS,
    )
}

fn run_f4_bounded(
    arm: &str,
    equations: &[Pol],
    n: usize,
    expected_positive: bool,
    max_degree: u32,
    budget_secs: u64,
) -> F4Run {
    let input: Vec<F4Poly> = equations.iter().map(Pol::to_f4).collect();
    let report = f4_fp::f4(
        &input,
        n,
        P,
        &F4Options::new(Ordering::Grevlex, max_degree)
            .with_budget(Duration::from_secs(budget_secs)),
    );
    let complete =
        !report.timed_out && report.pairs_above_bound == 0 && report.staircase_at_stop.is_none();
    F4Run {
        arm: arm.into(),
        variables: n,
        equations: equations.len(),
        input_max_degree: equations.iter().map(Pol::degree).max().unwrap_or(0),
        inconsistent: report.inconsistent,
        complete,
        correct: complete && ((!report.inconsistent) == expected_positive),
        timed_out: report.timed_out,
        pairs_above_bound: report.pairs_above_bound,
        solving_degree_max: report.solving_degree_max,
        max_cols_to_solution: report.max_cols_to_solution,
        field_ops: report.field_ops,
    }
}

fn coefficients(field: &mut CountedField, u: u64, t: u64, b: u64) -> (u64, u64, u64) {
    let t2 = field.square(t);
    let u2 = field.square(u);
    let a2 = field.square(A);
    let au = field.mul(A, u);
    let qa = field.square(subm(t, u));
    let mut qb_inner = field.mul(t2, u);
    qb_inner = addm(qb_inner, field.mul(t, u2));
    qb_inner = addm(qb_inner, field.mul(t, A));
    qb_inner = addm(qb_inner, au);
    qb_inner = addm(qb_inner, addm(b, b));
    let qb = negm(addm(qb_inner, qb_inner));
    let mut qc = field.mul(t2, u2);
    let t_inner = field.mul(t, addm(au, addm(b, b)));
    qc = subm(qc, addm(t_inner, t_inner));
    qc = addm(qc, a2);
    let bu = field.mul(b, u);
    qc = subm(qc, addm(addm(bu, bu), addm(bu, bu)));
    (qa, qb, qc)
}

fn solve_quadratic(field: &mut CountedField, qa: u64, qb: u64, qc: u64) -> QuadraticRoots {
    if qa == 0 {
        if qb == 0 {
            return QuadraticRoots {
                roots: [0; 2],
                len: 0,
            };
        }
        let inverse_qb = field.pow(qb, P - 2);
        return QuadraticRoots {
            roots: [field.mul(negm(qc), inverse_qb), 0],
            len: 1,
        };
    }
    let qb2 = field.square(qb);
    let ac = field.mul(qa, qc);
    let discriminant = subm(qb2, addm(addm(ac, ac), addm(ac, ac)));
    let sqrt = field.pow(discriminant, (P + 1) / 4);
    if field.square(sqrt) != discriminant {
        return QuadraticRoots {
            roots: [0; 2],
            len: 0,
        };
    }
    let denominator_inverse = field.pow(addm(qa, qa), P - 2);
    let first = field.mul(addm(negm(qb), sqrt), denominator_inverse);
    let second = field.mul(subm(negm(qb), sqrt), denominator_inverse);
    if first == second {
        QuadraticRoots {
            roots: [first, 0],
            len: 1,
        }
    } else {
        QuadraticRoots {
            roots: [first, second],
            len: 2,
        }
    }
}

fn image(
    left: &[u64],
    right: &[u64],
    b: u64,
    field: &mut CountedField,
    solves: &mut u64,
    roots_returned: &mut u64,
) -> BTreeSet<u64> {
    let mut out = BTreeSet::new();
    for &u in left {
        for &v in right {
            *solves += 1;
            let (qa, qb, qc) = coefficients(field, u, v, b);
            let roots = solve_quadratic(field, qa, qb, qc);
            *roots_returned += roots.len as u64;
            out.extend(roots.roots[..roots.len].iter().copied());
        }
    }
    out
}

fn staged_run(leaves: [&[u64]; 4], target_x: u64, b: u64, expected_positive: bool) -> StagedRun {
    let mut field = CountedField::default();
    let mut pair_image_solves = 0;
    let mut roots_returned = 0;
    let left_image = image(
        leaves[0],
        leaves[1],
        b,
        &mut field,
        &mut pair_image_solves,
        &mut roots_returned,
    );
    let right_image = image(
        leaves[2],
        leaves[3],
        b,
        &mut field,
        &mut pair_image_solves,
        &mut roots_returned,
    );
    let mut join_solves = 0;
    let mut index_lookups = 0;
    let mut witness = None;
    for &left in &left_image {
        join_solves += 1;
        let (qa, qb, qc) = coefficients(&mut field, left, target_x, b);
        let roots = solve_quadratic(&mut field, qa, qb, qc);
        roots_returned += roots.len as u64;
        for &right in &roots.roots[..roots.len] {
            index_lookups += 1;
            if witness.is_none() && right_image.contains(&right) {
                witness = Some([left, right]);
            }
        }
    }
    let mut stream = String::new();
    for value in &left_image {
        writeln!(stream, "L,{value}").expect("writing to String cannot fail");
    }
    for value in &right_image {
        writeln!(stream, "R,{value}").expect("writing to String cannot fail");
    }
    if let Some([left, right]) = witness {
        writeln!(stream, "W,{left},{right}").expect("writing to String cannot fail");
    }
    let positive = witness.is_some();
    StagedRun {
        positive,
        correct: positive == expected_positive,
        pair_image_solves,
        join_solves,
        roots_returned,
        index_lookups,
        field_multiplications: field.multiplications,
        maximum_local_degree: 2,
        left_image_size: left_image.len(),
        right_image_size: right_image.len(),
        witness,
        image_sha256: hex::encode(sha256(stream.as_bytes())),
    }
}

fn run_terminal(
    terminal: u64,
    b: u64,
    square_roots: &[Vec<u64>],
) -> Result<TerminalResult, String> {
    let boundaries = roots(DEPTH - 1, terminal);
    if boundaries.len() != 1usize << (DEPTH - 1) {
        return Err(format!(
            "terminal {terminal}: expected 16 boundaries, found {}",
            boundaries.len()
        ));
    }
    let mut x_by_boundary = BTreeMap::new();
    let mut points_by_boundary = BTreeMap::new();
    let mut boundary_records = Vec::new();
    for boundary in boundaries {
        let xs = liftable_x(boundary, b, square_roots);
        let points = signed_points(&xs, b, square_roots);
        boundary_records.push(BoundaryRecord {
            boundary,
            liftable_x: xs.clone(),
            signed_points: points.len(),
        });
        if !xs.is_empty() {
            x_by_boundary.insert(boundary, xs);
            points_by_boundary.insert(boundary, points);
        }
    }
    let nonempty: Vec<(u64, Vec<Affine>)> = points_by_boundary
        .iter()
        .map(|(&boundary, points)| (boundary, points.clone()))
        .collect();
    let (planted_source_boundaries, target) = planted_target(&nonempty)?;
    let ids: Vec<u64> = x_by_boundary.keys().copied().collect();
    let mut positive_cells = Vec::new();
    let mut negative_cells = Vec::new();
    'outer: for &b1 in &ids {
        for &b2 in &ids {
            for &b3 in &ids {
                for &b4 in &ids {
                    let tuple = [b1, b2, b3, b4];
                    let (positive, signed_tuples) = classify(
                        [
                            &points_by_boundary[&b1],
                            &points_by_boundary[&b2],
                            &points_by_boundary[&b3],
                            &points_by_boundary[&b4],
                        ],
                        target,
                    );
                    let entry = (tuple, signed_tuples);
                    if positive && positive_cells.len() < CELLS_PER_CLASS {
                        positive_cells.push(entry);
                    } else if !positive && negative_cells.len() < CELLS_PER_CLASS {
                        negative_cells.push(entry);
                    }
                    if positive_cells.len() == CELLS_PER_CLASS
                        && negative_cells.len() == CELLS_PER_CLASS
                    {
                        break 'outer;
                    }
                }
            }
        }
    }
    if positive_cells.len() != CELLS_PER_CLASS || negative_cells.len() != CELLS_PER_CLASS {
        return Err(format!(
            "terminal {terminal}: selected {} positive and {} negative cells",
            positive_cells.len(),
            negative_cells.len()
        ));
    }
    let mut selected = Vec::with_capacity(2 * CELLS_PER_CLASS);
    selected.extend(positive_cells.into_iter().map(|entry| (true, entry)));
    selected.extend(negative_cells.into_iter().map(|entry| (false, entry)));
    let mut cells = Vec::new();
    for (expected_positive, (boundaries, exhaustive_signed_tuples)) in selected {
        let leaves = [
            x_by_boundary[&boundaries[0]].as_slice(),
            x_by_boundary[&boundaries[1]].as_slice(),
            x_by_boundary[&boundaries[2]].as_slice(),
            x_by_boundary[&boundaries[3]].as_slice(),
        ];
        let (baseline_equations, baseline_n) = baseline_system(boundaries, target.x, b);
        let baseline = run_f4(
            "curve-lift-baseline",
            &baseline_equations,
            baseline_n,
            expected_positive,
        );
        let (specialised_equations, specialised_n) = specialised_system(leaves, target.x, b)?;
        let specialised = run_f4(
            "liftability-specialised",
            &specialised_equations,
            specialised_n,
            expected_positive,
        );
        let staged = staged_run(leaves, target.x, b, expected_positive);
        cells.push(CellResult {
            boundaries,
            expected_positive,
            exhaustive_signed_tuples,
            baseline,
            specialised,
            staged,
        });
    }
    Ok(TerminalResult {
        terminal,
        boundaries: boundary_records,
        planted_source_boundaries,
        target: [target.x, target.y],
        selected_positive_cells: CELLS_PER_CLASS,
        selected_negative_cells: CELLS_PER_CLASS,
        cells,
    })
}

fn s3_value(x: u64, y: u64, z: u64, b: u64) -> u64 {
    let q = mulm(x, y);
    let sum = addm(x, y);
    let difference = subm(x, y);
    let first = mulm(mulm(z, z), mulm(difference, difference));
    let inner = addm(mulm(sum, addm(q, A)), mulm(2, b));
    let second = negm(mulm(mulm(2, z), inner));
    let third = subm(mulm(subm(q, A), subm(q, A)), mulm(mulm(4, b), sum));
    addm(addm(first, second), third)
}

fn atom_reference(
    xs: [u64; 4],
    target: Affine,
    b: u64,
    square_roots: &[Vec<u64>],
) -> AtomReference {
    let point_sets: Vec<Vec<Affine>> = xs
        .into_iter()
        .map(|x| {
            square_roots[curve_rhs(x, b) as usize]
                .iter()
                .map(|&y| Affine { x, y })
                .collect()
        })
        .collect();
    let mut affine = false;
    let mut identity = false;
    let mut signed_tuples = 0;
    for &p1 in &point_sets[0] {
        for &p2 in &point_sets[1] {
            for &p3 in &point_sets[2] {
                for &p4 in &point_sets[3] {
                    signed_tuples += 1;
                    let left = add_points(Point::Affine(p1), Point::Affine(p2));
                    let right = add_points(Point::Affine(p3), Point::Affine(p4));
                    let sum = add_points(left, right);
                    if sum != Point::Affine(target) && sum != Point::Affine(negate(target)) {
                        continue;
                    }
                    if left == Point::Infinity || right == Point::Infinity {
                        identity = true;
                    } else {
                        affine = true;
                    }
                }
            }
        }
    }
    AtomReference {
        affine,
        identity,
        full: affine || identity,
        signed_tuples,
    }
}

fn identity_chart(xs: [u64; 4], target_x: u64, b: u64) -> (bool, u64) {
    let mut positive = false;
    let mut checks = 0;
    if xs[0] == xs[1] {
        checks += 1;
        positive |= s3_value(xs[2], xs[3], target_x, b) == 0;
    }
    if xs[2] == xs[3] {
        checks += 1;
        positive |= s3_value(xs[0], xs[1], target_x, b) == 0;
    }
    (positive, checks)
}

fn atomized_staged_run(
    leaves: [&[u64]; 4],
    target_x: u64,
    b: u64,
    expected_positive: bool,
) -> AtomizedStagedRun {
    let mut field = CountedField::default();
    let mut pair_image_solves = 0;
    let mut roots_returned = 0;
    let left_image = image(
        leaves[0],
        leaves[1],
        b,
        &mut field,
        &mut pair_image_solves,
        &mut roots_returned,
    );
    let right_image = image(
        leaves[2],
        leaves[3],
        b,
        &mut field,
        &mut pair_image_solves,
        &mut roots_returned,
    );
    let left_identity = leaves[0].iter().any(|x| leaves[1].contains(x));
    let right_identity = leaves[2].iter().any(|x| leaves[3].contains(x));
    let mut join_solves = 0;
    let mut index_lookups = 0;
    let mut witness_chart = None;
    for &left in &left_image {
        join_solves += 1;
        let (qa, qb, qc) = coefficients(&mut field, left, target_x, b);
        let roots = solve_quadratic(&mut field, qa, qb, qc);
        roots_returned += roots.len as u64;
        for &right in &roots.roots[..roots.len] {
            index_lookups += 1;
            if witness_chart.is_none() && right_image.contains(&right) {
                witness_chart = Some("affine-affine".into());
            }
        }
    }
    let mut identity_index_checks = 0;
    if left_identity {
        identity_index_checks += 1;
        if witness_chart.is_none() && right_image.contains(&target_x) {
            witness_chart = Some("left-identity".into());
        }
    }
    if right_identity {
        identity_index_checks += 1;
        if witness_chart.is_none() && left_image.contains(&target_x) {
            witness_chart = Some("right-identity".into());
        }
    }
    let positive = witness_chart.is_some();
    AtomizedStagedRun {
        positive,
        correct: positive == expected_positive,
        pair_image_solves,
        join_solves,
        roots_returned,
        index_lookups,
        field_multiplications: field.multiplications,
        left_identity,
        right_identity,
        identity_index_checks,
        witness_chart,
    }
}

fn run_atomized_terminal(
    terminal: &TerminalResult,
    b: u64,
    square_roots: &[Vec<u64>],
) -> Result<AtomizedTerminalResult, String> {
    let x_by_boundary: BTreeMap<u64, Vec<u64>> = terminal
        .boundaries
        .iter()
        .filter(|record| !record.liftable_x.is_empty())
        .map(|record| (record.boundary, record.liftable_x.clone()))
        .collect();
    let target = Affine {
        x: terminal.target[0],
        y: terminal.target[1],
    };
    let mut atomized_subcomponents = 0;
    let mut terminal_signed_tuples = 0;
    let mut cells = Vec::new();
    for parent in &terminal.cells {
        let roots = [
            &x_by_boundary[&parent.boundaries[0]],
            &x_by_boundary[&parent.boundaries[1]],
            &x_by_boundary[&parent.boundaries[2]],
            &x_by_boundary[&parent.boundaries[3]],
        ];
        let mut atom_runs = Vec::new();
        let mut affine_positive = false;
        let mut identity_positive = false;
        let mut signed_tuples = 0;
        for &x1 in roots[0] {
            for &x2 in roots[1] {
                for &x3 in roots[2] {
                    for &x4 in roots[3] {
                        let xs = [x1, x2, x3, x4];
                        let reference = atom_reference(xs, target, b, square_roots);
                        let (identity_chart_positive, identity_s3_checks) =
                            identity_chart(xs, target.x, b);
                        let singleton = [[x1], [x2], [x3], [x4]];
                        let leaves = [
                            singleton[0].as_slice(),
                            singleton[1].as_slice(),
                            singleton[2].as_slice(),
                            singleton[3].as_slice(),
                        ];
                        let (equations, n) = specialised_system(leaves, target.x, b)?;
                        let affine_f4 = run_f4("atomized-affine", &equations, n, reference.affine);
                        affine_positive |= !affine_f4.inconsistent;
                        identity_positive |= identity_chart_positive;
                        signed_tuples += reference.signed_tuples;
                        atom_runs.push(AtomRun {
                            x: xs,
                            identity_chart_correct: identity_chart_positive == reference.identity,
                            identity_chart_positive,
                            identity_s3_checks,
                            reference,
                            affine_f4,
                        });
                    }
                }
            }
        }
        let combined_positive = affine_positive || identity_positive;
        let leaves = [
            roots[0].as_slice(),
            roots[1].as_slice(),
            roots[2].as_slice(),
            roots[3].as_slice(),
        ];
        let staged = atomized_staged_run(leaves, target.x, b, parent.expected_positive);
        atomized_subcomponents += atom_runs.len();
        terminal_signed_tuples += signed_tuples;
        cells.push(AtomizedCellResult {
            boundaries: parent.boundaries,
            expected_positive: parent.expected_positive,
            atoms: atom_runs.len(),
            signed_tuples,
            affine_positive,
            identity_positive,
            combined_positive,
            false_negative: parent.expected_positive && !combined_positive,
            false_positive: !parent.expected_positive && combined_positive,
            atom_runs,
            staged,
        });
    }
    Ok(AtomizedTerminalResult {
        terminal: terminal.terminal,
        target: terminal.target,
        parent_cells: terminal.cells.len(),
        atomized_subcomponents,
        signed_tuples: terminal_signed_tuples,
        cells,
    })
}

fn run_atomized(round9_path: &PathBuf, out: Option<PathBuf>) -> Result<(), String> {
    let bytes = std::fs::read(round9_path).map_err(|error| error.to_string())?;
    let round9_sha256 = hex::encode(sha256(&bytes));
    if round9_sha256 != ROUND9_SHA256 {
        return Err(format!(
            "round-9 SHA-256 mismatch: expected {ROUND9_SHA256}, got {round9_sha256}"
        ));
    }
    let round9: ExperimentResult =
        serde_json::from_slice(&bytes).map_err(|error| error.to_string())?;
    if round9.schema != "p256.dickson_balanced_s5_degree/v1"
        || round9.curve != CURVE_SLUG
        || round9.prime != P
        || round9.a != A
        || round9.depth != DEPTH
        || round9.terminals != TERMINALS
        || round9.terminal_results.len() != TERMINALS.len()
        || round9
            .terminal_results
            .iter()
            .any(|terminal| terminal.cells.len() != 2 * CELLS_PER_CLASS)
    {
        return Err("round-9 header or frozen cell count mismatch".into());
    }
    let square_roots = square_root_table();
    let mut terminal_results = Vec::new();
    for terminal in &round9.terminal_results {
        terminal_results.push(run_atomized_terminal(terminal, round9.b, &square_roots)?);
    }
    let result = AtomizedExperimentResult {
        schema: "p256.dickson_atomized_s5_degree/v1".into(),
        curve: CURVE_SLUG.into(),
        prime: P,
        a: A,
        b: round9.b,
        round9_sha256,
        max_degree: MAX_DEGREE,
        budget_seconds_per_atom: BUDGET_SECS,
        terminal_results,
    };
    let text = serde_json::to_string_pretty(&result).map_err(|error| error.to_string())? + "\n";
    match out {
        Some(path) => std::fs::write(path, text).map_err(|error| error.to_string())?,
        None => print!("{text}"),
    }
    Ok(())
}

fn affine_pair_image_reference(
    left_x: u64,
    right_x: u64,
    b: u64,
    square_roots: &[Vec<u64>],
) -> BTreeSet<u64> {
    let mut out = BTreeSet::new();
    for &left_y in &square_roots[curve_rhs(left_x, b) as usize] {
        for &right_y in &square_roots[curve_rhs(right_x, b) as usize] {
            let sum = add_points(
                Point::Affine(Affine {
                    x: left_x,
                    y: left_y,
                }),
                Point::Affine(Affine {
                    x: right_x,
                    y: right_y,
                }),
            );
            if let Point::Affine(point) = sum {
                out.insert(point.x);
            }
        }
    }
    out
}

fn counted_pair_image(left_x: u64, right_x: u64, b: u64) -> (BTreeSet<u64>, u64, u64, u64) {
    let mut field = CountedField::default();
    let mut solves = 0;
    let mut roots_returned = 0;
    let image = image(
        &[left_x],
        &[right_x],
        b,
        &mut field,
        &mut solves,
        &mut roots_returned,
    );
    (image, solves, roots_returned, field.multiplications)
}

fn quadratic_polynomial(qa: u64, qb: u64, qc: u64) -> Pol {
    let variable = Pol::var(1, 0);
    variable
        .square()
        .scale(qa)
        .add(&variable.scale(qb))
        .add(&Pol::constant(1, qc))
}

fn run_image_terminal(
    terminal: &AtomizedTerminalResult,
    b: u64,
    square_roots: &[Vec<u64>],
) -> Result<ImageTerminalResult, String> {
    let target_x = terminal.target[0];
    let mut terminal_leaf_atoms = 0;
    let mut terminal_image_atoms = 0;
    let mut cells = Vec::new();
    for parent in &terminal.cells {
        let mut leaf_runs = Vec::new();
        let mut affine_positive = false;
        let mut identity_positive = false;
        let mut image_atoms = 0;
        for atom in &parent.atom_runs {
            let (left_image, left_solves, left_roots, left_muls) =
                counted_pair_image(atom.x[0], atom.x[1], b);
            let (right_image, right_solves, right_roots, right_muls) =
                counted_pair_image(atom.x[2], atom.x[3], b);
            let left_reference = affine_pair_image_reference(atom.x[0], atom.x[1], b, square_roots);
            let right_reference =
                affine_pair_image_reference(atom.x[2], atom.x[3], b, square_roots);
            if left_image != left_reference || right_image != right_reference {
                return Err(format!(
                    "pair image mismatch at terminal {}, leaf atom {:?}",
                    terminal.terminal, atom.x
                ));
            }
            let right_values: Vec<u64> = right_image.iter().copied().collect();
            if right_values.is_empty() || right_values.len() > 2 {
                return Err(format!(
                    "right pair image at {:?} has {} affine roots",
                    atom.x,
                    right_values.len()
                ));
            }
            let right_membership = leaf_polynomial(1, 0, &right_values)?;
            let mut atom_image_runs = Vec::new();
            for &left in &left_image {
                let mut coefficient_field = CountedField::default();
                let (qa, qb, qc) = coefficients(&mut coefficient_field, left, target_x, b);
                let target_polynomial = quadratic_polynomial(qa, qb, qc);
                let expected_positive = right_values
                    .iter()
                    .any(|&right| s3_value(left, right, target_x, b) == 0);
                let f4 = run_f4_bounded(
                    "image-atomized-final-join",
                    &[target_polynomial, right_membership.clone()],
                    1,
                    expected_positive,
                    IMAGE_MAX_DEGREE,
                    IMAGE_BUDGET_SECS,
                );
                affine_positive |= !f4.inconsistent;
                atom_image_runs.push(ImageAtomRun {
                    left_image_x: left,
                    expected_positive,
                    coefficient_multiplications: coefficient_field.multiplications,
                    f4,
                });
            }
            identity_positive |= atom.identity_chart_positive;
            image_atoms += atom_image_runs.len();
            leaf_runs.push(ImageLeafRun {
                leaf_x: atom.x,
                left_affine_image: left_image.iter().copied().collect(),
                right_affine_image: right_values,
                image_build_solves: left_solves + right_solves,
                image_build_roots: left_roots + right_roots,
                image_build_multiplications: left_muls + right_muls,
                identity_chart_positive: atom.identity_chart_positive,
                identity_s3_checks: atom.identity_s3_checks,
                image_atoms: atom_image_runs,
            });
        }
        let combined_positive = affine_positive || identity_positive;
        terminal_leaf_atoms += leaf_runs.len();
        terminal_image_atoms += image_atoms;
        cells.push(ImageCellResult {
            boundaries: parent.boundaries,
            expected_positive: parent.expected_positive,
            leaf_atoms: leaf_runs.len(),
            image_atoms,
            signed_tuples: parent.signed_tuples,
            affine_positive,
            identity_positive,
            combined_positive,
            false_negative: parent.expected_positive && !combined_positive,
            false_positive: !parent.expected_positive && combined_positive,
            leaf_runs,
        });
    }
    Ok(ImageTerminalResult {
        terminal: terminal.terminal,
        target: terminal.target,
        parent_cells: terminal.cells.len(),
        leaf_atoms: terminal_leaf_atoms,
        image_atoms: terminal_image_atoms,
        signed_tuples: terminal.signed_tuples,
        cells,
    })
}

fn run_image_atomized(round10_path: &PathBuf, out: Option<PathBuf>) -> Result<(), String> {
    let bytes = std::fs::read(round10_path).map_err(|error| error.to_string())?;
    let round10_sha256 = hex::encode(sha256(&bytes));
    if round10_sha256 != ROUND10_SHA256 {
        return Err(format!(
            "round-10 SHA-256 mismatch: expected {ROUND10_SHA256}, got {round10_sha256}"
        ));
    }
    let round10: AtomizedExperimentResult =
        serde_json::from_slice(&bytes).map_err(|error| error.to_string())?;
    let total_atoms: usize = round10
        .terminal_results
        .iter()
        .map(|terminal| terminal.atomized_subcomponents)
        .sum();
    let total_signed: usize = round10
        .terminal_results
        .iter()
        .map(|terminal| terminal.signed_tuples)
        .sum();
    if round10.schema != "p256.dickson_atomized_s5_degree/v1"
        || round10.curve != CURVE_SLUG
        || round10.prime != P
        || round10.a != A
        || round10.terminal_results.len() != TERMINALS.len()
        || total_atoms != 104
        || total_signed != 1664
        || round10
            .terminal_results
            .iter()
            .any(|terminal| terminal.cells.len() != 2 * CELLS_PER_CLASS)
    {
        return Err("round-10 header or frozen evidence count mismatch".into());
    }
    let square_roots = square_root_table();
    let mut terminal_results = Vec::new();
    for terminal in &round10.terminal_results {
        terminal_results.push(run_image_terminal(terminal, round10.b, &square_roots)?);
    }
    let result = ImageExperimentResult {
        schema: "p256.dickson_image_atomized_s5_degree/v1".into(),
        curve: CURVE_SLUG.into(),
        prime: P,
        a: A,
        b: round10.b,
        round10_sha256,
        max_degree: IMAGE_MAX_DEGREE,
        budget_seconds_per_image_atom: IMAGE_BUDGET_SECS,
        terminal_results,
    };
    let text = serde_json::to_string_pretty(&result).map_err(|error| error.to_string())? + "\n";
    match out {
        Some(path) => std::fs::write(path, text).map_err(|error| error.to_string())?,
        None => print!("{text}"),
    }
    Ok(())
}

fn sum_sets(point_sets: &[Vec<Affine>]) -> BTreeSet<Point> {
    let mut sums = BTreeSet::from([Point::Infinity]);
    for points in point_sets {
        let mut next = BTreeSet::new();
        for &sum in &sums {
            for &point in points {
                next.insert(add_points(sum, Point::Affine(point)));
            }
        }
        sums = next;
    }
    sums
}

fn first_finite_sum(point_sets: &[Vec<Affine>]) -> Option<Affine> {
    fn visit(index: usize, point_sets: &[Vec<Affine>], sum: Point) -> Option<Affine> {
        if index == point_sets.len() {
            return match sum {
                Point::Infinity => None,
                Point::Affine(point) => Some(point),
            };
        }
        for &point in &point_sets[index] {
            if let Some(result) =
                visit(index + 1, point_sets, add_points(sum, Point::Affine(point)))
            {
                return Some(result);
            }
        }
        None
    }
    visit(0, point_sets, Point::Infinity)
}

fn boundary_cell_reference(
    boundaries: [u64; 8],
    points_by_boundary: &BTreeMap<u64, Vec<Affine>>,
    target: Affine,
) -> (bool, usize) {
    let sets: Vec<Vec<Affine>> = boundaries
        .iter()
        .map(|boundary| points_by_boundary[boundary].clone())
        .collect();
    let signed_tuples = sets.iter().map(Vec::len).product();
    let sums = sum_sets(&sets);
    (
        sums.contains(&Point::Affine(target)) || sums.contains(&Point::Affine(negate(target))),
        signed_tuples,
    )
}

fn enumerate_x_atoms(leaves: [&[u64]; 8]) -> Vec<[u64; 8]> {
    fn visit(index: usize, leaves: &[&[u64]; 8], current: &mut [u64; 8], out: &mut Vec<[u64; 8]>) {
        if index == 8 {
            out.push(*current);
            return;
        }
        for &value in leaves[index] {
            current[index] = value;
            visit(index + 1, leaves, current, out);
        }
    }
    let mut out = Vec::new();
    visit(0, &leaves, &mut [0; 8], &mut out);
    out
}

fn reference_image(xs: &[u64], b: u64, square_roots: &[Vec<u64>]) -> ImageState {
    let point_sets: Vec<Vec<Affine>> = xs
        .iter()
        .map(|&x| {
            square_roots[curve_rhs(x, b) as usize]
                .iter()
                .map(|&y| Affine { x, y })
                .collect()
        })
        .collect();
    let sums = sum_sets(&point_sets);
    ImageState {
        affine: sums
            .iter()
            .filter_map(|point| match point {
                Point::Infinity => None,
                Point::Affine(point) => Some(point.x),
            })
            .collect(),
        identity: sums.contains(&Point::Infinity),
    }
}

fn pair_image_state(
    left: u64,
    right: u64,
    b: u64,
    field: &mut CountedField,
    solves: &mut u64,
    roots_returned: &mut u64,
) -> ImageState {
    *solves += 1;
    let (qa, qb, qc) = coefficients(field, left, right, b);
    let roots = solve_quadratic(field, qa, qb, qc);
    *roots_returned += roots.len as u64;
    ImageState {
        affine: roots.roots[..roots.len].iter().copied().collect(),
        identity: left == right,
    }
}

fn compose_image_states(
    left: &ImageState,
    right: &ImageState,
    b: u64,
    field: &mut CountedField,
    solves: &mut u64,
    roots_returned: &mut u64,
) -> ImageState {
    let mut affine = BTreeSet::new();
    if left.identity {
        affine.extend(right.affine.iter().copied());
    }
    if right.identity {
        affine.extend(left.affine.iter().copied());
    }
    for &u in &left.affine {
        for &v in &right.affine {
            *solves += 1;
            let (qa, qb, qc) = coefficients(field, u, v, b);
            let roots = solve_quadratic(field, qa, qb, qc);
            *roots_returned += roots.len as u64;
            affine.extend(roots.roots[..roots.len].iter().copied());
        }
    }
    let shared_affine = left.affine.iter().any(|x| right.affine.contains(x));
    ImageState {
        affine,
        identity: (left.identity && right.identity) || shared_affine,
    }
}

#[allow(clippy::too_many_arguments)]
fn run_s9_cell(
    kind: &str,
    boundaries: [u64; 8],
    expected_positive: bool,
    signed_tuples: usize,
    target: Affine,
    x_by_boundary: &BTreeMap<u64, Vec<u64>>,
    b: u64,
    square_roots: &[Vec<u64>],
) -> Result<S9CellResult, String> {
    let leaves = [
        x_by_boundary[&boundaries[0]].as_slice(),
        x_by_boundary[&boundaries[1]].as_slice(),
        x_by_boundary[&boundaries[2]].as_slice(),
        x_by_boundary[&boundaries[3]].as_slice(),
        x_by_boundary[&boundaries[4]].as_slice(),
        x_by_boundary[&boundaries[5]].as_slice(),
        x_by_boundary[&boundaries[6]].as_slice(),
        x_by_boundary[&boundaries[7]].as_slice(),
    ];
    let atoms = enumerate_x_atoms(leaves);
    let mut candidate_positive = false;
    let mut quadratic_solves = 0;
    let mut roots_returned = 0;
    let mut final_index_lookups = 0;
    let mut field = CountedField::default();
    let mut maximum_two_leaf_affine_image = 0;
    let mut maximum_four_leaf_affine_image = 0;
    let mut four_leaf_identity_images = 0;
    let mut stream = String::new();
    for xs in &atoms {
        let pairs = [
            pair_image_state(
                xs[0],
                xs[1],
                b,
                &mut field,
                &mut quadratic_solves,
                &mut roots_returned,
            ),
            pair_image_state(
                xs[2],
                xs[3],
                b,
                &mut field,
                &mut quadratic_solves,
                &mut roots_returned,
            ),
            pair_image_state(
                xs[4],
                xs[5],
                b,
                &mut field,
                &mut quadratic_solves,
                &mut roots_returned,
            ),
            pair_image_state(
                xs[6],
                xs[7],
                b,
                &mut field,
                &mut quadratic_solves,
                &mut roots_returned,
            ),
        ];
        maximum_two_leaf_affine_image = maximum_two_leaf_affine_image.max(
            pairs
                .iter()
                .map(|image| image.affine.len())
                .max()
                .unwrap_or(0),
        );
        let left = compose_image_states(
            &pairs[0],
            &pairs[1],
            b,
            &mut field,
            &mut quadratic_solves,
            &mut roots_returned,
        );
        let right = compose_image_states(
            &pairs[2],
            &pairs[3],
            b,
            &mut field,
            &mut quadratic_solves,
            &mut roots_returned,
        );
        let left_reference = reference_image(&xs[..4], b, square_roots);
        let right_reference = reference_image(&xs[4..], b, square_roots);
        if left != left_reference || right != right_reference {
            return Err(format!(
                "four-leaf image mismatch in {kind} at x atom {xs:?}"
            ));
        }
        maximum_four_leaf_affine_image = maximum_four_leaf_affine_image
            .max(left.affine.len())
            .max(right.affine.len());
        four_leaf_identity_images += u64::from(left.identity) + u64::from(right.identity);
        let mut atom_positive = false;
        for &left_x in &left.affine {
            quadratic_solves += 1;
            let (qa, qb, qc) = coefficients(&mut field, left_x, target.x, b);
            let roots = solve_quadratic(&mut field, qa, qb, qc);
            roots_returned += roots.len as u64;
            for &right_x in &roots.roots[..roots.len] {
                final_index_lookups += 1;
                atom_positive |= right.affine.contains(&right_x);
            }
        }
        atom_positive |= left.identity && right.affine.contains(&target.x);
        atom_positive |= right.identity && left.affine.contains(&target.x);
        let reference = reference_image(xs, b, square_roots);
        let reference_positive = reference.affine.contains(&target.x);
        if atom_positive != reference_positive {
            return Err(format!(
                "final eight-leaf classification mismatch in {kind} at x atom {xs:?}"
            ));
        }
        candidate_positive |= atom_positive;
        writeln!(
            stream,
            "{:?}|{:?}:{}|{:?}:{}|{}|{}",
            xs,
            left.affine,
            u8::from(left.identity),
            right.affine,
            u8::from(right.identity),
            u8::from(atom_positive),
            u8::from(reference_positive)
        )
        .expect("writing to String cannot fail");
    }
    Ok(S9CellResult {
        kind: kind.into(),
        boundaries,
        expected_positive,
        candidate_positive,
        x_atoms: atoms.len(),
        signed_tuples,
        quadratic_solves,
        roots_returned,
        final_index_lookups,
        field_multiplications: field.multiplications,
        maximum_two_leaf_affine_image,
        maximum_four_leaf_affine_image,
        four_leaf_identity_images,
        false_negative: expected_positive && !candidate_positive,
        false_positive: !expected_positive && candidate_positive,
        image_stream_sha256: hex::encode(sha256(stream.as_bytes())),
    })
}

fn frozen_nonempty_inventory(terminal: u64) -> Vec<(u64, Vec<u64>)> {
    match terminal {
        0 => vec![
            (25, vec![519]),
            (41, vec![1061]),
            (262, vec![372]),
            (342, vec![583]),
            (384, vec![409]),
            (433, vec![678]),
            (718, vec![408, 743]),
            (767, vec![725]),
            (809, vec![1040]),
            (894, vec![969]),
            (1110, vec![320]),
            (1126, vec![235]),
        ],
        369 => vec![
            (94, vec![163, 988]),
            (166, vec![390]),
            (288, vec![574, 577]),
            (507, vec![462]),
            (539, vec![328]),
            (614, vec![935]),
            (626, vec![138]),
            (644, vec![229]),
            (704, vec![330]),
            (863, vec![273, 878]),
            (985, vec![318]),
            (1057, vec![470, 681]),
        ],
        _ => Vec::new(),
    }
}

fn run_s9_terminal(
    terminal: u64,
    b: u64,
    square_roots: &[Vec<u64>],
) -> Result<S9TerminalResult, String> {
    let mut inventory = Vec::new();
    let mut x_by_boundary = BTreeMap::new();
    let mut points_by_boundary = BTreeMap::new();
    for boundary in roots(DEPTH - 1, terminal) {
        let xs = liftable_x(boundary, b, square_roots);
        if xs.is_empty() {
            continue;
        }
        let points = signed_points(&xs, b, square_roots);
        inventory.push((boundary, xs.clone()));
        x_by_boundary.insert(boundary, xs);
        points_by_boundary.insert(boundary, points);
    }
    if inventory != frozen_nonempty_inventory(terminal) {
        return Err(format!("terminal {terminal} factor-base inventory changed"));
    }
    let selected_boundaries: [u64; 8] = inventory[..8]
        .iter()
        .map(|entry| entry.0)
        .collect::<Vec<_>>()
        .try_into()
        .expect("eight selected boundaries");
    let selected_point_sets: Vec<Vec<Affine>> = selected_boundaries
        .iter()
        .map(|boundary| points_by_boundary[boundary].clone())
        .collect();
    let target = first_finite_sum(&selected_point_sets)
        .ok_or_else(|| format!("terminal {terminal} planted target stayed infinity"))?;
    let reverse: [u64; 8] = selected_boundaries
        .iter()
        .rev()
        .copied()
        .collect::<Vec<_>>()
        .try_into()
        .expect("eight reversed boundaries");
    let mut selected_cells = Vec::new();
    for (kind, boundaries) in [
        ("planted-forward", selected_boundaries),
        ("planted-reverse", reverse),
    ] {
        let (positive, signed_tuples) =
            boundary_cell_reference(boundaries, &points_by_boundary, target);
        if !positive {
            return Err(format!("terminal {terminal} {kind} cell is not positive"));
        }
        selected_cells.push((kind.to_string(), boundaries, positive, signed_tuples));
    }
    let pool = [inventory[0].0, inventory[1].0];
    for code in 0u16..256 {
        let mut boundaries = [0u64; 8];
        for (index, boundary) in boundaries.iter_mut().enumerate() {
            *boundary = pool[((code >> (7 - index)) & 1) as usize];
        }
        let (positive, signed_tuples) =
            boundary_cell_reference(boundaries, &points_by_boundary, target);
        if !positive {
            let negative_index = selected_cells.len() - 2;
            selected_cells.push((
                format!("negative-{negative_index}"),
                boundaries,
                false,
                signed_tuples,
            ));
            if selected_cells.len() == 4 {
                break;
            }
        }
    }
    if selected_cells.len() != 4 {
        return Err(format!(
            "terminal {terminal} did not yield two negative controls"
        ));
    }
    let mut cells = Vec::new();
    for (kind, boundaries, expected_positive, signed_tuples) in selected_cells {
        cells.push(run_s9_cell(
            &kind,
            boundaries,
            expected_positive,
            signed_tuples,
            target,
            &x_by_boundary,
            b,
            square_roots,
        )?);
    }
    Ok(S9TerminalResult {
        terminal,
        target: [target.x, target.y],
        nonempty_boundaries: inventory.len(),
        selected_boundaries,
        cells,
    })
}

fn run_s9(round11_path: &PathBuf, out: Option<PathBuf>) -> Result<(), String> {
    let bytes = std::fs::read(round11_path).map_err(|error| error.to_string())?;
    let round11_sha256 = hex::encode(sha256(&bytes));
    if round11_sha256 != ROUND11_SHA256 {
        return Err(format!(
            "round-11 SHA-256 mismatch: expected {ROUND11_SHA256}, got {round11_sha256}"
        ));
    }
    let p256 = CurveParams::p256();
    let b = (&p256.b % P).to_u64().ok_or("P-256 b reduction failed")?;
    let square_roots = square_root_table();
    let mut terminal_results = Vec::new();
    for terminal in TERMINALS {
        terminal_results.push(run_s9_terminal(terminal, b, &square_roots)?);
    }
    let result = S9ExperimentResult {
        schema: "p256.dickson_s9_image_growth/v1".into(),
        curve: CURVE_SLUG.into(),
        prime: P,
        a: A,
        b,
        round11_sha256,
        local_maximum_degree: 2,
        terminal_results,
    };
    let text = serde_json::to_string_pretty(&result).map_err(|error| error.to_string())? + "\n";
    match out {
        Some(path) => std::fs::write(path, text).map_err(|error| error.to_string())?,
        None => print!("{text}"),
    }
    Ok(())
}

fn run_round9(out: Option<PathBuf>) -> Result<(), String> {
    let p256 = CurveParams::p256();
    let b = (&p256.b % P).to_u64().ok_or("P-256 b reduction failed")?;
    let square_roots = square_root_table();
    let mut terminal_results = Vec::new();
    for terminal in TERMINALS {
        terminal_results.push(run_terminal(terminal, b, &square_roots)?);
    }
    let result = ExperimentResult {
        schema: "p256.dickson_balanced_s5_degree/v1".into(),
        curve: CURVE_SLUG.into(),
        prime: P,
        a: A,
        b,
        depth: DEPTH,
        terminals: TERMINALS.to_vec(),
        cells_per_class: CELLS_PER_CLASS,
        max_degree: MAX_DEGREE,
        budget_seconds_per_arm_cell: BUDGET_SECS,
        terminal_results,
    };
    let text = serde_json::to_string_pretty(&result).map_err(|error| error.to_string())? + "\n";
    match out {
        Some(path) => std::fs::write(path, text).map_err(|error| error.to_string())?,
        None => print!("{text}"),
    }
    Ok(())
}

fn run(cli: Cli) -> Result<(), String> {
    match (
        cli.atomize_round9,
        cli.image_atomize_round10,
        cli.scale_eight_round11,
    ) {
        (Some(_), Some(_), _) | (Some(_), _, Some(_)) | (_, Some(_), Some(_)) => {
            Err("choose only one continuation mode".into())
        }
        (Some(path), None, None) => run_atomized(&path, cli.out),
        (None, Some(path), None) => run_image_atomized(&path, cli.out),
        (None, None, Some(path)) => run_s9(&path, cli.out),
        (None, None, None) => run_round9(cli.out),
    }
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

    fn direct_s3(x: u64, y: u64, z: u64, b: u64) -> u64 {
        let q = mulm(x, y);
        let sum = addm(x, y);
        let difference = subm(x, y);
        let first = mulm(mulm(z, z), mulm(difference, difference));
        let inner = addm(mulm(sum, addm(q, A)), mulm(2, b));
        let second = negm(mulm(mulm(2, z), inner));
        let third = subm(mulm(subm(q, A), subm(q, A)), mulm(mulm(4, b), sum));
        addm(addm(first, second), third)
    }

    #[test]
    fn variable_s3_quadraticisation_equals_direct_polynomial() {
        let n = 7;
        let x = Pol::var(n, 0);
        let y = Pol::var(n, 1);
        let z = Pol::var(n, 2);
        let s = Pol::var(n, 3);
        let q = Pol::var(n, 4);
        let h = Pol::var(n, 5);
        let k = Pol::var(n, 6);
        let b = 37;
        let mut equations = Vec::new();
        variable_s3(&mut equations, &x, &y, &z, &s, &q, &h, &k, b);
        assert_eq!(equations.len(), 5);
        assert_eq!(equations.iter().map(Pol::degree).max(), Some(2));
        for (xv, yv, zv) in [(1, 2, 3), (17, 301, 500), (1149, 42, 91)] {
            let values = [
                xv,
                yv,
                zv,
                addm(xv, yv),
                mulm(xv, yv),
                mulm(zv, zv),
                mulm(zv, addm(xv, yv)),
            ];
            for equation in &equations[..4] {
                assert_eq!(equation.evaluate(&values), 0);
            }
            assert_eq!(equations[4].evaluate(&values), direct_s3(xv, yv, zv, b));
        }
    }

    #[test]
    fn specialised_leaf_polynomial_has_exact_roots() {
        for roots in [&[17u64][..], &[17u64, 301][..]] {
            let polynomial = leaf_polynomial(1, 0, roots).unwrap();
            let actual: Vec<u64> = (0..P)
                .filter(|&value| polynomial.evaluate(&[value]) == 0)
                .collect();
            assert_eq!(actual, roots);
            assert!(polynomial.degree() <= 2);
        }
    }

    #[test]
    fn mismatching_positive_cells_require_an_identity_chart() {
        let p256 = CurveParams::p256();
        let b = (&p256.b % P).to_u64().unwrap();
        let square_roots = square_root_table();
        let cases = [
            (0, vec![[25, 25, 767, 1110], [25, 25, 1110, 767]]),
            (369, vec![[94, 94, 288, 1057]]),
        ];
        for (terminal, mismatches) in cases {
            let mut point_sets = BTreeMap::new();
            for boundary in roots(DEPTH - 1, terminal) {
                let xs = liftable_x(boundary, b, &square_roots);
                let points = signed_points(&xs, b, &square_roots);
                if !points.is_empty() {
                    point_sets.insert(boundary, points);
                }
            }
            let nonempty: Vec<(u64, Vec<Affine>)> = point_sets
                .iter()
                .map(|(&boundary, points)| (boundary, points.clone()))
                .collect();
            let (_, target) = planted_target(&nonempty).unwrap();
            for boundaries in mismatches {
                let sets = [
                    &point_sets[&boundaries[0]],
                    &point_sets[&boundaries[1]],
                    &point_sets[&boundaries[2]],
                    &point_sets[&boundaries[3]],
                ];
                let mut witnesses = 0;
                for &p1 in sets[0] {
                    for &p2 in sets[1] {
                        for &p3 in sets[2] {
                            for &p4 in sets[3] {
                                let sum = add_four([p1, p2, p3, p4]);
                                if sum != Point::Affine(target)
                                    && sum != Point::Affine(negate(target))
                                {
                                    continue;
                                }
                                witnesses += 1;
                                let left = add_points(Point::Affine(p1), Point::Affine(p2));
                                let right = add_points(Point::Affine(p3), Point::Affine(p4));
                                assert!(
                                    left == Point::Infinity || right == Point::Infinity,
                                    "cell {boundaries:?} had an affine-intermediate witness"
                                );
                            }
                        }
                    }
                }
                assert!(
                    witnesses > 0,
                    "cell {boundaries:?} had no group-law witness"
                );
            }
        }
    }
}

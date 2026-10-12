//! Direct quadraticization of Semaev S3 on frozen Dickson fibres.

use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;
use std::process::ExitCode;
use std::time::Duration;

use clap::Parser;
use crypto_lib::cryptanalysis::f4_fp::{self, F4Options, Ordering, Poly};
use crypto_lib::ecc::curve::CurveParams;
use num_traits::ToPrimitive;
use serde::Serialize;

const P: u64 = 1151;
const MAX_DEGREE: u32 = 12;
const BUDGET_SECS: u64 = 30;
const TARGETS: usize = 3;

#[derive(Parser)]
#[command(about = "Measure direct quadraticized S3 on frozen P-256 Dickson analogues")]
struct Cli {
    /// Compact raw JSON output; stdout when omitted.
    #[arg(long)]
    out: Option<PathBuf>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize)]
#[serde(rename_all = "kebab-case")]
enum Layout {
    Grouped,
    AuxFirst,
    LevelInterleaved,
    ReverseBlocks,
}

const LAYOUTS: [Layout; 4] = [
    Layout::Grouped,
    Layout::AuxFirst,
    Layout::LevelInterleaved,
    Layout::ReverseBlocks,
];

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

    fn degree(&self) -> u32 {
        self.terms
            .keys()
            .map(|exponent| exponent.iter().sum())
            .max()
            .unwrap_or(0)
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

#[derive(Clone, Debug)]
struct Indices {
    chain1: Vec<usize>,
    y1: usize,
    u1: usize,
    chain2: Vec<usize>,
    y2: usize,
    u2: usize,
    q: usize,
    n: usize,
}

impl Indices {
    fn new(depth: u32, layout: Layout) -> Self {
        let t = depth as usize;
        match layout {
            Layout::Grouped => Self {
                chain1: (0..t).collect(),
                y1: t,
                u1: t + 1,
                chain2: (t + 2..2 * t + 2).collect(),
                y2: 2 * t + 2,
                u2: 2 * t + 3,
                q: 2 * t + 4,
                n: 2 * t + 5,
            },
            Layout::AuxFirst => Self {
                q: 0,
                u1: 1,
                u2: 2,
                y1: 3,
                y2: 4,
                chain1: (5..5 + t).collect(),
                chain2: (5 + t..5 + 2 * t).collect(),
                n: 2 * t + 5,
            },
            Layout::LevelInterleaved => Self {
                chain1: (0..t).map(|level| 2 * level).collect(),
                chain2: (0..t).map(|level| 2 * level + 1).collect(),
                y1: 2 * t,
                y2: 2 * t + 1,
                u1: 2 * t + 2,
                u2: 2 * t + 3,
                q: 2 * t + 4,
                n: 2 * t + 5,
            },
            Layout::ReverseBlocks => Self {
                q: 0,
                u2: 1,
                y2: 2,
                chain2: (0..t).map(|level| 3 + t - 1 - level).collect(),
                u1: 3 + t,
                y1: 4 + t,
                chain1: (0..t).map(|level| 5 + 2 * t - 1 - level).collect(),
                n: 2 * t + 5,
            },
        }
    }
}

#[derive(Clone, Debug, Serialize)]
struct Metric {
    target: [u64; 2],
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
    milliseconds: f64,
}

#[derive(Clone, Debug, Serialize)]
struct LayoutResult {
    layout: Layout,
    runs: Vec<Metric>,
}

#[derive(Clone, Debug, Serialize)]
struct Cell {
    depth: u32,
    terminal: u64,
    liftable_roots: usize,
    common_targets_available: usize,
    selected_targets: Vec<[u64; 2]>,
    quartic_s3: Vec<Metric>,
    quadratic_s3: Vec<LayoutResult>,
}

#[derive(Clone, Debug, Serialize)]
struct ExperimentResult {
    schema: String,
    curve: String,
    prime: u64,
    a: u64,
    b: u64,
    max_degree: u32,
    budget_seconds: u64,
    cells: Vec<Cell>,
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

fn square_roots(value: u64) -> Vec<u64> {
    (0..P)
        .filter(|candidate| mulm(*candidate, *candidate) == value)
        .collect()
}

fn curve_rhs(x: u64, a: u64, b: u64) -> u64 {
    addm(addm(mulm(mulm(x, x), x), mulm(a, x)), b)
}

fn curve_points(a: u64, b: u64) -> Vec<Affine> {
    let mut points = Vec::new();
    for x in 0..P {
        for y in square_roots(curve_rhs(x, a, b)) {
            points.push(Affine { x, y });
        }
    }
    points
}

fn add_points(left: Affine, right: Affine, a: u64) -> Point {
    let slope = if left.x == right.x {
        if addm(left.y, right.y) == 0 || left.y == 0 {
            return Point::Infinity;
        }
        mulm(
            addm(mulm(3, mulm(left.x, left.x)), a),
            inverse(mulm(2, left.y)),
        )
    } else {
        mulm(subm(right.y, left.y), inverse(subm(right.x, left.x)))
    };
    let x = subm(subm(mulm(slope, slope), left.x), right.x);
    let y = subm(mulm(slope, subm(left.x, x)), left.y);
    Point::Affine(Affine { x, y })
}

fn dickson_roots(depth: u32, terminal: u64) -> Vec<u64> {
    (0..P)
        .filter(|&x| {
            let mut value = x;
            for _ in 0..depth {
                value = subm(mulm(value, value), 2);
            }
            value == terminal
        })
        .collect()
}

fn signed_points(xs: &[u64], a: u64, b: u64) -> Vec<Affine> {
    xs.iter()
        .flat_map(|&x| {
            square_roots(curve_rhs(x, a, b))
                .into_iter()
                .map(move |y| Affine { x, y })
        })
        .collect()
}

fn decompositions(points: &[Affine], a: u64) -> BTreeSet<Affine> {
    let mut out = BTreeSet::new();
    for &left in points {
        for &right in points {
            if let Point::Affine(sum) = add_points(left, right, a) {
                out.insert(sum);
            }
        }
    }
    out
}

fn add_chain(equations: &mut Vec<Pol>, n: usize, chain: &[usize], terminal: u64) -> Pol {
    for (level, &index) in chain.iter().enumerate() {
        let current = Pol::var(n, index);
        let mut equation = current.square().sub(&Pol::constant(n, 2));
        if level + 1 < chain.len() {
            equation = equation.sub(&Pol::var(n, chain[level + 1]));
        } else {
            equation = equation.sub(&Pol::constant(n, terminal));
        }
        equations.push(equation);
    }
    Pol::var(n, chain[0])
}

fn s3(x1: &Pol, x2: &Pol, x3: &Pol, a: u64, b: u64) -> Pol {
    let n = x1.n;
    let difference = x1.sub(x2);
    let first = difference.square().mul(&x3.square());
    let sum = x1.add(x2);
    let product = x1.mul(x2);
    let inner = sum
        .mul(&product.add(&Pol::constant(n, a)))
        .add(&Pol::constant(n, mulm(2, b)));
    let second = inner.mul(x3).scale(P - 2);
    let shifted = product.sub(&Pol::constant(n, a));
    let third = shifted.square().sub(&sum.scale(mulm(4, b)));
    first.add(&second).add(&third)
}

fn quartic_system(depth: u32, terminal: u64, target_x: u64, a: u64, b: u64) -> (Vec<Pol>, usize) {
    let t = depth as usize;
    let block = t + 1;
    let n = 2 * block;
    let mut equations = Vec::new();
    let chain1: Vec<usize> = (0..t).collect();
    let chain2: Vec<usize> = (block..block + t).collect();
    let x1 = add_chain(&mut equations, n, &chain1, terminal);
    let x2 = add_chain(&mut equations, n, &chain2, terminal);
    for (y_index, x) in [(t, &x1), (block + t, &x2)] {
        let y = Pol::var(n, y_index);
        let rhs = x.square().mul(x).add(&x.scale(a)).add(&Pol::constant(n, b));
        equations.push(y.square().sub(&rhs));
    }
    equations.push(s3(&x1, &x2, &Pol::constant(n, target_x), a, b));
    (equations, n)
}

fn quadratic_system(
    depth: u32,
    terminal: u64,
    target_x: u64,
    a: u64,
    b: u64,
    layout: Layout,
) -> (Vec<Pol>, usize) {
    let indices = Indices::new(depth, layout);
    let n = indices.n;
    let mut equations = Vec::new();
    let x1 = add_chain(&mut equations, n, &indices.chain1, terminal);
    let x2 = add_chain(&mut equations, n, &indices.chain2, terminal);
    for (x, y_index, u_index) in [(&x1, indices.y1, indices.u1), (&x2, indices.y2, indices.u2)] {
        let y = Pol::var(n, y_index);
        let u = Pol::var(n, u_index);
        equations.push(u.sub(&x.square()));
        equations.push(
            y.square()
                .sub(&u.mul(x))
                .sub(&x.scale(a))
                .sub(&Pol::constant(n, b)),
        );
    }
    let q = Pol::var(n, indices.q);
    equations.push(q.sub(&x1.mul(&x2)));
    let sum = x1.add(&x2);
    let difference = x1.sub(&x2);
    let first = difference.square().scale(mulm(target_x, target_x));
    let inner = sum
        .mul(&q.add(&Pol::constant(n, a)))
        .add(&Pol::constant(n, mulm(2, b)));
    let second = inner.scale(subm(0, mulm(2, target_x)));
    let third = q
        .sub(&Pol::constant(n, a))
        .square()
        .sub(&sum.scale(mulm(4, b)));
    equations.push(first.add(&second).add(&third));
    (equations, n)
}

fn measure(equations: Vec<Pol>, n: usize, target: Affine) -> Metric {
    let input_max_degree = equations.iter().map(Pol::degree).max().unwrap_or(0);
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
    Metric {
        target: [target.x, target.y],
        variables: n,
        equations: equations.len(),
        input_max_degree,
        inconsistent: report.inconsistent,
        complete,
        correct: complete && !report.inconsistent,
        timed_out: report.timed_out,
        pairs_above_bound: report.pairs_above_bound,
        solving_degree_max: report.solving_degree_max,
        max_cols_to_solution: report.max_cols_to_solution,
        field_ops: report.field_ops,
        milliseconds: report.ms,
    }
}

fn run_cell(depth: u32, terminal: u64, a: u64, b: u64, all_points: &[Affine]) -> Cell {
    let zero_x: Vec<u64> = dickson_roots(depth, 0)
        .into_iter()
        .filter(|&x| !square_roots(curve_rhs(x, a, b)).is_empty())
        .collect();
    let candidate_x: Vec<u64> = dickson_roots(depth, terminal)
        .into_iter()
        .filter(|&x| !square_roots(curve_rhs(x, a, b)).is_empty())
        .collect();
    let zero_decompositions = decompositions(&signed_points(&zero_x, a, b), a);
    let candidate_decompositions = decompositions(&signed_points(&candidate_x, a, b), a);
    let common: Vec<Affine> = all_points
        .iter()
        .copied()
        .filter(|target| {
            zero_decompositions.contains(target) && candidate_decompositions.contains(target)
        })
        .collect();
    let targets: Vec<Affine> = common.iter().copied().take(TARGETS).collect();
    let quartic_s3 = targets
        .iter()
        .copied()
        .map(|target| {
            let (equations, n) = quartic_system(depth, terminal, target.x, a, b);
            measure(equations, n, target)
        })
        .collect();
    let quadratic_s3 = LAYOUTS
        .into_iter()
        .map(|layout| LayoutResult {
            layout,
            runs: targets
                .iter()
                .copied()
                .map(|target| {
                    let (equations, n) = quadratic_system(depth, terminal, target.x, a, b, layout);
                    measure(equations, n, target)
                })
                .collect(),
        })
        .collect();
    Cell {
        depth,
        terminal,
        liftable_roots: candidate_x.len(),
        common_targets_available: common.len(),
        selected_targets: targets.iter().map(|target| [target.x, target.y]).collect(),
        quartic_s3,
        quadratic_s3,
    }
}

fn run(cli: Cli) -> Result<(), String> {
    let p256 = CurveParams::p256();
    let b = (&p256.b % P).to_u64().ok_or("P-256 b reduction failed")?;
    let a = P - 3;
    let all_points = curve_points(a, b);
    let cells = [(4, 0), (4, 782), (5, 0), (5, 369)]
        .into_iter()
        .map(|(depth, terminal)| run_cell(depth, terminal, a, b, &all_points))
        .collect::<Vec<_>>();
    if cells
        .iter()
        .flat_map(|cell| {
            cell.quartic_s3.iter().chain(
                cell.quadratic_s3
                    .iter()
                    .flat_map(|layout| layout.runs.iter()),
            )
        })
        .any(|metric| metric.complete && !metric.correct)
    {
        return Err("a certified F4 verdict contradicted an exhaustive positive target".into());
    }
    let result = ExperimentResult {
        schema: "p256.dickson_quadratic_s3/v1".into(),
        curve: crypto_lib::cryptanalysis::p256_dickson_factor_base::CURVE_SLUG.into(),
        prime: P,
        a,
        b,
        max_degree: MAX_DEGREE,
        budget_seconds: BUDGET_SECS,
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
    fn layouts_are_permutations() {
        for depth in [4, 5] {
            for layout in LAYOUTS {
                let indices = Indices::new(depth, layout);
                let mut values = indices.chain1.clone();
                values.extend([indices.y1, indices.u1, indices.y2, indices.u2, indices.q]);
                values.extend(indices.chain2);
                values.sort_unstable();
                assert_eq!(values, (0..indices.n).collect::<Vec<_>>());
            }
        }
    }

    #[test]
    fn direct_s3_quadraticization_has_degree_two() {
        let p256 = CurveParams::p256();
        let b = (&p256.b % P).to_u64().unwrap();
        let target = curve_points(P - 3, b)[0];
        for layout in LAYOUTS {
            let (system, _) = quadratic_system(4, 782, target.x, P - 3, b, layout);
            assert_eq!(system.iter().map(Pol::degree).max(), Some(2));
        }
    }
}

//! Component-wise quadratic S3 on branched Dickson fibres.

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
const DEPTH: u32 = 5;
const MAX_DEGREE: u32 = 12;
const BUDGET_SECS: u64 = 30;
const TARGETS: usize = 3;

#[derive(Parser)]
#[command(about = "Measure branched quadratic S3 on frozen P-256 Dickson analogues")]
struct Cli {
    /// Compact raw JSON output; stdout when omitted.
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

#[derive(Clone, Debug, Serialize)]
struct ComponentRun {
    target: [u64; 2],
    boundary1: u64,
    boundary2: u64,
    expected_positive: bool,
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
struct BranchResult {
    branch_depth: u32,
    boundaries: Vec<u64>,
    components_per_target: usize,
    runs: Vec<ComponentRun>,
}

#[derive(Clone, Debug, Serialize)]
struct TerminalResult {
    terminal: u64,
    liftable_roots: usize,
    common_targets_available: usize,
    selected_targets: Vec<[u64; 2]>,
    branches: Vec<BranchResult>,
}

#[derive(Clone, Debug, Serialize)]
struct ExperimentResult {
    schema: String,
    curve: String,
    prime: u64,
    a: u64,
    b: u64,
    depth: u32,
    max_degree: u32,
    budget_seconds_per_component: u64,
    terminals: Vec<TerminalResult>,
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

fn signed_points(xs: &[u64], a: u64, b: u64) -> Vec<Affine> {
    xs.iter()
        .flat_map(|&x| {
            square_roots(curve_rhs(x, a, b))
                .into_iter()
                .map(move |y| Affine { x, y })
        })
        .collect()
}

fn sums(points1: &[Affine], points2: &[Affine], a: u64) -> BTreeSet<Affine> {
    let mut out = BTreeSet::new();
    for &left in points1 {
        for &right in points2 {
            if let Point::Affine(sum) = add_points(left, right, a) {
                out.insert(sum);
            }
        }
    }
    out
}

fn add_chain(equations: &mut Vec<Pol>, n: usize, start: usize, depth: u32, terminal: u64) -> Pol {
    for level in 0..depth as usize {
        let current = Pol::var(n, start + level);
        let mut equation = current.square().sub(&Pol::constant(n, 2));
        if level + 1 < depth as usize {
            equation = equation.sub(&Pol::var(n, start + level + 1));
        } else {
            equation = equation.sub(&Pol::constant(n, terminal));
        }
        equations.push(equation);
    }
    Pol::var(n, start)
}

fn system(
    depth: u32,
    boundary1: u64,
    boundary2: u64,
    target_x: u64,
    a: u64,
    b: u64,
) -> (Vec<Pol>, usize) {
    let t = depth as usize;
    let block = t + 2;
    let n = 2 * block + 1;
    let mut equations = Vec::new();
    let x1 = add_chain(&mut equations, n, 0, depth, boundary1);
    let x2 = add_chain(&mut equations, n, block, depth, boundary2);
    let y1 = Pol::var(n, t);
    let u1 = Pol::var(n, t + 1);
    let y2 = Pol::var(n, block + t);
    let u2 = Pol::var(n, block + t + 1);
    for (x, y, u) in [(&x1, &y1, &u1), (&x2, &y2, &u2)] {
        equations.push(u.sub(&x.square()));
        equations.push(
            y.square()
                .sub(&u.mul(x))
                .sub(&x.scale(a))
                .sub(&Pol::constant(n, b)),
        );
    }
    let q = Pol::var(n, 2 * block);
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

fn measure(
    depth: u32,
    boundary1: u64,
    boundary2: u64,
    target: Affine,
    expected_positive: bool,
    a: u64,
    b: u64,
) -> ComponentRun {
    let (equations, n) = system(depth, boundary1, boundary2, target.x, a, b);
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
    ComponentRun {
        target: [target.x, target.y],
        boundary1,
        boundary2,
        expected_positive,
        variables: n,
        equations: equations.len(),
        input_max_degree,
        inconsistent: report.inconsistent,
        complete,
        correct: complete && ((!report.inconsistent) == expected_positive),
        timed_out: report.timed_out,
        pairs_above_bound: report.pairs_above_bound,
        solving_degree_max: report.solving_degree_max,
        max_cols_to_solution: report.max_cols_to_solution,
        field_ops: report.field_ops,
        milliseconds: report.ms,
    }
}

fn target_set(terminal: u64, a: u64, b: u64, all_points: &[Affine]) -> (Vec<Affine>, usize) {
    let zero_x: Vec<u64> = roots(DEPTH, 0)
        .into_iter()
        .filter(|&x| !square_roots(curve_rhs(x, a, b)).is_empty())
        .collect();
    let candidate_x: Vec<u64> = roots(DEPTH, terminal)
        .into_iter()
        .filter(|&x| !square_roots(curve_rhs(x, a, b)).is_empty())
        .collect();
    let zero_points = signed_points(&zero_x, a, b);
    let candidate_points = signed_points(&candidate_x, a, b);
    let zero_sums = sums(&zero_points, &zero_points, a);
    let candidate_sums = sums(&candidate_points, &candidate_points, a);
    let common: Vec<Affine> = all_points
        .iter()
        .copied()
        .filter(|point| zero_sums.contains(point) && candidate_sums.contains(point))
        .collect();
    (common.iter().copied().take(TARGETS).collect(), common.len())
}

fn run_terminal(terminal: u64, a: u64, b: u64, all_points: &[Affine]) -> TerminalResult {
    let full_x: Vec<u64> = roots(DEPTH, terminal)
        .into_iter()
        .filter(|&x| !square_roots(curve_rhs(x, a, b)).is_empty())
        .collect();
    let (targets, available) = target_set(terminal, a, b, all_points);
    let branches = [0, 1, 2]
        .into_iter()
        .map(|branch_depth| {
            let boundaries = roots(branch_depth, terminal);
            let remaining = DEPTH - branch_depth;
            let point_sets: BTreeMap<u64, Vec<Affine>> = boundaries
                .iter()
                .map(|&boundary| {
                    let xs: Vec<u64> = roots(remaining, boundary)
                        .into_iter()
                        .filter(|&x| !square_roots(curve_rhs(x, a, b)).is_empty())
                        .collect();
                    (boundary, signed_points(&xs, a, b))
                })
                .collect();
            let mut runs = Vec::new();
            for &target in &targets {
                for &boundary1 in &boundaries {
                    for &boundary2 in &boundaries {
                        let expected_positive =
                            sums(&point_sets[&boundary1], &point_sets[&boundary2], a)
                                .contains(&target);
                        runs.push(measure(
                            remaining,
                            boundary1,
                            boundary2,
                            target,
                            expected_positive,
                            a,
                            b,
                        ));
                    }
                }
            }
            BranchResult {
                branch_depth,
                components_per_target: boundaries.len() * boundaries.len(),
                boundaries,
                runs,
            }
        })
        .collect();
    TerminalResult {
        terminal,
        liftable_roots: full_x.len(),
        common_targets_available: available,
        selected_targets: targets.iter().map(|target| [target.x, target.y]).collect(),
        branches,
    }
}

fn run(cli: Cli) -> Result<(), String> {
    let p256 = CurveParams::p256();
    let b = (&p256.b % P).to_u64().ok_or("P-256 b reduction failed")?;
    let a = P - 3;
    let all_points = curve_points(a, b);
    let result = ExperimentResult {
        schema: "p256.dickson_branched_s3/v1".into(),
        curve: crypto_lib::cryptanalysis::p256_dickson_factor_base::CURVE_SLUG.into(),
        prime: P,
        a,
        b,
        depth: DEPTH,
        max_degree: MAX_DEGREE,
        budget_seconds_per_component: BUDGET_SECS,
        terminals: [0, 369]
            .into_iter()
            .map(|terminal| run_terminal(terminal, a, b, &all_points))
            .collect(),
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
    fn frozen_branches_have_power_of_two_boundaries() {
        for terminal in [0, 369] {
            assert_eq!(roots(0, terminal), vec![terminal]);
            assert_eq!(roots(1, terminal).len(), 2);
            assert_eq!(roots(2, terminal).len(), 4);
        }
    }

    #[test]
    fn branched_system_remains_quadratic() {
        let p256 = CurveParams::p256();
        let b = (&p256.b % P).to_u64().unwrap();
        let (equations, _) = system(3, 1, 2, 3, P - 3, b);
        assert_eq!(equations.iter().map(Pol::degree).max(), Some(2));
    }
}

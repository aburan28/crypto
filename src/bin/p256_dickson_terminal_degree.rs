//! Sweep full Dickson terminal fibres on the frozen P-256 small-prime analogue.

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
const TARGETS_PER_TERMINAL: usize = 3;

#[derive(Parser)]
#[command(about = "Sweep Dickson terminal fibres for paired S3 solving degree")]
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
struct SolveMetric {
    target: [u64; 2],
    variables: usize,
    equations: usize,
    input_max_degree: u32,
    inconsistent: bool,
    complete: bool,
    correct: bool,
    timed_out: bool,
    pairs_above_bound: usize,
    degree_reached: u32,
    solving_degree_max: u32,
    max_cols_to_solution: usize,
    field_ops: u64,
    milliseconds: f64,
}

#[derive(Clone, Debug, Serialize)]
struct PairedRun {
    candidate: SolveMetric,
    terminal_zero: SolveMetric,
}

#[derive(Clone, Debug, Serialize)]
struct TerminalResult {
    terminal: u64,
    full_roots: usize,
    liftable_roots: usize,
    common_targets_available: usize,
    selected_targets: Vec<[u64; 2]>,
    paired_runs: Vec<PairedRun>,
}

#[derive(Clone, Debug, Serialize)]
struct DepthResult {
    depth: u32,
    terminal_zero_liftable_roots: usize,
    eligible_terminals: usize,
    winner_terminal: Option<u64>,
    terminals: Vec<TerminalResult>,
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
    targets_per_terminal: usize,
    depths: Vec<DepthResult>,
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
    base %= P;
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
    assert_ne!(value % P, 0);
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

fn signed_points(abscissae: &[u64], a: u64, b: u64) -> Vec<Affine> {
    abscissae
        .iter()
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

fn dickson_fibres(depth: u32) -> Vec<Vec<u64>> {
    let mut fibres = vec![Vec::new(); P as usize];
    for x in 0..P {
        let mut terminal = x;
        for _ in 0..depth {
            terminal = subm(mulm(terminal, terminal), 2);
        }
        fibres[terminal as usize].push(x);
    }
    fibres
}

fn add_domain(equations: &mut Vec<Pol>, n: usize, base: usize, depth: u32, terminal: u64) -> Pol {
    for index in 0..depth as usize {
        let current = Pol::var(n, base + index);
        let mut equation = current.square().sub(&Pol::constant(n, 2));
        if index + 1 < depth as usize {
            equation = equation.sub(&Pol::var(n, base + index + 1));
        } else {
            equation = equation.sub(&Pol::constant(n, terminal));
        }
        equations.push(equation);
    }
    Pol::var(n, base)
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

fn build_system(depth: u32, terminal: u64, target_x: u64, a: u64, b: u64) -> (Vec<Pol>, usize) {
    let block = depth as usize + 1;
    let n = 2 * block;
    let mut equations = Vec::new();
    let x1 = add_domain(&mut equations, n, 0, depth, terminal);
    let x2 = add_domain(&mut equations, n, block, depth, terminal);
    for (base, x) in [(0, &x1), (block, &x2)] {
        let y = Pol::var(n, base + depth as usize);
        let rhs = x.square().mul(x).add(&x.scale(a)).add(&Pol::constant(n, b));
        equations.push(y.square().sub(&rhs));
    }
    equations.push(s3(&x1, &x2, &Pol::constant(n, target_x), a, b));
    (equations, n)
}

fn solve(depth: u32, terminal: u64, target: Affine, a: u64, b: u64) -> SolveMetric {
    let (equations, n) = build_system(depth, terminal, target.x, a, b);
    let input_max_degree = equations.iter().map(Pol::degree).max().unwrap_or(0);
    let input: Vec<Poly> = equations.iter().map(Pol::to_f4).collect();
    let options =
        F4Options::new(Ordering::Grevlex, MAX_DEGREE).with_budget(Duration::from_secs(BUDGET_SECS));
    let report = f4_fp::f4(&input, n, P, &options);
    let complete =
        !report.timed_out && report.pairs_above_bound == 0 && report.staircase_at_stop.is_none();
    SolveMetric {
        target: [target.x, target.y],
        variables: n,
        equations: equations.len(),
        input_max_degree,
        inconsistent: report.inconsistent,
        complete,
        correct: complete && !report.inconsistent,
        timed_out: report.timed_out,
        pairs_above_bound: report.pairs_above_bound,
        degree_reached: report.degree_reached,
        solving_degree_max: report.solving_degree_max,
        max_cols_to_solution: report.max_cols_to_solution,
        field_ops: report.field_ops,
        milliseconds: report.ms,
    }
}

fn median_u32(values: impl Iterator<Item = u32>) -> Option<u32> {
    let mut values: Vec<u32> = values.collect();
    values.sort_unstable();
    (!values.is_empty()).then(|| values[values.len() / 2])
}

fn median_usize(values: impl Iterator<Item = usize>) -> Option<usize> {
    let mut values: Vec<usize> = values.collect();
    values.sort_unstable();
    (!values.is_empty()).then(|| values[values.len() / 2])
}

fn median_u64(values: impl Iterator<Item = u64>) -> Option<u64> {
    let mut values: Vec<u64> = values.collect();
    values.sort_unstable();
    (!values.is_empty()).then(|| values[values.len() / 2])
}

fn winner(terminals: &[TerminalResult]) -> Option<u64> {
    terminals
        .iter()
        .filter(|terminal| terminal.terminal != 0)
        .filter_map(|terminal| {
            let completed: Vec<&SolveMetric> = terminal
                .paired_runs
                .iter()
                .map(|run| &run.candidate)
                .filter(|metric| metric.complete && metric.correct)
                .collect();
            (completed.len() == TARGETS_PER_TERMINAL).then(|| {
                (
                    median_u32(completed.iter().map(|metric| metric.solving_degree_max)).unwrap(),
                    median_usize(completed.iter().map(|metric| metric.max_cols_to_solution))
                        .unwrap(),
                    median_u64(completed.iter().map(|metric| metric.field_ops)).unwrap(),
                    usize::MAX - terminal.liftable_roots,
                    terminal.terminal,
                )
            })
        })
        .min()
        .map(|key| key.4)
}

fn run_depth(depth: u32, a: u64, b: u64, all_points: &[Affine]) -> Result<DepthResult, String> {
    let fibres = dickson_fibres(depth);
    let width = 1usize << depth;
    let baseline_roots = fibres[0].clone();
    if baseline_roots.len() != width {
        return Err(format!(
            "terminal zero has {} roots at depth {depth}, expected {width}",
            baseline_roots.len()
        ));
    }
    let baseline_x: Vec<u64> = baseline_roots
        .iter()
        .copied()
        .filter(|&x| !square_roots(curve_rhs(x, a, b)).is_empty())
        .collect();
    let baseline_points = signed_points(&baseline_x, a, b);
    let baseline_decompositions = decompositions(&baseline_points, a);
    let baseline_count = baseline_x.len();

    let mut terminals = Vec::new();
    for terminal in 0..P {
        if terminal == 2 || terminal == P - 2 || fibres[terminal as usize].len() != width {
            continue;
        }
        let xs: Vec<u64> = fibres[terminal as usize]
            .iter()
            .copied()
            .filter(|&x| !square_roots(curve_rhs(x, a, b)).is_empty())
            .collect();
        if xs.len() < baseline_count {
            continue;
        }
        let points = signed_points(&xs, a, b);
        let candidate_decompositions = decompositions(&points, a);
        let common: Vec<Affine> = all_points
            .iter()
            .copied()
            .filter(|target| {
                baseline_decompositions.contains(target)
                    && candidate_decompositions.contains(target)
            })
            .collect();
        let selected: Vec<Affine> = common.iter().copied().take(TARGETS_PER_TERMINAL).collect();
        let paired_runs = selected
            .iter()
            .copied()
            .map(|target| PairedRun {
                candidate: solve(depth, terminal, target, a, b),
                terminal_zero: solve(depth, 0, target, a, b),
            })
            .collect();
        terminals.push(TerminalResult {
            terminal,
            full_roots: width,
            liftable_roots: xs.len(),
            common_targets_available: common.len(),
            selected_targets: selected.iter().map(|target| [target.x, target.y]).collect(),
            paired_runs,
        });
    }
    let winner_terminal = winner(&terminals);
    Ok(DepthResult {
        depth,
        terminal_zero_liftable_roots: baseline_count,
        eligible_terminals: terminals.len(),
        winner_terminal,
        terminals,
    })
}

fn run(cli: Cli) -> Result<(), String> {
    let p256 = CurveParams::p256();
    let b = (&p256.b % P).to_u64().ok_or("P-256 b reduction failed")?;
    let a = P - 3;
    let all_points = curve_points(a, b);
    let result = ExperimentResult {
        schema: "p256.dickson_terminal_degree/v1".into(),
        curve: crypto_lib::cryptanalysis::p256_dickson_factor_base::CURVE_SLUG.into(),
        prime: P,
        a,
        b,
        max_degree: MAX_DEGREE,
        budget_seconds: BUDGET_SECS,
        targets_per_terminal: TARGETS_PER_TERMINAL,
        depths: [4, 5]
            .into_iter()
            .map(|depth| run_depth(depth, a, b, &all_points))
            .collect::<Result<Vec<_>, _>>()?,
    };
    if result
        .depths
        .iter()
        .flat_map(|depth| &depth.terminals)
        .flat_map(|terminal| &terminal.paired_runs)
        .any(|run| {
            (run.candidate.complete && !run.candidate.correct)
                || (run.terminal_zero.complete && !run.terminal_zero.correct)
        })
    {
        return Err("a certified F4 verdict contradicted an exhaustive positive target".into());
    }
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
    fn frozen_prime_has_full_zero_fibres() {
        for depth in [4, 5] {
            assert_eq!(dickson_fibres(depth)[0].len(), 1usize << depth);
        }
    }

    #[test]
    fn s3_system_has_degree_four() {
        let p256 = CurveParams::p256();
        let b = (&p256.b % P).to_u64().unwrap();
        let target = curve_points(P - 3, b)[0];
        let (system, _) = build_system(4, 0, target.x, P - 3, b);
        assert_eq!(system.iter().map(Pol::degree).max(), Some(4));
    }
}

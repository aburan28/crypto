//! Residual-depth scaling for quadratic S3 on branched Dickson fibres.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::path::PathBuf;
use std::process::ExitCode;
use std::time::Duration;

use clap::Parser;
use crypto_lib::cryptanalysis::f4_fp::{self, F4Options, Ordering, Poly};
use crypto_lib::ecc::curve::CurveParams;
use crypto_lib::hash::sha256;
use num_traits::ToPrimitive;
use serde::Serialize;

const P: u64 = 7681;
const DEPTHS: [u32; 4] = [5, 6, 7, 8];
const RESIDUAL_DEPTHS: [u32; 3] = [1, 2, 3];
const MAX_DEPTH: u32 = 8;
const MAX_DEGREE: u32 = 12;
const BUDGET_SECS: u64 = 5;
const MIN_LIFTABLE_ROOTS: usize = 96;
const TERMINALS: usize = 2;

#[derive(Parser)]
#[command(about = "Scale residual-depth quadratic S3 across frozen Dickson fibres")]
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
struct ComponentRecord {
    boundary1: u64,
    boundary2: u64,
    expected_positive: bool,
    inconsistent: bool,
    complete: bool,
    correct: bool,
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
    residual_depth: u32,
    branch_depth: u32,
    target_kind: String,
    target: [u64; 2],
    variables: usize,
    equations: usize,
    input_max_degree: u32,
    boundaries: usize,
    components: usize,
    expected_positive_components: usize,
    consistent_components: usize,
    complete_components: usize,
    correct_components: usize,
    timed_out_components: usize,
    pairs_above_bound_total: usize,
    degree_histogram: BTreeMap<u32, usize>,
    max_solving_degree: u32,
    max_cols_to_solution: usize,
    total_field_ops: u64,
    component_stream_sha256: String,
    positive_components: Vec<ComponentRecord>,
    exceptional_components: Vec<ComponentRecord>,
    max_degree_example: Option<ComponentRecord>,
    max_columns_example: Option<ComponentRecord>,
    max_ops_example: Option<ComponentRecord>,
}

#[derive(Clone, Debug, Serialize)]
struct TerminalSelection {
    terminal: u64,
    depth8_roots: usize,
    depth8_liftable_x: usize,
    depth8_signed_points: usize,
}

#[derive(Clone, Debug, Serialize)]
struct ExperimentResult {
    schema: String,
    curve: String,
    prime: u64,
    a: u64,
    b: u64,
    depths: Vec<u32>,
    residual_depths: Vec<u32>,
    max_degree: u32,
    budget_seconds_per_component: u64,
    minimum_liftable_depth8_roots: usize,
    terminal_selection: Vec<TerminalSelection>,
    total_curve_points: usize,
    cells: Vec<CellResult>,
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

fn square_root_table() -> Vec<Vec<u64>> {
    let mut table = vec![Vec::new(); P as usize];
    for y in 0..P {
        table[mulm(y, y) as usize].push(y);
    }
    table
}

fn curve_rhs(x: u64, a: u64, b: u64) -> u64 {
    addm(addm(mulm(mulm(x, x), x), mulm(a, x)), b)
}

fn curve_points(a: u64, b: u64, roots: &[Vec<u64>]) -> Vec<Affine> {
    let mut points = Vec::new();
    for x in 0..P {
        for &y in &roots[curve_rhs(x, a, b) as usize] {
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

fn signed_points(xs: &[u64], a: u64, b: u64, square_roots: &[Vec<u64>]) -> Vec<Affine> {
    xs.iter()
        .flat_map(|&x| {
            square_roots[curve_rhs(x, a, b) as usize]
                .iter()
                .copied()
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
    residual_depth: u32,
    boundary1: u64,
    boundary2: u64,
    target_x: u64,
    a: u64,
    b: u64,
) -> (Vec<Pol>, usize) {
    let t = residual_depth as usize;
    let block = t + 2;
    let n = 2 * block + 1;
    let mut equations = Vec::new();
    let x1 = add_chain(&mut equations, n, 0, residual_depth, boundary1);
    let x2 = add_chain(&mut equations, n, block, residual_depth, boundary2);
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
    residual_depth: u32,
    boundary1: u64,
    boundary2: u64,
    target: Affine,
    expected_positive: bool,
    a: u64,
    b: u64,
) -> ComponentRecord {
    let (equations, n) = system(residual_depth, boundary1, boundary2, target.x, a, b);
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
    ComponentRecord {
        boundary1,
        boundary2,
        expected_positive,
        inconsistent: report.inconsistent,
        complete,
        correct: complete && (report.inconsistent != expected_positive),
        timed_out: report.timed_out,
        pairs_above_bound: report.pairs_above_bound,
        solving_degree_max: report.solving_degree_max,
        max_cols_to_solution: report.max_cols_to_solution,
        field_ops: report.field_ops,
    }
}

fn digest_hex(bytes: &[u8]) -> String {
    hex::encode(sha256(bytes))
}

fn choose_terminals(a: u64, b: u64, square_roots: &[Vec<u64>]) -> Vec<TerminalSelection> {
    let mut selected = Vec::new();
    for terminal in 0..P {
        let xs = roots(MAX_DEPTH, terminal);
        if xs.len() != 1usize << MAX_DEPTH {
            continue;
        }
        let liftable_x = xs
            .iter()
            .filter(|&&x| !square_roots[curve_rhs(x, a, b) as usize].is_empty())
            .count();
        if liftable_x < MIN_LIFTABLE_ROOTS {
            continue;
        }
        let signed = signed_points(&xs, a, b, square_roots).len();
        selected.push(TerminalSelection {
            terminal,
            depth8_roots: xs.len(),
            depth8_liftable_x: liftable_x,
            depth8_signed_points: signed,
        });
        if selected.len() == TERMINALS {
            break;
        }
    }
    selected
}

fn choose_targets(
    terminal: u64,
    depth: u32,
    a: u64,
    b: u64,
    square_roots: &[Vec<u64>],
    all_points: &[Affine],
) -> Result<(Affine, Affine), String> {
    let xs: Vec<u64> = roots(depth, terminal)
        .into_iter()
        .filter(|&x| !square_roots[curve_rhs(x, a, b) as usize].is_empty())
        .collect();
    let points = signed_points(&xs, a, b, square_roots);
    let sum_set = sums(&points, &points, a);
    let positive = all_points
        .iter()
        .copied()
        .find(|point| sum_set.contains(point))
        .ok_or_else(|| format!("no positive target for terminal {terminal}, depth {depth}"))?;
    let negative = all_points
        .iter()
        .copied()
        .find(|point| !sum_set.contains(point))
        .ok_or_else(|| format!("no negative target for terminal {terminal}, depth {depth}"))?;
    Ok((positive, negative))
}

fn canonical_record(stream: &mut String, record: &ComponentRecord) {
    writeln!(
        stream,
        "{},{},{},{},{},{},{},{},{},{},{}",
        record.boundary1,
        record.boundary2,
        u8::from(record.expected_positive),
        u8::from(record.inconsistent),
        u8::from(record.complete),
        u8::from(record.correct),
        u8::from(record.timed_out),
        record.pairs_above_bound,
        record.solving_degree_max,
        record.max_cols_to_solution,
        record.field_ops,
    )
    .expect("writing to String cannot fail");
}

fn better_example<F>(
    current: &Option<ComponentRecord>,
    candidate: &ComponentRecord,
    value: F,
) -> bool
where
    F: Fn(&ComponentRecord) -> u64,
{
    current
        .as_ref()
        .is_none_or(|existing| value(candidate) > value(existing))
}

#[allow(clippy::too_many_arguments)]
fn run_cell(
    terminal: u64,
    depth: u32,
    residual_depth: u32,
    target_kind: &str,
    target: Affine,
    a: u64,
    b: u64,
    square_roots: &[Vec<u64>],
) -> Result<CellResult, String> {
    let branch_depth = depth - residual_depth;
    let boundaries = roots(branch_depth, terminal);
    let expected_boundaries = 1usize << branch_depth;
    if boundaries.len() != expected_boundaries {
        return Err(format!(
            "terminal {terminal}, depth {depth}, residual {residual_depth}: expected {expected_boundaries} boundaries, got {}",
            boundaries.len()
        ));
    }
    let point_sets: BTreeMap<u64, Vec<Affine>> = boundaries
        .iter()
        .map(|&boundary| {
            let xs: Vec<u64> = roots(residual_depth, boundary)
                .into_iter()
                .filter(|&x| !square_roots[curve_rhs(x, a, b) as usize].is_empty())
                .collect();
            (boundary, signed_points(&xs, a, b, square_roots))
        })
        .collect();

    let mut stream = String::new();
    let mut degree_histogram = BTreeMap::new();
    let mut expected_positive_components = 0;
    let mut consistent_components = 0;
    let mut complete_components = 0;
    let mut correct_components = 0;
    let mut timed_out_components = 0;
    let mut pairs_above_bound_total = 0;
    let mut max_solving_degree = 0;
    let mut max_cols_to_solution = 0;
    let mut total_field_ops = 0;
    let mut positive_components = Vec::new();
    let mut exceptional_components = Vec::new();
    let mut max_degree_example = None;
    let mut max_columns_example = None;
    let mut max_ops_example = None;

    for &boundary1 in &boundaries {
        for &boundary2 in &boundaries {
            let expected_positive =
                sums(&point_sets[&boundary1], &point_sets[&boundary2], a).contains(&target);
            let record = measure(
                residual_depth,
                boundary1,
                boundary2,
                target,
                expected_positive,
                a,
                b,
            );
            canonical_record(&mut stream, &record);
            *degree_histogram
                .entry(record.solving_degree_max)
                .or_insert(0) += 1;
            expected_positive_components += usize::from(record.expected_positive);
            consistent_components += usize::from(!record.inconsistent);
            complete_components += usize::from(record.complete);
            correct_components += usize::from(record.correct);
            timed_out_components += usize::from(record.timed_out);
            pairs_above_bound_total += record.pairs_above_bound;
            max_solving_degree = max_solving_degree.max(record.solving_degree_max);
            max_cols_to_solution = max_cols_to_solution.max(record.max_cols_to_solution);
            total_field_ops += record.field_ops;
            if record.expected_positive {
                positive_components.push(record.clone());
            }
            if !record.correct || !record.complete {
                exceptional_components.push(record.clone());
            }
            if better_example(&max_degree_example, &record, |entry| {
                u64::from(entry.solving_degree_max)
            }) {
                max_degree_example = Some(record.clone());
            }
            if better_example(&max_columns_example, &record, |entry| {
                entry.max_cols_to_solution as u64
            }) {
                max_columns_example = Some(record.clone());
            }
            if better_example(&max_ops_example, &record, |entry| entry.field_ops) {
                max_ops_example = Some(record);
            }
        }
    }

    let (template, variables) = system(residual_depth, 0, 0, target.x, a, b);
    Ok(CellResult {
        terminal,
        depth,
        residual_depth,
        branch_depth,
        target_kind: target_kind.into(),
        target: [target.x, target.y],
        variables,
        equations: template.len(),
        input_max_degree: template.iter().map(Pol::degree).max().unwrap_or(0),
        boundaries: boundaries.len(),
        components: boundaries.len() * boundaries.len(),
        expected_positive_components,
        consistent_components,
        complete_components,
        correct_components,
        timed_out_components,
        pairs_above_bound_total,
        degree_histogram,
        max_solving_degree,
        max_cols_to_solution,
        total_field_ops,
        component_stream_sha256: digest_hex(stream.as_bytes()),
        positive_components,
        exceptional_components,
        max_degree_example,
        max_columns_example,
        max_ops_example,
    })
}

fn run(cli: Cli) -> Result<(), String> {
    let p256 = CurveParams::p256();
    let b = (&p256.b % P).to_u64().ok_or("P-256 b reduction failed")?;
    let a = P - 3;
    let square_roots = square_root_table();
    let all_points = curve_points(a, b, &square_roots);
    let terminal_selection = choose_terminals(a, b, &square_roots);
    if terminal_selection.len() != TERMINALS {
        return Err(format!(
            "expected {TERMINALS} eligible terminals, found {}",
            terminal_selection.len()
        ));
    }

    let mut cells = Vec::new();
    for selection in &terminal_selection {
        for depth in DEPTHS {
            let (positive, negative) =
                choose_targets(selection.terminal, depth, a, b, &square_roots, &all_points)?;
            for residual_depth in RESIDUAL_DEPTHS {
                for (target_kind, target) in [("positive", positive), ("negative", negative)] {
                    cells.push(run_cell(
                        selection.terminal,
                        depth,
                        residual_depth,
                        target_kind,
                        target,
                        a,
                        b,
                        &square_roots,
                    )?);
                }
            }
        }
    }

    let result = ExperimentResult {
        schema: "p256.dickson_residual_scaling/v1".into(),
        curve: crypto_lib::cryptanalysis::p256_dickson_factor_base::CURVE_SLUG.into(),
        prime: P,
        a,
        b,
        depths: DEPTHS.to_vec(),
        residual_depths: RESIDUAL_DEPTHS.to_vec(),
        max_degree: MAX_DEGREE,
        budget_seconds_per_component: BUDGET_SECS,
        minimum_liftable_depth8_roots: MIN_LIFTABLE_ROOTS,
        terminal_selection,
        total_curve_points: all_points.len(),
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
    fn selected_terminals_have_full_depth8_fibres() {
        let p256 = CurveParams::p256();
        let b = (&p256.b % P).to_u64().unwrap();
        let square_roots = square_root_table();
        let selected = choose_terminals(P - 3, b, &square_roots);
        assert_eq!(selected.len(), TERMINALS);
        for terminal in selected {
            assert_eq!(terminal.depth8_roots, 1 << MAX_DEPTH);
            assert!(terminal.depth8_liftable_x >= MIN_LIFTABLE_ROOTS);
        }
    }

    #[test]
    fn residual_system_remains_quadratic() {
        let p256 = CurveParams::p256();
        let b = (&p256.b % P).to_u64().unwrap();
        for residual_depth in RESIDUAL_DEPTHS {
            let (equations, _) = system(residual_depth, 1, 2, 3, P - 3, b);
            assert_eq!(equations.iter().map(Pol::degree).max(), Some(2));
        }
    }
}

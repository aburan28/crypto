//! Paired small-prime solving-degree experiment for P-256 factor-base shapes.

use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;
use std::process::ExitCode;
use std::time::Duration;

use clap::Parser;
use crypto_lib::cryptanalysis::f4_fp::{self, F4Options, Ordering, Poly};
use crypto_lib::ecc::curve::CurveParams;
use num_traits::ToPrimitive;
use serde::Serialize;

const P: u64 = 127;
const MAX_DEGREE: u32 = 12;
const BUDGET_SECS: u64 = 30;
const TARGETS_PER_CLASS: usize = 3;

#[derive(Parser)]
#[command(about = "Measure Dickson and affine-bitbox summation-system solving degree")]
struct Cli {
    /// Compact raw JSON output; stdout when omitted.
    #[arg(long)]
    out: Option<PathBuf>,
    /// Run only one frozen depth (2, 3, 4, or 5).
    #[arg(long)]
    depth: Option<u32>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize)]
#[serde(rename_all = "kebab-case")]
enum Family {
    Dickson,
    AffineBitbox,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize)]
#[serde(rename_all = "kebab-case")]
enum Formulation {
    S3,
    QuadraticIncidence,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "kebab-case")]
enum Expected {
    Positive,
    Negative,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
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
struct FactorBaseSummary {
    family: Family,
    offset: Option<u64>,
    abscissae: Vec<u64>,
    signed_points: usize,
}

#[derive(Clone, Debug, Serialize)]
struct TargetSummary {
    positive_available: usize,
    negative_available: usize,
    selected_positive: Vec<[u64; 2]>,
    selected_negative: Vec<[u64; 2]>,
}

#[derive(Clone, Debug, Serialize)]
struct Run {
    depth: u32,
    formulation: Formulation,
    family: Family,
    expected: Expected,
    target: [u64; 2],
    variables: usize,
    equations: usize,
    input_max_degree: u32,
    inconsistent: bool,
    complete: bool,
    certified_correct: bool,
    certified_mismatch: bool,
    timed_out: bool,
    pairs_above_bound: usize,
    degree_reached: u32,
    solving_degree_max: u32,
    max_cols_to_solution: usize,
    field_ops: u64,
    milliseconds: f64,
}

#[derive(Clone, Debug, Serialize)]
struct DepthResult {
    depth: u32,
    dickson: FactorBaseSummary,
    affine_bitbox: FactorBaseSummary,
    targets: TargetSummary,
    runs: Vec<Run>,
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
    targets_per_class: usize,
    depths: Vec<DepthResult>,
}

fn addm(left: u64, right: u64) -> u64 {
    (left + right) % P
}

fn subm(left: u64, right: u64) -> u64 {
    (left + P - right % P) % P
}

fn mulm(left: u64, right: u64) -> u64 {
    (left * right) % P
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

fn dickson_abscissae(depth: u32, a: u64, b: u64) -> Vec<u64> {
    (0..P)
        .filter(|&x| {
            let mut value = x;
            for _ in 0..depth {
                value = subm(mulm(value, value), 2);
            }
            value == 0 && !square_roots(curve_rhs(x, a, b)).is_empty()
        })
        .collect()
}

fn bitbox_abscissae(depth: u32, offset: u64, a: u64, b: u64) -> Vec<u64> {
    let width = 1u64 << depth;
    (offset..offset + width)
        .filter(|&x| !square_roots(curve_rhs(x, a, b)).is_empty())
        .collect()
}

fn select_bitbox(depth: u32, wanted: usize, a: u64, b: u64) -> (u64, Vec<u64>) {
    let width = 1u64 << depth;
    (0..=P - width)
        .map(|offset| {
            let xs = bitbox_abscissae(depth, offset, a, b);
            let difference = xs.len().abs_diff(wanted);
            (difference, offset, xs)
        })
        .min_by_key(|(difference, offset, _)| (*difference, *offset))
        .map(|(_, offset, xs)| (offset, xs))
        .expect("the frozen toy interval scan is non-empty")
}

fn decomposition_sets(points: &[Affine], a: u64) -> (BTreeSet<Affine>, BTreeSet<Affine>) {
    let mut all = BTreeSet::new();
    let mut generic = BTreeSet::new();
    for &left in points {
        for &right in points {
            if let Point::Affine(sum) = add_points(left, right, a) {
                all.insert(sum);
                if left.x != right.x {
                    generic.insert(sum);
                }
            }
        }
    }
    (all, generic)
}

fn targets(
    all_points: &[Affine],
    dickson_points: &[Affine],
    bitbox_points: &[Affine],
    a: u64,
) -> (Vec<Affine>, Vec<Affine>, TargetSummary) {
    let (dickson_all, dickson_generic) = decomposition_sets(dickson_points, a);
    let (bitbox_all, bitbox_generic) = decomposition_sets(bitbox_points, a);
    let positives: Vec<Affine> = all_points
        .iter()
        .copied()
        .filter(|target| dickson_generic.contains(target) && bitbox_generic.contains(target))
        .collect();
    let negatives: Vec<Affine> = all_points
        .iter()
        .copied()
        .filter(|target| !dickson_all.contains(target) && !bitbox_all.contains(target))
        .collect();
    let selected_positive: Vec<Affine> =
        positives.iter().copied().take(TARGETS_PER_CLASS).collect();
    let selected_negative: Vec<Affine> =
        negatives.iter().copied().take(TARGETS_PER_CLASS).collect();
    let summary = TargetSummary {
        positive_available: positives.len(),
        negative_available: negatives.len(),
        selected_positive: selected_positive.iter().map(|p| [p.x, p.y]).collect(),
        selected_negative: selected_negative.iter().map(|p| [p.x, p.y]).collect(),
    };
    (selected_positive, selected_negative, summary)
}

fn linear_bitbox_x(n: usize, base: usize, depth: u32, offset: u64) -> Pol {
    let mut x = Pol::constant(n, offset);
    let mut weight = 1;
    for index in 0..depth as usize {
        x = x.add(&Pol::var(n, base + index).scale(weight));
        weight = addm(weight, weight);
    }
    x
}

fn add_domain(
    equations: &mut Vec<Pol>,
    n: usize,
    base: usize,
    depth: u32,
    family: Family,
    offset: u64,
) -> Pol {
    match family {
        Family::Dickson => {
            for index in 0..depth as usize {
                let current = Pol::var(n, base + index);
                let mut equation = current.square().sub(&Pol::constant(n, 2));
                if index + 1 < depth as usize {
                    equation = equation.sub(&Pol::var(n, base + index + 1));
                }
                equations.push(equation);
            }
            Pol::var(n, base)
        }
        Family::AffineBitbox => {
            for index in 0..depth as usize {
                let bit = Pol::var(n, base + index);
                equations.push(bit.square().sub(&bit));
            }
            linear_bitbox_x(n, base, depth, offset)
        }
    }
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

fn build_system(
    depth: u32,
    family: Family,
    formulation: Formulation,
    offset: u64,
    target: Affine,
    a: u64,
    b: u64,
) -> (Vec<Pol>, usize) {
    match formulation {
        Formulation::S3 => {
            let block = depth as usize + 1;
            let n = 2 * block;
            let mut equations = Vec::new();
            let x1 = add_domain(&mut equations, n, 0, depth, family, offset);
            let x2 = add_domain(&mut equations, n, block, depth, family, offset);
            for (base, x) in [(0, &x1), (block, &x2)] {
                let y = Pol::var(n, base + depth as usize);
                let rhs = x.square().mul(x).add(&x.scale(a)).add(&Pol::constant(n, b));
                equations.push(y.square().sub(&rhs));
            }
            equations.push(s3(&x1, &x2, &Pol::constant(n, target.x), a, b));
            (equations, n)
        }
        Formulation::QuadraticIncidence => {
            let block = depth as usize + 2;
            let lambda_index = 2 * block;
            let inverse_index = lambda_index + 1;
            let n = inverse_index + 1;
            let mut equations = Vec::new();
            let x1 = add_domain(&mut equations, n, 0, depth, family, offset);
            let x2 = add_domain(&mut equations, n, block, depth, family, offset);
            let y1 = Pol::var(n, depth as usize);
            let y2 = Pol::var(n, block + depth as usize);
            for (base, x, y) in [(0, &x1, &y1), (block, &x2, &y2)] {
                let square = Pol::var(n, base + depth as usize + 1);
                equations.push(square.sub(&x.square()));
                equations.push(
                    y.square()
                        .sub(&square.mul(x))
                        .sub(&x.scale(a))
                        .sub(&Pol::constant(n, b)),
                );
            }
            let lambda = Pol::var(n, lambda_index);
            let inverse = Pol::var(n, inverse_index);
            let denominator = x2.sub(&x1);
            equations.push(denominator.mul(&inverse).sub(&Pol::constant(n, 1)));
            equations.push(lambda.mul(&denominator).sub(&y2.sub(&y1)));
            equations.push(
                lambda
                    .square()
                    .sub(&x1)
                    .sub(&x2)
                    .sub(&Pol::constant(n, target.x)),
            );
            equations.push(
                lambda
                    .mul(&x1.sub(&Pol::constant(n, target.x)))
                    .sub(&y1)
                    .sub(&Pol::constant(n, target.y)),
            );
            (equations, n)
        }
    }
}

fn measure(
    depth: u32,
    formulation: Formulation,
    family: Family,
    expected: Expected,
    target: Affine,
    offset: u64,
    a: u64,
    b: u64,
) -> Run {
    let (equations, n) = build_system(depth, family, formulation, offset, target, a, b);
    let input_max_degree = equations.iter().map(Pol::degree).max().unwrap_or(0);
    let input: Vec<Poly> = equations.iter().map(Pol::to_f4).collect();
    let options =
        F4Options::new(Ordering::Grevlex, MAX_DEGREE).with_budget(Duration::from_secs(BUDGET_SECS));
    let report = f4_fp::f4(&input, n, P, &options);
    let complete =
        !report.timed_out && report.pairs_above_bound == 0 && report.staircase_at_stop.is_none();
    let observed_positive = !report.inconsistent;
    let expected_positive = expected == Expected::Positive;
    Run {
        depth,
        formulation,
        family,
        expected,
        target: [target.x, target.y],
        variables: n,
        equations: equations.len(),
        input_max_degree,
        inconsistent: report.inconsistent,
        complete,
        certified_correct: complete && observed_positive == expected_positive,
        certified_mismatch: complete && observed_positive != expected_positive,
        timed_out: report.timed_out,
        pairs_above_bound: report.pairs_above_bound,
        degree_reached: report.degree_reached,
        solving_degree_max: report.solving_degree_max,
        max_cols_to_solution: report.max_cols_to_solution,
        field_ops: report.field_ops,
        milliseconds: report.ms,
    }
}

fn run_depth(depth: u32, a: u64, b: u64, all_points: &[Affine]) -> DepthResult {
    let dickson_x = dickson_abscissae(depth, a, b);
    let (offset, bitbox_x) = select_bitbox(depth, dickson_x.len(), a, b);
    let dickson_points = signed_points(&dickson_x, a, b);
    let bitbox_points = signed_points(&bitbox_x, a, b);
    let (positives, negatives, target_summary) =
        targets(all_points, &dickson_points, &bitbox_points, a);
    let mut runs = Vec::new();
    for formulation in [Formulation::S3, Formulation::QuadraticIncidence] {
        for family in [Family::Dickson, Family::AffineBitbox] {
            let family_offset = if family == Family::AffineBitbox {
                offset
            } else {
                0
            };
            for &target in &positives {
                runs.push(measure(
                    depth,
                    formulation,
                    family,
                    Expected::Positive,
                    target,
                    family_offset,
                    a,
                    b,
                ));
            }
            for &target in &negatives {
                runs.push(measure(
                    depth,
                    formulation,
                    family,
                    Expected::Negative,
                    target,
                    family_offset,
                    a,
                    b,
                ));
            }
        }
    }
    DepthResult {
        depth,
        dickson: FactorBaseSummary {
            family: Family::Dickson,
            offset: None,
            abscissae: dickson_x,
            signed_points: dickson_points.len(),
        },
        affine_bitbox: FactorBaseSummary {
            family: Family::AffineBitbox,
            offset: Some(offset),
            abscissae: bitbox_x,
            signed_points: bitbox_points.len(),
        },
        targets: target_summary,
        runs,
    }
}

fn run(cli: Cli) -> Result<(), String> {
    let p256 = CurveParams::p256();
    let b = (&p256.b % P).to_u64().ok_or("P-256 b reduction failed")?;
    let a = P - 3;
    let all_points = curve_points(a, b);
    let depths: Vec<u32> = match cli.depth {
        Some(depth) if (2..=5).contains(&depth) => vec![depth],
        Some(depth) => return Err(format!("depth {depth} is outside the frozen set 2..=5")),
        None => (2..=5).collect(),
    };
    let result = ExperimentResult {
        schema: "p256.factor_base_degree/v1".into(),
        curve: crypto_lib::cryptanalysis::p256_dickson_factor_base::CURVE_SLUG.into(),
        prime: P,
        a,
        b,
        max_degree: MAX_DEGREE,
        budget_seconds: BUDGET_SECS,
        targets_per_class: TARGETS_PER_CLASS,
        depths: depths
            .into_iter()
            .map(|depth| run_depth(depth, a, b, &all_points))
            .collect(),
    };
    let text = serde_json::to_string_pretty(&result).map_err(|error| error.to_string())? + "\n";
    match cli.out {
        Some(path) => std::fs::write(path, text).map_err(|error| error.to_string())?,
        None => print!("{text}"),
    }
    if result
        .depths
        .iter()
        .flat_map(|depth| &depth.runs)
        .any(|run| run.certified_mismatch)
    {
        return Err("a certified F4 verdict disagreed with exhaustive enumeration".into());
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
    fn p127_curve_and_dickson_fibres_exist() {
        let p256 = CurveParams::p256();
        let b = (&p256.b % P).to_u64().unwrap();
        let points = curve_points(P - 3, b);
        assert!(!points.is_empty());
        for depth in 2..=5 {
            assert!(!dickson_abscissae(depth, P - 3, b).is_empty());
        }
    }

    #[test]
    fn formulations_have_frozen_input_degrees() {
        let p256 = CurveParams::p256();
        let b = (&p256.b % P).to_u64().unwrap();
        let target = curve_points(P - 3, b)[0];
        for family in [Family::Dickson, Family::AffineBitbox] {
            let (s3_system, _) = build_system(2, family, Formulation::S3, 0, target, P - 3, b);
            assert_eq!(s3_system.iter().map(Pol::degree).max(), Some(4));
            let (incidence, _) = build_system(
                2,
                family,
                Formulation::QuadraticIncidence,
                0,
                target,
                P - 3,
                b,
            );
            assert_eq!(incidence.iter().map(Pol::degree).max(), Some(2));
        }
    }
}

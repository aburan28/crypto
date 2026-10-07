//! Round 296: exact Dickson-aware divisor/Riemann--Roch selector screen.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::{Duration, Instant};

use clap::Parser;
use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use crypto_lib::cryptanalysis::f4_fp::{self, F4Options, Ordering, Poly as F4Poly};
use crypto_lib::cryptanalysis::p256_dickson_factor_base::{
    self, CURVE_SLUG, WideFactorBaseDump, WideFbPoint,
};
use crypto_lib::ecc::{CurveParams, Point};
use crypto_lib::hash::sha256;
use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};
use serde::Serialize;

const ROUND295_SHA256: &str = "9550f2bbaceb9e297fb35c12e79480581ca80bd58b03a8493e86ea7c1eda91a5";
const FACTOR_BASE_SHA256: &str = "d27516ca40a612ecf3ebabfa8ae04776084c7948110e20f221da438e1e64d8f8";
const FB_ID: &str = "FB1hc72514a2a8d3";
const P: u64 = 1151;
const A: u64 = P - 3;
const MAX_DEGREE: u32 = 7;
const BUDGET_SECS: u64 = 20;
const TOY_SIZES: [(usize, usize); 4] = [(3, 5), (5, 8), (7, 11), (9, 14)];

#[derive(Parser)]
#[command(about = "Evaluate the exact P-256 Dickson/RR selector formulation")]
struct Cli {
    #[arg(long)]
    round295: PathBuf,
    #[arg(long)]
    factor_base: PathBuf,
    #[arg(long)]
    out: PathBuf,
    #[arg(long)]
    telemetry_out: PathBuf,
}

#[derive(Clone, Debug, Serialize)]
struct Dependency {
    path: String,
    bytes: u64,
    sha256: String,
}

#[derive(Clone, Debug, Serialize)]
struct MacaulayDimension {
    degree: u32,
    columns_at_most_degree: String,
    generic_rows_at_most_degree: String,
    dense_bytes_at_32_bytes_per_entry: String,
}

#[derive(Clone, Debug, Serialize)]
struct SelectorShape {
    arity: usize,
    columns: usize,
    complement_degree: usize,
    projective_variables: usize,
    projective_equations: usize,
    normalized_variables: usize,
    normalized_equations: usize,
    maximum_input_degree: u32,
    input_monomials: usize,
    macaulay: Vec<MacaulayDimension>,
}

#[derive(Clone, Debug, Serialize)]
struct P256Certificate {
    fb_id: String,
    columns: usize,
    signed_points: usize,
    depth8_columns: usize,
    depth6_columns: usize,
    component_overlap: usize,
    common_chain_rows_checked: usize,
    common_chain_exact: bool,
    overbroad_depth8_columns: usize,
    overbroad_union_columns: usize,
    overbroad_extra_columns: usize,
    support_degree: usize,
    support_coefficients_sha256: String,
    support_roots_checked: usize,
    support_derivatives_nonzero: usize,
    support_replay_byte_identical: bool,
    public_scalar: String,
    public_target: [String; 2],
    selector: SelectorShape,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
struct Affine {
    x: u64,
    y: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ToyPoint {
    Infinity,
    Affine(Affine),
}

#[derive(Clone, Debug, Serialize)]
struct TargetReceipt {
    kind: String,
    target: [u64; 2],
    planted_columns: Option<Vec<usize>>,
    planted_negative_mask: Option<u64>,
    candidates: u64,
    relations: u64,
    replayed_relations: u64,
    replay_failures: u64,
    witness_sha256: String,
    planted_equation_witness_verified: Option<bool>,
    variables: usize,
    equations: usize,
    input_monomials: usize,
    input_maximum_degree: u32,
    f4_inconsistent: bool,
    f4_classification: Option<bool>,
    f4_complete: bool,
    f4_timed_out: bool,
    f4_pairs_above_bound: usize,
    solving_degree_max: u32,
    max_rows: usize,
    max_cols: usize,
    max_cols_to_solution: usize,
    field_operations: u64,
    false_positive: bool,
    false_negative: bool,
}

#[derive(Clone, Debug, Serialize)]
struct TargetTelemetry {
    arity: usize,
    columns: usize,
    kind: String,
    wall_ms: f64,
}

#[derive(Clone, Debug, Serialize)]
struct TelemetryReceipt {
    schema: String,
    targets: Vec<TargetTelemetry>,
    total_wall_ms: f64,
    time_fit_exponent_in_variables: Option<f64>,
}

#[derive(Clone, Debug, Serialize)]
struct ToyCase {
    arity: usize,
    columns: usize,
    support_x: Vec<u64>,
    support_sha256: String,
    signed_domain: u64,
    planted: TargetReceipt,
    unplanted: TargetReceipt,
}

#[derive(Clone, Debug, Serialize)]
struct Gates {
    common_chain_exact: bool,
    planted_witnesses_exact: bool,
    zero_false_positives_and_false_negatives_on_complete_instances: bool,
    largest_three_complete_sizes_available: bool,
    structured_residual_degree_at_most_5: bool,
    non_increasing_degree_slope: bool,
    projected_materialized_storage_below_2_50: bool,
    projected_cost_per_usable_relation_below_2_103: bool,
    projected_collection_below_2_120: bool,
    promoted: bool,
}

#[derive(Serialize)]
struct ResultReceipt {
    schema: String,
    curve: String,
    screening_round: u32,
    execution_status: String,
    round295: Dependency,
    factor_base: Dependency,
    p256: P256Certificate,
    toy_prime: u64,
    toy_a: u64,
    toy_b: u64,
    toy_support_pool: Vec<u64>,
    toy_support_pool_sha256: String,
    cases: Vec<ToyCase>,
    complete_targets: usize,
    censored_targets: usize,
    false_positives: u64,
    false_negatives: u64,
    degree_fit_slope: Option<f64>,
    field_ops_fit_exponent_in_variables: Option<f64>,
    width_fit_exponent_in_variables: Option<f64>,
    p256_relation_attempted: bool,
    p256_relations_reported: u64,
    gates: Gates,
    classification: String,
    dominant_obstruction: String,
    decision: String,
    semantic_evidence_sha256: String,
    result_json_bytes: u64,
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
                let exponent = left_exp
                    .iter()
                    .zip(right_exp)
                    .map(|(left, right)| left + right)
                    .collect::<Vec<_>>();
                let value = out.terms.entry(exponent).or_insert(0);
                *value = addm(*value, mulm(*left_coefficient, *right_coefficient));
            }
        }
        out.terms.retain(|_, coefficient| *coefficient != 0);
        out
    }

    fn degree(&self) -> u32 {
        self.terms
            .keys()
            .map(|exponent| exponent.iter().sum())
            .max()
            .unwrap_or(0)
    }

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
        let terms = self
            .terms
            .iter()
            .map(|(exponent, coefficient)| (exponent.clone(), *coefficient))
            .collect::<Vec<_>>();
        f4_fp::normalise(&terms, P, Ordering::Grevlex)
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

fn invm(value: u64) -> u64 {
    powm(value, P - 2)
}

fn negm(value: u64) -> u64 {
    if value == 0 { 0 } else { P - value }
}

fn poly_mul_u64(left: &[u64], right: &[u64]) -> Vec<u64> {
    let mut out = vec![0; left.len() + right.len() - 1];
    for (i, &a) in left.iter().enumerate() {
        for (j, &b) in right.iter().enumerate() {
            out[i + j] = addm(out[i + j], mulm(a, b));
        }
    }
    out
}

fn poly_sub_u64(left: &[u64], right: &[u64]) -> Vec<u64> {
    let len = left.len().max(right.len());
    (0..len)
        .map(|i| subm(*left.get(i).unwrap_or(&0), *right.get(i).unwrap_or(&0)))
        .collect()
}

fn product_roots_u64(roots: &[u64]) -> Vec<u64> {
    roots
        .iter()
        .fold(vec![1], |poly, &root| poly_mul_u64(&poly, &[negm(root), 1]))
}

fn parse_big(value: &str) -> Result<BigUint, String> {
    let (digits, radix) = value
        .strip_prefix("0x")
        .map_or((value, 10), |digits| (digits, 16));
    BigUint::parse_bytes(digits.as_bytes(), radix)
        .ok_or_else(|| format!("invalid integer `{value}`"))
}

fn lower_hex(value: &BigUint) -> String {
    format!("0x{}", value.to_str_radix(16))
}

fn checked_dependency(path: &Path, expected: &str) -> Result<(Dependency, Vec<u8>), String> {
    let bytes = fs::read(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let digest = sha256_hex(&bytes);
    if digest != expected {
        return Err(format!(
            "{} hash mismatch: expected {expected}, got {digest}",
            path.display()
        ));
    }
    Ok((
        Dependency {
            path: path.display().to_string(),
            bytes: bytes.len() as u64,
            sha256: digest,
        },
        bytes,
    ))
}

fn big_poly_from_roots(roots: &[BigUint], p: &BigUint) -> Vec<BigUint> {
    roots.iter().fold(vec![BigUint::one()], |poly, root| {
        let neg = if root.is_zero() {
            BigUint::zero()
        } else {
            p - root
        };
        let mut out = vec![BigUint::zero(); poly.len() + 1];
        for (i, coefficient) in poly.iter().enumerate() {
            out[i] = (&out[i] + coefficient * &neg) % p;
            out[i + 1] = (&out[i + 1] + coefficient) % p;
        }
        out
    })
}

fn big_poly_eval(poly: &[BigUint], x: &BigUint, p: &BigUint) -> BigUint {
    poly.iter()
        .rev()
        .fold(BigUint::zero(), |value, coefficient| {
            (value * x + coefficient) % p
        })
}

fn big_derivative(poly: &[BigUint], p: &BigUint) -> Vec<BigUint> {
    poly.iter()
        .enumerate()
        .skip(1)
        .map(|(i, coefficient)| (coefficient * BigUint::from(i)) % p)
        .collect()
}

fn coefficient_stream(poly: &[BigUint]) -> Result<Vec<u8>, String> {
    let mut out = Vec::with_capacity(poly.len() * 32);
    for coefficient in poly {
        let raw = coefficient.to_bytes_be();
        if raw.len() > 32 {
            return Err("P-256 coefficient exceeded 32 bytes".into());
        }
        out.extend(std::iter::repeat_n(0, 32 - raw.len()));
        out.extend(raw);
    }
    Ok(out)
}

fn dump_columns(points: &[WideFbPoint]) -> Result<Vec<BigUint>, String> {
    if !points.len().is_multiple_of(2) {
        return Err("factor base has an odd signed-point count".into());
    }
    let curve = CurveParams::p256();
    let mut out = Vec::with_capacity(points.len() / 2);
    for (column, pair) in points.chunks_exact(2).enumerate() {
        if pair[0].col != column as u64
            || pair[1].col != column as u64
            || pair[0].x != pair[1].x
            || pair[0].coef != "1"
            || pair[1].coef != (&curve.n - BigUint::one()).to_string()
        {
            return Err(format!("invalid signed pair at column {column}"));
        }
        let x = parse_big(&pair[0].x)?;
        let y0 = parse_big(&pair[0].y)?;
        let y1 = parse_big(&pair[1].y)?;
        let rhs = (((&x * &x) % &curve.p) * &x + &curve.a * &x + &curve.b) % &curve.p;
        if (&y0 + &y1) % &curve.p != BigUint::zero()
            || (&y0 * &y0) % &curve.p != rhs
            || (&y1 * &y1) % &curve.p != rhs
        {
            return Err(format!("invalid P-256 point pair at column {column}"));
        }
        out.push(x);
    }
    if !out.windows(2).all(|pair| pair[0] < pair[1]) {
        return Err("factor-base columns are not strictly sorted".into());
    }
    Ok(out)
}

fn build_component(
    depth: u32,
    root: &str,
) -> Result<p256_dickson_factor_base::WideBuildResult, String> {
    let built = p256_dickson_factor_base::build(&format!(
        "dickson-torus:depth={depth},root_exponent={root}"
    ))?;
    let replay = p256_dickson_factor_base::verify(&built.dump)?;
    if replay != built {
        return Err(format!("depth-{depth} component replay mismatch"));
    }
    Ok(built)
}

fn component_x(built: &p256_dickson_factor_base::WideBuildResult) -> Result<Vec<BigUint>, String> {
    dump_columns(&built.dump.points)
}

fn dickson_chain(x: &BigUint, depth: usize, p: &BigUint) -> Vec<BigUint> {
    let mut value = x.clone();
    let mut out = Vec::with_capacity(depth);
    for _ in 0..depth {
        value = if (&value * &value) % p >= BigUint::from(2u8) {
            ((&value * &value) % p) - BigUint::from(2u8)
        } else {
            ((&value * &value) % p) + p - BigUint::from(2u8)
        };
        out.push(value.clone());
    }
    out
}

fn set_union(left: &[BigUint], right: &[BigUint]) -> BTreeSet<BigUint> {
    left.iter().chain(right).cloned().collect()
}

fn binomial(n: usize, k: usize) -> BigUint {
    if k > n {
        return BigUint::zero();
    }
    let k = k.min(n - k);
    let mut value = BigUint::one();
    for i in 0..k {
        value *= BigUint::from(n - i);
        value /= BigUint::from(i + 1);
    }
    value
}

fn p256_input_monomials(support: &[BigUint], m: usize) -> usize {
    let columns = support.len() - 1;
    let complement = columns - m;
    let d = m.div_ceil(2);
    let mut total = 0usize;
    for k in 0..columns {
        for i in 0..=m {
            if k >= i {
                let j = k - i;
                if j <= complement {
                    total += 1;
                }
            }
        }
        if !support[k].is_zero() {
            total += 1;
        }
    }
    for k in 0..=m {
        for i in 0..=d {
            let Some(j) = k.checked_sub(i) else {
                continue;
            };
            if i <= j && j <= d && !(i == d && j == d) {
                total += 1;
            }
        }
        for shift in [0usize, 1, 3] {
            let Some(sum) = k.checked_sub(shift) else {
                continue;
            };
            for i in 0..d.saturating_sub(1) {
                let Some(j) = sum.checked_sub(i) else {
                    continue;
                };
                if i <= j && j < d - 1 {
                    total += 1;
                }
            }
        }
        if k >= 1 && k - 1 <= m {
            total += 1;
        }
        if k <= m {
            total += 1;
        }
    }
    total
}

fn selector_shape(columns: usize, m: usize, support: &[BigUint]) -> SelectorShape {
    let normalized_variables = columns + m;
    let normalized_equations = columns + m + 1;
    let mut macaulay = Vec::new();
    for degree in 2..=7usize {
        let columns_at_degree = binomial(normalized_variables + degree, degree);
        let rows = BigUint::from(normalized_equations)
            * binomial(normalized_variables + degree - 2, degree - 2);
        let bytes = &columns_at_degree * &rows * BigUint::from(32u8);
        macaulay.push(MacaulayDimension {
            degree: degree as u32,
            columns_at_most_degree: columns_at_degree.to_string(),
            generic_rows_at_most_degree: rows.to_string(),
            dense_bytes_at_32_bytes_per_entry: bytes.to_string(),
        });
    }
    SelectorShape {
        arity: m,
        columns,
        complement_degree: columns - m,
        projective_variables: columns + m + 2,
        projective_equations: columns + m + 2,
        normalized_variables,
        normalized_equations,
        maximum_input_degree: 2,
        input_monomials: p256_input_monomials(support, m),
        macaulay,
    }
}

fn p256_certificate(dump: &WideFactorBaseDump) -> Result<P256Certificate, String> {
    if dump.schema != "ecbench.factor_base_dump/v1-wide"
        || dump.curve.slug != CURVE_SLUG
        || dump.factor_base.fb_id != FB_ID
        || dump.factor_base.columns != 164
        || dump.factor_base.signed_points != 328
    {
        return Err("Round-295 factor-base identity mismatch".into());
    }
    let params = &dump.factor_base.params;
    let root8 = params
        .get("component_0_root_exponent")
        .ok_or("missing depth-8 root")?;
    let root6 = params
        .get("component_1_root_exponent")
        .ok_or("missing depth-6 root")?;
    let terminal8 = parse_big(
        params
            .get("component_0_terminal")
            .ok_or("missing depth-8 terminal")?,
    )?;
    let terminal6 = parse_big(
        params
            .get("component_1_terminal")
            .ok_or("missing depth-6 terminal")?,
    )?;
    let depth8 = build_component(8, root8)?;
    let depth6 = build_component(6, root6)?;
    if depth8.dump.factor_base.fb_id != params["component_0_fb_id"]
        || depth6.dump.factor_base.fb_id != params["component_1_fb_id"]
    {
        return Err("component FB1 mismatch".into());
    }
    let x8 = component_x(&depth8)?;
    let x6 = component_x(&depth6)?;
    let selected = dump_columns(&dump.points)?;
    let exact_union = set_union(&x8, &x6);
    if exact_union.iter().cloned().collect::<Vec<_>>() != selected {
        return Err("rebuilt component union differs from selected factor base".into());
    }
    let intersection = x8.iter().filter(|x| x6.binary_search(x).is_ok()).count();
    let curve = CurveParams::p256();
    let mut common_rows = 0usize;
    for x in &selected {
        let chain = dickson_chain(x, 8, &curve.p);
        if chain[7] != terminal8 && chain[5] != terminal6 {
            return Err("stored column failed the common-chain terminal product".into());
        }
        common_rows += 1;
    }

    let overbroad = build_component(8, root6)?;
    let overbroad_x = component_x(&overbroad)?;
    let overbroad_union = set_union(&x8, &overbroad_x);
    let overbroad_extra = overbroad_union.difference(&exact_union).count();
    if overbroad_extra == 0 || !x6.iter().all(|x| overbroad_x.binary_search(x).is_ok()) {
        return Err(
            "overbroad depth-8 control did not contain the depth-6 component plus extras".into(),
        );
    }

    let support = big_poly_from_roots(&selected, &curve.p);
    let support_replay = big_poly_from_roots(&selected, &curve.p);
    let stream = coefficient_stream(&support)?;
    let replay_stream = coefficient_stream(&support_replay)?;
    let derivative = big_derivative(&support, &curve.p);
    let roots_checked = selected
        .iter()
        .filter(|x| big_poly_eval(&support, x, &curve.p).is_zero())
        .count();
    let derivative_nonzero = selected
        .iter()
        .filter(|x| !big_poly_eval(&derivative, x, &curve.p).is_zero())
        .count();
    if roots_checked != selected.len() || derivative_nonzero != selected.len() {
        return Err("support polynomial root or square-free check failed".into());
    }

    let mut scalar = BigUint::from_bytes_be(&sha256(
        format!("{CURVE_SLUG}/rr-selector-round296/public").as_bytes(),
    )) % &curve.n;
    if scalar.is_zero() {
        scalar = BigUint::one();
    }
    let target = curve.generator().scalar_mul_vartime(&scalar, &curve.a_fe());
    let (target_x, target_y) = match target {
        Point::Affine { x, y } => (x.value, y.value),
        Point::Infinity => return Err("P-256 public target is infinity".into()),
    };
    if selected.binary_search(&target_x).is_ok() {
        return Err("P-256 public target abscissa is in the factor base".into());
    }

    Ok(P256Certificate {
        fb_id: FB_ID.into(),
        columns: selected.len(),
        signed_points: dump.points.len(),
        depth8_columns: x8.len(),
        depth6_columns: x6.len(),
        component_overlap: intersection,
        common_chain_rows_checked: common_rows,
        common_chain_exact: exact_union.len() == selected.len(),
        overbroad_depth8_columns: overbroad_x.len(),
        overbroad_union_columns: overbroad_union.len(),
        overbroad_extra_columns: overbroad_extra,
        support_degree: support.len() - 1,
        support_coefficients_sha256: sha256_hex(&stream),
        support_roots_checked: roots_checked,
        support_derivatives_nonzero: derivative_nonzero,
        support_replay_byte_identical: stream == replay_stream,
        public_scalar: lower_hex(&scalar),
        public_target: [lower_hex(&target_x), lower_hex(&target_y)],
        selector: selector_shape(selected.len(), 109, &support),
    })
}

fn toy_rhs(x: u64, b: u64) -> u64 {
    addm(addm(mulm(mulm(x, x), x), mulm(A, x)), b)
}

fn square_roots() -> Vec<Vec<u64>> {
    let mut table = vec![Vec::new(); P as usize];
    for value in 0..P {
        table[mulm(value, value) as usize].push(value);
    }
    table
}

fn toy_iterate(mut value: u64, depth: u32) -> u64 {
    for _ in 0..depth {
        value = subm(mulm(value, value), 2);
    }
    value
}

fn toy_support_pool(b: u64, square_roots: &[Vec<u64>]) -> Vec<u64> {
    let mut values = (0..P)
        .filter(|&x| {
            matches!(toy_iterate(x, 5), 0 | 369) && !square_roots[toy_rhs(x, b) as usize].is_empty()
        })
        .collect::<Vec<_>>();
    values.sort_by(|left, right| {
        let left_hash = sha256(format!("rr-selector-round296/support/{left}").as_bytes());
        let right_hash = sha256(format!("rr-selector-round296/support/{right}").as_bytes());
        left_hash.cmp(&right_hash).then(left.cmp(right))
    });
    values
}

fn add_points(left: ToyPoint, right: ToyPoint) -> ToyPoint {
    let (left, right) = match (left, right) {
        (ToyPoint::Infinity, value) | (value, ToyPoint::Infinity) => return value,
        (ToyPoint::Affine(left), ToyPoint::Affine(right)) => (left, right),
    };
    let slope = if left.x == right.x {
        if addm(left.y, right.y) == 0 || left.y == 0 {
            return ToyPoint::Infinity;
        }
        mulm(
            addm(mulm(3, mulm(left.x, left.x)), A),
            invm(mulm(2, left.y)),
        )
    } else {
        mulm(subm(right.y, left.y), invm(subm(right.x, left.x)))
    };
    let x = subm(subm(mulm(slope, slope), left.x), right.x);
    let y = subm(mulm(slope, subm(left.x, x)), left.y);
    ToyPoint::Affine(Affine { x, y })
}

fn negate(point: Affine) -> Affine {
    Affine {
        x: point.x,
        y: negm(point.y),
    }
}

fn sum_points(points: impl IntoIterator<Item = Affine>) -> ToyPoint {
    points.into_iter().fold(ToyPoint::Infinity, |sum, point| {
        add_points(sum, ToyPoint::Affine(point))
    })
}

fn low_point(x: u64, b: u64, square_roots: &[Vec<u64>]) -> Result<Affine, String> {
    let roots = &square_roots[toy_rhs(x, b) as usize];
    let y = roots
        .iter()
        .copied()
        .min()
        .ok_or("support x is not liftable")?;
    Ok(Affine { x, y })
}

fn next_combination(combination: &mut [usize], limit: usize) -> bool {
    let width = combination.len();
    for position in (0..width).rev() {
        if combination[position] < limit - width + position {
            combination[position] += 1;
            for next in position + 1..width {
                combination[next] = combination[next - 1] + 1;
            }
            return true;
        }
    }
    false
}

fn combination_list(width: usize, choose: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    let mut combination = (0..choose).collect::<Vec<_>>();
    loop {
        out.push(combination.clone());
        if !next_combination(&mut combination, width) {
            break;
        }
    }
    out
}

fn planted_target(
    support: &[u64],
    m: usize,
    b: u64,
    square_roots: &[Vec<u64>],
) -> Result<(Affine, Vec<usize>, u64), String> {
    let support_set = support.iter().copied().collect::<BTreeSet<_>>();
    let mut combinations = combination_list(support.len(), m);
    combinations.sort_by(|left, right| {
        let left_key = format!(
            "rr-selector-round296/planted/{m}/{}/{}",
            support.len(),
            left.iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(",")
        );
        let right_key = format!(
            "rr-selector-round296/planted/{m}/{}/{}",
            support.len(),
            right
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(",")
        );
        sha256(left_key.as_bytes())
            .cmp(&sha256(right_key.as_bytes()))
            .then(left.cmp(right))
    });
    for (rank, columns) in combinations.into_iter().enumerate() {
        let sign_hash = sha256(
            format!(
                "rr-selector-round296/planted/{m}/{}/signs/{rank}",
                support.len()
            )
            .as_bytes(),
        );
        let mut mask = 0u64;
        let mut points = Vec::with_capacity(m);
        for (position, &column) in columns.iter().enumerate() {
            let mut point = low_point(support[column], b, square_roots)?;
            if sign_hash[position / 8] >> (position % 8) & 1 == 1 {
                point = negate(point);
                mask |= 1 << position;
            }
            points.push(point);
        }
        if let ToyPoint::Affine(target) = sum_points(points) {
            if !support_set.contains(&target.x) {
                return Ok((target, columns, mask));
            }
        }
    }
    Err(format!("no admissible planted target for m={m}"))
}

fn unplanted_target(
    support: &[u64],
    m: usize,
    b: u64,
    square_roots: &[Vec<u64>],
) -> Result<Affine, String> {
    let support_set = support.iter().copied().collect::<BTreeSet<_>>();
    let mut points = (0..P)
        .filter(|x| !support_set.contains(x))
        .flat_map(|x| {
            square_roots[toy_rhs(x, b) as usize]
                .iter()
                .copied()
                .map(move |y| Affine { x, y })
        })
        .collect::<Vec<_>>();
    points.sort_by(|left, right| {
        let left_hash = sha256(
            format!(
                "rr-selector-round296/unplanted/{m}/{}/{},{}",
                support.len(),
                left.x,
                left.y
            )
            .as_bytes(),
        );
        let right_hash = sha256(
            format!(
                "rr-selector-round296/unplanted/{m}/{}/{},{}",
                support.len(),
                right.x,
                right.y
            )
            .as_bytes(),
        );
        left_hash.cmp(&right_hash).then(left.cmp(right))
    });
    points
        .into_iter()
        .next()
        .ok_or_else(|| "toy curve has no target".to_string())
}

#[derive(Clone, Debug)]
struct Reference {
    candidates: u64,
    witnesses: Vec<(Vec<usize>, u64)>,
    digest: String,
    replay_failures: u64,
}

fn reference_relations(
    support: &[u64],
    m: usize,
    target: Affine,
    b: u64,
    square_roots: &[Vec<u64>],
) -> Result<Reference, String> {
    let base_points = support
        .iter()
        .map(|&x| low_point(x, b, square_roots))
        .collect::<Result<Vec<_>, _>>()?;
    let negative_target = negate(target);
    let mut candidates = 0u64;
    let mut witnesses = Vec::new();
    let mut stream = String::new();
    let mut combination = (0..m).collect::<Vec<_>>();
    loop {
        for mask in 0..1u64 << m {
            candidates += 1;
            let points = combination.iter().enumerate().map(|(position, &column)| {
                let point = base_points[column];
                if mask >> position & 1 == 1 {
                    negate(point)
                } else {
                    point
                }
            });
            if matches!(sum_points(points), ToyPoint::Affine(sum) if sum == target || sum == negative_target)
            {
                use std::fmt::Write as _;
                writeln!(
                    stream,
                    "{}:{mask}",
                    combination
                        .iter()
                        .map(usize::to_string)
                        .collect::<Vec<_>>()
                        .join(",")
                )
                .map_err(|error| error.to_string())?;
                witnesses.push((combination.clone(), mask));
            }
        }
        if !next_combination(&mut combination, support.len()) {
            break;
        }
    }
    let replay_failures = witnesses
        .iter()
        .filter(|(columns, mask)| {
            let points = columns.iter().enumerate().map(|(position, &column)| {
                let point = base_points[column];
                if mask >> position & 1 == 1 {
                    negate(point)
                } else {
                    point
                }
            });
            !matches!(sum_points(points), ToyPoint::Affine(sum) if sum == target || sum == negative_target)
        })
        .count() as u64;
    Ok(Reference {
        candidates,
        witnesses,
        digest: sha256_hex(stream.as_bytes()),
        replay_failures,
    })
}

fn null_vector(mut matrix: Vec<Vec<u64>>) -> Option<Vec<u64>> {
    let rows = matrix.len();
    let columns = matrix.first()?.len();
    let mut pivots = Vec::new();
    let mut row = 0usize;
    for column in 0..columns {
        let Some(pivot) = (row..rows).find(|&candidate| matrix[candidate][column] != 0) else {
            continue;
        };
        matrix.swap(row, pivot);
        let inverse = invm(matrix[row][column]);
        for value in &mut matrix[row] {
            *value = mulm(*value, inverse);
        }
        let pivot_row = matrix[row].clone();
        for other in 0..rows {
            if other == row || matrix[other][column] == 0 {
                continue;
            }
            let factor = matrix[other][column];
            for (value, &pivot_value) in matrix[other].iter_mut().zip(&pivot_row) {
                *value = subm(*value, mulm(factor, pivot_value));
            }
        }
        pivots.push((row, column));
        row += 1;
        if row == rows {
            break;
        }
    }
    let pivot_columns = pivots
        .iter()
        .map(|(_, column)| *column)
        .collect::<BTreeSet<_>>();
    let free = (0..columns).find(|column| !pivot_columns.contains(column))?;
    let mut vector = vec![0; columns];
    vector[free] = 1;
    for &(pivot_row, pivot_column) in pivots.iter().rev() {
        let sum = ((pivot_column + 1)..columns).fold(0, |value, column| {
            addm(value, mulm(matrix[pivot_row][column], vector[column]))
        });
        vector[pivot_column] = negm(sum);
    }
    matrix
        .iter()
        .all(|row| {
            row.iter()
                .zip(&vector)
                .fold(0, |sum, (&a, &b)| addm(sum, mulm(a, b)))
                == 0
        })
        .then_some(vector)
}

#[derive(Clone, Debug)]
struct EquationWitness {
    values: Vec<u64>,
    equations_verified: bool,
}

fn rr_witness(
    support: &[u64],
    columns: &[usize],
    negative_mask: u64,
    target: Affine,
    b: u64,
    square_roots: &[Vec<u64>],
) -> Result<EquationWitness, String> {
    let m = columns.len();
    let d = m.div_ceil(2);
    let mut points = Vec::with_capacity(m);
    for (position, &column) in columns.iter().enumerate() {
        let mut point = low_point(support[column], b, square_roots)?;
        if negative_mask >> position & 1 == 1 {
            point = negate(point);
        }
        points.push(point);
    }
    if sum_points(points.iter().copied()) != ToyPoint::Affine(target) {
        return Err("planted points do not sum to the planted target".into());
    }
    points.push(negate(target));
    let matrix = points
        .iter()
        .map(|point| {
            let mut row = Vec::with_capacity(m + 1);
            let mut power = 1u64;
            for _ in 0..=d {
                row.push(power);
                power = mulm(power, point.x);
            }
            power = 1;
            for _ in 0..d - 1 {
                row.push(mulm(point.y, power));
                power = mulm(power, point.x);
            }
            row
        })
        .collect::<Vec<_>>();
    let mut function = null_vector(matrix).ok_or("planted RR interpolation has no null vector")?;
    let leading = function[d];
    if leading == 0 {
        return Err("planted RR interpolation has zero leading A coefficient".into());
    }
    let scale = invm(leading);
    for value in &mut function {
        *value = mulm(*value, scale);
    }
    let a_poly = function[..=d].to_vec();
    let b_poly = function[d + 1..].to_vec();
    let selected_x = columns
        .iter()
        .map(|&column| support[column])
        .collect::<Vec<_>>();
    let selected_set = columns.iter().copied().collect::<BTreeSet<_>>();
    let complement_x = support
        .iter()
        .enumerate()
        .filter(|(column, _)| !selected_set.contains(column))
        .map(|(_, &x)| x)
        .collect::<Vec<_>>();
    let g = product_roots_u64(&selected_x);
    let q = product_roots_u64(&complement_x);
    let support_poly = product_roots_u64(support);
    if poly_mul_u64(&g, &q) != support_poly {
        return Err("planted support factorization mismatch".into());
    }
    let a_squared = poly_mul_u64(&a_poly, &a_poly);
    let b_squared = poly_mul_u64(&b_poly, &b_poly);
    let curve_b_squared = poly_mul_u64(&[b, A, 0, 1], &b_squared);
    let norm = poly_sub_u64(&a_squared, &curve_b_squared);
    let rhs = poly_mul_u64(&[negm(target.x), 1], &g);
    if norm != rhs {
        return Err("planted normalized norm identity mismatch".into());
    }
    let mut values = Vec::with_capacity(support.len() + m);
    values.extend_from_slice(&g[..m]);
    values.extend_from_slice(&q[..support.len() - m]);
    values.extend_from_slice(&a_poly[..d]);
    values.extend_from_slice(&b_poly);
    let (equations, variables) = rr_system(support, m, target.x, b)?;
    if values.len() != variables {
        return Err("planted witness variable count mismatch".into());
    }
    let equations_verified = equations
        .iter()
        .all(|equation| equation.evaluate(&values) == 0);
    Ok(EquationWitness {
        values,
        equations_verified,
    })
}

fn poly_convolution(left: &[Pol], right: &[Pol]) -> Vec<Pol> {
    let n = left[0].n;
    let mut out = vec![Pol::zero(n); left.len() + right.len() - 1];
    for (i, a) in left.iter().enumerate() {
        for (j, b) in right.iter().enumerate() {
            out[i + j] = out[i + j].add(&a.mul(b));
        }
    }
    out
}

fn rr_system(
    support: &[u64],
    m: usize,
    target_x: u64,
    curve_b: u64,
) -> Result<(Vec<Pol>, usize), String> {
    let columns = support.len();
    if m % 2 != 1 || m >= columns {
        return Err("RR selector requires odd m below the support width".into());
    }
    let complement = columns - m;
    let d = m.div_ceil(2);
    let variables = columns + m;
    let mut next = 0usize;
    let mut g = (0..m)
        .map(|_| {
            let value = Pol::var(variables, next);
            next += 1;
            value
        })
        .collect::<Vec<_>>();
    g.push(Pol::constant(variables, 1));
    let mut q = (0..complement)
        .map(|_| {
            let value = Pol::var(variables, next);
            next += 1;
            value
        })
        .collect::<Vec<_>>();
    q.push(Pol::constant(variables, 1));
    let mut a_poly = (0..d)
        .map(|_| {
            let value = Pol::var(variables, next);
            next += 1;
            value
        })
        .collect::<Vec<_>>();
    a_poly.push(Pol::constant(variables, 1));
    let b_poly = (0..d - 1)
        .map(|_| {
            let value = Pol::var(variables, next);
            next += 1;
            value
        })
        .collect::<Vec<_>>();
    if next != variables {
        return Err("RR selector variable layout mismatch".into());
    }

    let support_poly = product_roots_u64(support);
    let factor = poly_convolution(&g, &q);
    let mut equations = (0..columns)
        .map(|coefficient| {
            factor[coefficient].sub(&Pol::constant(variables, support_poly[coefficient]))
        })
        .collect::<Vec<_>>();
    let a_squared = poly_convolution(&a_poly, &a_poly);
    let b_squared = poly_convolution(&b_poly, &b_poly);
    let curve = [
        Pol::constant(variables, curve_b),
        Pol::constant(variables, A),
        Pol::zero(variables),
        Pol::constant(variables, 1),
    ];
    let curve_b_squared = poly_convolution(&curve, &b_squared);
    let target = [
        Pol::constant(variables, negm(target_x)),
        Pol::constant(variables, 1),
    ];
    let rhs = poly_convolution(&target, &g);
    for coefficient in 0..=m {
        equations.push(
            a_squared[coefficient]
                .sub(&curve_b_squared[coefficient])
                .sub(&rhs[coefficient]),
        );
    }
    if equations.len() != columns + m + 1 || equations.iter().any(|equation| equation.degree() > 2)
    {
        return Err("RR selector equation layout mismatch".into());
    }
    Ok((equations, variables))
}

fn target_receipt(
    kind: &str,
    support: &[u64],
    m: usize,
    target: Affine,
    planted: Option<(&[usize], u64)>,
    b: u64,
    square_roots: &[Vec<u64>],
) -> Result<(TargetReceipt, f64), String> {
    let reference = reference_relations(support, m, target, b, square_roots)?;
    let witness = planted
        .map(|(columns, mask)| rr_witness(support, columns, mask, target, b, square_roots))
        .transpose()?;
    let (equations, variables) = rr_system(support, m, target.x, b)?;
    if let Some(witness) = &witness {
        if witness.values.len() != variables || !witness.equations_verified {
            return Err("planted coefficient witness failed".into());
        }
    }
    let input_monomials = equations.iter().map(|equation| equation.terms.len()).sum();
    let input = equations.iter().map(Pol::to_f4).collect::<Vec<_>>();
    let started = Instant::now();
    let report = f4_fp::f4(
        &input,
        variables,
        P,
        &F4Options::new(Ordering::Grevlex, MAX_DEGREE)
            .with_budget(Duration::from_secs(BUDGET_SECS)),
    );
    let wall_ms = started.elapsed().as_secs_f64() * 1000.0;
    let complete =
        !report.timed_out && report.pairs_above_bound == 0 && report.staircase_at_stop.is_none();
    let algebraic_classification = complete.then_some(!report.inconsistent);
    let expected = !reference.witnesses.is_empty();
    let false_positive = algebraic_classification == Some(true) && !expected;
    let false_negative = algebraic_classification == Some(false) && expected;
    Ok((
        TargetReceipt {
            kind: kind.into(),
            target: [target.x, target.y],
            planted_columns: planted.map(|(columns, _)| columns.to_vec()),
            planted_negative_mask: planted.map(|(_, mask)| mask),
            candidates: reference.candidates,
            relations: reference.witnesses.len() as u64,
            replayed_relations: reference.witnesses.len() as u64,
            replay_failures: reference.replay_failures,
            witness_sha256: reference.digest,
            planted_equation_witness_verified: witness.map(|value| value.equations_verified),
            variables,
            equations: equations.len(),
            input_monomials,
            input_maximum_degree: equations.iter().map(Pol::degree).max().unwrap_or(0),
            f4_inconsistent: report.inconsistent,
            f4_classification: algebraic_classification,
            f4_complete: complete,
            f4_timed_out: report.timed_out,
            f4_pairs_above_bound: report.pairs_above_bound,
            solving_degree_max: report.solving_degree_max,
            max_rows: report.max_rows,
            max_cols: report.max_cols,
            max_cols_to_solution: report.max_cols_to_solution,
            field_operations: report.field_ops,
            false_positive,
            false_negative,
        },
        wall_ms,
    ))
}

fn fit_slope(points: &[(f64, f64)], logarithmic: bool) -> Option<f64> {
    if points.len() < 3 || points.iter().any(|&(x, y)| x <= 0.0 || y <= 0.0) {
        return None;
    }
    let transformed = points
        .iter()
        .map(|&(x, y)| {
            if logarithmic {
                (x.ln(), y.ln())
            } else {
                (x, y)
            }
        })
        .collect::<Vec<_>>();
    let mean_x = transformed.iter().map(|(x, _)| x).sum::<f64>() / transformed.len() as f64;
    let mean_y = transformed.iter().map(|(_, y)| y).sum::<f64>() / transformed.len() as f64;
    let numerator = transformed
        .iter()
        .map(|(x, y)| (x - mean_x) * (y - mean_y))
        .sum::<f64>();
    let denominator = transformed
        .iter()
        .map(|(x, _)| (x - mean_x).powi(2))
        .sum::<f64>();
    (denominator > 0.0).then_some(numerator / denominator)
}

fn stream_u64(values: &[u64]) -> Vec<u8> {
    values
        .iter()
        .flat_map(|value| value.to_be_bytes())
        .collect()
}

fn write_result(mut result: ResultReceipt, path: &Path) -> Result<(), String> {
    let semantic = serde_json::json!({
        "schema": result.schema,
        "curve": result.curve,
        "p256": result.p256,
        "cases": result.cases,
        "gates": result.gates,
        "classification": result.classification,
        "dominant_obstruction": result.dominant_obstruction,
        "decision": result.decision,
    });
    result.semantic_evidence_sha256 =
        sha256_hex(&serde_json::to_vec(&semantic).map_err(|error| error.to_string())?);
    loop {
        let text = serde_json::to_string_pretty(&result).map_err(|error| error.to_string())? + "\n";
        let bytes = text.len() as u64;
        if bytes == result.result_json_bytes {
            fs::write(path, text).map_err(|error| format!("{}: {error}", path.display()))?;
            break;
        }
        result.result_json_bytes = bytes;
    }
    Ok(())
}

fn run(cli: Cli) -> Result<(), String> {
    let (round295, _) = checked_dependency(&cli.round295, ROUND295_SHA256)?;
    let (factor_base, dump_bytes) = checked_dependency(&cli.factor_base, FACTOR_BASE_SHA256)?;
    let dump: WideFactorBaseDump =
        serde_json::from_slice(&dump_bytes).map_err(|error| error.to_string())?;
    let p256 = p256_certificate(&dump)?;

    let curve = CurveParams::p256();
    let toy_b = (&curve.b % P).to_u64().ok_or("P-256 b reduction failed")?;
    let roots = square_roots();
    let pool = toy_support_pool(toy_b, &roots);
    if pool.len() < TOY_SIZES.last().expect("nonempty sizes").1 {
        return Err(format!("toy support pool has only {} columns", pool.len()));
    }
    let pool_sha256 = sha256_hex(&stream_u64(&pool));
    let mut cases = Vec::new();
    let mut telemetry_targets = Vec::new();
    for (m, columns) in TOY_SIZES {
        let support = pool[..columns].to_vec();
        let support_sha256 = sha256_hex(&stream_u64(&support));
        let (planted_point, planted_columns, planted_mask) =
            planted_target(&support, m, toy_b, &roots)?;
        let (planted, planted_ms) = target_receipt(
            "planted",
            &support,
            m,
            planted_point,
            Some((&planted_columns, planted_mask)),
            toy_b,
            &roots,
        )?;
        telemetry_targets.push(TargetTelemetry {
            arity: m,
            columns,
            kind: "planted".into(),
            wall_ms: planted_ms,
        });
        let unplanted_point = unplanted_target(&support, m, toy_b, &roots)?;
        let (unplanted, unplanted_ms) = target_receipt(
            "unplanted",
            &support,
            m,
            unplanted_point,
            None,
            toy_b,
            &roots,
        )?;
        telemetry_targets.push(TargetTelemetry {
            arity: m,
            columns,
            kind: "unplanted".into(),
            wall_ms: unplanted_ms,
        });
        let signed_domain = binomial(columns, m)
            .to_u64()
            .ok_or("toy binomial does not fit u64")?
            * (1u64 << m);
        if planted.candidates != signed_domain || unplanted.candidates != signed_domain {
            return Err("exhaustive reference candidate count mismatch".into());
        }
        cases.push(ToyCase {
            arity: m,
            columns,
            support_x: support,
            support_sha256,
            signed_domain,
            planted,
            unplanted,
        });
    }

    let completed = cases
        .iter()
        .filter(|case| case.planted.f4_complete && case.unplanted.f4_complete)
        .collect::<Vec<_>>();
    let complete_targets = cases
        .iter()
        .flat_map(|case| [&case.planted, &case.unplanted])
        .filter(|target| target.f4_complete)
        .count();
    let censored_targets = 2 * cases.len() - complete_targets;
    let false_positives = cases
        .iter()
        .flat_map(|case| [&case.planted, &case.unplanted])
        .filter(|target| target.false_positive)
        .count() as u64;
    let false_negatives = cases
        .iter()
        .flat_map(|case| [&case.planted, &case.unplanted])
        .filter(|target| target.false_negative)
        .count() as u64;
    let degree_points = completed
        .iter()
        .map(|case| {
            (
                case.planted.variables as f64,
                case.planted
                    .solving_degree_max
                    .max(case.unplanted.solving_degree_max) as f64,
            )
        })
        .collect::<Vec<_>>();
    let field_points = completed
        .iter()
        .map(|case| {
            (
                case.planted.variables as f64,
                case.planted
                    .field_operations
                    .max(case.unplanted.field_operations) as f64,
            )
        })
        .collect::<Vec<_>>();
    let width_points = completed
        .iter()
        .map(|case| {
            (
                case.planted.variables as f64,
                case.planted.max_cols.max(case.unplanted.max_cols) as f64,
            )
        })
        .collect::<Vec<_>>();
    let time_points = completed
        .iter()
        .filter_map(|case| {
            let times = telemetry_targets
                .iter()
                .filter(|target| target.arity == case.arity && target.columns == case.columns)
                .map(|target| target.wall_ms)
                .collect::<Vec<_>>();
            (!times.is_empty()).then_some((
                case.planted.variables as f64,
                times.into_iter().fold(0.0, f64::max),
            ))
        })
        .collect::<Vec<_>>();
    let degree_fit_slope = fit_slope(&degree_points, false);
    let field_ops_fit = fit_slope(&field_points, true);
    let width_fit = fit_slope(&width_points, true);
    let time_fit = fit_slope(&time_points, true);
    let largest_three_available = completed.len() >= 3;
    let largest_three = completed.iter().rev().take(3).copied().collect::<Vec<_>>();
    let degree_at_most_5 = largest_three_available
        && largest_three.iter().all(|case| {
            case.planted.solving_degree_max <= 5 && case.unplanted.solving_degree_max <= 5
        });
    let non_increasing_degree =
        largest_three_available && degree_fit_slope.is_some_and(|slope| slope <= 0.0);
    let planted_exact = cases.iter().all(|case| {
        case.planted.planted_equation_witness_verified == Some(true)
            && case.planted.replay_failures == 0
    });
    let no_classification_errors = false_positives == 0
        && false_negatives == 0
        && cases
            .iter()
            .flat_map(|case| [&case.planted, &case.unplanted])
            .filter(|target| target.f4_complete)
            .all(|target| target.replay_failures == 0);
    let gates = Gates {
        common_chain_exact: p256.common_chain_exact,
        planted_witnesses_exact: planted_exact,
        zero_false_positives_and_false_negatives_on_complete_instances: no_classification_errors,
        largest_three_complete_sizes_available: largest_three_available,
        structured_residual_degree_at_most_5: degree_at_most_5,
        non_increasing_degree_slope: non_increasing_degree,
        projected_materialized_storage_below_2_50: false,
        projected_cost_per_usable_relation_below_2_103: false,
        projected_collection_below_2_120: false,
        promoted: false,
    };
    let dominant_obstruction = if !largest_three_available {
        "the normalized quadratic selector did not complete on three sizes; censored branches cannot support a P-256 projection"
    } else if !degree_at_most_5 || !non_increasing_degree {
        "global solving degree grows despite quadratic Dickson and norm generators"
    } else {
        "toy degree alone does not price a 273-variable structured solver below the storage and cost gates"
    };
    let result = ResultReceipt {
        schema: "p256.dickson_rr_selector/v1".into(),
        curve: CURVE_SLUG.into(),
        screening_round: 296,
        execution_status: "complete".into(),
        round295,
        factor_base,
        p256,
        toy_prime: P,
        toy_a: A,
        toy_b,
        toy_support_pool: pool,
        toy_support_pool_sha256: pool_sha256,
        cases,
        complete_targets,
        censored_targets,
        false_positives,
        false_negatives,
        degree_fit_slope,
        field_ops_fit_exponent_in_variables: field_ops_fit,
        width_fit_exponent_in_variables: width_fit,
        p256_relation_attempted: false,
        p256_relations_reported: 0,
        gates,
        classification: "structured-selector experiment; promotion rejected".into(),
        dominant_obstruction: dominant_obstruction.into(),
        decision: "retain the exact common-chain and normalized RR formulation, but do not promote or attempt an unplanted P-256 relation unless a specialized solver clears the unset global storage and cost gates".into(),
        semantic_evidence_sha256: String::new(),
        result_json_bytes: 0,
    };
    write_result(result, &cli.out)?;
    let telemetry = TelemetryReceipt {
        schema: "p256.dickson_rr_selector_telemetry/v1".into(),
        total_wall_ms: telemetry_targets.iter().map(|target| target.wall_ms).sum(),
        targets: telemetry_targets,
        time_fit_exponent_in_variables: time_fit,
    };
    let telemetry_text =
        serde_json::to_string_pretty(&telemetry).map_err(|error| error.to_string())? + "\n";
    fs::write(&cli.telemetry_out, telemetry_text)
        .map_err(|error| format!("{}: {error}", cli.telemetry_out.display()))?;
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
    fn support_factorization_and_normalized_witness_hold() {
        let b = (&CurveParams::p256().b % P).to_u64().unwrap();
        let roots = square_roots();
        let pool = toy_support_pool(b, &roots);
        let support = &pool[..5];
        let (target, columns, mask) = planted_target(support, 3, b, &roots).unwrap();
        let witness = rr_witness(support, &columns, mask, target, b, &roots).unwrap();
        assert_eq!(witness.values.len(), 8);
        assert!(witness.equations_verified);
    }

    #[test]
    fn rr_system_has_registered_shape() {
        let b = (&CurveParams::p256().b % P).to_u64().unwrap();
        let roots = square_roots();
        let pool = toy_support_pool(b, &roots);
        for (m, columns) in TOY_SIZES {
            let (equations, variables) = rr_system(&pool[..columns], m, 17, b).unwrap();
            assert_eq!(variables, columns + m);
            assert_eq!(equations.len(), columns + m + 1);
            assert!(equations.iter().all(|equation| equation.degree() <= 2));
        }
    }

    #[test]
    fn reference_candidate_count_is_exact() {
        let b = (&CurveParams::p256().b % P).to_u64().unwrap();
        let roots = square_roots();
        let pool = toy_support_pool(b, &roots);
        let support = &pool[..5];
        let target = unplanted_target(support, 3, b, &roots).unwrap();
        let reference = reference_relations(support, 3, target, b, &roots).unwrap();
        assert_eq!(reference.candidates, 80);
        assert_eq!(reference.replay_failures, 0);
    }
}

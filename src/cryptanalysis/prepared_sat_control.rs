//! Prepared one-target SAT control on the already disclosed synthetic n17 instance.
//! Native execution/audit transport is supplied separately. No imported scalar.
use std::{collections::BTreeMap, time::Instant};

use num_bigint::BigUint;
use rand::{rngs::StdRng, Rng, SeedableRng};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};

use super::koblitz_fast::{FastCurve, FastPoint};
use crate::binary_ecc::{BinaryCurve, BinaryPoint, F2mElement, IrreduciblePoly};

pub const PREPARATION_SHA256: &str =
    "91856ab78550436d3f668367f9aebd9e2c0604bd64b1472d9d19ec318e2b144e";
pub const STATE_SHA256: &str = "edbff76da6442b9f2e5e8235682c9bf1052465f310c765ba8a37d018a60bf107";
pub const TARGET: [u64; 2] = [52411, 72106];
pub const ORDER: u64 = 65587;
pub const PHASES: [&str; 5] = [
    "target_query",
    "target_pdp",
    "target_relation_check",
    "target_descent",
    "target_recovery_check",
];

fn require(condition: bool, message: &str) -> Result<(), String> {
    if condition {
        Ok(())
    } else {
        Err(message.into())
    }
}
pub fn sha256(bytes: &[u8]) -> String {
    hex::encode(crate::hash::sha256::sha256(bytes))
}
pub fn canonical_sha(value: &Value) -> Result<String, String> {
    fn check(v: &Value) -> Result<(), String> {
        match v {
            Value::Number(n) if !n.is_u64() && !n.is_i64() => {
                Err("identity contains a float".into())
            }
            Value::Array(a) => a.iter().try_for_each(check),
            Value::Object(o) => o.values().try_for_each(check),
            _ => Ok(()),
        }
    }
    check(value)?;
    Ok(sha256(value.to_string().as_bytes()))
}
fn integer(value: &Value) -> Result<u64, String> {
    value
        .as_u64()
        .or_else(|| value.as_str().and_then(|s| s.parse().ok()))
        .ok_or("invalid integer".into())
}
fn pair(value: &Value) -> Result<[u64; 2], String> {
    let values = value
        .as_array()
        .filter(|v| v.len() == 2)
        .ok_or("invalid point pair")?;
    let p = [integer(&values[0])?, integer(&values[1])?];
    require(p.iter().all(|&v| v < 1 << 17), "point outside n17 field")?;
    Ok(p)
}
fn point(p: [u64; 2]) -> BinaryPoint {
    BinaryPoint::Affine {
        x: F2mElement::from_biguint(&BigUint::from(p[0]), 17),
        y: F2mElement::from_biguint(&BigUint::from(p[1]), 17),
    }
}
pub fn encoded(p: FastPoint) -> Option<[u64; 2]> {
    (!p.infinity).then_some([p.x, p.y])
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ControlConfig {
    pub schema_version: u32,
    pub question: String,
    pub target: [u64; 2],
    pub descent_query_seed: u64,
    pub export_nonce: u64,
    pub conflict_budget: u64,
    pub max_queries: usize,
    pub exporter_timeout_ms: u64,
    pub solver_timeout_ms: u64,
    pub controller_timeout_ms: u64,
}
impl ControlConfig {
    pub fn validate(&self) -> Result<(), String> {
        require(
            self.schema_version == 1
                && self.question == "native-prepared-known-input-control-v1"
                && self.target == TARGET,
            "only the disclosed one-target n17 control is admitted",
        )?;
        require(
            (1..=8).contains(&self.max_queries)
                && (100_000..=1_000_000).contains(&self.conflict_budget)
                && (1..=30_000).contains(&self.exporter_timeout_ms)
                && (1..=60_000).contains(&self.solver_timeout_ms)
                && (30_000..=900_000).contains(&self.controller_timeout_ms),
            "control resource/budget envelope differs",
        )
    }
}

pub struct PreparedState {
    pub curve: FastCurve,
    pub generator: FastPoint,
    pub geometry: Vec<FastPoint>,
    pub logs: Vec<u64>,
    pub projection: Vec<Option<(usize, u64)>>,
}
impl PreparedState {
    pub fn load(document: &Value) -> Result<Self, String> {
        require(
            canonical_sha(document)? == PREPARATION_SHA256,
            "preparation seal differs",
        )?;
        let fixture = &document["certificate"]["inputs"]["fixture"];
        require(
            fixture["degree"] == 17
                && fixture["curve_a"] == 1
                && integer(&fixture["subgroup_order"])? == ORDER
                && integer(&fixture["cofactor"])? == 2
                && fixture["irreducible"]["low_terms"] == json!([0, 3])
                && pair(&fixture["generator"])? == [43693, 23339]
                && fixture["targets"] == json!([])
                && fixture["target_scalar_constructed"] == false
                && document["record_sha256"] == STATE_SHA256,
            "preparation is not the exact target-independent toy state",
        )?;
        let curve = BinaryCurve {
            m: 17,
            irreducible: IrreduciblePoly {
                degree: 17,
                low_terms: vec![0, 3],
            },
            a: F2mElement::one(17),
            b: F2mElement::one(17),
            generator: point([43693, 23339]),
            order: BigUint::from(ORDER),
            cofactor: BigUint::from(2u64),
        };
        let fast = FastCurve::new(&curve).ok_or("native curve backend unavailable")?;
        let g = fast.lift(&curve.generator);
        require(
            curve.is_on_curve(&curve.generator) && fast.mul_u64(g, ORDER).infinity,
            "invalid generator",
        )?;
        let geometry = document["certificate"]["inputs"]["base"]
            .as_array()
            .filter(|a| a.len() == 63)
            .ok_or("geometric base count differs")?
            .iter()
            .map(|p| {
                let p = point(pair(p)?);
                require(curve.is_on_curve(&p), "geometric point is not on curve")?;
                Ok(fast.lift(&p))
            })
            .collect::<Result<Vec<_>, String>>()?;
        let columns = document["record"]["column_logs"]
            .as_array()
            .filter(|a| a.len() == 29)
            .ok_or("column count differs")?;
        let mut logs = Vec::new();
        let mut points = Vec::new();
        for column in columns {
            let p = fast.lift(&point(pair(&column["point"])?));
            let log = integer(&column["log"])?;
            require(
                log < ORDER && fast.mul_u64(g, log) == p,
                "prepared column log fails scalar replay",
            )?;
            points.push(p);
            logs.push(log);
        }
        let locations = document["record"]["projection"]
            .as_array()
            .filter(|a| a.len() == 63)
            .ok_or("projection count differs")?;
        let mut projection = Vec::new();
        for (&p, location) in geometry.iter().zip(locations) {
            let image = fast.mul_u64(p, 2);
            if location.is_null() {
                require(image.infinity, "projection dropped a usable point")?;
                projection.push(None);
            } else {
                let [column, coefficient] = pair(location)?;
                require(
                    column < 29
                        && coefficient < ORDER
                        && fast.mul_u64(points[column as usize], coefficient) == image,
                    "projection coefficient fails native group replay",
                )?;
                projection.push(Some((column as usize, coefficient)));
            }
        }
        Ok(Self {
            curve: fast,
            generator: g,
            geometry,
            logs,
            projection,
        })
    }
    pub fn target(&self, p: [u64; 2]) -> Result<FastPoint, String> {
        require(p == TARGET, "target differs from disclosed control")?;
        let p = FastPoint::affine(p[0], p[1]);
        require(
            self.curve.mul_u64(p, ORDER).infinity,
            "target is outside subgroup",
        )?;
        Ok(p)
    }
    pub fn recover(&self, a: u64, b: u64, indices: &[usize; 3]) -> Result<u64, String> {
        require(
            a < ORDER && b > 0 && b < ORDER && indices.iter().all(|&i| i < 63),
            "invalid descent relation",
        )?;
        let projected = indices
            .iter()
            .filter_map(|&i| self.projection[i])
            .fold(0, |sum, (column, coefficient)| {
                (sum + coefficient * self.logs[column] % ORDER) % ORDER
            });
        Ok((projected + ORDER - 2 * a % ORDER) % ORDER * pow(2 * b % ORDER, ORDER - 2) % ORDER)
    }
    pub fn lift(&self, model: &[bool], public: FastPoint) -> Option<[usize; 3]> {
        if model.len() < 51 {
            return None;
        }
        let xs: [u64; 3] =
            std::array::from_fn(|i| (0..6).fold(0, |v, b| v | ((model[6 * i + b] as u64) << b)));
        let candidates = xs.map(|x| {
            self.geometry
                .iter()
                .enumerate()
                .filter(|(_, p)| p.x == x)
                .map(|(i, _)| i)
                .collect::<Vec<_>>()
        });
        for &a in &candidates[0] {
            for &b in &candidates[1] {
                for &c in &candidates[2] {
                    if self.curve.add(
                        self.curve.add(self.geometry[a], self.geometry[b]),
                        self.geometry[c],
                    ) == public
                    {
                        return Some([a, b, c]);
                    }
                }
            }
        }
        None
    }
}
fn pow(mut a: u64, mut exponent: u64) -> u64 {
    let mut answer = 1;
    while exponent != 0 {
        if exponent & 1 == 1 {
            answer = answer * a % ORDER;
        }
        a = a * a % ORDER;
        exponent >>= 1;
    }
    answer
}

/// Every nanosecond of the target-dependent interval has one phase owner.
pub struct OnlineClock {
    start: Instant,
    last: Instant,
    owner: usize,
    costs: [u64; 5],
}
impl OnlineClock {
    pub fn start() -> Self {
        let start = Instant::now();
        Self {
            start,
            last: start,
            owner: 0,
            costs: [0; 5],
        }
    }
    pub fn phase(&mut self, owner: usize) {
        assert!(owner < 5);
        let now = Instant::now();
        self.costs[self.owner] += now.duration_since(self.last).as_nanos() as u64;
        self.last = now;
        self.owner = owner;
    }
    pub fn finish(mut self) -> (u64, BTreeMap<String, u64>) {
        let stop = Instant::now();
        self.costs[self.owner] += stop.duration_since(self.last).as_nanos() as u64;
        let total = stop.duration_since(self.start).as_nanos() as u64;
        assert_eq!(self.costs.iter().sum::<u64>(), total);
        (
            total,
            PHASES
                .into_iter()
                .zip(self.costs)
                .map(|(k, v)| (k.into(), v))
                .collect(),
        )
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NativeOutput {
    pub exit_code: Option<i32>,
    pub timed_out: bool,
    pub stdout: String,
    pub receipt: Value,
}
#[derive(Clone, Debug)]
pub struct QueryOutput {
    pub native: NativeOutput,
    pub manifest: Value,
    pub anf: String,
    pub cnf: String,
    pub source_receipt: Value,
}
pub trait QueryBackend {
    fn query(
        &mut self,
        trial: usize,
        public: [u64; 2],
        config: &ControlConfig,
        clock: &mut OnlineClock,
    ) -> Result<QueryOutput, String>;
    fn progress(&mut self, attempt: &Value) -> Result<(), String>;
}

/// Complete producer flow. Callback-only controls never establish real native admission.
pub fn solve<B: QueryBackend>(
    preparation: &Value,
    config: &ControlConfig,
    backend: &mut B,
) -> Result<Value, String> {
    config.validate()?;
    let prepared = PreparedState::load(preparation)?;
    let target = prepared.target(config.target)?;
    let mut rng = StdRng::seed_from_u64(config.descent_query_seed ^ 0x44455343454e5400);
    let mut clock = OnlineClock::start();
    let mut attempts = Vec::new();
    let mut scalar = None;
    for trial in 0..config.max_queries {
        clock.phase(0);
        let (a, b) = (rng.gen_range(1..ORDER), rng.gen_range(1..ORDER));
        let public = prepared.curve.add(
            prepared.curve.mul_u64(prepared.generator, a),
            prepared.curve.mul_u64(target, b),
        );
        clock.phase(1);
        let mut row = json!({"trial":trial,"a":a,"b":b,"public_point":encoded(public),"status":"IDENTITY_QUERY",
            "source_model_valid":false,"witness_indices":null,"candidate_scalar":null});
        if let Some(p) = encoded(public) {
            match backend.query(trial, p, config, &mut clock) {
                Err(reason) => {
                    row["status"] = json!("QUERY_TRANSPORT_FAILURE");
                    row["reason"] = json!(reason);
                }
                Ok(output) => {
                    row["native"] = json!(output.native);
                    row["source_receipt"] = output.source_receipt;
                    clock.phase(2);
                    let status = native_status(&output.native);
                    row["status"] = json!(status);
                    if status == "SAT_MODEL" {
                        let validated = (|| {
                            validate_manifest(&output.manifest, p)?;
                            let count = output.manifest["exports"]["cryptominisat_xor_dimacs"]
                                ["variables"]
                                .as_u64()
                                .ok_or("missing expanded variable count")?
                                as usize;
                            let model = parse_model(&output.native.stdout, count)?;
                            verify_anf(&output.anf, &model[..51])?;
                            verify_cnf(&output.cnf, &model)?;
                            let indices = prepared
                                .lift(&model, public)
                                .ok_or("source model does not lift to a group witness")?;
                            Ok::<_, String>((model, indices))
                        })();
                        match validated {
                            Err(reason) => {
                                row["status"] = json!("INVALID_OR_NONLIFTING_SOURCE_MODEL");
                                row["reason"] = json!(reason);
                            }
                            Ok((model, indices)) => {
                                row["source_model_valid"] = json!(true);
                                row["witness_indices"] = json!(indices);
                                row["source_model_sha256"] = json!(sha256(
                                    &model[..51].iter().map(|&b| u8::from(b)).collect::<Vec<_>>()
                                ));
                                row["expanded_model_sha256"] = json!(sha256(
                                    &model.iter().map(|&b| u8::from(b)).collect::<Vec<_>>()
                                ));
                                clock.phase(3);
                                let candidate = prepared.recover(a, b, &indices)?;
                                clock.phase(4);
                                require(
                                    prepared.curve.mul_u64(prepared.generator, candidate) == target,
                                    "target scalar fails native group replay",
                                )?;
                                scalar = Some(candidate);
                                row["status"] = json!("VALID_POINT_WITNESS");
                                row["candidate_scalar"] = json!(candidate);
                            }
                        }
                    }
                }
            }
        }
        if scalar.is_some() {
            attempts.push(row);
            break;
        }
        clock.phase(1);
        backend.progress(&row)?;
        attempts.push(row);
    }
    let (attempted, phases) = clock.finish();
    Ok(
        json!({"schema_version":1,"method":"native-prepared-sat-control-v1","status":if scalar.is_some(){"complete"}else{"incomplete"},
        "target_input":config.target,"config":config,"preparation_sha256":PREPARATION_SHA256,"mathematical_state_sha256":STATE_SHA256,
        "target_attempts":attempts,"recovered_scalar":scalar,"scalar_verified":scalar.is_some(),
        "online_attempt_wall_ns":attempted,"online_wall_ns":scalar.map(|_|attempted),"online_phases_ns":phases,
        "ordinary_queries_executed":0,"fresh_paired_qualification":false,"headline_eligible":false,
        "source_bound_execution_admitted":false,"promotion_eligible":false,"online_speedup":null}),
    )
}

pub fn native_status(output: &NativeOutput) -> &'static str {
    if output.timed_out {
        return "TIMEOUT";
    }
    let statuses = output
        .stdout
        .lines()
        .filter(|line| line.starts_with("s "))
        .collect::<Vec<_>>();
    match (output.exit_code, statuses.as_slice()) {
        (Some(10), ["s SATISFIABLE"]) => "SAT_MODEL",
        (Some(20), ["s UNSATISFIABLE"]) => "SOURCE_UNSAT",
        (Some(15), ["s INDETERMINATE"]) => "CONFLICT_BUDGET_INCONCLUSIVE",
        (Some(0), ["s UNKNOWN"]) => "UNKNOWN_INCONCLUSIVE",
        _ => "NATIVE_ERROR",
    }
}

pub fn validate_manifest(manifest: &Value, point: [u64; 2]) -> Result<(), String> {
    require(
        manifest["n"] == 17
            && manifest["curve_a"] == 1
            && manifest["ell"] == 6
            && manifest["m"] == 3
            && manifest["source_variables"] == 51
            && manifest["source_equations"] == 50
            && manifest["representation"] == "symmetrised_s4"
            && manifest["factor_base_basis_bitmasks"] == json!(["1", "2", "4", "8", "16", "32"])
            && integer(&manifest["target"]["x"])? == point[0]
            && integer(&manifest["target"]["y"])? == point[1],
        "source instance differs from registered query",
    )
}
pub fn parse_model(stdout: &str, count: usize) -> Result<Vec<bool>, String> {
    require(
        (51..=4096).contains(&count),
        "model width outside n17 envelope",
    )?;
    let mut model = vec![None; count];
    let mut ended = false;
    for line in stdout.lines().filter_map(|l| l.strip_prefix("v ")) {
        for token in line.split_whitespace() {
            require(!ended, "model has tokens after terminator")?;
            let literal = token.parse::<i32>().map_err(|_| "invalid model literal")?;
            if literal == 0 {
                ended = true;
                continue;
            }
            let index = literal.unsigned_abs() as usize;
            require(
                index > 0 && index <= count && model[index - 1].is_none(),
                "out-of-range or duplicated model variable",
            )?;
            model[index - 1] = Some(literal > 0);
        }
    }
    require(ended, "model is unterminated")?;
    model
        .into_iter()
        .map(|v| v.ok_or_else(|| "model is incomplete".into()))
        .collect()
}

fn lines(text: &str) -> impl Iterator<Item = &str> {
    text.lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('c'))
}
fn header(line: &str, n: usize) -> Result<usize, String> {
    let t = line.split_whitespace().collect::<Vec<_>>();
    require(
        t.len() == 4 && t[0] == "p" && t[1] == "cnf" && t[2].parse::<usize>().ok() == Some(n),
        "source header differs",
    )?;
    let rows = t[3]
        .parse::<usize>()
        .map_err(|_| "invalid source row count")?;
    require(
        rows > 0 && rows <= 100_000,
        "source row count outside envelope",
    )?;
    Ok(rows)
}
fn var(token: &str, n: usize) -> Result<usize, String> {
    let v = token
        .parse::<usize>()
        .map_err(|_| "invalid source variable")?;
    require(v > 0 && v <= n, "source variable outside model")?;
    Ok(v - 1)
}
pub fn verify_anf(text: &str, model: &[bool]) -> Result<(), String> {
    check_anf(text, model, true)
}
fn check_anf(text: &str, model: &[bool], enforce_model: bool) -> Result<(), String> {
    require(model.len() == 51, "source model is not 51 bits")?;
    let mut input = lines(text);
    let expected = header(input.next().ok_or("missing ANF header")?, 51)?;
    require(expected == 50, "ANF equation count differs")?;
    let mut count = 0;
    for line in input {
        let t = line.split_whitespace().collect::<Vec<_>>();
        require(
            t.len() >= 2 && t[0] == "x" && t.last() == Some(&"0"),
            "malformed ANF equation",
        )?;
        let (mut i, mut value, mut marker) = (1, false, false);
        while i < t.len() - 1 {
            if t[i] == "T" {
                require(!marker, "duplicate constant marker")?;
                marker = true;
                i += 1;
            } else if let Some(d) = t[i].strip_prefix('.') {
                let d = d.parse::<usize>().map_err(|_| "bad monomial degree")?;
                require(
                    d > 0 && d <= 51 && i + d < t.len() - 1,
                    "truncated monomial",
                )?;
                let mut product = true;
                for token in &t[i + 1..=i + d] {
                    product &= model[var(token, 51)?];
                }
                value ^= product;
                i += d + 1;
            } else {
                value ^= model[var(t[i], 51)?];
                i += 1;
            }
        }
        require(
            !enforce_model || value ^ marker,
            "source ANF equation fails",
        )?;
        count += 1;
    }
    require(count == expected, "source ANF row count differs")
}
pub fn verify_cnf(text: &str, model: &[bool]) -> Result<(), String> {
    check_cnf(text, model, true)
}
fn check_cnf(text: &str, model: &[bool], enforce_model: bool) -> Result<(), String> {
    let mut input = lines(text);
    let expected = header(input.next().ok_or("missing CNF header")?, model.len())?;
    let mut count = 0;
    for line in input {
        let t = line.split_whitespace().collect::<Vec<_>>();
        let xor = t.first() == Some(&"x");
        require(
            !t.is_empty() && t.last() == Some(&"0"),
            "unterminated source clause",
        )?;
        let (mut parity, mut satisfied) = (false, false);
        for token in &t[usize::from(xor)..t.len() - 1] {
            let literal = token.parse::<i32>().map_err(|_| "invalid source literal")?;
            let index = literal.unsigned_abs() as usize;
            require(
                index > 0 && index <= model.len(),
                "source literal outside model",
            )?;
            let value = model[index - 1] ^ (literal < 0);
            parity ^= value;
            satisfied |= value;
        }
        require(
            !enforce_model || if xor { parity } else { satisfied },
            "source CNF/XOR clause fails",
        )?;
        count += 1;
    }
    require(count == expected, "source CNF row count differs")
}
pub fn validate_source_syntax(anf: &str, cnf: &str, width: usize) -> Result<(), String> {
    require(
        (51..=4096).contains(&width),
        "source width outside envelope",
    )?;
    check_anf(anf, &[false; 51], false)?;
    check_cnf(cnf, &vec![false; width], false)
}

#[cfg(test)]
mod tests {
    use super::*;
    const PREP: &str = include_str!("../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-ic-state-v1/sat-preparation.json");
    const MANIFEST: &str = include_str!("../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-report-contract-v1/sat-source-check/manifest.json");
    const ANF: &str = include_str!("../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-report-contract-v1/sat-source-check/instance.anf");
    const CNF: &str = include_str!("../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-report-contract-v1/sat-source-check/instance.xor.cnf");
    const MODEL: &str = include_str!("../../research/ic_candidate_tournament_20260915/goal_20260924/prepared-report-contract-v1/sat-source-check/result.json");

    fn config() -> ControlConfig {
        ControlConfig {
            schema_version: 1,
            question: "native-prepared-known-input-control-v1".into(),
            target: TARGET,
            descent_query_seed: 2026100102,
            export_nonce: 2026100103,
            conflict_budget: 100_000,
            max_queries: 2,
            exporter_timeout_ms: 30_000,
            solver_timeout_ms: 60_000,
            controller_timeout_ms: 180_000,
        }
    }
    fn retained_model() -> Vec<bool> {
        serde_json::from_str::<Value>(MODEL).unwrap()["cnf_assignment"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_bool().unwrap())
            .collect()
    }
    fn native(code: i32, status: &str) -> NativeOutput {
        NativeOutput {
            exit_code: Some(code),
            timed_out: false,
            stdout: format!("s {status}\n"),
            receipt: json!({"test_fixture":true}),
        }
    }
    fn source(model: &[bool]) -> QueryOutput {
        let mut out = native(10, "SATISFIABLE");
        out.stdout.push_str("v ");
        for (i, &b) in model.iter().enumerate() {
            out.stdout.push_str(&format!(
                "{} ",
                if b { i as i32 + 1 } else { -(i as i32 + 1) }
            ));
        }
        out.stdout.push_str("0\n");
        QueryOutput {
            native: out,
            manifest: serde_json::from_str(MANIFEST).unwrap(),
            anf: ANF.into(),
            cnf: CNF.into(),
            source_receipt: json!({"retained_fixture":true}),
        }
    }
    struct Retained {
        model: Vec<bool>,
        failed: Vec<Value>,
        calls: usize,
    }
    impl QueryBackend for Retained {
        fn query(
            &mut self,
            trial: usize,
            p: [u64; 2],
            _: &ControlConfig,
            _: &mut OnlineClock,
        ) -> Result<QueryOutput, String> {
            self.calls += 1;
            if trial == 0 {
                return Ok(QueryOutput {
                    native: native(15, "INDETERMINATE"),
                    manifest: Value::Null,
                    anf: String::new(),
                    cnf: String::new(),
                    source_receipt: json!({"test_fixture":true}),
                });
            }
            assert_eq!(p, [62577, 27783]);
            Ok(source(&self.model))
        }
        fn progress(&mut self, row: &Value) -> Result<(), String> {
            self.failed.push(row.clone());
            Ok(())
        }
    }
    #[test]
    fn complete_retained_control_recovers_without_scalar_input() {
        let prep = serde_json::from_str(PREP).unwrap();
        let mut backend = Retained {
            model: retained_model(),
            failed: Vec::new(),
            calls: 0,
        };
        let report = solve(&prep, &config(), &mut backend).unwrap();
        assert_eq!(report["status"], "complete");
        assert_eq!(report["recovered_scalar"], 24886);
        assert_eq!(backend.failed.len(), 1);
        assert_eq!(backend.calls, 2);
        assert_eq!(
            report["target_attempts"][0]["status"],
            "CONFLICT_BUDGET_INCONCLUSIVE"
        );
        assert_eq!(
            report["target_attempts"][1]["witness_indices"],
            json!([29, 51, 2])
        );
        assert_eq!(
            report["online_phases_ns"]
                .as_object()
                .unwrap()
                .values()
                .map(|v| v.as_u64().unwrap())
                .sum::<u64>(),
            report["online_wall_ns"].as_u64().unwrap()
        );
        assert_eq!(report["headline_eligible"], false);
        assert_eq!(report["source_bound_execution_admitted"], false);
    }
    #[test]
    fn all_expanded_model_bit_corruptions_reject() {
        let model = retained_model();
        verify_anf(ANF, &model[..51]).unwrap();
        verify_cnf(CNF, &model).unwrap();
        for bit in 0..model.len() {
            let mut bad = model.clone();
            bad[bit] = !bad[bit];
            assert!(
                verify_anf(ANF, &bad[..51]).is_err() || verify_cnf(CNF, &bad).is_err(),
                "accepted corruption {bit}"
            );
        }
    }
    #[test]
    fn invalid_model_retained_as_failure_and_never_complete() {
        let prep = serde_json::from_str(PREP).unwrap();
        let mut model = retained_model();
        model[0] = !model[0];
        let mut backend = Retained {
            model,
            failed: Vec::new(),
            calls: 0,
        };
        let report = solve(&prep, &config(), &mut backend).unwrap();
        assert_eq!(report["status"], "incomplete");
        assert!(report["recovered_scalar"].is_null());
        assert!(report["online_wall_ns"].is_null());
        assert_eq!(backend.failed.len(), 2);
        assert_eq!(
            report["target_attempts"][1]["status"],
            "INVALID_OR_NONLIFTING_SOURCE_MODEL"
        );
    }
    #[test]
    fn statuses_require_exit_and_unique_native_status() {
        for (code, status, want) in [
            (10, "SATISFIABLE", "SAT_MODEL"),
            (20, "UNSATISFIABLE", "SOURCE_UNSAT"),
            (15, "INDETERMINATE", "CONFLICT_BUDGET_INCONCLUSIVE"),
            (0, "UNKNOWN", "UNKNOWN_INCONCLUSIVE"),
            (0, "SATISFIABLE", "NATIVE_ERROR"),
            (10, "UNKNOWN", "NATIVE_ERROR"),
        ] {
            let mut out = native(code, status);
            assert_eq!(native_status(&out), want);
            out.stdout.push_str("s UNKNOWN\n");
            assert_eq!(native_status(&out), "NATIVE_ERROR");
            out.timed_out = true;
            assert_eq!(native_status(&out), "TIMEOUT");
        }
    }
    #[test]
    fn rejects_partial_duplicate_unterminated_and_oversized_models() {
        for text in [
            "v 1 0\n",
            "v 1 1 0\n",
            "v 1\n",
            "v 2147483647 0\n",
            "v -2147483648 0\n",
            "v 0 1\n",
        ] {
            assert!(parse_model(text, 51).is_err());
        }
        assert!(parse_model("v 0\n", 50).is_err());
        assert!(parse_model("v 0\n", 4097).is_err());
    }
    #[test]
    fn seal_config_and_source_layout_fail_closed() {
        let mut prep: Value = serde_json::from_str(PREP).unwrap();
        prep["record"]["column_logs"][0]["log"] = json!(0);
        assert!(PreparedState::load(&prep).is_err());
        let mut cfg = config();
        cfg.target = [1, 2];
        assert!(cfg.validate().is_err());
        let mut manifest: Value = serde_json::from_str(MANIFEST).unwrap();
        validate_manifest(&manifest, [62577, 27783]).unwrap();
        assert!(validate_manifest(&manifest, [62577, 27784]).is_err());
        manifest["source_variables"] = json!(52);
        assert!(validate_manifest(&manifest, [62577, 27783]).is_err());
    }
}

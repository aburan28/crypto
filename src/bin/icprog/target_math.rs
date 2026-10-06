//! Independent data-only audit for the new disclosed n17 F5 target worker.
//! No producer arithmetic, search, worker invocation or runtime admission.
use super::{
    oracle::Curve, ordinary_preparation, prepared_sat, sat_control::native, sat_query_law::QueryLaw,
};
use crypto_lib::cryptanalysis::prepared_sat_control::{canonical_sha, sha256};
use serde::{
    de::{self, MapAccess, SeqAccess, Visitor},
    Deserialize, Deserializer, Serialize,
};
use serde_json::{json, Map, Value};
use std::{collections::BTreeSet, fmt, fs, path::Path, time::Instant};

const SCOPE: &str = "disclosed-n17-native-one-target-f5-v1";
const RECORD_SCOPE: &str = "disclosed-n17-target-attempt-files-v1";
const PHASES: [&str; 14] = [
    "setup",
    "factor_base",
    "precompute",
    "queries",
    "pdp",
    "relation_check",
    "matrix_build",
    "relation_la",
    "target_query",
    "target_pdp",
    "target_relation_check",
    "target_descent",
    "recovery_check",
    "rho_solve",
];
const ONLINE: [&str; 6] = [
    "target_query",
    "target_pdp",
    "target_relation_check",
    "target_descent",
    "recovery_check",
    "rho_solve",
];

/// Reject duplicate keys for both record types. Only the source-bound ordinary
/// preparation producer may carry finite diagnostic floats outside its math.
struct StrictJson<const ALLOW_FLOAT: bool>(Value);
impl<'de, const ALLOW_FLOAT: bool> Deserialize<'de> for StrictJson<ALLOW_FLOAT> {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        struct V<const ALLOW_FLOAT: bool>;
        impl<'de, const ALLOW_FLOAT: bool> Visitor<'de> for V<ALLOW_FLOAT> {
            type Value = StrictJson<ALLOW_FLOAT>;
            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                if ALLOW_FLOAT {
                    f.write_str("JSON with unique keys and finite numbers")
                } else {
                    f.write_str("JSON with unique keys and integer numbers")
                }
            }
            fn visit_bool<E: de::Error>(self, v: bool) -> Result<Self::Value, E> {
                Ok(StrictJson(json!(v)))
            }
            fn visit_u64<E: de::Error>(self, v: u64) -> Result<Self::Value, E> {
                Ok(StrictJson(json!(v)))
            }
            fn visit_i64<E: de::Error>(self, v: i64) -> Result<Self::Value, E> {
                Ok(StrictJson(json!(v)))
            }
            fn visit_f64<E: de::Error>(self, v: f64) -> Result<Self::Value, E> {
                if !ALLOW_FLOAT {
                    return Err(de::Error::custom(
                        "floating point is forbidden in target records",
                    ));
                }
                let n = serde_json::Number::from_f64(v)
                    .ok_or_else(|| de::Error::custom("nonfinite diagnostic number"))?;
                Ok(StrictJson(Value::Number(n)))
            }
            fn visit_str<E: de::Error>(self, v: &str) -> Result<Self::Value, E> {
                Ok(StrictJson(json!(v)))
            }
            fn visit_string<E: de::Error>(self, v: String) -> Result<Self::Value, E> {
                Ok(StrictJson(json!(v)))
            }
            fn visit_none<E: de::Error>(self) -> Result<Self::Value, E> {
                Ok(StrictJson(Value::Null))
            }
            fn visit_unit<E: de::Error>(self) -> Result<Self::Value, E> {
                Ok(StrictJson(Value::Null))
            }
            fn visit_seq<A: SeqAccess<'de>>(self, mut a: A) -> Result<Self::Value, A::Error> {
                let mut out = Vec::new();
                while let Some(v) = a.next_element::<StrictJson<ALLOW_FLOAT>>()? {
                    out.push(v.0);
                }
                Ok(StrictJson(Value::Array(out)))
            }
            fn visit_map<A: MapAccess<'de>>(self, mut a: A) -> Result<Self::Value, A::Error> {
                let mut out = Map::new();
                while let Some(k) = a.next_key::<String>()? {
                    if out.contains_key(&k) {
                        return Err(de::Error::custom("duplicate JSON key"));
                    }
                    out.insert(k, a.next_value::<StrictJson<ALLOW_FLOAT>>()?.0);
                }
                Ok(StrictJson(Value::Object(out)))
            }
        }
        d.deserialize_any(V::<ALLOW_FLOAT>)
    }
}
pub(super) fn parse(bytes: &[u8]) -> Result<Value, String> {
    serde_json::from_slice::<StrictJson<false>>(bytes)
        .map(|v| v.0)
        .map_err(|e| e.to_string())
}
pub(super) fn parse_ordinary_producer(bytes: &[u8]) -> Result<Value, String> {
    serde_json::from_slice::<StrictJson<true>>(bytes)
        .map(|v| v.0)
        .map_err(|e| e.to_string())
}
fn keys(v: &Value, fields: &[&str]) -> Result<(), String> {
    native::require(
        v.as_object()
            .is_some_and(|o| o.len() == fields.len() && fields.iter().all(|f| o.contains_key(*f))),
        "undeclared or missing audit fields",
    )
}
fn number(v: &Value) -> Result<u64, String> {
    v.as_u64().ok_or("missing unsigned integer".into())
}
fn digest(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn power(mut a: u64, mut n: u64, r: u64) -> u64 {
    let mut out = 1;
    while n != 0 {
        if n & 1 != 0 {
            out = out * a % r;
        }
        a = a * a % r;
        n >>= 1;
    }
    out
}
#[derive(Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Plan {
    schema_version: u32,
    question: String,
    algorithm_seed: u64,
    max_queries: usize,
}
#[derive(Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Config {
    schema_version: u32,
    question: String,
    target: [u64; 2],
    plan: Plan,
    worker_timeout_ms: u64,
}
impl Config {
    fn validate(&self) -> Result<(), String> {
        native::require(
            self.schema_version == 1
                && self.question == SCOPE
                && self.target.iter().all(|&v| v < 1 << 17)
                && self.plan.schema_version == 1
                && self.plan.question == "prepared-public-target-n17-f5-v1"
                && (1..=8).contains(&self.plan.max_queries)
                && (30_000..=900_000).contains(&self.worker_timeout_ms),
            "outside exact n17 target audit scope",
        )
    }
}
fn decode(curve: &Curve, v: &Value) -> Result<super::oracle::Point, String> {
    curve.decode(&super::json::parse(&v.to_string())?)
}

/// Independently replay the producer report. No claim about who executed it.
pub(super) fn verify(preparation: &Value, cfg: &Config, producer: &Value) -> Result<Value, String> {
    cfg.validate()?;
    let prep = ordinary_preparation::audit(preparation)?;
    native::require(
        prep["mathematical_preparation_complete"] == true
            && prep["declared_family"] == "matrix_f5"
            && prep["rank"] == 29
            && prep["usable_points"] == 62,
        "target lacks complete independently reconstructed F5 preparation",
    )?;
    keys(
        producer,
        &[
            "schema_version",
            "question",
            "plan",
            "target",
            "target_count",
            "preparation_family",
            "preparation_math_sha256",
            "actual_usable_base_points",
            "geometric_base_points",
            "effective_columns",
            "dispatch",
            "reusable_validation_ns",
            "solver_setup_ns",
            "report",
            "verified_scalar",
            "independent_recovery_certificate",
            "failure",
            "interruption",
            "costs",
            "exclusive_ic_online_phases_ns",
            "timing_class",
            "candidate_id",
            "workload_id",
            "run_id",
            "source_bound_execution_admitted",
            "fresh_paired_qualification",
            "headline_eligible",
            "promotion_eligible",
            "full_goal_complete",
            "online_speedup",
            "operation_counts",
        ],
    )?;
    native::require(
        producer["schema_version"] == 1
            && producer["question"] == "prepared-public-target-n17-f5-producer-v1"
            && producer["plan"] == json!(cfg.plan)
            && producer["target"] == json!(cfg.target)
            && producer["target_count"] == 1
            && producer["preparation_family"] == "matrix_f5"
            && producer["preparation_math_sha256"] == canonical_sha(preparation)?
            && producer["actual_usable_base_points"] == 62
            && producer["geometric_base_points"] == 63
            && producer["effective_columns"] == 29
            && producer["timing_class"] == "ordinary-host-diagnostic"
            && producer["interruption"].is_null(),
        "target producer scope/preparation/count differs or solve is interrupted",
    )?;
    for f in [
        "source_bound_execution_admitted",
        "fresh_paired_qualification",
        "headline_eligible",
        "promotion_eligible",
        "full_goal_complete",
    ] {
        native::require(
            producer[f] == false,
            "producer cannot admit its own runtime or comparison",
        )?;
    }
    for f in [
        "candidate_id",
        "workload_id",
        "run_id",
        "online_speedup",
        "operation_counts",
    ] {
        native::require(
            producer[f].is_null(),
            "producer adds unregistered identity or performance claim",
        )?;
    }
    number(&producer["reusable_validation_ns"])?;
    number(&producer["solver_setup_ns"])?;
    let dispatch = &producer["dispatch"];
    native::require(
        dispatch
            == &json!({"strategy":"Groebner", "field_kernel":dispatch["field_kernel"], "pair_table":false,
        "query_rule":"seeded-sample", "summands":3, "direct_collision":false, "engine":{"MatrixF5":{"max_degree":3}}, "node_budget":8192,
        "collapse_negation":true, "collapse_projected_orbits":false, "attempt_record_contract":"start-and-completed-pdp-before-recovery-v1"})
            && matches!(
                dispatch["field_kernel"].as_str(),
                Some("pmull" | "pclmulqdq" | "portable")
            ),
        "undeclared target dispatch",
    )?;
    let curve = Curve::new(&super::json::parse(&preparation["fixture"].to_string())?)?;
    let target = decode(&curve, &json!(cfg.target))?;
    native::require(
        target.is_some() && curve.mul(target, curve.r as u128).is_none(),
        "target outside nonidentity n17 subgroup",
    )?;
    let base = preparation["geometric_base"]
        .as_array()
        .ok_or("missing geometry")?
        .iter()
        .map(|p| decode(&curve, p))
        .collect::<Result<Vec<_>, _>>()?;
    let (_, projection) = prepared_sat::projected_columns(&curve, &base)?;
    let logs = prep["independently_recovered_column_logs"]
        .as_array()
        .ok_or("missing checked logs")?
        .iter()
        .map(number)
        .collect::<Result<Vec<_>, _>>()?;
    let pairs = prepared_sat::pair_sums(&curve, &base);
    let report = &producer["report"];
    keys(
        report,
        &["trials", "claimed_scalar", "relation", "attempts"],
    )?;
    let attempts = report["attempts"]
        .as_array()
        .filter(|a| !a.is_empty() && a.len() <= cfg.plan.max_queries)
        .ok_or("missing or over-cap target attempts")?;
    native::require(
        report["trials"] == attempts.len(),
        "target trial count differs",
    )?;
    let mut law = QueryLaw::new(cfg.plan.algorithm_seed);
    let mut recovered = None;
    let mut relation = Value::Null;
    let mut audited = Vec::new();
    let mut mix = std::collections::BTreeMap::<String, usize>::new();
    for (trial, attempt) in attempts.iter().enumerate() {
        native::require(
            recovered.is_none(),
            "query continued after recoverable IC witness",
        )?;
        keys(attempt, &["trial", "a", "b", "pdp"])?;
        let (a, b) = (law.scalar(), law.scalar());
        native::require(
            attempt["trial"] == trial && attempt["a"] == a && attempt["b"] == b,
            "target query differs from independent seeded law",
        )?;
        let query = curve.add(
            curve.mul(Some(curve.g), a as u128),
            curve.mul(target, b as u128),
        );
        let pdp = &attempt["pdp"];
        keys(pdp, &["outcome", "points", "stats"])?;
        let outcome = pdp["outcome"].as_str().ok_or("missing outcome")?;
        let stats = &pdp["stats"]["stats"];
        if outcome == "identity" {
            native::require(
                query.is_none()
                    && pdp["points"].is_null()
                    && pdp["stats"] == json!({"family":"none"}),
                "false identity or hidden identity work",
            )?;
            // The frozen IC source skips a direct collision; never count it as IC recovery.
        } else {
            native::require(
                query.is_some(),
                "nonidentity PDP disposition on identity query",
            )?;
            keys(&pdp["stats"], &["family", "engine", "stats"])?;
            keys(
                stats,
                &[
                    "reductions",
                    "infeasible_branches",
                    "propagations",
                    "splits",
                    "exhausted",
                    "unsupported",
                    "max_degree_built",
                    "oversize",
                    "eliminated",
                ],
            )?;
            native::require(
                pdp["stats"]["family"] == "groebner"
                    && pdp["stats"]["engine"] == json!({"MatrixF5":{"max_degree":3}})
                    && number(&stats["reductions"])? <= 8192
                    && number(&stats["max_degree_built"])? <= 3
                    && stats["exhausted"].is_boolean()
                    && stats["unsupported"].is_boolean(),
                "target engine, budget or status differs",
            )?;
            for k in [
                "infeasible_branches",
                "propagations",
                "splits",
                "oversize",
                "eliminated",
            ] {
                number(&stats[k])?;
            }
            match outcome {
                "witness" => {
                    native::require(
                        stats["unsupported"] == false,
                        "unsupported frontend claims witness",
                    )?;
                    let indices = pdp["points"]
                        .as_array()
                        .filter(|p| p.len() == 3)
                        .ok_or("missing three-point witness")?
                        .iter()
                        .map(|v| {
                            number(v).and_then(|v| {
                                usize::try_from(v)
                                    .ok()
                                    .filter(|&v| v < base.len())
                                    .ok_or("witness index outside geometry".into())
                            })
                        })
                        .collect::<Result<Vec<_>, _>>()?;
                    native::require(
                        indices.iter().fold(None, |s, &i| curve.add(s, base[i])) == query,
                        "target witness fails independent full-point addition",
                    )?;
                    let sum = indices
                        .iter()
                        .filter_map(|&i| projection[i])
                        .fold(0, |s, (c, k)| (s + k * logs[c] % curve.r) % curve.r);
                    let denominator = curve.h * b % curve.r;
                    native::require(denominator != 0, "noninvertible projected recovery")?;
                    let d = (sum + curve.r - curve.h * a % curve.r) % curve.r
                        * power(denominator, curve.r - 2, curve.r)
                        % curve.r;
                    native::require(
                        curve.mul(Some(curve.g), d as u128) == target,
                        "projected IC recovery fails independent scalar replay",
                    )?;
                    recovered = Some(d);
                    relation = json!({"a":a,"b":b,"points":indices});
                }
                "proved_unsat" => native::require(
                    pdp["points"].is_null()
                        && stats["exhausted"] == false
                        && stats["unsupported"] == false
                        && !prepared_sat::has_three_sum(&curve, &base, &pairs, query),
                    "incomplete or false geometric negative",
                )?,
                "incomplete" | "unsupported" => native::require(
                    pdp["points"].is_null()
                        && stats["exhausted"] == true
                        && stats["unsupported"] == (outcome == "unsupported"),
                    "inconclusive target disposition differs",
                )?,
                _ => {
                    return Err(
                        "unaccepted target outcome; original failure must remain retained".into(),
                    )
                }
            }
        }
        *mix.entry(outcome.into()).or_default() += 1;
        audited.push(json!({"trial":trial,"a":a,"b":b,"public_query":query,"outcome":outcome}));
    }
    let complete = recovered.is_some();
    native::require(
        report["relation"] == relation
            && report["claimed_scalar"]
                == recovered
                    .map(|d| json!(d.to_string()))
                    .unwrap_or(Value::Null)
            && producer["verified_scalar"] == json!(recovered),
        "answer is not the recorded IC recovery",
    )?;
    if let Some(d) = recovered {
        native::require(
            producer["failure"].is_null()
                && producer["independent_recovery_certificate"]
                    == json!({"check":"independent-shift-and-add-scalar-replay",
            "generator":[curve.g.0,curve.g.1],"scalar":d,"target":cfg.target,"inside_online_interval":true}),
            "online independent recovery certificate differs",
        )?;
    } else {
        native::require(
            attempts.len() == cfg.plan.max_queries
                && producer["failure"] == "query budget exhausted without scalar recovery"
                && producer["independent_recovery_certificate"].is_null(),
            "unexplained early stop or failed solve claims recovery",
        )?;
    }
    let timing = verify_costs(producer)?;
    Ok(
        json!({"schema_version":1,"status":"PASS_TARGET_MATHEMATICS_RUNTIME_NOT_ADMITTED", "target":cfg.target,
        "preparation_math_sha256":canonical_sha(preparation)?,"preparation":prep,"audited_attempts":audited,"outcome_mix":mix,
        "verified_recovery":complete,"verified_scalar":recovered,"timing":timing,"solver_calls_by_auditor":0,
        "source_bound_execution_admitted":false,"fresh_paired_qualification":false,"headline_eligible":false,"promotion_eligible":false,
        "full_goal_complete":false,"online_wall_ns":null,"online_speedup":null}),
    )
}

pub(super) fn verify_costs(producer: &Value) -> Result<Value, String> {
    let clock = &producer["costs"];
    keys(
        clock,
        &[
            "schema_version",
            "phases_ns",
            "observed_wall_ns",
            "online_phases_ns",
            "online_wall_ns",
        ],
    )?;
    native::require(clock["schema_version"] == 1, "timing schema differs")?;
    keys(&clock["phases_ns"], &PHASES)?;
    keys(&clock["online_phases_ns"], &ONLINE)?;
    let observed = number(&clock["observed_wall_ns"])?;
    let online = number(&clock["online_wall_ns"])?;
    let mut raw_sum = 0u64;
    let mut online_sum = 0u64;
    for (name, v) in clock["phases_ns"].as_object().ok_or("missing phases")? {
        if !v.is_null() {
            raw_sum = raw_sum
                .checked_add(number(v)?)
                .ok_or("raw phase overflow")?;
        }
        if ONLINE.contains(&name.as_str()) {
            native::require(
                clock["online_phases_ns"][name] == *v,
                "target-only raw and online phase differ",
            )?;
            if !v.is_null() {
                online_sum = online_sum
                    .checked_add(number(v)?)
                    .ok_or("online phase overflow")?;
            }
        } else if name != "setup" {
            native::require(v.is_null(), "reusable work inside target-only clock")?;
        }
    }
    native::require(
        raw_sum == observed
            && online_sum == online
            && online > 0
            && online <= observed
            && clock["online_phases_ns"]["rho_solve"].is_null(),
        "target timing ledger does not close",
    )?;
    let mut exclusive = Map::new();
    let mut missing = Vec::new();
    for (name, exported) in [
        ("target_query", "target_query"),
        ("target_pdp", "target_pdp"),
        ("target_relation_check", "target_relation_check"),
        ("target_descent", "target_descent"),
        ("recovery_check", "target_recovery_check"),
    ] {
        let v = &clock["online_phases_ns"][name];
        if v.is_null() {
            missing.push(exported);
        }
        exclusive.insert(exported.into(), v.clone());
    }
    native::require(
        producer["exclusive_ic_online_phases_ns"] == Value::Object(exclusive.clone()),
        "five-phase export differs from raw clock",
    )?;
    Ok(
        json!({"producer_recorded_online_interval_ns":online,"five_phase_recorded_total_ns":if missing.is_empty(){Some(online)}else{None},
        "exclusive_ic_online_phases_ns":exclusive,"missing_costs":missing,"source_timing_attested":false}),
    )
}

fn record_inventory(root: &Path, cap: usize) -> Result<Value, String> {
    let mut names = BTreeSet::new();
    for entry in fs::read_dir(root).map_err(|e| e.to_string())? {
        let e = entry.map_err(|e| e.to_string())?;
        native::require(
            e.file_type().map_err(|e| e.to_string())?.is_file(),
            "target journal member not regular",
        )?;
        names.insert(
            e.file_name()
                .into_string()
                .map_err(|_| "non-UTF8 target member")?,
        );
        native::require(names.len() <= 1 + 2 * cap, "target journal exceeds cap")?;
    }
    let files = names
        .iter()
        .map(|name| {
            let bytes = native::read(
                &root.join(name),
                if name == "binding.json" {
                    65536
                } else {
                    8 * 1024 * 1024
                },
            )?;
            Ok(json!({"path":name,"bytes":bytes.len(),"sha256":sha256(&bytes)}))
        })
        .collect::<Result<Vec<Value>, String>>()?;
    Ok(json!({"names":names,"files":files}))
}

/// Compare every original durable envelope with the mathematically audited report.
pub(super) fn verify_records(
    execution: &Path,
    cfg: &Config,
    producer: &Value,
    seal: &str,
    worker: &str,
) -> Result<Value, String> {
    native::require(
        digest(seal) && digest(worker),
        "explicit target record seal/worker pin required",
    )?;
    let root = execution.join("attempts");
    let before = record_inventory(&root, cfg.plan.max_queries)?;
    let binding = parse(&native::read(&root.join("binding.json"), 65536)?)?;
    native::require(
        binding
            == json!({"schema_version":1,"question":RECORD_SCOPE,"registration_sha256":seal,"worker_sha256":worker,
        "family":"matrix_f5","target":cfg.target,"algorithm_seed":cfg.plan.algorithm_seed,"max_queries":cfg.plan.max_queries}),
        "durable target binding differs",
    )?;
    let attempts = producer["report"]["attempts"]
        .as_array()
        .filter(|a| !a.is_empty() && a.len() <= cfg.plan.max_queries)
        .ok_or("missing bounded target attempts")?;
    let mut expected = BTreeSet::from(["binding.json".to_string()]);
    let mut law = QueryLaw::new(cfg.plan.algorithm_seed);
    let binding_sha = canonical_sha(&binding)?;
    for (trial, attempt) in attempts.iter().enumerate() {
        let (a, b) = (law.scalar(), law.scalar());
        let header = json!({"trial":trial,"a":a,"b":b,"target":cfg.target});
        let start_name = format!("start-{trial:03}.json");
        let done_name = format!("completed-{trial:03}.json");
        expected.insert(start_name.clone());
        expected.insert(done_name.clone());
        let start = parse(&native::read(&root.join(start_name), 8 * 1024 * 1024)?)?;
        let done = parse(&native::read(&root.join(done_name), 8 * 1024 * 1024)?)?;
        let start_expected = json!({"schema_version":1,"binding_sha256":binding_sha,"header":header,"start_sha256":null,
            "body":{"trial":trial,"a":a,"b":b,"target":cfg.target,"query_rule":"seeded-sample-aG-plus-bQ","decomposition_started":false}});
        native::require(
            start == start_expected
                && done
                    == json!({"schema_version":1,"binding_sha256":binding_sha,"header":header,
            "start_sha256":canonical_sha(&start_expected)?,"body":{"target":cfg.target,"attempt":attempt,"record_kind":"completed-pdp-before-scalar-recovery"}}),
            "durable target records differ from independent law or producer attempts",
        )?;
    }
    native::require(
        before["names"] == json!(expected)
            && record_inventory(&root, cfg.plan.max_queries)? == before,
        "target record set has gaps/pending/unknown members or changed during audit",
    )?;
    Ok(before)
}

pub(super) fn run(
    preparation: &Path,
    config: &Path,
    producer: &Path,
    execution: &Path,
    seal: &str,
    worker: &str,
    out: &Path,
) -> Result<String, String> {
    let started = Instant::now();
    let execution = execution.canonicalize().map_err(|e| e.to_string())?;
    let output_parent = out
        .parent()
        .ok_or("audit output lacks parent")?
        .canonicalize()
        .map_err(|e| e.to_string())?;
    native::require(
        !output_parent.starts_with(&execution),
        "audit output must be outside original target evidence",
    )?;
    let inputs = [preparation, config, producer]
        .map(|p| native::read(p, 16 * 1024 * 1024))
        .into_iter()
        .collect::<Result<Vec<_>, _>>()?;
    let before = record_inventory(&execution.join("attempts"), 8)?;
    let result = (|| {
        let cfg: Config = serde_json::from_value(parse(&inputs[1])?).map_err(|e| e.to_string())?;
        let report = parse(&inputs[2])?;
        let mut result = verify(&parse(&inputs[0])?, &cfg, &report)?;
        result["attempt_files"] = verify_records(&execution, &cfg, &report, seal, worker)?;
        Ok::<Value, String>(result)
    })();
    for (path, bytes) in [preparation, config, producer].into_iter().zip(&inputs) {
        native::require(
            native::read(path, 16 * 1024 * 1024)? == *bytes,
            "target audit input changed during replay",
        )?;
    }
    native::require(
        record_inventory(&execution.join("attempts"), 8)? == before,
        "target journal changed during replay",
    )?;
    let error = result.as_ref().err().cloned();
    let mut receipt = result.unwrap_or_else(|e| json!({"schema_version":1,"status":"TARGET_MATHEMATICS_NOT_VERIFIED","error":e,
        "source_bound_execution_admitted":false,"verified_recovery":false,"fresh_paired_qualification":false,"headline_eligible":false,
        "promotion_eligible":false,"full_goal_complete":false,"online_wall_ns":null,"online_speedup":null,"solver_calls_by_auditor":0}));
    receipt["input_sha256"] = json!({"preparation":sha256(&inputs[0]),"config":sha256(&inputs[1]),"producer":sha256(&inputs[2])});
    // Retain byte custody of a rejected transcript too, without admitting it.
    receipt["attempt_files"] = before;
    receipt["registration_sha256"] = json!(seal);
    receipt["worker_sha256"] = json!(worker);
    receipt["independent_audit_wall_ns"] = json!(started.elapsed().as_nanos().to_string());
    receipt["audit_scope"] = json!("original target data and mathematics only; source/build/runtime/timing custody remain unverified");
    native::save(out, &receipt)?;
    if let Some(error) = error {
        return Err(format!(
            "target mathematical audit failed; receipt retained: {error}"
        ));
    }
    serde_json::to_string_pretty(&receipt).map_err(|e| e.to_string())
}

#[cfg(test)]
#[path = "target_math_tests.rs"]
mod tests;

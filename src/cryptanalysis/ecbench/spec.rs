//! The experiment spec (`ecbench.spec/v1`) and the plan it expands to.
//!
//! A spec names *what* is measured: curves, how many independent
//! single-target workloads per curve, the arms (methods), and how many
//! interleaved rounds.  *Where* it runs (CPUs, lock, thresholds) is a
//! property of the host and a flag of `ecbench run`, recorded in the
//! session, so the same spec runs unchanged on a laptop, a lab box and an
//! isolab worker.
//!
//! Expansion is deterministic: the same spec gives the same workloads,
//! the same execution order and the same algorithm seeds everywhere.

use serde::{Deserialize, Serialize};
use serde_json::{json, Value};

use crate::cryptanalysis::ecbench::canonical::{derive_u64, short_id};
use crate::cryptanalysis::ecbench::methods::{decl, resolve, Applies, MethodSpec, ResolvedMethod};
use crate::cryptanalysis::ecbench::workload::{CurveSpec, Instance, Workload};

pub const SPEC_SCHEMA: &str = "ecbench.spec/v1";

/// What an arm is in the comparison.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Role {
    /// The method the boundary is measured against (normally rho).
    Reference,
    /// The unmodified method a candidate changes.
    Baseline,
    Candidate,
    /// A repeat of another arm: the noise floor (A/A).
    Control,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ArmSpec {
    pub name: String,
    pub role: Role,
    pub method: MethodSpec,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WorkloadPlan {
    pub curves: Vec<CurveSpec>,
    /// Independent single-target workloads per curve.  Each is its own
    /// row; they are never pooled into a multi-target run.
    pub targets_per_curve: u64,
    pub target_seed: u64,
}

/// Order of the arms inside a round.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Order {
    /// Spec order on even rounds, reversed on odd ones: a linear drift
    /// (thermal, a slowly filling page cache) cancels over a pair.
    Alternate,
    /// A seeded permutation per round.
    Shuffle,
}

/// The isolation level a wall-clock figure needs (`docs/ecbench/README.md`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum Level {
    L0,
    L1,
    L2,
    L3,
}

impl Level {
    pub fn name(self) -> &'static str {
        match self {
            Level::L0 => "L0",
            Level::L1 => "L1",
            Level::L2 => "L2",
            Level::L3 => "L3",
        }
    }
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "L0" => Some(Level::L0),
            "L1" => Some(Level::L1),
            "L2" => Some(Level::L2),
            "L3" => Some(Level::L3),
            _ => None,
        }
    }
}

fn default_warmup() -> u32 {
    1
}
fn default_order() -> Order {
    Order::Alternate
}
fn default_level() -> Level {
    Level::L2
}
fn default_timeout() -> u64 {
    600
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Measurement {
    /// Measured rounds: every arm once per workload per round.
    pub rounds: u32,
    /// Unmeasured rounds first, to take one-time costs off round one.
    #[serde(default = "default_warmup")]
    pub warmup: u32,
    #[serde(default = "default_order")]
    pub order: Order,
    /// Master seed of the algorithm seeds.  Each (workload, round) gets
    /// one seed shared by every arm, so arms are paired.
    pub seed: u64,
    /// Below this level a wall-clock figure is descriptive, never a result.
    #[serde(default = "default_level")]
    pub isolation_required: Level,
    #[serde(default = "default_timeout")]
    pub timeout_seconds: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Spec {
    pub schema: String,
    pub label: String,
    #[serde(default)]
    pub description: String,
    pub workloads: WorkloadPlan,
    pub arms: Vec<ArmSpec>,
    pub measurement: Measurement,
}

/// One arm with its method resolved.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Arm {
    pub name: String,
    pub role: Role,
    pub method: ResolvedMethod,
}

/// One planned execution.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Execution {
    pub seq: u64,
    pub round: u32,
    pub warmup: bool,
    pub arm: usize,
    pub workload: usize,
    pub algorithm_seed: u64,
}

/// A spec expanded.
pub struct Plan {
    pub spec: Spec,
    /// `ECS1h` + 12 hex of the spec without its label and description.
    pub spec_id: String,
    pub spec_sha256: String,
    pub arms: Vec<Arm>,
    pub workloads: Vec<Workload>,
    /// One built curve per spec curve, in spec order.
    pub instances: Vec<Instance>,
    /// The index into `instances` of each workload's curve.
    pub curve_of: Vec<usize>,
    pub executions: Vec<Execution>,
}

impl Plan {
    /// The built curve workload `w` runs on.
    pub fn instance(&self, w: usize) -> &Instance {
        &self.instances[self.curve_of[w]]
    }
}

impl Spec {
    pub fn from_json(text: &str) -> Result<Self, String> {
        let spec: Spec = serde_json::from_str(text).map_err(|e| format!("spec: {e}"))?;
        if spec.schema != SPEC_SCHEMA {
            return Err(format!(
                "spec schema is `{}`, expected `{SPEC_SCHEMA}`",
                spec.schema
            ));
        }
        Ok(spec)
    }

    fn identity_view(&self) -> Value {
        json!({
            "schema": self.schema,
            "workloads": serde_json::to_value(&self.workloads).unwrap_or(Value::Null),
            "arms": serde_json::to_value(&self.arms).unwrap_or(Value::Null),
            "measurement": serde_json::to_value(&self.measurement).unwrap_or(Value::Null),
        })
    }
}

/// Validate and expand.  Every refusal happens here, before anything runs.
pub fn plan(spec: Spec) -> Result<Plan, String> {
    if spec.arms.is_empty() {
        return Err("a spec needs at least one arm".into());
    }
    if spec.measurement.rounds == 0 {
        return Err("rounds must be at least 1".into());
    }
    if spec.workloads.curves.is_empty() || spec.workloads.targets_per_curve == 0 {
        return Err("a spec needs at least one curve and one target per curve".into());
    }
    let mut names = std::collections::BTreeSet::new();
    let mut arms = Vec::new();
    for a in &spec.arms {
        if !names.insert(a.name.clone()) {
            return Err(format!("arm name `{}` is used twice", a.name));
        }
        if a.name.is_empty()
            || !a
                .name
                .chars()
                .all(|c| c.is_ascii_alphanumeric() || "-_.".contains(c))
        {
            return Err(format!("arm name `{}` must be [A-Za-z0-9._-]+", a.name));
        }
        arms.push(Arm {
            name: a.name.clone(),
            role: a.role,
            method: resolve(&a.method)?,
        });
    }
    let mut workloads = Vec::new();
    let mut instances = Vec::new();
    let mut curve_of = Vec::new();
    for c in &spec.workloads.curves {
        // One build per curve: a prime search counts points in O(p) per
        // candidate, so building per target would multiply that.
        let inst = c.build()?;
        for i in 0..spec.workloads.targets_per_curve {
            let w = Workload::on(c, &inst, spec.workloads.target_seed, i)?;
            for a in &arms {
                let d = decl(&a.method.id).expect("resolved");
                if d.applies == Applies::KoblitzOnly && w.curve.family != "koblitz" {
                    return Err(format!(
                        "arm `{}` ({}) runs on Koblitz curves only, and {} is {}; put it in a spec of its own",
                        a.name, a.method.id, w.curve.slug, w.curve.family
                    ));
                }
            }
            workloads.push(w);
            curve_of.push(instances.len());
        }
        instances.push(inst);
    }
    let m = &spec.measurement;
    let mut executions = Vec::new();
    let mut seq = 0u64;
    for round in 0..(m.warmup + m.rounds) {
        let warmup = round < m.warmup;
        for (wi, w) in workloads.iter().enumerate() {
            // One seed per (workload, round), shared by every arm.
            let wid = u64::from_str_radix(&w.workload_id[1..], 16).unwrap_or(0);
            let algorithm_seed = derive_u64(
                if warmup {
                    "ecbench.warmup_seed"
                } else {
                    "ecbench.algorithm_seed"
                },
                &[m.seed, wid, round as u64],
            );
            let mut order: Vec<usize> = (0..arms.len()).collect();
            match m.order {
                Order::Alternate => {
                    if round % 2 == 1 {
                        order.reverse();
                    }
                }
                Order::Shuffle => {
                    let mut s = derive_u64("ecbench.order", &[m.seed, wid, round as u64]);
                    for i in (1..order.len()).rev() {
                        s = crate::cryptanalysis::ecbench::canonical::splitmix64(s);
                        order.swap(i, (s % (i as u64 + 1)) as usize);
                    }
                }
            }
            for arm in order {
                executions.push(Execution {
                    seq,
                    round,
                    warmup,
                    arm,
                    workload: wi,
                    algorithm_seed,
                });
                seq += 1;
            }
        }
    }
    let (spec_id, spec_sha256) = short_id("ECS1", &spec.identity_view())?;
    Ok(Plan {
        spec,
        spec_id,
        spec_sha256,
        arms,
        workloads,
        instances,
        curve_of,
        executions,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spec(arms: &[(&str, &str)]) -> Spec {
        Spec::from_json(
            &json!({
                "schema": SPEC_SCHEMA,
                "label": "t",
                "workloads": {"curves": [{"kind": "prime_search", "bits": 14, "seed": 59297}],
                              "targets_per_curve": 2, "target_seed": 1},
                "arms": arms.iter().map(|(n, m)| json!({"name": n, "role": "candidate", "method": {"id": m}})).collect::<Vec<_>>(),
                "measurement": {"rounds": 3, "seed": 9},
            })
            .to_string(),
        )
        .unwrap()
    }

    #[test]
    fn plan_pairs_seeds_and_alternates() {
        let p = plan(spec(&[("a", "rho.negation"), ("b", "bsgs.textbook")])).unwrap();
        // (1 warm-up + 3 rounds) × 2 workloads × 2 arms.
        assert_eq!(p.executions.len(), 16);
        let r1: Vec<_> = p
            .executions
            .iter()
            .filter(|e| e.round == 1 && e.workload == 0)
            .collect();
        assert_eq!(r1[0].algorithm_seed, r1[1].algorithm_seed);
        assert_eq!((r1[0].arm, r1[1].arm), (1, 0));
        assert!(p.spec_id.starts_with("ECS1h"));
    }

    #[test]
    fn label_does_not_change_the_spec_id() {
        let a = plan(spec(&[("a", "rho.plain")])).unwrap();
        let mut s = spec(&[("a", "rho.plain")]);
        s.label = "other".into();
        assert_eq!(plan(s).unwrap().spec_id, a.spec_id);
    }

    #[test]
    fn koblitz_only_arm_on_a_prime_curve_is_refused() {
        assert!(plan(spec(&[("f", "rho.signed_frobenius")])).is_err());
        assert!(plan(spec(&[("a", "rho.plain"), ("a", "rho.plain")])).is_err());
    }
}

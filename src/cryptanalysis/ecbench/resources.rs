//! Cold resource vectors beside the counted-operation table.
//!
//! SAT conflicts, hash-table probes and bytes of memory are not constant
//! multiples of a curve addition. Preserve their native units and the
//! whole-solve PMU measurements instead of inventing a conversion to GAE.

use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

use serde::Serialize;

use super::record::Record;
use super::runner::{read_records, read_session};

pub const RESOURCE_SCHEMA: &str = "ecbench.resources/v1";

#[derive(Debug, Serialize)]
pub struct ArmResources {
    pub arm: String,
    pub method: String,
    /// All measured attempts are retained, including failures.
    pub measured: usize,
    pub verified: usize,
    pub statuses: BTreeMap<String, usize>,
    pub lower_bound: bool,
    pub unpriced_units: Vec<String>,
    /// The sum exists only if every measured attempt has a GAE total.
    pub total_gae_sum: Option<f64>,
    /// Fork-to-reap time and CPU time cover all measured attempts.
    pub process_wall_ns_sum: u64,
    pub user_ns_sum: u64,
    pub sys_ns_sum: u64,
    /// In-process solve time covers every method phase, if all runs report it.
    pub solve_wall_ns_sum: Option<u64>,
    /// Peak resident memory of any one measured child, not a sum.
    pub peak_rss_kib: u64,
    /// User-space PMU counts across the solve, present only when every
    /// measured attempt had the corresponding counter available.
    pub instructions_sum: Option<u64>,
    pub cycles_sum: Option<u64>,
    pub pmu_complete_runs: usize,
    pub pmu_errors: Vec<String>,
    /// Cold charged work by named phase; the sum is the arm's charged GAE.
    pub phase_gae_totals: BTreeMap<String, f64>,
    /// Wall totals are descriptive until the session earns its required level.
    /// Coverage states how many measured attempts reported each phase clock.
    pub phase_wall_ns_totals: BTreeMap<String, u64>,
    pub phase_wall_covered_runs: BTreeMap<String, usize>,
    pub online_wall_ns_sum: Option<u64>,
    pub online_complete_runs: usize,
    pub online_phase_ns_totals: BTreeMap<String, u64>,
    /// Integer counters have their original units. Phase names prevent a
    /// `lookups` counter in two phases from being mistaken for one unit.
    pub phase_native_totals: BTreeMap<String, u64>,
    pub method_counter_totals: BTreeMap<String, u64>,
}

#[derive(Debug, Serialize)]
pub struct Resources {
    pub schema: &'static str,
    pub session_id: String,
    pub records_sha256: Option<String>,
    pub scope: &'static str,
    pub arms: Vec<ArmResources>,
}

fn add(target: &mut u64, value: u64, label: &str) -> Result<(), String> {
    *target = target
        .checked_add(value)
        .ok_or_else(|| format!("resource total overflow: {label}"))?;
    Ok(())
}

fn add_map(map: &mut BTreeMap<String, u64>, key: String, value: u64) -> Result<(), String> {
    add(map.entry(key.clone()).or_default(), value, &key)
}

fn add_optional(sum: &mut Option<u64>, value: Option<u64>, label: &str) -> Result<(), String> {
    if let Some(acc) = sum {
        if let Some(value) = value {
            add(acc, value, label)?;
        } else {
            *sum = None;
        }
    }
    Ok(())
}

fn summarize(arm: &str, records: &[&Record]) -> Result<ArmResources, String> {
    let mut out = ArmResources {
        arm: arm.into(),
        method: records[0].method.id.clone(),
        measured: records.len(),
        verified: 0,
        statuses: BTreeMap::new(),
        lower_bound: false,
        unpriced_units: Vec::new(),
        total_gae_sum: Some(0.0),
        process_wall_ns_sum: 0,
        user_ns_sum: 0,
        sys_ns_sum: 0,
        solve_wall_ns_sum: Some(0),
        peak_rss_kib: 0,
        instructions_sum: Some(0),
        cycles_sum: Some(0),
        pmu_complete_runs: 0,
        pmu_errors: Vec::new(),
        phase_gae_totals: BTreeMap::new(),
        phase_wall_ns_totals: BTreeMap::new(),
        phase_wall_covered_runs: BTreeMap::new(),
        online_wall_ns_sum: Some(0),
        online_complete_runs: 0,
        online_phase_ns_totals: BTreeMap::new(),
        phase_native_totals: BTreeMap::new(),
        method_counter_totals: BTreeMap::new(),
    };
    let mut unpriced = BTreeSet::new();
    let mut pmu_errors = BTreeSet::new();
    for r in records {
        out.verified += usize::from(r.counts());
        *out.statuses.entry(r.outcome.status.clone()).or_default() += 1;
        out.lower_bound |= r.cost.lower_bound;
        unpriced.extend(r.cost.unpriced.iter().cloned());
        out.total_gae_sum = out
            .total_gae_sum
            .zip(r.cost.total_gae)
            .map(|(sum, value)| sum + value);
        add(
            &mut out.process_wall_ns_sum,
            r.time.process_wall_ns,
            "process_wall_ns",
        )?;
        add(&mut out.user_ns_sum, r.time.user_ns, "user_ns")?;
        add(&mut out.sys_ns_sum, r.time.sys_ns, "sys_ns")?;
        add_optional(
            &mut out.solve_wall_ns_sum,
            r.time.solve_wall_ns,
            "solve_wall_ns",
        )?;
        out.peak_rss_kib = out.peak_rss_kib.max(r.time.max_rss_kib);
        let hw = r.time.hw_solve.as_ref();
        if hw.and_then(|h| h.instructions).is_some() && hw.and_then(|h| h.cycles).is_some() {
            out.pmu_complete_runs += 1;
        }
        if let Some(error) = hw.and_then(|h| h.error.as_ref()) {
            pmu_errors.insert(error.clone());
        }
        add_optional(
            &mut out.instructions_sum,
            hw.and_then(|h| h.instructions),
            "instructions",
        )?;
        add_optional(&mut out.cycles_sum, hw.and_then(|h| h.cycles), "cycles")?;
        if let Some(window) = &r.online {
            out.online_complete_runs += 1;
            add_optional(
                &mut out.online_wall_ns_sum,
                Some(window.wall_ns),
                "online_wall_ns",
            )?;
            for (name, &value) in &window.phases_ns {
                add_map(&mut out.online_phase_ns_totals, name.clone(), value)?;
            }
        } else {
            out.online_wall_ns_sum = None;
        }
        for phase in &r.phases {
            *out.phase_gae_totals.entry(phase.name.clone()).or_default() += phase.gae;
            if let Some(value) = phase.wall_ns {
                add_map(&mut out.phase_wall_ns_totals, phase.name.clone(), value)?;
                *out.phase_wall_covered_runs
                    .entry(phase.name.clone())
                    .or_default() += 1;
            }
            for (name, &value) in &phase.native {
                add_map(
                    &mut out.phase_native_totals,
                    format!("{}.{}", phase.name, name),
                    value,
                )?;
            }
        }
        for (name, &value) in &r.counters {
            add_map(&mut out.method_counter_totals, name.clone(), value)?;
        }
    }
    out.unpriced_units = unpriced.into_iter().collect();
    out.pmu_errors = pmu_errors.into_iter().collect();
    Ok(out)
}

pub fn session(dir: &Path) -> Result<Resources, String> {
    let session = read_session(dir)?;
    if session.status != "complete" {
        return Err(format!(
            "resource table requires a complete session: {}",
            session.status
        ));
    }
    let records = read_records(dir)?;
    let mut by_arm: BTreeMap<&str, Vec<&Record>> = BTreeMap::new();
    for r in &records {
        if !r.warmup {
            by_arm.entry(r.arm.as_str()).or_default().push(r);
        }
    }
    let mut arms = Vec::with_capacity(by_arm.len());
    for (arm, recs) in by_arm {
        arms.push(summarize(arm, &recs)?);
    }
    Ok(Resources {
        schema: RESOURCE_SCHEMA,
        session_id: session.session_id,
        records_sha256: session.records_sha256,
        scope: "all measured attempts; solve PMU excludes process startup",
        arms,
    })
}

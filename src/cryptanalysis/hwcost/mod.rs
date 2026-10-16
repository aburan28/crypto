//! `hwcost`: a calibrated analytical model of how long a calculation takes
//! on given hardware and what it costs.
//!
//! A [`Hardware`] profile and a [`Workload`] go in; a [`Prediction`] with a
//! time and cost breakdown comes out.  Per phase of the workload:
//!
//! - **Compute.** `iterations × Σ_op (count_op × cycles_op) / clock`, the
//!   cycles per operation taken from the profile (measured on the hot loop by
//!   [`measure`], or entered from instruction tables, llvm-mca or uiCA with
//!   the source named).
//! - **Memory.** `iterations × bytes_per_iteration / bandwidth` of the
//!   first memory level that holds the per-core working set.  The phase takes
//!   `max(compute, memory)`: the roofline, applied to a CPU.
//! - **Parallelism.** Amdahl with a measured scaling efficiency:
//!   `T₁ · ((1 − f) + f / (cores · efficiency(cores)))`.  Efficiency between
//!   measured core counts is interpolated; beyond them it is held at the
//!   last measured value and the prediction says so.
//! - **Transfers.** `messages × latency + bytes / bandwidth` per data
//!   movement; transfers that overlap compute take a `max` with it (per
//!   channel they still add up), the rest add to the wall time.
//! - **Cost.** Wall hours × machines × $/hour, plus egress and storage.
//!
//! Every constant carries its source.  A prediction lists the constants it
//! used that are not measured (sources beginning `assumed` or
//! `illustrative`), so a reader can see which numbers the answer rests on.
//! An operation the workload names but the profile does not price is an
//! error, never a zero.

pub mod fit;
pub mod kernel;
pub mod measure;
pub mod rho;

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

pub const HARDWARE_SCHEMA: &str = "hwcost.hardware.v1";
pub const WORKLOAD_SCHEMA: &str = "hwcost.workload.v1";
pub const PREDICTION_SCHEMA: &str = "hwcost.prediction.v1";

// ── Inputs ──────────────────────────────────────────────────────────

/// The cost of one operation on one core.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct OpCost {
    /// Cycles per operation used by the model (reciprocal throughput when
    /// the hot loop has independent work in flight).
    pub cycles: f64,
    /// Cycles per operation in a dependent chain, when known.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub latency_cycles: Option<f64>,
    pub source: String,
}

/// One level of the memory hierarchy.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct MemLevel {
    pub name: String,
    /// Capacity in bytes; `None` for main memory (unbounded).
    pub capacity_bytes: Option<f64>,
    /// A level shared by all cores of the machine (its capacity is divided
    /// among the cores in use when deciding where a working set fits).
    pub shared: bool,
    /// Read bandwidth one core gets alone, GB/s (10⁹ bytes per second).
    pub read_gbps_one_core: f64,
    /// Aggregate read bandwidth with every core of the profile reading.
    pub read_gbps_all_cores: f64,
    pub source: String,
}

/// Throughput relative to perfect scaling at a core count.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct ScalePoint {
    pub cores: u32,
    pub efficiency: f64,
}

/// A data path: disk, PCIe, network.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Channel {
    pub latency_s: f64,
    pub bandwidth_bytes_per_s: f64,
    pub source: String,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Price {
    pub usd_per_hour: f64,
    #[serde(default)]
    pub egress_usd_per_gb: f64,
    #[serde(default)]
    pub storage_usd_per_gb_month: f64,
    pub source: String,
}

/// A machine type.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Hardware {
    pub schema: String,
    pub name: String,
    /// `cpu` or `gpu`; a GPU profile counts SMs as cores and prices
    /// instruction classes per SM-cycle.
    pub kind: String,
    pub clock_ghz: f64,
    pub clock_source: String,
    pub cores: u32,
    #[serde(default)]
    pub isa: Vec<String>,
    pub ops: BTreeMap<String, OpCost>,
    pub memory: Vec<MemLevel>,
    pub scaling: Vec<ScalePoint>,
    pub scaling_source: String,
    #[serde(default)]
    pub channels: BTreeMap<String, Channel>,
    pub price: Price,
    #[serde(default)]
    pub notes: Vec<String>,
}

/// A count given either plainly or as a power of two (`{"log2": 128}`),
/// so ECDLP-sized work stays readable.  (Deserialised by hand: the crate's
/// `serde_json` has `arbitrary_precision`, which untagged enums mishandle.)
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Count {
    Plain(f64),
    Log2 { log2: f64 },
}

impl Serialize for Count {
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        match *self {
            Count::Plain(v) => s.serialize_f64(v),
            Count::Log2 { log2 } => {
                use serde::ser::SerializeMap;
                let mut m = s.serialize_map(Some(1))?;
                m.serialize_entry("log2", &log2)?;
                m.end()
            }
        }
    }
}

impl<'de> Deserialize<'de> for Count {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        use serde::de::Error;
        let v = serde_json::Value::deserialize(d)?;
        let num = |x: &serde_json::Value| {
            x.as_f64()
                .ok_or_else(|| D::Error::custom("count is not a number"))
        };
        match &v {
            serde_json::Value::Number(_) => Ok(Count::Plain(num(&v)?)),
            serde_json::Value::Object(m) if m.len() == 1 && m.contains_key("log2") => {
                Ok(Count::Log2 {
                    log2: num(&m["log2"])?,
                })
            }
            _ => Err(D::Error::custom("a count is a number or {\"log2\": x}")),
        }
    }
}

impl Count {
    pub fn value(self) -> f64 {
        match self {
            Count::Plain(v) => v,
            Count::Log2 { log2 } => log2.exp2(),
        }
    }
    pub fn log2(self) -> f64 {
        match self {
            Count::Plain(v) => v.log2(),
            Count::Log2 { log2 } => log2,
        }
    }
}

impl Default for Count {
    fn default() -> Self {
        Count::Plain(0.0)
    }
}

fn one() -> f64 {
    1.0
}

fn one_count() -> Count {
    Count::Plain(1.0)
}

/// A loop: `iterations` times the same work.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Phase {
    pub name: String,
    pub iterations: Count,
    /// Operation name (a key of [`Hardware::ops`]) → operations per
    /// iteration.
    pub ops_per_iteration: BTreeMap<String, f64>,
    /// Cycles per iteration not covered by the counted operations (loop
    /// control, partition lookups), when calibrated.
    #[serde(default)]
    pub overhead_cycles_per_iteration: f64,
    #[serde(default)]
    pub overhead_source: String,
    /// Bytes read plus bytes written per iteration.
    pub bytes_per_iteration: f64,
    /// Bytes each core keeps live (decides the memory level).
    pub working_set_bytes_per_core: f64,
    /// Fraction of the phase that runs in parallel (Amdahl).
    #[serde(default = "one")]
    pub parallel_fraction: f64,
}

/// A data movement.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Transfer {
    pub name: String,
    /// A key of [`Hardware::channels`].
    pub channel: String,
    pub bytes: Count,
    /// Number of separate messages or requests (each pays the latency).
    #[serde(default = "one_count")]
    pub messages: Count,
    /// Runs concurrently with compute (takes a `max`, not a sum).
    #[serde(default)]
    pub overlaps_compute: bool,
    /// Each machine moves its share over its own link.
    #[serde(default)]
    pub per_machine_links: bool,
}

/// Data kept at rest.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Storage {
    pub name: String,
    pub bytes: Count,
    /// Months stored; `None` stores it for the predicted wall time.
    #[serde(default)]
    pub months: Option<f64>,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Workload {
    pub schema: String,
    pub name: String,
    pub phases: Vec<Phase>,
    #[serde(default)]
    pub transfers: Vec<Transfer>,
    #[serde(default)]
    pub storage: Vec<Storage>,
    #[serde(default)]
    pub egress_bytes: Count,
    #[serde(default)]
    pub notes: Vec<String>,
}

/// How many machines of the profile, and how many cores of each, run it.
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq)]
pub struct Deployment {
    pub machines: f64,
    pub cores_per_machine: u32,
}

// ── Output ──────────────────────────────────────────────────────────

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct OpShare {
    pub op: String,
    pub per_iteration: f64,
    pub cycles_per_iteration: f64,
    pub share: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct PhaseResult {
    pub name: String,
    pub iterations_log2: f64,
    pub cycles_per_iteration: f64,
    pub ns_per_iteration_compute: f64,
    pub ns_per_iteration_memory: f64,
    pub memory_level: String,
    pub bound: String,
    pub ops: Vec<OpShare>,
    /// One core, alone.
    pub time_one_core_s: f64,
    pub efficiency: f64,
    pub compute_wall_s: f64,
    pub memory_wall_s: f64,
    pub wall_s: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct TransferResult {
    pub name: String,
    pub channel: String,
    pub bytes_log2: f64,
    pub time_s: f64,
    pub overlaps_compute: bool,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct CostResult {
    pub machine_hours: f64,
    pub compute_usd: f64,
    pub egress_usd: f64,
    pub storage_usd: f64,
    pub total_usd: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Prediction {
    pub schema: String,
    pub hardware: String,
    pub workload: String,
    pub deployment: Deployment,
    pub phases: Vec<PhaseResult>,
    pub transfers: Vec<TransferResult>,
    /// Phases plus non-overlapping transfers, or the overlapping
    /// transfers where those take longer.
    pub wall_s: f64,
    pub cost: CostResult,
    /// Unmeasured constants this answer rests on.
    pub assumptions: Vec<String>,
}

// ── The model ───────────────────────────────────────────────────────

fn unmeasured(source: &str) -> bool {
    let s = source.trim_start().to_ascii_lowercase();
    s.starts_with("assumed") || s.starts_with("illustrative")
}

/// Cycles of operation `op` and their source.  `name@latency` prices the
/// operation at its dependent-chain latency, for an operation on the
/// critical path (the one inversion of a batch); a bare name prices it at
/// [`OpCost::cycles`].
pub fn op_cycles<'a>(hw: &'a Hardware, op: &str) -> Result<(f64, &'a str), String> {
    let (base, latency) = match op.strip_suffix("@latency") {
        Some(b) => (b, true),
        None => (op, false),
    };
    let c = hw
        .ops
        .get(base)
        .ok_or_else(|| format!("{} prices no operation `{base}`", hw.name))?;
    if latency {
        let l = c
            .latency_cycles
            .ok_or_else(|| format!("{} has no latency for `{base}`", hw.name))?;
        Ok((l, c.source.as_str()))
    } else {
        Ok((c.cycles, c.source.as_str()))
    }
}

/// Scaling efficiency at `cores`, interpolated between measured points.
/// The flag is set when `cores` lies beyond the last measured point.
pub fn efficiency(hw: &Hardware, cores: u32) -> Result<(f64, bool), String> {
    let mut pts = hw.scaling.clone();
    pts.sort_by_key(|p| p.cores);
    let first = pts
        .first()
        .ok_or("hardware profile has no scaling points")?;
    if cores <= first.cores {
        return Ok((first.efficiency, false));
    }
    for w in pts.windows(2) {
        let (a, b) = (&w[0], &w[1]);
        if cores <= b.cores {
            let t = f64::from(cores - a.cores) / f64::from(b.cores - a.cores);
            return Ok((a.efficiency + t * (b.efficiency - a.efficiency), false));
        }
    }
    Ok((pts.last().unwrap().efficiency, true))
}

/// The first level that holds `working_set` bytes per core when `cores`
/// cores share the machine.
pub fn memory_level(hw: &Hardware, working_set: f64, cores: u32) -> Result<&MemLevel, String> {
    for lvl in &hw.memory {
        match lvl.capacity_bytes {
            None => return Ok(lvl),
            Some(c) => {
                let per_core = if lvl.shared {
                    c / f64::from(cores.max(1))
                } else {
                    c
                };
                if working_set <= per_core {
                    return Ok(lvl);
                }
            }
        }
    }
    Err(format!(
        "no memory level of {} holds {working_set} bytes per core",
        hw.name
    ))
}

/// Aggregate read bandwidth (bytes/s) of a level with `cores` cores
/// reading: linear in cores up to the measured all-core aggregate.
pub fn aggregate_bandwidth(lvl: &MemLevel, cores: u32) -> f64 {
    (lvl.read_gbps_one_core * f64::from(cores)).min(lvl.read_gbps_all_cores) * 1e9
}

pub fn predict(hw: &Hardware, wl: &Workload, dep: Deployment) -> Result<Prediction, String> {
    if hw.schema != HARDWARE_SCHEMA {
        return Err(format!(
            "hardware schema {} is not {HARDWARE_SCHEMA}",
            hw.schema
        ));
    }
    if wl.schema != WORKLOAD_SCHEMA {
        return Err(format!(
            "workload schema {} is not {WORKLOAD_SCHEMA}",
            wl.schema
        ));
    }
    if dep.machines < 1.0 || dep.cores_per_machine == 0 {
        return Err("a deployment needs at least one machine and one core".into());
    }
    if dep.cores_per_machine > hw.cores {
        return Err(format!(
            "{} has {} cores, not {}",
            hw.name, hw.cores, dep.cores_per_machine
        ));
    }
    let hz = hw.clock_ghz * 1e9;
    let mut assumptions: Vec<String> = Vec::new();
    let mut note = |s: String| {
        if !assumptions.contains(&s) {
            assumptions.push(s);
        }
    };
    if unmeasured(&hw.clock_source) {
        note(format!("clock {} GHz: {}", hw.clock_ghz, hw.clock_source));
    }
    let (eff, beyond) = efficiency(hw, dep.cores_per_machine)?;
    if unmeasured(&hw.scaling_source) {
        note(format!("scaling efficiency: {}", hw.scaling_source));
    }
    if beyond {
        note(format!(
            "assumed: efficiency at {} cores held at the last measured point ({eff:.3})",
            dep.cores_per_machine
        ));
    }
    if dep.machines > 1.0 {
        note("assumed: machines scale perfectly with each other (independent work; coordination priced only as transfers)".into());
    }
    let total_cores = dep.machines * f64::from(dep.cores_per_machine);

    let mut phases = Vec::new();
    for ph in &wl.phases {
        let mut cyc = ph.overhead_cycles_per_iteration;
        if ph.overhead_cycles_per_iteration != 0.0 && unmeasured(&ph.overhead_source) {
            note(format!("{} overhead: {}", ph.name, ph.overhead_source));
        }
        let mut ops = Vec::new();
        for (op, &n) in &ph.ops_per_iteration {
            let c = op_cycles(hw, op).map_err(|e| format!("{e} (phase {})", ph.name))?;
            if unmeasured(&c.1) {
                note(format!("{op} = {} cycles: {}", c.0, c.1));
            }
            cyc += n * c.0;
            ops.push(OpShare {
                op: op.clone(),
                per_iteration: n,
                cycles_per_iteration: n * c.0,
                share: 0.0,
            });
        }
        for o in &mut ops {
            o.share = if cyc > 0.0 {
                o.cycles_per_iteration / cyc
            } else {
                0.0
            };
        }
        let iters = ph.iterations.value();
        let lvl = memory_level(hw, ph.working_set_bytes_per_core, dep.cores_per_machine)?;
        if unmeasured(&lvl.source) {
            note(format!("{} bandwidth: {}", lvl.name, lvl.source));
        }
        let bw_one = lvl.read_gbps_one_core * 1e9;
        let bw_machine = aggregate_bandwidth(lvl, dep.cores_per_machine);
        let ns_comp = cyc / hz * 1e9;
        let ns_mem = ph.bytes_per_iteration / bw_one * 1e9;
        let t1 = iters * (ns_comp.max(ns_mem)) * 1e-9;
        let f = ph.parallel_fraction.clamp(0.0, 1.0);
        let t_comp1 = iters * ns_comp * 1e-9;
        let compute_wall = t_comp1 * ((1.0 - f) + f / (total_cores * eff));
        let memory_wall = iters * ph.bytes_per_iteration / (bw_machine * dep.machines);
        let (wall, bound) = if memory_wall > compute_wall {
            (memory_wall, "memory")
        } else {
            (compute_wall, "compute")
        };
        phases.push(PhaseResult {
            name: ph.name.clone(),
            iterations_log2: ph.iterations.log2(),
            cycles_per_iteration: cyc,
            ns_per_iteration_compute: ns_comp,
            ns_per_iteration_memory: ns_mem,
            memory_level: lvl.name.clone(),
            bound: bound.into(),
            ops,
            time_one_core_s: t1,
            efficiency: eff,
            compute_wall_s: compute_wall,
            memory_wall_s: memory_wall,
            wall_s: wall,
        });
    }
    let phase_s: f64 = phases.iter().map(|p| p.wall_s).sum();

    let mut transfers = Vec::new();
    let mut serial_s = 0.0;
    let mut overlapped: BTreeMap<String, f64> = BTreeMap::new();
    for t in &wl.transfers {
        let ch = hw.channels.get(&t.channel).ok_or_else(|| {
            format!(
                "{} has no channel `{}` (transfer {})",
                hw.name, t.channel, t.name
            )
        })?;
        if unmeasured(&ch.source) {
            note(format!("{} channel: {}", t.channel, ch.source));
        }
        let links = if t.per_machine_links {
            dep.machines
        } else {
            1.0
        };
        let time = (t.messages.value() * ch.latency_s + t.bytes.value() / ch.bandwidth_bytes_per_s)
            / links;
        if t.overlaps_compute {
            *overlapped.entry(t.channel.clone()).or_default() += time;
        } else {
            serial_s += time;
        }
        transfers.push(TransferResult {
            name: t.name.clone(),
            channel: t.channel.clone(),
            bytes_log2: t.bytes.log2(),
            time_s: time,
            overlaps_compute: t.overlaps_compute,
        });
    }
    let overlap_s = overlapped.values().cloned().fold(0.0, f64::max);
    let wall_s = phase_s.max(overlap_s) + serial_s;

    if unmeasured(&hw.price.source) {
        note(format!("price: {}", hw.price.source));
    }
    let machine_hours = wall_s / 3600.0 * dep.machines;
    let compute_usd = machine_hours * hw.price.usd_per_hour;
    let egress_usd = wl.egress_bytes.value() / 1e9 * hw.price.egress_usd_per_gb;
    let months = wall_s / (3600.0 * 24.0 * 365.25 / 12.0);
    let storage_usd: f64 = wl
        .storage
        .iter()
        .map(|s| {
            s.bytes.value() / 1e9 * s.months.unwrap_or(months) * hw.price.storage_usd_per_gb_month
        })
        .sum();

    Ok(Prediction {
        schema: PREDICTION_SCHEMA.into(),
        hardware: hw.name.clone(),
        workload: wl.name.clone(),
        deployment: dep,
        phases,
        transfers,
        wall_s,
        cost: CostResult {
            machine_hours,
            compute_usd,
            egress_usd,
            storage_usd,
            total_usd: compute_usd + egress_usd + storage_usd,
        },
        assumptions,
    })
}

/// A duration in the largest unit that keeps it readable.
pub fn human_duration(s: f64) -> String {
    const MIN: f64 = 60.0;
    const HOUR: f64 = 3600.0;
    const DAY: f64 = 86400.0;
    const YEAR: f64 = 365.25 * DAY;
    if !s.is_finite() {
        return format!("{s}");
    }
    if s < 1e-3 {
        format!("{:.3} µs", s * 1e6)
    } else if s < 1.0 {
        format!("{:.3} ms", s * 1e3)
    } else if s < MIN {
        format!("{s:.3} s")
    } else if s < HOUR {
        format!("{:.2} min", s / MIN)
    } else if s < DAY {
        format!("{:.2} h", s / HOUR)
    } else if s < YEAR {
        format!("{:.2} days", s / DAY)
    } else if s < 1e4 * YEAR {
        format!("{:.2} years", s / YEAR)
    } else {
        format!("2^{:.1} years", (s / YEAR).log2())
    }
}

/// Dollars, with powers of two once they stop being readable.
pub fn human_usd(x: f64) -> String {
    if x < 1e9 {
        format!("${x:.2}")
    } else {
        format!("$2^{:.1}", x.log2())
    }
}

/// The prediction as a text table.
pub fn render(p: &Prediction) -> String {
    use std::fmt::Write;
    let mut s = String::new();
    let _ = writeln!(
        s,
        "{} on {} × {} ({} cores each)",
        p.workload, p.deployment.machines, p.hardware, p.deployment.cores_per_machine
    );
    for ph in &p.phases {
        let _ = writeln!(
            s,
            "  phase {}: 2^{:.2} iterations, {:.1} cycles/iter ({:.2} ns compute, {:.2} ns {} memory per iter, {}-bound), efficiency {:.3}",
            ph.name,
            ph.iterations_log2,
            ph.cycles_per_iteration,
            ph.ns_per_iteration_compute,
            ph.ns_per_iteration_memory,
            ph.memory_level,
            ph.bound,
            ph.efficiency
        );
        for o in &ph.ops {
            let _ = writeln!(
                s,
                "      {:<14} {:>9.4}/iter {:>10.1} cycles {:>6.1}%",
                o.op,
                o.per_iteration,
                o.cycles_per_iteration,
                100.0 * o.share
            );
        }
        let _ = writeln!(
            s,
            "      one core {}; deployment compute {}, memory {}",
            human_duration(ph.time_one_core_s),
            human_duration(ph.compute_wall_s),
            human_duration(ph.memory_wall_s)
        );
    }
    for t in &p.transfers {
        let _ = writeln!(
            s,
            "  transfer {} over {}: 2^{:.1} bytes, {}{}",
            t.name,
            t.channel,
            t.bytes_log2,
            human_duration(t.time_s),
            if t.overlaps_compute {
                " (overlaps compute)"
            } else {
                ""
            }
        );
    }
    let _ = writeln!(s, "  wall time {}", human_duration(p.wall_s));
    let _ = writeln!(
        s,
        "  cost {} = compute {} ({:.3e} machine-hours) + egress {} + storage {}",
        human_usd(p.cost.total_usd),
        human_usd(p.cost.compute_usd),
        p.cost.machine_hours,
        human_usd(p.cost.egress_usd),
        human_usd(p.cost.storage_usd)
    );
    if !p.assumptions.is_empty() {
        let _ = writeln!(s, "  rests on unmeasured constants:");
        for a in &p.assumptions {
            let _ = writeln!(s, "    - {a}");
        }
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hw() -> Hardware {
        let mut ops = BTreeMap::new();
        ops.insert(
            "mul".into(),
            OpCost {
                cycles: 10.0,
                latency_cycles: None,
                source: "measured".into(),
            },
        );
        ops.insert(
            "inv".into(),
            OpCost {
                cycles: 1000.0,
                latency_cycles: None,
                source: "assumed".into(),
            },
        );
        let mut channels = BTreeMap::new();
        channels.insert(
            "network".into(),
            Channel {
                latency_s: 1e-3,
                bandwidth_bytes_per_s: 1e8,
                source: "assumed: 1 Gbit/s".into(),
            },
        );
        Hardware {
            schema: HARDWARE_SCHEMA.into(),
            name: "toy".into(),
            kind: "cpu".into(),
            clock_ghz: 1.0,
            clock_source: "measured".into(),
            cores: 4,
            isa: vec![],
            ops,
            memory: vec![
                MemLevel {
                    name: "L1".into(),
                    capacity_bytes: Some(1000.0),
                    shared: false,
                    read_gbps_one_core: 100.0,
                    read_gbps_all_cores: 400.0,
                    source: "measured".into(),
                },
                MemLevel {
                    name: "DRAM".into(),
                    capacity_bytes: None,
                    shared: true,
                    read_gbps_one_core: 1.0,
                    read_gbps_all_cores: 2.0,
                    source: "measured".into(),
                },
            ],
            scaling: vec![
                ScalePoint {
                    cores: 1,
                    efficiency: 1.0,
                },
                ScalePoint {
                    cores: 4,
                    efficiency: 0.5,
                },
            ],
            scaling_source: "measured".into(),
            channels,
            price: Price {
                usd_per_hour: 3.6,
                egress_usd_per_gb: 0.0,
                storage_usd_per_gb_month: 0.0,
                source: "illustrative".into(),
            },
            notes: vec![],
        }
    }

    fn wl(ws: f64) -> Workload {
        let mut ops = BTreeMap::new();
        ops.insert("mul".into(), 5.0);
        ops.insert("inv".into(), 0.01);
        Workload {
            schema: WORKLOAD_SCHEMA.into(),
            name: "w".into(),
            phases: vec![Phase {
                name: "walk".into(),
                iterations: Count::Plain(1e9),
                ops_per_iteration: ops,
                overhead_cycles_per_iteration: 0.0,
                overhead_source: String::new(),
                bytes_per_iteration: 100.0,
                working_set_bytes_per_core: ws,
                parallel_fraction: 1.0,
            }],
            transfers: vec![],
            storage: vec![],
            egress_bytes: Count::Plain(0.0),
            notes: vec![],
        }
    }

    #[test]
    fn compute_bound_phase_is_ops_times_cycles_over_clock() {
        let p = predict(
            &hw(),
            &wl(500.0),
            Deployment {
                machines: 1.0,
                cores_per_machine: 1,
            },
        )
        .unwrap();
        let ph = &p.phases[0];
        assert_eq!(ph.cycles_per_iteration, 60.0);
        assert_eq!(ph.bound, "compute");
        assert!((p.wall_s - 60.0).abs() < 1e-9);
        assert!((p.cost.compute_usd - 0.06).abs() < 1e-12);
        assert!(p.assumptions.iter().any(|a| a.starts_with("inv")));
    }

    #[test]
    fn working_set_outside_cache_turns_memory_bound() {
        let p = predict(
            &hw(),
            &wl(1e6),
            Deployment {
                machines: 1.0,
                cores_per_machine: 1,
            },
        )
        .unwrap();
        let ph = &p.phases[0];
        assert_eq!(ph.memory_level, "DRAM");
        assert_eq!(ph.bound, "memory");
        assert!((p.wall_s - 100.0).abs() < 1e-9);
    }

    #[test]
    fn scaling_interpolates_and_flags_extrapolation() {
        let h = hw();
        assert_eq!(efficiency(&h, 1).unwrap(), (1.0, false));
        let (e, beyond) = efficiency(&h, 2).unwrap();
        assert!((e - (1.0 - 0.5 / 3.0)).abs() < 1e-12 && !beyond);
        let p = predict(
            &h,
            &wl(500.0),
            Deployment {
                machines: 1.0,
                cores_per_machine: 4,
            },
        )
        .unwrap();
        assert!((p.wall_s - 60.0 / (4.0 * 0.5)).abs() < 1e-9);
    }

    #[test]
    fn overlapping_transfers_take_a_max_and_serial_ones_add() {
        let mut w = wl(500.0);
        w.transfers.push(Transfer {
            name: "up".into(),
            channel: "network".into(),
            bytes: Count::Log2 { log2: 30.0 },
            messages: Count::Plain(1000.0),
            overlaps_compute: true,
            per_machine_links: false,
        });
        let p = predict(
            &hw(),
            &w,
            Deployment {
                machines: 1.0,
                cores_per_machine: 1,
            },
        )
        .unwrap();
        let t = 1000.0 * 1e-3 + 2f64.powi(30) / 1e8;
        assert!((p.transfers[0].time_s - t).abs() < 1e-9);
        assert!(
            (p.wall_s - 60.0).abs() < 1e-9,
            "hidden behind 60 s of compute"
        );
        w.transfers[0].overlaps_compute = false;
        let p = predict(
            &hw(),
            &w,
            Deployment {
                machines: 1.0,
                cores_per_machine: 1,
            },
        )
        .unwrap();
        assert!((p.wall_s - (60.0 + t)).abs() < 1e-9);
    }

    #[test]
    fn unpriced_operation_is_an_error() {
        let mut w = wl(500.0);
        w.phases[0].ops_per_iteration.insert("sqr".into(), 1.0);
        let e = predict(
            &hw(),
            &w,
            Deployment {
                machines: 1.0,
                cores_per_machine: 1,
            },
        )
        .unwrap_err();
        assert!(e.contains("`sqr`"));
    }

    #[test]
    fn counts_parse_plain_or_log2() {
        let c: Count = serde_json::from_str(r#"{"log2": 10}"#).unwrap();
        assert_eq!(c.value(), 1024.0);
        let c: Count = serde_json::from_str("12.5").unwrap();
        assert_eq!(c.value(), 12.5);
    }
}

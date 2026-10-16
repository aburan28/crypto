//! From a [`Measurements`] run to a calibrated [`Hardware`] profile, and the
//! validation of the model against the rho kernel it was measured on.
//!
//! Three predictions of the kernel's ns per iteration are compared with
//! its measured median at every field and batch size:
//!
//! - **A, a priori** (no fitted parameter): the counted operations priced at
//!   their microbenchmarked throughput, the batch inversion and its chained
//!   multiplications at their latency, under the measured memory roofline;
//! - **L, all-latency** (descriptive): every operation at its
//!   dependent-chain latency, as if no two operations overlapped;
//! - **C, calibrated**: A plus one per-iteration residual per field (the
//!   uncounted loads, branches and loop control, less any serial latency A
//!   over-counts; it may be negative), fitted
//!   on the small batches (`B ≤ fit_max_batch`) and tested on the larger ones;
//! - **S, size-extrapolated** (largest field only): A and the residual with
//!   every per-operation cost fitted against the limb count on the smaller
//!   `FpWide` fields and extrapolated.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use super::kernel::{ELEMS_MOVED_PER_ITERATION, ELEMS_RESIDENT_PER_WALK, TABLE};
use super::measure::{field_spec, KernelCounts, Measurements, FIELDS};
use super::{Channel, Hardware, MemLevel, OpCost, Price, ScalePoint, HARDWARE_SCHEMA};

pub const VALIDATION_SCHEMA: &str = "hwcost.validation.v1";

/// The base operations of the kernel and how their cost is extrapolated in
/// the limb count `n`: `α + γ·n²` for the products and the inversion,
/// `α + β·n` for additions and the per-iteration residual.
pub const OPS: [&str; 4] = ["mul", "sqr", "add", "inv"];

fn quadratic(op: &str) -> bool {
    matches!(op, "mul" | "sqr" | "inv")
}

/// Least squares `y = a + b·x`.
pub fn linfit(xs: &[f64], ys: &[f64]) -> Option<(f64, f64)> {
    let n = xs.len() as f64;
    if xs.len() < 2 {
        return None;
    }
    let mx = xs.iter().sum::<f64>() / n;
    let my = ys.iter().sum::<f64>() / n;
    let sxx: f64 = xs.iter().map(|x| (x - mx) * (x - mx)).sum();
    if sxx == 0.0 {
        return None;
    }
    let sxy: f64 = xs.iter().zip(ys).map(|(x, y)| (x - mx) * (y - my)).sum();
    let b = sxy / sxx;
    Some((my - b * mx, b))
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Row {
    pub field: String,
    pub batch: usize,
    /// `fit` or `test` for the calibrated model.
    pub set: String,
    pub measured_ns: f64,
    pub measured_q1: f64,
    pub measured_q3: f64,
    pub working_set_bytes: f64,
    pub memory_level: String,
    pub compute_ns: f64,
    pub memory_ns: f64,
    pub bound: String,
    pub a_priori_ns: f64,
    pub latency_ns: f64,
    pub calibrated_ns: f64,
    pub size_extrapolated_ns: Option<f64>,
    pub err_a_priori: f64,
    pub err_latency: f64,
    pub err_calibrated: f64,
    pub err_size_extrapolated: Option<f64>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Claim {
    pub id: String,
    pub statement: String,
    pub value: f64,
    pub threshold: f64,
    pub verdict: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SizeFit {
    pub op: String,
    pub model: String,
    pub alpha: f64,
    pub slope: f64,
    pub extrapolated_ns: f64,
    pub measured_ns: Option<f64>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Validation {
    pub schema: String,
    pub fit_max_batch: usize,
    pub clock_ghz: f64,
    /// The `imul`-chain cross-check of the clock.
    pub clock_imul_ghz: Option<f64>,
    /// Fitted per-iteration residual (measured − A over the fit set), ns.
    pub residual_ns: BTreeMap<String, f64>,
    /// Measured efficiency per field and thread count.
    pub scaling: BTreeMap<String, BTreeMap<usize, f64>>,
    pub size_fits: Vec<SizeFit>,
    pub size_target: Option<String>,
    pub rows: Vec<Row>,
    pub claims: Vec<Claim>,
}

fn arm_ns(m: &Measurements, field: &str, op: &str, mode: &str) -> Result<f64, String> {
    m.median(&format!("{field}.{op}.{mode}"))
        .ok_or_else(|| format!("no samples for {field}.{op}.{mode}"))
}

pub fn clock_ghz(m: &Measurements) -> (f64, String) {
    match m.median("clock.add_chain") {
        Some(g) => (
            g,
            format!(
                "measured: dependent add chain, median of {} reps",
                m.params.reps
            ),
        ),
        None => (
            1.0,
            "assumed: no clock measurement; cycles are nanoseconds".into(),
        ),
    }
}

/// Memory levels from the cache topology and the bandwidth sweep: each
/// level takes the largest swept size that fills at most half of it
/// (per thread for private levels; summed over threads for shared ones).
pub fn memory_levels(m: &Measurements) -> Vec<MemLevel> {
    let all = *m.params.threads.last().unwrap_or(&1);
    let bw = |bytes: usize, t: usize| m.median(&format!("mem.{t}t.{bytes}"));
    let mut sizes = m.params.mem_sizes.clone();
    sizes.sort();
    let mut levels = Vec::new();
    let mut data_caches: Vec<_> = m
        .caches
        .iter()
        .filter(|c| c.kind != "Instruction")
        .collect();
    data_caches.sort_by_key(|c| c.level);
    for c in data_caches {
        let shared = !c.shared_cpu_list.is_empty()
            && (c.shared_cpu_list.contains('-') || c.shared_cpu_list.contains(','));
        let cap = c.size_bytes as f64;
        let one = sizes.iter().rev().find(|&&s| (s as f64) <= cap / 2.0);
        let per_thread_all = if shared {
            cap / 2.0 / all as f64
        } else {
            cap / 2.0
        };
        let multi = sizes.iter().rev().find(|&&s| (s as f64) <= per_thread_all);
        let (Some(&s1), Some(&sa)) = (one, multi) else {
            continue;
        };
        let (Some(b1), Some(ba)) = (bw(s1, 1), bw(sa, all)) else {
            continue;
        };
        levels.push(MemLevel {
            name: format!("L{}{}", c.level, if c.kind == "Data" { "d" } else { "" }),
            capacity_bytes: Some(cap),
            shared,
            read_gbps_one_core: b1,
            read_gbps_all_cores: ba,
            source: format!(
                "measured: sequential read, {s1} B on 1 thread / {sa} B per thread on {all} (median)"
            ),
        });
    }
    if let Some(&big) = sizes.last() {
        if let (Some(b1), Some(ba)) = (bw(big, 1), bw(big, all)) {
            levels.push(MemLevel {
                name: "DRAM".into(),
                capacity_bytes: None,
                shared: true,
                read_gbps_one_core: b1,
                read_gbps_all_cores: ba,
                source: format!("measured: sequential read, {big} B per thread (median)"),
            });
        }
    }
    levels
}

/// Efficiency `ns₁ / (T · ns_T)` per field and thread count.
pub fn scaling(m: &Measurements) -> BTreeMap<String, BTreeMap<usize, f64>> {
    let mut out = BTreeMap::new();
    for f in &m.params.fields {
        let Some(one) = m.median(&format!("scaling.{f}.1t")) else {
            continue;
        };
        let mut row = BTreeMap::new();
        for &t in &m.params.threads {
            if let Some(nt) = m.median(&format!("scaling.{f}.{t}t")) {
                row.insert(t, one / (t as f64 * nt));
            }
        }
        out.insert(f.clone(), row);
    }
    out
}

fn counts_of<'a>(m: &'a Measurements, f: &str, b: usize) -> Option<&'a KernelCounts> {
    m.counts.iter().find(|c| c.field == f && c.batch == b)
}

/// Compute ns per iteration from counts and per-op costs (ns): the
/// batch's serial chain (chained multiplications, the inversion) at
/// latency, the rest at throughput, or everything at latency.
fn compute_ns(
    c: &KernelCounts,
    thr: &dyn Fn(&str) -> f64,
    lat: &dyn Fn(&str) -> f64,
    all_latency: bool,
) -> f64 {
    let p = |op: &str| if all_latency { lat(op) } else { thr(op) };
    c.mul * p("mul")
        + c.mul_dep * lat("mul")
        + c.sqr * p("sqr")
        + c.add * p("add")
        + c.inv * lat("inv")
}

fn level_for(levels: &[MemLevel], ws: f64) -> Option<&MemLevel> {
    // Validation runs one thread, so a shared level is all its own.
    levels
        .iter()
        .find(|l| l.capacity_bytes.is_none_or(|c| ws <= c))
}

pub fn validate(m: &Measurements, fit_max_batch: usize) -> Result<(Validation, Hardware), String> {
    let (ghz, clock_src) = clock_ghz(m);
    let cores = *m.params.threads.last().unwrap_or(&1) as u32;
    let levels = memory_levels(m);
    if levels.is_empty() {
        return Err("no memory level could be priced from the bandwidth sweep".into());
    }
    let scal = scaling(m);

    // Per-op ns, throughput and latency.
    let mut thr: BTreeMap<(String, String), f64> = BTreeMap::new();
    let mut lat: BTreeMap<(String, String), f64> = BTreeMap::new();
    for f in &m.params.fields {
        for op in OPS {
            thr.insert((f.clone(), op.into()), arm_ns(m, f, op, "throughput")?);
            lat.insert((f.clone(), op.into()), arm_ns(m, f, op, "latency")?);
        }
    }

    // Model A, L and the per-field residual fit (model C).
    struct Pre {
        field: String,
        batch: usize,
        measured: (f64, f64, f64),
        ws: f64,
        level: String,
        comp: f64,
        mem: f64,
        lat: f64,
    }
    let mut pre = Vec::new();
    for f in &m.params.fields {
        let spec = field_spec(f)?;
        let e = spec.elem_bytes() as f64;
        for &b in &m.params.batches {
            let arm = format!("kernel.{f}.{b}");
            let Some(meas) = m.median(&arm) else { continue };
            let (q1, q3) = m.quartiles(&arm).unwrap();
            let c = counts_of(m, f, b).ok_or_else(|| format!("no counts for {f} B={b}"))?;
            let t = |op: &str| thr[&(f.clone(), op.to_string())];
            let l = |op: &str| lat[&(f.clone(), op.to_string())];
            let comp = compute_ns(c, &t, &l, false);
            let comp_lat = compute_ns(c, &t, &l, true);
            let ws = ELEMS_RESIDENT_PER_WALK as f64 * e * b as f64 + 2.0 * e * TABLE as f64;
            let lvl = level_for(&levels, ws).ok_or("no memory level")?;
            let mem = ELEMS_MOVED_PER_ITERATION as f64 * e / lvl.read_gbps_one_core;
            pre.push(Pre {
                field: f.clone(),
                batch: b,
                measured: (meas, q1, q3),
                ws,
                level: lvl.name.clone(),
                comp,
                mem,
                lat: comp_lat.max(mem),
            });
        }
    }
    let mut residual: BTreeMap<String, f64> = BTreeMap::new();
    for f in &m.params.fields {
        let res: Vec<f64> = pre
            .iter()
            .filter(|p| &p.field == f && p.batch <= fit_max_batch)
            .map(|p| p.measured.0 - p.comp.max(p.mem))
            .collect();
        if res.is_empty() {
            return Err(format!(
                "{f}: no batch ≤ {fit_max_batch} to fit the residual on"
            ));
        }
        residual.insert(f.clone(), res.iter().sum::<f64>() / res.len() as f64);
    }

    // Model S: per-op costs and residual against limbs on the FpWide fields
    // below the largest, extrapolated to it.
    let wide: Vec<_> = FIELDS
        .iter()
        .filter(|s| s.name.starts_with('w') && m.params.fields.iter().any(|f| f == s.name))
        .collect();
    let mut size_fits = Vec::new();
    let mut size_target = None;
    let mut s_ns: BTreeMap<usize, f64> = BTreeMap::new();
    if wide.len() >= 3 {
        let target = wide.last().unwrap();
        let train = &wide[..wide.len() - 1];
        let x_of = |n: usize, quad: bool| if quad { (n * n) as f64 } else { n as f64 };
        let mut ext_thr: BTreeMap<String, f64> = BTreeMap::new();
        let mut ext_lat: BTreeMap<String, f64> = BTreeMap::new();
        for op in OPS {
            for (mode, table, out) in [
                ("throughput", &thr, &mut ext_thr),
                ("latency", &lat, &mut ext_lat),
            ] {
                let xs: Vec<f64> = train.iter().map(|s| x_of(s.limbs, quadratic(op))).collect();
                let ys: Vec<f64> = train
                    .iter()
                    .map(|s| table[&(s.name.to_string(), op.to_string())])
                    .collect();
                let (a, g) = linfit(&xs, &ys).ok_or("degenerate size fit")?;
                let v = a + g * x_of(target.limbs, quadratic(op));
                out.insert(op.to_string(), v);
                size_fits.push(SizeFit {
                    op: format!("{op}.{mode}"),
                    model: if quadratic(op) {
                        "a + g*n^2"
                    } else {
                        "a + b*n"
                    }
                    .into(),
                    alpha: a,
                    slope: g,
                    extrapolated_ns: v,
                    measured_ns: table
                        .get(&(target.name.to_string(), op.to_string()))
                        .copied(),
                });
            }
        }
        let xs: Vec<f64> = train.iter().map(|s| s.limbs as f64).collect();
        let ys: Vec<f64> = train.iter().map(|s| residual[s.name]).collect();
        let (a, g) = linfit(&xs, &ys).ok_or("degenerate residual fit")?;
        let ov = a + g * target.limbs as f64;
        size_fits.push(SizeFit {
            op: "residual".into(),
            model: "a + b*n".into(),
            alpha: a,
            slope: g,
            extrapolated_ns: ov,
            measured_ns: residual.get(target.name).copied(),
        });
        for p in pre.iter().filter(|p| p.field == target.name) {
            let c = counts_of(m, &p.field, p.batch).unwrap();
            let comp = compute_ns(c, &|op| ext_thr[op], &|op| ext_lat[op], false);
            s_ns.insert(p.batch, comp.max(p.mem) + ov);
        }
        size_target = Some(target.name.to_string());
    }

    let rows: Vec<Row> = pre
        .iter()
        .map(|p| {
            let a = p.comp.max(p.mem);
            let c = a + residual[&p.field];
            let s = if size_target.as_deref() == Some(p.field.as_str()) {
                s_ns.get(&p.batch).copied()
            } else {
                None
            };
            let err = |x: f64| x / p.measured.0 - 1.0;
            Row {
                field: p.field.clone(),
                batch: p.batch,
                set: if p.batch <= fit_max_batch {
                    "fit"
                } else {
                    "test"
                }
                .into(),
                measured_ns: p.measured.0,
                measured_q1: p.measured.1,
                measured_q3: p.measured.2,
                working_set_bytes: p.ws,
                memory_level: p.level.clone(),
                compute_ns: p.comp,
                memory_ns: p.mem,
                bound: if p.mem > p.comp { "memory" } else { "compute" }.into(),
                a_priori_ns: a,
                latency_ns: p.lat,
                calibrated_ns: c,
                size_extrapolated_ns: s,
                err_a_priori: err(a),
                err_latency: err(p.lat),
                err_calibrated: err(c),
                err_size_extrapolated: s.map(err),
            }
        })
        .collect();

    let mut claims = claims(&rows, &m.params.fields, size_target.as_deref());
    if let Some(g) = m.median("clock.imul_chain") {
        let dev = (g / ghz - 1.0).abs();
        claims.insert(
            0,
            Claim {
                id: "C0".into(),
                statement: "clock: add-chain and imul-chain estimates agree, |relative difference|"
                    .into(),
                value: dev,
                threshold: 0.05,
                verdict: if dev <= 0.05 { "PASS" } else { "FAIL" }.into(),
            },
        );
    }

    // The profile.
    let mut ops = BTreeMap::new();
    let src = |what: &str| {
        format!(
            "measured: hwcost measure, {what}, median of {} reps, at {ghz:.3} GHz",
            m.params.reps
        )
    };
    for f in &m.params.fields {
        for op in OPS {
            let t = thr[&(f.clone(), op.to_string())];
            let l = lat[&(f.clone(), op.to_string())];
            ops.insert(
                format!("{f}.{op}"),
                OpCost {
                    cycles: t * ghz,
                    latency_cycles: Some(l * ghz),
                    source: src(&format!(
                        "{f}.{op} throughput (8 chains; inv 4) and latency"
                    )),
                },
            );
        }
        ops.insert(
            format!("{f}.rho_residual"),
            OpCost {
                cycles: residual[f] * ghz,
                latency_cycles: None,
                source: format!(
                    "measured: fitted per-iteration residual (measured − a priori) of the rho kernel over batches ≤ {fit_max_batch} (hwcost fit)"
                ),
            },
        );
    }
    let mut eff_by_t: BTreeMap<usize, Vec<f64>> = BTreeMap::new();
    for row in scal.values() {
        for (&t, &e) in row {
            eff_by_t.entry(t).or_default().push(e);
        }
    }
    let scaling_pts: Vec<ScalePoint> = eff_by_t
        .into_iter()
        .map(|(t, mut v)| ScalePoint {
            cores: t as u32,
            efficiency: super::measure::median(&mut v).unwrap_or(1.0),
        })
        .collect();
    let host = &m.host;
    let cpu = host
        .pointer("/stable/cpu_model")
        .and_then(|v| v.as_str())
        .unwrap_or("unknown CPU")
        .to_string();
    let isa: Vec<String> = host
        .pointer("/stable/features")
        .and_then(|v| v.as_array())
        .map(|a| {
            a.iter()
                .filter_map(|x| x.as_str().map(String::from))
                .collect()
        })
        .unwrap_or_default();
    let mut channels = BTreeMap::new();
    channels.insert(
        "network".into(),
        Channel {
            latency_s: 1e-3,
            bandwidth_bytes_per_s: 1.25e9,
            source: "assumed: 10 Gbit/s, 1 ms per message (not measured)".into(),
        },
    );
    channels.insert(
        "disk".into(),
        Channel {
            latency_s: 1e-3,
            bandwidth_bytes_per_s: 125e6,
            source: "assumed: block-storage baseline 125 MB/s, 1 ms per request (not measured)"
                .into(),
        },
    );
    let hw = Hardware {
        schema: HARDWARE_SCHEMA.into(),
        name: format!("{cpu}, {cores} logical CPUs"),
        kind: "cpu".into(),
        clock_ghz: ghz,
        clock_source: clock_src,
        cores,
        isa,
        ops,
        memory: levels,
        scaling: scaling_pts,
        scaling_source: format!(
            "measured: rho kernel, B = {} per thread, median efficiency over fields {:?}",
            m.params.scaling_batch, m.params.fields
        ),
        channels,
        price: Price {
            usd_per_hour: 0.1785,
            egress_usd_per_gb: 0.09,
            storage_usd_per_gb_month: 0.023,
            source: "illustrative: AWS us-east-1 on-demand c7i.xlarge (4 vCPU) $0.1785/h, internet egress $0.09/GB, S3 standard $0.023/GB-month; this host's own cost is not priced".into(),
        },
        notes: vec![
            format!("compiled target features: {:?}", m.compiled_features),
            "op costs are for the crate's fp_backend fields as compiled; another ISA build (target-cpu=native, ADX/MULX, AVX-512 IFMA) needs its own measurement".into(),
        ],
    };
    Ok((
        Validation {
            schema: VALIDATION_SCHEMA.into(),
            fit_max_batch,
            clock_ghz: ghz,
            clock_imul_ghz: m.median("clock.imul_chain"),
            residual_ns: residual,
            scaling: scal,
            size_fits,
            size_target,
            rows,
            claims,
        },
        hw,
    ))
}

fn median_abs(v: impl Iterator<Item = f64>) -> f64 {
    let mut v: Vec<f64> = v.map(f64::abs).collect();
    super::measure::median(&mut v).unwrap_or(f64::NAN)
}

/// The claims frozen in the study protocol.
fn claims(rows: &[Row], fields: &[String], size_target: Option<&str>) -> Vec<Claim> {
    let mut out = Vec::new();
    let le = |id: String, statement: String, value: f64, threshold: f64| Claim {
        id,
        statement,
        value,
        threshold,
        verdict: if value.is_finite() && value <= threshold {
            "PASS"
        } else {
            "FAIL"
        }
        .into(),
    };
    for f in fields {
        let test = rows.iter().filter(|r| &r.field == f && r.set == "test");
        out.push(le(
            format!("V1.{f}"),
            format!("{f}: calibrated model, test batches, median |error|"),
            median_abs(test.map(|r| r.err_calibrated)),
            0.30,
        ));
    }
    for f in fields {
        let all = rows.iter().filter(|r| &r.field == f);
        out.push(le(
            format!("V2.{f}"),
            format!("{f}: a priori model (no fitted parameter), all batches, median |error|"),
            median_abs(all.map(|r| r.err_a_priori)),
            0.30,
        ));
    }
    if let Some(t) = size_target {
        out.push(le(
            format!("V3.{t}"),
            format!("{t}: size-extrapolated model (fitted on smaller limb counts), all batches, median |error|"),
            median_abs(rows.iter().filter(|r| r.field == t).filter_map(|r| r.err_size_extrapolated)),
            0.30,
        ));
    }
    out
}

/// Validation rows as CSV.
pub fn rows_csv(v: &Validation) -> String {
    let mut s = String::from("field,batch,set,measured_ns,measured_q1,measured_q3,working_set_bytes,memory_level,compute_ns,memory_ns,bound,a_priori_ns,latency_ns,calibrated_ns,size_extrapolated_ns,err_a_priori,err_latency,err_calibrated,err_size_extrapolated\n");
    for r in &v.rows {
        let o = |x: Option<f64>| x.map(|v| format!("{v:.6}")).unwrap_or_default();
        s.push_str(&format!(
            "{},{},{},{:.4},{:.4},{:.4},{},{},{:.4},{:.4},{},{:.4},{:.4},{:.4},{},{:.6},{:.6},{:.6},{}\n",
            r.field,
            r.batch,
            r.set,
            r.measured_ns,
            r.measured_q1,
            r.measured_q3,
            r.working_set_bytes,
            r.memory_level,
            r.compute_ns,
            r.memory_ns,
            r.bound,
            r.a_priori_ns,
            r.latency_ns,
            r.calibrated_ns,
            o(r.size_extrapolated_ns),
            r.err_a_priori,
            r.err_latency,
            r.err_calibrated,
            o(r.err_size_extrapolated),
        ));
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linfit_recovers_a_line() {
        let (a, b) = linfit(&[1.0, 4.0, 9.0], &[3.0, 9.0, 19.0]).unwrap();
        assert!((a - 1.0).abs() < 1e-12 && (b - 2.0).abs() < 1e-12);
        assert!(linfit(&[2.0, 2.0], &[1.0, 3.0]).is_none());
    }
}

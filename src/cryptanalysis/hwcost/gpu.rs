//! GPU kernels through the model: an SM is a core, an issue pipe is a
//! [`super::OpCost::pipe`], an operation costs `1 / (lane-ops per SM-clock)`
//! SM-cycles per lane-instruction, and the pipes of an update overlap
//! ([`super::Phase::pipes_overlap`]): the busiest pipe sets the SM-clocks
//! per update, as in `ecc2k130/roofline.py`.
//!
//! Validated against the ECC2K-130 build sweep
//! (`ecc2k130/benchmarks/roofline/calibration.json`): each row holds a
//! build's dynamic lane-instruction counts per update (ALU, FMA, CLMAD,
//! issue, LSU) and its measured rate on one GPU.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use super::fit::{median_abs, Claim};
use super::{
    predict, Count, Deployment, Hardware, MemLevel, OpCost, Phase, Price, ScalePoint, Workload,
    HARDWARE_SCHEMA, WORKLOAD_SCHEMA,
};

/// One build of the sweep.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct GpuRow {
    pub gpu: String,
    pub arm: String,
    /// Measured updates per second, in billions.
    pub bps: f64,
    #[serde(rename = "smClocksPerUpdate")]
    pub sm_clocks_per_update: f64,
    pub alu: f64,
    pub fma: f64,
    pub clmad: f64,
    pub issue: f64,
    #[serde(default)]
    pub lsu: Option<f64>,
}

#[derive(Deserialize)]
struct CalibrationFile {
    rows: Vec<GpuRow>,
}

pub fn load_rows(path: &std::path::Path) -> Result<Vec<GpuRow>, String> {
    let s = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let f: CalibrationFile =
        serde_json::from_str(&s).map_err(|e| format!("{}: {e}", path.display()))?;
    Ok(f.rows)
}

/// How ALU and FMA work share the integer datapath.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PipeModel {
    /// ALU and FMA are separate pipes (roofline model M4).
    Separate,
    /// ALU and FMA share one datapath measured on a mixed stream (M5).
    SharedInt,
}

impl PipeModel {
    pub fn name(self) -> &'static str {
        match self {
            PipeModel::Separate => "separate",
            PipeModel::SharedInt => "shared-int",
        }
    }
}

struct Ceiling {
    rate: f64,
    source: &'static str,
}

struct GpuSpec {
    id: &'static str,
    name: &'static str,
    sms: u32,
    clock_ghz: f64,
    clock_source: &'static str,
    alu: Ceiling,
    fma: Ceiling,
    int_shared: Ceiling,
    clmad: Ceiling,
    lsu: Ceiling,
    issue: Ceiling,
    usd_per_hour: f64,
}

const ISSUE: Ceiling = Ceiling {
    rate: 128.0,
    source:
        "assumed: 4 schedulers × 32 lanes per SM-clock (architectural; ecc2k130/roofline.py:77-104)",
};

const SPECS: [GpuSpec; 2] = [
    GpuSpec {
        id: "rtx-pro-6000",
        name: "RTX PRO 6000 Blackwell Server Edition (sm_120), 188 SMs",
        sms: 188,
        clock_ghz: 2.415,
        clock_source: "measured: median SM clock sampled after the 20b runs (ecc2k130/benchmarks/two-chains/summary.json:43-45; the clock calibration.json's smClocksPerUpdate uses)",
        alu: Ceiling { rate: 62.1, source: "measured: LOP3 stream, ecc2k130/benchmarks/hardware-limits/result.json (THROUGHPUT-CEILING.md:27-42)" },
        fma: Ceiling { rate: 62.3, source: "measured: IMAD stream, ecc2k130/benchmarks/hardware-limits/result.json (THROUGHPUT-CEILING.md:27-42)" },
        int_shared: Ceiling { rate: 69.3, source: "measured: mixed LOP3+IMAD stream (ecc2k130/THROUGHPUT-CEILING.md:41-42)" },
        clmad: Ceiling { rate: 1.99, source: "measured: CLMAD lo+hi, ecc2k130/benchmarks/fast-clmad/probe/summary.json:29-34" },
        lsu: Ceiling { rate: 15.9, source: "measured: LDS stream, ecc2k130/benchmarks/hardware-limits/result.json (THROUGHPUT-CEILING.md:27-42)" },
        issue: ISSUE,
        usd_per_hour: 3.031,
    },
    GpuSpec {
        id: "b200",
        name: "B200 (sm_100), 148 SMs",
        sms: 148,
        clock_ghz: 1.965,
        clock_source: "measured: SM clock under load (ecc2k130/benchmarks/fast-clmad/probe/summary.json:68-73; ecc2k130/AUTOSWEEP.md:20-21)",
        alu: Ceiling { rate: 63.7, source: "measured: LOP3, ecc2k130/benchmarks/fast-clmad/probe/summary.json:68-73" },
        fma: Ceiling { rate: 64.0, source: "assumed: not measured on B200; architectural 64 (ecc2k130/roofline.py)" },
        int_shared: Ceiling { rate: 69.0, source: "assumed: the shared-int width of ecc2k130/roofline.py's B200 model (calibration.json models.b200)" },
        clmad: Ceiling { rate: 29.1, source: "measured: CLMAD lo+hi, ecc2k130/benchmarks/fast-clmad/probe/summary.json:68-73" },
        lsu: Ceiling { rate: 16.0, source: "assumed: not measured on B200; architectural 16" },
        issue: ISSUE,
        usd_per_hour: 6.25,
    },
];

/// Expected ECC2K-130 work, log₂ iterations (paper/ecc2k130/sec-solvetime.tex:18-21).
pub const ECC2K130_LOG2_WORK: f64 = 60.9;

/// Dollars for the expected ECC2K-130 work at `bps` billion updates per
/// second and `usd_per_hour` (paper/ecc2k130's solvetime arithmetic).
pub fn usd_per_solve(bps: f64, usd_per_hour: f64) -> f64 {
    ECC2K130_LOG2_WORK.exp2() / (bps * 1e9) / 3600.0 * usd_per_hour
}

pub fn usd_per_hour(gpu: &str) -> Result<f64, String> {
    SPECS
        .iter()
        .find(|s| s.id == gpu)
        .map(|s| s.usd_per_hour)
        .ok_or_else(|| format!("no GPU spec `{gpu}`"))
}

pub fn gpu_ids() -> Vec<&'static str> {
    SPECS.iter().map(|s| s.id).collect()
}

/// The profile of a GPU under a pipe model.
pub fn profile(gpu: &str, model: PipeModel) -> Result<Hardware, String> {
    let s = SPECS
        .iter()
        .find(|s| s.id == gpu)
        .ok_or_else(|| format!("no GPU spec `{gpu}` (have {:?})", gpu_ids()))?;
    let mut ops = BTreeMap::new();
    let mut op = |name: &str, c: &Ceiling, pipe: &str| {
        ops.insert(
            name.to_string(),
            OpCost {
                cycles: 1.0 / c.rate,
                latency_cycles: None,
                pipe: Some(pipe.to_string()),
                source: format!("{} lane-ops per SM-clock: {}", c.rate, c.source),
            },
        );
    };
    match model {
        PipeModel::Separate => {
            op("alu", &s.alu, "alu");
            op("fma", &s.fma, "fma");
        }
        PipeModel::SharedInt => {
            op("alu", &s.int_shared, "int");
            op("fma", &s.int_shared, "int");
        }
    }
    op("clmad", &s.clmad, "clmad");
    op("lsu", &s.lsu, "lsu");
    op("issue", &s.issue, "issue");
    Ok(Hardware {
        schema: HARDWARE_SCHEMA.into(),
        name: format!("{} [{}]", s.name, model.name()),
        kind: "gpu".into(),
        clock_ghz: s.clock_ghz,
        clock_source: s.clock_source.into(),
        cores: s.sms,
        isa: vec![],
        ops,
        memory: vec![MemLevel {
            name: "registers".into(),
            capacity_bytes: None,
            shared: false,
            read_gbps_one_core: 1e6,
            read_gbps_all_cores: 1e9,
            source: "assumed: the update is register-resident; no memory term".into(),
        }],
        scaling: vec![
            ScalePoint { cores: 1, efficiency: 1.0 },
            ScalePoint { cores: s.sms, efficiency: 1.0 },
        ],
        scaling_source: "assumed: SMs are independent (efficiency 1; occupancy and latency are what a calibrated η absorbs)".into(),
        channels: BTreeMap::new(),
        price: Price {
            usd_per_hour: s.usd_per_hour,
            egress_usd_per_gb: 0.0,
            storage_usd_per_gb_month: 0.0,
            source: "measured list price: Modal, 2026-09-19 (ecc2k130/benchmarks/modal-gpus/SURVEY.md:111-122)".into(),
        },
        notes: vec![format!("pipe model: {}", model.name())],
    })
}

/// `iterations` updates of a sweep row's kernel.
pub fn workload(row: &GpuRow, iterations: Count) -> Workload {
    let mut ops = BTreeMap::new();
    ops.insert("alu".to_string(), row.alu);
    ops.insert("fma".to_string(), row.fma);
    ops.insert("clmad".to_string(), row.clmad);
    ops.insert("issue".to_string(), row.issue);
    if let Some(l) = row.lsu {
        ops.insert("lsu".to_string(), l);
    }
    Workload {
        schema: WORKLOAD_SCHEMA.into(),
        name: format!("{} {}", row.gpu, row.arm),
        phases: vec![Phase {
            name: "update".into(),
            iterations,
            ops_per_iteration: ops,
            overhead_cycles_per_iteration: 0.0,
            overhead_source: String::new(),
            bytes_per_iteration: 0.0,
            working_set_bytes_per_core: 0.0,
            parallel_fraction: 1.0,
            pipes_overlap: true,
        }],
        transfers: vec![],
        storage: vec![],
        egress_bytes: Count::Plain(0.0),
        notes: vec!["dynamic lane-instruction counts per update from ecc2k130/benchmarks/roofline/calibration.json".into()],
    }
}

/// Predicted billions of updates per second on a whole GPU, and the
/// binding pipe.
pub fn predicted_bps(hw: &Hardware, row: &GpuRow) -> Result<(f64, String), String> {
    let n = 1e12;
    let p = predict(
        hw,
        &workload(row, Count::Plain(n)),
        Deployment {
            machines: 1.0,
            cores_per_machine: hw.cores,
        },
    )?;
    let pipe = p.phases[0].binding_pipe.clone().unwrap_or_default();
    Ok((n / p.wall_s / 1e9, pipe))
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CheckRow {
    pub gpu: String,
    pub arm: String,
    pub model: String,
    pub measured_bps: f64,
    pub a_priori_bps: f64,
    pub binding_pipe: String,
    /// measured / a priori: the share of the binding pipe's ceiling reached.
    pub eta: f64,
    /// `same-gpu` (η from this GPU's fit arm) or `cross-gpu` (η from the
    /// reference GPU's fit arm); `fit` for the fit arm itself.
    pub set: String,
    pub calibrated_bps: f64,
    pub err_a_priori: f64,
    pub err_calibrated: f64,
    /// Expected ECC2K-130 work at the measured and calibrated rates.
    pub usd_measured: f64,
    pub usd_calibrated: f64,
}

/// Every row under both pipe models.  η is fitted on the row whose arm
/// starts with `fit_arm` on each GPU (`same-gpu` rows) and, for every other
/// GPU, also taken from `reference_gpu`'s fit row (`cross-gpu` rows).
pub fn check(rows: &[GpuRow], fit_arm: &str, reference_gpu: &str) -> Result<Vec<CheckRow>, String> {
    let mut out = Vec::new();
    for model in [PipeModel::Separate, PipeModel::SharedInt] {
        let mut eta_fit: BTreeMap<String, f64> = BTreeMap::new();
        for r in rows.iter().filter(|r| r.arm.starts_with(fit_arm)) {
            let (pred, _) = predicted_bps(&profile(&r.gpu, model)?, r)?;
            eta_fit.entry(r.gpu.clone()).or_insert(r.bps / pred);
        }
        let reference = *eta_fit
            .get(reference_gpu)
            .ok_or_else(|| format!("no `{fit_arm}` row for {reference_gpu}"))?;
        for r in rows {
            let hw = profile(&r.gpu, model)?;
            let (pred, pipe) = predicted_bps(&hw, r)?;
            let is_fit = r.arm.starts_with(fit_arm);
            let mut sets = vec![(
                if is_fit { "fit" } else { "same-gpu" },
                *eta_fit
                    .get(&r.gpu)
                    .ok_or_else(|| format!("no `{fit_arm}` row for {}", r.gpu))?,
            )];
            if r.gpu != reference_gpu {
                sets.push(("cross-gpu", reference));
            }
            for (set, eta) in sets {
                let cal = pred * eta;
                out.push(CheckRow {
                    gpu: r.gpu.clone(),
                    arm: r.arm.clone(),
                    model: model.name().into(),
                    measured_bps: r.bps,
                    a_priori_bps: pred,
                    binding_pipe: pipe.clone(),
                    eta: r.bps / pred,
                    set: set.into(),
                    calibrated_bps: cal,
                    err_a_priori: pred / r.bps - 1.0,
                    err_calibrated: cal / r.bps - 1.0,
                    usd_measured: usd_per_solve(r.bps, hw.price.usd_per_hour),
                    usd_calibrated: usd_per_solve(cal, hw.price.usd_per_hour),
                });
            }
        }
    }
    Ok(out)
}

/// The protocol's claims G1–G3 (separate pipes) and the shared-int
/// model's same quantities (reported).
pub fn claims(rows: &[CheckRow]) -> Vec<Claim> {
    let mut out = Vec::new();
    for model in ["separate", "shared-int"] {
        let m = |r: &&CheckRow| r.model == model;
        let own: Vec<&CheckRow> = rows
            .iter()
            .filter(m)
            .filter(|r| r.set != "cross-gpu")
            .collect();
        let same: Vec<&CheckRow> = rows
            .iter()
            .filter(m)
            .filter(|r| r.set == "same-gpu")
            .collect();
        let cross: Vec<&CheckRow> = rows
            .iter()
            .filter(m)
            .filter(|r| r.set == "cross-gpu")
            .collect();
        let g1 = median_abs(own.iter().map(|r| r.err_a_priori));
        let g2 = median_abs(same.iter().map(|r| r.err_calibrated));
        let g3 = median_abs(cross.iter().map(|r| r.err_calibrated));
        let over = own.iter().filter(|r| r.eta > 1.0).count() as f64;
        if model == "separate" {
            out.push(Claim::at_most(
                "G1",
                &format!("a priori, all {} builds, median |error|", own.len()),
                g1,
                0.30,
            ));
            out.push(Claim::report(
                "G1.eta>1",
                "builds measured above the a priori ceiling (η > 1)",
                over,
            ));
            out.push(Claim::at_most(
                "G2",
                &format!(
                    "η from each GPU's base build, {} other builds, median |error|",
                    same.len()
                ),
                g2,
                0.30,
            ));
            out.push(Claim::at_most("G2.10", "the same, against 10%", g2, 0.10));
            out.push(Claim::at_most(
                "G3",
                &format!(
                    "η from the RTX PRO 6000 base build carried to {} B200 builds, median |error|",
                    cross.len()
                ),
                g3,
                0.30,
            ));
        } else {
            out.push(Claim::report(
                "G1.shared",
                "shared-int pipes: a priori, median |error|",
                g1,
            ));
            out.push(Claim::report(
                "G1.shared.eta>1",
                "shared-int pipes: builds with η > 1",
                over,
            ));
            out.push(Claim::report(
                "G2.shared",
                "shared-int pipes: η per GPU, median |error|",
                g2,
            ));
            out.push(Claim::report(
                "G3.shared",
                "shared-int pipes: η carried across GPUs, median |error|",
                g3,
            ));
        }
    }
    out
}

pub fn rows_csv(rows: &[CheckRow]) -> String {
    let mut s = String::from("gpu,arm,model,set,measured_bps,a_priori_bps,binding_pipe,eta,calibrated_bps,err_a_priori,err_calibrated,usd_measured,usd_calibrated\n");
    for r in rows {
        s.push_str(&format!(
            "{},\"{}\",{},{},{:.6},{:.6},{},{:.6},{:.6},{:.6},{:.6},{:.0},{:.0}\n",
            r.gpu,
            r.arm,
            r.model,
            r.set,
            r.measured_bps,
            r.a_priori_bps,
            r.binding_pipe,
            r.eta,
            r.calibrated_bps,
            r.err_a_priori,
            r.err_calibrated,
            r.usd_measured,
            r.usd_calibrated
        ));
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(alu: f64, clmad: f64) -> GpuRow {
        GpuRow {
            gpu: "rtx-pro-6000".into(),
            arm: "base".into(),
            bps: 20.0,
            sm_clocks_per_update: 0.0,
            alu,
            fma: 100.0,
            clmad,
            issue: 1000.0,
            lsu: None,
        }
    }

    #[test]
    fn busiest_pipe_sets_the_rate() {
        let hw = profile("rtx-pro-6000", PipeModel::Separate).unwrap();
        let (bps, pipe) = predicted_bps(&hw, &row(1319.6, 33.125)).unwrap();
        assert_eq!(pipe, "alu");
        let want = 188.0 * 2.415e9 / (1319.6 / 62.1) / 1e9;
        assert!((bps - want).abs() < 1e-9 * want);
        let (_, pipe) = predicted_bps(&hw, &row(100.0, 60.0)).unwrap();
        assert_eq!(pipe, "clmad");
    }

    #[test]
    fn shared_int_adds_alu_and_fma() {
        let hw = profile("rtx-pro-6000", PipeModel::SharedInt).unwrap();
        let (bps, pipe) = predicted_bps(&hw, &row(1319.6, 33.125)).unwrap();
        assert_eq!(pipe, "int");
        let want = 188.0 * 2.415e9 / ((1319.6 + 100.0) / 69.3) / 1e9;
        assert!((bps - want).abs() < 1e-9 * want);
    }
}

//! `hwcost` — a calibrated analytical model of how long a calculation
//! takes on given hardware and what it costs.
//!
//! ```text
//! hwcost measure --out DIR [--quick]           # calibrate this host
//! hwcost fit --run DIR [--fit-max-batch 1024]  # profile + validation
//! hwcost predict --hardware H.json --workload W.json [--machines M] [--cores C]
//! hwcost rho --hardware H.json --group-log2 64,128,256 [--machines M] [--cores C]
//! ```
//!
//! The model and its schemas: `crypto_lib::cryptanalysis::hwcost`.

use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::{Parser, Subcommand};

use crypto_lib::cryptanalysis::hwcost::fit::{self, rows_csv};
use crypto_lib::cryptanalysis::hwcost::measure::{self, Measurements};
use crypto_lib::cryptanalysis::hwcost::rho::{self, RhoParams};
use crypto_lib::cryptanalysis::hwcost::{self as hw, Deployment, Hardware, Workload};
use crypto_lib::cryptanalysis::hwcost::{compare, gpu, ic};

#[derive(Parser)]
#[command(
    name = "hwcost",
    about = "Calibrated analytical time-and-cost model (compute, memory roofline, parallel scaling, transfers, $)"
)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Measure this host: clock, field-op costs, memory bandwidth, the rho
    /// kernel over fields and batch sizes, and its thread scaling.
    Measure {
        #[arg(long)]
        out: PathBuf,
        /// A short run (2 fields, 3 batch sizes, 2 reps) for smoke tests.
        #[arg(long)]
        quick: bool,
        #[arg(long)]
        reps: Option<usize>,
        #[arg(long)]
        target_ms: Option<f64>,
        /// Measure only these fields (comma-separated; default all).
        #[arg(long, value_delimiter = ',')]
        fields: Vec<String>,
        /// Also time `prime_fast::rho_parallel` on a curve of this many bits.
        #[arg(long)]
        prime_fast_rho_bits: Option<u32>,
    },
    /// Fit a hardware profile from a measurement run and validate the model
    /// against the measured kernel.
    Fit {
        #[arg(long)]
        run: PathBuf,
        /// Batch sizes up to this are the fit set; larger ones are tested.
        #[arg(long, default_value_t = 1024)]
        fit_max_batch: usize,
        /// Where to write profile.json, validation.json, validation.csv
        /// (default: the run directory).
        #[arg(long)]
        out: Option<PathBuf>,
        /// Take model C's residuals from this fitted run instead of fitting
        /// them here (calibration transfer across builds or hosts).
        #[arg(long)]
        residual_from: Option<PathBuf>,
    },
    /// Predict time and cost of a workload on a hardware profile.
    Predict {
        #[arg(long)]
        hardware: PathBuf,
        #[arg(long)]
        workload: PathBuf,
        #[arg(long, default_value_t = 1.0)]
        machines: f64,
        /// Cores per machine (default: all of the profile's).
        #[arg(long)]
        cores: Option<u32>,
        #[arg(long)]
        json: bool,
    },
    /// Pollard rho with distinguished points on a prime-order curve:
    /// generate the workload from the kernel's exact op counts and predict.
    Rho {
        #[arg(long)]
        hardware: PathBuf,
        /// Comma-separated log2 group orders.
        #[arg(long, value_delimiter = ',')]
        group_log2: Vec<f64>,
        /// A measured field, or `auto` for the smallest that holds the group.
        #[arg(long, default_value = "auto")]
        field: String,
        #[arg(long, default_value_t = 1024)]
        batch: usize,
        /// Comma-separated machine counts.
        #[arg(long, value_delimiter = ',', default_value = "1")]
        machines: Vec<f64>,
        #[arg(long)]
        cores: Option<u32>,
        #[arg(long)]
        negation_map: bool,
        /// Write each generated workload into this directory.
        #[arg(long)]
        write_workloads: Option<PathBuf>,
        /// Print one CSV row per prediction instead of tables.
        #[arg(long)]
        csv: bool,
    },
    /// Compare two builds measured in alternation (A1, B1, A2, B2):
    /// per-arm ratios, drift and classification, and the fp64bf / fp64
    /// kernel ratio in every run.
    Compare {
        /// Baseline runs A1,A2.
        #[arg(long, value_delimiter = ',')]
        base: Vec<PathBuf>,
        /// Candidate runs B1,B2.
        #[arg(long, value_delimiter = ',')]
        cand: Vec<PathBuf>,
        /// Smallest margin an effect must clear.
        #[arg(long, default_value_t = 0.02)]
        floor: f64,
        /// Claim `ID:PREFIX:MIN_BATCH:T`: max(r1, r2) over the arms named
        /// `PREFIX…` (with batch ≥ MIN_BATCH where the arm has one) is ≤ T.
        #[arg(long = "max-ratio")]
        max_ratio: Vec<String>,
        /// Claim `ID:PREFIX:MIN_BATCH:T`: max |r_k − 1| over those arms is ≤ T.
        #[arg(long = "max-dev")]
        max_dev: Vec<String>,
        #[arg(long)]
        out: PathBuf,
    },
    /// Kernel ratio of two fields measured in the same runs, at every
    /// batch size, and the claim that it stays at or below a threshold.
    Pair {
        #[arg(long, value_delimiter = ',')]
        runs: Vec<PathBuf>,
        #[arg(long)]
        a: String,
        #[arg(long)]
        b: String,
        #[arg(long, default_value_t = 4096)]
        min_batch: usize,
        #[arg(long, default_value_t = 0.85)]
        threshold: f64,
        #[arg(long)]
        id: String,
        #[arg(long)]
        out: PathBuf,
    },
    /// Validate the GPU pipe model against a build sweep's recorded
    /// per-pipe counts and measured rates.
    GpuCheck {
        /// `ecc2k130/benchmarks/roofline/calibration.json`.
        #[arg(long)]
        calibration: PathBuf,
        /// Arm name prefix that fits each GPU's efficiency η.
        #[arg(long, default_value = "base")]
        fit_arm: String,
        /// GPU whose η is also carried to the others.
        #[arg(long, default_value = "rtx-pro-6000")]
        reference_gpu: String,
        /// Writes gpu_check.csv and the profiles used.
        /// Also print the expected ECC2K-130 cost at a measured rate,
        /// `GPU:BPS` (repeatable), as a check of the cost arithmetic.
        #[arg(long)]
        check_rate: Vec<String>,
        #[arg(long)]
        out: PathBuf,
    },
    /// Price an IC pipeline and its matched rho in wall time and dollars
    /// from a frozen ecbench campaign (sessions + §5a band file).
    IcPrice {
        #[arg(long)]
        band: PathBuf,
        #[arg(long)]
        sessions: PathBuf,
        #[arg(long)]
        ic_arm: String,
        #[arg(long)]
        rho_arm: String,
        /// log2 r range "lo,hi" the exponents are fitted over.
        #[arg(long, value_delimiter = ',')]
        tier: Vec<f64>,
        /// Comma-separated log2 r to extrapolate to.
        #[arg(long, value_delimiter = ',')]
        targets: Vec<f64>,
        /// Price of one core-hour (default: one vCPU of the illustrative
        /// c7i.xlarge, $0.1785/h ÷ 4).
        #[arg(long, default_value_t = 0.1785 / 4.0)]
        usd_per_core_hour: f64,
        #[arg(long)]
        out: PathBuf,
    },
}

fn read_json<T: serde::de::DeserializeOwned>(p: &Path) -> Result<T, String> {
    let s = std::fs::read_to_string(p).map_err(|e| format!("{}: {e}", p.display()))?;
    serde_json::from_str(&s).map_err(|e| format!("{}: {e}", p.display()))
}

fn write(p: &Path, s: &str) -> Result<(), String> {
    std::fs::write(p, s).map_err(|e| format!("{}: {e}", p.display()))
}

fn json<T: serde::Serialize>(v: &T) -> String {
    serde_json::to_string_pretty(v).unwrap() + "\n"
}

fn run() -> Result<(), String> {
    match Cli::parse().cmd {
        Cmd::Measure {
            out,
            quick,
            reps,
            target_ms,
            fields,
            prime_fast_rho_bits,
        } => {
            let mut params = measure::default_params(quick);
            params.prime_fast_rho_bits = prime_fast_rho_bits;
            if let Some(r) = reps {
                params.reps = r;
            }
            if let Some(t) = target_ms {
                params.target_ms = t;
            }
            if !fields.is_empty() {
                for f in &fields {
                    measure::field_spec(f)?;
                }
                params.fields = fields;
            }
            std::fs::create_dir_all(&out).map_err(|e| e.to_string())?;
            let m = measure::run(params, &mut |s| eprintln!("{s}"))?;
            write(&out.join("measurements.json"), &json(&m))?;
            let mut csv = String::from("group,arm,rep,value,unit\n");
            for s in &m.samples {
                csv.push_str(&format!(
                    "{},{},{},{},{}\n",
                    s.group, s.arm, s.rep, s.value, s.unit
                ));
            }
            write(&out.join("samples.csv"), &csv)?;
            let mut med = String::from("arm,median\n");
            for (a, v) in measure::medians(&m) {
                med.push_str(&format!("{a},{v}\n"));
            }
            write(&out.join("medians.csv"), &med)?;
            eprintln!("wrote {}", out.display());
        }
        Cmd::Fit {
            run,
            fit_max_batch,
            out,
            residual_from,
        } => {
            let m: Measurements = read_json(&run.join("measurements.json"))?;
            let other: Option<fit::Validation> = match &residual_from {
                Some(d) => Some(read_json(&d.join("validation.json"))?),
                None => None,
            };
            let (v, profile) =
                fit::validate_with(&m, fit_max_batch, other.as_ref().map(|o| &o.residual_ns))?;
            let out = out.unwrap_or(run);
            std::fs::create_dir_all(&out).map_err(|e| e.to_string())?;
            write(&out.join("profile.json"), &json(&profile))?;
            write(&out.join("validation.json"), &json(&v))?;
            write(&out.join("validation.csv"), &rows_csv(&v))?;
            println!(
                "clock {:.3} GHz; residual ns/iter {:?}",
                v.clock_ghz, v.residual_ns
            );
            println!("scaling {:?}", v.scaling);
            for c in &v.claims {
                println!(
                    "{:<8} {:<4} {:.4} (threshold {}) {}",
                    c.id, c.verdict, c.value, c.threshold, c.statement
                );
            }
        }
        Cmd::Predict {
            hardware,
            workload,
            machines,
            cores,
            json: as_json,
        } => {
            let h: Hardware = read_json(&hardware)?;
            let w: Workload = read_json(&workload)?;
            let dep = Deployment {
                machines,
                cores_per_machine: cores.unwrap_or(h.cores),
            };
            let p = hw::predict(&h, &w, dep)?;
            if as_json {
                print!("{}", json(&p));
            } else {
                print!("{}", hw::render(&p));
            }
        }
        Cmd::Rho {
            hardware,
            group_log2,
            field,
            batch,
            machines,
            cores,
            negation_map,
            write_workloads,
            csv,
        } => {
            let h: Hardware = read_json(&hardware)?;
            let cores = cores.unwrap_or(h.cores);
            if csv {
                println!("group_log2,field,batch,machines,cores,negation_map,iterations_log2,dp_bits,dp_log2,ns_per_iter_compute,bound,memory_level,efficiency,wall_s,wall_human,machine_hours,compute_usd,storage_usd,total_usd,assumptions");
            }
            if let Some(d) = &write_workloads {
                std::fs::create_dir_all(d).map_err(|e| e.to_string())?;
            }
            for &g in &group_log2 {
                let f = if field == "auto" {
                    rho::field_for(g)?.to_string()
                } else {
                    field.clone()
                };
                let counts = measure::counts_for(&f, batch, 20261016)?;
                for &mach in &machines {
                    let p = RhoParams {
                        field: f.clone(),
                        group_log2: g,
                        batch,
                        machines: mach,
                        cores_per_machine: cores,
                        negation_map,
                        dp_per_message: 1024.0,
                    };
                    let (w, sh) = rho::workload(&p, &counts)?;
                    if let Some(d) = &write_workloads {
                        write(
                            &d.join(format!("rho-{f}-n{g}-B{batch}-m{mach}.json")),
                            &json(&w),
                        )?;
                    }
                    let pr = hw::predict(
                        &h,
                        &w,
                        Deployment {
                            machines: mach,
                            cores_per_machine: cores,
                        },
                    )?;
                    if csv {
                        let walk = pr.phases.iter().find(|p| p.name == "walk").unwrap();
                        println!(
                            "{g},{f},{batch},{mach},{cores},{negation_map},{:.4},{},{:.4},{:.4},{},{},{:.4},{:.6e},{},{:.6e},{:.6e},{:.6e},{:.6e},\"{}\"",
                            sh.iterations_log2,
                            sh.dp_bits,
                            sh.distinguished_points_log2,
                            walk.ns_per_iteration_compute,
                            walk.bound,
                            walk.memory_level,
                            walk.efficiency,
                            pr.wall_s,
                            hw::human_duration(pr.wall_s),
                            pr.cost.machine_hours,
                            pr.cost.compute_usd,
                            pr.cost.storage_usd,
                            pr.cost.total_usd,
                            pr.assumptions.join("; ").replace('"', "'"),
                        );
                    } else {
                        print!("{}", hw::render(&pr));
                    }
                }
            }
        }
        Cmd::Compare {
            base,
            cand,
            floor,
            max_ratio,
            max_dev,
            out,
        } => {
            if base.len() != 2 || cand.len() != 2 {
                return Err("--base and --cand take two runs each".into());
            }
            let load = |p: &PathBuf| read_json::<Measurements>(&p.join("measurements.json"));
            let (a1, a2, b1, b2) = (
                load(&base[0])?,
                load(&base[1])?,
                load(&cand[0])?,
                load(&cand[1])?,
            );
            let rows = compare::compare(&a1, &a2, &b1, &b2, floor);
            let name = |p: &PathBuf| {
                p.file_name()
                    .map(|s| s.to_string_lossy().into_owned())
                    .unwrap_or_default()
            };
            let runs = vec![
                (name(&base[0]), &a1),
                (name(&cand[0]), &b1),
                (name(&base[1]), &a2),
                (name(&cand[1]), &b2),
            ];
            let pairs = compare::fp64_pairs(&runs);
            std::fs::create_dir_all(&out).map_err(|e| e.to_string())?;
            write(&out.join("compare.csv"), &compare::compare_csv(&rows))?;
            write(
                &out.join("fp64_pairs.csv"),
                &compare::pairs_csv(&pairs, "fp64", "fp64bf"),
            )?;
            let mut counts: std::collections::BTreeMap<(String, String), usize> =
                Default::default();
            for r in &rows {
                let group = r.arm.split('.').next().unwrap_or("").to_string();
                *counts.entry((group, r.class.clone())).or_default() += 1;
            }
            for ((g, c), n) in &counts {
                println!("{g:<8} {c:<20} {n}");
            }
            let base_names = [name(&base[0]), name(&base[1])];
            let f1: Vec<&compare::PairRow> = pairs
                .iter()
                .filter(|p| p.batch >= 4096 && base_names.contains(&p.run))
                .collect();
            let worst = f1.iter().map(|p| p.ratio).fold(f64::MIN, f64::max);
            let mut claims = Vec::new();
            // F1 (round 2) only where the runs measured `fp64bf`.
            if !f1.is_empty() {
                claims.push(fit::Claim::at_most(
                    "F1",
                    &format!(
                        "worst fp64bf/fp64 kernel ratio at B ≥ 4096 in the baseline runs ({} cells)",
                        f1.len()
                    ),
                    worst,
                    0.85,
                ));
            }
            for (spec, dev) in max_ratio
                .iter()
                .map(|s| (s, false))
                .chain(max_dev.iter().map(|s| (s, true)))
            {
                claims.push(compare::ratio_claim(&rows, spec, dev)?);
            }
            for c in &claims {
                println!(
                    "{} {}: {:.4} (threshold {}) {}",
                    c.id, c.verdict, c.value, c.threshold, c.statement
                );
            }
            let classes: Vec<serde_json::Value> = counts
                .iter()
                .map(|((g, c), n)| serde_json::json!({"group": g, "class": c, "arms": n}))
                .collect();
            write(
                &out.join("summary.json"),
                &json(&serde_json::json!({"classes": classes, "claims": claims})),
            )?;
        }
        Cmd::Pair {
            runs,
            a,
            b,
            min_batch,
            threshold,
            id,
            out,
        } => {
            let ms: Vec<(String, Measurements)> = runs
                .iter()
                .map(|p| {
                    let n = p
                        .file_name()
                        .map(|s| s.to_string_lossy().into_owned())
                        .unwrap_or_default();
                    read_json::<Measurements>(&p.join("measurements.json")).map(|m| (n, m))
                })
                .collect::<Result<_, _>>()?;
            let refs: Vec<(String, &Measurements)> =
                ms.iter().map(|(n, m)| (n.clone(), m)).collect();
            let rows = compare::field_pairs(&refs, &a, &b);
            let cells: Vec<&compare::PairRow> =
                rows.iter().filter(|r| r.batch >= min_batch).collect();
            let worst = cells.iter().map(|r| r.ratio).fold(f64::MIN, f64::max);
            let c = fit::Claim::at_most(
                &id,
                &format!(
                    "worst {b}/{a} kernel ratio at B ≥ {min_batch} over {} runs ({} cells)",
                    ms.len(),
                    cells.len()
                ),
                if cells.is_empty() { f64::NAN } else { worst },
                threshold,
            );
            std::fs::create_dir_all(&out).map_err(|e| e.to_string())?;
            write(&out.join("pairs.csv"), &compare::pairs_csv(&rows, &a, &b))?;
            write(&out.join("claims.json"), &json(&vec![c.clone()]))?;
            println!(
                "{} {}: {:.4} (threshold {})",
                c.id, c.verdict, c.value, threshold
            );
        }
        Cmd::GpuCheck {
            calibration,
            fit_arm,
            reference_gpu,
            check_rate,
            out,
        } => {
            let rows = gpu::load_rows(&calibration)?;
            let check = gpu::check(&rows, &fit_arm, &reference_gpu)?;
            std::fs::create_dir_all(&out).map_err(|e| e.to_string())?;
            write(&out.join("gpu_check.csv"), &gpu::rows_csv(&check))?;
            write(&out.join("gpu_check.json"), &json(&check))?;
            let claims = gpu::claims(&check);
            write(&out.join("claims.json"), &json(&claims))?;
            for c in &claims {
                println!(
                    "{:<18} {:<6} {:.4} {}",
                    c.id, c.verdict, c.value, c.statement
                );
            }
            for id in gpu::gpu_ids() {
                for m in [gpu::PipeModel::Separate, gpu::PipeModel::SharedInt] {
                    write(
                        &out.join(format!("profile-{id}-{}.json", m.name())),
                        &json(&gpu::profile(id, m)?),
                    )?;
                }
            }
            println!("{} rows checked against {} builds", check.len(), rows.len());
            for cr in &check_rate {
                let (g, b) = cr.split_once(':').ok_or("--check-rate takes GPU:BPS")?;
                let bps: f64 = b.parse().map_err(|_| format!("bad rate `{b}`"))?;
                let usd = gpu::usd_per_hour(g)?;
                println!(
                    "expected ECC2K-130 work 2^{} at {bps} B/s on {g} (${usd}/h): {:.1} h, ${:.0}",
                    gpu::ECC2K130_LOG2_WORK,
                    gpu::ECC2K130_LOG2_WORK.exp2() / (bps * 1e9) / 3600.0,
                    gpu::usd_per_solve(bps, usd)
                );
            }
        }
        Cmd::IcPrice {
            band,
            sessions,
            ic_arm,
            rho_arm,
            tier,
            targets,
            usd_per_core_hour,
            out,
        } => {
            if tier.len() != 2 {
                return Err("--tier takes lo,hi".into());
            }
            let pairs = ic::pairs(&band, &sessions, &ic_arm, &rho_arm, usd_per_core_hour)?;
            let fits = ic::fits(&pairs, (tier[0], tier[1]))?;
            let proj = ic::project(&fits, &targets, usd_per_core_hour);
            std::fs::create_dir_all(&out).map_err(|e| e.to_string())?;
            write(&out.join("pairs.csv"), &ic::pairs_csv(&pairs))?;
            write(&out.join("pairs.json"), &json(&pairs))?;
            write(&out.join("fits.json"), &json(&fits))?;
            write(&out.join("projections.json"), &json(&proj))?;
            let largest = pairs.last().ok_or("no paired curve")?;
            let gp = ic::project_gae(&fits, largest.ns_per_gae, &targets, usd_per_core_hour);
            // Post hoc (protocol addendum 1): P(w) fits × ns per GAE of the
            // matched rho at its largest measured size, named in the last two
            // columns.
            let mut gcsv = String::from("log2_r,w,ns_per_gae,ic_core_years,rho_core_years,ic_usd,rho_usd,ratio,ns_per_gae_arm,ns_per_gae_log2_r\n");
            for g in &gp {
                gcsv.push_str(&format!(
                    "{},{},{:.3},{:.6e},{:.6e},{:.6e},{:.6e},{:.4},{},{:.2}\n",
                    g.log2_r,
                    g.w,
                    g.ns_per_gae,
                    g.ic_core_years,
                    g.rho_core_years,
                    g.ic_usd,
                    g.rho_usd,
                    g.ratio,
                    largest.rho.arm,
                    largest.log2_r
                ));
            }
            write(&out.join("projections_gae.csv"), &gcsv)?;
            let summary = ic::summary(&pairs);
            write(&out.join("summary.json"), &json(&summary))?;
            for c in &summary {
                println!("{:<18} {:.4} {}", c.id, c.value, c.statement);
            }
            let mut csv = String::from("log2_r,ic_core_years,rho_core_years,ic_usd,rho_usd,wall_ratio,band_ratio_w0,band_ratio_w01,band_ratio_w1,rank\n");
            for p in &proj {
                csv.push_str(&format!(
                    "{},{:.6e},{:.6e},{:.6e},{:.6e},{:.4},{:.4},{:.4},{:.4},{}\n",
                    p.log2_r,
                    p.ic_wall_core_years,
                    p.rho_wall_core_years,
                    p.ic_usd,
                    p.rho_usd,
                    p.wall_ratio,
                    p.band_ratio[0],
                    p.band_ratio[1],
                    p.band_ratio[2],
                    p.rank
                ));
            }
            write(&out.join("projections.csv"), &csv)?;
            println!(
                "{} curves paired ({ic_arm} vs {rho_arm}); fits over log2 r in [{}, {}]",
                pairs.len(),
                tier[0],
                tier[1]
            );
        }
    }
    Ok(())
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("hwcost: {e}");
            ExitCode::FAILURE
        }
    }
}

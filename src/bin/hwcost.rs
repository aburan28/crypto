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

use crypto_lib::cryptanalysis::hwcost::fit::{rows_csv, validate};
use crypto_lib::cryptanalysis::hwcost::measure::{self, Measurements};
use crypto_lib::cryptanalysis::hwcost::rho::{self, RhoParams};
use crypto_lib::cryptanalysis::hwcost::{self as hw, Deployment, Hardware, Workload};

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
        } => {
            let mut params = measure::default_params(quick);
            if let Some(r) = reps {
                params.reps = r;
            }
            if let Some(t) = target_ms {
                params.target_ms = t;
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
        } => {
            let m: Measurements = read_json(&run.join("measurements.json"))?;
            let (v, profile) = validate(&m, fit_max_batch)?;
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

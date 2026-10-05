//! Walk the `F_p`-isogeny class of a prime-field curve and record every
//! curve it reaches in the repository's curve formats.
//!
//! ```text
//! isogeny_walk class  --curve p256 --max-ell 61
//! isogeny_walk walk   --curve p256 --max-ell 13 --max-curves 1000 --out DIR
//! isogeny_walk verify --curve p256 --dir DIR
//! isogeny_walk traits --curve p256 --dir DIR --shard 0 --of 4 --out SHARD0
//! isogeny_walk collect --out TRAITS SHARD0 SHARD1 SHARD2 SHARD3
//! isogeny_walk plan   --curve p256 --max-ell 61 --max-curves 20000 \
//!                     --commit <sha> --shards 16 --out SPECS
//! ```
//!
//! `walk` also runs the class audits (`ecc_safety`, the structural report,
//! the PKM signals) and every curve detector
//! (`src/cryptanalysis/isogeny_walk/traits.rs`).  `traits` runs the
//! detectors on one shard of an existing walk; `plan` writes one `taskq`
//! task spec per shard (`taskq/README.md`, submit each with
//! `taskq submit --spec FILE`), and `collect` merges and checks the
//! shards.
//!
//! `walk` writes `curves.yaml` (the `docs/curves/ic/curves.yaml` format),
//! `isogeny_routes.json` (nodes, kernel-certified edges and `IW1` routes)
//! and `walk.json` (configuration, class invariants, counts).  `verify`
//! replays every node identity, order audit and edge certificate from
//! `isogeny_routes.json` alone.  See `src/cryptanalysis/isogeny_walk/`.

use std::path::PathBuf;
use std::process::ExitCode;

use clap::{Args, Parser, Subcommand};
use crypto_lib::cryptanalysis::isogeny_walk::queue::{self, PlanWalk};
use crypto_lib::cryptanalysis::isogeny_walk::record::{sha256_hex, V};
use crypto_lib::cryptanalysis::isogeny_walk::walk::{
    self, ClassInfo, EllKind, RouteIndex, StartCurve, Walk, WalkConfig,
};
use crypto_lib::ecc::curve::CurveParams;
use num_bigint::BigUint;

#[derive(Parser)]
#[command(
    about = "Isogeny walks over a prime-field curve's isogeny class, recorded as ICV1/EC1 curves and IW1 routes"
)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Print the class invariants and the splitting of each ℓ.
    Class {
        #[command(flatten)]
        curve: CurveArgs,
        #[command(flatten)]
        primes: PrimeArgs,
    },
    /// Walk the ℓ-isogeny graphs breadth first and write the records.
    Walk {
        #[command(flatten)]
        curve: CurveArgs,
        #[command(flatten)]
        primes: PrimeArgs,
        /// Stop expanding once this many curves are known.
        #[arg(long, default_value_t = 256)]
        max_curves: usize,
        /// Do not expand curves this many steps from the root.
        #[arg(long)]
        max_depth: Option<usize>,
        /// Seed for root splitting (results do not depend on it).
        #[arg(long, default_value_t = 1)]
        seed: u64,
        /// Points per order audit.
        #[arg(long, default_value_t = 1)]
        audit_points: usize,
        /// Worker threads (default: all cores).
        #[arg(long)]
        threads: Option<usize>,
        /// Also walk Atkin primes (they have no rational ℓ-isogeny).
        #[arg(long)]
        include_atkin: bool,
        /// Walked curves on which the class audits are re-run to check that
        /// their verdicts are the root's.
        #[arg(long, default_value_t = 8)]
        class_audit_sample: usize,
        /// Skip the class audits.
        #[arg(long)]
        no_class_audits: bool,
        /// Output directory.
        #[arg(long)]
        out: PathBuf,
    },
    /// Run the curve detectors on one shard of a walk directory.
    Traits {
        #[command(flatten)]
        curve: CurveArgs,
        /// A walk directory (its isogeny_routes.json is read).
        #[arg(long)]
        dir: PathBuf,
        #[arg(long, default_value_t = 0)]
        shard: usize,
        #[arg(long, default_value_t = 1)]
        of: usize,
        /// Walked curves shard 0 re-audits with the class audits.
        #[arg(long, default_value_t = 8)]
        class_audit_sample: usize,
        /// Output directory; default $TASKQ_OUTPUT_DIR.
        #[arg(long)]
        out: Option<PathBuf>,
    },
    /// Merge trait shard directories, refusing gaps, duplicates and mixed walks.
    Collect {
        #[arg(long)]
        out: PathBuf,
        dirs: Vec<PathBuf>,
    },
    /// Write one taskq task spec per trait shard of a walk.
    Plan {
        #[command(flatten)]
        curve: CurveArgs,
        #[command(flatten)]
        primes: PrimeArgs,
        #[arg(long, default_value_t = 256)]
        max_curves: usize,
        #[arg(long)]
        max_depth: Option<usize>,
        #[arg(long)]
        include_atkin: bool,
        /// Full commit sha the workers check out (must be pushed).
        #[arg(long)]
        commit: String,
        #[arg(long, default_value = "cpu")]
        queue: String,
        #[arg(long, default_value_t = 4)]
        shards: usize,
        /// Per-shard run limit in seconds.
        #[arg(long, default_value_t = 3600)]
        timeout: u64,
        /// Directory for the spec files.
        #[arg(long)]
        out: PathBuf,
    },
    /// Re-verify a walk directory from its isogeny_routes.json.
    Verify {
        #[command(flatten)]
        curve: CurveArgs,
        #[arg(long)]
        dir: PathBuf,
        #[arg(long, default_value_t = 2)]
        audit_points: usize,
    },
}

#[derive(Args)]
struct CurveArgs {
    /// p256, p224, or custom (with --p --a --b --order --gx --gy).
    #[arg(long, default_value = "p256")]
    curve: String,
    #[arg(long)]
    p: Option<String>,
    #[arg(long)]
    a: Option<String>,
    #[arg(long)]
    b: Option<String>,
    /// The subgroup order r of the generator; the group order is r·h.
    #[arg(long)]
    order: Option<String>,
    #[arg(long, default_value_t = 1)]
    cofactor: u32,
    #[arg(long)]
    gx: Option<String>,
    #[arg(long)]
    gy: Option<String>,
    /// Label for a custom curve.
    #[arg(long, default_value = "custom")]
    name: String,
}

#[derive(Args)]
struct PrimeArgs {
    /// Comma-separated odd primes ℓ.
    #[arg(long, value_delimiter = ',')]
    primes: Option<Vec<u64>>,
    /// Every odd prime ℓ up to this bound (ignored with --primes).
    #[arg(long, default_value_t = 13)]
    max_ell: u64,
}

fn int(s: &Option<String>, what: &str) -> Result<BigUint, String> {
    let s = s
        .as_deref()
        .ok_or(format!("--{what} is required for a custom curve"))?;
    let (digits, radix) = match s.strip_prefix("0x") {
        Some(h) => (h, 16),
        None => (s, 10),
    };
    BigUint::parse_bytes(digits.as_bytes(), radix).ok_or(format!("--{what}: not an integer"))
}

impl CurveArgs {
    fn start(&self) -> Result<StartCurve, String> {
        match self.curve.as_str() {
            "p256" | "P-256" => Ok(StartCurve::p256()),
            "p224" | "P-224" => Ok(StartCurve::p224()),
            "custom" => {
                let name: &'static str = Box::leak(self.name.clone().into_boxed_str());
                let params = CurveParams {
                    name,
                    p: int(&self.p, "p")?,
                    a: int(&self.a, "a")?,
                    b: int(&self.b, "b")?,
                    gx: int(&self.gx, "gx")?,
                    gy: int(&self.gy, "gy")?,
                    n: int(&self.order, "order")?,
                    h: self.cofactor,
                };
                Ok(StartCurve::from_params(&params, false))
            }
            other => Err(format!("unknown curve {other}; use p256, p224 or custom")),
        }
    }
}

impl CurveArgs {
    /// The arguments that select this curve, for a planned command.
    fn argv(&self) -> Vec<String> {
        let mut v = vec!["--curve".to_string(), self.curve.clone()];
        if self.curve == "custom" {
            for (k, x) in [
                ("--p", &self.p),
                ("--a", &self.a),
                ("--b", &self.b),
                ("--order", &self.order),
                ("--gx", &self.gx),
                ("--gy", &self.gy),
            ] {
                if let Some(x) = x {
                    v.extend([k.to_string(), x.clone()]);
                }
            }
            v.extend([
                "--cofactor".into(),
                self.cofactor.to_string(),
                "--name".into(),
                self.name.clone(),
            ]);
        }
        v
    }
}

impl PrimeArgs {
    fn list(&self) -> Vec<u64> {
        match &self.primes {
            Some(p) => p.clone(),
            None => (3..=self.max_ell)
                .step_by(2)
                .filter(|&l| (2..l).take_while(|d| d * d <= l).all(|d| l % d != 0))
                .collect(),
        }
    }
}

fn write(path: &PathBuf, text: &str) -> Result<String, String> {
    std::fs::write(path, text).map_err(|e| format!("{}: {e}", path.display()))?;
    Ok(sha256_hex(text))
}

fn run(cli: Cli) -> Result<(), String> {
    match cli.cmd {
        Cmd::Class { curve, primes } => {
            let start = curve.start()?;
            let class = ClassInfo::compute(&start, &primes.list());
            print!("{}", class.record().json());
            Ok(())
        }
        Cmd::Walk {
            curve,
            primes,
            max_curves,
            max_depth,
            seed,
            audit_points,
            threads,
            include_atkin,
            class_audit_sample,
            no_class_audits,
            out,
        } => {
            if let Some(t) = threads {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(t)
                    .build_global()
                    .map_err(|e| e.to_string())?;
            }
            let start = curve.start()?;
            let requested = primes.list();
            let class = ClassInfo::compute(&start, &requested);
            let mut walked = Vec::new();
            let mut skipped = Vec::new();
            for &(l, kind, _) in &class.ells {
                if kind == EllKind::Atkin && !include_atkin {
                    skipped.push((
                        l,
                        "atkin: (t^2-4p / ell) = -1, no rational ell-isogeny".to_string(),
                    ));
                } else {
                    walked.push(l);
                }
            }
            if walked.is_empty() {
                return Err("no walkable prime among those requested".into());
            }
            let config = WalkConfig {
                primes: walked,
                max_curves,
                max_depth: max_depth.unwrap_or(usize::MAX),
                seed,
                audit_points,
                skipped,
            };
            let mut w = Walk::new(start, config)?;
            eprintln!(
                "isogeny_walk: {} from {}; ℓ ∈ {:?}; Φ_ℓ built in {:?} ms",
                w.start.name,
                w.start.tag,
                w.config.primes,
                w.stats.phi_ms.iter().map(|x| x.1).collect::<Vec<_>>()
            );
            w.run();
            if !no_class_audits {
                w.run_class_audits(class_audit_sample);
            }
            let ids: Vec<_> = (0..w.nodes.len()).map(|i| w.node_ids(i)).collect();
            let routes = RouteIndex::build(&w, &ids);
            std::fs::create_dir_all(&out).map_err(|e| e.to_string())?;
            let yaml_sha = write(&out.join("curves.yaml"), &w.curves_yaml(&ids, &routes))?;
            let routes_sha = write(
                &out.join("isogeny_routes.json"),
                &w.routes_json(&ids, &routes).json(),
            )?;
            let mut summary = match w.summary(&ids) {
                V::Map(m) => m,
                _ => unreachable!(),
            };
            summary.push((
                "outputs".into(),
                V::map(vec![
                    ("curves.yaml", V::s(yaml_sha)),
                    ("isogeny_routes.json", V::s(routes_sha)),
                ]),
            ));
            let summary = V::Map(summary);
            write(&out.join("walk.json"), &summary.json())?;
            eprintln!(
                "isogeny_walk: {} curves ({} expanded), {} verified edges, {} failures in {} ms; wrote {}",
                w.nodes.len(),
                w.nodes.iter().filter(|n| n.expanded).count(),
                w.edges.len(),
                w.failures.len(),
                w.stats.walk_ms,
                out.display()
            );
            if w.failures.is_empty() {
                Ok(())
            } else {
                Err(format!(
                    "{} roots did not yield a certified edge (see walk.json)",
                    w.failures.len()
                ))
            }
        }
        Cmd::Traits {
            curve,
            dir,
            shard,
            of,
            class_audit_sample,
            out,
        } => {
            let start = curve.start()?;
            let out = out
                .or_else(|| std::env::var_os("TASKQ_OUTPUT_DIR").map(PathBuf::from))
                .ok_or("--out or TASKQ_OUTPUT_DIR is required")?;
            let path = dir.join("isogeny_routes.json");
            let text =
                std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
            let shard_out = queue::run_shard(&text, &start, shard, of, class_audit_sample)?;
            queue::write_shard(&out, &shard_out)?;
            eprintln!(
                "isogeny_walk: trait shard {shard}/{of} of {} written to {}",
                path.display(),
                out.display()
            );
            Ok(())
        }
        Cmd::Collect { out, dirs } => {
            let (merged, summary, class) = queue::collect(&dirs)?;
            std::fs::create_dir_all(&out).map_err(|e| e.to_string())?;
            write(&out.join(queue::TRAITS_FILE), &merged)?;
            write(&out.join("collect.json"), &summary.json())?;
            if let Some(c) = class {
                write(&out.join(queue::CLASS_FILE), &c)?;
            }
            print!("{}", summary.json());
            Ok(())
        }
        Cmd::Plan {
            curve,
            primes,
            max_curves,
            max_depth,
            include_atkin,
            commit,
            queue: queue_name,
            shards,
            timeout,
            out,
        } => {
            curve.start()?;
            let mut walk_args = vec![
                "--primes".to_string(),
                primes
                    .list()
                    .iter()
                    .map(u64::to_string)
                    .collect::<Vec<_>>()
                    .join(","),
                "--max-curves".into(),
                max_curves.to_string(),
            ];
            if let Some(d) = max_depth {
                walk_args.extend(["--max-depth".into(), d.to_string()]);
            }
            if include_atkin {
                walk_args.push("--include-atkin".into());
            }
            let plan_walk = PlanWalk {
                curve_args: curve.argv(),
                walk_args,
            };
            let specs = queue::plan(&plan_walk, &commit, &queue_name, shards, timeout)?;
            std::fs::create_dir_all(&out).map_err(|e| e.to_string())?;
            for (i, spec) in specs.iter().enumerate() {
                write(
                    &out.join(format!("shard-{i:04}-of-{shards:04}.json")),
                    &spec.json(),
                )?;
            }
            eprintln!(
                "isogeny_walk: {shards} taskq specs in {}; submit each with `taskq submit --spec FILE`",
                out.display()
            );
            Ok(())
        }
        Cmd::Verify {
            curve,
            dir,
            audit_points,
        } => {
            let start = curve.start()?;
            let path = dir.join("isogeny_routes.json");
            let text =
                std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
            let routes: serde_json::Value =
                serde_json::from_str(&text).map_err(|e| e.to_string())?;
            let (nodes, edges) = walk::verify_routes(&routes, &start, audit_points)?;
            let receipt = V::map(vec![
                ("schema", V::s("isogeny-walk-verify/v1")),
                ("isogeny_routes_sha256", V::s(sha256_hex(&text))),
                ("curves_verified", V::int(nodes)),
                ("edges_verified", V::int(edges)),
                ("audit_points", V::int(audit_points)),
                ("result", V::s("pass")),
            ]);
            print!("{}", receipt.json());
            Ok(())
        }
    }
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("isogeny_walk: {e}");
            ExitCode::FAILURE
        }
    }
}

//! Live progress for the fixed public Koblitz index-calculus example.
use clap::{Parser, ValueEnum};
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    koblitz_index_calculus_dlp_with_progress, DecompositionStrategy, KoblitzCurve, KoblitzIcEvent,
    KoblitzIcOptions,
};
use num_bigint::BigUint;
use std::{
    io::{self, Write},
    process::ExitCode,
    time::Instant,
};

#[derive(Parser)]
#[command(
    name = "ic",
    version,
    about = "Index calculus with live stage progress"
)]
#[command(
    long_about = "Run the full index-calculus pipeline on the fixed public K_0 / GF(2^9) example with known logarithm 53. Displays factor-base construction, relation collection, linear algebra, and verification. This demo does not accept Certicom challenges or external targets."
)]
struct Cli {
    /// Decomposition engine for the same fixed example.
    #[arg(long, value_enum, default_value_t = Solver::Groebner)]
    solver: Solver,
}
#[derive(Clone, Copy, Debug, ValueEnum)]
enum Solver {
    Groebner,
    Sat,
    Enumerate,
}
impl Solver {
    fn strategy(self) -> DecompositionStrategy {
        match self {
            Self::Groebner => DecompositionStrategy::Groebner,
            Self::Sat => DecompositionStrategy::Sat,
            Self::Enumerate => DecompositionStrategy::Enumerate,
        }
    }
}
struct Progress {
    since: Instant,
    update: Instant,
    collected: usize,
}
impl Progress {
    fn new() -> Self {
        Self {
            since: Instant::now(),
            update: Instant::now(),
            collected: 0,
        }
    }
    fn start(&mut self, number: u8, label: &str) {
        self.since = Instant::now();
        println!("[{number}/4] {label} ...");
    }
    fn event(&mut self, event: KoblitzIcEvent) {
        match event {
            KoblitzIcEvent::FactorBaseStarted => self.start(1, "Factor base"),
            KoblitzIcEvent::FactorBaseReady { points, orbits } => {
                println!("      {points} points; {orbits} Frobenius orbit columns");
                println!("      done in {:.3}s", self.since.elapsed().as_secs_f64());
            }
            KoblitzIcEvent::RelationCollectionStarted { wanted } => {
                self.start(2, "Relation collection");
                println!("      target: {wanted} relations (including surplus)");
            }
            KoblitzIcEvent::RelationProgress {
                collected,
                wanted,
                trials,
            } => {
                if collected != self.collected || self.update.elapsed().as_secs() >= 1 {
                    println!(
                        "      {collected}/{wanted} relations; {trials} trials; {:.3}s",
                        self.since.elapsed().as_secs_f64()
                    );
                    self.collected = collected;
                    self.update = Instant::now();
                }
            }
            KoblitzIcEvent::RelationCollectionFinished { collected, trials } => {
                println!("      collection stopped: {collected} relations from {trials} trials in {:.3}s", self.since.elapsed().as_secs_f64());
            }
            KoblitzIcEvent::LinearAlgebraStarted { rows, columns } => {
                self.start(3, "Linear algebra");
                println!(
                    "      {rows} x {columns} matrix; Gaussian elimination modulo subgroup order"
                );
            }
            KoblitzIcEvent::LinearAlgebraFinished => println!(
                "      system solved in {:.3}s",
                self.since.elapsed().as_secs_f64()
            ),
            KoblitzIcEvent::LinearAlgebraIncomplete => {
                println!("      current matrix did not yield a candidate");
            }
            KoblitzIcEvent::LinearAlgebraSkipped => {
                println!("[3/4] Linear algebra skipped: a direct relation supplied the candidate")
            }
            KoblitzIcEvent::VerificationStarted => self.start(4, "Verification"),
            KoblitzIcEvent::VerificationFinished { verified } => println!(
                "      [d]G = Q: {}; {:.3}s",
                if verified { "PASS" } else { "FAIL" },
                self.since.elapsed().as_secs_f64()
            ),
        }
        // Keep output live when redirected to a log too.
        let _ = io::stdout().flush();
    }
}
fn run(cli: Cli) -> Result<(), String> {
    let started = Instant::now();
    println!("ic — end-to-end index calculus");
    println!("Fixture: public synthetic K_0 / GF(2^9); known logarithm 53");
    println!("Engine: {:?}\n", cli.solver);
    let curve = KoblitzCurve::new(0, 9).ok_or("could not construct the demo curve")?;
    let expected = BigUint::from(53u32);
    let target = curve.mul(curve.generator(), &expected);
    let options = KoblitzIcOptions {
        strategy: cli.solver.strategy(),
        collapse_negation: false,
        stop_on_verified_rank: false,
        allow_direct_relation: false,
        ..KoblitzIcOptions::default()
    };
    let mut progress = Progress::new();
    let report =
        koblitz_index_calculus_dlp_with_progress(&curve, &target, &options, &mut |event| {
            progress.event(event)
        })
        .ok_or("pipeline could not complete; see the last active stage")?;
    if report.sat_invalid_models != 0 {
        return Err(format!("invalid SAT models: {}", report.sat_invalid_models));
    }
    let recovered = report.log.ok_or_else(|| {
        format!(
            "no verified result: {} relations, {} trials, {} inconclusive SAT attempts",
            report.relations, report.trials, report.sat_unknowns,
        )
    })?;
    if recovered != expected {
        return Err(format!(
            "fixture mismatch: expected {expected}, recovered {recovered}"
        ));
    }
    println!("\nComplete: recovered {recovered}; expected {expected}; fixture check PASS");
    println!(
        "Work: {} F4 reductions; {} SAT calls; {} inconclusive SAT attempts",
        report.reductions, report.sat_calls, report.sat_unknowns
    );
    println!("Total elapsed: {:.3}s", started.elapsed().as_secs_f64());
    Ok(())
}
fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("ic: ERROR: {error}");
            ExitCode::FAILURE
        }
    }
}

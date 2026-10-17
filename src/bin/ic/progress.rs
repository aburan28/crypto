//! Optional live text only; never retained in JSON reports or replay evidence.
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    KoblitzIcEvent, KoblitzRelationAttemptDisposition, PdpOutcome,
};
use std::{fmt, io::Write, sync::OnceLock, time::Instant};

static START: OnceLock<Instant> = OnceLock::new();

pub fn init(verbose: bool) {
    if verbose {
        let _ = START.set(Instant::now());
    }
}

pub fn enabled() -> bool {
    START.get().is_some()
}

/// One lock and flush per line keeps parallel workers' messages intact and
/// makes redirected output visible immediately, including through `tee`.
pub fn line(message: fmt::Arguments<'_>) {
    if let Some(start) = START.get() {
        let mut out = std::io::stdout().lock();
        let _ = writeln!(
            out,
            "[{:9.3}s] ic: {message}",
            start.elapsed().as_secs_f64()
        );
        let _ = out.flush();
    }
}

pub fn outcome(outcome: PdpOutcome) -> &'static str {
    match outcome {
        PdpOutcome::Witness => "decomposition found (pending relation verification)",
        PdpOutcome::ProvedUnsat => "proved unsatisfiable",
        PdpOutcome::Incomplete => "budget exhausted; unresolved",
        PdpOutcome::Unsupported => "unsupported",
        PdpOutcome::InvalidModel => "invalid model",
        PdpOutcome::Unresolved => "no decomposition found; unresolved",
        PdpOutcome::Identity => "identity; no decomposition attempted",
    }
}

pub fn event(event: &KoblitzIcEvent, max_trials: usize) {
    if !enabled() {
        return;
    }
    use KoblitzIcEvent::*;
    match event {
        FactorBaseStarted => line(format_args!("[1/4] factor base: building")),
        FactorBaseReady { points, orbits } => line(format_args!(
            "[1/4] factor base: ready; {points} points, {orbits} columns"
        )),
        PairTableStarted => line(format_args!("[1/4] pair table: building")),
        PairTableReady { entries } => {
            line(format_args!("[1/4] pair table: ready; {entries} entries"))
        }
        RelationCollectionStarted { wanted } => line(format_args!(
            "[2/4] collection: started; {wanted} relations requested, trial budget {max_trials}"
        )),
        RelationBatchStarted {
            first_trial,
            trials,
        } => line(format_args!(
            "[2/4] batch: starting trials {first_trial}..={}",
            first_trial + trials - 1
        )),
        RelationAttemptFinished {
            trial,
            disposition,
            conflicts,
            collected,
        } => {
            let state = match disposition {
                KoblitzRelationAttemptDisposition::RelationFound => "relation found",
                KoblitzRelationAttemptDisposition::Refuted => "refuted",
                KoblitzRelationAttemptDisposition::Unsupported => "unsupported",
                KoblitzRelationAttemptDisposition::Unknown => "unresolved",
                KoblitzRelationAttemptDisposition::InvalidModel => "invalid model",
                KoblitzRelationAttemptDisposition::DirectSkipped => "direct relation skipped",
                KoblitzRelationAttemptDisposition::DirectSolved => {
                    "direct relation solved and verified"
                }
            };
            line(format_args!("[2/4] trial {trial}/{max_trials}: {state}; {collected} relations; {conflicts} conflicts"));
        }
        RelationProgress {
            collected,
            wanted,
            trials,
        } => line(format_args!(
            "[2/4] collection: {collected}/{wanted} relations after {trials} trials"
        )),
        RelationCollectionFinished { collected, trials } => line(format_args!(
            "[2/4] collection: stopped; {collected} relations, {trials} trials"
        )),
        LinearAlgebraStarted { rows, columns } => line(format_args!(
            "[3/4] linear algebra: solving {rows} x {columns} relation matrix"
        )),
        MatrixRank {
            rows,
            columns,
            rank,
            candidate_produced,
        } => line(format_args!(
            "[3/4] matrix: {rows} rows, {columns} columns, rank {rank}; candidate {}",
            if *candidate_produced {
                "produced"
            } else {
                "not determined"
            }
        )),
        LinearAlgebraFinished => line(format_args!("[3/4] linear algebra: complete")),
        LinearAlgebraIncomplete => line(format_args!("[3/4] linear algebra: incomplete")),
        LinearAlgebraSkipped => line(format_args!(
            "[3/4] linear algebra: skipped; direct relation"
        )),
        VerificationStarted => line(format_args!("[4/4] verification: checking [d]G = Q")),
        VerificationFinished { verified } => line(format_args!(
            "[4/4] verification: {}",
            if *verified { "passed" } else { "failed" }
        )),
    }
}

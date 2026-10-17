//! Optional live text only; never retained in JSON reports or replay evidence.
use crypto_lib::binary_ecc::BinaryPoint;
use crypto_lib::cryptanalysis::koblitz_index_calculus::{
    FrobeniusFactorBase, KoblitzIcEvent, PdpOutcome, QueryTrace,
};
use std::{fmt, io::Write, sync::OnceLock, time::Instant};

static START: OnceLock<(Instant, u8)> = OnceLock::new();

pub fn init(verbosity: u8) {
    if verbosity > 0 {
        let _ = START.set((Instant::now(), verbosity));
    }
}

pub fn enabled() -> bool {
    START.get().is_some()
}

pub fn level() -> u8 {
    START.get().map_or(0, |(_, level)| *level)
}

/// One lock and flush per line keeps parallel workers' messages intact and
/// makes redirected output visible immediately, including through `tee`.
pub fn line(message: fmt::Arguments<'_>) {
    if let Some((start, _)) = START.get() {
        let mut out = std::io::stdout().lock();
        let _ = writeln!(
            out,
            "[{:9.3}s] ic: {message}",
            start.elapsed().as_secs_f64()
        );
        let _ = out.flush();
    }
}

fn point(point: &BinaryPoint) -> String {
    match point {
        BinaryPoint::Infinity => "O".into(),
        BinaryPoint::Affine { x, y } => format!(
            "(0x{}, 0x{})",
            x.to_biguint().to_str_radix(16),
            y.to_biguint().to_str_radix(16)
        ),
    }
}

/// Print the same coefficients/witness the producer just used. Never recover
/// coefficients from a seed, rerun an oracle, or claim the scalar is verified.
pub fn query(
    context: fmt::Arguments<'_>,
    trace: QueryTrace<'_>,
    max_trials: u64,
    fb: &FrobeniusFactorBase,
) {
    if !enabled() || (trace.outcome != PdpOutcome::Witness && level() < 3) {
        return;
    }
    let equation = match trace
        .points
        .filter(|_| trace.outcome == PdpOutcome::Witness)
    {
        Some(points) => points
            .iter()
            .map(|i| format!("F[{i}]"))
            .collect::<Vec<_>>()
            .join(" + "),
        None => String::new(),
    };
    let b = trace.b.to_string();
    let target_term = if b == "0" {
        String::new()
    } else {
        format!(" + [{b}]Q")
    };
    line(format_args!(
        "{context} {}/{}: {}; R = [{}]G{}{}{}",
        trace.trial + 1,
        max_trials,
        if trace.outcome == PdpOutcome::Witness {
            "relation found (pending verification)"
        } else {
            outcome(trace.outcome)
        },
        trace.a,
        target_term,
        if equation.is_empty() { "" } else { " = " },
        equation
    ));
    if level() >= 2 {
        if let Some(points) = trace
            .points
            .filter(|_| trace.outcome == PdpOutcome::Witness)
        {
            for &index in points {
                if let Some(value) = fb.points.get(index) {
                    line(format_args!(
                        "{context} {}: F[{index}] = {}",
                        trace.trial + 1,
                        point(value)
                    ));
                }
            }
        }
    }
}

pub fn inputs(context: fmt::Arguments<'_>, generator: &BinaryPoint, target: Option<&BinaryPoint>) {
    if level() >= 2 {
        line(format_args!("{context}: G = {}", point(generator)));
        if let Some(target) = target {
            line(format_args!("{context}: Q = {}", point(target)));
        }
    }
}

pub fn pipeline<E: fmt::Debug>(
    event: crypto_lib::cryptanalysis::ic_framework::PipelineEvent<'_, E>,
    max_trials: u64,
    cofactor: u64,
    modulus: u64,
) {
    use crypto_lib::cryptanalysis::{ic_boundary::Relation, ic_framework::PipelineEvent::*};
    match event {
        FactorBaseStarted => line(format_args!("bench: factor base: building")),
        FactorBaseReady { points, columns } => {
            line(format_args!("bench: factor base: ready; {} points, {columns} columns", points.len()));
            if level() >= 2 {
                for (i, point) in points.iter().enumerate() {
                    line(format_args!("bench: F[{i}] = {point:?}"));
                }
            }
        }
        PrecomputeStarted => line(format_args!("bench: oracle precomputation: building")),
        PrecomputeFinished => line(format_args!("bench: oracle precomputation: ready; collection and incremental linear algebra started")),
        Query { trial, a, b, point, relations } => {
            if relations.is_empty() {
                if level() >= 3 {
                    // A frontend without a completeness certificate has only
                    // failed to return a relation. Do not print "refuted".
                    line(format_args!("bench: query {} (trial budget {max_trials}): no relation returned; R = [{a}]G + [{b}]Q = {point:?}", trial + 1));
                }
            } else {
                for relation in relations {
                    match relation {
                        Relation::Summands(points) => {
                            let terms = points.iter().map(|i| format!("F[{i}]")).collect::<Vec<_>>().join(" + ");
                            line(format_args!("bench: query {} (trial budget {max_trials}): relation found (pending verification); R = [{a}]G + [{b}]Q = {terms}", trial + 1));
                        }
                        Relation::Row { columns, a, b } => {
                            let terms = columns.iter().map(|(i, c)| format!("{c}*L[{i}]")).collect::<Vec<_>>().join(" + ");
                            // A merged large-prime row carries its OWN
                            // coefficients and projected-column logarithms.
                            line(format_args!("bench: query {} (trial budget {max_trials}): merged relation row; {cofactor}*{a} + {cofactor}*{b}*d = {terms} (mod {modulus}); L[i] = log_G(projected column i)", trial + 1));
                        }
                    }
                }
                if level() >= 2 {
                    line(format_args!("bench: query {}: R = {point:?}", trial + 1));
                }
            }
        }
        Finished { trials, relations, rank, verified } => line(format_args!("bench: collection, matrix and verification finished; {trials} trials, {relations} relations, rank {rank}; verified={verified}")),
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
        } if level() >= 3 => line(format_args!(
            "[2/4] batch: starting trials {first_trial}..={}",
            first_trial + trials - 1
        )),
        RelationAttemptFinished { trial, conflicts, collected, .. } if level() >= 3 => line(format_args!("[2/4] trial {trial} details: {collected} relations collected; {conflicts} solver conflicts")),
        RelationBatchStarted { .. } | RelationAttemptFinished { .. } => {}
        RelationProgress {
            collected,
            wanted,
            trials,
        } if level() >= 2 => line(format_args!(
            "[2/4] collection: {collected}/{wanted} relations after {trials} trials"
        )),
        RelationProgress { .. } => {}
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

//! Exact modular target-span diagnostic on the frozen compact-orbit holdout.
//! Usage: koblitz_target_span_gate <repository-root> <output.json>
//! No producer runs or target selection occur here; see the preregistered
//! experiments/koblitz-target-span-gate-20261006/PROTOCOL.md.

use serde_json::{json, Value};
use std::fs;
use std::path::{Path, PathBuf};

fn one_json(path: &Path) -> Value {
    let body = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut lines = body.lines().filter(|line| !line.is_empty());
    let first = lines.next().expect("empty JSONL file");
    assert!(
        lines.next().is_none(),
        "expected one row: {}",
        path.display()
    );
    serde_json::from_str(first).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

fn json_lines(path: &Path) -> Vec<Value> {
    let body = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    body.lines()
        .filter(|line| !line.is_empty())
        .map(|line| {
            serde_json::from_str(line).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        })
        .collect()
}

fn u64_at(row: &Value, key: &str) -> u64 {
    row[key]
        .as_u64()
        .unwrap_or_else(|| panic!("missing integer {key}"))
}

fn str_at<'a>(row: &'a Value, key: &str) -> &'a str {
    row[key]
        .as_str()
        .unwrap_or_else(|| panic!("missing string {key}"))
}

fn mul(a: u64, b: u64, modulus: u64) -> u64 {
    ((a as u128 * b as u128) % modulus as u128) as u64
}

fn sub(a: u64, b: u64, modulus: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        modulus - (b - a)
    }
}

fn power(mut base: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut value = 1;
    while exponent != 0 {
        if exponent & 1 != 0 {
            value = mul(value, base, modulus);
        }
        base = mul(base, base, modulus);
        exponent >>= 1;
    }
    value
}

struct RowSpan {
    modulus: u64,
    columns: usize,
    pivots: Vec<Option<Vec<u64>>>,
    rank: usize,
}

impl RowSpan {
    fn new(modulus: u64, columns: usize) -> Self {
        Self {
            modulus,
            columns,
            pivots: vec![None; columns],
            rank: 0,
        }
    }

    fn insert(&mut self, mut row: Vec<u64>) -> bool {
        assert_eq!(row.len(), self.columns + 1);
        assert!(row.iter().all(|&value| value < self.modulus));
        for column in 0..self.columns {
            if row[column] == 0 {
                continue;
            }
            if let Some(pivot) = &self.pivots[column] {
                let coefficient = row[column];
                for (value, &p) in row.iter_mut().zip(pivot).skip(column) {
                    *value = sub(*value, mul(coefficient, p, self.modulus), self.modulus);
                }
            } else {
                let inverse = power(row[column], self.modulus - 2, self.modulus);
                assert_eq!(mul(row[column], inverse, self.modulus), 1);
                for value in row.iter_mut().skip(column) {
                    *value = mul(*value, inverse, self.modulus);
                }
                self.pivots[column] = Some(row);
                self.rank += 1;
                return true;
            }
        }
        assert_eq!(row[self.columns], 0, "inconsistent rank equation");
        false
    }

    // If the target vector is in the row span, its RHS is the unique
    // combination of recorded scalar RHS values implied by these rows.
    fn target_log(&self, target: &[u64]) -> Option<u64> {
        assert_eq!(target.len(), self.columns);
        let mut remainder = target.to_vec();
        let mut recovered = 0;
        for column in 0..self.columns {
            let coefficient = remainder[column];
            if coefficient == 0 {
                continue;
            }
            let pivot = self.pivots[column].as_ref()?;
            recovered =
                (recovered + mul(coefficient, pivot[self.columns], self.modulus)) % self.modulus;
            for (value, &p) in remainder.iter_mut().zip(pivot).skip(column) {
                *value = sub(*value, mul(coefficient, p, self.modulus), self.modulus);
            }
        }
        assert!(remainder.iter().all(|&value| value == 0));
        Some(recovered)
    }
}

fn read_vec(row: &Value, key: &str) -> Vec<u64> {
    serde_json::from_value(row[key].clone()).unwrap_or_else(|_| panic!("invalid {key}"))
}

fn analyze_run(root: &Path, n: u64, columns: usize, repeat: u64) -> Value {
    let dir = root
        .join("experiments/koblitz-base-size-cold-panel-20261004/runs/holdout")
        .join(format!("n{n}/k{columns}/r{repeat}"));
    let base = one_json(&dir.join("base.jsonl"));
    let target = one_json(&dir.join("ic_targets.jsonl"));
    let replay = one_json(&dir.join("replay.jsonl"));
    let status = one_json(&dir.join("status.jsonl"));
    let trace = json_lines(&dir.join("rank.jsonl"));
    assert!(trace.len() >= 3, "incomplete rank trace: {}", dir.display());
    let header = &trace[0];
    let solution = trace.last().unwrap();
    assert_eq!(str_at(header, "kind"), "compact_orbit_rank_header");
    assert_eq!(str_at(solution, "kind"), "compact_orbit_rank_solution");
    assert_eq!(u64_at(header, "n"), n);
    assert_eq!(u64_at(header, "orbit_columns") as usize, columns);
    assert_eq!(u64_at(&base, "n"), n);
    assert_eq!(u64_at(&base, "orbit_columns") as usize, columns);
    assert_eq!(u64_at(&target, "n"), n);
    assert_eq!(u64_at(&replay, "n"), n);
    assert_eq!(u64_at(&status, "n"), n);
    assert_eq!(u64_at(&status, "K") as usize, columns);
    assert_eq!(u64_at(&status, "repeat"), repeat);
    for field in ["rho_status", "ic_status", "replay_status"] {
        assert_eq!(u64_at(&status, field), 0);
    }
    assert_eq!(replay["verified"], true);
    assert_eq!(replay["target_relation_replayed"], true);
    assert_eq!(target["group_verified"], true);
    assert_eq!(u64_at(&target, "exit_code"), 0);
    assert_eq!(str_at(header, "base_hash"), str_at(&base, "base_hash"));
    assert_eq!(u64_at(&replay, "orbit_columns") as usize, columns);
    let modulus = u64_at(header, "subgroup_order");
    assert_eq!(u64_at(&replay, "subgroup_order"), modulus);
    let labels: Vec<[u64; 2]> =
        serde_json::from_value(base["factor_base_point_labels"].clone()).expect("base labels");
    assert_eq!(labels.len() as u64, u64_at(header, "factor_base_points"));
    let target_indices = read_vec(&target, "point_indices");
    assert_eq!(target_indices.len(), 4);
    let mut target_row = vec![0; columns];
    for index in target_indices {
        let [column, coefficient] = labels[index as usize];
        assert!((column as usize) < columns && coefficient < modulus);
        target_row[column as usize] = (target_row[column as usize] + coefficient) % modulus;
    }
    let producer_scalar = u64_at(&target, "recovered_scalar");
    assert_eq!(u64_at(&replay, "target_scalar"), producer_scalar);
    let logs = read_vec(solution, "logs");
    assert_eq!(logs.len(), columns);
    assert_eq!(
        target_row
            .iter()
            .zip(&logs)
            .fold(0, |sum, (&a, &b)| (sum + mul(a, b, modulus)) % modulus),
        producer_scalar
    );

    let mut span = RowSpan::new(modulus, columns);
    let mut first_span_rank = None;
    let mut probes_total = 0;
    let mut probes_at_first_span = None;
    let mut attempt_count = 0;
    for attempt in &trace[1..trace.len() - 1] {
        assert_eq!(str_at(attempt, "kind"), "compact_orbit_rank_attempt");
        assert_eq!(u64_at(attempt, "attempt_index"), attempt_count);
        assert_eq!(u64_at(attempt, "rank_before") as usize, span.rank);
        attempt_count += 1;
        assert_eq!(attempt["found"], true, "unexpected frozen rank miss");
        let row = read_vec(attempt, "row");
        assert_eq!(row[columns], u64_at(attempt, "scalar"));
        assert_eq!(
            row[..columns]
                .iter()
                .zip(&logs)
                .fold(0, |sum, (&a, &b)| (sum + mul(a, b, modulus)) % modulus),
            row[columns]
        );
        let gained = span.insert(row);
        assert_eq!(attempt["gained"], gained);
        assert_eq!(u64_at(attempt, "rank_after") as usize, span.rank);
        probes_total += u64_at(attempt, "probes");
        if first_span_rank.is_none() {
            if let Some(scalar) = span.target_log(&target_row) {
                assert_eq!(scalar, producer_scalar);
                first_span_rank = Some(span.rank);
                probes_at_first_span = Some(probes_total);
            }
        }
    }
    assert_eq!(attempt_count, u64_at(solution, "attempts"));
    assert_eq!(span.rank, columns);
    assert_eq!(u64_at(solution, "rank") as usize, columns);
    assert_eq!(u64_at(solution, "failures"), 0);
    assert_eq!(u64_at(solution, "relations"), columns as u64);
    assert_eq!(span.target_log(&target_row), Some(producer_scalar));
    let first = first_span_rank.expect("full rank did not span target");
    let prefix_probes = probes_at_first_span.unwrap();
    json!({
        "n":n, "K":columns, "repeat":repeat,
        "target":target["target"], "base_hash":str_at(&base,"base_hash"),
        "actual_base_points":labels.len(),
        "first_target_span_rank":first,
        "full_rank":span.rank,
        "proper_prefix":first < columns,
        "rank_rows_saved_counterfactual":columns - first,
        "rank_probes_total":probes_total,
        "rank_probes_prefix":prefix_probes,
        "rank_probes_suffix_counterfactual":probes_total - prefix_probes,
        "recovered_scalar":producer_scalar,
        "replay_verified":true
    })
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert_eq!(
        args.len(),
        3,
        "usage: koblitz_target_span_gate <repo-root> <output.json>"
    );
    let root = PathBuf::from(&args[1]);
    let mut runs = Vec::new();
    let mut cells = Vec::new();
    for (n, choices) in [(41, [64, 85]), (53, [160, 220])] {
        for columns in choices {
            let cell_runs: Vec<Value> = (1..=6)
                .map(|repeat| analyze_run(&root, n, columns, repeat))
                .collect();
            let reference = &cell_runs[0];
            for row in &cell_runs[1..] {
                for field in [
                    "target",
                    "base_hash",
                    "actual_base_points",
                    "recovered_scalar",
                    "rank_probes_total",
                ] {
                    assert_eq!(row[field], reference[field], "repeat changed {field}");
                }
            }
            let threshold = columns * 4 / 5;
            let admitted = cell_runs
                .iter()
                .all(|row| u64_at(row, "first_target_span_rank") as usize <= threshold);
            let first_ranks: Vec<u64> = cell_runs
                .iter()
                .map(|row| u64_at(row, "first_target_span_rank"))
                .collect();
            cells.push(json!({
                "n":n, "K":columns, "distinct_targets":1,
                "repeats":6, "admission_prefix_limit":threshold,
                "first_target_span_ranks":first_ranks,
                "admit_new_producer":admitted
            }));
            runs.extend(cell_runs);
        }
    }
    let admitted = cells.iter().any(|cell| cell["admit_new_producer"] == true);
    let output = json!({
        "kind":"compact_orbit_frozen_target_span_gate",
        "schema_version":1,
        "input":"koblitz-base-size-cold-panel-20261004 held-out rank and target transcripts",
        "run_count":runs.len(), "distinct_point_base_cells":cells.len(),
        "admit_new_producer":admitted,
        "cells":cells, "runs":runs
    });
    fs::write(&args[2], format!("{output}\n")).expect("write result");
    println!(
        "admit_new_producer={admitted}; checked {} frozen runs",
        runs.len()
    );
}

#[cfg(test)]
mod tests {
    use super::RowSpan;

    #[test]
    fn detects_a_proper_prefix_and_rejects_a_missing_pivot() {
        let mut span = RowSpan::new(7, 3);
        assert_eq!(span.target_log(&[2, 1, 0]), None);
        assert!(span.insert(vec![1, 0, 0, 3]));
        assert_eq!(span.target_log(&[2, 0, 0]), Some(6));
        assert_eq!(span.target_log(&[2, 1, 0]), None);
        assert!(span.insert(vec![0, 1, 0, 4]));
        assert_eq!(span.target_log(&[2, 1, 0]), Some(3));
        assert_eq!(span.target_log(&[0, 0, 1]), None);
    }
}

//! Complete user-space instruction count for one `ecbench exec` solve.
//!
//! The opt-in child emits a reset marker immediately before `methods::solve`
//! and another immediately after it, including an error return. Callgrind
//! writes one part for each marker, while `ic_measurement` may write more
//! parts inside the solve. Sum precisely the parts after the first marker
//! through the second; the first part contains input and workload setup.
//! `Ir` is Callgrind's simulated instruction count, not a PMU retired-
//! instruction counter or a wall-time speedup.

use std::path::{Path, PathBuf};

use serde::Serialize;

const BEFORE: &str = "ecbench_before_solve";
const AFTER: &str = "ecbench_solve";
const BEFORE_ONLINE: &str = "ecbench_before_online";
const AFTER_ONLINE: &str = "ecbench_online";

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct Part {
    pub file: String,
    pub trigger: String,
    pub ir: u64,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct SolveIr {
    pub schema: &'static str,
    pub unit: &'static str,
    pub solve_ir: u64,
    pub prework_ir: u64,
    pub solve_parts: usize,
    pub parts: Vec<Part>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct OnlineIr {
    pub schema: &'static str,
    pub unit: &'static str,
    pub pre_online_ir: u64,
    pub online_ir: u64,
    pub post_online_ir: u64,
    pub solve_ir: u64,
    pub prework_ir: u64,
    pub parts: Vec<Part>,
}

fn sum_ir(parts: &[Part]) -> Result<u64, String> {
    parts.iter().try_fold(0u64, |sum, part| {
        sum.checked_add(part.ir)
            .ok_or_else(|| "Callgrind Ir total overflowed u64".to_string())
    })
}

fn parse_part(path: &Path, content: &str) -> Result<Part, String> {
    let mut trigger = None;
    let mut events: Vec<&str> = Vec::new();
    let mut totals: Option<Vec<u64>> = None;
    for line in content.lines() {
        if let Some(value) = line.strip_prefix("desc: Trigger:") {
            trigger = Some(
                value
                    .trim()
                    .strip_prefix("Client Request: ")
                    .unwrap_or(value.trim())
                    .to_string(),
            );
        } else if let Some(value) = line.strip_prefix("events:") {
            events = value.split_whitespace().collect();
        } else if let Some(value) = line
            .strip_prefix("summary:")
            .or_else(|| line.strip_prefix("totals:"))
        {
            totals = Some(
                value
                    .split_whitespace()
                    .map(|s| {
                        s.parse::<u64>().map_err(|_| {
                            format!("{} has invalid event count `{s}`", path.display())
                        })
                    })
                    .collect::<Result<Vec<_>, _>>()?,
            );
        }
    }
    let index = events
        .iter()
        .position(|event| *event == "Ir")
        .ok_or_else(|| format!("{} has no Ir event", path.display()))?;
    let ir = totals
        .and_then(|values| values.get(index).copied())
        .ok_or_else(|| format!("{} has no Ir total", path.display()))?;
    Ok(Part {
        file: path
            .file_name()
            .ok_or_else(|| format!("{} has no file name", path.display()))?
            .to_string_lossy()
            .into_owned(),
        trigger: trigger.ok_or_else(|| format!("{} has no trigger", path.display()))?,
        ir,
    })
}

fn summarize(parts: Vec<Part>) -> Result<SolveIr, String> {
    let before: Vec<_> = parts
        .iter()
        .enumerate()
        .filter_map(|(i, part)| (part.trigger == BEFORE).then_some(i))
        .collect();
    let after: Vec<_> = parts
        .iter()
        .enumerate()
        .filter_map(|(i, part)| (part.trigger == AFTER).then_some(i))
        .collect();
    if before.len() != 1 || after.len() != 1 || before[0] != 0 || before[0] >= after[0] {
        return Err(format!(
            "expected exactly one ordered {BEFORE}/{AFTER} marker; found {before:?}/{after:?}"
        ));
    }
    let solve_parts = after[0] - before[0];
    let solve_ir = sum_ir(&parts[before[0] + 1..=after[0]])?;
    Ok(SolveIr {
        schema: "ecbench.callgrind_ir/v1",
        unit: "callgrind.Ir",
        solve_ir,
        prework_ir: parts[before[0]].ir,
        solve_parts,
        parts,
    })
}

fn online_from_solve(complete: SolveIr) -> Result<OnlineIr, String> {
    let parts = &complete.parts;
    let before: Vec<_> = parts
        .iter()
        .enumerate()
        .filter_map(|(i, part)| (part.trigger == BEFORE_ONLINE).then_some(i))
        .collect();
    let after: Vec<_> = parts
        .iter()
        .enumerate()
        .filter_map(|(i, part)| (part.trigger == AFTER_ONLINE).then_some(i))
        .collect();
    let solve_end = parts
        .iter()
        .position(|part| part.trigger == AFTER)
        .ok_or("missing outer solve marker")?;
    if before.len() != 1
        || after.len() != 1
        || before[0] == 0
        || before[0] >= after[0]
        || after[0] >= solve_end
    {
        return Err(format!(
            "expected exactly one ordered online interval inside the solve; found {before:?}/{after:?}, solve end {solve_end}"
        ));
    }
    let pre_online_ir = sum_ir(&parts[1..=before[0]])?;
    let online_ir = sum_ir(&parts[before[0] + 1..=after[0]])?;
    let post_online_ir = sum_ir(&parts[after[0] + 1..=solve_end])?;
    if pre_online_ir
        .checked_add(online_ir)
        .and_then(|sum| sum.checked_add(post_online_ir))
        != Some(complete.solve_ir)
    {
        return Err("online instruction phases do not sum to complete solve".into());
    }
    Ok(OnlineIr {
        schema: "ecbench.callgrind_online_ir/v1",
        unit: "callgrind.Ir",
        pre_online_ir,
        online_ir,
        post_online_ir,
        solve_ir: complete.solve_ir,
        prework_ir: complete.prework_ir,
        parts: complete.parts,
    })
}

/// Read `<prefix>.1`, `<prefix>.2`, ... in numeric order. The terminal
/// `<prefix>` file is deliberately excluded because it follows the solve.
pub fn solve_ir(prefix: &Path) -> Result<SolveIr, String> {
    let dir = prefix
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let base = prefix
        .file_name()
        .ok_or_else(|| "Callgrind prefix needs a file name".to_string())?
        .to_string_lossy();
    let mut numbered: Vec<(u64, PathBuf)> = Vec::new();
    for entry in std::fs::read_dir(dir).map_err(|e| format!("{}: {e}", dir.display()))? {
        let path = entry.map_err(|e| e.to_string())?.path();
        let Some(name) = path.file_name().and_then(|n| n.to_str()) else {
            continue;
        };
        let Some(suffix) = name.strip_prefix(&format!("{base}.")) else {
            continue;
        };
        if !suffix.is_empty() && suffix.bytes().all(|b| b.is_ascii_digit()) {
            let n = suffix.parse::<u64>().map_err(|e| format!("{name}: {e}"))?;
            numbered.push((n, path));
        }
    }
    numbered.sort_by_key(|(n, _)| *n);
    if numbered.is_empty() {
        return Err(format!(
            "{} has no numbered Callgrind parts",
            prefix.display()
        ));
    }
    if numbered[0].0 != 1 {
        return Err("first Callgrind part number is not 1".into());
    }
    for pair in numbered.windows(2) {
        if pair[0].0 == pair[1].0 {
            return Err("duplicate Callgrind part number".into());
        }
        if pair[0].0 + 1 != pair[1].0 {
            return Err("missing Callgrind part number".into());
        }
    }
    let parts = numbered
        .into_iter()
        .map(|(_, path)| {
            let content =
                std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
            parse_part(&path, &content)
        })
        .collect::<Result<Vec<_>, _>>()?;
    summarize(parts)
}

/// Extract the one-target online interval from a complete solve profile.
/// Missing, duplicate, reversed or out-of-solve markers are errors.
pub fn online_ir(prefix: &Path) -> Result<OnlineIr, String> {
    online_from_solve(solve_ir(prefix)?)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn part(name: &str, trigger: &str, ir: u64) -> Part {
        Part {
            file: name.into(),
            trigger: trigger.into(),
            ir,
        }
    }

    #[test]
    fn counts_only_the_bracketed_parts_including_internal_phase_dumps() {
        let result = summarize(vec![
            part("p.1", BEFORE, 100),
            part("p.2", "generic_ic_setup", 7),
            part("p.3", "target_query", 11),
            part("p.4", AFTER, 13),
        ])
        .unwrap();
        assert_eq!(result.prework_ir, 100);
        assert_eq!(result.solve_ir, 31);
        assert_eq!(result.solve_parts, 3);
    }

    #[test]
    fn separates_reusable_and_one_target_instructions_without_losing_work() {
        let complete = summarize(vec![
            part("p.1", BEFORE, 100),
            part("p.2", BEFORE_ONLINE, 10),
            part("p.3", "generic_ic_target_query", 7),
            part("p.4", AFTER_ONLINE, 13),
            part("p.5", AFTER, 3),
        ])
        .unwrap();
        let online = online_from_solve(complete).unwrap();
        assert_eq!(online.pre_online_ir, 10);
        assert_eq!(online.online_ir, 20);
        assert_eq!(online.post_online_ir, 3);
        assert_eq!(online.solve_ir, 33);
        assert_eq!(online.prework_ir, 100);
    }

    #[test]
    fn refuses_missing_duplicate_reversed_or_outside_online_markers() {
        let base = vec![part("p.1", BEFORE, 1), part("p.2", AFTER, 2)];
        assert!(online_from_solve(summarize(base).unwrap()).is_err());
        let reversed = vec![
            part("p.1", BEFORE, 1),
            part("p.2", AFTER_ONLINE, 2),
            part("p.3", BEFORE_ONLINE, 3),
            part("p.4", AFTER, 4),
        ];
        assert!(online_from_solve(summarize(reversed).unwrap()).is_err());
        let duplicate = vec![
            part("p.1", BEFORE, 1),
            part("p.2", BEFORE_ONLINE, 2),
            part("p.3", BEFORE_ONLINE, 3),
            part("p.4", AFTER_ONLINE, 4),
            part("p.5", AFTER, 5),
        ];
        assert!(online_from_solve(summarize(duplicate).unwrap()).is_err());
        let outside = vec![
            part("p.1", BEFORE, 1),
            part("p.2", BEFORE_ONLINE, 2),
            part("p.3", AFTER, 3),
            part("p.4", AFTER_ONLINE, 4),
        ];
        assert!(online_from_solve(summarize(outside).unwrap()).is_err());
    }

    #[test]
    fn refuses_an_incomplete_or_reversed_profile() {
        assert!(summarize(vec![part("p.1", BEFORE, 100)]).is_err());
        assert!(summarize(vec![part("p.1", AFTER, 3), part("p.2", BEFORE, 4)]).is_err());
        assert!(summarize(vec![
            part("p.1", BEFORE, 1),
            part("p.2", AFTER, 2),
            part("p.3", AFTER, 3),
        ])
        .is_err());
    }

    #[test]
    fn parses_the_ir_column_of_a_callgrind_part() {
        let input =
            "desc: Trigger: Client Request: ecbench_solve\nevents: Ir Dr Dw\nsummary: 42 9 8\n";
        assert_eq!(parse_part(Path::new("p.2"), input).unwrap().ir, 42);
        assert!(parse_part(Path::new("p.2"), "events: Dr\nsummary: 9\n").is_err());
    }
}

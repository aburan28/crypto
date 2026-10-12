//! A callgrind profile split by measurement phase: the retired
//! `harness/callgrind_phases.py`, native, with the same output byte for
//! byte.
//!
//! `ic_measurement` issues a valgrind `DUMP_STATS_AT` request at every
//! phase boundary, labelled with the phase that just ended, so a run leaves
//! one part file per segment (`<out>.1`, `<out>.2`, ...).  Each part's
//! trigger names its phase and its `summary:` or `totals:` line holds its
//! counts.  This sums the parts by label.  The first `generic_ic_setup`
//! part is also annotated by function, exclusive costs.  A run with no
//! parts (`ic workflow` dumps none) is read from its base file as one part
//! and annotated whole.
//!
//! The annotation is `callgrind_annotate`'s (valgrind's own tool); only its
//! text is read here.

use std::path::{Path, PathBuf};
use std::process::Command;

use super::json::J;

/// A part's trigger and its counts, event by event.
fn read_part(path: &Path) -> Result<(String, Vec<(String, i128)>), String> {
    let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let text = String::from_utf8_lossy(&bytes);
    let (mut events, mut totals, mut trigger) = (Vec::new(), None, "unknown".to_string());
    // `^key:\s+(.*)$`: at least one space after the key, then the rest.
    let after = |line: &str, key: &str| -> Option<String> {
        let rest = line.strip_prefix(key)?;
        rest.starts_with(char::is_whitespace)
            .then(|| rest.trim_start().to_string())
    };
    for line in text.lines() {
        if let Some(rest) = after(line, "desc: Trigger:") {
            let t = rest.trim();
            trigger = t.strip_prefix("Client Request: ").unwrap_or(t).to_string();
        } else if let Some(rest) = after(line, "events:") {
            events = rest.split_whitespace().map(str::to_string).collect();
        } else if let Some(rest) = after(line, "summary:").or_else(|| after(line, "totals:")) {
            totals = Some(
                rest.split_whitespace()
                    .map(|x| {
                        x.parse::<i128>()
                            .map_err(|_| format!("{}: `{x}` is not a count", path.display()))
                    })
                    .collect::<Result<Vec<_>, _>>()?,
            );
        }
    }
    let totals = totals.unwrap_or_default();
    let counts = events
        .into_iter()
        .enumerate()
        .map(|(i, e)| (e, totals.get(i).copied().unwrap_or(0)))
        .collect();
    Ok((trigger, counts))
}

/// The digits after a part file's last dot, as a number (0 for none).
fn part_number(path: &Path) -> u64 {
    path.file_name()
        .and_then(|n| n.to_str())
        .and_then(|n| n.rsplit_once('.'))
        .and_then(|(_, s)| s.parse().ok())
        .unwrap_or(0)
}

fn is_digits(s: &str) -> bool {
    !s.is_empty() && s.chars().all(|c| c.is_ascii_digit())
}

/// `re.sub(r"\(\s*[\d.]+%\)", " ", line)`: every percentage in
/// parentheses becomes one space.
fn strip_percentages(line: &str) -> String {
    let chars: Vec<char> = line.chars().collect();
    let mut out = String::with_capacity(line.len());
    let mut i = 0;
    while i < chars.len() {
        if chars[i] == '(' {
            let mut j = i + 1;
            while j < chars.len() && chars[j].is_whitespace() {
                j += 1;
            }
            let digits = j;
            while j < chars.len() && (chars[j].is_ascii_digit() || chars[j] == '.') {
                j += 1;
            }
            if j > digits && chars.get(j) == Some(&'%') && chars.get(j + 1) == Some(&')') {
                out.push(' ');
                i = j + 2;
                continue;
            }
        }
        out.push(chars[i]);
        i += 1;
    }
    out
}

/// `re.sub(r"\s*\[.*\]$", "", name)`: from the first `[`, with the spaces
/// before it, to the end, when the name ends in `]`.  (Inside a name such
/// as `<[T] as Trait>` that cuts at the name's own bracket, as the script
/// did.)
fn strip_object(name: &str) -> String {
    match name.find('[') {
        Some(at) if name.ends_with(']') => name[..at].trim_end().to_string(),
        _ => name.to_string(),
    }
}

/// `callgrind_annotate`'s function table, read as the script read it: the
/// first `top` rows with a count for every event shown.
pub fn parse_annotate(text: &str, top: usize) -> Result<Vec<J>, String> {
    let (mut rows, mut events, mut in_table) = (Vec::new(), Vec::<String>::new(), false);
    for line in text.lines() {
        if let Some((_, shown)) = line
            .starts_with("Events shown:")
            .then(|| line.split_once(':'))
            .flatten()
        {
            events = shown.split_whitespace().map(str::to_string).collect();
        }
        if line.contains("file:function") {
            in_table = true;
            continue;
        }
        if !in_table || line.trim().is_empty() || line.starts_with('-') {
            continue;
        }
        let stripped = strip_percentages(line);
        let cells: Vec<&str> = stripped.split_whitespace().collect();
        let mut numbers: Vec<i128> = Vec::new();
        let mut i = 0;
        while i < cells.len() && numbers.len() < events.len() {
            let c = cells[i];
            if c == "." {
                numbers.push(0);
            } else if !c.is_empty() && c.chars().all(|ch| ch.is_ascii_digit() || ch == ',') {
                let digits = c.replace(',', "");
                numbers.push(
                    digits
                        .parse()
                        .map_err(|_| format!("`{c}` in callgrind_annotate's table"))?,
                );
            } else {
                break;
            }
            i += 1;
        }
        if numbers.len() < events.len() {
            continue;
        }
        let name = strip_object(&cells[i..].join(" "));
        let function = match name.strip_prefix("???:") {
            Some(rest) => rest.to_string(),
            None => name,
        };
        let mut kv = vec![("function".to_string(), J::Str(function))];
        kv.extend(
            events
                .iter()
                .zip(&numbers)
                .map(|(e, &v)| (e.clone(), J::Int(v))),
        );
        rows.push(J::Obj(kv));
        if rows.len() >= top {
            break;
        }
    }
    Ok(rows)
}

fn annotate(part: &Path, top: usize) -> Result<Vec<J>, String> {
    let out = Command::new("callgrind_annotate")
        .args(["--inclusive=no", "--threshold=99"])
        .arg(part)
        .output()
        .map_err(|e| format!("callgrind_annotate: {e}"))?;
    parse_annotate(&String::from_utf8_lossy(&out.stdout), top)
}

fn counts_obj(counts: &[(String, i128)]) -> J {
    J::Obj(
        counts
            .iter()
            .map(|(e, v)| (e.clone(), J::Int(*v)))
            .collect(),
    )
}

/// The phases of the profile `<dir>/<stem>`, as `callgrind_phases.py`
/// printed them.
pub fn phases(dir: &Path, stem: &str, top: usize) -> Result<J, String> {
    let prefix = format!("{stem}.");
    let mut parts: Vec<PathBuf> = std::fs::read_dir(dir)
        .map_err(|e| format!("{}: {e}", dir.display()))?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| {
            p.file_name().and_then(|n| n.to_str()).is_some_and(|n| {
                n.starts_with(&prefix) && n.rsplit_once('.').is_some_and(|(_, s)| is_digits(s))
            })
        })
        .collect();
    parts.sort_by_key(|p| part_number(p));
    let mut whole = None;
    if parts.is_empty() && dir.join(stem).exists() {
        // A command without phase dumps (`ic workflow`): the base file is
        // the whole run, as one part.
        whole = Some(dir.join(stem));
        parts.push(dir.join(stem));
    }
    // Label → (parts, counts), in the order labels first appear.
    let mut by_label: Vec<(String, i128, Vec<(String, i128)>)> = Vec::new();
    let mut first_setup: Option<PathBuf> = None;
    for p in &parts {
        let (label, counts) = read_part(p)?;
        let slot = match by_label.iter().position(|(l, _, _)| *l == label) {
            Some(i) => &mut by_label[i],
            None => {
                by_label.push((label.clone(), 0, Vec::new()));
                by_label.last_mut().expect("just pushed")
            }
        };
        slot.1 += 1;
        for (e, v) in counts {
            match slot.2.iter_mut().find(|(k, _)| *k == e) {
                Some((_, sum)) => *sum += v,
                None => slot.2.push((e, v)),
            }
        }
        if first_setup.is_none() && label == "generic_ic_setup" {
            first_setup = Some(p.clone());
        }
    }
    let mut total: Vec<(String, i128)> = Vec::new();
    for (_, _, counts) in &by_label {
        for (e, v) in counts {
            match total.iter_mut().find(|(k, _)| k == e) {
                Some((_, sum)) => *sum += v,
                None => total.push((e.clone(), *v)),
            }
        }
    }
    let ir = |counts: &[(String, i128)]| {
        counts
            .iter()
            .find(|(e, _)| e == "Ir")
            .map_or(0, |(_, v)| *v)
    };
    // A stable sort by descending instructions, as the script's.
    by_label.sort_by_key(|(_, _, counts)| std::cmp::Reverse(ir(counts)));
    let by_phase = J::Obj(
        by_label
            .iter()
            .map(|(label, n, counts)| {
                (
                    label.clone(),
                    J::Obj(vec![
                        ("parts".into(), J::Int(*n)),
                        ("counts".into(), counts_obj(counts)),
                    ]),
                )
            })
            .collect(),
    );
    let name = |p: &Path| {
        J::Str(
            p.file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_default(),
        )
    };
    let mut doc = vec![
        ("stem".to_string(), J::Str(stem.into())),
        ("parts".into(), J::Int(parts.len() as i128)),
        ("total".into(), counts_obj(&total)),
        ("by_phase".into(), by_phase),
        (
            "first_setup_part".into(),
            first_setup.as_deref().map_or(J::Null, name),
        ),
        (
            "first_setup_part_counts".into(),
            match &first_setup {
                Some(p) => counts_obj(&read_part(p)?.1),
                None => J::Null,
            },
        ),
        (
            "first_setup_part_top_functions".into(),
            J::Arr(match &first_setup {
                Some(p) => annotate(p, top)?,
                None => Vec::new(),
            }),
        ),
    ];
    if let Some(w) = &whole {
        doc.push(("whole_run_top_functions".into(), J::Arr(annotate(w, top)?)));
    }
    Ok(J::Obj(doc))
}

// ── R02b's callgrind control ────────────────────────────────────────

/// The arithmetic kernels whose instructions must not move: every function
/// of `koblitz_fast`, `Gf2`'s field operations, and the bulk key.
pub const KERNEL_MARKS: [&str; 3] = [
    "koblitz_fast::",
    "semaev_decomp::Gf2::",
    "PairSumTable::keys_of",
];
/// The whole run's instructions must lie this close to the base's.
pub const TOTAL_TOLERANCE: f64 = 5e-3;

struct Profile {
    total_ir: J,
    status: J,
    recovered: Vec<J>,
    /// Function → exclusive instructions, in the profile's order.
    functions: Vec<(String, J)>,
}

impl Profile {
    fn function(&self, name: &str) -> J {
        self.functions
            .iter()
            .find(|(f, _)| f == name)
            .map_or(J::Null, |(_, v)| v.clone())
    }
}

/// R02b's callgrind control (its PROTOCOL.md) from a run tree's
/// `callgrind/`: both arms' instructions per row, the arithmetic kernels
/// function by function, and every other function whose count differs,
/// named.  It is R02b's declared `analyse.py` `callgrind()`, natively.
pub fn control(dir: &Path) -> Result<J, String> {
    let mut tags: Vec<String> = match std::fs::read_dir(dir) {
        Ok(entries) => entries
            .filter_map(|e| e.ok())
            .filter_map(|e| e.file_name().to_str().map(str::to_string))
            .filter(|n| n.ends_with(".phases.json"))
            .collect(),
        Err(_) => Vec::new(),
    };
    tags.sort();
    let mut profiles: Vec<(String, Profile)> = Vec::new();
    for name in tags {
        let tag = name.trim_end_matches(".phases.json").to_string();
        let text = std::fs::read_to_string(dir.join(&name)).map_err(|e| format!("{name}: {e}"))?;
        let doc = if text.is_empty() {
            J::Obj(Vec::new())
        } else {
            super::json::parse(&text).map_err(|e| format!("{name}: {e}"))?
        };
        let wf = super::json::read_opt(&dir.join(format!("{tag}.workflow.json")))?
            .unwrap_or(J::Obj(Vec::new()));
        let mut functions: Vec<(String, J)> = Vec::new();
        for e in doc
            .get("whole_run_top_functions")
            .and_then(J::as_arr)
            .unwrap_or(&[])
        {
            let f = e.at("function")?.as_str().ok_or("`function`")?.to_string();
            let ir = e.at("Ir")?.clone();
            match functions.iter_mut().find(|(k, _)| *k == f) {
                Some((_, v)) => *v = ir,
                None => functions.push((f, ir)),
            }
        }
        profiles.push((
            tag,
            Profile {
                total_ir: doc
                    .get("total")
                    .and_then(|t| t.get("Ir"))
                    .cloned()
                    .unwrap_or(J::Null),
                status: wf.get("status").cloned().unwrap_or(J::Null),
                recovered: wf
                    .get("solutions")
                    .and_then(|s| s.get("items"))
                    .and_then(J::as_arr)
                    .unwrap_or(&[])
                    .iter()
                    .map(|s| s.get("recovered").cloned().unwrap_or(J::Null))
                    .collect(),
                functions,
            },
        ));
    }
    let mut pairs = Vec::new();
    for (tag, base) in &profiles {
        let Some(row) = tag.strip_prefix("base-") else {
            continue;
        };
        let Some((_, cand)) = profiles.iter().find(|(t, _)| *t == format!("cand-{row}")) else {
            continue;
        };
        let (Some(b), Some(c)) = (base.total_ir.as_f64(), cand.total_ir.as_f64()) else {
            continue;
        };
        if !(base.total_ir.truthy() && cand.total_ir.truthy()) {
            continue;
        }
        let mut names: Vec<&String> = base
            .functions
            .iter()
            .chain(&cand.functions)
            .map(|(f, _)| f)
            .collect();
        names.sort();
        names.dedup();
        let kernels: Vec<&String> = names
            .iter()
            .copied()
            .filter(|f| KERNEL_MARKS.iter().any(|m| f.contains(m)))
            .collect();
        let differing: Vec<(String, J)> = names
            .iter()
            .filter_map(|f| {
                let (x, y) = (base.function(f), cand.function(f));
                (!super::json::py_eq(&x, &y)).then(|| ((*f).clone(), J::Arr(vec![x, y])))
            })
            .collect();
        let ratio = b / c;
        let complete = J::Str("complete".into());
        pairs.push((
            row.to_string(),
            J::Obj(vec![
                ("ir_base_over_cand".into(), J::Float(ratio)),
                (
                    "kernels_identical".into(),
                    J::Bool(
                        kernels
                            .iter()
                            .all(|f| !differing.iter().any(|(d, _)| d == *f)),
                    ),
                ),
                (
                    "kernels".into(),
                    J::Arr(kernels.iter().map(|f| J::Str((*f).clone())).collect()),
                ),
                ("differing_functions".into(), J::Obj(differing)),
                (
                    "total_within_tolerance".into(),
                    J::Bool((ratio - 1.0).abs() <= TOTAL_TOLERANCE),
                ),
                (
                    "same_log".into(),
                    J::Bool(
                        !base.recovered.is_empty()
                            && super::json::py_eq(
                                &J::Arr(base.recovered.clone()),
                                &J::Arr(cand.recovered.clone()),
                            )
                            && super::json::py_eq(&base.status, &cand.status)
                            && super::json::py_eq(&cand.status, &complete),
                    ),
                ),
            ]),
        ));
    }
    let flag = |p: &J, k: &str| p.get(k).is_some_and(J::truthy);
    let held = !pairs.is_empty()
        && pairs.iter().all(|(_, p)| {
            flag(p, "kernels_identical") && flag(p, "total_within_tolerance") && flag(p, "same_log")
        });
    Ok(J::Obj(vec![
        ("pairs".into(), J::Obj(pairs)),
        ("control_held".into(), J::Bool(held)),
        (
            "rule".into(),
            J::Str(format!(
                "every arithmetic kernel ({}) identical to the instruction, the total within {} of 1, \
                 the same logarithm (PROTOCOL.md)",
                KERNEL_MARKS.join(", "),
                super::json::py_float(TOTAL_TOLERANCE)
            )),
        ),
    ]))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn percentages_become_spaces_and_objects_are_cut_as_the_script_cut_them() {
        assert_eq!(
            strip_percentages("1,028 (26.93%) 16 ( 8.64%) (x)"),
            "1,028   16   (x)"
        );
        assert_eq!(strip_object("???:f [/bin/ic]"), "???:f");
        assert_eq!(strip_object("f<[u8]> [/bin/ic]"), "f<");
        assert_eq!(strip_object("f [unclosed"), "f [unclosed");
    }

    #[test]
    fn a_table_row_reads_its_counts_and_its_function() {
        let text = "Events shown:     Ir Dr\n\
                    Ir          Dr        file:function\n\
                    --------------------------------\n\
                    2,000 (90.00%)  .   ???:a::b [/bin/ic]\n\
                    100 ( 4.50%)   7 (1.0%)  ./x.S:memset [/lib/libc.so.6]\n\
                      ./not/a/row.c\n";
        let rows = parse_annotate(text, 10).unwrap();
        assert_eq!(rows.len(), 2);
        assert_eq!(
            super::super::json::dumps_line(&rows[0], false),
            r#"{"function": "a::b", "Ir": 2000, "Dr": 0}"#
        );
        assert_eq!(
            super::super::json::dumps_line(&rows[1], false),
            r#"{"function": "./x.S:memset", "Ir": 100, "Dr": 7}"#
        );
        assert_eq!(parse_annotate(text, 1).unwrap().len(), 1);
    }
}

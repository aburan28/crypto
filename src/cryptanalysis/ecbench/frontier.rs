//! The frontier: the bounds of a domain nobody has beaten.
//!
//! Built from committed bound records, never written by hand.  Bounds are
//! grouped by domain and compared on the axes a bound carries (`ops`, the
//! ratio to the generic floor; `memory`, table entries per `√r`; and, when
//! asked, `uncharged`, the work the unit counts but does not price).  A
//! bound is *dominated* when another bound of the domain is at least as
//! good on every axis both know and clearly better on one, where clearly
//! better means the intervals do not overlap: two bounds measured on
//! different sessions cannot be paired, so the frontier is conservative
//! and leaves ties standing.  The paired answer is a challenge
//! ([`challenge`](super::challenge)).
//!
//! Unknown is not zero: an axis one side did not report is left out of the
//! comparison and the entry says so.  Inadmissible bounds are listed and
//! never on the frontier.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::cryptanalysis::ecbench::bounds::{seal_document, Bound, Domain, FIELD_AXES};

pub const FRONTIER_SCHEMA: &str = "ecbench.frontier/v1";
pub const FRONTIER_PREFIX: &str = "ECFR1h";
/// The axes dominance reads by default.  `uncharged` is reported on every
/// entry and used as an axis only when asked.
pub const DEFAULT_AXES: &[&str] = &["ops", "memory"];
/// The axes every entry reports whether or not they decide.
const REPORTED_AXES: &[&str] = &["ops", "memory", "uncharged"];

/// Every axis name a frontier or a challenge may read: the three every
/// bound carries and the three field-operation axes a bound carries when
/// its runs counted them.
pub fn is_axis(name: &str) -> bool {
    REPORTED_AXES.contains(&name) || FIELD_AXES.contains(&name)
}

/// The axis vocabulary, for error messages.
pub fn axis_names() -> String {
    REPORTED_AXES
        .iter()
        .chain(FIELD_AXES.iter())
        .copied()
        .collect::<Vec<_>>()
        .join(", ")
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Axis {
    pub known: bool,
    pub value: Option<f64>,
    pub ci95: Option<(f64, f64)>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Entry {
    pub bound_id: String,
    pub label: String,
    pub method_id: String,
    pub method: String,
    pub params: BTreeMap<String, String>,
    pub level: String,
    pub log2_r: Vec<f64>,
    pub runs: u64,
    pub axes: BTreeMap<String, Axis>,
    pub alpha: Option<f64>,
    pub alpha_ci95: Option<(f64, f64)>,
    pub declared_alpha: Option<f64>,
    pub scaling_claim: bool,
    pub bounded: bool,
    pub is_frontier: bool,
    pub dominated_by: Vec<String>,
    /// Why, per dominating bound; and which axes were left out as unknown.
    pub dominance_notes: Vec<String>,
    pub improves_on: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DomainFrontier {
    pub domain: Domain,
    pub domain_id: String,
    pub axes: Vec<String>,
    pub floor_model: String,
    pub entries: Vec<Entry>,
    /// The frontier entry with the least ops value, and the least memory.
    pub ops_leader: Option<String>,
    pub memory_leader: Option<String>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Rejected {
    pub bound_id: String,
    pub label: String,
    pub reasons: Vec<String>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Frontier {
    pub schema: String,
    pub frontier_id: String,
    pub axes: Vec<String>,
    /// Every bound read, admissible or not, sorted.
    pub built_from: Vec<String>,
    pub domains: Vec<DomainFrontier>,
    pub inadmissible: Vec<Rejected>,
}

impl Frontier {
    pub fn seal(&mut self) -> Result<String, String> {
        self.frontier_id = String::new();
        let (id, text) = seal_document(self, "frontier_id", FRONTIER_PREFIX)?;
        self.frontier_id = id;
        Ok(text)
    }
}

/// Every `*.json` bound under `paths` (files or directories), in path
/// order.  A file that is not a sealed bound is an error, not skipped: a
/// frontier built from a directory with a broken record is wrong, not
/// smaller.
pub fn load_bounds(paths: &[PathBuf]) -> Result<Vec<Bound>, String> {
    let mut files: Vec<PathBuf> = Vec::new();
    for p in paths {
        if p.is_dir() {
            let mut inner: Vec<PathBuf> = std::fs::read_dir(p)
                .map_err(|e| format!("{}: {e}", p.display()))?
                .filter_map(|e| e.ok().map(|e| e.path()))
                .filter(|q| q.extension().and_then(|x| x.to_str()) == Some("json"))
                .collect();
            inner.sort();
            files.extend(inner);
        } else {
            files.push(p.clone());
        }
    }
    let mut out = Vec::new();
    let mut seen: BTreeSet<String> = BTreeSet::new();
    for f in files {
        let b = Bound::read(&f)?;
        if !seen.insert(b.bound_id.clone()) {
            return Err(format!(
                "{}: bound {} appears twice",
                f.display(),
                b.bound_id
            ));
        }
        out.push(b);
    }
    Ok(out)
}

fn axis_of(b: &Bound, name: &str) -> Axis {
    match b.dimensions.get(name) {
        Some(d) if d.known => Axis {
            known: true,
            value: d.value,
            ci95: d.ci95,
        },
        _ => Axis {
            known: false,
            value: None,
            ci95: None,
        },
    }
}

fn entry_of(b: &Bound, axes: &[String]) -> Entry {
    let mut ax = BTreeMap::new();
    // A field-operation axis rides on an entry only when its bound
    // carries the dimension — the bound's runs counted it — whether or
    // not dominance was asked to read it; an entry without the key says
    // nothing, and `dominates` leaves the axis out as unknown.  Bounds
    // fitted before the axes existed therefore give the entries, and the
    // page, they always gave.
    for name in axes
        .iter()
        .map(String::as_str)
        .chain(REPORTED_AXES.iter().copied())
        .chain(FIELD_AXES.iter().copied())
    {
        if FIELD_AXES.contains(&name) && !b.dimensions.contains_key(name) {
            continue;
        }
        ax.entry(name.to_string())
            .or_insert_with(|| axis_of(b, name));
    }
    Entry {
        bound_id: b.bound_id.clone(),
        label: b.label.clone(),
        method_id: b.method.method_id.clone(),
        method: b.method.id.clone(),
        params: b.method.params.clone(),
        level: b.level.clone(),
        log2_r: b.sizes.iter().map(|s| s.log2_r).collect(),
        runs: b.provenance.verified,
        axes: ax,
        alpha: b.fit.alpha,
        alpha_ci95: b.fit.alpha_ci95,
        declared_alpha: b.fit.declared_alpha,
        scaling_claim: b.fit.scaling_claim,
        bounded: b.admissibility.bounded,
        is_frontier: false,
        dominated_by: vec![],
        dominance_notes: vec![],
        improves_on: b.improves_on.clone(),
    }
}

/// How `a` compares with `b` on one axis, lower being better.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Order {
    /// `a`'s interval lies wholly below `b`'s.
    ClearlyBetter,
    ClearlyWorse,
    /// The intervals overlap (or only points are known and they differ
    /// within nothing an interval could say).
    Indistinguishable,
    Unknown,
}

/// Compare two axes conservatively: with intervals on both sides, clearly
/// better needs disjoint intervals; with a point on either side, only the
/// points can be read, and a strict inequality is taken as clear.
pub fn order(a: &Axis, b: &Axis) -> Order {
    if !a.known || !b.known {
        return Order::Unknown;
    }
    let (Some(va), Some(vb)) = (a.value, b.value) else {
        return Order::Unknown;
    };
    match (a.ci95, b.ci95) {
        (Some((alo, ahi)), Some((blo, bhi))) => {
            if ahi < blo {
                Order::ClearlyBetter
            } else if alo > bhi {
                Order::ClearlyWorse
            } else {
                Order::Indistinguishable
            }
        }
        _ => {
            if va < vb {
                Order::ClearlyBetter
            } else if va > vb {
                Order::ClearlyWorse
            } else {
                Order::Indistinguishable
            }
        }
    }
}

/// `a` dominates `b` on `axes`: on every axis both know, `a`'s point
/// estimate is no larger and `a` is not clearly worse; on at least one, `a`
/// is clearly better.  The point condition keeps a wide interval from
/// knocking a precise, better entry off the frontier through another
/// axis.  Returns the axes it wins on and the axes left out as unknown.
pub fn dominates(a: &Entry, b: &Entry, axes: &[String]) -> (bool, Vec<String>, Vec<String>) {
    let mut wins = Vec::new();
    let mut unknown = Vec::new();
    for name in axes {
        let (Some(x), Some(y)) = (a.axes.get(name), b.axes.get(name)) else {
            unknown.push(name.clone());
            continue;
        };
        match order(x, y) {
            Order::ClearlyBetter => wins.push(name.clone()),
            Order::ClearlyWorse => return (false, vec![], unknown),
            Order::Indistinguishable => {
                if let (Some(vx), Some(vy)) = (x.value, y.value) {
                    if vx > vy {
                        return (false, vec![], unknown);
                    }
                }
            }
            Order::Unknown => unknown.push(name.clone()),
        }
    }
    (!wins.is_empty(), wins, unknown)
}

/// Build the frontier of every domain in `bounds`.
pub fn build(bounds: &[Bound], axes: &[String]) -> Result<Frontier, String> {
    for a in axes {
        if !is_axis(a) {
            return Err(format!("unknown axis `{a}`; axes are {}", axis_names()));
        }
    }
    let mut by_domain: BTreeMap<String, (Domain, Vec<&Bound>)> = BTreeMap::new();
    let mut inadmissible = Vec::new();
    let mut built_from: Vec<String> = bounds.iter().map(|b| b.bound_id.clone()).collect();
    built_from.sort();
    for b in bounds {
        if b.domain.id()? != b.domain_id {
            return Err(format!(
                "bound {} carries domain id {} but its domain hashes to {}",
                b.bound_id,
                b.domain_id,
                b.domain.id()?
            ));
        }
        if !b.admissible() {
            inadmissible.push(Rejected {
                bound_id: b.bound_id.clone(),
                label: b.label.clone(),
                reasons: b.admissibility.reasons.clone(),
            });
            continue;
        }
        by_domain
            .entry(b.domain_id.clone())
            .or_insert_with(|| (b.domain.clone(), Vec::new()))
            .1
            .push(b);
    }
    let mut domains = Vec::new();
    for (domain_id, (domain, bs)) in by_domain {
        let mut entries: Vec<Entry> = bs.iter().map(|b| entry_of(b, axes)).collect();
        let snapshot = entries.clone();
        for e in entries.iter_mut() {
            let mut all_unknown: BTreeSet<String> = BTreeSet::new();
            for other in snapshot.iter().filter(|o| o.bound_id != e.bound_id) {
                let (dom, wins, unknown) = dominates(other, e, axes);
                all_unknown.extend(unknown.iter().cloned());
                if dom {
                    e.dominated_by.push(other.bound_id.clone());
                    e.dominance_notes.push(format!(
                        "{} ({}) is clearly better on {}{}",
                        other.bound_id,
                        other.method,
                        wins.join(", "),
                        if unknown.is_empty() {
                            String::new()
                        } else {
                            format!("; {} unknown on one side and left out", unknown.join(", "))
                        }
                    ));
                }
            }
            e.is_frontier = e.dominated_by.is_empty();
            if e.is_frontier && !all_unknown.is_empty() {
                e.dominance_notes.push(format!(
                    "on the frontier with {} unknown against some entries",
                    all_unknown.into_iter().collect::<Vec<_>>().join(", ")
                ));
            }
        }
        let leader = |axis: &str| {
            entries
                .iter()
                .filter(|e| e.is_frontier)
                .filter_map(|e| {
                    e.axes
                        .get(axis)
                        .and_then(|a| a.value)
                        .map(|v| (v, &e.bound_id))
                })
                .min_by(|a, b| a.0.total_cmp(&b.0))
                .map(|(_, id)| id.clone())
        };
        let (ops_leader, memory_leader) = (leader("ops"), leader("memory"));
        // Frontier first, then by ops.
        entries.sort_by(|a, b| {
            b.is_frontier.cmp(&a.is_frontier).then_with(|| {
                let va = a
                    .axes
                    .get("ops")
                    .and_then(|x| x.value)
                    .unwrap_or(f64::INFINITY);
                let vb = b
                    .axes
                    .get("ops")
                    .and_then(|x| x.value)
                    .unwrap_or(f64::INFINITY);
                va.total_cmp(&vb)
            })
        });
        domains.push(DomainFrontier {
            domain,
            domain_id,
            axes: axes.to_vec(),
            floor_model: "S_floor = sqrt(pi / 2A), A the usable automorphisms; ops = S / S_floor"
                .into(),
            entries,
            ops_leader,
            memory_leader,
        });
    }
    let mut f = Frontier {
        schema: FRONTIER_SCHEMA.into(),
        frontier_id: String::new(),
        axes: axes.to_vec(),
        built_from,
        domains,
        inadmissible,
    };
    f.seal()?;
    Ok(f)
}

fn fmt(v: Option<f64>) -> String {
    v.map(|x| format!("{x:.3}")).unwrap_or_else(|| "–".into())
}

fn fmt_ci(ci: Option<(f64, f64)>) -> String {
    ci.map(|(l, h)| format!("[{l:.3}, {h:.3}]"))
        .unwrap_or_else(|| "–".into())
}

/// The frontier as a page: one table per domain, frontier rows first.
pub fn render_markdown(f: &Frontier, records_dir: &str) -> String {
    let mut out = String::new();
    out.push_str("# Measured bounds and the frontier\n\n");
    out.push_str(&format!(
        "Generated by `ecbench frontier build` from the bound records in `{records_dir}`; do not edit. Frontier `{}`, built from {} bound(s), dominance on `{}`. The protocol is [`README.md`](README.md).\n\n",
        f.frontier_id,
        f.built_from.len(),
        f.axes.join("`, `")
    ));
    out.push_str("A row is **on the frontier** when no other admissible bound of its domain is at least as good on every axis both report and clearly better (disjoint 95 % intervals) on one. Ties stand. `ops` is the mean ratio of `S = gae / √r` to the generic floor `√(π / 2A)`; `memory` is table entries per `√r`; `uncharged` is counted-but-unpriced work per `√r`. `α` is the fitted exponent of `gae = C · r^α` over the sizes listed, a scaling claim only with four or more sizes. Toy stays toy: nothing here speaks about cryptographic sizes.\n\n");
    // The field-operation columns appear in a domain only when some entry
    // there carries the axis, so a page built from bounds that predate
    // them is the page it was.
    let field_cols = |d: &DomainFrontier| -> Vec<&'static str> {
        FIELD_AXES
            .iter()
            .copied()
            .filter(|a| d.entries.iter().any(|e| e.axes.contains_key(*a)))
            .collect()
    };
    if f.domains.iter().any(|d| !field_cols(d).is_empty()) {
        out.push_str("`field muls`, `field sqrs` and `field invs` are the modular multiplications, squarings and inversions behind the group operations, per `√r`: the primitive level, counted today for the prime-field group law the generic walks and tables run on. A row showing `unknown` there did not count them; unknown is not zero.\n\n");
    }
    for d in &f.domains {
        out.push_str(&format!(
            "## {} · {} · {} targets · {} · tier {}\n\n",
            d.domain.problem, d.domain.family, d.domain.target_kind, d.domain.unit, d.domain.tier
        ));
        out.push_str(&format!(
            "Domain `{}`. Floor: {}. Ops leader: `{}`; memory leader: `{}`.\n\n",
            d.domain_id,
            d.floor_model,
            d.ops_leader.as_deref().unwrap_or("–"),
            d.memory_leader.as_deref().unwrap_or("–")
        ));
        let cols = field_cols(d);
        let col_title = |a: &str| format!(" {} (/√r) |", a.replace('_', " "));
        out.push_str("| frontier | method | params | ops (× floor) | 95 % | memory (/√r) | uncharged (/√r) |");
        for a in &cols {
            out.push_str(&col_title(a));
        }
        out.push_str(
            " α | 95 % | declared α | sizes (log₂ r) | runs | bounded | bound | dominated by |\n",
        );
        out.push_str("|:--|:--|:--|--:|:--|--:|--:|");
        for _ in &cols {
            out.push_str("--:|");
        }
        out.push_str("--:|:--|--:|:--|--:|:--|:--|:--|\n");
        for e in &d.entries {
            let ax = |n: &str| e.axes.get(n);
            let params = if e.params.is_empty() {
                "–".to_string()
            } else {
                e.params
                    .iter()
                    .map(|(k, v)| format!("{k}={v}"))
                    .collect::<Vec<_>>()
                    .join(", ")
            };
            let field_cells: String = cols
                .iter()
                .map(|a| {
                    format!(
                        " {} |",
                        ax(a)
                            .filter(|x| x.known)
                            .map(|x| fmt(x.value))
                            .unwrap_or_else(|| "unknown".into())
                    )
                })
                .collect();
            out.push_str(&format!(
                "| {} | `{}` | {} | {} | {} | {} | {} |{} {} | {} | {} | {} | {} | {} | `{}` | {} |\n",
                if e.is_frontier { "**yes**" } else { "no" },
                e.method,
                params,
                fmt(ax("ops").and_then(|a| a.value)),
                fmt_ci(ax("ops").and_then(|a| a.ci95)),
                ax("memory")
                    .filter(|a| a.known)
                    .map(|a| fmt(a.value))
                    .unwrap_or_else(|| "unknown".into()),
                fmt(ax("uncharged").and_then(|a| a.value)),
                field_cells,
                fmt(e.alpha),
                fmt_ci(e.alpha_ci95),
                fmt(e.declared_alpha),
                e.log2_r
                    .iter()
                    .map(|x| format!("{x:.1}"))
                    .collect::<Vec<_>>()
                    .join(", "),
                e.runs,
                if e.bounded { "yes" } else { "no" },
                e.bound_id,
                if e.dominated_by.is_empty() {
                    "–".to_string()
                } else {
                    e.dominated_by
                        .iter()
                        .map(|x| format!("`{x}`"))
                        .collect::<Vec<_>>()
                        .join(", ")
                }
            ));
        }
        out.push('\n');
        let notes: Vec<String> = d
            .entries
            .iter()
            .flat_map(|e| {
                e.dominance_notes
                    .iter()
                    .map(move |n| format!("- `{}` ({}): {n}", e.bound_id, e.method))
            })
            .collect();
        if !notes.is_empty() {
            out.push_str(&notes.join("\n"));
            out.push_str("\n\n");
        }
    }
    if !f.inadmissible.is_empty() {
        out.push_str("## Inadmissible bounds (listed, never on a frontier)\n\n");
        for r in &f.inadmissible {
            out.push_str(&format!(
                "- `{}` {}: {}\n",
                r.bound_id,
                r.label,
                r.reasons.join("; ")
            ));
        }
        out.push('\n');
    }
    out
}

/// Whether `path` holds exactly the rendered page.
pub fn check_markdown(f: &Frontier, records_dir: &str, path: &Path) -> Result<(), String> {
    let on_disk = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let fresh = render_markdown(f, records_dir);
    if on_disk != fresh {
        return Err(format!(
            "{} is stale: rebuild it with `ecbench frontier build`",
            path.display()
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn axis(v: f64, ci: Option<(f64, f64)>) -> Axis {
        Axis {
            known: true,
            value: Some(v),
            ci95: ci,
        }
    }

    fn entry(id: &str, ops: Axis, memory: Axis) -> Entry {
        let mut axes = BTreeMap::new();
        axes.insert("ops".to_string(), ops);
        axes.insert("memory".to_string(), memory);
        axes.insert(
            "uncharged".to_string(),
            Axis {
                known: true,
                value: Some(0.0),
                ci95: None,
            },
        );
        Entry {
            bound_id: id.into(),
            label: id.into(),
            method_id: id.into(),
            method: id.into(),
            params: BTreeMap::new(),
            level: "constant".into(),
            log2_r: vec![16.0],
            runs: 8,
            axes,
            alpha: None,
            alpha_ci95: None,
            declared_alpha: None,
            scaling_claim: false,
            bounded: false,
            is_frontier: false,
            dominated_by: vec![],
            dominance_notes: vec![],
            improves_on: vec![],
        }
    }

    #[test]
    fn disjoint_intervals_decide_and_overlapping_ones_tie() {
        let a = axis(1.0, Some((0.9, 1.1)));
        let b = axis(1.5, Some((1.3, 1.7)));
        let c = axis(1.15, Some((1.0, 1.3)));
        assert_eq!(order(&a, &b), Order::ClearlyBetter);
        assert_eq!(order(&b, &a), Order::ClearlyWorse);
        assert_eq!(order(&a, &c), Order::Indistinguishable);
        let unknown = Axis {
            known: false,
            value: None,
            ci95: None,
        };
        assert_eq!(order(&a, &unknown), Order::Unknown);
        // Points only: strict inequality.
        assert_eq!(
            order(&axis(1.0, None), &axis(2.0, None)),
            Order::ClearlyBetter
        );
    }

    #[test]
    fn dominance_needs_a_clear_win_and_no_clear_loss() {
        let axes: Vec<String> = DEFAULT_AXES.iter().map(|s| s.to_string()).collect();
        let rho = entry(
            "rho",
            axis(1.0, Some((0.9, 1.1))),
            axis(0.08, Some((0.07, 0.09))),
        );
        let bsgs = entry(
            "bsgs",
            axis(1.7, Some((1.6, 1.8))),
            axis(1.0, Some((0.98, 1.02))),
        );
        let (d, wins, _) = dominates(&rho, &bsgs, &axes);
        assert!(d);
        assert_eq!(wins, vec!["ops".to_string(), "memory".to_string()]);
        assert!(!dominates(&bsgs, &rho, &axes).0);
        // Better on ops, worse on memory: a trade, not dominance.
        let table = entry(
            "table",
            axis(0.5, Some((0.45, 0.55))),
            axis(3.0, Some((2.9, 3.1))),
        );
        assert!(!dominates(&table, &rho, &axes).0);
        assert!(!dominates(&rho, &table, &axes).0);
        // A wide interval with the worse point estimate on ops cannot
        // dominate through memory.
        let wide = entry(
            "wide",
            axis(1.6, Some((0.5, 3.0))),
            axis(0.01, Some((0.0, 0.02))),
        );
        assert!(!dominates(&wide, &rho, &axes).0);
        assert!(!dominates(&rho, &wide, &axes).0);
        // Memory unknown on one side: ops decides, and the axis is reported.
        let mut ic = entry("ic", axis(5.0, Some((4.0, 6.0))), axis(0.0, None));
        ic.axes.get_mut("memory").unwrap().known = false;
        let (d, wins, unknown) = dominates(&rho, &ic, &axes);
        assert!(d);
        assert_eq!(wins, vec!["ops".to_string()]);
        assert_eq!(unknown, vec!["memory".to_string()]);
    }
}

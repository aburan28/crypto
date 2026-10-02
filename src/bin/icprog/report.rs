//! A round's tables for its README, rendered from its `analysis.json`
//! only, so that every figure the README quotes is the analysis's own
//! (AGENTS.md §7: the page cites and never computes).

use super::json::J;

fn ci(v: Option<&J>) -> String {
    let f = |k: &str| v.and_then(|c| c.get(k)).and_then(J::as_f64);
    match (f("geomean"), f("lo"), f("hi")) {
        (Some(g), Some(lo), Some(hi)) => format!("{g:.3} [{lo:.3}, {hi:.3}]"),
        (Some(g), _, _) => format!("{g:.3}"),
        _ => "—".into(),
    }
}

fn pair(a: Option<f64>, b: Option<f64>) -> String {
    match (a, b) {
        (Some(a), Some(b)) => format!("{a:.2} → {b:.2}"),
        _ => "—".into(),
    }
}

fn num(row: &J, key: &str) -> Option<f64> {
    row.get(key).and_then(J::as_f64)
}

fn int(row: &J, key: &str) -> String {
    row.get(key)
        .and_then(J::as_i128)
        .map_or("—".into(), |v| v.to_string())
}

fn slug(row: &J) -> String {
    format!("`{}`", row.get("slug").and_then(J::as_str).unwrap_or("?"))
}

fn units(row: &J) -> String {
    let u = row.get("collect_units_per_summand");
    pair(
        u.and_then(|u| u.get("base")).and_then(J::as_f64),
        u.and_then(|u| u.get("cand")).and_then(J::as_f64),
    )
}

/// The suite and holdout tables of a round whose analysis has R05's shape.
pub fn r05(doc: &J) -> Result<String, String> {
    let mut out = String::new();
    out.push_str(
        "| curve | log₂ r | rows | rounds | cold [95%] | online | collection stage | rho online | `S` cold, base → candidate | collection, units a summand | R01's A/A, cold |\n\
         |:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|\n",
    );
    for row in doc.at("suite")?.as_arr().ok_or("`suite` is not a list")? {
        out.push_str(&format!(
            "| {} | {:.1} | {} | {} | {} | {} | {} | {} | {} | {} | {} |\n",
            slug(row),
            num(row, "log2_r").unwrap_or(f64::NAN),
            int(row, "rows"),
            int(row, "rounds"),
            ci(row.get("cold")),
            ci(row.get("online")),
            ci(row.get("collect_stage")),
            ci(row.get("rho_online")),
            pair(num(row, "s_cold_base"), num(row, "s_cold_cand")),
            units(row),
            ci(row.get("aa_cold")),
        ));
    }
    out.push('\n');
    out.push_str(
        "| curve | log₂ r | rows | rounds | pairs | cold [95%] | online | collection stage | rho online | `S` cold, base → candidate |\n\
         |:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|\n",
    );
    for row in doc
        .at("holdouts")?
        .as_arr()
        .ok_or("`holdouts` is not a list")?
    {
        let pairs = row
            .get("cold")
            .and_then(|c| c.get("n"))
            .and_then(J::as_i128)
            .map_or("—".into(), |v| v.to_string());
        out.push_str(&format!(
            "| {} | {:.1} | {} | {} | {} | {} | {} | {} | {} | {} |\n",
            slug(row),
            num(row, "log2_r").unwrap_or(f64::NAN),
            int(row, "rows"),
            int(row, "rounds"),
            pairs,
            ci(row.get("cold")),
            ci(row.get("online")),
            ci(row.get("collect_stage")),
            ci(row.get("rho_online")),
            pair(num(row, "s_cold_base"), num(row, "s_cold_cand")),
        ));
    }
    Ok(out)
}

fn rounded(x: f64, places: usize) -> J {
    J::Float(
        format!("{x:.places$}")
            .parse()
            .expect("a formatted float parses"),
    )
}

/// A host manifest's identifying fields, as the ledger carries them.
fn host_fields(host: &J) -> J {
    J::Obj(
        [
            "cpu_model",
            "logical_cores",
            "memory",
            "os",
            "rustc",
            "transparent_hugepage",
        ]
        .iter()
        .filter_map(|k| host.get(k).map(|v| (k.to_string(), v.clone())))
        .collect(),
    )
}

/// A new baseline's ledger entry (IC_TOOL_PROGRAM.md §7), from the round's
/// analysis and run tree, the ledger so far, and `meta`: the entry's
/// identity (`baseline`, `round`, `class`, `commit`, `built_from`, `base`,
/// `src_tree`, `binary_sha256`, `source`, `what_changed`, `measured`) and
/// `previous`, the baseline whose sizes carry each curve's EC1 identity.
pub fn baseline_entry(
    analysis: &J,
    ledger: &J,
    host: &J,
    host_resumed: Option<&J>,
    meta: &J,
) -> Result<J, String> {
    let baselines = ledger
        .at("baselines")?
        .as_arr()
        .ok_or("`baselines` is not a list")?;
    let find = |name: &str| {
        baselines
            .iter()
            .find(|b| b.get("baseline").and_then(J::as_str) == Some(name))
            .ok_or_else(|| format!("the ledger has no {name}"))
    };
    let v0 = find("v0")?;
    let previous = find(meta.at("previous")?.as_str().ok_or("`previous`")?)?;
    let size_of = |b: &J, alias: &str| -> Option<J> {
        b.get("sizes")?
            .as_arr()?
            .iter()
            .find(|s| s.get("curve_alias").and_then(J::as_str) == Some(alias))
            .cloned()
    };
    let mut sizes = Vec::new();
    for row in analysis
        .at("suite")?
        .as_arr()
        .ok_or("`suite` is not a list")?
    {
        let alias = row.at("slug")?.as_str().ok_or("`slug`")?;
        let prev = size_of(previous, alias)
            .ok_or_else(|| format!("{alias} is not in the previous baseline"))?;
        let unit = size_of(v0, alias)
            .and_then(|s| s.get("unit_ns").cloned())
            .ok_or_else(|| format!("{alias} has no v0 unit"))?;
        let cold = row.at("cold")?;
        let f = |k: &str| -> Result<f64, String> {
            cold.at(k)?.as_f64().ok_or_else(|| format!("cold.{k}"))
        };
        let s = |k: &str| -> Result<f64, String> {
            row.at(k)?
                .as_f64()
                .ok_or_else(|| format!("{k} is not a number"))
        };
        sizes.push(J::Obj(vec![
            ("curve_alias".into(), J::Str(alias.into())),
            ("ec1".into(), prev.at("ec1")?.clone()),
            ("curve_uid".into(), prev.at("curve_uid")?.clone()),
            ("log2_r".into(), row.at("log2_r")?.clone()),
            ("unit_ns".into(), unit),
            ("rows_measured".into(), row.at("rows")?.clone()),
            ("rounds".into(), row.at("rounds")?.clone()),
            ("s_cold".into(), rounded(s("s_cold_cand")?, 3)),
            (
                "s_cold_base_same_runs".into(),
                rounded(s("s_cold_base")?, 3),
            ),
            (
                "cold_ratio_over_base".into(),
                J::Obj(vec![
                    ("geomean".into(), rounded(f("geomean")?, 4)),
                    ("lo".into(), rounded(f("lo")?, 4)),
                    ("hi".into(), rounded(f("hi")?, 4)),
                    ("n".into(), cold.at("n")?.clone()),
                ]),
            ),
        ]));
    }
    let mut kv = Vec::new();
    for k in [
        "baseline",
        "round",
        "class",
        "commit",
        "built_from",
        "base",
        "src_tree",
        "binary_sha256",
    ] {
        kv.push((k.to_string(), meta.at(k)?.clone()));
    }
    kv.push(("host".into(), host_fields(host)));
    if let Some(h) = host_resumed {
        kv.push(("host_resumed".into(), host_fields(h)));
    }
    for k in ["source", "what_changed", "measured"] {
        kv.push((k.to_string(), meta.at(k)?.clone()));
    }
    kv.push(("sizes".into(), J::Arr(sizes)));
    Ok(J::Obj(kv))
}

#[cfg(test)]
mod tests {
    use super::super::json;
    use super::*;

    #[test]
    fn a_row_renders_its_interval_and_units() {
        let doc = json::parse(
            r#"{"suite": [{"slug": "s", "log2_r": 44.3, "rows": 8, "rounds": 10,
                "cold": {"n": 80, "geomean": 1.29, "lo": 1.25, "hi": 1.33},
                "s_cold_base": 8.04, "s_cold_cand": 6.2,
                "collect_units_per_summand": {"base": 2.3, "cand": 1.1}}],
                "holdouts": []}"#,
        )
        .unwrap();
        let t = r05(&doc).unwrap();
        assert!(t.contains("| `s` | 44.3 | 8 | 10 | 1.290 [1.250, 1.330] | — | — | — | 8.04 → 6.20 | 2.30 → 1.10 | — |"));
    }
}

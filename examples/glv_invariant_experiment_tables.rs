//! Print the tables of `RESEARCH_GLV_INVARIANT_FACTOR_BASES.md` §6 and §8
//! (E1–E15) from the frozen experiment files.
//!
//! ```text
//! cargo run --release --example glv_invariant_experiment_tables -- experiments/23_glv_invariant_e*.json
//! ```
//!
//! Every number in the note's §6 and §8 is printed by this binary; the
//! note cites, it does not compute.  A file's `rows[*].experiment` selects
//! its table.  The E1–E11 printers are the native replacement of the
//! retired Python script `scripts/glv_invariant_experiment_tables.py` (AGENTS.md, "no
//! Python") and reproduce its output byte for byte on the merged files,
//! Python's float formatting included.

use serde_json::Value;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::env;
use std::fmt::Write as _;
use std::fs;

static NULL: Value = Value::Null;

// ── access ─────────────────────────────────────────────────────────

/// `v[a][b]...` along a dotted path; a missing key is `null`.
fn at<'a>(v: &'a Value, path: &str) -> &'a Value {
    path.split('.')
        .fold(v, |acc, k| acc.get(k).unwrap_or(&NULL))
}

fn num(v: &Value) -> f64 {
    v.as_f64()
        .unwrap_or_else(|| panic!("expected a number, found {v}"))
}

fn fl(v: &Value, path: &str) -> f64 {
    num(at(v, path))
}

fn int(v: &Value) -> i128 {
    v.as_i64()
        .map(i128::from)
        .or_else(|| v.as_u64().map(i128::from))
        .unwrap_or_else(|| panic!("expected an integer, found {v}"))
}

fn it(v: &Value, path: &str) -> i128 {
    int(at(v, path))
}

/// Python truthiness.
fn truthy(v: &Value) -> bool {
    match v {
        Value::Null => false,
        Value::Bool(b) => *b,
        Value::Number(n) => n.as_f64().is_some_and(|x| x != 0.0),
        Value::String(s) => !s.is_empty(),
        Value::Array(a) => !a.is_empty(),
        Value::Object(o) => !o.is_empty(),
    }
}

fn bool_at(v: &Value, path: &str) -> bool {
    at(v, path)
        .as_bool()
        .unwrap_or_else(|| panic!("expected a bool at {path}"))
}

// ── Python-compatible formatting ───────────────────────────────────

/// `repr(float)`.
fn py_float(x: f64) -> String {
    let s = format!("{x:?}");
    match s.split_once('e') {
        Some((m, e)) => {
            let (sign, digits) = e.strip_prefix('-').map_or(("+", e), |d| ("-", d));
            format!("{m}e{sign}{digits:0>2}")
        }
        None => s,
    }
}

/// `str(x)` for a JSON value.
fn py_str(v: &Value) -> String {
    match v {
        Value::Null => "None".into(),
        Value::Bool(true) => "True".into(),
        Value::Bool(false) => "False".into(),
        Value::Number(n) if n.is_f64() => py_float(n.as_f64().unwrap()),
        Value::Number(n) => n.to_string(),
        Value::String(s) => s.clone(),
        other => other.to_string(),
    }
}

fn py_bool(b: bool) -> &'static str {
    if b {
        "True"
    } else {
        "False"
    }
}

fn yes(b: bool) -> &'static str {
    if b {
        "yes"
    } else {
        "NO"
    }
}

/// The script's `fmt`: `—` for `None`, two decimals for a float, `str`
/// otherwise.
fn fmt(v: &Value) -> String {
    match v {
        Value::Null => "—".into(),
        Value::Number(n) if n.is_f64() => format!("{:.2}", n.as_f64().unwrap()),
        other => py_str(other),
    }
}

fn fmt_f(x: Option<f64>) -> String {
    x.map_or("—".into(), |x| format!("{x:.2}"))
}

/// `round(x, n)`: the correctly rounded decimal, read back.
fn py_round(x: f64, n: usize) -> f64 {
    format!("{x:.n$}").parse().unwrap()
}

/// `format(x, '.3g')`.
fn g3(x: f64) -> String {
    if x == 0.0 {
        return "0".into();
    }
    let sci = format!("{x:.2e}");
    let (mant, exp) = sci.split_once('e').unwrap();
    let exp: i32 = exp.parse().unwrap();
    if (-4..3).contains(&exp) {
        let prec = (2 - exp) as usize;
        let s = format!("{x:.prec$}");
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            s
        }
    } else {
        let m = if mant.contains('.') {
            mant.trim_end_matches('0').trim_end_matches('.')
        } else {
            mant
        };
        let sign = if exp < 0 { '-' } else { '+' };
        format!("{m}e{sign}{:02}", exp.abs())
    }
}

fn mean(xs: &[f64]) -> f64 {
    xs.iter().sum::<f64>() / xs.len() as f64
}

fn min(xs: &[f64]) -> f64 {
    xs.iter().copied().fold(f64::INFINITY, f64::min)
}

fn max(xs: &[f64]) -> f64 {
    xs.iter().copied().fold(f64::NEG_INFINITY, f64::max)
}

/// `fmt(mean(xs)) (fmt(min(xs))–fmt(max(xs)))`, `—` for an empty list.
fn mean_range(xs: &[f64]) -> String {
    if xs.is_empty() {
        "— (—–—)".into()
    } else {
        format!("{:.2} ({:.2}–{:.2})", mean(xs), min(xs), max(xs))
    }
}

fn mean_or_dash(xs: &[f64]) -> String {
    if xs.is_empty() {
        "—".into()
    } else {
        format!("{:.2}", mean(xs))
    }
}

/// The sorted distinct values of `xs`, joined as Python's `str`.
fn joined_set(mut xs: Vec<f64>) -> String {
    xs.sort_by(f64::total_cmp);
    xs.dedup();
    xs.iter()
        .map(|&x| py_float(x))
        .collect::<Vec<_>>()
        .join(", ")
}

// ── sorting and grouping, as Python's tuples ───────────────────────

#[derive(Clone, Debug)]
enum K {
    S(String),
    I(i128),
    F(f64),
}

impl PartialEq for K {
    fn eq(&self, o: &Self) -> bool {
        self.cmp(o) == Ordering::Equal
    }
}
impl Eq for K {}
impl PartialOrd for K {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for K {
    fn cmp(&self, o: &Self) -> Ordering {
        match (self, o) {
            (K::S(a), K::S(b)) => a.cmp(b),
            (K::I(a), K::I(b)) => a.cmp(b),
            (K::F(a), K::F(b)) => a.total_cmp(b),
            (K::I(a), K::F(b)) => (*a as f64).total_cmp(b),
            (K::F(a), K::I(b)) => a.total_cmp(&(*b as f64)),
            _ => panic!("mixed key types"),
        }
    }
}

fn key(r: &Value, fields: &[&str]) -> Vec<K> {
    fields
        .iter()
        .map(|f| match at(r, f) {
            Value::String(s) => K::S(s.clone()),
            Value::Number(n) if n.is_f64() => K::F(n.as_f64().unwrap()),
            v @ Value::Number(_) => K::I(int(v)),
            v => panic!("unsortable key {f}: {v}"),
        })
        .collect()
}

fn sorted<'a>(rows: &[&'a Value], fields: &[&str]) -> Vec<&'a Value> {
    let mut v = rows.to_vec();
    v.sort_by_cached_key(|r| key(r, fields));
    v
}

fn grouped<'a>(rows: &[&'a Value], fields: &[&str]) -> BTreeMap<Vec<K>, Vec<&'a Value>> {
    let mut g: BTreeMap<Vec<K>, Vec<&Value>> = BTreeMap::new();
    for r in rows {
        g.entry(key(r, fields)).or_default().push(r);
    }
    g
}

fn ks(k: &K) -> String {
    match k {
        K::S(s) => s.clone(),
        K::I(i) => i.to_string(),
        K::F(f) => py_float(*f),
    }
}

/// The values at `path` that are not `null`, as floats.
fn present(rows: &[&Value], path: &str) -> Vec<f64> {
    rows.iter()
        .map(|r| at(r, path))
        .filter(|v| !v.is_null())
        .map(num)
        .collect()
}

// ── E1–E7 ──────────────────────────────────────────────────────────

fn both_verified(s: &Value) -> bool {
    bool_at(s, "folded.verified") && bool_at(s, "control.verified")
}

fn e1(o: &mut String, rows: &[&Value]) {
    o.push_str("### E1 — automorphism fold at the square point, one stream, both arms\n\n");
    o.push_str("| family | log2 r | h | oracle | seed | cols fold | cols control | deficiency fold / control | full-rank rel fold | full-rank rel control | ratio | full-rank rel / cols fold | full-rank rel / cols control | rank fraction at k = cols, fold / control | first-pin rel fold / control | trials | hit rate | correct |\n");
    o.push_str("|:--|--:|--:|:--|--:|--:|--:|:--|--:|--:|--:|--:|--:|:--|:--|--:|--:|:--|\n");
    for r in sorted(rows, &["family", "oracle", "log2_r", "seed"]) {
        let s = at(r, "stream");
        let (f, c) = (at(s, "folded"), at(s, "control"));
        let _ = writeln!(
            o,
            "| {} | {:.1} | {} | {} | {} | {} | {} | {} / {} | {} | {} | {} | {} | {} | {} / {} | {} / {} | {} | {:.4} | {} |",
            py_str(at(r, "family")),
            fl(r, "log2_r"),
            py_str(at(r, "cofactor")),
            py_str(at(r, "oracle")),
            py_str(at(r, "seed")),
            py_str(at(f, "columns")),
            py_str(at(c, "columns")),
            py_str(at(f, "deficiency_total")),
            py_str(at(c, "deficiency_total")),
            fmt(at(f, "square_relations")),
            fmt(at(c, "square_relations")),
            fmt(at(s, "square_ratio")),
            fmt(at(f, "square_over_columns")),
            fmt(at(c, "square_over_columns")),
            fmt(at(f, "rank_fraction_at_columns")),
            fmt(at(c, "rank_fraction_at_columns")),
            fmt(at(f, "first_pin_relations")),
            fmt(at(c, "first_pin_relations")),
            py_str(at(s, "trials")),
            fl(s, "hit_rate"),
            yes(both_verified(s)),
        );
    }
    o.push_str("\n#### E1 summary by family, oracle and size (means over seeds)\n\n");
    o.push_str("| family | oracle | log2 r (mean) | seeds | column ratio | full-rank ratio mean (min–max) | full-rank rel / cols, fold | full-rank rel / cols, control | rank fraction at k = cols, fold | control | first-pin ratio mean | all correct |\n");
    o.push_str("|:--|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|\n");
    for (k, rs) in grouped(rows, &["family", "oracle", "bits"]) {
        let sq = present(&rs, "stream.square_ratio");
        let fp = present(&rs, "stream.first_pin_ratio");
        let fo = present(&rs, "stream.folded.square_over_columns");
        let co = present(&rs, "stream.control.square_over_columns");
        let ok = rs.iter().all(|r| both_verified(at(r, "stream")));
        let col = rs
            .iter()
            .map(|r| py_round(fl(r, "stream.column_ratio"), 2))
            .collect();
        let rf = present(&rs, "stream.folded.rank_fraction_at_columns");
        let rc = present(&rs, "stream.control.rank_fraction_at_columns");
        let lr: Vec<f64> = rs.iter().map(|r| fl(r, "log2_r")).collect();
        let _ = writeln!(
            o,
            "| {} | {} | {:.1} | {} | {} | {} | {} | {} | {} | {} | {} | {} |",
            ks(&k[0]),
            ks(&k[1]),
            mean(&lr),
            rs.len(),
            joined_set(col),
            mean_range(&sq),
            mean_or_dash(&fo),
            mean_or_dash(&co),
            mean_or_dash(&rf),
            mean_or_dash(&rc),
            mean_or_dash(&fp),
            yes(ok),
        );
    }
}

fn e2(o: &mut String, rows: &[&Value]) {
    o.push_str(
        "### E2 — GLS line (`ψ`, order 4) against Koblitz orbit (`τ`, order n), one driver\n\n",
    );
    o.push_str("| family | size | log2 r | h | eigenvalue order | oracle | seed | points | cols fold | cols control | pts/col | column ratio | deficiency fold / control | square rel fold | square rel control | ratio | square/cols fold | trials | hit rate | oracle F_p muls / call | agreement with subtract | correct |\n");
    o.push_str("|:--|:--|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|:--|:--|\n");
    let per_call = |r: &Value| {
        let (m, c) = (at(r, "oracle_fp_muls"), at(r, "oracle_calls"));
        (truthy(m) && truthy(c)).then(|| num(m) / num(c))
    };
    let agr_str =
        |agree: i128, targets: i128, dis: i128| format!("{agree}/{targets} ({dis} disagree)");
    for r in sorted(rows, &["family", "log2_r", "seed"]) {
        let s = at(r, "stream");
        let (f, c) = (at(s, "folded"), at(s, "control"));
        let size = if r.get("p_bits").is_some() {
            format!("p=2^{}", py_str(at(r, "p_bits")))
        } else {
            format!("n={}", py_str(at(r, "n")))
        };
        let per = per_call(r).map_or("—".into(), |x| format!("{x:.0}"));
        let a = at(r, "agreement_with_subtract");
        let agr = if truthy(a) {
            agr_str(it(a, "agree"), it(a, "targets"), it(a, "disagree"))
        } else {
            "—".into()
        };
        let _ = writeln!(
            o,
            "| {} | {} | {:.1} | {} | {} | {} | {} | {} | {} | {} | {:.1} | {:.2} | {} / {} | {} | {} | {} | {} | {} | {:.4} | {} | {} | {} |",
            py_str(at(r, "family")),
            size,
            fl(r, "log2_r"),
            py_str(at(r, "cofactor")),
            py_str(at(r, "eigenvalue_order")),
            py_str(at(r, "oracle")),
            py_str(at(r, "seed")),
            py_str(at(f, "points")),
            py_str(at(f, "columns")),
            py_str(at(c, "columns")),
            fl(f, "points_per_column"),
            fl(s, "column_ratio"),
            py_str(at(f, "deficiency_total")),
            py_str(at(c, "deficiency_total")),
            fmt(at(f, "square_relations")),
            fmt(at(c, "square_relations")),
            fmt(at(s, "square_ratio")),
            fmt(at(f, "square_over_columns")),
            py_str(at(s, "trials")),
            fl(s, "hit_rate"),
            per,
            agr,
            yes(both_verified(s)),
        );
    }
    o.push_str("\n#### E2 summary\n\n");
    o.push_str("| family | eigenvalue order | instances | log2 r | points per column | column ratio | full-rank ratio mean (min–max) | rank fraction at k = cols, fold / control | oracle F_p muls / call | agreement with subtract | all correct |\n");
    o.push_str("|:--|--:|--:|:--|--:|--:|--:|--:|--:|:--|:--|\n");
    for (k, rs) in grouped(rows, &["family", "eigenvalue_order"]) {
        let sq = present(&rs, "stream.square_ratio");
        let ok = rs.iter().all(|r| both_verified(at(r, "stream")));
        let ppc = rs
            .iter()
            .map(|r| py_round(fl(r, "stream.folded.points_per_column"), 1))
            .collect();
        let col = rs
            .iter()
            .map(|r| py_round(fl(r, "stream.column_ratio"), 2))
            .collect();
        let rf = present(&rs, "stream.folded.rank_fraction_at_columns");
        let rc = present(&rs, "stream.control.rank_fraction_at_columns");
        let lr: Vec<f64> = rs.iter().map(|r| fl(r, "log2_r")).collect();
        let per: Vec<f64> = rs.iter().filter_map(|r| per_call(r)).collect();
        let agr: Vec<&Value> = rs
            .iter()
            .map(|r| at(r, "agreement_with_subtract"))
            .filter(|a| truthy(a))
            .collect();
        let agr_s = if agr.is_empty() {
            "—".into()
        } else {
            agr_str(
                agr.iter().map(|a| it(a, "agree")).sum(),
                agr.iter().map(|a| it(a, "targets")).sum(),
                agr.iter().map(|a| it(a, "disagree")).sum(),
            )
        };
        let _ = writeln!(
            o,
            "| {} | {} | {} | {:.1}–{:.1} | {} | {} | {} | {} / {} | {} | {} | {} |",
            ks(&k[0]),
            ks(&k[1]),
            rs.len(),
            min(&lr),
            max(&lr),
            joined_set(ppc),
            joined_set(col),
            mean_range(&sq),
            mean_or_dash(&rf),
            mean_or_dash(&rc),
            if per.is_empty() {
                "—".into()
            } else {
                format!("{:.0}", mean(&per))
            },
            agr_s,
            yes(ok),
        );
    }
}

fn e3(o: &mut String, rows: &[&Value]) {
    o.push_str("### E3 — composite groups on twisted CM curves\n\n");
    o.push_str("| family | p | log2 r | h | ord λ_ψ | ord λ_aut | aut = ±ψ on ⟨G⟩ | pts/col negation | pts/col ψ | pts/col aut | pts/col ψ+aut | cols ψ+aut | square ratio ψ+aut vs negation | square ratio ψ vs negation | correct |\n");
    o.push_str("|:--|--:|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|:--|\n");
    for r in sorted(rows, &["family", "p_bits", "seed"]) {
        let f = at(r, "folds");
        let sb = at(r, "stream_psi_aut_vs_negation");
        let sp = at(r, "stream_psi_vs_negation");
        let ok = if sb.is_null() && sp.is_null() {
            "folds only (r < 16·columns)".to_string()
        } else {
            yes((sb.is_null() || both_verified(sb)) && (sp.is_null() || both_verified(sp)))
                .to_string()
        };
        let ratio = |s: &Value| {
            if s.is_null() {
                "—".into()
            } else {
                fmt(at(s, "square_ratio"))
            }
        };
        let _ = writeln!(
            o,
            "| {} | 2^{} | {:.1} | {} | {} | {} | {} | {:.1} | {:.1} | {:.1} | {:.1} | {} | {} | {} | {} |",
            py_str(at(r, "family")),
            py_str(at(r, "p_bits")),
            fl(r, "log2_r"),
            py_str(at(r, "cofactor")),
            py_str(at(r, "psi_eigenvalue_order")),
            py_str(at(r, "aut_eigenvalue_order")),
            py_str(at(r, "aut_is_plus_minus_psi_on_subgroup")),
            fl(f, "negation.points_per_column"),
            fl(f, "psi.points_per_column"),
            fl(f, "aut.points_per_column"),
            num(&f["psi+aut"]["points_per_column"]),
            py_str(&f["psi+aut"]["columns"]),
            ratio(sb),
            ratio(sp),
            ok,
        );
    }
}

fn e4(o: &mut String, rows: &[&Value]) {
    o.push_str("### E4 — type C: degree-2 and degree-3 CM endomorphisms against a base\n\n");
    o.push_str("| family | degree | log2 r | h | seed | maps found | map | trace | ord_r(λ) | images in base | base points | chance fraction | verified |\n");
    o.push_str("|:--|--:|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|:--|\n");
    let maps = |r: &Value| r["maps"].as_array().cloned().unwrap_or_default();
    for r in sorted(rows, &["family", "degree", "log2_r", "seed"]) {
        let ms = maps(r);
        if ms.is_empty() {
            let _ = writeln!(
                o,
                "| {} | {} | {:.1} | {} | {} | 0 | — | — | — | — | {} | — | — |",
                py_str(at(r, "family")),
                py_str(at(r, "degree")),
                fl(r, "log2_r"),
                py_str(at(r, "cofactor")),
                py_str(at(r, "seed")),
                py_str(at(r, "base_points")),
            );
        }
        for m in &ms {
            let ov = at(m, "overlap");
            let name = m["name"].as_str().unwrap();
            let trace = name
                .split("trace=")
                .nth(1)
                .map_or("—".to_string(), |t| t.trim_end_matches(']').to_string());
            let _ = writeln!(
                o,
                "| {} | {} | {:.1} | {} | {} | {} | `{}` | {} | {} | {} | {} | {:.5} | {} |",
                py_str(at(r, "family")),
                py_str(at(r, "degree")),
                fl(r, "log2_r"),
                py_str(at(r, "cofactor")),
                py_str(at(r, "seed")),
                ms.len(),
                name,
                trace,
                py_str(at(ov, "eigenvalue_order")),
                py_str(at(ov, "images_in_base")),
                py_str(at(ov, "base_points")),
                fl(ov, "chance_fraction"),
                py_str(at(m, "verified")),
            );
        }
    }
    o.push_str("\n#### E4 summary\n\n");
    o.push_str("| family | degree | instances | maps per instance | ord_r(λ) min–max | images in base, total | base points, total | all verified |\n");
    o.push_str("|:--|--:|--:|--:|--:|--:|--:|:--|\n");
    for (k, rs) in grouped(rows, &["family", "degree"]) {
        let all: Vec<Value> = rs.iter().flat_map(|r| maps(r)).collect();
        let orders: Vec<i128> = all
            .iter()
            .map(|m| it(m, "overlap.eigenvalue_order"))
            .collect();
        let inside: i128 = all.iter().map(|m| it(m, "overlap.images_in_base")).sum();
        let pts: i128 = all.iter().map(|m| it(m, "overlap.base_points")).sum();
        let ok =
            all.iter().all(|m| bool_at(m, "verified")) && rs.iter().all(|r| !maps(r).is_empty());
        let mut counts: Vec<usize> = rs.iter().map(|r| maps(r).len()).collect();
        counts.sort();
        counts.dedup();
        let (lo, hi) = match (orders.iter().min(), orders.iter().max()) {
            (Some(a), Some(b)) => (a.to_string(), b.to_string()),
            _ => ("—".into(), "—".into()),
        };
        let _ = writeln!(
            o,
            "| {} | {} | {} | {} | {}–{} | {} | {} | {} |",
            ks(&k[0]),
            ks(&k[1]),
            rs.len(),
            counts
                .iter()
                .map(|c| c.to_string())
                .collect::<Vec<_>>()
                .join(", "),
            lo,
            hi,
            inside,
            pts,
            yes(ok),
        );
    }
}

fn e5(o: &mut String, rows: &[&Value]) {
    o.push_str("### E5 — subfield curves on E(F_{p³}): the Frobenius line\n\n");
    o.push_str("| family | p | log2 r | #E(F_p) | h / #E(F_p) | ord λ_π | ord λ_ζ | ζ ∈ ⟨π⟩ on ⟨G⟩ | pts/col negation | pts/col π | pts/col ζ | pts/col π+ζ | cols π | cols negation | full-rank rel π | full-rank rel negation | ratio | rank fraction at k = cols, π / negation | trials | hit rate | oracle muls/call | degenerate | correct |\n");
    o.push_str("|:--|--:|--:|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|--:|--:|--:|--:|:--|\n");
    for r in sorted(rows, &["family", "p_bits", "seed"]) {
        let f = at(r, "folds");
        let st = at(r, "relation_stream");
        let tail = if truthy(st) {
            let s = at(st, "stream");
            let (fo, co) = (at(s, "folded"), at(s, "control"));
            let per = if truthy(at(st, "oracle_calls")) {
                format!("{:.0}", fl(st, "oracle_fp_muls") / fl(st, "oracle_calls"))
            } else {
                "—".into()
            };
            format!(
                "{} | {} | {} | {} | {} | {} / {} | {} | {:.5} | {} | {} | {} |",
                py_str(at(fo, "columns")),
                py_str(at(co, "columns")),
                fmt(at(fo, "square_relations")),
                fmt(at(co, "square_relations")),
                fmt(at(s, "square_ratio")),
                fmt(at(fo, "rank_fraction_at_columns")),
                fmt(at(co, "rank_fraction_at_columns")),
                py_str(at(s, "trials")),
                fl(s, "hit_rate"),
                per,
                py_str(at(st, "oracle_degenerate")),
                yes(both_verified(s)),
            )
        } else {
            format!(
                "{} | {} | — | — | — | — | — | — | — | — | not an index calculus: the base is the subgroup |",
                py_str(at(f, "frobenius.columns")),
                py_str(at(f, "negation.columns")),
            )
        };
        let raw = |k: &str| {
            f.get(k)
                .map_or("—".to_string(), |x| py_str(at(x, "points_per_column")))
        };
        let _ = writeln!(
            o,
            "| {} | 2^{} | {:.1} | {} | {} | {} | {} | {} | {:.1} | {:.1} | {} | {} | {}",
            py_str(at(r, "family")),
            py_str(at(r, "p_bits")),
            fl(r, "log2_r"),
            py_str(at(r, "base_order")),
            it(r, "cofactor").div_euclid(it(r, "base_order")),
            py_str(at(r, "frobenius_eigenvalue_order")),
            fmt(at(r, "zeta_eigenvalue_order")),
            fmt(at(r, "zeta_is_frobenius_power")),
            fl(f, "negation.points_per_column"),
            fl(f, "frobenius.points_per_column"),
            raw("zeta"),
            raw("frobenius+zeta"),
            tail,
        );
    }
    o.push_str("\n#### E5 summary\n\n");
    o.push_str("| family | instances | log2 r | pts/col negation | pts/col π | pts/col ζ | pts/col π+ζ | ζ ∈ ⟨π⟩ on ⟨G⟩ | streams | full-rank ratio mean (min–max) | rank fraction at k = cols, π / negation | oracle muls/call | degenerate | all correct |\n");
    o.push_str("|:--|--:|:--|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|:--|\n");
    for (k, rs) in grouped(rows, &["family"]) {
        let ppc = |name: &str| {
            let xs: Vec<f64> = rs
                .iter()
                .filter_map(|r| r["folds"].get(name))
                .map(|x| py_round(fl(x, "points_per_column"), 1))
                .collect();
            if xs.is_empty() {
                "—".into()
            } else {
                joined_set(xs)
            }
        };
        let zf: Vec<&Value> = rs
            .iter()
            .map(|r| at(r, "zeta_is_frobenius_power"))
            .collect();
        let zf_s = if zf.iter().all(|v| v.is_null()) {
            "—"
        } else if zf.iter().all(|v| v.as_bool() == Some(true)) {
            "all"
        } else {
            "NO"
        };
        let sts: Vec<&Value> = rs
            .iter()
            .map(|r| at(r, "relation_stream"))
            .filter(|s| truthy(s))
            .collect();
        let sq = present(&sts, "stream.square_ratio");
        let rf = present(&sts, "stream.folded.rank_fraction_at_columns");
        let rc = present(&sts, "stream.control.rank_fraction_at_columns");
        let per: Vec<f64> = sts
            .iter()
            .filter(|st| truthy(at(st, "oracle_calls")))
            .map(|st| fl(st, "oracle_fp_muls") / fl(st, "oracle_calls"))
            .collect();
        let deg: i128 = sts.iter().map(|st| it(st, "oracle_degenerate")).sum();
        let ok = sts.iter().all(|st| both_verified(at(st, "stream")));
        let lr: Vec<f64> = rs.iter().map(|r| fl(r, "log2_r")).collect();
        let _ =
            writeln!(
            o,
            "| {} | {} | {:.1}–{:.1} | {} | {} | {} | {} | {} | {} | {} | {} / {} | {} | {} | {} |",
            ks(&k[0]),
            rs.len(),
            min(&lr),
            max(&lr),
            ppc("negation"),
            ppc("frobenius"),
            ppc("zeta"),
            ppc("frobenius+zeta"),
            zf_s,
            sts.len(),
            mean_range(&sq),
            mean_or_dash(&rf),
            mean_or_dash(&rc),
            if per.is_empty() {
                "—".into()
            } else {
                format!("{:.0}", mean(&per))
            },
            if sts.is_empty() {
                "—".into()
            } else {
                deg.to_string()
            },
            if sts.is_empty() { "folds only" } else { yes(ok) },
        );
    }
}

fn e6(o: &mut String, rows: &[&Value]) {
    o.push_str(
        "### E6 — the matched rho: folded by the automorphism group against the negation walk\n\n",
    );
    o.push_str("| family | log2 r | h | seed | walks | S negation (mean) | S folded (mean) | S ratio | steps ratio | expected √(A/2) | all verified |\n");
    o.push_str("|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|\n");
    for r in sorted(rows, &["family", "log2_r", "seed"]) {
        let _ = writeln!(
            o,
            "| {} | {:.1} | {} | {} | {} | {:.3} | {:.3} | {:.2} | {:.2} | {:.2} | {} |",
            py_str(at(r, "family")),
            fl(r, "log2_r"),
            py_str(at(r, "cofactor")),
            py_str(at(r, "seed")),
            py_str(at(r, "rho_runs")),
            fl(r, "negation_s_mean"),
            fl(r, "folded_s_mean"),
            fl(r, "s_ratio"),
            fl(r, "steps_ratio"),
            fl(r, "expected_ratio"),
            py_str(at(r, "all_verified")),
        );
    }
    o.push_str("\n#### E6 summary\n\n");
    o.push_str("| family | instances | walks | S ratio mean (min–max) | steps ratio mean | expected | all verified |\n");
    o.push_str("|:--|--:|--:|--:|--:|--:|:--|\n");
    for (k, rs) in grouped(rows, &["family"]) {
        let sr: Vec<f64> = rs.iter().map(|r| fl(r, "s_ratio")).collect();
        let st: Vec<f64> = rs.iter().map(|r| fl(r, "steps_ratio")).collect();
        let _ = writeln!(
            o,
            "| {} | {} | {} | {:.2} ({:.2}–{:.2}) | {:.2} | {:.2} | {} |",
            ks(&k[0]),
            rs.len(),
            rs.iter().map(|r| it(r, "rho_runs")).sum::<i128>(),
            mean(&sr),
            min(&sr),
            max(&sr),
            mean(&st),
            fl(rs[0], "expected_ratio"),
            py_bool(rs.iter().all(|r| bool_at(r, "all_verified"))),
        );
    }
}

fn e7(o: &mut String, rows: &[&Value]) {
    o.push_str("### E7 — three summands and orbit duplicates on j = 0\n\n");
    o.push_str("| log2 r | seed | m | group order | seed abscissae | points | cols fold | cols control | square rel fold | square rel control | ratio | zero-support rows fold / control | single-column rows fold / control | orbit duplicates | targets | hit rate | correct |\n");
    o.push_str("|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|--:|--:|--:|:--|\n");
    for r in sorted(rows, &["log2_r", "seed", "summands"]) {
        let s = at(r, "stream");
        let (f, c) = (at(s, "folded"), at(s, "control"));
        let _ = writeln!(
            o,
            "| {:.1} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} / {} | {} / {} | {} | {} | {:.4} | {} |",
            fl(r, "log2_r"),
            py_str(at(r, "seed")),
            py_str(at(r, "summands")),
            py_str(at(r, "group_order")),
            py_str(at(r, "seed_abscissae")),
            py_str(at(f, "points")),
            py_str(at(f, "columns")),
            py_str(at(c, "columns")),
            fmt(at(f, "square_relations")),
            fmt(at(c, "square_relations")),
            fmt(at(s, "square_ratio")),
            py_str(at(f, "zero_support_rows")),
            py_str(at(c, "zero_support_rows")),
            py_str(at(f, "single_column_rows")),
            py_str(at(c, "single_column_rows")),
            py_str(at(s, "orbit_duplicates")),
            py_str(at(s, "trials")),
            fl(s, "hit_rate"),
            yes(both_verified(s)),
        );
    }
}

// ── E8, E11: every phase priced ────────────────────────────────────

/// Least-squares slope of `log2 y` against `log2 x`.
fn fit_exponent(points: &[(f64, f64)]) -> Option<f64> {
    let pts: Vec<(f64, f64)> = points
        .iter()
        .filter(|(x, y)| *x > 0.0 && *y > 0.0)
        .map(|(x, y)| (x.log2(), y.log2()))
        .collect();
    if pts.len() < 2 {
        return None;
    }
    let mx = mean(&pts.iter().map(|p| p.0).collect::<Vec<_>>());
    let my = mean(&pts.iter().map(|p| p.1).collect::<Vec<_>>());
    let sxx: f64 = pts.iter().map(|p| (p.0 - mx).powi(2)).sum();
    if sxx == 0.0 {
        return None;
    }
    Some(pts.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum::<f64>() / sxx)
}

#[derive(Clone, Copy, Debug)]
struct Arm {
    build: f64,
    table: f64,
    stream: f64,
    la: f64,
    total: f64,
    s: f64,
}

#[derive(Clone, Copy, Debug)]
struct Walk {
    total: f64,
    s: f64,
    per_op: f64,
    verified: bool,
}

#[derive(Clone, Copy, Debug)]
struct Costs {
    folded: Arm,
    control: Arm,
    rho_negation: Walk,
    rho_folded: Walk,
    muls_per_add: f64,
}

/// What one stream costs, in `F_p` multiplications, before it is shared
/// between the arms.
struct Priced<'a> {
    /// The stream report (`folded`, `control`, `trials`, `group_ops`).
    stream: &'a Value,
    /// The stream's own field operations plus anything its oracle counts
    /// itself (a solver's multiplications, a key's uncounted ones).
    stream_total: f64,
    /// The pair table, built once and charged to both arms.
    table: f64,
    /// Base build per arm.
    build_folded: f64,
    build_control: f64,
}

fn price_arms(r: &Value, p: &Priced) -> (Arm, Arm) {
    let sqrt_r = fl(r, "r").sqrt();
    let s = p.stream;
    let arm = |name: &str, build: f64| {
        let a = at(s, name);
        let (sqt, trials) = (at(a, "square_trials"), at(s, "trials"));
        let share = if truthy(sqt) && truthy(trials) {
            num(sqt) / num(trials)
        } else {
            1.0
        };
        let stream = p.stream_total * share;
        let la = fl(a, "row_ops");
        let total = build + p.table + stream + la;
        Arm {
            build,
            table: p.table,
            stream,
            la,
            total,
            s: total / sqrt_r,
        }
    };
    (
        arm("folded", p.build_folded),
        arm("control", p.build_control),
    )
}

/// The matched rho walks of a row, in `F_p` multiplications.
fn price_walks(r: &Value, k_inv: f64) -> (Walk, Walk) {
    (
        price_walk(r, k_inv, "negation"),
        price_walk(r, k_inv, "folded"),
    )
}

/// One kind of matched walk of a row (`negation` or `folded`).
fn price_walk(r: &Value, k_inv: f64, name: &str) -> Walk {
    let sqrt_r = fl(r, "r").sqrt();
    let walks = r["walks"].as_array().unwrap();
    let tot: Vec<f64> = walks
        .iter()
        .map(|w| fl(w, &format!("{name}_fp_muls")) + k_inv * fl(w, &format!("{name}_fp_invs")))
        .collect();
    let gae: Vec<f64> = walks.iter().map(|w| fl(&w[name], "gae")).collect();
    Walk {
        total: mean(&tot),
        s: mean(&tot) / sqrt_r,
        per_op: mean(&tot) / mean(&gae).max(1.0),
        verified: walks.iter().all(|w| bool_at(&w[name], "verified")),
    }
}

fn unit(v: &Value, k_inv: f64) -> f64 {
    fl(v, "fp_muls") + k_inv * fl(v, "fp_invs")
}

/// Every phase of one E8/E11 row in `F_p` multiplications (inversions
/// priced by the row's measured conversion; a row operation of the linear
/// algebra is one `Z/rZ` multiplication, priced as one).
fn e8_costs(r: &Value) -> Costs {
    let k = fl(r, "inversion_in_multiplications");
    let s = at(r, "stream");
    let mut stream_total = unit(at(r, "stream_cost"), k);
    // The algebraic oracle counts its own multiplications (E11); they are
    // stream cost, shared by the arms in the same proportion.
    if r.get("solver").is_some() {
        stream_total += fl(r, "solver.fp_muls");
    }
    let table = if r.get("pair_table").is_some() {
        unit(at(r, "pair_table"), k)
    } else {
        0.0
    };
    let (folded, control) = price_arms(
        r,
        &Priced {
            stream: s,
            stream_total,
            table,
            build_folded: unit(at(r, "base_build.folded"), k),
            build_control: unit(at(r, "base_build.control"), k),
        },
    );
    let (rho_negation, rho_folded) = price_walks(r, k);
    Costs {
        folded,
        control,
        rho_negation,
        rho_folded,
        muls_per_add: fl(r, "stream_cost.fp_muls")
            / (fl(s, "group_ops.adds") + fl(s, "group_ops.doubles")).max(1.0),
    }
}

fn e8(o: &mut String, rows: &[&Value], label: &str, oracle: &str) {
    let _ = writeln!(o, "### {label} — three summands on the Frobenius line of a subfield curve, {oracle}: relations\n");
    o.push_str("| p | log2 r | h / #E(F_p) | cols π | cols negation | full-rank rel π | full-rank rel negation | ratio | rank fraction at k = cols, π / negation | single-column rows π / negation | orbit duplicates | targets | hit rate | correct |\n");
    o.push_str("|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|--:|--:|--:|:--|\n");
    let rows = sorted(rows, &["log2_r", "seed"]);
    for r in &rows {
        let s = at(r, "stream");
        let (f, c) = (at(s, "folded"), at(s, "control"));
        let _ = writeln!(
            o,
            "| 2^{} | {:.1} | {} | {} | {} | {} | {} | {} | {} / {} | {} / {} | {} | {} | {:.4} | {} |",
            py_str(at(r, "p_bits")),
            fl(r, "log2_r"),
            it(r, "cofactor").div_euclid(it(r, "base_order")),
            py_str(at(f, "columns")),
            py_str(at(c, "columns")),
            fmt(at(f, "square_relations")),
            fmt(at(c, "square_relations")),
            fmt(at(s, "square_ratio")),
            fmt(at(f, "rank_fraction_at_columns")),
            fmt(at(c, "rank_fraction_at_columns")),
            py_str(at(f, "single_column_rows")),
            py_str(at(c, "single_column_rows")),
            py_str(at(s, "orbit_duplicates")),
            py_str(at(s, "trials")),
            fl(s, "hit_rate"),
            yes(both_verified(s)),
        );
    }
    let _ = writeln!(o, "\n#### {label} — every phase priced, one unit (F_p multiplications; inversion = measured factor; LA row op = 1); rows that reached full rank on both arms\n");
    o.push_str("| p | log2 r | inv / mul (measured) | muls per addition | base build π / neg | pair table | stream (incl. oracle) π / neg | LA π / neg | total π | total negation | S π | S negation | S neg / S π | rho S negation | rho S folded | rho muls per group op, negation / folded | rho steps ratio (expected) | S π / rho S folded | all verified |\n");
    o.push_str("|--:|--:|--:|--:|:--|--:|:--|:--|--:|--:|--:|--:|--:|--:|--:|:--|:--|--:|:--|\n");
    let mut fits: BTreeMap<&str, Vec<(f64, f64)>> = BTreeMap::new();
    let at_full = |r: &Value| {
        truthy(at(r, "stream.folded.square_relations"))
            && truthy(at(r, "stream.control.square_relations"))
    };
    let full: Vec<&Value> = rows.iter().copied().filter(|r| at_full(r)).collect();
    for r in &rows {
        if !at_full(r) {
            let _ = writeln!(
                o,
                "| 2^{} | {:.1} | — | — | not full rank on both arms: the stream exhausted the group's targets | | | | | | | | | | | | | | |",
                py_str(at(r, "p_bits")),
                fl(r, "log2_r"),
            );
            continue;
        }
        let c = e8_costs(r);
        let (f, n) = (c.folded, c.control);
        let (rn, rf) = (c.rho_negation, c.rho_folded);
        let ok = both_verified(at(r, "stream")) && rn.verified && rf.verified;
        let rr = fl(r, "r");
        for (key, val) in [
            ("folded", f.total),
            ("control", n.total),
            ("rho_folded", rf.total),
            ("rho_negation", rn.total),
        ] {
            fits.entry(key).or_default().push((rr, val));
        }
        let _ = writeln!(
            o,
            "| 2^{} | {:.1} | {:.1} | {:.1} | {} / {} | {} | {} / {} | {} / {} | {} | {} | {:.1} | {:.1} | {:.2} | {:.1} | {:.1} | {:.0} / {:.0} | {:.2} ({:.2}) | {:.0} | {} |",
            py_str(at(r, "p_bits")),
            fl(r, "log2_r"),
            fl(r, "inversion_in_multiplications"),
            c.muls_per_add,
            g3(f.build),
            g3(n.build),
            g3(f.table),
            g3(f.stream),
            g3(n.stream),
            g3(f.la),
            g3(n.la),
            g3(f.total),
            g3(n.total),
            f.s,
            n.s,
            n.s / f.s,
            rn.s,
            rf.s,
            rn.per_op,
            rf.per_op,
            fl(r, "rho_steps_ratio"),
            fl(r, "rho_expected_ratio"),
            f.s / rf.s,
            yes(ok),
        );
    }
    let _ = writeln!(
        o,
        "\n#### {label} summary (rows at full rank on both arms)\n"
    );
    o.push_str("| arms | instances (of run) | log2 r | column ratio | full-rank ratio mean (min–max) | rank fraction at k = cols, π / negation | S neg / S π mean (min–max) | S π / rho S folded, min–max | rho steps ratio mean (expected) | rho muls per group op, negation / folded | rho S negation / rho S folded, mean | fitted exponent of total: π, negation, rho folded, rho negation (rho: 0.50) | all verified |\n");
    o.push_str("|:--|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|\n");
    if full.is_empty() {
        let _ = writeln!(o, "| π-line fold vs negation, m = 3, {oracle} | 0 ({}) | no row reached full rank on both arms | | | | | | | | | | |", rows.len());
        return;
    }
    let sq = present(&full, "stream.square_ratio");
    let rf_ = present(&full, "stream.folded.rank_fraction_at_columns");
    let rc_ = present(&full, "stream.control.rank_fraction_at_columns");
    let costs: Vec<Costs> = full.iter().map(|r| e8_costs(r)).collect();
    let s_ratio: Vec<f64> = costs.iter().map(|c| c.control.s / c.folded.s).collect();
    let vs_rho: Vec<f64> = costs.iter().map(|c| c.folded.s / c.rho_folded.s).collect();
    let rho_ratio: Vec<f64> = costs
        .iter()
        .map(|c| c.rho_negation.s / c.rho_folded.s)
        .collect();
    let ok = rows
        .iter()
        .all(|r| both_verified(at(r, "stream")) && bool_at(r, "rho_all_verified"));
    let lr: Vec<f64> = full.iter().map(|r| fl(r, "log2_r")).collect();
    let col = full
        .iter()
        .map(|r| py_round(fl(r, "stream.column_ratio"), 2))
        .collect();
    let exps = ["folded", "control", "rho_folded", "rho_negation"]
        .iter()
        .map(|k| fmt_f(fit_exponent(fits.get(k).map_or(&[][..], |v| v))))
        .collect::<Vec<_>>()
        .join(", ");
    let steps: Vec<f64> = full.iter().map(|r| fl(r, "rho_steps_ratio")).collect();
    let _ = writeln!(
        o,
        "| π-line fold vs negation, m = 3, {oracle} | {} ({}) | {:.1}–{:.1} | {} | {} | {} / {} | {:.2} ({:.2}–{:.2}) | {:.0}–{:.0} | {:.2} ({:.2}) | {:.0} / {:.0} | {:.2} | {} | {} |",
        full.len(),
        rows.len(),
        min(&lr),
        max(&lr),
        joined_set(col),
        mean_range(&sq),
        mean_or_dash(&rf_),
        mean_or_dash(&rc_),
        mean(&s_ratio),
        min(&s_ratio),
        max(&s_ratio),
        min(&vs_rho),
        max(&vs_rho),
        mean(&steps),
        fl(full[0], "rho_expected_ratio"),
        mean(&costs.iter().map(|c| c.rho_negation.per_op).collect::<Vec<_>>()),
        mean(&costs.iter().map(|c| c.rho_folded.per_op).collect::<Vec<_>>()),
        mean(&rho_ratio),
        exps,
        yes(ok),
    );
}

fn e9(o: &mut String, rows: &[&Value]) {
    o.push_str("### E9 — the matched folded rho on the F_{p²} and F_{p³} groups\n\n");
    o.push_str("| family | p | log2 r | generators | A | walks | S negation (mean) | S folded (mean) | S ratio | steps ratio | expected √(A/2) | all verified |\n");
    o.push_str("|:--|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|:--|\n");
    for r in sorted(rows, &["family", "log2_r", "seed"]) {
        let gens: Vec<String> = r["generators"]
            .as_array()
            .unwrap()
            .iter()
            .map(py_str)
            .collect();
        let _ = writeln!(
            o,
            "| {} | 2^{} | {:.1} | {} | {} | {} | {:.3} | {:.3} | {:.2} | {:.2} | {:.2} | {} |",
            py_str(at(r, "family")),
            py_str(at(r, "p_bits")),
            fl(r, "log2_r"),
            gens.join(", "),
            py_str(at(r, "group_order_of_fold")),
            py_str(at(r, "rho_runs")),
            fl(r, "negation_s_mean"),
            fl(r, "folded_s_mean"),
            fl(r, "s_ratio"),
            fl(r, "steps_ratio"),
            fl(r, "expected_ratio"),
            py_str(at(r, "all_verified")),
        );
    }
    o.push_str("\n#### E9 summary\n\n");
    o.push_str("| family | A | instances | log2 r | walks | S ratio mean (min–max) | steps ratio mean | expected | all verified |\n");
    o.push_str("|:--|--:|--:|:--|--:|--:|--:|--:|:--|\n");
    for (k, rs) in grouped(rows, &["family", "group_order_of_fold"]) {
        let sr: Vec<f64> = rs.iter().map(|r| fl(r, "s_ratio")).collect();
        let lr: Vec<f64> = rs.iter().map(|r| fl(r, "log2_r")).collect();
        let st: Vec<f64> = rs.iter().map(|r| fl(r, "steps_ratio")).collect();
        let _ = writeln!(
            o,
            "| {} | {} | {} | {:.1}–{:.1} | {} | {:.2} ({:.2}–{:.2}) | {:.2} | {:.2} | {} |",
            ks(&k[0]),
            ks(&k[1]),
            rs.len(),
            min(&lr),
            max(&lr),
            rs.iter().map(|r| it(r, "rho_runs")).sum::<i128>() * 2,
            mean(&sr),
            min(&sr),
            max(&sr),
            mean(&st),
            fl(rs[0], "expected_ratio"),
            py_bool(rs.iter().all(|r| bool_at(r, "all_verified"))),
        );
    }
}

fn e11(o: &mut String, rows: &[&Value]) {
    e8(o, rows, "E11", "the algebraic S₄ oracle");
    o.push_str("\n#### E11 — the solver, per instance\n\n");
    o.push_str("| p | log2 r | solver calls | F_p muls per call | Macaulay muls per call | unsolved (border unreachable) | retried at 11 / 12 / 13 | root triples | unliftable systems | agreement with the pair table (S₄ only / pair table only) | wall s per call |\n");
    o.push_str("|--:|--:|--:|--:|--:|:--|:--|--:|--:|:--|--:|\n");
    for r in sorted(rows, &["log2_r", "seed"]) {
        let sv = at(r, "solver");
        let calls = fl(sv, "calls").max(1.0);
        let a = at(r, "agreement_with_mitm");
        let agr = if truthy(a) {
            format!(
                "{}/{} ({} / {})",
                py_str(at(a, "agree")),
                py_str(at(a, "targets")),
                py_str(at(a, "s4_only")),
                py_str(at(a, "mitm_only"))
            )
        } else {
            "—".into()
        };
        let _ = writeln!(
            o,
            "| 2^{} | {:.1} | {} | {:.0} | {:.0} | {} ({}) | {} / {} / {} | {} | {} | {} | {:.4} |",
            py_str(at(r, "p_bits")),
            fl(r, "log2_r"),
            py_str(at(sv, "calls")),
            fl(sv, "fp_muls") / calls,
            fl(sv, "macaulay_muls") / calls,
            py_str(at(sv, "unsolved")),
            py_str(at(sv, "border_unreachable")),
            py_str(at(sv, "retried_at_degree_11")),
            py_str(at(sv, "retried_at_degree_12")),
            py_str(at(sv, "retried_at_degree_13")),
            py_str(at(sv, "root_triples")),
            py_str(at(sv, "unliftable_systems")),
            agr,
            fl(sv, "wall_ns") / calls / 1e9,
        );
    }
    o.push_str("\n#### E11 — solver summary\n\n");
    let calls: i128 = rows.iter().map(|r| it(r, "solver.calls")).sum();
    let muls: i128 = rows.iter().map(|r| it(r, "solver.fp_muls")).sum();
    let uns: i128 = rows.iter().map(|r| it(r, "solver.unsolved")).sum();
    let agr: Vec<&Value> = rows
        .iter()
        .map(|r| at(r, "agreement_with_mitm"))
        .filter(|a| truthy(a))
        .collect();
    let sum = |k: &str| agr.iter().map(|a| it(a, k)).sum::<i128>();
    o.push_str("| calls | F_p muls per call (mean) | unsolved | unsolved fraction | agreement with the pair table | S₄ only | pair table only |\n");
    o.push_str("|--:|--:|--:|--:|:--|--:|--:|\n");
    let _ = writeln!(
        o,
        "| {} | {:.0} | {} | {:.4} | {}/{} | {} | {} |",
        calls,
        muls as f64 / calls.max(1) as f64,
        uns,
        uns as f64 / calls.max(1) as f64,
        sum("agree"),
        sum("targets"),
        sum("s4_only"),
        sum("mitm_only"),
    );
}

// ── E12: the pair table over orbit representatives ────────────────

/// One E12 table variant (`negation_table` or `orbit_table`) of a
/// subfield row, priced in `F_p` multiplications.
struct E12Variant {
    pair_sums: f64,
    table: f64,
    stream_per_target: f64,
    folded: Arm,
    control: Arm,
}

fn e12_variant(r: &Value, name: &str) -> E12Variant {
    let k = fl(r, "inversion_in_multiplications");
    let v = at(r, name);
    let uncounted = |key: &str| {
        let x = at(v, key);
        if x.is_null() {
            0.0
        } else {
            num(x)
        }
    };
    let table = unit(at(v, "table"), k) + uncounted("build_uncounted_muls");
    let stream_total = unit(at(v, "stream_cost"), k) + uncounted("stream_uncounted_muls");
    let s = at(v, "stream");
    let (folded, control) = price_arms(
        r,
        &Priced {
            stream: s,
            stream_total,
            table,
            build_folded: unit(at(r, "base_build.folded"), k),
            build_control: unit(at(r, "base_build.control"), k),
        },
    );
    E12Variant {
        pair_sums: fl(v, "table.group_ops.adds"),
        table,
        stream_per_target: stream_total / fl(s, "trials"),
        folded,
        control,
    }
}

/// One E12 variant of an `F_p` row, priced in group additions.
fn e12p_variant(r: &Value, name: &str) -> E12Variant {
    let mul = fl(r, "mul_in_adds");
    let sqrt = fl(r, "sqrt_in_adds");
    let v = at(r, name);
    let gae = |g: &Value| fl(g, "adds") + fl(g, "doubles");
    let uncounted = |key: &str| {
        let x = at(v, key);
        if x.is_null() {
            0.0
        } else {
            num(x)
        }
    };
    let table = gae(at(v, "table_group_ops")) + mul * uncounted("build_uncounted_muls");
    let s = at(v, "stream");
    let stream_total = gae(at(s, "group_ops")) + mul * uncounted("stream_uncounted_muls");
    let build = |arm: &str| {
        let b = at(r, &format!("base_build.{arm}"));
        gae(at(b, "group_ops")) + sqrt * fl(b, "sqrt_solves") + mul * fl(b, "endomorphism_maps")
    };
    let sqrt_r = fl(r, "r").sqrt();
    let arm = |a: &str, build: f64| {
        let x = at(s, a);
        let (sqt, trials) = (at(x, "square_trials"), at(s, "trials"));
        let share = if truthy(sqt) && truthy(trials) {
            num(sqt) / num(trials)
        } else {
            1.0
        };
        let stream = stream_total * share;
        let la = mul * fl(x, "row_ops");
        let total = build + table + stream + la;
        Arm {
            build,
            table,
            stream,
            la,
            total,
            s: total / sqrt_r,
        }
    };
    E12Variant {
        pair_sums: fl(v, "table_group_ops.adds"),
        table,
        stream_per_target: stream_total / fl(s, "trials"),
        folded: arm("folded", build("folded")),
        control: arm("control", build("control")),
    }
}

fn agreement_e12(a: &Value) -> String {
    if truthy(a) {
        format!(
            "{}/{} ({} / {})",
            py_str(at(a, "agree")),
            py_str(at(a, "targets")),
            py_str(at(a, "orbit_only")),
            py_str(at(a, "negation_only"))
        )
    } else {
        "—".into()
    }
}

fn e12_ok(r: &Value) -> bool {
    both_verified(at(r, "negation_table.stream"))
        && both_verified(at(r, "orbit_table.stream"))
        && bool_at(r, "rho_all_verified")
        && it(r, "orbit_table.probe_stats.image_mismatches") == 0
}

fn e12_full(r: &Value) -> bool {
    ["negation_table", "orbit_table"].iter().all(|v| {
        truthy(at(r, &format!("{v}.stream.folded.square_relations")))
            && truthy(at(r, &format!("{v}.stream.control.square_relations")))
    })
}

fn rho_s(r: &Value, walk: &str, k_inv: Option<f64>) -> f64 {
    match k_inv {
        Some(k) => {
            let (n, f) = price_walks(r, k);
            if walk == "negation" {
                n.s
            } else {
                f.s
            }
        }
        None => mean(
            &r["walks"]
                .as_array()
                .unwrap()
                .iter()
                .map(|w| fl(&w[walk], "s"))
                .collect::<Vec<_>>(),
        ),
    }
}

/// E12 on the subfield line (`e12`, `F_p` multiplications) or on `F_p`
/// (`e12p`, group additions).
fn e12(o: &mut String, rows: &[&Value], prime: bool) {
    let (title, unit_s, size_key) = if prime {
        (
            "E12 — the orbit pair table on F_p, two summands (j = 0: key x³, w = 6; j = 1728: key x², w = 4)",
            "group additions (F_p multiplication and square root at the row's measured factors)",
            "bits",
        )
    } else {
        (
            "E12 — the orbit pair table on the Frobenius line, three summands (key: the minimal polynomial of x, w = 6)",
            "F_p multiplications (inversion = measured factor; the key's uncounted multiplications added; LA row op = 1)",
            "p_bits",
        )
    };
    let _ = writeln!(o, "### {title}\n");
    let _ = writeln!(o, "Unit: {unit_s}.  \"Negation\" is the full pair table of E8 (pairs over `⟨−1⟩`-representatives); \"orbit\" stores one entry per `H`-orbit of pairs.\n");
    o.push_str("| family | p | log2 r | pair sums built, negation / orbit | ÷ | orbit entries | table cost, negation / orbit | ÷ | stream cost per target, negation / orbit | orbit / negation (falsifies \"unchanged\" above 1.20) | keys / hits / maps | key collisions / image mismatches | full-rank rel π, negation table / orbit table | S π, negation table / orbit table | S π ratio | S negation arm, negation table / orbit table | rho S folded | S π (orbit) / rho S folded | agreement (orbit only / negation only) | correct |\n");
    o.push_str(
        "|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|:--|:--|--:|:--|--:|--:|:--|:--|\n",
    );
    let fields: [&str; 4] = ["family", size_key, "log2_r", "seed"];
    let rows = sorted(rows, &fields);
    for r in &rows {
        let size = if prime {
            format!("{}-bit", py_str(at(r, "bits")))
        } else {
            format!("2^{}", py_str(at(r, "p_bits")))
        };
        let price = |name| {
            if prime {
                e12p_variant(r, name)
            } else {
                e12_variant(r, name)
            }
        };
        let (n, b) = (price("negation_table"), price("orbit_table"));
        let ps = at(r, "orbit_table.probe_stats");
        let rho = rho_s(
            r,
            "folded",
            (!prime).then(|| fl(r, "inversion_in_multiplications")),
        );
        let full = e12_full(r);
        let s_cell = |x: f64| {
            if full {
                format!("{x:.1}")
            } else {
                "—".into()
            }
        };
        let _ = writeln!(
            o,
            "| {} | {} | {:.1} | {:.0} / {:.0} | {:.2} | {} | {} / {} | {:.2} | {} / {} | {:.3} | {} / {} / {} | {} / {} | {} / {} | {} / {} | {} | {} / {} | {:.2} | {} | {} | {} |",
            py_str(at(r, "family")),
            size,
            fl(r, "log2_r"),
            n.pair_sums,
            b.pair_sums,
            n.pair_sums / b.pair_sums,
            py_str(at(r, "orbit_table.entries")),
            g3(n.table),
            g3(b.table),
            n.table / b.table,
            g3(n.stream_per_target),
            g3(b.stream_per_target),
            b.stream_per_target / n.stream_per_target,
            py_str(at(ps, "keys")),
            py_str(at(ps, "hits")),
            py_str(at(ps, "maps")),
            py_str(at(ps, "key_collisions")),
            py_str(at(ps, "image_mismatches")),
            fmt(at(r, "negation_table.stream.folded.square_relations")),
            fmt(at(r, "orbit_table.stream.folded.square_relations")),
            s_cell(n.folded.s),
            s_cell(b.folded.s),
            if full { format!("{:.3}", n.folded.s / b.folded.s) } else { "—".into() },
            s_cell(n.control.s),
            s_cell(b.control.s),
            rho,
            if full { format!("{:.0}", b.folded.s / rho) } else { "—".into() },
            agreement_e12(at(r, "agreement")),
            yes(e12_ok(r)),
        );
    }
    let _ = writeln!(
        o,
        "\n#### {} summary\n",
        if prime {
            "E12 (F_p)"
        } else {
            "E12 (Frobenius line)"
        }
    );
    o.push_str("| family | instances (full rank on both tables) | log2 r | w | pair sums ÷ mean (min–max) | table cost ÷ mean (min–max) | stream cost per target, orbit / negation, mean (max) | table share of total π, negation table, mean (max) | S π ratio, negation / orbit table, mean (min–max) | S π (orbit) / rho S folded, min–max | key collisions | image mismatches | agreement (orbit only / negation only) | all correct |\n");
    o.push_str("|:--|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|\n");
    for (k, rs) in grouped(&rows, &["family"]) {
        let price = |r: &Value, name| {
            if prime {
                e12p_variant(r, name)
            } else {
                e12_variant(r, name)
            }
        };
        let pairs: Vec<(E12Variant, E12Variant)> = rs
            .iter()
            .map(|r| (price(r, "negation_table"), price(r, "orbit_table")))
            .collect();
        let full: Vec<usize> = (0..rs.len()).filter(|&i| e12_full(rs[i])).collect();
        let sums: Vec<f64> = pairs
            .iter()
            .map(|(n, b)| n.pair_sums / b.pair_sums)
            .collect();
        let tables: Vec<f64> = pairs.iter().map(|(n, b)| n.table / b.table).collect();
        let probe: Vec<f64> = pairs
            .iter()
            .map(|(n, b)| b.stream_per_target / n.stream_per_target)
            .collect();
        let share: Vec<f64> = full
            .iter()
            .map(|&i| pairs[i].0.folded.table / pairs[i].0.folded.total)
            .collect();
        let s_ratio: Vec<f64> = full
            .iter()
            .map(|&i| pairs[i].0.folded.s / pairs[i].1.folded.s)
            .collect();
        let vs_rho: Vec<f64> = full
            .iter()
            .map(|&i| {
                let r = rs[i];
                pairs[i].1.folded.s
                    / rho_s(
                        r,
                        "folded",
                        (!prime).then(|| fl(r, "inversion_in_multiplications")),
                    )
            })
            .collect();
        let lr: Vec<f64> = rs.iter().map(|r| fl(r, "log2_r")).collect();
        let mut ws: Vec<i128> = rs.iter().map(|r| it(r, "group_order_of_fold")).collect();
        ws.sort();
        ws.dedup();
        let tot = |p: &str| rs.iter().map(|r| it(r, p)).sum::<i128>();
        let agr: Vec<&Value> = rs
            .iter()
            .map(|r| at(r, "agreement"))
            .filter(|a| truthy(a))
            .collect();
        let asum = |k: &str| agr.iter().map(|a| it(a, k)).sum::<i128>();
        let _ = writeln!(
            o,
            "| {} | {} ({}) | {:.1}–{:.1} | {} | {} | {} | {:.3} ({:.3}) | {} | {} | {} | {} | {} | {}/{} ({} / {}) | {} |",
            ks(&k[0]),
            full.len(),
            rs.len(),
            min(&lr),
            max(&lr),
            ws.iter().map(|w| w.to_string()).collect::<Vec<_>>().join(", "),
            mean_range(&sums),
            mean_range(&tables),
            mean(&probe),
            max(&probe),
            if share.is_empty() { "—".into() } else { format!("{:.3} ({:.3})", mean(&share), max(&share)) },
            if s_ratio.is_empty() { "—".into() } else { format!("{:.3} ({:.3}–{:.3})", mean(&s_ratio), min(&s_ratio), max(&s_ratio)) },
            if vs_rho.is_empty() { "—".into() } else { format!("{:.0}–{:.0}", min(&vs_rho), max(&vs_rho)) },
            tot("orbit_table.probe_stats.key_collisions"),
            tot("orbit_table.probe_stats.image_mismatches"),
            asum("agree"),
            asum("targets"),
            asum("orbit_only"),
            asum("negation_only"),
            yes(rs.iter().all(|r| e12_ok(r))),
        );
    }
}

// ── E13: the 2-torsion symmetry of the system with the fold ────────

/// One E13 arm (a base and an oracle), priced in `F_p` multiplications.
/// The control (negation) arm of each stream never reaches full rank (its
/// rank ceiling is structural, §8.7), so only the folded arm is priced.
struct E13Arm {
    /// Phases: base build, the polynomials' set-up, the stream's group
    /// arithmetic, the solver, the linear algebra.
    phases: [f64; 5],
    columns: f64,
    square: Option<f64>,
    calls: f64,
    muls_per_call: f64,
    total: f64,
    s: f64,
}

fn e13_arm(r: &Value, stream_key: &str, base: &str) -> Option<E13Arm> {
    let v = at(r, stream_key);
    if v.is_null() {
        return None;
    }
    let k = fl(r, "inversion_in_multiplications");
    let s = at(v, "stream");
    let f = at(s, "folded");
    let share = {
        let (sqt, trials) = (at(f, "square_trials"), at(s, "trials"));
        if truthy(sqt) && truthy(trials) {
            num(sqt) / num(trials)
        } else {
            1.0
        }
    };
    let phases = [
        unit(at(r, &format!("base_build.{base}")), k),
        fl(r, "polynomials.setup_fp_muls"),
        unit(at(v, "stream_cost"), k) * share,
        fl(v, "solver.fp_muls") * share,
        fl(f, "row_ops"),
    ];
    let total: f64 = phases.iter().sum();
    let calls = fl(v, "solver.calls");
    let sq = at(f, "square_relations");
    Some(E13Arm {
        phases,
        columns: fl(f, "columns"),
        square: (!sq.is_null()).then(|| num(sq)),
        calls,
        muls_per_call: fl(v, "solver.fp_muls") / calls.max(1.0),
        total,
        s: total / fl(r, "r").sqrt(),
    })
}

/// E16: the phases of the best arm (fold 12 + D₃, E13's rows at every
/// size) against `p` and against `r`, beside the matched folded rho.
fn e16(o: &mut String, rows: &[&Value]) {
    o.push_str("### E16 — the best arm's phases against p and r, beside the matched rho (fold 12 + D₃, E13's driver)\n\n");
    o.push_str("| p bits | #E(F_p) | log2 r | log2 r − 2·log2 #E(F_p) | cols | full-rank rel | solver | group arithmetic | linear algebra | total | LA share | S | rho S folded | S / rho S |\n");
    o.push_str("|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|\n");
    let rows = sorted(rows, &["log2_r", "seed"]);
    let mut fits: BTreeMap<&str, (Vec<(f64, f64)>, Vec<(f64, f64)>)> = BTreeMap::new();
    let mut gap: Vec<(f64, f64)> = Vec::new();
    for r in &rows {
        let Some(a) = e13_arm(r, "d3_fold12", "fold12") else {
            continue;
        };
        if a.square.is_none() {
            continue;
        }
        let rr = fl(r, "r");
        let ep = fl(r, "base_order");
        let k = fl(r, "inversion_in_multiplications");
        let (_, rho) = price_walks(r, k);
        for (name, v) in [
            ("solver", a.phases[3]),
            ("group arithmetic", a.phases[2]),
            ("linear algebra", a.phases[4]),
            ("total", a.total),
            ("rho folded (total)", rho.total),
        ] {
            let e = fits.entry(name).or_default();
            e.0.push((ep, v));
            e.1.push((rr, v));
        }
        gap.push((rr, a.s / rho.s));
        let _ = writeln!(
            o,
            "| {} | {} | {:.1} | {:+.2} | {} | {:.0} | {} | {} | {} | {} | {:.4} | {:.0} | {:.1} | {:.0} |",
            py_str(at(r, "p_bits")),
            py_str(at(r, "base_order")),
            fl(r, "log2_r"),
            fl(r, "log2_r") - 2.0 * ep.log2(),
            a.columns,
            a.square.unwrap(),
            g3(a.phases[3]),
            g3(a.phases[2]),
            g3(a.phases[4]),
            g3(a.total),
            a.phases[4] / a.total,
            a.s,
            rho.s,
            a.s / rho.s,
        );
    }
    o.push_str("\n#### E16 — fitted exponents (least squares, log–log), rows at full rank\n\n");
    o.push_str("| quantity | against #E(F_p) ≈ p | against r | rows |\n");
    o.push_str("|:--|--:|--:|--:|\n");
    for (name, (by_p, by_r)) in &fits {
        let _ = writeln!(
            o,
            "| {} | {} | {} | {} |",
            name,
            fmt_f(fit_exponent(by_p)),
            fmt_f(fit_exponent(by_r)),
            by_p.len()
        );
    }
    let _ = writeln!(
        o,
        "| S / rho S | — | {} | {} |",
        fmt_f(fit_exponent(&gap)),
        gap.len()
    );
}

fn e17_ok(r: &Value) -> bool {
    ["d3", "s3"]
        .iter()
        .map(|k| at(r, k))
        .filter(|v| !v.is_null())
        .all(|v| bool_at(v, "stream.folded.verified") && bool_at(v, "stream.control.verified"))
        && bool_at(r, "rho_all_verified")
}

/// E17 (goal G2): the D₃ oracle and the `τ_T` fold on the full group
/// `E(F_{p³})`, `r ≈ p³/h`: the solve constant against S₃, agreement
/// with the pair table, every phase of the `D₃` arm against the matched
/// negation rho, the fitted exponents, and the model `S/rho S = A·r^a +
/// B·r^b` (relation part, linear algebra) the fits imply.
fn e17(o: &mut String, rows: &[&Value]) {
    let rows = sorted(rows, &["log2_r", "seed"]);
    o.push_str("### E17 — the D₃ oracle and the τ_T fold on the full group E(F_{p³}): solve constant and agreement\n\n");
    o.push_str("| p | log2 r | h | cols ⟨−1, τ_T⟩ / ⟨−1⟩ | ratio | D₃ calls | D₃ muls per call | R(q₃) degree max | D₃ unsolved | S₃ muls per call | S₃ / D₃ | agreement D₃ / S₃ disagreements, all (non-degenerate) of targets; pair-table hits, cancelling | correct |\n");
    o.push_str("|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|:--|:--|\n");
    let mut ratios = Vec::new();
    for r in &rows {
        let d3 = e13_arm(r, "d3", "fold4").unwrap();
        let s3 = e13_arm(r, "s3", "fold4");
        let ag = at(r, "agreement");
        let agr = if truthy(ag) {
            format!(
                "{} ({}) / {} ({}) of {}; {}, {}",
                py_str(at(ag, "d3_disagreements")),
                py_str(at(ag, "d3_disagreements_nondegenerate")),
                py_str(at(ag, "s3_disagreements")),
                py_str(at(ag, "s3_disagreements_nondegenerate")),
                py_str(at(ag, "targets")),
                py_str(at(ag, "pair_table_hits")),
                py_str(at(ag, "pair_table_cancelling")),
            )
        } else {
            "—".into()
        };
        let cols = (fl(r, "columns.fold4"), fl(r, "columns.control"));
        if let Some(s) = &s3 {
            ratios.push(s.muls_per_call / d3.muls_per_call);
        }
        let _ = writeln!(
            o,
            "| {} | {:.1} | {} | {:.0} / {:.0} | {:.2} | {:.0} | {:.0} | {} | {} | {} | {} | {} | {} |",
            py_str(at(r, "p")),
            fl(r, "log2_r"),
            py_str(at(r, "cofactor")),
            cols.0,
            cols.1,
            cols.1 / cols.0,
            d3.calls,
            d3.muls_per_call,
            py_str(at(r, "d3.solver.resultant_degree_max")),
            py_str(at(r, "d3.solver.unsolved")),
            s3.as_ref().map_or("—".into(), |s| format!("{:.0}", s.muls_per_call)),
            s3.as_ref().map_or("—".into(), |s| format!("{:.1}", s.muls_per_call / d3.muls_per_call)),
            agr,
            yes(e17_ok(r)),
        );
    }
    if !ratios.is_empty() {
        let _ = writeln!(
            o,
            "\nS₃ / D₃ muls per call over {} rows: {}.\n",
            ratios.len(),
            mean_range(&ratios)
        );
    }
    o.push_str("### E17 — every phase of the D₃ arm (base ⟨−1, τ_T⟩) and of the S₃ arm on the same base, beside the matched negation rho (F_p multiplications)\n\n");
    o.push_str("| p | log2 r | h | cols | rel to pinned | set-up | base | group arithmetic | solver | linear algebra | total D₃ | S D₃ | rho S | S / rho S, D₃ | S / rho S, S₃ |\n");
    o.push_str("|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|\n");
    let mut fits: BTreeMap<&str, Vec<(f64, f64)>> = BTreeMap::new();
    let (mut rel_gap, mut la_gap) = (Vec::new(), Vec::new());
    for r in &rows {
        let d3 = e13_arm(r, "d3", "fold4").unwrap();
        let Some(sq) = d3.square else {
            continue;
        };
        let s3 = e13_arm(r, "s3", "fold4");
        let k = fl(r, "inversion_in_multiplications");
        let rho = price_walk(r, k, "negation");
        let rr = fl(r, "r");
        for (name, v) in [
            ("solver", d3.phases[3]),
            ("group arithmetic", d3.phases[2]),
            ("linear algebra", d3.phases[4]),
            ("total", d3.total),
            ("rho (total)", rho.total),
            ("S / rho S, D₃", d3.s / rho.s),
        ] {
            fits.entry(name).or_default().push((rr, v));
        }
        if let Some(s) = &s3 {
            fits.entry("S / rho S, S₃")
                .or_default()
                .push((rr, s.s / rho.s));
        }
        rel_gap.push((rr, (d3.total - d3.phases[4]) / rho.total));
        la_gap.push((rr, d3.phases[4] / rho.total));
        let _ = writeln!(
            o,
            "| {} | {:.1} | {} | {:.0} | {:.0} | {} | {} | {} | {} | {} | {} | {:.0} | {:.2} | {:.1} | {} |",
            py_str(at(r, "p")),
            fl(r, "log2_r"),
            py_str(at(r, "cofactor")),
            d3.columns,
            sq,
            g3(d3.phases[1]),
            g3(d3.phases[0]),
            g3(d3.phases[2]),
            g3(d3.phases[3]),
            g3(d3.phases[4]),
            g3(d3.total),
            d3.s,
            rho.s,
            d3.s / rho.s,
            s3.as_ref().map_or("—".into(), |s| format!("{:.0}", s.s / rho.s)),
        );
    }
    o.push_str("\n#### E17 — fitted exponents against r (least squares, log–log), rows with the logarithm pinned\n\n");
    o.push_str("| quantity | exponent | rows |\n|:--|--:|--:|\n");
    for (name, pts) in &fits {
        let _ = writeln!(
            o,
            "| {} | {} | {} |",
            name,
            fmt_f(fit_exponent(pts)),
            pts.len()
        );
    }
    // The model the fits imply: S/rho S = A·r^a (everything but the linear
    // algebra) + B·r^b (the linear algebra), each fitted on its own.
    let fit_line = |pts: &[(f64, f64)]| -> Option<(f64, f64)> {
        let a = fit_exponent(pts)?;
        let lx: Vec<f64> = pts.iter().map(|p| p.0.log2()).collect();
        let ly: Vec<f64> = pts.iter().map(|p| p.1.log2()).collect();
        Some((a, mean(&ly) - a * mean(&lx)))
    };
    if let (Some((a, ca)), Some((b, cb))) = (fit_line(&rel_gap), fit_line(&la_gap)) {
        let _ = writeln!(
            o,
            "\nModel (extrapolation, not measurement): S / rho S = 2^{{{ca:.2}}}·r^{{{a:.3}}} (relation part) + 2^{{{cb:.2}}}·r^{{{b:.3}}} (linear algebra)."
        );
        if a < 0.0 && b > 0.0 {
            // d/dlog r of A r^a + B r^b = 0  ⇒  log2 r* = (log2(−aA) − log2(bB)) / (b − a).
            let l_star = ((-a).log2() + ca - b.log2() - cb) / (b - a);
            let s_star = 2f64.powf(ca + a * l_star) + 2f64.powf(cb + b * l_star);
            let _ = writeln!(
                o,
                "Its minimum: S / rho S = {s_star:.3} at r ≈ 2^{l_star:.1}; the relation part alone reaches 1 at r ≈ 2^{:.1}.",
                -ca / a
            );
        }
    }
}

/// One E18 row's walks of one fold: mean `F_p` multiplications, mean
/// steps, mean `S = muls / √n`, and whether every walk recovered `d`.
fn e18_walk(r: &Value, fold: &str) -> Option<(f64, f64, Vec<f64>, bool)> {
    let w: Vec<&Value> = r["rho"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|w| w["fold"].as_str() == Some(fold))
        .collect();
    if w.is_empty() {
        return None;
    }
    let sq = fl(r, "n").sqrt();
    let muls: Vec<f64> = w.iter().map(|w| fl(w, "fp_muls")).collect();
    let steps: Vec<f64> = w.iter().map(|w| fl(w, "steps")).collect();
    Some((
        mean(&muls),
        mean(&steps),
        muls.iter().map(|m| m / sq).collect(),
        w.iter().all(|w| bool_at(w, "correct")),
    ))
}

/// E18 (goal G2, note §9.1): folds of the base on the `k = 4` linear
/// algebra.  `r = LA·16 / rho` in `F_p` multiplications (a mod-`n`
/// multiplication is `16`, §11.16), per curve against that curve's walks,
/// and pooled the way §11.19 pools rho: `S` per fold over every walk of
/// the family, `r = LA·16 / (S·√n)`.
fn e18(o: &mut String, rows: &[&Value]) {
    let rows = sorted(rows, &["family", "bits", "seed"]);
    o.push_str(
        "### E18 — folds of the base on the k = 4 linear algebra, against the matched rho\n\n",
    );
    o.push_str("| family | p | log2 n | h | base | cols control / folded | rate | targets with > 1 decomposition | control core (φ) | folded core (φ) | duplicates skipped (folded) | LA control | LA folded | LA ratio | rho steps unfolded / negation / ψ | correct (control, folded, walks) |\n");
    o.push_str("|:--|--:|--:|--:|--:|:--|--:|--:|:--|:--|--:|--:|--:|--:|:--|:--|\n");
    let mut pooled: BTreeMap<(String, String), Vec<f64>> = BTreeMap::new();
    for r in &rows {
        let fam = py_str(at(r, "family"));
        let hist: Vec<f64> = r["decompositions_per_target"]
            .as_array()
            .unwrap()
            .iter()
            .map(num)
            .collect();
        let multi: f64 = hist[2..].iter().sum();
        let walks: Vec<_> = ["None", "Negation", "Psi"]
            .iter()
            .map(|f| e18_walk(r, f))
            .collect();
        for (f, w) in ["None", "Negation", "Psi"].iter().zip(&walks) {
            if let Some((_, _, s, _)) = w {
                pooled
                    .entry((fam.clone(), f.to_string()))
                    .or_default()
                    .extend(s);
            }
        }
        let steps = walks
            .iter()
            .map(|w| w.as_ref().map_or("—".into(), |w| format!("{:.0}", w.1)))
            .collect::<Vec<_>>()
            .join(" / ");
        let walks_ok = walks.iter().flatten().all(|w| w.3);
        let _ = writeln!(
            o,
            "| {} | {} | {:.1} | {} | {} | {} / {} | {:.4} | {:.0} of {:.0} | {} ({:.3}) | {} ({:.3}) | {} | {} | {} | {:.2} | {} | {}, {}, {} |",
            fam,
            py_str(at(r, "p")),
            fl(r, "bits"),
            py_str(at(r, "cofactor")),
            py_str(at(r, "base_points")),
            py_str(at(r, "control.columns")),
            py_str(at(r, "folded.columns")),
            fl(r, "decomposition_rate"),
            multi,
            fl(r, "targets_decomposed"),
            py_str(at(r, "control.unknowns")),
            fl(r, "control.phi"),
            py_str(at(r, "folded.unknowns")),
            fl(r, "folded.phi"),
            py_str(at(r, "folded.duplicate_rows")),
            g3(fl(r, "control.la_ops")),
            g3(fl(r, "folded.la_ops")),
            fl(r, "control.la_ops") / fl(r, "folded.la_ops"),
            steps,
            yes(bool_at(r, "control.correct")),
            yes(bool_at(r, "folded.correct")),
            yes(walks_ok),
        );
    }
    o.push_str("\n#### E18 — rho's S per fold, pooled over every walk of the family (F_p multiplications / √n)\n\n");
    o.push_str(
        "| family | walk | S (mean ± s.e.) | walks | unfolded / this |\n|:--|:--|--:|--:|--:|\n",
    );
    let stat = |v: &[f64]| -> (f64, f64) {
        let m = mean(v);
        let sd =
            (v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (v.len().max(2) - 1) as f64).sqrt();
        (m, sd / (v.len() as f64).sqrt())
    };
    for ((fam, fold), v) in &pooled {
        let (m, se) = stat(v);
        let base = stat(&pooled[&(fam.clone(), "None".to_string())]).0;
        let _ = writeln!(
            o,
            "| {fam} | {fold} | {m:.3} ± {se:.3} | {} | {:.3} |",
            v.len(),
            base / m
        );
    }
    o.push_str("\n#### E18 — r = LA / rho, against the pooled S of each walk\n\n");
    o.push_str("| family | p | log2 n | r control, unfolded | r control, negation | r control, ψ | r folded, unfolded | r folded, negation | r folded, ψ |\n");
    o.push_str("|:--|--:|--:|--:|--:|--:|--:|--:|--:|\n");
    let mut fits: BTreeMap<String, Vec<(f64, f64)>> = BTreeMap::new();
    let mut means: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    for r in &rows {
        let fam = py_str(at(r, "family"));
        let sq = fl(r, "n").sqrt();
        let cell = |arm: &str, fold: &str| -> Option<f64> {
            let v = pooled.get(&(fam.clone(), fold.to_string()))?;
            Some(fl(r, &format!("{arm}.la_ops")) * 16.0 / (mean(v) * sq))
        };
        let mut cells = Vec::new();
        for arm in ["control", "folded"] {
            for fold in ["None", "Negation", "Psi"] {
                let v = cell(arm, fold);
                if let Some(x) = v {
                    let key = format!("{fam} | {arm} | {fold}");
                    fits.entry(key.clone()).or_default().push((fl(r, "n"), x));
                    means.entry(key).or_default().push(x);
                }
                cells.push(v.map_or("—".into(), |x| format!("{x:.3}")));
            }
        }
        let _ = writeln!(
            o,
            "| {} | {} | {:.1} | {} |",
            fam,
            py_str(at(r, "p")),
            fl(r, "bits"),
            cells.join(" | ")
        );
    }
    o.push_str("\n| family | arm | walk | r, mean ± s.e. over curves | fitted exponent in n | curves |\n|:--|:--|:--|--:|--:|--:|\n");
    for (key, v) in &means {
        let (m, se) = stat(v);
        let _ = writeln!(
            o,
            "| {key} | {m:.3} ± {se:.3} | {} | {} |",
            fmt_f(fit_exponent(&fits[key])),
            v.len()
        );
    }
}

fn e13_ok(r: &Value) -> bool {
    ["d3_fold12", "d3_fold6", "s3_fold12", "s3_fold6"]
        .iter()
        .map(|k| at(r, k))
        .filter(|v| !v.is_null())
        .all(|v| bool_at(v, "stream.folded.verified") && bool_at(v, "stream.control.verified"))
        && bool_at(r, "rho_all_verified")
        && it(r, "d3_fold12.solver.unsolved") == 0
        && it(r, "d3_fold6.solver.unsolved") == 0
}

fn e13(o: &mut String, rows: &[&Value]) {
    o.push_str("### E13 — FGHR's 2-torsion symmetry and the fold, on the Y-line of a subfield curve: relations and solvers\n\n");
    o.push_str("Arms: base fold 12 = `⟨−1, π, τ_T⟩`, fold 6 = `⟨−1, π⟩`, negation = `⟨−1⟩`; oracle D₃ = the system symmetrised by even sign changes and permutations (16 solutions, conic resultant), S₃ = by permutations only (64, Macaulay).\n\n");
    o.push_str("| p | log2 r | h | cols fold 12 / fold 6 / negation | full-rank rel fold 12 / fold 6 (D₃) | full-rank rel fold 12 / fold 6 (S₃) | negation rank at stop, D₃ streams (of cols) | D₃ calls fold 12 / fold 6 | D₃ muls per call | R(q₃) degree max | D₃ unsolved | S₃ muls per call | S₃ quotient per call | S₃ / D₃ muls per call | agreement with the pair table: D₃ / S₃ disagreements (targets, pair-table hits) | correct |\n");
    o.push_str("|--:|--:|--:|:--|:--|:--|:--|:--|--:|--:|--:|--:|--:|--:|:--|:--|\n");
    let rows = sorted(rows, &["log2_r", "seed"]);
    for r in &rows {
        let a = e13_arm(r, "d3_fold12", "fold12").unwrap();
        let b = e13_arm(r, "d3_fold6", "fold6").unwrap();
        let c = e13_arm(r, "s3_fold12", "fold12");
        let d = e13_arm(r, "s3_fold6", "fold6");
        let sq = |x: &E13Arm| x.square.map_or("—".to_string(), |v| format!("{v:.0}"));
        let s3pair = match (&c, &d) {
            (Some(c), Some(d)) => format!("{} / {}", sq(c), sq(d)),
            (Some(c), None) => format!("{} / —", sq(c)),
            _ => "—".into(),
        };
        let ag = at(r, "agreement");
        let agr = if truthy(ag) {
            format!(
                "{} / {} ({}, {})",
                py_str(at(ag, "d3_disagreements")),
                py_str(at(ag, "s3_disagreements")),
                py_str(at(ag, "targets")),
                py_str(at(ag, "pair_table_hits"))
            )
        } else {
            "—".into()
        };
        let s3_quot = at(r, "s3_fold12.solver.quotient_dim_total");
        let _ = writeln!(
            o,
            "| 2^{} | {:.1} | {} | {} / {} / {} | {} / {} | {} | {} / {} ({}) | {:.0} / {:.0} | {:.0} | {} | {} | {} | {} | {} | {} | {} |",
            py_str(at(r, "p_bits")),
            fl(r, "log2_r"),
            py_str(at(r, "cofactor")),
            py_str(at(r, "columns.fold12")),
            py_str(at(r, "columns.fold6")),
            py_str(at(r, "columns.control")),
            sq(&a),
            sq(&b),
            s3pair,
            py_str(at(r, "d3_fold12.stream.control.rank_final")),
            py_str(at(r, "d3_fold6.stream.control.rank_final")),
            py_str(at(r, "columns.control")),
            a.calls,
            b.calls,
            (fl(r, "d3_fold12.solver.fp_muls") + fl(r, "d3_fold6.solver.fp_muls"))
                / (a.calls + b.calls).max(1.0),
            py_str(at(r, "d3_fold12.solver.resultant_degree_max")),
            it(r, "d3_fold12.solver.unsolved") + it(r, "d3_fold6.solver.unsolved"),
            c.as_ref().map_or("—".into(), |c| format!("{:.0}", c.muls_per_call)),
            if s3_quot.is_null() {
                "—".into()
            } else {
                format!("{:.1}", num(s3_quot) / c.as_ref().unwrap().calls.max(1.0))
            },
            c.as_ref()
                .map_or("—".into(), |c| format!("{:.1}", c.muls_per_call / a.muls_per_call)),
            agr,
            yes(e13_ok(r)),
        );
    }
    o.push_str("\n#### E13 — every phase priced (F_p multiplications; inversion = measured factor; each solver's own count; the polynomials' set-up charged to every arm; LA row op = 1); folded arms only\n\n");
    o.push_str("| p | log2 r | S fold 12 + D₃ | S fold 6 + D₃ | S fold 12 + S₃ | S fold 6 + S₃ | base lever: S(6, D₃) / S(12, D₃) | system lever: S(12, S₃) / S(12, D₃) | product | combined: S(6, S₃) / S(12, D₃) | combined / product (falsifies \"multiply\" below 0.80) | rho S folded | S(12, D₃) / rho S folded |\n");
    o.push_str("|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|\n");
    let mut fits: BTreeMap<&str, Vec<(f64, f64)>> = BTreeMap::new();
    let mut summary: Vec<(f64, f64, Option<(f64, f64, f64)>, f64)> = Vec::new();
    for r in &rows {
        let a = e13_arm(r, "d3_fold12", "fold12").unwrap();
        let b = e13_arm(r, "d3_fold6", "fold6").unwrap();
        let c = e13_arm(r, "s3_fold12", "fold12");
        let d = e13_arm(r, "s3_fold6", "fold6");
        let rho = rho_s(r, "folded", Some(fl(r, "inversion_in_multiplications")));
        let full_ab = a.square.is_some() && b.square.is_some();
        let rr = fl(r, "r");
        let cell = |x: Option<&E13Arm>| match x {
            Some(x) if x.square.is_some() => format!("{:.1}", x.s),
            Some(_) => "not full rank".into(),
            None => "—".into(),
        };
        let lever = full_ab.then(|| b.s / a.s);
        let sys = c
            .as_ref()
            .filter(|c| c.square.is_some() && a.square.is_some())
            .map(|c| c.s / a.s);
        let comb = d
            .as_ref()
            .filter(|d| d.square.is_some() && a.square.is_some())
            .map(|d| d.s / a.s);
        if a.square.is_some() {
            fits.entry("fold12_d3").or_default().push((rr, a.total));
            fits.entry("rho_folded")
                .or_default()
                .push((rr, rho * rr.sqrt()));
        }
        if b.square.is_some() {
            fits.entry("fold6_d3").or_default().push((rr, b.total));
        }
        if let Some(l) = lever {
            summary.push((
                l,
                a.s / rho,
                match (sys, comb) {
                    (Some(s), Some(c)) => Some((s, c, c / (l * s))),
                    _ => None,
                },
                a.columns / b.columns,
            ));
        }
        let _ = writeln!(
            o,
            "| 2^{} | {:.1} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {:.1} | {} |",
            py_str(at(r, "p_bits")),
            fl(r, "log2_r"),
            cell(Some(&a)),
            cell(Some(&b)),
            cell(c.as_ref()),
            cell(d.as_ref()),
            fmt_f(lever),
            fmt_f(sys),
            fmt_f(lever.zip(sys).map(|(l, s)| l * s)),
            fmt_f(comb),
            fmt_f(lever.zip(sys).zip(comb).map(|((l, s), c)| c / (l * s))),
            rho,
            if a.square.is_some() {
                format!("{:.0}", a.s / rho)
            } else {
                "—".into()
            },
        );
    }
    o.push_str("\n#### E13 — the phases of fold 12 + D₃, rows at full rank (F_p multiplications; share of the total)\n\n");
    o.push_str("| p | log2 r | base build | polynomial set-up | stream group arithmetic | D₃ solver | linear algebra | total | dominant phase |\n");
    o.push_str("|--:|--:|--:|--:|--:|--:|--:|--:|:--|\n");
    let names = [
        "base build",
        "set-up",
        "group arithmetic",
        "solver",
        "linear algebra",
    ];
    let mut by_p: Vec<(f64, f64)> = Vec::new();
    for r in &rows {
        let a = e13_arm(r, "d3_fold12", "fold12").unwrap();
        if a.square.is_none() {
            continue;
        }
        by_p.push((2f64.powf(fl(r, "p_bits")), a.total));
        let cell = |x: f64| format!("{} ({:.2})", g3(x), x / a.total);
        let top = (0..5)
            .max_by(|&i, &j| a.phases[i].total_cmp(&a.phases[j]))
            .unwrap();
        let _ = writeln!(
            o,
            "| 2^{} | {:.1} | {} | {} | {} | {} | {} | {} | {} |",
            py_str(at(r, "p_bits")),
            fl(r, "log2_r"),
            cell(a.phases[0]),
            cell(a.phases[1]),
            cell(a.phases[2]),
            cell(a.phases[3]),
            cell(a.phases[4]),
            g3(a.total),
            names[top],
        );
    }
    let _ = writeln!(
        o,
        "\nFitted exponent of the fold 12 + D₃ total against p (the line's field size, over the {} rows above; the stream is p·ln p solves at a cost that grows with log p, so 1 + o(1) is expected): {}.",
        by_p.len(),
        fmt_f(fit_exponent(&by_p))
    );
    o.push_str("\n#### E13 summary (rows where both D₃ arms reached full rank)\n\n");
    o.push_str("| instances (of run) | log2 r | column ratio fold 12 / fold 6 | base lever S(6, D₃) / S(12, D₃), mean (min–max) | rows with S₃ on both bases | system lever S(12, S₃) / S(12, D₃), mean (min–max) | combined S(6, S₃) / S(12, D₃), mean (min–max) | combined / product, min | S(12, D₃) / rho S folded, min–max | fitted exponent of total against r: fold 12 + D₃, fold 6 + D₃, rho folded (rho: 0.50) | all correct |\n");
    o.push_str("|--:|:--|--:|--:|--:|--:|--:|--:|--:|:--|:--|\n");
    let lr: Vec<f64> = rows.iter().map(|r| fl(r, "log2_r")).collect();
    let lev: Vec<f64> = summary.iter().map(|x| x.0).collect();
    let vs: Vec<f64> = summary.iter().map(|x| x.1).collect();
    let both: Vec<(f64, f64, f64)> = summary.iter().filter_map(|x| x.2).collect();
    let sys: Vec<f64> = both.iter().map(|x| x.0).collect();
    let comb: Vec<f64> = both.iter().map(|x| x.1).collect();
    let mult: Vec<f64> = both.iter().map(|x| x.2).collect();
    let cols: Vec<f64> = summary.iter().map(|x| py_round(x.3, 2)).collect();
    let exps = ["fold12_d3", "fold6_d3", "rho_folded"]
        .iter()
        .map(|k| fmt_f(fit_exponent(fits.get(k).map_or(&[][..], |v| v))))
        .collect::<Vec<_>>()
        .join(", ");
    let _ = writeln!(
        o,
        "| {} ({}) | {:.1}–{:.1} | {} | {} | {} | {} | {} | {} | {} | {} | {} |",
        summary.len(),
        rows.len(),
        min(&lr),
        max(&lr),
        joined_set(cols),
        mean_range(&lev),
        both.len(),
        mean_range(&sys),
        mean_range(&comb),
        if mult.is_empty() {
            "—".into()
        } else {
            format!("{:.2}", min(&mult))
        },
        if vs.is_empty() {
            "—".into()
        } else {
            format!("{:.0}–{:.0}", min(&vs), max(&vs))
        },
        exps,
        yes(rows.iter().all(|r| e13_ok(r))),
    );
}

// ── E14: Q-curves of degree 2 and 3 over F_{p²} ────────────────────

fn e14(o: &mut String, rows: &[&Value]) {
    o.push_str("### E14 — Q-curves of degree 2 and 3 over F_{p²}: ψ = π ∘ ι ∘ φ against the F_p-line base\n\n");
    o.push_str("| d | p | log2 r | h | twist | seed | map | λ² ≡ | ord_r(λ) | images in base | base points | chance fraction | expected at chance | verified |\n");
    o.push_str("|--:|--:|--:|--:|:--|--:|:--|:--|--:|--:|--:|--:|--:|:--|\n");
    let sign = |m: &Value| match at(m, "square_is_plus_d") {
        Value::Bool(true) => "+d",
        Value::Bool(false) => "−d",
        _ => "neither",
    };
    let sorted_rows = sorted(rows, &["degree", "log2_r", "seed"]);
    for r in &sorted_rows {
        for m in r["maps"].as_array().unwrap() {
            let ov = at(m, "overlap");
            let _ = writeln!(
                o,
                "| {} | {} | {:.1} | {} | {} | {} | `{}` | {} | {} | {} | {} | {:.5} | {:.2} | {} |",
                py_str(at(r, "degree")),
                py_str(at(r, "p")),
                fl(r, "log2_r"),
                py_str(at(r, "cofactor")),
                if bool_at(r, "twisted") { "yes" } else { "no" },
                py_str(at(r, "seed")),
                py_str(at(m, "name")),
                sign(m),
                py_str(at(ov, "eigenvalue_order")),
                py_str(at(ov, "images_in_base")),
                py_str(at(ov, "base_points")),
                fl(ov, "chance_fraction"),
                fl(ov, "chance_fraction") * fl(ov, "base_points"),
                py_str(at(m, "verified")),
            );
        }
    }
    o.push_str("\n#### E14 summary\n\n");
    o.push_str("| d | instances | maps | log2 r | λ² ≡ +d / −d | ord_r(λ) min–max | ord_r(λ) / r, min | images in base, total | expected at chance, total | base points, total | all verified |\n");
    o.push_str("|--:|--:|--:|:--|:--|--:|--:|--:|--:|--:|:--|\n");
    for (k, rs) in grouped(&sorted_rows, &["degree"]) {
        let maps: Vec<(&Value, &Value)> = rs
            .iter()
            .flat_map(|r| r["maps"].as_array().unwrap().iter().map(move |m| (*r, m)))
            .collect();
        let ords: Vec<i128> = maps
            .iter()
            .map(|(_, m)| it(m, "overlap.eigenvalue_order"))
            .collect();
        let rel: Vec<f64> = maps
            .iter()
            .map(|(r, m)| fl(m, "overlap.eigenvalue_order") / fl(r, "r"))
            .collect();
        let inside: i128 = maps
            .iter()
            .map(|(_, m)| it(m, "overlap.images_in_base"))
            .sum();
        let expected: f64 = maps
            .iter()
            .map(|(_, m)| fl(m, "overlap.chance_fraction") * fl(m, "overlap.base_points"))
            .sum();
        let pts: i128 = maps.iter().map(|(_, m)| it(m, "overlap.base_points")).sum();
        let plus = maps.iter().filter(|(_, m)| sign(m) == "+d").count();
        let minus = maps.iter().filter(|(_, m)| sign(m) == "−d").count();
        let lr: Vec<f64> = rs.iter().map(|r| fl(r, "log2_r")).collect();
        let _ = writeln!(
            o,
            "| {} | {} | {} | {:.1}–{:.1} | {} / {} | {}–{} | {:.3} | {} | {:.1} | {} | {} |",
            ks(&k[0]),
            rs.len(),
            maps.len(),
            min(&lr),
            max(&lr),
            plus,
            minus,
            ords.iter().min().copied().unwrap_or(0),
            ords.iter().max().copied().unwrap_or(0),
            min(&rel),
            inside,
            expected,
            pts,
            yes(maps.iter().all(|(_, m)| bool_at(m, "verified"))),
        );
    }
}

// ── E15: the ECC2K-130 family, two summands ────────────────────────

fn e15(o: &mut String, rows: &[&Value]) {
    o.push_str("### E15 — two summands on E_0: y² + xy = x³ + 1 over GF(2^m), the Frobenius-orbit base against the abscissa control\n\n");
    o.push_str("| curve | m | ord_m(2) | intermediate subfields | log2 r | h | h = #E_0(F_2) | subspace dim | seed | cols fold / control | D fold / control | single-column rows fold / control (of relations) | chance, 1 / cols fold | full-rank rel fold / control | ratio | column ratio | rank fraction at k = cols, fold / control | correct |\n");
    o.push_str("|:--|--:|--:|:--|--:|:--|:--|--:|--:|:--|:--|:--|--:|:--|--:|--:|:--|:--|\n");
    let factors = |r: &Value| {
        r["cofactor_factors"]
            .as_array()
            .unwrap()
            .iter()
            .map(py_str)
            .collect::<Vec<_>>()
            .join("·")
    };
    let subfields = |r: &Value| {
        let v: Vec<String> = r["intermediate_subfield_degrees"]
            .as_array()
            .unwrap()
            .iter()
            .map(|d| format!("GF(2^{})", py_str(d)))
            .collect();
        if v.is_empty() {
            "none".to_string()
        } else {
            v.join(", ")
        }
    };
    let sorted_rows = sorted(rows, &["n", "seed"]);
    for r in &sorted_rows {
        let s = at(r, "stream");
        let (f, c) = (at(s, "folded"), at(s, "control"));
        let rel = fl(s, "relations").max(1.0);
        let _ = writeln!(
            o,
            "| `{}` | {} | {} | {} | {:.1} | {} = {} | {} | {} | {} | {} / {} | {} / {} | {} ({:.3}) / {} ({:.3}) | {:.3} | {} / {} | {} | {:.2} | {} / {} | {} |",
            py_str(at(r, "instance")),
            py_str(at(r, "n")),
            py_str(at(r, "ord_n_2")),
            subfields(r),
            fl(r, "log2_r"),
            py_str(at(r, "cofactor")),
            factors(r),
            if bool_at(r, "cofactor_is_e0_f2") { "yes" } else { "no" },
            py_str(at(r, "subspace_dimension")),
            py_str(at(r, "seed")),
            py_str(at(f, "columns")),
            py_str(at(c, "columns")),
            py_str(at(f, "deficiency_total")),
            py_str(at(c, "deficiency_total")),
            py_str(at(f, "single_column_rows")),
            fl(f, "single_column_rows") / rel,
            py_str(at(c, "single_column_rows")),
            fl(c, "single_column_rows") / rel,
            1.0 / fl(f, "columns"),
            fmt(at(f, "square_relations")),
            fmt(at(c, "square_relations")),
            fmt(at(s, "square_ratio")),
            fl(s, "column_ratio"),
            fmt(at(f, "rank_fraction_at_columns")),
            fmt(at(c, "rank_fraction_at_columns")),
            yes(both_verified(s)),
        );
    }
    o.push_str("\n#### E15 summary by degree\n\n");
    o.push_str("| curve | m | h | h = #E_0(F_2) | instances | D fold / control | single-column fraction, fold, mean (min–max) | against chance | full-rank ratio mean (min–max) | column ratio | all correct |\n");
    o.push_str("|:--|--:|:--|:--|--:|:--|--:|--:|--:|--:|:--|\n");
    for (_, rs) in grouped(&sorted_rows, &["n"]) {
        let r0 = rs[0];
        let frac: Vec<f64> = rs
            .iter()
            .map(|r| fl(r, "stream.folded.single_column_rows") / fl(r, "stream.relations").max(1.0))
            .collect();
        let chance = 1.0 / fl(r0, "stream.folded.columns");
        let sq = present(&rs, "stream.square_ratio");
        let mut ds: Vec<String> = rs
            .iter()
            .map(|r| {
                format!(
                    "{} / {}",
                    py_str(at(r, "stream.folded.deficiency_total")),
                    py_str(at(r, "stream.control.deficiency_total"))
                )
            })
            .collect();
        ds.sort();
        ds.dedup();
        let _ = writeln!(
            o,
            "| `{}` | {} | {} = {} | {} | {} | {} | {:.3} ({:.3}–{:.3}) | {:.1}× | {} | {:.2} | {} |",
            py_str(at(r0, "instance")),
            py_str(at(r0, "n")),
            py_str(at(r0, "cofactor")),
            factors(r0),
            if bool_at(r0, "cofactor_is_e0_f2") { "yes" } else { "no" },
            rs.len(),
            ds.join(", "),
            mean(&frac),
            min(&frac),
            max(&frac),
            mean(&frac) / chance,
            mean_range(&sq),
            fl(r0, "stream.column_ratio"),
            yes(rs.iter().all(|r| both_verified(at(r, "stream")))),
        );
    }
}

fn main() {
    let mut rows: Vec<Value> = Vec::new();
    for path in env::args().skip(1) {
        let text = fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
        let v: Value = serde_json::from_str(&text).unwrap_or_else(|e| panic!("{path}: {e}"));
        // `{"rows": [...]}` (this binary's runners) or a bare array (`quartic_folds`).
        rows.extend(
            v["rows"]
                .as_array()
                .or(v.as_array())
                .cloned()
                .unwrap_or_default(),
        );
    }
    let mut by: BTreeMap<String, Vec<&Value>> = BTreeMap::new();
    for r in &rows {
        by.entry(py_str(at(r, "experiment"))).or_default().push(r);
    }
    let mut o = String::new();
    type Printer = fn(&mut String, &[&Value]);
    let printers: [(&str, Printer); 18] = [
        ("e1", e1),
        ("e2", e2),
        ("e3", e3),
        ("e4", e4),
        ("e5", e5),
        ("e6", e6),
        ("e7", e7),
        ("e8", |o, r| e8(o, r, "E8", "the pair-table oracle")),
        ("e9", e9),
        ("e11", e11),
        ("e12", |o, r| e12(o, r, false)),
        ("e12p", |o, r| e12(o, r, true)),
        ("e13", e13),
        ("e13", e16),
        ("e17", e17),
        ("e18", e18),
        ("e14", e14),
        ("e15", e15),
    ];
    for (exp, print) in printers {
        if let Some(rs) = by.get(exp) {
            print(&mut o, rs);
            o.push('\n');
        }
    }
    print!("{o}");
}

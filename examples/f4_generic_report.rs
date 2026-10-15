//! Analysis for `research/f4_generic_field_20261015/` (the F4 tower engine
//! over a generic prime field). Every number in the study's tables comes
//! from here; `plot.py` only lays out the CSVs this writes.
//!
//! ```bash
//! # the pre-registered speed claims (PROTOCOL.md §4)
//! cargo run --release --example f4_generic_report -- speed runs/speed OUTDIR
//! # the measured cells against the recorded references (PROTOCOL.md §5–6)
//! cargo run --release --example f4_generic_report -- cells OUTDIR \
//!     --ref ref/*.jsonl --rows runs/depth/*.jsonl
//! ```
//!
//! `speed` reads `CELL.ARM.repN.jsonl` (one row each), checks that the arms
//! of a cell agree on every field but the clock and the backend label (gate
//! G3), and writes `raw.csv`, `summary.csv`, `claims.csv` and `tables.md`.
//!
//! `cells` groups the rows by (prime, kind, m, control, N), lists `D`,
//! width and steps per target, the multiply-adds and the F4 time, the
//! verification outcome, and the reference values of the same
//! (kind family, m, control, N) from the recorded rows; it applies the
//! decision rule of PROTOCOL §6 per series and writes `cells.tsv`,
//! `series.md`, `verdicts.md` and `plot_data.csv`.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

use serde_json::Value;

fn read_rows(path: &str) -> Vec<Value> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    text.lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .map(|l| serde_json::from_str::<Value>(l).unwrap_or_else(|e| panic!("{path}: {e}")))
        .filter(|v| v.get("kind").is_some())
        .collect()
}

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    assert!(n > 0, "median of nothing");
    if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

fn u(r: &Value, k: &str) -> u64 {
    r[k].as_u64().unwrap_or_else(|| panic!("field {k}"))
}

fn st(r: &Value, k: &str) -> String {
    match &r[k] {
        Value::String(s) => s.clone(),
        v => v.to_string(),
    }
}

// ── speed ──────────────────────────────────────────────────────────

struct Claim {
    id: &'static str,
    cell: &'static str,
    base: &'static str,
    cand: &'static str,
    threshold: f64,
    what: &'static str,
}

const CLAIMS: [Claim; 5] = [
    Claim {
        id: "S1a",
        cell: "C1",
        base: "u32",
        cand: "fp32",
        threshold: 0.95,
        what: "split Dickson p1, m=2 N=16: generic Fp32 no regression vs f4_fp_tower",
    },
    Claim {
        id: "S1b",
        cell: "C2",
        base: "u32",
        cand: "fp32",
        threshold: 0.95,
        what: "torus Dickson T31, m=2 N=16: generic Fp32 no regression vs f4_fp_tower",
    },
    Claim {
        id: "S2a",
        cell: "C3",
        base: "u32",
        cand: "fp32",
        threshold: 1.5,
        what: "Kummer p2 (> 2^31), m=2 N=16: generic Fp32 >= 1.5x f4_fp_tower",
    },
    Claim {
        id: "S2b",
        cell: "C4",
        base: "u32",
        cand: "fp32",
        threshold: 1.5,
        what: "Kummer p2 (> 2^31), m=3 N=12: generic Fp32 >= 1.5x f4_fp_tower",
    },
    Claim {
        id: "S3",
        cell: "C5",
        base: "fp64",
        cand: "fp52",
        threshold: 1.8,
        what: "Kummer v40 (45-bit), m=2 N=16: Fp52 IFMA >= 1.8x Fp64",
    },
];

fn speed(dir: &str, out: &str) {
    std::fs::create_dir_all(out).unwrap();
    // (cell, arm) -> rep -> row
    let mut rows: BTreeMap<(String, String), BTreeMap<u32, Value>> = BTreeMap::new();
    for e in std::fs::read_dir(dir).unwrap() {
        let path = e.unwrap().path();
        let name = path.file_name().unwrap().to_string_lossy().to_string();
        let Some(stem) = name.strip_suffix(".jsonl") else {
            continue;
        };
        let parts: Vec<&str> = stem.split('.').collect();
        if parts.len() != 3 {
            continue;
        }
        let rep: u32 = parts[2].trim_start_matches("rep").parse().unwrap();
        let mut rs = read_rows(path.to_str().unwrap());
        assert_eq!(rs.len(), 1, "{name}: one row per run");
        rows.entry((parts[0].into(), parts[1].into()))
            .or_default()
            .insert(rep, rs.remove(0));
    }
    // G3: the arms of a cell agree on every field but the clock and labels.
    let strip = |r: &Value| -> Value {
        let mut r = r.clone();
        let o = r.as_object_mut().unwrap();
        for k in [
            "ms",
            "wall_s",
            "engine",
            "backend",
            "builder",
            "sig",
            "sig_order",
            "sig_rewrite",
            "sig_steps",
        ] {
            o.remove(k);
        }
        r
    };
    let mut gate = String::new();
    let cells: BTreeSet<String> = rows.keys().map(|k| k.0.clone()).collect();
    let mut g3_ok = true;
    for c in &cells {
        let mut reference: Option<Value> = None;
        for ((cc, arm), reps) in &rows {
            if cc != c {
                continue;
            }
            for (rep, r) in reps {
                assert_eq!(
                    r["timed_out"],
                    Value::Bool(false),
                    "{c} {arm} rep {rep} stopped"
                );
                let s = strip(r);
                match &reference {
                    None => reference = Some(s),
                    Some(x) => {
                        if *x != s {
                            g3_ok = false;
                            let _ = writeln!(gate, "G3 FAIL: {c} {arm} rep {rep} differs");
                        }
                    }
                }
            }
        }
    }
    let _ = writeln!(
        gate,
        "G3 (arms agree on every field but clock/labels): {}",
        if g3_ok { "PASS" } else { "FAIL" }
    );

    let mut raw = String::from("cell,arm,rep,ms,muladds,max_cols_to_solution,solving_degree_max,steps_to_solution,p,kind,m,N,backend\n");
    for ((c, arm), reps) in &rows {
        for (rep, r) in reps {
            let _ = writeln!(
                raw,
                "{c},{arm},{rep},{:.3},{},{},{},{},{},{},{},{},{}",
                r["ms"].as_f64().unwrap(),
                u(r, "muladds"),
                u(r, "max_cols_to_solution"),
                u(r, "solving_degree_max"),
                u(r, "steps_to_solution"),
                st(r, "p"),
                st(r, "kind"),
                u(r, "m"),
                u(r, "N"),
                r.get("backend").map_or("f4_fp_tower".into(), |b| b
                    .as_str()
                    .unwrap_or("")
                    .to_string()),
            );
        }
    }
    std::fs::write(format!("{out}/raw.csv"), raw).unwrap();

    // Per-arm medians and ns per multiply-add.
    let mut summary = String::from("cell,arm,reps,median_ms,min_ms,max_ms,muladds,ns_per_muladd\n");
    let mut ns_of: BTreeMap<(String, String), f64> = BTreeMap::new();
    let mut tables = String::new();
    let _ = writeln!(tables, "{gate}");
    let _ = writeln!(tables, "| cell | arm | p | kind | m | N | D / width / steps | muladds | median ms | min–max ms | ns / muladd |");
    let _ = writeln!(tables, "|:--|:--|:--|:--|--:|--:|:--|--:|--:|:--|--:|");
    for ((c, arm), reps) in &rows {
        let mut ms: Vec<f64> = reps.values().map(|r| r["ms"].as_f64().unwrap()).collect();
        let (mn, mx) = (
            ms.iter().cloned().fold(f64::INFINITY, f64::min),
            ms.iter().cloned().fold(0.0, f64::max),
        );
        let med = median(&mut ms);
        let r0 = reps.values().next().unwrap();
        let ma = u(r0, "muladds");
        let ns = med * 1e6 / ma as f64;
        ns_of.insert((c.clone(), arm.clone()), ns);
        let _ = writeln!(
            summary,
            "{c},{arm},{},{med:.3},{mn:.3},{mx:.3},{ma},{ns:.4}",
            reps.len()
        );
        let _ = writeln!(
            tables,
            "| {c} | {arm} | {} | {} | {} | {} | {} / {} / {} | {ma} | {med:.0} | {mn:.0}–{mx:.0} | {ns:.3} |",
            st(r0, "p"),
            st(r0, "kind"),
            u(r0, "m"),
            u(r0, "N"),
            u(r0, "solving_degree_max"),
            u(r0, "max_cols_to_solution"),
            u(r0, "steps_to_solution"),
        );
    }
    std::fs::write(format!("{out}/summary.csv"), summary).unwrap();

    let mut claims = String::from("id,claim,cell,baseline,candidate,baseline_median_ms,candidate_median_ms,ratio_median_paired,ratio_min,ratio_max,threshold,verdict\n");
    let _ = writeln!(tables, "\n| id | claim | baseline ms | candidate ms | ratio (median of paired) | paired min–max | threshold | verdict |");
    let _ = writeln!(tables, "|:--|:--|--:|--:|--:|:--|--:|:--|");
    let mut pass = 0;
    for cl in &CLAIMS {
        let b = rows.get(&(cl.cell.into(), cl.base.into()));
        let c = rows.get(&(cl.cell.into(), cl.cand.into()));
        let (Some(b), Some(c)) = (b, c) else {
            let _ = writeln!(
                claims,
                "{},{},{},{},{},,,,,,{},MISSING",
                cl.id, cl.what, cl.cell, cl.base, cl.cand, cl.threshold
            );
            continue;
        };
        let mut ratios: Vec<f64> = b
            .iter()
            .filter_map(|(rep, rb)| {
                c.get(rep)
                    .map(|rc| rb["ms"].as_f64().unwrap() / rc["ms"].as_f64().unwrap())
            })
            .collect();
        let mut bm: Vec<f64> = b.values().map(|r| r["ms"].as_f64().unwrap()).collect();
        let mut cm: Vec<f64> = c.values().map(|r| r["ms"].as_f64().unwrap()).collect();
        let (lo, hi) = (
            ratios.iter().cloned().fold(f64::INFINITY, f64::min),
            ratios.iter().cloned().fold(0.0, f64::max),
        );
        let med = median(&mut ratios);
        let ok = g3_ok && med >= cl.threshold;
        pass += usize::from(ok);
        let verdict = if ok { "PASS" } else { "FAIL" };
        let (bmed, cmed) = (median(&mut bm), median(&mut cm));
        let _ = writeln!(
            claims,
            "{},{},{},{},{},{bmed:.3},{cmed:.3},{med:.4},{lo:.4},{hi:.4},{},{verdict}",
            cl.id, cl.what, cl.cell, cl.base, cl.cand, cl.threshold
        );
        let _ = writeln!(
            tables,
            "| {} | {} | {bmed:.0} | {cmed:.0} | {med:.3} | {lo:.3}–{hi:.3} | {} | {verdict} |",
            cl.id, cl.what, cl.threshold
        );
    }
    let _ = writeln!(
        tables,
        "\n{pass} PASS of {} pre-registered claims.",
        CLAIMS.len()
    );
    // Descriptive: cost per multiply-add relative to generic Fp32 at p1 (C1).
    if let Some(&base) = ns_of.get(&("C1".to_string(), "fp32".to_string())) {
        let _ = writeln!(
            tables,
            "\n| cell | arm | ns / muladd | × generic Fp32 (C1) |\n|:--|:--|--:|--:|"
        );
        for ((c, arm), ns) in &ns_of {
            let _ = writeln!(tables, "| {c} | {arm} | {ns:.3} | {:.2} |", ns / base);
        }
    }
    std::fs::write(format!("{out}/claims.csv"), claims).unwrap();
    std::fs::write(format!("{out}/tables.md"), &tables).unwrap();
    print!("{tables}");
}

// ── cells ──────────────────────────────────────────────────────────

/// The kind family a reference row stands for: the torus and the split
/// Dickson share `D` (round 8), so both are references for a torus cell.
fn family(kind: &str) -> &'static str {
    match kind {
        "kummer" => "kummer",
        "dickson" | "dickson-torus" => "dickson",
        _ => "other",
    }
}

fn label(r: &Value) -> String {
    r.get("prime")
        .and_then(Value::as_str)
        .map(String::from)
        .unwrap_or_else(|| st(r, "p"))
}

type Key = (String, String, u64, String, u64); // prime, kind, m, control, N

fn cells(out: &str, refs: &[String], files: &[String]) {
    std::fs::create_dir_all(out).unwrap();
    // reference: (family, m, control, N) -> (D set, width set, steps set, sources)
    let mut reference: BTreeMap<
        (String, u64, String, u64),
        (
            BTreeSet<u64>,
            BTreeSet<u64>,
            BTreeSet<u64>,
            BTreeSet<String>,
        ),
    > = BTreeMap::new();
    for f in refs {
        for r in read_rows(f) {
            if r["timed_out"] == Value::Bool(true) || st(&r, "target") != "random" {
                continue;
            }
            let k = (
                family(&st(&r, "kind")).to_string(),
                u(&r, "m"),
                st(&r, "control"),
                u(&r, "N"),
            );
            let e = reference.entry(k).or_default();
            e.0.insert(u(&r, "solving_degree_max"));
            e.1.insert(u(&r, "max_cols_to_solution"));
            e.2.insert(u(&r, "steps_to_solution"));
            let src = std::path::Path::new(f)
                .file_stem()
                .unwrap()
                .to_string_lossy()
                .to_string();
            e.3.insert(format!("{src}@{}", st(&r, "p")));
        }
    }
    let mut groups: BTreeMap<Key, Vec<Value>> = BTreeMap::new();
    for f in files {
        for r in read_rows(f) {
            let k = (
                label(&r),
                st(&r, "kind"),
                u(&r, "m"),
                st(&r, "control"),
                u(&r, "N"),
            );
            groups.entry(k).or_default().push(r);
        }
    }
    let mut tsv = String::from("prime\tp\tkind\tm\tcontrol\tN\ttarget\ttarget_index\tD\twidth\tsteps\tmuladds\tms\tbackend\tstopped\tinconsistent\tverify_ok\tplanted_ok\tmax_dense_entries\n");
    let mut plot = String::from(
        "prime,kind,m,control,N,D,width,steps,muladds_median,ms_median,backend,ref_D\n",
    );
    let mut md = String::new();
    // series verdicts: (prime, kind, m, control) -> (equal, differ, unvalued)
    let mut series: BTreeMap<(String, String, u64, String), (u32, u32, u32, Vec<String>)> =
        BTreeMap::new();
    let _ = writeln!(md, "| prime | kind | m | control | N | targets | D | width | steps | muladds (median) | F4 time (median) | verify | reference D / width / steps | D vs reference |");
    let _ = writeln!(
        md,
        "|:--|:--|--:|:--|--:|:--|:--|:--|:--|--:|--:|:--|:--|:--|"
    );
    for ((prime, kind, m, control, n), rs) in &groups {
        let mut ds = BTreeSet::new();
        let mut ws = BTreeSet::new();
        let mut ss = BTreeSet::new();
        let mut mas = Vec::new();
        let mut mss = Vec::new();
        let mut verify = Vec::new();
        let mut kinds_t = Vec::new();
        let mut stopped = 0;
        let mut backend = String::new();
        for r in rs {
            let stop = r["timed_out"] == Value::Bool(true);
            let v = r
                .get("verify")
                .and_then(|v| v.get("ok"))
                .and_then(Value::as_bool);
            let _ = writeln!(
                tsv,
                "{prime}\t{}\t{kind}\t{m}\t{control}\t{n}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{}\t{stop}\t{}\t{}\t{}\t{}",
                st(r, "p"),
                st(r, "target"),
                u(r, "target_index"),
                u(r, "solving_degree_max"),
                u(r, "max_cols_to_solution"),
                u(r, "steps_to_solution"),
                u(r, "muladds"),
                r["ms"].as_f64().unwrap(),
                st(r, "backend"),
                r["inconsistent"],
                v.map_or("-".into(), |b| b.to_string()),
                r.get("planted_ok").map_or("-".into(), |x| x.to_string()),
                r.get("max_dense_entries").map_or("-".into(), |x| x.to_string()),
            );
            backend = st(r, "backend");
            kinds_t.push(st(r, "target"));
            if stop {
                stopped += 1;
                continue;
            }
            ds.insert(u(r, "solving_degree_max"));
            ws.insert(u(r, "max_cols_to_solution"));
            ss.insert(u(r, "steps_to_solution"));
            mas.push(u(r, "muladds") as f64);
            mss.push(r["ms"].as_f64().unwrap());
            if let Some(b) = v {
                verify.push(b);
            }
        }
        let fmt = |s: &BTreeSet<u64>| s.iter().map(u64::to_string).collect::<Vec<_>>().join("/");
        let rk = (family(kind).to_string(), *m, control.clone(), *n);
        let rf = reference.get(&rk);
        let sk = (prime.clone(), kind.clone(), *m, control.clone());
        let e = series.entry(sk).or_default();
        let cmp = if ds.is_empty() {
            e.2 += 1;
            "not valued (stopped)".to_string()
        } else if let Some((rd, _, _, _)) = rf {
            if ds == *rd {
                e.0 += 1;
                "equal".into()
            } else {
                e.1 += 1;
                e.3.push(format!("N={n}: {} vs {}", fmt(&ds), fmt(rd)));
                "DIFFERS".into()
            }
        } else {
            e.2 += 1;
            "no reference".into()
        };
        let all_targets: BTreeSet<&str> = kinds_t.iter().map(String::as_str).collect();
        let vtxt = if verify.is_empty() {
            "-".to_string()
        } else {
            format!(
                "{}/{} ok",
                verify.iter().filter(|&&b| b).count(),
                verify.len()
            )
        };
        let (mam, msm) = if mas.is_empty() {
            (f64::NAN, f64::NAN)
        } else {
            (median(&mut mas), median(&mut mss))
        };
        let reftxt = rf.map_or("—".to_string(), |(d, w, s, src)| {
            format!(
                "{} / {} / {} ({})",
                fmt(d),
                fmt(w),
                fmt(s),
                src.iter().cloned().collect::<Vec<_>>().join(", ")
            )
        });
        let _ = writeln!(
            md,
            "| {prime} | {kind} | {m} | {control} | {n} | {} × {}{} | {} | {} | {} | {:.3e} | {:.1} s | {vtxt} | {reftxt} | {cmp} |",
            rs.len(),
            all_targets.iter().cloned().collect::<Vec<_>>().join("+"),
            if stopped > 0 { format!(" ({stopped} stopped)") } else { String::new() },
            fmt(&ds),
            fmt(&ws),
            fmt(&ss),
            mam,
            msm / 1e3,
        );
        if !ds.is_empty() {
            let _ = writeln!(
                plot,
                "{prime},{kind},{m},{control},{n},{},{},{},{mam},{msm},{backend},{}",
                ds.iter().max().unwrap(),
                ws.iter().max().unwrap(),
                ss.iter().max().unwrap(),
                rf.map_or(String::new(), |(d, _, _, _)| fmt(d)),
            );
        }
    }
    let mut verdicts = String::from("| prime | kind | m | control | cells equal | cells differing | not valued | verdict | differences |\n|:--|:--|--:|:--|--:|--:|--:|:--|:--|\n");
    for ((prime, kind, m, control), (eq, df, nv, diffs)) in &series {
        let verdict = if *df > 0 {
            "H_same refuted"
        } else if *eq > 0 {
            "H_same supported"
        } else {
            "not valued"
        };
        let _ = writeln!(
            verdicts,
            "| {prime} | {kind} | {m} | {control} | {eq} | {df} | {nv} | {verdict} | {} |",
            diffs.join("; ")
        );
    }
    std::fs::write(format!("{out}/cells.tsv"), tsv).unwrap();
    std::fs::write(format!("{out}/series.md"), &md).unwrap();
    std::fs::write(format!("{out}/verdicts.md"), &verdicts).unwrap();
    std::fs::write(format!("{out}/plot_data.csv"), plot).unwrap();
    print!("{md}\n{verdicts}");
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("speed") => speed(&args[1], &args[2]),
        Some("cells") => {
            let out = &args[1];
            let mut refs = Vec::new();
            let mut rows = Vec::new();
            let mut cur: Option<&mut Vec<String>> = None;
            for a in &args[2..] {
                match a.as_str() {
                    "--ref" => cur = Some(&mut refs),
                    "--rows" => cur = Some(&mut rows),
                    _ => cur.as_mut().expect("--ref or --rows first").push(a.clone()),
                }
            }
            cells(out, &refs, &rows);
        }
        _ => {
            eprintln!(
                "usage: f4_generic_report speed DIR OUT | cells OUT --ref FILES… --rows FILES…"
            );
            std::process::exit(2);
        }
    }
}

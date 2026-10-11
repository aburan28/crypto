//! Summarise a benchmark JSONL file (V1/V2 groups) into Markdown tables:
//! `report results/baseline.jsonl > results/baseline.md`. The native replacement of the earlier
//! `scripts/report.py`; it prints the same tables. Records without a `group` are skipped.
use isogeny_algos::json::Json;
use std::collections::HashMap;

fn fmt_ns(x: f64) -> String {
    if x < 1e3 {
        format!("{x:.0} ns")
    } else if x < 1e6 {
        format!("{:.1} us", x / 1e3)
    } else if x < 1e9 {
        format!("{:.2} ms", x / 1e6)
    } else {
        format!("{:.2} s", x / 1e9)
    }
}

/// Python's str() of a parsed value.
fn py_str(v: Option<&Json>) -> String {
    match v {
        None | Some(Json::Null) => "None".into(),
        Some(Json::Str(s)) => s.clone(),
        Some(Json::Raw(r)) => r.clone(),
        Some(Json::Bool(b)) => (if *b { "True" } else { "False" }).into(),
        Some(Json::Num(n)) => n.to_string(),
        Some(Json::Float(x)) => x.to_string(),
        Some(other) => py_dumps(other),
    }
}

/// Python's json.dumps with its default separators and ensure_ascii.
fn py_dumps(v: &Json) -> String {
    fn esc(s: &str) -> String {
        let mut o = String::from("\"");
        for c in s.chars() {
            match c {
                '"' => o.push_str("\\\""),
                '\\' => o.push_str("\\\\"),
                '\n' => o.push_str("\\n"),
                '\r' => o.push_str("\\r"),
                '\t' => o.push_str("\\t"),
                c if (c as u32) < 0x20 || (c as u32) > 0x7e => {
                    let mut buf = [0u16; 2];
                    for u in c.encode_utf16(&mut buf) {
                        o.push_str(&format!("\\u{:04x}", u));
                    }
                }
                c => o.push(c),
            }
        }
        o.push('"');
        o
    }
    match v {
        Json::Null => "null".into(),
        Json::Bool(b) => (if *b { "true" } else { "false" }).into(),
        Json::Num(n) => n.to_string(),
        Json::Float(x) => x.to_string(),
        Json::Raw(r) => r.clone(),
        Json::Str(s) => esc(s),
        Json::Arr(a) => format!(
            "[{}]",
            a.iter().map(py_dumps).collect::<Vec<_>>().join(", ")
        ),
        Json::Obj(kv) => format!(
            "{{{}}}",
            kv.iter()
                .map(|(k, v)| format!("{}: {}", esc(k), py_dumps(v)))
                .collect::<Vec<_>>()
                .join(", ")
        ),
    }
}

/// statistics.median
fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

fn num(r: &Json, k: &str) -> f64 {
    r.get(k).and_then(Json::as_f64).unwrap_or(f64::NAN)
}
fn text(r: &Json, k: &str) -> String {
    py_str(r.get(k))
}
fn int_of(s: &str) -> i64 {
    s.trim().parse::<f64>().map(|x| x as i64).unwrap_or(0)
}
fn verified(r: &Json) -> bool {
    matches!(r.get("verified"), Some(Json::Bool(true)))
}

struct Out(Vec<String>);
impl Out {
    fn table(&mut self, title: &str, header: &[String], rows: &[Vec<String>]) {
        self.0.push(format!("\n### {title}\n"));
        self.0.push(format!("| {} |", header.join(" | ")));
        self.0
            .push(format!("|{}|", vec!["---"; header.len()].join("|")));
        for r in rows {
            self.0.push(format!("| {} |", r.join(" | ")));
        }
    }
    #[allow(clippy::too_many_arguments)]
    fn pivot(&mut self, title: &str, recs: &[&Json], rowkey: &str, colkey: &str, rowname: &str) {
        let mut rows: Vec<String> = vec![];
        for r in recs {
            let v = text(r, rowkey);
            if !rows.contains(&v) {
                rows.push(v);
            }
        }
        rows.sort_by_key(|v| int_of(v));
        let mut cols: Vec<String> = vec![];
        for r in recs {
            let c = text(r, colkey);
            if !cols.contains(&c) {
                cols.push(c);
            }
        }
        let body: Vec<Vec<String>> = rows
            .iter()
            .map(|rw| {
                let mut line = vec![rw.clone()];
                for c in &cols {
                    let m = recs
                        .iter()
                        .find(|r| text(r, rowkey) == *rw && text(r, colkey) == *c);
                    line.push(m.map_or("-".into(), |r| fmt_ns(num(r, "median_ns"))));
                }
                line
            })
            .collect();
        let mut header = vec![rowname.to_string()];
        header.extend(cols);
        self.table(title, &header, &body);
    }
}

fn s(v: &str) -> String {
    v.to_string()
}

fn main() {
    let path = std::env::args().nth(1).unwrap_or_else(|| {
        eprintln!("usage: report <results.jsonl> > <out.md>");
        std::process::exit(2);
    });
    let data = std::fs::read_to_string(&path).unwrap_or_else(|e| {
        eprintln!("{path}: {e}");
        std::process::exit(2);
    });
    let recs: Vec<Json> = data
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| Json::parse(l).unwrap_or_else(|e| panic!("{path}: {e}")))
        .filter(|r| r.get("group").is_some())
        .collect();
    let group = |g: &str| {
        recs.iter()
            .filter(|r| text(r, "group") == g)
            .collect::<Vec<&Json>>()
    };
    let mut out = Out(vec![]);
    let uniq_sorted = |v: Vec<String>| {
        let mut u: Vec<String> = vec![];
        for x in v {
            if !u.contains(&x) {
                u.push(x);
            }
        }
        u.sort_by_key(|x| int_of(x));
        u
    };

    // kernel -> isogeny
    let k = group("kernel");
    if !k.is_empty() {
        let akey = |r: &Json| {
            let t = text(r, "task");
            text(r, "algo")
                + if t == "codomain_given_h" {
                    "(h given)"
                } else if t.contains("poly+formulas") {
                    "(poly+formulas)"
                } else {
                    ""
                }
        };
        for bits in uniq_sorted(k.iter().map(|r| text(r, "p_bits")).collect()) {
            for task in ["codomain_from_point", "codomain_given_h", "eval_point"] {
                let sel: Vec<&Json> = k
                    .iter()
                    .copied()
                    .filter(|r| {
                        text(r, "p_bits") == bits
                            && (text(r, "task") == task
                                || (task == "codomain_from_point"
                                    && text(r, "task").starts_with("codomain_from_point")))
                    })
                    .collect();
                let mut algos: Vec<String> = vec![];
                for r in &sel {
                    let key = akey(r);
                    if !algos.contains(&key) {
                        algos.push(key);
                    }
                }
                let mut ells: Vec<i64> = sel.iter().map(|r| int_of(&text(r, "ell"))).collect();
                ells.sort();
                ells.dedup();
                let rows: Vec<Vec<String>> = ells
                    .iter()
                    .map(|&l| {
                        let mut row = vec![l.to_string()];
                        for a in &algos {
                            let m = sel
                                .iter()
                                .find(|r| int_of(&text(r, "ell")) == l && akey(r) == *a);
                            row.push(m.map_or("-".into(), |r| fmt_ns(num(r, "median_ns"))));
                        }
                        row
                    })
                    .collect();
                if !rows.is_empty() {
                    let mut header = vec![s("l")];
                    header.extend(algos.clone());
                    out.table(
                        &format!("Kernel->isogeny, {bits}-bit p, task `{task}` (median)"),
                        &header,
                        &rows,
                    );
                }
            }
        }
    }

    // find
    let f = group("find");
    if !f.is_empty() {
        for bits in uniq_sorted(f.iter().map(|r| text(r, "p_bits")).collect()) {
            let sel: Vec<&Json> = f
                .iter()
                .copied()
                .filter(|r| text(r, "p_bits") == bits)
                .collect();
            let mut ells: Vec<i64> = sel.iter().map(|r| int_of(&text(r, "ell"))).collect();
            ells.sort();
            ells.dedup();
            let rows: Vec<Vec<String>> = ells
                .iter()
                .map(|&l| {
                    let mut d: HashMap<String, &Json> = HashMap::new();
                    for r in sel.iter().filter(|r| int_of(&text(r, "ell")) == l) {
                        d.insert(text(r, "algo"), r);
                    }
                    let g = |a: &str| d.get(a).map_or("-".into(), |r| fmt_ns(num(r, "median_ns")));
                    let kf = d
                        .get("divpoly_factor")
                        .and_then(|r| r.get("kernels_found"))
                        .map_or("-".into(), |v| py_str(Some(v)));
                    vec![
                        l.to_string(),
                        kf,
                        g("divpoly_factor"),
                        g("elkies_phi_bmss"),
                        g("phi_roots_only"),
                        g("phi_setup"),
                    ]
                })
                .collect();
            let header: Vec<String> = [
                "l",
                "#kernels",
                "divpoly factor",
                "Phi roots+Elkies+BMSS",
                "Phi roots only",
                "Phi_l setup (once)",
            ]
            .iter()
            .map(|x| s(x))
            .collect();
            out.table(
                &format!("Kernel finding, {bits}-bit p (median; Phi setup is a single run)"),
                &header,
                &rows,
            );
        }
    }

    // path
    let p = group("path");
    if !p.is_empty() {
        let mut keys: Vec<(String, String)> = vec![];
        for r in &p {
            let kk = (text(r, "algo"), text(r, "p_bits"));
            if !keys.contains(&kk) {
                keys.push(kk);
            }
        }
        keys.sort_by_key(|a| (a.0.clone(), int_of(&a.1)));
        let rows: Vec<Vec<String>> = keys
            .iter()
            .map(|(algo, bits)| {
                let rs: Vec<&Json> = p
                    .iter()
                    .copied()
                    .filter(|r| text(r, "algo") == *algo && text(r, "p_bits") == *bits)
                    .collect();
                let okr: Vec<&Json> = rs.iter().copied().filter(|r| verified(r)).collect();
                let med = if okr.is_empty() {
                    s("-")
                } else {
                    fmt_ns(median(okr.iter().map(|r| num(r, "ns")).collect()))
                };
                let mm = |key: &str| {
                    let v: Vec<f64> = okr
                        .iter()
                        .filter_map(|r| r.get(key).and_then(Json::as_f64))
                        .collect();
                    if v.is_empty() {
                        s("-")
                    } else {
                        format!("{:.0}", median(v))
                    }
                };
                let work: Vec<String> = [
                    "nodes_expanded",
                    "steps",
                    "nodes",
                    "walk_steps",
                    "bfs_nodes",
                ]
                .iter()
                .filter(|k| okr.iter().any(|r| r.get(k).is_some()))
                .map(|k| format!("{k}={}", mm(k)))
                .collect();
                vec![
                    algo.clone(),
                    bits.clone(),
                    format!("{}/{}", okr.len(), rs.len()),
                    med,
                    mm("path_len"),
                    work.join(" "),
                ]
            })
            .collect();
        let header: Vec<String> = [
            "algorithm",
            "p bits",
            "verified",
            "median time",
            "path len",
            "work counters (median)",
        ]
        .iter()
        .map(|x| s(x))
        .collect();
        out.table(
            "Path finding (median over verified instances)",
            &header,
            &rows,
        );
        let bad: Vec<&Json> = p.iter().copied().filter(|r| !verified(r)).collect();
        if !bad.is_empty() {
            out.0.push(s(
                "\nUnverified / failed instances (kept in the raw data):\n",
            ));
            for r in bad {
                out.0.push(format!(
                    "- {} p_bits={} instance={} ns={} note: {}",
                    text(r, "algo"),
                    text(r, "p_bits"),
                    text(r, "instance"),
                    text(r, "ns"),
                    text(r, "note")
                ));
            }
        }
    }

    // V2 groups
    let k2 = group("kernel2");
    if !k2.is_empty() {
        out.pivot("V2 kernel -> isogeny, odd degree, 32-bit p (median): general Velu from P, x-only Velu from x(P), Kohel with h given", &k2, "ell", "algo", "l");
    }
    let k2e = group("kernel2_even");
    if !k2e.is_empty() {
        out.pivot(
            "V2 kernel -> isogeny, even / composite cyclic degree n, 32-bit p (median)",
            &k2e,
            "degree",
            "algo",
            "n",
        );
    }
    let k2m = group("kernel2_mont");
    if !k2m.is_empty() {
        out.pivot(
            "V2 Montgomery x-only Velu vs Weierstrass Velu for the same CSIDH kernel (median)",
            &k2m,
            "ell",
            "algo",
            "l",
        );
    }
    let mut ch = group("chain");
    if !ch.is_empty() {
        ch.sort_by(|a, b| {
            let key = |r: &Json| {
                (
                    int_of(&text(r, "ell")),
                    int_of(&text(r, "e")),
                    text(r, "p_bits"),
                    text(r, "strategy"),
                )
            };
            key(a).cmp(&key(b))
        });
        let rows: Vec<Vec<String>> = ch
            .iter()
            .map(|r| {
                vec![
                    text(r, "ell"),
                    text(r, "e"),
                    text(r, "p_bits"),
                    text(r, "strategy"),
                    fmt_ns(num(r, "median_ns")),
                    (num(r, "l_mults") as i64).to_string(),
                    (num(r, "evals") as i64).to_string(),
                    (num(r, "builds") as i64).to_string(),
                ]
            })
            .collect();
        let header: Vec<String> = [
            "l",
            "e",
            "p bits",
            "strategy",
            "time",
            "l-mults",
            "point evals",
            "isogeny builds",
        ]
        .iter()
        .map(|x| s(x))
        .collect();
        out.table(
            "V2 l^e-kernel chains over F_p^2 (median); identical codomain for all strategies",
            &header,
            &rows,
        );
    }
    let bmss = group("bmss");
    let bm: Vec<&Json> = bmss
        .iter()
        .copied()
        .filter(|r| {
            !["end_to_end_via_phi", "sigma_from_phi", "dual_isogeny"]
                .contains(&text(r, "algo").as_str())
        })
        .collect();
    if !bm.is_empty() {
        out.pivot(
            "V2 (E, E~) -> isogeny, 40-bit p (median; sigma supplied where required)",
            &bm,
            "ell",
            "algo",
            "l",
        );
    }
    let ee: Vec<&Json> = bmss
        .iter()
        .copied()
        .filter(|r| ["sigma_from_phi", "dual_isogeny"].contains(&text(r, "algo").as_str()))
        .collect();
    if !ee.is_empty() {
        out.pivot(
            "V2 sigma from Phi and dual isogeny (median)",
            &ee,
            "ell",
            "algo",
            "l",
        );
    }
    let e2e: Vec<Json> = bmss
        .iter()
        .filter(|r| text(r, "algo") == "end_to_end_via_phi")
        .map(|r| {
            let mut r = (*r).clone();
            if let Json::Obj(kv) = &mut r {
                let m = format!(
                    "e2e_{}",
                    py_str(kv.iter().find(|(k, _)| k == "method").map(|(_, v)| v))
                );
                kv.push((s("algo2"), Json::Str(m)));
            }
            r
        })
        .collect();
    if !e2e.is_empty() {
        let refs: Vec<&Json> = e2e.iter().collect();
        out.pivot(
            "V2 end-to-end from (E, l) with Phi_l given (median)",
            &refs,
            "ell",
            "algo2",
            "l",
        );
    }
    let cs = group("csidh");
    if !cs.is_empty() {
        let rows: Vec<Vec<String>> = cs
            .iter()
            .map(|r| {
                let Json::Obj(kv) = r else { unreachable!() };
                let extra: Vec<(String, Json)> = kv
                    .iter()
                    .filter(|(k, _)| {
                        ![
                            "group",
                            "algo",
                            "verified",
                            "note",
                            "median_ns",
                            "min_ns",
                            "reps",
                        ]
                        .contains(&k.as_str())
                    })
                    .cloned()
                    .collect();
                let t = r
                    .get("median_ns")
                    .or_else(|| r.get("ns"))
                    .and_then(Json::as_f64);
                vec![
                    text(r, "algo"),
                    py_dumps(&Json::Obj(extra)),
                    t.filter(|&x| x != 0.0).map_or(s("-"), fmt_ns),
                    s(if verified(r) { "yes" } else { "NO" }),
                ]
            })
            .collect();
        let header: Vec<String> = ["what", "parameters / counters", "time", "verified"]
            .iter()
            .map(|x| s(x))
            .collect();
        out.table("V2 CSIDH-style action", &header, &rows);
    }
    let en = group("endo");
    if !en.is_empty() {
        let rows: Vec<Vec<String>> = en
            .iter()
            .map(|r| {
                vec![
                    text(r, "p_bits"),
                    text(r, "instance"),
                    text(r, "height_3"),
                    text(r, "level_3"),
                    fmt_ns(num(r, "median_ns")),
                    s(if verified(r) { "yes" } else { "NO" }),
                ]
            })
            .collect();
        let header: Vec<String> = [
            "p bits",
            "instance",
            "height of 3-volcano",
            "level of E",
            "time",
            "verified",
        ]
        .iter()
        .map(|x| s(x))
        .collect();
        out.table("V2 Kohel End(E) conductor (median)", &header, &rows);
    }
    let wk = group("walks");
    if !wk.is_empty() {
        let mut keys: Vec<(String, String)> = vec![];
        for r in &wk {
            let kk = (text(r, "algo"), text(r, "p_bits"));
            if !keys.contains(&kk) {
                keys.push(kk);
            }
        }
        keys.sort_by_key(|a| (int_of(&a.1), a.0.clone()));
        let rows: Vec<Vec<String>> = keys
            .iter()
            .map(|(algo, bits)| {
                let rs: Vec<&Json> = wk
                    .iter()
                    .copied()
                    .filter(|r| text(r, "algo") == *algo && text(r, "p_bits") == *bits)
                    .collect();
                let ok: Vec<&Json> = rs.iter().copied().filter(|r| verified(r)).collect();
                let med = if ok.is_empty() {
                    s("-")
                } else {
                    fmt_ns(median(ok.iter().map(|r| num(r, "ns")).collect()))
                };
                let steps: Vec<f64> = ok
                    .iter()
                    .filter_map(|r| {
                        r.get("steps")
                            .or_else(|| r.get("nodes_expanded"))
                            .and_then(Json::as_f64)
                    })
                    .collect();
                vec![
                    algo.clone(),
                    bits.clone(),
                    format!("{}/{}", ok.len(), rs.len()),
                    med,
                    if steps.is_empty() {
                        s("-")
                    } else {
                        format!("{:.0}", median(steps))
                    },
                ]
            })
            .collect();
        let header: Vec<String> = [
            "algorithm",
            "p bits",
            "verified",
            "median time",
            "median steps / nodes",
        ]
        .iter()
        .map(|x| s(x))
        .collect();
        out.table(
            "V2 walks on identical instances (median over verified instances)",
            &header,
            &rows,
        );
    }
    println!("{}", out.0.join("\n"));
}

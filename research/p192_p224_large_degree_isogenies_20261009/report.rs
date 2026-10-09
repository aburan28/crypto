//! Deterministic report, vector visual, native PDF, and append-only model registration.
use num_bigint::BigUint;
use serde_json::{json, value::RawValue, Value};
use std::{fmt::Write, fs, path::Path};
#[allow(dead_code)]
#[path = "../../src/hash/sha256.rs"]
mod hash_sha256;
fn hash(b: &[u8]) -> String {
    hash_sha256::sha256(b)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}
fn load(p: &Path) -> Value {
    serde_json::from_slice(&fs::read(p).unwrap()).unwrap()
}
fn number(s: &str) -> Value {
    serde_json::from_str(s).unwrap()
}
fn big(v: &Value) -> BigUint {
    BigUint::parse_bytes(v.as_str().unwrap().as_bytes(), 10).unwrap()
}
fn pretty(v: &Value) -> String {
    let mut bytes = Vec::new();
    let mut ser = serde_json::Serializer::with_formatter(
        &mut bytes,
        serde_json::ser::PrettyFormatter::with_indent(b" "),
    );
    serde::Serialize::serialize(v, &mut ser).unwrap();
    String::from_utf8(bytes).unwrap()
}
fn register(registry: &Path, curves: &[Value], source: &str) {
    let text = fs::read_to_string(registry).unwrap();
    let start = text.find("\"curves\": [").unwrap() + "\"curves\": ".len();
    let mut stream =
        serde_json::Deserializer::from_str(&text[start..]).into_iter::<Vec<Box<RawValue>>>();
    let mut entries = stream.next().unwrap().unwrap();
    let end = start + stream.byte_offset();
    let mut held: Vec<String> = entries
        .iter()
        .map(|r| {
            serde_json::from_str::<Value>(r.get()).unwrap()["slug"]
                .as_str()
                .unwrap()
                .to_owned()
        })
        .collect();
    for c in curves {
        let slug = c["name"].as_str().unwrap();
        if held.iter().any(|s| s == slug) {
            let existing: Value = serde_json::from_str(
                entries
                    .iter()
                    .find(|r| serde_json::from_str::<Value>(r.get()).unwrap()["slug"] == slug)
                    .unwrap()
                    .get(),
            )
            .unwrap();
            assert!(existing["representations"].as_array().unwrap().iter().any(|r|r["ec1"]==c["target_ec1"] && r["curve_uid"]==c["target_curve_uid"]),
                "existing target model needs its mapped-generator representation registered");
            continue;
        }
        let p = big(&c["p"]);
        let n = big(&c["group_order"]);
        let trace = &p + BigUint::from(1u8) - n;
        let field = json!({"characteristic":number(&p.to_string()),"degree":1,"representation":"prime","element_encoding":"hex integer modulo the characteristic"});
        let representation = json!({"model":"short Weierstrass","coefficients":[0,0,0,number(c["a"].as_str().unwrap()),number(c["b"].as_str().unwrap())],
            "subgroup_order":c["subgroup_order"],"cofactor":1,"generator":c["generator"],"target_group":"prime-order subgroup"});
        let e = json!({"slug":slug,"icv1":c["icv1"],"family":"prime","params":{"p":c["p"],"a":c["a"],"b":c["b"]},
            "model_json":c["model_json"],"trace":number(&trace.to_string()),"order":c["group_order"],"j":c["j"],"end":"unk",
            "aliases":[slug],"standard_names":[],"sources":[source],"representations":[{"ec1":c["target_ec1"],"curve_uid":c["target_curve_uid"],
            "field_sha256":hash(serde_json::to_string(&field).unwrap().as_bytes()),"field":field,"curve":representation,"sources":[source]}]});
        let item = pretty(&e)
            .lines()
            .enumerate()
            .map(|(i, l)| {
                if i == 0 {
                    l.to_owned()
                } else {
                    format!("  {l}")
                }
            })
            .collect::<Vec<_>>()
            .join("\n");
        entries.push(RawValue::from_string(item).unwrap());
        held.push(slug.to_owned());
    }
    // Preserve the existing registry's order and raw records. Append only the
    // newly replayed models; sorting here would reorder unrelated standards.
    let output = format!(
        "{}[\n  {}\n ]{}",
        &text[..start],
        entries
            .iter()
            .map(|e| e.get())
            .collect::<Vec<_>>()
            .join(",\n  "),
        &text[end..]
    );
    let _: Value = serde_json::from_str(&output).unwrap();
    fs::write(registry, output).unwrap();
}

#[derive(Clone)]
enum Op {
    Text(f64, f64, f64, String),
    Rect(f64, f64, f64, f64, String),
    Line(f64, f64, f64, f64, String),
    Dot(f64, f64, f64, String),
}
fn xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}
fn svg(ops: &[Op]) -> String {
    let mut s=String::from("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1100\" height=\"760\" viewBox=\"0 0 1100 760\"><rect width=\"1100\" height=\"760\" fill=\"white\"/><g font-family=\"Arial,Helvetica,sans-serif\" fill=\"#172637\">");
    for op in ops {
        match op {
        Op::Text(x,y,f,t)=>write!(s,"<text x=\"{x}\" y=\"{y}\" font-size=\"{f}\">{}</text>",xml(t)).unwrap(),
        Op::Rect(x,y,w,h,c)=>write!(s,"<rect x=\"{x}\" y=\"{y}\" width=\"{w}\" height=\"{h}\" rx=\"9\" fill=\"{c}\"/>").unwrap(),
        Op::Line(x,y,a,b,c)=>write!(s,"<line x1=\"{x}\" y1=\"{y}\" x2=\"{a}\" y2=\"{b}\" stroke=\"{c}\" stroke-width=\"2\"/>").unwrap(),
        Op::Dot(x,y,r,c)=>write!(s,"<circle cx=\"{x}\" cy=\"{y}\" r=\"{r}\" fill=\"{c}\"/>").unwrap(),
    }
    }
    s.push_str("</g></svg>\n");
    s
}
fn escape(s: &str) -> String {
    s.replace('\\', "\\\\")
        .replace('(', "\\(")
        .replace(')', "\\)")
}
fn text(s: &mut String, x: f64, y: f64, size: f64, t: &str) {
    writeln!(
        s,
        "0.08 0.14 0.2 rg BT /F1 {size} Tf {x} {y} Td ({}) Tj ET",
        escape(t)
    )
    .unwrap();
}
fn wrapped(s: &mut String, y: &mut f64, t: &str, size: f64) {
    let mut line = String::new();
    for word in t.split_whitespace() {
        if line.len() + word.len() > 96 {
            text(s, 42., *y, size, &line);
            *y -= 16.;
            line.clear();
        }
        if !line.is_empty() {
            line.push(' ');
        }
        line.push_str(word);
    }
    if !line.is_empty() {
        text(s, 42., *y, size, &line);
        *y -= 16.;
    }
    *y -= 8.;
}
fn rgb(c: &str) -> (f64, f64, f64) {
    let n = u32::from_str_radix(&c[1..], 16).unwrap();
    (
        ((n >> 16) & 255) as f64 / 255.,
        ((n >> 8) & 255) as f64 / 255.,
        (n & 255) as f64 / 255.,
    )
}

fn emit_pdf(path: &Path, pages: Vec<(u32, u32, String)>) {
    let mut objects: Vec<String> = vec![
        String::new(),
        "<< /Type /Catalog /Pages 2 0 R >>".into(),
        String::new(),
        "<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>".into(),
    ];
    let mut kids = vec![];
    for (w, h, body) in pages {
        let n = objects.len();
        kids.push(format!("{n} 0 R"));
        objects.push(format!("<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {w} {h}] /Resources << /Font << /F1 3 0 R >> >> /Contents {} 0 R >>",n+1));
        objects.push(format!(
            "<< /Length {} >>\nstream\n{}endstream",
            body.len(),
            body
        ));
    }
    objects[2] = format!(
        "<< /Type /Pages /Count {} /Kids [{}] >>",
        kids.len(),
        kids.join(" ")
    );
    let mut bytes = String::from("%PDF-1.4\n");
    let mut offsets = vec![0];
    for (i, o) in objects.iter().enumerate().skip(1) {
        offsets.push(bytes.len());
        writeln!(bytes, "{i} 0 obj\n{o}\nendobj").unwrap();
    }
    let xref = bytes.len();
    writeln!(bytes, "xref\n0 {}\n0000000000 65535 f ", objects.len()).unwrap();
    for offset in offsets.iter().skip(1) {
        writeln!(bytes, "{offset:010} 00000 n ").unwrap();
    }
    writeln!(
        bytes,
        "trailer\n<< /Size {} /Root 1 0 R >>\nstartxref\n{xref}\n%%EOF",
        objects.len()
    )
    .unwrap();
    fs::write(path, bytes).unwrap();
}
fn vector_pdf(ops: &[Op]) -> String {
    let mut s = String::from("q 0.71 0 0 -0.71 26 568 cm\n");
    for op in ops {
        match op {
            Op::Text(x, y, size, t) => writeln!(
                s,
                "0.08 0.14 0.2 rg BT /F1 {size} Tf 1 0 0 -1 {x} {y} Tm ({}) Tj ET",
                escape(t)
            )
            .unwrap(),
            Op::Rect(x, y, w, h, c) => {
                let (r, g, b) = rgb(c);
                writeln!(s, "{r} {g} {b} rg {x} {y} {w} {h} re f").unwrap();
            }
            Op::Line(x, y, a, b, c) => {
                let (r, g, bb) = rgb(c);
                writeln!(s, "{r} {g} {bb} RG 2 w {x} {y} m {a} {b} l S").unwrap();
            }
            Op::Dot(x, y, r, c) => {
                let (rr, g, b) = rgb(c);
                writeln!(
                    s,
                    "{rr} {g} {b} rg {} {} {} {} re f",
                    x - r,
                    y - r,
                    2. * r,
                    2. * r
                )
                .unwrap();
            }
        }
    }
    s.push_str("Q\n");
    s
}
fn main() {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() == 3 && args[2] == "--check-manifest" {
        let study = Path::new(&args[0]);
        let m = load(&study.join("MANIFEST.json"));
        let mut count = 0;
        for group in ["files", "sources", "canonical_outputs", "executables"] {
            for r in m[group].as_array().unwrap() {
                let p = Path::new(r["path"].as_str().unwrap());
                let p = if group == "files" {
                    study.join(p)
                } else {
                    p.to_path_buf()
                };
                let b = fs::read(&p).unwrap();
                assert_eq!(
                    b.len() as u64,
                    r["bytes"].as_u64().unwrap(),
                    "{}",
                    p.display()
                );
                assert_eq!(hash(&b), r["sha256"], "{}", p.display());
                count += 1;
            }
        }
        println!("Continuation manifest verified: {count} bound files.");
        return;
    }
    assert_eq!(args.len(), 2);
    let study = Path::new(&args[0]);
    let mut accepted = vec![];
    let mut models = vec![];
    for dir in ["kernel-p192-10453", "kernel-p224-1471"] {
        let replay = load(&study.join(dir).join("replay.json"));
        assert_eq!(replay["status"], "PASS");
        assert_eq!(replay["degree_coverage"], "PARTIAL");
        let records = replay["records"].as_array().unwrap();
        assert_eq!(records.len(), 1);
        let record = records[0].clone();
        assert_eq!(record["kernel_check"], "PASS");
        assert_eq!(record["exact_rational_map_check"], "PASS");
        assert_eq!(record["public_scalar_transport_checks"], 20);
        accepted.push(record);
        models.extend(
            load(&study.join(dir).join("curves.json"))["curves"]
                .as_array()
                .unwrap()
                .iter()
                .cloned(),
        );
    }
    fs::write(study.join("VERIFIED_MODELS.json"),format!("{}\n",serde_json::to_string_pretty(&json!({"schema":"verified-single-eigenline-models/v1","degree_coverage":"PARTIAL","records":accepted,"curves":models})).unwrap())).unwrap();
    register(
        Path::new(&args[1]),
        &models,
        "research/p192_p224_large_degree_isogenies_20261009/VERIFIED_MODELS.json",
    );
    let mut attempts = vec![];
    for dir in [
        "run-p192-509",
        "run-p224-521",
        "run-p192-1021",
        "run-p224-1031",
        "kernel-p224-1471",
        "kernel-p192-10453",
    ] {
        let summary = load(&study.join(dir).join("search.json"));
        let mut row = summary["attempts"][0].clone();
        row["directory"] = json!(dir);
        row["budget_seconds"] = summary["timeout_seconds"].clone();
        row["method"] = json!(if dir.starts_with("kernel") {
            "kernel-first, one eigenline"
        } else {
            "full Hecke/Newton"
        });
        attempts.push(row);
    }
    let mut md=String::from("# Verified high-prime isogeny continuation\n\n2026-10-09. Continuing the requested P-192/P-224 search produced explicit, independently verified rational maps beyond degree 1009 on both sources. Each new map covers one selected Frobenius eigenline; complete enumeration of the two eigenlines remains unresolved. The previous 22 maps and all earlier failures remain preserved.\n\n| Source | Degree | Extension degree | Kernel degree | Numerator / denominator coefficients | Verified maps | Degree coverage |\n|---|---:|---:|---:|---:|---:|---|\n");
    for (i, r) in accepted.iter().enumerate() {
        writeln!(
            md,
            "| {} | {} | {} | {} | {} / {} | 1 | PARTIAL, one of two eigenlines |",
            if i == 0 { "P-192" } else { "P-224" },
            r["ell"],
            if i == 0 { 3 } else { 5 },
            r["kernel_degree"],
            r["map_numerator_coefficients"],
            r["map_denominator_coefficients"]
        )
        .unwrap();
    }
    md.push_str("\nThe [kernel-first protocol](KERNEL_PROTOCOL.md) fixes the source traces, seed, extension validation, eigenlines, budgets and unchanged mathematical map checks. Native finite extensions pass Rabin irreducibility checks. Nonzero prime-degree torsion and the selected Frobenius eigenvalue are checked, and every kernel coefficient must descend to the source prime field. The existing native Kohel API constructs the map. Independent replay proves squarefreeness, division-polynomial torsion, subgroup closure and the Velu codomain, checks exact rational-map substitution, and verifies 20 fresh public scalar transports plus the mapped generator's published subgroup order.\n\n");
    md.push_str("![Verified maps and completed probe outcomes](SEARCH.svg)\n\n## Preserved full-Hecke outcomes\n\n| Source | Degree | Method | Budget seconds | Outcome | Receipt |\n|---|---:|---|---:|---|---|\n");
    for r in &attempts {
        let name = r["curve"].as_str().unwrap();
        let ell = r["ell"].as_u64().unwrap();
        writeln!(
            md,
            "| {} | {} | {} | {} | {} | [receipt]({}/{}/ell-{}.receipt.json) |",
            name,
            ell,
            r["method"].as_str().unwrap(),
            r["budget_seconds"],
            r["receipt"]["status"].as_str().unwrap(),
            r["directory"].as_str().unwrap(),
            name,
            ell
        )
        .unwrap();
    }
    md.push_str("\nAll four full-Hecke jobs retained their two-map acceptance rule and completed with bounded timeouts. The separately preregistered kernel-first route certifies individual maps and records degree coverage PARTIAL; it makes no full-enumeration claim. Process elapsed values and RSS samples are operational receipts on a shared host. No matched speed ratio or ECDLP-cost improvement was measured.\n\n## Implementation and validation\n\nThe modular-series setup now uses a field-valued divisor sieve and avoids u64 cube/sum overflow. The standalone release suite passed 114 tests; the supervisor passed both structural tests; extension arithmetic passed four tests including a separately enumerated point-count fixture. Independent verification uses the original walker field and exact kernel conditions. High-degree polynomial multiplication, reciprocal reduction, memoized division recurrences and block modular composition are checked against the unchanged schoolbook/Horner/all-index routines. All 12 verifier tests pass, and the known-answer single-map control passes the full replay. No mathematical verification is omitted.\n\nP-224 degree-1471 replay passed with the original all-index kernel routine. The initial P-192 eager replay was stopped before it produced a certificate; its construction bytes and logs remain preserved. The validated memoized verifier then completed degree-10453 replay with every original kernel and exact-map check.\n\nSources and reproduction are in [COMMANDS.md](COMMANDS.md) and [execution notes](EXECUTION_NOTES.md). [P-192 replay](kernel-p192-10453/replay.json), [P-224 replay](kernel-p224-1471/replay.json), exact registered models and mapped generators in each run's curves.json, and [MANIFEST.json](MANIFEST.json) bind the evidence. The [PDF](SEARCH.pdf) and vector share editable [report.rs](report.rs) source. Native catalogue/aliases/rosters/browser refresh accompanies these two models; existing IC measurement bytes are preserved. The root-library and stale traits-catalogue publication gates are recorded in the execution notes and are not promoted as passing checks.\n");
    let mut ops=vec![Op::Text(30.,40.,25.,"Verified prime-degree isogenies beyond 1009".into()),Op::Text(30.,72.,15.,"One independently verified Frobenius eigenline per source; second eigenline unresolved.".into())];
    for (i, r) in accepted.iter().enumerate() {
        let y = 105. + i as f64 * 190.;
        for (side, prefix, x, color) in [
            ("Source", "source", 30., "#edf4fb"),
            ("Verified target", "target", 620., "#e7f4ee"),
        ] {
            ops.push(Op::Rect(
                x,
                y,
                if x == 30. { 465. } else { 450. },
                145.,
                color.into(),
            ));
            ops.push(Op::Text(
                x + 15.,
                y + 25.,
                18.,
                format!("{}: {}", if i == 0 { "P-192" } else { "P-224" }, side),
            ));
            for (key, yy, size) in [("icv1", 49., 11.), ("ec1", 70., 12.)] {
                ops.push(Op::Text(
                    x + 15.,
                    y + yy,
                    size,
                    r[format!("{prefix}_{key}")].as_str().unwrap().into(),
                ));
            }
            let uid = r[format!("{prefix}_curve_uid")].as_str().unwrap();
            ops.push(Op::Text(x + 15., y + 93., 10., uid[..53].into()));
            ops.push(Op::Text(x + 15., y + 108., 10., uid[53..].into()));
            ops.push(Op::Text(
                x + 15.,
                y + 130.,
                12.,
                "Kernel, exact map and subgroup transport: PASS".into(),
            ));
        }
        ops.push(Op::Line(495., y + 77., 620., y + 77., "#136e4e".into()));
        ops.push(Op::Line(620., y + 77., 608., y + 71., "#136e4e".into()));
        ops.push(Op::Line(620., y + 77., 608., y + 83., "#136e4e".into()));
        ops.push(Op::Text(
            501.,
            y + 60.,
            12.,
            format!("degree {} ->", r["ell"]),
        ));
    }
    ops.push(Op::Text(
        30.,
        510.,
        20.,
        "Completed explicit probes (logarithmic prime-degree axis)".into(),
    ));
    ops.push(Op::Text(30.,535.,13.,"Population: six preregistered construction attempts. Green: verified map. Orange: timeout.".into()));
    let xpos =
        |n: u64| 90. + ((n as f64).ln() - 509f64.ln()) / (10453f64.ln() - 509f64.ln()) * 940.;
    for (i, name) in ["p192", "p224"].iter().enumerate() {
        let y = 580. + i as f64 * 68.;
        ops.push(Op::Text(30., y + 5., 14., name.to_uppercase()));
        for r in attempts.iter().filter(|r| r["curve"] == *name) {
            let ell = r["ell"].as_u64().unwrap();
            ops.push(Op::Dot(
                xpos(ell),
                y,
                5.,
                if r["receipt"]["status"] == "PASS" {
                    "#136e4e"
                } else {
                    "#cc7c14"
                }
                .into(),
            ));
        }
    }
    for n in [509, 1021, 1471, 4093, 10453] {
        ops.push(Op::Line(xpos(n), 668., xpos(n), 675., "#333333".into()));
        ops.push(Op::Text(xpos(n) - 14., 694., 12., n.to_string()));
    }
    ops.push(Op::Text(30.,730.,13.,"Verified individual maps retain PARTIAL degree coverage; the other eigenline is not drawn.".into()));
    let mut p1 = String::new();
    text(&mut p1, 42., 790., 22., "High-prime isogeny continuation");
    text(
        &mut p1,
        42.,
        763.,
        13.,
        "P-192 / P-224 | 2026-10-09 | independently replayed maps",
    );
    let mut y = 728.;
    wrapped(&mut p1,&mut y,"Result: P-192 degree 10453 and P-224 degree 1471. Each is an explicit rational map over its source prime field, verified independently by the existing walker arithmetic. One of the two eigenlines is constructed at each degree; full enumeration remains unresolved.",11.);
    for (i, r) in accepted.iter().enumerate() {
        wrapped(&mut p1,&mut y,&format!("{}: kernel degree {}; numerator {} coefficients, denominator {} coefficients. Kernel, codomain, exact rational-map substitution and 20 fresh public scalar transports pass; the mapped nonidentity generator satisfies the published prime subgroup order.",if i==0{"P-192"}else{"P-224"},r["kernel_degree"],r["map_numerator_coefficients"],r["map_denominator_coefficients"]),11.);
    }
    wrapped(&mut p1,&mut y,"Method: choose a recorded low-order Frobenius eigenline, construct prime-degree torsion over a validated degree-three or degree-five finite extension, verify its order/eigenvalue, descend its cyclic kernel to the base field, and use the native Kohel map API. All extension and map checks remain explicit.",11.);
    wrapped(&mut p1,&mut y,"Preserved full-Hecke attempts: P-192 degrees 509 and 1021, P-224 degrees 521 and 1031. They timed out under the preregistered 600- or 1800-second budgets. They still required two maps. Kernel-first records separately require one map and retain PARTIAL coverage.",11.);
    wrapped(&mut p1,&mut y,"Bounds: one construction at a time; seed 1; sampled RSS cap 8 GiB, 25 ms samples, at most 100 ms transient monitor-unavailability grace. Kernel cases each had 1800 seconds. Elapsed receipts are operational records on a shared host, without a matched speed or ECDLP-cost comparison.",11.);
    wrapped(&mut p1,&mut y,"Evidence: KERNEL_PROTOCOL.md, PROTOCOL.md, COMMANDS.md, kernel-p192-10453/replay.json, kernel-p224-1471/replay.json and each run's raw CLI record and receipt. The prior 2026-10-08 report and all failures remain preserved.",11.);
    let mut p3 = String::new();
    text(&mut p3, 42., 790., 22., "Validation and provenance");
    let mut y = 748.;
    for t in ["114 standalone release tests pass. The supervisor's two structural/order tests and all four extension tests pass. Reducible moduli are rejected; source-sized extension inverses and Frobenius relations pass; an exact small-field point census independently checks the trace recurrence.","The field-valued divisor sieve removes redundant divisor scans and u64 overflow from modular-series setup. Optional stderr stage traces preserve the later construction bottlenecks. This stage correction is not reported as an end-to-end speed measurement.","The independent verifier retains squarefreeness, division-polynomial torsion, subgroup closure, the Velu codomain, normalization/coprimality, exact curve-equation substitution and fresh public subgroup transport. Exact Karatsuba, reciprocal reduction, lazy recurrences and block composition agree with the unchanged reference routines; 12 verifier tests pass.","P-224's certificate completed on the original all-index recurrence. The initial eager P-192 replay was stopped before producing a certificate, and its raw construction and logs remain preserved. The memoized exact verifier then completed all checks for degree 10453.","Existing registry records and IC measurements are preserved. New target models are registered with their mapped generators and exact EC1/UIDs, then native covers, aliases, leaderboard rosters and browser identities are refreshed. Root-library and dependent traits export gates are recorded honestly in EXECUTION_NOTES.md.","Sources: the native algorithm API and prior walker modules, KERNEL_PROTOCOL.md and the committed source snapshots. report.rs is the editable report/vector/PDF source. MANIFEST.json binds source, executables, corpus bytes and canonical outputs. PDF pages are rendered with Poppler and visually checked."]{wrapped(&mut p3,&mut y,t,11.);}
    fs::write(study.join("RESULTS.md"), md).unwrap();
    fs::write(study.join("SEARCH.svg"), svg(&ops)).unwrap();
    emit_pdf(
        &study.join("SEARCH.pdf"),
        vec![(595, 842, p1), (842, 595, vector_pdf(&ops)), (595, 842, p3)],
    );
    fn visit(root: &Path, dir: &Path, out: &mut Vec<Value>) {
        for item in fs::read_dir(dir).unwrap() {
            let p = item.unwrap().path();
            if p.is_dir() {
                if p.file_name().unwrap() != "target" {
                    visit(root, &p, out);
                }
            } else if p.file_name().unwrap() != "MANIFEST.json" {
                let b = fs::read(&p).unwrap();
                out.push(json!({"path":p.strip_prefix(root).unwrap().to_string_lossy(),"bytes":b.len(),"sha256":hash(&b)}));
            }
        }
    }
    let mut files = vec![];
    visit(study, study, &mut files);
    files.sort_by_key(|v| v["path"].as_str().unwrap().to_owned());
    let mut canonical = vec![];
    for p in [
        "docs/curves/registry.json",
        "docs/curves/covers.json",
        "docs/curves/cover-links.yaml",
        "docs/curves/traits.json",
        "src/cryptanalysis/curve_aliases.json",
        "docs/ic/leaderboard.json",
        "docs/ic/LEADERBOARD.md",
        "docs/ic-leaderboard.html",
        "docs/browser/data.json",
    ] {
        let b = fs::read(p).unwrap();
        canonical.push(json!({"path":p,"bytes":b.len(),"sha256":hash(&b)}));
    }
    let mut sources = vec![];
    visit(Path::new(""), Path::new("isogeny_algos/src"), &mut sources);
    for p in [
        "isogeny_algos/Cargo.toml",
        "isogeny_algos/Cargo.lock",
        "src/cryptanalysis/isogeny_walk/field.rs",
        "src/cryptanalysis/isogeny_walk/poly.rs",
        "src/cryptanalysis/isogeny_walk/kernel.rs",
        "src/cryptanalysis/isogeny_walk/curve.rs",
        "src/cryptanalysis/isogeny_walk/modpoly.rs",
        "src/hash/sha256.rs",
    ] {
        let b = fs::read(p).unwrap();
        sources.push(json!({"path":p,"bytes":b.len(),"sha256":hash(&b)}));
    }
    let mut executables = vec![];
    for p in [
        "/private/tmp/p192-p224-continuation.f9sck4/bin/isogeny-algos-sieve-97c80540b",
        "/private/tmp/p192-p224-continuation.f9sck4/kernel-first-v3",
        "/private/tmp/p192-p224-continuation.f9sck4/supervisor",
        "/private/tmp/p192-p224-continuation.f9sck4/supervisor-single",
        "/private/tmp/p192-p224-continuation.f9sck4/replay-one-eager",
        "/private/tmp/p192-p224-continuation.f9sck4/replay-one-memoized",
        "/private/tmp/p192-p224-continuation.f9sck4/replay/release/continuation-report",
    ] {
        let b = fs::read(p).unwrap();
        executables
            .push(json!({"path":p,"bytes":b.len(),"sha256":hash(&b),"included_in_git":false}));
    }
    let rev = std::process::Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .unwrap();
    assert!(rev.status.success());
    let rev = String::from_utf8(rev.stdout).unwrap();
    fs::write(study.join("MANIFEST.json"),format!("{}\n",serde_json::to_string_pretty(&json!({"schema":"high-prime-isogeny-continuation/v1","date":"2026-10-09","source_snapshot_revision":rev.trim(),"algorithm_source_revision":"97c80540b","kernel_helper_source_revision":"2fe16c12b","original_study":"research/p192_p224_large_degree_isogenies_20261008","scope":"individual certified maps; PARTIAL two-eigenline enumeration","files":files,"sources":sources,"executables":executables,"canonical_outputs":canonical})).unwrap())).unwrap();
    println!("Verified continuation report, diagram, PDF, registry and manifest written.");
}

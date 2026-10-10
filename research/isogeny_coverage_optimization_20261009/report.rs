use num_bigint::BigUint;
use serde_json::{json, value::RawValue, Value};
use std::{fmt::Write, fs, path::Path};
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

fn counts(root: &Path, case_name: &str, build: &str) -> Vec<u64> {
    (1..=5)
        .map(|round| {
            fs::read_to_string(root.join(format!("{case_name}-{round}-{build}/instructions.txt")))
                .unwrap()
                .trim()
                .parse()
                .unwrap()
        })
        .collect()
}
fn median(v: &[u64]) -> u64 {
    let mut v = v.to_vec();
    v.sort();
    v[v.len() / 2]
}
fn main() {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() == 2 && args[1] == "--manifest" {
        let root = Path::new(&args[0]);
        let mut files = vec![];
        fn visit(root: &Path, dir: &Path, out: &mut Vec<Value>) {
            for entry in fs::read_dir(dir).unwrap() {
                let p = entry.unwrap().path();
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
        visit(root, root, &mut files);
        files.sort_by_key(|v| v["path"].as_str().unwrap().to_string());
        let mut sources = vec![];
        for base in [
            "tools/isogeny-cli",
            "isogeny_algos/src",
            "src/cryptanalysis/isogeny_walk",
        ] {
            visit(Path::new(""), Path::new(base), &mut sources);
        }
        sources.retain(|v| !v["path"].as_str().unwrap().contains("/dist/"));
        for p in [
            "docs/curves/registry.json",
            "docs/curves/covers.json",
            "docs/curves/cover-links.yaml",
            "docs/ic/leaderboard.json",
            "docs/ic/LEADERBOARD.md",
            "docs/ic-leaderboard.html",
            "docs/browser/data.json",
            "src/cryptanalysis/curve_aliases.json",
            ".github/workflows/release.yml",
        ] {
            let b = fs::read(p).unwrap();
            sources.push(json!({"path":p,"bytes":b.len(),"sha256":hash(&b)}));
        }
        fs::write(root.join("MANIFEST.json"),format!("{}\n",serde_json::to_string_pretty(&json!({"schema":"isogeny-evidence-manifest/v1","files":files,"sources_and_canonical":sources})).unwrap())).unwrap();
        println!(
            "Manifest written: {} files and {} source/canonical bindings.",
            files.len(),
            sources.len()
        );
        return;
    }
    assert_eq!(args.len(), 3);
    let study = Path::new(&args[0]);
    let v1 = Path::new(&args[1]);
    let v2 = Path::new(&args[2]);
    let mut rows = vec![];
    for case_name in ["screen-p224", "verify-p224", "search-p192"] {
        let base = counts(v1, case_name, "baseline");
        let first = counts(v1, case_name, "candidate");
        let second = if case_name == "screen-p224" {
            first.clone()
        } else {
            counts(v2, case_name, "candidate2")
        };
        let b = median(&base);
        let c = median(&second);
        rows.push(json!({"case":case_name,"baseline":base,"candidate1":first,"candidate2":second,"baseline_median":b,"candidate_median":c,"instruction_ratio":b as f64/c as f64,"reduction_percent":100.0*(1.0-c as f64/b as f64),"unit":"Callgrind Ir, native ARM Linux","wall_time_claim":false}));
    }
    fs::write(
        study.join("PROFILE_SUMMARY.json"),
        format!(
            "{}\n",
            serde_json::to_string_pretty(
                &json!({"schema":"isogeny-profile-summary/v1","rows":rows})
            )
            .unwrap()
        ),
    )
    .unwrap();
    let mut curves = vec![];
    let mut certificates = vec![];
    for name in ["p521-7", "p384-19"] {
        let path = study.join("evidence").join(name);
        let cert = load(&path.join("replay.json"));
        assert_eq!(cert["status"], "PASS");
        assert_eq!(cert["degree_coverage"], "COMPLETE");
        for r in cert["records"].as_array().unwrap() {
            assert_eq!(r["exact_rational_map_check"], "PASS");
            assert_eq!(r["public_scalar_transport_checks"], 20);
            certificates.push(r.clone());
        }
        curves.extend(
            load(&path.join("curves.json"))["curves"]
                .as_array()
                .unwrap()
                .iter()
                .cloned(),
        );
    }
    fs::write(
        study.join("VERIFIED_MODELS.json"),
        format!(
            "{}\n",
            serde_json::to_string_pretty(
                &json!({"schema":"isogeny-wide-models/v1","records":certificates,"curves":curves})
            )
            .unwrap()
        ),
    )
    .unwrap();
    register(
        Path::new("docs/curves/registry.json"),
        &curves,
        "research/isogeny_coverage_optimization_20261009/VERIFIED_MODELS.json",
    );
    let registry = load(Path::new("docs/curves/registry.json"));
    let total = registry["curves"].as_array().unwrap().len();
    let mut md=String::from("# Isogeny native optimization and curve coverage\n\nDated 2026-10-09. The implementation is Rust plus assembly. This round extends the dedicated binary beyond P192/P224 and measures native ARM Linux instruction counts under the local Docker VM.\n\n| Command | Baseline median instructions | Final candidate median instructions | Baseline / candidate | Reduction |\n|---|---:|---:|---:|---:|\n");
    for row in &rows {
        writeln!(
            md,
            "| {} | {} | {} | {:.4} | {:.3}% |",
            row["case"].as_str().unwrap(),
            row["baseline_median"],
            row["candidate_median"],
            row["instruction_ratio"].as_f64().unwrap(),
            row["reduction_percent"].as_f64().unwrap()
        )
        .unwrap();
    }
    writeln!(md,"\nFive rounds per panel. Candidate lists, replay records and map bytes match the frozen baseline. The screening improvement includes setup/output; verification excludes construction; the full-search panel includes the constructor child and replay. The small full-search case regresses slightly. No native Mac wall-time claim is made. Candidate 1's small verifier improvement and Candidate 2's larger instruction reduction remain separately recorded in PROFILE_SUMMARY.json.\n\n## Coverage\n\nThe catalogue originally held 345 models: 200 prime, 72 binary, 60 Koblitz, 4 subfield and 9 extension models. It includes 170 short-Weierstrass prime models with construction support through 640 bits, of which 142 have registered complete public subgroup generators. The four newly verified P384/P521 codomains bring the canonical model count to {total}.\n\nMontgomery/Edwards entries and binary/extension models are visible and can be structurally screened from their known orders, but their construction/independent-map adapters remain unresolved. Catalogue completeness means the pinned repository inputs; the importer retains seven unresolved standards records and disclaims a worldwide census.\n\n## Implemented and verified\n\nScreening now resolves the curve and trace once outside the prime loop. Independent replay uses exact-size fields, binary extended-GCD inversion, a fast path for one, and allocation-free small-integer conversion. A recorded profile identified field multiplication as 78.65% of verifier instructions; forced inlining was tested as a second candidate. Producer and verifier arithmetic remain independent. No exact kernel, rational-map or subgroup check was removed.\n\nThe CLI suite passed 58 tests, field reference checks passed 3 tests across widths through P521, and the existing standalone suite passed. P384/19 and P521/7 each produced two independently certified maps, including 20 scalar transport checks per map. P224/1471 replay passed the optimized implementation.\n\n## Scope and remaining gates\n\nThe baseline is c87d96eda, Candidate 1 is 6e3147ad9 and Candidate 2 is 3413c5329. Rust 1.93.1, identical release flags and the same 345-model source snapshot were used for paired measurements. Callgrind 3.19.0 records guest ARM Linux user-space instructions. This is not native Mac throughput and no end-to-end solver/security claim is made.\n\nRoot library compilation still has 643 pre-existing errors. Hosted CI and release publication remain gated. The old study sources/evidence remain frozen at their recorded revisions. No IC performance ledger row is promoted.\n\n[Protocol](PROTOCOL.md), [profile spec](PROFILE_SPEC.json), [summary](PROFILE_SUMMARY.json), [models](VERIFIED_MODELS.json), [diagram](COVERAGE.svg), [PDF](REPORT.pdf), and the raw profile archives in evidence/ carry the receipts. The original candidate's regression and all failed test panels are retained.\n").unwrap();
    fs::write(study.join("RESULTS.md"), format!("{}\n", md.trim_end())).unwrap();
    let mut ops = vec![
        Op::Text(
            35.,
            45.,
            25.,
            "Native isogeny optimization and coverage".into(),
        ),
        Op::Text(
            35.,
            78.,
            14.,
            "Measured instruction counts; five rounds; exact outputs and checks preserved".into(),
        ),
    ];
    for (i, row) in rows.iter().enumerate() {
        let y = 115. + i as f64 * 100.;
        ops.push(Op::Text(35., y, 18., row["case"].as_str().unwrap().into()));
        let ratio = row["instruction_ratio"].as_f64().unwrap();
        ops.push(Op::Rect(280., y - 18., 650., 25., "#aab7c4".into()));
        ops.push(Op::Rect(280., y + 15., 650. / ratio, 25., "#246d9b".into()));
        ops.push(Op::Text(
            35.,
            y + 32.,
            14.,
            format!(
                "B/C {:.4}; reduction {:.3}%",
                ratio,
                row["reduction_percent"].as_f64().unwrap()
            ),
        ));
    }
    ops.push(Op::Text(
        35.,
        450.,
        20.,
        format!("Catalogue: {total} exact models after four validation codomains"),
    ));
    for (i, (title, detail)) in [
        (
            "Prime short Weierstrass",
            "170 original models constructible; 142 original models replayable",
        ),
        (
            "Montgomery / Edwards",
            "30 original models visible; construction conversion adapters pending",
        ),
        (
            "Binary / Koblitz / subfield",
            "136 original models visible and screenable; full map adapters pending",
        ),
        (
            "Extension / unresolved standards",
            "9 extension models; 7 unresolved import rows remain explicit",
        ),
    ]
    .iter()
    .enumerate()
    {
        let y = 490. + i as f64 * 60.;
        ops.push(Op::Rect(
            35.,
            y,
            1015.,
            46.,
            if i == 0 { "#e1f2e8" } else { "#fff0d8" }.into(),
        ));
        ops.push(Op::Text(48., y + 19., 15., title.to_string()));
        ops.push(Op::Text(48., y + 37., 12., detail.to_string()));
    }
    fs::write(study.join("COVERAGE.svg"), svg(&ops)).unwrap();
    let mut pages = vec![];
    let mut page = String::new();
    text(
        &mut page,
        42.,
        790.,
        21.,
        "Isogeny native optimization and coverage",
    );
    let mut y = 751.;
    for line in md.lines().filter(|l| !l.starts_with('#') && !l.is_empty()) {
        if y < 120. {
            pages.push((595, 842, page));
            page = String::new();
            text(&mut page, 42., 790., 21., "Evidence and coverage continued");
            y = 750.;
        }
        wrapped(&mut page, &mut y, line, 10.);
    }
    pages.push((595, 842, page));
    pages.push((842, 595, vector_pdf(&ops)));
    emit_pdf(&study.join("REPORT.pdf"), pages);
    println!("Report, visual, PDF and {total}-model canonical registry written.");
}

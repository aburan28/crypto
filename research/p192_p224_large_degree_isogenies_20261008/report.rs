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
    entries.sort_by_key(|r| {
        let v: Value = serde_json::from_str(r.get()).unwrap();
        let rank = match v["family"].as_str().unwrap() {
            "koblitz" => 0,
            "subfield" => 1,
            "binary" => 2,
            "prime" => 3,
            _ => 4,
        };
        let degree = v["params"]["n"]
            .as_u64()
            .or_else(|| v["params"]["m"].as_u64())
            .unwrap_or(0);
        let p = v["params"]["p"]
            .as_str()
            .and_then(|s| BigUint::parse_bytes(s.as_bytes(), 10))
            .unwrap_or_default();
        let a = if rank == 0 {
            v["params"]["a"].as_u64().unwrap_or(0)
        } else {
            0
        };
        (rank, degree, p, a, v["slug"].as_str().unwrap().to_owned())
    });
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

struct Stats {
    name: String,
    source: Value,
    screen: Value,
    attempts: Vec<Value>,
    split: usize,
    inert: usize,
    repeated: usize,
    pass: usize,
    timeouts: usize,
    fail: usize,
    max: Option<u64>,
    maps: usize,
}
#[derive(Clone)]
enum Op {
    Text(f64, f64, f64, String),
    Rect(f64, f64, f64, f64, String),
    Line(f64, f64, f64, f64, String),
    Dot(f64, f64, f64, String),
}
fn scene(stats: &[Stats], replay: &Value) -> Vec<Op> {
    let mut o = vec![
        Op::Text(
            30.,
            40.,
            25.,
            "P-192 / P-224: large prime-degree search".into(),
        ),
        Op::Text(
            30.,
            72.,
            15.,
            "Explicit maps and structural screening have different verification obligations."
                .into(),
        ),
    ];
    for (i, s) in stats.iter().enumerate() {
        let y = 105. + i as f64 * 190.;
        let records = replay["records"].as_array().unwrap();
        let root = records.iter().find(|r| {
            r["ell"].as_u64() == s.max
                && r["source_icv1"].as_str().unwrap().contains(if i == 0 {
                    "fp192"
                } else {
                    "fp224"
                })
        });
        o.push(Op::Rect(30., y, 465., 145., "#edf4fb".into()));
        o.push(Op::Text(
            45.,
            y + 25.,
            18.,
            format!("{}: source model", s.name),
        ));
        o.push(Op::Text(
            45.,
            y + 49.,
            11.,
            s.source["source_icv1"].as_str().unwrap().into(),
        ));
        o.push(Op::Text(
            45.,
            y + 70.,
            12.,
            s.source["source_ec1"].as_str().unwrap().into(),
        ));
        let uid = s.source["source_curve_uid"].as_str().unwrap();
        o.push(Op::Text(45., y + 93., 10., uid[..53].into()));
        o.push(Op::Text(45., y + 108., 10., uid[53..].into()));
        o.push(Op::Text(
            45.,
            y + 130.,
            12.,
            format!(
                "{} certified maps; largest degree {}",
                s.maps,
                s.max.map_or("none".into(), |l| l.to_string())
            ),
        ));
        if let Some(root) = root {
            o.push(Op::Rect(620., y, 450., 145., "#e7f4ee".into()));
            o.push(Op::Text(
                635.,
                y + 25.,
                18.,
                "One independently certified target".into(),
            ));
            o.push(Op::Text(
                635.,
                y + 49.,
                11.,
                root["target_icv1"].as_str().unwrap().into(),
            ));
            o.push(Op::Text(
                635.,
                y + 70.,
                11.,
                root["target_ec1"].as_str().unwrap().into(),
            ));
            let uid = root["target_curve_uid"].as_str().unwrap();
            o.push(Op::Text(635., y + 93., 10., uid[..53].into()));
            o.push(Op::Text(635., y + 108., 10., uid[53..].into()));
            o.push(Op::Text(
                635.,
                y + 130.,
                12.,
                "Kernel, codomain, and public subgroup transport: PASS".into(),
            ));
            o.push(Op::Line(495., y + 77., 620., y + 77., "#136e4e".into()));
            o.push(Op::Line(620., y + 77., 608., y + 71., "#136e4e".into()));
            o.push(Op::Line(620., y + 77., 608., y + 83., "#136e4e".into()));
            o.push(Op::Text(
                507.,
                y + 64.,
                13.,
                format!("degree {} ->", s.max.unwrap()),
            ));
        } else {
            o.push(Op::Rect(620., y, 450., 145., "#fff1df".into()));
            o.push(Op::Text(
                635.,
                y + 25.,
                18.,
                "Explicit map remains unresolved".into(),
            ));
            o.push(Op::Text(
                635.,
                y + 57.,
                13.,
                format!("{} split candidate degrees in the frozen screen", s.split),
            ));
            o.push(Op::Text(
                635.,
                y + 82.,
                13.,
                format!(
                    "{} construction attempts; {} timeouts",
                    s.attempts.len(),
                    s.timeouts
                ),
            ));
            o.push(Op::Text(
                635.,
                y + 110.,
                12.,
                "No target edge is asserted from a candidate degree.".into(),
            ));
        }
    }
    o.push(Op::Text(
        30.,
        510.,
        20.,
        "Coverage by prime degree (logarithmic degree axis)".into(),
    ));
    o.push(Op::Text(
        30.,
        535.,
        13.,
        "Population: every prime in [67, 4093], separately for each source. Counts are exact."
            .into(),
    ));
    let x = |l: u64| 90. + ((l as f64).ln() - 67f64.ln()) / (4093f64.ln() - 67f64.ln()) * 940.;
    for (i, s) in stats.iter().enumerate() {
        let y = 580. + i as f64 * 68.;
        o.push(Op::Text(30., y + 5., 14., s.name.clone()));
        for r in s.screen["screen"].as_array().unwrap() {
            let ell = r["ell"].as_u64().unwrap();
            let state = s
                .attempts
                .iter()
                .find(|a| a["ell"].as_u64() == Some(ell))
                .map(|a| a["receipt"]["status"].as_str().unwrap());
            let color = match state {
                Some("PASS") => "#136e4e",
                Some("TIMEOUT") => "#cc7c14",
                Some(_) => "#ad3333",
                None => {
                    if r["status"] == "split" {
                        "#819ec5"
                    } else {
                        "#d6dbe3"
                    }
                }
            };
            o.push(Op::Dot(
                x(ell),
                y,
                if state.is_some() { 4. } else { 1.7 },
                color.into(),
            ));
        }
    }
    for l in [67, 127, 257, 509, 1009, 2003, 4093] {
        o.push(Op::Line(x(l), 668., x(l), 675., "#333333".into()));
        o.push(Op::Text(x(l) - 14., 694., 12., l.to_string()));
    }
    o.push(Op::Text(30.,730.,13.,"Green: replayed degree. Orange: timed out. Red: failed. Blue: split candidate only. Gray: other screen status.".into()));
    o
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
fn pdf(path: &Path, stats: &[Stats], ops: &[Op]) {
    let mut p1 = String::new();
    text(&mut p1, 42., 790., 22., "Large prime-degree isogenies");
    text(
        &mut p1,
        42.,
        763.,
        13.,
        "P-192 and P-224 | bounded search | 2026-10-08",
    );
    let mut y = 728.;
    wrapped(&mut p1,&mut y,"Question: construct rational prime-degree maps using the new Hecke/Newton modular-polynomial and BMSS fastElkies' algorithms, including degrees above 1009.",11.);
    for s in stats {
        wrapped(&mut p1,&mut y,&format!("{}: {} prime degrees screened; {} split, {} inert, {} repeated eigenvalue. {} explicit construction attempts; {} passed, {} timed out, {} failed. {} independently certified maps. Largest certified degree: {}.",s.name,s.screen["screen"].as_array().unwrap().len(),s.split,s.inert,s.repeated,s.attempts.len(),s.pass,s.timeouts,s.fail,s.maps,s.max.map_or("none".into(), |l| l.to_string())),11.);
    }
    text(
        &mut p1,
        42.,
        y,
        13.,
        "Construction probes at or beyond degree 1009",
    );
    y -= 24.;
    for s in stats {
        for a in s
            .attempts
            .iter()
            .filter(|a| a["ell"].as_u64().unwrap() >= 1009)
        {
            text(
                &mut p1,
                42.,
                y,
                10.,
                &format!(
                    "{} | degree {} | {} | receipt: evidence-v4/{}/ell-{}.receipt.json",
                    s.name,
                    a["ell"],
                    a["receipt"]["status"].as_str().unwrap(),
                    a["curve"].as_str().unwrap(),
                    a["ell"]
                ),
            );
            y -= 16.;
        }
    }
    y -= 8.;
    wrapped(&mut p1,&mut y,"Bounds: all primes from 67 through 4093 were screened. Every split degree through 257 was attempted, plus the first split prime at or above 509, 1009, 2003 and 4000 for each source. Each construction had a 180-second budget and a sampled 8 GiB RSS cap. Other split primes remain candidates.",11.);
    wrapped(&mut p1,&mut y,"Verification: the existing walker's independent implementation checks squarefreeness, torsion, subgroup closure, and the Velu codomain. Exact rational-map substitution and 20 fresh public scalar-transport checks pass for every map; each mapped nonidentity generator satisfies the published prime order.",11.);
    wrapped(&mut p1,&mut y,"Interpretation boundary: a split Frobenius polynomial supplies candidate eigenlines. Only completed, independently replayed kernels and maps count as certified constructions. A timeout is a resource-bounded outcome, not a nonexistence result. ECDLP-cost changes were not measured.",11.);
    text(&mut p1, 42., y - 8., 13., "Sources and evidence");
    y -= 36.;
    for t in [
        "BMSS: https://arxiv.org/abs/cs/0609020",
        "New native algorithms: isogeny_algos/docs/USAGE.md",
        "Frozen protocol: PROTOCOL.md; raw evidence: evidence-v4/",
        "Search: evidence-v4/search.json; independent certification: evidence-v4/replay.json",
        "Exact target models and generators: evidence-v4/curves.json",
    ] {
        wrapped(&mut p1, &mut y, t, 10.);
    }
    let mut p2 = String::from("q 0.71 0 0 -0.71 26 568 cm\n");
    for op in ops {
        match op {
            Op::Text(x, y, size, t) => writeln!(
                p2,
                "0.08 0.14 0.2 rg BT /F1 {size} Tf 1 0 0 -1 {x} {y} Tm ({}) Tj ET",
                escape(t)
            )
            .unwrap(),
            Op::Rect(x, y, w, h, c) => {
                let (r, g, b) = rgb(c);
                writeln!(p2, "{r} {g} {b} rg {x} {y} {w} {h} re f").unwrap()
            }
            Op::Line(x, y, a, b, c) => {
                let (r, g, bb) = rgb(c);
                writeln!(p2, "{r} {g} {bb} RG 2 w {x} {y} m {a} {b} l S").unwrap()
            }
            Op::Dot(x, y, r, c) => {
                let (rr, g, b) = rgb(c);
                writeln!(
                    p2,
                    "{rr} {g} {b} rg {} {} {} {} re f",
                    x - r,
                    y - r,
                    2. * r,
                    2. * r
                )
                .unwrap()
            }
        }
    }
    p2.push_str("Q\n");
    let mut p3 = String::new();
    text(&mut p3, 42., 790., 21., "Validation and reproduction");
    let mut y = 751.;
    for t in ["Two native fixes were required: P-192 now selects exactly three Montgomery limbs, and P-224's Hecke precondition compares the full characteristic against a small integer rather than its low word (1). The final standalone release suite passes 111 tests, including the quadratic-wrapper consistency check; the supervisor's two tests also pass.",
        "Evidence history: evidence/ retains the initial P-192 preflight failure. evidence-v2/ retains the first 40-degree search and all 16 original P-224 assertion failures. It also preserves an interrupted P-192 degree-149 invocation without a receipt, followed by the sealed retry. Those extra invocations are separate from the final distinct-degree counts.",
        "Final evidence-v4/ reuses the 24 sealed P-192 trials byte for byte and executes fresh P-224 trials. evidence-v3/ retains a preflight RSS-monitor shutdown race. The corrected monitor samples every 25 ms and allows at most 100 ms of unavailable samples while the process exits; persistent unavailability fails closed.",
        "The mandated root cargo test --release --lib failed to compile with 643 pre-existing errors both before and after the rebase. Independent replay and standalone tests pass separately. The repository-wide validation and pre-push publication gate remains unresolved; see validation/lib-test-after-rebase.log.",
        "Environment: Apple M4 Pro, 14 logical cores, 48 GiB RAM, macOS arm64. Elapsed values in process receipts are operational time-budget records on a shared host. No timing or speedup comparison is made.",
        "Reproduce from the feature branch: build isogeny_algos in release mode; run its p192_p224_large_degree_search example with the CLI path and a NEW output directory; build this study's replay binary; replay the directory against docs/curves/registry.json; then run the report binary.",
        "Canonical views: registered target models, native cover certificates, aliases, leaderboard rosters and browser identities are refreshed from the replay. Existing IC measurement bytes, the boundary ledger and progress timeline retain their performance results. See RESULTS.md for validation and delivery status.",
        "Source: report.rs. Final progress: SEARCH-final.log. MANIFEST.json binds protocol, source, executable versions and evidence bytes. SEARCH.svg is editable. Failures and timeouts remain alongside the completed maps."] {wrapped(&mut p3,&mut y,t,11.);}
    let pages = vec![(595, 842, p1), (842, 595, p2), (595, 842, p3)];
    let mut objects = vec![
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
    objects[2] = format!("<< /Type /Pages /Count 3 /Kids [{}] >>", kids.join(" "));
    let mut bytes = String::from("%PDF-1.4\n");
    let mut offsets = vec![0];
    for (i, obj) in objects.iter().enumerate().skip(1) {
        offsets.push(bytes.len());
        writeln!(bytes, "{i} 0 obj\n{obj}\nendobj").unwrap();
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
fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    assert_eq!(args.len(), 2, "usage: report STUDY_DIR REGISTRY_JSON");
    let study = Path::new(&args[0]);
    let evidence = study.join("evidence-v4");
    let search = load(&evidence.join("search.json"));
    assert_eq!(
        search["protocol_sha256"],
        hash(&fs::read(study.join("PROTOCOL.md")).unwrap()),
        "the frozen protocol changed after the search"
    );
    let replay = load(&evidence.join("replay.json"));
    assert_eq!(replay["status"], "PASS");
    for record in replay["records"].as_array().unwrap() {
        assert_eq!(record["status"], "PASS");
        assert_eq!(record["exact_rational_map_check"], "PASS");
    }
    let curves = load(&evidence.join("curves.json"));
    let registry = load(Path::new(&args[1]));
    let mut stats = vec![];
    for (preset, name) in [("p192", "P-192"), ("p224", "P-224")] {
        let screen = load(&evidence.join(format!("{preset}-screen.json")));
        let rows = screen["screen"].as_array().unwrap();
        let split = rows.iter().filter(|r| r["status"] == "split").count();
        let inert = rows.iter().filter(|r| r["status"] == "inert").count();
        let repeated = rows.len() - split - inert;
        let preflight = load(&evidence.join(preset).join("order.json"));
        let slug = &preflight["result"]["icv1"]["slug"];
        let entry = registry["curves"]
            .as_array()
            .unwrap()
            .iter()
            .find(|c| &c["slug"] == slug)
            .unwrap();
        let standard = entry["representations"]
            .as_array()
            .unwrap()
            .iter()
            .find(|r| r["curve"]["cofactor"] == 1)
            .unwrap();
        let source = json!({"source_icv1":slug,"source_ec1":standard["ec1"],"source_curve_uid":standard["curve_uid"]});
        let attempts: Vec<Value> = search["attempts"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|a| a["curve"] == preset)
            .cloned()
            .collect();
        let pass = attempts
            .iter()
            .filter(|a| a["receipt"]["status"] == "PASS")
            .count();
        let timeouts = attempts
            .iter()
            .filter(|a| a["receipt"]["status"] == "TIMEOUT")
            .count();
        let fail = attempts.len() - pass - timeouts;
        let certified: Vec<_> = replay["records"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|r| &r["source_icv1"] == slug)
            .collect();
        assert_eq!(
            certified.len(),
            pass * 2,
            "every completed construction requires two independently certified maps"
        );
        for attempt in attempts.iter().filter(|a| a["receipt"]["status"] == "PASS") {
            let degree_records: Vec<_> = certified
                .iter()
                .filter(|r| r["ell"] == attempt["ell"])
                .collect();
            assert_eq!(degree_records.len(), 2);
            assert!(degree_records.iter().any(|r| r["index"] == 0));
            assert!(degree_records.iter().any(|r| r["index"] == 1));
        }
        let max = certified.iter().map(|r| r["ell"].as_u64().unwrap()).max();
        stats.push(Stats {
            name: name.into(),
            source,
            screen,
            attempts,
            split,
            inert,
            repeated,
            pass,
            timeouts,
            fail,
            max,
            maps: certified.len(),
        });
    }
    register(
        Path::new(&args[1]),
        curves["curves"].as_array().unwrap(),
        "research/p192_p224_large_degree_isogenies_20261008/evidence-v4/curves.json",
    );
    let mut md=String::from("# P-192 / P-224 large prime-degree isogeny search\n\nDated 2026-10-08. The requested search used the new native algorithms and extended prime-degree screening beyond 1009, through 4093. Construction and structural screening are reported separately.\n\n| Source | Prime degrees screened | Split | Inert | Repeated eigenvalue | Attempts | Completed degrees | Timeouts | Failed | Certified maps | Largest certified degree |\n|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n");
    for s in &stats {
        writeln!(
            md,
            "| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |",
            s.name,
            s.screen["screen"].as_array().unwrap().len(),
            s.split,
            s.inert,
            s.repeated,
            s.attempts.len(),
            s.pass,
            s.timeouts,
            s.fail,
            s.maps,
            s.max.map_or("none".into(), |l| l.to_string())
        )
        .unwrap();
    }
    md.push_str("\nThe frozen [protocol](PROTOCOL.md) attempted every split prime from 67 through 257 plus the first split prime at or above 509, 1009, 2003 and 4000 for each source. Each child process had a 180-second construction budget. All remaining split primes are candidate degrees without constructed maps. Each accepted degree supplied two maps, with explicit monic kernel polynomials, codomains and rational x-maps.\n\nThe [independent replay](evidence-v4/replay.json) uses the existing walker's separate field, polynomial, division-polynomial torsion, subgroup-closure and Velu codomain implementations. Each map also passed 20 public scalar-transport checks on fresh deterministic points. The standard generator's image was nonidentity and satisfied the published prime group order. Exact registered target models, generators and EC1 representations are in [curves.json](evidence-v4/curves.json). This replay is independent by implementation on the same host.\n\n![Certified examples and coverage](SEARCH.svg)\n\n## Construction beyond 1009\n\n| Source | Degree | Outcome | Evidence |\n|---|---:|---|---|\n");
    for s in &stats {
        for a in &s.attempts {
            let ell = a["ell"].as_u64().unwrap();
            if ell > 1009 {
                let preset = a["curve"].as_str().unwrap();
                writeln!(
                    md,
                    "| {} | {} | {} | [receipt](evidence-v4/{}/ell-{}.receipt.json) |",
                    s.name,
                    ell,
                    a["receipt"]["status"].as_str().unwrap(),
                    preset,
                    ell
                )
                .unwrap();
            }
        }
    }
    md.push_str(r#"
A split degree is a structural candidate. A failed or timed-out construction leaves its explicit map unresolved and does not prove nonexistence. Construction beyond 1009 remains unresolved where the individual receipts time out. The search measured no ECDLP-cost change. Process elapsed values are operational time-budget records on a shared host, without isolated benchmark measurements.

## Validation, source and delivery status

Two native implementation fixes were required. The initial P-192 preflight exposed a limb-dispatch bug; the failed launch remains in [evidence/](evidence/) and [SEARCH.log](SEARCH.log). The dispatcher now chooses the exact Montgomery limb count. P-224's Hecke modular-polynomial precondition used the low-word characteristic accessor, whose value is 1 for this 224-bit prime. The precondition now uses an exact characteristic comparison, which also forwards through the existing quadratic-field wrapper. Its regression compares the Hecke polynomial with the separate linear-algebra route and requires two verified maps at a split degree. The final [standalone release suite](validation/algorithms-tests-final.log) passed all 111 tests. The [search example](validation/example-tests-final.log) passed both tests. The [exact map-identity tests](validation/exact-map-tests.log) passed, including rejection of a changed numerator and changed target coefficient.

The final table counts distinct planned degrees. [evidence-v2/search.json](evidence-v2/search.json) preserves the first complete 40-degree search, including the 16 P-224 assertion failures before its fix. It also retains one interrupted P-192 degree-149 invocation without a completed receipt, separately from the sealed retry. Final [evidence-v4/](evidence-v4/) reuses all 24 sealed P-192 trials and their screen byte for byte after command/digest checks; the P-224 trials were executed fresh with the corrected CLI. The original CLI is preserved locally as `isogeny-algos-pre-p224-fix`, with its digest in the manifest.

[evidence-v3/](evidence-v3/) preserves a P-224 preflight RSS-monitor shutdown race without reclassifying its receipt. The corrected supervisor retains the 180-second deadline and sampled 8 GiB RSS cap, samples every 25 ms, and allows at most 100 ms of transient unavailable samples while the child exits. Persistent unavailability fails closed. The final receipts include missing-sample counts and the grace interval. Brief memory growth between samples remains possible.

The mandatory root release-library check failed to compile with 643 existing errors before and after rebasing onto upstream `2fe5cec8a4a9d8d55c8e0abec9fe22852f3a3726`. Its [original log](validation/lib-test.log) and [post-rebase log](validation/lib-test-after-rebase.log) are retained. The independent replay executable imports the original five walker modules directly; every completed map passes kernel checks, exact polynomial substitution, and fresh public subgroup transport. The repository-wide validation and pre-push publication gate remains unresolved until the baseline compiles. See [execution notes](EXECUTION_NOTES.md).

The target registry is updated natively from frozen replayed models. Native cover certificates and the linked cover graph, aliases, leaderboard curve rosters, and browser curve identities are refreshed together. The catalogue tools reproduce the existing outputs exactly before extension. Existing IC measurement bytes are preserved. The boundary ledger, performance scoreboard, progress timeline, and older walker guide receive no new performance row: the experiment supplies maps and exact models without an IC/rho ratio or boundary promotion.

Sources: [BMSS](https://arxiv.org/abs/cs/0609020), [native CLI/API](../../isogeny_algos/docs/USAGE.md), [existing walker](../../docs/isogeny-walk/README.md). The vector and PDF share editable [report.rs](report.rs) source; [SEARCH.pdf](SEARCH.pdf) includes the report and visual.

## Reproduction

```sh
cargo build --manifest-path isogeny_algos/Cargo.toml --locked --release --bin isogeny-algos --example p192_p224_large_degree_search
isogeny_algos/target/release/examples/p192_p224_large_degree_search isogeny_algos/target/release/isogeny-algos NEW_EVIDENCE_DIR c70c32d486a3ac7531fe27f7193d9f09caa58344
cargo build --manifest-path research/p192_p224_large_degree_isogenies_20261008/Cargo.toml --locked --release --offline --bins
research/p192_p224_large_degree_isogenies_20261008/target/release/replay NEW_EVIDENCE_DIR docs/curves/registry.json
```

[COMMANDS.md](COMMANDS.md) records the actual multi-phase execution and catalogue closeout. [MANIFEST.json](MANIFEST.json) binds the frozen protocol, source files, executable versions, canonical outputs, and evidence bytes. Failed builds and preflight launches remain preserved in their original directories.
"#);
    fs::write(study.join("RESULTS.md"), md).unwrap();
    let ops = scene(&stats, &replay);
    fs::write(study.join("SEARCH.svg"), svg(&ops)).unwrap();
    pdf(&study.join("SEARCH.pdf"), &stats, &ops);
    let mut files = vec![];
    fn visit(root: &Path, dir: &Path, out: &mut Vec<Value>) {
        for e in fs::read_dir(dir).unwrap() {
            let p = e.unwrap().path();
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
    visit(study, study, &mut files);
    files.sort_by_key(|f| f["path"].as_str().unwrap().to_owned());
    let mut sources = vec![];
    for p in [
        "isogeny_algos/Cargo.toml",
        "isogeny_algos/Cargo.lock",
        "isogeny_algos/tests/api.rs",
        "isogeny_algos/src/api.rs",
        "isogeny_algos/examples/p192_p224_large_degree_search.rs",
        "src/cryptanalysis/isogeny_walk/field.rs",
        "src/cryptanalysis/isogeny_walk/poly.rs",
        "src/cryptanalysis/isogeny_walk/curve.rs",
        "src/cryptanalysis/isogeny_walk/kernel.rs",
        "src/cryptanalysis/isogeny_walk/modpoly.rs",
        "src/hash/sha256.rs",
        "src/bin/curve_cover_check/checker.rs",
        "src/bin/curve_cover_check/models.rs",
        "src/bin/curve_cover_check/links.rs",
        "src/binary_ecc/f2m.rs",
        "src/cryptanalysis/ecc2k130_guard.rs",
        "src/ct_bignum.rs",
        "src/utils/mod.rs",
        "src/utils/encoding.rs",
        "src/utils/random.rs",
    ] {
        let b = fs::read(p).unwrap();
        sources.push(json!({"path":p,"bytes":b.len(),"sha256":hash(&b)}));
    }
    visit(Path::new(""), Path::new("isogeny_algos/src"), &mut sources);
    sources.sort_by_key(|s| s["path"].as_str().unwrap().to_owned());
    sources.dedup_by(|a, b| a["path"] == b["path"]);
    let mut canonical_outputs = vec![];
    for path in [
        "docs/curves/registry.json",
        "docs/curves/covers.json",
        "docs/curves/cover-links.yaml",
        "src/cryptanalysis/curve_aliases.json",
        "docs/ic/leaderboard.json",
        "docs/ic/LEADERBOARD.md",
        "docs/ic-leaderboard.html",
        "docs/browser/data.json",
    ] {
        let b = fs::read(path).unwrap();
        canonical_outputs.push(json!({"path":path,"bytes":b.len(),"sha256":hash(&b)}));
    }
    let mut executables = vec![];
    for p in [
        "isogeny_algos/target/release/isogeny-algos",
        "isogeny_algos/target/release/isogeny-algos-pre-p224-fix",
        "isogeny_algos/target/release/large-degree-search-supervised",
        "isogeny_algos/target/release/large-degree-search-resumable",
        "isogeny_algos/target/release/large-degree-search-monitor-fixed",
        "research/p192_p224_large_degree_isogenies_20261008/target/release/replay",
        "research/p192_p224_large_degree_isogenies_20261008/target/release/report",
        "research/p192_p224_large_degree_isogenies_20261008/target/release/catalogue-covers",
        "research/p192_p224_large_degree_isogenies_20261008/target/release/catalogue-views",
    ] {
        let b = fs::read(p).unwrap();
        executables.push(json!({"path":p,"bytes":b.len(),"sha256":hash(&b),"optimization":"release, level 3","included_in_git":false}));
    }
    let source_revision = std::process::Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .unwrap();
    assert!(source_revision.status.success());
    let source_revision = String::from_utf8(source_revision.stdout).unwrap();
    let execution_phases = json!([
        {"evidence":"evidence/","purpose":"initial P-192 preflight failure","construction_source_revision":"7819ee8c096b9c99617a8ff14ff001d9f4efd4e9","executable_digest":"not retained"},
        {"evidence":"evidence-v2/","purpose":"first complete 40-degree search; additional interrupted degree-149 invocation retained separately","construction_executable":"isogeny_algos/target/release/isogeny-algos-pre-p224-fix","construction_source_revision":"ddf0799aa95e1edab8ab626ce9d31b75f12f3aa4","supervisor_executables":["isogeny_algos/target/release/large-degree-search-supervised","isogeny_algos/target/release/large-degree-search-resumable"]},
        {"evidence":"evidence-v3/","purpose":"corrected P-224 preflight stopped by the preserved monitor race","construction_executable":"isogeny_algos/target/release/isogeny-algos","construction_source_revision":"631d857e8414aed95c1b0254f0365bc0ff56f3c5","p192":"sealed evidence-v2 files reused"},
        {"evidence":"evidence-v4/","purpose":"final degree outcomes and strongest whole-search independent replay","construction_executable":"isogeny_algos/target/release/isogeny-algos","construction_source_revision":"631d857e8414aed95c1b0254f0365bc0ff56f3c5","supervisor_executable":"isogeny_algos/target/release/large-degree-search-monitor-fixed","p192":"24 sealed evidence-v2 trials copied byte for byte; not new construction invocations","p224":"16 fresh construction trials"}
    ]);
    fs::write(study.join("MANIFEST.json"),format!("{}\n",serde_json::to_string_pretty(&json!({"schema":"large-degree-isogeny-manifest/v1","base_revision":search["source_commit"],"integrated_upstream_revision":"2fe5cec8a4a9d8d55c8e0abec9fe22852f3a3726","source_snapshot_revision":source_revision.trim(),"protocol_sha256":search["protocol_sha256"],"execution_phases":execution_phases,"host":{"cpu":"Apple M4 Pro","logical_cores":14,"ram_gib":48,"os":"macOS arm64","benchmark_isolation":"not a timing benchmark"},"sources":sources,"executables":executables,"canonical_outputs":canonical_outputs,"files":files})).unwrap())).unwrap();
    println!("Report, vector visual, PDF, manifest and model registration written.");
}

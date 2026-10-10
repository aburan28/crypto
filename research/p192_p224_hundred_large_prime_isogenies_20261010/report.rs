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

include!("../isogeny_catalogue_setup_20261009/visuals.rs");
fn files(root: &Path, dir: &Path, rows: &mut Vec<Value>) {
    let mut paths: Vec<_> = fs::read_dir(dir)
        .unwrap()
        .map(|r| r.unwrap().path())
        .collect();
    paths.sort();
    for p in paths {
        if p.is_dir() {
            files(root, &p, rows);
        } else if !matches!(
            p.extension().and_then(|s| s.to_str()),
            Some("pid" | "pending")
        ) {
            let b = fs::read(&p).unwrap();
            rows.push(json!({"path":p.strip_prefix(root).unwrap().to_string_lossy(),"bytes":b.len(),"sha256":hash(&b)}));
        }
    }
}
fn main() {
    let repo = std::env::current_dir().unwrap();
    let study = repo.join("research/p192_p224_hundred_large_prime_isogenies_20261010");
    let mut curves = vec![];
    let mut rows = vec![];
    let mut totals = vec![];
    let mut keys = std::collections::HashSet::new();
    for name in ["p192", "p224"] {
        let path = study.join("evidence").join(name);
        let m = load(&path.join("walk.json"));
        assert_eq!(m["status"], "PASS");
        assert_eq!(m["verified_maps"], 50);
        let (mut construct, mut verify, mut invocation) = (0u64, 0u64, 0u64);
        for step in m["steps"].as_array().unwrap() {
            let dir = path.join(step["directory"].as_str().unwrap());
            let replay = load(&dir.join("replay.json"));
            assert_eq!(replay["status"], "PASS");
            let r = &replay["records"][0];
            assert!(keys.insert(format!("{}:{}", r["source_icv1"], r["certificate_sha256"])));
            assert_eq!(
                step["replay_sha256"],
                hash(&fs::read(dir.join("replay.json")).unwrap())
            );
            assert_eq!(
                step["curves_sha256"],
                hash(&fs::read(dir.join("curves.json")).unwrap())
            );
            construct += step["construction_elapsed_ms_operational_only"]
                .as_u64()
                .unwrap();
            verify += step["verification_elapsed_ms"].as_u64().unwrap();
            invocation += step["invocation_elapsed_ms"].as_u64().unwrap();
            rows.push(r.clone());
            curves.extend(
                load(&dir.join("curves.json"))["curves"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .cloned(),
            );
        }
        totals.push(json!({"curve":name,"ell":m["ell"],"maps":50,"construction_elapsed_ms_operational_only":construct,
            "verification_elapsed_ms":verify,"summed_invocation_elapsed_ms":invocation}));
    }
    assert_eq!(keys.len(), 100);
    register(
        &repo.join("docs/curves/registry.json"),
        &curves,
        "research/p192_p224_hundred_large_prime_isogenies_20261010/RESULTS.md",
    );
    let summary = json!({"schema":"hundred-large-prime-maps/v1","status":"PASS","new_distinct_maps":100,
        "distinct_prime_degrees":2,"degree_coverage":"PARTIAL","records":rows,"totals":totals});
    fs::write(
        study.join("SUMMARY.json"),
        format!("{}\n", serde_json::to_string_pretty(&summary).unwrap()),
    )
    .unwrap();
    let mut md=String::from("# One hundred large prime-degree isogeny maps\n\nThe campaign independently verified 100 distinct maps: 50 consecutive edges in each public P192/P224 isogeny class. There are two distinct prime degrees, 2341 and 1117. The first maps start on the published sources; later maps start on their independently verified codomains.\n\n");
    md.push_str("The [frozen protocol](PROTOCOL.md) defines map-level uniqueness and retains one-eigenline coverage PARTIAL at every step. Existing P192/10453 and P224/1471 maps and repeated replays are excluded. This meets a hundred-map target with repeated prime degrees allowed; a hundred different prime degrees remains a separate task.\n\n| Public class | Prime degree | New maps | Constructor seconds | Independent replay seconds | Summed invocation seconds |\n|---|---:|---:|---:|---:|---:|\n");
    for t in &totals {
        writeln!(
            md,
            "| {} | {} | 50 | {:.3} | {:.3} | {:.3} |",
            t["curve"].as_str().unwrap(),
            t["ell"],
            t["construction_elapsed_ms_operational_only"]
                .as_u64()
                .unwrap() as f64
                / 1000.,
            t["verification_elapsed_ms"].as_u64().unwrap() as f64 / 1000.,
            t["summed_invocation_elapsed_ms"].as_u64().unwrap() as f64 / 1000.
        )
        .unwrap();
    }
    md.push_str("\nThese are complete operational intervals on a shared Apple M4 Pro host, with two concurrent walks. Constructors have an 1800-second limit and 8 GiB RSS cap. Independent replay is outside that constructor cap. These measurements do not establish a matched speedup. All 100 exact kernel, codomain and rational curve-equation checks passed, alongside subgroup-order checks and 20 public scalar-transport checks per map (2000 total). Replay uses separate field and polynomial implementations on the same host.\n\n[SUMMARY.json](SUMMARY.json) lists all exact source/destination ICV1 and EC1 records, degrees and map hashes. [P192 walk](evidence/p192/walk.json) and [P224 walk](evidence/p224/walk.json) link every raw kernel and rational map, resource and stderr receipt, mapped public generator and independent replay. The embedded source revision was clean; both walks retain the binary digest.\n\n![Verified routes and coverage](SEARCH.svg)\n\nThe screen covered every prime from 1010 through 1000000 on both sources. Orders at most 8 selected 2 P192 and 3 P224 candidates; orders at most 64 selected 40 and 38. These are structural candidates. The generalized constructor supports checked extension degrees 2 through 16. All enabled degrees passed field-law tests, and the original degree-3/5 selections remain unchanged.\n\n## Validation and delivery\n\nThe CLI passed 62 tests plus 3 limb-field reference tests. Standalone algorithm suites passed 114 tests. The required root release library check retains 643 pre-existing compilation errors. An additional strict CLI Clippy diagnostic retains 35 errors from inherited/generated modules; that check is not presented as passing. The required algorithm Clippy job is checked separately. PR #1606 remains open and the release remains unpublished.\n\nCanonical updates append the 100 verified models and mapped generators, rebuild covers and links, and refresh aliases, leaderboard rosters and browser joins. The complete route graph is in [ROUTES.mmd](ROUTES.mmd). The IC boundary ledger, performance scoreboard and timeline retain their existing measurements because this round does not measure an IC improvement. Earlier frozen studies retain their historical findings.\n\nRun the actual executable with: isogeny walk --curve p224 --ell 1117 --count 50 --out ./p224-walk. Replay it with: isogeny verify-walk --input ./p224-walk. Both arithmetic implementations and pinned public anchors are embedded. Runtime Cargo and a checkout are unnecessary. See [tool instructions](../../tools/isogeny-cli/README.md), [commands](COMMANDS.md) and [manifest](MANIFEST.json).\n");
    md = md.replace("The required root release library check retains 643 pre-existing compilation errors.",
        "The initial required root compilation failed with 643 errors. Three existing upstream type-repair commits restored compilation. The first complete runtime attempt exposed catalogue parser gaps, a serde-feature-dependent fixture, one sparse fixture omission and a one MiB digest buffer overflowing the default test stack. Their original diagnostics are retained, and the parser, fixture reference and buffer allocation were repaired without changing map acceptance checks.");
    let log_path = study.join("evidence/validation/root-final-fixed.log");
    if log_path.exists() {
        let log = fs::read_to_string(&log_path).unwrap();
        let outcome = log
            .lines()
            .rev()
            .find(|l| {
                l.contains("test result:")
                    || l.contains("could not compile")
                    || l.contains("fatal runtime error")
            })
            .unwrap_or("IN_PROGRESS: final required library test is still running");
        writeln!(md,"\nFinal required library result: {}. See [raw validation](evidence/validation/root-final-fixed.log).",outcome).unwrap();
    } else {
        md.push_str("\nThe final required library suite is running. Push/publication remains gated on its result.\n");
    }
    let gate_prefix = study.join("evidence/validation/root-complete-prefix.log");
    if gate_prefix.exists() {
        let log = fs::read_to_string(&gate_prefix).unwrap();
        let failed = log
            .lines()
            .filter(|l| l.starts_with("test ") && l.ends_with("FAILED"))
            .count();
        writeln!(md,"\nThe unsandboxed full library run has at least {failed} observed failures outside the dedicated map verification path. Its retained [log prefix](evidence/validation/root-complete-prefix.log) is a checkpoint, not a completed passing suite. PR push is awaiting the owner's explicit choice to continue broader repairs or publish with merge/release still blocked.").unwrap();
    }
    md.push_str("\nThe rebuilt 449-model executable passed 62 CLI tests and 3 independent limb-field reference tests. It also passed a two-step saved-chain replay and complete two-map P384/19 construction/replay from a directory without a checkout. The [Apple silicon executable archive](evidence/isogeny-aarch64-apple-darwin-3c39c5f23.tar.gz) includes BUILD.json and source commit 3c39c5f23e0bd1590d220fb32386ad9ed79bf399. Its source state was clean.\n");
    fs::write(study.join("RESULTS.md"), &md).unwrap();
    let mut ops = vec![
        Op::Text(35., 45., 26., "100 verified large prime-degree maps".into()),
        Op::Text(
            35.,
            77.,
            15.,
            "50 per class; 2 distinct degrees; one of two eigenlines per edge (PARTIAL)".into(),
        ),
    ];
    for (i, class) in ["p192", "p224"].iter().enumerate() {
        let yy = 112. + 178. * i as f64;
        let col = if i == 0 { "#185d91" } else { "#228357" };
        let first = &rows[i * 50];
        let last = &rows[i * 50 + 49];
        ops.push(Op::Rect(35., yy, 1030., 150., "#edf3f8".into()));
        ops.push(Op::Text(
            52.,
            yy + 26.,
            17.,
            format!(
                "{} class: prime degree {} on every edge",
                class.to_uppercase(),
                first["ell"]
            ),
        ));
        ops.push(Op::Text(
            52.,
            yy + 51.,
            12.,
            first["source_icv1"].as_str().unwrap().into(),
        ));
        ops.push(Op::Text(
            52.,
            yy + 73.,
            12.,
            first["target_icv1"].as_str().unwrap().into(),
        ));
        ops.push(Op::Line(66., yy + 89., 1000., yy + 89., col.into()));
        ops.push(Op::Text(
            80.,
            yy + 109.,
            13.,
            "Verified edge 1 -> ... -> verified edge 50; full identities in SUMMARY.json".into(),
        ));
        ops.push(Op::Text(
            52.,
            yy + 133.,
            12.,
            last["target_icv1"].as_str().unwrap().into(),
        ));
        let y = 525. + i as f64 * 55.;
        ops.push(Op::Text(40., y + 20., 16., class.to_uppercase()));
        ops.push(Op::Rect(135., y, 850., 31., col.into()));
        ops.push(Op::Text(1000., y + 22., 17., "50".into()));
    }
    ops.push(Op::Text(
        35.,
        495.,
        19.,
        "New distinct maps independently replayed (count)".into(),
    ));
    for tick in 0..=5 {
        let x = 135. + 170. * tick as f64;
        ops.push(Op::Line(x, 625., x, 632., "#64748b".into()));
        ops.push(Op::Text(x - 5., 650., 12., (tick * 10).to_string()));
    }
    ops.push(Op::Text(990., 671., 12., "maps".into()));
    ops.push(Op::Text(
        35.,
        701.,
        13.,
        "Population: 100 maps, two public curve classes; exact checks, no sampling uncertainty."
            .into(),
    ));
    ops.push(Op::Text(35.,727.,13.,"Kernel + exact rational identity + subgroup preservation + 2000 public scalar checks: PASS".into()));
    fs::write(study.join("SEARCH.svg"), svg(&ops)).unwrap();
    let mut pages = vec![];
    let mut s = String::new();
    let mut y = 750.;
    for line in md.lines() {
        if line.starts_with("![") || line.starts_with("|---") {
            continue;
        }
        if line.starts_with("| ") {
            if y < 95. {
                pages.push((612, 792, s));
                s = String::new();
                y = 750.;
            }
            let cells: Vec<_> = line.split('|').skip(1).take(6).map(str::trim).collect();
            let headings = [
                "Class",
                "Degree",
                "Maps",
                "Construct (s)",
                "Replay (s)",
                "Total (s)",
            ];
            for (i, x) in [42., 116., 190., 245., 360., 462.].iter().enumerate() {
                let value = if cells[0] == "Public class" {
                    headings[i]
                } else {
                    cells[i]
                };
                text(&mut s, *x, y, 9., value);
            }
            writeln!(s, "0.8 0.85 0.9 RG 42 {} m 570 {} l S", y - 5., y - 5.).unwrap();
            y -= 24.;
            continue;
        }
        let heading = line.starts_with('#');
        let mut clean = line.trim_start_matches('#').trim().to_string();
        let mut cursor = 0;
        while let Some(start) = clean[cursor..].find('[').map(|n| n + cursor) {
            let Some(mid) = clean[start..].find("](").map(|n| n + start) else {
                break;
            };
            let Some(end) = clean[mid + 2..].find(')').map(|n| n + mid + 2) else {
                break;
            };
            let replacement = format!("{} ({})", &clean[start + 1..mid], &clean[mid + 2..end]);
            clean.replace_range(start..=end, &replacement);
            cursor = start + replacement.len();
        }
        if heading {
            if y < 120. {
                pages.push((612, 792, s));
                s = String::new();
                y = 750.;
            }
            text(
                &mut s,
                42.,
                y,
                if line.starts_with("# ") { 20. } else { 15. },
                &clean,
            );
            y -= 32.;
            continue;
        }
        let mut segment = String::new();
        for word in clean.split_whitespace() {
            if segment.len() + word.len() > 90 {
                if y < 70. {
                    pages.push((612, 792, s));
                    s = String::new();
                    y = 750.;
                }
                text(&mut s, 42., y, 10., &segment);
                y -= 16.;
                segment.clear();
            }
            if !segment.is_empty() {
                segment.push(' ');
            }
            segment.push_str(word);
        }
        if y < 70. {
            pages.push((612, 792, s));
            s = String::new();
            y = 750.;
        }
        if !segment.is_empty() {
            text(&mut s, 42., y, 10., &segment);
            y -= 16.;
        }
        y -= 8.;
    }
    pages.push((612, 792, s));
    pages.push((1100, 760, vector_pdf(&ops)));
    emit_pdf(&study.join("REPORT.pdf"), pages);
    let mut graph = String::from("flowchart LR\n");
    for (i, r) in rows.iter().enumerate() {
        let class = i / 50;
        let index = i % 50;
        if index == 0 {
            writeln!(
                graph,
                "  c{class}_0[\"{}\"]",
                r["source_icv1"].as_str().unwrap()
            )
            .unwrap();
        }
        writeln!(
            graph,
            "  c{class}_{index} -->|\"prime degree {}; verified PARTIAL\"| c{class}_{}[\"{}\"]",
            r["ell"],
            index + 1,
            r["target_icv1"].as_str().unwrap()
        )
        .unwrap();
    }
    fs::write(study.join("ROUTES.mmd"), &graph).unwrap();
    fs::write(
        repo.join("docs/curves/large-prime-walks-20261010.mmd"),
        &graph,
    )
    .unwrap();
    fs::write(
        repo.join("docs/curves/large-prime-walks-20261010.json"),
        serde_json::to_vec_pretty(&summary).unwrap(),
    )
    .unwrap();
    let mut evidence = vec![];
    files(&study, &study.join("evidence"), &mut evidence);
    for p in [
        "PROTOCOL.md",
        "SUMMARY.json",
        "RESULTS.md",
        "SEARCH.svg",
        "ROUTES.mmd",
        "REPORT.pdf",
        "report.rs",
        "Cargo.toml",
        "COMMANDS.md",
        "ROOT_VALIDATION.md",
        "walk_controls.rs",
        "audit.rs",
        "Cargo.lock",
        "p192-screen-order8.json",
        "p192-screen-order64.json",
        "p224-screen-order8.json",
        "p224-screen-order64.json",
    ] {
        let b = fs::read(study.join(p)).unwrap();
        evidence.push(json!({"path":p,"bytes":b.len(),"sha256":hash(&b)}));
    }
    let mut canonical = vec![];
    for p in [
        "docs/curves/registry.json",
        "docs/curves/covers.json",
        "docs/curves/cover-links.yaml",
        "docs/curves/large-prime-walks-20261010.json",
        "docs/curves/large-prime-walks-20261010.mmd",
        "docs/ic/leaderboard.json",
        "docs/ic/LEADERBOARD.md",
        "docs/ic-leaderboard.html",
        "docs/browser/data.json",
        "src/cryptanalysis/curve_aliases.json",
    ] {
        let b = fs::read(repo.join(p)).unwrap();
        canonical.push(json!({"path":p,"bytes":b.len(),"sha256":hash(&b)}));
    }
    fs::write(study.join("MANIFEST.json"),format!("{}\n",serde_json::to_string_pretty(&json!({
        "schema":"hundred-large-prime-manifest/v1","status":"PASS","files":evidence,
        "research_status":"PASS","publication_status":"PENDING_REQUIRED_CHECKS",
        "canonical_outputs":canonical,
        "source_revision":"49a280dc60972c9d16b2205b4652cc6c5d87284e",
        "host":{"cpu":"Apple M4 Pro","logical_cores":14,"ram_gib":48,"os":"macOS arm64","isolation":"shared host; operational timings"}})).unwrap())).unwrap();
}

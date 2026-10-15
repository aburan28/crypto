//! Render the frozen fourth-rung control; no timing result is inferred.
use isogeny_algos::{json::Json, sha256::sha256_hex};
use std::{fmt::Write as _, path::Path};
fn s<'a>(j: &'a Json, k: &str) -> &'a str {
    j.get(k).unwrap().as_str().unwrap()
}
fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert_eq!(args.len(), 3, "make-rung4-report FROZEN_RUN FRESH_DIR");
    let root = Path::new(&args[1]);
    let raw = std::fs::read(root.join("INCREMENTAL_COMPARISON.json")).unwrap();
    let comparison = Json::parse(std::str::from_utf8(&raw).unwrap()).unwrap();
    assert_eq!(s(&comparison, "status"), "PASS");
    let Json::Arr(cases) = comparison.get("cases").unwrap() else {
        panic!("missing cases")
    };
    assert_eq!(cases.len(), 96);
    let mut csv = String::from("bits,ell,algorithm,operation,previous,candidate\n");
    for row in cases {
        for op in ["add_sub", "mul", "inv", "neg"] {
            let counts = row.get("operations").unwrap().get(op).unwrap();
            assert_eq!(s(counts, "baseline"), s(counts, "candidate"));
            writeln!(
                csv,
                "{},{},{},{op},{},{}",
                s(row, "bits"),
                s(row, "ell"),
                s(row, "algorithm"),
                s(counts, "baseline"),
                s(counts, "candidate")
            )
            .unwrap();
        }
    }
    let cumulative_raw = std::fs::read(root.join("CUMULATIVE_COMPARISON.json")).unwrap();
    let cumulative = Json::parse(std::str::from_utf8(&cumulative_raw).unwrap()).unwrap();
    assert_eq!(s(&cumulative, "status"), "PASS");
    assert_eq!(
        cumulative.get("candidate_run_sha256"),
        comparison.get("candidate_run_sha256")
    );
    assert_eq!(
        s(&comparison, "candidate_run_sha256"),
        sha256_hex(&std::fs::read(root.join("counts.jsonl")).unwrap())
    );
    let revision = std::fs::read_to_string(root.join("SOURCE_REVISION.txt")).unwrap();
    let binary_hash = sha256_hex(&std::fs::read(root.join("bin/isogeny-perf")).unwrap());
    let svg = r##"<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="180"><rect width="1000" height="180" fill="white"/><g font-family="Arial,sans-serif" fill="#17324d"><text x="20" y="30" font-size="22">Fourth BMSS candidate: exact controls, timing pending</text><g fill="#edf5fa" stroke="#7293a7"><rect x="20" y="55" width="280" height="75"/><rect x="360" y="55" width="280" height="75"/><rect x="700" y="55" width="280" height="75"/></g><text x="35" y="85" font-size="17">Original / previous / candidate</text><text x="35" y="112" font-size="15">8 methods · 4 widths · 3 degrees</text><text x="320" y="100" font-size="26">→</text><text x="375" y="85" font-size="17">96 exact output matches</text><text x="375" y="112" font-size="15">All four counters unchanged</text><text x="660" y="100" font-size="26">→</text><text x="715" y="85" font-size="17">Separate CPU campaigns</text><text x="715" y="112" font-size="15">Cumulative and incremental</text><text x="20" y="165" font-size="16">Shared Mac: zero qualified wall samples. An unchanged call count does not establish a timing gain.</text></g></svg>"##;
    let body = "The small-precision Weierstrass recurrence now takes a dedicated direct-division path. The batch-inverted path removes per-iteration Option dispatch and unused denominator conversion. This timing hypothesis follows the retained 256-bit/degree-7 Atkin regressions in the preceding AMD campaign; its timing benefit remains unverified.\n\nAll eight BMSS constructors at field widths 192/224/256/511 and prime degrees 7/31/101 reproduce the original and previous-rung full coefficients, codomain and kernel. Each case retains independent Kohel and public point-transport checks. All four generic Field counters are identical to the previous rung in all 96 cells. Inversion internals remain opaque. Setup, sigma, reference generation and verification stay outside constructor counts. These are synthetic supersingular field-width fixtures. Named P192/P224 whole-search and replay costs remain open.\n\nThe standalone release library, binary, integration and example tests pass; strict all-target Clippy and formatting pass. The original tested crate and continuation checkout compare byte for byte. Both comparison run hashes are checked before generating this report. The initial accidental selection of an older header-only baseline was rejected by the strict reader and remains documented in the frozen run record.\n\nThe next native wave retains separate original/candidate and previous/candidate twenty-command campaigns, each with five A/A and five interleaved A/B pairs. No timing ratio is established by this control run.";
    let typ = format!("#set page(paper: \"a4\", margin: 17mm)\n#set text(font: \"Arial\", size: 10pt)\n= Fourth BMSS candidate control\n2026-10-15 · 96/96 exact outputs; 384/384 counters unchanged\n\n#image(\"CONTROL.svg\", width: 100%)\n\n{body}\n\nSource revision:\n#text(size: 7pt)[{}]\nBinary SHA256:\n#text(size: 7pt)[{binary_hash}]\nIncremental comparison SHA256:\n#text(size: 7pt)[{}]\n", revision.trim(),sha256_hex(&raw));
    let out = Path::new(&args[2]);
    std::fs::create_dir(out).unwrap();
    for (file, bytes) in [("REPORT.typ",typ),("CONTROL.svg",svg.into()),("COUNTS.csv",csv),("REPORT.md",format!("# Fourth BMSS candidate control\n\n{body}\n\nSource `{}`; binary SHA256 `{binary_hash}`; incremental comparison SHA256 `{}`. Full raw runs, source, executable and validation logs are in `../RUNG4_RUNS.tar.gz`.\n",revision.trim(),sha256_hex(&raw))),("control.typ","#set page(width: 1000pt, height: 180pt, margin: 0pt)\n#image(\"CONTROL.svg\", width: 100%)\n".into())] { std::fs::write(out.join(file),bytes).unwrap(); }
}

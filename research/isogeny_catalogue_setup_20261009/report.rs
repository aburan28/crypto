use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{fmt::Write, fs, path::Path};
include!("visuals.rs");

fn load(path: &Path) -> Value {
    serde_json::from_slice(&fs::read(path).unwrap()).unwrap()
}
fn digest(bytes: &[u8]) -> String {
    Sha256::digest(bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}
fn entries(root: &Path, base: &Path, output: &mut Vec<Value>) {
    let mut paths: Vec<_> = fs::read_dir(root)
        .unwrap()
        .map(|p| p.unwrap().path())
        .collect();
    paths.sort();
    for path in paths {
        if path.is_dir() {
            if path.file_name().unwrap() != "target" && path.file_name().unwrap() != "dist" {
                entries(&path, base, output);
            }
        } else if path.file_name().unwrap() != "MANIFEST.json" {
            let bytes = fs::read(&path).unwrap();
            output.push(json!({"path":path.strip_prefix(base).unwrap().to_string_lossy(),"bytes":bytes.len(),"sha256":digest(&bytes)}));
        }
    }
}
fn counts(root: &Path, case: &str, build: &str) -> Vec<u64> {
    (1..=5)
        .map(|round| {
            let directory = root.join(format!("{case}-{round}-{build}"));
            assert_eq!(
                fs::read_to_string(directory.join("exit-status.txt"))
                    .unwrap()
                    .trim(),
                "0"
            );
            fs::read_to_string(directory.join("instructions.txt"))
                .unwrap()
                .trim()
                .parse()
                .unwrap()
        })
        .collect()
}
fn median(values: &[u64]) -> u64 {
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    sorted[sorted.len() / 2]
}
fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert_eq!(
        args.len(),
        3,
        "setup-report STUDY EVIDENCE or STUDY --manifest"
    );
    let study = Path::new(&args[1]);
    if args[2] == "--manifest" {
        let mut files = vec![];
        entries(study, study, &mut files);
        let mut sources = vec![];
        entries(Path::new("tools/isogeny-cli"), Path::new("."), &mut sources);
        let manifest =
            json!({"schema":"isogeny-setup-manifest/v1","files":files,"tool_sources":sources});
        fs::write(
            study.join("MANIFEST.json"),
            serde_json::to_vec_pretty(&manifest).unwrap(),
        )
        .unwrap();
        println!(
            "Manifest: {} study files, {} tool-source bindings",
            files.len(),
            sources.len()
        );
        return;
    }
    let evidence = Path::new(&args[2]);
    let reference_catalogue = load(&evidence.join("baseline-catalogue.json"));
    assert_eq!(
        reference_catalogue,
        load(&evidence.join("candidate-catalogue.json"))
    );
    assert_eq!(reference_catalogue["count"], 349);
    let mut rows = vec![];
    for case in ["screen-p224", "verify-p224", "search-p192"] {
        for round in 1..=5 {
            let a = evidence.join(format!("{case}-{round}-baseline"));
            let b = evidence.join(format!("{case}-{round}-candidate"));
            let aa = load(&a.join("stdout.json"));
            let bb = load(&b.join("stdout.json"));
            assert_eq!(aa["status"], "PASS");
            assert_eq!(bb["status"], "PASS");
            if case == "screen-p224" {
                assert!(aa["candidates"].is_array());
                assert_eq!(aa["candidates"], bb["candidates"]);
            } else {
                assert!(aa["certificate"]["records"].is_array());
                assert_eq!(aa["certificate"]["records"], bb["certificate"]["records"]);
                assert_eq!(
                    aa["certificate"]["degree_coverage"],
                    bb["certificate"]["degree_coverage"]
                );
                for record in bb["certificate"]["records"].as_array().unwrap() {
                    assert_eq!(record["exact_rational_map_check"], "PASS");
                    assert_eq!(record["kernel_check"], "PASS");
                    assert_eq!(record["public_scalar_transport_checks"], 20);
                }
                if case == "search-p192" {
                    assert_eq!(aa["degree_coverage"], "COMPLETE");
                    assert_eq!(
                        fs::read(a.join("results/p192/ell-73.json")).unwrap(),
                        fs::read(b.join("results/p192/ell-73.json")).unwrap()
                    );
                }
            }
        }
        let reference = counts(evidence, case, "baseline");
        let candidate = counts(evidence, case, "candidate");
        let am = median(&reference);
        let bm = median(&candidate);
        rows.push(json!({"case":case,"reference":reference,"candidate":candidate,"reference_median":am,"candidate_median":bm,
            "instruction_ratio":am as f64/bm as f64,"reduction_percent":100.0*(1.0-bm as f64/am as f64),
            "unit":"Callgrind Ir, ARM Linux, all traced processes","rounds":5,"correctness":"PASS","wall_time_claim":false}));
    }
    let primary = &rows[2];
    let acceptance = if primary["reduction_percent"].as_f64().unwrap() >= 1.0 {
        "PASS"
    } else {
        "NOT_MET"
    };
    let summary = json!({"schema":"isogeny-setup-profile/v1","reference":load(&evidence.join("baseline-version.json")),
        "candidate":load(&evidence.join("candidate-version.json")),"primary_acceptance":acceptance,"rows":rows});
    fs::write(
        study.join("PROFILE_SUMMARY.json"),
        serde_json::to_vec_pretty(&summary).unwrap(),
    )
    .unwrap();
    let mut report = String::from("# Native isogeny command setup and receipt hashing\n\nDated 2026-10-09. Matched 349-model reference/candidate builds, five interleaved rounds per panel. Values are guest ARM Linux instructions under the local Docker VM; native macOS throughput remains unmeasured.\n\n| Command | Reference median instructions | Candidate median instructions | Reference / candidate | Reduction |\n| --- | ---: | ---: | ---: | ---: |\n");
    for row in &rows {
        writeln!(
            report,
            "| {} | {} | {} | {:.6} | {:.6}% |",
            row["case"].as_str().unwrap(),
            row["reference_median"],
            row["candidate_median"],
            row["instruction_ratio"].as_f64().unwrap(),
            row["reduction_percent"].as_f64().unwrap()
        )
        .unwrap();
    }
    writeln!(report,"\nThe preregistered complete-search threshold was at least 1% fewer instructions; outcome: **{acceptance}**. Screening includes setup/output. Verification excludes construction. Full search includes the constructor child, executable/receipt hashing, output and independent replay. No stage ratio is an end-to-end solver gain.\n").unwrap();
    report.push_str("All catalogue JSON, normalized candidate lists, independent certificate records and complete-search map bytes match. Runtime SHA-256 matches the unchanged repository implementation at padding/block boundaries and for a 1 MiB input. The build-generated catalogue index preserves curve parameters, subgroup data and operation capabilities; the complete public inventory remains embedded. The candidate uses pinned RustCrypto SHA-256 with runtime CPU detection and a portable fallback. Constructor and independent-verifier arithmetic are unchanged.\n\nThe 349-model inventory, 174 prime-field construction models and 146 replay models remain. Montgomery/Edwards, binary and extension map adapters are still open requirements. The earlier high-degree maps and their coverage labels remain frozen. No new curve/map finding changes canonical graphs; the catalogue, cover graph and IC performance ledger remain unchanged.\n\nThe root gate still reports 643 compiler errors. Hosted CI was queued without an assigned runner at the start of this round. Merge and ReleaseMe publication require their applicable gates; local passing tests do not satisfy hosted checks.\n\n[Protocol](PROTOCOL.md), [summary](PROFILE_SUMMARY.json), [diagram](COMMAND_SETUP.svg), [PDF](REPORT.pdf), validation logs, and compact evidence archives carry the receipts. Source/executable hashes and toolchain versions identify the paired builds.\n");
    fs::write(study.join("RESULTS.md"), &report).unwrap();
    let mut ops = vec![
        Op::Text(
            38.,
            45.,
            26.,
            "Native isogeny command setup and receipt hashing".into(),
        ),
        Op::Text(
            38.,
            78.,
            16.,
            "349-model matched inputs | 5 interleaved rounds | ARM Linux Callgrind Ir".into(),
        ),
        Op::Rect(38., 106., 260., 88., "#e6edf4".into()),
        Op::Text(54., 136., 17., "Complete public catalogue".into()),
        Op::Text(54., 165., 15., "349 models; exact identities".into()),
        Op::Line(308., 150., 347., 150., "#60788d".into()),
        Op::Rect(357., 106., 300., 88., "#e4f1eb".into()),
        Op::Text(373., 136., 17., "Compact execution index".into()),
        Op::Text(373., 165., 15., "Native SHA-256; identical hashes".into()),
        Op::Line(667., 150., 706., 150., "#60788d".into()),
        Op::Rect(716., 106., 343., 88., "#e6edf4".into()),
        Op::Text(732., 136., 17., "Constructor + independent replay".into()),
        Op::Text(732., 165., 15., "Exact checks and maps preserved".into()),
        Op::Text(
            38.,
            231.,
            17.,
            "Executed instructions relative to each command's matched reference".into(),
        ),
    ];
    for (i, row) in rows.iter().enumerate() {
        let y = 279. + i as f64 * 128.;
        let fraction = row["candidate_median"].as_u64().unwrap() as f64
            / row["reference_median"].as_u64().unwrap() as f64;
        ops.push(Op::Text(38., y, 19., row["case"].as_str().unwrap().into()));
        ops.push(Op::Rect(220., y - 17., 590., 25., "#c8d3df".into()));
        ops.push(Op::Text(
            825.,
            y + 2.,
            16.,
            format!("Reference {}", row["reference_median"]),
        ));
        ops.push(Op::Rect(
            220.,
            y + 18.,
            590. * fraction,
            25.,
            "#438d75".into(),
        ));
        ops.push(Op::Text(
            825.,
            y + 37.,
            16.,
            format!("Candidate {}", row["candidate_median"]),
        ));
        ops.push(Op::Text(
            220.,
            y + 66.,
            15.,
            format!(
                "Ratio {:.6}; reduction {:.6}%",
                row["instruction_ratio"].as_f64().unwrap(),
                row["reduction_percent"].as_f64().unwrap()
            ),
        ));
    }
    ops.push(Op::Text(
        38.,
        697.,
        17.,
        format!("Complete-search acceptance: {acceptance}. Exact output comparisons: PASS."),
    ));
    ops.push(Op::Text(
        38.,
        732.,
        15.,
        "Catalogue adapters, root compilation, hosted CI and publication gaps remain explicit."
            .into(),
    ));
    fs::write(study.join("COMMAND_SETUP.svg"), svg(&ops)).unwrap();
    let mut page = String::new();
    let mut y = 746.;
    text(
        &mut page,
        42.,
        y,
        19.,
        "Isogeny setup and hashing continuation",
    );
    y -= 34.;
    for line in ["Matched 349-model inputs; Rust 1.93.1; portable release flags; native ARM Linux under the local Docker VM.",
        "Five interleaved rounds. Callgrind Ir is summed across parent and constructor child. Native macOS throughput is unmeasured."] {
        wrapped(&mut page,&mut y,line,11.);
    }
    for row in &rows {
        wrapped(
            &mut page,
            &mut y,
            &format!(
                "{}: reference {} instructions; candidate {}; ratio {:.6}; reduction {:.6}%.",
                row["case"].as_str().unwrap(),
                row["reference_median"],
                row["candidate_median"],
                row["instruction_ratio"].as_f64().unwrap(),
                row["reduction_percent"].as_f64().unwrap()
            ),
            11.,
        );
    }
    for line in [format!("Preregistered full-search acceptance: {acceptance}. The threshold is at least 1% fewer complete-command instructions."),
        "Catalogue JSON, candidate records, independent certificates and full-search map bytes match in every round. No exact mathematical check was removed.".into(),
        "Receipt hashes are checked against the unchanged SHA-256 implementation. Runtime CPU detection preserves a portable fallback. The compact index retains exact curve/subgroup inputs.".into(),
        "Inventory and support remain 349 models, 174 prime construction models and 146 replay models. Montgomery/Edwards, binary and extension map adapters remain open.".into(),
        "Root compilation retains 643 errors; hosted CI, merge and ReleaseMe publication remain pending. No new map changes canonical graphs or the IC performance ledger.".into(),
        "Receipts: PROTOCOL.md, PROFILE_SUMMARY.json, MANIFEST.json, validation logs and compact raw evidence archives. The next page contains the editable vector diagram.".into()] {
        wrapped(&mut page,&mut y,&line,11.);
    }
    emit_pdf(
        &study.join("REPORT.pdf"),
        vec![(595, 842, page), (842, 595, vector_pdf(&ops))],
    );
    println!("All 15 paired output comparisons PASS; primary acceptance {acceptance}.");
}

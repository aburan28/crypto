//! Render only frozen, native-reader-validated kernel baseline records.
use isogeny_algos::{json::Json, sha256::sha256_hex};
use std::{fmt::Write as _, path::Path};

fn field<'a>(j: &'a Json, k: &str) -> Result<&'a str, String> {
    j.get(k).and_then(Json::as_str).ok_or_else(|| format!("missing {k}"))
}
fn main() -> Result<(), String> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 5 {
        return Err("make-kernel-report BASELINE.jsonl REPEAT.jsonl COMPARISON.json FRESH_DIR".into());
    }
    let baseline = std::fs::read(&args[1]).map_err(|e| e.to_string())?;
    let repeat = std::fs::read(&args[2]).map_err(|e| e.to_string())?;
    let comparison_bytes = std::fs::read(&args[3]).map_err(|e| e.to_string())?;
    let comparison = Json::parse(std::str::from_utf8(&comparison_bytes).map_err(|e| e.to_string())?)?;
    if field(&comparison, "status")? != "PASS"
        || comparison.get("exact_output_matches").and_then(Json::as_f64) != Some(84.0)
        || field(&comparison, "baseline_run_sha256")? != sha256_hex(&baseline)
        || field(&comparison, "candidate_run_sha256")? != sha256_hex(&repeat)
    {
        return Err("missing complete, hash-linked 84-case comparison".into());
    }
    let Json::Arr(cases) = comparison.get("cases").ok_or("missing cases")? else {
        return Err("invalid cases".into());
    };
    for case in cases {
        for op in ["add_sub", "mul", "inv", "neg"] {
            let values = case.get("operations").and_then(|j| j.get(op)).ok_or("missing counts")?;
            if field(values, "baseline")? != field(values, "candidate")? {
                return Err("repeat counts differ".into());
            }
        }
    }
    let header = Json::parse(std::str::from_utf8(&baseline).map_err(|e| e.to_string())?.lines().next().ok_or("empty run")?)?;
    let source = field(&header, "source_revision_supplied_by_launcher")?;
    let binary = field(&header, "binary_sha256")?;
    if source.len() != 40 || !source.bytes().all(|b| b.is_ascii_hexdigit()) {
        return Err("invalid source".into());
    }
    let names = ["velu_from_point", "velu_cyclic", "velu_xonly", "velu_xonly_fast", "sqrt_velu", "sqrt_velu_fast", "kohel_given_kernel"];
    let mut selected = vec![];
    for name in names {
        let row = cases.iter().find(|r| {
            r.get("bits").and_then(Json::as_str) == Some("224")
                && r.get("ell").and_then(Json::as_str) == Some("101")
                && r.get("algorithm").and_then(Json::as_str) == Some(name)
        }).ok_or_else(|| format!("missing selected cell {name}"))?;
        let count = |op| -> Result<u64, String> {
            field(row.get("operations").and_then(|j| j.get(op)).ok_or("missing operation")?, "baseline")?
                .parse().map_err(|_| "invalid integer count".into())
        };
        let input = match name {
            "kohel_given_kernel" => "kernel polynomial",
            "velu_xonly" | "velu_xonly_fast" | "sqrt_velu_fast" => "generator x",
            _ => "affine generator",
        };
        selected.push((name, input, count("mul")?, count("inv")?));
    }
    let mut svg = String::from("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1100\" height=\"430\" viewBox=\"0 0 1100 430\"><rect width=\"1100\" height=\"430\" fill=\"white\"/><g font-family=\"Arial,sans-serif\" fill=\"#17324d\"><text x=\"20\" y=\"30\" font-size=\"21\">Kernel constructor baseline · 224 field bits / degree 101</text><text x=\"20\" y=\"55\" font-size=\"14\">Generic Field calls; inversion internals opaque; n=2 identical recorded counts; shared Apple M4 Pro</text><text x=\"310\" y=\"83\" font-size=\"16\">Multiplication calls</text><text x=\"755\" y=\"83\" font-size=\"16\">Inversion calls</text>");
    let max_mul = selected.iter().map(|r| r.2).max().unwrap().max(1) as f64;
    let max_inv = selected.iter().map(|r| r.3).max().unwrap().max(1) as f64;
    let mut csv = String::from("method,supplied_input,multiplication_calls,inversion_calls,repeated_runs\n");
    let mut table = String::from("| Method | Supplied input | Multiplication calls | Inversion calls |\n|---|---|---:|---:|\n");
    let mut typst_rows = String::new();
    for (index, (name, input, mul, inv)) in selected.iter().enumerate() {
        let y = 104 + 39 * index;
        let color = if *name == "kohel_given_kernel" { "#ad742b" } else { "#2b7188" };
        writeln!(svg, "<text x=\"20\" y=\"{}\" font-size=\"15\">{name}</text><text x=\"20\" y=\"{}\" font-size=\"12\">{input}</text><rect x=\"310\" y=\"{y}\" width=\"{:.2}\" height=\"18\" fill=\"{color}\"/><text x=\"690\" y=\"{}\" text-anchor=\"end\" font-size=\"14\">{mul}</text><rect x=\"755\" y=\"{y}\" width=\"{:.2}\" height=\"18\" fill=\"{color}\"/><text x=\"1080\" y=\"{}\" text-anchor=\"end\" font-size=\"14\">{inv}</text>", y + 12, y + 27, 300.0 * *mul as f64 / max_mul, y + 14, 240.0 * *inv as f64 / max_inv, y + 14).unwrap();
        writeln!(csv, "{name},{input},{mul},{inv},2").unwrap();
        writeln!(table, "| {name} | {input} | {mul} | {inv} |").unwrap();
        writeln!(typst_rows, "[{name}], [{input}], [{mul}], [{inv}],").unwrap();
    }
    svg.push_str("<text x=\"310\" y=\"392\" font-size=\"12\">Bars start at zero; each panel uses its own linear count scale.</text><text x=\"20\" y=\"418\" font-size=\"13\">Amber: supplied kernel polynomial. Different input boundaries; no end-to-end or wall-time ratio is admitted.</text></g></svg>\n");
    let pipeline = "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1100\" height=\"160\"><rect width=\"1100\" height=\"160\" fill=\"white\"/><g font-family=\"Arial,sans-serif\" fill=\"#17324d\"><rect x=\"10\" y=\"15\" width=\"300\" height=\"90\" fill=\"#edf5fa\" stroke=\"#6c91ac\"/><text x=\"25\" y=\"45\" font-size=\"18\">Supplied inputs kept distinct</text><text x=\"25\" y=\"74\" font-size=\"15\">Affine generator / x / kernel polynomial</text><text x=\"332\" y=\"65\" font-size=\"24\">→</text><rect x=\"370\" y=\"15\" width=\"310\" height=\"90\" fill=\"#edf5fa\" stroke=\"#6c91ac\"/><text x=\"385\" y=\"45\" font-size=\"18\">7 constructors · 4 widths · 3 degrees</text><text x=\"385\" y=\"74\" font-size=\"15\">Save counts before verification</text><text x=\"702\" y=\"65\" font-size=\"24\">→</text><rect x=\"740\" y=\"15\" width=\"350\" height=\"90\" fill=\"#edf5fa\" stroke=\"#6c91ac\"/><text x=\"755\" y=\"45\" font-size=\"18\">Kohel + Vélu · 2ell+1 point transports</text><text x=\"755\" y=\"74\" font-size=\"15\">84/84 repeated fingerprints and counts</text><text x=\"15\" y=\"140\" font-size=\"14\">Fixture/reference construction and verification excluded. Native representation allocation included. No wall-time admission.</text></g></svg>\n";
    let report = format!("# Seven kernel constructor baselines\n\n2026-10-15. Source `{source}`. Two complete 84-case runs: seven methods, four field widths (192/224/256/511), degrees 7/31/101. The native reader accepted both runs and 84 paired output fingerprints. All four counted operations repeat exactly in every cell.\n\n## Selected measured counts\n\n224 field bits, degree 101; generic Field API calls per construction, inversion internals opaque.\n\n{table}\n![Counts](COUNTS.svg)\n\n![Boundary](PIPELINE.svg)\n\nThese are baseline counts, not optimization gains. Supplied-kernel and supplied-generator rows have distinct boundaries. Field/fixture generation, references, fingerprints and verification are excluded. Owned result allocation is included in the constructor.\n\nCorrectness: degree, codomain and full point transports at 2*ell+1 distinct non-kernel abscissae agree with both Kohel and point-enumerating Vélu. Counted/native outputs agree. These fixture fingerprints are distinct from the BMSS series' full coefficient comparison. The curves are synthetic field-width fixtures, not named NIST models.\n\nShared Apple M4 Pro, macOS ARM64, 14 logical CPUs, 48 GiB. No timing samples requested; zero qualified wall samples. Other physical CPU baselines and kernel optimization rounds remain open.\n\nBinary SHA256: `{binary}`. Comparison SHA256: `{}`. Full raw runs, source and executable are retained in the neighboring evidence archive; see KERNEL_BASELINE_PROTOCOL.md.\n", sha256_hex(&comparison_bytes));
    let typst = format!("#set page(paper: \"a4\", margin: 16mm)\n#set text(font: \"Arial\", size: 8.5pt)\n#set heading(numbering: none)\n= Seven kernel constructor baselines\n2026-10-15 · Two repeated native 84-case runs\n\nSeven methods at 192/224/256/511 field bits, degrees 7/31/101. All 84 paired output fingerprints and all four field-call counters repeat exactly. Synthetic supersingular fixtures; these are not the named NIST curve models.\n\n== Measured calls at 224 bits / degree 101\n#table(columns: (1.6fr, 1.3fr, 1fr, 0.8fr), inset: 4pt, [Method], [Supplied input], [Mul calls], [Inv calls],\n{typst_rows})\n#image(\"COUNTS.svg\", width: 100%)\n#image(\"PIPELINE.svg\", width: 100%)\n\nGeneric Field API calls, inversion internals opaque; n=2 identical recorded counts. Different supplied-input boundaries remain separate. Reference/fixture generation, fingerprints and verification are excluded; owned representation allocation is included. These baseline counts establish no kernel speedup.\n\nDegree, codomain and full point transports at 2ell+1 distinct non-kernel abscissae agree with both Kohel and point-enumerating Vélu. Counts are saved before verification. These fixture fingerprints are distinct from the BMSS coefficient comparison.\n\nShared Apple M4 Pro, 14 logical CPUs, 48 GiB, macOS ARM64. Zero timing samples requested or qualified. Other physical CPU baselines and optimization rounds remain open.\n\nSource: #text(size: 7pt)[{source}]\nBinary SHA256: #text(size: 6.5pt)[{binary}]\nComparison SHA256: #text(size: 6.5pt)[{}]\n\nRead REPORT.md, COUNTS.csv and KERNEL_BASELINE_PROTOCOL.md with the frozen raw records and source/binary archive. Both complete matrices, including all other widths and degrees, remain in the raw comparison.\n", sha256_hex(&comparison_bytes));
    let output_note = "Kohel's result contains all rational-map coefficients. The other constructor results store compact evaluation representations. Their build counts include different output materialization costs.";
    let report = format!("{report}\n{output_note}\n");
    let typst = format!("{typst}\n\n{output_note}\n");
    let dir = Path::new(&args[4]);
    std::fs::create_dir(dir).map_err(|e| e.to_string())?;
    for (name, bytes) in [
        ("COUNTS.svg", svg), ("PIPELINE.svg", pipeline.into()), ("COUNTS.csv", csv),
        ("REPORT.md", report), ("REPORT.typ", typst),
        ("counts.typ", "#set page(width: 1100pt, height: 430pt, margin: 0pt)\n#image(\"COUNTS.svg\", width: 1100pt)\n".into()),
        ("pipeline.typ", "#set page(width: 1100pt, height: 160pt, margin: 0pt)\n#image(\"PIPELINE.svg\", width: 1100pt)\n".into()),
    ] {
        std::fs::write(dir.join(name), bytes).map_err(|e| e.to_string())?;
    }
    Ok(())
}

//! Render only a complete native artifact audit and its hash-linked timing analysis.
use isogeny_algos::{json::Json, sha256::sha256_hex};
use std::{fmt::Write as _, path::Path};
fn s<'a>(j: &'a Json, k: &str) -> Result<&'a str, String> {
    j.get(k)
        .and_then(Json::as_str)
        .ok_or_else(|| format!("missing {k}"))
}
fn n(j: &Json, k: &str) -> Result<f64, String> {
    j.get(k)
        .and_then(Json::as_f64)
        .filter(|x| x.is_finite() && *x > 0.0)
        .ok_or_else(|| format!("invalid {k}"))
}
fn read(path: &Path) -> Result<Json, String> {
    Json::parse(&std::fs::read_to_string(path).map_err(|e| e.to_string())?)
}
fn main() -> Result<(), String> {
    let a: Vec<_> = std::env::args().collect();
    if a.len() != 3 {
        return Err("make-cpu-report ARTIFACT_ROOT FRESH_DIR".into());
    }
    let root = Path::new(&a[1]);
    let audit = read(&root.join("ubuntu-24.04-audit/AUDIT.json"))?;
    let bytes = std::fs::read(root.join("ubuntu-24.04-audit/CPU_ANALYSIS.json"))
        .map_err(|e| e.to_string())?;
    let analysis = Json::parse(std::str::from_utf8(&bytes).map_err(|e| e.to_string())?)?;
    if s(&audit, "status")? != "PASS"
        || audit.get("wall_speedup_eligible") != Some(&Json::Bool(true))
        || audit.get("qualified_timing_runs").and_then(Json::as_f64) != Some(20.0)
        || s(&analysis, "status")? != "PASS"
        || s(&audit, "summary_sha256")? != s(&analysis, "source_summary_sha256")?
    {
        return Err("missing complete, independently audited timing admission".into());
    }
    let Some(Json::Arr(rows)) = analysis.get("qualified_wall_results") else {
        return Err("missing results".into());
    };
    if rows.len() != 64 {
        return Err("missing full timing matrix".into());
    }
    let mut csv = String::from("bits,ell,algorithm,paired_median,paired_min,paired_max,aa_min,aa_max,observed_noise_envelope,median_exceeds_aa_envelope\n");
    for r in rows {
        writeln!(
            csv,
            "{},{},{},{},{},{},{},{},{},{}",
            s(r, "bits")?,
            s(r, "ell")?,
            s(r, "algorithm")?,
            n(r, "ab_median")?,
            n(r, "ab_min")?,
            n(r, "ab_max")?,
            n(r, "aa_min")?,
            n(r, "aa_max")?,
            n(r, "observed_noise_envelope")?,
            r.get("median_exceeds_aa_envelope")
                .ok_or("missing noise flag")?
                .dump()
        )
        .unwrap();
    }
    let names = [
        "linear_algebra",
        "stark1972",
        "atkin1992",
        "atkin_modcomp",
        "elkies1992",
        "elkies1998",
        "fast_elkies",
        "fast_elkies_prime",
    ];
    let mut table = String::from("| Method | Paired median ratio | Five-round range | Observed A/A range |\n|---|---:|---:|---:|\n");
    let mut typst_rows = String::new();
    let x = |ratio: f64| 220.0 + (ratio.log2() + 1.0) * 700.0 / 6.0;
    let mut svg = String::from("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1100\" height=\"475\"><rect width=\"1100\" height=\"475\" fill=\"white\"/><g font-family=\"Arial,sans-serif\" fill=\"#17324d\"><text x=\"20\" y=\"28\" font-size=\"21\">AMD EPYC 7763 · 224-bit fixture · prime degree 31</text><text x=\"20\" y=\"53\" font-size=\"14\">Median of five paired round ratios; whiskers = observed min/max; grey = A/A range; log2 ratio axis</text>");
    for tick in [0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0] {
        writeln!(svg,"<line x1=\"{:.2}\" y1=\"75\" x2=\"{:.2}\" y2=\"395\" stroke=\"#d9e0e5\"/><text x=\"{:.2}\" y=\"420\" text-anchor=\"middle\" font-size=\"14\">{tick}</text>",x(tick),x(tick),x(tick)).unwrap();
    }
    for (i, name) in names.iter().enumerate() {
        let r = rows
            .iter()
            .find(|r| {
                r.get("bits").and_then(Json::as_str) == Some("224")
                    && r.get("ell").and_then(Json::as_str) == Some("31")
                    && r.get("algorithm").and_then(Json::as_str) == Some(*name)
            })
            .ok_or("missing selected cell")?;
        let (m, lo, hi, aa_lo, aa_hi) = (
            n(r, "ab_median")?,
            n(r, "ab_min")?,
            n(r, "ab_max")?,
            n(r, "aa_min")?,
            n(r, "aa_max")?,
        );
        let y = 94 + 39 * i;
        writeln!(svg,"<text x=\"20\" y=\"{}\" font-size=\"15\">{name}</text><rect x=\"{:.2}\" y=\"{}\" width=\"{:.2}\" height=\"24\" fill=\"#bfc7ce\"/><line x1=\"{:.2}\" y1=\"{y}\" x2=\"{:.2}\" y2=\"{y}\" stroke=\"#226e88\" stroke-width=\"3\"/><circle cx=\"{:.2}\" cy=\"{y}\" r=\"5\" fill=\"#226e88\"/><text x=\"1085\" y=\"{}\" text-anchor=\"end\" font-size=\"14\">{m:.3} [{lo:.3}, {hi:.3}]</text>",y+5,x(aa_lo),y-12,(x(aa_hi)-x(aa_lo)).max(1.0),x(lo),x(hi),x(m),y+5).unwrap();
        writeln!(
            table,
            "| {name} | {m:.4} | {lo:.4}–{hi:.4} | {aa_lo:.4}–{aa_hi:.4} |"
        )
        .unwrap();
        writeln!(
            typst_rows,
            "[{name}], [{m:.3}], [{lo:.3}–{hi:.3}], [{aa_lo:.3}–{aa_hi:.3}],"
        )
        .unwrap();
    }
    svg.push_str("<text x=\"20\" y=\"451\" font-size=\"14\">Unit: baseline/candidate constructor sample-median ratio. Setup, sigma and verification excluded.</text><text x=\"20\" y=\"471\" font-size=\"13\">Five A/A + five interleaved A/B pairs; 20/20 native guest-isolated runs. Other 56 cells, including slower rows, retained.</text></g></svg>\n");
    let pipeline="<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1100\" height=\"150\"><rect width=\"1100\" height=\"150\" fill=\"white\"/><g font-family=\"Arial,sans-serif\" fill=\"#17324d\"><rect x=\"10\" y=\"10\" width=\"300\" height=\"80\" fill=\"#edf5fa\" stroke=\"#6c91ac\"/><text x=\"25\" y=\"39\" font-size=\"17\">Frozen curves, codomains and sigma</text><text x=\"25\" y=\"67\" font-size=\"15\">Setup / independent reference excluded</text><text x=\"330\" y=\"60\" font-size=\"24\">→</text><rect x=\"370\" y=\"10\" width=\"310\" height=\"80\" fill=\"#edf5fa\" stroke=\"#6c91ac\"/><text x=\"385\" y=\"39\" font-size=\"17\">Native constructor · 5 samples per run</text><text x=\"385\" y=\"67\" font-size=\"15\">1 warm call; reserved guest CPU</text><text x=\"700\" y=\"60\" font-size=\"24\">→</text><rect x=\"740\" y=\"10\" width=\"350\" height=\"80\" fill=\"#edf5fa\" stroke=\"#6c91ac\"/><text x=\"755\" y=\"39\" font-size=\"17\">Exact map / kernel + point checks</text><text x=\"755\" y=\"67\" font-size=\"15\">Verify outside timed interval; retain failures</text><text x=\"20\" y=\"127\" font-size=\"14\">Every raw pair median, binary/source identity and isolation command rechecked. VM host neighbours/frequency unobserved.</text></g></svg>\n";
    let platforms="| Platform / recorded CPU | Full output matrix | Repeated valid runs | Qualified wall runs |\n|---|---:|---:|---:|\n| Linux x86-64 / AMD EPYC 7763 | 96/96 | 20/20 | 20/20 |\n| Linux ARM64 / implementer 0x41, part 0xd49 | 96/96 | 19/20 | CENSORED |\n| macOS ARM64 / Apple M1 (Virtual) | 96/96 | 20/20 | 0 |\n| macOS x86-64 / Intel Core i7-8700B | 96/96 | 20/20 | 0 |\n";
    let scope="Constructor-only measurements on synthetic supersingular field-width fixtures. Preparing sigma, fixture/reference generation, fingerprints and verification are excluded. These are not named P192/P224 whole-search or replay speedups. The AMD guest has four logical CPUs; the timed child is reserved and pinned by the native tool. Guest contention checks pass for all twenty commands. VM host neighbours and frequency remain unobserved. A/A envelopes and min/max ranges are descriptive, not confidence intervals.";
    let regressions="All 64 cells remain in RESULTS.csv and the raw analysis. At 256 bits/degree 7, Atkin 1992 has a paired median ratio 0.984707 and Atkin modular composition 0.970765; their A/A ranges are 0.994689–1.001896 and 0.993324–1.016337. These slower rows are retained. Several other small-case changes are inside the observed A/A envelope. No blanket per-method or cross-CPU win is inferred.";
    let arm="The ARM first A/A run was refused before execution: runner housekeeping consumed 0.29 CPU-seconds in a two-second settling window (limit 0.20), with PSI below its limit. The producer completed nineteen subsequent admitted runs. The requested complete campaign remains censored; no partial ARM timing ratio is promoted. Both macOS profiles have exact repeated counts and outputs but no admitted wall samples.";
    let report=format!("# First native CPU campaign\n\n2026-10-15. GitHub run [38766895864](https://github.com/aburan28/crypto/actions/runs/38766895864), baseline `79006eb85a400474e83f6a76fac0dca6b60a8b5d`, candidate `a21306ca20418a4c4bfcc8448de9b8851b8449ed` (first three BMSS rungs). Native compiler receipts and actual processors are retained.\n\n{platforms}\n\n## Admitted AMD 224-bit / degree-31 ratios\n\nMedian of five baseline/candidate paired round ratios. Each arm's round median uses five constructor samples.\n\n{table}\n![Ratios](RATIOS.svg)\n\n![Boundary](PIPELINE.svg)\n\n{scope}\n\n{regressions}\n\n{arm}\n\nThe native artifact auditor checks actual retained binary hashes, supplied source revisions, complete input/output matrices, repeated counts, isolation command/argument binding and every admitted sample median. It writes separate recheck records; originals are unchanged. The analysis SHA256 is `{}`. Linux legacy-family outputs remain exploratory: they lack complete result fingerprints and A/A controls, and keep every verification flag.\n",sha256_hex(&bytes));
    let typst=format!("#set page(paper: \"a4\", margin: 16mm)\n#set text(font: \"Arial\", size: 8.5pt)\n#set heading(numbering: none)\n= First native CPU timing campaign\n2026-10-15 · Run 38766895864 · AMD EPYC 7763\n\nEight BMSS methods, four field widths (192/224/256/511), degrees 7/31: 64 timed cells. Five A/A pairs and five interleaved A/B pairs, five constructor samples per arm per round. All twenty AMD commands pass native guest-isolation checks and independent artifact audit. Full 96-case coefficient comparisons pass on all four CPU profiles at degrees 7/31/101.\n\n== 224-bit fixture / degree 31\n#table(columns: (1.7fr, 0.8fr, 1.3fr, 1.3fr), inset: 3pt, [Method], [Median ratio], [Five-round range], [A/A range],\n{typst_rows})\n#image(\"RATIOS.svg\", width: 100%)\n#image(\"PIPELINE.svg\", width: 100%)\n\n{scope}\n\n{regressions}\n\n{arm}\n\nBaseline 79006eb85a40 → candidate a21306ca2041. Raw ZIPs, all 64 rows, medians, source/binary/isolation audits and failures remain in the linked evidence. Analysis SHA256:\n#text(size: 6.5pt)[{}]\n",sha256_hex(&bytes));
    let out = Path::new(&a[2]);
    std::fs::create_dir(out).map_err(|e| e.to_string())?;
    for (name,body) in [("REPORT.md",report),("REPORT.typ",typst),("RATIOS.svg",svg),("PIPELINE.svg",pipeline.into()),("RESULTS.csv",csv),("ratios.typ","#set page(width: 1100pt, height: 475pt, margin: 0pt)\n#image(\"RATIOS.svg\", width: 1100pt)\n".into()),("pipeline.typ","#set page(width: 1100pt, height: 150pt, margin: 0pt)\n#image(\"PIPELINE.svg\", width: 1100pt)\n".into())] {std::fs::write(out.join(name),body).map_err(|e| e.to_string())?;}
    Ok(())
}

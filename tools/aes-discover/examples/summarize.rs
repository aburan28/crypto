//! Render recorded observations. This never launches a scientific trial or
//! computes a new significance test. Original JSON remains the source of truth.
use aes_discover::digest;
use serde_json::{json, Value};
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::Path;

fn write_new(path: &Path, content: &str) -> Result<(), Box<dyn std::error::Error>> {
    let mut f = OpenOptions::new().write(true).create_new(true).open(path)?;
    f.write_all(content.as_bytes())?;
    Ok(())
}
fn xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}
fn bars(title: &str, subtitle: &str, rows: &[(String, u64)], maximum: u64) -> String {
    let height = 110 + rows.len() * 30;
    let mut s = format!(
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="820" height="{height}" viewBox="0 0 820 {height}"><rect width="820" height="{height}" fill="#f8fafc"/><g font-family="sans-serif" fill="#0f172a"><text x="22" y="28" font-size="20" font-weight="bold">{}</text><text x="22" y="51" font-size="12">{}</text>"##,
        xml(title),
        xml(subtitle)
    );
    for (i, (label, n)) in rows.iter().enumerate() {
        let y = 72 + i * 30;
        let width = 450.0 * (*n as f64) / (maximum.max(1) as f64);
        s.push_str(&format!(r##"<text x="22" y="{}" font-size="12">{}</text><rect x="270" y="{y}" width="{width:.2}" height="18" rx="3" fill="#2563eb"/><text x="{}" y="{}" font-size="12">{n}</text>"##,y+14,xml(label),278.0+width,y+14));
    }
    s.push_str("</g></svg>");
    s
}

fn run() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args_os().collect();
    if args.len() != 3 {
        return Err("usage: summarize RUN-DIRECTORY NEW-REPORT-DIRECTORY".into());
    }
    let input = Path::new(&args[1]);
    let output = Path::new(&args[2]);
    fs::create_dir(output)?;
    let mut paths: Vec<_> = fs::read_dir(input)?
        .map(|r| r.map(|x| x.path()))
        .collect::<Result<_, _>>()?;
    paths.sort();
    let mut families = Vec::new();
    let mut integrals = Vec::new();
    let mut bindings = Vec::new();
    let mut differential_bars = Vec::new();
    let mut integral_bars = Vec::new();
    let mut frozen_event_max = 0;
    let mut structure_max = 0;
    for path in &paths {
        let name = path.file_name().unwrap().to_string_lossy();
        if name.ends_with("-validation.json")
            || (name.starts_with("integral") && name.ends_with(".json"))
        {
            let bytes = fs::read(path)?;
            let v: Value = serde_json::from_slice(&bytes)?;
            bindings.push(json!({"file":name,"sha256":digest(&bytes)}));
            if name.ends_with("-validation.json") {
                let count = v["validated_candidates"]
                    .as_array()
                    .ok_or("missing candidate list")?
                    .len() as u64;
                let label = name.trim_end_matches("-validation.json").to_string();
                let spec_name = format!("{label}-discovery.json");
                let spec_bytes = fs::read(input.join(&spec_name))?;
                if v["spec_sha256"].as_str() != Some(digest(&spec_bytes).as_str()) {
                    return Err(format!("{name}: discovery artifact digest mismatch").into());
                }
                let spec: Value = serde_json::from_slice(&spec_bytes)?;
                let frozen = spec["candidates"]
                    .as_array()
                    .ok_or("missing frozen events")?
                    .len();
                frozen_event_max = frozen_event_max.max(frozen as u64);
                bindings.push(json!({"file":spec_name,"sha256":digest(&spec_bytes)}));
                differential_bars.push((label.clone(), count));
                families.push(json!({"label":label,"cipher":v["cipher"],"keys":v["validation_plan"]["keys"],
                    "pairs_per_event_per_key":v["validation_plan"]["samples_per_key_per_event"],
                    "alpha":v["validation_plan"]["alpha"],"primary_oracle":v["validation_plan"]["primary_oracle"],
                    "frozen_events":frozen,"searched_input_differences":spec["searched_input_differences"],
                    "discovery_keys":spec["discovery_key_sha256"].as_array().ok_or("missing discovery keys")?.len(),
                    "validated_events":count,"control_pass":v["random_control_passed"]}));
            } else {
                let rows = v["rows"].as_array().ok_or("missing integral rows")?;
                let positive = rows
                    .iter()
                    .filter(|r| r["balanced_output_bytes"] == 16)
                    .count() as u64;
                let null_count = rows
                    .iter()
                    .filter(|r| r["permutation_balanced_output_bytes"] == 16)
                    .count();
                let label = name.trim_end_matches(".json").to_string();
                structure_max = structure_max.max(rows.len() as u64);
                integral_bars.push((label.clone(), positive));
                integrals.push(
                    json!({"label":label,"structures":rows.len(),"all_bytes_balanced":positive,
                    "all_bytes_balanced_in_permutation":null_count,"cipher":v["cipher"]}),
                );
            }
        }
    }
    let qb = fs::read(input.join("quartets.json"))?;
    let quartet: Value = serde_json::from_slice(&qb)?;
    bindings.push(json!({"file":"quartets.json","sha256":digest(&qb)}));
    let accepted = quartet["accepted_nondegenerate_quartets"]
        .as_u64()
        .ok_or("missing quartet count")?;
    let quartet_bars = vec![
        ("Accepted nondegenerate quartets".into(), accepted),
        (
            "Inversion relation matches".into(),
            quartet["key_wrapped_inversion_matches"].as_u64().unwrap(),
        ),
        (
            "Raw AES S-box relation matches".into(),
            quartet["raw_aes_sbox_matches"].as_u64().unwrap(),
        ),
        (
            "Output affine undone matches".into(),
            quartet["aes_output_linear_affine_undone_matches"]
                .as_u64()
                .unwrap(),
        ),
        (
            "Random byte bijection matches".into(),
            quartet["random_bijection_matches"].as_u64().unwrap(),
        ),
    ];
    let status = fs::read_to_string(input.join("status.tsv"))?;
    let rows: Vec<_> = status.lines().collect();
    let successful = rows
        .iter()
        .filter(|line| line.split('\t').nth(1) == Some("0"))
        .count();
    let environment = fs::read_to_string(input.join("environment.txt"))?;
    let source_revision = environment
        .lines()
        .find(|line| line.len() == 40 && line.bytes().all(|c| c.is_ascii_hexdigit()))
        .ok_or("missing source revision in environment receipt")?;
    bindings.push(json!({"file":"environment.txt","sha256":digest(environment.as_bytes())}));
    bindings.push(json!({"file":"status.tsv","sha256":digest(status.as_bytes())}));
    let all_controls_passed =
        !families.is_empty() && families.iter().all(|f| f["control_pass"] == true);
    let alpha_sum: f64 = families
        .iter()
        .map(|f| f["alpha"].as_f64().ok_or("missing family alpha"))
        .collect::<Result<Vec<_>, _>>()?
        .iter()
        .sum();
    let positive_labels: Vec<_> = families
        .iter()
        .filter(|f| f["validated_events"].as_u64().unwrap_or(0) > 0)
        .map(|f| f["label"].as_str().unwrap_or("unlabeled"))
        .collect();
    let observations = format!("{} recorded differential families; {} families contain an event passing every tested key. Families with validated events: {}. All matched controls passed: {all_controls_passed}. The sum of recorded per-family alpha allocations is {alpha_sum:.6}. This allocation is meaningful only when frozen before execution; retrospective collection does not create a campaign error guarantee.",
        families.len(), positive_labels.len(), if positive_labels.is_empty() { "none".to_string() } else { positive_labels.join(", ") });
    let summary = json!({"schema":"aes-discover/recorded-summary/v1","run_directory":input.to_string_lossy(),
        "source_revision":source_revision,"all_controls_passed":all_controls_passed,"recorded_alpha_sum":alpha_sum,
        "recorded_processes":rows.len(),"successful_processes":successful,"families":families,"integrals":integrals,
        "quartet_component":quartet,"source_bindings":bindings,
        "scope":"descriptive counts derived from recorded JSON only; no new trials, statistical tests, scientific promotion, speed comparison or full-AES claim"});
    write_new(
        &output.join("summary.json"),
        &serde_json::to_string_pretty(&summary)?,
    )?;
    write_new(
        &output.join("differential-events.svg"),
        &bars(
            "Frozen differential events passing every tested key",
            "Recorded per-family counts; frozen event, key, sample and alpha settings are in summary.json",
            &differential_bars,
            frozen_event_max,
        ),
    )?;
    write_new(
        &output.join("integral-structures.svg"),
        &bars(
            "Structures with all 16 output bytes XOR-balanced",
            "256-text structures; per-configuration structure counts are in summary.json",
            &integral_bars,
            structure_max,
        ),
    )?;
    write_new(
        &output.join("quartet-component.svg"),
        &bars(
            "Cross-ratio: component identity and affine-layer ablation",
            &format!(
                "{} generated quartets; pole cases excluded; raw match counts",
                quartet["generated_quartets"]
            ),
            &quartet_bars,
            accepted,
        ),
    )?;
    let mut table=String::from("| Family | Keys | Pairs/event/key | Validated events | PRP control passed |\n| --- | ---: | ---: | ---: | --- |\n");
    for f in summary["families"].as_array().unwrap() {
        table.push_str(&format!(
            "| {} | {} | {} | {} | {} |\n",
            f["label"].as_str().unwrap(),
            f["keys"],
            f["pairs_per_event_per_key"],
            f["validated_events"],
            f["control_pass"]
        ));
    }
    let component_observations = format!("Of {accepted} accepted nondegenerate quartets, inversion matches: {}; raw AES S-box matches: {}; affine-undone matches: {}; random byte-bijection matches: {}. These are local component observations, not an invariant through an AES round or a reproduction of a complete published attack.",
        quartet["key_wrapped_inversion_matches"], quartet["raw_aes_sbox_matches"],
        quartet["aes_output_linear_affine_undone_matches"], quartet["random_bijection_matches"]);
    let md=format!("# Native AES laboratory: recorded suite\n\nStatus: {successful}/{} recorded commands exited successfully. Raw directory: `{}`.\n\n## Differential endpoints\n\n{table}\n{observations}\n\nThese observations establish neither security nor a novel attack. Per-family frozen event counts, discovery-key counts, input-difference counts and validation plans are in `summary.json`. Raw reports retain every key/event and the conservative Chernoff/PRP-coupling/Bonferroni bounds. A failed control withholds the affected family's candidates.\n\n![Recorded differential counts](differential-events.svg)\n\n## Component controls\n\n![Recorded integral counts](integral-structures.svg)\n\n![Recorded quartet counts](quartet-component.svg)\n\n{component_observations}\n\n## Provenance and limits\n\nThe recorded source revision is `{source_revision}`. The raw directory retains executable/source/protocol hashes, source dirty state, exact argv, seeds, timestamps, stdout/stderr and statuses. `summary.json` binds the summarized JSON, discovery artifacts, environment and statuses by SHA-256. This renderer derives descriptive observations only; it does not verify that every recorded claim or receipt is truthful. Deterministic replay is not independent scientific review.\n\nNo key recovery, hardware performance comparison, solver certificate, universal key-independence, security proof or accepted scientific-goal status follows. Source for the component experiment: https://eprint.iacr.org/2026/1549 . Its complete attack is not reproduced.\n",rows.len(),input.display());
    write_new(&output.join("REPORT.md"), &md)?;
    let typ = format!(
        r##"#set page(paper: "a4", margin: 18mm)
#set text(font: "DejaVu Sans", size: 10pt)
= AES differential discovery laboratory
Descriptive rendering of recorded synthetic fixtures

{successful}/{} recorded commands exited successfully. Recorded source revision:
{source_revision}. The raw directory retains
exact commands, seeds, executable/source/protocol hashes, logs and statuses.

== Differential endpoint observations
#text({observations_typst})

Per-family frozen event counts, keys, sample counts, input differences and alpha
are in summary.json. A failed matched control withholds that family's candidates.

#image("differential-events.svg", width: 100%)

This does not establish security, absence of a different distinguisher, novelty,
key recovery, or universal key independence. The native software's conservative
bounds and all per-key counts remain in the source JSON.

#pagebreak()
== Integral and component observations
#image("integral-structures.svg", width: 100%)
#image("quartet-component.svg", width: 100%)

#text({component_typst})

Source for the component experiment:
#link("https://eprint.iacr.org/2026/1549").

== Validation and limits
The native summarizer derives only descriptive counts from recorded results.
It does not independently verify the truth of every receipt. Summary JSON retains
source-artifact hashes. Deterministic replay is not independent scientific review.
No throughput, attack-cost gain or scientific campaign promotion follows.
"##,
        rows.len(),
        observations_typst = serde_json::to_string(&observations)?,
        component_typst = serde_json::to_string(&component_observations)?
    );
    write_new(&output.join("report.typ"), &typ)?;
    println!("Rendered recorded observations into {}", output.display());
    Ok(())
}
fn main() {
    if let Err(e) = run() {
        eprintln!("summarize: {e}");
        std::process::exit(2);
    }
}

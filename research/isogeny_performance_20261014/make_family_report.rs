//! Native coverage report; never derives performance ratios from the legacy harness.
use isogeny_algos::{json::Json, sha256::sha256_hex};
use std::{collections::BTreeMap, fmt::Write as _, path::Path};
fn s<'a>(j: &'a Json, k: &str) -> &'a str {
    j.get(k).unwrap().as_str().unwrap()
}
fn main() {
    let a: Vec<_> = std::env::args().collect();
    assert_eq!(a.len(), 3);
    let root = Path::new(&a[1]);
    let out = Path::new(&a[2]);
    std::fs::create_dir(out).unwrap();
    let labels = ["ubuntu-24.04", "ubuntu-24.04-arm", "apple-m4-pro"];
    let mut cells: BTreeMap<String, Vec<(u64, u64)>> = BTreeMap::new();
    let mut hashes = String::new();
    for (i, label) in labels.iter().enumerate() {
        let raw = std::fs::read(root.join(format!("final-{label}.json"))).unwrap();
        let j = Json::parse(std::str::from_utf8(&raw).unwrap()).unwrap();
        assert_eq!(s(&j, "status"), "EXPLORATORY_ONLY");
        assert_eq!(
            j.get("qualified_optimization_ratio"),
            Some(&Json::Bool(false))
        );
        let Json::Arr(groups) = j.get("groups").unwrap() else {
            panic!("missing groups")
        };
        let Json::Arr(failed) = j.get("verified_false_rows").unwrap() else {
            panic!("missing failed rows")
        };
        assert_eq!(s(&j, "rows"), "395");
        assert_eq!(failed.len(), 1);
        assert_eq!(
            s(failed[0].get("record").unwrap(), "algo"),
            "ghs_with_kohel_volcano"
        );
        for g in groups {
            let entry = cells.entry(s(g, "group").into()).or_insert(vec![(0, 0); 3]);
            entry[i] = (
                s(g, "rows").parse().unwrap(),
                s(g, "verified_true").parse().unwrap(),
            );
        }
        writeln!(hashes, "{label}: {}", sha256_hex(&raw)).unwrap();
    }
    assert_eq!(cells.len(), 24);
    let mut table=String::from("| Reported group | AMD true / rows | ARM true / rows | M4 true / rows |\n|---|---:|---:|---:|\n");
    let mut typ_rows = String::new();
    let mut csv = String::from("group,amd_rows,amd_true,arm_rows,arm_true,m4_rows,m4_true\n");
    for (group, values) in &cells {
        writeln!(
            table,
            "| {group} | {} / {} | {} / {} | {} / {} |",
            values[0].1, values[0].0, values[1].1, values[1].0, values[2].1, values[2].0
        )
        .unwrap();
        writeln!(
            typ_rows,
            "[{group}], [{} / {}], [{} / {}], [{} / {}],",
            values[0].1, values[0].0, values[1].1, values[1].0, values[2].1, values[2].0
        )
        .unwrap();
        writeln!(
            csv,
            "{group},{},{},{},{},{},{}",
            values[0].0, values[0].1, values[1].0, values[1].1, values[2].0, values[2].1
        )
        .unwrap();
    }
    let svg = r##"<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="190"><rect width="1000" height="190" fill="white"/><g font-family="Arial,sans-serif" fill="#17324d"><text x="20" y="28" font-size="21">Legacy coverage census: each profile retains 395 records</text><text x="20" y="62" font-size="17">AMD EPYC 7763</text><text x="330" y="62" font-size="17">ARM64 0x41 / 0xd49</text><text x="675" y="62" font-size="17">Apple M4 Pro</text><g fill="#eaf3f8" stroke="#6d91a6"><rect x="20" y="78" width="285" height="62"/><rect x="330" y="78" width="320" height="62"/><rect x="675" y="78" width="305" height="62"/></g><text x="35" y="104" font-size="17">394 true flags · 1 false flag</text><text x="345" y="104" font-size="17">394 true flags · 1 false flag</text><text x="690" y="104" font-size="17">394 true flags · 1 false flag</text><text x="35" y="128" font-size="14">Native outer isolation admitted</text><text x="345" y="128" font-size="14">Native outer isolation admitted</text><text x="690" y="128" font-size="14">Shared host; unqualified wall data</text><text x="20" y="173" font-size="15">No exact output fingerprints or matched A/A controls. No qualified optimization ratio on any profile.</text></g></svg>"##;
    let body="The 24 explicit command groups emit these 24 reported labels; aliases and split subgroups prevent a one-to-one command mapping. All rows, raw timings and verification flags remain preserved. The one false flag on every profile is ghs_with_kohel_volcano on a 17-bit fixture; its note describes split primes that may not generate the needed class-group orbit. A false flag is retained, not converted into a successful path.\n\nAMD and ARM use original source 79006eb85a40 from hosted run 38766895864. Their outer native isolation receipts pass, but those artifacts did not retain the legacy bench executable. The new M4 run uses fourth-candidate source 9a55ce6d22b9473c1d2afbf398bfdf51ae22b544 and retains executable, source, compiler and shell/resource receipts. These different revisions are not pooled or compared for speed.\n\nEach record is a legacy stage measurement, not an exact-output attestation. Neither complete input/output fingerprints nor matched A/A and interleaved A/B controls exist for this harness. The shared M4 supplies no qualified wall timing. Rows marked verified refer only to their existing harness-specific check, not independent replay. Named P192/P224 whole-search cost and the remaining exact family adapters stay open.";
    let typ=format!("#set page(paper: \"a4\", margin: 15mm)\n#set text(font: \"Arial\", size: 8pt)\n= Exploratory native family coverage\n2026-10-15 · Three native profiles · 395 records each\n\n#image(\"COVERAGE.svg\", width: 100%)\n\n#table(columns: (1.5fr,1fr,1fr,1fr), inset: 2pt, [Reported group], [AMD true / rows], [ARM true / rows], [M4 true / rows],\n{typ_rows})\n\n{body}\n\nCensus SHA256 receipts:\n#text(size: 5.5pt)[{}]\n",hashes.replace('\n',"\\\n"));
    for (file,data)in[("REPORT.md",format!("# Exploratory native family coverage\n\n![Coverage](COVERAGE.svg)\n\n{table}\n\n{body}\n\nCensus SHA256:\n\n```text\n{hashes}```\n")),("REPORT.typ",typ),("COVERAGE.svg",svg.into()),("COVERAGE.csv",csv),("coverage.typ","#set page(width: 1000pt, height: 190pt, margin: 0pt)\n#image(\"COVERAGE.svg\", width: 100%)\n".into())]{std::fs::write(out.join(file),data).unwrap();}
}

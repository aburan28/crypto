//! Re-derive the selected pair's full-fiber audit summary and diagram from receipts.

use std::fs;
use std::path::Path;

use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use serde_json::{json, Value};

fn number(row: &Value, key: &str) -> u64 {
    row[key].as_u64().unwrap_or_else(|| panic!("missing {key}"))
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    assert_eq!(args.len(), 5, "usage: fiber checker summary.json graph.svg");
    assert!(!Path::new(&args[3]).exists(), "refuse summary overwrite");
    assert!(!Path::new(&args[4]).exists(), "refuse graph overwrite");
    let fiber_bytes = fs::read(&args[1]).expect("fiber receipt");
    let checker_bytes = fs::read(&args[2]).expect("checker receipt");
    let fiber: Value = serde_json::from_slice(&fiber_bytes).expect("fiber JSON");
    let checker: Value = serde_json::from_slice(&checker_bytes).expect("checker JSON");
    assert_eq!(fiber["schema"], "different_j_order_fiber/v1");
    assert_eq!(checker["schema"], "different_j_order_fiber_check/v1");
    assert_eq!(fiber["all_checks_passed"], true);
    assert_eq!(checker["status"], "pass");
    for key in [
        "source_points",
        "image_points",
        "subgroup_points",
        "subgroup_pair_checks",
        "full_map_sha256",
        "subgroup_map_sha256",
    ] {
        assert_eq!(fiber[key], checker[key], "checker disagreement: {key}");
    }
    let source = number(&fiber, "source_points");
    let image = number(&fiber, "image_points");
    let subgroup = number(&fiber, "subgroup_points");
    let pair_checks = number(&fiber, "subgroup_pair_checks");
    assert_eq!(
        (source, image, subgroup, pair_checks),
        (44655, 14885, 229, 52441)
    );
    assert_eq!(source, 3 * image);
    assert_eq!(checker["three_point_fibers"], image);
    assert_eq!(fiber["fiber_histogram"]["3"], image);
    assert_eq!(pair_checks, subgroup * subgroup);
    assert_eq!(fiber["identity_fiber"].as_array().unwrap().len(), 3);
    let summary = json!({
        "schema":"different_j_order_fiber_summary/v1",
        "fiber_receipt_sha256":sha256_hex(&fiber_bytes),
        "checker_receipt_sha256":sha256_hex(&checker_bytes),
        "p":number(&fiber,"p"),"r":number(&fiber,"r"),
        "source_points":source,"image_points":image,"fiber_size":3,
        "subgroup_points":subgroup,"subgroup_pair_checks":pair_checks,
        "full_map_sha256":fiber["full_map_sha256"],
        "subgroup_map_sha256":fiber["subgroup_map_sha256"],
        "status":"replayed_full_fiber_audit",
        "dlp_recovery":Value::Null,"portable_counted_S":Value::Null
    });
    fs::write(
        &args[3],
        format!("{}\n", serde_json::to_string_pretty(&summary).unwrap()),
    )
    .unwrap();

    let image_width = 660.0 * (image as f64 / source as f64);
    let svg = format!(
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="440" viewBox="0 0 1000 440">
<rect width="1000" height="440" fill="#f8fafc"/>
<text x="48" y="37" font-family="Arial,sans-serif" font-size="22" font-weight="700" fill="#0f172a">p=211 full rational isogeny-fiber audit</text>
<text x="48" y="61" font-family="Arial,sans-serif" font-size="13" fill="#475569">Direct enumeration and independent replay · points over F₂₁₁² · run FIBER_V1</text>
<text x="48" y="113" font-family="Arial,sans-serif" font-size="16" fill="#0f172a">Source curve, enumerated</text>
<rect x="294" y="87" width="660" height="38" rx="4" fill="#0e7490"/>
<text x="316" y="113" font-family="Arial,sans-serif" font-size="17" font-weight="700" fill="white">{source} points</text>
<text x="48" y="176" font-family="Arial,sans-serif" font-size="16" fill="#0f172a">Distinct image points</text>
<rect x="294" y="150" width="{image_width:.1}" height="38" rx="4" fill="#ea580c"/>
<text x="316" y="176" font-family="Arial,sans-serif" font-size="17" font-weight="700" fill="white">{image} points</text>
<text x="48" y="228" font-family="Arial,sans-serif" font-size="14" fill="#334155">Every image point has exactly 3 rational preimages: {image} × 3 = {source}.</text>
<line x1="48" y1="249" x2="954" y2="249" stroke="#cbd5e1"/>
<rect x="48" y="273" width="392" height="123" rx="9" fill="#e0f2fe" stroke="#0e7490"/>
<text x="68" y="303" font-family="Arial,sans-serif" font-size="16" font-weight="700" fill="#0c4a6e">Kernel fiber</text>
<text x="68" y="334" font-family="Arial,sans-serif" font-size="15" fill="#0f172a">φ⁻¹(O) = &#123;O, (3,1), (3,−1)&#125;</text>
<text x="68" y="363" font-family="Arial,sans-serif" font-size="13" fill="#334155">Full map digest checked independently.</text>
<rect x="460" y="273" width="494" height="123" rx="9" fill="#ffedd5" stroke="#ea580c"/>
<text x="480" y="303" font-family="Arial,sans-serif" font-size="16" font-weight="700" fill="#9a3412">Order-{subgroup} subgroup transport</text>
<text x="480" y="334" font-family="Arial,sans-serif" font-size="15" fill="#0f172a">{subgroup} source and {subgroup} image points</text>
<text x="480" y="363" font-family="Arial,sans-serif" font-size="13" fill="#334155">{pair_checks} / {pair_checks} ordered-pair homomorphism checks</text>
<text x="48" y="426" font-family="Arial,sans-serif" font-size="12" fill="#64748b">Map/fiber correctness audit; relation yield, DLP cost and speed are outside this run.</text>
</svg>
"##
    );
    fs::write(&args[4], svg).unwrap();
    println!(
        "{}",
        json!({"status":"ok","source":source,"image":image,"pairs":pair_checks})
    );
}

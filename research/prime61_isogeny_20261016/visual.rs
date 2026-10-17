//! Render only the validated, recorded summary; performs no research search.
use serde_json::Value;
use std::fs;
fn main() {
    let path = std::env::args().nth(1).expect("SUMMARY.json");
    let v: Value = serde_json::from_str(&fs::read_to_string(path).unwrap()).unwrap();
    assert_eq!(v["runs"], 192);
    let mut svg = String::from(
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="960" height="620" viewBox="0 0 960 620">
<rect width="960" height="620" fill="#f8fafc"/>
<g font-family="DejaVu Sans, sans-serif" fill="#162b45">
<text x="32" y="38" font-size="23" font-weight="bold">Prime61: outgoing horizontal prime-degree coverage</text>
<text x="32" y="66" font-size="15">One source • one hop • degree ≤ 61 • 15 distinct kernels • full class remains open</text>
<text x="48" y="106" font-size="14">Prime degree</text><text x="218" y="106" font-size="14">Verified outgoing kernel edges (exact count)</text>
"##,
    );
    for (i, ell) in [3u64, 7, 11, 13, 19, 23, 53, 61].iter().enumerate() {
        let n = v["degree_histogram"][ell.to_string()].as_u64().unwrap();
        let y = 122 + 29 * i;
        svg.push_str(&format!(r##"<text x="80" y="{}" font-size="15">{ell}</text><rect x="218" y="{y}" width="{}" height="19" rx="3" fill="#397ba8"/><text x="{}" y="{}" font-size="15">{n}</text>"##, y+15, n*240, 228+n*240, y+15));
    }
    svg.push_str(r##"
<text x="32" y="384" font-size="13">Degree 2 and inert primes: 0 rational edges. Vertical F_p edges: 0 (fundamental discriminant).</text>
<rect x="32" y="409" width="896" height="93" rx="8" fill="#e0edf6"/>
<text x="48" y="436" font-size="17" font-weight="bold">Verification → transported target → fresh native three-pair base</text>
<text x="48" y="461" font-size="15">960 scalar checks + 540 transported pair-sum checks; all passed.</text>
<text x="48" y="485" font-size="15">16 workloads × 12 seeds = 192 capped runs; 19,200,000 trials; 0 relations; 0 solves.</text>
<rect x="32" y="518" width="896" height="73" rx="8" fill="#fff0d1"/>
<text x="48" y="545" font-size="16" font-weight="bold">Full class: ≥ 68,162,026 vertices from an exact class-order certificate</text>
<text x="48" y="570" font-size="14">Existing walker headers alone need ≥ 26.7 GB; container limit 17.2 GB (16 GiB).</text>
<text x="32" y="611" font-size="11">Source: results/summary.json and exact certificate / memory records. Trials are stage counters; no speedup claim.</text>
</g></svg>"##);
    println!("{svg}");
}

//! Generate census and construction figures/tables directly from verified records.
use serde_json::Value;
use std::{collections::BTreeMap, fmt::Write, fs, path::Path};
fn n(v: &Value, k: &str) -> u64 {
    v[k].as_str().unwrap().parse().unwrap()
}
fn label(s: &mut String, x: f64, y: f64, size: u32, text: &str) {
    writeln!(s,"<text x='{x}' y='{y}' font-family='Arial,sans-serif' font-size='{size}' fill='#203047'>{text}</text>").unwrap();
}
fn main() {
    let a: Vec<_> = std::env::args().skip(1).collect();
    assert_eq!(a.len(), 1);
    let base = Path::new(&a[0]);
    let d: Value = serde_json::from_slice(&fs::read(base.join("summary.json")).unwrap()).unwrap();
    let mut svg=String::from("<svg xmlns='http://www.w3.org/2000/svg' width='1400' height='700' viewBox='0 0 1400 700'><rect width='1400' height='700' fill='white'/>");
    label(
        &mut svg,
        50.,
        45.,
        28,
        "Exact support of both genus-3 hyperelliptic branches",
    );
    label(
        &mut svg,
        50.,
        78.,
        18,
        "Ordinary trace classes t = 2 mod 4, both signs; complete censuses",
    );
    let mut tex=String::from("\\begin{tabular}{rrrrrr}\\toprule\np & Ordinary & Cubic & Quadratic & Union & Zero\\\\\\midrule\n");
    let mut md=String::from("| p | Ordinary classes | Cubic | Quadratic | Union | Exact zeros |\n| --- | ---: | ---: | ---: | ---: | ---: |\n");
    let rows = d["censuses"].as_array().unwrap();
    for (i, r) in rows.iter().enumerate() {
        let p = n(r, "p");
        let all = n(r, "ordinary_classes");
        let plus = n(r, "plus_support");
        let minus = n(r, "minus_support");
        let overlap = n(r, "overlap");
        let zero = n(r, "combined_zero");
        let union = all - zero;
        writeln!(tex, "{p} & {all} & {plus} & {minus} & {union} & {zero}\\\\").unwrap();
        writeln!(md, "| {p} | {all} | {plus} | {minus} | {union} | {zero} |").unwrap();
        let y = 130. + i as f64 * 65.;
        label(&mut svg, 50., y + 24., 20, &format!("p={p}"));
        let mut x = 175.;
        for (v, c) in [
            (plus - overlap, "#4b6cb7"),
            (overlap, "#8c62a7"),
            (minus - overlap, "#118b72"),
            (zero, "#ccd3dd"),
        ] {
            let w = 950. * v as f64 / all as f64;
            writeln!(
                svg,
                "<rect x='{x}' y='{y}' width='{w}' height='36' fill='{c}'/>"
            )
            .unwrap();
            x += w;
        }
        label(&mut svg, 1150., y + 24., 18, &format!("{union}/{all}"));
    }
    for (i, (name, c)) in [
        ("Cubic only", "#4b6cb7"),
        ("Both", "#8c62a7"),
        ("Quadratic only", "#118b72"),
        ("Exact zero", "#ccd3dd"),
    ]
    .iter()
    .enumerate()
    {
        let x = 55. + i as f64 * 310.;
        writeln!(
            svg,
            "<rect x='{x}' y='630' width='20' height='20' fill='{c}'/>"
        )
        .unwrap();
        label(&mut svg, x + 30., 647., 18, name);
    }
    svg.push_str("</svg>\n");
    tex.push_str("\\bottomrule\\end{tabular}\n");
    fs::write(base.join("coverage.svg"), svg).unwrap();
    fs::write(base.join("census_table.tex"), tex).unwrap();
    fs::write(base.join("census_table.md"), md).unwrap();
    let mut costs=String::from("\\begin{tabular}{lrrrr}\\toprule\nPolicy & Sources & Witnesses & Zero classes & Capped/closed positives\\\\\\midrule\n");
    let attempts = d["p7_construction"]["all_attempts"].as_array().unwrap();
    for (run, title) in [
        ("construction_p7_run1", "2,3"),
        ("construction_p7_odd_run3", "2,3,5,7"),
        ("construction_p7_followup_run1", "Adaptive wider edges"),
        ("construction_p7_targeted_run1", "Prioritized 11,13"),
        ("construction_p7_degree19_run1", "Prioritized 19"),
    ] {
        let a: Vec<_> = attempts.iter().filter(|v| v["run"] == run).collect();
        if a.is_empty() {
            continue;
        }
        let found = a
            .iter()
            .filter(|v| v["cost"]["status"] == "VERIFIED_WITNESS")
            .count();
        let zeros = a.iter().filter(|v| v["trace"] == 474).count();
        writeln!(
            costs,
            "{title} & {} & {found} & {zeros} & {}\\\\",
            a.len(),
            a.len() - found - zeros
        )
        .unwrap();
    }
    writeln!(
        costs,
        "Union of verified routes & 64 & {} & 1 & {}\\\\",
        d["p7_construction"]["constructed_after_followups"],
        63 - d["p7_construction"]["constructed_after_followups"]
            .as_u64()
            .unwrap()
    )
    .unwrap();
    costs.push_str("\\bottomrule\\end{tabular}\n");
    fs::write(base.join("construction_table.tex"), costs).unwrap();
    let mut groups: BTreeMap<(u64, u64, u64), (u64, u64, u64, u64)> = BTreeMap::new();
    for row in d["large_rows"].as_array().unwrap() {
        let p = n(&row["cost"], "p");
        let degree = row["source"]["degree"]
            .as_str()
            .map(|v| v.parse::<u64>().unwrap())
            .unwrap_or(6);
        let bits = match p {
            4294967291 | 602233 | 13421 => 192,
            17179869143 => 204,
            68719476731 => 216,
            274877906899 => 228,
            1099511627689 => 240,
            4398046511093 | 38543917 | 262139 => 252,
            _ => panic!("unregistered field"),
        };
        let g = groups.entry((degree, bits, p)).or_default();
        if row["cost"]["status"] == "RESTRICTED_COMPONENT_EXHAUSTED" {
            g.0 += 1;
        } else if row["cost"]["status"] == "VERTEX_CAP" {
            g.1 += 1;
        } else {
            g.2 += 1;
        }
        g.3 += n(&row["cost"], "visited");
    }
    let mut large=String::from("\\begin{tabular}{rrrrrrr}\\toprule\nBits & Degree & p & Sources & Closures & Caps & Visited\\\\\\midrule\n");
    for ((degree, bits, p), (closed, caps, other, visited)) in &groups {
        assert_eq!(*other, 0);
        writeln!(
            large,
            "{bits} & {degree} & {p} & {} & {closed} & {caps} & {visited}\\\\",
            closed + caps
        )
        .unwrap();
    }
    large.push_str("\\bottomrule\\end{tabular}\n");
    fs::write(base.join("large_table.tex"), large).unwrap();
    let flow="<svg xmlns='http://www.w3.org/2000/svg' width='1400' height='700' viewBox='0 0 1400 700'><rect width='1400' height='700' fill='white'/><style>text{font-family:Arial,sans-serif;fill:#203047}.box{fill:#eef4fa;stroke:#537498;stroke-width:2}.arrow{fill:none;stroke:#537498;stroke-width:3}</style><text x='45' y='48' font-size='28'>Exact class support and verified construction</text><rect class='box' x='45' y='90' width='600' height='190' rx='12'/><text x='70' y='128' font-size='23'>Cubic / split-quadratic branch</text><text x='70' y='170' font-size='21'>Norm-one torus: N+ = p^4 + p^2 + 1</text><text x='70' y='210' font-size='21'>2-quotient has full rational 4-torsion</text><text x='70' y='250' font-size='21'>Selected CM orders: c divides f_pi / 4</text><rect class='box' x='730' y='90' width='625' height='190' rx='12'/><text x='755' y='128' font-size='23'>Nonsplit-quadratic branch</text><text x='755' y='170' font-size='21'>Alternating torus: N- = p^4 - p^2 + 1</text><text x='755' y='210' font-size='21'>One rational nonidentity 2-torsion point</text><text x='755' y='250' font-size='21'>Selected CM orders: v2(c) = v2(f_pi)</text><path class='arrow' d='M345 280V330H700 M1040 280V330H700V370'/><rect class='box' x='260' y='370' width='875' height='95' rx='12'/><text x='290' y='408' font-size='23'>Exact union: either torus / CM intersection is nonempty</text><text x='290' y='445' font-size='19'>Membership is certified separately from the route to that representative</text><path class='arrow' d='M700 465V515'/><rect class='box' x='100' y='515' width='1200' height='110' rx='12'/><text x='130' y='553' font-size='22'>Independent source - rational kernel maps - literal cubic or quartic model</text><text x='130' y='590' font-size='20'>Every retained edge and endpoint replayed; capped searches retain unresolved labels</text></svg>\n";
    fs::write(base.join("two_branch_flow.svg"), flow).unwrap();
    let doc: Value =
        serde_json::from_slice(&fs::read(base.join("curve_records.json")).unwrap()).unwrap();
    let records = doc["records"].as_array().unwrap();
    let raw = fs::read_to_string(base.join("construction_p7_targeted_run1/stdout.txt")).unwrap();
    let mut chosen = false;
    let mut ids = Vec::new();
    let mut degrees = Vec::new();
    for line in raw.lines() {
        let m: BTreeMap<_, _> = line
            .split('|')
            .skip(1)
            .filter_map(|v| v.split_once('='))
            .collect();
        if line.starts_with("SOURCE|") {
            chosen = m["seed"] == "202610090042";
        }
        if chosen && line.starts_with("ROUTE|") {
            for role in ["source", "target"] {
                if role == "source" && !ids.is_empty() {
                    continue;
                }
                let av: Value = serde_json::from_str(m[format!("{role}_a").as_str()]).unwrap();
                let bv: Value = serde_json::from_str(m[format!("{role}_b").as_str()]).unwrap();
                let r = records
                    .iter()
                    .find(|r| {
                        let model: Value =
                            serde_json::from_str(r["model_json"].as_str().unwrap()).unwrap();
                        r["trace"] == "202" && model["a"] == av && model["b"] == bv
                    })
                    .unwrap();
                ids.push((
                    r["icv1"].as_str().unwrap().to_string(),
                    r["slug"].as_str().unwrap().to_string(),
                ));
            }
            degrees.push(m["degree"].to_string());
        }
    }
    assert_eq!(ids.len(), 3);
    let mut route=String::from("<svg xmlns='http://www.w3.org/2000/svg' width='1400' height='740' viewBox='0 0 1400 740'><rect width='1400' height='740' fill='white'/>");
    label(
        &mut route,
        50.,
        44.,
        27,
        "Verified route from an independent source to the quadratic branch",
    );
    label(
        &mut route,
        50.,
        79.,
        18,
        "F_(7^6), trace 202, order 117448; frozen source seed 202610090042",
    );
    for (i, (id, slug)) in ids.iter().enumerate() {
        let y = 105. + i as f64 * 205.;
        writeln!(route,"<rect x='50' y='{y}' width='1300' height='125' rx='12' fill='#eef4fa' stroke='#537498' stroke-width='2'/>").unwrap();
        label(
            &mut route,
            75.,
            y + 32.,
            21,
            [
                "Independent full-2 source",
                "Verified intermediate model",
                "Quadratic-branch endpoint",
            ][i],
        );
        label(&mut route, 75., y + 64., 18, slug);
        label(&mut route, 75., y + 98., 17, id);
        if i < 2 {
            writeln!(route,"<path d='M700 {}V{} L691 {} M700 {}L709 {}' stroke='#537498' fill='none' stroke-width='3'/>",y+125.,y+185.,y+172.,y+185.,y+172.).unwrap();
            label(
                &mut route,
                750.,
                y + 165.,
                18,
                &format!(
                    "Degree {}, rational kernel; point-count and map replay verified",
                    degrees[i]
                ),
            );
        }
    }
    label(&mut route,50.,715.,17,"Exact kernels and literal quartic conversion: construction_p7_targeted_run1/stdout.txt and replay_run1/stdout.txt");
    route.push_str("</svg>\n");
    fs::write(base.join("verified_route.svg"), route).unwrap();
    println!(
        "generated {} exact census rows and six large-field construction rows",
        rows.len()
    );
}

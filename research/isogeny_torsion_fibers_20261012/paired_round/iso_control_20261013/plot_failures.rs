//! Visualize verified, failed, and not-reached rows in bounded v1/v2 runs.
use std::collections::BTreeMap;
use std::env;
use std::fs;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    assert_eq!(
        args.len(),
        5,
        "plot_failures RAW.tsv VERSION COUNTS.tsv GRAPH.svg"
    );
    let raw = fs::read_to_string(&args[1])?;
    let mut lines = raw.lines();
    let header: Vec<&str> = lines.next().ok_or("empty raw TSV")?.split('\t').collect();
    let col = |name: &str| {
        header
            .iter()
            .position(|&v| v == name)
            .expect("required column")
    };
    let mut counts: BTreeMap<(u64, String), (u64, u64)> = BTreeMap::new();
    let mut target_status: BTreeMap<(u64, String, u64), bool> = BTreeMap::new();
    for (i, line) in lines.enumerate() {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), header.len(), "row {} width", i + 2);
        let p = fields[col("p")].parse::<u64>()?;
        let arm = fields[col("arm")].to_owned();
        let target_index = fields[col("target_index")].parse::<u64>()?;
        assert!(target_index < 8);
        let verified = fields[col("status")] == "verified";
        assert!(target_status
            .insert((p, arm.clone(), target_index), verified)
            .is_none());
        let row = counts.entry((p, arm)).or_default();
        if verified {
            row.0 += 1;
        } else {
            row.1 += 1;
        }
    }
    let arms = [
        "source_relations",
        "transport_relations",
        "isomorphism_relations",
        "rho_negation",
        "rho_six_class",
    ];
    let mut table = String::from("version\tp\tarm\tverified\tfailed\tnot_reached\n");
    let mut svg = format!("<svg xmlns='http://www.w3.org/2000/svg' width='980' height='650' viewBox='0 0 980 650'>\n<rect width='980' height='650' fill='#fff'/><style>text{{font-family:Arial,sans-serif;fill:#172033}}.title{{font-size:20px;font-weight:bold}}.label{{font-size:14px}}.small{{font-size:12px}}</style>\n<text x='30' y='35' class='title'>{} bounded run: arm completion on frozen targets</text>\n<text x='30' y='56' class='small'>Each row represents eight intended source targets. Gray means execution did not reach the target.</text>\n", args[2]);
    for (side, p) in [13, 199].iter().enumerate() {
        let top = 105.0 + side as f64 * 265.0;
        svg.push_str(&format!(
            "<text x='30' y='{:.1}' class='label'>F_{p}, subgroup order {}</text>\n",
            top - 20.0,
            if *p == 13 { 19 } else { 211 }
        ));
        for (j, arm) in arms.iter().enumerate() {
            let (ok, failed) = counts
                .get(&(*p, (*arm).into()))
                .copied()
                .unwrap_or_default();
            assert!(ok + failed <= 8);
            let absent = 8 - ok - failed;
            table.push_str(&format!(
                "{}\t{p}\t{arm}\t{ok}\t{failed}\t{absent}\n",
                args[2]
            ));
            let y = top + j as f64 * 43.0;
            svg.push_str(&format!(
                "<text x='35' y='{:.1}' class='small'>{}</text>\n",
                y + 17.0,
                arm.replace('_', " ")
            ));
            let x0 = 260.0;
            let width = 73.0;
            for k in 0..8 {
                let color = match target_status.get(&(*p, (*arm).into(), k)) {
                    Some(true) => "#0f766e",
                    Some(false) => "#d97706",
                    None => "#cbd5e1",
                };
                let x = x0 + k as f64 * 82.0;
                svg.push_str(&format!("<rect x='{x:.1}' y='{y:.1}' width='{width:.1}' height='25' rx='3' fill='{color}'/><text x='{:.1}' y='{:.1}' class='small'>{}</text>\n", x+30.0, y+17.0, k));
            }
            svg.push_str(&format!(
                "<text x='925' y='{:.1}' class='small'>{ok}/8</text>\n",
                y + 18.0
            ));
        }
    }
    svg.push_str("<rect x='30' y='595' width='15' height='15' fill='#0f766e'/><text x='52' y='608' class='small'>verified</text><rect x='155' y='595' width='15' height='15' fill='#d97706'/><text x='177' y='608' class='small'>failed</text><rect x='260' y='595' width='15' height='15' fill='#cbd5e1'/><text x='282' y='608' class='small'>not reached</text></svg>\n");
    fs::write(&args[3], table)?;
    fs::write(&args[4], svg)?;
    Ok(())
}

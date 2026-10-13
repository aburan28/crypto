//! Derive fixed-run summaries and an SVG from the verified raw TSV.
use std::collections::BTreeMap;
use std::env;
use std::fs;
use std::path::Path;

#[derive(Default)]
struct Arm {
    ops: Vec<u64>,
    setup: u64,
    collection: u64,
    descent: u64,
    verification: u64,
    rho: u64,
    relation_trials: u64,
    descent_trials: u64,
    row_mul_subs: u64,
    fb_lifts: u64,
}
fn col(header: &[&str], name: &str) -> usize {
    header.iter().position(|&v| v == name).expect("column")
}
fn num(row: &[&str], header: &[&str], name: &str) -> u64 {
    row[col(header, name)].parse().expect("integer")
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    assert_eq!(args.len(), 5, "summarize RAW SUMMARY PHASES SVG");
    let input = fs::read_to_string(&args[1])?;
    let mut lines = input.lines();
    let header: Vec<_> = lines.next().ok_or("empty raw file")?.split('\t').collect();
    let mut groups = BTreeMap::<String, Arm>::new();
    let mut targets = BTreeMap::<(usize, String), u64>::new();
    for line in lines {
        let row: Vec<_> = line.split('\t').collect();
        assert_eq!(row.len(), header.len());
        assert_eq!(row[col(&header, "status")], "verified");
        let name = row[col(&header, "arm")].to_string();
        let index = num(&row, &header, "target_index") as usize;
        let value = num(&row, &header, "total_gae");
        assert!(targets.insert((index, name.clone()), value).is_none());
        let arm = groups.entry(name).or_default();
        arm.ops.push(num(&row, &header, "total_gae"));
        arm.setup += num(&row, &header, "setup_gae");
        arm.collection += num(&row, &header, "collection_gae");
        arm.descent += num(&row, &header, "descent_gae");
        arm.verification += num(&row, &header, "verification_gae");
        arm.rho += num(&row, &header, "rho_gae");
        arm.relation_trials += num(&row, &header, "relation_trials");
        arm.descent_trials += num(&row, &header, "descent_trials");
        arm.row_mul_subs += num(&row, &header, "row_mul_subs");
        arm.fb_lifts += num(&row, &header, "fb_lifts");
    }
    assert_eq!(groups.len(), 3);
    let order = ["rho_negation", "source_relations", "transport_relations"];
    let rho_total: u64 = groups["rho_negation"].ops.iter().sum();
    let root = 103f64.sqrt();
    let floor = (std::f64::consts::PI / 4.0).sqrt();
    let mut summary = String::from("arm\tn\ttotal_gae\tmean_gae\tmin_gae\tmax_gae\tmean_s\tmin_s\tmax_s\ttotal_over_rho\tmean_s_over_A2_floor\trelation_trials\tdescent_trials\n");
    let mut phases = String::from("arm\tn\tsetup_gae\tcollection_gae\tdescent_gae\tverification_gae\trho_gae\trow_mul_subs\tfb_lifts\n");
    let version = if args[1].contains("v1") { "V1" } else { "V2" };
    let mut svg = format!("<svg xmlns='http://www.w3.org/2000/svg' width='1100' height='530' viewBox='0 0 1100 530'>\n<rect width='1100' height='530' fill='#fff'/><style>text{{font-family:Arial,sans-serif;fill:#172033}}.title{{font-size:22px;font-weight:bold}}.label{{font-size:16px}}.small{{font-size:13px}}</style>\n<text x='35' y='38' class='title'>Non-isomorphic F617 pair: complete counted group-operation S</text>\n<text x='35' y='63' class='small'>8 identical source targets, subgroup order 103; mean bar and observed minimum–maximum whisker. {version}, 24/24 verified arms.</text>\n");
    let colors = ["#64748b", "#0284c7", "#be185d"];
    let labels = [
        "negation rho",
        "source relations",
        "degree-three destination relations",
    ];
    let max_s = groups
        .values()
        .flat_map(|arm| arm.ops.iter())
        .max()
        .copied()
        .unwrap() as f64
        / root
        * 1.15;
    let bottom = 430f64;
    let height = 300f64;
    for tick in 0..=4 {
        let y = bottom - height * tick as f64 / 4.0;
        svg.push_str(&format!("<line x1='80' y1='{y:.1}' x2='1050' y2='{y:.1}' stroke='#e2e8f0'/><text x='34' y='{:.1}' class='small'>{:.0}</text>\n", y+4.0, max_s*tick as f64/4.0));
    }
    for (j, name) in order.iter().enumerate() {
        let arm = &groups[*name];
        assert_eq!(arm.ops.len(), 8);
        let total: u64 = arm.ops.iter().sum();
        assert_eq!(
            total,
            arm.setup + arm.collection + arm.descent + arm.verification + arm.rho
        );
        let mean = total as f64 / 8.0;
        let min = *arm.ops.iter().min().unwrap();
        let max = *arm.ops.iter().max().unwrap();
        summary.push_str(&format!("{name}\t8\t{total}\t{mean:.6}\t{min}\t{max}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\n", mean/root, min as f64/root, max as f64/root, total as f64/rho_total as f64, mean/root/floor, arm.relation_trials, arm.descent_trials));
        phases.push_str(&format!(
            "{name}\t8\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            arm.setup,
            arm.collection,
            arm.descent,
            arm.verification,
            arm.rho,
            arm.row_mul_subs,
            arm.fb_lifts
        ));
        let x = 140.0 + j as f64 * 320.0;
        let y = bottom - (mean / root) / max_s * height;
        let ymin = bottom - (min as f64 / root) / max_s * height;
        let ymax = bottom - (max as f64 / root) / max_s * height;
        svg.push_str(&format!("<rect x='{x:.1}' y='{y:.1}' width='100' height='{:.1}' fill='{}'/><line x1='{:.1}' y1='{ymax:.1}' x2='{:.1}' y2='{ymin:.1}' stroke='#111827' stroke-width='3'/><line x1='{:.1}' y1='{ymax:.1}' x2='{:.1}' y2='{ymax:.1}' stroke='#111827' stroke-width='3'/><line x1='{:.1}' y1='{ymin:.1}' x2='{:.1}' y2='{ymin:.1}' stroke='#111827' stroke-width='3'/><text x='{:.1}' y='{:.1}' class='label'>{:.3} S</text><text x='{:.1}' y='460' class='label'>{}</text>\n", bottom-y, colors[j], x+50.0, x+50.0, x+28.0, x+72.0, x+28.0, x+72.0, x, ymax-14.0, mean/root, x-45.0, labels[j]));
    }
    svg.push_str("<text x='35' y='505' class='small'>S = (group additions + doublings)/sqrt(103), including cold factor-base projection, collection, recovery, descent, and scalar verification.</text></svg>\n");
    fs::write(&args[2], summary)?;
    fs::write(&args[3], phases)?;
    fs::write(&args[4], svg)?;
    let mut paired = String::from("target_index\trho_gae\tsource_gae\tdestination_gae\tdestination_minus_source\tdestination_over_source\n");
    for i in 0..8 {
        let rho = targets[&(i, "rho_negation".into())];
        let src = targets[&(i, "source_relations".into())];
        let dst = targets[&(i, "transport_relations".into())];
        paired.push_str(&format!(
            "{i}\t{rho}\t{src}\t{dst}\t{}\t{:.6}\n",
            dst as i64 - src as i64,
            dst as f64 / src as f64
        ));
    }
    let summary_name = Path::new(&args[2]).file_name().unwrap().to_string_lossy();
    fs::write(
        Path::new(&args[2]).with_file_name(format!("paired_{summary_name}")),
        paired,
    )?;
    Ok(())
}

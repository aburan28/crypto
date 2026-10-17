//! Derive summaries and SVGs from the frozen v3 TSV; no solver work is rerun.
use std::collections::BTreeMap;
use std::env;
use std::fs;
use std::path::Path;

#[derive(Clone, Copy)]
struct Stats {
    mean: f64,
    min: f64,
    max: f64,
}

fn stats(values: &[u64]) -> Stats {
    assert_eq!(values.len(), 8);
    Stats {
        mean: values.iter().sum::<u64>() as f64 / values.len() as f64,
        min: *values.iter().min().unwrap() as f64,
        max: *values.iter().max().unwrap() as f64,
    }
}

fn esc(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    assert_eq!(
        args.len(),
        5,
        "plot RAW.tsv SUMMARY.tsv RELATION.svg RHO.svg"
    );
    let raw = fs::read_to_string(&args[1])?;
    let mut lines = raw.lines();
    let header: Vec<&str> = lines.next().ok_or("empty TSV")?.split('\t').collect();
    let col = |name: &str| {
        header
            .iter()
            .position(|&v| v == name)
            .expect("required column")
    };
    let phase_names = [
        "setup_gae",
        "collection_gae",
        "descent_gae",
        "verification_gae",
        "rho_gae",
        "field_muls",
        "field_sqrs",
        "field_invs",
        "fb_lifts",
        "pair_entries",
        "relation_trials",
        "relation_lookups",
        "independent_rows",
        "row_mul_subs",
        "row_inversions",
        "descent_trials",
        "descent_lookups",
        "map_evaluations",
        "map_muls",
        "map_sqrs",
        "map_invs",
        "rho_steps",
        "rho_walks",
        "rho_dp",
        "rho_unpriced_canonicalisations",
        "rho_canon_field_muls",
        "rho_canon_point_negations",
        "rho_attempts",
        "rho_class_setup_field_muls",
    ];
    let mut groups: BTreeMap<(u64, String), (u64, Vec<u64>)> = BTreeMap::new();
    let mut phases: BTreeMap<(u64, String), Vec<Vec<u64>>> = BTreeMap::new();
    for (n, line) in lines.enumerate() {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), header.len(), "row {} width", n + 2);
        assert_eq!(fields[col("status")], "verified", "row {} status", n + 2);
        let p = fields[col("p")].parse::<u64>()?;
        let r = fields[col("r")].parse::<u64>()?;
        let arm = fields[col("arm")].to_owned();
        groups
            .entry((p, arm.clone()))
            .or_insert((r, Vec::new()))
            .1
            .push(fields[col("total_gae")].parse()?);
        phases.entry((p, arm)).or_default().push(
            phase_names
                .iter()
                .map(|&name| fields[col(name)].parse::<u64>().unwrap())
                .collect(),
        );
    }
    assert_eq!(groups.len(), 10);
    let arms = [
        "rho_negation",
        "rho_six_class",
        "source_relations",
        "transport_relations",
        "isomorphism_relations",
    ];
    let mut summary = String::from("p\tr\tarm\tn\tmean_gae\tmin_gae\tmax_gae\tmean_s\tmin_s\tmax_s\tmean_s_over_rho_negation\tmean_s_over_A2_floor\tmean_s_over_A6_floor\n");
    for p in [13, 199] {
        let rho = stats(&groups[&(p, "rho_negation".into())].1);
        for arm in arms {
            let (r, values) = &groups[&(p, arm.into())];
            let s = stats(values);
            let root = (*r as f64).sqrt();
            let mean_s = s.mean / root;
            summary.push_str(&format!("{p}\t{r}\t{arm}\t8\t{:.6}\t{:.0}\t{:.0}\t{mean_s:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\n", s.mean, s.min, s.max, s.min/root, s.max/root, s.mean/rho.mean, mean_s/(std::f64::consts::PI/4.0).sqrt(), mean_s/(std::f64::consts::PI/12.0).sqrt()));
        }
    }
    fs::write(&args[2], summary)?;
    let mut phase = format!(
        "p\tarm\tn\t{}\n",
        phase_names
            .iter()
            .map(|name| format!("mean_{name}"))
            .collect::<Vec<_>>()
            .join("\t")
    );
    for ((p, arm), rows) in &phases {
        phase.push_str(&format!("{p}\t{arm}\t{}", rows.len()));
        for j in 0..phase_names.len() {
            let mean = rows.iter().map(|row| row[j] as f64).sum::<f64>() / rows.len() as f64;
            phase.push_str(&format!("\t{mean:.6}"));
        }
        phase.push('\n');
    }
    let filename = Path::new(&args[2]).file_name().unwrap().to_string_lossy();
    fs::write(
        Path::new(&args[2]).with_file_name(format!("phase_{filename}")),
        phase,
    )?;

    let mut svg = String::from("<svg xmlns='http://www.w3.org/2000/svg' width='1120' height='495' viewBox='0 0 1120 495'>\n<rect width='1120' height='495' fill='#fff'/>\n<style>text{font-family:Arial,sans-serif;fill:#172033}.title{font-size:20px;font-weight:bold}.label{font-size:14px}.small{font-size:12px}</style>\n<text x='35' y='35' class='title'>Same-target relation and rho controls: counted group-operation S</text>\n<text x='35' y='57' class='small'>v3: eight public targets per order; whiskers show observed minimum and maximum. Map field work and native table work are separate.</text>\n");
    let colors = ["#64748b", "#0284c7", "#be185d", "#0f766e"];
    let labels = ["rho A2", "source", "degree-3", "scale iso"];
    let graph_arms = [
        "rho_negation",
        "source_relations",
        "transport_relations",
        "isomorphism_relations",
    ];
    for (side, p) in [13, 199].iter().enumerate() {
        let left = 50.0 + side as f64 * 550.0;
        let top = 110.0;
        let bottom = 390.0;
        let r = groups[&(*p, "rho_negation".into())].0;
        let vals: Vec<Stats> = graph_arms
            .iter()
            .map(|arm| stats(&groups[&(*p, (*arm).into())].1))
            .collect();
        let max_s = vals
            .iter()
            .map(|v| v.max / (r as f64).sqrt())
            .fold(0.0, f64::max)
            * 1.12;
        svg.push_str(&format!(
            "<text x='{left:.1}' y='89' class='label'>F_{p}, subgroup order {r}</text>\n"
        ));
        for t in 0..=4 {
            let y = bottom - t as f64 * (bottom - top) / 4.0;
            let value = max_s * t as f64 / 4.0;
            svg.push_str(&format!("<line x1='{left:.1}' y1='{y:.1}' x2='{:.1}' y2='{y:.1}' stroke='#e5e7eb'/><text x='{:.1}' y='{:.1}' class='small'>{value:.0}</text>\n", left+500.0, left+2.0, y-4.0));
        }
        for (j, v) in vals.iter().enumerate() {
            let x = left + 70.0 + j as f64 * 105.0;
            let mean_s = v.mean / (r as f64).sqrt();
            let y = bottom - mean_s / max_s * (bottom - top);
            let ymin = bottom - v.min / (r as f64).sqrt() / max_s * (bottom - top);
            let ymax = bottom - v.max / (r as f64).sqrt() / max_s * (bottom - top);
            svg.push_str(&format!("<rect x='{x:.1}' y='{y:.1}' width='47' height='{:.1}' fill='{}'/><line x1='{:.1}' y1='{ymin:.1}' x2='{:.1}' y2='{ymax:.1}' stroke='#111827' stroke-width='2'/><line x1='{:.1}' y1='{ymin:.1}' x2='{:.1}' y2='{ymin:.1}' stroke='#111827' stroke-width='2'/><line x1='{:.1}' y1='{ymax:.1}' x2='{:.1}' y2='{ymax:.1}' stroke='#111827' stroke-width='2'/><text x='{:.1}' y='{:.1}' class='small'>{mean_s:.2}</text><text x='{:.1}' y='412' class='small'>{}</text>\n", (bottom-y).max(0.0), colors[j], x+23.5, x+23.5, x+15.0, x+32.0, x+15.0, x+32.0, x, ymax-7.0, x-8.0, esc(labels[j])));
        }
    }
    svg.push_str("<text x='35' y='468' class='small'>S = (group additions + doublings) / sqrt(r). Destination arms include source-point cross-verification; ratios use paired totals.</text>\n</svg>\n");
    fs::write(&args[3], svg)?;

    let mut rho_svg = String::from("<svg xmlns='http://www.w3.org/2000/svg' width='1000' height='480' viewBox='0 0 1000 480'>\n<rect width='1000' height='480' fill='#fff'/>\n<style>text{font-family:Arial,sans-serif;fill:#172033}.title{font-size:20px;font-weight:bold}.label{font-size:15px}.small{font-size:12px}</style>\n<text x='35' y='35' class='title'>Exact-target rho controls: counted group-operation S</text>\n<text x='35' y='57' class='small'>v3: eight targets per order, observed ranges, logarithmic horizontal scale. Failed six-class attempts are charged.</text>\n");
    let x_of = |v: f64| {
        160.0
            + v.log10()
                .clamp(-0.3, 5.0)
                .mul_add(770.0 / 5.3, 0.3 * 770.0 / 5.3)
    };
    for power in 0..=5 {
        let x = x_of(10f64.powi(power));
        rho_svg.push_str(&format!("<line x1='{x:.1}' y1='88' x2='{x:.1}' y2='407' stroke='#e5e7eb'/><text x='{:.1}' y='431' class='small'>10^{power}</text>\n", x-12.0));
    }
    for (j, (p, arm, label, color)) in [
        (13, "rho_negation", "F13 rho A2", "#64748b"),
        (13, "rho_six_class", "F13 rho A6", "#d97706"),
        (199, "rho_negation", "F199 rho A2", "#64748b"),
        (199, "rho_six_class", "F199 rho A6", "#d97706"),
    ]
    .iter()
    .enumerate()
    {
        let (r, values) = &groups[&(*p, (*arm).into())];
        let v = stats(values);
        let root = (*r as f64).sqrt();
        let y = 125.0 + j as f64 * 78.0;
        let (x, xmin, xmax) = (x_of(v.mean / root), x_of(v.min / root), x_of(v.max / root));
        rho_svg.push_str(&format!("<text x='35' y='{:.1}' class='label'>{}</text><line x1='{xmin:.1}' y1='{y:.1}' x2='{xmax:.1}' y2='{y:.1}' stroke='{color}' stroke-width='3'/><line x1='{xmin:.1}' y1='{:.1}' x2='{xmin:.1}' y2='{:.1}' stroke='{color}' stroke-width='2'/><line x1='{xmax:.1}' y1='{:.1}' x2='{xmax:.1}' y2='{:.1}' stroke='{color}' stroke-width='2'/><circle cx='{x:.1}' cy='{y:.1}' r='6' fill='{color}'/><text x='{:.1}' y='{:.1}' class='small'>{:.3} S</text>\n", y+5.0, esc(label), y-7.0, y+7.0, y-7.0, y+7.0, (x+10.0).min(875.0), y-10.0, v.mean/root));
    }
    rho_svg.push_str("<text x='35' y='465' class='small'>A2 uses negation; A6 uses cubic action plus negation. These are fixed implementation controls, not a best-rho calibration.</text>\n</svg>\n");
    fs::write(&args[4], rho_svg)?;
    Ok(())
}

//! Render the recorded counted-group-operation data; no solver work is rerun.
use std::collections::BTreeMap;
use std::env;
use std::fs;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    assert_eq!(args.len(), 4, "plot RAW.tsv SUMMARY.tsv GRAPH.svg");
    let raw = fs::read_to_string(&args[1])?;
    let mut lines = raw.lines();
    let header: Vec<&str> = lines.next().ok_or("empty raw file")?.split('\t').collect();
    let column = |name: &str| -> Result<usize, String> {
        header
            .iter()
            .position(|&s| s == name)
            .ok_or(format!("missing {name}"))
    };
    let p_col = column("p")?;
    let r_col = column("r")?;
    let arm_col = column("arm")?;
    let status_col = column("status")?;
    let ops_col = column("total_gae")?;
    let phase_names = [
        "setup_gae",
        "collection_gae",
        "descent_gae",
        "verification_gae",
        "rho_gae",
        "fb_lifts",
        "pair_entries",
        "relation_trials",
        "relation_lookups",
        "independent_rows",
        "row_mul_subs",
        "row_inversions",
        "descent_trials",
        "descent_lookups",
        "field_muls",
        "field_sqrs",
        "field_invs",
        "map_evaluations",
        "map_muls",
        "map_sqrs",
        "map_invs",
    ];
    let phase_cols: Vec<usize> = phase_names
        .iter()
        .map(|name| column(name))
        .collect::<Result<_, _>>()?;
    let mut groups: BTreeMap<(u64, String), (u64, Vec<f64>)> = BTreeMap::new();
    let mut phases: BTreeMap<(u64, String), Vec<Vec<u64>>> = BTreeMap::new();
    for (index, line) in lines.enumerate() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != header.len() {
            return Err(format!(
                "row {} has {} columns, expected {}",
                index + 2,
                fields.len(),
                header.len()
            )
            .into());
        }
        let p = fields[p_col].parse::<u64>()?;
        let r = fields[r_col].parse::<u64>()?;
        let arm = fields[arm_col].to_owned();
        if fields[status_col] != "verified" {
            return Err(format!("row {} is not verified", index + 2).into());
        }
        let counts = phase_cols
            .iter()
            .map(|&col| fields[col].parse::<u64>())
            .collect::<Result<Vec<_>, _>>()?;
        phases.entry((p, arm.clone())).or_default().push(counts);
        groups
            .entry((p, arm))
            .or_insert_with(|| (r, Vec::new()))
            .1
            .push(fields[ops_col].parse::<f64>()?);
    }
    assert_eq!(groups.len(), 6, "expected two curves times three arms");
    let mut summary = String::from("p\tr\tarm\tn\tmean_gae\tmin_gae\tmax_gae\tmean_s\tmin_s\tmax_s\tmean_s_over_rho\tmean_s_over_floor\n");
    let mut panels = Vec::new();
    for p in [13, 199] {
        let mut panel = Vec::new();
        let rho_values = &groups
            .get(&(p, "rho_negation".to_string()))
            .ok_or("missing rho")?
            .1;
        let rho_mean = rho_values.iter().sum::<f64>() / rho_values.len() as f64;
        let floor_s = (std::f64::consts::PI / 4.0).sqrt();
        for arm in ["rho_negation", "source_relations", "transport_relations"] {
            let (r, values) = groups.get(&(p, arm.to_string())).ok_or("missing arm")?;
            assert_eq!(values.len(), 8, "expected eight targets in each arm");
            let mean = values.iter().sum::<f64>() / values.len() as f64;
            let min = values.iter().copied().fold(f64::INFINITY, f64::min);
            let max = values.iter().copied().fold(0.0, f64::max);
            let root = (*r as f64).sqrt();
            summary.push_str(&format!("{p}\t{r}\t{arm}\t8\t{mean:.6}\t{min:.0}\t{max:.0}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\n", mean / root, min / root, max / root, mean / rho_mean, (mean / root) / floor_s));
            panel.push((arm, mean / root, min / root, max / root));
        }
        let r = groups.get(&(p, "rho_negation".to_string())).unwrap().0;
        panels.push((p, r, panel));
    }
    fs::write(&args[2], summary)?;
    let mut phase_summary = format!(
        "p\tarm\tn\t{}\n",
        phase_names
            .iter()
            .map(|name| format!("mean_{name}"))
            .collect::<Vec<_>>()
            .join("\t")
    );
    for ((p, arm), rows) in &phases {
        phase_summary.push_str(&format!("{p}\t{arm}\t{}", rows.len()));
        for col in 0..phase_names.len() {
            let mean = rows.iter().map(|row| row[col] as f64).sum::<f64>() / rows.len() as f64;
            phase_summary.push_str(&format!("\t{mean:.6}"));
        }
        phase_summary.push('\n');
    }
    let summary_name = std::path::Path::new(&args[2])
        .file_name()
        .ok_or("summary filename")?
        .to_string_lossy();
    let phase_path = std::path::Path::new(&args[2]).with_file_name(format!("phase_{summary_name}"));
    fs::write(phase_path, phase_summary)?;

    let mut svg = String::from("<svg xmlns='http://www.w3.org/2000/svg' width='1000' height='480' viewBox='0 0 1000 480'>\n<rect width='1000' height='480' fill='#fff'/>\n<style>text{font-family:Arial,sans-serif;fill:#1f2937} .title{font-size:20px;font-weight:bold} .label{font-size:14px} .small{font-size:12px}</style>\n<text x='35' y='34' class='title'>Same-target complete solves: counted group-operation S</text>\n<text x='35' y='56' class='small'>Mean over 8 public targets; whisker = observed minimum to maximum. Other native costs are counted separately.</text>\n");
    let colors = ["#64748b", "#0284c7", "#be185d"];
    let labels = ["matched rho", "source relations", "transport relations"];
    for (side, (p, r, rows)) in panels.iter().enumerate() {
        let left = 45.0 + side as f64 * 490.0;
        let top = 105.0;
        let bottom = 375.0;
        let max = rows.iter().map(|v| v.3).fold(0.0, f64::max) * 1.12;
        svg.push_str(&format!(
            "<text x='{left:.1}' y='87' class='label'>F_{p}, subgroup order {r}</text>\n"
        ));
        for t in 0..=4 {
            let y = bottom - t as f64 * (bottom - top) / 4.0;
            let val = max * t as f64 / 4.0;
            svg.push_str(&format!("<line x1='{left:.1}' y1='{y:.1}' x2='{:.1}' y2='{y:.1}' stroke='#e5e7eb'/><text x='{:.1}' y='{:.1}' class='small'>{val:.0}</text>\n", left + 445.0, left + 2.0, y - 4.0));
        }
        for (j, row) in rows.iter().enumerate() {
            let x = left + 75.0 + j as f64 * 128.0;
            let y = bottom - (row.1 / max) * (bottom - top);
            let ymin = bottom - (row.2 / max) * (bottom - top);
            let ymax = bottom - (row.3 / max) * (bottom - top);
            svg.push_str(&format!("<rect x='{x:.1}' y='{y:.1}' width='54' height='{:.1}' fill='{}'/><line x1='{:.1}' y1='{ymin:.1}' x2='{:.1}' y2='{ymax:.1}' stroke='#111827' stroke-width='2'/><line x1='{:.1}' y1='{ymin:.1}' x2='{:.1}' y2='{ymin:.1}' stroke='#111827' stroke-width='2'/><line x1='{:.1}' y1='{ymax:.1}' x2='{:.1}' y2='{ymax:.1}' stroke='#111827' stroke-width='2'/><text x='{:.1}' y='{:.1}' class='small'>{:.2}</text><text x='{:.1}' y='398' class='small'>{}</text>\n", (bottom-y).max(0.0), colors[j], x+27.0, x+27.0, x+18.0, x+36.0, x+18.0, x+36.0, x, ymax-7.0, row.1, x-20.0, labels[j]));
        }
    }
    svg.push_str("<text x='35' y='445' class='small'>S = (group additions + doublings) / sqrt(r). The raw TSV carries map field work, lifts, lookups, and row arithmetic.</text>\n</svg>\n");
    fs::write(&args[3], svg)?;
    Ok(())
}

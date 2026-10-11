//! Freeze independent coordinate replay for every constructed higher-degree control.
use std::{collections::BTreeMap, fmt::Write, fs, path::Path};
fn kv(s: &str) -> BTreeMap<&str, &str> {
    s.split('|')
        .skip(1)
        .filter_map(|v| v.split_once('='))
        .collect()
}
fn main() {
    let a: Vec<_> = std::env::args().skip(1).collect();
    assert_eq!(a.len(), 1);
    let base = Path::new(&a[0]);
    let raw = fs::read_to_string(base.join("general_controls_run1/stdout.txt")).unwrap();
    assert!(fs::read(base.join("general_controls_run1/stderr.txt"))
        .unwrap()
        .is_empty());
    let mut gp = String::new();
    let mut p = "";
    let mut k = "";
    let mut m = "";
    let mut order = "";
    let mut ea = "";
    let mut eb = "";
    let mut controls = 0;
    let mut noncoprime = 0;
    for line in raw.lines() {
        let v = kv(line);
        if line.starts_with("GENERAL_CONTROL|") {
            p = v["p"];
            k = v["degree"];
            m = v["modulus"];
            order = v["order"];
            ea = v["short_a"];
            eb = v["short_b"];
            controls += 1;
            noncoprime += usize::from(v["norm_section_gcd"] != "1");
            writeln!(gp, "check_model({p},{k},{m},{ea},{eb},{order});").unwrap();
        }
        if line.starts_with("LITERAL_MINUS|") {
            writeln!(
                gp,
                "check_minus({p},{k},{m},{ea},{eb},{},{},{},{},{},{},{order});",
                v["alpha"], v["d"], v["twist"], v["h"], v["c3"], v["a2"]
            )
            .unwrap();
        }
    }
    assert_eq!((controls, noncoprime), (16, 8));
    gp.push_str("print(\"ALGEBRA_COMPLETE|higher_degree_controls=16|noncoprime_sections=8\");\n");
    fs::write(base.join("general_replay_inputs.gp"), gp).unwrap();
    println!("16 independent literal replays, including eight noncoprime sections");
}

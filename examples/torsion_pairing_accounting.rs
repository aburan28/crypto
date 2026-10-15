//! Native, fail-closed accounting over the frozen pairing-transfer run records.
//! This program never participates in the target solve.

use std::fs;
use std::path::Path;

use serde_json::{json, Value};

const FLOOR: f64 = 0.886_226_925_452_758;

fn read(path: &Path) -> Value {
    serde_json::from_str(&fs::read_to_string(path).unwrap()).unwrap()
}

fn u(v: &Value, path: &[&str]) -> u64 {
    let mut at = v;
    for key in path {
        at = &at[*key];
    }
    at.as_u64()
        .unwrap_or_else(|| panic!("missing integer: {path:?}"))
}

fn f(v: &Value, path: &[&str]) -> f64 {
    let mut at = v;
    for key in path {
        at = &at[*key];
    }
    at.as_f64()
        .unwrap_or_else(|| panic!("missing number: {path:?}"))
}

fn median(mut values: Vec<f64>) -> Value {
    values.sort_by(f64::total_cmp);
    json!({"min":values[0],"median":values[values.len()/2],
        "max":values[values.len()-1],"n":values.len()})
}

fn main() {
    let root = std::env::args().nth(1).expect("study directory");
    let root = Path::new(&root);
    let runs = root.join("runs");
    let status = fs::read_to_string(runs.join("RUN_STATUS.tsv")).unwrap();
    assert_eq!(status.lines().count(), 31, "expected 30 completed commands");
    assert!(status.lines().skip(1).all(|line| line.ends_with("\t0")));
    let fixtures =
        fs::read_to_string(root.join("../extension_field_20261014/FIXTURES.jsonl")).unwrap();
    let fixture_rows: Vec<Value> = fixtures
        .lines()
        .take(2)
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    let mut all = Vec::new();
    let mut groups = Vec::new();

    for cell in 0..2usize {
        let expected = u(&fixture_rows[cell], &["answer_scalar"]);
        let r = u(&fixture_rows[cell], &["r"]);
        let p = u(&fixture_rows[cell], &["p"]);
        let mut pair_s = Vec::new();
        let mut rho_s = Vec::new();
        let mut rho_aa_s = Vec::new();
        let mut pair_to_rho = Vec::new();
        let mut pair_cpu = Vec::new();
        let mut rho_cpu = Vec::new();
        let mut rho_aa_cpu = Vec::new();
        let mut pair_muls = Vec::new();
        let mut rho_muls = Vec::new();
        for round in 0..3usize {
            let name = |suffix: &str| runs.join(format!("c{cell}_r{round}_{suffix}.json"));
            let pair = read(&name("pairing"));
            let verify = read(&name("verify"));
            let bilinear = read(&name("bilinear"));
            let rho = read(&name("rho"));
            let rho_aa = read(&name("rho_aa"));
            for v in [&pair, &verify, &bilinear, &rho, &rho_aa] {
                assert_eq!(u(v, &["cell"]), cell as u64);
                assert_eq!(u(v, &["p"]), p);
                assert_eq!(u(v, &["r"]), r);
            }
            assert_eq!(u(&pair, &["recovered"]), expected);
            assert_eq!(u(&verify, &["recovered"]), expected);
            assert_eq!(verify["verified"], true);
            for check in bilinear["checks"].as_array().unwrap() {
                assert_eq!(check["passed"], true);
            }
            for control in [&rho, &rho_aa] {
                assert_eq!(u(control, &["rho", "recovered"]), expected);
                assert_eq!(control["rho"]["verified"], true);
                assert_eq!(
                    u(control, &["seed"]),
                    0x5041_4952_0000_1000u64 + 256 * cell as u64 + round as u64
                );
            }
            assert_eq!(rho["rho"]["gae"], rho_aa["rho"]["gae"]);
            assert_eq!(rho["fp_muls"], rho_aa["fp_muls"]);

            let setup = u(&pair, &["setup", "cpu_ns"]);
            let selection = u(&pair, &["selection_and_pairings", "cpu_ns"]);
            let field_dlp = u(&pair, &["multiplicative_dlp", "cpu_ns"]);
            let cold = u(&pair, &["cold_cpu_ns"]);
            let residual = cold.checked_sub(setup + selection + field_dlp).unwrap();
            let replay = u(&verify, &["cpu_ns"]);
            let bilinear_cpu = u(&bilinear, &["cpu_ns"]);
            let total = cold + replay + bilinear_cpu;
            let rho_ns = u(&rho, &["cpu_ns"]);
            let rho_aa_ns = u(&rho_aa, &["cpu_ns"]);
            let calibration = f(&rho, &["ns_per_source_add_calibration"]);
            let divisor = calibration * (r as f64).sqrt();
            let pair_scaled = total as f64 / divisor;
            let rho_scaled = rho_ns as f64 / divisor;
            let rho_aa_scaled = rho_aa_ns as f64 / divisor;
            let ratio = total as f64 / rho_ns as f64;
            let pair_fp_muls = u(&pair, &["setup", "fp_muls"])
                + u(&pair, &["selection_and_pairings", "fp_muls"])
                + u(&pair, &["multiplicative_dlp", "fp_muls"])
                + u(&verify, &["fp_muls"])
                + u(&bilinear, &["fp_muls"]);
            let pair_fp_invs = u(&pair, &["setup", "fp_invs"])
                + u(&pair, &["selection_and_pairings", "fp_invs"])
                + u(&pair, &["multiplicative_dlp", "fp_invs"])
                + u(&verify, &["fp_invs"])
                + u(&bilinear, &["fp_invs"]);
            let k_mul_setup = u(&pair, &["setup", "k_cumulative", "extension_muls"]);
            let k_mul_select = u(
                &pair,
                &["selection_and_pairings", "k_cumulative", "extension_muls"],
            );
            let k_mul_dlp = u(
                &pair,
                &["multiplicative_dlp", "k_cumulative", "extension_muls"],
            );
            let k_mul_bilinear = u(&bilinear, &["field_counts", "extension_muls"]);
            let peak = [&pair, &verify, &bilinear]
                .iter()
                .map(|v| u(v, &["peak_rss_raw"]))
                .max()
                .unwrap();
            assert_eq!(pair["peak_rss_unit"], "bytes");
            assert_eq!(verify["peak_rss_unit"], "bytes");
            assert_eq!(bilinear["peak_rss_unit"], "bytes");
            pair_s.push(pair_scaled);
            rho_s.push(rho_scaled);
            rho_aa_s.push(rho_aa_scaled);
            pair_to_rho.push(ratio);
            pair_cpu.push(total as f64);
            rho_cpu.push(rho_ns as f64);
            rho_aa_cpu.push(rho_aa_ns as f64);
            pair_muls.push(pair_fp_muls as f64);
            rho_muls.push(u(&rho, &["fp_muls"]) as f64);
            all.push(json!({"cell":cell,"round":round,"p":p,"r":r,
                "recovered":expected,"verified":true,
                "pairing_cpu_phases_ns":{"setup":setup,"selection_and_pairings":selection,
                    "field_dlp":field_dlp,"cold_residual":residual,
                    "independent_scalar_replay":replay,"bilinearity_checks":bilinear_cpu,
                    "full_total":total},
                "pairing_fp_muls":pair_fp_muls,"pairing_fp_invs":pair_fp_invs,
                "pairing_extension_muls":{"setup":k_mul_setup,
                    "selection_and_pairings":k_mul_select-k_mul_setup,
                    "field_dlp":k_mul_dlp-k_mul_select,
                    "bilinearity_checks":k_mul_bilinear},
                "pairing_peak_rss_bytes":peak,
                "rho_cpu_ns":rho_ns,"rho_aa_cpu_ns":rho_aa_ns,
                "rho_fp_muls":u(&rho,&["fp_muls"]),
                "rho_fp_invs":u(&rho,&["fp_invs"]),
                "rho_gae":f(&rho,&["rho","gae"]),
                "rho_counted_s":f(&rho,&["rho","s"]),
                "ns_per_source_add_calibration":calibration,
                "s_cpu_gae":{"pairing":pair_scaled,"rho":rho_scaled,
                    "rho_aa":rho_aa_scaled},
                "ratio_to_floor":{"pairing":pair_scaled/FLOOR,
                    "rho":rho_scaled/FLOOR,"rho_aa":rho_aa_scaled/FLOOR},
                "pairing_over_rho_cpu_ratio":ratio}));
        }
        groups.push(json!({"cell":cell,"p":p,"r":r,"answer":expected,
            "pairing_cpu_ns":median(pair_cpu),"rho_cpu_ns":median(rho_cpu),
            "rho_aa_cpu_ns":median(rho_aa_cpu),
            "pairing_s_cpu_gae":median(pair_s),"rho_s_cpu_gae":median(rho_s),
            "rho_aa_s_cpu_gae":median(rho_aa_s),
            "pairing_over_rho_cpu_ratio":median(pair_to_rho),
            "pairing_fp_muls":median(pair_muls),"rho_fp_muls":median(rho_muls)}));
    }
    println!("{}", serde_json::to_string_pretty(&json!({
        "schema":"torsion_pairing_accounting/v1",
        "run_commit":"47deb8dccaff455c41500a00b99bf208c5ab2a0c",
        "floor_counted_gae_per_sqrt_r":FLOOR,
        "unit":"current-thread CPU ns converted by the paired rho source-add calibration, divided by sqrt(r)",
        "floor_caveat":"the counting floor is not a hardware time lower bound",
        "pairing_total":"cold pairing run plus independent scalar replay plus post-solve bilinearity check; each is a separate process",
        "rows":all,"by_cell":groups
    })).unwrap());
}

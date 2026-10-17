//! Posthoc within-rho endpoint diagnostics; no IC/rho speed-score admission.
use super::*;
use std::collections::BTreeMap;

pub fn run() {
    let root = Path::new(ROOT);
    let mut endpoints = BTreeMap::new();
    let mut public_targets = BTreeMap::new();
    let mut sources = Vec::new();
    for shard in 0..8 {
        let dir = root.join(format!("sessions/primary-{shard}"));
        let seal = read(dir.join("session.json"));
        assert_eq!(seal["status"], "complete");
        let bytes = fs::read(dir.join("records.jsonl")).unwrap();
        let digest = sha256_hex(&bytes);
        assert_eq!(seal["records_sha256"], digest);
        sources.push(json!({"session":shard,"records_sha256":digest}));
        for line in std::str::from_utf8(&bytes).unwrap().lines() {
            let v: Value = serde_json::from_str(line).unwrap();
            assert_eq!(v["outcome"]["status"], "verified");
            let c = &v["workload"]["curve_spec"];
            let key = (
                c["p"].as_u64().unwrap(),
                c["a"].as_u64().unwrap(),
                c["b"].as_u64().unwrap(),
            );
            let public = json!([v["workload"]["workload_id"], v["workload"]["target"]]);
            if let Some(previous) = public_targets.insert(key, public.clone()) {
                assert_eq!(previous, public, "same-endpoint arms use different targets");
            }
            if v["arm"] == "rho" {
                let cpu = (v["time"]["user_ns"].as_u64().unwrap() as f64
                    + v["time"]["sys_ns"].as_u64().unwrap() as f64)
                    / 1e9;
                assert!(endpoints.insert(key, json!({"process_cpu_seconds":cpu,"charged_group_work_lower_bound":v["cost"]["total_gae"],"unpriced":v["cost"]["unpriced"],"workload":v["workload"]["workload_id"]})).is_none());
            }
        }
    }
    assert_eq!(endpoints.len(), 200);
    assert_eq!(public_targets.len(), 200);
    let inputs = read(root.join("inputs/pairs.json"));
    let mut pairs = Vec::new();
    for pair in inputs.as_array().unwrap() {
        let p = pair["p"].as_u64().unwrap();
        pairs.push(json!({"p":p,"ell":pair["ell"],"crater":endpoints[&(p,pair["a0"].as_u64().unwrap(),pair["b0"].as_u64().unwrap())],"floor":endpoints[&(p,pair["a1"].as_u64().unwrap(),pair["b1"].as_u64().unwrap())]}));
    }
    assert_eq!(pairs.len(), 100);
    let mut out = json!({"schema":"conductor-gap.rho-diagnostics/v1","analysis":"posthoc diagnostic of the completed matched-target rho panel","verified_attempts":200,"paired_gaps":100,"same_public_target_across_methods":true,"sources":sources,"pairs":pairs,"scope":"L0; endpoint scheduling waves confound CPU. Charged group work omits canonicalisation; its ratio is a lower-bound diagnostic, not an eligible complete-work or IC/rho speed score.","bootstrap":"10000 whole-gap paired resamples; xorshift seed 20261021; ratio of means"});
    for metric in ["process_cpu_seconds", "charged_group_work_lower_bound"] {
        let xs: Vec<_> = pairs
            .iter()
            .map(|v| {
                (
                    v["crater"][metric].as_f64().unwrap(),
                    v["floor"][metric].as_f64().unwrap(),
                )
            })
            .collect();
        let ratio = |indices: &[usize]| {
            indices.iter().map(|&i| xs[i].1).sum::<f64>()
                / indices.iter().map(|&i| xs[i].0).sum::<f64>()
        };
        let mut seed = 20261021u64;
        let mut draws = Vec::new();
        for _ in 0..10000 {
            let indices: Vec<_> = (0..xs.len())
                .map(|_| {
                    seed ^= seed << 13;
                    seed ^= seed >> 7;
                    seed ^= seed << 17;
                    (seed % xs.len() as u64) as usize
                })
                .collect();
            draws.push(ratio(&indices));
        }
        draws.sort_by(f64::total_cmp);
        out[format!("floor_over_crater_{metric}")] = json!({"ratio_of_means":ratio(&(0..xs.len()).collect::<Vec<_>>()),"paired_percentile_95":[draws[250],draws[9749]],"crater_mean":xs.iter().map(|v|v.0).sum::<f64>()/xs.len() as f64,"floor_mean":xs.iter().map(|v|v.1).sum::<f64>()/xs.len() as f64});
    }
    write(root.join("results/paired-rho.json"), &out);
}

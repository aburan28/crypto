//! Audit a retained native integration-control packet and catalogue plan.
//! Build with rustc and the dedicated CLI's locked serde_json dependency.
use serde_json::{json, Value};
use std::{fs, path::Path};
fn read(p: &Path) -> Value {
    serde_json::from_slice(&fs::read(p).unwrap()).unwrap()
}
fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert_eq!(
        args.len(),
        4,
        "audit CONTROL_DIR POLICY_AND_PLAN_DIR REPORT_DIR"
    );
    let controls = Path::new(&args[1]);
    let data = Path::new(&args[2]);
    let out = Path::new(&args[3]);
    fs::create_dir_all(out).unwrap();
    let policy = read(&data.join("policy.json"));
    let plan = read(&data.join("plan.json"));
    let freeze = read(&data.join("freeze.json"));
    assert_eq!(policy["tool"]["git_commit"], freeze["source_commit"]);
    assert_eq!(policy["tool"]["git_dirty"], "false");
    assert_eq!(policy["season"], plan["season"]);
    let rows = policy["curves"].as_array().unwrap();
    let planned = plan["curves"].as_array().unwrap();
    assert_eq!(rows.len(), planned.len());
    let mut construction = 0;
    let mut replay = 0;
    let mut candidates = 0;
    let mut roster = vec![];
    for (row, p) in rows.iter().zip(planned) {
        assert!(!row["standard_names"].as_array().unwrap().is_empty());
        assert_eq!(row["source_icv1"], p["source_icv1"]);
        construction += usize::from(row["capabilities"]["construction"] == true);
        replay += usize::from(row["capabilities"]["independent_replay"] == true);
        let c = p["candidates"].as_array().unwrap();
        candidates += c.len();
        for pair in c.windows(2) {
            assert!(degree(&pair[0]) > degree(&pair[1]));
        }
        for c in c {
            assert_eq!(c["proof_status"], "CANDIDATE_ONLY");
            assert!(degree(c) > 10453);
        }
        roster.push(json!({"source_icv1":row["source_icv1"], "standard_names":row["standard_names"],
            "construction":row["capabilities"]["construction"], "independent_replay":row["capabilities"]["independent_replay"],
            "adapter_detail":row["capabilities"]["detail"],"screen_status":p["screen_status"],"candidate_count":c.len(),
            "highest_candidate":c.first().map(|v| &v["ell"])}));
    }
    let mut negatives = vec![];
    for label in [
        "composite",
        "inflated-degree",
        "source",
        "kernel",
        "map",
        "compact-only",
        "traversal",
        "missing-hash",
        "duplicate-kernel",
    ] {
        if label == "duplicate-kernel"
            && !controls.join("judged-duplicate-kernel/entry.json").exists()
        {
            continue;
        }
        let entry = read(&controls.join(format!("judged-{label}/entry.json")));
        assert_eq!(entry["status"], "REJECTED", "{label}");
        assert!(entry["score"].is_null());
        negatives.push(json!({"control":label,"status":entry["status"],"reason":entry["reason"]}));
    }
    let first = read(&controls.join("p192/RESULT.json"));
    let second = read(&controls.join("p521/RESULT.json"));
    assert_eq!(first["score"], "13");
    assert_eq!(second["score"], "7");
    assert_ne!(first["source_icv1"], second["source_icv1"]);
    let board = read(&controls.join("board/leaderboard.json"));
    assert_eq!(board["curves"].as_array().unwrap().len(), rows.len());
    assert_eq!(board["attempts"][2]["status"], "REJECTED");
    let winners: Vec<_> = board["curves"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|v| v["status"] == "VERIFIED")
        .collect();
    assert_eq!(winners.len(), 2);
    let mut commands = vec![];
    for n in 0..1000 {
        let f = controls.join(format!("command-{n:03}.json"));
        if !f.exists() {
            break;
        }
        let c = read(&f);
        assert!(!c["exit_code"].is_null());
        commands.push(json!({"file":f.file_name().unwrap().to_str(),"args":c["args"],"exit_code":c["exit_code"]}));
    }
    assert!(commands.len() >= 15);
    let mut stress = vec![];
    if data.join("resource-stress").is_dir() {
        for n in 0..20 {
            let entry = read(&data.join(format!("resource-stress/{n}/entry.json")));
            assert!(entry["status"] == "ACCEPTED" || entry["status"] == "CENSORED");
            if entry["status"] == "ACCEPTED" {
                assert!(
                    entry["verification"]["sampled_peak_resident_bytes"]
                        .as_str()
                        .unwrap()
                        .parse::<u64>()
                        .unwrap()
                        > 0
                );
            } else {
                assert!(entry["score"].is_null());
            }
            stress.push(entry);
        }
    }
    let stress_accepted = stress.iter().filter(|s| s["status"] == "ACCEPTED").count();
    let mut starting = vec![];
    let mut anchor_attempts = vec![];
    if data.join("anchor-board/leaderboard.json").exists() {
        let board = read(&data.join("anchor-board/leaderboard.json"));
        assert_eq!(board["season"], policy["season"]);
        assert_eq!(board["curves"].as_array().unwrap().len(), rows.len());
        for (i, entry) in board["attempts"].as_array().unwrap().iter().enumerate() {
            if entry["status"] == "ACCEPTED" {
                assert_eq!(entry["tool"]["git_commit"], freeze["source_commit"]);
                assert_eq!(entry["judge_binary_sha256"], freeze["binary_sha256"]);
                let verified =
                    read(&data.join(format!("anchor-board/entry-{i}/verification.json")));
                let maps = verified["fresh_certificate"]["records"].as_array().unwrap();
                assert_eq!(verified["status"], "PASS");
                for map in maps {
                    assert_eq!(map["status"], "PASS");
                    assert_eq!(map["source_icv1"], entry["source_icv1"]);
                    assert_eq!(
                        map["ell"].as_u64().unwrap().to_string(),
                        entry["prime_degree"].as_str().unwrap()
                    );
                    assert_eq!(map["public_scalar_transport_checks"], 20);
                }
                starting.push(json!({"source_icv1":entry["source_icv1"],"prime_degree":entry["prime_degree"],
                    "coverage":entry["degree_coverage"],"map_count":maps.len(),"map_sha256":entry["map_sha256"],
                    "targets":maps.iter().map(|m|m["target_icv1"].clone()).collect::<Vec<_>>()}));
            } else {
                assert!(entry["score"].is_null());
            }
            anchor_attempts.push(entry.clone());
        }
    }
    let short_budget = if data.join("anchor-short-budget/entry.json").exists() {
        let entry = read(&data.join("anchor-short-budget/entry.json"));
        if entry["status"] == "CENSORED" {
            assert!(entry["score"].is_null());
        }
        entry
    } else {
        Value::Null
    };
    let audit = json!({"schema":"isogeny-degree-challenge-audit/v1","status":"PASS","season":policy["season"],
        "registered_prime_models":rows.len(),"construction_models":construction,"replay_models":replay,
        "unavailable_replay_models":rows.len()-replay,"shortlisted_candidates":candidates,
        "candidate_degree_interval":[plan["lower"].clone(),plan["upper"].clone()],
        "known_answer_controls":[first,second],"rejected_mutations":negatives,
        "forged_score_rejected":true,"verified_source_rows":winners.len(),"commands":commands,"roster":roster,
        "resource_stress_count":stress.len(),"resource_stress":stress,
        "resource_stress_accepted":stress_accepted,"resource_stress_censored":stress.len()-stress_accepted,
        "starting_curve_rows":starting,"anchor_attempts":anchor_attempts,"short_budget_control":short_budget,
        "timing_claim":null,"new_discovery_count":0});
    fs::write(
        out.join("audit.full.json"),
        serde_json::to_vec_pretty(&audit).unwrap(),
    )
    .unwrap();
    let mut summary = audit.clone();
    for key in ["commands", "roster", "resource_stress", "anchor_attempts"] {
        summary.as_object_mut().unwrap().remove(key);
    }
    summary["command_count"] = json!(commands.len());
    summary["source_commit"] = freeze["source_commit"].clone();
    summary["binary_sha256"] = freeze["binary_sha256"].clone();
    summary["known_answer_controls"] = json!([&first, &second]
        .iter()
        .map(|c| json!({
        "source_icv1":c["source_icv1"],"prime_degree":c["prime_degree"],"status":c["status"]}))
        .collect::<Vec<_>>());
    summary["anchor_attempt_statuses"] = json!(anchor_attempts
        .iter()
        .map(|a| a["status"].clone())
        .collect::<Vec<_>>());
    summary["short_budget_control"] = if short_budget.is_null() {
        Value::Null
    } else {
        json!({"status":short_budget["status"],"reason":short_budget["reason"],"score":short_budget["score"]})
    };
    fs::write(
        out.join("audit.json"),
        serde_json::to_vec_pretty(&summary).unwrap(),
    )
    .unwrap();
    let mut report = format!("# Native admission observations\n\nSource `{}`; run started `{}`. Exact native counts from the pinned policy, plan, raw control entries and fresh certificates.\n\n| Observation | Count |\n|---|---:|\n| Standard prime models | {} |\n| Construction adapters | {} |\n| Replay adapters | {} |\n| Candidate curve-degree pairs | {} |\n| Corrupted-artifact controls rejected | {} |\n| Resource stress accepted | {} / {} |\n| Resource stress censored | {} |\n| Newly discovered maps | 0 |\n\nP192/13 and P521/7 were admitted as separate known-answer controls. All listed altered artifacts and the forged saved-score submission were rejected. Every leaderboard build performed fresh replay.\n",
        freeze["source_commit"].as_str().unwrap(),freeze["started_utc"].as_str().unwrap(),rows.len(),construction,replay,candidates,negatives.len(),stress_accepted,stress.len(),stress.len()-stress_accepted);
    if !starting.is_empty() {
        report.push_str("\nThe starting leaderboard re-admitted these historical maps. All passed twenty public scalar transports per map.\n\n| Exact source | Degree | Exact target | Coverage |\n|---|---:|---|---|\n");
        for row in &starting {
            for target in row["targets"].as_array().unwrap() {
                report.push_str(&format!(
                    "| `{}` | {} | `{}` | {} |\n",
                    row["source_icv1"].as_str().unwrap(),
                    row["prime_degree"].as_str().unwrap(),
                    target.as_str().unwrap(),
                    row["coverage"].as_str().unwrap()
                ));
            }
        }
        report.push_str(&format!("\nThe one-second P192/10453 control returned `{}` / `{}` with no score. The full-budget leaderboard used a separate fresh replay.\n",short_budget["status"].as_str().unwrap(),short_budget["reason"].as_str().unwrap()));
    }
    report.push_str("\n[Native audit summary](audit.json), [complete audit](audit.full.json.gz), [policy](policy.json.gz), [candidate plan](plan.json.gz), [raw packet](raw.tar.gz), [diagram](admission.svg), [visual report](report.pdf). Compressed JSON files preserve their original bytes; the raw archive also contains every input, map, certificate and attempt receipt. Raw receipts are shared-host operational records; no speed comparison is made.\n");
    fs::write(out.join("observations.md"), report).unwrap();
    // Diagram of admission: unlike a timing plot, these counts describe coverage only.
    let mut svg = String::from("<svg xmlns='http://www.w3.org/2000/svg' width='900' height='355' viewBox='0 0 900 355'><rect width='900' height='355' fill='#f5f7fa'/><g font-family='Arial' fill='#172b4d'><text x='30' y='34' font-size='21'>Prime Degree Frontier: native admission controls</text>");
    let boxes = [
        (30, 70, format!("{} standard prime models", rows.len())),
        (330, 70, format!("{} replay-capable models", replay)),
        (630, 70, "Fresh exact replay".into()),
        (330, 190, "P192 / 13 accepted".into()),
        (630, 190, "P521 / 7 accepted".into()),
    ];
    for (x, y, label) in boxes {
        svg.push_str(&format!("<rect x='{x}' y='{y}' width='240' height='60' rx='9' fill='white' stroke='#2463a0'/><text x='{}' y='{}' font-size='17'>{label}</text>",x+12,y+36));
    }
    svg.push_str("<path d='M270 100H330 M570 100H630 M750 130V190 M750 160H450V190' fill='none' stroke='#2463a0' stroke-width='2'/><text x='30' y='293' font-size='17'>{} altered-artifact controls rejected; forged leaderboard score rejected.</text><text x='30' y='324' font-size='15'>Public fixtures; separate source-curve scores. Candidate screens earn no score.</text></g></svg>");
    let svg = svg.replace("{}", &negatives.len().to_string());
    fs::write(out.join("admission.svg"), svg).unwrap();
    println!("{}",serde_json::to_string(&json!({"models":rows.len(),"construction":construction,"replay":replay,"candidates":candidates,"commands":commands.len()})).unwrap());
}
fn degree(v: &Value) -> u64 {
    v["ell"].as_str().unwrap().parse().unwrap()
}

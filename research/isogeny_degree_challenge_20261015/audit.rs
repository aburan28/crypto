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
    assert_eq!(policy["season"], plan["season"]);
    let rows = policy["curves"].as_array().unwrap();
    let planned = plan["curves"].as_array().unwrap();
    assert_eq!(rows.len(), planned.len());
    let mut construction = 0;
    let mut replay = 0;
    let mut candidates = 0;
    let mut roster = vec![];
    for (row, p) in rows.iter().zip(planned) {
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
    ] {
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
    let audit = json!({"schema":"isogeny-degree-challenge-audit/v1","status":"PASS","season":policy["season"],
        "registered_prime_models":rows.len(),"construction_models":construction,"replay_models":replay,
        "unavailable_replay_models":rows.len()-replay,"shortlisted_candidates":candidates,
        "candidate_degree_interval":[plan["lower"].clone(),plan["upper"].clone()],
        "known_answer_controls":[first,second],"rejected_mutations":negatives,
        "forged_score_rejected":true,"verified_source_rows":winners.len(),"commands":commands,"roster":roster,
        "timing_claim":null,"new_discovery_count":0});
    fs::write(
        out.join("audit.json"),
        serde_json::to_vec_pretty(&audit).unwrap(),
    )
    .unwrap();
    // Diagram of admission: unlike a timing plot, these counts describe coverage only.
    let mut svg = String::from("<svg xmlns='http://www.w3.org/2000/svg' width='900' height='355' viewBox='0 0 900 355'><rect width='900' height='355' fill='#f5f7fa'/><g font-family='Arial' fill='#172b4d'><text x='30' y='34' font-size='21'>Prime Degree Frontier: native admission controls</text>");
    let boxes = [
        (30, 70, format!("{} registered prime models", rows.len())),
        (330, 70, format!("{} replay-capable models", replay)),
        (630, 70, "Fresh exact replay".into()),
        (330, 190, "P192 / 13 accepted".into()),
        (630, 190, "P521 / 7 accepted".into()),
    ];
    for (x, y, label) in boxes {
        svg.push_str(&format!("<rect x='{x}' y='{y}' width='240' height='60' rx='9' fill='white' stroke='#2463a0'/><text x='{}' y='{}' font-size='17'>{label}</text>",x+12,y+36));
    }
    svg.push_str("<path d='M270 100H330 M570 100H630 M750 130V190 M750 160H450V190' fill='none' stroke='#2463a0' stroke-width='2'/><text x='30' y='293' font-size='17'>8 altered-artifact controls rejected; forged leaderboard score rejected.</text><text x='30' y='324' font-size='15'>Public fixtures; separate source-curve scores. Candidate screens earn no score.</text></g></svg>");
    fs::write(out.join("admission.svg"), svg).unwrap();
    println!("{}",serde_json::to_string(&json!({"models":rows.len(),"construction":construction,"replay":replay,"candidates":candidates,"commands":commands.len()})).unwrap());
}
fn degree(v: &Value) -> u64 {
    v["ell"].as_str().unwrap().parse().unwrap()
}

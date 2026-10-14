//! Freeze, execute structural controls, and analyze the ten-question campaign.
//! Measurement is exclusively through the existing native `ecbench` binary.
#[path = "ecbench_open_questions/analyze.rs"]
mod analyze;
#[path = "ecbench_open_questions/controls.rs"]
mod controls;

use serde_json::{json, Value};
use std::{fs, path::Path};

const SEED: u64 = 20261014;

fn save(path: &Path, value: &Value) -> Result<(), String> {
    if path.exists() {
        return Err(format!("refusing to overwrite {}", path.display()));
    }
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    fs::write(
        path,
        format!("{}\n", serde_json::to_string_pretty(value).unwrap()),
    )
    .map_err(|e| e.to_string())
}

fn arm(name: &str, role: &str, id: &str, params: Value) -> Value {
    json!({"name":name,"role":role,"method":{"id":id,"params":params}})
}

fn rho(koblitz: bool) -> Vec<Value> {
    let id = if koblitz {
        "rho.signed_frobenius_strong"
    } else {
        "rho.negation"
    };
    vec![
        arm("rho", "reference", id, json!({})),
        arm("rho-aa", "control", id, json!({})),
    ]
}

fn ic(name: &str, base: &str, oracle: &str, trials: u64, linalg: &str) -> Value {
    arm(
        name,
        "candidate",
        "ic.pipeline_counted",
        json!({
            "factor_base":base,"oracle":oracle,"max_trials":trials.to_string(),
            "targets":"walk","linalg":linalg,"solver_budget_seconds":"0"
        }),
    )
}

fn spec(label: &str, curves: Value, arms: Vec<Value>, targets: u64, rounds: u64) -> Value {
    json!({
        "schema":"ecbench.spec/v1", "label":label,
        "description":"PROTOCOL-v1.md. Public single-target diagnostics; full cold counts and unpriced counters retained. No wall-time or asymptotic claim from L0 measurements.",
        "workloads":{"curves":curves,"targets_per_curve":targets,"target_seed":SEED,"target_kind":"public"},
        "arms":arms,
        "measurement":{"rounds":rounds,"warmup":0,"order":"alternate","seed":SEED,"isolation_required":"L2","timeout_seconds":180}
    })
}

fn freeze(root: &Path) -> Result<(), String> {
    let dir = root.join("specs");
    let mut a = rho(false);
    for size in [32, 64] {
        a.push(ic(
            &format!("abs{size}-subtract2"),
            &format!("prime-abscissa:size={size}"),
            "rootfind-subtract:m=2",
            100000,
            "incremental-gauss",
        ));
    }
    a.push(ic(
        "abs64-mitm2",
        "prime-abscissa:size=64",
        "mitm:m=2,negation_folded=1",
        100000,
        "incremental-gauss",
    ));
    a.push(ic(
        "abs64-subtract3",
        "prime-abscissa:size=64",
        "rootfind-subtract:m=3",
        100000,
        "incremental-gauss",
    ));
    save(
        &dir.join("q01-prime.json"),
        &spec(
            "Q1/Q10 prime factor-base and complete-pipeline controls",
            json!([16, 18, 20, 22]
                .map(|bits| json!({"kind":"prime_search","bits":bits,"seed":59297}))),
            a,
            4,
            3,
        ),
    )?;

    let mut a = rho(false);
    for d in [3, 4] {
        a.push(ic(
            &format!("sub{d}-mitm3"),
            &format!("binary-subspace:dimension={d}"),
            "mitm:m=3",
            20000,
            "incremental-gauss",
        ));
    }
    save(
        &dir.join("q02-binary.json"),
        &spec(
            "Q2/Q6/Q10 ordinary binary baseline",
            json!([13, 17, 19, 23]
                .map(|n| json!({"kind":"binary_random","n":n,"seed":59297,"max_cofactor":8}))),
            a,
            4,
            2,
        ),
    )?;

    let mut a = rho(true);
    for engine in ["inherited-f4", "f6-ic"] {
        a.push(ic(
            engine,
            "koblitz-standard-subspace:dimension=3",
            &format!("pdp3-koblitz:m=3,engine={engine},degree=4,node_budget=4096"),
            8,
            "incremental-gauss",
        ));
    }
    a.push(ic(
        "rootfind-control",
        "koblitz-standard-subspace:dimension=3",
        "pdp3-koblitz:m=3,engine=rootfind",
        8,
        "incremental-gauss",
    ));
    save(&dir.join("q02-solvers.json"), &spec("Q2 solver diagnostics at fixed summand dimension3; intermediate field variables retained", json!([13,19,23,31].map(|n|json!({"kind":"koblitz","a":0,"n":n}))), a, 2, 1))?;

    let mut a = rho(true);
    for oracle in ["mitm:m=3,negation_folded=1", "mitm-frobenius:m=3"] {
        a.push(ic(
            if oracle.starts_with("mitm-f") {
                "orbit-folded"
            } else {
                "direct"
            },
            "compact-orbit-scan:columns=4,raw_x_cap=4096",
            oracle,
            256,
            "incremental-gauss",
        ));
    }
    save(
        &dir.join("q03-frobenius.json"),
        &spec(
            "Q3 same-point direct versus Frobenius orbit addition",
            json!([19, 23, 31, 37].map(|n| json!({"kind":"koblitz","a":0,"n":n}))),
            a,
            4,
            2,
        ),
    )?;

    let mut a = rho(true);
    for w in [1, 2, 3] {
        a.push(ic(
            &format!("gamma-weight-{w}"),
            &format!("cyclotomic:w={w},columns=8"),
            "mitm-frobenius:m=3",
            256,
            "incremental-gauss",
        ));
    }
    a.push(ic(
        "compact-control",
        "compact-orbit-scan:columns=8,raw_x_cap=4096",
        "mitm-frobenius:m=3",
        256,
        "incremental-gauss",
    ));
    save(
        &dir.join("q04-nonlinear.json"),
        &spec(
            "Q4 nonlinear gamma-weight bases; gamma weight is not a normal-basis Hamming ideal",
            json!([11, 23].map(|n| json!({"kind":"koblitz","a":1,"n":n}))),
            a,
            4,
            2,
        ),
    )?;

    let mut a = rho(true);
    for (m, lp) in [("3", 0), ("4", 0), ("6", 0), ("n-1", 0), ("3", 1), ("3", 2)] {
        a.push(arm(
            &format!("m{m}-lp{lp}"),
            "candidate",
            "ic.large_prime",
            json!({
                "small_dimension":"3","envelope_dimension":"4","summands":m,
                "large_primes":lp.to_string(),"max_trials":"256","max_states":"100000"
            }),
        ));
    }
    save(
        &dir.join("q06-many-summand.json"),
        &spec(
            "Q6/Q7 exact multi-summand and large-prime policies; n-1 is literal ten",
            json!([{"kind":"koblitz","a":1,"n":11}]),
            a,
            4,
            2,
        ),
    )?;

    let mut a = rho(true);
    for (name, stop) in [
        ("target-pin", "incremental-gauss"),
        ("full-rank", "incremental-gauss-full-rank"),
        ("checked-rank", "incremental-gauss-full-rank-checked"),
    ] {
        a.push(ic(
            name,
            "compact-orbit-scan:columns=4,raw_x_cap=4096",
            "mitm-frobenius:m=3",
            256,
            stop,
        ));
    }
    a.push(arm("shared-rank", "candidate", "ic.shared_rank", json!({"factor_base":"compact-orbit-scan:columns=4,raw_x_cap=4096","rank_seed":SEED.to_string(),"rank_max_trials":"256","target_max_attempts":"256"})));
    save(
        &dir.join("q07-rank.json"),
        &spec(
            "Q7 useful rank and target-blind reusable setup",
            json!([19, 23].map(|n| json!({"kind":"koblitz","a":0,"n":n}))),
            a,
            4,
            2,
        ),
    )?;

    let mut a = rho(false);
    a.push(ic(
        "abscissa-mitm",
        "prime-abscissa:size=64",
        "mitm:m=2,negation_folded=1",
        100000,
        "incremental-gauss",
    ));
    a.push(ic(
        "abscissa-subtract",
        "prime-abscissa:size=64",
        "rootfind-subtract:m=2",
        100000,
        "incremental-gauss",
    ));
    save(
        &dir.join("q08-isogeny.json"),
        &spec(
            "Q8 crater/floor diagnostic without fabricated isogeny transport",
            json!(["crater", "floor"].map(
                |level| json!({"kind":"conductor_gap","ell":263,"t":241,"level":level,"seed":59297})
            )),
            a,
            4,
            3,
        ),
    )?;

    // This input is deliberately retained as an applicability gate. It must
    // fail in plan/dispatch rather than masquerade as a degree-83 IC run.
    save(
        &dir.join("q10-n83-gate.json"),
        &spec(
            "Q10 mandatory N83 IC applicability gate",
            json!([{"kind":"koblitz","a":0,"n":83}]),
            vec![ic(
                "ic-n83",
                "compact-orbit-scan:columns=4,raw_x_cap=4096",
                "mitm-frobenius:m=3",
                256,
                "incremental-gauss",
            )],
            1,
            1,
        ),
    )?;
    Ok(())
}

fn main() -> Result<(), String> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let (command, root) = match args.as_slice() {
        [command, root] => (command.as_str(), Path::new(root)),
        _ => return Err("usage: ecbench_open_questions freeze|controls|analyze STUDY_DIR".into()),
    };
    match command {
        "freeze" => freeze(root),
        "controls" => save(&root.join("controls.json"), &controls::run()?),
        "analyze" => analyze::run(root),
        _ => Err(format!("unknown command {command}")),
    }
}

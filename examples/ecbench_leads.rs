//! Freeze the follow-up and analyze only ecbench's sealed records.
#[path = "ecbench_open_questions/analyze.rs"]
mod analyze;

use serde_json::{json, Value};
use std::{collections::BTreeMap, fmt::Write as _, fs, path::Path};

fn save(path: &Path, value: &Value) -> Result<(), String> {
    if path.exists() {
        return Err(format!("refusing to overwrite {}", path.display()));
    }
    fs::create_dir_all(path.parent().ok_or("missing parent")?).map_err(|e| e.to_string())?;
    fs::write(
        path,
        format!(
            "{}\n",
            serde_json::to_string_pretty(value).map_err(|e| e.to_string())?
        ),
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

fn pipeline(name: &str, base: &str, oracle: &str, cap: u64) -> Value {
    arm(
        name,
        "candidate",
        "ic.pipeline_counted",
        json!({"factor_base":base,"oracle":oracle,
        "linalg":"incremental-gauss","targets":"walk","max_trials":cap.to_string(),"solver_budget_seconds":"0"}),
    )
}

fn orbit(name: &str, columns: u64, mode: &str, cap: u64) -> Value {
    arm(
        name,
        "candidate",
        "ic.orbit_rank_v1",
        json!({"columns":columns.to_string(),
        "raw_x_cap":"4096","max_trials":cap.to_string(),"table_mode":mode,"max_pair_visits":"1000000"}),
    )
}

fn spec(label: &str, curves: Value, arms: Vec<Value>, targets: u64) -> Value {
    json!({"schema":"ecbench.spec/v1","label":label,"description":"PROTOCOL-v1.md; fresh public single-target workloads, cold bounded counts; no wall or asymptotic claim",
        "workloads":{"curves":curves,"targets_per_curve":targets,"target_seed":20261015,"target_kind":"public"},
        "measurement":{"rounds":1,"warmup":0,"order":"alternate","seed":20261015,"isolation_required":"L2","timeout_seconds":180},"arms":arms})
}

fn freeze(root: &Path) -> Result<(), String> {
    let dir = root.join("specs");
    let mut arms = rho(false);
    for size in [16, 32, 64, 128] {
        arms.push(pipeline(
            &format!("abs{size}-mitm2"),
            &format!("prime-abscissa:size={size}"),
            "mitm:m=2,negation_folded=1",
            100000,
        ));
    }
    arms.push(pipeline(
        "abs64-subtract3",
        "prime-abscissa:size=64",
        "rootfind-subtract:m=3",
        100000,
    ));
    save(
        &dir.join("prime-balance.json"),
        &spec(
            "Fresh-target prime factor-base balance",
            json!([16, 18, 20, 22]
                .map(|bits| json!({"kind":"prime_search","bits":bits,"seed":59297}))),
            arms,
            8,
        ),
    )?;
    let mut arms = rho(true);
    for engine in ["inherited-f4", "f6-ic"] {
        for (version, oracle) in [("v1", "pdp3-koblitz"), ("v2", "pdp3-koblitz-v2")] {
            arms.push(pipeline(
                &format!("{engine}-{version}"),
                "koblitz-standard-subspace:dimension=3",
                &format!("{oracle}:m=3,engine={engine},degree=4,node_budget=4096"),
                8,
            ));
        }
    }
    save(
        &dir.join("binary-shape.json"),
        &spec(
            "Versioned chained-system shape comparison",
            json!([13, 19, 23, 31].map(|n| json!({"kind":"koblitz","a":0,"n":n}))),
            arms,
            2,
        ),
    )?;
    let mut arms = rho(true);
    arms.push(orbit("exact-c4", 4, "direct", 2048));
    arms.push(orbit("orbit-c4", 4, "orbit", 2048));
    arms.push(pipeline(
        "existing-orbit-c4",
        "compact-orbit-scan:columns=4,raw_x_cap=4096",
        "mitm-frobenius:m=3",
        2048,
    ));
    save(
        &dir.join("orbit-ladder.json"),
        &spec(
            "Width-generic exact and quotient rank ladder",
            json!([13, 19, 23, 31].map(|n| json!({"kind":"koblitz","a":0,"n":n}))),
            arms,
            4,
        ),
    )?;
    let mut arms = rho(true);
    arms.push(orbit("orbit-c4", 4, "orbit", 2048));
    arms.push(orbit("orbit-c8", 8, "orbit", 2048));
    save(
        &dir.join("orbit-heldout.json"),
        &spec(
            "Held-out degree37 rank/yield diagnostic",
            json!([{"kind":"koblitz","a":0,"n":37}]),
            arms,
            4,
        ),
    )?;
    save(
        &dir.join("wide-n83.json"),
        &spec(
            "Exact E0 N83 ordinary-IC bounded diagnostic",
            json!([{"kind":"koblitz","a":0,"n":83}]),
            vec![
                orbit("orbit-c4", 4, "orbit", 256),
                orbit("orbit-c16", 16, "orbit", 256),
            ],
            2,
        ),
    )?;
    save(
        &dir.join("wide-legacy-refusal.json"),
        &spec(
            "Preserve previous wide pipeline applicability control",
            json!([{"kind":"koblitz","a":0,"n":83}]),
            vec![pipeline(
                "legacy-ic",
                "compact-orbit-scan:columns=4,raw_x_cap=4096",
                "mitm-frobenius:m=3",
                256,
            )],
            1,
        ),
    )?;
    Ok(())
}

fn shapes(root: &Path) -> Result<(), String> {
    use crypto_lib::binary_ecc::F2mElement;
    use crypto_lib::cryptanalysis::{
        ic_boundary::koblitz_instance,
        koblitz_groebner::{build_decomposition_system, FieldStructure},
    };
    let mut rows = Vec::new();
    for n in [13, 19, 23, 31] {
        let inst = koblitz_instance(0, n).ok_or("missing curve")?;
        let kc = inst.koblitz.as_ref().ok_or("missing Koblitz data")?;
        let st = FieldStructure::new(n, &kc.curve.irreducible);
        let basis: Vec<_> = (0..3)
            .map(|i| F2mElement::from_bit_positions(&[i], n))
            .collect();
        for m in [2, 3] {
            let crypto_lib::binary_ecc::BinaryPoint::Affine { x, .. } = kc.generator() else {
                return Err("infinite generator".into());
            };
            let sys = build_decomposition_system(&basis, x, &kc.curve.b, m, &st)
                .ok_or("presentation refused")?;
            let expected_vars = m * 3 + (m - 2) * n as usize;
            let expected_eq = (m - 1) * n as usize;
            if sys.n_vars != expected_vars || sys.equations.len() != expected_eq {
                return Err("actual presentation mismatch".into());
            }
            rows.push(
                json!({"curve":inst.name,"n":n,"ell":3,"m":m,"constructed_vars":sys.n_vars,
                "constructed_equations":sys.equations.len(),"legacy_vars":m*3,"legacy_equations":n,
                "v2_matches_actual":true}),
            );
        }
    }
    save(
        &root.join("shape-controls.json"),
        &json!({"schema":"ecbench-leads-shapes/v1","rows":rows,"ok":true}),
    )
}

fn read(path: impl AsRef<Path>) -> Result<Value, String> {
    serde_json::from_slice(&fs::read(path).map_err(|e| e.to_string())?).map_err(|e| e.to_string())
}

fn register_gap_models(root: &Path) -> Result<(), String> {
    use crypto_lib::cryptanalysis::ecbench::{
        canonical::{canonical, sha256_hex},
        conductor_gap,
    };
    let catalog_text =
        fs::read_to_string("docs/curves/registry.json").map_err(|e| e.to_string())?;
    let mut catalog: Value = serde_json::from_str(&catalog_text).map_err(|e| e.to_string())?;
    let mut aliases = read("src/cryptanalysis/curve_aliases.json")?;
    let p = conductor_gap::params(263, 241)?;
    let mut additions = Vec::new();
    let mut new_rows = Vec::new();
    for (level, inst) in [
        ("crater", conductor_gap::crater(p)?),
        ("floor", conductor_gap::floor(p, 59297)?),
    ] {
        let id = inst.curve_id();
        let parts: Vec<_> = id.icv1.split(':').collect();
        let field = json!({"characteristic":inst.curve.p,"degree":1,"representation":"prime",
            "element_encoding":"hex integer modulo the characteristic"});
        let curve = json!({"model":"short Weierstrass","coefficients":[0,0,0,inst.curve.a,inst.curve.b],
            "subgroup_order":inst.r.to_string(),"cofactor":inst.cofactor,
            "generator":[format!("0x{:x}",inst.generator.0),format!("0x{:x}",inst.generator.1)],"target_group":"prime-order subgroup"});
        let hash = sha256_hex(canonical(&json!({"field":field,"curve":curve}))?.as_bytes());
        let ec1 = format!(
            "EC1P{}Cecbenchgaph{}",
            64 - inst.curve.p.leading_zeros(),
            &hash[..12]
        );
        let row = json!({"slug":id.slug,"icv1":id.icv1,"family":"prime","model_json":id.model_json,
            "params":{"p":inst.curve.p.to_string(),"a":inst.curve.a.to_string(),"b":inst.curve.b.to_string(),
                "generator_call":format!("conductor_gap::{level}(params(263,241), seed59297 where applicable)")},
            "trace":parts[2].parse::<i64>().map_err(|e|e.to_string())?,"order":parts[3],"j":parts[4],"end":parts[5],
            "aliases":[ec1],"standard_names":[],"sources":["docs/curves/sources/ecbench-ten-questions.json"],
            "representations":[{"ec1":ec1,"curve_uid":format!("urn:ec-record:1:sha256:{hash}"),
                "field_sha256":sha256_hex(canonical(&field)?.as_bytes()),"field":field,"curve":curve,
                "sources":["docs/curves/sources/ecbench-ten-questions.json"]}]});
        if let Some(old) = catalog["curves"]
            .as_array()
            .unwrap()
            .iter()
            .find(|r| r["slug"] == row["slug"])
        {
            if old["model_json"] != row["model_json"] {
                return Err("registry model disagreement".into());
            }
        } else {
            catalog["curves"].as_array_mut().unwrap().push(row.clone());
            new_rows.push(row.clone());
            for name in [&id.slug, &id.icv1, &ec1] {
                aliases["aliases"][name.to_lowercase()] = json!(id.slug);
            }
        }
        additions.push(row);
    }
    save(
        &root.join("registered-gap-models.json"),
        &json!({"curves":additions}),
    )?;
    save(
        Path::new("docs/curves/sources/ecbench-ten-questions.json"),
        &json!({"schema":"ecbench-native-registry-additions/v1",
        "generated_by":"examples/ecbench_leads.rs register","curves":additions}),
    )?;
    let pretty_one = |v: &Value| -> Result<Vec<u8>, String> {
        let mut bytes = Vec::new();
        let mut serializer = serde_json::Serializer::with_formatter(
            &mut bytes,
            serde_json::ser::PrettyFormatter::with_indent(b" "),
        );
        serde::Serialize::serialize(v, &mut serializer).map_err(|e| e.to_string())?;
        bytes.push(b'\n');
        Ok(bytes)
    };
    // Append rows in the existing one-space catalog layout; preserve every
    // historical row's raw bytes and key order, not just its parsed value.
    if !new_rows.is_empty() {
        let start = catalog_text
            .find("\"curves\":")
            .ok_or("catalog curves key")?
            + "\"curves\":".len();
        let mut stream = serde_json::Deserializer::from_str(&catalog_text[start..])
            .into_iter::<Box<serde_json::value::RawValue>>();
        let _array = stream
            .next()
            .ok_or("catalog array")?
            .map_err(|e| e.to_string())?;
        let end = start + stream.byte_offset();
        let insert = catalog_text[..end - 1].trim_end().len();
        let rows = new_rows
            .iter()
            .map(|r| {
                pretty_one(r)
                    .and_then(|bytes| String::from_utf8(bytes).map_err(|e| e.to_string()))
                    .map(|s| {
                        s.trim_end()
                            .lines()
                            .map(|line| format!("  {line}"))
                            .collect::<Vec<_>>()
                            .join("\n")
                    })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mut output = catalog_text;
        output.insert_str(insert, &format!(",\n{}", rows.join(",\n")));
        if serde_json::from_str::<Value>(&output).map_err(|e| e.to_string())? != catalog {
            return Err("registry append changed data".into());
        }
        fs::write("docs/curves/registry.json", output).map_err(|e| e.to_string())?;
    }
    let mut bytes = Vec::new();
    let mut serializer = serde_json::Serializer::with_formatter(
        &mut bytes,
        serde_json::ser::PrettyFormatter::with_indent(b""),
    );
    serde::Serialize::serialize(&aliases["aliases"], &mut serializer).map_err(|e| e.to_string())?;
    let text = format!("{{\n\"schema_version\": {},\n\"generated_from\": {},\n\"generated_by\": {},\n\"aliases\": {}\n}}\n",
        aliases["schema_version"],aliases["generated_from"],aliases["generated_by"],
        String::from_utf8(bytes).map_err(|e|e.to_string())?);
    if serde_json::from_str::<Value>(&text).map_err(|e| e.to_string())? != aliases {
        return Err("alias formatting changed metadata".into());
    }
    fs::write("src/cryptanalysis/curve_aliases.json", text).map_err(|e| e.to_string())?;
    Ok(())
}

fn publish(root: &Path) -> Result<(), String> {
    let summary = read(root.join("summary.json"))?;
    let shapes = read(root.join("shape-controls.json"))?;
    let failed = read(root.join("failed-replays-v2.json"))?;
    if shapes["ok"] != true || failed["ok"] != true {
        return Err("controls/replays not accepted".into());
    }
    let mut statuses = BTreeMap::<String, u64>::new();
    let mut targets = Vec::new();
    let mut renders = Vec::new();
    let mut shape_pairs = BTreeMap::<(String, String, String), BTreeMap<String, Value>>::new();
    let mut shape_checks = Vec::new();
    for s in summary["sessions"].as_array().ok_or("missing sessions")? {
        if s["audit_ok"] != true {
            return Err("failed session audit".into());
        }
        let name = s["name"].as_str().ok_or("missing session name")?;
        let session = read(root.join("sessions").join(name).join("session.json"))?;
        let text = fs::read_to_string(root.join("sessions").join(name).join("records.jsonl"))
            .map_err(|e| e.to_string())?;
        for line in text.lines() {
            let r: Value = serde_json::from_str(line).map_err(|e| e.to_string())?;
            if r["warmup"] != false {
                continue;
            }
            let status = r["outcome"]["status"].as_str().ok_or("missing status")?;
            *statuses.entry(status.into()).or_default() += 1;
            targets.push(json!({"session":name,"session_id":r["session_id"],"run_id":r["run_id"],"record_id":r["record_id"],
                "curve":r["workload"]["curve"],"workload_id":r["workload"]["workload_id"],"target":r["workload"]["target"],
                "arm":r["arm"],"status":status,"cost":r["cost"],"counters":r["counters"],"factor_base":r["factor_base"],
                "solver":r["solver"],"fully_priced_ratio":null,"online_wall_speedup":null}));
            if name == "binary-shape" && !r["arm"].as_str().unwrap().starts_with("rho") {
                let arm = r["arm"].as_str().unwrap();
                let (engine, version) = arm.rsplit_once('-').ok_or("shape arm")?;
                shape_pairs
                    .entry((
                        r["workload"]["workload_id"].to_string(),
                        r["round"].to_string(),
                        engine.into(),
                    ))
                    .or_default()
                    .insert(version.into(), r);
            }
        }
        let rows: Vec<_> = summary["rows"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|r| r["session"] == name)
            .cloned()
            .collect();
        renders.push(
            json!({"session":name,"ledger":s,"metadata":session,"rows":rows,
            "path":format!("reports/run-20261015T{:04}Z-{name}",
                session["started_unix_ms"].as_u64().unwrap_or(0)/60000%1440/60*100
                +session["started_unix_ms"].as_u64().unwrap_or(0)/60000%60)}),
        );
    }
    for ((workload, round, engine), versions) in shape_pairs {
        let (a, b) = (
            versions.get("v1").ok_or("missing v1")?,
            versions.get("v2").ok_or("missing v2")?,
        );
        let stable = a["cost"]["total_gae"] == b["cost"]["total_gae"]
            && a["outcome"] == b["outcome"]
            && a["counters"] == b["counters"];
        let n = b["workload"]["curve"]["field_degree"]
            .as_u64()
            .ok_or("missing n")?;
        let actual =
            b["solver"]["n_vars"] == json!(9 + n) && b["solver"]["n_equations"] == json!(2 * n);
        if !stable || !actual {
            return Err(format!("shape correction changed stable results or has incorrect metadata: {workload} {engine}"));
        }
        shape_checks.push(json!({"workload":workload,"round":round,"engine":engine,"stable_outcome_counts":stable,
            "legacy_variables":a["solver"]["n_vars"],"actual_variables":b["solver"]["n_vars"],
            "legacy_equations":a["solver"]["n_equations"],"actual_equations":b["solver"]["n_equations"],"actual_shape_verified":actual}));
    }
    let attempts: u64 = statuses.values().sum();
    let verified = statuses.get("verified").copied().unwrap_or(0);
    let results = json!({"schema":"ecbench-leads-results/v1","attempts":attempts,"verified":verified,"statuses":statuses,
        "all_audits_ok":true,"failed_replays":failed["replays"].as_array().unwrap().len(),"shape_pair_checks":shape_checks,
        "fully_priced_speedup":null,"asymptotic_claim":false,"original_questions_resolved":0});
    save(&root.join("results.json"), &results)?;
    save(&root.join("per-target.json"), &json!(targets))?;
    save(
        &root.join("render-data.json"),
        &json!({"results":results,"sessions":renders,"shapes":shapes,"rows":summary["rows"]}),
    )?;
    let mut report=format!("# Pursuing the ecbench leads\n\n{attempts} fresh public single-target attempts; {verified} externally verified recoveries. Status counts: `{statuses:?}`. Every session audit passed. The native supplement separately replayed {} unsuccessful outcomes. All ten original mathematical questions remain partial.\n\n[Frozen protocol](PROTOCOL-v1.md), [source-linked parent campaign](../ecbench_open_questions_20261014/REPORT.md), [summary](summary.json), [all per-target rows](per-target.json), [shape construction controls](shape-controls.json), [results and paired metadata checks](results.json), [failure replays](failed-replays-v2.json), [exact command manifest](execution-manifest.json).\n\nPublic target seed20261015; constructor seed59297 for the unchanged prime models. One target per workload, cold setup per execution, no shared cross-target state, L0. Wall times are inadmissible. S is counted GAE/sqrt(r); unpriced native work remains explicit, and count quotients do not bound full-cost ratios. Ranges and means below use verified runs only, never failures as zero. Ratios are null whenever a curve/arm is censored. The degree83 panel has no matched rho cost and emits no ratio.\n\n## Observed lead tests\n\nThe versioned PDP3 adapter reports the actual chained coordinates, including intermediate field variables and all equation blocks. Eight independently constructed presentations and {} paired engine/target checks establish this metadata correction on the frozen domains. They do not establish solving-degree growth.\n\nThe new `ic.orbit_rank_v1` collector uses unknown base logs and the same relation/rank code at both integer widths. Its orbit table stores finite enumerated relative-phase descriptors and lifts hits to full point equations before modular elimination. It is not a low-degree quotient correspondence. The exact E0 degree83 field/subgroup is retained. A bounded completed execution is distinct from a verified degree83 recovery.\n\nThe prime panel tests the setup/yield balance on eight fresh single-target workloads per model. Base selection from this panel is descriptive; it is not a held-out speedup or exponent result.\n\n",failed["replays"].as_array().unwrap().len(),shape_checks.len());
    report.push_str("| Curve | Arm | Verified / attempted | Mean counted S | Count quotient / rho | Status |\n|---|---|---:|---:|---:|---|\n");
    for r in summary["rows"].as_array().unwrap() {
        let number = |v: &Value| v.as_f64().map(|x| format!("{x:.6}")).unwrap_or("—".into());
        writeln!(
            report,
            "| {} | {} | {} / {} | {} | {} | {} |",
            r["curve"].as_str().unwrap(),
            r["arm"].as_str().unwrap(),
            r["verified"],
            r["records"],
            number(&r["mean_s"]),
            number(&r["paired_complete_counted_ratio_candidate_over_rho"]),
            r["statuses"]
        )
        .unwrap();
    }
    report.push_str("\n## Original question status\n\n| Q | Follow-up evidence | Original obligation still open |\n|---|---|---|\n|1|Fresh-target prime base-size panel|Scalable cheap decomposition and full-cost complexity|\n|2|Versioned system shape, constructed systems and paired old/v2 controls|Solving-degree bounds at growing dimensions|\n|3|Direct/quotient point-table controls and width-generic exact lifting|Compact algebraic correspondence without enumerated phase search|\n|4|Parent nonlinear-base controls retained|Scalable nonlinear membership and favorable full cost|\n|5|Parent exact quasi-subfield census retained|Useful general family and a complete IC pipeline|\n|6|Actual chained variable/equation counts|Higher-arity lifting/solving advantage|\n|7|Modular rank and target pinning at both widths|General useful-rank law and reliable paid descent|\n|8|Parent same-order curve diagnostics retained|Certified paid isogeny/Hecke transport|\n|9|Parent degree83 pullbacks and small generic recovery retained|General norm and cheaper elliptic-component solver|\n|10|Fresh paired ladder and exact degree83 bounded IC executions|Successful N83 recovery, full costs, independent/L2/hardware evidence and growing-size crossover|\n\nSources and paper-specific qualifications remain in the [29-reference literature record](../ecbench_open_questions_20261014/literature/REPORT.md). In particular, first-fall degree does not establish solving degree; enumerated quotient descriptors do not provide a compact algebraic system; a subgroup-preserving cover does not price a destination solver. No new mathematical theorem is claimed.\n\n## Run evidence\n\n");
    for s in &renders {
        let name = s["session"].as_str().unwrap();
        let path = s["path"].as_str().unwrap();
        writeln!(report,"- [{name}]({path}/REPORT.md): [sealed session](sessions/{name}/session.json), [audit]({name}.audit.json), [resources]({name}.resources.json).").unwrap();
        fs::create_dir_all(root.join(path)).map_err(|e| e.to_string())?;
        save(&root.join(path).join("data.json"), s)?;
        let prose=format!("# {name}: ecbench follow-up\n\nSession `{}`; source `{}`; binary `{}`; host `{}`.\n\n{} attempts, {} verified recoveries; {} exact verified replays. [Spec](../../specs/{name}.json), [sealed records](../../sessions/{name}/records.jsonl), [audit](../../{name}.audit.json), [resource vector](../../{name}.resources.json), [supplemental failure replay](../../failed-replays-v2.json).\n\nFresh public targets, seed20261015, cold setup repeated, no planted log. L0; no admitted wall time. Figures are counted diagnostics. Successful-only means/ranges do not represent censored solve costs; incomplete count quotients do not bound full cost. The wide panel has no matched rho run.\n\n![Recorded counted costs](cost.svg)\n",s["ledger"]["session_id"],s["metadata"]["git_commit"],s["metadata"]["binary_sha256"],s["metadata"]["env_class_id"],s["ledger"]["records"],s["ledger"]["verified_records"],s["ledger"]["replays"]);
        fs::write(root.join(path).join("REPORT.md"), prose).map_err(|e| e.to_string())?;
    }
    report.push_str("\n## Remaining gates and publication\n\nIndependent-host/L2/GPU access remains unavailable. The PR is unmerged and subject to CI/review. Frozen parent sessions were not changed. Canonical diagnostics are refreshed from sealed data; theoretical barrier curves, admitted speedup graphs and cover/isogeny diagrams retain their prior values because this iteration proves no new barrier crossing or transport map. See the [resolution problem](../../problems/ecbench-ten-question-resolution-gates.md).\n");
    fs::write(root.join("REPORT.md"), report).map_err(|e| e.to_string())?;
    Ok(())
}

fn main() -> Result<(), String> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let [cmd, root] = args.as_slice() else {
        return Err("usage: ecbench_leads freeze|shapes|analyze STUDY".into());
    };
    match cmd.as_str() {
        "freeze" => freeze(Path::new(root)),
        "shapes" => shapes(Path::new(root)),
        "analyze" => analyze::run(Path::new(root)),
        "publish" => publish(Path::new(root)),
        "register" => register_gap_models(Path::new(root)),
        _ => Err("unknown command".into()),
    }
}

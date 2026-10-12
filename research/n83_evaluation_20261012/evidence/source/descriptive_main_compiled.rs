//! Post-panel descriptive export. This is separate from the frozen primary estimator.
//! No solving, sampling, bootstrap, new confidence intervals, or timing occurs here.
use serde_json::{json,Value};
use std::collections::{BTreeMap,BTreeSet};
use std::fs::{self,File,OpenOptions};
use std::io::{BufRead,BufReader,Write};
use std::os::unix::fs::PermissionsExt;
use std::path::{Path,PathBuf};
use std::process::{Command,Stdio};
use std::time::{SystemTime,UNIX_EPOCH};

type Result<T>=std::result::Result<T,Box<dyn std::error::Error>>;
fn error(s:impl Into<String>)->Box<dyn std::error::Error>{s.into().into()}
fn hash(path:&Path)->Result<String>{
    let out=Command::new("/usr/bin/sha256sum").arg("--").arg(path).output()?;
    if !out.status.success(){return Err(error(format!("hash failed: {}",path.display())));}
    Ok(String::from_utf8(out.stdout)?.split_whitespace().next().ok_or_else(||error("empty hash"))?.into())
}
fn hash_bytes(bytes:&[u8])->Result<String>{
    let mut child=Command::new("/usr/bin/sha256sum").stdin(Stdio::piped()).stdout(Stdio::piped()).spawn()?;
    child.stdin.take().ok_or_else(||error("missing hash stdin"))?.write_all(bytes)?;
    let out=child.wait_with_output()?;
    if !out.status.success(){return Err(error("compiled-source hash failed"));}
    Ok(String::from_utf8(out.stdout)?.split_whitespace().next().ok_or_else(||error("empty source hash"))?.into())
}
fn json_file(path:&Path)->Result<Value>{Ok(serde_json::from_reader(File::open(path)?)?)}
fn jsonl(path:&Path)->Result<Vec<Value>>{
    let mut out=Vec::new();for (i,line) in BufReader::new(File::open(path)?).lines().enumerate(){
        let line=line?;if line.trim().is_empty(){continue;}
        out.push(serde_json::from_str(&line).map_err(|e|error(format!("{}:{}: {e}",path.display(),i+1)))?);
    }Ok(out)
}
fn nonnegative(v:&Value)->Option<f64>{v.as_f64().filter(|x|x.is_finite()&&*x>=0.)}
fn successful(r:&Value)->bool{r["status"]=="COMPLETED_VERIFIED"&&r["semantic"]["valid"]==true&&r["semantic"]["online_boundary_qualified"]==true}
fn positive(v:&Value)->Option<f64>{v.as_f64().filter(|x|x.is_finite()&&*x>0.)}
fn reciprocal_presentation(matched:&Value,mean_key:&str,ci_key:&str)->Value{
    let n=matched["n"].as_u64().unwrap_or(0);
    let ratio=positive(&matched[mean_key]);
    let inverse=ratio.map(|x|1./x).filter(|x|x.is_finite()&&*x>0.);
    let inverse_ci=if n>=2{
        matched[ci_key].as_array().filter(|a|a.len()==2).and_then(|a|{
            let lower=positive(&a[0])?;let upper=positive(&a[1])?;
            if lower>upper{return None;}
            let reversed=[1./upper,1./lower];
            if reversed.iter().all(|x|x.is_finite()&&*x>0.){Some(reversed)}else{None}
        })
    }else{None};
    json!({"presentation_class":"ALGEBRAIC_RECIPROCAL_OF_FROZEN_PRIMARY_ESTIMATOR","new_estimator":false,"new_resampling":false,
        "n":n,"primary_mean_key":mean_key,"primary_ci_key":ci_key,"strong_over_ic_primary_ratio":matched[mean_key],
        "ic_over_strong_primary_ratio":inverse,"strong_over_ic_primary_95_ci":matched[ci_key],"ic_over_strong_primary_95_ci":inverse_ci,
        "ratio_identity":"IC/strong = 1 / (strong/IC)","interval_identity":"IC/strong interval = [1 / original_upper, 1 / original_lower]",
        "confidence_interval_available":inverse_ci.is_some(),"no_confidence_interval_fabricated_for_n_below_two":true,
        "missing_nonpositive_nonfinite_or_invalid_values_remain_null":true})
}
fn describe(rows:&[Value],key:&str,successful_only:bool)->Value{
    let selected:Vec<_>=rows.iter().filter(|r|!successful_only||r["successful_verified"]==true).collect();
    let values:Vec<_>=selected.iter().filter_map(|r|nonnegative(&r[key]).map(|v|(r["pair_id"].clone(),r["status"].clone(),v))).collect();
    // Charged timeout/failure wall, CPU and RSS are observations, not solved-DLP costs.
    let n=values.len();let statistics=if n==0{json!({"minimum":null,"maximum":null,"arithmetic_mean":null})}else{
        let mut sum=0.;let mut compensation=0.;let mut min=f64::INFINITY;let mut max=f64::NEG_INFINITY;
        for (_,_,v) in &values {let y=*v-compensation;let t=sum+y;compensation=(t-sum)-y;sum=t;min=min.min(*v);max=max.max(*v);}
        json!({"minimum":min,"maximum":max,"arithmetic_mean":sum/n as f64})
    };
    json!({"metric":key,"population":if successful_only{"successful_verified_formal_runs_only"}else{"all_formal_attempts_with_observed_charged_values_including_failures_and_caps"},
        "planned_formal_rows":rows.len(),"selected_rows":selected.len(),"observed_values":n,"missing_values":selected.len()-n,
        "statistics":statistics,"raw_values":values.iter().map(|(id,status,v)|json!({"pair_id":id,"status":status,"value":v})).collect::<Vec<_>>(),
        "zero_imputation":false,"confidence_interval":null,"secondary_descriptive_only":true})
}
fn raw_file(r:&Value,path_key:&str,hash_key:&str)->Result<Value>{
    let Some(s)=r[path_key].as_str()else{return Ok(json!({"path":null,"recorded_sha256":r[hash_key],"observed_sha256":null,"available":false,"hash_match":null}));};
    let p=Path::new(s);if !p.is_file(){return Ok(json!({"path":s,"recorded_sha256":r[hash_key],"observed_sha256":null,"available":false,"hash_match":null}));}
    let observed=hash(p)?;let matches=r[hash_key].as_str()==Some(observed.as_str());
    Ok(json!({"path":s,"recorded_sha256":r[hash_key],"observed_sha256":observed,"available":true,"hash_match":matches}))
}
fn method_row(pair:&Value,arm:&str,path:&Path)->Result<Value>{
    if !path.is_file(){return Ok(json!({"pair_id":pair["pair_id"],"arm":arm,"status":"NOT_ATTEMPTED_OR_RECEIPT_UNAVAILABLE","successful_verified":false,"receipt_path":path,
        "receipt_sha256":null,"online_ms":null,"preparation_ms":null,"whole_process_wall_s":null,"cpu_user_s":null,"cpu_system_s":null,"cpu_total_s":null,"peak_rss_bytes":null,
        "seed":{"walk_seed":null,"rank_seed":null},"raw_receipt":null,"missing_receipt":true}));}
    let r=json_file(path)?;
    if r["schema"]!="n83.process-receipt/v1"{return Err(error(format!("unexpected receipt schema: {}",path.display())));}
    if let Some(recorded_arm)=r["arm"].as_str(){
        if recorded_arm!=arm{return Err(error(format!("receipt arm mismatch: {}",path.display())));}
    }else if r["status"]!="LAUNCH_FAILED"{return Err(error(format!("receipt arm missing: {}",path.display())));}
    let row=&r["semantic"]["selected_row"];
    let walk_seed=if arm=="rho" {row.get("walk_seed").cloned().unwrap_or_else(||r["argv"][2].clone())}else{Value::Null};
    let rank_seed=if arm=="ic" {r["argv"][3].clone()}else{Value::Null};
    Ok(json!({"pair_id":pair["pair_id"],"arm":arm,"status":r["status"],"successful_verified":successful(&r),"receipt_path":path,"receipt_sha256":hash(path)?,
        "online_ms":r["semantic"]["producer_online_ms"],"preparation_ms":r["semantic"]["producer_preparation_ms"],
        "whole_process_wall_s":r["elapsed_process_wall_s"],"cpu_user_s":r["cpu_user_s"],"cpu_system_s":r["cpu_system_s"],"cpu_total_s":r["cpu_total_s"],"peak_rss_bytes":r["peak_rss_bytes"],
        "seed":{"walk_seed":walk_seed,"rank_seed":rank_seed},"timeout":r["timeout"],"return_code":r["return_code"],"termination_signal":r["termination_signal"],"online_incomplete_or_censored":!successful(&r),
        "producer_phase_sum_ms":r["semantic"]["producer_phase_sum_ms"],"producer_phase_gap_ms":r["semantic"]["producer_phase_gap_ms"],
        "five_phase_exact_conformance":r["semantic"]["five_phase_exact_conformance"],"resource_envelope":r["resource_envelope"],"target_sha256":r["target_sha256"],
        "raw_file_integrity":{"stdout":raw_file(&r,"stdout_path","stdout_sha256")?,"stderr":raw_file(&r,"stderr_path","stderr_sha256")?,"semantic":raw_file(&r,"semantic_path","semantic_sha256")?},
        "raw_receipt":r,"missing_receipt":false}))
}
fn run()->Result<()>{
    let args:Vec<_>=std::env::args().skip(1).collect();if args.len()%2!=0{return Err(error("all arguments need values"));}
    let mut opts=BTreeMap::new();for p in args.chunks(2){if !["--pairs","--primary-statistics","--out"].contains(&p[0].as_str()){return Err(error(format!("unknown option {}",p[0])));}if opts.insert(p[0].as_str(),p[1].as_str()).is_some(){return Err(error("duplicate option"));}}
    let required=|key:&str|opts.get(key).map(|s|PathBuf::from(*s)).ok_or_else(||error(format!("missing {key}")));
    let pairs_path=required("--pairs")?;let primary_path=required("--primary-statistics")?;let out=required("--out")?;
    let pairs=jsonl(&pairs_path)?;let primary=json_file(&primary_path)?;let pairs_hash=hash(&pairs_path)?;
    if primary["schema"]!="n83.paired-statistics/v1"{return Err(error("frozen primary statistics schema mismatch"));}
    if primary["pairs_sha256"].as_str()!=Some(pairs_hash.as_str()){return Err(error("primary analysis pairs hash does not match supplied pairs"));}
    let formal:Vec<_>=pairs.iter().filter(|p|p["kind"]=="matched").collect();
    if formal.len()!=5{return Err(error(format!("expected frozen five formal matched pairs, found {}",formal.len())));}
    let mut ids=BTreeSet::new();for p in &formal{let id=p["pair_id"].as_str().ok_or_else(||error("pair_id absent"))?;if !ids.insert(id){return Err(error("duplicate formal pair_id"));}}
    let mut method_rows:BTreeMap<String,Vec<Value>>=BTreeMap::new();
    for p in formal{
        for (arm,key) in [("rho","baseline_receipt"),("ic","candidate_receipt")]{
            let path=PathBuf::from(p[key].as_str().ok_or_else(||error(format!("{key} missing")))?);
            method_rows.entry(arm.into()).or_default().push(method_row(p,arm,&path)?);
        }
    }
    let mut methods=json!({});
    for (arm,rows) in &method_rows{
        let mut counts:BTreeMap<String,usize>=BTreeMap::new();for r in rows{*counts.entry(r["status"].as_str().unwrap_or("STATUS_UNAVAILABLE").into()).or_default()+=1;}
        let mut solved=json!({});let mut attempts=json!({});
        for metric in ["online_ms","preparation_ms","whole_process_wall_s","cpu_user_s","cpu_system_s","cpu_total_s","peak_rss_bytes"] {solved[metric]=describe(rows,metric,true);}
        for metric in ["whole_process_wall_s","cpu_user_s","cpu_system_s","cpu_total_s","peak_rss_bytes"] {attempts[metric]=describe(rows,metric,false);}
        methods[arm]=json!({"planned_formal_runs":rows.len(),"status_counts":counts,"successful_verified_count":rows.iter().filter(|r|r["successful_verified"]==true).count(),
            "observed_receipt_count":rows.iter().filter(|r|r["missing_receipt"]==false).count(),"successful_run_descriptions":solved,"all_attempt_charge_descriptions":attempts,"rows":rows});
    }
    let protocol_path=pairs_path.parent().ok_or_else(||error("pairs file has no parent"))?.join("protocol.json");
    let protocol=if protocol_path.is_file(){Some(json_file(&protocol_path)?)}else{None};
    let executable=std::env::current_exe()?;
    let reciprocal=json!({"scope":"Algebraic presentation of the unchanged frozen matched-pair primary outputs; no new estimator, bootstrap or confidence interval.",
        "online":reciprocal_presentation(&primary["matched"],"geometric_mean_paired_online_ratio","paired_percentile_bootstrap_95_ci"),
        "whole_process":reciprocal_presentation(&primary["matched"],"geometric_mean_paired_whole_process_ratio","whole_process_bootstrap_95_ci")});
    let summary=json!({"schema":"n83.secondary-descriptive-report-data/v1","created_epoch_ms":SystemTime::now().duration_since(UNIX_EPOCH)?.as_millis(),
        "analysis_class":"SECONDARY_DESCRIPTIVE_NOT_PREREGISTERED_PRIMARY_ESTIMATOR","primary_estimator_unchanged":true,
        "statistical_scope":"Five timing repeats of one public synthetic point on one host. Descriptive min/max/arithmetic means only; no new CI, population inference, or primary-estimator replacement.",
        "accounting":"Successful online and preparation observations are separate from whole-process wall/CPU/RSS. Failed and capped process costs remain observed attempt charges; missing quantities stay null, never zero.",
        "claim_boundary":"Diagnostic export only; it does not establish admitted vs_rho, full ecbench/common-operation-equivalent cost, independent-host replay, or deployed-curve conclusions.",
        "units":{"online_ms":"milliseconds","preparation_ms":"milliseconds","whole_process_wall_s":"seconds","cpu_user_s":"seconds","cpu_system_s":"seconds","cpu_total_s":"seconds","peak_rss_bytes":"bytes"},
        "source_provenance":{"language":"Rust","compiled_source_main_rs_sha256":hash_bytes(include_bytes!("main.rs"))?,"compiled_manifest_sha256":hash_bytes(include_bytes!("../Cargo.toml"))?,"executable_path":executable,"executable_sha256":hash(&executable)?,"compiled_source_revision":option_env!("N83_SOURCE_COMMIT"),"hash_tool":"/usr/bin/sha256sum (native, post-panel)"},
        "inputs":{"pairs_path":pairs_path,"pairs_sha256":pairs_hash,"primary_statistics_path":primary_path,"primary_statistics_sha256":hash(&primary_path)?,"protocol_path":protocol_path,"protocol_sha256":if protocol_path.is_file(){Some(hash(&protocol_path)?)}else{None}},
        "primary_reciprocal_presentation":reciprocal,"frozen_primary_statistics_raw":primary,"frozen_protocol_raw":protocol,"pairs_raw":pairs,"methods":methods});
    let mut f=OpenOptions::new().write(true).create_new(true).open(&out)?;serde_json::to_writer_pretty(&mut f,&summary)?;f.write_all(b"\n")?;f.sync_all()?;
    fs::set_permissions(&out,fs::Permissions::from_mode(0o444))?;
    println!("{}",json!({"report_data_path":out,"report_data_sha256":hash(&out)?,"analysis_class":"SECONDARY_DESCRIPTIVE","formal_pairs":5,"primary_estimator_unchanged":true}));Ok(())
}
fn main(){if let Err(e)=run(){eprintln!("descriptive export error: {e}\nusage: n83-descriptive-summary --pairs pairs.jsonl --primary-statistics statistics.json --out NEW_REPORT_DATA_JSON");std::process::exit(2);}}

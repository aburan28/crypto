//! Post-hoc analysis supporting the frozen protocol's wall-censored attempts.
#[path = "analysis_censored.rs"]
mod analysis;
use crypto_lib::cryptanalysis::ecbench::canonical::sha256_hex;
use rand::{rngs::StdRng, SeedableRng};
use serde_json::{json, Value};
use std::{fs, path::Path};
const ROOT: &str = "research/conductor_gap_followups_20261016";
fn read(p: impl AsRef<Path>) -> Value {
    serde_json::from_slice(&fs::read(p).unwrap()).unwrap()
}
fn write(p: impl AsRef<Path>, v: &Value) {
    let p = p.as_ref();
    let tmp = p.with_extension(format!("tmp-{}", std::process::id()));
    fs::write(&tmp, serde_json::to_string_pretty(v).unwrap() + "\n").unwrap();
    fs::File::open(&tmp).unwrap().sync_all().unwrap();
    fs::rename(tmp, p).unwrap();
}
fn u(v: &Value, k: &str) -> u64 {
    v[k].as_u64().unwrap()
}
fn compare_prefix(a: &Value, b: &Value, ae: &[Value], be: &[Value]) -> Value {
    let qa = u(a, "queries");
    let qb = u(b, "queries");
    let common = qa.min(qb);
    let ap: Vec<_> = ae.iter().filter(|e| u(e, "query") < common).collect();
    let bp: Vec<_> = be.iter().filter(|e| u(e, "query") < common).collect();
    assert_eq!(ap, bp, "same-crater common query prefix differs");
    let full = qa == qb;
    if full {
        assert_eq!(ae, be);
        assert_eq!(a["verified"], b["verified"]);
    } else {
        assert!(
            a["native"]["status"] == "TIME_LIMIT" || b["native"]["status"] == "TIME_LIMIT",
            "unequal query counts require a wall censor"
        );
    }
    if a["verified"] == true && b["verified"] == true {
        assert_eq!(
            qa, qb,
            "two fully recovered same-stream runs must stop at the same query"
        );
        assert_eq!(a["native"]["recovered"], b["native"]["recovered"]);
    }
    json!({"common_queries":common,"common_events":ap.len(),"complete_streams_equal":full,"wall_censor_present":a["native"]["status"]=="TIME_LIMIT"||b["native"]["status"]=="TIME_LIMIT"})
}
fn main() {
    match std::env::args().nth(1).as_deref() {
        Some("progress") => analysis::progress(),
        Some("analyze") => analysis::analyze(),
        _ => panic!("progress | analyze"),
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    fn row(q: u64, status: &str) -> Value {
        json!({"queries":q,"verified":false,"native":{"status":status}})
    }
    #[test]
    fn censored_prefix_allows_later_events_but_rejects_a_changed_common_event() {
        let a = row(10, "TIME_LIMIT");
        let b = row(20, "QUERY_LIMIT");
        let e = json!({"query":5,"alpha":100,"rank":1});
        let later = json!({"query":12,"alpha":107,"rank":2});
        let c = compare_prefix(&a, &b, &[e.clone()], &[e.clone(), later.clone()]);
        assert_eq!(c["common_events"], 1);
        assert_eq!(c["complete_streams_equal"], false);
        let corrupt = json!({"query":5,"alpha":101,"rank":1});
        assert!(std::panic::catch_unwind(|| compare_prefix(
            &a,
            &b,
            &[e.clone()],
            &[corrupt.clone(), later.clone()]
        ))
        .is_err());
    }
    #[test]
    fn unequal_uncensored_counts_are_rejected() {
        assert!(std::panic::catch_unwind(|| compare_prefix(
            &row(10, "QUERY_LIMIT"),
            &row(20, "QUERY_LIMIT"),
            &[],
            &[]
        ))
        .is_err());
    }
    #[test]
    fn recovered_streams_must_match_completely() {
        let mut a = row(10, "PASS");
        a["verified"] = json!(true);
        a["native"]["recovered"] = json!(123);
        let b = a.clone();
        assert_eq!(
            compare_prefix(&a, &b, &[], &[])["complete_streams_equal"],
            true
        );
        let mut c = a.clone();
        c["queries"] = json!(20);
        c["native"]["status"] = json!("TIME_LIMIT");
        assert!(std::panic::catch_unwind(|| compare_prefix(&a, &c, &[], &[])).is_err());
    }
}

use rusqlite::{Connection, OpenFlags};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{env, error::Error, fs, path::Path};

type Result<T> = std::result::Result<T, Box<dyn Error>>;

const DB: &str = "../weak_curve_database_20261012/weak_curves.sqlite";
const CASES: &str = "../weak_curve_case_index_20261012/published_cases.sqlite";
const CURVE_REGISTRY: &str = "../../docs/curves/registry.json";
const TRANSFER_REGISTRY: &str = "../../docs/curves/jacobian-transfers.json";

fn file_hash(path: &str) -> Result<String> {
    let bytes = fs::read(path)?;
    Ok(format!("{:x}", Sha256::digest(bytes)))
}

fn connect(path: &str) -> Result<Connection> {
    Ok(Connection::open_with_flags(
        path,
        OpenFlags::SQLITE_OPEN_READ_ONLY,
    )?)
}

fn models(db: &Connection) -> Result<Vec<Value>> {
    let mut stmt = db.prepare(
        r#"SELECT m.icv1,m.slug,m.characteristic,m.extension_degree,m.field_bits,m.trace,m.curve_order,
         m.cubic_model_label,m.quadratic_model_label,m.subgroup_advantage,m.inventory_law,m.model_sha256,
         json_extract(p.family_checks_json,'$."WCF1:g3-cubic".status'),
         json_extract(p.family_checks_json,'$."WCF1:g3-quadratic".status')
         FROM curve_model m JOIN curve_property p USING(icv1)
         ORDER BY m.field_bits,m.extension_degree,m.icv1"#,
    )?;
    let rows = stmt.query_map([], |row| {
        let cubic: String = row.get(7)?;
        let quadratic: String = row.get(8)?;
        let signal = if cubic == "yes" || quadratic == "yes" {
            "exact_covering_model"
        } else if cubic == "outside_domain" {
            "outside_covering_index"
        } else if cubic == "unclassified_in_this_index" {
            "unclassified_covering_model"
        } else {
            "no_indexed_covering_model"
        };
        Ok(json!({
            "id": row.get::<_,String>(0)?, "slug": row.get::<_,String>(1)?,
            "p": row.get::<_,String>(2)?, "degree": row.get::<_,i64>(3)?,
            "bits": row.get::<_,i64>(4)?, "trace": row.get::<_,String>(5)?,
            "order": row.get::<_,String>(6)?, "cubic": cubic,
            "quadratic": quadratic, "subgroup_advantage": row.get::<_,String>(9)?,
            "inventory_law": row.get::<_,String>(10)?, "sha256": row.get::<_,String>(11)?,
            "cubic_screen": row.get::<_,Option<String>>(12)?,
            "quadratic_screen": row.get::<_,Option<String>>(13)?,
            "signal": signal, "attack_ready": false,
            "source": "local_corpus"
        }))
    })?;
    rows.collect::<rusqlite::Result<Vec<_>>>()
        .map_err(Into::into)
}

fn cases(db: &Connection) -> Result<Vec<Value>> {
    let mut stmt = db.prepare(
        "SELECT case_no,d,n,e,cover_genus,coefficient,printed_page \
         FROM published_case ORDER BY case_no",
    )?;
    let rows = stmt.query_map([], |row| {
        Ok(json!({
            "case_no": row.get::<_,i64>(0)?, "degree": row.get::<_,i64>(1)?,
            "rank": row.get::<_,i64>(2)?, "excess_genus": row.get::<_,i64>(3)?,
            "cover_genus": row.get::<_,i64>(4)?, "coefficient": row.get::<_,String>(5)?,
            "printed_page": row.get::<_,i64>(6)?,
            "source_url": "https://eprint.iacr.org/2015/805.pdf"
        }))
    })?;
    rows.collect::<rusqlite::Result<Vec<_>>>()
        .map_err(Into::into)
}

fn finding(
    subject_id: &str,
    code: &str,
    criterion: &str,
    scope: &str,
    selector: &str,
    grade: &str,
    statement: &str,
) -> Value {
    // Canonical compact JSON is sorted by serde_json's map implementation.
    // Evidence grades, source URLs and measurements are deliberately outside
    // the identity so a finding keeps its ID as its verification improves.
    let key = json!({"v":1,"subject_id":subject_id,"mechanism":code,"criterion":criterion,
        "scope":scope,"selector":selector});
    let digest = format!("{:x}", Sha256::digest(serde_json::to_vec(&key).unwrap()));
    json!({"id":format!("WID1:{code}:{}", &digest[..12]),"criterion_id":criterion,
        "mechanism":code,"scope":scope,"selector":selector,"grade":grade,
        "statement":statement})
}

fn attach_findings(records: &mut [Value]) {
    for record in records {
        let id = record["id"].as_str().unwrap();
        let mut findings = Vec::new();
        match record["signal"].as_str().unwrap() {
            "anomalous_prime_field" => findings.push(finding(
                id,
                "ANOM",
                "anomalous_prime",
                "selected_subgroup",
                "r=101",
                "exact_condition",
                "Trace one: #E(F_101)=101 for the selected order-101 point.",
            )),
            "smooth_selected_group" => findings.push(finding(
                id,
                "PH",
                "smooth_selected_order",
                "selected_subgroup",
                "r=10",
                "exact_condition",
                "The selected cyclic group has order 10=2·5.",
            )),
            "low_embedding_degree_candidate" => findings.push(finding(
                id,
                "MOV",
                "mov_weil_pairing",
                "selected_subgroup",
                "r=13;k=2",
                "candidate",
                "The selected order has embedding degree two; pairing transfer remains open.",
            )),
            "exact_covering_model" => {
                if record["cubic"] == "yes" {
                    findings.push(finding(
                        id,
                        "GHS3C",
                        "cubic_norm_one",
                        "model",
                        "genus3-cubic",
                        "covering_model",
                        "The exact model has the recorded cubic genus-three covering label.",
                    ));
                }
                if record["quadratic"] == "yes" {
                    findings.push(finding(
                        id,
                        "GHS3Q",
                        "quadratic_torus",
                        "model",
                        "genus3-quadratic",
                        "covering_model",
                        "The exact model has the recorded quadratic genus-three covering label.",
                    ));
                }
            }
            _ => {}
        }
        record["findings"] = json!(findings);
    }
}

fn named_curves() -> Result<Vec<Value>> {
    let registry: Value = serde_json::from_slice(&fs::read(CURVE_REGISTRY)?)?;
    let transfers: Value = serde_json::from_slice(&fs::read(TRANSFER_REGISTRY)?)?;
    let registry_curves = registry["curves"]
        .as_array()
        .ok_or("curve registry has no curves")?;
    let transfer_curves = transfers["curves"]
        .as_array()
        .ok_or("transfer registry has no curves")?;
    let mut records = Vec::new();
    for c in registry_curves {
        let names = c["standard_names"]
            .as_array()
            .ok_or("missing standard_names")?;
        if names.is_empty() {
            continue;
        }
        let id = c["icv1"].as_str().ok_or("missing ICV1")?;
        let slug = c["slug"].as_str().ok_or("missing slug")?;
        let degree = if c["family"] == "prime" {
            1
        } else {
            c["params"]["m"]
                .as_i64()
                .or_else(|| c["params"]["n"].as_i64())
                .ok_or("missing binary degree")?
        };
        let p = if degree == 1 {
            c["params"]["p"].as_str().ok_or("missing prime modulus")?
        } else {
            "2"
        };
        let bits = if degree == 1 {
            slug.split('-')
                .nth(1)
                .ok_or("missing field tag")?
                .strip_prefix("fp")
                .ok_or("bad prime field tag")?
                .parse::<i64>()?
        } else {
            degree
        };
        let trace = id.split(':').nth(2).ok_or("missing ICV1 trace")?;
        let name = names[0].as_str().ok_or("invalid standard name")?;
        let source_url = c["standards_provenance"]
            .as_array()
            .and_then(|provenance| {
                provenance
                    .iter()
                    .flat_map(|p| p["sources"].as_array().into_iter().flatten())
                    .find_map(|s| s["url"].as_str())
            })
            .unwrap_or("https://github.com/aburan28/crypto/blob/main/docs/curves/registry.json");
        let mut record = json!({"id":id,"slug":slug,"name":name,"aliases":names,
            "p":p,"degree":degree,"bits":bits,"trace":trace,"order":c["order"],
            "field_modulus":c["params"]["modulus"],"model_a":c["params"]["a"],
            "model_b":if c["family"] == "koblitz" {json!("0x1")} else {c["params"]["b"].clone()},
            "signal":"named_curve_unassessed","attack_ready":false,
            "source":"named_curve_registry","source_url":source_url,
            "subgroup_advantage":"not_evaluated_in_this_snapshot","findings":[]});
        if name == "Oakley Group 3" || name == "Oakley Group 4" {
            let transfer = transfer_curves
                .iter()
                .find(|t| t["icv1"] == id)
                .ok_or("Oakley transfer record absent")?;
            let best = transfer["best_route"]
                .as_u64()
                .ok_or("Oakley best route absent")? as usize;
            let route = &transfer["routes"][best];
            if route["kind"] != "ghs_weil_descent" || route["evidence"] != "exact_invariant" {
                return Err("Oakley best route is not an exact GHS invariant".into());
            }
            let (
                expected_degree,
                expected_b,
                expected_order,
                r,
                cofactor,
                generator_multiplier,
                generator,
            ) = if name == "Oakley Group 3" {
                (
                    155,
                    "0x7338f",
                    "45671926166590716193865565914344635196769237316",
                    "3805993847215893016155463826195386266397436443",
                    12,
                    4,
                    [123, 456],
                )
            } else {
                (
                    185,
                    "0x1ee9",
                    "49039857307708443467467104857652682248052385001045053116",
                    "12259964326927110866866776214413170562013096250261263279",
                    4,
                    2,
                    [24, 13],
                )
            };
            if degree != expected_degree
                || c["params"]["b"] != expected_b
                || c["order"] != expected_order
                || transfer["subgroup"]["order"] != r
                || route["invariants"]["relative_degree"] != 5
                || route["invariants"]["magic_number"] != 5
                || route["destination"]["genus"] != "15"
                || route["cost"]["destination"]["enge_gaudry_hypothesis_holds"] != false
            {
                return Err(
                    "Oakley source or route drifted from checked RFC and transfer values".into(),
                );
            }
            let base_degree = route["invariants"]["base_degree"].as_i64().unwrap();
            let selector = format!("r={r};tower=5x{base_degree};g=15");
            let finding = finding(id,"GHS","ghs_binary","selected_subgroup",&selector,
                "exact_invariant_candidate",
                "The degree-five binary tower has GHS magic number five and predicted genus 15; a cover and subgroup transport are still needed.");
            record["signal"] = json!("binary_ghs_invariant");
            record["r"] = json!(r);
            record["full_order_cofactor"] = json!(cofactor);
            record["published_generator_order_multiplier"] = json!(generator_multiplier);
            record["published_generator"] = json!(generator);
            record["extension_factorization"] = json!([5, base_degree]);
            record["ghs_magic_number"] = json!(5);
            record["ghs_genus_formula"] = json!(15);
            record["ghs_base_degree"] = json!(base_degree);
            record["ghs_cover_constructed"] = json!(false);
            record["ghs_subgroup_map_certified"] = json!(false);
            record["projection"] = json!({"model":route["cost"]["destination"]["solver"],
                "log2_total":route["cost"]["log2_total"],
                "log2_rho_reference":transfer["reference_log2"],
                "hypothesis_satisfied":false,"status":"formula_only"});
            record["source_url"] = json!("https://www.rfc-editor.org/rfc/rfc2412.html");
            record["source_urls"] = json!([
                "https://www.rfc-editor.org/rfc/rfc2412.html",
                "https://www.rfc-editor.org/rfc/rfc2409.html",
                "https://github.com/aburan28/crypto/blob/main/docs/curves/JACOBIAN_TRANSFERS.md"
            ]);
            record["subgroup_advantage"] =
                json!("GHS invariant exact; cover, transport and qualified cost open");
            record["findings"] = json!([finding]);
        }
        records.push(record);
    }
    records.sort_by(|a, b| a["id"].as_str().cmp(&b["id"].as_str()));
    if records.len() != 103 {
        return Err("named curve registry count drifted".into());
    }
    Ok(records)
}

fn mul(a: u64, b: u64, p: u64) -> u64 {
    (a * b) % p
}
fn rhs(x: u64, a: u64, b: u64, p: u64) -> u64 {
    (mul(mul(x, x, p), x, p) + mul(a, x, p) + b) % p
}
fn nonsingular(a: u64, b: u64, p: u64) -> bool {
    !(4 * mul(mul(a, a, p), a, p) + 27 * mul(b, b, p)).is_multiple_of(p)
}
fn root_counts(p: u64) -> Vec<u64> {
    let mut counts = vec![0; p as usize];
    for y in 0..p {
        counts[mul(y, y, p) as usize] += 1;
    }
    counts
}
fn point_count(a: u64, b: u64, p: u64, roots: &[u64]) -> u64 {
    1 + (0..p).map(|x| roots[rhs(x, a, b, p) as usize]).sum::<u64>()
}
type Point = Option<(u64, u64)>;
fn inv(x: u64, p: u64) -> u64 {
    // The calling fields have prime cardinality; Fermat inverse.
    let mut base = x;
    let mut exp = p - 2;
    let mut out = 1;
    while exp > 0 {
        if exp & 1 == 1 {
            out = mul(out, base, p);
        }
        base = mul(base, base, p);
        exp >>= 1;
    }
    out
}
fn add(a: Point, b: Point, curve_a: u64, p: u64) -> Point {
    let (x1, y1) = match a {
        Some(v) => v,
        None => return b,
    };
    let (x2, y2) = match b {
        Some(v) => v,
        None => return a,
    };
    if x1 == x2 && (y1 + y2) % p == 0 {
        return None;
    }
    let (num, den) = if x1 == x2 {
        ((3 * mul(x1, x1, p) + curve_a) % p, (2 * y1) % p)
    } else {
        ((y2 + p - y1) % p, (x2 + p - x1) % p)
    };
    let slope = mul(num, inv(den, p), p);
    let x3 = (mul(slope, slope, p) + 2 * p - x1 - x2) % p;
    let y3 = (mul(slope, (x1 + p - x3) % p, p) + p - y1) % p;
    Some((x3, y3))
}
fn point_order(point: Point, a: u64, p: u64, bound: u64) -> u64 {
    let mut acc = None;
    for k in 1..=bound {
        acc = add(acc, point, a, p);
        if acc.is_none() {
            return k;
        }
    }
    panic!("point order exceeds curve order")
}
fn first_point(a: u64, b: u64, p: u64) -> Point {
    for x in 0..p {
        for y in 0..p {
            if mul(y, y, p) == rhs(x, a, b, p) {
                return Some((x, y));
            }
        }
    }
    None
}
fn prime_examples() -> Vec<Value> {
    let p = 101;
    let roots = root_counts(p);
    let mut anomalous = None;
    let mut smooth = None;
    for a in 0..p {
        for b in 0..p {
            if !nonsingular(a, b, p) {
                continue;
            }
            let order = point_count(a, b, p, &roots);
            if order == p && anomalous.is_none() {
                let point = first_point(a, b, p).unwrap();
                assert_eq!(point_order(Some(point), a, p, order), p);
                assert_eq!((a, b), (1, 32));
                anomalous = Some(json!({"id":"ICV1:fp-101:1:101:99:unk:unk:r:da77b88029ba",
                    "slug":"icv1-fp7-t1-da77b880",
                    "legacy_id":format!("EX1:F101:a{a}:b{b}:anomalous"),
                    "p":"101","degree":1,"bits":7,"a":a,"b":b,"order":order,
                    "trace":1,"generator":point,"generator_order":p,
                    "signal":"anomalous_prime_field","attack_ready":true,
                    "proof":"exhaustive_point_count","subgroup_advantage":"published_trace_one_reduction"}));
            }
            if order == 100 && smooth.is_none() {
                let mut best = (None, 0);
                for x in 0..p {
                    for y in 0..p {
                        if mul(y, y, p) == rhs(x, a, b, p) {
                            let point = Some((x, y));
                            let ord = point_order(point, a, p, order);
                            if ord > best.1 {
                                best = (point, ord);
                            }
                        }
                    }
                }
                assert_eq!(best.1, 10);
                assert_eq!((a, b), (1, 0));
                smooth = Some(json!({"id":"ICV1:fp-101:2:100:11:unk:unk:r:6861181a950c",
                    "slug":"icv1-fp7-t2-6861181a",
                    "legacy_id":format!("EX1:F101:a{a}:b{b}:smooth"),
                    "p":"101","degree":1,"bits":7,"a":a,"b":b,"order":order,
                    "trace":2,"generator":best.0,"generator_order":best.1,
                    "curve_order_factors":{"2":2,"5":2},
                    "generator_order_factors":{"2":1,"5":1},
                    "group_exponent":best.1,"group_structure":[10,10],
                    "signal":"smooth_selected_group",
                    "attack_ready":true,"proof":"exhaustive_point_count_and_group_order",
                    "subgroup_advantage":"pohlig_hellman_applicable"}));
            }
            if anomalous.is_some() && smooth.is_some() {
                break;
            }
        }
        if anomalous.is_some() && smooth.is_some() {
            break;
        }
    }
    let p2 = 103;
    let a2 = 1;
    let b2 = 0;
    let order2 = point_count(a2, b2, p2, &root_counts(p2));
    assert_eq!(order2, 104);
    let mut pair_point = None;
    for x in 0..p2 {
        for y in 0..p2 {
            if mul(y, y, p2) == rhs(x, a2, b2, p2)
                && point_order(Some((x, y)), a2, p2, order2) == 13
            {
                pair_point = Some((x, y));
                break;
            }
        }
        if pair_point.is_some() {
            break;
        }
    }
    assert!(pair_point.is_some());
    vec![
        anomalous.unwrap(),
        smooth.unwrap(),
        json!({
            "id":"ICV1:fp-103:0:104:80:unk:unk:r:8c4948d207a4",
            "slug":"icv1-fp7-t0-8c4948d2",
            "legacy_id":"EX1:F103:a1:b0:pairing_candidate",
            "p":"103","degree":1,
            "bits":7,"a":a2,"b":b2,"order":order2,"trace":0,
            "generator":pair_point,"generator_order":13,"embedding_degree":2,
            "signal":"low_embedding_degree_candidate","attack_ready":false,
            "proof":"exhaustive_point_count_group_order_and_modular_order",
            "subgroup_advantage":"pairing_map_and_destination_cost_open"
        }),
    ]
}

fn snapshot() -> Result<Value> {
    let db = connect(DB)?;
    let case_db = connect(CASES)?;
    let mut models = models(&db)?;
    let cases = cases(&case_db)?;
    let mut examples = prime_examples();
    attach_findings(&mut models);
    attach_findings(&mut examples);
    let named = named_curves()?;
    let mut curve_ids = std::collections::HashSet::new();
    let mut finding_ids = std::collections::HashSet::new();
    for record in models.iter().chain(&examples).chain(&named) {
        let id = record["id"].as_str().ok_or("curve ID missing")?;
        if !curve_ids.insert(id) {
            return Err(format!("duplicate curve ID: {id}").into());
        }
        for entry in record["findings"]
            .as_array()
            .ok_or("finding list missing")?
        {
            let id = entry["id"].as_str().ok_or("finding ID missing")?;
            if !finding_ids.insert(id) {
                return Err(format!("duplicate finding ID: {id}").into());
            }
        }
    }
    if finding_ids.len() != 29 {
        return Err("finding count drifted".into());
    }
    if curve_ids.len() != 920 {
        return Err("curve identity count drifted".into());
    }
    let positive = models
        .iter()
        .filter(|m| m["signal"] == "exact_covering_model")
        .count();
    let large = models
        .iter()
        .filter(|m| m["bits"].as_i64().unwrap() >= 192)
        .count();
    let large_survivors = models
        .iter()
        .filter(|m| {
            m["bits"].as_i64().unwrap() >= 192 && m["cubic_screen"] == "necessary_condition_passed"
        })
        .count();
    let large_depth_excluded = models
        .iter()
        .filter(|m| {
            m["bits"].as_i64().unwrap() >= 192 && m["cubic_screen"] == "excluded_by_depth_condition"
        })
        .count();
    let large_outside = models
        .iter()
        .filter(|m| {
            m["bits"].as_i64().unwrap() >= 192 && m["cubic_screen"] == "outside_recorded_domain"
        })
        .count();
    assert_eq!(
        (
            models.len(),
            positive,
            large,
            large_survivors,
            large_depth_excluded,
            large_outside,
            cases.len()
        ),
        (814, 24, 537, 247, 137, 153, 36)
    );
    Ok(json!({"schema":"weakcurves.snapshot.v2", "sources": {
        "model_database_sha256": file_hash(DB)?, "case_database_sha256": file_hash(CASES)?,
        "curve_registry_sha256":file_hash(CURVE_REGISTRY)?,
        "transfer_registry_sha256":file_hash(TRANSFER_REGISTRY)?},
        "counts":{"exact_models":models.len(),"exact_covering_models":positive,
            "large_models":large,"large_cubic_screen_survivors":large_survivors,
            "large_cubic_depth_excluded":large_depth_excluded,"large_outside_domain":large_outside,
            "published_cases":cases.len(),"prime_field_examples":examples.len(),
            "named_curves":named.len(),"oakley_ghs_candidates":2,
            "findings":finding_ids.len()},
        "models":models,"published_cases":cases,"prime_examples":examples,
        "named_curves":named}))
}

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 || !matches!(args[1].as_str(), "build" | "check") {
        return Err("usage: weak-curve-leaderboard <build|check> <output.json>".into());
    }
    let data = serde_json::to_vec_pretty(&snapshot()?)?;
    if args[1] == "build" {
        fs::write(&args[2], data)?;
    } else if fs::read(Path::new(&args[2]))? != data {
        return Err("snapshot differs from frozen source databases".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn prime_field_examples_are_exact() {
        let records = prime_examples();
        assert_eq!(records.len(), 3);
        for r in records {
            let p = r["p"].as_str().unwrap().parse::<u64>().unwrap();
            let a = r["a"].as_u64().unwrap();
            let b = r["b"].as_u64().unwrap();
            assert!(nonsingular(a, b, p));
            assert_eq!(
                point_count(a, b, p, &root_counts(p)),
                r["order"].as_u64().unwrap()
            );
            let point: (u64, u64) = serde_json::from_value(r["generator"].clone()).unwrap();
            assert_eq!(
                point_order(Some(point), a, p, r["order"].as_u64().unwrap()),
                r["generator_order"].as_u64().unwrap()
            );
        }
    }
    #[test]
    fn pairing_candidate_has_exact_degree_two() {
        assert_eq!(103 % 13, 12);
        assert_eq!((103 * 103) % 13, 1);
    }
    #[test]
    fn named_curve_findings_preserve_exact_subjects() {
        let records = named_curves().unwrap();
        assert_eq!(records.len(), 103);
        let oakley: Vec<_> = records
            .iter()
            .filter(|c| c["name"].as_str().unwrap().starts_with("Oakley Group "))
            .collect();
        assert_eq!(oakley.len(), 2);
        let mut ids = std::collections::HashSet::new();
        for c in oakley {
            assert_eq!(c["ghs_magic_number"], 5);
            assert_eq!(c["ghs_genus_formula"], 15);
            assert_eq!(c["projection"]["hypothesis_satisfied"], false);
            let f = &c["findings"][0];
            assert!(f["id"].as_str().unwrap().starts_with("WID1:GHS:"));
            assert!(ids.insert(f["id"].as_str().unwrap()));
        }
    }
    #[test]
    fn finding_id_is_stable_across_grade_changes() {
        let a = finding(
            "ICV1:example",
            "GHS",
            "ghs_binary",
            "selected_subgroup",
            "r=13;tower=5x31;g=15",
            "candidate",
            "first statement",
        );
        let b = finding(
            "ICV1:example",
            "GHS",
            "ghs_binary",
            "selected_subgroup",
            "r=13;tower=5x31;g=15",
            "verified_map",
            "updated statement",
        );
        assert_eq!(a["id"], b["id"]);
        assert_ne!(
            a["id"],
            finding(
                "ICV1:other",
                "GHS",
                "ghs_binary",
                "selected_subgroup",
                "r=13;tower=5x31;g=15",
                "candidate",
                "first statement"
            )["id"]
        );
    }
}

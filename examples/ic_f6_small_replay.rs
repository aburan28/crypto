//! Independent group replay and phase-ledger check of the frozen n9 cold runs.
use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use num_bigint::BigUint;
use serde_json::{json, Value};

fn decimal(value: &Value) -> BigUint {
    value
        .as_str()
        .expect("decimal string")
        .parse()
        .expect("valid decimal")
}

fn point(value: &Value) -> BinaryPoint {
    BinaryPoint::Affine {
        x: F2mElement::from_biguint(&decimal(&value[0]), 9),
        y: F2mElement::from_biguint(&decimal(&value[1]), 9),
    }
}

fn phase_sum(value: &Value) -> u64 {
    value
        .as_object()
        .expect("phase map")
        .values()
        .filter_map(Value::as_u64)
        .sum()
}

fn main() {
    let curve = KoblitzCurve::new(1, 9).expect("frozen n9 curve");
    let mut rows = Vec::new();
    for path in std::env::args().skip(1) {
        let raw: Value =
            serde_json::from_slice(&std::fs::read(&path).expect("run file")).expect("valid JSON");
        let target = point(&raw["fixture"]["targets"][0]);
        let generator = point(&raw["fixture"]["generator"]);
        let scalar = decimal(&raw["solutions"][0]["recovered"]);
        let group_replay = raw["status"] == "complete"
            && scalar < curve.subgroup_order
            && generator == curve.generator().clone()
            && curve.mul(&generator, &scalar) == target;
        let online_ns = raw["online_wall_ns"].as_u64().expect("online interval");
        let timing = &raw["generic_phase_timing"];
        let phase_closure = phase_sum(&timing["online_phases_ns"]) == online_ns
            && phase_sum(&timing["phases_ns"])
                == timing["observed_wall_ns"].as_u64().expect("cold interval");
        let relation_certified =
            raw["log_table_report"]["verified"] == true && raw["rejected_relations"] == 0;
        let base: Vec<_> = raw["factor_base"]
            .as_array()
            .expect("factor base")
            .iter()
            .map(point)
            .collect();
        let relations_replayed = raw["relations"]
            .as_array()
            .expect("ordinary relations")
            .iter()
            .all(|relation| {
                let mut sum = BinaryPoint::Infinity;
                for index in relation["points"].as_array().expect("point indices") {
                    let point = &base[index.as_u64().expect("index") as usize];
                    sum = curve.add(&sum, point);
                }
                let a = BigUint::from(relation["a"].as_u64().expect("relation scalar"));
                curve.mul(&generator, &a) == sum
            });
        let column_logs_replayed = raw["column_logs"]
            .as_array()
            .expect("certified column logs")
            .iter()
            .all(|column| {
                curve.mul(&generator, &decimal(&column["log"])) == point(&column["point"])
            });
        assert!(
            group_replay
                && phase_closure
                && relation_certified
                && relations_replayed
                && column_logs_replayed,
            "{path}"
        );
        rows.push(json!({
            "run_file":path,
            "recovered":scalar.to_string(),
            "target":raw["fixture"]["targets"][0],
            "group_replay":group_replay,
            "phase_closure":phase_closure,
            "relation_certified":relation_certified,
            "relations_replayed":relations_replayed,
            "column_logs_replayed":column_logs_replayed,
            "cold_inside_worker_ns":timing["observed_wall_ns"],
            "online_ns":online_ns
        }));
    }
    assert_eq!(rows.len(), 3, "three frozen arms");
    println!(
        "{}",
        json!({"schema_version":1,"status":"verified","runs":rows})
    );
}

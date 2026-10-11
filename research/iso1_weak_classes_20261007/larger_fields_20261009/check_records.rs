//! Documentary ICV1 and raw-receipt consistency; curve arithmetic is performed by GP.
use crypto_lib::hash::sha256::Sha256;
use serde_json::Value;
use std::{collections::BTreeSet, fs, path::Path};
fn sha(bytes: &[u8]) -> String {
    let mut h = Sha256::new();
    h.update(bytes);
    h.finalize().iter().map(|b| format!("{b:02x}")).collect()
}
fn check(path: &Path) {
    let data: Value = serde_json::from_slice(&fs::read(path).unwrap()).unwrap();
    assert_eq!(data["paired_counts"].as_u64(), Some(17));
    assert_eq!(data["geometry_controls"].as_u64(), Some(49));
    let mut unique = BTreeSet::new();
    for record in data["records"].as_array().unwrap() {
        let p = record["characteristic"].as_u64().unwrap();
        let k = record["field_degree"].as_u64().unwrap();
        let modulus: Vec<_> = record["modulus_low_coefficients"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_u64().unwrap().to_string())
            .collect();
        let field = format!(
            "fpk-{p}-{k}-{}",
            &sha(format!("fpk-modulus:{p}:{}", modulus.join(",")).as_bytes())[..8]
        );
        for which in ["source", "target"] {
            let c = &record[which];
            let canonical = c["model_json"].as_str().unwrap();
            let parsed: Value = serde_json::from_str(canonical).unwrap();
            assert_eq!(serde_json::to_string(&parsed).unwrap(), canonical);
            let h = sha(canonical.as_bytes());
            assert_eq!(c["model_sha256"].as_str().unwrap(), h);
            assert_eq!(parsed["field"].as_str().unwrap(), field);
            let j: Vec<_> = c["j_coefficients"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_str().unwrap())
                .collect();
            let trace = c["trace"].as_str().unwrap();
            let order = c["group_order"].as_str().unwrap();
            let id = format!(
                "ICV1:{field}:{trace}:{order}:{}:unk:unk:r:{}",
                j.join(","),
                &h[..12]
            );
            assert_eq!(c["icv1"].as_str().unwrap(), id);
            let signed = if let Some(t) = trace.strip_prefix('-') {
                format!("tm{t}")
            } else {
                format!("t{trace}")
            };
            let slug = format!("icv1-fp{}k{k}-{signed}-{}", 64 - p.leading_zeros(), &h[..8]);
            assert_eq!(c["slug"].as_str().unwrap(), slug);
            assert!(unique.insert(id));
            assert!(c["subgroup_order"].is_null() && c["volcano_level"].is_null());
        }
    }
    assert_eq!(unique.len(), 34);
    let base = path.parent().unwrap();
    let receipts = fs::read_to_string(base.join("evidence_run1/receipts.tsv")).unwrap();
    for row in receipts.lines().skip(1) {
        let columns: Vec<_> = row.split('\t').collect();
        assert_eq!(columns[7], "COMPLETE");
        for (kind, col) in [("stdout", 10), ("stderr", 11)] {
            let file = base.join(format!(
                "evidence_run1/p{}_n{}_{}.{}",
                columns[0], columns[2], columns[3], kind
            ));
            assert_eq!(sha(&fs::read(file).unwrap()), columns[col]);
        }
    }
    println!("PASS: 34 canonical ICV1 identities and 68 raw output hashes");
}
fn main() {
    check(Path::new(
        &std::env::args().nth(1).expect("curve_records.json"),
    ));
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn retained_models_and_receipts_agree() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
        check(&root.join("larger_fields_20261009/curve_records.json"));
    }
}

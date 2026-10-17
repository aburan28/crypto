//! Incremental native catalogue refresh. Never recompute or alter baseline rows.
//! Full debug replay currently hits an existing binary-subfield trace assertion.
use crypto_lib::cryptanalysis::curve_traits::{self, Factoring, TraitFile};
use serde_json::json;
use std::{collections::BTreeMap, fs};

fn digest(bytes: &[u8]) -> String {
    crypto_lib::hash::sha256::sha256(bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}
fn main() {
    let text = fs::read_to_string("docs/curves/traits.json").unwrap();
    let mut file: TraitFile = serde_json::from_str(&text).unwrap();
    let old_hash = digest(text.as_bytes());
    assert_eq!(
        old_hash, "7395b479aea3e90f3fd2010dd8d4f570c60dc4fb2baf83bda97165b0e659b844",
        "start from the recorded baseline traits file"
    );
    let registry = fs::read_to_string("docs/curves/registry.json").unwrap();
    // Reusing old rows is legitimate only if every other input is unchanged.
    for input in &file.inputs {
        if input.path != "docs/curves/registry.json" {
            assert_eq!(
                input.sha256,
                digest(&fs::read(&input.path).unwrap()),
                "changed baseline input {}",
                input.path
            );
        }
    }
    let before: BTreeMap<_, _> = file
        .curves
        .iter()
        .map(|c| (c.slug.clone(), c.clone()))
        .collect();
    let hints =
        curve_traits::read_hints(&fs::read_to_string("docs/curves/factor-hints.txt").unwrap())
            .unwrap();
    let fx = Factoring::new(file.factor_budget, hints);
    let mut rows = vec![];
    let mut added = vec![];
    let mut remaining = before.clone();
    for c in curve_traits::read_registry(&registry).unwrap() {
        if let Some(row) = remaining.remove(&c.slug) {
            rows.push(row);
        } else {
            assert!(c.slug.starts_with("icv1-fp57-") || c.slug.starts_with("icv1-fp59-"));
            added.push(c.slug.clone());
            rows.push(curve_traits::compute(&c, &fx).unwrap());
        }
    }
    assert_eq!(added.len(), 6);
    for row in &file.curves {
        if remaining.contains_key(&row.slug) {
            rows.push(row.clone());
        }
    }
    for row in &rows {
        if let Some(original) = before.get(&row.slug) {
            assert_eq!(row, original);
        }
    }
    let registry_hash = digest(registry.as_bytes());
    file.registry_sha256 = registry_hash.clone();
    let input = file
        .inputs
        .iter_mut()
        .find(|i| i.path == "docs/curves/registry.json")
        .unwrap();
    input.sha256 = registry_hash;
    file.curves = rows;
    let refreshed = file.to_json();
    fs::write("docs/curves/traits.json", &refreshed).unwrap();
    println!(
        "{}",
        json!({"status":"INCREMENTAL_NATIVE_REFRESH","baseline_sha256":old_hash,"output_sha256":digest(refreshed.as_bytes()),"added":added,"baseline_rows_preserved":before.len(),"other_input_hashes_unchanged":true,"full_debug_replay":"BLOCKED_BY_EXISTING_BINARY_SUBFIELD_ASSERTION; no assertion disabled"})
    );
}

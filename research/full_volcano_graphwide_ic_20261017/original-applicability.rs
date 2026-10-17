//! Native arithmetic inventory, not a pairing execution or DLP measurement.
use crypto_lib::cryptanalysis::curve_traits::arith;
use serde_json::{json, Value};
use std::{env, fs, path::Path};

fn main() {
    let a: Vec<_> = env::args().collect();
    assert_eq!(
        a.len(),
        3,
        "original-applicability INPUT_DIRECTORY FRESH_JSON"
    );
    let out = Path::new(&a[2]);
    assert!(!out.exists());
    let mut rows = vec![];
    for label in ["A", "B", "C", "D", "E"] {
        let path = Path::new(&a[1]).join(format!("{label}.json"));
        let v: Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
        let p = v["p"].as_u64().unwrap();
        let n = v["order"]
            .as_u64()
            .or_else(|| v["group_order"].as_u64())
            .unwrap();
        let factors = arith::factor(&n.into(), 1_000_000);
        assert!(
            factors.complete(),
            "incomplete order factorization for {label}"
        );
        let r: u64 = factors
            .largest_prime()
            .unwrap()
            .to_string()
            .parse()
            .unwrap();
        let r1 = arith::factor(&(r - 1).into(), 1_000_000);
        assert!(r1.complete());
        let k: u64 = arith::multiplicative_order(&p.into(), &r.into(), &r1)
            .to_string()
            .parse()
            .unwrap();
        let coefficient_bytes = (u64::BITS - p.leading_zeros()).div_ceil(8);
        rows.push(json!({"case":label,"input":v,"order_factorization":factors.primes.iter().map(|(q,e)|json!([q.to_string(),e])).collect::<Vec<_>>(),
            "largest_prime_subgroup":r,"cofactor":n/r,"r_divides_order_once":(n/r)%r!=0,
            "r_minus_1_factorization":r1.primes.iter().map(|(q,e)|json!([q.to_string(),e])).collect::<Vec<_>>(),
            "embedding_degree":k,"destination":"F_(p^k)",
            "dense_packed_coefficient_bytes_per_element":(k as u128 * coefficient_bytes as u128).to_string(),
            "dense_u64_coefficient_bytes_per_element":(k as u128 * 8).to_string(),
            "representation_bound_scope":"Specified dense coefficient layouts only; not a bound on other representations or algorithms",
            "pairing_execution":"NOT_RUN","finite_field_index_calculus":"NOT_RUN"}));
    }
    fs::write(out, serde_json::to_vec_pretty(&json!({"schema":"original-volcano-ff-applicability/v1","scope":"Exact original inputs, largest prime subgroups; arithmetic inventory only","rows":rows})).unwrap()).unwrap();
}

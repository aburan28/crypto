//! The exact field and curve records of the Koblitz curves the index-calculus
//! pipeline builds, for EC1 curve identities (docs/curve-identities.md).
//!
//! `KoblitzCurve::new(a, n)` fixes the field's defining polynomial, the
//! prime-order subgroup and its generator.  A parameter file names only
//! `(a, n)`, so its curve's identity is read off the library here rather
//! than inferred from the degree.  The records use the encoding of the
//! repository's existing EC1 records for `E_0`: modulus exponents from
//! the degree down, field elements as `0x`-prefixed polynomial-basis
//! integers, and `[a1, a2, a3, a4, a6]` for `y² + a1·xy + a3·y = x³ +
//! a2·x² + a4·x + a6`.  It performs no certification: `tools/curve_identity.py`
//! hashes the metadata, and nothing more.
//!
//!     cargo run --release --example koblitz_curve_records -- 1:19 0:37 ...
use crypto_lib::binary_ecc::BinaryPoint;
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use serde_json::{json, Value};

fn record(a: u8, n: u32) -> Value {
    let kc = KoblitzCurve::new(a, n).expect("a Koblitz curve");
    let c = &kc.curve;
    assert_eq!(c.m, n);
    assert_eq!(
        c.order, kc.subgroup_order,
        "the generator's order is the subgroup's"
    );
    assert_eq!(c.cofactor, kc.cofactor);
    let mut low = c.irreducible.low_terms.clone();
    low.sort_unstable_by(|x, y| y.cmp(x));
    let mut exponents = vec![n];
    exponents.extend(low);
    let (x, y) = match &c.generator {
        BinaryPoint::Affine { x, y } => (x.to_biguint(), y.to_biguint()),
        BinaryPoint::Infinity => panic!("the generator is the point at infinity"),
    };
    let cofactor = u64::try_from(&kc.cofactor).expect("a small cofactor");
    json!({
        "a": a,
        "n": n,
        "field": {
            "characteristic": 2,
            "degree": n,
            "representation": "polynomial",
            "modulus_exponents": exponents,
            "element_encoding": "hex polynomial coefficient bitset, least significant bit is constant",
        },
        "curve": {
            "model": "binary Weierstrass",
            "coefficients": [1, a, 0, 0, 1],
            "subgroup_order": kc.subgroup_order.to_string(),
            "cofactor": cofactor,
            "generator": [format!("{x:#x}"), format!("{y:#x}")],
            "target_group": "prime-order subgroup",
        },
        "group_order": kc.group_order.to_string(),
    })
}

fn main() {
    let curves: Vec<Value> = std::env::args()
        .skip(1)
        .map(|arg| {
            let (a, n) = arg.split_once(':').expect("a:n");
            record(a.parse().expect("a"), n.parse().expect("n"))
        })
        .collect();
    println!(
        "{}",
        serde_json::to_string_pretty(&json!({ "curves": curves })).expect("json")
    );
}

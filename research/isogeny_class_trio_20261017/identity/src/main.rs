//! Prints the ICV1 identity of each trio pick with the repository's native
//! producer (src/cryptanalysis/curve_id.rs). Inputs are the recorded picks
//! in ../data/pick-*.json; the values are passed explicitly below.

use crypto_lib::cryptanalysis::curve_id;
use num_bigint::BigUint;

fn big(s: &str) -> BigUint {
    if let Some(h) = s.strip_prefix("0x") {
        BigUint::parse_bytes(h.as_bytes(), 16).unwrap()
    } else {
        BigUint::parse_bytes(s.as_bytes(), 10).unwrap()
    }
}

fn show(name: &str, id: Option<curve_id::CurveId>) {
    let id = id.expect("identity rejected");
    println!(
        "{{\"pick\":\"{name}\",\"icv1\":\"{}\",\"slug\":\"{}\",\"model_json\":{:?}}}",
        id.icv1, id.slug, id.model_json
    );
}

fn main() {
    // K-41: x^41 + x^3 + 1, a = 0, b = 1, #E = 2199025563772
    show("K-41", curve_id::koblitz(0, 41, &big("0x20000000009"), &big("2199025563772")));
    // B-47: x^47 + x^5 + 1, a2 = 1, b = 0x52e3db551eb3, #E = 140737484079278;
    // End(E) = Z[pi] (root on the 31-volcano floor, f = 31^2): disc = Delta.
    show(
        "B-47",
        curve_id::binary(47, &big("0x800000000021"), &big("1"), &big("0x52e3db551eb3"), &big("140737484079278"), Some(-544665341266711)),
    );
    // P-47: p = 2^47 - 115, a = -3, b = 0xda8a554f2ba, #E = 140737497140563
    let p = big("140737488355213");
    show("P-47", curve_id::prime(&p, &(&p - 3u32), &big("0xda8a554f2ba"), &big("140737497140563")));
    // Alternates (61-bit): B-61 x^61 + x^5 + x^2 + x + 1, a2 = 1; P-61 over 2^61 - 1; K-61 a = 0.
    show("K-61", curve_id::koblitz(0, 61, &big("0x2000000000000027"), &big("2305843009055095052")));
    show(
        "B-61",
        curve_id::binary(61, &big("0x2000000000000027"), &big("1"), &big("0x9c60f836d06fc8"), &big("2305843011871087922"), Some(-2161629330377202847)),
    );
    let p61 = big("2305843009213693951");
    show("P-61", curve_id::prime(&p61, &(&p61 - 3u32), &big("0x54aa66b458b77ec"), &big("2305843010027268529")));
}

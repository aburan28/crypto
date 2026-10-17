//! Recompute exact ICV1 model and EC1 subgroup-representation identities.
use crypto_lib::cryptanalysis::{curve_id, isogeny_walk::record};
use num_bigint::BigUint;

fn main() {
    let p = BigUint::from(617u32);
    let zero = BigUint::from(0u32);
    let order = BigUint::from(618u32);
    let r = BigUint::from(103u32);
    let h = BigUint::from(6u32);
    for (name, b, x, y) in [
        ("source", 1u32, 66u32, 264u32),
        ("destination", 590, 483, 13),
    ] {
        let b = BigUint::from(b);
        let x = BigUint::from(x);
        let y = BigUint::from(y);
        let id = curve_id::prime(&p, &zero, &b, &order).expect("valid curve identity");
        let (field, model) = record::ec1_records(&p, &zero, &b, &r, &h, (&x, &y));
        let (ec1, uid) = record::ec1_identity(&p, &field, &model, "fp");
        println!(
            "{name}\ticv1={}\tslug={}\tec1={}\tuid={}\tmodel_json={}",
            id.icv1, id.slug, ec1, uid, id.model_json
        );
    }
}

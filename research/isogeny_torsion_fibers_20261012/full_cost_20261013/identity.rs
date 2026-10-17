//! Emit repository ICV1 and EC1 identities for the four verified pairs.
//! Link to the release crypto_lib and num_bigint rlibs built from this tree.
use crypto_lib::cryptanalysis::curve_id;
use crypto_lib::cryptanalysis::isogeny_walk::record::{ec1_identity, ec1_records};
use num_bigint::BigUint;

fn main() {
    println!("p\tr\tarm\tb\tgx\tgy\ticv1\tslug\tec1\tcurve_uid");
    let rows = [
        (617u64, 103u64, [(1, 66, 264), (590, 483, 13)]),
        (2633, 439, [(1, 240, 728), (2606, 2453, 643)]),
        (9677, 1613, [(1, 3188, 1881), (9650, 4515, 3563)]),
        (38561, 6427, [(1, 38503, 10141), (38534, 35110, 19035)]),
    ];
    for (p, r, pair) in rows {
        let p = BigUint::from(p);
        let r = BigUint::from(r);
        let order = &r * BigUint::from(6u8);
        for (which, (b, gx, gy)) in pair.into_iter().enumerate() {
            let a = BigUint::from(0u8);
            let b = BigUint::from(b as u64);
            let gx = BigUint::from(gx as u64);
            let gy = BigUint::from(gy as u64);
            let id = curve_id::prime(&p, &a, &b, &order).expect("valid model and order");
            let (field, curve) = ec1_records(&p, &a, &b, &r, &BigUint::from(6u8), (&gx, &gy));
            let (ec1, uid) = ec1_identity(&p, &field, &curve, "j0");
            println!(
                "{p}\t{r}\t{}\t{b}\t{gx}\t{gy}\t{}\t{}\t{ec1}\t{uid}",
                if which == 0 { "source" } else { "destination" },
                id.icv1,
                id.slug
            );
        }
    }
}

//! Independent exact multiset support check for the frozen N83 panel sizes.
//! Uses only standard-library integer arithmetic; no group or solver claims.

fn multiset_count(points: u128, summands: u128) -> u128 {
    let mut count = 1u128;
    for i in 1..=summands {
        count = count
            .checked_mul(points + i - 1)
            .expect("count product fits u128")
            / i;
    }
    count
}

fn main() {
    let curves = [
        ("a0", 2_417_851_639_230_796_216_685_689u128),
        ("a1", 8_569_786_107_849_059u128),
    ];
    println!("curve,K,m,points,multisets,nonidentity_order,hit_ceiling");
    for (curve, order) in curves {
        for k in [64u128, 256, 600] {
            let points = 166 * k;
            for m in 2..=4 {
                let count = multiset_count(points, m);
                let ceiling = ((count as f64) / ((order - 1) as f64)).min(1.0);
                println!("{curve},{k},{m},{points},{count},{},{ceiling:.9e}", order - 1);
            }
        }
    }
    let structured = multiset_count(332_166, 4);
    println!(
        "a0_structured_d12,2001,4,332166,{structured},{},{:.9e}",
        curves[0].1 - 1,
        (structured as f64) / ((curves[0].1 - 1) as f64)
    );
}

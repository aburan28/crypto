//! Bounded prime-field baseline using the existing native relation collector.
//! This does not reproduce an unavailable isogeny/candidate-selection driver.
//! Usage: prime61_factor_base [bits=61] [trials=100000] [runs=12] [base=3]
use crypto_lib::cryptanalysis::{
    curve_id,
    prime_fast::{collect_relations, find_a3_curve, solve_relations, FastFactorBase},
};
use num_bigint::BigUint;
use serde_json::json;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let arg = |i: usize, default: u64| -> u64 {
        args.get(i)
            .map(|s| s.parse().expect("unsigned integer"))
            .unwrap_or(default)
    };
    let bits = u32::try_from(arg(1, 61)).unwrap();
    let trials = arg(2, 100_000);
    let runs = arg(3, 12);
    let size = usize::try_from(arg(4, 3)).unwrap();
    assert!((8..=61).contains(&bits) && trials > 0 && runs > 0 && size > 0);
    let (curve, g) = find_a3_curve(bits, 3).expect("prime-order curve search");
    if bits == 61 {
        assert_eq!(curve.f.p, (1u64 << 61) - 1);
    }
    let a = curve.f.from_mont(curve.a);
    let b = curve.f.from_mont(curve.b);
    let id = curve_id::prime(
        &BigUint::from(curve.f.p),
        &BigUint::from(a),
        &BigUint::from(b),
        &BigUint::from(curve.n),
    )
    .unwrap();
    let fb = FastFactorBase::build(&curve, size);
    assert_eq!(fb.len(), size);
    for &point in &fb.pts {
        assert!(curve.is_on_curve(point));
        assert!(curve.mul(Some(point), curve.n).is_none());
    }
    // A synthetic known-answer target. Its scalar is never passed to collection.
    let secret = 0x20261016u64 % (curve.n - 1) + 1;
    let q = curve.mul(Some(g), secret).unwrap();
    println!(
        "{}",
        json!({"kind":"manifest", "schema":"prime61.factor-base-baseline/v1",
        "curve":{"p":curve.f.p,"a":a,"b":b,"group_order":curve.n,"r":curve.n,
            "cofactor":1,"gx":curve.canonical(g).0,"gy":curve.canonical(g).1,
            "icv1":id.icv1,"slug":id.slug},
        "target":curve.canonical(q),"target_kind":"synthetic-known-answer",
        "target_scalar_for_replay_only":secret,"runs":runs,"max_targets_per_run":trials,
        "inverse_pairs":fb.len(),"factor_base":fb.pts.iter().map(|&p|curve.canonical(p)).collect::<Vec<_>>(),
        "summands":2,"candidate_search_performed":false,
        "isogeny_transport_performed":false,"wall_time_claim":null,
        "total_operation_cost":null,"speedup":null})
    );
    for run in 0..runs {
        let seed = 2026101600 + run;
        let needed = size + 17;
        let (relations, stats) = collect_relations(&curve, g, q, &fb, None, needed, trials, seed);
        // Independently replay every emitted relation using scalar multiplication
        // and ordinary additions, rather than the batched decomposition path.
        for relation in &relations {
            let lhs = curve.combine(g, q, relation.a, relation.b);
            let mut rhs = None;
            for &(index, sign) in &relation.entries {
                let p = fb.pts[index as usize];
                let p = if sign < 0 { curve.neg(p) } else { p };
                rhs = curve.add(rhs, curve.mul(Some(p), sign.unsigned_abs() as u64));
            }
            assert_eq!(lhs, rhs, "relation replay");
        }
        assert_eq!(stats.rejected, 0, "invalid decomposition");
        let recovered = solve_relations(&curve, size, &relations);
        let verified = recovered.is_some_and(|k| curve.mul(Some(g), k) == Some(q));
        if let Some(k) = recovered {
            assert!(verified && k == secret, "recovered scalar replay");
        }
        println!(
            "{}",
            json!({"kind":"run","run":run,"seed":seed,
            "targets":stats.targets,"relations":stats.relations,"needed":needed,
            "rows_one":stats.rows_one,"rows_two":stats.rows_two,"rows_three":stats.rows_three,
            "rejected":stats.rejected,"relation_replays":relations.len(),
            "scalar":recovered,"verified":verified,
            "status":if verified {"VERIFIED_SOLVE"} else {"EXHAUSTED_TRIAL_CAP"}})
        );
    }
}

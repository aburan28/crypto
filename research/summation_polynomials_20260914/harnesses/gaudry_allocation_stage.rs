//! Matched residual corpus for representation-only Gaudry optimizations.
use crypto_lib::cryptanalysis::gaudry_cubic::*;
use rand::{rngs::StdRng, SeedableRng};
use serde_json::json;
use std::{env, fs, time::Instant};

fn main() {
    let args: Vec<String> = env::args().collect();
    assert_eq!(args.len(), 5, "p seed count output.json");
    let p = args[1].parse().unwrap();
    let seed: u64 = args[2].parse().unwrap();
    let count: u64 = args[3].parse().unwrap();
    let start = Instant::now();
    let inst = generate_instance3(p, seed);
    let mut rng = StdRng::seed_from_u64(seed);
    let base = SubspaceBase::build(&inst, &mut rng);
    let pre = SymmetrisedS4::precompute(&inst.curve);
    let mut stats = SolveStats::default();
    let mut rows = Vec::new();
    for k in 1..=count {
        let target = inst.curve.mul(&inst.curve.g, k * 104_729 + 3);
        let actual = subspace_triple_oracle_groebner(
            &inst, &base, &pre, &target, &mut rng, &mut stats,
        );
        let (reference, _) = subspace_triple_oracle(&inst, &base, &target, &mut rng);
        // MITM also has two-term solutions, outside the S4 triple contract.
        let triples = |v: &Vec<Vec<(usize,i64)>>| {
            v.iter().filter(|t| t.len()!=2).cloned().collect::<Vec<_>>()
        };
        assert_eq!(triples(&actual), triples(&reference), "residual {k}");
        for terms in &actual {
            let mut residual = target;
            for &(i, coefficient) in terms {
                residual = inst.curve.sub(&residual, &inst.curve.mul_signed(&base.points[i], coefficient));
            }
            assert!(residual.inf);
        }
        rows.push(json!({"k":k,"target":target,"solutions":actual}));
    }
    fs::write(&args[4], serde_json::to_string_pretty(&json!({
        "p":p,"seed":seed,"count":count,"stats":stats,"rows":rows,
        "wall_s":start.elapsed().as_secs_f64()
    })).unwrap()).unwrap();
}

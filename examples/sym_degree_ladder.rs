//! Refutation-degree ladder for the torsion-symmetrised `m = 3` system.
//!
//! One process per cell `(K_a, n, ℓ)`.  Each draw picks a random subspace
//! `V ∋ 1` of dimension `ℓ` in the `u = 1/(x + 1)` frame (the frame in which
//! translation by the 2-torsion point is `u ↦ u + 1`), builds the symmetrised
//! factor base over it, and a uniform target in the prime-order subgroup.  The
//! symmetrised `S₄` system in the translation invariants has `3(ℓ − 1) + 1`
//! Boolean unknowns.  Its roots are counted exactly by brute force over the
//! cube, independently of every Macaulay code path; only a system with no
//! root is measured, by `solving_degree` up to `--d-max`.
//!
//! Output: one JSON line per draw.  The metric is the refutation degree `D`
//! (resolved, or a lower bound `≥ d_max + 1`); wall time is recorded, never the
//! metric.
//!
//!     sym_degree_ladder --a 1 --n 17 --ell 4 --unsat 4 --max-draws 256 \
//!         --d-max 8 --seed 20260929 --out cell.jsonl
use crypto_lib::binary_ecc::{BinaryPoint, F2mElement};
use crypto_lib::cryptanalysis::koblitz_bench::random_subspace_basis;
use crypto_lib::cryptanalysis::koblitz_groebner::{solving_degree, system_degree, FieldStructure};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use crypto_lib::cryptanalysis::koblitz_symmetrised::{
    build_symmetrised_factor_base_from_basis, build_symmetrised_system,
};
use crypto_lib::cryptanalysis::pq_groebner_f2::F2BoolPoly;
use num_bigint::BigUint;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::io::Write;
use std::time::Instant;

const M: usize = 3;

fn eval(p: &F2BoolPoly, a: u64) -> bool {
    p.terms
        .iter()
        .fold(false, |acc, t| acc ^ (t.mask & a == t.mask))
}

fn roots(polys: &[F2BoolPoly], n_vars: usize) -> u64 {
    (0u64..1 << n_vars)
        .filter(|&a| polys.iter().all(|p| !eval(p, a)))
        .count() as u64
}

fn word(e: &F2mElement) -> u64 {
    e.raw_bits().first().copied().unwrap_or(0)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let get = |name: &str| {
        args.iter()
            .position(|a| a == name)
            .and_then(|i| args.get(i + 1))
            .cloned()
    };
    let num = |name: &str, default: u64| -> u64 {
        get(name).map_or(default, |v| {
            v.parse().unwrap_or_else(|_| panic!("bad {name}"))
        })
    };
    let a = num("--a", 1) as u8;
    let n = num("--n", 17) as u32;
    let ell = num("--ell", 3) as usize;
    let want_unsat = num("--unsat", 4) as u32;
    let max_draws = num("--max-draws", 256) as u32;
    let d_max = num("--d-max", 8) as u32;
    let seed = num("--seed", 20260929);
    assert!((2..=8).contains(&ell), "need 2 ≤ ℓ ≤ 8");
    let mut out = get("--out").map(|p| {
        assert!(
            !std::path::Path::new(&p).exists(),
            "--out exists; never overwritten"
        );
        std::fs::File::create(p).expect("create --out")
    });

    let kc = KoblitzCurve::new(a, n).unwrap_or_else(|| panic!("no curve K_{a}/2^{n}"));
    let r64: u64 = kc.subgroup_order.iter_u64_digits().next().unwrap_or(0);
    assert!(
        kc.subgroup_order.bits() <= 63 && r64 > 2,
        "subgroup order out of range"
    );
    let st = FieldStructure::new(n, &kc.curve.irreducible);
    let cell_seed = seed ^ (u64::from(a) << 48) ^ (u64::from(n) << 40) ^ ((ell as u64) << 32);
    let mut rng = StdRng::seed_from_u64(cell_seed);
    let cell = format!("K{a}n{n}l{ell}");
    let g = kc.generator().clone();
    let (mut unsat, mut sat, mut skipped) = (0u32, 0u32, 0u32);

    for draw in 0..max_draws {
        if unsat >= want_unsat {
            break;
        }
        // V ∋ 1: 1 first, then ℓ − 1 random elements; a dependent draw, or one
        // whose Artin–Schreier image is dependent, is redrawn (counted).
        let mut v_basis = vec![F2mElement::one(n)];
        v_basis.extend(random_subspace_basis(n, ell - 1, &mut rng));
        let k = 1 + rng.gen::<u64>() % (r64 - 1);
        let Some(fb) = build_symmetrised_factor_base_from_basis(&kc, v_basis.clone()) else {
            skipped += 1;
            continue;
        };
        let target = kc.mul(&g, &BigUint::from(k));
        let Some(sys) = build_symmetrised_system(&kc, &fb, &target, M, &st) else {
            skipped += 1;
            continue;
        };
        let started = Instant::now();
        let count = roots(&sys.equations, sys.n_vars);
        let x_r = match &target {
            BinaryPoint::Affine { x, .. } => word(x),
            BinaryPoint::Infinity => unreachable!("k in [1, r-1]"),
        };
        let mut line = format!(
            r#"{{"cell":"{cell}","a":{a},"n":{n},"ell":{ell},"m":{M},"seed":{seed},"draw":{draw},"v_basis":{:?},"k":{k},"x_r":{x_r},"n_vars":{},"n_eqs":{},"system_degree":{},"roots":{count}"#,
            v_basis.iter().map(word).collect::<Vec<_>>(),
            sys.n_vars,
            sys.equations.len(),
            system_degree(&sys.equations),
        );
        if count > 0 {
            sat += 1;
            line.push_str(r#","outcome":{"kind":"satisfiable"}"#);
        } else {
            unsat += 1;
            let (d, profs) = solving_degree(&sys.equations, sys.n_vars, d_max);
            let built = profs.last().map(|p| p.degree);
            let refuted = profs.last().map(|p| p.refuted).unwrap_or(false);
            line.push_str(&match d {
                Some(degree) => format!(
                    r#","outcome":{{"kind":"resolved","degree":{degree},"refuted":{refuted}}}"#
                ),
                None if built == Some(d_max) => {
                    format!(r#","outcome":{{"kind":"at_least","degree":{}}}"#, d_max + 1)
                }
                None => format!(
                    r#","outcome":{{"kind":"caps","built":{}}}"#,
                    built.map_or("null".to_string(), |b| b.to_string())
                ),
            });
        }
        line.push_str(&format!(r#","secs":{}}}"#, started.elapsed().as_secs_f64()));
        println!("{line}");
        if let Some(f) = out.as_mut() {
            writeln!(f, "{line}").expect("write --out");
            f.flush().ok();
        }
    }
    eprintln!("{cell}: unsat {unsat}, satisfiable {sat}, redrawn {skipped}");
}

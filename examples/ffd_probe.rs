//! First-fall-degree probe for the Crossbred AutoLab X5 beats.
//!
//! ```bash
//! cargo run --release --example ffd_probe -- 9:3 9:4 15:3 15:4
//! FFD_TRIALS=16 cargo run --release --example ffd_probe -- 9:4 15:4
//! cargo run --release --example ffd_probe -- --sym --a 0 7:4 9:4 15:4
//! ```
//!
//! Chained-`x` uses `ffd_summary` with factor index 0, `d_max = 3`,
//! seed `0x5EED` (scaling-target convention).
//!
//! `--sym` prices the chained symmetrised `S₃`
//! (`build_chained_symmetrised_system`): divisor
//! `divisor_for_dimension(n, (n+1)/m)`, `d_max = 4` (bilinear, so a
//! fall at 4 would be invisible at the x-arm's `d_max = 3`), same seed.
//! That is not the unchained `S₄`.

use crypto_lib::binary_ecc::BinaryPoint;
use crypto_lib::cryptanalysis::koblitz_bench::{ffd_summary, format_ffd_table, FfdSummary};
use crypto_lib::cryptanalysis::koblitz_groebner::{
    first_fall_degree, system_degree, FieldStructure,
};
use crypto_lib::cryptanalysis::koblitz_index_calculus::KoblitzCurve;
use crypto_lib::cryptanalysis::koblitz_symmetrised::{
    build_chained_symmetrised_system, build_symmetrised_factor_base, divisor_for_dimension,
};
use num_bigint::BigUint;

fn parse_cases(args: &[String]) -> Vec<(u32, usize)> {
    args.iter()
        .filter_map(|a| {
            let (n, m) = a.split_once(':')?;
            Some((n.parse().ok()?, m.parse().ok()?))
        })
        .collect()
}

fn ffd_chained_sym(
    a: u8,
    n: u32,
    m: usize,
    d_max: u32,
    seed: u64,
    trials: usize,
) -> Option<FfdSummary> {
    let kc = KoblitzCurve::new(a, n).or_else(|| KoblitzCurve::new(1 - a, n))?;
    let target_dim = (n + 1).div_ceil(m as u32);
    let div = divisor_for_dimension(n, target_dim)?;
    let fb = build_symmetrised_factor_base(&kc, &div)?;
    let st = FieldStructure::new(n, &kc.curve.irreducible);
    let mut falls: Vec<Option<u32>> = Vec::with_capacity(trials);
    let mut syz = Vec::with_capacity(trials);
    let mut first: Option<(usize, usize, u32, u32)> = None;
    let mut t = 0u64;
    while falls.len() < trials && t < trials as u64 * 32 {
        t += 1;
        let scalar = seed.wrapping_add(t).saturating_mul(7).wrapping_add(1);
        let point = kc.mul(kc.generator(), &BigUint::from(scalar));
        let BinaryPoint::Affine { .. } = point else {
            continue;
        };
        let Some(sys) = build_chained_symmetrised_system(&kc, &fb, &point, m, &st) else {
            continue;
        };
        let (fall, macaulay) = first_fall_degree(&sys.equations, sys.n_vars, d_max);
        falls.push(fall);
        if let Some(q) = macaulay.iter().find(|q| q.degree == 2) {
            syz.push(q.syzygies() as f64);
        }
        if first.is_none() {
            first = Some((
                sys.n_vars,
                sys.equations.len(),
                system_degree(&sys.equations),
                fb.ell as u32,
            ));
        }
    }
    let (n_vars, n_eqs, degree, ell) = first?;
    let seen: Vec<u32> = falls.iter().flatten().copied().collect();
    Some(FfdSummary {
        n,
        ell,
        m,
        n_vars,
        n_eqs,
        degree,
        trials: falls.len(),
        fall_min: seen.iter().copied().min(),
        fall_max: seen.iter().copied().max(),
        no_fall: falls.iter().filter(|f| f.is_none()).count(),
        mean_syzygies_d2: if syz.is_empty() {
            0.0
        } else {
            syz.iter().sum::<f64>() / syz.len() as f64
        },
    })
}

fn main() {
    let raw: Vec<String> = std::env::args().skip(1).collect();
    let mut sym = false;
    let mut a: u8 = 0;
    let mut rest = Vec::new();
    let mut i = 0;
    while i < raw.len() {
        match raw[i].as_str() {
            "--sym" => sym = true,
            "--a" => {
                i += 1;
                a = raw.get(i).and_then(|s| s.parse().ok()).unwrap_or(0);
            }
            other => rest.push(other.to_string()),
        }
        i += 1;
    }
    let cases: Vec<(u32, usize)> = if rest.is_empty() {
        if sym {
            vec![(7, 4), (9, 4), (15, 4)]
        } else {
            vec![(9, 3), (9, 4), (15, 3), (15, 4)]
        }
    } else {
        parse_cases(&rest)
    };
    let trials: usize = std::env::var("FFD_TRIALS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(4);
    let d_max: u32 = std::env::var("FFD_DMAX")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(if sym { 4 } else { 3 });
    let seed = 0x5EED_u64;
    if sym {
        eprintln!(
            "frame=chained-sym a={a} d_max={d_max} trials={trials} seed={seed:#x} divisor=(n+1)/m"
        );
    } else {
        eprintln!("frame=chained-x factor_index=0 d_max={d_max} trials={trials} seed={seed:#x}");
    }
    let mut rows = Vec::new();
    for (n, m) in cases {
        let row = if sym {
            ffd_chained_sym(a, n, m, d_max, seed, trials)
        } else {
            ffd_summary(n, 0, m, d_max, seed, trials)
        };
        match row {
            Some(s) => rows.push(s),
            None => eprintln!("skip n={n} m={m}: no system under 64 vars"),
        }
    }
    print!("{}", format_ffd_table(&rows));
}

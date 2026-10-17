//! The entry points for the kernel formulas, chains, BMSS variants, duals, CSIDH, path finding,
//! KLPT and Richelot: every check passes on known inputs, and the methods agree with each other.
use isogeny_algos::api::{
    self, BmssMethod, ChainStrategy, PathMethod, PrimeCurve, Status, VeluMethod,
};
use isogeny_algos::int::Int;
use std::process::Command;

const P61: u64 = 2305843009213693951;

fn curve(p: u64, a: i64, b: i64) -> PrimeCurve {
    PrimeCurve::new(Int::from(p), Int::from(a), Int::from(b)).unwrap()
}

fn all_pass(cs: &[api::Check]) -> bool {
    cs.iter().all(|c| c.status == Status::Pass)
}

#[test]
fn kernel_formulas_agree_with_kohel() {
    // #E = 5 * 59 * 73 * ..., so E has a rational point of order 73
    let c = curve(P61, 7, 7);
    let ch = api::ell_power_chain(&c, 73, 1, None, None, ChainStrategy::Optimal, 1).unwrap();
    assert!(all_pass(&ch.checks));
    let x = &ch.kernel_point.0;
    let mut cods = vec![];
    for m in VeluMethod::ALL {
        let r = api::isogeny_from_x(&c, x, 73, m, 2).unwrap();
        assert!(all_pass(&r.checks), "{}", m.name());
        cods.push(r.codomain);
    }
    assert!(cods.windows(2).all(|w| w[0] == w[1]));
    assert_eq!(cods[0], ch.codomain);
}

#[test]
fn chains_by_every_strategy_reach_the_same_codomain() {
    // #E = 3^3 * 137 * ...
    let c = curve(P61, 4, 7);
    let n = api::count_points(&c, 1).unwrap().order;
    let mut cods = vec![];
    for st in [
        ChainStrategy::Naive,
        ChainStrategy::Balanced,
        ChainStrategy::Optimal,
    ] {
        let r = api::ell_power_chain(&c, 3, 3, None, Some(n.clone()), st, 5).unwrap();
        assert!(all_pass(&r.checks));
        assert_eq!(r.path.len(), 4);
        cods.push((r.kernel_point.clone(), r.codomain));
    }
    assert!(cods.windows(2).all(|w| w[0] == w[1]));
    // no point of order 2^5: the 2-part of the group is not cyclic enough
    assert!(api::ell_power_chain(
        &curve(P61, 1, 7),
        2,
        6,
        None,
        None,
        ChainStrategy::Optimal,
        1
    )
    .is_err());
}

#[test]
fn every_bmss_method_finds_the_same_kernels() {
    let c = curve(P61, 7, 7);
    let reference = api::isogenies(&c, 73, 1).unwrap();
    assert_eq!(reference.len(), 2);
    for m in BmssMethod::ALL {
        let rs = api::isogenies_by(&c, 73, m, 1).unwrap();
        assert_eq!(rs.len(), 2, "{}", m.name());
        for (r, s) in rs.iter().zip(&reference) {
            assert!(r.verified(), "{}", m.name());
            assert_eq!(r.kernel, s.kernel, "{}", m.name());
        }
    }
}

#[test]
fn duals_compose_to_multiplication() {
    let c = curve(P61, 7, 7);
    let ds = api::duals(&c, 73, 1).unwrap();
    assert_eq!(ds.len(), 2);
    for d in &ds {
        assert!(d.isogeny.verified() && all_pass(&d.checks));
        assert_eq!(d.dual_kernel.len(), 37);
    }
}

#[test]
fn csidh_action_small_and_512() {
    let primes = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 59];
    let exps = [1, -1, 2, 0, 0, 1, -3, 0, 1, 0, 0, 2];
    let r = api::csidh_action(&primes, &Int::zero(), &exps, 1).unwrap();
    assert!(all_pass(&r.checks));
    // the action is a group action: undoing it returns to A = 0
    let back: Vec<i32> = exps.iter().map(|e| -e).collect();
    let r2 = api::csidh_action(&primes, &r.a_out, &back, 2).unwrap();
    assert!(all_pass(&r2.checks));
    assert!(r2.a_out.is_zero());
    let ps = api::csidh512_primes();
    assert_eq!(ps.len(), 74);
    let mut e = vec![0i32; 74];
    e[0] = 1;
    e[73] = -1;
    let r = api::csidh_action(&ps, &Int::zero(), &e, 1).unwrap();
    assert!(all_pass(&r.checks));
    // a curve that is not supersingular is refused
    assert!(api::csidh_action(&primes, &Int::from(5i64), &exps, 1).is_err());
}

#[test]
fn path_methods_find_verified_paths() {
    for m in [
        PathMethod::Galbraith,
        PathMethod::Ghs,
        PathMethod::GalbraithStolbunov,
    ] {
        for seed in 1..4 {
            let r = api::isogeny_path(
                &Int::from(1000003u64),
                &Int::from(499965u64),
                &Int::from(620467u64),
                &[2, 3, 5, 7],
                m,
                100_000,
                seed,
            )
            .unwrap();
            assert!(all_pass(&r.checks), "{} seed {seed}", m.name());
            let mut seen = r.js.clone();
            seen.sort();
            seen.dedup();
            assert_eq!(seen.len(), r.js.len(), "no repeated vertex");
        }
    }
}

#[test]
fn klpt_prime_and_composite_norms() {
    let p = Int::from(1152921504606847067u64);
    for n in [Int::from(1000003u64), Int::from(1215u64)] {
        let r = api::klpt_random(&p, &n, 2, 1).unwrap();
        assert!(all_pass(&r.checks));
        assert_eq!(r.ideal_out.len(), 4);
    }
    assert!(api::klpt_random(&Int::from(1000033u64), &Int::from(7u64), 2, 1).is_err());
}

#[test]
fn richelot_small_and_large_fields() {
    let roots: Vec<Int> = (1..=6).map(|v| Int::from(v as u64)).collect();
    let small = api::richelot_isogenies(&Int::from(211u64), &roots, &Int::from(3u64), 1).unwrap();
    assert_eq!(small.len(), 15);
    assert!(small.iter().all(|r| all_pass(&r.checks)));
    let big = api::richelot_isogenies(&Int::from(P61), &roots, &Int::from(3u64), 1).unwrap();
    assert!(big
        .iter()
        .filter(|r| r.jacobian.is_some())
        .all(|r| all_pass(&r.checks)));
}

#[test]
fn cli_commands_exit_zero() {
    let bin = env!("CARGO_BIN_EXE_isogeny-algos");
    let p = P61.to_string();
    let runs: Vec<Vec<&str>> = vec![
        vec![
            "chain", "--p", &p, "--a", "4", "--b", "7", "--ell", "3", "--n", "3",
        ],
        vec![
            "isogenies",
            "--p",
            &p,
            "--a",
            "7",
            "--b",
            "7",
            "--ell",
            "73",
            "--method",
            "atkin1992",
        ],
        vec!["dual", "--p", &p, "--a", "7", "--b", "7", "--ell", "73"],
        vec![
            "csidh",
            "--primes",
            "3,5,7,11,13,17,19,23,29,31,37,59",
            "--exps",
            "1,-1,2,0,0,1,-3,0,1,0,0,2",
        ],
        vec![
            "path", "--p", "1000003", "--j1", "499965", "--j2", "620467", "--method", "ghs",
        ],
        vec!["klpt", "--p", "1152921504606847067", "--norm", "1000003"],
        vec![
            "richelot",
            "--p",
            "211",
            "--roots",
            "1,2,3,4,5,6",
            "--lc",
            "3",
        ],
    ];
    for args in runs {
        let out = Command::new(bin).args(&args).output().unwrap();
        let text = String::from_utf8_lossy(&out.stdout);
        assert_eq!(out.status.code(), Some(0), "{args:?}: {text}");
        assert!(text.contains("\"status\":\"PASS\""), "{args:?}");
    }
    let x = "2302858932999476951";
    for m in ["kohel", "velu", "xonly", "sqrtvelu"] {
        let out = Command::new(bin)
            .args([
                "velu",
                "--p",
                &p,
                "--a",
                "7",
                "--b",
                "7",
                "--ell",
                "73",
                "--kernel-x",
                x,
                "--method",
                m,
            ])
            .output()
            .unwrap();
        assert_eq!(out.status.code(), Some(0), "velu {m}");
    }
}

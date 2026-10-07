use isogeny_algos::curve::*;
use isogeny_algos::field::*;
use isogeny_algos::find::bmss::{self, Method};
use isogeny_algos::find::modpoly::Phi;
use isogeny_algos::kernel::kohel::kohel;
use isogeny_algos::testdata::*;

#[test]
fn bmss_family_agrees_with_kohel() {
    let mut rng = Rng::new(51);
    let fp = Zp::new(next_prime(1u64 << 30));
    let mut failures = vec![];
    for &ell in &[3u64, 5, 7, 11, 13, 17, 19, 23] {
        let phi = Phi::compute(&fp, ell as usize);
        for _trial in 0..3 {
            let (e, p) = loop {
                let (e, p) = curve_with_point(&fp, ell, &mut rng);
                let j = jinv(&fp, &e);
                if j != 0 && j != 1728 {
                    break (e, p);
                }
            };
            let (g, _) = kernel_poly_from_point(&fp, &e, &p, ell);
            let d = g.len() - 1;
            let sigma_true = if d >= 1 {
                fp.neg(fp.mul(2, g[d - 1]))
            } else {
                0
            };
            let kh = kohel(&fp, &e, &g, ell);
            let et = kh.cod;
            // sigma from Phi (needs the normalised codomain)
            let sphi = bmss::sigma_from_phi(&fp, &phi, &e, &et);
            if sphi != Some(sigma_true) {
                failures.push(format!(
                    "sigma_from_phi l={ell}: got {:?} want {}",
                    sphi, sigma_true
                ));
            }
            for m in Method::ALL {
                let sigma = if m.needs_sigma() {
                    Some(sigma_true)
                } else {
                    None
                };
                match bmss::isogeny(&fp, m, &e, &et, ell as usize, sigma) {
                    None => failures.push(format!("{} l={ell}: returned None", m.name())),
                    Some(iso) => {
                        if iso.ker != g || iso.num != kh.num {
                            failures.push(format!("{} l={ell}: wrong isogeny", m.name()));
                        }
                    }
                }
            }
        }
    }
    for f in &failures {
        eprintln!("FAIL {f}");
    }
    assert!(failures.is_empty(), "{} failures", failures.len());
}

#[test]
fn bmss_family_larger_l_no_phi() {
    let mut rng = Rng::new(52);
    for &(bits, ell) in &[
        (34u32, 29u64),
        (34, 31),
        (36, 37),
        (36, 41),
        (40, 53),
        (40, 61),
        (44, 101),
    ] {
        let fp = Zp::new(next_prime(1u64 << bits));
        for _trial in 0..2 {
            let (e, p) = loop {
                let (e, p) = curve_with_point(&fp, ell, &mut rng);
                let j = jinv(&fp, &e);
                if j != 0 && j != 1728 {
                    break (e, p);
                }
            };
            let (g, _) = kernel_poly_from_point(&fp, &e, &p, ell);
            let d = g.len() - 1;
            let sigma_true = fp.neg(fp.mul(2, g[d - 1]));
            let kh = kohel(&fp, &e, &g, ell);
            for m in Method::ALL {
                let sigma = if m.needs_sigma() {
                    Some(sigma_true)
                } else {
                    None
                };
                let iso = bmss::isogeny(&fp, m, &e, &kh.cod, ell as usize, sigma);
                assert!(iso.is_some(), "{} failed at l={ell}", m.name());
                let iso = iso.unwrap();
                assert!(
                    iso.ker == g && iso.num == kh.num,
                    "{} wrong at l={ell}",
                    m.name()
                );
            }
        }
    }
}

/// The algorithms must also work when the isogeny's kernel is not rational pointwise
/// (kernel polynomial irreducible), i.e. from the curves alone.
#[test]
fn bmss_family_on_galois_stable_kernels() {
    use isogeny_algos::find::divpoly;
    let mut rng = Rng::new(53);
    let fp = Zp::new(next_prime(1u64 << 30));
    let mut tested = 0;
    for _ in 0..300 {
        let c = Curve::new(fp.random(&mut rng), fp.random(&mut rng));
        if !is_smooth(&fp, &c) || jinv(&fp, &c) == 0 || jinv(&fp, &c) == 1728 {
            continue;
        }
        for ell in [5u64, 7] {
            for g in divpoly::kernel_polys(&fp, &c, ell, &mut rng) {
                let d = g.len() - 1;
                let kh = kohel(&fp, &c, &g, ell);
                let sigma = fp.neg(fp.mul(2, g[d - 1]));
                for m in Method::ALL {
                    let s = if m.needs_sigma() { Some(sigma) } else { None };
                    let iso = bmss::isogeny(&fp, m, &c, &kh.cod, ell as usize, s)
                        .unwrap_or_else(|| panic!("{} failed", m.name()));
                    assert_eq!(iso.ker, g);
                }
                tested += 1;
            }
        }
        if tested >= 12 {
            break;
        }
    }
    assert!(tested >= 6);
}

#[test]
fn end_to_end_from_phi_matches_v1() {
    use isogeny_algos::find::elkies;
    let mut rng = Rng::new(54);
    let fp = Zp::new(next_prime(1u64 << 30));
    for &ell in &[5usize, 7, 11, 13] {
        let phi = Phi::compute(&fp, ell);
        let mut checked = 0;
        while checked < 3 {
            let (e, _) = curve_with_point(&fp, ell as u64, &mut rng);
            let j = jinv(&fp, &e);
            if j == 0 || j == 1728 {
                continue;
            }
            let mut want: Vec<_> = elkies::elkies_isogenies(&fp, &phi, &e, &mut rng)
                .into_iter()
                .map(|i| i.ker)
                .collect();
            want.sort();
            assert!(!want.is_empty());
            for m in Method::ALL {
                let mut got: Vec<_> = bmss::isogenies_via_phi(&fp, &phi, &e, m, &mut rng)
                    .into_iter()
                    .map(|i| i.ker)
                    .collect();
                got.sort();
                assert_eq!(got, want, "{} l={ell}", m.name());
            }
            checked += 1;
        }
    }
}

use isogeny_algos::bigint::Big;
use isogeny_algos::curve::*;
use isogeny_algos::field::*;
use isogeny_algos::find::sea::*;
use isogeny_algos::fpm::FpM;
use std::collections::HashMap;

#[test]
fn atkin_candidate_sets_partition() {
    // every t mod l with a non-zero ratio falls in exactly one order class
    for l in [3u64, 5, 7, 11, 13] {
        for q in 1..l {
            let mut total = 0;
            for r in 1..=(l as usize + 1) {
                total += atkin_candidates(l, q, r).len();
            }
            assert!(total <= l as usize);
        }
    }
}

#[test]
fn sea_matches_bsgs_at_40_and_61_bits() {
    let mut rng = Rng::new(1100);
    for bits in [40u32, 61] {
        let fp = Zp::new(next_prime((1u64 << (bits - 1)) + 12345));
        let mut phis = HashMap::new();
        for _ in 0..4 {
            let e = loop {
                let e = Curve::new(fp.random(&mut rng), fp.random(&mut rng));
                let j = jinv(&fp, &e);
                if is_smooth(&fp, &e) && j != 0 && j != 1728 {
                    break e;
                }
            };
            let (n, st) = sea(&fp, &e, 43, &mut phis, &mut rng).expect("SEA");
            let reference = order(&fp, &e, &mut rng);
            assert_eq!(n, Big::from_u64(reference), "bits = {bits}, stats {st:?}");
        }
    }
}

#[test]
fn sea_at_128_bits() {
    let mut rng = Rng::new(1101);
    // 2^127 - 1
    let f = FpM::<2>::from_dec("170141183460469231731687303715884105727");
    let mut phis = HashMap::new();
    let e = loop {
        let e = Curve::new(f.random(&mut rng), f.random(&mut rng));
        let j = jinv(&f, &e);
        if is_smooth(&f, &e) && j != f.zero() && j != f.from_u64(1728) {
            break e;
        }
    };
    let (n, st) = sea(&f, &e, 61, &mut phis, &mut rng).expect("SEA");
    for _ in 0..5 {
        let p = random_point_f(&f, &e, &mut rng);
        assert_eq!(pmul_big(&f, &e, &p, &n), Pt::Inf, "stats {st:?}");
    }
    eprintln!("128-bit SEA: elkies {:?}, atkin {:?}, candidates {}", st.elkies, st.atkin, st.candidates);
}

/// Isogeny cycles: the eigenvalue mod l^k on the rational cyclic l^k-subgroup gives t mod l^k,
/// checked against the BSGS group order.
#[test]
fn isogeny_cycles_give_t_mod_l_power() {
    use isogeny_algos::find::bmss;
    use isogeny_algos::find::elkies::elkies_codomain;
    use isogeny_algos::find::modpoly::Phi;
    let mut rng = Rng::new(1200);
    // q = 2 mod 3: for q = 1 mod 3 the eigenvalues mod 3 (product q) cannot be distinct
    let mut p0 = next_prime((1u64 << 39) + 777);
    while p0 % 3 != 2 {
        p0 = next_prime(p0 + 1);
    }
    let fp = Zp::new(p0);
    let q = fp.p as i128;
    let mut reached = std::collections::HashMap::new();
    for &(ell, kmax) in &[(3u64, 4u32), (5, 3), (7, 2), (11, 2)] {
        let phi = Phi::compute(&fp, ell as usize);
        let mut done = 0;
        let mut tries = 0;
        while done < 3 && tries < 200 {
            tries += 1;
            let e = Curve::new(fp.random(&mut rng), fp.random(&mut rng));
            let j = jinv(&fp, &e);
            if !is_smooth(&fp, &e) || j == 0 || j == 1728 {
                continue;
            }
            let roots = phi.neighbors(&fp, j, &mut rng);
            if roots.len() != 2 {
                continue; // Elkies prime with two rational l-isogenies (distinct eigenvalues)
            }
            let Some(et) = elkies_codomain(&fp, &phi, &e, roots[0]) else { continue };
            let Some(iso) = bmss::isogeny(&fp, bmss::Method::FastElkiesPrime, &e, &et, ell as usize, None) else { continue };
            let Some(lam1) = elkies_eigenvalue(&fp, &e, &iso.ker, ell) else { continue };
            let (lam, k) = eigenvalue_cycle(&fp, &phi, &e, &iso, ell, lam1, kmax, &mut rng);
            let m = (ell as i128).pow(k);
            let t = q + 1 - order(&fp, &e, &mut rng) as i128;
            // t = lam + q / lam mod l^k
            let lam_inv = (1..m).find(|&v| (v * lam as i128) % m == 1).unwrap();
            let tk = (lam as i128 + (q % m) * lam_inv) % m;
            assert_eq!(tk, t.rem_euclid(m), "l = {ell}, k = {k}");
            *reached.entry(ell).or_insert(0) = (*reached.get(&ell).unwrap_or(&0)).max(k);
            done += 1;
        }
    }
    eprintln!("max k reached per l: {reached:?}");
    assert!(reached.get(&3).copied().unwrap_or(0) >= 3);
    assert!(reached.get(&5).copied().unwrap_or(0) >= 2);
}

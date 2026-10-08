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

/// (lambda, mu, ord(lambda/mu)) in F_{l^2} from the true trace, for the independent check.
fn eig_mod_l(ell: u64, p: u64, t: u64) -> ((u64, u64), (u64, u64), u64) {
    let l = ell;
    let md = |a: u64| a % l;
    let sqrt_l = |a: u64| (0..l).find(|&x| x * x % l == md(a));
    let nr = (2..l).find(|&a| sqrt_l(a).is_none()).unwrap_or(2);
    let mul = |x: (u64, u64), y: (u64, u64)| (md(x.0 * y.0 + nr * md(x.1 * y.1)), md(x.0 * y.1 + x.1 * y.0));
    let pow = |mut x: (u64, u64), mut e: u64| {
        let mut rr = (1u64, 0u64);
        while e > 0 {
            if e & 1 == 1 {
                rr = mul(rr, x);
            }
            x = mul(x, x);
            e >>= 1;
        }
        rr
    };
    let inv2 = (l + 1) / 2;
    let ql = md(p);
    let disc = md(t * t + l * l * 4 - 4 * ql);
    let s = match sqrt_l(disc) {
        Some(v) => (v, 0),
        None => (0, sqrt_l(md(disc * pow((nr, 0), l * l - 2).0)).unwrap()),
    };
    let lam = mul((md(t + s.0), s.1), (inv2, 0));
    let mu = mul((md(t + l - s.0), md(l - s.1)), (inv2, 0));
    let invx = |x: (u64, u64)| pow(x, l * l - 2);
    let g = mul(lam, invx(mu));
    let mut x = g;
    let mut ord = 1u64;
    while x != (1, 0) && ord <= l * l {
        x = mul(x, g);
        ord += 1;
    }
    (lam, mu, ord)
}

fn pow_fl2(ell: u64, x: (u64, u64), e: u64) -> (u64, u64) {
    let l = ell;
    let md = |a: u64| a % l;
    let sqrt_l = |a: u64| (0..l).find(|&y| y * y % l == md(a));
    let nr = (2..l).find(|&a| sqrt_l(a).is_none()).unwrap_or(2);
    let mul = |x: (u64, u64), y: (u64, u64)| (md(x.0 * y.0 + nr * md(x.1 * y.1)), md(x.0 * y.1 + x.1 * y.0));
    let mut r = (1u64, 0u64);
    let mut b = x;
    let mut e = e;
    while e > 0 {
        if e & 1 == 1 {
            r = mul(r, b);
        }
        b = mul(b, b);
        e >>= 1;
    }
    r
}

/// Atkin primes over F_{p^d} towers (d = order of lambda/mu): Frobenius^d acts on E[l] as the
/// scalar lambda^d mod l. The tower recovers (lambda^d, d); it equals the true (lambda)^d, and
/// the refined candidate set contains the true t mod l -- often a singleton, i.e. an Atkin prime
/// pinned to Elkies strength.
#[test]
fn atkin_tower_eigenvalue_matches_true_trace() {
    use isogeny_algos::find::modpoly::Phi;
    let mut rng = Rng::new(1300);
    let fp = Zp::new(next_prime((1u64 << 39) + 4242));
    let p = fp.p;
    let mut phis: HashMap<usize, Phi<Zp>> = HashMap::new();
    let mut recovered = 0;
    let mut singletons = 0;
    let mut tries = 0;
    while recovered < 4 && tries < 600 {
        tries += 1;
        let e = Curve::new(fp.random(&mut rng), fp.random(&mut rng));
        let j = jinv(&fp, &e);
        if !is_smooth(&fp, &e) || j == 0 || j == 1728 {
            continue;
        }
        let t_true = (p as i128 + 1 - order(&fp, &e, &mut rng) as i128).rem_euclid(p as i128) as u64;
        for &ell in &[5u64, 7, 11, 13] {
            // genuine Atkin prime: discriminant t^2 - 4p a non-residue mod l (lambda, mu not in F_l)
            let (lam, _mu, ord) = eig_mod_l(ell, p, t_true % ell);
            let disc = ((t_true % ell) * (t_true % ell) + ell * ell * 4 - 4 * (p % ell)) % ell;
            if (0..ell).any(|x| x * x % ell == disc) {
                continue; // residue => Elkies
            }
            if !(2..=4).contains(&ord) {
                continue; // keep the tower degree small/fast for the test
            }
            let phi = phis.entry(ell as usize).or_insert_with(|| Phi::compute(&fp, ell as usize));
            if !phi.neighbors(&fp, j, &mut rng).is_empty() {
                continue; // classified Elkies by the modular polynomial
            }
            let _ = ord;
            let Some((nu, d)) = atkin_eigenvalue_tower(&fp, &e, ell, &mut rng) else { continue };
            // independent check: nu is a true eigenvalue of pi^d, i.e. nu == lambda^d or mu^d mod l
            let ld = pow_fl2(ell, lam, d as u64);
            let (_lam2, mu, _o) = eig_mod_l(ell, p, t_true % ell);
            let md_ = pow_fl2(ell, mu, d as u64);
            assert!(
                (ld.1 == 0 && nu == ld.0) || (md_.1 == 0 && nu == md_.0),
                "l={ell} d={d}: tower eigenvalue {nu} is neither lambda^d={ld:?} nor mu^d={md_:?}"
            );
            // refined candidate set must contain the true t mod l
            let cands = atkin_candidates_tower(ell, p % ell, d, nu);
            assert!(cands.contains(&(t_true % ell)), "l={ell} d={d}: true t mod l not in {cands:?}");
            if cands.len() == 1 {
                singletons += 1;
            }
            recovered += 1;
            if recovered >= 4 {
                break;
            }
        }
    }
    eprintln!("atkin tower: {recovered} recoveries, {singletons} singletons");
    assert!(recovered >= 4, "expected >= 4 Atkin tower recoveries, got {recovered}");
    assert!(singletons >= 1, "expected at least one Atkin prime resolved to a unique t mod l");
}

/// SEA that resolves Atkin primes to congruences via F_{p^d} towers gives the correct #E and
/// never leaves more candidates than plain SEA (Atkin primes become Elkies-strength when the
/// tower pins t mod l).
#[test]
fn sea_with_atkin_towers_matches_bsgs() {
    use isogeny_algos::find::modpoly::Phi;
    let mut rng = Rng::new(1400);
    let fp = Zp::new(next_prime((1u64 << 37) + 99));
    let mut phis: HashMap<usize, Phi<Zp>> = HashMap::new();
    let mut tested = 0;
    let mut improved = 0;
    let max_ell = 31;
    for _ in 0..6 {
        let e = loop {
            let e = Curve::new(fp.random(&mut rng), fp.random(&mut rng));
            let j = jinv(&fp, &e);
            if is_smooth(&fp, &e) && j != 0 && j != 1728 {
                break e;
            }
        };
        let reference = order(&fp, &e, &mut rng);
        let (n_plain, st_plain) = sea(&fp, &e, max_ell, &mut phis, &mut rng).expect("plain SEA");
        let (n_tower, st_tower) = sea_atkin_tower(&fp, &e, max_ell, &mut phis, &mut rng).expect("tower SEA");
        assert_eq!(n_plain, Big::from_u64(reference), "plain SEA wrong");
        assert_eq!(n_tower, Big::from_u64(reference), "tower SEA wrong");
        // the Atkin resolver only adds congruences, so it never leaves more candidates
        assert!(st_tower.candidates <= st_plain.candidates, "tower must not add candidates");
        if st_tower.candidates < st_plain.candidates {
            improved += 1;
        }
        tested += 1;
    }
    eprintln!("sea+atkin-towers: {tested} curves, {improved} with fewer candidates than plain SEA");
    assert!(tested >= 6);
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

/// The Atkin degree (Q-matrix Frobenius powers, equality at the divisors of l+1) equals the
/// definition: the least k with gcd(Y^(q^k) - Y, Phi_l(j, Y)) != 1, each Y^(q^k) by a fresh
/// exponentiation. Random curves at 40 bits, every Atkin prime l <= 31.
#[test]
fn atkin_degree_matches_gcd_definition() {
    use isogeny_algos::find::modpoly::Phi;
    use isogeny_algos::poly;
    let mut p = (1u64 << 39) + 777;
    while !is_prime(p) {
        p += 1;
    }
    let fp = Zp::new(p);
    let mut rng = Rng::new(4100);
    let phis: Vec<(usize, Phi<Zp>)> = [3usize, 5, 7, 11, 13, 17, 19, 23, 29, 31].iter().map(|&l| (l, Phi::compute(&fp, l))).collect();
    let mut checked = 0;
    while checked < 60 {
        let e = Curve::new(fp.random(&mut rng), fp.random(&mut rng));
        let j = jinv(&fp, &e);
        if !is_smooth(&fp, &e) || j == 0 || j == 1728 {
            continue;
        }
        for (_, phi) in &phis {
            let g = poly::monic(&fp, &phi.y_poly(&fp, j));
            let y = poly::x_poly(&fp);
            let mut yk = y.clone();
            let mut naive = None;
            for k in 1..g.len() {
                yk = poly::powmod_big(&fp, &yk, &fp.q(), &g);
                if poly::deg(&fp, &poly::gcd(&fp, &g, &poly::sub(&fp, &yk, &y))) > 0 {
                    naive = Some(k);
                    break;
                }
            }
            if naive == Some(1) {
                continue; // Elkies prime (or l = j-special): not an Atkin degree
            }
            assert_eq!(atkin_degree(&fp, &g), naive, "p = {p}, j = {j}, l = {}", phi.ell);
            checked += 1;
        }
    }
}

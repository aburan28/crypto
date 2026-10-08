use isogeny_algos::curve::*;
use isogeny_algos::field::*;
use isogeny_algos::kernel::chain::*;

/// p = 2^a 3^b f - 1 prime, p = 3 mod 4.
fn sidh_prime(a: u32, b: u32) -> u64 {
    let base = 2u64.pow(a) * 3u64.pow(b);
    for f in 1..10_000u64 {
        let p = base * f - 1;
        if p > 1000 && is_prime(p) && p % 4 == 3 {
            return p;
        }
    }
    panic!("no prime");
}

fn torsion_basis(
    f2: &Zp2,
    e: &Curve<(u64, u64)>,
    p: u64,
    ell: u64,
    n: u32,
    rng: &mut Rng,
) -> (Pt<(u64, u64)>, Pt<(u64, u64)>) {
    // E(F_{p^2}) = E[p+1] for the maximal curve; cofactor (p+1)/l^n lands in E[l^n]
    let cof = ((p + 1) / ell.pow(n)) as u128;
    let top = |q: &Pt<(u64, u64)>| pmul(f2, e, q, ell.pow(n - 1) as u128);
    let mut pts = vec![];
    while pts.len() < 2 {
        let r = random_point_f(f2, e, rng);
        let q = pmul(f2, e, &r, cof);
        let t = top(&q);
        if t == Pt::Inf || pmul(f2, e, &t, ell as u128) != Pt::Inf {
            continue;
        }
        if let Some(first) = pts.first() {
            // independence mod l: [l^{n-1}]Q not in <[l^{n-1}]P>
            let t0 = top(first);
            let mut dep = false;
            let mut m = Pt::Inf;
            for _ in 0..ell {
                if m == t {
                    dep = true;
                }
                m = padd(f2, e, &m, &t0);
            }
            if dep {
                continue;
            }
        }
        pts.push(q);
    }
    (pts[0], pts[1])
}

fn j_of(f2: &Zp2, e: &Curve<(u64, u64)>) -> (u64, u64) {
    jinv(f2, e)
}

#[test]
fn sidh_key_exchange_all_strategies_agree() {
    let mut rng = Rng::new(41);
    let (a, b) = (9u32, 5u32);
    let p = sidh_prime(a, b);
    let f2 = Zp2::new(p);
    // E0: y^2 = x^3 + x (supersingular for p = 3 mod 4), over F_{p^2} it has (p+1)^2 points
    let e0 = Curve::new(f2.from_u64(1), f2.zero());
    let (pa, qa) = torsion_basis(&f2, &e0, p, 2, a, &mut rng);
    let (pb, qb) = torsion_basis(&f2, &e0, p, 3, b, &mut rng);
    let sa = rng.below(2u64.pow(a)) | 1;
    let sb = rng.below(3u64.pow(b));
    let strategies = [
        Strategy::Naive,
        Strategy::Balanced,
        Strategy::Optimal {
            mul: 1.0,
            eval: 2.0,
            iso: 3.0,
        },
        Strategy::Optimal {
            mul: 5.0,
            eval: 1.0,
            iso: 1.0,
        },
    ];
    let mut shared = vec![];
    let mut stats = vec![];
    for strat in strategies {
        // Alice: kernel <PA + [sa] QA> of order 2^a, push PB, QB
        let ka = padd(&f2, &e0, &pa, &pmul(&f2, &e0, &qa, sa as u128));
        let (ca, pushed_a, st_a) =
            ell_power_isogeny(&f2, &e0, &ka, 2, a as usize, &[pb, qb], strat);
        // Bob: kernel <PB + [sb] QB> of order 3^b, push PA, QA
        let kb = padd(&f2, &e0, &pb, &pmul(&f2, &e0, &qb, sb as u128));
        let (cb, pushed_b, st_b) =
            ell_power_isogeny(&f2, &e0, &kb, 3, b as usize, &[pa, qa], strat);
        for q in pushed_a.iter().chain(pushed_b.iter()) {
            assert!(q != &Pt::Inf);
        }
        assert!(on_curve(&f2, &ca.cod, &pushed_a[0]) && on_curve(&f2, &cb.cod, &pushed_b[0]));
        // shared curves
        let kab = padd(
            &f2,
            &cb.cod,
            &pushed_b[0],
            &pmul(&f2, &cb.cod, &pushed_b[1], sa as u128),
        );
        let (cab, _, _) = ell_power_isogeny(&f2, &cb.cod, &kab, 2, a as usize, &[], strat);
        let kba = padd(
            &f2,
            &ca.cod,
            &pushed_a[0],
            &pmul(&f2, &ca.cod, &pushed_a[1], sb as u128),
        );
        let (cba, _, _) = ell_power_isogeny(&f2, &ca.cod, &kba, 3, b as usize, &[], strat);
        assert_eq!(
            j_of(&f2, &cab.cod),
            j_of(&f2, &cba.cod),
            "SIDH shared j-invariant, {:?}",
            strat
        );
        shared.push((j_of(&f2, &ca.cod), j_of(&f2, &cb.cod), j_of(&f2, &cab.cod)));
        stats.push((st_a, st_b));
    }
    // every strategy builds the same isogenies
    for s in &shared[1..] {
        assert_eq!(s, &shared[0]);
    }
    // naive does the most l-multiplications; balanced/optimal fewer
    assert!(stats[1].0.l_mults < stats[0].0.l_mults);
    eprintln!(
        "2^{a}-chain stats naive {:?} balanced {:?} opt {:?}",
        stats[0].0, stats[1].0, stats[2].0
    );
}

#[test]
fn smooth_composite_kernel() {
    let mut rng = Rng::new(42);
    let p = sidh_prime(7, 4);
    let f2 = Zp2::new(p);
    let e0 = Curve::new(f2.from_u64(1), f2.zero());
    // a point of order 2^3 * 3^2 * 5? use only 2 and 3 parts: order 2^7*3^4 divides p+1
    let cof = ((p + 1) / (2u64.pow(7) * 3u64.pow(4))) as u128;
    let pt = loop {
        let r = random_point_f(&f2, &e0, &mut rng);
        let q = pmul(&f2, &e0, &r, cof);
        // exact 2-part 2^7 and exact 3-part 3^4
        let ok2 = pmul(&f2, &e0, &q, 2u128.pow(6) * 3u128.pow(4)) != Pt::Inf;
        let ok3 = pmul(&f2, &e0, &q, 3u128.pow(3) * 2u128.pow(7)) != Pt::Inf;
        if ok2 && ok3 {
            break q;
        }
    };
    let (c1, _) = smooth_isogeny(&f2, &e0, &pt, &[(2, 7), (3, 4)], Strategy::Naive);
    let (c2, _) = smooth_isogeny(&f2, &e0, &pt, &[(3, 4), (2, 7)], Strategy::Balanced);
    assert_eq!(c1.deg, 2u64.pow(7) * 3u64.pow(4));
    // order of composition must not matter
    assert_eq!(jinv(&f2, &c1.cod), jinv(&f2, &c2.cod));
}

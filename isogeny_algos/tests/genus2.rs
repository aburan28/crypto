use isogeny_algos::field::*;
use isogeny_algos::genus2::*;
use isogeny_algos::poly;

fn distinct(fp: &Zp, k: usize, rng: &mut Rng) -> Vec<u64> {
    let mut v: Vec<u64> = vec![];
    while v.len() < k {
        let x = fp.random(rng);
        if !v.contains(&x) {
            v.push(x);
        }
    }
    v
}

#[test]
fn richelot_preserves_l_polynomial() {
    let mut rng = Rng::new(1000);
    for p in [101u64, 211, 503] {
        let fp = Zp::new(p);
        for _ in 0..3 {
            let r = distinct(&fp, 6, &mut rng);
            let lc = fp.random(&mut rng) | 1;
            let f = poly::scale(&fp, &poly::from_roots(&fp, &r), lc);
            let l0 = lpoly(&fp, &f);
            let mut checked = 0;
            for (a, b, c) in matchings6() {
                let g = [
                    poly::scale(&fp, &poly::from_roots(&fp, &[r[a.0], r[a.1]]), lc),
                    poly::from_roots(&fp, &[r[b.0], r[b.1]]),
                    poly::from_roots(&fp, &[r[c.0], r[c.1]]),
                ];
                let rl = richelot(&fp, &g);
                if rl.delta != 0 {
                    let fc = rl.codomain(&fp);
                    assert!(fc.len() == 7 || fc.len() == 6, "quintic or sextic codomain");
                    assert_eq!(lpoly(&fp, &fc), l0, "p = {p}");
                    checked += 1;
                }
            }
            assert!(checked >= 10);
        }
    }
}

#[test]
fn degenerate_richelot_splits() {
    let mut rng = Rng::new(1001);
    for p in [101u64, 211, 503] {
        let fp = Zp::new(p);
        for _ in 0..4 {
            // F(z) = c prod (z^2 - rho_i^2): the matching {rho_i, -rho_i} has Delta = 0
            let rho = distinct(&fp, 3, &mut rng);
            if rho.contains(&0)
                || rho
                    .iter()
                    .enumerate()
                    .any(|(i, &x)| rho.iter().skip(i + 1).any(|&y| y == fp.neg(x)))
            {
                continue;
            }
            let c = fp.random(&mut rng) | 1;
            let roots = [
                rho[0],
                fp.neg(rho[0]),
                rho[1],
                fp.neg(rho[1]),
                rho[2],
                fp.neg(rho[2]),
            ];
            let f = poly::scale(&fp, &poly::from_roots(&fp, &roots), c);
            // then move it by a Moebius transformation x -> (2x + 3)/(x + 5) to hide the symmetry
            let (ma, mb, mc, md) = (2u64, 3u64, 1u64, 5u64);
            let mut g: Vec<u64> = vec![0; 7];
            for k in 0..7 {
                // sum f_k (a x + b)^k (c x + d)^(6-k)
                let mut t = vec![f[k]];
                for _ in 0..k {
                    t = poly::mul(&fp, &t, &vec![mb, ma]);
                }
                for _ in k..6 {
                    t = poly::mul(&fp, &t, &vec![md, mc]);
                }
                let s = poly::add(&fp, &g, &t);
                g = s;
                g.resize(7, 0);
            }
            if g[6] == 0 {
                continue;
            }
            // roots move by the inverse map
            let inv = |y: u64| fp.div(fp.sub(fp.mul(md, y), mb), fp.sub(ma, fp.mul(mc, y)));
            let nr: Vec<u64> = roots.iter().map(|&y| inv(y)).collect();
            let lcg = g[6];
            let gs = [
                poly::scale(&fp, &poly::from_roots(&fp, &[nr[0], nr[1]]), lcg),
                poly::from_roots(&fp, &[nr[2], nr[3]]),
                poly::from_roots(&fp, &[nr[4], nr[5]]),
            ];
            assert_eq!(
                poly::scale(
                    &fp,
                    &poly::mul(&fp, &poly::mul(&fp, &gs[0], &gs[1]), &gs[2]),
                    1
                ),
                g
            );
            let rl = richelot(&fp, &gs);
            assert_eq!(rl.delta, 0, "symmetric matching is degenerate");
            let Some((e1, e2)) = split(&fp, &g, &gs) else {
                continue;
            };
            let (a1, a2) = (elliptic_trace(&fp, &e1), elliptic_trace(&fp, &e2));
            assert_eq!(lpoly(&fp, &g), LPoly::product(a1, a2, p as i64), "p = {p}");
        }
    }
}

#[test]
fn gluing_matches_product_l_polynomial() {
    let mut rng = Rng::new(1002);
    for p in [101u64, 103, 211, 499, 503, 1009, 1019] {
        let fp = Zp::new(p);
        let mut done = 0;
        for _ in 0..20 {
            let a = distinct(&fp, 3, &mut rng);
            let b = distinct(&fp, 3, &mut rng);
            let Some(c) = glue(&fp, [a[0], a[1], a[2]], [b[0], b[1], b[2]]) else {
                continue;
            };
            if c.len() != 7 {
                continue;
            }
            let ta = elliptic_trace(&fp, &poly::from_roots(&fp, &a));
            let tb = elliptic_trace(&fp, &poly::from_roots(&fp, &b));
            let l = lpoly(&fp, &c);
            assert_eq!(
                l,
                LPoly::product(ta, tb, p as i64),
                "p = {p}: glued Jacobian ~ E1 x E2"
            );
            done += 1;
        }
        assert!(done >= 5);
    }
}

/// Ibukiyama–Katsura–Oort: number of superspecial genus-2 curves over F_p-bar, times 2880:
/// (p-1)(p^2+25p+166) - 90{1-(-1/p)} + 360{1-(-2/p)} + 160{1-(-3/p)} + (2304 if p = 4 mod 5).
fn iko_count(p: i64) -> i64 {
    let m1 = if p % 4 == 1 { 1 } else { -1 };
    let m2 = if p % 8 == 1 || p % 8 == 3 { 1 } else { -1 };
    let m3 = if p % 3 == 1 { 1 } else { -1 };
    let t = (p - 1) * (p * p + 25 * p + 166) - 90 * (1 - m1)
        + 360 * (1 - m2)
        + 160 * (1 - m3)
        + if p % 5 == 4 { 2304 } else { 0 };
    assert_eq!(t % 2880, 0);
    t / 2880
}

#[test]
fn superspecial_richelot_graph() {
    let mut rng = Rng::new(1003);
    for p in [11u64, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83] {
        let g = superspecial_graph(p, 100_000, &mut rng);
        // supersingular j over F_{p^2}: h = floor(p/12) + {0,1,1,2}
        let h = (p / 12 + [0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 2][(p % 12) as usize]) as usize;
        assert_eq!(g.products, h * (h + 1) / 2, "p = {p}: products E1 x E2");
        assert_eq!(
            g.jacobians as i64,
            iko_count(p as i64),
            "p = {p}: superspecial Jacobians (Ibukiyama-Katsura-Oort)"
        );
        eprintln!(
            "p = {p}: {} Jacobians, {} products, {} edges found",
            g.jacobians, g.products, g.edges
        );
    }
}

#[test]
fn richelot_correspondence_maps_points_to_codomain() {
    let mut rng = Rng::new(1004);
    for p in [103u64, 211, 503] {
        let fp = Zp::new(p);
        let f2 = Zp2::new(p);
        let r = distinct(&fp, 6, &mut rng);
        let lc = fp.random(&mut rng) | 1;
        let f = poly::scale(&fp, &poly::from_roots(&fp, &r), lc);
        let mut images = 0;
        for (a, b, c) in matchings6().into_iter().take(5) {
            let g = [
                poly::scale(&fp, &poly::from_roots(&fp, &[r[a.0], r[a.1]]), lc),
                poly::from_roots(&fp, &[r[b.0], r[b.1]]),
                poly::from_roots(&fp, &[r[c.0], r[c.1]]),
            ];
            let rl = richelot(&fp, &g);
            if rl.delta == 0 {
                continue;
            }
            let cod: Vec<(u64, u64)> = rl.codomain(&fp).iter().map(|&c| (c, 0)).collect();
            for _ in 0..20 {
                let x = fp.random(&mut rng);
                let Some(y) = fp.sqrt(poly::eval(&fp, &f, x)) else {
                    continue;
                };
                if y == 0 {
                    continue;
                }
                let img = rl.image_point(&fp, &g, &f2, |c| (c, 0), x, y, &mut rng);
                assert!(img.len() <= 2);
                for (xp, yp) in img {
                    assert_eq!(
                        f2.mul(yp, yp),
                        poly::eval(&f2, &cod, xp),
                        "image on the codomain, p = {p}"
                    );
                    images += 1;
                }
            }
        }
        assert!(images > 20);
    }
}

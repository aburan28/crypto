//! Integer Coppersmith (Coron bivariate, trivariate with resultant elimination) on planted
//! norm-equation instances.
use isogeny_algos::field::Rng;
use isogeny_algos::fixed_degree::coppersmith::*;
use isogeny_algos::int::Int;
use isogeny_algos::mpoly::MPoly;

fn int(v: i64) -> Int {
    Int::from(v)
}

fn rnd_signed(rng: &mut Rng, bound: i64) -> Int {
    int((rng.next() % (2 * bound as u64 + 1)) as i64 - bound)
}

/// Random symmetric positive-definite integral Gram matrix with entries of about `bits` bits:
/// `G = B^T B` with `B` having entries of about `bits / 2` bits plus a diagonal shift.
fn random_gram(rng: &mut Rng, bits: u32) -> [[Int; 4]; 4] {
    let h = bits / 2;
    let b: Vec<Vec<Int>> = (0..4)
        .map(|_| {
            (0..4)
                .map(|_| {
                    let v = Int::from(rng.next() >> (64 - h));
                    if rng.next() & 1 == 1 {
                        -v
                    } else {
                        v
                    }
                })
                .collect()
        })
        .collect();
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            let mut s = Int::zero();
            for k in 0..4 {
                s = &s + &(&b[k][i] * &b[k][j]);
            }
            if i == j {
                s = &s + &Int::from(1u64 << h);
            }
            s
        })
    })
}

fn planted(rng: &mut Rng, xb: i64) -> [Int; 4] {
    std::array::from_fn(|_| rnd_signed(rng, xb))
}

#[test]
fn bivariate_planted() {
    let mut rng = Rng::new(1);
    // (Gram bits, |x_i| bound, k, shape), all inside the measured success region
    for (bits, xb, k, shape) in [
        (20u32, 12i64, 2u32, Shape::Box),
        (30, 30, 2, Shape::Box),
        (40, 60, 2, Shape::Box),
        (20, 12, 2, Shape::Simplex),
        (30, 60, 2, Shape::Simplex),
        (40, 500, 3, Shape::Simplex),
    ] {
        for _ in 0..2 {
            let g = random_gram(&mut rng, bits);
            let x = planted(&mut rng, xb);
            let d = quad_form(&g, &x);
            let f = norm_eq_poly(&g, &d, &x[2..]);
            let params = CoppersmithParams::new(k, shape);
            let xbi = int(xb + 1);
            let (roots, rep) = coron_bivariate_with(&f, &xbi, &xbi, &params);
            for (a, b) in &roots {
                assert!(f.eval(&[a.clone(), b.clone()]).is_zero());
            }
            assert!(
                roots.contains(&(x[0].clone(), x[1].clone())),
                "bits={bits} X={xb} k={k} {shape:?}: missing planted root; {rep:?}"
            );
            assert!(rep.complete);
        }
    }
    // the plain entry point (Coron box lattice)
    let g = random_gram(&mut rng, 24);
    let x = planted(&mut rng, 20);
    let d = quad_form(&g, &x);
    let f = norm_eq_poly(&g, &d, &x[2..]);
    let r = coron_bivariate(&f, &int(21), &int(21), 2);
    assert!(r.contains(&(x[0].clone(), x[1].clone())));
}

#[test]
fn bivariate_zero_constant_and_generic() {
    // p(0,0) = 0 forces the shift path; generic (non-norm) polynomial with degree 2 in each var
    let mut rng = Rng::new(2);
    for _ in 0..2 {
        let (x0, y0) = (rnd_signed(&mut rng, 15), rnd_signed(&mut rng, 15));
        let mut f = MPoly::zero(2);
        for e in [
            [2u32, 2],
            [2, 1],
            [1, 2],
            [2, 0],
            [0, 2],
            [1, 1],
            [1, 0],
            [0, 1],
        ] {
            f.add_term(&e, &Int::from(rng.next() >> 34));
        }
        let v = f.eval(&[x0.clone(), y0.clone()]);
        f = &f - &MPoly::constant(2, &v);
        let (r, rep) = coron_bivariate_with(
            &f,
            &int(16),
            &int(16),
            &CoppersmithParams::new(3, Shape::Box),
        );
        assert!(r.contains(&(x0.clone(), y0.clone())), "{r:?} {rep:?}");
        for (a, b) in &r {
            assert!(f.eval(&[a.clone(), b.clone()]).is_zero());
        }
    }
    // f = A xy + B x + C y + D x^2 vanishes at (0,0): exercises the shift path (p(0) = 0).
    // (Coron needs XY < W^(1/3-eps), so the coefficients must be large compared with X, Y.)
    let x = MPoly::var(2, 0);
    let y = MPoly::var(2, 1);
    let co: Vec<Int> = (0..4).map(|_| Int::from(rng.next() >> 30)).collect();
    let f = &(&(&(&x * &y).scale(&co[0]) + &x.scale(&co[1])) + &y.scale(&co[2]))
        + &(&x * &x).scale(&co[3]);
    let (r, rep) =
        coron_bivariate_with(&f, &int(6), &int(6), &CoppersmithParams::new(2, Shape::Box));
    assert!(rep.complete, "{rep:?}");
    for (a, b) in &r {
        assert!(f.eval(&[a.clone(), b.clone()]).is_zero());
    }
    let mut brute = vec![];
    for a in -5i64..=5 {
        for b in -5i64..=5 {
            if f.eval(&[int(a), int(b)]).is_zero() {
                brute.push((int(a), int(b)));
            }
        }
    }
    assert_eq!(r, brute);
    assert!(r.contains(&(int(0), int(0))));
}

#[test]
fn trivariate_planted() {
    let mut rng = Rng::new(3);
    for (bits, xb, k) in [(30u32, 5i64, 2u32), (40, 8, 1), (40, 12, 2)] {
        for _ in 0..2 {
            let g = random_gram(&mut rng, bits);
            let x = planted(&mut rng, xb);
            let d = quad_form(&g, &x);
            let f = norm_eq_poly(&g, &d, &x[3..]);
            let params = CoppersmithParams::new(k, Shape::Simplex);
            let b = int(xb + 1);
            let (roots, rep) =
                trivariate_small_roots_report(&f, &[b.clone(), b.clone(), b.clone()], &params);
            for r in &roots {
                assert!(f.eval(r).is_zero());
            }
            assert!(
                roots.contains(&[x[0].clone(), x[1].clone(), x[2].clone()]),
                "bits={bits} X={xb} k={k}: missing planted root; {rep:?}"
            );
        }
    }
}

#[test]
fn negative_no_small_root() {
    let mut rng = Rng::new(4);
    let mut tested_bi = 0;
    let mut tested_tri = 0;
    let xb = 12i64;
    while tested_bi < 3 || tested_tri < 2 {
        let g = random_gram(&mut rng, 24);
        let x = planted(&mut rng, xb);
        // perturb d so that (very likely) no solution exists; confirm by brute force
        let d = &quad_form(&g, &x) + &int(1);
        if tested_bi < 3 {
            let f = norm_eq_poly(&g, &d, &x[2..]);
            let mut any = false;
            for a in -xb..=xb {
                for b in -xb..=xb {
                    any |= f.eval(&[int(a), int(b)]).is_zero();
                }
            }
            if !any {
                let r = coron_bivariate(&f, &int(xb + 1), &int(xb + 1), 2);
                assert!(r.is_empty(), "false roots {r:?}");
                tested_bi += 1;
            }
        }
        if tested_tri < 2 {
            let f = norm_eq_poly(&g, &d, &x[3..]);
            let mut any = false;
            for a in -xb..=xb {
                for b in -xb..=xb {
                    for c in -xb..=xb {
                        any |= f.eval(&[int(a), int(b), int(c)]).is_zero();
                    }
                }
            }
            if !any {
                let b = int(xb + 1);
                let r = trivariate_small_roots(
                    &f,
                    &[b.clone(), b.clone(), b],
                    &CoppersmithParams::new(2, Shape::Simplex),
                );
                assert!(r.is_empty(), "false roots {r:?}");
                tested_tri += 1;
            }
        }
    }
    // and the full solver returns None when there is no small solution
    let g = random_gram(&mut rng, 24);
    let x = planted(&mut rng, 3);
    let d = &quad_form(&g, &x) + &int(1);
    let mut any = false;
    for a in -3i64..=3 {
        for b in -3i64..=3 {
            for c in -3i64..=3 {
                for e in -3i64..=3 {
                    any |= quad_form(&g, &[int(a), int(b), int(c), int(e)]) == d;
                }
            }
        }
    }
    if !any {
        assert!(solve_norm_eq_coppersmith(&g, &d, &int(3), NormEqMode::Bivariate).is_none());
        assert!(solve_norm_eq_coppersmith(&g, &d, &int(3), NormEqMode::Trivariate).is_none());
    }
}

#[test]
fn norm_equation_solver() {
    let mut rng = Rng::new(5);
    for mode in [NormEqMode::Bivariate, NormEqMode::Trivariate] {
        for _ in 0..2 {
            let g = random_gram(
                &mut rng,
                if mode == NormEqMode::Bivariate {
                    30
                } else {
                    40
                },
            );
            let x = planted(&mut rng, 6);
            let d = quad_form(&g, &x);
            let t = std::time::Instant::now();
            let sol =
                solve_norm_eq_coppersmith(&g, &d, &int(6), mode).expect("planted solution exists");
            assert_eq!(quad_form(&g, &sol), d);
            assert!(sol.iter().all(|v| v.abs() <= int(6)));
            eprintln!(
                "solve_norm_eq {mode:?}: planted {x:?} found {sol:?} in {:?}",
                t.elapsed()
            );
        }
    }
    // non-symmetric Gram convention: Q(x) = sum G[i][j] x_i x_j
    let mut g = random_gram(&mut rng, 20);
    let t = &g[0][1] + &int(5);
    g[0][1] = t;
    let t = &g[1][0] - &int(5);
    g[1][0] = t;
    let x = [int(2), int(-1), int(3), int(1)];
    let d = quad_form(&g, &x);
    let sol = solve_norm_eq_coppersmith(&g, &d, &int(3), NormEqMode::Bivariate).unwrap();
    assert_eq!(quad_form(&g, &sol), d);
}

/// Print (not assert) the largest X (powers of 2) for which the planted root is recovered.
#[test]
fn measure_success_bounds() {
    let bits = 30u32;
    eprintln!("--- measured success bounds (Gram entries ~2^{bits}; W ~ 2^{bits} X^2) ---");
    let cases: Vec<(&str, usize, Shape, u32)> = vec![
        ("bivariate box", 2, Shape::Box, 1),
        ("bivariate box", 2, Shape::Box, 2),
        ("bivariate simplex", 2, Shape::Simplex, 1),
        ("bivariate simplex", 2, Shape::Simplex, 2),
        ("bivariate simplex", 2, Shape::Simplex, 3),
        ("bivariate simplex", 2, Shape::Simplex, 4),
        ("trivariate simplex", 3, Shape::Simplex, 1),
        ("trivariate simplex", 3, Shape::Simplex, 2),
    ];
    let mut rng = Rng::new(6);
    let g = random_gram(&mut rng, bits);
    for (name, nv, shape, k) in cases {
        let params = CoppersmithParams::new(k, shape);
        let mut best = 0u32;
        let mut pred = 0u32;
        let mut dim = 0;
        for e in 1..=20u32 {
            let xb = (1i64 << e) - 1;
            let x: [Int; 4] = std::array::from_fn(|i| {
                // planted root near the bound (worst case)
                let v = xb - (i as i64);
                if i % 2 == 0 {
                    int(v)
                } else {
                    int(-v)
                }
            });
            let d = quad_form(&g, &x);
            let f = norm_eq_poly(&g, &d, &x[nv..]);
            let bnd = vec![int(xb + 1); nv];
            if lattice_margin_bits(&f, &bnd, &params) > 0.0 {
                pred = e;
            }
            let ok = if nv == 2 {
                let (r, rep) = coron_bivariate_with(&f, &bnd[0], &bnd[1], &params);
                dim = rep.dim;
                r.contains(&(x[0].clone(), x[1].clone()))
            } else {
                let (r, rep) = trivariate_small_roots_report(
                    &f,
                    &[bnd[0].clone(), bnd[1].clone(), bnd[2].clone()],
                    &params,
                );
                dim = rep.dim;
                r.contains(&[x[0].clone(), x[1].clone(), x[2].clone()])
            };
            if ok {
                best = e;
            } else if e > best + 2 {
                break;
            }
        }
        // asymptotic (k -> inf), with W ~ 2^bits X^2 for these instances:
        //   box (Coron, Newton polygon of the conic):  XY < W^(1/3)   ->  X < 2^(bits/4)
        //   simplex, 2 variables:                       XY < W^(1/2)   ->  X < 2^(bits/2)
        //   simplex, 3 variables:                       XYZ < W^(1/2)  ->  X < 2^(bits/4)
        let asym = match (nv, shape) {
            (2, Shape::Box) => bits as f64 / 4.0,
            (2, Shape::Simplex) => bits as f64 / 2.0,
            _ => bits as f64 / 4.0,
        };
        eprintln!(
            "{name:20} k={k} dim={dim:3}: recovered up to X = 2^{best:2}  (det heuristic 2^{pred:2}, asymptotic X < 2^{asym})"
        );
    }
}

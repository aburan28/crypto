//! lit2026 W11 (b): x-only isogeny evaluation on ordinary binary curves
//! (after Banegas–Sarkis–Smith, ePrint 2026/704), counted and timed.
//!
//! ```sh
//! cargo run --release --example lit2026_w11_binary_xonly -- sect163 <steps> <reps> <rounds>
//! cargo run --release --example lit2026_w11_binary_xonly -- torsion19 <reps> <rounds>
//! ```
//!
//! **`sect163`.**  Over `F_{2^163}` (the sect163 polynomial
//! `z^163 + z^7 + z^6 + z^3 + 1`), for `ℓ = 3, 5, 7`: take the first
//! `b = elt(v·0x9E3779B9 ⊕ 0x5bd1e995)` whose `y² + xy = x³ + x² + b` has an
//! x-rational ℓ-kernel, and walk `steps` ℓ-isogenies, carrying two
//! abscissae.  Each step is evaluated by every evaluator of
//! `binary_velu` Part E plus the module's original `velu_x_map`; all must
//! agree, the codomain is `velu_codomain`, and (for the carried point that
//! has an ordinate) `velu_point_map_binary` must give the same abscissa.
//! Counts are per step (two points).  Timing: `reps` evaluations per arm,
//! arms rotated over `rounds` rounds.
//!
//! **`torsion19`.**  The `(n, a₂) = (19, 1)` Koblitz class walk's kernels
//! (`ℓ = 457`, `E[ℓ] ⊂ E(F_{q^8})`, `s = 228`): the walk step's two
//! abscissa evaluations per kernel by `velu_x_from_abscissae` (twice) and
//! by `velu_x_pair_fraction`, with `F_{q^8}` operation counts and time.

use std::time::Instant;

use num_bigint::BigUint;
use rand::rngs::StdRng;
use rand::SeedableRng;

use crypto_lib::binary_ecc::{BinaryPoint, F2mElement, IrreduciblePoly};
use crypto_lib::cryptanalysis::binary_torsion_walk::{
    all_kernels, ext_ops, ext_ops_reset, mult_order, order_over_extension, scalar_eigenvalue,
    torsion_basis, velu_x_from_abscissae, velu_x_pair_fraction, Ext, ExtCurve, Torsion,
};
use crypto_lib::cryptanalysis::binary_velu::{
    derivative, elt, kernel_poly_counted, rational_kernels_wide, velu_codomain,
    velu_point_map_binary, velu_x_map, xeval_batch_inverse, xeval_fraction, xeval_fraction_multi,
    xeval_kernel_poly, xeval_product, xeval_projective, y_for_x_odd, Curve, F2mCounts, F2mCtx,
};
use crypto_lib::cryptanalysis::koblitz_isogeny_cost::field_for;

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

fn hex(e: &F2mElement) -> String {
    format!("{:x}", e.to_biguint())
}

fn counts_json(c: F2mCounts) -> String {
    format!(
        "{{\"mul\":{},\"sqr\":{},\"inv\":{},\"add\":{}}}",
        c.mul, c.sqr, c.inv, c.add
    )
}

const ARMS: [&str; 8] = [
    "velu_x_map_original",
    "velu_point_map_binary_original",
    "batch_inverse",
    "kernel_poly_given_h",
    "fraction",
    "fraction_multi",
    "product_normalised",
    "projective",
];

fn sect163(steps: usize, reps: usize, rounds: usize) {
    let irr = IrreduciblePoly::deg_163();
    let n = 163u32;
    let one = F2mElement::one(n);
    for ell in [3u32, 5, 7] {
        let s = ((ell - 1) / 2) as usize;
        let (v0, mut b) = (1u64..4000)
            .map(|v| (v, elt(v.wrapping_mul(0x9E37_79B9) ^ 0x5bd1_e995, n)))
            .find(|(_, b)| !rational_kernels_wide(b, ell, &irr).is_empty())
            .expect("an x-rational kernel");
        // carried points: one on E with an ordinate, one arbitrary abscissa
        let mut px = (2u64..)
            .map(|v| elt(v * 7919 + 13, n))
            .find(|x| y_for_x_odd(x, &one, &b, &irr).is_some())
            .unwrap();
        let mut py = y_for_x_odd(&px, &one, &b, &irr).unwrap();
        let mut qx = elt(0xC0FFEE_1234_5678, n);
        // projective carriers (start at Z = 1, never normalised in the walk)
        let mut pproj = (px.clone(), one.clone());
        let mut qproj = (qx.clone(), one.clone());
        let mut prev_b: Option<F2mElement> = None;
        for step in 0..steps {
            let t_k = Instant::now();
            let kernels = rational_kernels_wide(&b, ell, &irr);
            let kernel_find_s = t_k.elapsed().as_secs_f64();
            // prefer a kernel whose codomain is not the previous curve
            let mut pick = None;
            for (i, xs) in kernels.iter().enumerate() {
                let c = F2mCtx::new(&irr);
                let h = kernel_poly_counted(&c, xs);
                let b2 = velu_codomain(&b, &h, n, &irr);
                if pick.is_none() || prev_b.as_ref() != Some(&b2) {
                    pick = Some((i, h, b2));
                    if prev_b.as_ref() != Some(&pick.as_ref().unwrap().2) {
                        break;
                    }
                }
            }
            let (ki, h, b_next) = pick.expect("kernel");
            let xs = &kernels[ki];
            let backtrack = prev_b.as_ref() == Some(&b_next);
            let dh = derivative(&h, n);
            let sqrt_b = F2mCtx::new(&irr).sqrt(&b);

            // ── counted evaluation of both carried points per arm ──
            let mut counts: Vec<(&str, F2mCounts)> = Vec::new();
            let c = F2mCtx::new(&irr);
            let want_p = velu_x_map(&px, &h, n, &irr).expect("not a kernel point");
            let want_q = velu_x_map(&qx, &h, n, &irr).expect("not a kernel point");
            let full = velu_point_map_binary(
                &BinaryPoint::Affine {
                    x: px.clone(),
                    y: py.clone(),
                },
                &h,
                n,
                &irr,
            )
            .unwrap();
            let BinaryPoint::Affine { x: fx, y: fy } = full else {
                panic!("image O")
            };
            assert_eq!(fx, want_p, "velu_point_map_binary abscissa");
            // image on the codomain y² + xy = x³ + x² + b'
            let lhs = fy.square(&irr).add(&fx.mul(&fy, &irr));
            let fx2 = fx.square(&irr);
            assert_eq!(
                lhs,
                fx2.mul(&fx, &irr).add(&fx2).add(&b_next),
                "image on codomain"
            );

            c.reset();
            let a = xeval_batch_inverse(&c, &px, xs).unwrap();
            let bq = xeval_batch_inverse(&c, &qx, xs).unwrap();
            assert_eq!((a, bq), (want_p.clone(), want_q.clone()));
            counts.push(("batch_inverse", c.counts()));

            c.reset();
            let hh = kernel_poly_counted(&c, xs);
            let build_h = c.counts();
            assert!(hh.eq_poly(&h));
            c.reset();
            let a = xeval_kernel_poly(&c, &px, &h, &dh).unwrap();
            let bq = xeval_kernel_poly(&c, &qx, &h, &dh).unwrap();
            assert_eq!((a, bq), (want_p.clone(), want_q.clone()));
            counts.push(("kernel_poly_given_h", c.counts()));

            c.reset();
            let a = xeval_fraction(&c, &px, xs).unwrap();
            let bq = xeval_fraction(&c, &qx, xs).unwrap();
            assert_eq!((a, bq), (want_p.clone(), want_q.clone()));
            counts.push(("fraction", c.counts()));

            c.reset();
            let m = xeval_fraction_multi(&c, &[px.clone(), qx.clone()], xs).unwrap();
            assert_eq!(m, vec![want_p.clone(), want_q.clone()]);
            counts.push(("fraction_multi", c.counts()));

            // product form: per-kernel constant u^{-2}, u = Π x_i, counted
            c.reset();
            let mut u = xs[0].clone();
            for x in &xs[1..] {
                u = c.mul(&u, x);
            }
            let u2inv = c.inv(&c.sqr(&u)).unwrap();
            let norm_cost = c.counts();
            c.reset();
            let a = xeval_product(&c, &px, xs, &sqrt_b, Some(&u2inv)).unwrap();
            let bq = xeval_product(&c, &qx, xs, &sqrt_b, Some(&u2inv)).unwrap();
            assert_eq!((a, bq), (want_p.clone(), want_q.clone()));
            counts.push(("product_normalised", c.counts()));

            // projective: carried (X:Z) from previous steps
            c.reset();
            let np = xeval_projective(&c, (&pproj.0, &pproj.1), xs).unwrap();
            let nq = xeval_projective(&c, (&qproj.0, &qproj.1), xs).unwrap();
            counts.push(("projective", c.counts()));
            let norm =
                |p: &(F2mElement, F2mElement)| p.0.mul(&p.1.flt_inverse(&irr).unwrap(), &irr);
            assert_eq!(norm(&np), want_p, "projective P");
            assert_eq!(norm(&nq), want_q, "projective Q");

            // ── timing, arms rotated per round ──
            let mut times: Vec<Vec<f64>> = vec![Vec::new(); ARMS.len()];
            let pt = BinaryPoint::Affine {
                x: px.clone(),
                y: py.clone(),
            };
            for round in 0..rounds {
                for j in 0..ARMS.len() {
                    let arm = (j + round) % ARMS.len();
                    let cc = F2mCtx::new(&irr);
                    let t0 = Instant::now();
                    for _ in 0..reps {
                        match arm {
                            0 => {
                                std::hint::black_box(velu_x_map(&px, &h, n, &irr));
                                std::hint::black_box(velu_x_map(&qx, &h, n, &irr));
                            }
                            1 => {
                                std::hint::black_box(velu_point_map_binary(&pt, &h, n, &irr));
                            }
                            2 => {
                                std::hint::black_box(xeval_batch_inverse(&cc, &px, xs));
                                std::hint::black_box(xeval_batch_inverse(&cc, &qx, xs));
                            }
                            3 => {
                                std::hint::black_box(xeval_kernel_poly(&cc, &px, &h, &dh));
                                std::hint::black_box(xeval_kernel_poly(&cc, &qx, &h, &dh));
                            }
                            4 => {
                                std::hint::black_box(xeval_fraction(&cc, &px, xs));
                                std::hint::black_box(xeval_fraction(&cc, &qx, xs));
                            }
                            5 => {
                                std::hint::black_box(xeval_fraction_multi(
                                    &cc,
                                    &[px.clone(), qx.clone()],
                                    xs,
                                ));
                            }
                            6 => {
                                std::hint::black_box(xeval_product(
                                    &cc,
                                    &px,
                                    xs,
                                    &sqrt_b,
                                    Some(&u2inv),
                                ));
                                std::hint::black_box(xeval_product(
                                    &cc,
                                    &qx,
                                    xs,
                                    &sqrt_b,
                                    Some(&u2inv),
                                ));
                            }
                            _ => {
                                std::hint::black_box(xeval_projective(
                                    &cc,
                                    (&pproj.0, &pproj.1),
                                    xs,
                                ));
                                std::hint::black_box(xeval_projective(
                                    &cc,
                                    (&qproj.0, &qproj.1),
                                    xs,
                                ));
                            }
                        }
                    }
                    times[arm].push(t0.elapsed().as_secs_f64() * 1e6 / reps as f64);
                }
            }
            let med: Vec<String> = ARMS
                .iter()
                .zip(times.iter_mut())
                .map(|(name, t)| {
                    let lo = t.iter().cloned().fold(f64::INFINITY, f64::min);
                    let hi = t.iter().cloned().fold(0.0, f64::max);
                    format!(
                        "\"{name}\":{{\"median_us\":{:.3},\"min_us\":{lo:.3},\"max_us\":{hi:.3}}}",
                        median(t)
                    )
                })
                .collect();
            let cj: Vec<String> = counts
                .iter()
                .map(|(k, v)| format!("\"{k}\":{}", counts_json(*v)))
                .collect();
            println!(
                "{{\"mode\":\"sect163\",\"field\":\"z^163+z^7+z^6+z^3+1\",\"a2\":1,\"ell\":{ell},\"s\":{s},\"scan_v\":{v0},\"step\":{step},\"b\":\"{}\",\"b_next\":\"{}\",\"j_next\":\"{}\",\"kernels_available\":{},\"kernel_index\":{ki},\"backtrack\":{backtrack},\"kernel_find_s\":{kernel_find_s:.3},\"x_phi_p\":\"{}\",\"x_phi_q\":\"{}\",\"all_evaluators_agree\":true,\"build_h_counts\":{},\"product_norm_counts\":{},\"counts_two_points\":{{{}}},\"time_two_points\":{{{}}},\"reps\":{reps},\"rounds\":{rounds}}}",
                hex(&b),
                hex(&b_next),
                hex(&b_next.flt_inverse(&irr).unwrap()),
                kernels.len(),
                hex(&want_p),
                hex(&want_q),
                counts_json(build_h),
                counts_json(norm_cost),
                cj.join(","),
                med.join(",")
            );
            // advance the walk
            prev_b = Some(b.clone());
            b = b_next;
            px = fx;
            py = fy;
            qx = want_q;
            pproj = np;
            qproj = nq;
        }
    }
}

fn torsion19(reps: usize, rounds: usize) {
    let (n, a2, ell) = (19u32, 1u8, 457u64);
    let irr = field_for(n).expect("field");
    let q = 1u64 << n;
    let base = Curve::new(n, &irr, a2, 1).expect("K_1");
    let t = (q + 1) as i64 - base.order as i64;
    let lam = scalar_eigenvalue(q, t, ell).expect("ℓ | c");
    let m = mult_order(lam, ell) as usize;
    let ext = Ext::new(&irr, m, 17);
    let om: BigUint = order_over_extension(q, t, m);
    let ec = ExtCurve::new(&ext, a2 as u64, 1);
    let mut rng = StdRng::seed_from_u64(23);
    let Torsion::Rank2(p1, p2) = torsion_basis(&ec, &om, ell, 16, &mut rng) else {
        panic!("rank 2 expected at the crater")
    };
    let kernels = all_kernels(&ec, &p1, &p2, ell);
    let (xp, xq) = (0x1_2345 % q, 0x5_4321 % q);
    // identity + counts over every kernel
    let (mut old_ops, mut new_ops) = ([0u64; 3], [0u64; 3]);
    let mut agree = 0usize;
    for k in &kernels {
        ext_ops_reset();
        let want = [
            velu_x_from_abscissae(&ext, &k.abscissae, xp),
            velu_x_from_abscissae(&ext, &k.abscissae, xq),
        ];
        let o = ext_ops();
        ext_ops_reset();
        let got = velu_x_pair_fraction(&ext, &k.abscissae, [xp, xq]);
        let nn = ext_ops();
        assert_eq!(got, want, "pair fraction vs batch inversion");
        agree += 1;
        for i in 0..3 {
            old_ops[i] += o[i];
            new_ops[i] += nn[i];
        }
    }
    let mut t_old = Vec::new();
    let mut t_new = Vec::new();
    for round in 0..rounds {
        for j in 0..2 {
            let arm = (j + round) % 2;
            let t0 = Instant::now();
            for _ in 0..reps {
                for k in &kernels {
                    if arm == 0 {
                        std::hint::black_box(velu_x_from_abscissae(&ext, &k.abscissae, xp));
                        std::hint::black_box(velu_x_from_abscissae(&ext, &k.abscissae, xq));
                    } else {
                        std::hint::black_box(velu_x_pair_fraction(&ext, &k.abscissae, [xp, xq]));
                    }
                }
            }
            let us = t0.elapsed().as_secs_f64() * 1e6 / (reps * kernels.len()) as f64;
            if arm == 0 {
                t_old.push(us)
            } else {
                t_new.push(us)
            }
        }
    }
    let (lo_o, hi_o) = (
        t_old.iter().cloned().fold(f64::INFINITY, f64::min),
        t_old.iter().cloned().fold(0.0, f64::max),
    );
    let (lo_n, hi_n) = (
        t_new.iter().cloned().fold(f64::INFINITY, f64::min),
        t_new.iter().cloned().fold(0.0, f64::max),
    );
    let ratios: Vec<f64> = t_old.iter().zip(&t_new).map(|(a, b)| a / b).collect();
    println!(
        "{{\"mode\":\"torsion19\",\"n\":{n},\"a2\":{a2},\"ell\":{ell},\"ext_degree\":{m},\"s\":{},\"kernels\":{},\"identical\":{agree},\"ext_ops_batch_two_calls\":{{\"mul\":{},\"sqr\":{},\"inv\":{}}},\"ext_ops_pair_fraction\":{{\"mul\":{},\"sqr\":{},\"inv\":{}}},\"us_per_kernel_batch\":{{\"median\":{:.3},\"min\":{lo_o:.3},\"max\":{hi_o:.3}}},\"us_per_kernel_fraction\":{{\"median\":{:.3},\"min\":{lo_n:.3},\"max\":{hi_n:.3}}},\"round_ratios\":{:?},\"reps\":{reps},\"rounds\":{rounds}}}",
        (ell - 1) / 2,
        kernels.len(),
        old_ops[0],
        old_ops[1],
        old_ops[2],
        new_ops[0],
        new_ops[1],
        new_ops[2],
        median(&mut t_old.clone()),
        median(&mut t_new.clone()),
        ratios
    );
}

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let num = |i: usize, d: usize| -> usize {
        a.get(i).map(|s| s.parse().expect("integer")).unwrap_or(d)
    };
    match a.first().map(String::as_str).unwrap_or("sect163") {
        "torsion19" => torsion19(num(1, 3), num(2, 5)),
        _ => sect163(num(1, 8), num(2, 200), num(3, 5)),
    }
}

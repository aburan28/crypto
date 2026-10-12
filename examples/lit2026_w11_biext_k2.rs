//! lit2026 W11 (c): the embedding-degree-2 Tate pairing on the repository's
//! `k = 2` path (supersingular `y² = x³ + x`, `p ≡ 3 mod 4`, distortion map)
//! by the repo's affine Miller loop, a denominator-free Jacobian Miller loop,
//! and a cubical (biextension) ladder — after Zheng–Lin–Zhao, ePrint
//! 2025/1652 — with counted `F_p` operations and identical reduced values.
//!
//! ```sh
//! cargo run --release --example lit2026_w11_biext_k2 -- <pairs_per_size> <reps> <rounds>
//! ```
//!
//! Instances: for each `r_bits` in 16..=52 step 4, the prime `r` just below
//! `2^{r_bits}` (with `r ≡ 3 mod 4`... any odd prime) and the least `c ≡ 0
//! (mod 4)` making `p = c·r − 1` a prime `< 2^62` with `p ≡ 3 (mod 4)` —
//! `p = r·c − 1` with `4 | c`, so `r | p + 1 = #E`.  `P` is the least-`x`
//! point of order `r`, `R = [k]P` for deterministic `k`.  At the toy size
//! `(p, r) = (1283, 107)` the values are also checked against
//! `eds_mov::modified_tate`.

use std::time::Instant;

use crypto_lib::cryptanalysis::biext_k2::{
    compare, cubical_ladder, ec_mul, miller_affine, miller_jacobian, point_of_order, Fp,
};
use crypto_lib::cryptanalysis::challenge_ladder::is_prime_u64;
use crypto_lib::cryptanalysis::eds_mov;

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let num = |i: usize, d: usize| -> usize {
        a.get(i).map(|s| s.parse().expect("integer")).unwrap_or(d)
    };
    let (pairs, reps, rounds) = (num(0, 4), num(1, 200), num(2, 5));

    // toy check against eds_mov
    {
        let (p, r) = (1283u64, 107u64);
        let f = Fp::new(p);
        let pp = point_of_order(&f, r).unwrap();
        for k in 1..=5u128 {
            let rr = ec_mul(&f, Some(pp), k * 7).unwrap();
            let want = eds_mov::modified_tate(pp, rr, r, 1, p).unwrap();
            let c = compare(p, r, pp, rr).unwrap();
            assert_eq!(c.value_affine, want);
            assert_eq!(c.value_jacobian, want);
            assert_eq!(c.value_cubical, want);
        }
        println!("{{\"toy_check\":\"p=1283 r=107: affine, Jacobian and cubical equal eds_mov::modified_tate on 5 pairs\"}}");
    }

    let mut r_bits = 16u32;
    while r_bits <= 52 {
        // r: largest prime below 2^r_bits
        let mut r = (1u64 << r_bits) - 1;
        while !is_prime_u64(r) {
            r -= 2;
        }
        // p = c·r − 1 prime, c ≡ 0 mod 4, p < 2^62
        let mut c = 4u64;
        let p = loop {
            let p = c as u128 * r as u128 - 1;
            assert!(p < (1u128 << 62), "no p below 2^62");
            if is_prime_u64(p as u64) {
                break p as u64;
            }
            c += 4;
        };
        assert_eq!(p % 4, 3);
        let f = Fp::new(p);
        let pp = point_of_order(&f, r).expect("order-r point");
        for j in 0..pairs {
            let k = 0x1234_5678_9abc_u128.wrapping_mul(j as u128 + 3) % r as u128;
            let rr = ec_mul(&f, Some(pp), k.max(2)).unwrap();
            let cmp = compare(p, r, pp, rr).expect("pairing");
            let identical = cmp.value_affine == cmp.value_jacobian
                && cmp.value_jacobian == cmp.value_cubical
                && cmp.value_cubical != (1, 0);
            assert!(identical, "values must agree and be non-trivial");
            // timing (loop only; final exponentiation excluded)
            let mut times = [Vec::new(), Vec::new(), Vec::new()];
            for round in 0..rounds {
                for jj in 0..3 {
                    let arm = (jj + round) % 3;
                    let t0 = Instant::now();
                    for _ in 0..reps {
                        let g = Fp::new(p);
                        match arm {
                            0 => {
                                std::hint::black_box(miller_affine(&g, pp, rr, r));
                            }
                            1 => {
                                std::hint::black_box(miller_jacobian(&g, pp, rr, r));
                            }
                            _ => {
                                std::hint::black_box(cubical_ladder(&g, pp, rr, r));
                            }
                        }
                    }
                    times[arm].push(t0.elapsed().as_secs_f64() * 1e6 / reps as f64);
                }
            }
            let fmt = |t: &mut Vec<f64>| {
                let lo = t.iter().cloned().fold(f64::INFINITY, f64::min);
                let hi = t.iter().cloned().fold(0.0, f64::max);
                format!(
                    "{{\"median_us\":{:.3},\"min_us\":{lo:.3},\"max_us\":{hi:.3}}}",
                    median(t)
                )
            };
            let [ta, tj, tc] = &mut times;
            // An inversion is a Fermat exponentiation: log2(p−2) squarings
            // and popcount(p−2) multiplications.
            let inv_cost = (63 - (p - 2).leading_zeros() + (p - 2).count_ones()) as u64;
            let o = |c: &crypto_lib::cryptanalysis::biext_k2::FpCounts| {
                format!(
                    "{{\"mul\":{},\"sqr\":{},\"inv\":{},\"add\":{},\"m_plus_s\":{},\"m_equiv\":{}}}",
                    c.mul,
                    c.sqr,
                    c.inv,
                    c.add,
                    c.m_plus_s(),
                    c.m_plus_s() + c.inv * inv_cost
                )
            };
            println!(
                "{{\"r_bits\":{r_bits},\"r\":{r},\"p\":{p},\"p_bits\":{},\"k\":{k},\"identical\":{identical},\"value\":[{},{}],\"ops_miller_affine\":{},\"ops_miller_jacobian\":{},\"ops_cubical\":{},\"ops_unsquare\":{},\"raw_cubical_is_square\":{},\"cubical_over_jacobian_m_plus_s\":{:.4},\"cubical_plus_unsquare_over_jacobian_m_plus_s\":{:.4},\"cubical_plus_unsquare_over_jacobian_m_equiv\":{:.4},\"time_miller_affine\":{},\"time_miller_jacobian\":{},\"time_cubical\":{},\"reps\":{reps},\"rounds\":{rounds}}}",
                64 - p.leading_zeros(),
                cmp.value_cubical.0,
                cmp.value_cubical.1,
                o(&cmp.ops_affine),
                o(&cmp.ops_jacobian),
                o(&cmp.ops_cubical),
                o(&cmp.ops_unsquare),
                cmp.value_cubical_raw == f.f2_mul(cmp.value_affine, cmp.value_affine),
                cmp.ops_cubical.m_plus_s() as f64 / cmp.ops_jacobian.m_plus_s() as f64,
                (cmp.ops_cubical.m_plus_s() + cmp.ops_unsquare.m_plus_s()) as f64
                    / cmp.ops_jacobian.m_plus_s() as f64,
                (cmp.ops_cubical.m_plus_s() + cmp.ops_unsquare.m_plus_s() + cmp.ops_cubical.inv * inv_cost) as f64
                    / (cmp.ops_jacobian.m_plus_s() + cmp.ops_jacobian.inv * inv_cost) as f64,
                fmt(ta),
                fmt(tj),
                fmt(tc)
            );
        }
        r_bits += 4;
    }
}

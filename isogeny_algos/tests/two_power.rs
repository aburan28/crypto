use isogeny_algos::bigint::Big;
use isogeny_algos::curve::*;
use isogeny_algos::field::{is_prime, Field, Rng, Zp, Zp2};
use isogeny_algos::fp2::Fp2;
use isogeny_algos::kernel::chain::{ell_power_isogeny, Strategy};
use isogeny_algos::kernel::montgomery::{self, ladder_p, proj24, XZ};
use isogeny_algos::kernel::two_power::*;

fn affine_a<F: Field>(f: &F, k: (F::E, F::E)) -> F::E {
    montgomery::affine_a(f, k)
}

/// x of a point of order exactly 2^e on y^2 = x^3 + 6 x^2 + x over F_{p^2}, p = 2^e c - 1,
/// with 2^(e-1) R != (0, 0) (so that every 4-isogeny step is non-degenerate).
fn kernel_point<F: Field>(f: &F, e: usize, cof: &Big, rng: &mut Rng) -> (F::E, F::E) {
    let a = f.from_u64(6);
    let k = proj24(f, a);
    loop {
        let x = f.random(rng);
        let rhs = f.mul(x, f.add(f.mul(x, f.add(x, a)), f.one()));
        let Some(y) = f.sqrt(rhs) else { continue };
        let r = ladder_p(f, k, (x, f.one()), cof);
        let h = xdbl_e(f, k, r, e - 1);
        if f.is_zero(h.1) || f.is_zero(h.0) {
            continue;
        }
        if !f.is_zero(xdbl_e(f, k, h, 1).1) {
            continue;
        }
        let _ = y;
        return (f.div(r.0, r.1), f.one());
    }
}

#[test]
fn four_and_two_chains_agree_with_weierstrass_velu() {
    let (e, b) = (20usize, 6u32);
    let base = (1u64 << e) * 3u64.pow(b);
    let p = (1..)
        .map(|c| base * c - 1)
        .find(|&p| is_prime(p) && p % 4 == 3)
        .unwrap();
    let cof = Big::from_u64((p + 1) >> e);
    let f = Zp2::new(p);
    let mut rng = Rng::new(31);
    let a0 = f.from_u64(6);
    for trial in 0..4 {
        let (rx, _) = kernel_point(&f, e, &cof, &mut rng);
        let r: XZ<(u64, u64)> = (rx, f.one());
        let extra = (f.random(&mut rng), f.one());
        let n4 = e / 2;
        let mut j_ref = None;
        for (name, splits) in [
            ("optimal", optimal_splits(n4, 1.0, 1.0)),
            (
                "naive",
                (0..=n4).map(|m| m.saturating_sub(1)).collect::<Vec<_>>(),
            ),
        ] {
            let mut st = TwoPowerStats::default();
            let mut pts = vec![extra];
            let cod =
                four_chain(&f, proj24(&f, a0), r, e, &mut pts, &splits, &mut st).expect("4-chain");
            let j = montgomery::j_invariant(&f, affine_a(&f, cod));
            assert_eq!(st.isogenies, n4, "{name}");
            if let Some(jr) = j_ref {
                assert_eq!(j, jr, "4-chain {name}");
            }
            j_ref = Some(j);
        }
        let j4 = j_ref.unwrap();
        // the same kernel as e two-isogenies
        let mut st = TwoPowerStats::default();
        let cod2 = two_chain(
            &f,
            proj24(&f, a0),
            r,
            e,
            &mut vec![],
            &optimal_splits(e, 1.0, 1.0),
            &mut st,
        )
        .expect("2-chain");
        assert_eq!(
            montgomery::j_invariant(&f, affine_a(&f, cod2)),
            j4,
            "2-chain"
        );
        // Weierstrass Velu chain (y from the curve equation)
        let w = montgomery::to_weierstrass(&f, a0);
        let xw = f.add(rx, f.div(a0, f.from_u64(3)));
        let yw = f.sqrt(rhs(&f, &w, xw)).unwrap();
        let (ch, _, _) = ell_power_isogeny(&f, &w, &Pt::Aff(xw, yw), 2, e, &[], Strategy::Balanced);
        assert_eq!(jinv(&f, &ch.cod), j4, "trial {trial}");
    }
}

#[test]
fn fp2_generic_matches_zp2() {
    let p = 1_000_000_007u64; // = 3 mod 4
    let a = Zp2::new(p);
    let g = Fp2::new(Zp::new(p));
    let mut rng = Rng::new(5);
    for _ in 0..200 {
        let (x, y) = (a.random(&mut rng), a.random(&mut rng));
        assert_eq!(a.mul(x, y), g.mul(x, y));
        assert_eq!(a.sq(x), g.sq(x));
        if x != a.zero() {
            assert_eq!(a.inv(x), g.inv(x));
        }
        let s = g.sqrt(g.sq(x)).unwrap();
        assert!(s == x || s == g.neg(x));
    }
}

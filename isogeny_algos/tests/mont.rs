use isogeny_algos::curve::*;
use isogeny_algos::field::*;
use isogeny_algos::kernel::{montgomery::*, velu::velu_from_point};

const PRIMES: [u64; 7] = [3, 5, 7, 11, 13, 17, 19];

fn csidh_p(n: usize) -> u64 {
    4 * PRIMES[..n].iter().product::<u64>() - 1
}

#[test]
fn montgomery_velu_matches_weierstrass_along_a_chain() {
    let mut rng = Rng::new(31);
    let p = csidh_p(6); // 4*3*5*7*11*13*17 - 1 = 1021019
    assert!(is_prime(p));
    let fp = Zp::new(p);
    let mut a = 0u64; // E_0: y^2 = x^3 + x is supersingular for p = 3 mod 4
    for step in 0..12 {
        let ell = PRIMES[step % 6];
        let e24 = a24(&fp, a);
        let w = to_weierstrass(&fp, a);
        let rhs = |u: u64| fp.add(fp.mul(fp.add(fp.mul(u, u), fp.mul(a, u)), u), u);
        // a point of order ell on E_a(F_p) (rational points: rhs square), cofactor (p+1)/ell
        let (kx, ky) = loop {
            let x = fp.random(&mut rng);
            if fp.sqrt(rhs(x)).is_none() {
                continue;
            }
            let k = ladder(&fp, e24, x, ((p + 1) / ell) as u128);
            if k.1 == 0 {
                continue;
            }
            let kx = fp.div(k.0, k.1);
            break (kx, fp.sqrt(rhs(kx)).unwrap());
        };
        let d = ((ell - 1) / 2) as usize;
        let ms = multiples(&fp, e24, (kx, 1), d);
        assert_eq!(
            ladder(&fp, e24, kx, ell as u128).1,
            0,
            "kernel point must have order l"
        );
        let a_new = velu_codomain(&fp, a, &ms, ell);
        // Weierstrass Velu on the converted curve
        let pw = Pt::Aff(fp.add(kx, fp.div(a, 3)), ky);
        let v = velu_from_point(&fp, &w, &pw, ell);
        assert_eq!(
            j_invariant(&fp, a_new),
            jinv(&fp, &v.cod),
            "step {step}, l = {ell}"
        );
        // x-map sends points of E_a to points of E_a' (square rhs stays square: B' = B pi^2)
        let rhs_new = |u: u64| fp.add(fp.mul(fp.add(fp.mul(u, u), fp.mul(a_new, u)), u), u);
        for _ in 0..20 {
            let u = fp.random(&mut rng);
            if fp.sqrt(rhs(u)).is_some() {
                let u2 = isog_x(&fp, &ms, u);
                assert!(
                    fp.sqrt(rhs_new(u2)).is_some(),
                    "image of a rational point must be rational"
                );
            }
        }
        // the codomain is again supersingular with p+1 points: [p+1]P' = O for a random x
        let e24n = a24(&fp, a_new);
        let x = fp.random(&mut rng);
        assert_eq!(
            ladder(&fp, e24n, x, (p + 1) as u128).1,
            0,
            "codomain must have order divisible by p+1"
        );
        a = a_new;
    }
}

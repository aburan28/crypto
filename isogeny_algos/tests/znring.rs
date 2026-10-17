//! Z/nZ products and dot products against schoolbook u128 arithmetic, for odd (Montgomery
//! kernels) and even (direct accumulation) moduli, prime and composite.
use isogeny_algos::field::{Field, Rng};
use isogeny_algos::strassen::ZnRing;

#[test]
fn znring_conv_and_dot_match_schoolbook() {
    let mut rng = Rng::new(7);
    for n in [
        2u64,
        12,
        1 << 61,
        (1 << 62) - 2,
        15,
        1_000_003 * 999_983,
        (1 << 61) - 1,
    ] {
        let r = ZnRing::new(n);
        for (la, lb) in [(0usize, 5usize), (1, 1), (3, 7), (40, 33), (200, 150)] {
            let a: Vec<u64> = (0..la).map(|_| r.random(&mut rng)).collect();
            let b: Vec<u64> = (0..lb).map(|_| r.random(&mut rng)).collect();
            let want: Vec<u64> = (0..(la + lb).saturating_sub(1))
                .map(|t| {
                    (0..la)
                        .filter(|&i| t >= i && t - i < lb)
                        .fold(0u128, |s, i| {
                            (s + a[i] as u128 * b[t - i] as u128) % n as u128
                        }) as u64
                })
                .collect();
            assert_eq!(r.conv(&a, &b), want, "n = {n}, {la} x {lb}");
            let k = want.len() / 2;
            assert_eq!(r.conv_trunc(&a, &b, k), want[..k], "n = {n}");
            let m = la.min(lb);
            let dot = (0..m).fold(0u128, |s, i| {
                (s + a[i] as u128 * b[lb - 1 - i] as u128) % n as u128
            });
            assert_eq!(r.dot_rev(&a, &b) as u128, dot, "n = {n}");
        }
    }
}

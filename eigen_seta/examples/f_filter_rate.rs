//! How many candidate secret degrees f survive key generation's number-
//! theoretic filters?  These filters are public (an attacker who guesses f
//! can apply them), so their acceptance rate lowers the effective size of
//! the f search space.  Replicates step 1 of `keygen_attempt` (kem.rs):
//!   f odd, gcd(f, 21) = 1, 2f != 2^E (mod d), m = 2^E - 7 d^2 f^2 a sum of
//!   two squares found by `two_squares`, and the schedule condition
//!   max(v2(c' - lamS), v2(c' + lamS)) = 2 for some sign choice.
//! The later keygen rejections (alpha search, gluing check, factor
//! identification, self-test) are not counted here.
//!
//!   cargo run --release --example f_filter_rate -- <Toy16|Toy32|L1|L3|L5> <accepts>
use eigen_seta::ntheory::{modpos, random_below, sqrt_mod_2power, two_squares, v2};
use eigen_seta::params::*;
use num_bigint::{BigInt, BigUint};
use num_traits::{One, Zero};
use rand_core::OsRng;

fn go<P: ParamSet<N>, const N: usize>(want: usize) {
    let mut rng = OsRng;
    let e = P::E;
    let d = BigUint::from(3u32).pow(P::F3 as u32);
    let n2 = BigUint::one() << e;
    let s_d = sqrt_mod_2power(&((BigUint::one() << (e + 2)) - 7u32), e as u64 + 2);
    let s_dm = &s_d % &n2;
    let lo = BigUint::parse_bytes(P::F_LO_HEX.as_bytes(), 16).unwrap();
    let (mut drawn, mut parity_gcd, mut two_f, mut sq, mut sched, mut acc) =
        (0u64, 0u64, 0u64, 0u64, 0u64, 0usize);
    while acc < want {
        drawn += 1;
        let f = &lo + random_below(&mut rng, &lo);
        let f = if (&f % 2u32).is_zero() { f + 1u32 } else { f };
        if (&f % 3u32).is_zero() || (&f % 7u32).is_zero() {
            parity_gcd += 1;
            continue;
        }
        let tf = BigInt::from(f.clone()) * 2 - BigInt::from(n2.clone());
        if modpos(&tf, &d).is_zero() {
            two_f += 1;
            continue;
        }
        let m = &n2 - &d * &d * &f * &f * 7u32;
        let Some((x, y)) = two_squares(&m, &mut rng, 1 << 14) else {
            sq += 1;
            continue;
        };
        let (c0, cp0) = if (&x % 2u32).is_zero() {
            (x, y)
        } else {
            (y, x)
        };
        let s0 = (&f * &s_dm) % &n2;
        let ds0 = BigInt::from((&d * &s0) % &n2);
        let mut ok = false;
        for sc in [1i32, -1] {
            for scp in [1i32, -1] {
                let c = BigInt::from(c0.clone()) * sc;
                let cp = BigInt::from(cp0.clone()) * scp;
                let lam_s = &ds0 + &c;
                ok |= v2(&(&cp - &lam_s)).max(v2(&(&cp + &lam_s))) == 2;
            }
        }
        if !ok {
            sched += 1;
            continue;
        }
        acc += 1;
    }
    let rate = acc as f64 / drawn as f64;
    println!(
        "{} drawn {} accepted {} rate {:.5} (log2 {:.2}); rejected: gcd(f,21) {} | 2f=2^E mod d {} | two-squares {} | schedule {}",
        P::NAME, drawn, acc, rate, rate.log2(), parity_gcd, two_f, sq, sched
    );
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let want: usize = a.get(2).and_then(|s| s.parse().ok()).unwrap_or(100);
    match a.get(1).map(|s| s.as_str()).unwrap_or("Toy16") {
        "Toy16" => go::<Toy16, 4>(want),
        "Toy32" => go::<Toy32, 6>(want),
        "L1" => go::<L1, 21>(want),
        "L3" => go::<L3, 32>(want),
        "L5" => go::<L5, 42>(want),
        _ => panic!("set"),
    }
}

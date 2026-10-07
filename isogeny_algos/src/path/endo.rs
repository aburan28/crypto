//! Kohel's algorithm (thesis 1996) for the endomorphism ring of an ordinary curve over F_p:
//! End(E) is the order of conductor f_E = prod l^{level_l(E)} in the maximal order O_K, where
//! level_l is the depth of E in the l-isogeny volcano of the Frobenius order Z[pi]
//! (height h_l = v_l(f_pi)); levels are found by comparing lengths of non-backtracking walks.
use super::graph::*;
use super::volcano::*;
use crate::field::{Field, Rng};

#[derive(Debug, Clone)]
pub struct EndRing {
    pub frobenius_disc: i128,
    pub fundamental_disc: i128,
    pub conductor_frobenius: u128,
    /// (l, height h_l, level of E)
    pub per_prime: Vec<(u64, u32, u32)>,
    pub conductor: u128,
    pub disc: i128,
}

fn factor_small(mut n: u128) -> (Vec<(u64, u32)>, u128) {
    let mut out = vec![];
    let mut q = 2u64;
    while (q as u128) * (q as u128) <= n && q < 2_000_000 {
        if n % q as u128 == 0 {
            let mut v = 0;
            while n % q as u128 == 0 {
                n /= q as u128;
                v += 1;
            }
            out.push((q, v));
        }
        q += if q == 2 { 1 } else { 2 };
    }
    (out, n)
}

/// End(E) for E/F_p with Frobenius trace `trace` and j-invariant j. All primes dividing f_pi must be odd
/// and have Phi_l in `cache`.
pub fn endomorphism_ring<F: Field>(
    f: &F,
    cache: &PhiCache,
    p: u64,
    trace: i64,
    j: F::E,
    rng: &mut Rng,
) -> Result<EndRing, String> {
    let d: i128 = (trace as i128) * (trace as i128) - 4 * (p as i128);
    if d >= 0 {
        return Err("not ordinary".into());
    }
    let (fac, rest) = factor_small(d.unsigned_abs());
    let mut fpi: u128 = 1;
    for &(q, v) in &fac {
        fpi *= (q as u128).pow(v / 2);
    }
    if rest > 1 {
        let s = (rest as f64).sqrt() as u128;
        if (s.saturating_sub(1)..=s + 1).any(|t| t * t == rest && t > 1) {
            return Err("square cofactor with a large prime in the conductor".into());
        }
    }
    while (d / ((fpi * fpi) as i128)).rem_euclid(4) > 1 {
        fpi /= 2;
    }
    let dk = d / ((fpi * fpi) as i128);
    let mut per_prime = vec![];
    let mut conductor: u128 = 1;
    for &(q, _) in &fac {
        // height of the l-volcano = v_l(f_pi)
        let h_actual = (0..)
            .take_while(|k| fpi % (q as u128).pow(k + 1) == 0)
            .count() as u32;
        if h_actual == 0 {
            continue;
        }
        if q == 2 {
            return Err("2 divides the conductor of Z[pi]: not supported".into());
        }
        if !cache.ells().contains(&(q as usize)) {
            return Err(format!("Phi_{q} not available"));
        }
        let (level, _) = level_and_up(f, cache, q as usize, h_actual, j, rng);
        per_prime.push((q, h_actual, level));
        conductor *= (q as u128).pow(level);
    }
    Ok(EndRing {
        frobenius_disc: d,
        fundamental_disc: dk,
        conductor_frobenius: fpi,
        per_prime,
        conductor,
        disc: dk * ((conductor * conductor) as i128),
    })
}

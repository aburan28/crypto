//! Route 1: the modular curves that carry the Weil restriction.
//!
//! Let `f` be the newform of a CM elliptic curve over `Q` with CM by
//! `Q(√−7)` reducing to `E_a` mod 2 — the conductor-49 curve for `a = 1`, its
//! twist by `χ_{−3}` (conductor 441) for `a = 0`, since `χ_{−3}(2) = −1` flips
//! the trace.  Let `ℓ = 2n + 1` be a prime with `⟨2, −1⟩ = (Z/ℓ)^*`
//! (`ℓ = 7, 11, 263` for `n = 3, 5, 131`), so 2 is inert in `Q(ζ_ℓ)^+` with
//! residue degree `n`, and let `χ` be an even character of conductor `ℓ` and
//! order `n`.  The newform `f ⊗ χ` has coefficient field `Q(ζ_n)`, nebentypus
//! `χ²` and level `N = cond(f)·ℓ²` (no change when `ℓ | cond(f)`); its abelian
//! variety has dimension `n − 1` and, by Eichler–Shimura, Frobenius at 2 with
//! eigenvalues `χ(2)^σ·{α, ᾱ}` — exactly the Weil numbers of `A_n(E_a)`.  It is a
//! factor of `J_H(N)` with `±H = {d : d ≡ ±1 mod ℓ}`, the kernel of `χ²`, so the
//! reduction mod 2 of `X_H(N)` is an explicit curve over `F_2` whose Jacobian
//! contains `A_n(E_a)`.  This module computes its genus.

use num_integer::Integer;
use serde::Serialize;

#[derive(Clone, Debug, Serialize)]
pub struct CosetGenus {
    pub level: u64,
    /// `None` for `Γ_0(N)`; `Some(ℓ)` for `Γ_H(N)` with `±H = {d ≡ ±1 mod ℓ}`.
    pub ell: Option<u64>,
    pub index_mu: u64,
    pub e2: u64,
    pub e3: u64,
    pub cusps: u64,
    pub genus: i64,
}

/// Genus of `Γ_0(N)` or `Γ_H(N)` by enumerating the cosets of the group in
/// `SL_2(Z)` — bottom rows `(c, d) mod N`, primitive, modulo `±H` — and the
/// permutation action of `S`, `T` and `ST` on them.
pub fn coset_genus(n: u64, ell: Option<u64>) -> CosetGenus {
    let nn = n as usize;
    let units: Vec<u64> = (1..n)
        .filter(|&u| u.gcd(&n) == 1)
        .filter(|&u| match ell {
            None => true,
            Some(l) => u % l == 1 || u % l == l - 1,
        })
        .collect();
    let mut id = vec![u32::MAX; nn * nn];
    let mut classes: Vec<(u64, u64)> = Vec::new();
    for c in 0..n {
        for d in 0..n {
            if c.gcd(&d).gcd(&n) != 1 || id[(c * n + d) as usize] != u32::MAX {
                continue;
            }
            let k = classes.len() as u32;
            for &u in &units {
                id[((u * c % n) * n + u * d % n) as usize] = k;
            }
            classes.push((c, d));
        }
    }
    let idx = |c: u64, d: u64| id[(c * n + d) as usize] as usize;
    let mu = classes.len();
    let mut e2 = 0u64;
    let mut e3 = 0u64;
    let mut t_perm = vec![0usize; mu];
    for (k, &(c, d)) in classes.iter().enumerate() {
        if idx(d, (n - c) % n) == k {
            e2 += 1;
        }
        if idx(d, (d + n - c) % n) == k {
            e3 += 1;
        }
        t_perm[k] = idx(c, (c + d) % n);
    }
    let mut seen = vec![false; mu];
    let mut cusps = 0u64;
    for k in 0..mu {
        if seen[k] {
            continue;
        }
        cusps += 1;
        let mut j = k;
        while !seen[j] {
            seen[j] = true;
            j = t_perm[j];
        }
    }
    let twelve_g = 12 + mu as i64 - 3 * e2 as i64 - 4 * e3 as i64 - 6 * cusps as i64;
    assert_eq!(twelve_g % 12, 0, "genus formula must give an integer");
    CosetGenus {
        level: n,
        ell,
        index_mu: mu as u64,
        e2,
        e3,
        cusps,
        genus: twelve_g / 12,
    }
}

pub fn factorize(mut n: u64) -> Vec<(u64, u32)> {
    let mut out = Vec::new();
    let mut p = 2;
    while p * p <= n {
        if n.is_multiple_of(p) {
            let mut k = 0;
            while n.is_multiple_of(p) {
                n /= p;
                k += 1;
            }
            out.push((p, k));
        }
        p += 1;
    }
    if n > 1 {
        out.push((n, 1));
    }
    out
}

fn kron_minus1(p: u64) -> i64 {
    if p == 2 {
        0
    } else if p % 4 == 1 {
        1
    } else {
        -1
    }
}

fn kron_minus3(p: u64) -> i64 {
    if p == 3 {
        0
    } else if p % 3 == 1 {
        1
    } else {
        -1
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct Gamma0Data {
    pub level: u64,
    pub mu0: u64,
    pub e2: u64,
    pub e3: u64,
    pub cusps0: u64,
    pub genus0: i64,
}

/// Index, elliptic points, cusps and genus of `Γ_0(N)` from the standard formulas.
pub fn gamma0(n: u64) -> Gamma0Data {
    let fac = factorize(n);
    let mu0: u64 = fac.iter().map(|&(p, k)| p.pow(k - 1) * (p + 1)).product();
    let e2 = if n.is_multiple_of(4) {
        0
    } else {
        fac.iter()
            .map(|&(p, _)| (1 + kron_minus1(p)) as u64)
            .product()
    };
    let e3 = if n.is_multiple_of(9) {
        0
    } else {
        fac.iter()
            .map(|&(p, _)| (1 + kron_minus3(p)) as u64)
            .product()
    };
    let cusps0: u64 = fac
        .iter()
        .map(|&(p, k)| {
            (0..=k)
                .map(|j| {
                    let m = j.min(k - j);
                    if m == 0 {
                        1
                    } else {
                        p.pow(m - 1) * (p - 1)
                    }
                })
                .sum::<u64>()
        })
        .product();
    let twelve_g = 12 + mu0 as i64 - 3 * e2 as i64 - 4 * e3 as i64 - 6 * cusps0 as i64;
    assert_eq!(twelve_g % 12, 0);
    Gamma0Data {
        level: n,
        mu0,
        e2,
        e3,
        cusps0,
        genus0: twelve_g / 12,
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct ModularRow {
    pub n: u32,
    pub koblitz_a: u8,
    pub ell: u64,
    pub newform_level: u64,
    pub level: u64,
    pub index_over_gamma0: u64,
    pub gamma0: Gamma0Data,
    /// `Γ_0(N)` has no elliptic points, so `X_H → X_0` can ramify only at cusps.
    pub elliptic_free: bool,
    /// Genus if every cusp splits — the cover is then étale, `g = n(g_0 − 1) + 1`.
    /// It does: the cover is cut out by a character of conductor `ℓ`, and `ℓ`
    /// divides `N/gcd(c, N/c)` for every cusp denominator `c`, so every cusp
    /// stabiliser lies in `±H` (validated by enumeration at small levels).
    pub genus_split_cusps: Option<i64>,
    /// Genus if every cusp ramified: the upper bound, `+ (n − 1)·c_0/2`.
    pub genus_upper_bound: Option<i64>,
    pub genus_exact_by_enumeration: Option<i64>,
    pub log2_genus: f64,
    pub dim_target: u32,
}

/// The modular curve carrying `A_n(E_a)`.
pub fn modular_row(n: u32, a: u8, ell: u64, enumerate_if_level_at_most: u64) -> ModularRow {
    let newform_level = if a == 1 { 49 } else { 441 };
    let level = if newform_level % ell == 0 {
        newform_level
    } else {
        newform_level * ell * ell
    };
    let g0 = gamma0(level);
    let idx = n as u64;
    let elliptic_free = g0.e2 == 0 && g0.e3 == 0;
    let split = elliptic_free.then(|| idx as i64 * (g0.genus0 - 1) + 1);
    let upper = split.map(|s| s + (idx as i64 - 1) * g0.cusps0 as i64 / 2);
    let exact = (level <= enumerate_if_level_at_most).then(|| coset_genus(level, Some(ell)).genus);
    assert!(
        exact.is_some() || elliptic_free,
        "levels with elliptic points must be enumerated"
    );
    ModularRow {
        n,
        koblitz_a: a,
        ell,
        newform_level,
        level,
        index_over_gamma0: idx,
        gamma0: g0,
        elliptic_free,
        genus_split_cusps: split,
        genus_upper_bound: upper,
        genus_exact_by_enumeration: exact,
        log2_genus: (exact.or(split).expect("one of the two") as f64).log2(),
        dim_target: n - 1,
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct CuspValidation {
    pub level: u64,
    pub ell: u64,
    pub enumerated: CosetGenus,
    pub gamma0_formula: Gamma0Data,
    pub gamma0_enumerated_genus: i64,
    pub cusps_split_completely: bool,
}

/// At levels small enough to enumerate, check (i) the `Γ_0(N)` formulas and
/// (ii) that every cusp of `X_0(N)` splits completely in `X_H(N)`.
pub fn validate_cusps(level: u64, ell: u64) -> CuspValidation {
    let g0 = gamma0(level);
    let enum0 = coset_genus(level, None);
    let enumerated = coset_genus(level, Some(ell));
    let idx = enumerated.index_mu / g0.mu0;
    CuspValidation {
        level,
        ell,
        cusps_split_completely: enumerated.cusps == idx * g0.cusps0,
        gamma0_enumerated_genus: enum0.genus,
        gamma0_formula: g0,
        enumerated,
    }
}

/// Facts the construction rests on at `n = 131`.
#[derive(Clone, Debug, Serialize)]
pub struct WeilMatch {
    pub a2_conductor49_curve: i64,
    pub kronecker_minus3_at_2: i64,
    pub a2_after_twist: i64,
    pub koblitz_trace_e0: i64,
    pub order_of_2_mod_263_up_to_sign: u64,
    pub two_inert_in_real_subfield: bool,
    /// `log₂ ∏_ζ |T² + ζT + 2ζ²|` (ζ primitive 131st roots) against
    /// `log₂ P_A(T)` at `T = 1, 2, 3`.
    pub product_vs_trace_zero_poly: Vec<(i64, f64, f64)>,
}

pub fn weil_match(pa: &[num_bigint::BigInt]) -> WeilMatch {
    use super::weil::{eval, log2_bigint};
    use num_bigint::BigInt;
    // a_2 of y² + xy = x³ − x² − 2x − 1: its reduction y² + xy = x³ + x² + 1 has
    // the points O and (0, 1) over F_2.
    let mut pts = 1i64;
    for x in 0..2i64 {
        for y in 0..2i64 {
            if (y * y + x * y - (x * x * x + x * x + 1)).rem_euclid(2) == 0 {
                pts += 1;
            }
        }
    }
    let a2 = 3 - pts;
    // (−3 | 2) = −1 since −3 ≡ 5 mod 8
    let kron = if (-3i64).rem_euclid(8) == 1 || (-3i64).rem_euclid(8) == 7 {
        1
    } else {
        -1
    };
    let mut k = 1u64;
    let mut x = 2u64;
    while x != 1 && x != 262 {
        x = x * 2 % 263;
        k += 1;
    }
    let mut checks = Vec::new();
    for t in [1i64, 2, 3] {
        let mut acc = 0.0f64;
        for j in 1..131u32 {
            let ang = 2.0 * std::f64::consts::PI * j as f64 / 131.0;
            let (zr, zi) = (ang.cos(), ang.sin());
            let (z2r, z2i) = ((2.0 * ang).cos(), (2.0 * ang).sin());
            let re = (t * t) as f64 + t as f64 * zr + 2.0 * z2r;
            let im = t as f64 * zi + 2.0 * z2i;
            acc += 0.5 * (re * re + im * im).log2();
        }
        checks.push((t, acc, log2_bigint(&eval(pa, &BigInt::from(t)))));
    }
    WeilMatch {
        a2_conductor49_curve: a2,
        kronecker_minus3_at_2: kron,
        a2_after_twist: a2 * kron,
        koblitz_trace_e0: -1,
        order_of_2_mod_263_up_to_sign: k,
        two_inert_in_real_subfield: k == 131,
        product_vs_trace_zero_poly: checks,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn x0_49_has_genus_one_and_x7_genus_three() {
        assert_eq!(coset_genus(49, None).genus, 1);
        assert_eq!(gamma0(49).genus0, 1);
        let x7 = coset_genus(49, Some(7));
        assert_eq!((x7.index_mu, x7.cusps, x7.genus), (168, 24, 3));
    }

    #[test]
    fn x0_121_formula_agrees() {
        assert_eq!(gamma0(121).genus0, 6);
        assert_eq!(coset_genus(121, None).genus, 6);
    }
}

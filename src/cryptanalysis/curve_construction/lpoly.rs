//! Characteristic polynomials of Frobenius from point counts.
//!
//! For a genus-`g` curve over `F_q`, the counts `N_1 … N_g` determine
//! `L(T) = ∏(1 − α_i T)` through Newton's identities
//! `k·c_k = Σ_{i≤k} S_i c_{k−i}` with `S_k = N_k − q^k − 1`, and the functional
//! equation `c_{g+j} = q^j c_{g−j}`.  Everything here is returned as the monic
//! characteristic polynomial `P(T) = T^{2g} L(1/T) = ∏(T − α_i)`, the
//! convention in which `P(1) = #Jac(F_q)` and `P_A | P_C` means "the
//! Jacobian of `C` contains `A` up to isogeny".

/// Monic `P(T)`, little-endian, or `None` when the counts violate Newton
/// integrality (a singular or non-curve model).
pub fn char_poly_from_counts(counts: &[i64], q: i64) -> Option<Vec<i128>> {
    let g = counts.len();
    let q = q as i128;
    let s: Vec<i128> = counts
        .iter()
        .enumerate()
        .map(|(i, &n)| n as i128 - q.pow(i as u32 + 1) - 1)
        .collect();
    let mut c = vec![0i128; 2 * g + 1];
    c[0] = 1;
    for k in 1..=g {
        let num: i128 = (1..=k).map(|i| s[i - 1] * c[k - i]).sum();
        if num % k as i128 != 0 {
            return None;
        }
        c[k] = num / k as i128;
    }
    for j in 1..=g {
        c[g + j] = q.pow(j as u32) * c[g - j];
    }
    let mut p = vec![0i128; 2 * g + 1];
    for (i, ci) in c.iter().enumerate() {
        p[2 * g - i] = *ci;
    }
    Some(p)
}

/// Exact division by a monic polynomial; `None` if the remainder is non-zero.
pub fn divide_exact(num: &[i128], den: &[i128]) -> Option<Vec<i128>> {
    assert_eq!(
        *den.last().expect("non-empty"),
        1,
        "denominator must be monic"
    );
    if num.len() < den.len() {
        return None;
    }
    let mut rem = num.to_vec();
    let mut quo = vec![0i128; num.len() - den.len() + 1];
    for i in (0..quo.len()).rev() {
        let c = rem[i + den.len() - 1];
        if c != 0 {
            for (j, d) in den.iter().enumerate() {
                rem[i + j] -= c * d;
            }
        }
        quo[i] = c;
    }
    rem.iter().all(|&x| x == 0).then_some(quo)
}

/// Point counts `#C(F_{q^k})`, `k = 1..=kmax`, implied by a characteristic
/// polynomial: `q^k + 1 − Σ α_i^k`, power sums from Newton's identities.
pub fn counts_from_char_poly(p: &[i128], q: i64, kmax: usize) -> Vec<i128> {
    let deg = p.len() - 1;
    // e_i with P(T) = Σ (−1)^i e_i T^{deg−i}
    let e: Vec<i128> = (0..=deg)
        .map(|i| if i % 2 == 0 { p[deg - i] } else { -p[deg - i] })
        .collect();
    let mut pw = vec![0i128; kmax + 1];
    for k in 1..=kmax {
        let mut acc = 0i128;
        for i in 1..k.min(deg + 1) {
            let term = e[i] * pw[k - i];
            acc += if i % 2 == 1 { term } else { -term };
        }
        if k <= deg {
            let term = k as i128 * e[k];
            acc += if k % 2 == 1 { term } else { -term };
        }
        pw[k] = acc;
    }
    (1..=kmax)
        .map(|k| (q as i128).pow(k as u32) + 1 - pw[k])
        .collect()
}

#[derive(Clone, Copy, Debug)]
struct C64 {
    re: f64,
    im: f64,
}

impl C64 {
    fn add(self, o: C64) -> C64 {
        C64 {
            re: self.re + o.re,
            im: self.im + o.im,
        }
    }
    fn sub(self, o: C64) -> C64 {
        C64 {
            re: self.re - o.re,
            im: self.im - o.im,
        }
    }
    fn mul(self, o: C64) -> C64 {
        C64 {
            re: self.re * o.re - self.im * o.im,
            im: self.re * o.im + self.im * o.re,
        }
    }
    fn div(self, o: C64) -> C64 {
        let d = o.re * o.re + o.im * o.im;
        C64 {
            re: (self.re * o.re + self.im * o.im) / d,
            im: (self.im * o.re - self.re * o.im) / d,
        }
    }
    fn norm2(self) -> f64 {
        self.re * self.re + self.im * self.im
    }
}

/// Roots of a monic real polynomial by Durand–Kerner, as `|root|²` values.
fn root_norms(p: &[i128]) -> Vec<f64> {
    let deg = p.len() - 1;
    let coef: Vec<f64> = p.iter().map(|&c| c as f64).collect();
    let eval = |z: C64| {
        let mut acc = C64 { re: 0.0, im: 0.0 };
        for &c in coef.iter().rev() {
            acc = acc.mul(z).add(C64 { re: c, im: 0.0 });
        }
        acc
    };
    let mut roots: Vec<C64> = (0..deg)
        .map(|i| {
            let ang = 0.4 + 2.0 * std::f64::consts::PI * i as f64 / deg as f64;
            C64 {
                re: 1.4 * ang.cos(),
                im: 1.4 * ang.sin(),
            }
        })
        .collect();
    for _ in 0..4000 {
        let prev = roots.clone();
        for i in 0..deg {
            let mut den = C64 { re: 1.0, im: 0.0 };
            for j in 0..deg {
                if i != j {
                    den = den.mul(roots[i].sub(roots[j]));
                }
            }
            roots[i] = roots[i].sub(eval(roots[i]).div(den));
        }
        let moved: f64 = roots
            .iter()
            .zip(&prev)
            .map(|(a, b)| a.sub(*b).norm2())
            .sum();
        if moved < 1e-28 {
            break;
        }
    }
    roots.iter().map(|z| z.norm2()).collect()
}

/// `true` when every root has absolute value `√q` (to 1e-6), i.e. `p` is a
/// Weil `q`-polynomial.  Used as an independent check on computed curves.
pub fn is_weil_polynomial(p: &[i128], q: i64) -> bool {
    p.len() > 1
        && root_norms(p)
            .iter()
            .all(|&n| (n - q as f64).abs() < 1e-6 * q as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn elliptic_curve_round_trip() {
        // E_0 over F_2: 4 points, P(T) = T² + T + 2.
        let p = char_poly_from_counts(&[4], 2).unwrap();
        assert_eq!(p, vec![2, 1, 1]);
        assert!(is_weil_polynomial(&p, 2));
        assert_eq!(counts_from_char_poly(&p, 2, 3), vec![4, 8, 4]);
    }

    #[test]
    fn division_detects_factors() {
        // (T² + T + 2)(T² − T + 2) = T⁴ + 3T² + 4
        let prod = vec![4, 0, 3, 0, 1];
        assert_eq!(divide_exact(&prod, &[2, 1, 1]), Some(vec![2, -1, 1]));
        assert_eq!(divide_exact(&prod, &[2, 0, 1]), None);
    }
}

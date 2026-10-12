//! Hilbert class polynomials `H_D(X) ∈ Z[X]` by the complex-analytic method.
//!
//! `H_D(X) = prod (X - j(tau_f))` over the primitive reduced forms `f = (a, b, c)` of
//! discriminant `D < 0`, with `tau_f = (-b + i sqrt|D|) / (2a)`.  Each `j(tau)` is evaluated as
//! `j = E_4^3 / Delta` with
//!
//! * `1/q = exp(pi sqrt|D| / a) * (cos(pi b/a) + i sin(pi b/a))` computed directly, so no
//!   relative precision is lost to the small `|q|`,
//! * `E_4 = 1 + 240 sum sigma_3(n) q^n`,
//! * `Delta / q = (prod (1 - q^n))^24` by Euler's pentagonal-number series,
//!
//! all in fixed point with `P` fractional bits on [`Int`].  `P` is chosen from the size of the
//! roots, and the result is accepted only if every coefficient's real part lies within `2^-16`
//! of an integer and its imaginary part is below `2^-16`; otherwise an error names the worst
//! coefficient.  The rounding test is necessary, not sufficient: callers that use `H_D mod p`
//! should check what they build from it (a CM curve with the expected trace), which certifies
//! the polynomial for that use.
use crate::int::Int;

/// Primitive reduced positive-definite forms `(a, b, c)` with `b^2 - 4ac = d`:
/// `|b| <= a <= c`, `b >= 0` when `|b| = a` or `a = c`.  Their number is `h(d)`.
pub fn reduced_forms(d: i64) -> Vec<(i64, i64, i64)> {
    assert!(
        d < 0 && (d.rem_euclid(4) == 0 || d.rem_euclid(4) == 1),
        "not a discriminant"
    );
    let ad = -d;
    let mut out = vec![];
    let mut a = 1i64;
    while 3 * a * a <= ad {
        let mut b = -a + 1;
        while b <= a {
            if (b - d).rem_euclid(2) == 0 {
                let num = b * b - d;
                if num % (4 * a) == 0 {
                    let c = num / (4 * a);
                    if c >= a && !(b < 0 && (a == c || -b == a)) && gcd3(a, b, c) == 1 {
                        out.push((a, b, c));
                    }
                }
            }
            b += 1;
        }
        a += 1;
    }
    out
}

fn gcd(mut a: i64, mut b: i64) -> i64 {
    a = a.abs();
    b = b.abs();
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}
fn gcd3(a: i64, b: i64, c: i64) -> i64 {
    gcd(gcd(a, b), c)
}

/// Fixed-point arithmetic with `p` fractional bits.
#[derive(Clone, Copy)]
struct Fx {
    p: usize,
}

#[derive(Clone, Debug)]
struct C {
    re: Int,
    im: Int,
}

impl Fx {
    fn one(&self) -> Int {
        Int::one().shl(self.p)
    }
    fn mul(&self, a: &Int, b: &Int) -> Int {
        (a * b).shr(self.p)
    }
    /// a / b in fixed point.
    fn div(&self, a: &Int, b: &Int) -> Int {
        let (q, _) = a.shl(self.p).divrem_trunc(b);
        q
    }
    fn div_small(&self, a: &Int, k: i64) -> Int {
        let (q, _) = a.divrem_trunc(&Int::from(k));
        q
    }
    fn cmul(&self, x: &C, y: &C) -> C {
        let re = &self.mul(&x.re, &y.re) - &self.mul(&x.im, &y.im);
        let im = &self.mul(&x.re, &y.im) + &self.mul(&x.im, &y.re);
        C { re, im }
    }
    fn cdiv(&self, x: &C, y: &C) -> C {
        // x / y = x * conj(y) / |y|^2
        let n = &self.mul(&y.re, &y.re) + &self.mul(&y.im, &y.im);
        let t = self.cmul(
            x,
            &C {
                re: y.re.clone(),
                im: -&y.im,
            },
        );
        C {
            re: self.div(&t.re, &n),
            im: self.div(&t.im, &n),
        }
    }
    /// pi by Machin's formula: 16 atan(1/5) - 4 atan(1/239).
    fn pi(&self) -> Int {
        let a5 = self.atan_inv(5);
        let a239 = self.atan_inv(239);
        &(&a5 * &Int::from(16i64)) - &(&a239 * &Int::from(4i64))
    }
    fn atan_inv(&self, k: i64) -> Int {
        // atan(1/k) = sum (-1)^n / ((2n+1) k^(2n+1))
        let k2 = k * k;
        let mut term = self.div_small(&self.one(), k); // 1/k^(2n+1)
        let mut sum = Int::zero();
        let mut n = 0i64;
        while !term.is_zero() {
            let t = self.div_small(&term, 2 * n + 1);
            sum = if n % 2 == 0 { &sum + &t } else { &sum - &t };
            term = self.div_small(&term, k2);
            n += 1;
        }
        sum
    }
    fn sqrt_u64(&self, n: u64) -> Int {
        Int::from(n).shl(2 * self.p).isqrt()
    }
    /// exp(y) for real y (any sign) by halving to |y| < 1/2, Taylor, then squaring.
    fn exp(&self, y: &Int) -> Int {
        let mut k = 0usize;
        let mut r = y.clone();
        let half = self.one().shr(1);
        while r.abs() > half {
            r = r.shr(1);
            k += 1;
        }
        let mut sum = self.one();
        let mut term = self.one();
        let mut n = 1i64;
        loop {
            term = self.div_small(&self.mul(&term, &r), n);
            if term.is_zero() {
                break;
            }
            sum = &sum + &term;
            n += 1;
        }
        for _ in 0..k {
            sum = self.mul(&sum, &sum);
        }
        sum
    }
    /// (cos t, sin t) by halving to |t| < 1/2, Taylor, then double-angle steps.
    fn cos_sin(&self, t: &Int) -> (Int, Int) {
        let mut k = 0usize;
        let mut r = t.clone();
        let half = self.one().shr(1);
        while r.abs() > half {
            r = r.shr(1);
            k += 1;
        }
        let r2 = self.mul(&r, &r);
        let (mut c, mut s) = (self.one(), r.clone());
        let (mut tc, mut ts) = (self.one(), r.clone());
        let mut n = 1i64;
        loop {
            tc = -self.div_small(&self.mul(&tc, &r2), (2 * n - 1) * (2 * n));
            ts = -self.div_small(&self.mul(&ts, &r2), (2 * n) * (2 * n + 1));
            if tc.is_zero() && ts.is_zero() {
                break;
            }
            c = &c + &tc;
            s = &s + &ts;
            n += 1;
        }
        for _ in 0..k {
            let c2 = &self.mul(&c, &c) - &self.mul(&s, &s);
            let s2 = self.mul(&(&c + &c), &s);
            c = c2;
            s = s2;
        }
        (c, s)
    }
}

fn sigma3(n: u64) -> i64 {
    let mut s = 0i64;
    let mut d = 1u64;
    while d * d <= n {
        if n.is_multiple_of(d) {
            s += (d * d * d) as i64;
            let e = n / d;
            if e != d {
                s += (e * e * e) as i64;
            }
        }
        d += 1;
    }
    s
}

/// `j(tau)` for `tau = (-b + i sqrt|d|)/(2a)`, in fixed point (real, imaginary).
fn j_of_form(fx: &Fx, pi: &Int, sqrt_ad: &Int, a: i64, b: i64) -> C {
    // y = pi sqrt|d| / a,  theta = pi b / a
    let y = fx.div_small(&fx.mul(pi, sqrt_ad), a);
    let theta = fx.div_small(&(pi * &Int::from(b)), a);
    let (ct, st) = fx.cos_sin(&theta);
    let ey = fx.exp(&y);
    let eny = fx.div(&fx.one(), &ey);
    // q = e^{-y} (cos theta - i sin theta), 1/q = e^{y} (cos theta + i sin theta)
    let q = C {
        re: fx.mul(&eny, &ct),
        im: -fx.mul(&eny, &st),
    };
    let inv_q = C {
        re: fx.mul(&ey, &ct),
        im: fx.mul(&ey, &st),
    };
    // E4 = 1 + 240 sum sigma3(n) q^n ; eta-product prod(1 - q^n) by the pentagonal series
    let mut e4 = C {
        re: fx.one(),
        im: Int::zero(),
    };
    let mut qn = q.clone();
    let mut n = 1u64;
    while !(qn.re.is_zero() && qn.im.is_zero()) {
        let c = sigma3(n) * 240;
        let ci = Int::from(c);
        e4.re = &e4.re + &(&qn.re * &ci);
        e4.im = &e4.im + &(&qn.im * &ci);
        qn = fx.cmul(&qn, &q);
        n += 1;
    }
    // prod (1 - q^n) = sum_{k in Z} (-1)^k q^{k(3k-1)/2}
    let mut eta = C {
        re: fx.one(),
        im: Int::zero(),
    };
    let mut k = 1i64;
    loop {
        let mut any = false;
        for e in [k * (3 * k - 1) / 2, k * (3 * k + 1) / 2] {
            let qe = cpow(fx, &q, e as u64);
            if qe.re.is_zero() && qe.im.is_zero() {
                continue;
            }
            any = true;
            if k % 2 == 1 {
                eta.re = &eta.re - &qe.re;
                eta.im = &eta.im - &qe.im;
            } else {
                eta.re = &eta.re + &qe.re;
                eta.im = &eta.im + &qe.im;
            }
        }
        if !any {
            break;
        }
        k += 1;
    }
    let eta2 = fx.cmul(&eta, &eta);
    let eta4 = fx.cmul(&eta2, &eta2);
    let eta8 = fx.cmul(&eta4, &eta4);
    let eta16 = fx.cmul(&eta8, &eta8);
    let eta24 = fx.cmul(&eta16, &eta8);
    let e4_3 = fx.cmul(&fx.cmul(&e4, &e4), &e4);
    fx.cmul(&inv_q, &fx.cdiv(&e4_3, &eta24))
}

fn cpow(fx: &Fx, x: &C, mut e: u64) -> C {
    let mut r = C {
        re: fx.one(),
        im: Int::zero(),
    };
    let mut b = x.clone();
    while e > 0 {
        if e & 1 == 1 {
            r = fx.cmul(&r, &b);
        }
        e >>= 1;
        if e > 0 {
            b = fx.cmul(&b, &b);
        }
    }
    r
}

/// The result of [`hilbert_class_poly`].
#[derive(Clone, Debug)]
pub struct HilbertPoly {
    pub d: i64,
    /// coefficients low to high; monic of degree h(d)
    pub coeffs: Vec<Int>,
    pub forms: Vec<(i64, i64, i64)>,
    /// fractional bits used
    pub precision: usize,
    /// log2 of the worst distance of a coefficient from the integer it was rounded to
    /// (real part) or from zero (imaginary part)
    pub worst_error_log2: f64,
}

/// `H_d(X)` for a discriminant `d < 0` (fundamental or not; primitive forms are used).
pub fn hilbert_class_poly(d: i64) -> Result<HilbertPoly, String> {
    let forms = reduced_forms(d);
    let ad = (-d) as f64;
    // log2 |j(tau_f)| ~ pi sqrt|d| / (a ln 2); the coefficients are bounded by the
    // elementary symmetric functions of the roots
    let mut bits = 0f64;
    for &(a, _, _) in &forms {
        bits += (std::f64::consts::PI * ad.sqrt() / (a as f64) / std::f64::consts::LN_2).max(10.0)
            + 2.0;
    }
    let h = forms.len();
    let p = bits.ceil() as usize + 2 * (64 - (h as u64).leading_zeros() as usize) + 96;
    let fx = Fx { p };
    let pi = fx.pi();
    let sqrt_ad = fx.sqrt_u64((-d) as u64);
    let mut poly: Vec<C> = vec![C {
        re: fx.one(),
        im: Int::zero(),
    }];
    for &(a, b, _) in &forms {
        let j = j_of_form(&fx, &pi, &sqrt_ad, a, b);
        // poly *= (X - j)
        let mut next = vec![
            C {
                re: Int::zero(),
                im: Int::zero()
            };
            poly.len() + 1
        ];
        for (k, c) in poly.iter().enumerate() {
            next[k + 1].re = &next[k + 1].re + &c.re;
            next[k + 1].im = &next[k + 1].im + &c.im;
            let t = fx.cmul(c, &j);
            next[k].re = &next[k].re - &t.re;
            next[k].im = &next[k].im - &t.im;
        }
        poly = next;
    }
    let one = fx.one();
    let half = one.shr(1);
    let mut coeffs = vec![];
    let mut worst = f64::NEG_INFINITY;
    for (k, c) in poly.iter().enumerate() {
        // round re / 2^p to the nearest integer
        let shifted = &c.re + &(if c.re.is_neg() { -&half } else { half.clone() });
        let n = shifted.shr(p);
        let err_re = (&c.re - &n.shl(p)).abs();
        let err_im = c.im.abs();
        let e = err_re.max(err_im);
        let el = if e.is_zero() {
            -(p as f64)
        } else {
            e.bits() as f64 - p as f64
        };
        worst = worst.max(el);
        if el > -16.0 {
            return Err(format!(
                "H_{d}: coefficient {k} not within 2^-16 of an integer (log2 error {el:.1}, precision {p} bits)"
            ));
        }
        coeffs.push(n);
    }
    Ok(HilbertPoly {
        d,
        coeffs,
        forms,
        precision: p,
        worst_error_log2: worst,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn class_numbers_from_forms() {
        // h(D) for a few fundamental discriminants
        for (d, h) in [
            (-3, 1),
            (-4, 1),
            (-7, 1),
            (-23, 3),
            (-47, 5),
            (-71, 7),
            (-619, 5),
            (-3019, 7),
        ] {
            assert_eq!(reduced_forms(d).len(), h, "h({d})");
        }
    }

    #[test]
    fn small_hilbert_polys_match_known_values() {
        // H_{-3} = X, H_{-4} = X - 1728, H_{-7} = X + 3375, H_{-8} = X - 8000
        for (d, c0) in [
            (-3i64, 0i64),
            (-4, -1728),
            (-7, 3375),
            (-8, -8000),
            (-11, 32768),
            (-163, 262537412640768000),
        ] {
            let hp = hilbert_class_poly(d).unwrap();
            assert_eq!(hp.coeffs.len(), 2);
            assert_eq!(hp.coeffs[1], Int::one());
            assert_eq!(hp.coeffs[0], Int::from(c0), "H_{d}");
        }
        // H_{-15} = X^2 + 191025 X - 121287375
        let hp = hilbert_class_poly(-15).unwrap();
        let want = [-121287375i64, 191025, 1];
        for (c, w) in hp.coeffs.iter().zip(want) {
            assert_eq!(*c, Int::from(w));
        }
        // H_{-23} = X^3 + 3491750 X^2 - 5151296875 X + 12771880859375
        let hp = hilbert_class_poly(-23).unwrap();
        let want = [12771880859375i64, -5151296875, 3491750, 1];
        for (c, w) in hp.coeffs.iter().zip(want) {
            assert_eq!(*c, Int::from(w));
        }
    }
}

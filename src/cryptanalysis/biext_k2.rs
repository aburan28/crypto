//! Embedding-degree-2 pairings by a cubical (biextension) ladder, counted
//! against Miller, on the repository's `k = 2` path
//! ([`crate::cryptanalysis::eds_mov`]: supersingular `y² = x³ + x` over
//! `F_p`, `p ≡ 3 (mod 4)`, distortion map `φ(x, y) = (−x, iy)`).
//!
//! lit2026 W11, after Zheng–Lin–Zhao, ePrint 2025/1652 ("Pairings with
//! embedding degree two via biextensions"), whose abstract reports 12–18%
//! fewer `F_p` multiplications than Miller.  The paper's PDF could not be
//! retrieved in this session, so the ladder here is our own construction
//! from cubical arithmetic on the Kummer line (Robert's formulation: the
//! non-reduced Tate pairing is the ratio `λ_{P,Q}/λ_P` read off
//! `ℓP̃ = λ_P·0̃` and `(ℓP+Q)~ = λ_{P,Q}·Q̃`), specialised to `k = 2`.
//!
//! Three evaluators of the reduced Tate pairing `t_r(P, φ(R))`, all on one
//! counting field so their `F_p` operation counts are comparable:
//!
//! - [`miller_affine`]: the repository's affine Miller loop of
//!   `eds_mov::miller_fp2` (an inversion per step, verticals kept), ported
//!   to the counting field;
//! - [`miller_jacobian`]: denominator-free Miller with `T` in Jacobian
//!   coordinates.  `x(φ(R)) = −x(R) ∈ F_p`, so every vertical line and
//!   every `F_p` scale factor is killed by the final exponentiation
//!   `(p² − 1)/r`; this is the usual optimised baseline;
//! - [`cubical_ladder`]: the three-point Montgomery ladder on
//!   `(nP, (n+1)P, nP + Q)` with cubical differential additions (divide by
//!   the difference's coordinates instead of cross-multiplying).  `P`,
//!   `Q = φ(R)` have `F_p` abscissae, so only the `nP + Q` lane is in
//!   `F_{p²}`.  After the final exponentiation `Z(rP + Q)` gives the
//!   **square** of the reduced Tate pairing (measured on every instance), so
//!   [`unsquare`] raises it to `(r+1)/2`; its cost is counted separately.
//!
//! All three are final-exponentiated and compared with each other and, at
//! toy sizes, with [`crate::cryptanalysis::eds_mov::modified_tate`].
//! Multiplication counts are for the loop only (the final exponentiation is
//! the same for all and is reported separately).

use std::cell::Cell;

/// Counted operations in `F_p`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct FpCounts {
    pub mul: u64,
    pub sqr: u64,
    pub inv: u64,
    pub add: u64,
}

impl FpCounts {
    /// `M + S` (squarings at full price) — the "F_p mults" figure.
    pub fn m_plus_s(&self) -> u64 {
        self.mul + self.sqr
    }
}

/// `F_p` with `p < 2^63`, operation-counted.
pub struct Fp {
    pub p: u64,
    c: Cell<FpCounts>,
}

impl Fp {
    pub fn new(p: u64) -> Self {
        assert!(p % 4 == 3 && p < (1u64 << 63));
        Self {
            p,
            c: Cell::new(FpCounts::default()),
        }
    }
    pub fn counts(&self) -> FpCounts {
        self.c.get()
    }
    pub fn reset(&self) {
        self.c.set(FpCounts::default());
    }
    fn bump(&self, g: impl FnOnce(&mut FpCounts)) {
        let mut c = self.c.get();
        g(&mut c);
        self.c.set(c);
    }
    #[inline]
    fn raw_mul(&self, a: u64, b: u64) -> u64 {
        ((a as u128 * b as u128) % self.p as u128) as u64
    }
    pub fn mul(&self, a: u64, b: u64) -> u64 {
        self.bump(|c| c.mul += 1);
        self.raw_mul(a, b)
    }
    pub fn sqr(&self, a: u64) -> u64 {
        self.bump(|c| c.sqr += 1);
        self.raw_mul(a, a)
    }
    pub fn add(&self, a: u64, b: u64) -> u64 {
        self.bump(|c| c.add += 1);
        let s = a as u128 + b as u128;
        (s % self.p as u128) as u64
    }
    pub fn sub(&self, a: u64, b: u64) -> u64 {
        self.bump(|c| c.add += 1);
        if a >= b {
            a - b
        } else {
            self.p - (b - a)
        }
    }
    pub fn neg(&self, a: u64) -> u64 {
        if a == 0 {
            0
        } else {
            self.p - a
        }
    }
    fn raw_pow(&self, mut a: u64, mut e: u64) -> u64 {
        let mut r = 1u64;
        while e > 0 {
            if e & 1 == 1 {
                r = self.raw_mul(r, a);
            }
            a = self.raw_mul(a, a);
            e >>= 1;
        }
        r
    }
    pub fn inv(&self, a: u64) -> u64 {
        self.bump(|c| c.inv += 1);
        self.raw_pow(a, self.p - 2)
    }
}

/// `a + b·i`, `i² = −1`.
pub type Fp2 = (u64, u64);

impl Fp {
    /// Karatsuba: `3M`.
    pub fn f2_mul(&self, x: Fp2, y: Fp2) -> Fp2 {
        let ac = self.mul(x.0, y.0);
        let bd = self.mul(x.1, y.1);
        let t = self.mul(self.add(x.0, x.1), self.add(y.0, y.1));
        (self.sub(ac, bd), self.sub(self.sub(t, ac), bd))
    }
    /// Complex squaring: `2M`.
    pub fn f2_sqr(&self, x: Fp2) -> Fp2 {
        let re = self.mul(self.add(x.0, x.1), self.sub(x.0, x.1));
        let ab = self.mul(x.0, x.1);
        (re, self.add(ab, ab))
    }
    /// `F_p × F_{p²}`: `2M`.
    pub fn f2_scale(&self, s: u64, x: Fp2) -> Fp2 {
        (self.mul(s, x.0), self.mul(s, x.1))
    }
    pub fn f2_add(&self, x: Fp2, y: Fp2) -> Fp2 {
        (self.add(x.0, y.0), self.add(x.1, y.1))
    }
    pub fn f2_sub(&self, x: Fp2, y: Fp2) -> Fp2 {
        (self.sub(x.0, y.0), self.sub(x.1, y.1))
    }
    /// `1/x` by the norm: `2S + 2M + I`.
    pub fn f2_inv(&self, x: Fp2) -> Fp2 {
        let n = self.add(self.sqr(x.0), self.sqr(x.1));
        let ni = self.inv(n);
        (self.mul(x.0, ni), self.mul(self.neg(x.1), ni))
    }
    /// `x^e`, uncounted (the final exponentiation, common to all three).
    pub fn f2_pow_uncounted(&self, mut x: Fp2, mut e: u128) -> Fp2 {
        let mm = |a: Fp2, b: Fp2| -> Fp2 {
            let re = (self.raw_mul(a.0, b.0) as u128 + self.p as u128
                - self.raw_mul(a.1, b.1) as u128)
                % self.p as u128;
            let im =
                (self.raw_mul(a.0, b.1) as u128 + self.raw_mul(a.1, b.0) as u128) % self.p as u128;
            (re as u64, im as u64)
        };
        let mut r = (1u64, 0u64);
        while e > 0 {
            if e & 1 == 1 {
                r = mm(r, x);
            }
            x = mm(x, x);
            e >>= 1;
        }
        r
    }
    /// The reduced pairing value: `f^{(p²−1)/r}`.
    pub fn final_exp(&self, f: Fp2, r: u64) -> Fp2 {
        let e = (self.p as u128 * self.p as u128 - 1) / r as u128;
        self.f2_pow_uncounted(f, e)
    }
}

/// Affine point on `y² = x³ + x` over `F_p` (`None` = `O`), uncounted
/// helpers for set-up.
pub type Pt = Option<(u64, u64)>;

pub fn ec_add(f: &Fp, a: Pt, b: Pt) -> Pt {
    let (x1, y1) = a?;
    let _ = b?;
    let (x2, y2) = b.unwrap();
    let p = f.p;
    let lam = if x1 == x2 {
        if (y1 as u128 + y2 as u128).is_multiple_of(p as u128) {
            return None;
        }
        let num = (3 * f.raw_mul(x1, x1) as u128 + 1) % p as u128;
        f.raw_mul(num as u64, f.raw_pow(f.raw_mul(2, y1), p - 2))
    } else {
        let num = (y2 + p - y1) % p;
        let den = (x2 + p - x1) % p;
        f.raw_mul(num, f.raw_pow(den, p - 2))
    };
    let x3 = ((f.raw_mul(lam, lam) as u128 + 2 * p as u128 - x1 as u128 - x2 as u128) % p as u128)
        as u64;
    let y3 =
        ((f.raw_mul(lam, (x1 + p - x3) % p) as u128 + p as u128 - y1 as u128) % p as u128) as u64;
    Some((x3, y3))
}

fn ec_add_opt(f: &Fp, a: Pt, b: Pt) -> Pt {
    match (a, b) {
        (None, _) => b,
        (_, None) => a,
        _ => ec_add(f, a, b),
    }
}

pub fn ec_mul(f: &Fp, pt: Pt, mut k: u128) -> Pt {
    let mut acc = None;
    let mut base = pt;
    while k > 0 {
        if k & 1 == 1 {
            acc = ec_add_opt(f, acc, base);
        }
        base = ec_add_opt(f, base, base);
        k >>= 1;
    }
    acc
}

/// The deterministic point of order `r` with the smallest `x`:
/// `[(p+1)/r]·(x, √(x³+x))`.
pub fn point_of_order(f: &Fp, r: u64) -> Option<(u64, u64)> {
    let p = f.p;
    let c = (p as u128 + 1) / r as u128;
    for x in 1..p.min(1_000_000) {
        let rhs = (f.raw_mul(f.raw_mul(x, x), x) as u128 + x as u128) % p as u128;
        let rhs = rhs as u64;
        let y = f.raw_pow(rhs, (p + 1) / 4);
        if f.raw_mul(y, y) != rhs || y == 0 {
            continue;
        }
        if let Some(q) = ec_mul(f, Some((x, y)), c) {
            if ec_mul(f, Some(q), r as u128).is_none() {
                return Some(q);
            }
        }
    }
    None
}

/// The Miller loop of `eds_mov::miller_fp2`, counted: `f_{r,P}(φ(R))`
/// with affine `T`, one `F_p` inversion and one `F_{p²}` inversion per
/// step.
pub fn miller_affine(f: &Fp, pp: (u64, u64), rr: (u64, u64), r: u64) -> Option<Fp2> {
    let rx: Fp2 = (f.neg(rr.0), 0);
    let ry: Fp2 = (0, rr.1);
    let line = |t: (u64, u64), lam: u64, x3: u64| -> Option<Fp2> {
        let l = f.f2_sub(
            f.f2_sub(ry, (t.1, 0)),
            f.f2_scale(lam, f.f2_sub(rx, (t.0, 0))),
        );
        let v = f.f2_sub(rx, (x3, 0));
        if v == (0, 0) {
            return None;
        }
        Some(f.f2_mul(l, f.f2_inv(v)))
    };
    let dbl = |t: (u64, u64)| -> Option<(Fp2, Pt)> {
        let (x1, y1) = t;
        if y1 == 0 {
            return Some((f.f2_sub(rx, (x1, 0)), None));
        }
        let x1s = f.sqr(x1);
        let num = f.add(f.add(f.add(x1s, x1s), x1s), 1);
        let lam = f.mul(num, f.inv(f.add(y1, y1)));
        let x3 = f.sub(f.sub(f.sqr(lam), x1), x1);
        let y3 = f.sub(f.mul(lam, f.sub(x1, x3)), y1);
        Some((line(t, lam, x3)?, Some((x3, y3))))
    };
    let mut acc: Fp2 = (1, 0);
    let mut t: Pt = Some(pp);
    let nbits = 64 - r.leading_zeros();
    for i in (0..nbits - 1).rev() {
        let (g, t2) = dbl(t?)?;
        acc = f.f2_mul(f.f2_sqr(acc), g);
        t = t2;
        if (r >> i) & 1 == 1 {
            let (g, t3) = match t {
                None => ((1, 0), Some(pp)),
                Some((x1, y1)) => {
                    let (x2, y2) = pp;
                    if x1 == x2 {
                        if (y1 as u128 + y2 as u128).is_multiple_of(f.p as u128) {
                            (f.f2_sub(rx, (x2, 0)), None)
                        } else {
                            dbl((x1, y1))?
                        }
                    } else {
                        let lam = f.mul(f.sub(y2, y1), f.inv(f.sub(x2, x1)));
                        let x3 = f.sub(f.sub(f.sqr(lam), x1), x2);
                        let y3 = f.sub(f.mul(lam, f.sub(x1, x3)), y1);
                        (line((x1, y1), lam, x3)?, Some((x3, y3)))
                    }
                }
            };
            acc = f.f2_mul(acc, g);
            t = t3;
        }
    }
    Some(acc)
}

/// Denominator-free Miller: `T = (X:Y:Z)` Jacobian, verticals and `F_p`
/// factors dropped (the final exponentiation kills them because
/// `x(φ(R)) ∈ F_p`).
pub fn miller_jacobian(f: &Fp, pp: (u64, u64), rr: (u64, u64), r: u64) -> Fp2 {
    let xq = f.neg(rr.0); // x(φ(R)) ∈ F_p
    let yq = rr.1; // y(φ(R)) = yq·i
    let (xp, yp) = pp;
    let (mut x, mut y, mut z) = (xp, yp, 1u64);
    let mut inf = false;
    let mut acc: Fp2 = (1, 0);
    let nbits = 64 - r.leading_zeros();
    for i in (0..nbits - 1).rev() {
        // Doubling with tangent line.
        acc = f.f2_sqr(acc);
        if !inf {
            if y == 0 {
                inf = true; // vertical tangent: an F_p value, dropped
            } else {
                let xx = f.sqr(x);
                let yy = f.sqr(y);
                let yyyy = f.sqr(yy);
                let zz = f.sqr(z);
                // S = 4·X·Y² = 2((X + Y²)² − X² − Y⁴)
                let s = f.sub(f.sub(f.sqr(f.add(x, yy)), xx), yyyy);
                let s = f.add(s, s);
                // M = 3X² + a·Z⁴ with a = 1
                let m = f.add(f.add(f.add(xx, xx), xx), f.sqr(zz));
                let x3 = f.sub(f.sqr(m), f.add(s, s));
                let y8 = {
                    let t = f.add(yyyy, yyyy);
                    let t = f.add(t, t);
                    f.add(t, t)
                };
                let y3 = f.sub(f.mul(m, f.sub(s, x3)), y8);
                let z3 = f.sub(f.sub(f.sqr(f.add(y, z)), yy), zz);
                // l̃ = −Z3·Y − M·Z·(ZZ·xq − X) + (Z3·Z·ZZ·yq)·i
                let re = f.add(f.mul(z3, y), f.mul(f.mul(m, z), f.sub(f.mul(zz, xq), x)));
                let im = f.mul(f.mul(z3, f.mul(z, zz)), yq);
                let l = (f.neg(re), im);
                acc = f.f2_mul(acc, l);
                x = x3;
                y = y3;
                z = z3;
            }
        }
        if (r >> i) & 1 == 1 {
            if inf {
                // O + P = P, line is the constant 1.
                x = xp;
                y = yp;
                z = 1;
                inf = false;
                continue;
            }
            // Mixed addition T + P with the chord through them.
            let z1z1 = f.sqr(z);
            let u2 = f.mul(xp, z1z1);
            let s2 = f.mul(f.mul(yp, z), z1z1);
            let h = f.sub(u2, x);
            let rr_ = f.sub(s2, y);
            if h == 0 {
                if rr_ == 0 {
                    // T = P: tangent; cannot occur for r prime > 2 inside
                    // the loop, kept for completeness via a doubling.
                    unreachable!("T = P inside the Miller loop");
                }
                // T = −P: vertical line, F_p-valued, dropped; T + P = O.
                inf = true;
                continue;
            }
            let hh = f.sqr(h);
            let hhh = f.mul(h, hh);
            let v = f.mul(x, hh);
            let x3 = f.sub(f.sub(f.sqr(rr_), hhh), f.add(v, v));
            let y3 = f.sub(f.mul(rr_, f.sub(v, x3)), f.mul(y, hhh));
            let z3 = f.mul(z, h);
            // l̃ = Z3·yq·i − Z3·yP − r·(xq − xP)
            let re = f.add(f.mul(z3, yp), f.mul(rr_, f.sub(xq, xp)));
            let im = f.mul(z3, yq);
            acc = f.f2_mul(acc, (f.neg(re), im));
            x = x3;
            y = y3;
            z = z3;
        }
    }
    acc
}

/// One cubical differential addition on the `F_p` lane:
/// `x̃(A+B) = ((U+V)²/X_d, (U−V)²/Z_d)` with `U = (X_A−Z_A)(X_B+Z_B)`,
/// `V = (X_A+Z_A)(X_B−Z_B)`, given `1/X_d` and `Z_d = 1`.
fn cadd_fp(f: &Fp, a: (u64, u64), b: (u64, u64), inv_xd: u64) -> (u64, u64) {
    let u = f.mul(f.sub(a.0, a.1), f.add(b.0, b.1));
    let v = f.mul(f.add(a.0, a.1), f.sub(b.0, b.1));
    let x = f.mul(f.sqr(f.add(u, v)), inv_xd);
    let z = f.sqr(f.sub(u, v));
    (x, z)
}

/// Cubical doubling (the Montgomery `xDBL` with `A = 0`, unnormalised):
/// `X' = (X+Z)²(X−Z)²`, `Z' = 4XZ·((X−Z)² + 4XZ/2)`.
fn cdbl_fp(f: &Fp, a: (u64, u64), half: u64) -> (u64, u64) {
    let t1 = f.sqr(f.add(a.0, a.1));
    let t2 = f.sqr(f.sub(a.0, a.1));
    let t3 = f.sub(t1, t2);
    let x = f.mul(t1, t2);
    let z = f.mul(t3, f.add(t2, f.mul(half, t3)));
    (x, z)
}

/// The difference of a mixed cubical addition, by its inverted
/// coordinates: `(1/X_d, 1/Z_d)`.
#[derive(Clone, Copy)]
enum Diff {
    /// `X_d ∈ F_p`, `Z_d = 1`.
    Fp(u64),
    /// `X_d ∈ F_{p²}`, `Z_d ∈ F_p`.
    Fp2(Fp2, u64),
}

/// The `F_{p²}` lane: `x̃(A + B)` with `A` on the `F_p` lane and `B` on
/// the `F_{p²}` lane, `X = (U+V)²/X_d`, `Z = (U−V)²/Z_d`.
fn cadd_mixed(f: &Fp, a: (u64, u64), b: (Fp2, Fp2), d: Diff) -> (Fp2, Fp2) {
    let u = f.f2_scale(f.sub(a.0, a.1), f.f2_add(b.0, b.1));
    let v = f.f2_scale(f.add(a.0, a.1), f.f2_sub(b.0, b.1));
    let s = f.f2_sqr(f.f2_add(u, v));
    let t = f.f2_sqr(f.f2_sub(u, v));
    match d {
        Diff::Fp(ix) => (f.f2_scale(ix, s), t),
        Diff::Fp2(ix, iz) => (f.f2_mul(ix, s), f.f2_scale(iz, t)),
    }
}

/// The cubical ladder for `t_r(P, φ(R))`, before the final exponentiation.
///
/// Lifts `P̃ = (x_P : 1)`, `Q̃ = (x_Q : 1)`, `(P+Q)~ = (x_{P+Q} : 1)` with
/// `Q = φ(R)`; the lift of `P − Q` is then fixed by the cubical relation
/// `cADD(P̃, Q̃, (P−Q)~) = (P+Q)~`, i.e. `X_{P−Q} = (U+V)²/X_{P+Q}`,
/// `Z_{P−Q} = (U−V)²` (an `F_p` value).  The ladder keeps
/// `(nP̃, (n+1)P̃, (nP+Q)~)`:
///
/// - bit 0: `2nP + Q = nP + (nP+Q)` (difference `−Q`, `x̃(−Q) = Q̃`);
/// - bit 1: `(2n+1)P + Q = (n+1)P + (nP+Q)` (difference `P − Q`).
///
/// At the end `rP̃ = (λ_P : 0)` with `λ_P ∈ F_p` and
/// `(rP+Q)~ = λ_{P,Q}·Q̃`, so `λ_{P,Q} = Z(rP+Q)`, and the pairing
/// `λ_{P,Q}/λ_P` equals `Z(rP+Q)` up to `F_p^*`, which the final
/// exponentiation removes.  `None` when `rP̃` is not `(·:0)`.
pub fn cubical_ladder(f: &Fp, pp: (u64, u64), rr: (u64, u64), r: u64) -> Option<Fp2> {
    let xq = f.neg(rr.0);
    let (xp, yp) = pp;
    // x(P + Q) in F_{p²}, Q = (xq, yq·i): one affine chord.
    let dx = f.sub(xq, xp); // ∈ F_p, nonzero for independent P, Q
    let dy: Fp2 = (f.neg(yp), rr.1); // yq·i − yp
                                     // The cubical lift of P − Q: X_{P−Q} = (U+V)²/X_{P+Q}, Z_{P−Q} = (U−V)²,
                                     // so 1/X_{P−Q} = X_{P+Q}/(U+V)² needs only F_p inverses.
    let u = f.mul(f.sub(xp, 1), f.add(xq, 1));
    let v = f.mul(f.add(xp, 1), f.sub(xq, 1));
    let s_pmq = f.sqr(f.add(u, v));
    let z_pmq = f.sqr(f.sub(u, v));
    // One inversion for the five F_p denominators (Montgomery's trick).
    let dens = [dx, s_pmq, z_pmq, xq, xp];
    let mut pre = [0u64; 5];
    let mut acc = dens[0];
    pre[0] = acc;
    for i in 1..5 {
        acc = f.mul(acc, dens[i]);
        pre[i] = acc;
    }
    let mut inv_acc = f.inv(acc);
    let mut invs = [0u64; 5];
    for i in (1..5).rev() {
        invs[i] = f.mul(inv_acc, pre[i - 1]);
        inv_acc = f.mul(inv_acc, dens[i]);
    }
    invs[0] = inv_acc;
    let [inv_dx, inv_s, inv_z, inv_xq, inv_xp] = invs;
    let lam = f.f2_scale(inv_dx, dy);
    let xpq: Fp2 = f.f2_sub(f.f2_sub(f.f2_sqr(lam), (xp, 0)), (xq, 0));
    let d_pmq = Diff::Fp2(f.f2_scale(inv_s, xpq), inv_z);
    let d_q = Diff::Fp(inv_xq);
    let half = f.p.div_ceil(2); // 1/2 mod p, no inversion

    let mut a = (xp, 1u64); // nP
    let mut b = cdbl_fp(f, a, half); // (n+1)P, n = 1
    let mut c: (Fp2, Fp2) = (xpq, (1, 0)); // nP + Q
    let nbits = 64 - r.leading_zeros();
    for i in (0..nbits - 1).rev() {
        if (r >> i) & 1 == 0 {
            c = cadd_mixed(f, a, c, d_q);
            b = cadd_fp(f, b, a, inv_xp);
            a = cdbl_fp(f, a, half);
        } else {
            c = cadd_mixed(f, b, c, d_pmq);
            a = cadd_fp(f, b, a, inv_xp);
            b = cdbl_fp(f, b, half);
        }
    }
    if a.1 != 0 {
        return None;
    }
    Some(c.1)
}

/// `v^{(r+1)/2}` for `v` of norm 1 (a reduced pairing value): turns the
/// cubical ladder's `t_r(P,φ(R))²` into `t_r(P,φ(R))` (`r` odd).
/// Cyclotomic squaring `(a, b)² = (2a² − 1, 2ab)` (`1S + 1M`), general
/// multiplication `3M`; counted.
pub fn unsquare(f: &Fp, v: Fp2, r: u64) -> Fp2 {
    let e = r.div_ceil(2); // (r+1)/2, r odd
    let mut acc = v;
    let top = 63 - e.leading_zeros();
    for i in (0..top).rev() {
        let a2 = f.sqr(acc.0);
        let ab = f.mul(acc.0, acc.1);
        acc = (f.sub(f.add(a2, a2), 1), f.add(ab, ab));
        if (e >> i) & 1 == 1 {
            acc = f.f2_mul(acc, v);
        }
    }
    acc
}

/// One pairing evaluated three ways, with the counts.
#[derive(Clone, Debug)]
pub struct PairingComparison {
    pub p: u64,
    pub r: u64,
    pub r_bits: u32,
    pub value_affine: Fp2,
    pub value_jacobian: Fp2,
    pub value_cubical: Fp2,
    pub ops_affine: FpCounts,
    pub ops_jacobian: FpCounts,
    pub ops_cubical: FpCounts,
    /// The reduced cubical value is the **square** of the Tate pairing
    /// (measured: `value_cubical_raw = value_affine²` on every instance);
    /// this is the raw value before [`unsquare`].
    pub value_cubical_raw: Fp2,
    /// Cost of [`unsquare`] (cyclotomic exponentiation by `(r+1)/2`).
    pub ops_unsquare: FpCounts,
}

/// Evaluate `t_r(P, φ(R))` by all three loops and final-exponentiate.
pub fn compare(p: u64, r: u64, pp: (u64, u64), rr: (u64, u64)) -> Option<PairingComparison> {
    let f = Fp::new(p);
    f.reset();
    let ma = miller_affine(&f, pp, rr, r)?;
    let ops_affine = f.counts();
    f.reset();
    let mj = miller_jacobian(&f, pp, rr, r);
    let ops_jacobian = f.counts();
    f.reset();
    let cl = cubical_ladder(&f, pp, rr, r)?;
    let ops_cubical = f.counts();
    let raw = f.final_exp(cl, r);
    f.reset();
    let unsq = unsquare(&f, raw, r);
    let ops_unsquare = f.counts();
    Some(PairingComparison {
        p,
        r,
        r_bits: 64 - r.leading_zeros(),
        value_affine: f.final_exp(ma, r),
        value_jacobian: f.final_exp(mj, r),
        value_cubical: unsq,
        ops_affine,
        ops_jacobian,
        ops_cubical,
        value_cubical_raw: raw,
        ops_unsquare,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::eds_mov;

    #[test]
    fn three_loops_agree_with_eds_mov_at_toy_size() {
        let (p, r) = (1283u64, 107u64);
        let f = Fp::new(p);
        let pp = point_of_order(&f, r).unwrap();
        let want = eds_mov::modified_tate(pp, pp, r, 1, p).unwrap();
        let c = compare(p, r, pp, pp).unwrap();
        assert_eq!(c.value_affine, want, "affine Miller port");
        assert_eq!(c.value_jacobian, want, "Jacobian Miller");
        assert_eq!(c.value_cubical, want, "cubical ladder, unsquared");
        assert_eq!(
            c.value_cubical_raw,
            f.f2_mul(want, want),
            "raw cubical = t²"
        );
        // and on pairs (kP, φ(P)), (P, φ(kP))
        for k in [2u128, 5, 33, 106] {
            let kp = ec_mul(&f, Some(pp), k).unwrap();
            let want = eds_mov::modified_tate(kp, pp, r, 1, p).unwrap();
            let c = compare(p, r, kp, pp).unwrap();
            assert_eq!(
                (c.value_affine, c.value_jacobian, c.value_cubical),
                (want, want, want),
                "k={k}"
            );
            let want = eds_mov::modified_tate(pp, kp, r, 1, p).unwrap();
            let c = compare(p, r, pp, kp).unwrap();
            assert_eq!(
                (c.value_affine, c.value_jacobian, c.value_cubical),
                (want, want, want),
                "phi k={k}"
            );
        }
    }
}

//! Cyclotomic factorisation of the NIST P-192 / P-224 / P-256 primes and
//! what it does, and does not do, to the **point group** and to the
//! **heights of lifted points**.
//!
//! Three questions, answered natively and written as a Markdown report plus
//! a JSON record:
//!
//! 1. **Algebraic (cyclotomic) factorisation.**  Each prime is `f(2^32)` for
//!    a sparse Solinas polynomial `f`.  `f ± 1` factors in `Z[t]`, and every
//!    `Φ_d(t)` factor forces a small prime factor of `p ± 1`.  The identities
//!    are verified by exact polynomial expansion and sealed as
//!    `identity.certificate/v1` records.  The full integer factorisations of
//!    `p ± 1`, `n ± 1`, the twist order, the trace and the CM discriminant
//!    are frozen inputs (found with PARI 2.17, see the report) that this
//!    tool *verifies*: product check plus Miller–Rabin on every factor, and
//!    an independent native trial-division/rho pass on each value.
//! 2. **Point-group consequences.**  The embedding degree `ord_n(p)`, the
//!    Cheon / Brown–Gallant cost as a function of the number of published
//!    powers `[α^i]G` (which is where the factorisation of `n ∓ 1` acts on
//!    the points), and the cyclotomic-norm factorisation
//!    `#E(F_{p^k}) = ∏_{d|k} N(Φ_d(π))`.
//! 3. **Point heights on lifts to `Q`.**  Integer and small-denominator
//!    rational points on the canonical lift `y² = x³ − 3x + b`, the
//!    canonical height `ĥ` of the generator lifted xedni-style to a curve
//!    through it, and the height pairing / regulator of a two-point lift.
//!    `ĥ` is computed as `lim h_x([2^k]P) / 4^k` (the `x`-coordinate
//!    normalisation, twice Silverman's AEC normalisation), with the last
//!    two terms of the sequence reported as the convergence check.
//!
//! ```bash
//! cargo run --release --example nist_cyclotomic_point_heights -- \
//!     [--out research/nist_cyclotomic_point_heights_20261011] \
//!     [--scan-x 1000000] [--scan-d 40] [--scan-u 100000] [--doublings 6]
//! ```
//!
//! Public parameters only; nothing here is a key or an attack run.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

use crypto_lib::cryptanalysis::curve_traits::arith::{
    self, is_prime, log2_abs, lucas_v, multiplicative_order, Factored,
};
use crypto_lib::cryptanalysis::identity_certificate::{issue, Expr, IdentityStatement};
use num_bigint::{BigInt, BigUint};
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};
use serde_json::{json, Value};

const LN2: f64 = std::f64::consts::LN_2;

// ────────────────────────────── curve data ──────────────────────────────

struct Curve {
    name: &'static str,
    /// Solinas polynomial, coefficient of `t^i` at index `i`, `p = f(2^32)`.
    f: &'static [i128],
    /// Algebraic factorisation of `f − 1` and `f + 1`: `(label, coeffs, exponent)`.
    f_minus_1: &'static [(&'static str, &'static [i128], u32)],
    f_plus_1: &'static [(&'static str, &'static [i128], u32)],
    p: &'static str,
    n: &'static str,
    b_hex: &'static str,
    gx_hex: &'static str,
    gy_hex: &'static str,
    /// Frozen complete factorisations (decimal primes with exponents).
    frozen: &'static [(&'static str, &'static [(&'static str, u32)])],
}

// Φ_1 = t − 1, Φ_2 = t + 1, Φ_4 = t² + 1.
const PHI1: &[i128] = &[-1, 1];
const PHI2: &[i128] = &[1, 1];
const PHI4: &[i128] = &[1, 0, 1];
const T: &[i128] = &[0, 1];

const CURVES: &[Curve] = &[
    Curve {
        name: "P-192",
        f: &[-1, 0, -1, 0, 0, 0, 1], // t^6 − t^2 − 1
        f_minus_1: &[("t^6 − t^2 − 2 (irreducible)", &[-2, 0, -1, 0, 0, 0, 1], 1)],
        f_plus_1: &[
            ("t", T, 2),
            ("Φ_1(t)", PHI1, 1),
            ("Φ_2(t)", PHI2, 1),
            ("Φ_4(t)", PHI4, 1),
        ],
        p: "6277101735386680763835789423207666416083908700390324961279",
        n: "6277101735386680763835789423176059013767194773182842284081",
        b_hex: "64210519e59c80e70fa7e9ab72243049feb8deecc146b9b1",
        gx_hex: "188da80eb03090f67cbf20eb43a18800f4ff0afd82ff1012",
        gy_hex: "07192b95ffc8da78631011ed6b24cdd573f977a11e794811",
        frozen: &[
            (
                "p-1",
                &[
                    ("2", 1),
                    ("59", 1),
                    ("149309", 1),
                    ("11393611", 1),
                    ("108341181769254293", 1),
                    ("288626509448065367648032903", 1),
                ],
            ),
            (
                "p+1",
                &[
                    ("2", 64),
                    ("3", 1),
                    ("5", 1),
                    ("17", 1),
                    ("257", 1),
                    ("641", 1),
                    ("65537", 1),
                    ("274177", 1),
                    ("6700417", 1),
                    ("67280421310721", 1),
                ],
            ),
            (
                "n-1",
                &[
                    ("2", 4),
                    ("5", 1),
                    ("2389", 1),
                    ("9564682313913860059195669", 1),
                    ("3433859179316188682119986911", 1),
                ],
            ),
            (
                "n+1",
                &[
                    ("2", 1),
                    ("3", 1),
                    ("41", 1),
                    ("116696426742853568300653", 1),
                    ("218658569360579814643606926178439", 1),
                ],
            ),
            (
                "twist",
                &[
                    ("23", 1),
                    ("10864375060560251605900677743", 1),
                    ("25120401793443689936479125511", 1),
                ],
            ),
            (
                "t",
                &[
                    ("89", 1),
                    ("179", 1),
                    ("383", 1),
                    ("6196490257", 1),
                    ("835990267859", 1),
                ],
            ),
            (
                "4p-t^2",
                &[
                    ("5", 1),
                    ("11", 1),
                    ("31", 1),
                    (
                        "14140398275856956083603613626459809774163796198179979683",
                        1,
                    ),
                ],
            ),
        ],
    },
    Curve {
        name: "P-224",
        f: &[1, 0, 0, -1, 0, 0, 0, 1], // t^7 − t^3 + 1
        f_minus_1: &[
            ("t", T, 3),
            ("Φ_1(t)", PHI1, 1),
            ("Φ_2(t)", PHI2, 1),
            ("Φ_4(t)", PHI4, 1),
        ],
        f_plus_1: &[("t^7 − t^3 + 2 (irreducible)", &[2, 0, 0, -1, 0, 0, 0, 1], 1)],
        p: "26959946667150639794667015087019630673557916260026308143510066298881",
        n: "26959946667150639794667015087019625940457807714424391721682722368061",
        b_hex: "b4050a850c04b3abf54132565044b0b7d7bfd8ba270b39432355ffb4",
        gx_hex: "b70e0cbd6bb4bf7f321390b94a03c1d356c21122343280d6115c1d21",
        gy_hex: "bd376388b5f723fb4c22dfe6cd4375a05a07476444d5819985007e34",
        frozen: &[
            (
                "p-1",
                &[
                    ("2", 96),
                    ("3", 1),
                    ("5", 1),
                    ("17", 1),
                    ("257", 1),
                    ("641", 1),
                    ("65537", 1),
                    ("274177", 1),
                    ("6700417", 1),
                    ("67280421310721", 1),
                ],
            ),
            (
                "p+1",
                &[
                    ("2", 1),
                    ("41", 1),
                    ("2671", 1),
                    ("10637177", 1),
                    ("78987407", 1),
                    ("146503177923784692041805407165164772746256602529", 1),
                ],
            ),
            (
                "n-1",
                &[
                    ("2", 2),
                    ("3", 6),
                    ("5", 1),
                    ("17", 1),
                    ("2153", 1),
                    (
                        "50520606258875818707470860153287666700917696099933389351507",
                        1,
                    ),
                ],
            ),
            (
                "n+1",
                &[
                    ("2", 1),
                    ("7", 1),
                    ("19", 1),
                    ("641707280681", 1),
                    ("157943015469068292359817151928968170360745244072702747", 1),
                ],
            ),
            (
                "twist",
                &[
                    ("3", 2),
                    ("11", 1),
                    ("47", 1),
                    ("3015283", 1),
                    ("40375823", 1),
                    ("267983539294927", 1),
                    ("177594041488131583478651368420021457", 1),
                ],
            ),
            (
                "t",
                &[("17626237", 1), ("55872356921", 1), ("4806058943693873", 1)],
            ),
            (
                "4p-t^2",
                &[
                    ("3", 3),
                    ("29", 1),
                    ("79", 1),
                    ("7523", 1),
                    ("40927", 1),
                    ("11549194661", 1),
                    ("388425074903852603481408727235725774844022299", 1),
                ],
            ),
        ],
    },
    Curve {
        name: "P-256",
        f: &[-1, 0, 0, 1, 0, 0, 1, -1, 1], // t^8 − t^7 + t^6 + t^3 − 1
        f_minus_1: &[
            ("Φ_1(t)", PHI1, 1),
            ("Φ_2(t)", PHI2, 1),
            (
                "t^6 − t^5 + 2t^4 − t^3 + 2t^2 + 2 (irreducible)",
                &[2, 0, 2, -1, 2, -1, 1],
                1,
            ),
        ],
        f_plus_1: &[
            ("t", T, 3),
            ("Φ_4(t)", PHI4, 1),
            ("t^3 − t^2 + 1 (irreducible)", &[1, 0, -1, 1], 1),
        ],
        p: "115792089210356248762697446949407573530086143415290314195533631308867097853951",
        n: "115792089210356248762697446949407573529996955224135760342422259061068512044369",
        b_hex: "5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b",
        gx_hex: "6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296",
        gy_hex: "4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5",
        frozen: &[
            (
                "p-1",
                &[
                    ("2", 1),
                    ("3", 1),
                    ("5", 2),
                    ("17", 1),
                    ("257", 1),
                    ("641", 1),
                    ("1531", 1),
                    ("65537", 1),
                    ("490463", 1),
                    ("6700417", 1),
                    ("835945042244614951780389953367877943453916927241", 1),
                ],
            ),
            (
                "p+1",
                &[
                    ("2", 96),
                    ("7", 1),
                    ("274177", 1),
                    ("67280421310721", 1),
                    ("11318308927973941931404914103", 1),
                ],
            ),
            (
                "n-1",
                &[
                    ("2", 4),
                    ("3", 1),
                    ("71", 1),
                    ("131", 1),
                    ("373", 1),
                    ("3407", 1),
                    ("17449", 1),
                    ("38189", 1),
                    ("187019741", 1),
                    ("622491383", 1),
                    ("1002328039319", 1),
                    ("2624747550333869278416773953", 1),
                ],
            ),
            (
                "n+1",
                &[
                    ("2", 1),
                    ("5", 1),
                    ("1879", 1),
                    ("176337611", 1),
                    (
                        "34946779280882916835155272231706129710560967816144871596921775673",
                        1,
                    ),
                ],
            ),
            (
                "twist",
                &[
                    ("3", 1),
                    ("5", 1),
                    ("13", 1),
                    ("179", 1),
                    (
                        "3317349640749355357762425066592395746459685764401801118712075735758936647",
                        1,
                    ),
                ],
            ),
            (
                "t",
                &[
                    ("113", 1),
                    ("2521", 1),
                    ("12713", 1),
                    ("3450079", 1),
                    ("7138040046852502726273", 1),
                ],
            ),
            (
                "4p-t^2",
                &[
                    ("3", 1),
                    ("5", 1),
                    ("456597257999", 1),
                    ("1428624589419343516204097", 1),
                    ("46523541035814968339936406074986559003387", 1),
                ],
            ),
        ],
    },
];

// ────────────────────────────── small helpers ──────────────────────────────

fn big(s: &str) -> BigInt {
    BigInt::parse_bytes(s.as_bytes(), 10).expect("decimal")
}
fn bighex(s: &str) -> BigInt {
    BigInt::parse_bytes(s.as_bytes(), 16).expect("hex")
}
fn ubig(x: &BigInt) -> BigUint {
    x.to_biguint().expect("non-negative")
}
fn bits(x: &BigInt) -> u64 {
    x.bits()
}
fn ln_abs(x: &BigInt) -> f64 {
    log2_abs(x) * LN2
}
fn r2(x: f64) -> f64 {
    arith::round_to(x, 2)
}

fn poly_mul(a: &[i128], b: &[i128]) -> Vec<i128> {
    let mut out = vec![0i128; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    out
}
fn poly_pow(a: &[i128], e: u32) -> Vec<i128> {
    (0..e).fold(vec![1], |acc, _| poly_mul(&acc, a))
}
fn poly_eval(a: &[i128], t: &BigInt) -> BigInt {
    a.iter()
        .rev()
        .fold(BigInt::zero(), |acc, &c| acc * t + BigInt::from(c))
}
fn poly_str(a: &[i128]) -> String {
    let mut s = String::new();
    for (i, &c) in a.iter().enumerate().rev() {
        if c == 0 {
            continue;
        }
        let sign = if c < 0 {
            "-"
        } else if s.is_empty() {
            ""
        } else {
            "+"
        };
        let mag = c.abs();
        let coef = if mag == 1 && i > 0 {
            String::new()
        } else {
            mag.to_string()
        };
        let var = match i {
            0 => String::new(),
            1 => "t".into(),
            _ => format!("t^{i}"),
        };
        if !s.is_empty() {
            s.push(' ');
        }
        s.push_str(sign);
        if !s.is_empty() && !sign.is_empty() {
            s.push(' ');
        }
        s.push_str(&coef);
        s.push_str(&var);
    }
    if s.is_empty() {
        "0".into()
    } else {
        s
    }
}
fn poly_expr(a: &[i128]) -> Expr {
    let t = || Expr::var("t");
    let terms: Vec<Expr> = a
        .iter()
        .enumerate()
        .filter(|(_, &c)| c != 0)
        .map(|(i, &c)| {
            let mut f = vec![Expr::Const(c as i64)];
            f.extend((0..i).map(|_| t()));
            Expr::product(f)
        })
        .collect();
    Expr::sum(terms)
}

fn frozen_factored(list: &[(&str, u32)]) -> Factored {
    Factored {
        primes: list.iter().map(|(p, e)| (ubig(&big(p)), *e)).collect(),
        composites: vec![],
    }
}
fn factored_value(f: &Factored) -> BigUint {
    f.primes
        .iter()
        .fold(BigUint::one(), |acc, (p, e)| acc * p.pow(*e))
        * f.composites
            .iter()
            .fold(BigUint::one(), |acc, (c, e)| acc * c.pow(*e))
}
fn factored_str(f: &Factored) -> String {
    let mut parts: Vec<String> = f
        .primes
        .iter()
        .map(|(p, e)| {
            if *e > 1 {
                format!("{p}^{e}")
            } else {
                p.to_string()
            }
        })
        .collect();
    parts.extend(f.composites.iter().map(|(c, e)| {
        let s = format!("C{}", c.bits());
        if *e > 1 {
            format!("{s}^{e}")
        } else {
            s
        }
    }));
    parts.join(" · ")
}
fn divisors(f: &Factored) -> Vec<BigInt> {
    let mut ds = vec![BigInt::one()];
    for (p, e) in &f.primes {
        let p = BigInt::from(p.clone());
        let mut next = Vec::with_capacity(ds.len() * (*e as usize + 1));
        for d in &ds {
            let mut q = d.clone();
            for _ in 0..=*e {
                next.push(q.clone());
                q *= &p;
            }
        }
        ds = next;
    }
    ds.sort();
    ds
}

// ────────────────────────────── rationals and E(Q) ──────────────────────────────

#[derive(Clone, Debug, PartialEq, Eq)]
struct Rat {
    n: BigInt,
    d: BigInt,
}
impl Rat {
    fn new(n: BigInt, d: BigInt) -> Rat {
        assert!(!d.is_zero());
        let g = n.gcd(&d);
        let (mut n, mut d) = (n / &g, d / &g);
        if d.is_negative() {
            n = -n;
            d = -d;
        }
        Rat { n, d }
    }
    fn int(x: &BigInt) -> Rat {
        Rat {
            n: x.clone(),
            d: BigInt::one(),
        }
    }
    fn small(x: i64) -> Rat {
        Rat::int(&BigInt::from(x))
    }
    fn add(&self, o: &Rat) -> Rat {
        Rat::new(&self.n * &o.d + &o.n * &self.d, &self.d * &o.d)
    }
    fn sub(&self, o: &Rat) -> Rat {
        Rat::new(&self.n * &o.d - &o.n * &self.d, &self.d * &o.d)
    }
    fn mul(&self, o: &Rat) -> Rat {
        Rat::new(&self.n * &o.n, &self.d * &o.d)
    }
    fn div(&self, o: &Rat) -> Rat {
        Rat::new(&self.n * &o.d, &self.d * &o.n)
    }
    fn is_zero(&self) -> bool {
        self.n.is_zero()
    }
    fn ln_abs(&self) -> f64 {
        ln_abs(&self.n) - ln_abs(&self.d)
    }
}

/// `y² = x³ + A x + B` over `Q`.
struct CurveQ {
    a: Rat,
    b: Rat,
}
type Pt = Option<(Rat, Rat)>;

impl CurveQ {
    fn on_curve(&self, p: &Pt) -> bool {
        match p {
            None => true,
            Some((x, y)) => {
                let lhs = y.mul(y);
                let rhs = x.mul(x).mul(x).add(&self.a.mul(x)).add(&self.b);
                lhs == rhs
            }
        }
    }
    fn add(&self, p: &Pt, q: &Pt) -> Pt {
        let (x1, y1) = match p {
            None => return q.clone(),
            Some(v) => v,
        };
        let (x2, y2) = match q {
            None => return p.clone(),
            Some(v) => v,
        };
        let lambda = if x1 == x2 {
            if y1.add(y2).is_zero() {
                return None;
            }
            let num = Rat::small(3).mul(&x1.mul(x1)).add(&self.a);
            num.div(&Rat::small(2).mul(y1))
        } else {
            y2.sub(y1).div(&x2.sub(x1))
        };
        let x3 = lambda.mul(&lambda).sub(x1).sub(x2);
        let y3 = lambda.mul(&x1.sub(&x3)).sub(y1);
        Some((x3, y3))
    }
    /// Discriminant `−16(4A³ + 27B²)`.
    fn disc(&self) -> Rat {
        let a3 = self.a.mul(&self.a).mul(&self.a);
        let b2 = self.b.mul(&self.b);
        Rat::small(-16).mul(&Rat::small(4).mul(&a3).add(&Rat::small(27).mul(&b2)))
    }
    /// `ĥ(P) ≈ h_x([2^k]P) / 4^k` for `k = 0..=doublings`; returns the
    /// sequence so the reader sees the convergence.  The curve and point
    /// are first moved to an integral model `(x, y) ↦ (u²x, u³y)` (the
    /// canonical height is an isomorphism invariant and the naive-height
    /// limit does not depend on the model), then doubled in Jacobian
    /// coordinates over `Z`, with one gcd per evaluated height.
    fn height_sequence(&self, p: &Pt, doublings: u32) -> Vec<f64> {
        let (x, y) = match p {
            None => return vec![0.0; doublings as usize + 1],
            Some(v) => v,
        };
        let u: BigInt = self.a.d.lcm(&self.b.d).lcm(&x.d).lcm(&y.d);
        let u2 = &u * &u;
        let a_int = self.a.mul(&Rat::int(&(&u2 * &u2)));
        assert!(a_int.d.is_one());
        let a = a_int.n;
        let xs = x.mul(&Rat::int(&u2));
        let ys = y.mul(&Rat::int(&(&u2 * &u)));
        assert!(xs.d.is_one() && ys.d.is_one());
        let (mut xj, mut yj, mut zj) = (xs.n, ys.n, BigInt::one());
        let mut seq = Vec::new();
        for k in 0..=doublings {
            let z2 = &zj * &zj;
            let g = xj.gcd(&z2);
            let h = ln_abs(&(&xj / &g)).max(ln_abs(&(&z2 / &g)));
            seq.push(h / 4f64.powi(k as i32));
            if k < doublings {
                // Jacobian doubling: S = 4XY², M = 3X² + aZ⁴, X' = M² − 2S, Y' = M(S − X') − 8Y⁴, Z' = 2YZ.
                let y2 = &yj * &yj;
                let s = BigInt::from(4) * &xj * &y2;
                let m = BigInt::from(3) * &xj * &xj + &a * &z2 * &z2;
                let x3 = &m * &m - BigInt::from(2) * &s;
                let y3 = &m * (&s - &x3) - BigInt::from(8) * &y2 * &y2;
                let z3 = BigInt::from(2) * &yj * &zj;
                xj = x3;
                yj = y3;
                zj = z3;
            }
        }
        seq
    }
}

/// Affine arithmetic on `y² = x³ − 3x + b` over `F_p` for the lifts.
fn ec_add_fp(
    p: &BigInt,
    a: &BigInt,
    p1: &(BigInt, BigInt),
    p2: &(BigInt, BigInt),
) -> (BigInt, BigInt) {
    let modp = |x: BigInt| x.mod_floor(p);
    let inv = |x: &BigInt| modp(x.clone()).modpow(&(p - 2), p);
    let lambda = if p1.0 == p2.0 {
        modp((BigInt::from(3) * &p1.0 * &p1.0 + a) * inv(&(BigInt::from(2) * &p1.1)))
    } else {
        modp((&p2.1 - &p1.1) * inv(&(&p2.0 - &p1.0)))
    };
    let x3 = modp(&lambda * &lambda - &p1.0 - &p2.0);
    let y3 = modp(&lambda * (&p1.0 - &x3) - &p1.1);
    (x3, y3)
}
fn ec_mul_fp(p: &BigInt, a: &BigInt, k: &BigInt, g: &(BigInt, BigInt)) -> (BigInt, BigInt) {
    let mut r: Option<(BigInt, BigInt)> = None;
    for i in (0..k.bits()).rev() {
        if let Some(cur) = &r {
            r = Some(ec_add_fp(p, a, cur, cur));
        }
        if k.bit(i) {
            r = Some(match &r {
                None => g.clone(),
                Some(cur) => ec_add_fp(p, a, cur, g),
            });
        }
    }
    r.expect("k > 0")
}

// ────────────────────────────── main ──────────────────────────────

struct Opts {
    out: PathBuf,
    scan_x: i64,
    scan_d: i64,
    scan_u: i64,
    doublings: u32,
}

fn parse_opts() -> Opts {
    let mut o = Opts {
        out: PathBuf::from("research/nist_cyclotomic_point_heights_20261011"),
        scan_x: 1_000_000,
        scan_d: 40,
        scan_u: 100_000,
        doublings: 5,
    };
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < args.len() {
        let v = args.get(i + 1).map(|s| s.as_str()).unwrap_or("");
        match args[i].as_str() {
            "--out" => o.out = PathBuf::from(v),
            "--scan-x" => o.scan_x = v.parse().expect("--scan-x"),
            "--scan-d" => o.scan_d = v.parse().expect("--scan-d"),
            "--scan-u" => o.scan_u = v.parse().expect("--scan-u"),
            "--doublings" => o.doublings = v.parse().expect("--doublings"),
            other => panic!("unknown flag {other}"),
        }
        i += 2;
    }
    o
}

fn main() {
    let opts = parse_opts();
    fs::create_dir_all(&opts.out).expect("out dir");
    let t32 = BigInt::from(1u64) << 32;
    let mut md = String::new();
    let mut js = serde_json::Map::new();
    let started = Instant::now();
    let mut svg_series: Vec<(String, Vec<(u32, f64)>, f64)> = Vec::new();

    writeln!(
        md,
        "# NIST P-192 / P-224 / P-256: cyclotomic factorisation and point heights\n"
    )
    .unwrap();
    writeln!(md, "Generated by `cargo run --release --example nist_cyclotomic_point_heights` (native; every number below is recomputed or verified here). Heights are in nats, `ĥ := lim h_x([2^k]P)/4^k` (the `x`-coordinate normalisation, twice Silverman's AEC `ĥ`). Options: scan |x| ≤ {}, denominators d ≤ {}, |u| ≤ {}, {} doublings.\n", opts.scan_x, opts.scan_d, opts.scan_u, opts.doublings).unwrap();

    for c in CURVES {
        let p = big(c.p);
        let n = big(c.n);
        let b = bighex(c.b_hex);
        let gx = bighex(c.gx_hex);
        let gy = bighex(c.gy_hex);
        let t: BigInt = &p + 1 - &n;
        let twist: BigInt = &p + 1 + &t;
        let mut cj = serde_json::Map::new();
        writeln!(md, "\n## {}\n", c.name).unwrap();
        writeln!(
            md,
            "`p = f(2^32)`, `f = {}`; trace `t` has {} bits; `n` is prime: {}.\n",
            poly_str(c.f),
            bits(&t),
            is_prime(&ubig(&n))
        )
        .unwrap();
        assert_eq!(poly_eval(c.f, &t32), p, "{}: f(2^32) != p", c.name);
        assert!(is_prime(&ubig(&p)) && is_prime(&ubig(&n)));
        let on_curve: BigInt =
            (&gy * &gy - (&gx * &gx * &gx - BigInt::from(3) * &gx + &b)).mod_floor(&p);
        assert!(on_curve.is_zero(), "{}: G is not on the curve", c.name);

        // ── 1. algebraic factorisation of f ± 1 ──
        writeln!(md, "### 1. Algebraic factorisation of `f ± 1` in `Z[t]`\n").unwrap();
        writeln!(md, "| value | factorisation | factor sizes at `t = 2^32` | certificate |\n|:--|:--|:--|:--|").unwrap();
        let mut identity_records = Vec::new();
        for (label, target, factors) in
            [("f − 1", -1i128, c.f_minus_1), ("f + 1", 1i128, c.f_plus_1)]
        {
            let mut lhs = c.f.to_vec();
            lhs[0] += target;
            let mut prod = vec![1i128];
            for (_, g, e) in factors {
                prod = poly_mul(&prod, &poly_pow(g, *e));
            }
            while prod.len() > lhs.len() && *prod.last().unwrap() == 0 {
                prod.pop();
            }
            assert_eq!(
                prod, lhs,
                "{}: {} factorisation does not expand back",
                c.name, label
            );
            let fstr: Vec<String> = factors
                .iter()
                .map(|(l, g, e)| {
                    let base = if l.starts_with('Φ') || *l == "t" {
                        l.to_string()
                    } else {
                        format!("({})", poly_str(g))
                    };
                    if *e > 1 {
                        format!("{base}^{e}")
                    } else {
                        base
                    }
                })
                .collect();
            let sizes: Vec<String> = factors
                .iter()
                .map(|(l, g, e)| {
                    let v = poly_eval(g, &t32);
                    let tag = if l.starts_with('Φ') {
                        format!(" = {}", l.replace("(t)", "(2^32)"))
                    } else {
                        String::new()
                    };
                    format!(
                        "{} bits{}{}",
                        bits(&v),
                        if *e > 1 {
                            format!(" ×{e}")
                        } else {
                            String::new()
                        },
                        tag
                    )
                })
                .collect();
            let identity_ref = if factors.len() > 1 {
                let rhs_terms: Vec<Expr> = factors
                    .iter()
                    .flat_map(|(_, g, e)| (0..*e).map(move |_| poly_expr(g)))
                    .collect();
                let st = IdentityStatement {
                    name: format!("{} {}: f(t) {} 1 = {}", c.name, label, if target < 0 { "−" } else { "+" }, fstr.join(" · ")),
                    source: "examples/nist_cyclotomic_point_heights.rs (exact expansion in Z[t], then sealed)".into(),
                    variables: vec!["t".into()],
                    degree_bound: (lhs.len() - 1) as u32,
                    lhs: poly_expr(&lhs),
                    rhs: Expr::product(rhs_terms),
                };
                let record = issue(st, 20261011, 64, "issued and checked by the same tool run; the exact Z[t] expansion is the independent check", None)
                    .expect("identity holds");
                let fname = format!(
                    "identity_{}_{}.json",
                    c.name.to_lowercase().replace('-', ""),
                    if target < 0 { "fm1" } else { "fp1" }
                );
                fs::write(
                    opts.out.join(&fname),
                    serde_json::to_string_pretty(&record).unwrap(),
                )
                .unwrap();
                identity_records.push(json!({"file": fname, "id": record.id}));
                format!("`{}`", record.id)
            } else {
                "— (irreducible)".into()
            };
            writeln!(
                md,
                "| `{label}` | `{}` | {} | {} |",
                fstr.join(" · "),
                sizes.join("; "),
                identity_ref
            )
            .unwrap();
        }
        cj.insert(
            "identity_certificates".into(),
            Value::Array(identity_records),
        );

        // ── 2. frozen integer factorisations, verified ──
        writeln!(
            md,
            "\n### 2. Integer factorisations (frozen inputs, verified here)\n"
        )
        .unwrap();
        writeln!(md, "Discovery: PARI 2.17 `factor`. Verification here: product equals the value, every listed prime passes Miller–Rabin (bases 2…53), and a native trial-division + Brent-rho pass (budget 2^20 rho iterations per composite) is reported so the reader sees which factors the repository's own code reproduces unaided.\n").unwrap();
        writeln!(md, "| value | bits | frozen factorisation (verified) | largest prime (bits) | native pass |\n|:--|--:|:--|--:|:--|").unwrap();
        let mut facs: BTreeMap<&str, Factored> = BTreeMap::new();
        let values: Vec<(&str, BigInt)> = vec![
            ("p-1", &p - 1),
            ("p+1", &p + 1),
            ("n-1", &n - 1),
            ("n+1", &n + 1),
            ("twist", twist.clone()),
            ("t", t.abs()),
            ("4p-t^2", BigInt::from(4) * &p - &t * &t),
        ];
        let mut fj = serde_json::Map::new();
        for (label, value) in &values {
            let list = c
                .frozen
                .iter()
                .find(|(l, _)| l == label)
                .map(|(_, l)| *l)
                .expect("frozen list");
            let f = frozen_factored(list);
            assert_eq!(
                factored_value(&f),
                ubig(value),
                "{}: frozen factorisation of {label} is wrong",
                c.name
            );
            for (q, _) in &f.primes {
                assert!(is_prime(q), "{}: {label} factor {q} is not prime", c.name);
            }
            let largest = f.primes.last().unwrap().0.bits();
            let native = arith::factor(&ubig(value), 1 << 20);
            let native_str = if native.complete() {
                "complete".to_string()
            } else {
                format!(
                    "{} (left {} composite{})",
                    factored_str(&native),
                    native.composites.len(),
                    if native.composites.len() == 1 {
                        ""
                    } else {
                        "s"
                    }
                )
            };
            writeln!(
                md,
                "| `{label}` | {} | `{}` | {} | {} |",
                bits(value),
                factored_str(&f),
                largest,
                native_str
            )
            .unwrap();
            fj.insert(label.to_string(), json!({"value": value.to_string(), "factors": list.iter().map(|(q,e)| json!([q, e])).collect::<Vec<_>>(), "largest_prime_bits": largest, "native_complete": native.complete()}));
            facs.insert(label, f);
        }
        cj.insert("factorisations".into(), Value::Object(fj));

        // ── 3. point-group consequences ──
        writeln!(md, "\n### 3. Consequences for the point group `E(F_p)`\n").unwrap();
        let k = multiplicative_order(&ubig(&p.mod_floor(&n)), &ubig(&n), &facs["n-1"]);
        let cof = ubig(&(&n - 1)) / &k;
        writeln!(md, "- **Embedding degree** `ord_n(p) = (n − 1)/{cof}` ({} bits): the MOV/Frey–Rück transfer lands in `F_{{p^k}}` with `k ≈ 2^{}`; the cyclotomic structure of `p ± 1` is unreachable from the curve.", k.bits(), k.bits() - 1).unwrap();
        writeln!(md, "- **`#E(F_{{p^2}}) = n · n_twist`**, and the twist order `p + 1 + t = N(Φ_2(π))` is `{}` (largest prime {} bits).", factored_str(&facs["twist"]), facs["twist"].primes.last().unwrap().0.bits()).unwrap();
        let s1 = &t;
        assert_eq!(&n * &twist, &p * &p + 1 - lucas_v(s1, &p, 2));
        cj.insert(
            "embedding_degree".into(),
            json!({"k": k.to_string(), "cofactor": cof.to_string(), "bits": k.bits()}),
        );

        // Cheon / Brown–Gallant: d | n−1 needs d+1 powers, d | n+1 needs 2d powers; cost √(n/d) + √d.
        writeln!(md, "\n**Cheon / Brown–Gallant with published powers `[α^i]G`, `i ≤ Q`.** `d | n − 1` needs `Q ≥ d`; `d | n + 1` needs `Q ≥ 2d`. Cost is `√(n/d) + √d` group operations against `√n = 2^{}` generic.\n", r2(log2_abs(&n) / 2.0)).unwrap();
        writeln!(md, "| `log2 Q` | best `d \\| n−1` | cost (bits) | best `d \\| n+1` | cost (bits) | saving vs √n (bits) |\n|--:|--:|--:|--:|--:|--:|").unwrap();
        let d1 = divisors(&facs["n-1"]);
        let d2 = divisors(&facs["n+1"]);
        let smooth1 = log2_abs(&(&n - 1)) - facs["n-1"].primes.last().unwrap().0.bits() as f64;
        let smooth2 = log2_abs(&(&n + 1)) - facs["n+1"].primes.last().unwrap().0.bits() as f64;
        let mut cheon = Vec::new();
        let mut series = Vec::new();
        let cost = |d: &BigInt| log2_abs(&((&n / d).sqrt() + d.sqrt()));
        for lq in [16u32, 24, 32, 40, 48, 56, 64, 80, 96] {
            let q = BigInt::one() << lq;
            let b1 = d1.iter().rev().find(|d| **d < q).unwrap();
            let b2 = d2.iter().rev().find(|d| BigInt::from(2) * *d <= q).unwrap();
            let (c1, c2) = (cost(b1), cost(b2));
            let best = c1.min(c2);
            writeln!(
                md,
                "| {lq} | 2^{} | {} | 2^{} | {} | {} |",
                r2(log2_abs(b1)),
                r2(c1),
                r2(log2_abs(b2)),
                r2(c2),
                r2(log2_abs(&n) / 2.0 - best)
            )
            .unwrap();
            series.push((lq, best));
            cheon.push(json!({"log2_q": lq, "d_minus": b1.to_string(), "cost_minus_bits": r2(c1), "d_plus": b2.to_string(), "cost_plus_bits": r2(c2)}));
        }
        writeln!(md, "\nSmooth part of `n − 1` (everything but its largest prime): 2^{}; of `n + 1`: 2^{}. The saving is capped by the largest divisor the smooth part offers below `Q`.\n", r2(smooth1), r2(smooth2)).unwrap();
        svg_series.push((c.name.to_string(), series, log2_abs(&n) / 2.0));
        cj.insert("cheon".into(), json!({"rows": cheon, "n_minus_1_smooth_bits": r2(smooth1), "n_plus_1_smooth_bits": r2(smooth2)}));

        // Cyclotomic norms N_d = N(Φ_d(π)), π² = tπ − p.
        writeln!(md, "**Cyclotomic factorisation of the extension orders.** `#E(F_{{p^k}}) = ∏_{{d|k}} N_d`, `N_d = N_{{Q(π)/Q}}(Φ_d(π))`, verified against `p^k + 1 − V_k(t, p)` for `k ≤ 12`. Each `N_d` is trial-divided and rho'd with budget 2^16; `C` marks an unsplit composite, `P` a probable prime.\n").unwrap();
        writeln!(md, "| `d` | bits | factorisation of `N_d` |\n|--:|--:|:--|").unwrap();
        let norm = |a: &BigInt, bb: &BigInt| a * a + a * bb * &t + bb * bb * &p; // N(a + bπ)
        let mut norms = BTreeMap::new();
        for d in 1u32..=12 {
            let phi = cyclotomic(d);
            // Horner in Z[π]: value = a + bπ.
            let (mut a, mut bb) = (BigInt::zero(), BigInt::zero());
            for &coef in phi.iter().rev() {
                // (a + bπ)·π = aπ + b(tπ − p) = −bp + (a + bt)π
                let na = -&bb * &p + BigInt::from(coef);
                let nb = &a + &bb * &t;
                a = na;
                bb = nb;
            }
            let nd = norm(&a, &bb);
            assert!(nd.is_positive());
            let f = arith::factor(&ubig(&nd), 1 << 16);
            let parts: Vec<String> = f
                .primes
                .iter()
                .map(|(q, e)| {
                    let s = if q.bits() > 40 {
                        format!("P{}", q.bits())
                    } else {
                        q.to_string()
                    };
                    if *e > 1 {
                        format!("{s}^{e}")
                    } else {
                        s
                    }
                })
                .chain(f.composites.iter().map(|(q, e)| {
                    if *e > 1 {
                        format!("C{}^{e}", q.bits())
                    } else {
                        format!("C{}", q.bits())
                    }
                }))
                .collect();
            writeln!(md, "| {d} | {} | `{}` |", bits(&nd), parts.join(" · ")).unwrap();
            norms.insert(d, (nd, parts));
        }
        for k in 1u32..=12 {
            let prod: BigInt = norms
                .iter()
                .filter(|(d, _)| k % **d == 0)
                .map(|(_, (v, _))| v.clone())
                .product();
            let order = p.pow(k) + 1 - lucas_v(s1, &p, k);
            assert_eq!(prod, order, "{}: cyclotomic product != #E(F_p^{k})", c.name);
        }
        assert_eq!(norms[&1].0, n);
        assert_eq!(norms[&2].0, twist);
        cj.insert(
            "cyclotomic_norms".into(),
            Value::Object(
                norms
                    .iter()
                    .map(|(d, (v, parts))| {
                        (d.to_string(), json!({"bits": bits(v), "factors": parts}))
                    })
                    .collect(),
            ),
        );

        // ── 4. heights ──
        writeln!(md, "\n### 4. Point heights on lifts to `Q`\n").unwrap();
        let eb = CurveQ {
            a: Rat::small(-3),
            b: Rat::int(&b),
        };
        let disc_b = eb.disc();
        writeln!(md, "**Canonical lift** `E_b/Q: y² = x³ − 3x + b` with the integer `b ∈ [0, p)`: `log2 |Δ| = {}` (`Δ = −432 (b − 2)(b + 2)`).\n", r2(log2_abs(&disc_b.n))).unwrap();
        // integer-point census with residue filter (squares mod 64, 63, 65, 11)
        let sq_table = |m: u64| -> Vec<bool> {
            let mut v = vec![false; m as usize];
            for i in 0..m {
                v[((i * i) % m) as usize] = true;
            }
            v
        };
        let mods = [64u64, 63, 65, 11];
        let tables: Vec<Vec<bool>> = mods.iter().map(|&m| sq_table(m)).collect();
        let bmods: Vec<i128> = mods.iter().map(|&m| (&b % m).to_i128().unwrap()).collect();
        let passes = |u: i128, dpow4: i128, dpow6: i128| -> bool {
            mods.iter()
                .zip(&tables)
                .zip(&bmods)
                .all(|((&m, tab), &bm)| {
                    let m = m as i128;
                    let v = ((u % m) * (u % m) % m * (u % m) % m
                        - 3 * (u % m) % m * (dpow4 % m) % m
                        + bm * (dpow6 % m) % m)
                        .rem_euclid(m);
                    tab[v as usize]
                })
        };
        let t0 = Instant::now();
        let mut int_pts = Vec::new();
        let mut candidates = 0u64;
        for x in -opts.scan_x..=opts.scan_x {
            if !passes(x as i128, 1, 1) {
                continue;
            }
            candidates += 1;
            let xb = BigInt::from(x);
            let v = &xb * &xb * &xb - BigInt::from(3) * &xb + &b;
            if v.is_negative() {
                continue;
            }
            let s = v.sqrt();
            if &s * &s == v {
                int_pts.push((x, s.to_string()));
            }
        }
        let t_int = t0.elapsed().as_secs_f64();
        let t0 = Instant::now();
        let mut rat_pts = Vec::new();
        let mut rat_candidates = 0u64;
        for d in 2..=opts.scan_d {
            let (d4, d6) = ((d as i128).pow(4), (d as i128).pow(6));
            let (bd4, bd6) = (BigInt::from(d4), BigInt::from(d6));
            for u in -opts.scan_u..=opts.scan_u {
                if u.gcd(&d) != 1 || !passes(u as i128, d4, d6) {
                    continue;
                }
                rat_candidates += 1;
                let ub = BigInt::from(u);
                let w = &ub * &ub * &ub - BigInt::from(3) * &ub * &bd4 + &b * &bd6;
                if w.is_negative() {
                    continue;
                }
                let s = w.sqrt();
                if &s * &s == w {
                    rat_pts.push((u, d, s.to_string()));
                }
            }
        }
        let t_rat = t0.elapsed().as_secs_f64();
        let trials =
            (2 * opts.scan_x + 1) as f64 + (opts.scan_d - 1) as f64 * (2 * opts.scan_u + 1) as f64;
        let expected = trials * 2f64.powf(-log2_abs(&b) / 2.0);
        writeln!(md, "| census | trials | passed residue filter | points found | heuristic expectation |\n|:--|--:|--:|--:|--:|").unwrap();
        writeln!(
            md,
            "| integer `x`, `\\|x\\| ≤ {}` | {} | {} | {} | — |",
            opts.scan_x,
            2 * opts.scan_x + 1,
            candidates,
            int_pts.len()
        )
        .unwrap();
        writeln!(
            md,
            "| rational `x = u/d²`, `2 ≤ d ≤ {}`, `\\|u\\| ≤ {}` | {} | {} | {} | — |",
            opts.scan_d,
            opts.scan_u,
            (opts.scan_d - 1) * (2 * opts.scan_u + 1),
            rat_candidates,
            rat_pts.len()
        )
        .unwrap();
        writeln!(md, "| both | {:.3e} | | {} | {:.1e} (each trial is a square with probability ≈ `b^{{−1/2}}`) |", trials, int_pts.len() + rat_pts.len(), expected).unwrap();
        writeln!(
            md,
            "\nScan time {:.1} s + {:.1} s. Points found: `{:?}` / `{:?}`.\n",
            t_int, t_rat, int_pts, rat_pts
        )
        .unwrap();

        // xedni-style lift of G: curve through the integer point (gx, gy) with a = −3
        let b1 = &gy * &gy - &gx * &gx * &gx + BigInt::from(3) * &gx;
        let e1 = CurveQ {
            a: Rat::small(-3),
            b: Rat::int(&b1),
        };
        let g1: Pt = Some((Rat::int(&gx), Rat::int(&gy)));
        assert!(e1.on_curve(&g1));
        let seq = e1.height_sequence(&g1, opts.doublings);
        let hhat = *seq.last().unwrap();
        let conv = (seq[seq.len() - 1] - seq[seq.len() - 2]).abs();
        let ld1 = e1.disc().ln_abs();
        writeln!(md, "**One-point lift (xedni, `r = 1`).** `E_1: y² = x³ − 3x + B_1` through `G̃ = (G_x, G_y) ∈ Z²`, `B_1 = G_y² − G_x³ + 3G_x` ({} bits), `log |Δ(E_1)| = {:.1}` nats.\n", bits(&b1), ld1).unwrap();
        writeln!(md, "| `k` | `h_x([2^k]G̃) / 4^k` |\n|--:|--:|").unwrap();
        for (k, h) in seq.iter().enumerate() {
            writeln!(md, "| {k} | {:.4} |", h).unwrap();
        }
        writeln!(md, "\n`ĥ(G̃) = {:.3}` nats = {:.1} bits (last step moved {:.1e}); `ĥ / log p = {:.4}`; `ĥ / log|Δ(E_1)| = {:.4}`; `h_x(G̃) − ĥ(G̃) = {:.3}` nats.\n", hhat, hhat / LN2, conv, hhat / ln_abs(&p), hhat / ld1, ln_abs(&gx) - hhat).unwrap();

        // two-point lift through G̃ and R̃ = lift of [k]G
        let kscal = big("1311768467294899695"); // 0x1234567890abcdef
        let r = ec_mul_fp(&p, &BigInt::from(-3), &kscal, &(gx.clone(), gy.clone()));
        let (x1, y1, x2, y2) = (Rat::int(&gx), Rat::int(&gy), Rat::int(&r.0), Rat::int(&r.1));
        let c1 = y1.mul(&y1).sub(&x1.mul(&x1).mul(&x1));
        let c2 = y2.mul(&y2).sub(&x2.mul(&x2).mul(&x2));
        let a2 = c1.sub(&c2).div(&x1.sub(&x2));
        let b2 = c1.sub(&a2.mul(&x1));
        let e2 = CurveQ { a: a2, b: b2 };
        let pg: Pt = Some((x1, y1));
        let pr: Pt = Some((x2, y2));
        assert!(e2.on_curve(&pg) && e2.on_curve(&pr));
        let hg = *e2.height_sequence(&pg, opts.doublings).last().unwrap();
        let hr = *e2.height_sequence(&pr, opts.doublings).last().unwrap();
        let hgr = *e2
            .height_sequence(&e2.add(&pg, &pr), opts.doublings)
            .last()
            .unwrap();
        let pairing = (hgr - hg - hr) / 2.0;
        let reg = hg * hr - pairing * pairing;
        let cosang = pairing / (hg * hr).sqrt();
        let ld2 = e2.disc().ln_abs();
        writeln!(md, "**Two-point lift (`r = 2`).** `E_2: y² = x³ + A x + B` with `A, B ∈ Q` through `G̃` and `R̃`, `R = [k]G`, `k = 0x1234567890abcdef`; `log |Δ(E_2)| = {:.1}` nats (`A` has a {}-bit denominator).\n", ld2, bits(&e2.a.d)).unwrap();
        writeln!(md, "| `ĥ(G̃)` | `ĥ(R̃)` | `ĥ(G̃ + R̃)` | `⟨G̃, R̃⟩` | cos ∠ | regulator `det` | `ĥ / log p` |\n|--:|--:|--:|--:|--:|--:|--:|").unwrap();
        writeln!(
            md,
            "| {:.2} | {:.2} | {:.2} | {:.3} | {:.5} | {:.4e} | {:.3}, {:.3} |\n",
            hg,
            hr,
            hgr,
            pairing,
            cosang,
            reg,
            hg / ln_abs(&p),
            hr / ln_abs(&p)
        )
        .unwrap();
        writeln!(md, "Reading: the lifted points are independent (regulator ≈ `ĥ(G̃)·ĥ(R̃)`, pairing ≈ 0). On `E_1` the generator's height is `log p` to within one nat; on `E_2` it is larger because clearing the denominator of `A` scales the model. Lifting cannot make a uniformly random point of `E(F_p)` smaller than about `log p` nats: a point admits a lift of naive height ≤ `H` only if `|x| ≤ e^H` or `|x − p| ≤ e^H`, a fraction `2e^H/p` of the group, which is exactly the small-`x` factor base and nothing more.\n").unwrap();
        cj.insert("heights".into(), json!({
            "canonical_lift_log2_disc": r2(log2_abs(&disc_b.n)),
            "integer_points": int_pts, "rational_points": rat_pts, "trials": trials, "expected_points": expected,
            "r1": {"log2_B1": bits(&b1), "hhat_nats": hhat, "sequence": seq, "ln_disc": ld1, "hhat_over_ln_p": hhat / ln_abs(&p), "hhat_over_ln_disc": hhat / ld1},
            "r2": {"h_g": hg, "h_r": hr, "h_gr": hgr, "pairing": pairing, "cos": cosang, "regulator": reg, "ln_disc": ld2}
        }));
        js.insert(c.name.into(), Value::Object(cj));
    }

    fs::write(opts.out.join("cheon_cost.svg"), cheon_svg(&svg_series)).unwrap();
    writeln!(
        md,
        "\n---\nTotal run time {:.1} s.",
        started.elapsed().as_secs_f64()
    )
    .unwrap();
    fs::write(opts.out.join("RESULTS.md"), &md).unwrap();
    fs::write(
        opts.out.join("results.json"),
        serde_json::to_string_pretty(&Value::Object(js)).unwrap(),
    )
    .unwrap();
    eprintln!("wrote {}/RESULTS.md and results.json", opts.out.display());
}

/// Coefficients of `Φ_d(x)`, low degree first, for `d ≤ 12` via
/// `x^d − 1 = ∏_{e | d} Φ_e(x)` and exact polynomial division.
#[allow(clippy::manual_is_multiple_of)] // keeps the 1.87 toolchain in CI
fn cyclotomic(d: u32) -> Vec<i128> {
    let mut num = vec![0i128; d as usize + 1];
    num[0] = -1;
    num[d as usize] = 1;
    for e in 1..d {
        if d % e == 0 {
            num = poly_div_exact(&num, &cyclotomic(e));
        }
    }
    num
}
fn poly_div_exact(a: &[i128], b: &[i128]) -> Vec<i128> {
    let mut rem = a.to_vec();
    let db = b.len() - 1;
    let mut q = vec![0i128; a.len() - db];
    for i in (0..q.len()).rev() {
        let coef = rem[i + db] / b[db];
        q[i] = coef;
        for (j, &bc) in b.iter().enumerate() {
            rem[i + j] -= coef * bc;
        }
    }
    assert!(rem.iter().all(|&c| c == 0), "inexact division");
    q
}

/// Line chart: best Cheon / Brown–Gallant cost (bits) against `log2 Q`
/// published powers, one series per curve, dashed generic `√n` lines.
fn cheon_svg(series: &[(String, Vec<(u32, f64)>, f64)]) -> String {
    let (w, h, ml, mr, mt, mb) = (720.0, 440.0, 64.0, 24.0, 40.0, 76.0);
    let (x0, x1, y0, y1) = (16.0, 96.0, 40.0, 132.0);
    let px = |x: f64| ml + (x - x0) / (x1 - x0) * (w - ml - mr);
    let py = |y: f64| mt + (y1 - y) / (y1 - y0) * (h - mt - mb);
    let colours = ["#1b6ca8", "#c0392b", "#2e8b57"];
    let mut s = String::new();
    writeln!(s, r##"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}" font-family="Helvetica, Arial, sans-serif" font-size="12">"##).unwrap();
    writeln!(s, r##"<rect width="100%" height="100%" fill="white"/>"##).unwrap();
    writeln!(s, r##"<text x="{}" y="22" text-anchor="middle" font-size="14" font-weight="bold">Cheon / Brown–Gallant cost from the factorisation of n ∓ 1 (best of both sides)</text>"##, w / 2.0).unwrap();
    for gx in (x0 as u32..=x1 as u32).step_by(16) {
        writeln!(s, r##"<line x1="{0}" y1="{1}" x2="{0}" y2="{2}" stroke="#ddd"/><text x="{0}" y="{3}" text-anchor="middle">{4}</text>"##, px(gx as f64), py(y0), py(y1), py(y0) + 16.0, gx).unwrap();
    }
    for gy in (y0 as u32..=y1 as u32).step_by(16) {
        writeln!(s, r##"<line x1="{0}" y1="{1}" x2="{2}" y2="{1}" stroke="#ddd"/><text x="{3}" y="{4}" text-anchor="end">{5}</text>"##, px(x0), py(gy as f64), px(x1), px(x0) - 6.0, py(gy as f64) + 4.0, gy).unwrap();
    }
    writeln!(s, r##"<text x="{}" y="{}" text-anchor="middle">log2 Q — number of published powers [α^i]G</text>"##, (px(x0) + px(x1)) / 2.0, h - 34.0).unwrap();
    writeln!(s, r##"<text transform="translate(16,{}) rotate(-90)" text-anchor="middle">cost, log2 group operations</text>"##, (py(y0) + py(y1)) / 2.0).unwrap();
    for (i, (name, pts, generic)) in series.iter().enumerate() {
        let col = colours[i % colours.len()];
        writeln!(s, r##"<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{col}" stroke-dasharray="6 4" stroke-width="1.2"/>"##, px(x0), py(*generic), px(x1), py(*generic)).unwrap();
        let path: Vec<String> = pts
            .iter()
            .enumerate()
            .map(|(j, (lq, c))| {
                format!(
                    "{}{:.1},{:.1}",
                    if j == 0 { "M" } else { "L" },
                    px(*lq as f64),
                    py(*c)
                )
            })
            .collect();
        writeln!(
            s,
            r##"<path d="{}" fill="none" stroke="{col}" stroke-width="2.4"/>"##,
            path.join(" ")
        )
        .unwrap();
        for (lq, c) in pts {
            writeln!(
                s,
                r##"<circle cx="{:.1}" cy="{:.1}" r="3.2" fill="{col}"/>"##,
                px(*lq as f64),
                py(*c)
            )
            .unwrap();
        }
        let ly = py(y0) - 62.0 + 18.0 * i as f64;
        writeln!(s, r##"<line x1="{}" y1="{ly}" x2="{}" y2="{ly}" stroke="{col}" stroke-width="2.4"/><text x="{}" y="{}">{name} (dashed: generic √n = 2^{:.0})</text>"##, px(x0) + 12.0, px(x0) + 40.0, px(x0) + 48.0, ly + 4.0, generic).unwrap();
    }
    writeln!(s, r##"<text x="{}" y="{}" font-size="10" fill="#444">Source: results.json here. Exact divisors of n − 1, n + 1; cost √(n/d) + √d, d ≤ Q (resp. 2d ≤ Q). Three curves, exact values, no uncertainty.</text>"##, ml, h - 12.0).unwrap();
    s.push_str("</svg>\n");
    s
}

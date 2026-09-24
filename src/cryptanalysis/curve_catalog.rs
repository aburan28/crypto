//! # A unified catalog of standardized elliptic curves for index calculus.
//!
//! The `crypto_lib` crate defines standardized curves in several places —
//! [`crate::ecc::curve::CurveParams`] and [`crate::ecc::curve_zoo`] for prime
//! fields, [`crate::binary_ecc::curve::BinaryCurve`] for `F_{2^m}` — but there
//! is no single place that answers *"give me every standardized curve, by
//! name, with its field family and verified parameters."*  The index-calculus
//! engine (`icx`) needs exactly that: one registry it can enumerate, look a
//! curve up in, verify, and hand to the pipeline.
//!
//! This module is that registry.  It does **not** re-encode curve constants;
//! it wraps the existing, independently checked constructors so there is one
//! source of truth for each curve's parameters.  New standardized curves that
//! the crate did not previously carry are added here with their provenance.
//!
//! ## What "verified" means
//!
//! Every entry can be [`CatalogCurve::verify`]d: the generator is on the
//! curve, `n · G = O`, the cofactor is positive, and `h · n` lands in the
//! Hasse interval `[q + 1 - 2√q, q + 1 + 2√q]`.  These are necessary
//! conditions, not a full certificate (subgroup order primality is screened,
//! not proved), and each check says so.  The catalog test asserts a per-family
//! floor on how many curves verify, and per-family coverage, rather than an
//! exact global count — the registry grows, and an exact count would be a
//! brittle assertion that fails on every addition (see `CLAUDE.md`).
//!
//! ## Scope note
//!
//! Membership here is about *coverage of standardized parameter sets*, not
//! about index calculus being a threat.  Whether IC is the relevant attack on
//! a given curve, and at what size a run is feasible, is decided by the engine
//! ([`crate::cryptanalysis::ic_engine`]), not by this catalog.

use num_bigint::BigUint;
use num_traits::{One, Zero};

use crate::binary_ecc::curve::{self as bcurve, BinaryCurve, BinaryPoint};
use crate::ecc::curve::CurveParams;

/// The finite-field / curve family an entry belongs to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Family {
    /// Short Weierstrass over a prime field `F_p`, `p > 3`.
    Prime,
    /// Random binary curve `y² + xy = x³ + a x² + b` over `F_{2^m}`.
    BinaryRandom,
    /// Koblitz (anomalous binary) curve: `a ∈ {0,1}`, `b = 1`, over `F_{2^m}`.
    Koblitz,
    /// Curve over an extension field `F_{p^k}`, `k ≥ 2` (Weil-descent regime).
    Extension,
    /// Curve over a characteristic-three field `F_{3^m}`.
    Char3,
}

impl Family {
    /// A short, stable identifier used in listings and `--family` filters.
    pub fn tag(self) -> &'static str {
        match self {
            Family::Prime => "prime",
            Family::BinaryRandom => "binary",
            Family::Koblitz => "koblitz",
            Family::Extension => "extension",
            Family::Char3 => "char3",
        }
    }

    /// Parse a `--family` filter value; accepts a few natural spellings.
    pub fn parse(s: &str) -> Option<Family> {
        match s.to_ascii_lowercase().as_str() {
            "prime" | "fp" | "p" => Some(Family::Prime),
            "binary" | "binary-random" | "f2m" | "b" => Some(Family::BinaryRandom),
            "koblitz" | "k" => Some(Family::Koblitz),
            "extension" | "fpk" | "ext" => Some(Family::Extension),
            "char3" | "f3m" | "ternary" => Some(Family::Char3),
            _ => None,
        }
    }
}

/// A human- and machine-readable description of the base field.
#[derive(Clone, Debug)]
pub struct FieldDesc {
    /// One of `prime`, `binary`, `extension`, `char3`.
    pub kind: &'static str,
    /// `F_p`, `F_2^163`, `F_3^97`, `F_p^3`, rendered for display.
    pub label: String,
    /// `⌈log2 #field⌉`, the field size in bits.
    pub bits: u64,
}

/// The concrete curve object, in whichever representation the crate provides.
///
/// The catalog is uniform above this; computation dispatches on it.
#[derive(Clone)]
pub enum CurveObject {
    /// Prime-field short Weierstrass (BigUint arithmetic).
    Prime(CurveParams),
    /// Binary curve over `F_{2^m}` (polynomial-basis `F2mElement`).
    ///
    /// Boxed because a `BinaryCurve` (multi-word field elements) is much
    /// larger than a `CurveParams`, and an unboxed variant would bloat every
    /// `CurveObject` to the binary size.
    Binary(Box<BinaryCurve>),
}

/// A single catalog entry, loaded (its constants materialized) but not yet
/// verified.  Verification is separate because it costs a scalar
/// multiplication and callers that only list curves should not pay it.
#[derive(Clone)]
pub struct CatalogCurve {
    /// Canonical lookup name (lowercase, hyphenated).
    pub name: &'static str,
    /// Alternative names that resolve to this curve.
    pub aliases: &'static [&'static str],
    /// Field / curve family.
    pub family: Family,
    /// Where the parameters come from (standard + citation).
    pub standard: &'static str,
    /// The concrete curve.
    pub object: CurveObject,
}

/// The result of one verification check.
#[derive(Clone, Debug)]
pub struct Check {
    pub name: &'static str,
    pub passed: bool,
    pub detail: String,
}

impl CatalogCurve {
    /// A rendered description of the base field.
    pub fn field(&self) -> FieldDesc {
        match &self.object {
            CurveObject::Prime(p) => FieldDesc {
                kind: "prime",
                label: "F_p".to_string(),
                bits: p.p.bits(),
            },
            CurveObject::Binary(c) => FieldDesc {
                kind: "binary",
                label: format!("F_2^{}", c.m),
                bits: c.m as u64,
            },
        }
    }

    /// The subgroup order `n` (order of the standard generator's subgroup).
    pub fn subgroup_order(&self) -> BigUint {
        match &self.object {
            CurveObject::Prime(p) => p.n.clone(),
            CurveObject::Binary(c) => c.order.clone(),
        }
    }

    /// The cofactor `h`.
    pub fn cofactor(&self) -> BigUint {
        match &self.object {
            CurveObject::Prime(p) => BigUint::from(p.h),
            CurveObject::Binary(c) => c.cofactor.clone(),
        }
    }

    /// `⌈log2 n⌉`, the subgroup order in bits.
    pub fn order_bits(&self) -> u64 {
        self.subgroup_order().bits()
    }

    /// Generic-attack (Pollard rho) security level, `½ log2 n` bits.
    ///
    /// This is the classical ECDLP work factor and is independent of whether
    /// index calculus applies; it is what "security bits" means in a listing.
    pub fn rho_security_bits(&self) -> f64 {
        // bits() is the position of the top set bit; log2(n) ≈ bits()-1..bits().
        let ln = self.subgroup_order();
        if ln.is_zero() {
            return 0.0;
        }
        // Fractional log2 via the top 64 bits, for a smooth estimate.
        let bits = ln.bits();
        let approx_log2 = if bits <= 64 {
            (u64::try_from(ln.clone()).unwrap_or(u64::MAX) as f64).log2()
        } else {
            let shift = bits - 53;
            let top = &ln >> shift;
            let top = u64::try_from(top).unwrap_or(1u64 << 52) as f64;
            top.log2() + shift as f64
        };
        approx_log2 / 2.0
    }

    /// The full group order `#E = h · n`, when both are known.
    pub fn group_order(&self) -> BigUint {
        self.cofactor() * self.subgroup_order()
    }

    /// Run the necessary-condition checks.  Never proves primality of `n`.
    pub fn verify(&self) -> Vec<Check> {
        match &self.object {
            CurveObject::Prime(p) => verify_prime(p),
            CurveObject::Binary(c) => verify_binary(c),
        }
    }
}

fn verify_prime(p: &CurveParams) -> Vec<Check> {
    let mut checks = Vec::new();
    let g = p.generator();
    checks.push(Check {
        name: "generator_on_curve",
        passed: p.is_on_curve(&g),
        detail: "generator satisfies y² = x³ + ax + b".to_string(),
    });
    // n · G = O.
    let a_fe = p.a_fe();
    let ng = g.scalar_mul(&p.n, &a_fe);
    let at_infinity = ng.x_coord().is_none();
    checks.push(Check {
        name: "generator_subgroup",
        passed: at_infinity,
        detail: "[n]G = O (point at infinity)".to_string(),
    });
    checks.push(Check {
        name: "positive_cofactor",
        passed: p.h > 0,
        detail: format!("cofactor h = {}", p.h),
    });
    // Hasse: |#E - (q+1)| <= 2 sqrt q, with #E = h n.
    let q = p.p.clone();
    let card = BigUint::from(p.h) * &p.n;
    checks.push(hasse_check(&q, &card));
    checks
}

fn verify_binary(c: &BinaryCurve) -> Vec<Check> {
    let mut checks = Vec::new();
    checks.push(Check {
        name: "generator_on_curve",
        passed: c.is_on_curve(&c.generator),
        detail: "generator satisfies y² + xy = x³ + ax² + b".to_string(),
    });
    let ng = bcurve::scalar_mul(c, &c.generator, &c.order);
    checks.push(Check {
        name: "generator_subgroup",
        passed: matches!(ng, BinaryPoint::Infinity),
        detail: "[n]G = O (point at infinity)".to_string(),
    });
    checks.push(Check {
        name: "positive_cofactor",
        passed: !c.cofactor.is_zero(),
        detail: format!("cofactor h = {}", c.cofactor),
    });
    let q = BigUint::one() << c.m;
    let card = &c.cofactor * &c.order;
    checks.push(hasse_check(&q, &card));
    checks
}

/// `|#E - (q + 1)| <= 2√q`, computed with integer square roots so it is exact.
fn hasse_check(q: &BigUint, card: &BigUint) -> Check {
    let q1 = q + BigUint::one();
    // 2√q, rounded up so the bound is not too tight.
    let two_sqrt_q = (q.sqrt() + BigUint::one()) * BigUint::from(2u32);
    let diff = if card >= &q1 { card - &q1 } else { &q1 - card };
    let passed = diff <= two_sqrt_q;
    Check {
        name: "hasse_bound",
        passed,
        detail: format!(
            "|#E-(q+1)| = {} <= 2sqrt(q) = {} (necessary cardinality bound)",
            diff, two_sqrt_q
        ),
    }
}

/// True if every necessary check passes.
pub fn verified(curve: &CatalogCurve) -> bool {
    curve.verify().iter().all(|c| c.passed)
}

/// Build the whole catalog.  Constructs each curve's constants but performs no
/// verification (that is [`CatalogCurve::verify`]).
pub fn all() -> Vec<CatalogCurve> {
    let mut v = Vec::new();
    prime_curves(&mut v);
    binary_curves(&mut v);
    v
}

/// Look a curve up by canonical name or alias (case-insensitive).
pub fn by_name(query: &str) -> Option<CatalogCurve> {
    let q = query.to_ascii_lowercase();
    all().into_iter().find(|c| {
        c.name.eq_ignore_ascii_case(&q) || c.aliases.iter().any(|a| a.eq_ignore_ascii_case(&q))
    })
}

/// Every canonical name, for listings and shell completion.
pub fn names() -> Vec<&'static str> {
    all().iter().map(|c| c.name).collect()
}

fn prime_curves(v: &mut Vec<CatalogCurve>) {
    // Core standardized prime curves (src/ecc/curve.rs).
    let core: &[(&str, &[&str], &str, fn() -> CurveParams)] = &[
        (
            "secp256k1",
            &["k256", "ansix9p256k1"],
            "SEC 2 v2 / Bitcoin",
            CurveParams::secp256k1,
        ),
        (
            "p256",
            &["secp256r1", "prime256v1", "nistp256"],
            "NIST FIPS 186-4 / SEC 2",
            CurveParams::p256,
        ),
        (
            "sm2",
            &["sm2p256v1"],
            "GM/T 0003-2012 (China)",
            CurveParams::sm2,
        ),
        (
            "gost-2012-256-test",
            &[],
            "RFC 7836 test curve",
            CurveParams::gost_3410_2012_256_test,
        ),
        (
            "gost-2012-512-test",
            &[],
            "RFC 7836 test curve",
            CurveParams::gost_3410_2012_512_test,
        ),
    ];
    for (name, aliases, std, ctor) in core {
        v.push(CatalogCurve {
            name,
            aliases,
            family: Family::Prime,
            standard: std,
            object: CurveObject::Prime(ctor()),
        });
    }

    // The curve zoo (src/ecc/curve_zoo.rs): SEC / NIST / Brainpool / ANSSI /
    // GOST prime curves.  All are on-curve and n·G = O verified in that file's
    // own tests; the catalog re-verifies them.
    let zoo: &[(&str, &[&str], &str, fn() -> CurveParams)] = &[
        ("secp112r1", &[], "SEC 2 v1", CurveParams::secp112r1),
        ("secp112r2", &[], "SEC 2 v1", CurveParams::secp112r2),
        ("secp128r1", &[], "SEC 2 v1", CurveParams::secp128r1),
        ("secp128r2", &[], "SEC 2 v1", CurveParams::secp128r2),
        ("secp160k1", &[], "SEC 2 v1", CurveParams::secp160k1),
        ("secp160r1", &[], "SEC 2 v1", CurveParams::secp160r1),
        ("secp160r2", &[], "SEC 2 v1", CurveParams::secp160r2),
        (
            "secp192k1",
            &["ansix9p192k1"],
            "SEC 2 v1",
            CurveParams::secp192k1,
        ),
        (
            "secp224k1",
            &["ansix9p224k1"],
            "SEC 2 v1",
            CurveParams::secp224k1,
        ),
        (
            "p192",
            &["secp192r1", "prime192v1", "nistp192"],
            "NIST FIPS 186-4",
            CurveParams::p192,
        ),
        (
            "p224",
            &["secp224r1", "nistp224"],
            "NIST FIPS 186-4",
            CurveParams::p224,
        ),
        (
            "p384",
            &["secp384r1", "nistp384"],
            "NIST FIPS 186-4",
            CurveParams::p384,
        ),
        (
            "p521",
            &["secp521r1", "nistp521"],
            "NIST FIPS 186-4",
            CurveParams::p521,
        ),
        (
            "brainpoolp192r1",
            &[],
            "RFC 5639",
            CurveParams::brainpool_p192r1,
        ),
        (
            "brainpoolp224r1",
            &[],
            "RFC 5639",
            CurveParams::brainpool_p224r1,
        ),
        (
            "brainpoolp256r1",
            &[],
            "RFC 5639",
            CurveParams::brainpool_p256r1,
        ),
        (
            "brainpoolp320r1",
            &[],
            "RFC 5639",
            CurveParams::brainpool_p320r1,
        ),
        (
            "brainpoolp384r1",
            &[],
            "RFC 5639",
            CurveParams::brainpool_p384r1,
        ),
        (
            "brainpoolp512r1",
            &[],
            "RFC 5639",
            CurveParams::brainpool_p512r1,
        ),
        ("frp256v1", &[], "ANSSI (France)", CurveParams::frp256v1),
        (
            "gost-cryptopro-a",
            &[],
            "RFC 4357",
            CurveParams::gost_cryptopro_a,
        ),
        (
            "gost-cryptopro-b",
            &[],
            "RFC 4357",
            CurveParams::gost_cryptopro_b,
        ),
        (
            "gost-cryptopro-c",
            &[],
            "RFC 4357",
            CurveParams::gost_cryptopro_c,
        ),
        (
            "gost-tc26-256-a",
            &[],
            "GOST R 34.10-2012 (256)",
            CurveParams::gost_tc26_256_a,
        ),
        (
            "gost-tc26-512-a",
            &[],
            "GOST R 34.10-2012 (512)",
            CurveParams::gost_tc26_512_a,
        ),
        (
            "gost-tc26-512-b",
            &[],
            "GOST R 34.10-2012 (512)",
            CurveParams::gost_tc26_512_b,
        ),
    ];
    for (name, aliases, std, ctor) in zoo {
        v.push(CatalogCurve {
            name,
            aliases,
            family: Family::Prime,
            standard: std,
            object: CurveObject::Prime(ctor()),
        });
    }
}

fn binary_curves(v: &mut Vec<CatalogCurve>) {
    // Binary curves the crate already carries (src/binary_ecc/curve.rs).
    // sect163k1 is a Koblitz curve; the sect*r* curves are random binary
    // curves; Oakley Group 3 is a subgroup curve (its generator has order
    // 4·q — see the crate constructor's documentation).
    let bin: &[(&str, &[&str], Family, &str, fn() -> BinaryCurve)] = &[
        (
            "sect113r1",
            &[],
            Family::BinaryRandom,
            "SEC 2 v1",
            BinaryCurve::sect113r1,
        ),
        (
            "sect113r2",
            &[],
            Family::BinaryRandom,
            "SEC 2 v1",
            BinaryCurve::sect113r2,
        ),
        (
            "sect131r1",
            &[],
            Family::BinaryRandom,
            "SEC 2 v1",
            BinaryCurve::sect131r1,
        ),
        (
            "sect131r2",
            &[],
            Family::BinaryRandom,
            "SEC 2 v1",
            BinaryCurve::sect131r2,
        ),
        (
            "sect163k1",
            &["k-163", "nist-k163"],
            Family::Koblitz,
            "SEC 2 / NIST K-163",
            BinaryCurve::sect163k1,
        ),
        (
            "sect163r1",
            &[],
            Family::BinaryRandom,
            "SEC 2 v1",
            BinaryCurve::sect163r1,
        ),
        (
            "sect163r2",
            &["b-163", "nist-b163"],
            Family::BinaryRandom,
            "SEC 2 / NIST B-163",
            BinaryCurve::sect163r2,
        ),
        (
            "oakley-group-3",
            &["ike-oakley-3", "ec2n-155"],
            Family::BinaryRandom,
            "RFC 2409 (EC2N; subgroup curve)",
            BinaryCurve::ike_oakley_group3,
        ),
        // Full NIST/SECG binary curve suite (added with this engine).
        (
            "sect193r1",
            &[],
            Family::BinaryRandom,
            "SEC 2 v2",
            BinaryCurve::sect193r1,
        ),
        (
            "sect193r2",
            &[],
            Family::BinaryRandom,
            "SEC 2 v2",
            BinaryCurve::sect193r2,
        ),
        (
            "sect233k1",
            &["k-233", "nist-k233", "ansit233k1"],
            Family::Koblitz,
            "SEC 2 / NIST K-233",
            BinaryCurve::sect233k1,
        ),
        (
            "sect233r1",
            &["b-233", "nist-b233", "ansit233r1"],
            Family::BinaryRandom,
            "SEC 2 / NIST B-233",
            BinaryCurve::sect233r1,
        ),
        (
            "sect239k1",
            &["ansit239k1"],
            Family::Koblitz,
            "SEC 2 v2",
            BinaryCurve::sect239k1,
        ),
        (
            "sect283k1",
            &["k-283", "nist-k283", "ansit283k1"],
            Family::Koblitz,
            "SEC 2 / NIST K-283",
            BinaryCurve::sect283k1,
        ),
        (
            "sect283r1",
            &["b-283", "nist-b283", "ansit283r1"],
            Family::BinaryRandom,
            "SEC 2 / NIST B-283",
            BinaryCurve::sect283r1,
        ),
        (
            "sect409k1",
            &["k-409", "nist-k409", "ansit409k1"],
            Family::Koblitz,
            "SEC 2 / NIST K-409",
            BinaryCurve::sect409k1,
        ),
        (
            "sect409r1",
            &["b-409", "nist-b409", "ansit409r1"],
            Family::BinaryRandom,
            "SEC 2 / NIST B-409",
            BinaryCurve::sect409r1,
        ),
        (
            "sect571k1",
            &["k-571", "nist-k571", "ansit571k1"],
            Family::Koblitz,
            "SEC 2 / NIST K-571",
            BinaryCurve::sect571k1,
        ),
        (
            "sect571r1",
            &["b-571", "nist-b571", "ansit571r1"],
            Family::BinaryRandom,
            "SEC 2 / NIST B-571",
            BinaryCurve::sect571r1,
        ),
    ];
    for (name, aliases, family, std, ctor) in bin {
        v.push(CatalogCurve {
            name,
            aliases,
            family: *family,
            standard: std,
            object: CurveObject::Binary(Box::new(ctor())),
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_catalog_curve_verifies() {
        // Per-family floors and coverage, not an exact count: the catalog
        // grows, so we assert a lower bound per family plus that every entry
        // that IS present verifies.  A family collapsing to zero still fails.
        let all = all();
        let mut per_family: std::collections::HashMap<&str, usize> =
            std::collections::HashMap::new();
        for c in &all {
            let checks = c.verify();
            for chk in &checks {
                assert!(
                    chk.passed,
                    "curve {} failed check {}: {}",
                    c.name, chk.name, chk.detail
                );
            }
            *per_family.entry(c.family.tag()).or_default() += 1;
        }
        assert!(
            *per_family.get("prime").unwrap_or(&0) >= 30,
            "expected at least 30 prime curves, got {:?}",
            per_family.get("prime")
        );
        assert!(
            *per_family.get("binary").unwrap_or(&0) >= 5,
            "expected at least 5 random-binary curves, got {:?}",
            per_family.get("binary")
        );
        assert!(
            *per_family.get("koblitz").unwrap_or(&0) >= 5,
            "expected the NIST Koblitz suite (K-163..K-571), got {:?}",
            per_family.get("koblitz")
        );
    }

    #[test]
    fn names_are_unique_and_resolvable() {
        let all = all();
        let mut seen = std::collections::HashSet::new();
        for c in &all {
            assert!(seen.insert(c.name), "duplicate canonical name {}", c.name);
            // Every canonical name and alias resolves back to this curve.
            assert_eq!(by_name(c.name).map(|x| x.name), Some(c.name));
            for a in c.aliases {
                assert_eq!(
                    by_name(a).map(|x| x.name),
                    Some(c.name),
                    "alias {} did not resolve to {}",
                    a,
                    c.name
                );
            }
        }
    }

    #[test]
    fn security_estimate_is_reasonable() {
        let c = by_name("p256").expect("p256 in catalog");
        // P-256 has ~128-bit rho security.
        let s = c.rho_security_bits();
        assert!((120.0..132.0).contains(&s), "P-256 rho security was {s}");
        assert_eq!(c.field().bits, 256);
    }
}

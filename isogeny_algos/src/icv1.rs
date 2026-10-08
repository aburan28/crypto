//! ICV1 identity of a prime-field curve y^2 = x^3 + a x + b over GF(p), as the crypto repository
//! names curves (`docs/curves/ICV1.md` there; AGENTS.md section 11):
//!
//! ```text
//! ICV1:fp-<p>:<trace>:<order>:<j>:<end>:<level>:<path>:<model12>
//! slug: icv1-fp<bits of p>-t<trace>-<model8>   (negative trace written tm<|t|>)
//! ```
//!
//! `model12` / `model8` are prefixes of SHA-256 of the model JSON
//! `{"a":"<a>","b":"<b>","field":"fp-<p>","form":"y^2=x^3+a*x+b","p":"<p>","v":"1"}` (decimal,
//! reduced mod p, sorted keys, no whitespace). This mirrors `cryptanalysis::curve_id::prime` in
//! the crypto crate and is pinned to the same vectors (`tests/icv1.rs`). `end` and `level` are
//! `unk` and `path` is `r`: this module certifies neither the endomorphism ring nor a volcano
//! position. The identity names a curve *model*; joins across repositories also need the EC1
//! curve UID of the exact record (generator, subgroup), which this crate does not produce.
use crate::int::Int;

/// The identity strings of one curve model.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CurveId {
    /// the full record string, `ICV1:fp-...`
    pub icv1: String,
    /// the name used in text, `icv1-fp...`
    pub slug: String,
    /// the hashed model JSON
    pub model_json: String,
}

/// ICV1 identity of y^2 = x^3 + a x + b over GF(p) with #E(GF(p)) = `order` (the whole group).
/// None if the curve is singular, p is not an odd prime above 3, or `order` violates the Hasse
/// bound (a subgroup order passed as the group's).
pub fn prime(p: &Int, a: &Int, b: &Int, order: &Int) -> Option<CurveId> {
    if *p <= Int::from(3i64) || !p.is_probable_prime() {
        return None;
    }
    let (a, b) = (a.modulo(p), b.modulo(p));
    let a3 = a.pow_mod(&Int::from(3i64), p);
    let disc = (&(&Int::from(4i64) * &a3) + &(&Int::from(27i64) * &(&b * &b))).modulo(p);
    let trace = &(p + &Int::one()) - order;
    if disc.is_zero() || &trace * &trace > &Int::from(4i64) * p {
        return None;
    }
    let field = format!("fp-{p}");
    let model_json = format!(
        r#"{{"a":"{a}","b":"{b}","field":"{field}","form":"y^2=x^3+a*x+b","p":"{p}","v":"1"}}"#
    );
    let j = (&(&Int::from(1728i64 * 4) * &a3) * &disc.inv_mod(p)?).modulo(p);
    let model = crate::sha256::sha256_hex(model_json.as_bytes());
    let icv1 = format!(
        "ICV1:{field}:{trace}:{order}:{j}:unk:unk:r:{}",
        &model[..12]
    );
    let t = if trace.is_neg() {
        format!("tm{}", trace.abs())
    } else {
        format!("t{trace}")
    };
    let slug = format!("icv1-fp{}-{t}-{}", p.bits(), &model[..8]);
    Some(CurveId {
        icv1,
        slug,
        model_json,
    })
}

//! Exact same-field changes of coordinates; exceptional affine points are
//! interpreted on smooth projective curves, never dropped from the map.
use super::checker::{number, ModelError};
use crypto_lib::utils::mod_inverse;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use serde_json::{json, Value};

pub(super) fn normalize(v: &Value) -> Result<Option<(Value, Value)>, ModelError> {
    let form = v["form"].as_str().unwrap_or("");
    if ![
        "B*y^2=x^3+A*x^2+x",
        "a*x^2+y^2=1+d*x^2*y^2",
        "x^2+y^2=c^2*(1+d*x^2*y^2)",
    ]
    .contains(&form)
    {
        return Ok(None);
    }
    let get = |k: &str| number(v[k].as_str().ok_or_else(|| format!("missing {k}"))?);
    let p = get("p")?;
    if p <= BigUint::from(3u32) || p.bits() > 4096 {
        return Err(ModelError::Unsupported(
            "model conversion requires prime characteristic > 3".into(),
        ));
    }
    if v["field"] != format!("fp-{p}") {
        return Err("inconsistent field label".to_string().into());
    }
    let inv = |x: &BigUint| {
        mod_inverse(x, &p).ok_or_else(|| "noninvertible model conversion denominator".to_string())
    };
    let sub = |a: &BigUint, b: &BigUint| (a + &p - b % &p) % &p;
    let canonical = |x: &BigUint| {
        if x < &p {
            Ok(())
        } else {
            Err("noncanonical model coefficient".to_string())
        }
    };
    let (a, b, scale) = if form == "B*y^2=x^3+A*x^2+x" {
        let a = get("A")?;
        let b = get("B")?;
        canonical(&a)?;
        canonical(&b)?;
        (a, b, BigUint::one())
    } else {
        let d = get("d")?;
        canonical(&d)?;
        let (a, d, c) = if form.starts_with("a*") {
            let a = get("a")?;
            canonical(&a)?;
            (a, d, BigUint::one())
        } else {
            let c = get("c")?;
            canonical(&c)?;
            if c.is_zero() {
                return Err("singular Edwards scale".to_string().into());
            }
            (
                BigUint::one(),
                d * c.modpow(&BigUint::from(4u32), &p) % &p,
                c,
            )
        };
        if a.is_zero() || d.is_zero() || a == d {
            return Err("singular Edwards model".to_string().into());
        }
        let den = inv(&sub(&a, &d))?;
        (
            (BigUint::from(2u32) * (&a + &d) * &den) % &p,
            (BigUint::from(4u32) * den) % &p,
            c,
        )
    };
    if b.is_zero() || sub(&(&a * &a % &p), &BigUint::from(4u32)).is_zero() {
        return Err("singular Montgomery model".to_string().into());
    }
    let third = inv(&BigUint::from(3u32))?;
    let shift = &a * &third % &p;
    let bi = inv(&b)?;
    let sa = sub(&BigUint::one(), &(&a * &a * &third % &p)) * &bi * &bi % &p;
    let sb = sub(
        &(BigUint::from(2u32) * &shift * &shift * &shift % &p),
        &shift,
    ) * &bi
        * &bi
        * &bi
        % &p;
    // Coefficient identities after x_M=B*X-A/3, y_M=B*Y.
    if (&sa * &b * &b + &a * &a * &third) % &p != BigUint::one()
        || (&sb * &b * &b * &b + &shift) % &p != BigUint::from(2u32) * &shift * &shift * &shift % &p
    {
        return Err("model conversion coefficient identity failed"
            .to_string()
            .into());
    }
    let short = json!({"v":"1","form":"y^2=x^3+a*x+b","p":p.to_string(),"field":format!("fp-{p}"),"a":sa.to_string(),"b":sb.to_string()});
    let map = json!({"schema_version":"elliptic-model-map/v1","source_model":short,"target_model":v,
        "kind":if form.starts_with("B*"){"montgomery_from_short"}else{"edwards_from_short"},
        "montgomery_A":format!("0x{a:x}"),"montgomery_B":format!("0x{b:x}"),"shift_A_over_3":format!("0x{shift:x}"),"edwards_scale":format!("0x{scale:x}"),
        "forward":if form.starts_with("B*"){json!({"x":"B*X-A/3","y":"B*Y"})}else{json!({"x":"c*(B*X-A/3)/(B*Y)","y":"c*(B*X-A/3-1)/(B*X-A/3+1)"})},
        "inverse":if form.starts_with("B*"){json!({"X":"(x+A/3)/B","Y":"y/B"})}else{json!({"X":"((c+y)/(c-y)+A/3)/B","Y":"c*(c+y)/(B*x*(c-y))"})},
        "exceptional_points":"unique extension of the birational map to smooth projective models","degree":1,"status":"verified_algebraic_identity"});
    Ok(Some((short, map)))
}

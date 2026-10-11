use crypto_lib::{
    binary_ecc::{F2mElement, IrreduciblePoly},
    cryptanalysis::ecc2k130_guard::is_probable_prime,
    hash::sha256::sha256,
    utils::mod_inverse,
};
use num_bigint::BigUint;
use num_traits::{One, Zero};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::collections::HashSet;

type Poly = Vec<BigUint>;
type Check<T> = Result<T, String>;
pub(super) fn digest(bytes: &[u8]) -> String {
    sha256(bytes).iter().map(|b| format!("{b:02x}")).collect()
}
pub(super) fn number(s: &str) -> Check<BigUint> {
    if s.len() > 1300 {
        return Err("integer exceeds field-size limit".into());
    }
    let (s, radix) = s.strip_prefix("0x").map_or((s, 10), |s| (s, 16));
    BigUint::parse_bytes(s.as_bytes(), radix).ok_or_else(|| "invalid unsigned integer".into())
}
fn get_number(v: &Value, key: &str) -> Check<BigUint> {
    number(
        v[key]
            .as_str()
            .ok_or_else(|| format!("missing string {key}"))?,
    )
}
/// An extension's coefficient list: exactly `k` numbers, each below `p`.
fn digit_list(v: &Value, key: &str, p: &BigUint, k: usize) -> Check<Poly> {
    let list = v[key]
        .as_array()
        .filter(|list| list.len() == k)
        .ok_or_else(|| format!("{key} is not a list of {k} numbers"))?;
    list.iter()
        .map(|x| {
            let n = number(x.as_str().ok_or_else(|| format!("{key}: not a string"))?)?;
            if &n < p {
                Ok(n)
            } else {
                Err("noncanonical coefficient".to_string())
            }
        })
        .collect()
}
fn trim(mut p: Poly) -> Poly {
    while p.last().is_some_and(Zero::is_zero) {
        p.pop();
    }
    p
}
fn bit_rem(mut a: BigUint, b: &BigUint) -> BigUint {
    while !a.is_zero() && a.bits() >= b.bits() {
        let shift = (a.bits() - b.bits()) as usize;
        a ^= b << shift;
    }
    a
}
fn bit_gcd(mut a: BigUint, mut b: BigUint) -> BigUint {
    while !b.is_zero() {
        let r = bit_rem(a, &b);
        a = b;
        b = r;
    }
    a
}
/// Rabin's checkpoints for degree `m`: `m/l` for every prime `l | m`.
fn prime_quotients(m: usize) -> HashSet<usize> {
    let mut n = m;
    let mut checkpoints = HashSet::new();
    let mut l = 2;
    while l * l <= n {
        if n.is_multiple_of(l) {
            checkpoints.insert(m / l);
            while n.is_multiple_of(l) {
                n /= l;
            }
        }
        l += 1;
    }
    if n > 1 {
        checkpoints.insert(m / n);
    }
    checkpoints
}
fn irreducible(q: &BigUint, poly: &IrreduciblePoly) -> bool {
    // Rabin: z^(2^m)=z, and gcd(z^(2^(m/l))-z,q)=1 for every prime l|m.
    let m = poly.degree;
    let checkpoints = prime_quotients(m as usize);
    let z = bit_rem(BigUint::from(2u32), q);
    let mut t = F2mElement::from_biguint(&z, m);
    for i in 1..=m {
        t = t.square(poly);
        if checkpoints.contains(&(i as usize))
            && bit_gcd(t.to_biguint() ^ &z, q.clone()) != BigUint::one()
        {
            return false;
        }
    }
    t.to_biguint() == z
}
/// Rabin over `GF(p)`: monic `f = t^k + low(t)` is irreducible iff
/// `t^(p^k) = t` modulo `f`, and `gcd(t^(p^(k/l)) - t, f) = 1` for every
/// prime `l | k`.
pub(super) fn irreducible_over_prime(p: &BigUint, low: &[BigUint]) -> bool {
    let k = low.len();
    if k < 2 {
        return k == 1;
    }
    let field = Field::Extension {
        p: p.clone(),
        low: low.to_vec(),
    };
    let gf_p = Field::Prime(p.clone());
    let mut f = low.to_vec();
    f.push(BigUint::one());
    // The packed element `t` is the integer `p`.
    let t = p.clone();
    let checkpoints = prime_quotients(k);
    let mut x = t.clone();
    for i in 1..=k {
        x = field.pow(&x, p);
        if checkpoints.contains(&i) {
            let diff = trim(unpack(p, k, &field.add(&x, &field.neg(&t))));
            if !matches!(gf_p.gcd(f.clone(), diff), Ok(g) if g.len() == 1) {
                return false;
            }
        }
    }
    x == t
}

/// The base-`p` digits of a packed `GF(p^k)` element, `k` of them.
fn unpack(p: &BigUint, k: usize, a: &BigUint) -> Poly {
    let mut a = a.clone();
    (0..k)
        .map(|_| {
            let d = &a % p;
            a /= p;
            d
        })
        .collect()
}
fn pack(p: &BigUint, digits: &[BigUint]) -> BigUint {
    digits
        .iter()
        .rev()
        .fold(BigUint::zero(), |acc, d| acc * p + d)
}
/// `a*b` in `GF(p)[t]/(t^k + low(t))`, on digit vectors of length `k`.
fn ext_mul(p: &BigUint, low: &[BigUint], a: &[BigUint], b: &[BigUint]) -> Poly {
    let k = low.len();
    let mut c = vec![BigUint::zero(); 2 * k - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            c[i + j] = (&c[i + j] + x * y) % p;
        }
    }
    // t^k = -low(t), from the top down.
    for i in (k..2 * k - 1).rev() {
        let top = std::mem::take(&mut c[i]);
        for (j, l) in low.iter().enumerate() {
            c[i - k + j] = (&c[i - k + j] + p - &top * l % p) % p;
        }
    }
    c.truncate(k);
    c
}

pub(super) enum Field {
    Prime(BigUint),
    Binary(IrreduciblePoly),
    /// `GF(p^k) = GF(p)[t]/(t^k + low(t))`, ICV1's extension part. The
    /// element `e_0 + e_1*t + ... + e_{k-1}*t^(k-1)` is packed as the
    /// integer `e_0 + e_1*p + ... + e_{k-1}*p^(k-1)`, as a binary element
    /// packs its polynomial-basis bits.
    Extension {
        p: BigUint,
        low: Vec<BigUint>,
    },
}
impl Field {
    pub(super) fn add(&self, a: &BigUint, b: &BigUint) -> BigUint {
        match self {
            Self::Prime(p) => (a + b) % p,
            Self::Binary(_) => a ^ b,
            Self::Extension { p, low } => {
                let (a, b) = (unpack(p, low.len(), a), unpack(p, low.len(), b));
                let sum: Poly = a.iter().zip(&b).map(|(x, y)| (x + y) % p).collect();
                pack(p, &sum)
            }
        }
    }
    fn neg(&self, a: &BigUint) -> BigUint {
        match self {
            Self::Prime(p) => (p - a) % p,
            Self::Binary(_) => a.clone(),
            Self::Extension { p, low } => {
                let digits: Poly = unpack(p, low.len(), a)
                    .iter()
                    .map(|x| (p - x) % p)
                    .collect();
                pack(p, &digits)
            }
        }
    }
    pub(super) fn mul(&self, a: &BigUint, b: &BigUint) -> BigUint {
        match self {
            Self::Prime(p) => a * b % p,
            Self::Binary(p) => F2mElement::from_biguint(a, p.degree)
                .mul(&F2mElement::from_biguint(b, p.degree), p)
                .to_biguint(),
            Self::Extension { p, low } => {
                let k = low.len();
                pack(p, &ext_mul(p, low, &unpack(p, k, a), &unpack(p, k, b)))
            }
        }
    }
    /// The integer `n` in the prime subfield.
    pub(super) fn constant(&self, n: u32) -> BigUint {
        match self {
            Self::Prime(p) | Self::Extension { p, .. } => BigUint::from(n) % p,
            Self::Binary(_) => BigUint::from(n & 1),
        }
    }
    fn pow(&self, a: &BigUint, e: &BigUint) -> BigUint {
        let mut r = self.constant(1);
        for i in (0..e.bits()).rev() {
            r = self.mul(&r, &r);
            if e.bit(i) {
                r = self.mul(&r, a);
            }
        }
        r
    }
    fn inverse(&self, a: &BigUint) -> Check<BigUint> {
        match self {
            Self::Prime(p) => mod_inverse(a, p),
            Self::Binary(p) => F2mElement::from_biguint(a, p.degree)
                .flt_inverse(p)
                .map(|a| a.to_biguint()),
            // a^(p^k - 2), by Fermat in the multiplicative group.
            Self::Extension { p, low } => {
                (!a.is_zero()).then(|| self.pow(a, &(p.pow(low.len() as u32) - 2u32)))
            }
        }
        .ok_or_else(|| "noninvertible field element".into())
    }
    fn canonical(&self, a: &BigUint) -> bool {
        match self {
            Self::Prime(p) => a < p,
            Self::Binary(p) => a.bits() <= p.degree as u64,
            Self::Extension { p, low } => a < &p.pow(low.len() as u32),
        }
    }
    fn add_poly(&self, a: &[BigUint], b: &[BigUint]) -> Poly {
        let mut c = vec![BigUint::zero(); a.len().max(b.len())];
        for (i, x) in a.iter().enumerate() {
            c[i] = x.clone();
        }
        for (i, x) in b.iter().enumerate() {
            c[i] = self.add(&c[i], x);
        }
        trim(c)
    }
    fn sub_poly(&self, a: &[BigUint], b: &[BigUint]) -> Poly {
        self.add_poly(a, &b.iter().map(|v| self.neg(v)).collect::<Poly>())
    }
    fn mul_poly(&self, a: &[BigUint], b: &[BigUint]) -> Poly {
        if a.is_empty() || b.is_empty() {
            return vec![];
        }
        let mut c = vec![BigUint::zero(); a.len() + b.len() - 1];
        for (i, x) in a.iter().enumerate() {
            for (j, y) in b.iter().enumerate() {
                c[i + j] = self.add(&c[i + j], &self.mul(x, y));
            }
        }
        trim(c)
    }
    fn gcd(&self, mut a: Poly, mut b: Poly) -> Check<Poly> {
        while !b.is_empty() {
            let inv = self.inverse(b.last().unwrap())?;
            while !a.is_empty() && a.len() >= b.len() {
                let offset = a.len() - b.len();
                let c = self.mul(a.last().unwrap(), &inv);
                for (i, v) in b.iter().enumerate() {
                    a[i + offset] = self.add(&a[i + offset], &self.neg(&self.mul(&c, v)));
                }
                a = trim(a);
            }
            (a, b) = (b, a);
        }
        Ok(a)
    }
    pub(super) fn eval(&self, p: &[BigUint], u: &BigUint) -> BigUint {
        p.iter()
            .rev()
            .fold(BigUint::zero(), |acc, c| self.add(&self.mul(&acc, u), c))
    }
}

pub(super) struct Model {
    pub field: Field,
    pub a: BigUint,
    pub b: BigUint,
    pub normalization: Option<Value>,
}
pub(super) enum ModelError {
    Unsupported(String),
    Invalid(String),
}
impl From<String> for ModelError {
    fn from(s: String) -> Self {
        Self::Invalid(s)
    }
}
pub(super) fn model(v: &Value) -> Result<Model, ModelError> {
    if v["v"] != "1" {
        return Err(ModelError::Unsupported("unsupported model version".into()));
    }
    if let Some((short, map)) = super::models::normalize(v)? {
        let mut m = model(&short)?;
        m.normalization = Some(map);
        return Ok(m);
    }
    let form = v["form"].as_str().unwrap_or("");
    let field = match form {
        // ICV1's extension part: GF(p)[t]/(f), f = t^k + c_{k-1}t^(k-1) +
        // ... + c_0 monic, coefficients listed over the basis 1, t, ...
        "y^2=x^3+a*x+b" if v.get("k").is_some() => {
            let p = get_number(v, "p")?;
            let k: usize = v["k"]
                .as_str()
                .and_then(|k| k.parse().ok())
                .filter(|k| (2..=4096).contains(k))
                .ok_or_else(|| "invalid extension degree".to_string())?;
            if p < BigUint::from(2u32) || p.bits() * k as u64 > 4096 {
                return Err("invalid field size".to_string().into());
            }
            if p <= BigUint::from(3u32) {
                return Err(ModelError::Unsupported(
                    "prime characteristic 2 or 3".into(),
                ));
            }
            if !is_probable_prime(&p) {
                return Err("modulus failed probable-prime screen".to_string().into());
            }
            let low = digit_list(v, "modulus", &p, k)?;
            let listed: Vec<String> = low.iter().map(|c| c.to_string()).collect();
            let label = format!(
                "fpk-{p}-{k}-{}",
                &digest(format!("fpk-modulus:{p}:{}", listed.join(",")).as_bytes())[..8]
            );
            if v["field"] != label {
                return Err("inconsistent extension field label".to_string().into());
            }
            if !irreducible_over_prime(&p, &low) {
                return Err("extension modulus is reducible".to_string().into());
            }
            Field::Extension { p, low }
        }
        "y^2=x^3+a*x+b" => {
            let p = get_number(v, "p")?;
            if p.bits() > 4096 || p < BigUint::from(2u32) {
                return Err("invalid field size".to_string().into());
            }
            if p <= BigUint::from(3u32) {
                return Err(ModelError::Unsupported(
                    "prime characteristic 2 or 3".into(),
                ));
            }
            if !is_probable_prime(&p) {
                return Err("modulus failed probable-prime screen".to_string().into());
            }
            if v["field"] != format!("fp-{p}") {
                return Err("inconsistent prime field label".to_string().into());
            }
            Field::Prime(p)
        }
        "y^2+xy=x^3+a*x^2+b" => {
            let q = get_number(v, "modulus")?;
            if q.bits() < 2 || q.bits() > 4097 || !q.bit(0) {
                return Err("invalid binary modulus".to_string().into());
            }
            let m = (q.bits() - 1) as u32;
            let p = IrreduciblePoly {
                degree: m,
                low_terms: (0..m).filter(|i| q.bit(*i as u64)).collect(),
            };
            let label = format!(
                "f2m-{m}-{}",
                &digest(format!("f2m-modulus:0x{q:x}").as_bytes())[..8]
            );
            if v["field"] != label {
                return Err("inconsistent binary field label".to_string().into());
            }
            if !irreducible(&q, &p) {
                return Err("binary modulus is reducible".to_string().into());
            }
            Field::Binary(p)
        }
        _ => {
            return Err(ModelError::Unsupported(
                "unsupported elliptic model form".into(),
            ))
        }
    };
    let (a, b) = match &field {
        Field::Extension { p, low } => (
            pack(p, &digit_list(v, "a", p, low.len())?),
            pack(p, &digit_list(v, "b", p, low.len())?),
        ),
        _ => (get_number(v, "a")?, get_number(v, "b")?),
    };
    if !field.canonical(&a) || !field.canonical(&b) {
        return Err("noncanonical coefficient".to_string().into());
    }
    let nonsingular = match &field {
        Field::Prime(_) | Field::Extension { .. } => !field
            .add(
                &field.mul(&field.constant(4), &field.mul(&field.mul(&a, &a), &a)),
                &field.mul(&field.constant(27), &field.mul(&b, &b)),
            )
            .is_zero(),
        Field::Binary(_) => !b.is_zero(),
    };
    if !nonsingular {
        return Err("singular elliptic model".to_string().into());
    }
    Ok(Model {
        field,
        a,
        b,
        normalization: None,
    })
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Certificate {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub target_model_map: Option<Value>,
    pub construction: String,
    pub genus: u32,
    pub degree: u32,
    pub h: Vec<String>,
    pub f: Vec<String>,
    pub x: Vec<String>,
    pub y_v: Vec<String>,
    pub y_0: Vec<String>,
}
fn encoded(p: Poly) -> Vec<String> {
    trim(p).iter().map(|a| format!("0x{a:x}")).collect()
}
pub(super) fn decoded(p: &[String], k: &Field) -> Check<Poly> {
    if p.len() > 16 {
        return Err("certificate polynomial exceeds supported degree".into());
    }
    let p: Poly = p.iter().map(|a| number(a)).collect::<Check<_>>()?;
    if p.iter().any(|a| !k.canonical(a)) || p.last().is_some_and(Zero::is_zero) {
        return Err("noncanonical certificate polynomial".into());
    }
    Ok(p)
}
fn monomial(degree: usize) -> Poly {
    let mut p = vec![BigUint::zero(); degree + 1];
    p[degree] = BigUint::one();
    p
}
pub(super) fn construct(m: &Model) -> Certificate {
    let k = &m.field;
    let (name, g, degree, h, f, x, y_v, y_0) = match k {
        // Characteristic p > 3, over GF(p) or GF(p^k): the constants 0..3
        // are distinct field elements.
        Field::Prime(_) | Field::Extension { .. } => {
            let cubic = vec![m.b.clone(), m.a.clone(), BigUint::zero(), BigUint::one()];
            let c = (0..4u32)
                .map(BigUint::from)
                .find(|c| !k.eval(&cubic, c).is_zero())
                .unwrap();
            let x = vec![c, BigUint::zero(), BigUint::one()];
            let x2 = k.mul_poly(&x, &x);
            let f = k.add_poly(
                &k.add_poly(
                    &k.mul_poly(&x2, &x),
                    &k.mul_poly(std::slice::from_ref(&m.a), &x),
                ),
                std::slice::from_ref(&m.b),
            );
            (
                "prime_quadratic_pullback_v1",
                2,
                2,
                vec![],
                f,
                x,
                monomial(0),
                vec![],
            )
        }
        Field::Binary(p) => {
            let d = F2mElement::from_biguint(&m.b, p.degree)
                .square_k_times(p.degree - 1, p)
                .to_biguint();
            let mut f = monomial(7);
            f[1] = d.clone();
            f[4] = m.a.clone();
            (
                "binary_cubic_pullback_v1",
                3,
                3,
                monomial(2),
                f,
                monomial(3),
                monomial(1),
                vec![d],
            )
        }
    };
    Certificate {
        target_model_map: m.normalization.clone(),
        construction: name.into(),
        genus: g,
        degree,
        h: encoded(h),
        f: encoded(f),
        x: encoded(x),
        y_v: encoded(y_v),
        y_0: encoded(y_0),
    }
}

pub(super) fn verify(m: &Model, c: &Certificate) -> Check<()> {
    if c.target_model_map != m.normalization {
        return Err("target-model isomorphism certificate mismatch".into());
    }
    let k = &m.field;
    let h = decoded(&c.h, k)?;
    let f = decoded(&c.f, k)?;
    let x = decoded(&c.x, k)?;
    let a = decoded(&c.y_v, k)?;
    let b = decoded(&c.y_0, k)?;
    // Reduce the supplied map's target equation modulo v^2=f-h*v.
    let aa = k.mul_poly(&a, &a);
    let ab = k.mul_poly(&a, &b);
    let mut linear = k.sub_poly(&k.add_poly(&ab, &ab), &k.mul_poly(&aa, &h));
    let mut constant = k.add_poly(&k.mul_poly(&aa, &f), &k.mul_poly(&b, &b));
    let x2 = k.mul_poly(&x, &x);
    let mut rhs = k.add_poly(&k.mul_poly(&x2, &x), std::slice::from_ref(&m.b));
    match k {
        Field::Prime(_) | Field::Extension { .. } => {
            rhs = k.add_poly(&rhs, &k.mul_poly(std::slice::from_ref(&m.a), &x))
        }
        Field::Binary(_) => {
            linear = k.add_poly(&linear, &k.mul_poly(&x, &a));
            constant = k.add_poly(&constant, &k.mul_poly(&x, &b));
            rhs = k.add_poly(&rhs, &k.mul_poly(std::slice::from_ref(&m.a), &x2));
        }
    }
    if !linear.is_empty() || constant != rhs {
        return Err("map fails polynomial substitution".into());
    }
    // Equation identities alone do not certify genus, separability or map degree.
    match k {
        Field::Prime(p) | Field::Extension { p, .. } => {
            if c.construction != "prime_quadratic_pullback_v1"
                || c.genus != 2
                || c.degree != 2
                || !h.is_empty()
                || !b.is_empty()
                || a != monomial(0)
                || x.len() != 3
                || !x[1].is_zero()
                || x[2] != BigUint::one()
                || f.len() != 7
                || f[6] != BigUint::one()
                || f[0].is_zero()
            {
                return Err("prime genus/degree hypotheses fail".into());
            }
            let derivative = trim(
                (1..f.len())
                    .map(|i| k.mul(&f[i], &(BigUint::from(i) % p)))
                    .collect(),
            );
            if k.gcd(f, derivative)?.len() != 1 {
                return Err("cover sextic is not squarefree".into());
            }
            // Squarefree sextic => genus 2; quadratic coordinate-line map => degree 2.
        }
        Field::Binary(_) => {
            if c.construction != "binary_cubic_pullback_v1"
                || c.genus != 3
                || c.degree != 3
                || h != monomial(2)
                || x != monomial(3)
                || a != monomial(1)
                || b.len() != 1
                || b[0].is_zero()
                || k.mul(&b[0], &b[0]) != m.b
            {
                return Err("binary genus/degree hypotheses fail".into());
            }
            let mut expected = monomial(7);
            expected[1] = b[0].clone();
            expected[4] = m.a.clone();
            if f != expected {
                return Err("binary Artin-Schreier pole conditions fail".into());
            }
            // v/u^2 has AS RHS u^3+a+d/u^3: two odd poles of order 3 => genus 3.
            // The coordinate-line map u -> u^3 is separable of degree 3.
        }
    }
    Ok(())
}

pub(super) fn catalog(bytes: &[u8]) -> Check<Value> {
    let registry: Value = serde_json::from_slice(bytes).map_err(|e| e.to_string())?;
    if registry["schema_version"] != 1 {
        return Err("unsupported registry schema".into());
    }
    let rows = registry["curves"]
        .as_array()
        .ok_or("missing curves array")?;
    let mut seen = HashSet::new();
    let mut curves = vec![];
    let (mut verified, mut unsupported, mut invalid) = (0, 0, 0);
    for row in rows {
        let slug = row["slug"].as_str().ok_or("missing curve slug")?;
        if !seen.insert(slug) {
            return Err(format!("duplicate slug {slug}"));
        }
        let raw = row["model_json"]
            .as_str()
            .ok_or("missing exact model_json")?;
        let hash = digest(raw.as_bytes());
        let icv1 = row["icv1"].as_str().ok_or("missing ICV1 identity")?;
        if !slug.ends_with(&format!("-{}", &hash[..8]))
            || !icv1.ends_with(&format!(":{}", &hash[..12]))
        {
            return Err(format!("model identity mismatch for {slug}"));
        }
        let value: Value = serde_json::from_str(raw).map_err(|e| format!("{slug}: {e}"))?;
        let mut finding = json!({"slug":slug,"model_sha256":hash,"exists":null,"certificate":null,
            "map_direction":"H -> E","field_relation":null,"field_check":null,
            "minimal_genus":null,"minimal_degree":null,"subfield_descent":"not_tested",
            "subgroup_transfer":"not_tested","dlp_advantage":null});
        match model(&value) {
            Ok(m) => {
                let cert = construct(&m);
                verify(&m, &cert)
                    .map_err(|e| format!("internal certificate failure for {slug}: {e}"))?;
                finding["status"] = json!("verified_over_declared_field");
                finding["exists"] = json!(true);
                finding["certificate"] = json!(cert);
                finding["field_relation"] = json!("same_field");
                finding["field_check"] = json!(match m.field {
                    Field::Prime(_) =>
                        "20_fixed_base_probable_prime_screen; primality remains an input assumption",
                    Field::Binary(_) => "irreducibility_proved_by_Rabin",
                    Field::Extension { .. } =>
                        "20_fixed_base_probable_prime_screen of p, whose primality remains an input assumption; irreducibility of the modulus over GF(p) proved by Rabin",
                });
                finding["reason"]=json!("polynomial substitution and genus/degree conditions checked; smooth projective models understood");
                verified += 1;
            }
            Err(ModelError::Unsupported(reason)) => {
                finding["status"] = json!("unsupported");
                finding["reason"] = json!(reason);
                unsupported += 1;
            }
            Err(ModelError::Invalid(reason)) => {
                finding["status"] = json!("invalid_input");
                finding["reason"] = json!(reason);
                invalid += 1;
            }
        }
        curves.push(finding);
    }
    Ok(
        json!({"schema_version":"curve-covers/v1","generated_by":"cargo run --bin curve_cover_check --",
        "registry_sha256":digest(bytes),
        "checker_source_sha256":digest(concat!(include_str!("checker.rs"),include_str!("models.rs"),include_str!("links.rs"),include_str!("../curve_cover_check.rs")).as_bytes()),
        "scope":"same-field existence for the supplied elliptic model; no optimality, descent, subgroup-map or DLP-cost claim",
        "coefficient_encoding":"0x field elements, ascending powers of u; H: v^2+h(u)*v=f(u); x=x(u), y=y_v(u)*v+y_0(u); a GF(p^k) element e_0+e_1*t+...+e_{k-1}*t^(k-1) is the integer e_0+e_1*p+...+e_{k-1}*p^(k-1)",
        "proof":"docs/curves/COVERS.md","summary":{"verified":verified,"unsupported":unsupported,"invalid_input":invalid},"curves":curves}),
    )
}

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
fn number(s: &str) -> Check<BigUint> {
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
fn irreducible(q: &BigUint, poly: &IrreduciblePoly) -> bool {
    // Rabin: z^(2^m)=z, and gcd(z^(2^(m/l))-z,q)=1 for every prime l|m.
    let m = poly.degree;
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
    let z = bit_rem(BigUint::from(2u32), q);
    let mut t = F2mElement::from_biguint(&z, m);
    for i in 1..=m {
        t = t.square(poly);
        if checkpoints.contains(&i) && bit_gcd(t.to_biguint() ^ &z, q.clone()) != BigUint::one() {
            return false;
        }
    }
    t.to_biguint() == z
}

pub(super) enum Field {
    Prime(BigUint),
    Binary(IrreduciblePoly),
}
impl Field {
    pub(super) fn add(&self, a: &BigUint, b: &BigUint) -> BigUint {
        match self {
            Self::Prime(p) => (a + b) % p,
            Self::Binary(_) => a ^ b,
        }
    }
    fn neg(&self, a: &BigUint) -> BigUint {
        match self {
            Self::Prime(p) => (p - a) % p,
            Self::Binary(_) => a.clone(),
        }
    }
    pub(super) fn mul(&self, a: &BigUint, b: &BigUint) -> BigUint {
        match self {
            Self::Prime(p) => a * b % p,
            Self::Binary(p) => F2mElement::from_biguint(a, p.degree)
                .mul(&F2mElement::from_biguint(b, p.degree), p)
                .to_biguint(),
        }
    }
    fn inverse(&self, a: &BigUint) -> Check<BigUint> {
        match self {
            Self::Prime(p) => mod_inverse(a, p),
            Self::Binary(p) => F2mElement::from_biguint(a, p.degree)
                .flt_inverse(p)
                .map(|a| a.to_biguint()),
        }
        .ok_or_else(|| "noninvertible field element".into())
    }
    fn canonical(&self, a: &BigUint) -> bool {
        match self {
            Self::Prime(p) => a < p,
            Self::Binary(p) => a.bits() <= p.degree as u64,
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
    let form = v["form"].as_str().unwrap_or("");
    let field = match form {
        // ICV1's extension part (B5a): its coefficients are lists over
        // GF(p^k), and no construction here covers them.
        "y^2=x^3+a*x+b" if v.get("k").is_some() => {
            return Err(ModelError::Unsupported(
                "an extension field GF(p^k): no cover construction for it here".into(),
            ));
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
    let a = get_number(v, "a")?;
    let b = get_number(v, "b")?;
    if !field.canonical(&a) || !field.canonical(&b) {
        return Err("noncanonical coefficient".to_string().into());
    }
    let nonsingular = match &field {
        Field::Prime(_) => !field
            .add(
                &field.mul(&BigUint::from(4u32), &field.mul(&field.mul(&a, &a), &a)),
                &field.mul(&BigUint::from(27u32), &field.mul(&b, &b)),
            )
            .is_zero(),
        Field::Binary(_) => !b.is_zero(),
    };
    if !nonsingular {
        return Err("singular elliptic model".to_string().into());
    }
    Ok(Model { field, a, b })
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Certificate {
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
        Field::Prime(_) => {
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
        Field::Prime(_) => rhs = k.add_poly(&rhs, &k.mul_poly(std::slice::from_ref(&m.a), &x)),
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
        Field::Prime(p) => {
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
        "checker_source_sha256":digest(concat!(include_str!("checker.rs"),include_str!("../curve_cover_check.rs")).as_bytes()),
        "scope":"same-field existence for the supplied elliptic model; no optimality, descent, subgroup-map or DLP-cost claim",
        "coefficient_encoding":"0x field elements, ascending powers of u; H: v^2+h(u)*v=f(u); x=x(u), y=y_v(u)*v+y_0(u)",
        "proof":"docs/curves/COVERS.md","summary":{"verified":verified,"unsupported":unsupported,"invalid_input":invalid},"curves":curves}),
    )
}

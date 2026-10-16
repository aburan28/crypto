//! Compact known-endomorphism certificates. Expanded-map challenge scores are separate.
use crate::identity_certificate::{self as idc, Expr, IdentityCertificate, IdentityStatement};
use num_bigint::{BigInt, BigUint, Sign};
use num_traits::{One, Zero};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{collections::BTreeMap, fs, path::Path, time::Instant};

pub const MAX_EXPONENT: u32 = 65536;
pub const P: &str =
    "115792089237316195423570985008687907853269984665640564039457584007908834671663";
pub const N: &str =
    "115792089237316195423570985008687907852837564279074904382605163141518161494337";
pub const GX: &str =
    "55066263022277343669578718895168534326250603453777594175500187360389116729240";
pub const GY: &str =
    "32670510020758816978083085130507043184471273380659243275938904335757337482424";
pub type Point = Option<[String; 2]>;
type Affine = Option<(BigUint, BigUint)>;

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct Domain {
    pub name: String,
    pub p: String,
    pub a: String,
    pub b: String,
    pub n: String,
    pub h: String,
    pub generator: Point,
    pub source_icv1: String,
    pub primality_basis: String,
}
pub fn domain() -> Domain {
    Domain { name:"secp256k1".into(), p:P.into(), a:"0".into(), b:"7".into(), n:N.into(), h:"1".into(),
        generator:Some([GX.into(),GY.into()]),
        source_icv1:"icv1-fp256-t432420386565659656852420866390673177327-76dadd18".into(),
        primality_basis:"trusted SEC 2 v2 section 2.4.1 parameters; p and n primality are not independently certified here".into() }
}
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct LucasLehmer {
    pub exponent: u32,
    pub transitions: u32,
    pub terminal: String,
    pub states_sha256: String,
    pub theorem: String,
}
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct Cornacchia {
    pub root: String,
    pub x: String,
    pub y: String,
    pub divisions: u32,
    pub division_sha256: String,
}
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct Evaluation {
    pub input: Point,
    pub output: Point,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Certificate {
    pub schema: String,
    pub representation: String,
    pub domain: Domain,
    pub exponent: u32,
    pub degree: String,
    pub a: String,
    pub b: String,
    pub beta: String,
    pub lambda: String,
    pub differential: String,
    pub subgroup_scalar: String,
    pub subgroup_inverse: String,
    pub primality: LucasLehmer,
    pub cornacchia: Cornacchia,
    pub identity: IdentityCertificate,
    pub evaluations: Vec<Evaluation>,
    pub tool: Value,
}
pub fn unsigned(s: &str) -> Result<BigUint, String> {
    if s.is_empty()
        || s.len() > 20000
        || !s.bytes().all(|c| c.is_ascii_digit())
        || (s.len() > 1 && s.starts_with('0'))
    {
        return Err("noncanonical bounded decimal integer".into());
    }
    BigUint::parse_bytes(s.as_bytes(), 10).ok_or("invalid integer".into())
}
pub fn signed(s: &str) -> Result<BigInt, String> {
    if let Some(t) = s.strip_prefix('-') {
        let n = unsigned(t)?;
        if n.is_zero() {
            return Err("negative zero".into());
        }
        Ok(-BigInt::from(n))
    } else {
        Ok(BigInt::from(unsigned(s)?))
    }
}
pub fn exponent_prime(e: u32) -> bool {
    e >= 3
        && e <= MAX_EXPONENT
        && e % 2 == 1
        && (2..)
            .take_while(|d| u64::from(*d) * u64::from(*d) <= u64::from(e))
            .all(|d| e % d != 0)
}
pub fn state_hash(h: &mut Sha256, index: u32, decimal: &str) {
    h.update(index.to_le_bytes());
    h.update((decimal.len() as u64).to_le_bytes());
    h.update(decimal.as_bytes());
}
// Producer uses Mersenne folding. The verifier uses another integer implementation
// and general modular reduction, covering every recurrence transition.
pub fn lucas_lehmer(e: u32) -> Result<LucasLehmer, String> {
    if !exponent_prime(e) {
        return Err("exponent must be an odd prime in 3..=65536".into());
    }
    let m = (BigUint::one() << e) - BigUint::one();
    let mut s = BigUint::from(4u8);
    let mut hash = Sha256::new();
    state_hash(&mut hash, 0, &s.to_string());
    for i in 1..=e - 2 {
        let square = &s * &s;
        let mut t = if square >= BigUint::from(2u8) {
            square - BigUint::from(2u8)
        } else {
            square + &m - BigUint::from(2u8)
        };
        while t.bits() > u64::from(e) {
            t = (&t & &m) + (t >> e);
        }
        if t >= m {
            t -= &m;
        }
        s = t;
        state_hash(&mut hash, i, &s.to_string());
    }
    let receipt = LucasLehmer {
        exponent: e,
        transitions: e - 2,
        terminal: s.to_string(),
        states_sha256: hex::encode(hash.finalize()),
        theorem: "Keith Conrad, The Lucas-Lehmer Test, Theorem 2.2".into(),
    };
    if !s.is_zero() {
        return Err(format!(
            "Lucas-Lehmer rejects M{e}; terminal {}",
            receipt.terminal
        ));
    }
    Ok(receipt)
}
fn isqrt(n: &BigUint) -> BigUint {
    if n <= &BigUint::one() {
        return n.clone();
    }
    let mut x = BigUint::one() << n.bits().div_ceil(2);
    loop {
        let y = (&x + n / &x) >> 1;
        if y >= x {
            return x;
        }
        x = y;
    }
}
pub fn cornacchia(m: &BigUint) -> Result<Cornacchia, String> {
    let root = (m - BigUint::from(3u8)).modpow(&((m + BigUint::one()) >> 2), m);
    let root = root.clone().min(m - &root);
    if (&root * &root + BigUint::from(3u8)) % m != BigUint::zero() {
        return Err("square root of -3 failed".into());
    }
    let (mut a, mut b) = (m.clone(), root.clone());
    let mut h = Sha256::new();
    let mut divisions = 0;
    while &b * &b >= *m {
        if b.is_zero() {
            return Err("Cornacchia zero divisor".into());
        }
        let q = &a / &b;
        let r = &a % &b;
        h.update(format!("{q}:{r}\n").as_bytes());
        divisions += 1;
        a = b;
        b = r;
    }
    let t = m - &b * &b;
    if &t % BigUint::from(3u8) != BigUint::zero() {
        return Err("Cornacchia remainder not divisible by three".into());
    }
    let y = isqrt(&(&t / BigUint::from(3u8)));
    if &b * &b + BigUint::from(3u8) * &y * &y != *m {
        return Err("Cornacchia square test failed".into());
    }
    Ok(Cornacchia {
        root: root.to_string(),
        x: b.to_string(),
        y: y.to_string(),
        divisions,
        division_sha256: hex::encode(h.finalize()),
    })
}
pub fn norm_statement() -> IdentityStatement {
    let a = Expr::var("a");
    let b = Expr::var("b");
    let t = Expr::var("t");
    let mul = Expr::times;
    let add = Expr::plus;
    let neg = Expr::minus;
    let b2 = mul(b.clone(), b.clone());
    IdentityStatement {
        name: "Eisenstein endomorphism norm before specialization".into(),
        source: "MIT 18.783 Lecture 12 dual identities; independent review TASK-20261016-3b4ca9"
            .into(),
        variables: vec!["a".into(), "b".into(), "t".into()],
        degree_bound: 4,
        lhs: mul(
            add(a.clone(), mul(b.clone(), t.clone())),
            add(
                add(a.clone(), neg(b.clone())),
                neg(mul(b.clone(), t.clone())),
            ),
        ),
        rhs: add(
            add(add(mul(a.clone(), a.clone()), neg(mul(a, b))), b2.clone()),
            neg(mul(
                b2,
                add(add(mul(t.clone(), t.clone()), t), Expr::Const(1)),
            )),
        ),
    }
}
type Polynomial = BTreeMap<[u16; 3], BigInt>;
fn expand(e: &Expr, depth: u32) -> Result<Polynomial, String> {
    if depth > 32 {
        return Err("identity depth limit".into());
    }
    match e {
        Expr::Const(c) => Ok(BTreeMap::from([([0; 3], BigInt::from(*c))])),
        Expr::Var(v) => {
            let i = ["a", "b", "t"]
                .iter()
                .position(|x| *x == v)
                .ok_or("undeclared identity variable")?;
            let mut m = [0; 3];
            m[i] = 1;
            Ok(BTreeMap::from([(m, BigInt::one())]))
        }
        Expr::Neg(x) => Ok(expand(x, depth + 1)?
            .into_iter()
            .map(|(m, c)| (m, -c))
            .collect()),
        Expr::Add(xs) => {
            let mut p = Polynomial::new();
            for x in xs {
                for (m, c) in expand(x, depth + 1)? {
                    *p.entry(m).or_default() += c;
                }
            }
            p.retain(|_, c| !c.is_zero());
            Ok(p)
        }
        Expr::Mul(xs) => {
            let mut p = BTreeMap::from([([0; 3], BigInt::one())]);
            for x in xs {
                let q = expand(x, depth + 1)?;
                let mut r = Polynomial::new();
                for (m, c) in &p {
                    for (n, d) in &q {
                        let k = [m[0] + n[0], m[1] + n[1], m[2] + n[2]];
                        if k.iter().sum::<u16>() > 8 {
                            return Err("identity degree limit".into());
                        }
                        *r.entry(k).or_default() += c * d;
                    }
                }
                r.retain(|_, c| !c.is_zero());
                if r.len() > 256 {
                    return Err("identity size limit".into());
                }
                p = r;
            }
            Ok(p)
        }
    }
}
pub fn exact_identity(statement: &IdentityStatement) -> Result<(), String> {
    if statement.variables != ["a", "b", "t"] {
        return Err("identity variable order".into());
    }
    let normalize = |mut p: Polynomial| {
        p.retain(|_, c| !c.is_zero());
        p
    };
    if normalize(expand(&statement.lhs, 0)?) != normalize(expand(&statement.rhs, 0)?) {
        return Err("integer polynomial coefficients differ".into());
    }
    Ok(())
}
pub fn identity() -> Result<IdentityCertificate, String> {
    let s = norm_statement();
    exact_identity(&s)?;
    idc::issue(
        s,
        20261016,
        32,
        "native issuer; exact integer expansion separately checked; no formal proof pointer",
        None,
    )
    .map_err(|e| e.to_string())
}
fn pnum(s: &str) -> BigUint {
    unsigned(s).unwrap()
}
fn sub(a: &BigUint, b: &BigUint, p: &BigUint) -> BigUint {
    (a + p - b) % p
}
fn on_curve(q: &Affine, p: &BigUint) -> bool {
    q.as_ref().is_none_or(|(x, y)| {
        x < p && y < p && (y * y) % p == ((x * x % p) * x + BigUint::from(7u8)) % p
    })
}
fn add(q: &Affine, r: &Affine, p: &BigUint) -> Affine {
    let (Some((x, y)), Some((u, v))) = (q, r) else {
        return if q.is_none() { r.clone() } else { q.clone() };
    };
    if x == u && (y != v || y.is_zero()) {
        return None;
    }
    let (top, bottom) = if x == u {
        (BigUint::from(3u8) * x * x % p, BigUint::from(2u8) * y % p)
    } else {
        (sub(v, y, p), sub(u, x, p))
    };
    let slope = top * bottom.modpow(&(p - BigUint::from(2u8)), p) % p;
    let xx = sub(&sub(&(&slope * &slope % p), x, p), u, p);
    let yy = sub(&(&slope * sub(x, &xx, p) % p), y, p);
    Some((xx, yy))
}
// Projective Jacobian producer, one field inversion at the final conversion.
fn scalar(q: &Affine, k: &BigInt, p: &BigUint) -> Affine {
    let Some((x, y)) = q else {
        return None;
    };
    let mut y = y.clone();
    if k.sign() == Sign::Minus {
        y = sub(&BigUint::zero(), &y, p);
    }
    let (mut xx, mut yy, mut zz) = (BigUint::one(), BigUint::one(), BigUint::zero());
    for i in (0..k.magnitude().bits()).rev() {
        if !zz.is_zero() {
            if yy.is_zero() {
                zz = BigUint::zero();
            } else {
                let a = &xx * &xx % p;
                let b = &yy * &yy % p;
                let c = &b * &b % p;
                let d = BigUint::from(2u8)
                    * sub(
                        &sub(&((&xx + &b).modpow(&BigUint::from(2u8), p)), &a, p),
                        &c,
                        p,
                    )
                    % p;
                let e = BigUint::from(3u8) * a % p;
                let f = &e * &e % p;
                let nx = sub(&f, &(BigUint::from(2u8) * &d % p), p);
                let ny = sub(
                    &(&e * sub(&d, &nx, p) % p),
                    &(BigUint::from(8u8) * c % p),
                    p,
                );
                zz = BigUint::from(2u8) * &yy * &zz % p;
                xx = nx;
                yy = ny;
            }
        }
        if k.magnitude().bit(i) {
            if zz.is_zero() {
                xx = x.clone();
                yy = y.clone();
                zz = BigUint::one();
                continue;
            }
            let z2 = &zz * &zz % p;
            let u = x * &z2 % p;
            let s = &y * &zz * &z2 % p;
            let h = sub(&u, &xx, p);
            let r = sub(&s, &yy, p);
            if h.is_zero() {
                if r.is_zero() {
                    let affine = Some((x.clone(), y.clone()));
                    let doubled = add(&affine, &affine, p);
                    if let Some((a, b)) = doubled {
                        xx = a;
                        yy = b;
                        zz = BigUint::one();
                    } else {
                        zz = BigUint::zero();
                    }
                } else {
                    zz = BigUint::zero();
                }
                continue;
            }
            let h2 = &h * &h % p;
            let h3 = &h * &h2 % p;
            let v = &xx * &h2 % p;
            let nx = sub(
                &sub(&(&r * &r % p), &h3, p),
                &(BigUint::from(2u8) * &v % p),
                p,
            );
            let ny = sub(&(&r * sub(&v, &nx, p) % p), &(&yy * &h3 % p), p);
            zz = zz * &h % p;
            xx = nx;
            yy = ny;
        }
    }
    if zz.is_zero() {
        return None;
    }
    let zi = zz.modpow(&(p - BigUint::from(2u8)), p);
    let zi2 = &zi * &zi % p;
    Some((xx * &zi2 % p, yy * zi2 * zi % p))
}
fn encode(q: &Affine) -> Point {
    q.as_ref().map(|(x, y)| [x.to_string(), y.to_string()])
}
fn decode(q: &Point) -> Result<Affine, String> {
    q.as_ref()
        .map(|v| Ok((unsigned(&v[0])?, unsigned(&v[1])?)))
        .transpose()
}
pub fn geometric(q: &Point, a: &BigInt, b: &BigInt, beta: &BigUint) -> Result<Point, String> {
    let p = pnum(P);
    let q = decode(q)?;
    if !on_curve(&q, &p) {
        return Err("off-curve input".into());
    }
    let phi = q.as_ref().map(|(x, y)| (x * beta % &p, y.clone()));
    Ok(encode(&add(&scalar(&q, a, &p), &scalar(&phi, b, &p), &p)))
}
fn cube_root(p: &BigUint) -> BigUint {
    let e = (p - BigUint::one()) / BigUint::from(3u8);
    for seed in 2u32..100 {
        let r = BigUint::from(seed).modpow(&e, p);
        if r != BigUint::one() {
            return r.clone().min(&r * &r % p);
        }
    }
    panic!("no nontrivial cube root in fixed search")
}
pub fn construct(e: u32) -> Result<Certificate, String> {
    let primality = lucas_lehmer(e)?;
    let degree = (BigUint::one() << e) - BigUint::one();
    let cornacchia = cornacchia(&degree)?;
    let x = unsigned(&cornacchia.x)?;
    let y = unsigned(&cornacchia.y)?;
    let a = BigInt::from(&x + &y);
    let b = BigInt::from(BigUint::from(2u8) * &y);
    let p = pnum(P);
    let n = pnum(N);
    let beta = cube_root(&p);
    let mut lambda = cube_root(&n);
    let g = decode(&domain().generator)?;
    let phi = g.as_ref().map(|(x, y)| (x * &beta % &p, y.clone()));
    if scalar(&g, &BigInt::from(lambda.clone()), &p) != phi {
        lambda = &lambda * &lambda % &n;
    }
    if scalar(&g, &BigInt::from(lambda.clone()), &p) != phi {
        return Err("beta/lambda pairing failed".into());
    }
    let au = a.magnitude();
    let bu = b.magnitude();
    let differential = (au + bu * &beta) % &p;
    let k = (au + bu * &lambda) % &n;
    if differential.is_zero() || k.is_zero() || degree == p {
        return Err("separability or subgroup invertibility failed".into());
    }
    let mut points = vec![
        None,
        domain().generator.clone(),
        g.as_ref()
            .map(|(x, y)| [x.to_string(), sub(&BigUint::zero(), y, &p).to_string()]),
    ];
    for k in [2, 3, 17, 101, 65537] {
        points.push(encode(&scalar(&g, &BigInt::from(k), &p)));
    }
    let evaluations = points
        .into_iter()
        .map(|input| {
            Ok(Evaluation {
                output: geometric(&input, &a, &b, &beta)?,
                input,
            })
        })
        .collect::<Result<Vec<_>, String>>()?;
    Ok(Certificate {
        schema: "compact-cm-isogeny/v1".into(),
        representation:
            "compact known-endomorphism self-map; unreduced signed integer coefficients".into(),
        domain: domain(),
        exponent: e,
        degree: degree.to_string(),
        a: a.to_string(),
        b: b.to_string(),
        beta: beta.to_string(),
        lambda: lambda.to_string(),
        differential: differential.to_string(),
        subgroup_scalar: k.to_string(),
        subgroup_inverse: k.modpow(&(&n - BigUint::from(2u8)), &n).to_string(),
        primality,
        cornacchia,
        identity: identity()?,
        evaluations,
        tool: crate::provenance(),
    })
}
fn fresh(out: &Path) -> Result<(), String> {
    fs::create_dir(out).map_err(|e| format!("fresh output directory required: {e}"))
}
fn save(path: &Path, v: &Value) -> Result<(), String> {
    fs::write(
        path,
        format!("{}\n", serde_json::to_string_pretty(v).unwrap()),
    )
    .map_err(|e| e.to_string())
}
fn limited_bytes(path: &Path) -> Result<Vec<u8>, String> {
    use std::io::Read;
    let mut bytes = Vec::new();
    fs::File::open(path)
        .map_err(|e| e.to_string())?
        .take(2 * 1024 * 1024 + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| e.to_string())?;
    if bytes.len() > 2 * 1024 * 1024 {
        return Err("certificate exceeds 2 MiB".into());
    }
    Ok(bytes)
}
pub fn read(path: &Path) -> Result<Certificate, String> {
    serde_json::from_slice(&limited_bytes(path)?).map_err(|e| e.to_string())
}
fn execute_unprotected(args: &[String]) -> Result<(Value, i32), String> {
    match args.first().map(String::as_str) {
        Some("construct") => {
            let o=crate::Options::parse(&args[1..], &["exponent","out"])?;
            let exponent=u32::try_from(o.number("exponent",None)?).map_err(|_|"exponent overflow")?;
            let out=Path::new(o.need("out")?);
            fresh(out)?;
            let started=Instant::now();
            let (record,code)=match construct(exponent) {
                Ok(certificate) => {
                    save(&out.join("certificate.json"), &serde_json::to_value(certificate).unwrap())?;
                    (json!({"schema":"compact-cm-construction/v1","status":"CONSTRUCTED","independent_verification":"NOT_RUN","elapsed_ms_operational_only":started.elapsed().as_millis(),"tool":crate::provenance()}),0)
                }
                Err(error) => (json!({"status":"REJECTED","error":error,"elapsed_ms_operational_only":started.elapsed().as_millis()}),1),
            };
            save(&out.join("construction.json"), &record)?;
            Ok((record,code))
        }
        Some("verify") => {
            let o=crate::Options::parse(&args[1..], &["input","out"])?;
            let out=Path::new(o.need("out")?);
            fresh(out)?;
            let started=Instant::now();
            // Hash the same bounded snapshot that is parsed and checked.
            let bytes=limited_bytes(Path::new(o.need("input")?))?;
            let result=serde_json::from_slice::<Certificate>(&bytes).map_err(|e|e.to_string()).and_then(|c|crate::compact_cm_verify::verify(&c));
            let (mut record,code)=match result {
                Ok(record) => (record,0),
                Err(error) => (json!({"status":"REJECTED","error":error,"score":Value::Null}),1),
            };
            record["elapsed_ms_operational_only"]=json!(started.elapsed().as_millis());
            record["input_sha256"]=json!(crate::sha256_hex(&bytes));
            record["tool"]=crate::provenance();
            save(&out.join("verification.json"), &record)?;
            Ok((record,code))
        }
        Some("evaluate") => {
            let o=crate::Options::parse(&args[1..], &["input","x","y","mode"])?;
            let certificate=read(Path::new(o.need("input")?))?;
            crate::compact_cm_verify::verify(&certificate)?;
            let point=Some([o.need("x")?.into(),o.need("y")?.into()]);
            let result=match o.need("mode")? {
                "geometric" => geometric(&point,&signed(&certificate.a)?,&signed(&certificate.b)?,&unsigned(&certificate.beta)?)?,
                "subgroup" => crate::compact_cm_verify::subgroup_evaluate(&certificate,&point)?,
                _ => return Err("mode must be geometric or subgroup".into()),
            };
            Ok((json!({"status":"PASS","evaluation_semantics":o.need("mode")?,"point":result,"geometric_degree":certificate.degree,"degree_proof":"separate unreduced certificate"}),0))
        }
        _ => Err("use isogeny compact-cm construct --exponent E --out FRESH; verify --input certificate.json --out FRESH; evaluate --input certificate.json --x X --y Y --mode geometric|subgroup".into()),
    }
}

/// Public construction and verification run in the existing native resource
/// supervisor. Workers wait for a positive RSS sample before doing arithmetic.
pub fn execute(args: &[String]) -> Result<(Value, i32), String> {
    if args.first().is_some_and(|s| s == "__worker") {
        if args.len() < 4 || args[args.len() - 2] != "--rss-ready" {
            return Err("resource acknowledgement is required".into());
        }
        let ready = Path::new(&args[args.len() - 1]);
        let started = Instant::now();
        while !ready.is_file() {
            if started.elapsed().as_secs() >= 10 {
                return Ok((
                    json!({"status":"CENSORED","reason":"RESOURCE_HANDSHAKE_UNAVAILABLE"}),
                    3,
                ));
            }
            std::thread::sleep(std::time::Duration::from_millis(5));
        }
        let result = std::panic::catch_unwind(|| execute_unprotected(&args[1..args.len() - 2]));
        fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(ready.with_extension("done"))
            .map_err(|e| e.to_string())?;
        return match result {
            Ok(Ok((record, code))) => Ok((
                json!({"status":if code==0 {"PASS"}else{"FAIL"},"record":record}),
                code,
            )),
            Ok(Err(error)) => Ok((json!({"status":"REJECTED","error":error}), 1)),
            Err(_) => Ok((
                json!({"status":"FAIL","error":"native compact worker panicked"}),
                1,
            )),
        };
    }
    let operation = args.first().map(String::as_str).unwrap_or("");
    if operation == "evaluate" {
        return execute_unprotected(args);
    }
    let allowed = match operation {
        "construct" => vec!["exponent", "out", "timeout"],
        "verify" => vec!["input", "out", "timeout"],
        _ => return execute_unprotected(args),
    };
    let o = crate::Options::parse(&args[1..], &allowed)?;
    let maximum = if operation == "construct" { 1800 } else { 3600 };
    let timeout = o.number("timeout", Some(maximum))?;
    if timeout == 0 || timeout > maximum {
        return Err(format!("timeout must be 1..={maximum}"));
    }
    // Validate all required arguments before taking a fresh output directory.
    if operation == "construct" {
        o.number("exponent", None)?;
    } else {
        o.need("input")?;
    }
    let out = Path::new(o.need("out")?);
    fresh(out)?;
    let out = fs::canonicalize(out).map_err(|e| e.to_string())?;
    let artifact = out.join("artifact");
    let ready = out.join("resources-ready");
    let mut command = vec!["compact-cm".into(), "__worker".into(), operation.into()];
    for key in if operation == "construct" {
        vec!["exponent"]
    } else {
        vec!["input"]
    } {
        command.push(format!("--{key}"));
        command.push(if key == "input" {
            fs::canonicalize(o.need(key)?)
                .map_err(|e| e.to_string())?
                .to_string_lossy()
                .into_owned()
        } else {
            o.need(key)?.into()
        });
    }
    command.extend([
        "--out".into(),
        artifact.to_string_lossy().into_owned(),
        "--rss-ready".into(),
        ready.to_string_lossy().into_owned(),
    ]);
    let receipt = crate::supervisor::run(
        &std::env::current_exe().map_err(|e| e.to_string())?,
        &out,
        &command,
        "native",
        timeout,
    );
    let receipt: Value = serde_json::from_str(&receipt.dump()).map_err(|e| e.to_string())?;
    let worker: Value =
        serde_json::from_slice(&fs::read(out.join("native.json")).map_err(|e| e.to_string())?)
            .unwrap_or(Value::Null);
    let code = if receipt["status"] == "PASS" {
        0
    } else if receipt["status"] == "FAIL" {
        1
    } else {
        3
    };
    let mut record = if let Some(r) = worker.get("record") {
        r.clone()
    } else {
        json!({"status":if code==3{"CENSORED"}else{"REJECTED"},"error":worker.get("error").cloned().unwrap_or(receipt["status"].clone()),"score":Value::Null})
    };
    if code == 0 && operation == "construct" {
        let bytes = limited_bytes(&artifact.join("certificate.json"))?;
        fs::write(out.join("certificate.json"), &bytes).map_err(|e| e.to_string())?;
        record["certificate_sha256"] = json!(crate::sha256_hex(&bytes));
    }
    if code != 0 {
        record["status"] = json!(if code == 3 { "CENSORED" } else { "REJECTED" });
        record["score"] = Value::Null;
    }
    record["resources"] = receipt;
    save(
        &out.join(if operation == "construct" {
            "construction.json"
        } else {
            "verification.json"
        }),
        &record,
    )?;
    Ok((record, code))
}

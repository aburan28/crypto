//! Independent replay: repository limb integers, general remainders, and affine
//! fixed-limb field arithmetic. Does not call the producer's LL, Cornacchia,
//! scalar multiplication, point evaluator, or polynomial expansion.
use crate::compact_cm::{Certificate, Point};
use crate::identity_certificate::{self as idc, Expr, IdentityStatement};
use crate::verification_field::{Element, PrimeField};
use isogeny_algos::{bigint::Big, int::Int};
use num_bigint::BigUint;
use serde_json::{json, Value};
use sha2::{Digest, Sha256};

type Fe = Element<4>;
type Q = Option<(Fe, Fe)>;
fn z(s: &str) -> Result<Int, String> {
    let (negative, t) = s.strip_prefix('-').map_or((false, s), |t| (true, t));
    if t.is_empty()
        || t.len() > 20000
        || !t.bytes().all(|c| c.is_ascii_digit())
        || (t.len() > 1 && t.starts_with('0'))
        || (negative && t == "0")
    {
        return Err("invalid canonical signed decimal".into());
    }
    let a = Int::from_big(&Big::from_dec(t));
    Ok(if negative { -a } else { a })
}
fn bu(s: &str) -> Result<BigUint, String> {
    let a = z(s)?;
    if a.is_neg() {
        return Err("negative field coordinate".into());
    }
    BigUint::parse_bytes(s.as_bytes(), 10).ok_or("invalid field coordinate".into())
}
fn require(ok: bool, reason: &str) -> Result<(), String> {
    if ok {
        Ok(())
    } else {
        Err(reason.into())
    }
}
fn prime_exponent(e: u32) -> bool {
    e >= 3
        && e <= 65536
        && e % 2 == 1
        && (2..)
            .take_while(|d| u64::from(*d) * u64::from(*d) <= u64::from(e))
            .all(|d| e % d != 0)
}
fn ll(c: &Certificate) -> Result<Int, String> {
    let e = c.exponent;
    require(
        prime_exponent(e),
        "exponent is not an odd prime in 3..=65536",
    )?;
    let m = &Int::one().shl(e as usize) - &Int::one();
    require(z(&c.degree)? == m, "degree differs from 2^e-1")?;
    require(
        c.primality.exponent == e && c.primality.transitions == e - 2,
        "LL iteration metadata mismatch",
    )?;
    require(
        c.primality.theorem == "Keith Conrad, The Lucas-Lehmer Test, Theorem 2.2",
        "LL theorem mismatch",
    )?;
    let mut s = Int::from(4i64);
    let mut h = Sha256::new();
    for i in 0..=e - 2 {
        let text = s.to_string();
        h.update(i.to_le_bytes());
        h.update((text.len() as u64).to_le_bytes());
        h.update(text.as_bytes());
        if i < e - 2 {
            s = (&(&s * &s) - &Int::from(2i64)).modulo(&m);
        }
    }
    require(
        s.is_zero() && c.primality.terminal == "0",
        "LL terminal is not zero",
    )?;
    require(
        hex::encode(h.finalize()) == c.primality.states_sha256,
        "LL recurrence digest mismatch",
    )?;
    Ok(m)
}
fn cornacchia(c: &Certificate, m: &Int) -> Result<(), String> {
    let root = z(&c.cornacchia.root)?;
    require(
        !root.is_neg() && &root * &Int::from(2i64) <= *m,
        "Cornacchia root not normalized",
    )?;
    require(
        (&(&root * &root) + &Int::from(3i64)).modulo(m).is_zero(),
        "Cornacchia root not square root of -3",
    )?;
    let (mut a, mut b) = (m.clone(), root);
    let mut h = Sha256::new();
    let mut count = 0;
    while &b * &b >= *m {
        require(!b.is_zero(), "Cornacchia zero divisor")?;
        let (q, r) = a.divrem_trunc(&b);
        require(
            !r.is_neg() && r < b && &(&q * &b) + &r == a,
            "Cornacchia Euclidean equation failed",
        )?;
        h.update(format!("{q}:{r}\n").as_bytes());
        a = b;
        b = r;
        count += 1;
    }
    let x = z(&c.cornacchia.x)?;
    let y = z(&c.cornacchia.y)?;
    require(
        x == b && !x.is_neg() && !y.is_neg(),
        "Cornacchia x/y or stopping mismatch",
    )?;
    let residue = m - &(&x * &x);
    let (q, r) = residue.divrem_trunc(&Int::from(3i64));
    require(
        r.is_zero() && !q.is_neg() && q.isqrt() == y && &y * &y == q,
        "Cornacchia exact square test failed",
    )?;
    require(
        &(&x * &x) + &(&Int::from(3i64) * &(&y * &y)) == *m,
        "Cornacchia norm mismatch",
    )?;
    require(
        count == c.cornacchia.divisions
            && hex::encode(h.finalize()) == c.cornacchia.division_sha256,
        "Cornacchia division transcript mismatch",
    )?;
    require(
        z(&c.a)? == &x + &y && z(&c.b)? == &Int::from(2i64) * &y,
        "unreduced canonical coefficient transform failed",
    )?;
    Ok(())
}
// Dense independent integer coefficient expansion, degrees <=4 in a,b,t.
fn polynomial(e: &Expr) -> Result<Vec<Int>, String> {
    let mut out = vec![Int::zero(); 125];
    match e {
        Expr::Const(c) => out[0] = Int::from(*c),
        Expr::Var(v) => {
            out[match v.as_str() {
                "a" => 25,
                "b" => 5,
                "t" => 1,
                _ => return Err("undeclared variable".into()),
            }] = Int::one();
        }
        Expr::Neg(x) => {
            out = polynomial(x)?.into_iter().map(|a| -a).collect();
        }
        Expr::Add(xs) => {
            for x in xs {
                for (i, c) in polynomial(x)?.into_iter().enumerate() {
                    out[i] = &out[i] + &c;
                }
            }
        }
        Expr::Mul(xs) => {
            out[0] = Int::one();
            for x in xs {
                let q = polynomial(x)?;
                let mut r = vec![Int::zero(); 125];
                for (i, a) in out.iter().enumerate().filter(|(_, a)| !a.is_zero()) {
                    for (j, b) in q.iter().enumerate().filter(|(_, b)| !b.is_zero()) {
                        let aa = i / 25 + j / 25;
                        let bb = (i / 5) % 5 + (j / 5) % 5;
                        let tt = i % 5 + j % 5;
                        require(aa + bb + tt <= 4, "exact identity degree limit")?;
                        let k = aa * 25 + bb * 5 + tt;
                        r[k] = &r[k] + &(a * b);
                    }
                }
                out = r;
            }
        }
    }
    Ok(out)
}
pub fn exact_identity(s: &IdentityStatement) -> Result<(), String> {
    require(s.variables == ["a", "b", "t"], "identity variable mismatch")?;
    require(
        polynomial(&s.lhs)? == polynomial(&s.rhs)?,
        "exact Z identity coefficient mismatch",
    )
}
struct Engine {
    f: PrimeField<4>,
}
impl Engine {
    fn new() -> Self {
        Self {
            f: PrimeField::new(&BigUint::parse_bytes(crate::compact_cm::P.as_bytes(), 10).unwrap())
                .unwrap(),
        }
    }
    fn decode(&self, q: &Point) -> Result<Q, String> {
        q.as_ref()
            .map(|v| {
                let x = bu(&v[0])?;
                let y = bu(&v[1])?;
                let p = BigUint::parse_bytes(crate::compact_cm::P.as_bytes(), 10).unwrap();
                require(x < p && y < p, "coordinate out of range")?;
                let q = Some((self.f.from_big(&x), self.f.from_big(&y)));
                require(self.on_curve(&q), "off-curve input")?;
                Ok(q.unwrap())
            })
            .transpose()
    }
    fn encode(&self, q: &Q) -> Point {
        q.as_ref()
            .map(|(x, y)| [self.f.to_big(x).to_string(), self.f.to_big(y).to_string()])
    }
    fn on_curve(&self, q: &Q) -> bool {
        q.as_ref().is_none_or(|(x, y)| {
            self.f.sqr(y)
                == self
                    .f
                    .add(&self.f.mul(&self.f.sqr(x), x), &self.f.from_u64(7))
        })
    }
    fn neg(&self, q: &Q) -> Q {
        q.map(|(x, y)| (x, self.f.neg(&y)))
    }
    fn add(&self, q: &Q, r: &Q) -> Q {
        let (Some((x, y)), Some((u, v))) = (q, r) else {
            return if q.is_none() { *r } else { *q };
        };
        if x == u && (y != v || self.f.is_zero(y)) {
            return None;
        }
        let (top, bottom) = if x == u {
            (
                self.f.mul(&self.f.from_u64(3), &self.f.sqr(x)),
                self.f.mul(&self.f.from_u64(2), y),
            )
        } else {
            (self.f.sub(v, y), self.f.sub(u, x))
        };
        let slope = self.f.div(&top, &bottom).unwrap();
        let xx = self.f.sub(&self.f.sub(&self.f.sqr(&slope), x), u);
        let yy = self.f.sub(&self.f.mul(&slope, &self.f.sub(x, &xx)), y);
        Some((xx, yy))
    }
    fn scalar(&self, q: &Q, k: &Int) -> Q {
        let mut base = if k.is_neg() { self.neg(q) } else { *q };
        let mut acc = None;
        for i in 0..k.bits() {
            if k.mag().bit(i) {
                acc = self.add(&acc, &base);
            }
            base = self.add(&base, &base);
        }
        acc
    }
    fn alpha(&self, q: &Q, a: &Int, b: &Int, beta: &Fe) -> Q {
        let phi = q.map(|(x, y)| (self.f.mul(&x, beta), y));
        self.add(&self.scalar(q, a), &self.scalar(&phi, b))
    }
}
pub fn subgroup_evaluate(c: &Certificate, q: &Point) -> Result<Point, String> {
    let e = Engine::new();
    let q = e.decode(q)?;
    require(
        e.scalar(&q, &z(&c.domain.n)?).is_none(),
        "point outside subgroup",
    )?;
    Ok(e.encode(&e.scalar(&q, &z(&c.subgroup_scalar)?)))
}
pub fn verify(c: &Certificate) -> Result<Value, String> {
    require(c.schema == "compact-cm-isogeny/v1", "wrong schema")?;
    require(
        c.representation
            == "compact known-endomorphism self-map; unreduced signed integer coefficients",
        "wrong representation",
    )?;
    require(
        c.domain == crate::compact_cm::domain(),
        "source differs from trusted exact SEC2 tuple",
    )?;
    let p = z(&c.domain.p)?;
    let n = z(&c.domain.n)?;
    let ell = ll(c)?;
    require(ell != p, "degree equals characteristic")?;
    require(
        c.identity.statement == crate::compact_cm::norm_statement(),
        "unrecognized exact identity statement",
    )?;
    exact_identity(&c.identity.statement)?;
    require(
        c.identity.points == 32 && c.identity.seed == 20261016 && c.identity.proof.is_none(),
        "unexpected evaluated identity parameters",
    )?;
    idc::check_with_extra(&c.identity, Some((0x4d41502d434d, 16)))
        .map_err(|e| e.to_string())
        .and_then(|v| match v {
            idc::Verdict::Accept { .. } => Ok(()),
            _ => Err("evaluated IDC1 rejected".into()),
        })?;
    let a = z(&c.a)?;
    let b = z(&c.b)?;
    let norm = &(&(&a * &a) - &(&a * &b)) + &(&b * &b);
    require(
        norm == ell,
        "geometric norm from unreduced integers differs from degree",
    )?;
    cornacchia(c, &ell)?;
    let beta = z(&c.beta)?;
    let lambda = z(&c.lambda)?;
    require(
        !beta.is_neg()
            && beta < p
            && beta != Int::one()
            && (&(&beta * &beta) + &(&beta + &Int::one()))
                .modulo(&p)
                .is_zero(),
        "beta is not a primitive cube root",
    )?;
    require(
        !lambda.is_neg()
            && lambda < n
            && (&(&lambda * &lambda) + &(&lambda + &Int::one()))
                .modulo(&n)
                .is_zero(),
        "lambda is not a cube-root eigenvalue",
    )?;
    let derivative = (&a + &(&b * &beta)).modulo(&p);
    let derivative_dual = (&a + &(&b * &(&beta * &beta))).modulo(&p);
    require(
        !derivative.is_zero()
            && z(&c.differential)? == derivative
            && (&derivative * &derivative_dual).modulo(&p) == ell.modulo(&p),
        "differential/separability check failed",
    )?;
    let k = (&a + &(&b * &lambda)).modulo(&n);
    let kd = (&a + &(&b * &(&lambda * &lambda))).modulo(&n);
    let inverse = z(&c.subgroup_inverse)?;
    require(
        !k.is_zero()
            && k == z(&c.subgroup_scalar)?
            && !inverse.is_neg()
            && inverse < n
            && (&k * &inverse).modulo(&n) == Int::one()
            && (&k * &kd).modulo(&n) == ell.modulo(&n),
        "subgroup invertibility/product failed",
    )?;
    let e = Engine::new();
    let g = e.decode(&c.domain.generator)?;
    require(
        g.is_some() && e.scalar(&g, &n).is_none(),
        "generator/order check failed",
    )?;
    let width = (&Int::from(4i64) * &p).isqrt();
    let center = &p + &Int::one();
    let lower = &center - &width;
    let upper = &center + &width;
    require(
        lower > Int::zero() && lower <= n && upper >= n && upper < (&n * &Int::from(2i64)),
        "Hasse interval has ambiguous multiples",
    )?;
    let beta_f = e.f.from_big(&bu(&c.beta)?);
    let phi_g = g.map(|(x, y)| (e.f.mul(&x, &beta_f), y));
    require(
        e.scalar(&g, &lambda) == phi_g,
        "full-point beta/lambda pairing failed",
    )?;
    let mut expected = vec![None, g, e.neg(&g)];
    for seed in [2, 3, 17, 101, 65537] {
        expected.push(e.scalar(&g, &Int::from(seed as i64)));
    }
    require(
        c.evaluations.len() == expected.len(),
        "fixed point-set coverage mismatch",
    )?;
    let dual_a = &a - &b;
    let dual_b = -b.clone();
    let mut rows = vec![];
    for (record, q) in c.evaluations.iter().zip(&expected) {
        require(record.input == e.encode(q), "fixed public input changed")?;
        let output = e.alpha(q, &a, &b, &beta_f);
        require(
            e.on_curve(&output) && record.output == e.encode(&output),
            "independent geometric point replay failed",
        )?;
        require(
            e.scalar(q, &k) == output,
            "separate subgroup evaluator mismatch",
        )?;
        let dual = e.alpha(q, &dual_a, &dual_b, &beta_f);
        let ell_q = e.scalar(q, &ell);
        require(
            e.alpha(&output, &dual_a, &dual_b, &beta_f) == ell_q
                && e.alpha(&dual, &a, &b, &beta_f) == ell_q,
            "dual composition failed",
        )?;
        require(
            e.alpha(&e.add(q, &g), &a, &b, &beta_f)
                == e.add(&output, &e.alpha(&g, &a, &b, &beta_f)),
            "group-law compatibility failed",
        )?;
        rows.push(json!({"input":record.input,"output":record.output,"independent_replay":"PASS","dual_compositions":"PASS","subgroup_action":"PASS","group_law":"PASS"}));
    }
    Ok(
        json!({"schema":"compact-cm-verification/v1","status":"VERIFIED_COMPACT","representation":"compact known-endomorphism self-map","source_icv1":c.domain.source_icv1,"geometric_prime_degree":c.degree,"degree_bits":c.exponent,"scientific_promotion":false,"expanded_map_challenge_score":Value::Null,"checks":{"domain":"PASS_TRUSTED_SEC2_PARAMETERS","p_n_primality":"TRUSTED_STANDARD_NOT_SELF_CONTAINED_CERTIFICATES","generator_and_hasse_order":"PASS","exponent_primality":"PASS","lucas_lehmer":"PASS","cornacchia":"PASS","exact_integer_identity":"PASS","evaluated_identity_supplement":"PASS","unreduced_norm":"PASS","separability":"PASS","beta_lambda_pairing":"PASS","subgroup_invertibility":"PASS","independent_point_replay":"PASS"},"identity_id":c.identity.id,"differential":c.differential,"public_point_checks":rows}),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compact_cm as producer;
    #[test]
    fn exact_source_tuple_matches_the_registered_standard() {
        let row = crate::catalog::resolve("secp256k1").unwrap();
        let domain = producer::domain();
        assert_eq!(row["slug"], domain.source_icv1);
        for (key, value) in [("p", domain.p.as_str()), ("a", "0"), ("b", "7")] {
            assert_eq!(
                crate::catalog::integer(&row["params"][key]),
                z(value).unwrap()
            );
        }
        assert_eq!(
            crate::catalog::integer(&row["order"]),
            z(&domain.n).unwrap()
        );
        let representation = crate::catalog::representation(row).unwrap();
        for index in 0..2 {
            assert_eq!(
                crate::catalog::integer(&representation["curve"]["generator"][index]),
                z(&domain.generator.as_ref().unwrap()[index]).unwrap()
            );
        }
    }
    #[test]
    fn small_positive_rungs_and_distinct_composite_controls() {
        assert!(producer::exponent_prime(11));
        assert!(producer::lucas_lehmer(11)
            .unwrap_err()
            .contains("Lucas-Lehmer rejects"));
        assert!(!producer::exponent_prime(9));
        for (exponent, degree) in [(3, "7"), (5, "31"), (7, "127")] {
            let c = producer::construct(exponent).unwrap();
            assert_eq!(verify(&c).unwrap()["geometric_prime_degree"], degree);
        }
    }
    #[test]
    fn signed_zero_and_negative_geometric_expressions_replay_independently() {
        let c = producer::construct(3).unwrap();
        let e = Engine::new();
        let beta = e.f.from_big(&bu(&c.beta).unwrap());
        let g = e.decode(&c.domain.generator).unwrap();
        for (a, b) in [(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1), (2, -3)] {
            let actual = producer::geometric(
                &c.domain.generator,
                &a.into(),
                &b.into(),
                &bu(&c.beta).unwrap(),
            )
            .unwrap();
            assert_eq!(
                actual,
                e.encode(&e.alpha(&g, &Int::from(a as i64), &Int::from(b as i64), &beta))
            );
        }
        assert!(producer::geometric(
            &Some(["0".into(), "0".into()]),
            &1.into(),
            &0.into(),
            &bu(&c.beta).unwrap()
        )
        .is_err());
    }
    #[test]
    fn exact_integer_checker_catches_modulus_sized_error_that_evaluated_idc_accepts() {
        let mut statement = producer::norm_statement();
        statement.rhs = Expr::plus(statement.rhs, Expr::Const(idc::MODULUS as i64));
        let cert = idc::issue(
            statement.clone(),
            19,
            32,
            "deliberate modular alias control",
            None,
        )
        .unwrap();
        assert!(matches!(
            idc::check(&cert).unwrap(),
            idc::Verdict::Accept { .. }
        ));
        assert!(producer::exact_identity(&statement).is_err());
        assert!(exact_identity(&statement).is_err());
        let mut sign = producer::norm_statement();
        sign.rhs = Expr::minus(sign.rhs);
        assert!(exact_identity(&sign).is_err());
    }
    #[test]
    fn adding_subgroup_order_preserves_points_and_changes_geometric_degree() {
        let c = producer::construct(3).unwrap();
        let a = z(&c.a).unwrap();
        let b = z(&c.b).unwrap();
        let n = z(&c.domain.n).unwrap();
        let alias = &a + &n;
        let original = &(&(&a * &a) - &(&a * &b)) + &(&b * &b);
        let changed = &(&(&alias * &alias) - &(&alias * &b)) + &(&b * &b);
        assert_eq!(
            &changed - &original,
            &n * &(&(&Int::from(2i64) * &z(&c.cornacchia.x).unwrap()) + &n)
        );
        assert!(changed > original);
        for row in &c.evaluations {
            assert_eq!(
                producer::geometric(
                    &row.input,
                    &producer::signed(&alias.to_string()).unwrap(),
                    &producer::signed(&c.b).unwrap(),
                    &bu(&c.beta).unwrap()
                )
                .unwrap(),
                row.output
            );
        }
        let mut bad = c.clone();
        bad.a = alias.to_string();
        assert!(verify(&bad).is_err());
    }
    #[test]
    fn complete_control_and_other_primitive_root_are_accepted_but_wrong_pairing_is_rejected() {
        let mut c = producer::construct(3).unwrap();
        assert_eq!(verify(&c).unwrap()["status"], "VERIFIED_COMPACT");
        let n = z(&c.domain.n).unwrap();
        let p = z(&c.domain.p).unwrap();
        let lambda = z(&c.lambda).unwrap();
        let other_lambda = (&lambda * &lambda).modulo(&n);
        let a = z(&c.a).unwrap();
        let b = z(&c.b).unwrap();
        let k = (&a + &(&b * &other_lambda)).modulo(&n);
        c.lambda = other_lambda.to_string();
        c.subgroup_scalar = k.to_string();
        c.subgroup_inverse = k.inv_mod(&n).unwrap().to_string();
        assert!(verify(&c).unwrap_err().contains("pairing"));
        let beta = z(&c.beta).unwrap();
        let other_beta = (&beta * &beta).modulo(&p);
        c.beta = other_beta.to_string();
        c.differential = (&a + &(&b * &other_beta)).modulo(&p).to_string();
        for row in &mut c.evaluations {
            row.output = producer::geometric(
                &row.input,
                &producer::signed(&c.a).unwrap(),
                &producer::signed(&c.b).unwrap(),
                &bu(&c.beta).unwrap(),
            )
            .unwrap();
        }
        assert_eq!(verify(&c).unwrap()["status"], "VERIFIED_COMPACT");
    }
}

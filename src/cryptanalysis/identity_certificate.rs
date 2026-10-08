//! # Identity certificates: machine-checkable records for algebraic identity claims.
//!
//! ## Why this exists
//!
//! Alman and Vassilevska Williams (`arXiv:2610.06783v1`) refuted the 3SUM and
//! APSP hypotheses with an algorithm whose base identity was found by machine
//! and whose main results were verified in Lean 4. This repository already
//! verifies *runs* (sealed `ecbench` sessions, replay receipts) but had no
//! record for *identity* claims: "polynomial `L` equals polynomial `R`".
//! Such a claim usually lived in prose, sometimes beside a unit test whose
//! evaluation points the claim's author chose. A checker that trusts the
//! author's points is a checker the author can satisfy by accident or by
//! construction. The next machine-found identity to arrive here should be
//! accepted only through a record a stranger can replay without trusting
//! anyone.
//!
//! ## What a certificate is
//!
//! An [`IdentityCertificate`] is a sealed record (`identity.certificate/v1`,
//! id `IDC1h<16 hex>`) containing an [`IdentityStatement`] (two expression
//! trees over named variables, a declared total-degree bound, a source), the
//! modulus `p = 2^61 - 1`, the number of evaluation points, a seed, and the
//! SHA-256 digest of the evaluated values. The checker ([`check`]):
//!
//! 1. recomputes the id from the canonical statement bytes and refuses a
//!    record whose id does not match (the statement cannot be swapped under a
//!    sealed id);
//! 2. derives every evaluation point deterministically from the **id and the
//!    seed together**, so the author cannot choose points independently of
//!    the statement they are certifying;
//! 3. evaluates both sides at every point and rejects on the first mismatch,
//!    returning the point as a counterexample;
//! 4. recomputes the digest of the left-hand values and rejects on mismatch;
//! 5. optionally evaluates at extra points of the *checker's* choosing
//!    ([`check_with_extra`]), which the certificate never saw.
//!
//! The soundness argument is Schwartz-Zippel: a nonzero polynomial of total
//! degree at most `d` over `F_p` vanishes at a uniformly random point with
//! probability at most `d / p`, so `k` independent points accept a false
//! identity with probability at most `(d / p)^k`. [`false_accept_log2`]
//! reports `k * (log2 d - 61)`; for the ten-multiplication identity below
//! (`d = 3`, `k = 32`) that is below `2^{-1900}`.
//!
//! ## What a certificate does not prove
//!
//! It certifies the identity **over `F_p`**. An identity over `Z` whose two
//! sides differ by a polynomial every coefficient of which is divisible by
//! `p = 2^61 - 1` would pass. The degree bound is syntactically checked, not
//! trusted: [`syntactic_degree`] must not exceed the declared bound. A
//! certificate with `proof: Some(ProofRef)` additionally names a formal
//! proof artifact by path and SHA-256; this module records the reference and
//! does not run the prover, so the `proved` tier is a pointer the reader must
//! follow, exactly as the repository's citation rules require.
//!
//! ## The worked example
//!
//! [`schoenhage_identity`] is Lemma 6 of the paper: Schoenhage's ten-term
//! identity computing a `3x3` outer product and a `2x2` inner product
//! together, plus a harmless error polynomial. It is the load-bearing
//! identity under every speedup in that paper, which makes it the right
//! first thing to certify here.

use std::collections::BTreeMap;
use std::fmt;

use serde::{Deserialize, Serialize};

use crate::hash::sha256::sha256;

/// Schema tag written into every certificate.
pub const SCHEMA: &str = "identity.certificate/v1";

/// The evaluation field: the Mersenne prime `2^61 - 1`.
pub const MODULUS: u64 = (1u64 << 61) - 1;

/// Polynomial expression over named variables with integer constants.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "op", content = "args", rename_all = "snake_case")]
pub enum Expr {
    /// A named variable; must appear in the statement's variable list.
    Var(String),
    /// An integer constant.
    Const(i64),
    /// Additive inverse.
    Neg(Box<Expr>),
    /// Sum of the arguments (empty sum is zero).
    Add(Vec<Expr>),
    /// Product of the arguments (empty product is one).
    Mul(Vec<Expr>),
}

impl Expr {
    /// Variable reference.
    pub fn var(name: &str) -> Expr {
        Expr::Var(name.to_string())
    }
    /// Sum of two expressions.
    pub fn plus(a: Expr, b: Expr) -> Expr {
        Expr::Add(vec![a, b])
    }
    /// Product of two expressions.
    pub fn times(a: Expr, b: Expr) -> Expr {
        Expr::Mul(vec![a, b])
    }
    /// Sum of many expressions.
    pub fn sum(terms: Vec<Expr>) -> Expr {
        Expr::Add(terms)
    }
    /// Product of many expressions.
    pub fn product(terms: Vec<Expr>) -> Expr {
        Expr::Mul(terms)
    }
    /// Additive inverse.
    pub fn minus(a: Expr) -> Expr {
        Expr::Neg(Box::new(a))
    }
}

/// The claim: `lhs == rhs` as polynomials in `variables`.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct IdentityStatement {
    /// Human-readable name.
    pub name: String,
    /// Where the identity comes from (paper, lemma, or "machine search run X").
    pub source: String,
    /// The variables, in a fixed order; evaluation points are assigned in this order.
    pub variables: Vec<String>,
    /// Declared upper bound on the total degree of `lhs - rhs`.
    pub degree_bound: u32,
    /// Left-hand side.
    pub lhs: Expr,
    /// Right-hand side.
    pub rhs: Expr,
}

/// A pointer to a formal proof artifact. Recorded, never executed here.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ProofRef {
    /// Prover kind, e.g. `lean4`.
    pub kind: String,
    /// Repository-relative path or durable location of the proof source.
    pub path: String,
    /// SHA-256 (hex) of that artifact.
    pub sha256: String,
}

/// Sealed certificate record.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct IdentityCertificate {
    /// Always [`SCHEMA`].
    pub schema: String,
    /// `IDC1h` + first 16 hex of SHA-256 over the canonical statement bytes.
    pub id: String,
    /// The claim.
    pub statement: IdentityStatement,
    /// Always [`MODULUS`] in v1.
    pub modulus: u64,
    /// Number of evaluation points.
    pub points: u32,
    /// Author-supplied seed, mixed with the id before any point is drawn.
    pub seed: u64,
    /// SHA-256 (hex) over the little-endian `u64` left-hand values, in point order.
    pub value_digest: String,
    /// Who issued and who independently checked; free text, part of the record.
    pub independence: String,
    /// Optional formal proof pointer (the `proved` tier).
    pub proof: Option<ProofRef>,
}

/// Why a statement or certificate was refused before evaluation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Malformed {
    /// `schema` field is not [`SCHEMA`].
    Schema(String),
    /// `modulus` is not [`MODULUS`].
    Modulus(u64),
    /// `points == 0`.
    NoPoints,
    /// The expression names a variable not in `variables`.
    UnknownVariable(String),
    /// A variable name is listed twice.
    DuplicateVariable(String),
    /// Syntactic degree exceeds the declared bound.
    DegreeBound { declared: u32, syntactic: u32 },
    /// Recomputed id differs from the recorded one.
    IdMismatch {
        recorded: String,
        recomputed: String,
    },
    /// `value_digest` is not 64 hex characters.
    DigestShape,
}

impl fmt::Display for Malformed {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Malformed::Schema(s) => write!(f, "schema {s:?} is not {SCHEMA}"),
            Malformed::Modulus(m) => write!(f, "modulus {m} is not {MODULUS}"),
            Malformed::NoPoints => write!(f, "points must be positive"),
            Malformed::UnknownVariable(v) => write!(f, "unknown variable {v}"),
            Malformed::DuplicateVariable(v) => write!(f, "duplicate variable {v}"),
            Malformed::DegreeBound {
                declared,
                syntactic,
            } => write!(
                f,
                "declared degree bound {declared} below syntactic degree {syntactic}"
            ),
            Malformed::IdMismatch {
                recorded,
                recomputed,
            } => write!(f, "id {recorded} does not match statement ({recomputed})"),
            Malformed::DigestShape => write!(f, "value_digest is not 64 hex chars"),
        }
    }
}

impl std::error::Error for Malformed {}

/// A point at which the two sides disagree.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Counterexample {
    /// Index of the point in the derived sequence (or in the extra sequence).
    pub index: u32,
    /// Variable assignments, in statement order.
    pub assignment: Vec<u64>,
    /// Left-hand value.
    pub lhs: u64,
    /// Right-hand value.
    pub rhs: u64,
}

/// Outcome of [`check`].
#[derive(Clone, Debug, PartialEq)]
pub enum Verdict {
    /// Every derived point (and every extra point) agreed and the digest matched.
    Accept {
        /// Points evaluated from the certificate.
        points: u32,
        /// Extra checker-chosen points evaluated.
        extra_points: u32,
        /// `log2` of the Schwartz-Zippel false-accept bound over all points.
        false_accept_log2: f64,
    },
    /// The sides disagree at a point.
    Reject(Counterexample),
    /// All points agreed but the recorded digest does not match the values.
    DigestMismatch {
        /// Digest recorded in the certificate.
        recorded: String,
        /// Digest recomputed by the checker.
        recomputed: String,
    },
}

// ---------------------------------------------------------------------------
// Field arithmetic and point derivation
// ---------------------------------------------------------------------------

#[inline]
fn fadd(a: u64, b: u64) -> u64 {
    let s = a + b; // both < 2^61, no overflow
    if s >= MODULUS {
        s - MODULUS
    } else {
        s
    }
}

#[inline]
fn fneg(a: u64) -> u64 {
    if a == 0 {
        0
    } else {
        MODULUS - a
    }
}

#[inline]
fn fmul(a: u64, b: u64) -> u64 {
    let prod = (a as u128) * (b as u128);
    // 2^61 == 1 mod p, so fold the high part down once, then normalise.
    let lo = (prod & (MODULUS as u128)) as u64;
    let hi = (prod >> 61) as u64;
    let folded = fadd(lo, hi % MODULUS);
    folded % MODULUS
}

fn fconst(c: i64) -> u64 {
    let m = MODULUS as i128;
    let r = ((c as i128) % m + m) % m;
    r as u64
}

/// splitmix64 step, used only to expand a 256-bit hash into field elements.
fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Field element from a 64-bit draw. Rejection would be exact; folding is
/// uniform to within `2^-61` and keeps derivation branch-free.
fn draw_field(state: &mut u64) -> u64 {
    splitmix64(state) % MODULUS
}

/// Derive the point stream state from `(id, seed, domain)`.
fn stream_state(id: &str, seed: u64, domain: &[u8]) -> u64 {
    let mut buf = Vec::with_capacity(id.len() + 8 + domain.len());
    buf.extend_from_slice(id.as_bytes());
    buf.extend_from_slice(&seed.to_le_bytes());
    buf.extend_from_slice(domain);
    let h = sha256(&buf);
    u64::from_le_bytes(h[..8].try_into().expect("8 bytes"))
}

// ---------------------------------------------------------------------------
// Statement validation, canonical bytes, id
// ---------------------------------------------------------------------------

/// Syntactic total degree: `Var = 1`, `Const = 0`, `Neg = inner`, `Add = max`, `Mul = sum`.
pub fn syntactic_degree(e: &Expr) -> u32 {
    match e {
        Expr::Var(_) => 1,
        Expr::Const(_) => 0,
        Expr::Neg(a) => syntactic_degree(a),
        Expr::Add(ts) => ts.iter().map(syntactic_degree).max().unwrap_or(0),
        Expr::Mul(ts) => ts.iter().map(syntactic_degree).sum(),
    }
}

fn collect_vars<'a>(e: &'a Expr, out: &mut Vec<&'a str>) {
    match e {
        Expr::Var(v) => out.push(v),
        Expr::Const(_) => {}
        Expr::Neg(a) => collect_vars(a, out),
        Expr::Add(ts) | Expr::Mul(ts) => ts.iter().for_each(|t| collect_vars(t, out)),
    }
}

/// Validate a statement's variables and degree bound.
pub fn validate_statement(s: &IdentityStatement) -> Result<(), Malformed> {
    let mut seen = BTreeMap::new();
    for v in &s.variables {
        if seen.insert(v.as_str(), ()).is_some() {
            return Err(Malformed::DuplicateVariable(v.clone()));
        }
    }
    let mut used = Vec::new();
    collect_vars(&s.lhs, &mut used);
    collect_vars(&s.rhs, &mut used);
    for v in used {
        if !seen.contains_key(v) {
            return Err(Malformed::UnknownVariable(v.to_string()));
        }
    }
    let syntactic = syntactic_degree(&s.lhs).max(syntactic_degree(&s.rhs));
    if syntactic > s.degree_bound {
        return Err(Malformed::DegreeBound {
            declared: s.degree_bound,
            syntactic,
        });
    }
    Ok(())
}

/// Canonical bytes: `serde_json::Value` serialises maps with sorted keys and
/// no whitespace, which makes the id independent of field order and layout.
pub fn canonical_statement_bytes(s: &IdentityStatement) -> Vec<u8> {
    let value = serde_json::to_value(s).expect("statement serialises");
    serde_json::to_vec(&value).expect("value serialises")
}

/// `IDC1h` + first 16 hex characters of SHA-256 over the canonical bytes.
pub fn statement_id(s: &IdentityStatement) -> String {
    let h = sha256(&canonical_statement_bytes(s));
    format!("IDC1h{}", hex::encode(&h[..8]))
}

// ---------------------------------------------------------------------------
// Evaluation
// ---------------------------------------------------------------------------

fn eval(e: &Expr, env: &BTreeMap<&str, u64>) -> u64 {
    match e {
        Expr::Var(v) => *env.get(v.as_str()).expect("validated variable"),
        Expr::Const(c) => fconst(*c),
        Expr::Neg(a) => fneg(eval(a, env)),
        Expr::Add(ts) => ts.iter().fold(0u64, |acc, t| fadd(acc, eval(t, env))),
        Expr::Mul(ts) => ts.iter().fold(1u64, |acc, t| fmul(acc, eval(t, env))),
    }
}

fn evaluate_points(
    s: &IdentityStatement,
    id: &str,
    seed: u64,
    domain: &[u8],
    points: u32,
) -> Result<Vec<u64>, Counterexample> {
    let mut state = stream_state(id, seed, domain);
    let mut lhs_values = Vec::with_capacity(points as usize);
    for index in 0..points {
        let assignment: Vec<u64> = s.variables.iter().map(|_| draw_field(&mut state)).collect();
        let env: BTreeMap<&str, u64> = s
            .variables
            .iter()
            .map(String::as_str)
            .zip(assignment.iter().copied())
            .collect();
        let l = eval(&s.lhs, &env);
        let r = eval(&s.rhs, &env);
        if l != r {
            return Err(Counterexample {
                index,
                assignment,
                lhs: l,
                rhs: r,
            });
        }
        lhs_values.push(l);
    }
    Ok(lhs_values)
}

fn digest_values(values: &[u64]) -> String {
    let mut buf = Vec::with_capacity(values.len() * 8);
    for v in values {
        buf.extend_from_slice(&v.to_le_bytes());
    }
    hex::encode(sha256(&buf))
}

/// `log2` of the Schwartz-Zippel false-accept bound `(d / p)^k`.
pub fn false_accept_log2(degree_bound: u32, points: u32) -> f64 {
    if degree_bound == 0 {
        return f64::NEG_INFINITY;
    }
    (points as f64) * ((degree_bound as f64).log2() - 61.0)
}

// ---------------------------------------------------------------------------
// Issue and check
// ---------------------------------------------------------------------------

/// Issuance failure.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum IssueError {
    /// The statement is malformed.
    Malformed(Malformed),
    /// The identity is false at a derived point; no certificate is issued.
    False(Counterexample),
}

impl fmt::Display for IssueError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            IssueError::Malformed(m) => write!(f, "malformed: {m}"),
            IssueError::False(c) => write!(
                f,
                "identity false at point {} (lhs {} != rhs {})",
                c.index, c.lhs, c.rhs
            ),
        }
    }
}

impl std::error::Error for IssueError {}

/// Issue a certificate. Refuses if the statement is malformed or the identity
/// fails at any derived point — an issuer cannot seal a false identity.
pub fn issue(
    statement: IdentityStatement,
    seed: u64,
    points: u32,
    independence: &str,
    proof: Option<ProofRef>,
) -> Result<IdentityCertificate, IssueError> {
    if points == 0 {
        return Err(IssueError::Malformed(Malformed::NoPoints));
    }
    validate_statement(&statement).map_err(IssueError::Malformed)?;
    let id = statement_id(&statement);
    let values =
        evaluate_points(&statement, &id, seed, b"cert", points).map_err(IssueError::False)?;
    Ok(IdentityCertificate {
        schema: SCHEMA.to_string(),
        id,
        statement,
        modulus: MODULUS,
        points,
        seed,
        value_digest: digest_values(&values),
        independence: independence.to_string(),
        proof,
    })
}

/// Structural checks that precede any evaluation.
pub fn validate_certificate(c: &IdentityCertificate) -> Result<(), Malformed> {
    if c.schema != SCHEMA {
        return Err(Malformed::Schema(c.schema.clone()));
    }
    if c.modulus != MODULUS {
        return Err(Malformed::Modulus(c.modulus));
    }
    if c.points == 0 {
        return Err(Malformed::NoPoints);
    }
    if c.value_digest.len() != 64 || !c.value_digest.bytes().all(|b| b.is_ascii_hexdigit()) {
        return Err(Malformed::DigestShape);
    }
    validate_statement(&c.statement)?;
    let recomputed = statement_id(&c.statement);
    if recomputed != c.id {
        return Err(Malformed::IdMismatch {
            recorded: c.id.clone(),
            recomputed,
        });
    }
    Ok(())
}

/// Check a certificate by re-deriving every point and recomputing the digest.
pub fn check(c: &IdentityCertificate) -> Result<Verdict, Malformed> {
    check_with_extra(c, None)
}

/// [`check`], plus `extra` = `(checker_seed, count)` points the certificate
/// never saw. The extra stream is domain-separated from the certificate's.
pub fn check_with_extra(
    c: &IdentityCertificate,
    extra: Option<(u64, u32)>,
) -> Result<Verdict, Malformed> {
    validate_certificate(c)?;
    let values = match evaluate_points(&c.statement, &c.id, c.seed, b"cert", c.points) {
        Ok(v) => v,
        Err(cx) => return Ok(Verdict::Reject(cx)),
    };
    let recomputed = digest_values(&values);
    if recomputed != c.value_digest {
        return Ok(Verdict::DigestMismatch {
            recorded: c.value_digest.clone(),
            recomputed,
        });
    }
    let mut extra_points = 0;
    if let Some((checker_seed, count)) = extra {
        if let Err(cx) = evaluate_points(&c.statement, &c.id, checker_seed, b"extra", count) {
            return Ok(Verdict::Reject(cx));
        }
        extra_points = count;
    }
    Ok(Verdict::Accept {
        points: c.points,
        extra_points,
        false_accept_log2: false_accept_log2(c.statement.degree_bound, c.points + extra_points),
    })
}

// ---------------------------------------------------------------------------
// Worked example: Schoenhage's ten-multiplication identity (paper Lemma 6)
// ---------------------------------------------------------------------------

/// Lemma 6 of `arXiv:2610.06783v1`: with `phat`/`qhat` the `3x3` extensions
/// of the `2x2` inner-product matrices whose columns (rows) sum to zero,
///
/// ```text
/// sum_{i,j} (x_i + phat_ij)(y_j + qhat_ij)(z_ij + z_0)
///   - (x_1 + x_2 + x_3)(y_1 + y_2 + y_3) z_0
/// = G + E,
/// G = sum_{i,j} x_i y_j z_ij + (sum_{i,j<=2} p_ij q_ij) z_0,
/// E = sum_{i,j} (x_i qhat_ij + phat_ij y_j + phat_ij qhat_ij) z_ij.
/// ```
///
/// Twenty-four variables, total degree 3.
pub fn schoenhage_identity() -> IdentityStatement {
    use Expr::{Const, Var};
    let v = |s: &str| Var(s.to_string());
    let x = |i: usize| v(&format!("x{i}"));
    let y = |j: usize| v(&format!("y{j}"));
    let z = |i: usize, j: usize| v(&format!("z{i}{j}"));
    let p = |i: usize, j: usize| v(&format!("p{i}{j}"));
    let q = |i: usize, j: usize| v(&format!("q{i}{j}"));

    // phat: columns sum to zero; third column is zero.
    let phat = |i: usize, j: usize| -> Expr {
        match (i, j) {
            (1 | 2, 1 | 2) => p(i, j),
            (3, 1 | 2) => Expr::minus(Expr::plus(p(1, j), p(2, j))),
            _ => Const(0),
        }
    };
    // qhat: rows sum to zero; third row is zero.
    let qhat = |i: usize, j: usize| -> Expr {
        match (i, j) {
            (1 | 2, 1 | 2) => q(i, j),
            (1 | 2, 3) => Expr::minus(Expr::plus(q(i, 1), q(i, 2))),
            _ => Const(0),
        }
    };

    let mut lhs_terms = Vec::new();
    for i in 1..=3 {
        for j in 1..=3 {
            lhs_terms.push(Expr::product(vec![
                Expr::plus(x(i), phat(i, j)),
                Expr::plus(y(j), qhat(i, j)),
                Expr::plus(z(i, j), v("z0")),
            ]));
        }
    }
    lhs_terms.push(Expr::product(vec![
        Expr::minus(Expr::sum(vec![x(1), x(2), x(3)])),
        Expr::sum(vec![y(1), y(2), y(3)]),
        v("z0"),
    ]));
    let lhs = Expr::sum(lhs_terms);

    let mut g_terms = Vec::new();
    for i in 1..=3 {
        for j in 1..=3 {
            g_terms.push(Expr::product(vec![x(i), y(j), z(i, j)]));
        }
    }
    let inner = Expr::sum(
        (1..=2)
            .flat_map(|i| (1..=2).map(move |j| (i, j)))
            .map(|(i, j)| Expr::times(p(i, j), q(i, j)))
            .collect(),
    );
    g_terms.push(Expr::times(inner, v("z0")));

    let mut e_terms = Vec::new();
    for i in 1..=3 {
        for j in 1..=3 {
            e_terms.push(Expr::times(
                Expr::sum(vec![
                    Expr::times(x(i), qhat(i, j)),
                    Expr::times(phat(i, j), y(j)),
                    Expr::times(phat(i, j), qhat(i, j)),
                ]),
                z(i, j),
            ));
        }
    }
    let rhs = Expr::plus(Expr::sum(g_terms), Expr::sum(e_terms));

    let mut variables: Vec<String> = Vec::new();
    for i in 1..=3 {
        variables.push(format!("x{i}"));
    }
    for j in 1..=3 {
        variables.push(format!("y{j}"));
    }
    for (i, j) in [(1, 1), (1, 2), (2, 1), (2, 2)] {
        variables.push(format!("p{i}{j}"));
    }
    for (i, j) in [(1, 1), (1, 2), (2, 1), (2, 2)] {
        variables.push(format!("q{i}{j}"));
    }
    for i in 1..=3 {
        for j in 1..=3 {
            variables.push(format!("z{i}{j}"));
        }
    }
    variables.push("z0".to_string());

    IdentityStatement {
        name: "Schoenhage ten-multiplication identity (3x3 outer + 2x2 inner, with error E)"
            .to_string(),
        source: "Alman-Vassilevska Williams, arXiv:2610.06783v1, Lemma 6 (Schoenhage 1981, Sec. 6, eq. 6.3, k = n = 3)".to_string(),
        variables,
        degree_bound: 3,
        lhs,
        rhs,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const INDEP: &str =
        "issued by the module's own constructor; checked by replaying from the sealed id and seed";

    #[test]
    fn schoenhage_identity_certifies_and_checks() {
        let cert = issue(schoenhage_identity(), 7, 32, INDEP, None).expect("issues");
        assert!(cert.id.starts_with("IDC1h") && cert.id.len() == 5 + 16);
        match check(&cert).expect("well-formed") {
            Verdict::Accept {
                points,
                extra_points,
                false_accept_log2,
            } => {
                assert_eq!(points, 32);
                assert_eq!(extra_points, 0);
                assert!(false_accept_log2 < -1800.0);
            }
            other => panic!("expected accept, got {other:?}"),
        }
        match check_with_extra(&cert, Some((0xDEAD_BEEF, 16))).expect("well-formed") {
            Verdict::Accept { extra_points, .. } => assert_eq!(extra_points, 16),
            other => panic!("expected accept, got {other:?}"),
        }
    }

    #[test]
    fn the_identity_without_its_error_term_is_false_and_cannot_be_issued() {
        // Drop E: claim lhs == G alone. Lemma 6 says this is false.
        let mut s = schoenhage_identity();
        if let Expr::Add(ref mut parts) = s.rhs {
            parts.pop();
        } else {
            panic!("rhs shape");
        }
        s.name = "Schoenhage identity with E dropped (false)".to_string();
        match issue(s.clone(), 7, 32, INDEP, None) {
            Err(IssueError::False(cx)) => {
                assert_eq!(
                    cx.index, 0,
                    "a degree-3 false identity fails at the first point"
                );
                assert_ne!(cx.lhs, cx.rhs);
                assert_eq!(cx.assignment.len(), 24);
            }
            other => panic!("expected refusal, got {other:?}"),
        }
        // A forged certificate for the false statement is rejected on check.
        let forged = IdentityCertificate {
            schema: SCHEMA.to_string(),
            id: statement_id(&s),
            statement: s,
            modulus: MODULUS,
            points: 32,
            seed: 7,
            value_digest: "00".repeat(32),
            independence: "forged in test".to_string(),
            proof: None,
        };
        assert!(matches!(check(&forged), Ok(Verdict::Reject(_))));
    }

    #[test]
    fn swapping_the_statement_under_a_sealed_id_is_malformed() {
        let mut cert = issue(schoenhage_identity(), 7, 8, INDEP, None).expect("issues");
        cert.statement.name.push_str(" (edited)");
        assert!(matches!(check(&cert), Err(Malformed::IdMismatch { .. })));
    }

    #[test]
    fn tampering_with_the_digest_is_detected() {
        let mut cert = issue(schoenhage_identity(), 7, 8, INDEP, None).expect("issues");
        cert.value_digest = "ab".repeat(32);
        assert!(matches!(check(&cert), Ok(Verdict::DigestMismatch { .. })));
    }

    #[test]
    fn points_depend_on_the_id_so_the_author_cannot_pre_pick_them() {
        // Same seed, different statements: different point streams.
        let a = issue(schoenhage_identity(), 1, 4, INDEP, None).expect("issues");
        let mut s = schoenhage_identity();
        s.name.push_str(" v2");
        let b = issue(s, 1, 4, INDEP, None).expect("issues");
        assert_ne!(a.id, b.id);
        assert_ne!(a.value_digest, b.value_digest);
        // Same statement, same seed: deterministic.
        let c = issue(schoenhage_identity(), 1, 4, INDEP, None).expect("issues");
        assert_eq!(a, c);
    }

    #[test]
    fn degree_bound_is_checked_syntactically_not_trusted() {
        let mut s = schoenhage_identity();
        s.degree_bound = 2;
        assert_eq!(
            validate_statement(&s),
            Err(Malformed::DegreeBound {
                declared: 2,
                syntactic: 3
            })
        );
        let mut s = schoenhage_identity();
        s.lhs = Expr::plus(s.lhs, Expr::var("w"));
        assert_eq!(
            validate_statement(&s),
            Err(Malformed::UnknownVariable("w".to_string()))
        );
    }

    #[test]
    fn field_arithmetic_agrees_with_u128_reference() {
        let mut st = 0x1234_5678u64;
        for _ in 0..2000 {
            let a = draw_field(&mut st);
            let b = draw_field(&mut st);
            let want = ((a as u128 * b as u128) % MODULUS as u128) as u64;
            assert_eq!(fmul(a, b), want);
            assert_eq!(
                fadd(a, b),
                ((a as u128 + b as u128) % MODULUS as u128) as u64
            );
            assert_eq!(fadd(a, fneg(a)), 0);
        }
        assert_eq!(fmul(MODULUS - 1, MODULUS - 1), 1);
        assert_eq!(fconst(-1), MODULUS - 1);
        assert_eq!(
            fconst(i64::MIN),
            ((i64::MIN as i128).rem_euclid(MODULUS as i128)) as u64
        );
    }

    #[test]
    fn canonical_bytes_are_sorted_and_round_trip() {
        let s = schoenhage_identity();
        let bytes = canonical_statement_bytes(&s);
        let text = String::from_utf8(bytes.clone()).unwrap();
        assert!(text.starts_with("{\"degree_bound\":3,\"lhs\":"));
        let back: IdentityStatement = serde_json::from_slice(&bytes).unwrap();
        assert_eq!(back, s);
        assert_eq!(statement_id(&back), statement_id(&s));
    }

    /// The committed worked example must keep checking: this pins the file in
    /// `docs/identity-certificates/` the way conformance vectors pin an encoding.
    #[test]
    fn committed_schoenhage_certificate_still_checks() {
        let text = include_str!("../../docs/identity-certificates/schoenhage-10mult.json");
        let cert: IdentityCertificate = serde_json::from_str(text).expect("parses");
        assert_eq!(cert.statement, schoenhage_identity());
        assert!(matches!(check(&cert), Ok(Verdict::Accept { .. })));
    }
}

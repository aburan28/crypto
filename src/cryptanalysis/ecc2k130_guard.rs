//! Guards for the ECC2K-130 Pollard rho campaign: a certificate that the
//! shipping iteration has no fruitless cycle, and a stored, verifiable
//! campaign seed.
//!
//! # No fruitless cycles
//!
//! The campaign walks `R ← σʲ(R) + R` with `j = 3 + ((HW(x_R) / 2) mod 8)`,
//! the Hamming weight taken in normal-basis coordinates.  `σ` permutes those
//! coordinates and negation fixes `x`, so `j` is constant on each orbit
//! `{±σᵃ(R)}` and the iteration is equivariant: it descends to the orbit set
//! with no canonical representative inside the loop, so there is no sign for
//! a fruitless cycle to feed back through (`ecc2k130/Makefile`,
//! `check-cycles`).
//!
//! What equivariance leaves is arithmetic.  On the order-`ℓ` subgroup `σ`
//! acts as multiplication by `λ`, a step multiplies the walk's scalar by
//! `λʲ + 1`, and the multipliers commute, so a walk returns to its own orbit
//! after `k` steps exactly when the product over the multiset of exponents it
//! used satisfies
//!
//! ```text
//! ∏ (λ^{j_t} + 1) ≡ ±λⁱ  (mod ℓ)   for some i.
//! ```
//!
//! [`certify`] decides that for every multiset up to a length bound — an
//! exhaustive check of a finite space, so a clean result rules out every
//! cycle of that length or shorter — and checks the two degeneracies that
//! would break the addition itself (`λʲ = 1`, a doubling; `λʲ = −1`,
//! infinity).  It derives `ℓ` from the curve-order recurrence and `λ` from
//! `√−7`, then cross-checks both against the constants `codegen/gen.py`
//! verified on the challenge point itself (`[λ]P = σ(P)`), which the Python
//! check it replaces never compared.
//!
//! # Seeds
//!
//! A live walk seed is `run id (16) ‖ walk index (32) ‖ restart counter
//! (16)` and the client expands it with a splitmix64 finaliser into the 128
//! bits that pick the start point `Q + Σ cᵢ σⁱ(P)` ([`legacy_seed`],
//! [`legacy_start_bits`]).  That derivation has no entropy by design: every
//! distinguished point is replayed from its recorded 64-bit seed, which is
//! how a collision is solved.  [`generate_seed`] adds what it lacks: a
//! 256-bit master seed drawn from the operating system's CSPRNG (or a
//! recorded public beacon), bound to the published challenge instance
//! through BLAKE3's key-derivation mode, with every input stored so anyone
//! can re-derive it ([`verify_seed`]), and a keyed-BLAKE3 start derivation
//! with test vectors for the next corpus boundary ([`walk_start_bits`]).

use std::collections::HashMap;

use num_bigint::{BigInt, BigUint};
use num_traits::{One, ToPrimitive, Zero};
use serde::{Deserialize, Serialize};

/// Step exponents of the shipping walk, `j = 3 + ((HW / 2) mod 8)`.
pub const STEP_EXPONENTS: [u32; 8] = [3, 4, 5, 6, 7, 8, 9, 10];
/// Field degrees the client supports.
pub const CURVES: [u32; 5] = [23, 41, 83, 97, 131];
/// Hits a curve's chance expectation explains: at the toy degrees random
/// products land on an orbit scalar regularly, which says nothing about
/// `m = 131`.
const COINCIDENCE_EXPECTED: f64 = 0.01;

/// `ℓ` for ECC2K-130, as `codegen/gen.py` writes it (`ELL_DEC`).
pub const CHALLENGE_ELL: &str = "680564733841876926932320129493409985129";
/// `λ` with `σ(R) = [λ]R` on the challenge subgroup (`S_DEC`), which
/// `codegen/gen.py` selects by checking it against the challenge point.
pub const CHALLENGE_LAMBDA: &str = "196511074115861092422032515080945363956";

/// `#E(GF(2^m))` for `y² + xy = x³ + 1`: `2^m + 1 − V_m` with `V_0 = 2`,
/// `V_1 = −1`, `V_{k+1} = −V_k − 2V_{k−1}`.
pub fn koblitz_order(m: u32) -> BigUint {
    assert!(m >= 1, "a field degree is at least 1");
    let (mut v0, mut v1) = (BigInt::from(2), BigInt::from(-1));
    for _ in 1..m {
        let next = -&v1 - (&v0 << 1usize);
        v0 = std::mem::replace(&mut v1, next);
    }
    ((BigInt::one() << (m as usize)) + 1u32 - v1)
        .to_biguint()
        .expect("a group order is positive")
}

const SMALL_PRIMES: [u32; 20] = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
];

/// Miller–Rabin with the first twenty primes as bases.  A composite passes
/// with probability at most `4⁻²⁰` over random inputs; `m = 131` does not
/// rest on it, because its `ℓ` is also compared with [`CHALLENGE_ELL`].
pub fn is_probable_prime(n: &BigUint) -> bool {
    if *n < BigUint::from(2u32) {
        return false;
    }
    for &p in &SMALL_PRIMES {
        let p = BigUint::from(p);
        if *n == p {
            return true;
        }
        if (n % &p).is_zero() {
            return false;
        }
    }
    let one = BigUint::one();
    let n_minus_1 = n - &one;
    let mut d = n_minus_1.clone();
    let mut r = 0u32;
    while !d.bit(0) {
        d >>= 1usize;
        r += 1;
    }
    'witness: for &a in &SMALL_PRIMES {
        let mut x = BigUint::from(a).modpow(&d, n);
        if x == one || x == n_minus_1 {
            continue;
        }
        for _ in 1..r {
            x = &x * &x % n;
            if x == n_minus_1 {
                continue 'witness;
            }
        }
        return false;
    }
    true
}

/// A square root of `a` modulo the odd prime `p` (Tonelli–Shanks), or
/// `None` when `a` is a non-residue.
pub fn sqrt_mod(a: &BigUint, p: &BigUint) -> Option<BigUint> {
    let one = BigUint::one();
    let a = a % p;
    if a.is_zero() {
        return Some(a);
    }
    let p_minus_1 = p - &one;
    let half = &p_minus_1 >> 1usize;
    if a.modpow(&half, p) != one {
        return None;
    }
    let mut q = p_minus_1.clone();
    let mut s = 0u32;
    while !q.bit(0) {
        q >>= 1usize;
        s += 1;
    }
    let mut z = BigUint::from(2u32);
    while z.modpow(&half, p) != p_minus_1 {
        z += 1u32;
    }
    let mut m = s;
    let mut c = z.modpow(&q, p);
    let mut t = a.modpow(&q, p);
    let mut r = a.modpow(&((&q + &one) >> 1usize), p);
    while t != one {
        let mut i = 0u32;
        let mut tt = t.clone();
        while tt != one {
            tt = &tt * &tt % p;
            i += 1;
        }
        let b = c.modpow(&(BigUint::one() << ((m - i - 1) as usize)), p);
        m = i;
        c = &b * &b % p;
        t = &t * &c % p;
        r = &r * &b % p;
    }
    Some(r)
}

/// `λ` with `σ(R) = [λ]R` on the order-`ℓ` subgroup: the root of
/// `λ² + λ + 2 = 0` whose order divides `m`.
pub fn frobenius_eigenvalue(m: u32, ell: &BigUint) -> Option<BigUint> {
    let seven = BigUint::from(7u32);
    if *ell <= seven {
        return None;
    }
    let r = sqrt_mod(&(ell - &seven), ell)?;
    let inv2 = BigUint::from(2u32).modpow(&(ell - 2u32), ell);
    let minus_one = ell - 1u32;
    let exponent = BigUint::from(m);
    [&minus_one + &r, &minus_one + (ell - &r)]
        .into_iter()
        .map(|numerator| numerator % ell * &inv2 % ell)
        .find(|candidate| candidate.modpow(&exponent, ell).is_one())
}

/// A multiset of step exponents whose product lands on an orbit scalar.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct CycleHit {
    pub length: usize,
    pub exponents: Vec<u32>,
    pub orbit_scalar: String,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CycleVerdict {
    /// Not a curve of order `4ℓ` with a Frobenius eigenvalue of order `m`.
    Skipped,
    /// No multiset up to the bound lands on an orbit.
    Clean,
    /// Hits, but as many as chance predicts at this toy size.
    Coincidence,
    /// A hit chance does not explain: a walk can loop without ever reaching
    /// a distinguished point.
    Real,
    /// `λʲ = ±1` for a step exponent: the addition itself breaks.
    Degenerate,
}

impl CycleVerdict {
    pub fn fails(self) -> bool {
        matches!(self, CycleVerdict::Real | CycleVerdict::Degenerate)
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CurveCertificate {
    pub m: u32,
    pub verdict: CycleVerdict,
    pub skip_reason: Option<String>,
    pub group_order: String,
    pub ell: Option<String>,
    pub lambda: Option<String>,
    pub multisets_checked: u64,
    pub orbit_scalars: usize,
    pub hits: Vec<CycleHit>,
    /// Hits a uniformly random product would give over the same search.
    pub expected_by_chance: f64,
    pub degenerate: Vec<String>,
    /// `ℓ` and `λ` against the generated header, for `m = 131` only.
    pub matches_challenge_constants: Option<bool>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CycleCertificate {
    pub schema: String,
    pub iteration: String,
    pub step_exponents: Vec<u32>,
    pub max_length: usize,
    pub curves: Vec<CurveCertificate>,
    pub passed: bool,
}

struct Search<'a> {
    ell: &'a BigUint,
    steps: &'a [(u32, BigUint)],
    targets: &'a HashMap<BigUint, String>,
    max_length: usize,
    prefix: Vec<u32>,
    checked: u64,
    hits: Vec<CycleHit>,
}

impl Search<'_> {
    /// Every non-decreasing exponent sequence extending `prefix`: products
    /// commute, so a cycle is fixed by which exponents it uses and how
    /// often, not by their order.
    fn extend(&mut self, from: usize, product: &BigUint) {
        if self.prefix.len() == self.max_length {
            return;
        }
        for i in from..self.steps.len() {
            let (j, factor) = &self.steps[i];
            let next = product * factor % self.ell;
            self.prefix.push(*j);
            self.checked += 1;
            if let Some(label) = self.targets.get(&next) {
                self.hits.push(CycleHit {
                    length: self.prefix.len(),
                    exponents: self.prefix.clone(),
                    orbit_scalar: label.clone(),
                });
            }
            self.extend(i, &next);
            self.prefix.pop();
        }
    }
}

/// The orbit scalars `±λⁱ`, labelled.
fn orbit_scalars(lambda: &BigUint, m: u32, ell: &BigUint) -> HashMap<BigUint, String> {
    let mut out = HashMap::with_capacity(2 * m as usize);
    let mut power = BigUint::one();
    for i in 0..m {
        out.insert(ell - &power, format!("-lam^{i}"));
        out.insert(power.clone(), format!("+lam^{i}"));
        power = power * lambda % ell;
    }
    out
}

/// [`certify_curve_with`] on the shipping schedule.
pub fn certify_curve(m: u32, max_length: usize) -> CurveCertificate {
    certify_curve_with(m, max_length, &STEP_EXPONENTS)
}

/// The certificate for one curve and one step schedule.
pub fn certify_curve_with(m: u32, max_length: usize, exponents: &[u32]) -> CurveCertificate {
    let order = koblitz_order(m);
    let mut cert = CurveCertificate {
        m,
        verdict: CycleVerdict::Skipped,
        skip_reason: None,
        group_order: order.to_string(),
        ell: None,
        lambda: None,
        multisets_checked: 0,
        orbit_scalars: 0,
        hits: Vec::new(),
        expected_by_chance: 0.0,
        degenerate: Vec::new(),
        matches_challenge_constants: None,
    };
    let four = BigUint::from(4u32);
    if !(&order % &four).is_zero() {
        cert.skip_reason = Some("the group order is not divisible by 4".into());
        return cert;
    }
    let ell = &order / &four;
    if !is_probable_prime(&ell) {
        cert.skip_reason = Some("the group order over 4 is not prime".into());
        return cert;
    }
    let Some(lambda) = frobenius_eigenvalue(m, &ell) else {
        cert.skip_reason = Some("no Frobenius eigenvalue of order dividing m".into());
        return cert;
    };
    let minus_one = &ell - 1u32;
    let mut steps = Vec::with_capacity(exponents.len());
    for &j in exponents {
        let power = lambda.modpow(&BigUint::from(j), &ell);
        if power.is_one() {
            cert.degenerate.push(format!(
                "j = {j}: sigma^j is the identity, the step is a doubling"
            ));
        }
        if power == minus_one {
            cert.degenerate.push(format!(
                "j = {j}: sigma^j(R) = -R, the step lands on infinity"
            ));
        }
        steps.push((j, (power + 1u32) % &ell));
    }
    let targets = orbit_scalars(&lambda, m, &ell);
    let mut search = Search {
        ell: &ell,
        steps: &steps,
        targets: &targets,
        max_length,
        prefix: Vec::with_capacity(max_length),
        checked: 0,
        hits: Vec::new(),
    };
    search.extend(0, &BigUint::one());
    let ell_f = ell.to_f64().unwrap_or(f64::INFINITY);
    cert.expected_by_chance = search.checked as f64 * targets.len() as f64 / ell_f;
    cert.multisets_checked = search.checked;
    cert.orbit_scalars = targets.len();
    cert.hits = search.hits;
    cert.verdict = if !cert.degenerate.is_empty() {
        CycleVerdict::Degenerate
    } else if cert.hits.is_empty() {
        CycleVerdict::Clean
    } else if cert.expected_by_chance > COINCIDENCE_EXPECTED {
        CycleVerdict::Coincidence
    } else {
        CycleVerdict::Real
    };
    if m == 131 {
        cert.matches_challenge_constants =
            Some(ell.to_string() == CHALLENGE_ELL && lambda.to_string() == CHALLENGE_LAMBDA);
    }
    cert.ell = Some(ell.to_string());
    cert.lambda = Some(lambda.to_string());
    cert
}

/// The certificate over `curves`; it passes when no curve has a real hit
/// or a degeneracy and `m = 131`, if present, matches the challenge.
pub fn certify(curves: &[u32], max_length: usize) -> CycleCertificate {
    let curves: Vec<CurveCertificate> = curves
        .iter()
        .map(|&m| certify_curve(m, max_length))
        .collect();
    let passed = curves
        .iter()
        .all(|c| !c.verdict.fails() && c.matches_challenge_constants != Some(false));
    CycleCertificate {
        schema: "ecc2k130-no-fruitless-cycle-v1".into(),
        iteration: "R <- sigma^j(R) + R, j = 3 + ((HW(x_R) / 2) mod 8)".into(),
        step_exponents: STEP_EXPONENTS.to_vec(),
        max_length,
        curves,
        passed,
    }
}

// ── Seeds ──────────────────────────────────────────────────────────

/// The live client's seed for walk `walk` of run `run_id` after `restart`
/// restarts: `include/walk.h` `eccSeedFor` plus one per restart.
pub fn legacy_seed(run_id: u16, walk: u32, restart: u16) -> u64 {
    (u64::from(run_id) << 48) | (u64::from(walk) << 16) | u64::from(restart)
}

/// `(run id, walk index, restart counter)` of a live seed.
pub fn decode_legacy_seed(seed: u64) -> (u16, u32, u16) {
    ((seed >> 48) as u16, (seed >> 16) as u32, seed as u16)
}

/// `include/walk.h` `eccPrf`: a splitmix64 finaliser over
/// `seed + γ·(idx + 1)`.  A bijection of the seed, not a keyed PRF.
pub fn legacy_prf(seed: u64, idx: u64) -> u64 {
    let mut z = seed.wrapping_add(0x9E37_79B9_7F4A_7C15u64.wrapping_mul(idx + 1));
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// The 128 bits `c` the live client draws a start point `Q + Σ cᵢσⁱ(P)`
/// from: `[eccPrf(seed, 0), eccPrf(seed, 1)]`, bit `i` of the pair being
/// the coefficient of `σⁱ(P)`.
pub fn legacy_start_bits(seed: u64) -> [u64; 2] {
    [legacy_prf(seed, 0), legacy_prf(seed, 1)]
}

/// The published ECC2K-130 instance a campaign seed is bound to: the
/// Certicom challenge, which no one running the campaign chose.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ChallengeAnchor {
    pub campaign: String,
    pub curve: String,
    pub field_polynomial: String,
    pub ell: String,
    pub px: String,
    pub py: String,
    pub qx: String,
    pub qy: String,
}

impl ChallengeAnchor {
    /// The challenge as `codegen/gen.py` carries it (polynomial basis).
    pub fn ecc2k130(campaign: &str) -> Self {
        Self {
            campaign: campaign.into(),
            curve: "certicom-ecc2k-130: y^2 + xy = x^3 + 1 over GF(2^131)".into(),
            field_polynomial: "x^131 + x^13 + x^2 + x + 1".into(),
            ell: CHALLENGE_ELL.into(),
            px: "051c99bfa6f18de467c80c23b98c7994aa".into(),
            py: "042ea2d112ecec71fcf7e000d7efc978bd".into(),
            qx: "06c997f3e7f2c66a4a5d2fda13756a37b1".into(),
            qy: "04a38d11829d32d347bd0c0f584d546e9a".into(),
        }
    }

    fn canonical(&self) -> Vec<u8> {
        format!(
            "campaign={};curve={};field={};ell={};px={};py={};qx={};qy={}",
            self.campaign,
            self.curve,
            self.field_polynomial,
            self.ell,
            self.px,
            self.py,
            self.qx,
            self.qy
        )
        .into_bytes()
    }
}

/// BLAKE3 `derive_key` context: globally unique, fixed, never reused for
/// another purpose.
pub const SEED_KDF_CONTEXT: &str = "crypto ecc2k130 2026-09-25 rho campaign master seed v1";
const SEED_SCHEMA: &str = "ecc2k130-campaign-seed-v1";
const WALK_START_DERIVATION: &str = "first 16 bytes of BLAKE3-keyed(master seed, \"walk-start\" || run_id u16 LE || walk u32 LE || restart u16 LE); bit i selects sigma^i(P)";
const WALK_START_VECTORS: [(u16, u32, u16); 4] =
    [(1, 0, 0), (1, 1, 0), (196, 6_160_383, 7), (8000, 0, 0)];

/// A walk start under a master seed: the 128 bits the next corpus would
/// draw in place of [`legacy_start_bits`].  Same seed layout, so every
/// distinguished point still replays from its recorded 64-bit seed.
pub fn walk_start_bits(master: &[u8; 32], run_id: u16, walk: u32, restart: u16) -> [u8; 16] {
    let mut hasher = blake3::Hasher::new_keyed(master);
    hasher.update(b"walk-start");
    hasher.update(&run_id.to_le_bytes());
    hasher.update(&walk.to_le_bytes());
    hasher.update(&restart.to_le_bytes());
    let mut out = [0u8; 16];
    out.copy_from_slice(&hasher.finalize().as_bytes()[..16]);
    out
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct WalkStartVector {
    pub run_id: u16,
    pub walk: u32,
    pub restart: u16,
    pub legacy_seed_hex: String,
    pub start_bits_hex: String,
}

/// Everything a campaign seed was derived from, and the seed.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct SeedProvenance {
    pub schema: String,
    pub kdf: String,
    pub kdf_context: String,
    pub anchor: ChallengeAnchor,
    pub anchor_blake3: String,
    /// `getrandom: OS CSPRNG`, or `beacon: <label>` for a recorded public value.
    pub entropy_source: String,
    pub entropy_hex: String,
    pub master_seed_hex: String,
    pub created_unix: u64,
    pub walk_start_derivation: String,
    pub walk_start_vectors: Vec<WalkStartVector>,
    /// BLAKE3 of this record serialised with this field empty.
    pub commitment_blake3: String,
}

fn master_seed(anchor: &ChallengeAnchor, entropy: &[u8]) -> [u8; 32] {
    let mut material = anchor.canonical();
    material.push(0);
    material.extend_from_slice(entropy);
    blake3::derive_key(SEED_KDF_CONTEXT, &material)
}

fn commitment(p: &SeedProvenance) -> String {
    let mut bare = p.clone();
    bare.commitment_blake3 = String::new();
    let bytes = serde_json::to_vec(&bare).expect("a provenance record serialises");
    blake3::hash(&bytes).to_hex().to_string()
}

fn walk_start_vectors(master: &[u8; 32]) -> Vec<WalkStartVector> {
    WALK_START_VECTORS
        .iter()
        .map(|&(run_id, walk, restart)| WalkStartVector {
            run_id,
            walk,
            restart,
            legacy_seed_hex: format!("{:016x}", legacy_seed(run_id, walk, restart)),
            start_bits_hex: hex::encode(walk_start_bits(master, run_id, walk, restart)),
        })
        .collect()
}

/// 32 bytes from the operating system's CSPRNG (`getrandom(2)` on Linux).
pub fn os_entropy() -> Result<[u8; 32], String> {
    let mut buf = [0u8; 32];
    getrandom::getrandom(&mut buf).map_err(|e| format!("OS CSPRNG unavailable: {e}"))?;
    Ok(buf)
}

/// A campaign seed from `entropy`, bound to `anchor`.
pub fn generate_seed(
    anchor: ChallengeAnchor,
    entropy: &[u8],
    entropy_source: &str,
    created_unix: u64,
) -> SeedProvenance {
    let master = master_seed(&anchor, entropy);
    let mut p = SeedProvenance {
        schema: SEED_SCHEMA.into(),
        kdf: "BLAKE3 derive_key".into(),
        kdf_context: SEED_KDF_CONTEXT.into(),
        anchor_blake3: blake3::hash(&anchor.canonical()).to_hex().to_string(),
        anchor,
        entropy_source: entropy_source.into(),
        entropy_hex: hex::encode(entropy),
        master_seed_hex: hex::encode(master),
        created_unix,
        walk_start_derivation: WALK_START_DERIVATION.into(),
        walk_start_vectors: walk_start_vectors(&master),
        commitment_blake3: String::new(),
    };
    p.commitment_blake3 = commitment(&p);
    p
}

/// Re-derive every field of a seed record; each disagreement is reported.
pub fn verify_seed(p: &SeedProvenance) -> Result<(), Vec<String>> {
    let mut failures = Vec::new();
    if p.schema != SEED_SCHEMA {
        failures.push(format!("schema {:?}, expected {SEED_SCHEMA:?}", p.schema));
    }
    if p.kdf_context != SEED_KDF_CONTEXT {
        failures.push(format!(
            "KDF context {:?} is not this tool's",
            p.kdf_context
        ));
    }
    if p.anchor != ChallengeAnchor::ecc2k130(&p.anchor.campaign) {
        failures.push("the anchor is not the published ECC2K-130 challenge".into());
    }
    if p.anchor_blake3 != blake3::hash(&p.anchor.canonical()).to_hex().to_string() {
        failures.push("anchor_blake3 does not hash the recorded anchor".into());
    }
    let entropy = match hex::decode(&p.entropy_hex) {
        Ok(e) if e.len() >= 32 => Some(e),
        Ok(e) => {
            failures.push(format!(
                "{} bytes of entropy; at least 32 are required",
                e.len()
            ));
            None
        }
        Err(e) => {
            failures.push(format!("entropy_hex is not hex: {e}"));
            None
        }
    };
    if let Some(entropy) = entropy {
        let master = master_seed(&p.anchor, &entropy);
        if p.master_seed_hex != hex::encode(master) {
            failures.push("master_seed_hex does not re-derive from the recorded inputs".into());
        }
        if p.walk_start_vectors != walk_start_vectors(&master) {
            failures.push("walk start vectors do not re-derive from the master seed".into());
        }
    }
    if p.walk_start_derivation != WALK_START_DERIVATION {
        failures.push("walk_start_derivation is not this tool's".into());
    }
    if p.commitment_blake3 != commitment(p) {
        failures.push("commitment_blake3 does not match the record".into());
    }
    if failures.is_empty() {
        Ok(())
    } else {
        Err(failures)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn big(s: &str) -> BigUint {
        BigUint::parse_bytes(s.as_bytes(), 10).expect("decimal")
    }

    #[test]
    fn the_order_recurrence_gives_the_challenge_group_order() {
        assert_eq!(
            koblitz_order(131),
            big("2722258935367507707729280517973639940516")
        );
        assert_eq!(koblitz_order(131), big(CHALLENGE_ELL) * 4u32);
        assert_eq!(koblitz_order(23), BigUint::from(8_383_412u32));
    }

    #[test]
    fn miller_rabin_separates_primes_from_pseudoprimes() {
        for p in [2u64, 3, 71, 73, 2_095_853, 549_756_390_943] {
            assert!(is_probable_prime(&BigUint::from(p)), "{p}");
        }
        // 561 and 41041 are Carmichael numbers, 2047 a base-2 strong pseudoprime.
        for c in [0u64, 1, 4, 561, 2047, 41_041, 3_215_031_751, 2_095_853 * 3] {
            assert!(!is_probable_prime(&BigUint::from(c)), "{c}");
        }
        assert!(is_probable_prime(&big(CHALLENGE_ELL)));
    }

    #[test]
    fn square_roots_square_back() {
        let ell = big(CHALLENGE_ELL);
        let minus_seven = &ell - 7u32;
        let r = sqrt_mod(&minus_seven, &ell).expect("-7 is a square mod ell");
        assert_eq!(&r * &r % &ell, minus_seven);
        let p = BigUint::from(2_095_853u32);
        assert!(sqrt_mod(&BigUint::from(0u32), &p).unwrap().is_zero());
    }

    #[test]
    fn the_eigenvalue_is_the_one_the_generator_checked_on_the_point() {
        let ell = big(CHALLENGE_ELL);
        let lambda = frobenius_eigenvalue(131, &ell).expect("an eigenvalue of order 131");
        assert_eq!(lambda, big(CHALLENGE_LAMBDA));
        // λ² + λ + 2 ≡ 0, the Frobenius characteristic polynomial.
        assert!(((&lambda * &lambda + &lambda + 2u32) % &ell).is_zero());
    }

    /// The numbers `codegen/cycles.py --max-length 8` printed on 2026-09-25.
    #[test]
    fn every_curve_reproduces_the_python_check() {
        let reference = [
            (23, "2095853", "93194", 46, 0.2824501527540338),
            (
                41,
                "549756390943",
                "256851699273",
                82,
                1.9195011051893556e-06,
            ),
            (
                83,
                "2417851639230796216685689",
                "254512724090651164922414",
                166,
                8.835339461438658e-19,
            ),
            (
                97,
                "39614081257132074233778707191",
                "23986389105595106179003010550",
                194,
                6.302269094150751e-23,
            ),
            (
                131,
                CHALLENGE_ELL,
                CHALLENGE_LAMBDA,
                262,
                4.954235552239736e-33,
            ),
        ];
        for (m, ell, lambda, targets, expected) in reference {
            let c = certify_curve(m, 8);
            assert_eq!(c.verdict, CycleVerdict::Clean, "m = {m}");
            assert_eq!(c.ell.as_deref(), Some(ell), "m = {m}");
            assert_eq!(c.lambda.as_deref(), Some(lambda), "m = {m}");
            assert_eq!(c.multisets_checked, 12_869, "m = {m}");
            assert_eq!(c.orbit_scalars, targets, "m = {m}");
            assert!(c.hits.is_empty() && c.degenerate.is_empty(), "m = {m}");
            let rel = (c.expected_by_chance - expected).abs() / expected;
            assert!(
                rel < 1e-9,
                "m = {m}: {} vs {expected}",
                c.expected_by_chance
            );
        }
        assert_eq!(
            certify_curve(131, 8).matches_challenge_constants,
            Some(true)
        );
    }

    #[test]
    fn the_campaign_certificate_passes() {
        let cert = certify(&CURVES, 8);
        assert!(cert.passed);
        assert_eq!(cert.curves.len(), CURVES.len());
    }

    /// A positive control.  At `m = 19` two of the eight multipliers share a
    /// class (`WALK-CONSTANT.md` §4), so short cycles exist; the search must
    /// find exactly the ones `codegen/cycles.py` found, and call them what
    /// the chance expectation at that size says they are.
    #[test]
    fn the_toy_curve_with_real_cycles_is_found_and_called_a_coincidence() {
        let c = certify_curve(19, 8);
        assert_eq!(c.ell.as_deref(), Some("130873"));
        assert_eq!(c.lambda.as_deref(), Some("41811"));
        assert_eq!(c.verdict, CycleVerdict::Coincidence);
        let mut hits: Vec<(Vec<u32>, String)> = c
            .hits
            .iter()
            .map(|h| (h.exponents.clone(), h.orbit_scalar.clone()))
            .collect();
        hits.sort();
        assert_eq!(
            hits,
            vec![
                (vec![3, 5, 6, 7, 9], "-lam^16".to_string()),
                (vec![3, 5, 6, 7, 10], "-lam^7".to_string()),
            ]
        );
        assert!((c.expected_by_chance - 3.7366148861873727).abs() < 1e-9);
    }

    #[test]
    fn a_curve_that_is_not_four_times_a_prime_is_skipped() {
        let c = certify_curve(29, 8);
        assert_eq!(c.verdict, CycleVerdict::Skipped);
        assert!(c.skip_reason.is_some());
    }

    #[test]
    fn a_step_that_is_a_multiple_of_m_is_caught_as_a_doubling() {
        let c = certify_curve_with(131, 2, &[3, 131]);
        assert_eq!(c.verdict, CycleVerdict::Degenerate);
        assert!(c.degenerate[0].contains("doubling"), "{:?}", c.degenerate);
        assert!(c.verdict.fails());
    }

    #[test]
    fn a_product_on_an_orbit_is_reported_as_a_hit() {
        // A factor of λ is a step that returns to the orbit in one move: the
        // search must see it, with the exact multiset that closes it.
        let ell = big(CHALLENGE_ELL);
        let lambda = big(CHALLENGE_LAMBDA);
        let targets = orbit_scalars(&lambda, 131, &ell);
        let steps = vec![
            (
                3u32,
                (lambda.modpow(&BigUint::from(3u32), &ell) + 1u32) % &ell,
            ),
            (99u32, lambda.clone()),
        ];
        let mut search = Search {
            ell: &ell,
            steps: &steps,
            targets: &targets,
            max_length: 2,
            prefix: Vec::new(),
            checked: 0,
            hits: Vec::new(),
        };
        search.extend(0, &BigUint::one());
        assert_eq!(search.checked, 5);
        assert!(search
            .hits
            .iter()
            .any(|h| h.exponents == vec![99] && h.orbit_scalar == "+lam^1"));
        assert!(search
            .hits
            .iter()
            .any(|h| h.exponents == vec![99, 99] && h.orbit_scalar == "+lam^2"));
    }

    /// `include/walk.h` evaluated independently (Python) on 2026-09-25.
    #[test]
    fn legacy_seeds_match_the_client() {
        let vectors = [
            (
                1u16,
                0u32,
                0u16,
                0x0001_0000_0000_0000u64,
                0xa285_e7b0_deb6_3750u64,
                0x0ecf_0817_eb32_d65au64,
            ),
            (
                1,
                1,
                0,
                0x0001_0000_0001_0000,
                0xdb7b_f4d0_c7c2_7076,
                0xcbe4_16a3_8159_8dd5,
            ),
            (
                196,
                6_160_383,
                7,
                0x00c4_005d_ffff_0007,
                0x5013_806b_8b24_640c,
                0x0ab7_9eaa_cee7_4562,
            ),
            (
                8000,
                0,
                0,
                0x1f40_0000_0000_0000,
                0x81f5_1fa9_fbfd_e021,
                0xdfd7_3c00_bc19_8691,
            ),
        ];
        for (run, walk, restart, seed, c0, c1) in vectors {
            assert_eq!(legacy_seed(run, walk, restart), seed);
            assert_eq!(decode_legacy_seed(seed), (run, walk, restart));
            assert_eq!(legacy_start_bits(seed), [c0, c1]);
        }
    }

    fn fixed() -> SeedProvenance {
        generate_seed(
            ChallengeAnchor::ecc2k130("ecc2k-130"),
            &[7u8; 32],
            "beacon: test vector",
            1_790_000_000,
        )
    }

    #[test]
    fn a_seed_record_verifies_and_re_derives_identically() {
        let p = fixed();
        assert_eq!(verify_seed(&p), Ok(()));
        assert_eq!(p, fixed());
        assert_eq!(p.walk_start_vectors.len(), WALK_START_VECTORS.len());
        let json = serde_json::to_string(&p).unwrap();
        let back: SeedProvenance = serde_json::from_str(&json).unwrap();
        assert_eq!(verify_seed(&back), Ok(()));
    }

    #[test]
    fn every_tampered_field_fails_verification() {
        let tampers: Vec<(&str, Box<dyn Fn(&mut SeedProvenance)>)> = vec![
            (
                "entropy",
                Box::new(|p| p.entropy_hex = hex::encode([8u8; 32])),
            ),
            (
                "master",
                Box::new(|p| p.master_seed_hex = hex::encode([0u8; 32])),
            ),
            ("anchor", Box::new(|p| p.anchor.qx = "00".into())),
            (
                "vector",
                Box::new(|p| p.walk_start_vectors[0].start_bits_hex = "00".into()),
            ),
            ("context", Box::new(|p| p.kdf_context = "other".into())),
            ("source", Box::new(|p| p.entropy_source = "edited".into())),
            (
                "short",
                Box::new(|p| p.entropy_hex = hex::encode([7u8; 16])),
            ),
        ];
        for (what, tamper) in tampers {
            let mut p = fixed();
            tamper(&mut p);
            assert!(verify_seed(&p).is_err(), "{what} was not caught");
        }
    }

    #[test]
    fn os_entropy_gives_distinct_verifiable_seeds() {
        let anchor = ChallengeAnchor::ecc2k130("ecc2k-130");
        let a = generate_seed(
            anchor.clone(),
            &os_entropy().unwrap(),
            "getrandom: OS CSPRNG",
            0,
        );
        let b = generate_seed(anchor, &os_entropy().unwrap(), "getrandom: OS CSPRNG", 0);
        assert_ne!(a.master_seed_hex, b.master_seed_hex);
        assert_eq!(verify_seed(&a), Ok(()));
        assert_eq!(verify_seed(&b), Ok(()));
    }

    #[test]
    fn walk_starts_are_distinct_and_keyed() {
        let m1 = [1u8; 32];
        let m2 = [2u8; 32];
        let base = walk_start_bits(&m1, 1, 0, 0);
        assert_ne!(base, walk_start_bits(&m1, 1, 0, 1));
        assert_ne!(base, walk_start_bits(&m1, 1, 1, 0));
        assert_ne!(base, walk_start_bits(&m1, 2, 0, 0));
        assert_ne!(base, walk_start_bits(&m2, 1, 0, 0));
    }
}

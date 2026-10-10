//! **P-256 backdoor map** — what is confirmed, what is unresolved, and what
//! a curve-level trapdoor would have to be.
//!
//! Companion to `src/cryptanalysis/p256_speculation.rs`, which records
//! the four *public concern vectors* and their probes.  This module encodes
//! the *assessment* that sits on top of those probes as auditable records,
//! so the reasoning can be reviewed, extended, and falsified line by line
//! instead of living in chat transcripts:
//!
//! 1. [`dual_ec_drbg_record`] — the one **confirmed** NSA backdoor that
//!    involves P-256.  It lives in a random-number generator built on two
//!    P-256 points, not in the curve.
//! 2. [`verify_p256_seed_derivation`] — the ANSI X9.62 / FIPS 186 seed →
//!    `b` derivation, **re-executed** against the published seed and the
//!    real P-256 constants.  This is the part of "verifiably random" that
//!    actually verifies; [`seed_provenance`] records the part that does not
//!    (nothing constrains the seed itself).
//! 3. [`p256_trapdoor_hypotheses`] — the four ranked places a hidden
//!    curve-level weakness would have to live, each with the grind-cost
//!    model it implies and the public evidence against it.
//!    [`GrindBudget`] turns "could the generator have ground seeds to hit a
//!    weak class of density 2^−k?" into arithmetic on stated assumptions.
//! 4. [`p256_backdoor_assessment`] / [`format_backdoor_assessment`] — the
//!    overall verdict and a Markdown rendering for research notes.
//!
//! **Status: speculation with falsifiers.**  Every record carries a
//! [`ProbeVerdict`]; none claims a weakness was found.  The curve-level hypotheses are ranked by how
//! compatible they are with the public record and with 1999-era compute,
//! not by any measurement.  The numbers in [`GrindBudget::default_1999`]
//! are assumptions and are labelled as such.
//!
//! The narrative version is `RESEARCH_P256_BACKDOOR_SPECULATION.md`.

use num_bigint::BigUint;

use crate::ecc::curve::CurveParams;
use crate::hash::sha1::sha1;

/// How much public research or evidence exists for a hypothesis (mirrors
/// the enum of the same name in `p256_speculation.rs`, which is not
/// registered as a module on every branch; kept local so this module is
/// self-contained).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum PublicResearchStatus {
    /// A concrete public literature exists and is still being extended.
    Active,
    /// Adjacent attacks exist, but not against prime-field P-256 ECDLP.
    AdjacentOnly,
    /// No meaningful public research program appears to exist.
    NotPubliclyGrounded,
    /// The concern is inherently about unknown unpublished work.
    Conjectural,
}

/// What this codebase can say about P-256 for one hypothesis.  No variant
/// asserts that a weakness was found.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ProbeVerdict {
    /// The public, finite check found no anomaly.
    NoPublicAnomaly,
    /// A structural property exists, but no public EC attack follows.
    StructurePresentNoKnownAttack,
    /// Not meaningfully falsifiable by a public finite check.
    NotPubliclyTestable,
}

// ── Seed provenance: the X9.62 derivation, re-executed ────────────────

/// The published P-256 generation seed (FIPS 186-4 Appendix D.1.2.3,
/// SEC 2 §2.4.2, ANSI X9.62).
pub const P256_SEED: [u8; 20] = [
    0xc4, 0x9d, 0x36, 0x08, 0x86, 0xe7, 0x04, 0x93, 0x6a, 0x66, 0x78, 0xe1, 0x13, 0x9d, 0x26, 0xb7,
    0x81, 0x9f, 0x7e, 0x90,
];

/// The X9.62 intermediate `r` for P-256, pinned so a regression in the
/// derivation code (or in SHA-1) is caught independently of the final
/// `r·b² ≡ a³` check.
pub const P256_X962_R_HEX: &str =
    "7efba1662985be9403cb055c75d4f7e0ce8d84a9c5114abcaf3177680104fa0d";

/// Outcome of re-running the X9.62 A.3.3.1 derivation for one curve.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeedDerivation {
    pub seed: Vec<u8>,
    /// `s = ⌊(t − 1) / 160⌋` and `v = t − 160·s` for a `t`-bit prime.
    pub s: u32,
    pub v: u32,
    /// The derived integer `r` (the `W0 ‖ W1 ‖ … ‖ Ws` bit string).
    pub r: BigUint,
    /// `r·b² ≡ a³ (mod p)` with the curve's published `a`, `b`.
    pub holds: bool,
}

/// ANSI X9.62 A.3.3.1 candidate `r` from a 160-bit seed for a `t`-bit
/// prime field: `H = SHA-1(seed)`, `c0 =` rightmost `v` bits of `H`,
/// `W0 = c0` with its leftmost bit cleared, `W_i = SHA-1((seed + i) mod
/// 2^160)` for `i = 1..=s`, `r = W0 ‖ W1 ‖ … ‖ Ws`.
pub fn x962_candidate_r(seed: &[u8; 20], t_bits: u32) -> (u32, u32, BigUint) {
    let s = (t_bits - 1) / 160;
    let v = t_bits - 160 * s;
    let h = BigUint::from_bytes_be(&sha1(seed));
    let v_mask = (BigUint::from(1u8) << v) - BigUint::from(1u8);
    let c0 = &h & &v_mask;
    let w0_mask = (BigUint::from(1u8) << (v - 1)) - BigUint::from(1u8);
    let mut w = &c0 & &w0_mask;
    let z = BigUint::from_bytes_be(seed);
    let modulus = BigUint::from(1u8) << 160;
    for i in 1..=s {
        let zi: BigUint = (&z + BigUint::from(i)) % &modulus;
        let mut zi_bytes = zi.to_bytes_be();
        while zi_bytes.len() < 20 {
            zi_bytes.insert(0, 0);
        }
        let wi = BigUint::from_bytes_be(&sha1(&zi_bytes));
        w = (w << 160) | wi;
    }
    (s, v, w)
}

/// Re-run the X9.62 verification for any short-Weierstrass curve: does
/// `r(seed)·b² ≡ a³ (mod p)`?
pub fn verify_seed_derivation(curve: &CurveParams, seed: &[u8; 20]) -> SeedDerivation {
    let t_bits = curve.p.bits() as u32;
    let (s, v, r) = x962_candidate_r(seed, t_bits);
    let p = &curve.p;
    let lhs = (&r * &curve.b % p) * &curve.b % p;
    let rhs = (&curve.a * &curve.a % p) * &curve.a % p;
    SeedDerivation {
        seed: seed.to_vec(),
        s,
        v,
        r,
        holds: lhs == rhs,
    }
}

/// The P-256 instance of [`verify_seed_derivation`] with the published seed.
pub fn verify_p256_seed_derivation() -> SeedDerivation {
    verify_seed_derivation(&CurveParams::p256(), &P256_SEED)
}

/// What the seed story does and does not establish.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeedProvenance {
    pub seed_hex: String,
    /// `true` iff [`verify_p256_seed_derivation`] holds.
    pub derivation_verified: bool,
    pub what_verifies: &'static str,
    pub what_does_not: &'static str,
    pub grinding_model: &'static str,
    pub counter_argument: &'static str,
    pub mundane_account: &'static str,
    pub open_bounty: &'static str,
}

pub fn seed_provenance() -> SeedProvenance {
    SeedProvenance {
        seed_hex: hex::encode(P256_SEED),
        derivation_verified: verify_p256_seed_derivation().holds,
        what_verifies: "b is the X9.62 A.3.3.1 output of SHA-1 on the published seed with a = −3 and the Solinas prime fixed in advance (r·b² ≡ a³ mod p re-executed here).",
        what_does_not: "Nothing constrains the seed. The procedure is pseudo-random from the seed onward, so a generator who knew a secret weak class of curves could try seeds until SHA-1 landed in it (Bernstein–Lange, BADA55; Bernstein et al., 'How to manipulate curve standards').",
        grinding_model: "A weak class of density 2^−k costs ~2^k candidate seeds. Whether that was feasible in 1999 depends entirely on whether membership is testable from (p, a, b) directly (SHA-1 rate) or needs the group order (one SEA point count per candidate); see GrindBudget.",
        counter_argument: "Koblitz–Menezes, 'A riddle wrapped in an enigma' (2015): the weak class would have to stay unknown to the open community for 25 years, and NSA put P-256 and P-384 into Suite B for its own classified traffic, so a deliberately weak curve would have exposed US systems.",
        mundane_account: "A former NSA IAD technical director has said the seeds were SHA-1 of English phrases that the generator (Jerry Solinas) later lost.",
        open_bounty: "A public bounty for the seed preimages (announced 2023) was unclaimed as of this writing; a preimage would settle the question in the mundane direction.",
    }
}

// ── The confirmed backdoor ────────────────────────────────────────────

/// Dual_EC_DRBG, the one documented NSA backdoor that uses P-256.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DualEcRecord {
    pub standard: &'static str,
    pub mechanism: &'static str,
    pub shown_by: &'static str,
    pub attribution: &'static str,
    pub withdrawn: &'static str,
    pub field_incident: &'static str,
    /// Whether the backdoor says anything about ECDLP on the curve.
    pub implicates_curve_math: bool,
}

pub fn dual_ec_drbg_record() -> DualEcRecord {
    DualEcRecord {
        standard: "NIST SP 800-90A (2006), Dual Elliptic Curve Deterministic Random Bit Generator, instantiated on P-256 with two fixed points P and Q.",
        mechanism: "Output is x-coordinates of scalar multiples of P and Q. If Q = [d]P for a secret d, roughly 30 bytes of output reveal the internal state, after which every later output is predictable.",
        shown_by: "Shumow and Ferguson, CRYPTO 2007 rump session: the trapdoor is structurally possible and the constants' provenance was unexplained.",
        attribution: "2013 Snowden documents (Bullrun) and the reported RSA BSAFE default-on arrangement pointed to NSA as the author of the constants.",
        withdrawn: "NIST removed Dual_EC_DRBG from SP 800-90A in 2014.",
        field_incident: "Juniper ScreenOS (2015): an unknown party replaced Q with their own point, demonstrating the mechanism against deployed VPN traffic.",
        implicates_curve_math: false,
    }
}

// ── Curve-level trapdoor hypotheses ───────────────────────────────────

/// Where a hypothesised P-256 weakness would have to live.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum TrapdoorLocus {
    /// In `b`, reached by grinding the seed.
    CurveCoefficient,
    /// In the Solinas prime, which was chosen openly for speed.
    FieldPrime,
    /// In the endomorphism ring / isogeny class.
    EndomorphismOrIsogeny,
    /// A non-public ECDLP algorithm that applies to generic curves.
    UnknownAlgorithm,
}

/// What a seed grind for the hypothesis would have cost per candidate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GrindCostModel {
    /// Membership testable from `(p, a, b)` alone: one SHA-1 plus a cheap
    /// algebraic test per candidate.
    CoefficientTest,
    /// Membership needs the group order: one 256-bit SEA point count per
    /// candidate.
    PointCount,
    /// No grind needed — the parameter was chosen in the open.
    NoGrindNeeded,
    /// Not a grinding hypothesis at all.
    NotApplicable,
}

/// Stated assumptions for the grind-feasibility arithmetic.  These are
/// **assumptions**, not measurements; change them and re-read the table.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GrindBudget {
    /// Total compute the generator could have spent, in CPU-years.
    pub cpu_years: f64,
    /// Seconds per candidate under [`GrindCostModel::CoefficientTest`].
    pub coefficient_test_seconds: f64,
    /// Seconds per candidate under [`GrindCostModel::PointCount`].
    pub point_count_seconds: f64,
}

impl GrindBudget {
    /// A generous late-1990s budget: a thousand CPU-years, one microsecond
    /// per SHA-1-and-test candidate, ten minutes per 256-bit SEA count.
    pub fn default_1999() -> Self {
        GrindBudget {
            cpu_years: 1000.0,
            coefficient_test_seconds: 1e-6,
            point_count_seconds: 600.0,
        }
    }

    pub fn cpu_seconds(&self) -> f64 {
        self.cpu_years * 365.25 * 86_400.0
    }

    /// `log₂` of the number of candidates the budget buys under `model`,
    /// i.e. the `k` of the thinnest weak class (density `2^−k`) that could
    /// have been hit.  `None` when no grind is involved.
    pub fn log2_candidates(&self, model: GrindCostModel) -> Option<f64> {
        let per = match model {
            GrindCostModel::CoefficientTest => self.coefficient_test_seconds,
            GrindCostModel::PointCount => self.point_count_seconds,
            GrindCostModel::NoGrindNeeded | GrindCostModel::NotApplicable => return None,
        };
        Some((self.cpu_seconds() / per).log2())
    }
}

/// One ranked hypothesis for a curve-level trapdoor.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TrapdoorHypothesis {
    /// 1 = most compatible with the public record and 1999 compute.
    pub rank: u8,
    pub locus: TrapdoorLocus,
    pub name: &'static str,
    pub mechanism: &'static str,
    pub requires: &'static str,
    pub grind: GrindCostModel,
    pub public_status: PublicResearchStatus,
    pub evidence_against: &'static str,
    /// Probe in this crate that bears on the hypothesis.
    pub repo_probe: &'static str,
    pub verdict: ProbeVerdict,
}

/// The four hypotheses, ranked.
pub fn p256_trapdoor_hypotheses() -> Vec<TrapdoorHypothesis> {
    vec![
        TrapdoorHypothesis {
            rank: 1,
            locus: TrapdoorLocus::CurveCoefficient,
            name: "Weak class recognisable from (p, a, b) without point counting",
            mechanism: "Grind seeds until SHA-1 lands b in a class whose membership is an algebraic property of the coefficients alone; the published prime-order filter (a few hundred SEA counts per curve) hides a few extra bits of selection.",
            requires: "A secret, cheaply testable invariant of (p, a, b) that implies an ECDLP weakness, unknown to the public for 25 years.",
            grind: GrindCostModel::CoefficientTest,
            public_status: PublicResearchStatus::Conjectural,
            evidence_against: "Every public coefficient-level weak class (anomalous, small embedding degree, small CM discriminant, singular) is checked and absent. Residues b mod q over 10⁵ primes pass a Kolmogorov–Smirnov uniformity test, and Solinas-reduction micro-bit correlations are null, so b shows no detectable algebraic selection.",
            repo_probe: "b_seed_profile::run_b_seed_profile; solinas_correlations::run_correlation_experiment; p256_structural checks",
            verdict: ProbeVerdict::NoPublicAnomaly,
        },
        TrapdoorHypothesis {
            rank: 2,
            locus: TrapdoorLocus::FieldPrime,
            name: "The prime, not the curve",
            mechanism: "The Solinas shape p = 2²⁵⁶ − 2²²⁴ + 2¹⁹² + 2⁹⁶ − 1 was chosen openly for fast reduction. A weakness tied to special-form primes needs no seed grinding at all — only the knowledge that such primes are weak.",
            requires: "An ECDLP analogue of the finite-field special-prime trapdoor (hidden-SNFS primes, Fried–Gaudry–Heninger–Thomé 2016). None is known; the field-DLP trapdoor does not transfer to the curve group.",
            grind: GrindCostModel::NoGrindNeeded,
            public_status: PublicResearchStatus::AdjacentOnly,
            evidence_against: "The strongest structure found in this family is cyclotomic smoothness of p ± 1 (P-224's p − 1 is 46-bit smooth). It accelerates DLP in F_p^*, which P-256's huge embedding degree never reaches; summation-polynomial index calculus measured no gain from it.",
            repo_probe: "p256_speculation::p256_solinas_prime_profile; RESEARCH_NIST_SOLINAS_STRUCTURE.md; RESEARCH_NIST_SOLINAS_EXPERIMENTS.md experiment 4",
            verdict: ProbeVerdict::StructurePresentNoKnownAttack,
        },
        TrapdoorHypothesis {
            rank: 3,
            locus: TrapdoorLocus::EndomorphismOrIsogeny,
            name: "Hidden endomorphism or isogeny structure",
            mechanism: "Grind seeds for a curve whose endomorphism ring or isogeny class admits a transfer (class-group action, descent to a weak model, small-degree isogeny to a special curve).",
            requires: "A point count per candidate (the invariant is the Frobenius trace), plus a weak class of density far above the public ones.",
            grind: GrindCostModel::PointCount,
            public_status: PublicResearchStatus::Conjectural,
            evidence_against: "Embedding degree is huge, the curve is not anomalous, the CM discriminant and class number are enormous (no CRS/CSIDH-style action to exploit), and the isogeny-class walk found no weak curve at small degree. Point-count-gated grinding also caps the reachable class density near 2^−25 under the 1999 budget.",
            repo_probe: "p256_speculation::p256_embedding_degree_probe; p256_isogeny_cover / p256_isogeny_walk; secp256k1_cm_audit method applied to generic j",
            verdict: ProbeVerdict::NoPublicAnomaly,
        },
        TrapdoorHypothesis {
            rank: 4,
            locus: TrapdoorLocus::UnknownAlgorithm,
            name: "A non-public ECDLP algorithm for generic prime-field curves",
            mechanism: "No trapdoor in the parameters at all: a classified algorithm beats generic rho on ordinary curves over prime fields.",
            requires: "A result that 25 years of open research, including the prime-regime index-calculus ladder in this repository, has not reproduced; prime-field decomposition methods remain above rho at every measured size.",
            grind: GrindCostModel::NotApplicable,
            public_status: PublicResearchStatus::Conjectural,
            evidence_against: "Not falsifiable by a finite public check; recorded as the unfalsifiable residual.",
            repo_probe: "docs/ic/BOUNDARY_TARGETS.md regime C (vs_rho: not achieved); p256_speculation EllipticIndexCalculus row",
            verdict: ProbeVerdict::NotPubliclyTestable,
        },
    ]
}

// ── Overall assessment ────────────────────────────────────────────────

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AssessmentStatus {
    /// Reasoning over public evidence; nothing here is a measurement.
    Speculation,
}

#[derive(Clone, Debug, PartialEq)]
pub struct BackdoorAssessment {
    pub status: AssessmentStatus,
    pub confirmed: DualEcRecord,
    pub seed: SeedProvenance,
    pub budget: GrindBudget,
    pub hypotheses: Vec<TrapdoorHypothesis>,
    pub overall: &'static str,
    pub binary_curve_note: &'static str,
    pub what_would_change_it: &'static str,
}

pub fn p256_backdoor_assessment() -> BackdoorAssessment {
    BackdoorAssessment {
        status: AssessmentStatus::Speculation,
        confirmed: dual_ec_drbg_record(),
        seed: seed_provenance(),
        budget: GrindBudget::default_1999(),
        hypotheses: p256_trapdoor_hypotheses(),
        overall: "The P-256 curve is almost certainly clean; the documented sabotage was Dual_EC_DRBG plus pressure on implementations and standards. The lasting damage is that the NIST process cannot prove the curve clean, which is why rigid generation (Curve25519, Brainpool, the CFRG requirements) became the expectation.",
        binary_curve_note: "The NIST binary curves all use prime extension degrees (163, 233, 283, 409, 571), which excludes the Teske/GHS magic-number trapdoor implemented in ec_trapdoor.rs (it needs a composite degree). That reads as defensive design, not an attack surface.",
        what_would_change_it: "A seed preimage (closes the question in the mundane direction); a public ECDLP method that beats rho on generic prime-field curves (opens hypothesis 4); a non-uniform statistic on b or a coefficient-level invariant with an attack behind it (opens hypothesis 1).",
    }
}

/// Markdown rendering for research notes and CLIs.
pub fn format_backdoor_assessment(a: &BackdoorAssessment) -> String {
    let mut s = String::new();
    s.push_str("# P-256 backdoor map (status: speculation with falsifiers)\n\n");
    s.push_str("## Confirmed: Dual_EC_DRBG\n\n");
    for (k, v) in [
        ("Standard", a.confirmed.standard),
        ("Mechanism", a.confirmed.mechanism),
        ("Shown by", a.confirmed.shown_by),
        ("Attribution", a.confirmed.attribution),
        ("Withdrawn", a.confirmed.withdrawn),
        ("Field incident", a.confirmed.field_incident),
    ] {
        s.push_str(&format!("- **{k}.** {v}\n"));
    }
    s.push_str(&format!(
        "- **Implicates curve math:** {}\n\n",
        if a.confirmed.implicates_curve_math {
            "yes"
        } else {
            "no"
        }
    ));
    s.push_str("## Unresolved: the seed\n\n");
    s.push_str(&format!("- **Seed:** `{}`\n", a.seed.seed_hex));
    s.push_str(&format!(
        "- **X9.62 derivation re-executed:** {}\n",
        if a.seed.derivation_verified {
            "holds (r·b² ≡ a³ mod p)"
        } else {
            "FAILS"
        }
    ));
    for (k, v) in [
        ("What verifies", a.seed.what_verifies),
        ("What does not", a.seed.what_does_not),
        ("Grinding model", a.seed.grinding_model),
        ("Counter-argument", a.seed.counter_argument),
        ("Mundane account", a.seed.mundane_account),
        ("Open bounty", a.seed.open_bounty),
    ] {
        s.push_str(&format!("- **{k}.** {v}\n"));
    }
    s.push_str("\n## Grind budget (assumptions)\n\n");
    s.push_str(&format!(
        "{} CPU-years; {:.0e} s per coefficient test; {:.0} s per 256-bit point count.\n\n",
        a.budget.cpu_years, a.budget.coefficient_test_seconds, a.budget.point_count_seconds
    ));
    s.push_str(
        "| grind model | log₂ candidates | thinnest reachable weak class |\n|---|---:|---:|\n",
    );
    for model in [GrindCostModel::CoefficientTest, GrindCostModel::PointCount] {
        if let Some(k) = a.budget.log2_candidates(model) {
            s.push_str(&format!("| {model:?} | {k:.1} | 2^−{k:.0} |\n"));
        }
    }
    s.push_str("\n## Curve-level hypotheses, ranked\n\n");
    for h in &a.hypotheses {
        s.push_str(&format!("### {}. {} ({:?})\n\n", h.rank, h.name, h.locus));
        s.push_str(&format!("- **Mechanism.** {}\n", h.mechanism));
        s.push_str(&format!("- **Requires.** {}\n", h.requires));
        s.push_str(&format!("- **Grind model.** {:?}", h.grind));
        if let Some(k) = a.budget.log2_candidates(h.grind) {
            s.push_str(&format!(" (≈2^{k:.0} candidates under the budget)"));
        }
        s.push('\n');
        s.push_str(&format!("- **Public status.** {:?}\n", h.public_status));
        s.push_str(&format!("- **Evidence against.** {}\n", h.evidence_against));
        s.push_str(&format!("- **Repo probe.** {}\n", h.repo_probe));
        s.push_str(&format!("- **Verdict.** {:?}\n\n", h.verdict));
    }
    s.push_str("## Overall\n\n");
    s.push_str(&format!(
        "{}\n\n{}\n\n**What would change this:** {}\n",
        a.overall, a.binary_curve_note, a.what_would_change_it
    ));
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn p256_seed_derivation_holds_and_matches_pinned_r() {
        let d = verify_p256_seed_derivation();
        assert_eq!((d.s, d.v), (1, 96));
        assert!(
            d.holds,
            "r·b² ≡ a³ (mod p) must hold for the published seed"
        );
        assert_eq!(format!("{:x}", d.r), P256_X962_R_HEX);
    }

    #[test]
    fn tampered_seed_fails_the_derivation() {
        let mut seed = P256_SEED;
        seed[19] ^= 1;
        assert!(!verify_seed_derivation(&CurveParams::p256(), &seed).holds);
    }

    #[test]
    fn hypotheses_are_ranked_distinct_and_claim_no_weakness() {
        let hs = p256_trapdoor_hypotheses();
        let ranks: Vec<u8> = hs.iter().map(|h| h.rank).collect();
        assert_eq!(ranks, vec![1, 2, 3, 4]);
        let loci: std::collections::HashSet<_> = hs.iter().map(|h| h.locus).collect();
        assert_eq!(loci.len(), 4);
        for h in &hs {
            assert!(!h.repo_probe.is_empty() && !h.evidence_against.is_empty());
            assert!(matches!(
                h.verdict,
                ProbeVerdict::NoPublicAnomaly
                    | ProbeVerdict::StructurePresentNoKnownAttack
                    | ProbeVerdict::NotPubliclyTestable
            ));
        }
        // The prime-level hypothesis needs no grind; the algorithmic one is not a grind.
        assert_eq!(hs[1].grind, GrindCostModel::NoGrindNeeded);
        assert_eq!(hs[3].grind, GrindCostModel::NotApplicable);
    }

    #[test]
    fn grind_budget_orders_the_models() {
        let b = GrindBudget::default_1999();
        let coeff = b.log2_candidates(GrindCostModel::CoefficientTest).unwrap();
        let count = b.log2_candidates(GrindCostModel::PointCount).unwrap();
        assert!(
            coeff > count + 20.0,
            "coefficient tests buy ≫ more candidates than point counts"
        );
        assert!(
            (24.0..28.0).contains(&count),
            "≈2^25.6 point counts: {count}"
        );
        assert!(
            (54.0..56.0).contains(&coeff),
            "≈2^54.8 coefficient tests: {coeff}"
        );
        assert_eq!(b.log2_candidates(GrindCostModel::NoGrindNeeded), None);
    }

    #[test]
    fn assessment_renders_every_record() {
        let a = p256_backdoor_assessment();
        assert_eq!(a.status, AssessmentStatus::Speculation);
        assert!(a.seed.derivation_verified);
        assert!(!a.confirmed.implicates_curve_math);
        let md = format_backdoor_assessment(&a);
        for h in &a.hypotheses {
            assert!(md.contains(h.name), "missing {}", h.name);
        }
        assert!(
            md.contains("Dual_EC_DRBG")
                && md.contains(&a.seed.seed_hex)
                && md.contains("speculation")
        );
    }
}

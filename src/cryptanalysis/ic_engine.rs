//! # The `icx` engine: run-regime classification and honest cost estimates.
//!
//! This module decides, for a catalog curve, *what kind of index-calculus run
//! is meaningful* and *what it should be expected to cost*, without ever
//! overclaiming.  It is the piece that keeps "run `icx` against every curve"
//! honest: for each curve it produces exactly one of three outcomes, per
//! `docs/ic/ENGINE_AUDIT_AND_DESIGN.md` §1.
//!
//! - [`Regime::Attack`] — index calculus is the relevant sub-exponential
//!   attack on this family, and the instance is inside the size envelope this
//!   repository has actually solved end to end.  A run recovers a verified
//!   logarithm and prices it against Pollard rho.
//! - [`Regime::Scaled`] — the named curve is out of that envelope (or is a
//!   prime-field curve, where no sub-exponential IC is known).  A run executes
//!   the same pipeline on a same-family analogue at a feasible size and labels
//!   any statement about the named curve as an extrapolation.
//! - [`Regime::EstimateOnly`] — no run is attempted; the report carries the
//!   generic-attack cost and the family's asymptotic IC form, with heuristics
//!   named.
//!
//! Nothing here claims to break a curve, and the prime-field note is explicit:
//! index calculus is **not** a threat to prime-field ECDLP, and the CLI says
//! so per curve.

use crate::cryptanalysis::curve_catalog::{CatalogCurve, Family};
use crate::cryptanalysis::matmul_exponent::{
    self, Applicability, BinaryIcHeuristic, OmegaBound, BOUNDS,
};

/// Largest subgroup-order bit length this repository has solved end to end
/// with the framework (a 48-bit subgroup of `K_0/F_{2^61}`; see
/// `docs/ic/runs/`).  Above this, a real curve is run as a scaled analogue.
pub const ATTACK_ENVELOPE_BITS: u64 = 48;

/// What sort of run is meaningful for a curve.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Regime {
    /// IC is relevant and the instance is in-envelope: solve and price it.
    Attack,
    /// Run a same-family analogue at feasible size; extrapolate, labelled.
    Scaled,
    /// Estimate only; no solve attempted.
    EstimateOnly,
}

impl Regime {
    pub fn tag(self) -> &'static str {
        match self {
            Regime::Attack => "attack",
            Regime::Scaled => "scaled",
            Regime::EstimateOnly => "estimate",
        }
    }
}

/// True when Weil-descent / summation-polynomial index calculus is the
/// relevant sub-exponential attack on this family.  Prime-field curves are
/// excluded: no sub-exponential IC is known for prime-field ECDLP.
pub fn ic_relevant(family: Family) -> bool {
    matches!(
        family,
        Family::BinaryRandom | Family::Koblitz | Family::Extension | Family::Char3
    )
}

/// The regime classification plus a human-readable rationale.
#[derive(Clone, Debug)]
pub struct RegimeReport {
    pub regime: Regime,
    pub ic_relevant: bool,
    pub rationale: String,
}

/// Classify a curve for a `run`.  `envelope_bits` defaults to
/// [`ATTACK_ENVELOPE_BITS`]; callers may lower it (e.g. for a quick CI run).
pub fn classify(curve: &CatalogCurve, envelope_bits: u64) -> RegimeReport {
    let relevant = ic_relevant(curve.family);
    let order_bits = curve.order_bits();
    if relevant && order_bits <= envelope_bits {
        RegimeReport {
            regime: Regime::Attack,
            ic_relevant: true,
            rationale: format!(
                "index calculus is the relevant sub-exponential attack on the \
                 {} family, and the {}-bit subgroup is within the end-to-end \
                 envelope ({} bits) demonstrated in this repository",
                curve.family.tag(),
                order_bits,
                envelope_bits
            ),
        }
    } else if relevant {
        RegimeReport {
            regime: Regime::Scaled,
            ic_relevant: true,
            rationale: format!(
                "index calculus applies to the {} family, but the {}-bit \
                 subgroup exceeds the {}-bit end-to-end envelope; a run uses a \
                 same-family analogue and any statement about {} is an \
                 extrapolation",
                curve.family.tag(),
                order_bits,
                envelope_bits,
                curve.name
            ),
        }
    } else {
        // Prime field: IC is not a threat.  A run is still possible on a small
        // prime analogue for study, so this is Scaled, not EstimateOnly, but
        // the rationale is explicit that IC does not beat rho here.
        RegimeReport {
            regime: Regime::Scaled,
            ic_relevant: false,
            rationale: format!(
                "no sub-exponential index calculus is known for prime-field \
                 ECDLP; Pollard rho is the best attack on {}. A run executes a \
                 small prime-field IC study instance and does not threaten the \
                 named curve",
                curve.name
            ),
        }
    }
}

/// A cost estimate for a curve, in the repository's honest terms.
#[derive(Clone, Debug)]
pub struct CostEstimate {
    pub curve: String,
    pub family: &'static str,
    pub field_label: String,
    pub field_bits: u64,
    pub order_bits: u64,
    /// `log2` of the expected Pollard-rho operation count (the generic floor).
    pub rho_log2_ops: f64,
    /// Generic-attack security level in bits (`½ log2 n`, with family speedups
    /// noted separately in `notes`).
    pub rho_security_bits: f64,
    pub ic_relevant: bool,
    pub regime: Regime,
    /// Plain-language notes: applicability, family speedups, heuristics.
    pub notes: Vec<String>,
    /// The published heuristic exponent, per `ω` bound, for a curve over
    /// `F_{2^n}`; `None` for every other field.
    pub heuristic: Option<HeuristicIc>,
}

/// The Petit–Quisquater-type heuristic for a curve over `F_{2^n}`,
/// evaluated at every `ω` bound in [`matmul_exponent::BOUNDS`].  A model
/// resting on the first-fall-degree assumption, never a measurement, and no
/// row below Strassen's is an algorithm anyone can run.
#[derive(Clone, Debug)]
pub struct HeuristicIc {
    pub heuristic: BinaryIcHeuristic,
    pub field_degree: u32,
    pub rows: Vec<HeuristicRow>,
}

/// One `ω` bound's row of [`HeuristicIc`].
#[derive(Clone, Debug)]
pub struct HeuristicRow {
    pub bound: OmegaBound,
    /// What the bound is worth in characteristic 2.
    pub applicability: Applicability,
    /// `log2` of the heuristic cost at this `n`.
    pub log2_cost: f64,
    /// The field degree from which the heuristic stays below `2^{n/2}`.
    pub turning_point: u32,
    /// The exponent's ratio to its value at the tightest bound established
    /// in characteristic 2 (it is linear in `ω`).
    pub exponent_ratio_to_best_established: f64,
}

/// The heuristic table for `E(F_{2^n})`.
pub fn binary_heuristic(n: u32) -> HeuristicIc {
    let heuristic = BinaryIcHeuristic::KousidisWiemers2019;
    let best = matmul_exponent::best_established(2).value;
    let rows = BOUNDS
        .iter()
        .map(|b| HeuristicRow {
            bound: *b,
            applicability: b.in_characteristic(2),
            log2_cost: heuristic.log2_cost(n as f64, b.value),
            turning_point: heuristic.turning_point(b.value),
            exponent_ratio_to_best_established: b.value / best,
        })
        .collect();
    HeuristicIc {
        heuristic,
        field_degree: n,
        rows,
    }
}

/// Produce a cost estimate.  This is deliberately conservative: it reports the
/// generic-attack floor exactly and describes the IC picture qualitatively
/// with its heuristics named, rather than fabricating an IC operation count
/// for a size no one has run.
pub fn estimate(curve: &CatalogCurve, envelope_bits: u64) -> CostEstimate {
    let field = curve.field();
    let order_bits = curve.order_bits();
    // log2 n from the security estimate (which is ½ log2 n).
    let log2_n = curve.rho_security_bits() * 2.0;
    // Generic rho: ~sqrt(pi/4 * n) = 0.886 * sqrt(n) group operations.
    let rho_log2_ops = 0.5 * log2_n + (std::f64::consts::PI / 4.0).sqrt().log2();
    let cls = classify(curve, envelope_bits);

    let mut notes = Vec::new();
    match curve.family {
        Family::Prime => {
            notes.push(
                "Prime-field ECDLP has no known sub-exponential index calculus; \
                 rho is the best attack. IC is included here only as a scaled \
                 study instance."
                    .to_string(),
            );
        }
        Family::Koblitz => {
            notes.push(
                "Koblitz curves admit both a rho speedup (negation + Frobenius, \
                 ~1/sqrt(2n)) and Weil-descent index calculus (Gaudry, \
                 Semaev, GGMP). In this repository IC's whole-pipeline S has \
                 not dropped below rho at any measured size; see the scoreboard."
                    .to_string(),
            );
        }
        Family::BinaryRandom => {
            notes.push(
                "Binary curves admit Weil-descent / summation-polynomial index \
                 calculus (Diem, FGHR). Sub-exponential in theory; not below \
                 rho end-to-end at the sizes measured here."
                    .to_string(),
            );
        }
        Family::Extension => {
            notes.push(
                "Extension-field curves F_{p^k}, k>=3, are the Gaudry/Diem \
                 index-calculus regime; the crossover with rho depends on k and \
                 field size."
                    .to_string(),
            );
        }
        Family::Char3 => {
            notes.push(
                "Characteristic-three curves are covered for parameter \
                 completeness; descent-based IC applies at descent-feasible \
                 sizes only."
                    .to_string(),
            );
        }
    }
    notes.push(format!(
        "Run regime for this size: {} — {}",
        cls.regime.tag(),
        cls.rationale
    ));

    let heuristic = (field.kind == "binary").then(|| binary_heuristic(field.bits as u32));
    if let Some(h) = &heuristic {
        let at = |id: &str| {
            h.rows
                .iter()
                .find(|r| r.bound.id == id)
                .expect("bound listed")
        };
        let (strassen, best, nine_fourths) = (
            at(matmul_exponent::STRASSEN.id),
            at(matmul_exponent::ALMAN_ET_AL_2024.id),
            at(matmul_exponent::CORPUS_FAMILY_107.id),
        );
        notes.push(format!(
            "Heuristic only ({}; first-fall-degree assumption, which Kosters-Yeo \
             and Huang-Kosters-Yeo give evidence against): log2 T at n = {} is \
             {:.1} at w = log2 7 (turning point n = {}), {:.1} at w < 2.371339 \
             (published, galactic; n = {}), and {:.1} at w <= 9/4 (corpus family \
             107: existence only, holds in characteristic 2 if its unverified \
             characteristic-0 proof does; n = {}). w < 2.258 is undetermined in \
             characteristic 2. None of the sub-Strassen rows is an algorithm.",
            h.heuristic.id(),
            h.field_degree,
            strassen.log2_cost,
            strassen.turning_point,
            best.log2_cost,
            best.turning_point,
            nine_fourths.log2_cost,
            nine_fourths.turning_point,
        ));
    }

    CostEstimate {
        curve: curve.name.to_string(),
        family: curve.family.tag(),
        field_label: field.label,
        field_bits: field.bits,
        order_bits,
        rho_log2_ops,
        rho_security_bits: curve.rho_security_bits(),
        ic_relevant: cls.ic_relevant,
        regime: cls.regime,
        notes,
        heuristic,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::curve_catalog::by_name;

    #[test]
    fn prime_curves_are_not_ic_relevant() {
        let p256 = by_name("p256").unwrap();
        let r = classify(&p256, ATTACK_ENVELOPE_BITS);
        assert!(!r.ic_relevant);
        assert_eq!(r.regime, Regime::Scaled);
        let est = estimate(&p256, ATTACK_ENVELOPE_BITS);
        assert!(!est.ic_relevant);
        // ~128-bit rho security for P-256.
        assert!((120.0..132.0).contains(&est.rho_security_bits));
        // rho op count ~2^127.8.
        assert!((120.0..132.0).contains(&est.rho_log2_ops));
    }

    #[test]
    fn koblitz_is_ic_relevant_and_out_of_envelope_when_large() {
        let k163 = by_name("sect163k1").unwrap();
        let r = classify(&k163, ATTACK_ENVELOPE_BITS);
        assert!(r.ic_relevant);
        // 163-bit curve is far above the 48-bit envelope.
        assert_eq!(r.regime, Regime::Scaled);
    }

    #[test]
    fn binary_estimates_carry_the_omega_table_and_prime_ones_do_not() {
        assert!(estimate(&by_name("p256").unwrap(), ATTACK_ENVELOPE_BITS)
            .heuristic
            .is_none());
        let est = estimate(&by_name("sect283k1").unwrap(), ATTACK_ENVELOPE_BITS);
        let h = est.heuristic.expect("binary field");
        assert_eq!(h.field_degree, 283);
        assert_eq!(h.rows.len(), BOUNDS.len());
        let tags: Vec<&str> = h.rows.iter().map(|r| r.applicability.tag()).collect();
        assert_eq!(
            tags,
            [
                "established",
                "established",
                "established",
                "undetermined",
                "conditional"
            ]
        );
        // Below 2^{n/2} only at the two corpus values, and only one of those
        // is decided in characteristic 2.
        let below: Vec<&str> = h
            .rows
            .iter()
            .filter(|r| r.log2_cost < 283.0 / 2.0)
            .map(|r| r.bound.id)
            .collect();
        assert_eq!(below, ["corpus-2.258", "corpus-family-107-9/4"]);
        assert!(est.notes.iter().any(|n| n.contains("first-fall-degree")));
    }

    #[test]
    fn small_koblitz_analogue_is_attack_regime() {
        // Use a tiny synthetic envelope to exercise the Attack branch on a
        // real catalog curve: pretend the envelope is 200 bits.
        let k163 = by_name("sect163k1").unwrap();
        let r = classify(&k163, 200);
        assert_eq!(r.regime, Regime::Attack);
    }
}

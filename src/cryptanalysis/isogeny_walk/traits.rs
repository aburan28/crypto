//! Trait detection for the curves a walk finds.
//!
//! Two scopes, because an `F_p`-isogeny class shares `p`, `#E` and the
//! trace:
//!
//! - **Class audits** depend only on `(p, #E)` and are run once, on the
//!   root: the parameter audit [`crate::ecc_safety`], the structural report
//!   [`crate::cryptanalysis::p256_structural::curve_structural_report`] and
//!   the Petit–Kosters–Messeng signals [`crate::cryptanalysis::pkm_criterion`].
//!   [`class_audits`] re-runs the two that take a curve model on a sample of
//!   walked curves and records whether every verdict matches the root's.
//!   That is a check of the invariance, not an assumption of it.
//! - **Curve detectors** ([`Detector`]) read one curve's model and
//!   generator, and their values can differ between isogenous curves.
//!   Each writes one `trait_status` entry in `curves.yaml`.  Add a detector
//!   to [`default_detectors`] to have every walk record it.

use num_bigint::{BigInt, BigUint};
use rayon::prelude::*;

use super::curve::{self, Model};
use super::field::{Fe, Field};
use super::record::V;
use super::walk::StartCurve;
use crate::cryptanalysis::p256_structural::{self, CurveStructuralReport};
use crate::cryptanalysis::pkm_criterion::{self, PkmAuditReport};
use crate::ecc::curve::CurveParams;
use crate::ecc_safety::{self, AnalysisOptions, SafetyCheck, SafetyReport};

/// Trial-division bound for the structural report.
pub const STRUCTURAL_TRIAL_BOUND: u64 = 1 << 20;

/// One curve, as a detector sees it.
pub struct CurveCtx<'a> {
    pub field: &'a Field,
    pub model: &'a Model,
    pub generator: (Fe, Fe),
    pub start: &'a StartCurve,
}

/// A per-curve trait.  `detect` returns the value, and `status` its
/// `curves.schema.json` status (`proved` for a checked fact,
/// `derived_from_model` for a property of the recorded model).
pub trait Detector: Sync {
    fn name(&self) -> &'static str;
    fn status(&self) -> &'static str;
    fn detect(&self, c: &CurveCtx) -> V;
}

/// `4a³ + 27b² ≠ 0`.
pub struct NonSingular;
impl Detector for NonSingular {
    fn name(&self) -> &'static str {
        "non_singular"
    }
    fn status(&self) -> &'static str {
        "proved"
    }
    fn detect(&self, c: &CurveCtx) -> V {
        V::Bool(!c.field.is_zero(&c.model.discriminant(c.field)))
    }
}

/// The recorded generator is on the curve, not the identity, and
/// `[r]G = O`.
pub struct GeneratorValid;
impl Detector for GeneratorValid {
    fn name(&self) -> &'static str {
        "generator_valid"
    }
    fn status(&self) -> &'static str {
        "proved"
    }
    fn detect(&self, c: &CurveCtx) -> V {
        let (x, y) = &c.generator;
        let ok = c.model.on_curve(c.field, x, y)
            && curve::scalar_mul(c.field, c.model, x, y, &c.start.subgroup_order).is_none();
        V::Bool(ok)
    }
}

/// An `F_p`-isomorphic model with `a = −3` exists (the fast doubling
/// formulas NIST curves use).
pub struct AMinus3Model;
impl Detector for AMinus3Model {
    fn name(&self) -> &'static str {
        "a_minus_3_model"
    }
    fn status(&self) -> &'static str {
        "derived_from_model"
    }
    fn detect(&self, c: &CurveCtx) -> V {
        V::Bool(!curve::a_minus_3_models(c.field, c.model).is_empty())
    }
}

/// `#{x ∈ [0, 64) : x³ + ax + b is a nonzero square}` of the recorded
/// model: the prefix-residue statistic PR #1330's screen sized.
pub struct QrPrefix64;
impl Detector for QrPrefix64 {
    fn name(&self) -> &'static str {
        "qr_prefix_64"
    }
    fn status(&self) -> &'static str {
        "derived_from_model"
    }
    fn detect(&self, c: &CurveCtx) -> V {
        V::int(curve::qr_prefix_64(c.field, c.model))
    }
}

/// Bit lengths of the recorded `a` and `b`, each as the smaller of `v`
/// and `p − v` (a small signed coefficient is a special form).
pub struct CoefficientBits;
impl Detector for CoefficientBits {
    fn name(&self) -> &'static str {
        "coefficient_bits"
    }
    fn status(&self) -> &'static str {
        "derived_from_model"
    }
    fn detect(&self, c: &CurveCtx) -> V {
        let bits = |v: &Fe| {
            let x = c.field.to_big(v);
            let y = c.field.modulus() - &x;
            x.min(y).bits()
        };
        V::map(vec![
            ("a", V::int(bits(&c.model.a))),
            ("b", V::int(bits(&c.model.b))),
        ])
    }
}

/// The detectors every walk runs, in output order.
pub fn default_detectors() -> Vec<Box<dyn Detector>> {
    vec![
        Box::new(NonSingular),
        Box::new(GeneratorValid),
        Box::new(AMinus3Model),
        Box::new(QrPrefix64),
        Box::new(CoefficientBits),
    ]
}

fn safety_params(
    start: &StartCurve,
    a: &BigUint,
    b: &BigUint,
    g: (&BigUint, &BigUint),
) -> ecc_safety::CurveParams {
    ecc_safety::CurveParams {
        p: start.p.clone(),
        a: BigInt::from(a.clone()),
        b: BigInt::from(b.clone()),
        gx: g.0.clone(),
        gy: g.1.clone(),
        n: start.subgroup_order.clone(),
        h: start.cofactor.clone(),
    }
}

fn curve_params(
    start: &StartCurve,
    a: &BigUint,
    b: &BigUint,
    g: (&BigUint, &BigUint),
) -> CurveParams {
    let name: &'static str = Box::leak(start.name.clone().into_boxed_str());
    CurveParams {
        name,
        p: start.p.clone(),
        a: a.clone(),
        b: b.clone(),
        gx: g.0.clone(),
        gy: g.1.clone(),
        n: start.subgroup_order.clone(),
        h: start.cofactor.clone().try_into().unwrap_or(u32::MAX),
    }
}

fn check(c: &SafetyCheck) -> V {
    match c {
        SafetyCheck::Pass => V::map(vec![("verdict", V::s("pass"))]),
        SafetyCheck::Fail(why) => V::map(vec![
            ("verdict", V::s("fail")),
            ("detail", V::s(why.clone())),
        ]),
        SafetyCheck::Inconclusive(why) => V::map(vec![
            ("verdict", V::s("inconclusive")),
            ("detail", V::s(why.clone())),
        ]),
    }
}

fn safety_rows(r: &SafetyReport) -> Vec<(&'static str, &SafetyCheck)> {
    vec![
        ("order_size", &r.order_size),
        ("order_smoothness", &r.order_smoothness),
        ("effective_security", &r.effective_security),
        ("generator_on_curve", &r.generator_on_curve),
        ("non_singular", &r.non_singular),
        ("non_supersingular", &r.non_supersingular),
        ("non_anomalous", &r.non_anomalous),
        ("multi_target_margin", &r.multi_target_margin),
        ("no_weil_descent", &r.no_weil_descent),
        ("orbit_brittleness", &r.orbit_brittleness),
    ]
}

fn safety_record(r: &SafetyReport) -> V {
    V::map(
        safety_rows(r)
            .into_iter()
            .map(|(k, c)| (k, check(c)))
            .chain(std::iter::once(("all_pass", V::Bool(r.all_pass()))))
            .collect(),
    )
}

fn factors(fs: &[(u64, u32)]) -> V {
    V::Seq(
        fs.iter()
            .map(|(q, e)| V::Seq(vec![V::int(q), V::int(e)]))
            .collect(),
    )
}

fn score(x: f64) -> V {
    V::s(format!("{x:.6}"))
}

fn structural_record(r: &CurveStructuralReport) -> V {
    V::map(vec![
        ("cm_discriminant_abs_bits", V::int(r.cm_disc_bits)),
        ("cm_discriminant_trial_factors", factors(&r.cm_disc_factors)),
        (
            "cm_discriminant_residue_bits",
            V::int(r.cm_disc_residue.bits()),
        ),
        (
            "cm_discriminant_residue_is_prime",
            V::Bool(r.cm_disc_residue_is_prime),
        ),
        ("twist_order", V::big(&r.twist_order)),
        ("twist_trial_factors", factors(&r.twist_factors)),
        ("twist_residue_bits", V::int(r.twist_residue.bits())),
        ("twist_residue_is_prime", V::Bool(r.twist_residue_is_prime)),
        ("max_twist_leak_bits", score(r.max_twist_leak_bits)),
        (
            "extension_orders",
            V::Seq(
                r.extension_orders
                    .iter()
                    .map(|(k, o)| V::map(vec![("k", V::int(k)), ("order", V::big(o))]))
                    .collect(),
            ),
        ),
        (
            "phase10_target_parity",
            V::Seq(vec![V::int(r.target_parity.0), V::int(r.target_parity.1)]),
        ),
        ("phase10_blocked", V::Bool(r.blocked)),
        ("trial_bound", V::int(STRUCTURAL_TRIAL_BOUND)),
    ])
}

fn pkm_scores(r: &PkmAuditReport) -> Vec<(&'static str, f64)> {
    let mut v = vec![
        ("special_prime", r.special_prime.score),
        ("trace_factors", r.trace_factors.score),
        ("order_neighbourhood", r.order_neighbourhood.score),
        ("embedding_window", r.embedding_window.score),
        ("overall", r.overall_score),
    ];
    if let Some(ds) = &r.divisor_sets {
        v.push(("divisor_sets", ds.score));
    }
    v
}

fn pkm_record(r: &PkmAuditReport) -> V {
    V::map(vec![
        ("solinas_weight", V::int(r.special_prime.solinas_weight)),
        (
            "near_power_of_two",
            V::Bool(r.special_prime.near_power_of_two),
        ),
        (
            "trace_small_factors",
            factors(&r.trace_factors.small_factors),
        ),
        (
            "n_minus_1_small_factors",
            factors(&r.order_neighbourhood.small_factors),
        ),
        (
            "embedding_degree_within_cap",
            V::opt(r.embedding_window.embedding_degree.map(V::int)),
        ),
        ("embedding_cap", V::int(r.embedding_window.cap)),
        (
            "divisor_sets",
            V::s(if r.divisor_sets.is_some() {
                "run"
            } else {
                "not run: field above 20 bits"
            }),
        ),
        (
            "scores",
            V::Map(
                pkm_scores(r)
                    .into_iter()
                    .map(|(k, s)| (k.to_string(), score(s)))
                    .collect(),
            ),
        ),
    ])
}

/// The class audits on the root, and the same model-taking audits re-run
/// on `sample`, each a walked curve's `(a, b, generator)`.
pub fn class_audits(start: &StartCurve, sample: &[(BigUint, BigUint, BigUint, BigUint)]) -> V {
    let opts = AnalysisOptions::default();
    let g = (&start.gx, &start.gy);
    let safety = safety_params(start, &start.a, &start.b, g).analyse(&opts);
    let root_params = curve_params(start, &start.a, &start.b, g);
    let structural = p256_structural::curve_structural_report(&root_params, STRUCTURAL_TRIAL_BOUND);
    let pkm = pkm_criterion::audit(&root_params);
    let root_safety: Vec<String> = safety_rows(&safety)
        .iter()
        .map(|(_, c)| format!("{c:?}"))
        .collect();
    let root_pkm: Vec<String> = pkm_scores(&pkm)
        .iter()
        .map(|(_, s)| format!("{s:.6}"))
        .collect();
    let differences: Vec<V> = sample
        .par_iter()
        .filter_map(|(a, b, gx, gy)| {
            let s = safety_params(start, a, b, (gx, gy)).analyse(&opts);
            let s_rows: Vec<String> = safety_rows(&s)
                .iter()
                .map(|(_, c)| format!("{c:?}"))
                .collect();
            let k = pkm_criterion::audit(&curve_params(start, a, b, (gx, gy)));
            let k_rows: Vec<String> = pkm_scores(&k)
                .iter()
                .map(|(_, s)| format!("{s:.6}"))
                .collect();
            (s_rows != root_safety || k_rows != root_pkm).then(|| {
                V::map(vec![
                    ("a", V::big(a)),
                    ("b", V::big(b)),
                    ("safety", safety_record(&s)),
                    ("pkm", pkm_record(&k)),
                ])
            })
        })
        .collect();
    V::map(vec![
        (
            "scope",
            V::s("depend only on p and #E, so they hold for every curve in the class; run once on the root"),
        ),
        ("ecc_safety", safety_record(&safety)),
        ("structural", structural_record(&structural)),
        ("pkm", pkm_record(&pkm)),
        (
            "invariance_check",
            V::map(vec![
                ("audits", V::s("ecc_safety and pkm re-run on walked curves' own models and generators")),
                ("curves_sampled", V::int(sample.len())),
                ("identical_to_root", V::Bool(differences.is_empty())),
                ("differences", V::Seq(differences)),
            ]),
        ),
    ])
}

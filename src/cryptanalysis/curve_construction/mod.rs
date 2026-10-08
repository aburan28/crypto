//! Can a curve be built whose Jacobian carries the ECC2K-130 subgroup?
//!
//! `⟨G⟩ ⊂ E_0(F_{2^131})` is `A(F_2)` for the simple 130-dimensional trace-zero
//! part `A` of `Res_{F_{2^131}/F_2}(E_0)`; a curve `C/F_2` of genus 130…~300 with
//! `P_A | P_C` and an evaluable correspondence would put the discrete logarithm
//! below Pollard rho (`research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md`).
//! This module tests four ways of producing one:
//!
//! 1. lifting to characteristic 0 — the modular curve carrying the Weil
//!    restriction ([`modular`], with the Klein quartic as the `n = 3` case in
//!    [`klein`]);
//! 2. cyclic covers of `E_0` and of `P¹` ([`cyclic`], [`kummer`],
//!    [`stickelberger`]);
//! 3. direct search, exhaustive at toy size ([`enumerate`]);
//! 4. no curve at all — cited, not recomputed.
//!
//! Everything is derived or enumerated; index-calculus costs are model figures
//! from [`window`].  Protocol:
//! `research/notes/ecc2k130/curve_construction_20261006/PROTOCOL.md`.

pub mod cyclic;
pub mod enumerate;
pub mod gf2k;
pub mod klein;
pub mod kummer;
pub mod legacy;
pub mod lpoly;
pub mod modular;
pub mod stickelberger;
pub mod weil;
pub mod window;

use enumerate::{FamilyResult, Target};
use serde::Serialize;
use std::collections::BTreeSet;
use weil::{koblitz_trace, to_i128, trace_zero_poly, weil_restriction_poly};
use window::{l_half_bits, price, rational_model_points, Cell};

#[derive(Clone, Debug)]
pub struct RunConfig {
    pub ghs_samples: usize,
    pub seed: u64,
    pub max_prime: u64,
    pub genus5_hyperelliptic: bool,
}

impl Default for RunConfig {
    fn default() -> Self {
        RunConfig {
            ghs_samples: 400,
            seed: 20261006,
            max_prime: 500,
            genus5_hyperelliptic: true,
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct ToyTarget {
    pub name: String,
    pub n: u32,
    pub koblitz_a: u8,
    pub kind: String,
    pub dimension: usize,
    pub char_poly: Vec<i128>,
}

#[derive(Clone, Debug, Serialize)]
pub struct ToyLadder {
    pub targets: Vec<ToyTarget>,
    pub families: Vec<FamilyResult>,
    /// Every smooth plane quartic with `Jac ~ Res_{F_8/F_2}(E_1)` is the Klein
    /// quartic in other coordinates.
    pub quartic_hits_for_w3_e1: usize,
    pub quartic_hits_are_one_klein_orbit: bool,
}

#[derive(Clone, Debug, Serialize)]
pub struct Route1 {
    pub klein: klein::KleinReport,
    pub cusp_validation: Vec<modular::CuspValidation>,
    pub rows: Vec<modular::ModularRow>,
    pub weil_match_131: modular::WeilMatch,
}

#[derive(Clone, Debug, Serialize)]
pub struct Route2 {
    pub cyclic_bounds_e0: Vec<cyclic::CyclicBound>,
    /// Every base of genus ≤ 1 over F_2, at n = 131.
    pub cyclic_bounds_all_bases_131: Vec<cyclic::CyclicBound>,
    pub min_genus_all_bases_131: u64,
    pub e0_order_mod_131_over_f2d: Vec<(u32, u64)>,
    pub kummer: Vec<kummer::KummerFamily>,
    pub stickelberger: Vec<stickelberger::StickelbergerScan>,
}

#[derive(Clone, Debug, Serialize)]
pub struct Route3 {
    pub log2_hyperelliptic_models_genus_130: u32,
    pub dim_moduli_of_curves_genus_130: u32,
    pub dim_hyperelliptic_locus_genus_130: u32,
    pub note: String,
}

#[derive(Clone, Debug, Serialize)]
pub struct Route4 {
    pub verdict_cited: String,
    pub source: String,
}

#[derive(Clone, Debug, Serialize)]
pub struct CostCalibration {
    pub genus: usize,
    pub exact_model_bits: f64,
    pub l_half_bits: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct ConstructionRow {
    pub route: u8,
    pub construction: String,
    pub genus: String,
    pub log2_genus: f64,
    pub carries_a: String,
    pub cost_model: String,
    pub log2_ic_cost: Option<f64>,
    pub log2_ratio_to_rho: Option<f64>,
}

#[derive(Clone, Debug, Serialize)]
pub struct Report {
    pub schema: String,
    pub question: String,
    pub verdict: String,
    pub unit: String,
    pub protocol: String,
    pub legacy_replay: legacy::LegacyReplay,
    pub toy_ladder: ToyLadder,
    pub route1_modular: Route1,
    pub route2_cyclic: Route2,
    pub route3_search: Route3,
    pub route4_no_curve: Route4,
    pub cost_calibration: Vec<CostCalibration>,
    pub constructions: Vec<ConstructionRow>,
}

pub fn toy_targets() -> Vec<ToyTarget> {
    let mut out = Vec::new();
    for n in [3u32, 5] {
        for a in [0u8, 1] {
            let t = koblitz_trace(a);
            for (kind, p) in [
                ("trace-zero A", trace_zero_poly(t, n)),
                ("Weil restriction W", weil_restriction_poly(t, n)),
            ] {
                let short = if kind.starts_with("trace") { "A" } else { "W" };
                out.push(ToyTarget {
                    name: format!("{short}{n}(E{a})"),
                    n,
                    koblitz_a: a,
                    kind: kind.into(),
                    dimension: (p.len() - 1) / 2,
                    char_poly: to_i128(&p),
                });
            }
        }
    }
    out
}

pub fn run(cfg: &RunConfig) -> Report {
    let legacy = legacy::run(cfg.ghs_samples, cfg.seed);
    let rho = legacy.target.log2_rho_reference;

    // ---- toy ladder ---------------------------------------------------------
    let tt = toy_targets();
    let targets: Vec<Target> = tt
        .iter()
        .map(|t| Target {
            name: t.name.clone(),
            poly: t.char_poly.clone(),
        })
        .collect();
    let mut families = Vec::new();
    for g in 2..=4 {
        families.push(enumerate::hyperelliptic(g, &targets));
    }
    let (quartics, quartic_hits) = enumerate::plane_quartics(&targets);
    families.push(quartics);
    families.push(enumerate::canonical_genus4(&targets));
    if cfg.genus5_hyperelliptic {
        families.push(enumerate::hyperelliptic(5, &targets));
    }
    let w3e1 = targets
        .iter()
        .position(|t| t.name == "W3(E1)")
        .expect("target present");
    let monos = enumerate::quartic_monomials();
    let klein_form = klein::klein_form(&monos);
    let orbit: BTreeSet<u64> = enumerate::gl3_f2()
        .into_iter()
        .map(|m| enumerate::substitute_quartic(&monos, klein_form, m))
        .collect();
    let hits: BTreeSet<u64> = quartic_hits[w3e1].iter().copied().collect();
    let toy = ToyLadder {
        targets: tt,
        families,
        quartic_hits_for_w3_e1: hits.len(),
        quartic_hits_are_one_klein_orbit: hits == orbit,
    };

    // ---- route 1 ------------------------------------------------------------
    let cusp_validation = [(49u64, 7u64), (121, 11), (441, 7), (539, 11), (5929, 11)]
        .iter()
        .map(|&(n, l)| modular::validate_cusps(n, l))
        .collect();
    let mut rows = Vec::new();
    for (n, ell) in [(3u32, 7u64), (5, 11), (131, 263)] {
        for a in [1u8, 0] {
            rows.push(modular::modular_row(n, a, ell, 6000));
        }
    }
    let route1 = Route1 {
        klein: klein::run(cfg.max_prime),
        cusp_validation,
        rows: rows.clone(),
        weil_match_131: modular::weil_match(&trace_zero_poly(-1, 131)),
    };

    // ---- route 2 ------------------------------------------------------------
    let bounds: Vec<cyclic::CyclicBound> =
        [3u32, 5, 131].iter().map(|&n| cyclic::bound(n)).collect();
    let e0_mod: Vec<(u32, u64)> = (1..131u32)
        .filter(|d| 130 % d == 0)
        .map(|d| (d, cyclic::e0_order_mod(131, d)))
        .collect();
    let mut kummer_fams = Vec::new();
    for ell in [3u64, 5] {
        for r in 2..=4 {
            kummer_fams.push(kummer::family(ell, r));
        }
    }
    let all131 = cyclic::all_bases(131);
    let min_all = all131.iter().map(|b| b.min_genus).min().expect("non-empty");
    let route2 = Route2 {
        cyclic_bounds_e0: bounds.clone(),
        cyclic_bounds_all_bases_131: all131,
        min_genus_all_bases_131: min_all,
        e0_order_mod_131_over_f2d: e0_mod,
        kummer: kummer_fams,
        stickelberger: vec![stickelberger::scan(7), stickelberger::scan(7 * 263)],
    };

    // ---- costs beyond the legacy window -------------------------------------
    let exact_at = |g: usize| price(&rational_model_points(g), g, 200, "random-polynomial", rho);
    let fermat_g = 920usize;
    let cyclic_g = route2.min_genus_all_bases_131 as usize;
    let c_fermat = exact_at(fermat_g);
    let c_cyclic = exact_at(cyclic_g);
    let mut calibration: Vec<CostCalibration> = legacy
        .window
        .cells
        .iter()
        .chain([&c_fermat, &c_cyclic])
        .map(|c: &Cell| CostCalibration {
            genus: c.genus,
            exact_model_bits: c.log2_total,
            l_half_bits: l_half_bits(c.genus as f64),
        })
        .collect();
    calibration.sort_by_key(|c| c.genus);
    calibration.dedup_by_key(|c| c.genus);

    // ---- one table ------------------------------------------------------------
    let deep = route2
        .stickelberger
        .iter()
        .find(|s| s.m == 7 * 263)
        .expect("scan");
    let deep_matches: u64 = deep
        .matches_by_level
        .iter()
        .filter(|(lv, _)| **lv % 263 == 0)
        .map(|(_, c)| *c)
        .sum();
    let fermat_verdict = if deep_matches == 0 {
        "excluded: no factor of level 263 or 1841 has CM type induced from Q(sqrt -7)".to_string()
    } else {
        format!("{deep_matches} factors of level 263 or 1841 have CM type induced from Q(sqrt -7): not excluded")
    };
    let m1_genus = legacy.ghs.ecc2k130_descent_genus.clone();
    let ratio = |c: f64| Some(c - rho);
    let modular131 = rows
        .iter()
        .find(|r| r.n == 131 && r.koblitz_a == 0)
        .expect("row")
        .clone();
    let mod_g = modular131.genus_split_cusps.expect("elliptic-free level");
    let kl1 = bounds[2]
        .rows
        .iter()
        .find(|r| r.d == 1)
        .expect("d = 1")
        .min_genus;
    let kum = bounds[2]
        .rows
        .iter()
        .find(|r| r.d == 130)
        .expect("d = 130")
        .min_genus;
    let exact130 = &legacy.window.exact_genus_130_jac_a;
    let constructions = vec![
        ConstructionRow {
            route: 0,
            construction: "target: genus-130 curve over F_2 with Jac ~ A (hypothetical)".into(),
            genus: "130".into(),
            log2_genus: 130f64.log2(),
            carries_a: "by definition; existence unknown".into(),
            cost_model: "exact zeta-function model (stage model of a hypothetical curve)".into(),
            log2_ic_cost: Some(exact130.log2_total),
            log2_ratio_to_rho: ratio(exact130.log2_total),
        },
        ConstructionRow {
            route: 0,
            construction: "GHS / Hess descent to F_2, ECC2K-130 itself (magic number 1)".into(),
            genus: format!("{m1_genus} (2^(m-1) - 1, type I; 2^(m-1) = 1 in the generic formula)"),
            log2_genus: 0.0,
            carries_a: "no: the transfer is the zero map on <G>".into(),
            cost_model: "none".into(),
            log2_ic_cost: None,
            log2_ratio_to_rho: None,
        },
        ConstructionRow {
            route: 0,
            construction: "GHS / Hess descent to F_2, any curve over F_2^131 (magic 130 or 131)".into(),
            genus: "2^129 or 2^130".into(),
            log2_genus: 129.0,
            carries_a: "not established".into(),
            cost_model: "floor: bits to write one divisor".into(),
            log2_ic_cost: Some(129.0),
            log2_ratio_to_rho: ratio(129.0),
        },
        ConstructionRow {
            route: 1,
            construction: format!("modular curve X_H({}) mod 2 (CM newform of level 441 twisted by an order-131 character of conductor 263)", modular131.level),
            genus: mod_g.to_string(),
            log2_genus: (mod_g as f64).log2(),
            carries_a: "yes: A_g mod 2 has the Weil numbers of A (Eichler-Shimura)".into(),
            cost_model: "L_{2^g}(1/2, sqrt 2) extrapolation".into(),
            log2_ic_cost: Some(l_half_bits(mod_g as f64)),
            log2_ratio_to_rho: ratio(l_half_bits(mod_g as f64)),
        },
        ConstructionRow {
            route: 2,
            construction: "geometrically cyclic degree-131 cover of any curve of genus <= 1 over F_2 (best: P^1, d = 10)".into(),
            genus: format!(">= {cyclic_g}"),
            log2_genus: (cyclic_g as f64).log2(),
            carries_a: "necessary bound only (class field theory + mu_d-stability)".into(),
            cost_model: "exact zeta-function model, random-polynomial counts".into(),
            log2_ic_cost: Some(c_cyclic.log2_total),
            log2_ratio_to_rho: ratio(c_cyclic.log2_total),
        },
        ConstructionRow {
            route: 2,
            construction: "cyclic degree-131 cover of E_0 with F_2-rational group (Klein mechanism, d = 1)".into(),
            genus: format!(">= {kl1}"),
            log2_genus: (kl1 as f64).log2(),
            carries_a: "necessary bound only".into(),
            cost_model: "L_{2^g}(1/2, sqrt 2) extrapolation".into(),
            log2_ic_cost: Some(l_half_bits(kl1 as f64)),
            log2_ratio_to_rho: ratio(l_half_bits(kl1 as f64)),
        },
        ConstructionRow {
            route: 2,
            construction: "Kummer cover y^131 = f of E_0 (d = 130)".into(),
            genus: format!(">= {kum}"),
            log2_genus: (kum as f64).log2(),
            carries_a: "necessary bound only".into(),
            cost_model: "L_{2^g}(1/2, sqrt 2) extrapolation".into(),
            log2_ic_cost: Some(l_half_bits(kum as f64)),
            log2_ratio_to_rho: ratio(l_half_bits(kum as f64)),
        },
        ConstructionRow {
            route: 2,
            construction: "cyclic degree-131 cover of E_0 with two branch points".into(),
            genus: "131".into(),
            log2_genus: 131f64.log2(),
            carries_a: "excluded for every twist".into(),
            cost_model: "none".into(),
            log2_ic_cost: None,
            log2_ratio_to_rho: None,
        },
        ConstructionRow {
            route: 2,
            construction: "Fermat quotient y^1841 = x^a (1-x)^b".into(),
            genus: fermat_g.to_string(),
            log2_genus: (fermat_g as f64).log2(),
            carries_a: fermat_verdict,
            cost_model: "exact zeta-function model, random-polynomial counts".into(),
            log2_ic_cost: Some(c_fermat.log2_total),
            log2_ratio_to_rho: ratio(c_fermat.log2_total),
        },
    ];

    let verdict = format!(
        "No construction reaches the window [130, {}..{}]. Best explicit curve carrying A: the modular curve X_H({}) mod 2, genus {} (2^{:.2}); \
         geometrically cyclic degree-131 covers of any curve of genus <= 1 over F_2 need genus >= {}. Whether A is isogenous to a Jacobian of genus <= ~300 remains open.",
        legacy.window.crossover_between.0,
        legacy.window.crossover_between.1,
        modular131.level,
        mod_g,
        (mod_g as f64).log2(),
        cyclic_g
    );

    Report {
        schema: "ecc2k130_curve_construction/v1".into(),
        question: "can a curve over F_2 whose Jacobian carries the ECC2K-130 subgroup be produced, and at what genus?".into(),
        verdict,
        unit: "log2 operations (model); S = ops / sqrt(r)".into(),
        protocol: "research/notes/ecc2k130/curve_construction_20261006/PROTOCOL.md".into(),
        route3_search: Route3 {
            log2_hyperelliptic_models_genus_130: 3 * 130 + 5,
            dim_moduli_of_curves_genus_130: 3 * 130 - 3,
            dim_hyperelliptic_locus_genus_130: 2 * 130 - 1,
            note: "models y^2 + h y = f with deg h <= g+1, deg f <= 2g+2 number 2^(3g+5); the toy ladder is this search run to completion at genus <= 4".into(),
        },
        route4_no_curve: Route4 {
            verdict_cited: "ECC2K-130 itself: every route tried on the challenge family (decomposition targets, Weil-descent SAT, residual walks, m >= 4 algebra) is bounded above rho or closed; the m = 4 algebraic audit reads a cost slope of 0.985-1.06 bits per unit n against the 0.25 it would need".into(),
            source: "docs/index-calculus-scoreboard.html, best-of table row 'ECC2K-130 itself'".into(),
        },
        legacy_replay: legacy,
        toy_ladder: toy,
        route1_modular: route1,
        route2_cyclic: route2,
        cost_calibration: calibration,
        constructions,
    }
}

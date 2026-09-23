//! Local isogeny neighborhoods and evidence-backed ordinary volcano metadata.
//!
//! Graph distance from the starting vertex is not endomorphism-ring depth.
//! The surface is maximal at the chosen ell, not necessarily globally maximal.
//! Trace determines the Frobenius order and the maximum rational ell-depth;
//! it does not generally determine the curve's own endomorphism order.
//!
//! Neighbor enumeration here only finds pointwise rational kernels. For odd
//! ell it is also sampled, and can miss kernels even with rational generators.
//! It cannot find Frobenius-stable kernels without rational generators.
//! Consequently graph valency and truncated walks are not used as certificates
//! of a floor, horizontal edge, or crater size. No binary-field model is supported.

use super::cm::{cm_discriminant, CmData};
use super::velu::{velu_isogeny_2, velu_isogeny_odd, VeluIsogeny};
use super::SmallCurve;
use std::collections::{BTreeMap, BTreeSet, VecDeque};

/// One breadth-first distance bucket in the enumerated neighborhood.
#[derive(Clone, Debug, serde::Serialize)]
pub struct VolcanoLevel {
    /// Distance from the input vertex; zero does not imply the crater.
    pub bfs_distance: usize,
    /// j-invariants present at this level.
    pub j_invariants: Vec<u64>,
    /// Curves at this level (as `(p, a, b)` tuples).
    pub curves: Vec<(u64, u64, u64)>,
}

/// A bounded neighborhood, not a certified complete rational volcano.
#[derive(Clone, Debug, serde::Serialize)]
pub struct VolcanoMap {
    pub ell: u64,
    pub start_curve: (u64, u64, u64),
    pub start_cm: CmData,
    pub levels: Vec<VolcanoLevel>,
    /// Whether the start curve sits on the crater (level 0).
    pub start_on_crater: Option<bool>,
    /// Only ell=2 has a complete rational-kernel enumeration in this backend.
    pub rational_kernel_enumeration_complete: bool,
    /// Preserve enumerated self-loops and parallel kernels independently of vertices.
    pub edges: Vec<VolcanoEdge>,
}

/// A recorded kernel edge; no horizontal/vertical classification is inferred.
#[derive(Clone, Debug, serde::Serialize)]
pub struct VolcanoEdge {
    pub source_j: u64,
    pub target_j: u64,
    pub kernel_half: Vec<(u64, u64)>,
}

fn validate_degree(curve: &SmallCurve, ell: u64) {
    assert!(
        curve.p > 3,
        "Short Weierstrass volcano backend requires characteristic > 3"
    );
    assert!(
        ell >= 2 && ell != curve.p,
        "ell must be prime and different from the characteristic"
    );
    let mut divisor = 2;
    while divisor <= ell / divisor {
        assert!(ell % divisor != 0, "ell must be prime");
        divisor += 1;
    }
}

fn valuation(mut conductor: i64, ell: u64) -> u32 {
    let mut depth = 0;
    while conductor > 0 && conductor as u64 % ell == 0 {
        conductor /= ell as i64;
        depth += 1;
    }
    depth
}

fn certified_depth(cm: &CmData, ell: u64) -> Option<u32> {
    if cm.supersingular {
        None
    } else if let Some(f) = cm.endomorphism_conductor {
        Some(valuation(f, ell))
    } else if cm.frobenius_order_conductor as u64 % ell != 0 {
        // f_E divides f_pi, so this local statement is known even when f_E is not.
        Some(0)
    } else {
        None
    }
}

/// Compute the j-invariant of an elliptic curve `y² = x³ + ax + b`:
///
/// ```text
/// j = 1728 · 4a³ / (4a³ + 27 b²)
/// ```
pub fn j_invariant(curve: &SmallCurve) -> u64 {
    let p = curve.p as u128;
    let a = curve.a as u128 % p;
    let b = curve.b as u128 % p;
    let a3 = (a * a % p * a) % p;
    let b2 = (b * b) % p;
    let num = (4 * a3) % p * 1728 % p;
    let den = (4 * a3 + 27 * b2) % p;
    // For singular curves the denominator is zero; return 0.
    let den_inv = match crate::utils::mod_inverse(
        &num_bigint::BigUint::from(den as u64),
        &num_bigint::BigUint::from(p as u64),
    ) {
        Some(v) => v,
        None => return 0,
    };
    let inv_u = den_inv.iter_u64_digits().next().unwrap_or(0) as u128;
    ((num * inv_u) % p) as u64
}

/// Enumerate the `ℓ`-isogenous neighbours of `curve` via Vélu's
/// formulas for pointwise rational kernels only. Frobenius-stable kernels
/// without rational generators may be absent for odd ell.
pub fn neighbors_ell(curve: &SmallCurve, ell: u64) -> Vec<VeluIsogeny> {
    validate_degree(curve, ell);
    if ell == 2 {
        // 2-isogenies: kernel = order-2 subgroup = {O, T}, with
        // T = (x_T, 0).  Roots of x³ + a x + b = 0 over F_p.
        let mut out = Vec::new();
        for x in 0..curve.p {
            if curve.rhs(x) == 0 {
                if let Some(iso) = velu_isogeny_2(curve, x) {
                    out.push(iso);
                }
            }
        }
        out
    } else {
        velu_isogeny_odd(curve, ell)
    }
}

/// Enumerate a capped neighborhood. Buckets are BFS distances, never strata.
/// The vertex cap limits discovery, not edges between retained vertices.
/// Vertices at the depth boundary are retained but are not expanded.
/// `start_on_crater` is unknown unless separately justified by conductor data.
pub fn map_volcano(
    curve: &SmallCurve,
    ell: u64,
    max_depth: usize,
    max_vertices: usize,
) -> VolcanoMap {
    validate_degree(curve, ell);
    assert!(
        max_vertices > 0,
        "Neighborhood must have room for its input vertex"
    );
    let cm = cm_discriminant(curve);
    assert!(
        !cm.supersingular,
        "Ordinary volcano metadata does not apply to supersingular curves"
    );
    let start_on_crater = certified_depth(&cm, ell).map(|d| d == 0);
    let mut edges = Vec::new();

    // BFS by j-invariant.
    let start_j = j_invariant(curve);
    let mut levels: BTreeMap<usize, Vec<(u64, u64, u64, u64)>> = BTreeMap::new();
    levels.insert(0, vec![(start_j, curve.p, curve.a, curve.b)]);
    let mut seen = BTreeSet::<u64>::new();
    seen.insert(start_j);

    let mut frontier: VecDeque<(SmallCurve, usize)> = VecDeque::new();
    frontier.push_back((*curve, 0));
    let mut visited = 1usize;

    while let Some((c, d)) = frontier.pop_front() {
        if d >= max_depth {
            continue;
        }
        let neighbours = neighbors_ell(&c, ell);
        for iso in neighbours {
            let j = j_invariant(&iso.codomain);
            if !seen.contains(&j) && visited >= max_vertices {
                continue;
            }
            edges.push(VolcanoEdge {
                source_j: j_invariant(&c),
                target_j: j,
                kernel_half: iso.kernel_half.clone(),
            });
            if seen.insert(j) {
                visited += 1;
                levels.entry(d + 1).or_default().push((
                    j,
                    iso.codomain.p,
                    iso.codomain.a,
                    iso.codomain.b,
                ));
                frontier.push_back((iso.codomain, d + 1));
            }
        }
    }

    let level_vec: Vec<VolcanoLevel> = levels
        .into_iter()
        .map(|(level, entries)| VolcanoLevel {
            bfs_distance: level,
            j_invariants: entries.iter().map(|e| e.0).collect(),
            curves: entries.iter().map(|e| (e.1, e.2, e.3)).collect(),
        })
        .collect();

    VolcanoMap {
        ell,
        start_curve: (curve.p, curve.a, curve.b),
        start_cm: cm,
        levels: level_vec,
        start_on_crater,
        rational_kernel_enumeration_complete: ell == 2,
        edges,
    }
}

/// Maximum rational ell-volcano depth v_ell(f_pi), not a guessed walk length.
/// Returns None for supersingular curves or when the requested depth cap is too
/// small. The cap is never returned as though a traversal completed.
pub fn volcano_depth(curve: &SmallCurve, ell: u64, max_depth: usize) -> Option<usize> {
    validate_degree(curve, ell);
    let cm = cm_discriminant(curve);
    if cm.supersingular {
        return None;
    }
    let depth = valuation(cm.frobenius_order_conductor, ell) as usize;
    (depth <= max_depth).then_some(depth)
}

/// Only certify a singleton surface when its proven maximal endomorphism order
/// has class number one. Other cases remain unknown; equal traces and incomplete
/// rational-point kernel walks cannot establish horizontal cycles.
pub fn crater_size(curve: &SmallCurve, ell: u64, max_steps: usize) -> Option<usize> {
    validate_degree(curve, ell);
    let cm = cm_discriminant(curve);
    if max_steps == 0 || cm.supersingular || cm.endomorphism_conductor != Some(1) {
        return None;
    }
    matches!(
        cm.fundamental_disc,
        -3 | -4 | -7 | -8 | -11 | -19 | -43 | -67 | -163
    )
    .then_some(1)
}

/// Local ordinary volcano position. Unknown is distinct from zero/false.
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize)]
pub struct VolcanoPosition {
    pub on_crater: Option<bool>,
    /// v_ell(f_E), when certified; not BFS distance.
    pub depth: Option<u32>,
    /// v_ell(f_pi), the maximum rational depth of this component.
    pub max_depth: Option<u32>,
    pub crater_size: Option<u32>,
}

pub fn position(curve: &SmallCurve, ell: u64) -> VolcanoPosition {
    validate_degree(curve, ell);
    let cm = cm_discriminant(curve);
    let depth = certified_depth(&cm, ell);
    VolcanoPosition {
        on_crater: depth.map(|d| d == 0),
        depth,
        max_depth: (!cm.supersingular).then(|| valuation(cm.frobenius_order_conductor, ell)),
        crater_size: crater_size(curve, ell, 64).map(|n| n as u32),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::isogeny::{toy_curve_a, toy_curve_b, toy_curve_j0};

    #[test]
    fn j_invariant_in_range() {
        let j = j_invariant(&toy_curve_a());
        assert!(j < toy_curve_a().p);
    }

    #[test]
    fn j_invariant_j0_curve_is_zero() {
        // y² = x³ + 1 has j = 0 by construction (a = 0).
        let j = j_invariant(&toy_curve_j0());
        assert_eq!(j, 0);
    }

    #[test]
    fn volcano_smoke() {
        // Build a tiny volcano and just check it doesn't crash.
        let v = map_volcano(&toy_curve_b(), 2, 3, 16);
        assert!(!v.levels.is_empty());
    }

    #[test]
    fn local_depth_uses_the_curve_order_not_the_frobenius_order() {
        let pos = position(&toy_curve_j0(), 2);
        assert_eq!(pos.depth, Some(0));
        assert_eq!(pos.max_depth, Some(1));
        assert_eq!(pos.on_crater, Some(true));
        assert_eq!(pos.crater_size, Some(1));
        assert_eq!(volcano_depth(&toy_curve_j0(), 2, 0), None);
        assert_eq!(volcano_depth(&toy_curve_j0(), 2, 1), Some(1));
    }

    #[test]
    fn unknown_depth_and_crater_are_not_zero_or_one() {
        let curve = SmallCurve {
            name: "unknown",
            p: 103,
            a: 88,
            b: 22,
        };
        let pos = position(&curve, 2);
        assert_eq!(pos.depth, None);
        assert_eq!(pos.on_crater, None);
        assert_eq!(pos.max_depth, Some(1));
        assert_eq!(pos.crater_size, None);
        assert_eq!(crater_size(&curve, 2, 64), None);
        // Its Frobenius conductor is 2, proving local maximality at 3.
        assert_eq!(position(&curve, 3).depth, Some(0));
        let map = map_volcano(&curve, 2, 0, 1);
        assert_eq!(map.start_on_crater, None);
        assert_eq!(map.levels[0].bfs_distance, 0);
        let json = serde_json::to_value(map).unwrap();
        assert!(json["start_on_crater"].is_null());
        assert!(json["levels"][0].get("level").is_none());
    }

    #[test]
    fn self_loops_survive_vertex_deduplication() {
        let map = map_volcano(&toy_curve_j0(), 3, 1, 8);
        assert!(!map.rational_kernel_enumeration_complete);
        assert!(map.edges.iter().any(|edge| edge.source_j == edge.target_j));
    }

    #[test]
    fn parallel_kernels_survive_the_vertex_limit() {
        let map = map_volcano(&toy_curve_j0(), 2, 1, 2);
        assert_eq!(
            map.levels
                .iter()
                .map(|b| b.j_invariants.len())
                .sum::<usize>(),
            2
        );
        assert_eq!(map.edges.len(), 3);
        assert!(map
            .edges
            .iter()
            .all(|e| e.target_j == map.edges[0].target_j));
    }

    #[test]
    fn supersingular_positions_remain_unclassified() {
        let curve = SmallCurve {
            name: "ss",
            p: 7,
            a: 1,
            b: 0,
        };
        let pos = position(&curve, 2);
        assert_eq!(pos.depth, None);
        assert_eq!(pos.max_depth, None);
        assert_eq!(pos.on_crater, None);
        assert_eq!(pos.crater_size, None);
        assert_eq!(volcano_depth(&curve, 2, 8), None);
    }

    #[test]
    fn invalid_degrees_fail_before_valuation_or_enumeration() {
        for ell in [0, 1, 4, 103] {
            assert!(std::panic::catch_unwind(|| position(&toy_curve_j0(), ell)).is_err());
        }
    }
}

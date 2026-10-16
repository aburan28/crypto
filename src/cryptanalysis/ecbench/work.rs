//! Every ecbench record priced in `bitop/v1` (`crate::cryptanalysis::bitops`).
//!
//! A record's work is rebuilt from what it already carries — the curve, the
//! phases' group operations, the native counters and the solver block — so
//! every committed session reprices without re-running a solve, and a new
//! record carries the block itself.
//!
//! Three things keep the figure honest:
//!
//! - **Every counter has a rule or the record says so.**  A counter this
//!   table does not know is named in `unpriced` and makes the total a lower
//!   bound; it is never charged zero silently.
//! - **Each counter is priced once.**  Several counters repeat a phase's
//!   native value (`matrix_work` is the elimination's `row_ops`, the
//!   shared-rank `*_uncharged` totals are sums of phase natives, the rho
//!   `*_full` `natives` phase copies the top-level counters, `aux_*_gae` and
//!   `pair_table_build_adds` are already inside the phase's additions).
//!   Index-calculus records are priced from their phase natives; generic
//!   methods from their top-level counters.
//! - **Work no counter sees is listed.**  Each method family carries the
//!   work its code does with no counter at all (`known_uncounted`), from the
//!   classification of 2026-10-15.  The instruction audit
//!   (`ecbench audit-ir`) measures how large it is; a family leaves the list
//!   only when its code counts that work.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::cryptanalysis::bitops::{self, Field, Prim, Widths};
use crate::cryptanalysis::ecbench::methods::{PhaseRecord, SolverStats};
use crate::cryptanalysis::ecbench::workload::CurveFacts;

pub const WORK_SCHEMA: &str = "ecbench.work/v1";

/// One phase in bit operations.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct PhaseWork {
    pub name: String,
    /// Primitive counts, by `bitops::Prim` name.  Exact.
    pub prims: BTreeMap<String, u128>,
    /// Gates of the phase, as a decimal integer (exact past `2^53`).
    pub bitops: String,
    pub log2_bitops: f64,
}

/// A record's work in `bitop/v1`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct WorkBlock {
    pub schema: String,
    pub unit: String,
    pub widths: Widths,
    pub phases: Vec<PhaseWork>,
    /// Total gates, as a decimal integer.
    pub bitops: String,
    pub log2_bitops: f64,
    /// `bitops / √r`.
    pub bitops_per_sqrt_r: f64,
    /// The generic floor `√(πr/2A)` additions, in gates, with `A` the
    /// automorphisms available on the curve.
    pub floor_bitops: f64,
    /// `bitops / floor_bitops`: dimensionless, comparable across curves,
    /// fields and methods.
    pub ratio_to_floor: f64,
    /// Counters with no rule.  Nonempty makes `bitops` a lower bound.
    pub unpriced: Vec<String>,
    /// Work the method's code does that no counter sees.
    pub known_uncounted: Vec<String>,
    /// `unpriced` and `known_uncounted` are both empty.
    pub complete: bool,
}

/// How one counter is charged.
#[derive(Clone, Copy, Debug)]
enum Rule {
    /// Not work: a count of outcomes or sizes, work charged under another
    /// name, or a memory size (memory is reported, never priced).
    Shape,
    /// One increment is this sum of primitives.
    Prims(&'static [(Prim, u32)]),
}

use Prim::*;
use Rule::{Prims, Shape};

const HASH: Rule = Prims(&[(HashProbe, 1)]);
const ROW_OP: Rule = Prims(&[(ModrMul, 1), (ModrAdd, 1)]);

/// Native counters of the index-calculus methods (phase `native` maps).
fn ic_rule(counter: &str) -> Option<Rule> {
    Some(match counter {
        // Outcomes, sizes, parameters, flags.
        "trials"
        | "relations"
        | "s4_hits"
        | "frobfold_mismatches"
        | "lift_failures"
        | "unliftable_systems"
        | "direct_relations_skipped"
        | "repeated_targets_skipped"
        | "stride_steps"
        | "stride_restarts"
        | "stride_renewals"
        | "walk_steps"
        | "walk_restarts"
        | "walk_pool_restarts"
        | "walk_pool_offsets"
        | "walk_jumps"
        | "native_walk_calls"
        | "lane_count"
        | "lane_batch_inversion"
        | "lane_rounds"
        | "lane_prepared"
        | "lane_unused_prepared"
        | "rows"
        | "rank"
        | "columns"
        | "single_column_rows"
        | "two_column_rows"
        | "repeated_column_rows"
        | "pinned_by_repeated_row"
        | "first_target_pinned_rank"
        | "first_target_pinned_trial"
        | "first_target_pinned_relation"
        | "full_rank_reached"
        | "dependent_rows"
        | "solver_calls"
        | "solver_budget_exceeded"
        | "abscissae_scanned"
        | "abscissa_stride"
        | "derived_counts"
        | "abscissae_scanned_uncharged"
        | "point_lifts_uncharged"
        | "cofactor_projected_points"
        | "torsion_dropped"
        | "aux_curves_tried"
        | "aux_curves_singular"
        | "aux_curves_filtered_out"
        | "aux_orders_computed"
        | "aux_filter_false_positives"
        | "aux_orders_undetermined"
        | "aux_two_sylow_rejected"
        | "pair_table_representatives"
        | "pair_table_kind"
        | "pair_table_points"
        | "scan_window_len"
        | "window_div"
        | "abscissae"
        | "fiber_size"
        | "fiber_cofactor"
        | "fiber_complete"
        | "fiber_samples"
        | "fiber_gen_coset"
        | "fiber_scale"
        | "fiber_table_entries"
        | "hits"
        | "attempts"
        | "rank_trials"
        | "rank_hits"
        | "target_attempts"
        | "decompositions"
        | "accepted_relations"
        | "eliminated_full_relations"
        | "large_prime_pivots"
        | "base_points"
        | "small_columns"
        | "core_unknowns"
        | "targets_decomposed"
        | "relations_used"
        | "full_relations"
        | "la_attempts"
        | "solver_charged_trials" => Shape,
        c if c.starts_with("relations_with_") => Shape,
        // Work already inside the phase's additions and doublings.
        "base_log_columns_checked"
        | "relation_witnesses_checked"
        | "witnesses_checked"
        | "scalar_replays"
        | "walk_pool_ops"
        | "pair_table_build_adds"
        | "aux_search_gae"
        | "aux_count_gae"
        | "aux_tower_build_gae" => Shape,
        // Sub-totals of a priced counter of the same phase.
        "aux_legendre_symbols"
        | "aux_sqrt_solves"
        | "legacy_as_solves_uncharged"
        | "factor_point_inserts_uncharged" => Shape,
        // Proxies whose work is listed as uncounted (insert attempts).
        "pair_table_entries" | "aux_bsgs_table_entries" => Shape,
        // Memory, reported beside the unit.
        "lane_state_bytes_floor" | "factor_base_vec_bytes_floor_uncharged" => Shape,
        // The solver is priced from the record's solver block.
        "solver_ops" => Shape,
        // Hash probes and inserts.
        "lookups"
        | "target_guard_probes"
        | "hash_probes_uncharged"
        | "mitm_lookups_uncharged"
        | "combination_states_uncharged" => HASH,
        // y² = rhs: the right-hand side (two products, two additions) and a
        // Tonelli–Shanks root.
        "sqrt_solves" => Prims(&[(FieldSqrt, 1), (FieldMul, 2), (FieldAdd, 2)]),
        "legendre_symbols" => Prims(&[(Legendre, 1)]),
        "inversions" => Prims(&[(FieldInv, 1)]),
        "canonicalisations" => Prims(&[(Canon, 1)]),
        "frobenius_maps" => Prims(&[(FieldSqr, 2)]),
        "montgomery_reductions" => Prims(&[(FieldMul, 1)]),
        "row_ops" | "row_ops_uncharged" | "wiedemann_modr_mults" => ROW_OP,
        "lp_merge_ops_uncharged" | "merge_modr_mults" => {
            Prims(&[(ModrMul, 1), (ModrAdd, 1), (HashProbe, 1)])
        }
        // One `points_with_x`: 1/x², b·x⁻², the Artin–Schreier solve, y = x·u.
        "as_solves" => Prims(&[(FieldInv, 1), (FieldSqr, 1), (FieldMul, 2), (AsSolve, 1)]),
        "lift_field_squares_uncharged" | "min_key_field_squares_uncharged" => {
            Prims(&[(FieldSqr, 1)])
        }
        "lift_field_multiplies_uncharged" => Prims(&[(FieldMul, 1)]),
        // is_on_curve on y² + xy = x³ + ax² + b.
        "curve_equation_checks_uncharged" => Prims(&[(FieldMul, 3), (FieldSqr, 2), (FieldAdd, 3)]),
        "point_negations_uncharged" => Prims(&[(FieldAdd, 1)]),
        // A negation or a GLV automorphism, the scalar's image mod r, a probe.
        "endomorphism_maps" => Prims(&[(FieldMul, 2), (ModrMul, 1), (HashProbe, 1)]),
        // One Vélu image: an inversion, four products, three additions.
        "velu_evaluations_uncharged" => Prims(&[(FieldInv, 1), (FieldMul, 4), (FieldAdd, 3)]),
        // ic.k4_gaudry's oracle: real E(F_{p⁴}) additions it never charged.
        "oracle_table_adds_uncharged" | "oracle_probe_adds_uncharged" => Prims(&[(GroupAdd, 1)]),
        _ => return None,
    })
}

/// Counters of the generic methods (rho, BSGS, kangaroo, the claw).
fn generic_rule(method: &str, counter: &str) -> Option<Rule> {
    let strong = method.starts_with("rho.signed_frobenius_strong");
    let kangaroo = method.starts_with("kangaroo");
    Some(match counter {
        "walk_steps"
        | "walks"
        | "fruitless_walks"
        | "capped_walks"
        | "wasted_merges"
        | "table_entries"
        | "failed_collisions"
        | "fruitless_escapes"
        | "scalar_multiplications"
        | "jumps"
        | "jump_coefficient_bits"
        | "distinguished_point_bits"
        | "setup_additions"
        | "setup_doublings"
        | "start_additions"
        | "walk_operations"
        | "verification_operations"
        | "look_ahead_additions"
        | "cycle_escape_doublings"
        | "steps"
        | "cycles_length_2"
        | "cycles_length_4"
        | "cycles_other_length"
        | "walks_abandoned_in_a_cycle"
        | "walks_capped"
        | "useless_collisions"
        | "seed_orbits"
        | "table_target"
        | "base_points"
        | "table_descriptors_skipped"
        | "table_descriptors_visited"
        | "table_identity_sums"
        | "table_duplicate_classes"
        | "table_classes"
        | "queries"
        | "queries_visited"
        | "query_pairs_skipped"
        | "relations"
        | "baby_steps"
        | "giant_steps"
        | "jump_count"
        | "dp_bits"
        | "starts" => Shape,
        // Kangaroo's distinguished points are probed under lookups_uncharged.
        "distinguished_points" if kangaroo => Shape,
        "distinguished_points" => HASH,
        "table_queries_uncharged"
        | "table_inserts_uncharged"
        | "inserts_uncharged"
        | "lookups_uncharged" => HASH,
        // fold64 twice, a xor, a rotate, two multiplies, three shift-xors, a
        // reduction.
        "partition_hashes_uncharged" => Prims(&[(WordOp, 12)]),
        // x^{2^t}, y^{2^t} by linear power maps.
        "frobenius_maps_uncharged" => Prims(&[(AsSolve, 2)]),
        "canonicalisations_uncharged" if strong => {
            // canon(x), three more basis conversions, the Montgomery scalar
            // update.
            Prims(&[(Canon, 1), (AsSolve, 3), (ModrMul, 2), (ModrAdd, 1)])
        }
        "canonicalisations_uncharged" => match method {
            "claw.pair_table" => Prims(&[(Canon, 1)]),
            "rho.cubic" => Prims(&[(FieldMul, 2), (ModrMul, 3), (FieldAdd, 2)]),
            "rho.quartic" => Prims(&[(FieldMul, 1), (ModrMul, 2), (FieldAdd, 3)]),
            m if m.starts_with("rho.signed_frobenius") => {
                Prims(&[(Canon, 1), (AsSolve, 2), (FieldAdd, 1), (ModrMul, 2)])
            }
            // Negation: a field negation and the key comparison.
            _ => Prims(&[(FieldAdd, 2)]),
        },
        _ => return None,
    })
}

/// Methods whose search and verification phases carry their group work in
/// `gae` alone (`methods::rho_report` writes additions only for set-up).
fn gae_carries_group_ops(method: &str) -> bool {
    matches!(
        method,
        "rho.plain"
            | "rho.negation"
            | "rho.cubic"
            | "rho.quartic"
            | "rho.signed_frobenius"
            | "rho.signed_frobenius_full"
    )
}

fn is_ic(method: &str) -> bool {
    method.starts_with("ic.")
}

/// Work each family's code does with no counter (classification of
/// 2026-10-15; see `docs/units/BITOPS.md`).
fn known_uncounted(method: &str, solver: Option<&SolverStats>) -> Vec<String> {
    let mut out: Vec<&str> = Vec::new();
    if is_ic(method) {
        out.extend([
            "gauss_column_scans_and_pivot_inversions",
            "relation_row_construction",
            "walk_bookkeeping",
            "table_insert_attempts",
            "base_index_builds",
        ]);
        if method == "ic.k4_gaudry" {
            out.extend(["c4_modelled_algebraic_solve", "oracle_hash_and_sort_work"]);
        }
        if method == "ic.large_prime" {
            out.push("large_prime_merge_rebuilds");
        }
    } else if method.starts_with("rho.signed_frobenius_strong") {
        out.extend(["per_lane_dp_hash", "wasted_merge_inversions"]);
    } else if method.starts_with("kangaroo") {
        out.push("per_step_mix_hashes_and_scalar_updates");
    } else if method == "bsgs.negation" {
        out.push("class_key_negations");
    } else if method.starts_with("rho.") && method != "rho.frozen_reference" {
        out.push("useless_collision_inversions");
    }
    if let Some(s) = solver {
        if let Some(gap) = solver_gap(&s.op_unit) {
            out.push(gap);
        }
    }
    out.into_iter().map(str::to_string).collect()
}

/// What a solver unit's count leaves out, if anything.
fn solver_gap(op_unit: &str) -> Option<&'static str> {
    if op_unit.starts_with("word XORs (elimination only)") {
        Some("solver_f4_symbolic_preprocessing_pairs_and_m4ri_tables")
    } else if op_unit.starts_with("word XORs (elimination, specialisation") {
        Some("solver_macaulay_build_and_summation_polynomials")
    } else if op_unit.starts_with("monomial operations") {
        Some("solver_reducer_searches_criteria_and_enumeration")
    } else if op_unit.starts_with("word operations (partial") {
        Some("solver_per_point_linear_solves")
    } else if op_unit.starts_with("word XORs") || op_unit.starts_with("vector XORs") {
        None
    } else {
        Some("solver_engine_work")
    }
}

/// Primitives one solver unit stands for, if it has a price.
fn solver_prims(s: &SolverStats) -> Option<Vec<(Prim, u128)>> {
    let u = s.op_unit.as_str();
    if u.starts_with("word XORs") || u.starts_with("word operations") {
        // f4-f2 records also carry the XORs actually performed (the logical
        // count elides M4RI table reuse); charge what was done.
        let performed = s.extra.get("word_xors_performed").copied().unwrap_or(s.ops);
        let tables = s
            .extra
            .get("full_m4ri_table_word_xors")
            .copied()
            .unwrap_or(0);
        return Some(vec![(WordOp, (performed + tables) as u128)]);
    }
    if u.starts_with("monomial operations") || u == "monomial tests" {
        // A monomial over at most 64 variables is one word: a term touched
        // is a word operation.
        return Some(vec![(WordOp, s.ops as u128)]);
    }
    if let Some(lanes) = u
        .strip_prefix("vector XORs (Gray-code steps, ")
        .and_then(|r| r.split_whitespace().next())
        .and_then(|l| l.parse::<u128>().ok())
    {
        // A vector XOR of `lanes` 32-bit lanes.
        return Some(vec![(WordOp, s.ops as u128 * lanes / 2)]);
    }
    None
}

/// The field of a curve, from its facts.
pub fn widths_of(curve: &CurveFacts) -> Result<Widths, String> {
    let r_bits = 128 - curve.r.leading_zeros();
    let field = match curve.family.as_str() {
        "prime" => Field::Prime {
            bits: curve.field_bits.ok_or("prime curve without field_bits")?,
        },
        "binary" | "koblitz" => Field::Binary {
            n: curve
                .field_degree
                .ok_or("binary curve without field_degree")?,
        },
        "fp4" => Field::Fp4 {
            bits: curve.field_bits.ok_or("fp4 curve without field_bits")?,
        },
        f => return Err(format!("no bitop/v1 field for family `{f}`")),
    };
    Ok(Widths { field, r_bits })
}

fn add_prims(into: &mut BTreeMap<Prim, u128>, prims: &[(Prim, u32)], count: u128) {
    for (p, k) in prims {
        *into.entry(*p).or_default() += count * *k as u128;
    }
}

/// What a record must supply to be priced.
pub struct Priceable<'a> {
    pub method_id: &'a str,
    pub method_params: &'a BTreeMap<String, String>,
    pub curve: &'a CurveFacts,
    pub phases: &'a [PhaseRecord],
    pub counters: &'a BTreeMap<String, u64>,
    pub solver: Option<&'a SolverStats>,
}

/// Price a record in `bitop/v1`.
pub fn price(x: &Priceable) -> Result<WorkBlock, String> {
    let widths = widths_of(x.curve)?;
    let method = x.method_id;
    let ic = is_ic(method);
    let strong = method.starts_with("rho.signed_frobenius_strong");
    let mut unpriced: Vec<String> = Vec::new();
    let mut phases: Vec<(String, BTreeMap<Prim, u128>)> = Vec::new();

    for p in x.phases {
        // The `*_full` repricing appends a copy of the counters as a phase.
        if p.name == "natives" {
            continue;
        }
        let mut prims: BTreeMap<Prim, u128> = BTreeMap::new();
        if gae_carries_group_ops(method) && p.adds + p.doubles == 0 {
            let ops = p.gae.round() as u128;
            let doubles = if p.name == "search" {
                x.counters
                    .get("cycle_escape_doublings")
                    .copied()
                    .unwrap_or(0) as u128
            } else {
                0
            };
            add_prims(&mut prims, &[(GroupAdd, 1)], ops.saturating_sub(doubles));
            add_prims(&mut prims, &[(GroupDouble, 1)], doubles);
        } else {
            add_prims(&mut prims, &[(GroupAdd, 1)], p.adds as u128);
            add_prims(&mut prims, &[(GroupDouble, 1)], p.doubles as u128);
        }
        if strong {
            // 1.5·w additions per scalar multiplication counted as such,
            // rounded down.
            *prims.entry(GroupAdd).or_default() +=
                p.scalar_mults as u128 * 3 * widths.r_bits as u128 / 2;
        }
        if ic {
            for (name, &count) in p.native.iter().chain(p.work_native.iter()) {
                if count == 0 {
                    continue;
                }
                match ic_rule(name) {
                    Some(Shape) => {}
                    Some(Prims(ps)) => add_prims(&mut prims, ps, count as u128),
                    None => push_unique(&mut unpriced, name),
                }
            }
            if p.name == "relations" {
                if let Some(s) = x.solver {
                    match solver_prims(s) {
                        Some(ps) => {
                            for (prim, c) in ps {
                                *prims.entry(prim).or_default() += c;
                            }
                        }
                        None => push_unique(
                            &mut unpriced,
                            &format!("solver_{}", s.op_unit.replace(' ', "_")),
                        ),
                    }
                }
            }
        }
        if method == "ic.k4_gaudry" && p.name == "factor_base" {
            // The base's F_p products are carried as `gae − adds` at 97 a
            // group operation (`ecbench::fp4`); an F_{p⁴} product is 16.
            let fb_muls = ((p.gae - p.adds as f64) * 97.0).round().max(0.0) as u128;
            *prims.entry(FieldMul).or_default() += fb_muls / 16;
        }
        phases.push((p.name.clone(), prims));
    }

    if ic {
        // ic.k4_gaudry keeps its linear algebra and oracle in top-level
        // counters only; every other IC counter repeats a phase native or is
        // a shape figure.
        let k4 = method == "ic.k4_gaudry";
        for (name, &count) in x.counters {
            if count == 0 {
                continue;
            }
            let in_phase = x.phases.iter().any(|p| p.native.contains_key(name));
            if in_phase && !k4 {
                continue;
            }
            match (k4, name.as_str()) {
                (
                    _,
                    "targets_tried" | "relations_found" | "matrix_rows" | "matrix_rank"
                    | "matrix_dependent" | "matrix_work",
                ) => {}
                _ => match ic_rule(name).or_else(|| ic_rule(name.trim_end_matches("_uncharged"))) {
                    Some(Shape) => {}
                    Some(Prims(ps)) => {
                        let phase = match name.as_str() {
                            "wiedemann_modr_mults"
                            | "merge_modr_mults"
                            | "row_ops_uncharged"
                            | "lp_merge_ops_uncharged" => "linear_algebra",
                            _ => "relations",
                        };
                        let slot = phase_slot(&mut phases, phase);
                        add_prims(slot, ps, count as u128);
                    }
                    None => push_unique(&mut unpriced, name),
                },
            }
        }
        if k4 {
            // The algebraic solve is a parameter, c4 F_p products a trial;
            // charge it, and keep it listed as modelled, not counted.
            let c4: u128 = x
                .method_params
                .get("c4")
                .and_then(|v| v.parse::<f64>().ok())
                .map(|v| v.round() as u128)
                .unwrap_or(0);
            let trials = x
                .counters
                .get("solver_charged_trials")
                .copied()
                .unwrap_or(0) as u128;
            let slot = phase_slot(&mut phases, "relations");
            // F_p products, at the F_{p⁴} curve's base width: a sixteenth of
            // an F_{p⁴} product each.
            let fp_muls = c4 * trials;
            *slot.entry(FieldMul).or_default() += fp_muls / 16;
        }
    } else {
        let mut natives = BTreeMap::new();
        for (name, &count) in x.counters {
            if count == 0 {
                continue;
            }
            match generic_rule(method, name) {
                Some(Shape) => {}
                Some(Prims(ps)) => add_prims(&mut natives, ps, count as u128),
                None => push_unique(&mut unpriced, name),
            }
        }
        if !natives.is_empty() {
            phases.push(("natives".into(), natives));
        }
    }

    let mut out_phases = Vec::new();
    let mut total: u128 = 0;
    for (name, prims) in phases {
        let mut gates: u128 = 0;
        for (&p, &c) in &prims {
            match widths.gates(p) {
                Ok(g) => gates += g * c,
                Err(e) => push_unique(&mut unpriced, &e.to_string()),
            }
        }
        total += gates;
        out_phases.push(PhaseWork {
            name,
            prims: prims
                .iter()
                .map(|(p, c)| (p.name().to_string(), *c))
                .collect(),
            bitops: gates.to_string(),
            log2_bitops: bitops::log2(gates),
        });
    }
    let r = x.curve.r as f64;
    let add = widths.gates(GroupAdd).map_err(|e| e.to_string())? as f64;
    let floor =
        (std::f64::consts::PI * r / (2.0 * x.curve.automorphisms_available as f64)).sqrt() * add;
    let known = known_uncounted(method, x.solver);
    unpriced.sort();
    Ok(WorkBlock {
        schema: WORK_SCHEMA.into(),
        unit: bitops::UNIT.into(),
        widths,
        phases: out_phases,
        bitops: total.to_string(),
        log2_bitops: bitops::log2(total),
        bitops_per_sqrt_r: total as f64 / r.sqrt(),
        floor_bitops: floor,
        ratio_to_floor: total as f64 / floor,
        complete: unpriced.is_empty() && known.is_empty(),
        unpriced,
        known_uncounted: known,
    })
}

/// Price a committed record.
pub fn price_record(
    r: &crate::cryptanalysis::ecbench::record::Record,
) -> Result<WorkBlock, String> {
    price(&Priceable {
        method_id: &r.method.id,
        method_params: &r.method.params,
        curve: &r.workload.curve,
        phases: &r.phases,
        counters: &r.counters,
        solver: r.solver.as_ref(),
    })
}

/// The record as the bound fitter reads it in `bitop/v1`: the cost unit
/// becomes `bitop/v1`, the total and every phase are gates, `S` is gates
/// over `√r` times the gates of one affine addition on the curve (so
/// `S / √(π/2A)` keeps its meaning, the ratio to the generic floor), and a
/// record whose pricing is incomplete is a lower bound naming what it
/// misses.  Counts, identities and timings are unchanged.
pub fn bitop_view(
    r: &crate::cryptanalysis::ecbench::record::Record,
) -> Result<crate::cryptanalysis::ecbench::record::Record, String> {
    let w = match &r.work {
        Some(w) => w.clone(),
        None => price_record(r)?,
    };
    let mut v = r.clone();
    let add = w.widths.gates(Prim::GroupAdd).map_err(|e| e.to_string())? as f64;
    let total: f64 = w.bitops.parse().map_err(|_| "bitops is not an integer")?;
    let sqrt_r = (r.workload.curve.r as f64).sqrt();
    v.cost.unit = bitops::UNIT.into();
    if r.cost.total_gae.is_some() {
        v.cost.total_gae = Some(total);
        let s = total / (sqrt_r * add);
        v.cost.s = Some(s);
        v.boundaries.ratio_to_floor = Some(s / v.boundaries.floor_s);
    }
    v.cost.lower_bound = !w.complete;
    let mut missing = w.unpriced.clone();
    missing.extend(w.known_uncounted.iter().cloned());
    v.cost.unpriced = missing;
    for p in v.phases.iter_mut() {
        p.gae = w
            .phases
            .iter()
            .find(|q| q.name == p.name)
            .and_then(|q| q.bitops.parse::<f64>().ok())
            .unwrap_or(0.0);
    }
    // Native work generic methods count outside their phases is its own
    // phase in the view, so the stage table accounts for every gate.
    if let Some(n) = w.phases.iter().find(|q| q.name == "natives") {
        if !v.phases.iter().any(|p| p.name == "natives") {
            v.phases
                .push(crate::cryptanalysis::ecbench::methods::PhaseRecord {
                    name: "natives".into(),
                    gae: n.bitops.parse().unwrap_or(0.0),
                    ..Default::default()
                });
        }
    }
    v.work = Some(w);
    Ok(v)
}

/// One line of `ecbench work`: a record's identity and its work.
#[derive(Debug, Serialize)]
pub struct WorkLine {
    pub record_id: String,
    pub session_id: String,
    pub arm: String,
    pub method: String,
    pub params: BTreeMap<String, String>,
    pub curve: String,
    pub family: String,
    pub log2_r: f64,
    pub status: String,
    pub warmup: bool,
    pub total_gae: Option<f64>,
    pub work: Result<WorkBlock, String>,
}

pub fn work_line(r: &crate::cryptanalysis::ecbench::record::Record) -> WorkLine {
    WorkLine {
        record_id: r.record_id.clone(),
        session_id: r.session_id.clone(),
        arm: r.arm.clone(),
        method: r.method.id.clone(),
        params: r.method.params.clone(),
        curve: r.workload.curve.slug.clone(),
        family: r.workload.curve.family.clone(),
        log2_r: (r.workload.curve.r as f64).log2(),
        status: r.outcome.status.clone(),
        warmup: r.warmup,
        total_gae: r.cost.total_gae,
        work: price_record(r),
    }
}

fn push_unique(v: &mut Vec<String>, s: &str) {
    if !v.iter().any(|x| x == s) {
        v.push(s.to_string());
    }
}

fn phase_slot<'a>(
    phases: &'a mut Vec<(String, BTreeMap<Prim, u128>)>,
    name: &str,
) -> &'a mut BTreeMap<Prim, u128> {
    if let Some(i) = phases.iter().position(|(n, _)| n == name) {
        return &mut phases[i].1;
    }
    phases.push((name.to_string(), BTreeMap::new()));
    &mut phases.last_mut().expect("just pushed").1
}

#[cfg(test)]
mod tests {
    use super::*;

    fn curve(family: &str, bits: u32, r: u128, aut: u32) -> CurveFacts {
        CurveFacts {
            slug: "t".into(),
            icv1: "t".into(),
            family: family.into(),
            construction: String::new(),
            field_degree: (family != "prime").then_some(bits),
            field_bits: (family == "prime").then_some(bits),
            group_order: r,
            r,
            cofactor: 1,
            generator: [String::new(), String::new()],
            automorphisms_available: aut,
            registered: false,
            ec1: None,
            curve_uid: None,
        }
    }

    #[test]
    fn rho_search_work_comes_from_gae_and_doublings_are_split() {
        let c = curve("prime", 20, (1 << 19) + 21, 2);
        let phases = vec![
            PhaseRecord {
                name: "setup".into(),
                adds: 10,
                doubles: 5,
                gae: 15.0,
                ..Default::default()
            },
            PhaseRecord {
                name: "search".into(),
                gae: 100.0,
                ..Default::default()
            },
        ];
        let counters: BTreeMap<String, u64> = [
            ("cycle_escape_doublings".to_string(), 4),
            ("walks".into(), 3),
        ]
        .into_iter()
        .collect();
        let params = BTreeMap::new();
        let w = price(&Priceable {
            method_id: "rho.negation",
            method_params: &params,
            curve: &c,
            phases: &phases,
            counters: &counters,
            solver: None,
        })
        .unwrap();
        assert_eq!(w.phases[1].prims["group_add"], 96);
        assert_eq!(w.phases[1].prims["group_double"], 4);
        let ws = widths_of(&c).unwrap();
        let expect = 106 * ws.gates(GroupAdd).unwrap() + 9 * ws.gates(GroupDouble).unwrap();
        assert_eq!(w.bitops, expect.to_string());
        assert!(w.unpriced.is_empty());
    }

    #[test]
    fn an_unknown_native_is_named_never_zero() {
        let c = curve("koblitz", 23, (1 << 21) + 5, 46);
        let phases = vec![PhaseRecord {
            name: "relations".into(),
            adds: 7,
            native: [("lookups".to_string(), 2u64), ("mystery".into(), 1)]
                .into_iter()
                .collect(),
            ..Default::default()
        }];
        let params = BTreeMap::new();
        let counters = BTreeMap::new();
        let w = price(&Priceable {
            method_id: "ic.pipeline",
            method_params: &params,
            curve: &c,
            phases: &phases,
            counters: &counters,
            solver: None,
        })
        .unwrap();
        assert_eq!(w.unpriced, vec!["mystery".to_string()]);
        assert!(!w.complete);
        assert_eq!(w.phases[0].prims["hash_probe"], 2);
    }

    #[test]
    fn matrix_work_is_not_priced_twice() {
        let c = curve("prime", 24, (1 << 23) + 15, 2);
        let phases = vec![PhaseRecord {
            name: "linear_algebra".into(),
            native: [("row_ops".to_string(), 50u64)].into_iter().collect(),
            ..Default::default()
        }];
        let counters: BTreeMap<String, u64> =
            [("matrix_work".to_string(), 50)].into_iter().collect();
        let params = BTreeMap::new();
        let w = price(&Priceable {
            method_id: "ic.pipeline",
            method_params: &params,
            curve: &c,
            phases: &phases,
            counters: &counters,
            solver: None,
        })
        .unwrap();
        assert_eq!(w.phases[0].prims["modr_mul"], 50);
        assert_eq!(w.phases.len(), 1);
    }
}

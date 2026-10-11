//! lit2026 W3: audit of the repository's index-calculus (IC) and
//! precompute configurations against three generic lower bounds.
//!
//! - **Structured generic-group model** (Corrigan-Gibbs, Henzinger, Wu,
//!   ePrint 2026/384, Theorem 3.2): in a group of prime order `r`, if a
//!   `δ` fraction of labels is "constrained by ⋆" (Definition 3.1), every
//!   `T`-query algorithm has advantage
//!   `≤ δ(3T + 2) + (3T + 1)²/r + 1/r`, so constant success needs
//!   `T = Ω(min(√r, 1/δ))` group-operation queries.  ⋆ evaluations are free.
//! - **Non-adaptive preprocessing** (Dinur, Keller, Marmor, arXiv
//!   2505.00894, Theorem 1.1): success `≤ 3T²/N + 4 ln2 · S·T/N`.
//! - **Multi-instance** (Akshima et al., ePrint 2026/2253; Yun, EUROCRYPT
//!   2015): `k` generic logarithms cost `Θ(√(k·r))` in total.
//!
//! Sub-commands (all read committed records; no record is modified):
//!
//! - `table <repo-root>`: the δ / bound table for every configuration, as
//!   JSON on stdout.
//! - `fit <repo-root>`: the multi-target slope `d ln(total ops) / d ln k`
//!   on the committed Kuhn–Struik batch records
//!   (`research/ic_rho_koblitz_20260923/batch/*.json`).
//! - `measure [batches] [prime_bits]`: a fresh multi-target rho measurement
//!   at `k = 1, 2, 4, 8, 16` on a toy prime-field curve (`A = 2`) and on
//!   `K_0 / GF(2^41)` (`A = 82`), counted in group-addition equivalents,
//!   every logarithm checked against the planted one; JSON on stdout.
//!
//! Every derived number is computed here and nowhere else; the plotting
//! script only lays out what this program writes.

use crypto_lib::cryptanalysis::ecbench::conductor_gap;
use crypto_lib::cryptanalysis::ic_boundary::{
    koblitz_instance, rho_batch_with, rho_cap, signed_frobenius_rho_batch, BinaryGroup,
    CountedGroup, GroupOps, NegationClasses, RhoWalk,
};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{json, Value};
use std::fs;
use std::path::Path;

// ── The bounds ─────────────────────────────────────────────────────────

/// A recorded log2 r that differs from log2 of the exact r by more than
/// this is reported (one-decimal rounding stays inside it).
const LOG2_TOLERANCE: f64 = 0.05;

/// One-sided level for "relation rate above the counting ceiling",
/// Bonferroni-corrected for the ~100 rows that carry a measured rate.
const SIGNIFICANCE: f64 = 0.05 / 100.0;

/// `erfc(x)`, Numerical Recipes' Chebyshev fit (fractional error < 1.2e-7).
fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.5 * z);
    let ans = t
        * (-z * z - 1.26551223
            + t * (1.00002368
                + t * (0.37409196
                    + t * (0.09678418
                        + t * (-0.18628806
                            + t * (0.27886807
                                + t * (-1.13520398
                                    + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))))
            .exp();
    if x >= 0.0 {
        ans
    } else {
        2.0 - ans
    }
}

/// `ln C(f + j − 1, j)`: multisets of `j` factor-base points.
fn ln_multiset(f: f64, j: u32) -> f64 {
    let mut acc = 0.0;
    for i in 0..j {
        acc += (f + i as f64).ln() - (i as f64 + 1.0).ln();
    }
    acc
}

/// Upper bound on δ for an `m`-term decomposition oracle over a factor base
/// of `f` (signed) points, in a subgroup of order `r`, with `c` torsion
/// classes a relation may absorb:
/// `δ ≤ min(1, c · Σ_{j=1..m} C(f + j − 1, j) / r)`.
///
/// The constrained labels of Definition 3.1 are the factor base and every
/// `j`-sum of it (`j ≤ m`) — the inputs and outputs of a defined ⋆ — so
/// their number is at most the number of multisets.
pub fn delta_upper(f: f64, m: u32, r: f64, c: f64) -> f64 {
    let mut s = 0.0f64;
    for j in 1..=m {
        s += (ln_multiset(f, j) + c.ln() - r.ln()).exp();
    }
    s.min(1.0)
}

/// The asymptotic form requested by the program: `min(√(r/A), 1/δ)`.
pub fn bound_theta(r: f64, a: f64, delta: f64) -> f64 {
    (r / a).sqrt().min(1.0 / delta)
}

/// The same with each query read as `A` labels (one per automorphism image),
/// `min(√(r/A), 1/(A δ))`.  This is our extension, not a theorem of
/// 2026/384, which has no automorphisms; it is the weaker of the two.
pub fn bound_theta_aut(r: f64, a: f64, delta: f64) -> f64 {
    (r / a).sqrt().min(1.0 / (a * delta))
}

/// Theorem 3.2 with its constants, `A = 1`: the least `T` at which the
/// advantage bound `δ(3T+2) + (3T+1)²/r + 1/r` reaches `eps`.  Either
/// term must carry half of `eps − 1/r`, so
/// `T ≥ min(((ε'/2)/δ − 2)/3, (√(r ε'/2) − 1)/3)`, `ε' = ε − 1/r`.
pub fn bound_explicit(r: f64, delta: f64, eps: f64) -> f64 {
    let e = eps - 1.0 / r;
    if e <= 0.0 {
        return 0.0;
    }
    let t_delta = ((e / 2.0) / delta - 2.0) / 3.0;
    let t_sqrt = ((r * e / 2.0).sqrt() - 1.0) / 3.0;
    t_delta.min(t_sqrt).max(0.0)
}

/// Dinur–Keller–Marmor Theorem 1.1: least online `T` for success `eps` by a
/// non-adaptive algorithm with `s` advice bits in a group of prime order
/// `n`: `3T²/n + 4 ln2 · sT/n ≥ eps`.
pub fn dkm_nonadaptive_t(n: f64, s: f64, eps: f64) -> f64 {
    // 3T² + 4 ln2 s T − eps n = 0.
    let b = 4.0 * std::f64::consts::LN_2 * s;
    (-b + (b * b + 12.0 * eps * n).sqrt()) / 6.0
}

/// Corrigan-Gibbs–Kogan (adaptive, preprocessing): `S T² ≳ N`, so
/// `T ≳ √(N/S)`; constants dropped (Õ).
pub fn ck_adaptive_t(n: f64, s: f64) -> f64 {
    (n / s.max(1.0)).sqrt()
}

/// Kuhn–Struik batch law: per-target share at `k` targets relative to one.
pub fn batch_law(k: usize) -> f64 {
    let mut term = 1.0f64;
    let mut sum = 0.0f64;
    for i in 0..k {
        sum += term;
        term *= (2 * i + 1) as f64 / (2 * i + 2) as f64;
    }
    sum / k.max(1) as f64
}

// ── Reading records ────────────────────────────────────────────────────

fn load(root: &Path, rel: &str) -> Value {
    let p = root.join(rel);
    let s = fs::read_to_string(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display()));
    serde_json::from_str(&s).unwrap_or_else(|e| panic!("{}: {e}", p.display()))
}

fn load_jsonl_first(root: &Path, rel: &str) -> Value {
    let p = root.join(rel);
    let s = fs::read_to_string(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display()));
    let line = s.lines().find(|l| !l.trim().is_empty()).expect("empty jsonl");
    serde_json::from_str(line).unwrap()
}

fn num(v: &Value) -> Option<f64> {
    match v {
        Value::Number(n) => n.as_f64(),
        Value::String(s) => s.parse::<f64>().ok(),
        _ => None,
    }
}

/// `(r, cofactor)` of a registered slug, from `docs/curves/registry.json`.
fn registry_order(reg: &Value, slug: &str) -> Option<(f64, f64)> {
    let c = reg["curves"].as_array()?.iter().find(|c| c["slug"] == slug)?;
    let cur = &c["representations"][0]["curve"];
    Some((num(&cur["subgroup_order"])?, num(&cur["cofactor"]).unwrap_or(1.0)))
}

fn parse_m(variant: &str) -> Option<u32> {
    for m in 2..=5u32 {
        if variant.contains(&format!("m{m}")) {
            return Some(m);
        }
    }
    None
}

fn parse_n(name: &str) -> Option<u32> {
    // "K_0 / GF(2^41)" or "icv1-f2m41-..."
    if let Some(i) = name.find("2^") {
        let t: String = name[i + 2..].chars().take_while(|c| c.is_ascii_digit()).collect();
        return t.parse().ok();
    }
    if let Some(i) = name.find("f2m") {
        let t: String = name[i + 3..].chars().take_while(|c| c.is_ascii_digit()).collect();
        return t.parse().ok();
    }
    None
}

#[allow(clippy::too_many_arguments)]
fn row(
    id: String,
    family: &str,
    source: String,
    r: f64,
    r_source: &str,
    a: f64,
    m: u32,
    f: f64,
    torsion_classes: f64,
    delta_measured: Option<f64>,
    trials_total: Option<f64>,
    measured: f64,
    unit: &str,
    group_op_unit: bool,
    targets: u32,
    note: &str,
) -> Value {
    let d_hi = delta_upper(f, m, r, torsion_classes);
    // The bound uses the largest δ any record or count admits (the weakest
    // bound), so a flag is never an artefact of an optimistic δ.
    let d_bound = d_hi.max(delta_measured.unwrap_or(0.0));
    let th = bound_theta(r, a, d_bound);
    let th_aut = bound_theta_aut(r, a, d_bound);
    let ex = bound_explicit(r, d_bound, 0.5);
    let mut flags = Vec::new();
    if d_hi >= 1.0 {
        flags.push("sggm_vacuous_delta_1");
    }
    // Relations per trial is a mean count of representations, whose
    // expectation for a uniform point is at most the counting ceiling.  A
    // rate above it is tested as a Poisson excess over all trials.
    let mut excess = Value::Null;
    if let (Some(dm), Some(n)) = (delta_measured, trials_total) {
        if dm > d_hi {
            let expected = d_hi * n;
            let observed = dm * n;
            let z = (observed - 0.5 - expected) / expected.sqrt();
            let p = 0.5 * erfc(z / std::f64::consts::SQRT_2);
            excess = json!({"observed_relations": observed, "expected_at_ceiling": expected,
                "z": z, "p_one_sided": p});
            flags.push("measured_rate_above_counting_ceiling");
            if p < SIGNIFICANCE {
                flags.push("rate_above_ceiling_SIGNIFICANT_ACCOUNTING_CHECK");
            }
        }
    }
    if !group_op_unit {
        flags.push("unit_not_group_operations");
    }
    // One instance of a k-target batch costs at least the single bound;
    // the batch total is what is compared.
    if group_op_unit && measured < ex {
        flags.push("below_explicit_bound_ACCOUNTING_CHECK");
    }
    if group_op_unit && measured < th {
        flags.push("below_theta_bound");
    }
    // The strict count above allows every one of the c torsion shifts; if
    // sums of factor-base points land in the subgroup's coset uniformly
    // (1 in h), the expected constrained fraction is the c = 1 count.
    let d_uni = delta_upper(f, m, r, 1.0).max(delta_measured.unwrap_or(0.0));
    let th_uni = bound_theta(r, a, d_uni);
    if group_op_unit && measured < th_uni {
        flags.push("below_uniform_landing_bound");
    }
    json!({
        "id": id, "family": family, "source": source,
        "r": r, "log2_r": r.log2(), "r_source": r_source,
        "aut": a, "m": m, "factor_base_points": f, "torsion_classes": torsion_classes,
        "delta_measured": delta_measured, "trials_total": trials_total,
        "rate_above_ceiling_test": excess,
        "delta_counting_upper": d_hi, "delta_used": d_bound,
        "inv_delta": 1.0 / d_bound, "sqrt_r_over_aut": (r / a).sqrt(),
        "bound_theta": th, "bound_theta_aut": th_aut, "bound_explicit_eps_half": ex,
        "measured": measured, "measured_unit": unit, "measured_is_group_ops": group_op_unit,
        "targets": targets,
        "ratio_measured_over_theta": measured / th,
        "ratio_measured_over_theta_aut": measured / th_aut,
        "delta_uniform_landing": d_uni, "bound_theta_uniform_landing": th_uni,
        "ratio_measured_over_theta_uniform_landing": measured / th_uni,
        "flags": flags, "note": note,
    })
}

fn table(root: &Path) -> Value {
    let reg = load(root, "docs/curves/registry.json");
    let lb = load(root, "docs/ic/leaderboard.json");
    let bt = load(root, "docs/ic/boundary_targets.json");
    let mut rows = Vec::new();
    let mut checks = Vec::new();

    // 1. Leaderboard (Tables A, B, C): every recipe of every board entry.
    for e in lb["board"].as_array().unwrap() {
        let slug = e["slug"].as_str().unwrap();
        let regime = e["regime"].as_str().unwrap();
        let (r, h) = registry_order(&reg, slug).expect("slug in registry");
        let lb_log2 = num(&e["log2_r"]).unwrap();
        if (lb_log2 - r.log2()).abs() > LOG2_TOLERANCE {
            checks.push(json!({"check": "log2_r", "where": format!("docs/ic/leaderboard.json board[{slug}]"),
                "recorded": lb_log2, "from_registry_r": r.log2()}));
        }
        let targets = e["targets"].as_u64().unwrap_or(1) as u32;
        let a = match regime {
            "koblitz" | "koblitz_batch" => 2.0 * parse_n(slug).unwrap() as f64,
            _ => 2.0,
        };
        for rec in e["recipes"].as_array().unwrap() {
            let variant = rec["variant"].as_str().unwrap_or("?");
            let m = rec["m"].as_u64().unwrap() as u32;
            let f = num(&rec["base_points"]).unwrap();
            let s = num(&rec["s"]).unwrap();
            let trials = num(&rec["trials"]);
            let rels = num(&rec["relations"]);
            let dm = match (trials, rels) {
                (Some(t), Some(x)) if t > 0.0 => Some(x / t),
                _ => None,
            };
            // trials and relations are means over `repetitions` runs.
            let reps = num(&rec["repetitions"]).unwrap_or(1.0);
            let n_trials = trials.map(|t| t * reps);
            // S is per target: the k-target batch spent k·S·√r.
            let measured = s * r.sqrt() * targets as f64;
            rows.push(row(
                format!("{slug}/{variant}"),
                regime,
                format!("docs/ic/leaderboard.json board[slug={slug}].recipes[variant={variant}]"),
                r,
                "docs/curves/registry.json",
                a,
                m,
                f,
                h.max(1.0),
                dm,
                n_trials,
                measured,
                "group-addition equivalents (S·√r·k), whole process",
                true,
                targets,
                "delta_measured = relations / trials (fraction of tried points that decomposed)",
            ));
        }
    }

    // 2. Round-5 operation-counted Koblitz rows (small n), boundary_targets
    //    supplementary ledger.  r = 2^log2_r as recorded (4 decimals).
    let k5 = &bt["supplementary_charged_ledger"]["data"]["regimes"]["koblitz"]["records"]["vs_rho"]
        ["current"]["metrics"]["operation_counted_whole_process_2026_09_22_round5"]["by_instance"];
    for (i, x) in k5.as_array().unwrap().iter().enumerate() {
        let inst = x["instance"].as_str().unwrap();
        let variant = x["variant"].as_str().unwrap();
        let n = parse_n(inst).unwrap();
        let r = 2f64.powf(num(&x["log2_r"]).unwrap());
        let m = parse_m(variant).unwrap_or(2);
        let f = num(&x["signed_points"]).unwrap();
        let h = num(&x["cofactor"]).unwrap_or(1.0);
        let s = num(&x["S"]).unwrap();
        rows.push(row(
            format!("{inst}/{variant}"),
            "koblitz_round5",
            format!("docs/ic/boundary_targets.json supplementary_charged_ledger.data.regimes.koblitz.records.vs_rho.current.metrics.operation_counted_whole_process_2026_09_22_round5.by_instance[{i}]"),
            r,
            "2^log2_r (recorded to 4 decimals)",
            2.0 * n as f64,
            m,
            f,
            h.max(1.0),
            None,
            None,
            s * r.sqrt(),
            "group-addition equivalents (S·√r), whole process",
            true,
            1,
            "m parsed from the variant name",
        ));
    }

    // 3. Compact-orbit ladder n37..n83 (4-sum relations, B = 2nK points).
    //    Measured: native probes (S3 solve + root-table lookup), not group
    //    operations.
    struct Rung {
        n: u32,
        slug: &'static str,
        b_src: (&'static str, &'static str),
        target_probes: (&'static str, &'static str),
        rank_probes: (&'static str, &'static str),
        k_cols: (&'static str, &'static str),
        log2_rec: Option<(&'static str, &'static str)>,
        jsonl: bool,
    }
    let bt_path = "docs/ic/boundary_targets.json";
    let rungs = [
        Rung { n: 37, slug: "icv1-f2m37-tm534059-32aad96b",
            b_src: ("experiments/koblitz-n37-single-target-20261004/ic_summary.jsonl", "/factor_base_points"),
            target_probes: ("experiments/koblitz-n37-single-target-20261004/ic_targets.jsonl", "/probes"),
            rank_probes: ("experiments/koblitz-n37-single-target-20261004/ic_summary.jsonl", "/rank_probes_mean"),
            k_cols: ("experiments/koblitz-n37-single-target-20261004/ic_summary.jsonl", "/orbit_columns"),
            log2_rec: None, jsonl: true },
        Rung { n: 41, slug: "icv1-f2m41-tm2308219-7f48b14a",
            b_src: ("experiments/koblitz-n41-n53-single-target-20261004/runs/n41/r1/ic_summary.jsonl", "/factor_base_points"),
            target_probes: ("experiments/koblitz-n41-n53-single-target-20261004/runs/n41/r1/ic_targets.jsonl", "/probes"),
            rank_probes: ("experiments/koblitz-n41-n53-single-target-20261004/runs/n41/r1/ic_summary.jsonl", "/rank_probes_mean"),
            k_cols: ("experiments/koblitz-n41-n53-single-target-20261004/runs/n41/r1/ic_summary.jsonl", "/orbit_columns"),
            log2_rec: None, jsonl: true },
        Rung { n: 53, slug: "icv1-f2m53-tm56619371-dac20a85",
            b_src: ("experiments/koblitz-n41-n53-single-target-20261004/runs/n53/r1/ic_summary.jsonl", "/factor_base_points"),
            target_probes: ("experiments/koblitz-n41-n53-single-target-20261004/runs/n53/r1/ic_targets.jsonl", "/probes"),
            rank_probes: ("experiments/koblitz-n41-n53-single-target-20261004/runs/n53/r1/ic_summary.jsonl", "/rank_probes_mean"),
            k_cols: ("experiments/koblitz-n41-n53-single-target-20261004/runs/n53/r1/ic_summary.jsonl", "/orbit_columns"),
            log2_rec: None, jsonl: true },
        Rung { n: 61, slug: "icv1-f2m61-t158598901-ab42b6c5",
            b_src: (bt_path, "/regimes/koblitz/records/vs_rho/history/2/metrics/factor_base_points"),
            target_probes: ("experiments/koblitz-single-target-n61-20260925/claim_report_vs_rho.json", "/operation_accounting/ic_online_operations"),
            rank_probes: ("experiments/koblitz-single-target-n61-20260925/claim_report_vs_rho.json", "/operation_accounting/ic_rank_probes_mean_pooled"),
            k_cols: (bt_path, "/regimes/koblitz/records/vs_rho/history/2/metrics/orbit_columns"),
            log2_rec: Some((bt_path, "/regimes/koblitz/records/vs_rho/history/2/metrics/subgroup_log2_order")), jsonl: false },
        Rung { n: 71, slug: "icv1-f2m71-t48653080717-2cafaea3",
            b_src: (bt_path, "/regimes/koblitz/records/vs_rho/history/1/metrics/factor_base_points"),
            target_probes: ("experiments/koblitz-single-target-n71-20261002/claim_report_vs_rho.json", "/operation_accounting/ic_online_operations"),
            rank_probes: ("experiments/koblitz-single-target-n71-20261002/claim_report_vs_rho.json", "/operation_accounting/ic_rank_probes_mean_pooled"),
            k_cols: (bt_path, "/regimes/koblitz/records/vs_rho/history/1/metrics/orbit_columns"),
            log2_rec: Some((bt_path, "/regimes/koblitz/records/vs_rho/history/1/metrics/subgroup_log2_order")), jsonl: false },
        Rung { n: 73, slug: "icv1-f2m73-tm184271214331-9d25cc67",
            b_src: (bt_path, "/regimes/koblitz/records/vs_rho/history/0/metrics/factor_base_points"),
            target_probes: ("experiments/koblitz-single-target-n73-20261003/claim_report_vs_rho.json", "/operation_accounting/ic_online_operations"),
            rank_probes: ("experiments/koblitz-single-target-n73-20261003/claim_report_vs_rho.json", "/operation_accounting/ic_rank_probes_mean_pooled"),
            k_cols: (bt_path, "/regimes/koblitz/records/vs_rho/history/0/metrics/orbit_columns"),
            log2_rec: Some((bt_path, "/regimes/koblitz/records/vs_rho/history/0/metrics/subgroup_log2_order")), jsonl: false },
        Rung { n: 83, slug: "icv1-f2m83-t6151469093347-cdcc5432",
            b_src: (bt_path, "/regimes/koblitz/records/vs_rho/current/metrics/factor_base_points"),
            target_probes: ("experiments/koblitz-single-target-n83-20261006/claim_report_vs_rho.json", "/operation_accounting/ic_online_operations"),
            rank_probes: ("experiments/koblitz-single-target-n83-20261006/claim_report_vs_rho.json", "/operation_accounting/ic_rank_probes_mean_pooled"),
            k_cols: (bt_path, "/regimes/koblitz/records/vs_rho/current/metrics/orbit_columns"),
            log2_rec: Some((bt_path, "/regimes/koblitz/records/vs_rho/current/metrics/subgroup_log2_order")), jsonl: false },
    ];
    let get = |(file, ptr): (&str, &str), jsonl: bool| -> f64 {
        let v = if jsonl && file.ends_with(".jsonl") { load_jsonl_first(root, file) } else { load(root, file) };
        num(v.pointer(ptr).unwrap_or_else(|| panic!("{file} {ptr}"))).unwrap()
    };
    for g in &rungs {
        // r: the registry where the slug is registered, else the rung's
        // claim report (`subgroup_order`); the registry lacks n >= 71.
        let (r, r_source) = match registry_order(&reg, g.slug) {
            Some((r, _)) => (r, "docs/curves/registry.json".to_string()),
            None => (
                get((g.target_probes.0, "/subgroup_order"), false),
                format!("{} /subgroup_order", g.target_probes.0),
            ),
        };
        let a = 2.0 * g.n as f64;
        let b = get(g.b_src, g.jsonl);
        let tp = get(g.target_probes, g.jsonl);
        let rp = get(g.rank_probes, g.jsonl);
        let k = get(g.k_cols, g.jsonl);
        if let Some(src) = g.log2_rec {
            let rec = get(src, false);
            if (rec - r.log2()).abs() > LOG2_TOLERANCE {
                checks.push(json!({"check": "log2_r", "where": format!("{} {}", src.0, src.1),
                    "recorded": rec, "from_registry_r": r.log2(), "slug": g.slug}));
            }
        }
        // Excluded precompute in the ledger against its own components.
        if let Some(src) = g.log2_rec {
            let ledger = load(root, src.0);
            let ptr = src.1.replace("/subgroup_log2_order", "/ic_precomputation_ms_excluded_from_online");
            let claim = load(root, g.target_probes.0);
            if let (Some(ex), Some(d)) = (
                ledger.pointer(&ptr).and_then(num),
                claim.get("ic_precompute_detail").or_else(|| claim.get("precompute_detail")),
            ) {
                let part = |k: &str| d.get(k).and_then(num);
                let rank_max = d
                    .get("rank_stage_ms_range")
                    .and_then(|v| v.as_array())
                    .map(|a| a.iter().filter_map(num).fold(0.0, f64::max));
                let rank = part("rank_stage_ms").or(rank_max);
                let index = part("index_build_ms");
                let comps = [
                    ("base_load", part("base_load_ms")),
                    ("index_build", index),
                    ("rank_stage", rank),
                    ("linear_algebra", part("linear_algebra_ms")),
                ];
                let sum_known: f64 = comps.iter().filter_map(|c| c.1).sum();
                // The subset of recorded components whose sum is closest to
                // the recorded figure says what the figure contains.
                let known: Vec<(&str, f64)> = comps.iter().filter_map(|(k, v)| v.map(|v| (*k, v))).collect();
                let mut best = (f64::INFINITY, 0usize);
                for mask in 1..(1usize << known.len()) {
                    let s: f64 = (0..known.len()).filter(|i| mask >> i & 1 == 1).map(|i| known[i].1).sum();
                    let rel = (s - ex).abs() / ex;
                    if rel < best.0 {
                        best = (rel, mask);
                    }
                }
                let matched: Vec<&str> = (0..known.len()).filter(|i| best.1 >> i & 1 == 1).map(|i| known[i].0).collect();
                let includes_rank = matched.contains(&"rank_stage");
                let includes_index = index.map(|_| matched.contains(&"index_build"));
                if !includes_rank || includes_index == Some(false) {
                    checks.push(json!({"check": "excluded_precompute", "where": format!("{} {ptr}", src.0),
                        "slug": g.slug, "n": g.n, "recorded_ms": ex,
                        "components_ms": {"base_load": part("base_load_ms"), "index_build": index, "rank_stage": rank,
                            "linear_algebra": part("linear_algebra_ms")},
                        "components_source": format!("{} ic_precompute_detail|precompute_detail", g.target_probes.0),
                        "includes_rank_stage": includes_rank, "includes_index_build": includes_index,
                        "sum_of_recorded_components_ms": sum_known,
                        "closest_component_subset": matched, "closest_subset_relative_error": best.0}));
                }
            }
        }
        let total_probes = rp * k + tp;
        for (what, val) in [("target_probes_online", tp), ("rank_plus_target_probes", total_probes)] {
            let mut v = row(
                format!("compact-n{}/{what}", g.n),
                "koblitz_compact_orbit",
                format!("{} {} ; {} {}", g.target_probes.0, g.target_probes.1, g.rank_probes.0, g.rank_probes.1),
                r,
                &r_source,
                a,
                4,
                b,
                1.0,
                None,
                None,
                val,
                "probes (S3 quadratic solve + root-table lookup); not group operations",
                false,
                1,
                "4-sum relations over B = 2nK signed-orbit points; online = target probes with the index and logs precomputed",
            );
            v["orbit_columns"] = json!(k);
            v["probes_over_sqrt_r_over_aut"] = json!(val / (r / a).sqrt());
            // Advice: the S3 root index holds the K²n regular states.  The
            // adaptive preprocessing line (Corrigan-Gibbs–Kogan, also the
            // δ-free part of 2026/384 Thm B.1) is S·T² ≈ r/A.
            let s_entries = k * k * g.n as f64;
            v["advice_entries_k2n"] = json!(s_entries);
            v["ck_adaptive_t_floor"] = json!(ck_adaptive_t(r / a, s_entries));
            v["log2_s_t2_minus_log2_r_over_a"] = json!((s_entries * val * val).log2() - (r / a).log2());
            // Non-adaptive line (Dinur–Keller–Marmor, in entries): S·T vs r/A.
            v["s_t_over_r_over_a"] = json!(s_entries * val / (r / a));
            if what == "target_probes_online" {
                // The frozen target is one draw.  If first-hit probes were
                // exponential with the rank-stage mean, P(T <= observed):
                v["rank_probes_mean_per_relation"] = json!(rp);
                v["target_over_rank_mean"] = json!(tp / rp);
                v["p_exponential_at_most_observed"] = json!(1.0 - (-tp / rp).exp());
                v["s_t_over_r_over_a_at_rank_mean"] = json!(s_entries * rp / (r / a));
            }
            rows.push(v);
        }
    }

    // 4. Multiplier closure (Lead A1), research/transfer_aware_ic_20261010:
    //    a pair table of x(s_a + s_b) (advice) and online probes
    //    T = R − s_c − s_d with (c, d) drawn independently of R, stopped at
    //    the first table hit: a non-adaptive (S, T) algorithm in the sense
    //    of Dinur–Keller–Marmor.  Their Theorem 1.1 puts it on S·T ≳ N.
    let mut pre = Vec::new();
    for cell in ["c23", "c37", "c37b"] {
        let rel = format!("research/transfer_aware_ic_20261010/results/{cell}.json");
        let v = load(root, &rel);
        let r = num(&v["r"]).unwrap();
        let n = num(&v["n"]).unwrap();
        let a = 2.0 * n;
        for (i, arm) in v["arms"].as_array().unwrap().iter().enumerate() {
            let t = num(&arm["mean_first_hit_probes"]).unwrap();
            let pairs = num(&arm["pairs"]).unwrap();
            let canon = num(&arm["canonical"]).unwrap();
            let bits = r.log2().ceil();
            pre.push(json!({
                "id": format!("multiplier-closure/{cell}/{}-K{}-c{}", arm["kind"].as_str().unwrap_or("?"), arm["K"], arm["c"]),
                "source": format!("{rel} /arms/{i}"),
                "r": r, "n": n, "aut": a, "hits": arm["hits"],
                "pairs": pairs, "canonical_entries": canon, "mean_first_hit_probes": t,
                "pairs_times_t_over_r": pairs * t / r,
                "canonical_times_t_over_r_over_a": canon * t / (r / a),
                "dkm_nonadaptive_t_eps_half_canonical_bits": dkm_nonadaptive_t(r / a, canon * bits, 0.5),
                "dkm_nonadaptive_t_eps_half_pairs_bits": dkm_nonadaptive_t(r, pairs * bits, 0.5),
                "sqrt_r_over_a": (r / a).sqrt(),
                "logs_ok": arm["logs_ok"], "target_ok": arm["target_ok"],
            }));
        }
    }
    json!({
        "precompute_multiplier_closure": pre,
        "schema": "lit2026-w3-ggm-table/v1",
        "bounds": {
            "sggm_theorem_3_2": "Adv <= delta(3T+2) + (3T+1)^2/r + 1/r  (ePrint 2026/384, Thm 3.2, prime order r, T group-operation queries, star evaluations free)",
            "sggm_theorem_B_1": "Adv <= O~(S T^2/q + delta T) with S bits of advice",
            "dkm_theorem_1_1": "non-adaptive (S,T): Adv <= 3T^2/N + 4 ln2 S T/N  (arXiv 2505.00894)",
            "delta_definition": "fraction of subgroup labels constrained by star: the factor base and every j-sum (j<=m) of it; upper-bounded by c*sum_j C(|F|+j-1,j)/r with c torsion classes (cofactor h)",
            "delta_used": "max(counting upper bound, measured relations/trials): the weakest bound any record admits",
        },
        "rows": rows,
        "record_checks": checks,
    })
}

// ── Multi-target slope ─────────────────────────────────────────────────

/// OLS slope of `y` on `x`.
fn ols(x: &[f64], y: &[f64]) -> (f64, f64) {
    let n = x.len() as f64;
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let sxy: f64 = x.iter().zip(y).map(|(a, b)| (a - mx) * (b - my)).sum();
    let sxx: f64 = x.iter().map(|a| (a - mx).powi(2)).sum();
    let slope = sxy / sxx;
    (slope, my - slope * mx)
}

/// Slope of ln(mean batch total) on ln k, with a bootstrap over batches
/// within each k (2,000 resamples, fixed seed): `(slope, lo, hi)`.
fn slope_with_ci(by_k: &[(usize, Vec<f64>)], seed: u64) -> (f64, f64, f64) {
    let point = |sets: &[(usize, Vec<f64>)]| -> f64 {
        let x: Vec<f64> = sets.iter().map(|(k, _)| (*k as f64).ln()).collect();
        let y: Vec<f64> = sets
            .iter()
            .map(|(_, v)| (v.iter().sum::<f64>() / v.len() as f64).ln())
            .collect();
        ols(&x, &y).0
    };
    let s = point(by_k);
    let mut rng = StdRng::seed_from_u64(seed);
    let mut boots = Vec::with_capacity(2000);
    for _ in 0..2000 {
        let res: Vec<(usize, Vec<f64>)> = by_k
            .iter()
            .map(|(k, v)| (*k, (0..v.len()).map(|_| v[rng.gen_range(0..v.len())]).collect()))
            .collect();
        boots.push(point(&res));
    }
    boots.sort_by(|a, b| a.partial_cmp(b).unwrap());
    (s, boots[50], boots[1949])
}

fn law_slope(ks: &[usize]) -> f64 {
    let x: Vec<f64> = ks.iter().map(|&k| (k as f64).ln()).collect();
    let y: Vec<f64> = ks.iter().map(|&k| (k as f64 * batch_law(k)).ln()).collect();
    ols(&x, &y).0
}

fn summarize(curve: &str, r: f64, a: f64, by_k: Vec<(usize, Vec<f64>)>, seed: u64, source: &str) -> Value {
    let ks: Vec<usize> = by_k.iter().map(|(k, _)| *k).collect();
    let (s, lo, hi) = slope_with_ci(&by_k, seed);
    let per_k: Vec<Value> = by_k
        .iter()
        .map(|(k, v)| {
            let mean = v.iter().sum::<f64>() / v.len() as f64;
            let generic = ((*k as f64) * r / a).sqrt();
            json!({"k": k, "batches": v.len(), "mean_total_gae": mean,
                "mean_total_over_sqrt_k_r_over_A": mean / generic,
                "per_target_over_k1": None::<f64>,
                "batch_law": batch_law(*k)})
        })
        .collect();
    let k1 = by_k.iter().find(|(k, _)| *k == 1).map(|(_, v)| v.iter().sum::<f64>() / v.len() as f64);
    let per_k: Vec<Value> = per_k
        .into_iter()
        .map(|mut p| {
            if let Some(m1) = k1 {
                let k = p["k"].as_u64().unwrap() as f64;
                p["per_target_over_k1"] = json!(p["mean_total_gae"].as_f64().unwrap() / k / m1);
            }
            p
        })
        .collect();
    json!({"curve": curve, "r": r, "log2_r": r.log2(), "aut": a, "source": source,
        "ks": ks, "slope": s, "slope_ci95_bootstrap": [lo, hi],
        "kuhn_struik_law_slope_same_ks": law_slope(&ks), "asymptotic_slope": 0.5,
        "per_k": per_k})
}

fn fit(root: &Path) -> Value {
    let mut out = Vec::new();
    for f in ["k0n41", "k0n53", "k0n61"] {
        let rel = format!("research/ic_rho_koblitz_20260923/batch/{f}.json");
        let v = load(root, &rel);
        for c in v["curves"].as_array().unwrap() {
            let r = num(&c["r"]).unwrap();
            let a = num(&c["automorphisms"]).unwrap();
            let by_k: Vec<(usize, Vec<f64>)> = c["sizes"]
                .as_array()
                .unwrap()
                .iter()
                .map(|s| {
                    let k = s["k"].as_u64().unwrap() as usize;
                    let g: Vec<f64> = s["runs"].as_array().unwrap().iter().map(|x| num(&x["gae"]).unwrap()).collect();
                    (k, g)
                })
                .collect();
            out.push(summarize(c["instance"].as_str().unwrap(), r, a, by_k, 0x5733, &rel));
        }
    }
    json!({"schema": "lit2026-w3-multitarget-fit/v1", "unit": "group-addition equivalents (gae) per batch, whole batch",
        "fits": out})
}

// ── Fresh multi-target measurement ─────────────────────────────────────

fn derive(seed: u64, parts: &[u64]) -> u64 {
    let mut h = seed ^ 0x9E37_79B9_7F4A_7C15;
    for &p in parts {
        h = (h ^ p).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        h ^= h >> 31;
    }
    h
}

/// A prime-order curve over `GF(q)`, `q` about `bits` bits: the first
/// `(ℓ, t)` that `conductor_gap::params` accepts, then a random (`floor`)
/// curve of that prime order — ordinary `j`, `Aut = {±1}`.
fn toy_prime_instance(bits: u32) -> conductor_gap::GapParams {
    // q = (t² + 19ℓ²)/4 ≈ 2^bits: ℓ ≈ √(4·2^bits/19).
    let ell0 = ((4.0 * 2f64.powi(bits as i32) / 19.0).sqrt()) as u64 | 1;
    for ell in (ell0..ell0 + 200_000).step_by(2) {
        for t in (1..4000u64).step_by(2) {
            if let Ok(p) = conductor_gap::params(ell, t) {
                return p;
            }
        }
    }
    panic!("no parameters near {bits} bits");
}

fn measure(batches: usize, prime_bits: u32) -> Value {
    let ks = [1usize, 2, 4, 8, 16];
    let seed = 0x1A2026_0303u64;
    let mut curves = Vec::new();

    // (a) prime field, A = 2.
    let gp = toy_prime_instance(prime_bits);
    let inst = conductor_gap::floor(gp, seed).expect("floor curve");
    let g = inst.curve.clone();
    let gen = inst.generator_point();
    let r = inst.r;
    let mut by_k = Vec::new();
    let mut runs_out = Vec::new();
    let mut all_ok = true;
    for &k in &ks {
        let mut tot = Vec::new();
        for b in 0..batches {
            let mut rng = StdRng::seed_from_u64(derive(seed, &[1, k as u64, b as u64]));
            let planted: Vec<u64> = (0..k).map(|_| rng.gen_range(1..r)).collect();
            let mut ops = GroupOps::default();
            let targets: Vec<_> = planted.iter().map(|&d| g.mul(&mut ops, gen, d)).collect();
            let res = rho_batch_with(&g, &NegationClasses { r }, gen, &targets, r, derive(seed, &[2, k as u64, b as u64]),
                rho_cap(r, 64.0), RhoWalk::negation());
            let ok = res.per_target.iter().zip(&planted).all(|(t, &d)| t.verified && t.recovered == Some(d));
            all_ok &= ok;
            tot.push(res.gae);
            runs_out.push(json!({"k": k, "batch": b, "gae": res.gae, "all_verified": ok,
                "solved_on_an_earlier_trail": res.counters.get("solved_on_an_earlier_trail")}));
        }
        by_k.push((k, tot));
    }
    let mut s = summarize(&format!("prime-floor q={} N={} (ell={}, t={})", gp.q, gp.n, gp.ell, gp.t), r as f64, 2.0, by_k, 0x5734, "fresh: rho_batch_with, NegationClasses, RhoWalk::negation");
    s["all_verified"] = json!(all_ok);
    s["curve_params"] = json!({"p": g.p, "a": g.a, "b": g.b, "order": inst.group_order, "cofactor": inst.cofactor, "generator": inst.generator});
    s["runs"] = json!(runs_out);
    curves.push(s);

    // (b) K_0 / GF(2^41), A = 82.
    let kin = koblitz_instance(0, 41).expect("K_0/GF(2^41)");
    let bg = BinaryGroup(&kin.fast);
    let r = kin.r;
    let mut by_k = Vec::new();
    let mut runs_out = Vec::new();
    let mut all_ok = true;
    for &k in &ks {
        let mut tot = Vec::new();
        for b in 0..batches {
            let mut rng = StdRng::seed_from_u64(derive(seed, &[3, k as u64, b as u64]));
            let planted: Vec<u64> = (0..k).map(|_| rng.gen_range(1..r)).collect();
            let mut ops = GroupOps::default();
            let targets: Vec<_> = planted.iter().map(|&d| bg.mul(&mut ops, kin.generator, d)).collect();
            let res = signed_frobenius_rho_batch(&kin, &targets, derive(seed, &[4, k as u64, b as u64])).unwrap();
            let ok = res.per_target.iter().zip(&planted).all(|(t, &d)| t.verified && t.recovered == Some(d));
            all_ok &= ok;
            tot.push(res.gae);
            runs_out.push(json!({"k": k, "batch": b, "gae": res.gae, "all_verified": ok,
                "solved_on_an_earlier_trail": res.counters.get("solved_on_an_earlier_trail")}));
        }
        by_k.push((k, tot));
    }
    let mut s = summarize(&kin.name, r as f64, 82.0, by_k, 0x5735, "fresh: signed_frobenius_rho_batch");
    s["all_verified"] = json!(all_ok);
    s["runs"] = json!(runs_out);
    curves.push(s);

    json!({"schema": "lit2026-w3-multitarget-measure/v1", "batches_per_k": batches, "ks": ks, "seed": seed,
        "unit": "group-addition equivalents (gae) per batch: shared jump table, each target's setup, starts, walks, verification",
        "fits": curves})
}

fn sig(x: f64) -> String {
    if !x.is_finite() {
        return "—".into();
    }
    let a = x.abs();
    if a == 0.0 {
        "0".into()
    } else if !(1e-3..1e6).contains(&a) {
        format!("{x:.3e}")
    } else if a >= 100.0 {
        format!("{x:.0}")
    } else {
        format!("{x:.3}")
    }
}

fn f(v: &Value) -> f64 {
    num(v).unwrap_or(f64::NAN)
}

/// Markdown tables from the three JSON outputs (no new numbers).
fn report(t: &Value, fit: &Value, meas: &Value) -> String {
    let mut s = String::new();
    let rows = t["rows"].as_array().unwrap();
    s.push_str("### Per family (SGGM bound, group-operation units)\n\n");
    s.push_str("| family | rows | rows with δ_upper < 1 | min measured / min(√(r/A), 1/δ) | min measured / min(√(r/A), 1/δ), δ_upper < 1 only | min measured / √(r/A) | rows below any bound |\n|:--|--:|--:|--:|--:|--:|--:|\n");
    let mut fams: Vec<&str> = rows.iter().map(|r| r["family"].as_str().unwrap()).collect();
    fams.sort();
    fams.dedup();
    for fam in fams {
        let rs: Vec<&Value> = rows.iter().filter(|r| r["family"] == fam).collect();
        let nv: Vec<&&Value> = rs.iter().filter(|r| f(&r["delta_counting_upper"]) < 1.0).collect();
        let minr = rs.iter().map(|r| f(&r["ratio_measured_over_theta"])).fold(f64::INFINITY, f64::min);
        let minnv = nv.iter().map(|r| f(&r["ratio_measured_over_theta"])).fold(f64::INFINITY, f64::min);
        let minsq = rs.iter().map(|r| f(&r["measured"]) / f(&r["sqrt_r_over_aut"])).fold(f64::INFINITY, f64::min);
        let below = rs
            .iter()
            .filter(|r| r["flags"].as_array().unwrap().iter().any(|x| x.as_str().unwrap().starts_with("below")))
            .count();
        let unit = if rs[0]["measured_is_group_ops"] == true { "" } else { " (probes, not group ops)" };
        s.push_str(&format!("| {fam}{unit} | {} | {} | {} | {} | {} | {below} |\n", rs.len(), nv.len(), sig(minr), sig(minnv), sig(minsq)));
    }
    s.push_str("\n### Rows with δ_upper < 1 (the bound is not vacuous)\n\n| configuration | log2 r | A | m | abs F | δ measured | δ upper (count) | 1/δ | √(r/A) | bound | Thm 3.2, ε=½ | measured (GAE) | measured/bound | flags |\n|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|\n");
    for r in rows.iter().filter(|r| f(&r["delta_counting_upper"]) < 1.0) {
        s.push_str(&format!(
            "| {} | {:.2} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |\n",
            r["id"].as_str().unwrap(),
            f(&r["log2_r"]),
            f(&r["aut"]),
            r["m"],
            f(&r["factor_base_points"]),
            sig(f(&r["delta_measured"])),
            sig(f(&r["delta_counting_upper"])),
            sig(f(&r["inv_delta"])),
            sig(f(&r["sqrt_r_over_aut"])),
            sig(f(&r["bound_theta"])),
            sig(f(&r["bound_explicit_eps_half"])),
            sig(f(&r["measured"])),
            sig(f(&r["ratio_measured_over_theta"])),
            r["flags"].as_array().unwrap().iter().map(|x| x.as_str().unwrap()).collect::<Vec<_>>().join(", ")
        ));
    }
    s.push_str("\n### Compact-orbit ladder (4-sum relations; probes are S3 solve + root-table lookup, not group operations)\n\n| rung | log2 r | A | B = abs F | K | δ upper | measured | probes | probes / √(r/A) | S = K²n | log2(S·T²) − log2(r/A) | S·T/(r/A) | frozen / rank mean | P(T ≤ frozen) if exponential |\n|:--|--:|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|\n");
    for r in rows.iter().filter(|r| r["family"] == "koblitz_compact_orbit") {
        let id = r["id"].as_str().unwrap();
        let (rung, what) = id.split_once('/').unwrap();
        s.push_str(&format!(
            "| {rung} | {:.2} | {} | {} | {} | {} | {what} | {} | {} | {} | {:.1} | {} | {} | {} |\n",
            f(&r["log2_r"]),
            f(&r["aut"]),
            f(&r["factor_base_points"]),
            f(&r["orbit_columns"]),
            sig(f(&r["delta_counting_upper"])),
            sig(f(&r["measured"])),
            sig(f(&r["probes_over_sqrt_r_over_aut"])),
            sig(f(&r["advice_entries_k2n"])),
            f(&r["log2_s_t2_minus_log2_r_over_a"]),
            sig(f(&r["s_t_over_r_over_a"])),
            sig(f(&r["target_over_rank_mean"])),
            sig(f(&r["p_exponential_at_most_observed"])),
        ));
    }
    s.push_str("\nAt the rank-stage mean probes per relation, S·T/(r/A) is: ");
    let parts: Vec<String> = rows
        .iter()
        .filter(|r| r["id"].as_str().unwrap().ends_with("target_probes_online"))
        .map(|r| format!("{} {}", r["id"].as_str().unwrap().split('/').next().unwrap(), sig(f(&r["s_t_over_r_over_a_at_rank_mean"]))))
        .collect();
    s.push_str(&parts.join(", "));
    s.push_str(".\n\n### Multiplier closure (non-adaptive table arm)\n\n| arm | r | A | pairs | canonical entries | mean first-hit probes T | pairs·T/r | canonical·T/(r/A) | logs | target |\n|:--|--:|--:|--:|--:|--:|--:|--:|:--|:--|\n");
    for p in t["precompute_multiplier_closure"].as_array().unwrap() {
        s.push_str(&format!(
            "| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |\n",
            p["id"].as_str().unwrap(),
            f(&p["r"]),
            f(&p["aut"]),
            f(&p["pairs"]),
            f(&p["canonical_entries"]),
            sig(f(&p["mean_first_hit_probes"])),
            sig(f(&p["pairs_times_t_over_r"])),
            sig(f(&p["canonical_times_t_over_r_over_a"])),
            p["logs_ok"],
            p["target_ok"]
        ));
    }
    s.push_str("\n### Record checks\n\n");
    for c in t["record_checks"].as_array().unwrap() {
        if c["check"] == "log2_r" {
            s.push_str(&format!("- `{}`: recorded log2 r = {}, log2 of the exact r = {:.4}\n", c["where"].as_str().unwrap(), c["recorded"], f(&c["from_registry_r"])));
        } else {
            s.push_str(&format!(
                "- `{}` (n = {}): recorded excluded precompute {} ms; the rung's own components are {} (base load / index build / rank stage / linear algebra, ms); includes rank stage: {}; includes index build: {}\n",
                c["where"].as_str().unwrap(), c["n"], sig(f(&c["recorded_ms"])),
                ["base_load", "index_build", "rank_stage", "linear_algebra"].iter().map(|k| sig(f(&c["components_ms"][*k]))).collect::<Vec<_>>().join(" / "),
                c["includes_rank_stage"], c["includes_index_build"]));
        }
    }
    s.push_str("\n### Multi-target slope d ln(total gae) / d ln k\n\n| source | curve | log2 r | A | k | slope | 95% bootstrap | Kuhn–Struik law slope (same k) | all verified |\n|:--|:--|--:|--:|:--|--:|:--|--:|:--|\n");
    for (label, src) in [("committed records", fit), ("fresh measurement", meas)] {
        for c in src["fits"].as_array().unwrap() {
            let ci = c["slope_ci95_bootstrap"].as_array().unwrap();
            s.push_str(&format!(
                "| {label} | {} | {:.2} | {} | {:?} | {:.3} | [{:.3}, {:.3}] | {:.3} | {} |\n",
                c["curve"].as_str().unwrap(),
                f(&c["log2_r"]),
                f(&c["aut"]),
                c["ks"].as_array().unwrap().iter().map(|k| k.as_u64().unwrap()).collect::<Vec<_>>(),
                f(&c["slope"]),
                f(&ci[0]),
                f(&ci[1]),
                f(&c["kuhn_struik_law_slope_same_ks"]),
                c.get("all_verified").map(|v| v.to_string()).unwrap_or_else(|| "true (ledger §19.3)".into())
            ));
        }
    }
    s.push_str("\n| source | curve | k | batches | mean total gae | total / √(k r/A) | per target over k=1 | batch law |\n|:--|:--|--:|--:|--:|--:|--:|--:|\n");
    for (label, src) in [("records", fit), ("fresh", meas)] {
        for c in src["fits"].as_array().unwrap() {
            for p in c["per_k"].as_array().unwrap() {
                s.push_str(&format!(
                    "| {label} | {} | {} | {} | {} | {} | {} | {:.3} |\n",
                    c["curve"].as_str().unwrap(),
                    p["k"],
                    p["batches"],
                    sig(f(&p["mean_total_gae"])),
                    sig(f(&p["mean_total_over_sqrt_k_r_over_A"])),
                    sig(f(&p["per_target_over_k1"])),
                    f(&p["batch_law"])
                ));
            }
        }
    }
    s
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.get(1).map(String::as_str) == Some("report") {
        let rd = |i: usize| -> Value { serde_json::from_str(&fs::read_to_string(&args[i]).unwrap()).unwrap() };
        print!("{}", report(&rd(2), &rd(3), &rd(4)));
        return;
    }
    let out = match args.get(1).map(String::as_str) {
        Some("table") => table(Path::new(args.get(2).map(String::as_str).unwrap_or("."))),
        Some("fit") => fit(Path::new(args.get(2).map(String::as_str).unwrap_or("."))),
        Some("measure") => {
            let b = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(24);
            let bits = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(36);
            measure(b, bits)
        }
        _ => {
            eprintln!("usage: lit2026_w3_ggm_audit table|fit <repo-root> | measure [batches] [prime_bits]");
            std::process::exit(2);
        }
    };
    println!("{}", serde_json::to_string_pretty(&out).unwrap());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn delta_upper_counts_multisets() {
        // |F| = 4, m = 2, r = 1000: (4 + 10)/1000.
        let d = delta_upper(4.0, 2, 1000.0, 1.0);
        assert!((d - 0.014).abs() < 1e-12, "{d}");
        assert_eq!(delta_upper(1e6, 4, 1e9, 1.0), 1.0);
    }

    #[test]
    fn explicit_bound_inverts_theorem_3_2() {
        // At T = bound, one of the two terms equals (eps - 1/r)/2.
        let (r, d, eps) = (1e12, 1e-4, 0.5);
        let t = bound_explicit(r, d, eps);
        let e = eps - 1.0 / r;
        let lhs = (d * (3.0 * t + 2.0)).max((3.0 * t + 1.0).powi(2) / r);
        assert!((lhs - e / 2.0).abs() < 1e-9, "{lhs}");
        // sqrt branch when delta is tiny.
        let t2 = bound_explicit(1e12, 1e-12, 0.5);
        assert!((t2 - ((1e12f64 * (0.5 - 1e-12) / 2.0).sqrt() - 1.0) / 3.0).abs() < 1e-6);
    }

    #[test]
    fn dkm_root_satisfies_theorem_1_1() {
        let (n, s, eps) = (1e10, 1e3, 0.5);
        let t = dkm_nonadaptive_t(n, s, eps);
        let v = 3.0 * t * t / n + 4.0 * std::f64::consts::LN_2 * s * t / n;
        assert!((v - eps).abs() < 1e-9);
    }

    #[test]
    fn committed_records_audit() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"));
        let t = table(root);
        let rows = t["rows"].as_array().unwrap();
        assert_eq!(rows.len(), 208);
        // No configuration's recorded cost sits below its bound.
        for r in rows {
            let flags: Vec<&str> = r["flags"].as_array().unwrap().iter().map(|x| x.as_str().unwrap()).collect();
            assert!(!flags.iter().any(|f| f.starts_with("below") || f.contains("SIGNIFICANT")), "{}", r["id"]);
        }
        // The one recorded log2 r that disagrees with the exact order (n61).
        let checks = t["record_checks"].as_array().unwrap();
        let log2: Vec<&Value> = checks.iter().filter(|c| c["check"] == "log2_r").collect();
        assert_eq!(log2.len(), 1);
        assert_eq!(log2[0]["slug"], "icv1-f2m61-t158598901-ab42b6c5");
        // n61 and n71 record an excluded precompute without the rank stage.
        let ex: Vec<u64> = checks
            .iter()
            .filter(|c| c["check"] == "excluded_precompute")
            .map(|c| c["n"].as_u64().unwrap())
            .collect();
        assert_eq!(ex, vec![61, 71]);
        let fit = fit(root);
        assert_eq!(fit["fits"].as_array().unwrap().len(), 3);
    }

    #[test]
    fn slope_of_exact_square_root_is_half() {
        let by_k: Vec<(usize, Vec<f64>)> =
            [1usize, 2, 4, 8, 16].iter().map(|&k| (k, vec![1000.0 * (k as f64).sqrt(); 3])).collect();
        let (s, lo, hi) = slope_with_ci(&by_k, 1);
        assert!((s - 0.5).abs() < 1e-12 && (lo - 0.5).abs() < 1e-12 && (hi - 0.5).abs() < 1e-12);
        assert!((batch_law(1) - 1.0).abs() < 1e-15);
        assert!((batch_law(2) - 0.75).abs() < 1e-15);
    }
}

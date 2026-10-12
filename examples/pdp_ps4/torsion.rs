//! The 2-torsion half of PS-4.
//!
//! Translation by `T = (0, 1)` is `x ↦ 1/x`; in the repository's frame
//! `u = 1/(x + 1)` it is `u ↦ u + 1`, and in the ONB `1` is the all-ones word,
//! so `u ↦ u + 1` is the complement `wt(u) ↦ n − wt(u)`.  A decomposition
//! `R = P₁ + P₂ + P₃` stays a decomposition of the same `R` when two summands
//! are translated: `(P₁ + T) + (P₂ + T) + P₃ = R`, and the four thirds of the
//! pair `(P₁ + T, P₂ + T)` are those of `(P₁, P₂)`.  That is a symmetry of a
//! *fixed* target (unlike Frobenius), but only on a base that `T` preserves.
//!
//! This module (a) measures how much of each weight base `T` preserves, in
//! the `x`-frame and in the `u`-frame, and (b) tests the symmetry on the
//! smallest `T`-closed base that contains the weight base,
//! `B^± = {P : wt(u) ≤ w or wt(u) ≥ n − w}`: the same pair loop with and
//! without the `T` reduction against the pair-table oracle on that base.

use crate::base::Base;
use crate::curve::{lift, sorted3};
use crate::field::Onb;
use crate::kernel::{dedup, Second};
use crate::oracle::{o3, FxSet};
use crate::run::{both, both_json};
use crate::targets::{enumerate_weight, panel_targets};
use crate::{arg, arg_s, hex};
use std::io::Write;

/// Per-frame statistics of the weight base under `T`.
fn frame_stats(f: &Onb, w: usize) -> String {
    let n = f.n;
    let words: Vec<u64> = enumerate_weight(n, w)
        .into_iter()
        .map(|v| v as u64)
        .collect();
    // x-frame: base {x : wt(x) <= w}; T sends x != 0 to 1/x.
    let (mut xb, mut xb_pres, mut xb_all, mut xb_all_pres) = (0usize, 0usize, 0usize, 0usize);
    for &x in &words {
        if x == 0 {
            continue; // x = 0 is T itself; T + T = O, outside the finite base
        }
        xb_all += 1;
        let img = f.inv_c(x);
        let pres = img.count_ones() as usize <= w;
        xb_all_pres += usize::from(pres);
        if lift(f, x).is_some() {
            xb += 1;
            xb_pres += usize::from(pres);
        }
    }
    // u-frame: base {P : wt(u(P)) <= w}, u = 1/(x+1); x = 1/u + 1.
    let (mut ub, mut ub_pres, mut ub_identity_ok) = (0usize, 0usize, true);
    for &u in &words {
        if u == 0 {
            continue; // u = 0 is the point at infinity
        }
        let x = f.inv_c(u) ^ f.mask;
        if lift(f, x).is_none() {
            continue;
        }
        ub += 1;
        if x == 0 {
            continue; // T: its translate is O
        }
        let xi = f.inv_c(x); // x(P + T)
        let denom = xi ^ f.mask;
        if denom == 0 {
            continue; // x(P+T) = 1: u = infinity, not a finite frame value
        }
        let u2 = f.inv_c(denom);
        ub_identity_ok &= u2 == (u ^ f.mask);
        ub_pres += usize::from(u2.count_ones() as usize <= w);
    }
    format!(
        "{{\"x_frame\":{{\"base_abscissae_nonzero_liftable\":{xb},\"T_image_in_base\":{xb_pres},\"all_nonzero_words\":{xb_all},\"all_nonzero_words_T_image_light\":{xb_all_pres}}},\"u_frame\":{{\"base_abscissae_liftable\":{ub},\"T_image_in_base\":{ub_pres},\"translation_is_u_plus_1_on_all_tested\":{ub_identity_ok},\"complement_weight_exceeds_w\":{}}}}}",
        n - w > w
    )
}

pub fn cmd_torsion(args: &[String]) {
    let n: usize = arg(args, "--n").expect("--n");
    let w: usize = arg(args, "--w").expect("--w");
    let ntargets: usize = arg(args, "--targets").unwrap_or(6);
    let nplanted: usize = arg(args, "--planted").unwrap_or(2);
    let seed: u64 = arg(args, "--seed").unwrap_or(20261011);
    let only: Option<Vec<usize>> =
        arg_s(args, "--only").map(|s| s.split(',').filter_map(|v| v.parse().ok()).collect());
    let mut out: Box<dyn Write> = match arg_s(args, "--out") {
        Some(p) => Box::new(
            std::fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(p)
                .unwrap(),
        ),
        None => Box::new(std::io::stdout()),
    };
    let f = Onb::new(n).expect("type-II ONB");
    writeln!(
        out,
        "{{\"kind\":\"torsion_frames\",\"n\":{n},\"w\":{w},\"frames\":{}}}",
        frame_stats(&f, w)
    )
    .unwrap();

    // The T-closed base B^± in x-terms.
    let light: Vec<u64> = enumerate_weight(n, w)
        .into_iter()
        .map(|v| v as u64)
        .collect();
    let mut words: Vec<u64> = light.iter().flat_map(|&u| [u, u ^ f.mask]).collect();
    words.sort_unstable();
    words.dedup();
    let mut pairs: Vec<(u64, u64)> = Vec::new();
    for &u in &words {
        if u == 0 {
            continue;
        }
        let x = f.inv_c(u) ^ f.mask;
        if let Some(y) = lift(&f, x) {
            pairs.push((x, y));
        }
    }
    pairs.sort_unstable();
    let bt = Base {
        xs: pairs.iter().map(|p| p.0).collect(),
        ys: pairs.iter().map(|p| p.1).collect(),
        words: words.len() - 1,
    };
    let m = bt.xs.len();
    let set: FxSet = bt.xs.iter().copied().collect();
    // Closure: every non-zero abscissa has 1/x in the base; x = 0 is the boundary class.
    let closed = bt.xs.iter().all(|&x| x == 0 || set.contains(&f.inv_c(x)));
    let has_zero = bt.xs.first() == Some(&0);
    // Representatives of the classes {x, 1/x}: x < 1/x, and x = 0 alone.
    let mut reps: Vec<usize> = Vec::new();
    for (i, &x) in bt.xs.iter().enumerate() {
        if x == 0 || x < f.inv_c(x) {
            reps.push(i);
        }
    }
    let (m64, r64) = (m as u64, reps.len() as u64);
    writeln!(
        out,
        "{{\"kind\":\"torsion_base\",\"n\":{n},\"w\":{w},\"liftable_x\":{m},\"closed_under_x_to_1_over_x\":{closed},\"contains_T\":{has_zero},\"classes\":{},\"pairs_ordered\":{},\"pairs_reduced\":{},\"factor\":{:.4}}}",
        reps.len(),
        m64 * m64,
        r64 * m64,
        (m64 * m64) as f64 / (r64 * m64) as f64
    )
    .unwrap();

    let (targets, _) = panel_targets(n, w, 3, ntargets, nplanted, seed);
    let rs: Vec<(u64, u64)> = targets.iter().map(|t| (t.x, t.y)).collect();
    let oracle = o3(&f, &bt, &rs);
    let all_idx: Vec<usize> = (0..m).collect();
    let leaf = |t: u64| {
        f.note_lookup();
        set.contains(&t)
    };
    for (ti, tg) in targets.iter().enumerate() {
        if only.as_ref().is_some_and(|o| !o.contains(&ti)) {
            continue;
        }
        let r = (tg.x, tg.y);
        let sols = &oracle.sols[ti];
        let u = both(&f, &bt, r, &all_idx, Second::All, &leaf);
        let t = both(&f, &bt, r, &reps, Second::All, &leaf);
        // Closure of the reduced set: translate two non-T members by T.
        let mut clo: Vec<[u64; 3]> = Vec::new();
        for s in &t.full.distinct {
            clo.push(*s);
            for (i, j) in [(0usize, 1usize), (0, 2), (1, 2)] {
                if s[i] != 0 && s[j] != 0 {
                    let k = 3 - i - j;
                    clo.push(sorted3(f.inv_c(s[i]), f.inv_c(s[j]), s[k]));
                }
            }
        }
        let clo = dedup(clo);
        let u_eq = u.full.distinct == *sols;
        let t_subset = t
            .full
            .distinct
            .iter()
            .all(|s| sols.binary_search(s).is_ok());
        // Without the closure the reduced set is a subset; with it, it is the whole set.
        writeln!(
            out,
            "{{\"kind\":\"torsion_target\",\"n\":{n},\"w\":{w},\"idx\":{ti},\"stratum\":\"{}\",\"x_r\":\"{}\",\"base_liftable_x\":{m},\"oracle_solutions\":{},\"verdict_oracle\":{},\"loops\":{{\"U_ord_Tclosed\":{},\"T_rep\":{}}},\"checks\":{{\"U_ord_equals_oracle\":{u_eq},\"T_rep_subset_of_oracle\":{t_subset},\"T_rep_solutions\":{},\"closure_solutions\":{},\"closure_equals_oracle\":{},\"verdict_T_rep_equals_U\":{}}}}}",
            tg.stratum,
            hex(r.0),
            sols.len(),
            !sols.is_empty(),
            both_json(&u),
            both_json(&t),
            t.full.distinct.len(),
            clo.len(),
            clo == *sols,
            t.full.found == u.full.found
        )
        .unwrap();
        out.flush().unwrap();
    }
}

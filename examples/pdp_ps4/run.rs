//! `run --n N --w W`: the unreduced and reduced pair loops on one setting's
//! panel targets, checked against the pair-table oracle (and, when given, the
//! panel's own records).

use crate::base::{build_base, orbits, stable, Base};
use crate::curve::sorted3;
use crate::field::{Counts, Onb};
use crate::kernel::{dedup, run, LoopResult, Second};
use crate::oracle::{o3, O3};
use crate::targets::panel_targets;
use crate::{arg, arg_s, hex};
use std::fs::OpenOptions;
use std::io::Write;

pub struct Sum {
    pub pairs: u64,
    pub thirds: u64,
    pub exact_pairs: u64,
    pub rows: u64,
    pub hits: usize,
    pub found: bool,
    pub counts: Counts,
    pub ms: f64,
    pub distinct: Vec<[u64; 3]>,
}

pub fn summ(r: LoopResult) -> Sum {
    Sum {
        pairs: r.pairs,
        thirds: r.thirds,
        exact_pairs: r.exact_pairs,
        rows: r.rows,
        hits: r.hits.len(),
        found: r.found,
        counts: r.counts,
        ms: r.nanos as f64 / 1e6,
        distinct: dedup(r.hits),
    }
}

pub fn sum_json(s: &Sum) -> String {
    format!(
        "{{\"pairs\":{},\"thirds\":{},\"exact_pairs\":{},\"rows\":{},\"hits\":{},\"distinct_solutions\":{},\"found\":{},\"counts\":{},\"ms_unisolated\":{:.3}}}",
        s.pairs, s.thirds, s.exact_pairs, s.rows, s.hits, s.distinct.len(), s.found, s.counts.json(), s.ms
    )
}

/// A loop in full-scan and first-witness mode.  The first-witness run is
/// skipped (and the full run reused) when the full scan found nothing, since
/// the two then coincide.
pub struct Both {
    pub full: Sum,
    pub first: Sum,
    pub first_reused: bool,
}

pub fn both<L: Fn(u64) -> bool>(
    f: &Onb,
    b: &Base,
    r: (u64, u64),
    firsts: &[usize],
    second: Second,
    leaf: &L,
) -> Both {
    let full = summ(run(f, b, r, firsts, second, leaf, true));
    if full.found {
        let first = summ(run(f, b, r, firsts, second, leaf, false));
        Both {
            full,
            first,
            first_reused: false,
        }
    } else {
        let first = Sum {
            pairs: full.pairs,
            thirds: full.thirds,
            exact_pairs: full.exact_pairs,
            rows: full.rows,
            hits: 0,
            found: false,
            counts: full.counts,
            ms: full.ms,
            distinct: Vec::new(),
        };
        Both {
            full,
            first,
            first_reused: true,
        }
    }
}

pub fn both_json(x: &Both) -> String {
    format!(
        "{{\"full\":{},\"first_witness\":{},\"first_witness_reused_full\":{}}}",
        sum_json(&x.full),
        sum_json(&x.first),
        x.first_reused
    )
}

pub fn parse_panel(path: &str) -> (Option<u64>, Vec<serde_json::Value>) {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut header = None;
    let mut cells: Vec<serde_json::Value> = Vec::new();
    for line in text.lines() {
        let Ok(v) = serde_json::from_str::<serde_json::Value>(line) else {
            continue;
        };
        match v["kind"].as_str() {
            Some("header") => header = v["base_points"].as_u64(),
            Some("cell") => {
                let t = v["target"].as_u64().unwrap() as usize;
                if cells.len() <= t {
                    cells.resize(t + 1, serde_json::Value::Null);
                }
                if cells[t].is_null() {
                    cells[t] = v;
                }
            }
            _ => {}
        }
    }
    (header, cells)
}

pub fn cmd_run(args: &[String]) {
    let n: usize = arg(args, "--n").expect("--n");
    let w: usize = arg(args, "--w").expect("--w");
    let ntargets: usize = arg(args, "--targets").unwrap_or(6);
    let nplanted: usize = arg(args, "--planted").unwrap_or(2);
    let seed: u64 = arg(args, "--seed").unwrap_or(20261011);
    let skip_ord = args.iter().any(|a| a == "--no-ord");
    let skip_conj = args.iter().any(|a| a == "--no-conj");
    let out_path = arg_s(args, "--out");
    let mut out: Box<dyn Write> = match &out_path {
        Some(p) => Box::new(
            OpenOptions::new()
                .create(true)
                .append(true)
                .open(p)
                .unwrap(),
        ),
        None => Box::new(std::io::stdout()),
    };
    let panel = arg_s(args, "--panel").map(|p| parse_panel(&p));
    let panel_header: Option<u64> = panel.as_ref().and_then(|(h, _)| *h);

    let f = Onb::new(n).expect("type-II ONB");
    let (targets, panel_base_points) = panel_targets(n, w, 3, ntargets, nplanted, seed);
    let b = build_base(&f, w);
    let m = b.xs.len();
    assert!(stable(&f, &b), "the base must be Frobenius-stable");
    let base_points = 2 * m - usize::from(b.xs.first() == Some(&0));
    let orb = orbits(&f, &b);
    let reps_words: std::collections::HashSet<u64> = orb.reps.iter().map(|&i| b.xs[i]).collect();
    let short = orb.short_orbits(&b, n);
    let short_json: Vec<String> = short
        .iter()
        .map(|(x, s)| format!("{{\"x\":\"{}\",\"size\":{}}}", hex(*x), s))
        .collect();
    let (m64, r64) = (m as u64, orb.reps.len() as u64);
    writeln!(
        out,
        "{{\"kind\":\"setting\",\"n\":{n},\"w\":{w},\"words\":{},\"liftable_x\":{m},\"base_points\":{base_points},\"panel_base_points\":{panel_base_points},\"base_points_match_lib\":{},\"panel_header_base_points\":{},\"base_points_match_panel_header\":{},\"reps\":{},\"orbit_cover\":{},\"orbit_cover_equals_base\":{},\"short_orbits\":[{}],\"pairs_ordered\":{},\"pairs_unordered\":{},\"pairs_reduced\":{},\"factor_vs_ordered\":{:.4},\"factor_vs_unordered\":{:.4},\"n_over_1\":{n},\"n_over_2\":{:.1},\"targets\":{}}}",
        b.words,
        base_points == panel_base_points,
        panel_header.map_or("null".to_string(), |h| h.to_string()),
        opt_b(panel_header.map(|h| h as usize == base_points)),
        orb.reps.len(),
        orb.covered(),
        orb.covered() == m,
        short_json.join(","),
        m64 * m64,
        m64 * (m64 + 1) / 2,
        r64 * m64,
        (m64 * m64) as f64 / (r64 * m64) as f64,
        (m64 * (m64 + 1) / 2) as f64 / (r64 * m64) as f64,
        n as f64 / 2.0,
        targets.len()
    )
    .unwrap();

    // The oracle, one pass for every target.
    let rs: Vec<(u64, u64)> = targets.iter().map(|t| (t.x, t.y)).collect();
    let t0 = std::time::Instant::now();
    let oracle: O3 = o3(&f, &b, &rs);
    let o3_ms = t0.elapsed().as_secs_f64() * 1e3;
    eprintln!(
        "n={n} w={w}: base {m} abscissae, oracle {} keys, {} pair sums, {:.0} ms",
        oracle.keys, oracle.pair_sums, o3_ms
    );

    let all_idx: Vec<usize> = (0..m).collect();
    let leaf = |t: u64| f.weight_le(t, w as u32);
    for (ti, tg) in targets.iter().enumerate() {
        let r = (tg.x, tg.y);
        assert!(
            crate::curve::on_curve(&f, r.0, r.1),
            "target {ti} is not on the curve"
        );
        let x_in = r.0.count_ones() as usize <= w;
        let inv_in = r.0 != 0 && f.inv_c(r.0).count_ones() as usize <= w;
        let sols = &oracle.sols[ti];
        let verdict_o = !sols.is_empty();
        // Panel record, if given.
        let pj = panel.as_ref().and_then(|(_, cells)| cells.get(ti)).cloned();
        let (pm, pt, pc) = match &pj {
            Some(c) if !c.is_null() => (
                Some(c["x_r"].as_str() == Some(&format!("0x{:x}", r.0))),
                c["truth_decomposable"].as_bool(),
                c["ref_count"].as_u64(),
            ),
            _ => (None, None, None),
        };
        let u_unord = both(&f, &b, r, &all_idx, Second::FromFirst, &leaf);
        let u_ord = if skip_ord {
            None
        } else {
            Some(both(&f, &b, r, &all_idx, Second::All, &leaf))
        };
        let r_rep = both(&f, &b, r, &orb.reps, Second::All, &leaf);
        // Predicted reduced set: oracle triples with a representative member.
        let pred: Vec<[u64; 3]> = sols
            .iter()
            .filter(|t| t.iter().any(|x| reps_words.contains(x)))
            .copied()
            .collect();
        let eq_unord = u_unord.full.distinct == *sols;
        let eq_ord = u_ord.as_ref().map(|u| u.full.distinct == *sols);
        let rep_subset = r_rep
            .full
            .distinct
            .iter()
            .all(|t| sols.binary_search(t).is_ok());
        let rep_eq_pred = r_rep.full.distinct == pred;
        // Orbit closure: the reduced loop on every conjugate of the target.
        let conj_json = if skip_conj {
            "null".to_string()
        } else {
            let mut union: Vec<[u64; 3]> = Vec::new();
            let mut tot = Counts::default();
            let mut pairs = 0u64;
            let mut distinct_targets = 0usize;
            let mut ms = 0.0;
            let mut any = false;
            let mut seen: Vec<(u64, u64)> = Vec::new();
            for k in 0..n {
                let rk = (f.frob_raw(r.0, k), f.frob_raw(r.1, k));
                if seen.contains(&rk) {
                    break;
                }
                seen.push(rk);
                distinct_targets += 1;
                let s = summ(run(&f, &b, rk, &orb.reps, Second::All, &leaf, true));
                pairs += s.pairs;
                tot.plus(&s.counts);
                ms += s.ms;
                any |= s.found;
                for t in &s.distinct {
                    let back = n - k;
                    union.push(sorted3(
                        f.frob_raw(t[0], back),
                        f.frob_raw(t[1], back),
                        f.frob_raw(t[2], back),
                    ));
                }
            }
            let union = dedup(union);
            format!(
                "{{\"conjugates\":{distinct_targets},\"pairs\":{pairs},\"counts\":{},\"ms_unisolated\":{ms:.3},\"closure_solutions\":{},\"closure_equals_unreduced\":{},\"verdict_any_conjugate\":{}}}",
                tot.json(),
                union.len(),
                union == *sols,
                any
            )
        };
        writeln!(
            out,
            "{{\"kind\":\"target\",\"n\":{n},\"w\":{w},\"idx\":{ti},\"stratum\":\"{}\",\"x_r\":\"{}\",\"y_r\":\"{}\",\"x_r_weight_le_w\":{},\"inv_x_r_weight_le_w\":{},\"oracle\":{{\"decomposable\":{},\"solutions\":{},\"ref_count\":{},\"pair_sums\":{},\"keys\":{},\"ms_unisolated\":{:.1}}},\"panel\":{{\"x_r_match\":{},\"truth_decomposable\":{},\"ref_count\":{},\"verdict_match\":{},\"ref_count_match\":{}}},\"loops\":{{\"U_unord\":{},\"U_ord\":{},\"R_rep\":{}}},\"checks\":{{\"U_unord_equals_oracle\":{},\"U_ord_equals_oracle\":{},\"R_rep_subset_of_oracle\":{},\"R_rep_equals_predicted\":{},\"predicted_solutions\":{},\"verdict_oracle\":{},\"verdict_U_unord\":{},\"verdict_R_rep\":{},\"verdict_R_rep_equals_U\":{}}},\"conjugates\":{}}}",
            tg.stratum,
            hex(r.0),
            hex(r.1),
            x_in,
            inv_in,
            verdict_o,
            sols.len(),
            oracle.ref_count[ti],
            oracle.pair_sums,
            oracle.keys,
            o3_ms,
            opt_b(pm),
            opt_b(pt),
            pc.map_or("null".to_string(), |c| c.to_string()),
            opt_b(pt.map(|t| t == verdict_o)),
            opt_b(pc.map(|c| c == oracle.ref_count[ti])),
            both_json(&u_unord),
            u_ord.as_ref().map_or("null".to_string(), both_json),
            both_json(&r_rep),
            eq_unord,
            opt_b(eq_ord),
            rep_subset,
            rep_eq_pred,
            pred.len(),
            verdict_o,
            u_unord.full.found,
            r_rep.full.found,
            r_rep.full.found == u_unord.full.found,
            conj_json
        )
        .unwrap();
        out.flush().unwrap();
    }
}

fn opt_b(v: Option<bool>) -> String {
    v.map_or("null".to_string(), |b| b.to_string())
}

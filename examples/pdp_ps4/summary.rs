//! Markdown tables from the raw records.  Every number printed is read from a
//! record or is a ratio of two recorded numbers (the formula is in the table
//! header); nothing is measured here.

use crate::arg_s;
use crate::field::Counts;
use crate::model::{bitops, onb_costs};
use serde_json::Value;
use std::collections::BTreeMap;

fn read_lines(path: &str) -> Vec<Value> {
    match std::fs::read_to_string(path) {
        Ok(t) => t
            .lines()
            .filter_map(|l| serde_json::from_str::<Value>(l).ok())
            .collect(),
        Err(_) => Vec::new(),
    }
}

fn u(v: &Value) -> u64 {
    v.as_u64().unwrap_or(0)
}

fn fl(v: &Value) -> f64 {
    v.as_f64().unwrap_or(0.0)
}

fn counts_of(v: &Value) -> Counts {
    Counts {
        mul: u(&v["mul"]),
        sqr: u(&v["sqr"]),
        sqr_k: u(&v["sqr_k"]),
        inv: u(&v["inv"]),
        batch_calls: u(&v["batch_calls"]),
        batch_elems: u(&v["batch_elems"]),
        batch_mul: u(&v["batch_mul"]),
        add: u(&v["add"]),
        popcnt: u(&v["popcnt"]),
        lookup: u(&v["lookup"]),
    }
}

/// Sum of the counted operations, the pairs and thirds of one loop mode over targets.
#[derive(Default, Clone)]
struct Acc {
    n: usize,
    pairs: u64,
    thirds: u64,
    hits: u64,
    distinct: u64,
    counts: Counts,
    targets: u64,
}

impl Acc {
    fn add(&mut self, n: usize, s: &Value) {
        self.n = n;
        self.pairs += u(&s["pairs"]);
        self.thirds += u(&s["thirds"]);
        self.hits += u(&s["hits"]);
        self.distinct += u(&s["distinct_solutions"]);
        self.counts.plus(&counts_of(&s["counts"]));
        self.targets += 1;
    }
    fn mulq(&self) -> u64 {
        self.counts.mulq()
    }
    fn bits(&self) -> f64 {
        bitops(
            &self.counts,
            self.counts.popcnt + self.counts.lookup,
            &onb_costs(self.n),
        )
    }
}

fn f3(x: f64) -> String {
    if x.is_nan() {
        "n/a".to_string()
    } else {
        format!("{x:.3}")
    }
}

fn ratio(a: f64, b: f64) -> f64 {
    if b == 0.0 {
        f64::NAN
    } else {
        a / b
    }
}

pub fn cmd_summary(args: &[String]) {
    let raw = arg_s(args, "--raw").unwrap_or_else(|| "raw".to_string());
    let runs = read_lines(&format!("{raw}/runs.jsonl"));
    // Last record wins per (kind, n, w, idx).
    let mut settings: BTreeMap<(usize, usize), Value> = BTreeMap::new();
    let mut targets: BTreeMap<(usize, usize, u64), Value> = BTreeMap::new();
    for r in &runs {
        let n = u(&r["n"]) as usize;
        let w = u(&r["w"]) as usize;
        match r["kind"].as_str() {
            Some("setting") => {
                settings.insert((n, w), r.clone());
            }
            Some("target") => {
                targets.insert((n, w, u(&r["idx"])), r.clone());
            }
            _ => {}
        }
    }

    println!("### Table 1. The bases and their Frobenius orbits\n");
    println!("`m` = abscissae of weight `<= w` that lift (one point per abscissa; the base has `2m - 1` points, `x = 0` being `T`). `reps` = Frobenius orbits of the base. Pairs: ordered `m^2`, unordered `m(m+1)/2` (the PS-0/PS-2 loop), reduced `reps * m`.\n");
    println!("| n | w | words wt<=w | m | base points (= panel's) | reps | orbit sizes cover m | short orbits (x, size) | pairs ordered | pairs unordered | pairs reduced | factor vs ordered (m/reps) | n | factor vs unordered | n/2 | measured / predicted (ordered) |");
    println!("|--:|--:|--:|--:|--:|--:|:-:|:--|--:|--:|--:|--:|--:|--:|--:|--:|");
    for ((n, w), s) in &settings {
        let short: Vec<String> = s["short_orbits"]
            .as_array()
            .map(|a| {
                a.iter()
                    .map(|o| format!("({}, {})", o["x"].as_str().unwrap_or("?"), u(&o["size"])))
                    .collect()
            })
            .unwrap_or_default();
        println!(
            "| {n} | {w} | {} | {} | {} ({}) | {} | {} | {} | {} | {} | {} | {:.3} | {n} | {:.3} | {:.1} | {:.3} |",
            u(&s["words"]),
            u(&s["liftable_x"]),
            u(&s["base_points"]),
            match (s["base_points_match_lib"].as_bool(), s["base_points_match_panel_header"].as_bool()) {
                (Some(true), Some(true)) => "lib and panel header match",
                (Some(true), None) => "lib match, no panel run",
                _ => "MISMATCH",
            },
            u(&s["reps"]),
            if s["orbit_cover_equals_base"].as_bool() == Some(true) { "yes" } else { "NO" },
            short.join(" "),
            u(&s["pairs_ordered"]),
            u(&s["pairs_unordered"]),
            u(&s["pairs_reduced"]),
            fl(&s["factor_vs_ordered"]),
            fl(&s["factor_vs_unordered"]),
            fl(&s["n_over_2"]),
            fl(&s["factor_vs_ordered"]) / *n as f64
        );
    }

    // Per-setting accumulators.
    #[derive(Default)]
    struct Row {
        u_unord: Acc,
        u_ord: Acc,
        r_rep: Acc,
        u_unord_fw: Acc,
        u_ord_fw: Acc,
        r_rep_fw: Acc,
        conj_pairs: u64,
        conj_counts: Counts,
        conj_targets: u64,
        ntargets: u64,
        decomposable: u64,
        panel_n: u64,
        panel_verdict_ok: u64,
        panel_count_ok: u64,
        panel_x_ok: u64,
        u_unord_eq: u64,
        u_ord_eq: u64,
        u_ord_n: u64,
        rep_subset: u64,
        rep_pred: u64,
        verdict_same: u64,
        sols_total: u64,
        sols_rep: u64,
        closure_eq: u64,
        closure_n: u64,
        exact_pairs: u64,
        pairs_u_unord: u64,
        planted: u64,
        planted_lost: u64,
        uniform_dec: u64,
        uniform_lost: u64,
        degenerate: u64,
        // Verdict-only cost split by what the oracle says about the target.
        dec_u_fw: Acc,
        dec_r_fw: Acc,
        dec_r_says_no: u64,
        nodec_u: Acc,
        nodec_r: Acc,
    }
    let mut rows: BTreeMap<(usize, usize), Row> = BTreeMap::new();
    for ((n, w, _), t) in &targets {
        let row = rows.entry((*n, *w)).or_default();
        row.ntargets += 1;
        let l = &t["loops"];
        row.u_unord.add(*n, &l["U_unord"]["full"]);
        row.u_unord_fw.add(*n, &l["U_unord"]["first_witness"]);
        row.r_rep.add(*n, &l["R_rep"]["full"]);
        row.r_rep_fw.add(*n, &l["R_rep"]["first_witness"]);
        if !l["U_ord"].is_null() {
            row.u_ord.add(*n, &l["U_ord"]["full"]);
            row.u_ord_fw.add(*n, &l["U_ord"]["first_witness"]);
            row.u_ord_n += 1;
            row.u_ord_eq += u64::from(t["checks"]["U_ord_equals_oracle"].as_bool() == Some(true));
        }
        row.exact_pairs += u(&l["U_unord"]["full"]["exact_pairs"]);
        row.pairs_u_unord += u(&l["U_unord"]["full"]["pairs"]);
        let dec = t["oracle"]["decomposable"].as_bool() == Some(true);
        row.decomposable += u64::from(dec);
        if dec {
            row.dec_u_fw.add(*n, &l["U_unord"]["first_witness"]);
            row.dec_r_fw.add(*n, &l["R_rep"]["first_witness"]);
            row.dec_r_says_no += u64::from(l["R_rep"]["full"]["found"].as_bool() != Some(true));
        } else {
            row.nodec_u.add(*n, &l["U_unord"]["full"]);
            row.nodec_r.add(*n, &l["R_rep"]["full"]);
        }
        if !t["panel"]["verdict_match"].is_null() {
            row.panel_n += 1;
            row.panel_verdict_ok += u64::from(t["panel"]["verdict_match"].as_bool() == Some(true));
            row.panel_count_ok += u64::from(t["panel"]["ref_count_match"].as_bool() == Some(true));
            row.panel_x_ok += u64::from(t["panel"]["x_r_match"].as_bool() == Some(true));
        }
        let c = &t["checks"];
        row.u_unord_eq += u64::from(c["U_unord_equals_oracle"].as_bool() == Some(true));
        row.rep_subset += u64::from(c["R_rep_subset_of_oracle"].as_bool() == Some(true));
        row.rep_pred += u64::from(c["R_rep_equals_predicted"].as_bool() == Some(true));
        let same = c["verdict_R_rep_equals_U"].as_bool() == Some(true);
        row.verdict_same += u64::from(same);
        row.sols_total += u(&t["oracle"]["solutions"]);
        row.sols_rep += u(&l["R_rep"]["full"]["distinct_solutions"]);
        if !t["conjugates"].is_null() {
            row.closure_n += 1;
            row.closure_eq +=
                u64::from(t["conjugates"]["closure_equals_unreduced"].as_bool() == Some(true));
            row.conj_pairs += u(&t["conjugates"]["pairs"]);
            row.conj_counts.plus(&counts_of(&t["conjugates"]["counts"]));
            row.conj_targets += 1;
        }
        let planted = t["stratum"].as_str() == Some("planted");
        row.planted += u64::from(planted);
        row.planted_lost += u64::from(planted && !same);
        row.uniform_dec += u64::from(!planted && dec);
        row.uniform_lost += u64::from(!planted && !same);
        row.degenerate += u64::from(
            t["x_r_weight_le_w"].as_bool() == Some(true)
                || t["inv_x_r_weight_le_w"].as_bool() == Some(true),
        );
    }

    println!("\n### Table 2. Candidates and counted operations per target, reduced against unreduced (sums over the targets of the setting)\n");
    println!("Full-scan mode. `pairs` = `(P1, P2)` pairs visited; `thirds` = candidate thirds tested (4 per pair); `M` = counted field multiplications (`mul` + `batch_mul`, the library's accounting, inversion products included); `bits` = the bit-operation model of `model.rs` applied to the counted operations (a convention, not a measurement). `U` = unreduced, `R` = reduced (one representative per Frobenius orbit as first summand, every base point as second). Op factor / pair factor is the fraction of the candidate factor that the counted operations keep (the PS-4 success criterion asks for at least 0.8).\n");
    println!("| n | w | targets | pairs U-unord | pairs U-ord | pairs R | pair factor vs unord | pair factor vs ord | M U-unord | M U-ord | M R | M factor vs unord | M factor vs ord | M/pair U-unord | M/pair R | bits factor vs unord | bits factor vs ord | op/pair factor (unord) | op/pair factor (ord) |");
    println!("|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|");
    for ((n, w), r) in &rows {
        let pf_un = ratio(r.u_unord.pairs as f64, r.r_rep.pairs as f64);
        let pf_or = ratio(r.u_ord.pairs as f64, r.r_rep.pairs as f64);
        let mf_un = ratio(r.u_unord.mulq() as f64, r.r_rep.mulq() as f64);
        let mf_or = ratio(r.u_ord.mulq() as f64, r.r_rep.mulq() as f64);
        let bf_un = ratio(r.u_unord.bits(), r.r_rep.bits());
        let bf_or = ratio(r.u_ord.bits(), r.r_rep.bits());
        println!(
            "| {n} | {w} | {} | {} | {} | {} | {:.3} | {:.3} | {} | {} | {} | {:.3} | {:.3} | {:.3} | {:.3} | {:.3} | {:.3} | {:.3} | {:.3} |",
            r.ntargets,
            r.u_unord.pairs,
            r.u_ord.pairs,
            r.r_rep.pairs,
            pf_un,
            pf_or,
            r.u_unord.mulq(),
            r.u_ord.mulq(),
            r.r_rep.mulq(),
            mf_un,
            mf_or,
            ratio(r.u_unord.mulq() as f64, r.u_unord.pairs as f64),
            ratio(r.r_rep.mulq() as f64, r.r_rep.pairs as f64),
            bf_un,
            bf_or,
            mf_un / pf_un,
            mf_or / pf_or
        );
    }

    println!("\n### Table 2b. Counted field operations per pair, by component\n");
    println!("Full scans pooled over the targets of the setting. `mul` includes the products inside inversions (the library's accounting); `sqr`, `sqr_k` are squarings (free coordinate permutations here, counted so the operation vocabulary matches the control); `inv` is whole inversions (two per row for `U+-`, one for the row's Montgomery batch); `batch_mul` is `3` per inverted denominator; `add` is field additions; `popcnt` is the weight tests (4 per pair). Per candidate third, divide by 4.\n");
    println!("| n | w | loop | pairs | mul/pair | sqr/pair | sqr_k/pair | inv/pair | batch_mul/pair | M/pair | add/pair | popcnt/pair |");
    println!("|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|");
    for ((n, w), r) in &rows {
        for (name, a) in [
            ("U-unord", &r.u_unord),
            ("U-ord", &r.u_ord),
            ("R", &r.r_rep),
        ] {
            if a.pairs == 0 {
                continue;
            }
            let p = a.pairs as f64;
            println!(
                "| {n} | {w} | {name} | {} | {:.4} | {:.4} | {:.4} | {:.5} | {:.4} | {:.4} | {:.3} | {:.3} |",
                a.pairs,
                a.counts.mul as f64 / p,
                a.counts.sqr as f64 / p,
                a.counts.sqr_k as f64 / p,
                a.counts.inv as f64 / p,
                a.counts.batch_mul as f64 / p,
                a.mulq() as f64 / p,
                a.counts.add as f64 / p,
                a.counts.popcnt as f64 / p
            );
        }
    }

    println!("\n### Table 3. Verdict-only mode, split by what the oracle says\n");
    println!("*Decomposable targets*: pairs and `M` until the first witness, mean per target, for the unreduced unordered loop and the reduced loop; `R says no` counts decomposable targets on which the reduced loop found nothing in its full scan. *Refuting targets*: the full scan of both loops (the reduced loop's `no` refutes only the decompositions that contain a representative; it is not a refutation of `R`).\n");
    println!("| n | w | decomposable | mean pairs to witness U-unord | mean pairs to witness R | mean M U-unord | mean M R | M R / M U-unord | R says no | refuting | mean M U-unord (full) | mean M R (full) | M R / M U-unord |");
    println!("|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|");
    for ((n, w), r) in &rows {
        let d = r.dec_u_fw.targets.max(1) as f64;
        let nd = r.nodec_u.targets.max(1) as f64;
        println!(
            "| {n} | {w} | {} | {:.0} | {:.0} | {:.0} | {:.0} | {} | {} | {} | {:.0} | {:.0} | {} |",
            r.dec_u_fw.targets,
            r.dec_u_fw.pairs as f64 / d,
            r.dec_r_fw.pairs as f64 / d,
            r.dec_u_fw.mulq() as f64 / d,
            r.dec_r_fw.mulq() as f64 / d,
            f3(ratio(r.dec_r_fw.mulq() as f64, r.dec_u_fw.mulq() as f64)),
            r.dec_r_says_no,
            r.nodec_u.targets,
            r.nodec_u.mulq() as f64 / nd,
            r.nodec_r.mulq() as f64 / nd,
            f3(ratio(r.nodec_r.mulq() as f64, r.nodec_u.mulq() as f64))
        );
    }

    println!("\n### Table 4. Verdicts and solution sets\n");
    println!("`dec` = targets the pair-table oracle calls decomposable. Panel columns: the panel's own records (`examples/hamming_nb_sat_panel`) where they exist (x_R equal / verdict equal / `ref_count` equal, out of the targets compared). `U = oracle`: kernel solution set equals the oracle's, as a set of sorted abscissa triples. `R subset`, `R = predicted`: the reduced loop's set is contained in the oracle's, and equals the oracle's triples that have a representative member. `verdict kept` = the reduced loop's verdict equals the unreduced one (fixed target). `kept sols` = distinct solutions of the reduced loop over those of the oracle, pooled over the targets (in brackets, the fraction expected if the three members of a solution were independent draws from the base, `1 - (1 - reps/m)^3`; a derived expectation, not a measurement).`closure = U`: the union over the target's Frobenius conjugates of the reduced loop's solutions, mapped back, equals the unreduced set.\n");
    println!("| n | w | targets | dec | panel x_R | panel verdict | panel ref_count | U-unord = oracle | U-ord = oracle | R subset | R = predicted | verdict kept (fixed target) | uniform lost | planted lost | kept sols | closure = U |");
    println!("|--:|--:|--:|--:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|--:|:-:|");
    for ((n, w), r) in &rows {
        let pan = |k: u64| {
            if r.panel_n == 0 {
                "n/a".to_string()
            } else {
                format!("{k}/{}", r.panel_n)
            }
        };
        // Expectation if the three members of a solution were independent draws from the
        // base: a solution is kept when at least one member is a representative.
        let frac_rep = u(&settings[&(*n, *w)]["reps"]) as f64
            / u(&settings[&(*n, *w)]["liftable_x"]) as f64;
        let expect = 1.0 - (1.0 - frac_rep).powi(3);
        let kept = if r.sols_total == 0 {
            "n/a".to_string()
        } else {
            format!(
                "{}/{} = {:.3} (independence: {:.3})",
                r.sols_rep,
                r.sols_total,
                r.sols_rep as f64 / r.sols_total as f64,
                expect
            )
        };
        println!(
            "| {n} | {w} | {} | {} | {} | {} | {} | {}/{} | {}/{} | {}/{} | {}/{} | {}/{} | {}/{} | {}/{} | {} | {}/{} |",
            r.ntargets,
            r.decomposable,
            pan(r.panel_x_ok),
            pan(r.panel_verdict_ok),
            pan(r.panel_count_ok),
            r.u_unord_eq,
            r.ntargets,
            r.u_ord_eq,
            r.u_ord_n,
            r.rep_subset,
            r.ntargets,
            r.rep_pred,
            r.ntargets,
            r.verdict_same,
            r.ntargets,
            r.uniform_lost,
            r.ntargets - r.planted,
            r.planted_lost,
            r.planted,
            kept,
            r.closure_eq,
            r.closure_n
        );
    }

    println!("\n### Table 5. What the orbit-complete reading costs, and the yield per candidate\n");
    println!("`closure pairs` = pairs the reduced loop visits over all conjugates of the target (what the closure check ran); compared with the ordered unreduced loop on the one target. `hits/third` = raw hits per candidate third, pooled over the setting's targets (the stage yield per candidate).\n");
    println!("| n | w | closure pairs / U-ord pairs | closure M / U-ord M | closure M / U-unord M | hits U-unord | thirds U-unord | hits R | thirds R | hits per third U-unord | hits per third R | R / U-unord |");
    println!("|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|");
    for ((n, w), r) in &rows {
        let yu = ratio(r.u_unord.hits as f64, r.u_unord.thirds as f64);
        let yr = ratio(r.r_rep.hits as f64, r.r_rep.thirds as f64);
        println!(
            "| {n} | {w} | {:.3} | {:.3} | {:.3} | {} | {} | {} | {} | {:.3e} | {:.3e} | {:.3} |",
            ratio(r.conj_pairs as f64, r.u_ord.pairs as f64),
            ratio(r.conj_counts.mulq() as f64, r.u_ord.mulq() as f64),
            ratio(r.conj_counts.mulq() as f64, r.u_unord.mulq() as f64),
            r.u_unord.hits,
            r.u_unord.thirds,
            r.r_rep.hits,
            r.r_rep.thirds,
            yu,
            yr,
            ratio(yr, yu)
        );
    }

    println!("\n### Table 6. Exact-path pairs\n");
    println!("Pairs that left the batched formula (a zero denominator: `U = ±P2`, or `P1 = ±R`) and were done by the exact path, in the unordered full scans.\n");
    println!("| n | w | exact pairs | pairs | targets with `wt(x_R) <= w` or `wt(1/x_R) <= w` |");
    println!("|--:|--:|--:|--:|--:|");
    for ((n, w), r) in &rows {
        println!(
            "| {n} | {w} | {} | {} | {} |",
            r.exact_pairs, r.pairs_u_unord, r.degenerate
        );
    }

    // ---- the leaf comparison -------------------------------------------------
    let control = read_lines(&format!("{raw}/control.jsonl"));
    let mut rungs: BTreeMap<u64, Value> = BTreeMap::new();
    for c in &control {
        if c["kind"].as_str() == Some("control") {
            rungs.insert(u(&c["l"]), c.clone());
        }
    }
    println!("\n### Table 7. The subspace-leaf control (PS-2 group-law kernel on the PS-0 corpus, library counters)\n");
    println!("Per pair, full scan of `i1 <= i2`, the first 6 uniform and 2 planted corpus targets of each rung. Polynomial basis, `a = b = 1`; leaf `x3 < 2^l`.\n");
    println!("| l | n | targets | base m | pairs | mul/pair | sqr/pair | sqr_k/pair | inv/pair | batch_mul/pair | M/pair | model bits/pair | verdicts agree with corpus |");
    println!("|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|");
    for (l, c) in &rungs {
        let p = &c["per_pair"];
        println!(
            "| {l} | {} | {} | {} | {} | {:.4} | {:.4} | {:.4} | {:.5} | {:.4} | {:.4} | {:.1} | {}/{} |",
            u(&c["n"]),
            u(&c["targets"]),
            u(&c["base_liftable_x"]),
            u(&c["pairs"]),
            fl(&p["mul"]),
            fl(&p["sqr"]),
            fl(&p["sqr_k"]),
            fl(&p["inv"]),
            fl(&p["batch_mul"]),
            fl(&p["mulq"]),
            fl(&c["model"]["bitops_per_pair"]),
            u(&c["verdict_agrees_with_corpus"]),
            u(&c["targets"])
        );
    }

    println!(
        "\n### Table 8. Leaf cost per pair: weight base in the ONB against the subspace control\n"
    );
    println!("Weight side: the unordered full scans of Table 2, pooled over the setting's targets, per pair. Control: Table 7. `x` columns are weight / control. `n-near` = the control rung with `3l` nearest `n` (`dn` = its `n` minus this `n`); `m-near` = the rung whose base size is nearest in ratio. `bits x` prices both sides with the model of `model.rs` at each side's own field degree; `bits x (control at n)` scales the control's modelled bits by `(n / n_l)^2`, the leading term of the schoolbook model, to put it at this setting's degree. The stop condition of PS-4 is a ratio above 2. The model is a convention, not a measurement; the wall note is Table 9.\n");
    println!("| n | w | M/pair (ONB) | popcount/pair | model bits/pair (ONB) | n-near l | dn | M/pair x | bits x | bits x (control at n) | m-near l | M/pair x | bits x | range of bits x over l=8..12 |");
    println!("|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|");
    for ((n, w), r) in &rows {
        let a = &r.u_unord;
        let m_pair = ratio(a.mulq() as f64, a.pairs as f64);
        let pc_pair = ratio(a.counts.popcnt as f64, a.pairs as f64);
        let bits_pair = ratio(a.bits(), a.pairs as f64);
        let mset = u(&settings[&(*n, *w)]["liftable_x"]) as f64;
        let mut near_n: Option<u64> = None;
        let mut near_m: Option<u64> = None;
        for (l, c) in &rungs {
            let dn = (u(&c["n"]) as i64 - *n as i64).abs();
            if near_n.is_none_or(|b| dn < (u(&rungs[&b]["n"]) as i64 - *n as i64).abs()) {
                near_n = Some(*l);
            }
            let dm = (fl(&c["base_liftable_x"]) / mset).ln().abs();
            if near_m.is_none_or(|b| dm < (fl(&rungs[&b]["base_liftable_x"]) / mset).ln().abs()) {
                near_m = Some(*l);
            }
        }
        let cell = |l: Option<u64>| -> Option<(u64, f64, f64, f64, f64)> {
            let l = l?;
            let c = rungs.get(&l)?;
            let nl = fl(&c["n"]);
            let ctl_bits = fl(&c["model"]["bitops_per_pair"]);
            Some((
                l,
                m_pair / fl(&c["per_pair"]["mulq"]),
                bits_pair / ctl_bits,
                bits_pair / (ctl_bits * (*n as f64 / nl).powi(2)),
                nl - *n as f64,
            ))
        };
        let (nl, dn, nm, nb, nbs) = match cell(near_n) {
            Some((l, m, b, bs, dn)) => (
                l.to_string(),
                format!("{dn:+.0}"),
                format!("{m:.3}"),
                format!("{b:.3}"),
                format!("{bs:.3}"),
            ),
            None => ("-".into(), "-".into(), "-".into(), "-".into(), "-".into()),
        };
        let (ml, mm, mb) = match cell(near_m) {
            Some((l, m, b, _, _)) => (l.to_string(), format!("{m:.3}"), format!("{b:.3}")),
            None => ("-".into(), "-".into(), "-".into()),
        };
        let all: Vec<f64> = rungs
            .values()
            .map(|c| bits_pair / fl(&c["model"]["bitops_per_pair"]))
            .collect();
        let mx = all.iter().cloned().fold(f64::MIN, f64::max);
        let mn = all.iter().cloned().fold(f64::MAX, f64::min);
        println!("| {n} | {w} | {m_pair:.4} | {pc_pair:.4} | {bits_pair:.1} | {nl} | {dn} | {nm} | {nb} | {nbs} | {ml} | {mm} | {mb} | {mn:.3} to {mx:.3} |");
    }

    println!("\n### Table 9. Unisolated wall time per pair (a labelled note, not a result)\n");
    println!("`ms_unisolated` fields of the unordered full scans and of the control, on a shared container with no pinning. The weight kernel is portable Rust over a software cyclic convolution; the control uses `pclmulqdq`. Different code paths, so this is not a like-for-like measure of the two leaves.\n");
    println!("| setting | ns/pair (weight, ONB) |");
    println!("|:--|--:|");
    for ((n, w), r) in &rows {
        let mut ms = 0.0;
        let mut pairs = 0u64;
        for ((tn, tw, _), t) in &targets {
            if tn == n && tw == w {
                ms += fl(&t["loops"]["U_unord"]["full"]["ms_unisolated"]);
                pairs += u(&t["loops"]["U_unord"]["full"]["pairs"]);
            }
        }
        let _ = r;
        println!("| n={n} w={w} | {:.1} |", ms * 1e6 / pairs as f64);
    }
    for (l, c) in &rungs {
        println!(
            "| control l={l} (n={}) | {:.1} |",
            u(&c["n"]),
            fl(&c["ms_unisolated"]) * 1e6 / u(&c["pairs"]) as f64
        );
    }

    // ---- torsion ------------------------------------------------------------
    let tors = read_lines(&format!("{raw}/torsion.jsonl"));
    println!("\n### Table 10. Where the 2-torsion symmetry applies: how much of each weight base `T` preserves\n");
    println!("`x` frame: `T` is `x -> 1/x`. `u` frame (`u = 1/(x+1)`): `T` is `u -> u + 1`, the complement in the ONB (`1` is the all-ones word), of weight `n - wt(u)`. Counts are over abscissae that lift.\n");
    println!("| n | w | x-frame: base | x-frame: image in base | u-frame: base | u-frame: image in base | `u -> u+1` verified | n - w > w |");
    println!("|--:|--:|--:|--:|--:|--:|:-:|:-:|");
    let mut tf: BTreeMap<(usize, usize), &Value> = BTreeMap::new();
    for t in &tors {
        if t["kind"].as_str() == Some("torsion_frames") {
            tf.insert((u(&t["n"]) as usize, u(&t["w"]) as usize), t);
        }
    }
    for ((n, w), t) in &tf {
        let f = &t["frames"];
        println!(
            "| {n} | {w} | {} | {} | {} | {} | {} | {} |",
            u(&f["x_frame"]["base_abscissae_nonzero_liftable"]),
            u(&f["x_frame"]["T_image_in_base"]),
            u(&f["u_frame"]["base_abscissae_liftable"]),
            u(&f["u_frame"]["T_image_in_base"]),
            f["u_frame"]["translation_is_u_plus_1_on_all_tested"],
            f["u_frame"]["complement_weight_exceeds_w"]
        );
    }
    println!("\n### Table 11. The T reduction on the T-closed base `B^+- = {{wt(u) <= w or wt(u) >= n - w}}`\n");
    println!("Same panel targets, ordered loops on `B^+-` with the leaf a table probe (membership in `B^+-` by `x`; the `u`-weight test would need an inversion per candidate). `T-rep` = first summand one per class `{{x, 1/x}}` (plus `T` alone), second summand every base point. `closure` = the `T-rep` set closed under translating two non-`T` members by `T`.\n");
    println!("| n | w | targets | m (B^+-) | classes | pairs U-ord | pairs T-rep | pair factor | M factor | U = oracle | T-rep subset | closure = oracle | verdict kept |");
    println!("|--:|--:|--:|--:|--:|--:|--:|--:|--:|:-:|:-:|:-:|:-:|");
    let mut tb: BTreeMap<(usize, usize), &Value> = BTreeMap::new();
    for t in &tors {
        if t["kind"].as_str() == Some("torsion_base") {
            tb.insert((u(&t["n"]) as usize, u(&t["w"]) as usize), t);
        }
    }
    for ((n, w), b) in &tb {
        let (mut cnt, mut pu, mut pt, mut mu, mut mt) = (0u64, 0u64, 0u64, 0u64, 0u64);
        let (mut ueq, mut tsub, mut ceq, mut vk) = (0u64, 0u64, 0u64, 0u64);
        for t in &tors {
            if t["kind"].as_str() == Some("torsion_target")
                && u(&t["n"]) as usize == *n
                && u(&t["w"]) as usize == *w
            {
                cnt += 1;
                pu += u(&t["loops"]["U_ord_Tclosed"]["full"]["pairs"]);
                pt += u(&t["loops"]["T_rep"]["full"]["pairs"]);
                mu += counts_of(&t["loops"]["U_ord_Tclosed"]["full"]["counts"]).mulq();
                mt += counts_of(&t["loops"]["T_rep"]["full"]["counts"]).mulq();
                let c = &t["checks"];
                ueq += u64::from(c["U_ord_equals_oracle"].as_bool() == Some(true));
                tsub += u64::from(c["T_rep_subset_of_oracle"].as_bool() == Some(true));
                ceq += u64::from(c["closure_equals_oracle"].as_bool() == Some(true));
                vk += u64::from(c["verdict_T_rep_equals_U"].as_bool() == Some(true));
            }
        }
        println!(
            "| {n} | {w} | {cnt} | {} | {} | {pu} | {pt} | {:.4} | {:.4} | {ueq}/{cnt} | {tsub}/{cnt} | {ceq}/{cnt} | {vk}/{cnt} |",
            u(&b["liftable_x"]),
            u(&b["classes"]),
            ratio(pu as f64, pt as f64),
            ratio(mu as f64, mt as f64)
        );
    }

    // ---- checks ---------------------------------------------------------------
    let fc = read_lines(&format!("{raw}/field_checks.jsonl"));
    println!("\n### Table 12. The normal-basis field against the library's `Type2Onb`\n");
    println!("| n | exhaustive | products compared | product mismatches | squarings | mismatches | inversions | mismatches | inversion products counted / formula | lifts compared | lift mismatches |");
    println!("|--:|:-:|--:|--:|--:|--:|--:|--:|:-:|--:|--:|");
    for c in fc
        .iter()
        .filter(|c| c["kind"].as_str() == Some("field_check"))
    {
        println!(
            "| {} | {} | {} | {} | {} | {} | {} | {} | {}/{} | {} | {} |",
            u(&c["n"]),
            c["exhaustive"],
            u(&c["mul_pairs"]),
            u(&c["mul_mismatch"]),
            u(&c["sqr_n"]),
            u(&c["sqr_mismatch"]),
            u(&c["inv_n"]),
            u(&c["inv_mismatch"]),
            u(&c["inv_products_counted"]),
            u(&c["inv_products_formula"]),
            u(&c["lift_n"]),
            u(&c["lift_mismatch"])
        );
    }
    let cc = read_lines(&format!("{raw}/curve_checks.jsonl"));
    println!("\n### Table 13. The curve `K_0: y^2 + xy = x^3 + 1` (a = 0, b = 1): constants, point count, translation by `T`, `S3`\n");
    println!("| n | `T` on curve | `T+T = O` | points | points expected (Koblitz recurrence) | `T`-translation points tested (exhaustive?) | `x(P+T) = 1/x` mismatches (P, -P) | group law vs library: pairs / mismatches | `S3` zero on sums: sums / nonzero | `S3` on a random third: nonzero |");
    println!("|--:|:-:|:-:|--:|--:|:--|:-:|:-:|:-:|--:|");
    for c in cc
        .iter()
        .filter(|c| c["kind"].as_str() == Some("curve_check"))
    {
        println!(
            "| {} | {} | {} | {} | {} | {} ({}) | {} / {} | {} / {} | {} / {} | {} |",
            u(&c["n"]),
            c["T_on_curve"],
            c["T_plus_T_is_O"],
            c["points"],
            c["points_expected"],
            u(&c["T_translation_points"]),
            c["T_translation_exhaustive"],
            u(&c["T_translation_mismatch"]),
            u(&c["T_translation_neg_mismatch"]),
            u(&c["group_law_pairs"]),
            u(&c["group_law_mismatch_vs_lib"]),
            u(&c["S3_sums"]),
            u(&c["S3_nonzero_on_sums"]),
            u(&c["S3_nonzero_on_random_third"])
        );
    }
    let th = read_lines(&format!("{raw}/theorem_checks.jsonl"));
    println!("\n### Table 14. The exactness theorem of `ps1_algebra/README.md` section 3, re-checked for `K_0`\n");
    println!("`V` = the set the summands are drawn from; `S4` roots in `V^3` by the resultant of `S3(x1,x2,t)` and `S3(x3,x_R,t)`; decompositions = sorted triples of the exact group law. A *degenerate* target has `x_R in V` or `1/x_R in V`. Extra roots are roots that are not decompositions.\n");
    println!("| n | V | targets | degenerate | triples tested | S4 roots | decompositions | decompositions that are not roots | extra roots F1 | extra roots F2 | unexplained extras | extras at non-degenerate targets | targets with S4 roots = decompositions | targets with extras |");
    println!("|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|");
    for c in th
        .iter()
        .filter(|c| c["kind"].as_str() == Some("theorem_check"))
    {
        println!(
            "| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |",
            u(&c["n"]),
            c["V"].as_str().unwrap_or("?"),
            u(&c["targets_checked"]),
            u(&c["degenerate_targets"]),
            u(&c["triples_tested"]),
            u(&c["S4_roots_total"]),
            u(&c["rational_decompositions_total"]),
            u(&c["decompositions_that_are_not_roots"]),
            u(&c["extra_roots_F1"]),
            u(&c["extra_roots_F2"]),
            u(&c["extra_roots_unexplained"]),
            u(&c["extra_roots_at_nondegenerate_targets"]),
            u(&c["targets_where_roots_equal_decompositions"]),
            u(&c["targets_with_extra_roots"])
        );
    }

    // ---- the plotting data: derived here, laid out by plot_ps4.py ------------------
    if let Some(path) = arg_s(args, "--data") {
        let mut sj = Vec::new();
        let mut lj = Vec::new();
        for ((n, w), r) in &rows {
            let nt = r.ntargets as f64;
            sj.push(serde_json::json!({
                "n": n, "w": w, "targets": r.ntargets,
                "pair_factor_ord": ratio(r.u_ord.pairs as f64, r.r_rep.pairs as f64),
                "pair_factor_unord": ratio(r.u_unord.pairs as f64, r.r_rep.pairs as f64),
                "M_factor_ord": ratio(r.u_ord.mulq() as f64, r.r_rep.mulq() as f64),
                "M_factor_unord": ratio(r.u_unord.mulq() as f64, r.r_rep.mulq() as f64),
                "verdict_kept_frac": r.verdict_same as f64 / nt,
                "closure_eq_frac": r.closure_eq as f64 / (r.closure_n.max(1) as f64),
                "kept_sols_frac": ratio(r.sols_rep as f64, r.sols_total as f64),
                "dec_fw_ratio": ratio(r.dec_r_fw.mulq() as f64, r.dec_u_fw.mulq() as f64),
                "closure_M_over_ord": ratio(r.conj_counts.mulq() as f64, r.u_ord.mulq() as f64),
            }));
            let a = &r.u_unord;
            let m_pair = ratio(a.mulq() as f64, a.pairs as f64);
            let bits_pair = ratio(a.bits(), a.pairs as f64);
            let mut near: Option<u64> = None;
            for (l, c) in &rungs {
                let dn = (u(&c["n"]) as i64 - *n as i64).abs();
                if near.is_none_or(|b| dn < (u(&rungs[&b]["n"]) as i64 - *n as i64).abs()) {
                    near = Some(*l);
                }
            }
            if let Some(l) = near {
                let c = &rungs[&l];
                let nl = fl(&c["n"]);
                let ctl_bits = fl(&c["model"]["bitops_per_pair"]);
                lj.push(serde_json::json!({
                    "n": n, "w": w, "near_l": l,
                    "M_pair_onb": m_pair,
                    "M_pair_ctl": fl(&c["per_pair"]["mulq"]),
                    "bits_ratio": bits_pair / ctl_bits,
                    "bits_ratio_scaled": bits_pair / (ctl_bits * (*n as f64 / nl).powi(2)),
                }));
            }
        }
        let mut tj = Vec::new();
        for ((n, w), b) in &tb {
            let (mut pu, mut pt, mut mu, mut mt) = (0u64, 0u64, 0u64, 0u64);
            for t in &tors {
                if t["kind"].as_str() == Some("torsion_target")
                    && u(&t["n"]) as usize == *n
                    && u(&t["w"]) as usize == *w
                {
                    pu += u(&t["loops"]["U_ord_Tclosed"]["full"]["pairs"]);
                    pt += u(&t["loops"]["T_rep"]["full"]["pairs"]);
                    mu += counts_of(&t["loops"]["U_ord_Tclosed"]["full"]["counts"]).mulq();
                    mt += counts_of(&t["loops"]["T_rep"]["full"]["counts"]).mulq();
                }
            }
            tj.push(serde_json::json!({
                "n": n, "w": w, "m": u(&b["liftable_x"]),
                "pair_factor": ratio(pu as f64, pt as f64),
                "M_factor": ratio(mu as f64, mt as f64),
            }));
        }
        let doc = serde_json::json!({"settings": sj, "leaf": lj, "torsion": tj});
        std::fs::write(&path, serde_json::to_string_pretty(&doc).unwrap()).unwrap();
    }
}

//! Verdict analysis: reads the recorded campaign raw files, the replay and
//! control outputs, and writes results.json (all numbers derived here).
use crate::engine::*;
use crate::json::{dump, n as jn, obj, parse, st, V};
use std::collections::BTreeMap;

const B: usize = 10000;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: usize) -> usize { (self.next() % n as u64) as usize }
}

fn read(p: &str) -> V { parse(&std::fs::read_to_string(p).unwrap()) }
fn mean(x: &[f64]) -> f64 { if x.is_empty() { f64::NAN } else { x.iter().sum::<f64>() / x.len() as f64 } }
fn pct(v: &mut Vec<f64>, p: f64) -> f64 {
    v.retain(|x| x.is_finite());
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    if v.is_empty() { return f64::NAN; }
    let i = ((v.len() - 1) as f64 * p).round() as usize;
    v[i]
}
fn ci(v: &mut Vec<f64>) -> (f64, f64) { (pct(v, 0.025), pct(v, 0.975)) }

/// one cluster-bootstrap draw of the pooled mean over seed clusters
fn draw(cl: &[Vec<f64>], r: &mut Rng) -> f64 {
    let mut all = vec![];
    for _ in 0..cl.len() {
        let c = &cl[r.below(cl.len())];
        for _ in 0..c.len() { all.push(c[r.below(c.len())]); }
    }
    mean(&all)
}
fn pooled(cl: &[Vec<f64>]) -> Vec<f64> { cl.iter().flatten().cloned().collect() }

fn pearson(x: &[f64], y: &[f64]) -> f64 {
    let (mx, my) = (mean(x), mean(y));
    let mut sxy = 0.0; let mut sxx = 0.0; let mut syy = 0.0;
    for i in 0..x.len() { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx).powi(2); syy += (y[i] - my).powi(2); }
    sxy / (sxx * syy).sqrt()
}
fn ols_slope(x: &[f64], y: &[f64]) -> f64 {
    let (mx, my) = (mean(x), mean(y));
    let mut sxy = 0.0; let mut sxx = 0.0;
    for i in 0..x.len() { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx).powi(2); }
    sxy / sxx
}
fn civ(lo: f64, hi: f64) -> V { V::Arr(vec![jn(lo), jn(hi)]) }

#[derive(Default, Clone)]
struct Cell { seeds: Vec<Vec<f64>>, censored: usize, frac: Vec<f64>, l: usize, r: usize, wrel: Vec<f64> }

pub fn run(raw: &str, outdir: &str, out: &str) {
    let mut rng = Rng(20261015);
    let mut res: Vec<(String, V)> = vec![];

    // ── replay summary ──
    let mut rp = (0usize, 0usize, 0usize, 0usize, 0usize, 0usize, 0usize); let mut nshape = 0usize; // cells, structs, targets, matches, e1, meta_ok, unsat_total
    let mut replayed_cells = vec![];
    let mut rdir: Vec<_> = std::fs::read_dir(format!("{}/replay", outdir)).unwrap().filter_map(|e| e.ok()).map(|e| e.path()).collect();
    rdir.sort();
    for p in &rdir {
        let v = read(p.to_str().unwrap());
        rp.0 += 1;
        replayed_cells.push(st(p.file_stem().unwrap().to_str().unwrap()));
        for s in v.get("results").unwrap().arr().unwrap() {
            rp.1 += 1;
            rp.2 += s.get("replayed").unwrap().num().unwrap() as usize;
            rp.3 += s.get("replay_matches").unwrap().num().unwrap() as usize;
            rp.4 += s.get("e1_ok").unwrap().num().unwrap() as usize;
            rp.6 += s.get("unsat_targets").unwrap().num().unwrap() as usize;
            nshape += s.get("shape_recorded").and_then(|x| x.num()).unwrap_or(0.0) as usize;
            let ok = ["basis_match", "unsat_list_match", "fraction_match", "points_match"].iter()
                .all(|k| matches!(s.get(k), Some(V::Bool(true))));
            if ok { rp.5 += 1; }
        }
    }
    let mut delta_rep = vec![];
    if let Ok(rd) = std::fs::read_dir(format!("{}/delta", outdir)) {
        let mut ps: Vec<_> = rd.filter_map(|e| e.ok()).map(|e| e.path()).collect();
        ps.sort();
        for p in ps {
            let v = read(p.to_str().unwrap());
            let nn = v.get("field_degree").unwrap().num().unwrap();
            for r in v.get("rows").unwrap().arr().unwrap() {
                let ds: Vec<f64> = r.get("d_star").unwrap().arr().unwrap().iter().filter_map(|x| x.num()).collect();
                let mut bs: Vec<f64> = (0..B).map(|_| { let c = vec![ds.clone()]; draw(&c, &mut rng) }).collect();
                let (lo, hi) = ci(&mut bs);
                delta_rep.push(obj(vec![("n", jn(nn)), ("delta", r.get("delta").unwrap().clone()), ("dim", r.get("dim").unwrap().clone()),
                    ("d_star_mean", jn(mean(&ds))), ("ci95", civ(lo, hi)), ("unsat", jn(ds.len() as f64)),
                    ("mean_match_recorded", r.get("mean_match").unwrap().clone()),
                    ("fraction_match_recorded", r.get("fraction_match").unwrap().clone())]));
            }
        }
    }
    res.push(("replay".into(), obj(vec![
        ("cells", jn(rp.0 as f64)), ("structure_arms", jn(rp.1 as f64)), ("unsat_targets_in_cells", jn(rp.6 as f64)),
        ("targets_replayed", jn(rp.2 as f64)), ("targets_identical_on_recorded_fields", jn(rp.3 as f64)), ("targets_with_recorded_macaulay_shape", jn(nshape as f64)),
        ("arms_with_bases_targets_fraction_points_identical", jn(rp.5 as f64)),
        ("e1_invariants_hold", jn(rp.4 as f64)), ("replayed_cells", V::Arr(replayed_cells)),
        ("delta_replay", V::Arr(delta_rep)),
    ])));

    // ── recorded m=2 grid: build cells (n, structure, dim) ──
    let mut cells: BTreeMap<(usize, String, usize), Cell> = BTreeMap::new();
    let mut files: Vec<_> = std::fs::read_dir(raw).unwrap().filter_map(|e| e.ok()).map(|e| e.path())
        .filter(|p| { let s = p.file_name().unwrap().to_str().unwrap(); s.starts_with("m2_") && s.ends_with(".json") }).collect();
    files.sort();
    let mut fields: BTreeMap<usize, Field> = BTreeMap::new();
    for p in &files {
        let name = p.file_stem().unwrap().to_str().unwrap().to_string();
        let seed: u64 = name.rsplit("_s").next().unwrap().parse().unwrap();
        let v = read(p.to_str().unwrap());
        let n = v.get("field_degree").unwrap().num().unwrap() as usize;
        let dim = v.get("dim").unwrap().num().unwrap() as usize;
        let f = fields.entry(n).or_insert_with(|| Field::new(n)).clone();
        for r in v.get("results").unwrap().arr().unwrap() {
            if r.get("skipped").is_some() { continue; }
            let s = r.get("structure").unwrap().str().unwrap().to_string();
            let ds: Vec<f64> = r.get("targets").unwrap().arr().unwrap().iter().filter_map(|t| t.get("d_star").and_then(|x| x.num())).collect();
            let c = cells.entry((n, s.clone(), dim)).or_default();
            c.censored += r.get("censored_above_dmax").unwrap().num().unwrap() as usize;
            c.frac.push(r.get("decomposable_fraction").unwrap().num().unwrap());
            // exact invariants from the recorded bases
            let b1: Vec<u64> = r.get("basis").unwrap().arr().unwrap().iter().map(|x| u64::from_str_radix(x.str().unwrap().trim_start_matches("0x"), 16).unwrap()).collect();
            let b2: Vec<u64> = match r.get("basis2") { Some(V::Arr(a)) => a.iter().map(|x| u64::from_str_radix(x.str().unwrap().trim_start_matches("0x"), 16).unwrap()).collect(), _ => b1.clone() };
            let (rv, _) = product_dim(&f, &b1, &b2);
            let mut sm = b1.clone(); sm.extend_from_slice(&b2);
            c.l = rank_u64(&sm); c.r = c.r.max(rv);
            // work proxy at D* per target
            let fr = r.get("decomposable_fraction").unwrap().num().unwrap();
            for t in r.get("targets").unwrap().arr().unwrap() {
                if let (Some(d), Some(sh)) = (t.get("d_star").and_then(|x| x.num()), t.get("macaulay_shape")) {
                    let a = sh.get(&format!("{}", d as i64)).unwrap().arr().unwrap();
                    let w = a[0].num().unwrap() * a[1].num().unwrap();
                    if fr > 0.0 { c.wrel.push(w / fr); }
                }
            }
            let _ = seed;
            c.seeds.push(ds);
        }
    }
    let mut cell_rows = vec![];
    let mut cell_means: BTreeMap<(usize, String, usize), (f64, f64, f64)> = BTreeMap::new();
    for ((n, s, d), c) in &cells {
        let all = pooled(&c.seeds);
        let mut bs: Vec<f64> = (0..B).map(|_| draw(&c.seeds, &mut rng)).collect();
        let (lo, hi) = ci(&mut bs);
        cell_means.insert((*n, s.clone(), *d), (mean(&all), lo, hi));
        cell_rows.push(obj(vec![
            ("n", jn(*n as f64)), ("structure", st(s)), ("dim", jn(*d as f64)), ("nvars", jn(2.0 * *d as f64)),
            ("seeds", jn(c.seeds.len() as f64)), ("unsat_targets", jn(all.len() as f64)), ("censored", jn(c.censored as f64)),
            ("d_star_mean", jn(mean(&all))), ("ci95", civ(lo, hi)),
            ("decomposable_fraction_mean", jn(mean(&c.frac))),
            ("l_dim_V1_plus_V2", jn(c.l as f64)), ("r_dim_V1V2_max_over_seeds", jn(c.r as f64)),
            ("slack_n_minus_l_minus_r", jn(*n as f64 - c.l as f64 - c.r as f64)),
            ("wrel_mean", jn(mean(&c.wrel))),
        ]));
    }
    res.push(("cells".into(), V::Arr(cell_rows)));

    // ── H1 contrasts vs random at matched (n, dim) + N1 seed-to-seed null ──
    let mut h1 = vec![];
    let mut n1 = vec![];
    for ((n, s, d), c) in &cells {
        if s == "random" {
            let k = c.seeds.len();
            for i in 0..k { for j in (i + 1)..k {
                n1.push(mean(&c.seeds[i]) - mean(&c.seeds[j]));
            } }
            continue;
        }
        if let Some(cr) = cells.get(&(*n, "random".to_string(), *d)) {
            let diff = mean(&pooled(&c.seeds)) - mean(&pooled(&cr.seeds));
            let mut bs: Vec<f64> = (0..B).map(|_| draw(&c.seeds, &mut rng) - draw(&cr.seeds, &mut rng)).collect();
            let (lo, hi) = ci(&mut bs);
            h1.push(obj(vec![("n", jn(*n as f64)), ("dim", jn(*d as f64)), ("structure", st(s)),
                ("diff_vs_random", jn(diff)), ("ci95", civ(lo, hi)),
                ("r_struct", jn(c.r as f64)), ("r_random", jn(cr.r as f64)), ("l_struct", jn(c.l as f64)), ("l_random", jn(cr.l as f64))]));
        }
    }
    let mut n1s = n1.clone();
    let n1_abs95 = { let mut a: Vec<f64> = n1.iter().map(|x| x.abs()).collect(); pct(&mut a, 0.95) };
    res.push(("h1_contrasts".into(), V::Arr(h1)));
    res.push(("n1_random_seed_contrasts".into(), obj(vec![("count", jn(n1.len() as f64)), ("abs_q95", jn(n1_abs95)),
        ("q025", jn(pct(&mut n1s.clone(), 0.025))), ("q975", jn(pct(&mut n1s, 0.975))), ("values", V::Arr(n1.iter().map(|x| jn(*x)).collect()))])));

    // ── H2 slope of D* vs log2(dim) per (n, structure) ──
    let mut h2 = vec![];
    let mut groups: BTreeMap<(usize, String), Vec<(usize, Vec<Vec<f64>>)>> = BTreeMap::new();
    for ((n, s, d), c) in &cells { groups.entry((*n, s.clone())).or_default().push((*d, c.seeds.clone())); }
    for ((n, s), v) in &groups {
        if v.len() < 2 { continue; }
        let fit = |v: &Vec<(usize, Vec<f64>)>| {
            let mut x = vec![]; let mut y = vec![];
            for (d, ds) in v { for t in ds { x.push((*d as f64).log2()); y.push(*t); } }
            ols_slope(&x, &y)
        };
        let base: Vec<(usize, Vec<f64>)> = v.iter().map(|(d, cl)| (*d, pooled(cl))).collect();
        let sl = fit(&base);
        let mut bs: Vec<f64> = (0..2000).map(|_| {
            let bv: Vec<(usize, Vec<f64>)> = v.iter().map(|(d, cl)| {
                let p = pooled(cl);
                (*d, (0..p.len()).map(|_| p[rng.below(p.len())]).collect())
            }).collect();
            fit(&bv)
        }).collect();
        let (lo, hi) = ci(&mut bs);
        let dims: Vec<f64> = v.iter().map(|x| x.0 as f64).collect();
        h2.push(obj(vec![("n", jn(*n as f64)), ("structure", st(s)), ("dims", V::Arr(dims.iter().map(|x| jn(*x)).collect())),
            ("slope_per_doubling", jn(sl)), ("ci95", civ(lo, hi)), ("abandon_rule_fires", V::Bool(lo > 1.0))]));
    }
    res.push(("h2_slopes".into(), V::Arr(h2)));

    // P2 positive control: n=11 random dim5 - dim2
    if let (Some(a), Some(b)) = (cells.get(&(11, "random".into(), 5)), cells.get(&(11, "random".into(), 2))) {
        let diff = mean(&pooled(&a.seeds)) - mean(&pooled(&b.seeds));
        let mut bs: Vec<f64> = (0..B).map(|_| draw(&a.seeds, &mut rng) - draw(&b.seeds, &mut rng)).collect();
        let (lo, hi) = ci(&mut bs);
        res.push(("p2_dimension_control".into(), obj(vec![("diff_dim5_minus_dim2", jn(diff)), ("ci95", civ(lo, hi))])));
    }

    // ── H5 correlations over non-subfield cell means ──
    let mut xs_d = vec![]; let mut xs_n = vec![]; let mut ys = vec![]; let mut xs_k = vec![];
    for ((n, s, d), c) in &cells {
        if s == "subfield" { continue; }
        xs_d.push(*d as f64); xs_n.push(*n as f64); ys.push(mean(&pooled(&c.seeds)));
        xs_k.push((c.l + c.r) as f64 / *n as f64);
    }
    let m = ys.len();
    let boot_r = |x: &Vec<f64>, rng: &mut Rng| {
        let mut bs: Vec<f64> = (0..B).map(|_| {
            let idx: Vec<usize> = (0..m).map(|_| rng.below(m)).collect();
            pearson(&idx.iter().map(|&i| x[i]).collect::<Vec<_>>(), &idx.iter().map(|&i| ys[i]).collect::<Vec<_>>())
        }).collect();
        ci(&mut bs)
    };
    let (a1, b1) = boot_r(&xs_d, &mut rng);
    let (a2, b2) = boot_r(&xs_n, &mut rng);
    let (a3, b3) = boot_r(&xs_k, &mut rng);
    res.push(("h5_correlations".into(), obj(vec![
        ("cells", jn(m as f64)),
        ("r_dstar_dim", jn(pearson(&xs_d, &ys))), ("r_dstar_dim_ci95", civ(a1, b1)),
        ("r_dstar_n", jn(pearson(&xs_n, &ys))), ("r_dstar_n_ci95", civ(a2, b2)),
        ("r_dstar_(l+r)/n", jn(pearson(&xs_k, &ys))), ("r_dstar_(l+r)/n_ci95", civ(a3, b3)),
    ])));

    // ── H3 recorded delta sweeps ──
    let mut h3 = vec![];
    for nn in [11, 13, 15, 17, 19] {
        let p = format!("{}/delta_n{}_m2.json", raw, nn);
        if let Ok(t) = std::fs::read_to_string(&p) {
            let v = parse(&t);
            let rows: Vec<&V> = v.get("rows").unwrap().arr().unwrap().iter().filter(|r| r.get("skipped").is_none()).collect();
            let xs: Vec<f64> = rows.iter().map(|r| r.get("delta").unwrap().num().unwrap()).collect();
            let ys: Vec<f64> = rows.iter().map(|r| r.get("d_star_mean").and_then(|x| x.num()).unwrap_or(f64::NAN)).collect();
            let w0 = rows[0].get("work_per_relation_proxy").and_then(|x| x.num());
            let pts: Vec<V> = rows.iter().map(|r| obj(vec![
                ("delta", r.get("delta").unwrap().clone()), ("dim", r.get("dim").unwrap().clone()),
                ("d_star_mean", r.get("d_star_mean").cloned().unwrap_or(V::Null)),
                ("decomposable_fraction", r.get("decomposable_fraction").unwrap().clone()),
                ("unsat_measured", r.get("unsat_measured").unwrap().clone()),
                ("wrel", r.get("work_per_relation_proxy").cloned().unwrap_or(V::Null)),
                ("wrel_gain_vs_delta_m1", match (w0, r.get("work_per_relation_proxy").and_then(|x| x.num())) { (Some(a), Some(b)) => jn(a / b), _ => V::Null }),
            ])).collect();
            h3.push(obj(vec![("n", jn(nn as f64)), ("dstar_drop_per_plus2delta", jn(-2.0 * ols_slope(&xs, &ys))), ("rows", V::Arr(pts))]));
        }
    }
    res.push(("h3_delta".into(), V::Arr(h3)));

    // ── H4 sym4 recorded ──
    let mut h4 = vec![];
    for nm in ["sym4_n9_l2", "sym4_n15_l2", "sym4_n15_l3", "sym4_n19_l3", "sym4_n21_l3"] {
        if let Ok(t) = std::fs::read_to_string(format!("{}/{}.json", raw, nm)) {
            let v = parse(&t);
            let mut arms = vec![];
            for r in v.get("results").unwrap().arr().unwrap() {
                if r.get("skipped").is_some() { continue; }
                let ds: Vec<f64> = r.get("targets").map(|x| x.arr().unwrap().iter().filter_map(|t| t.get("d_star").and_then(|y| y.num())).collect()).unwrap_or_default();
                arms.push(obj(vec![("structure", r.get("structure").unwrap().clone()), ("nvars", r.get("nvars").cloned().unwrap_or(V::Null)),
                    ("d_star", V::Arr(ds.iter().map(|x| jn(*x)).collect())), ("d_star_mean", jn(mean(&ds))),
                    ("censored", r.get("censored_above_dmax").cloned().unwrap_or(V::Null)),
                    ("point_decomposable_fraction", r.get("point_decomposable_fraction").cloned().unwrap_or(V::Null))]));
            }
            h4.push(obj(vec![("cell", st(nm)), ("arms", V::Arr(arms))]));
        }
    }
    res.push(("h4_sym4".into(), V::Arr(h4)));

    // ── N2 invariant null vs Semaev (and exploratory symmetric-linear variant N2s) ──
    for (fname, key) in [("null.json", "n2_invariant_null"), ("null_sym.json", "n2s_exploratory_symmetric_linear_null")] {
    if let Ok(t) = std::fs::read_to_string(format!("{}/{}", outdir, fname)) {
        let v = parse(&t);
        let mut rows = vec![];
        let mut within = 0; let mut total = 0;
        for c in v.get("cells").unwrap().arr().unwrap() {
            let n = c.get("field_degree").unwrap().num().unwrap() as usize;
            let d = c.get("dim").unwrap().num().unwrap() as usize;
            let s = c.get("structure").unwrap().str().unwrap().to_string();
            let nd: Vec<f64> = c.get("d_star").unwrap().arr().unwrap().iter().filter_map(|x| x.num()).collect();
            let ncens = c.get("d_star").unwrap().arr().unwrap().iter().filter(|x| matches!(x, V::Null)).count();
            let Some(sc) = cells.get(&(n, s.clone(), d)) else { continue };
            let sm = mean(&pooled(&sc.seeds));
            let diff = sm - mean(&nd);
            let ncl = vec![nd.clone()];
            let mut bs: Vec<f64> = (0..B).map(|_| draw(&sc.seeds, &mut rng) - draw(&ncl, &mut rng)).collect();
            let (lo, hi) = ci(&mut bs);
            total += 1; if diff.abs() <= 0.25 { within += 1; }
            rows.push(obj(vec![("n", jn(n as f64)), ("dim", jn(d as f64)), ("structure", st(&s)),
                ("l", c.get("l").unwrap().clone()), ("q", c.get("q").unwrap().clone()), ("modal_share", c.get("modal_share").unwrap().clone()),
                ("semaev_mean", jn(sm)), ("null_mean", jn(mean(&nd))), ("null_censored", jn(ncens as f64)),
                ("diff", jn(diff)), ("ci95", civ(lo, hi))]));
        }
        let mut pq = vec![];
        for c in v.get("planted_q_control").unwrap().arr().unwrap() {
            let nd: Vec<f64> = c.get("d_star").unwrap().arr().unwrap().iter().filter_map(|x| x.num()).collect();
            let ncl = vec![nd.clone()];
            let mut bs: Vec<f64> = (0..B).map(|_| draw(&ncl, &mut rng)).collect();
            let (lo, hi) = ci(&mut bs);
            pq.push(obj(vec![("q", c.get("q").unwrap().clone()), ("d_star_mean", jn(mean(&nd))), ("ci95", civ(lo, hi)), ("count", jn(nd.len() as f64))]));
        }
        res.push((key.into(), obj(vec![("cells", jn(total as f64)), ("within_0.25", jn(within as f64)),
            ("share_within", jn(within as f64 / total.max(1) as f64)), ("hypothesis_I_pass_rule", st(">= 0.80")),
            ("rows", V::Arr(rows)), ("planted_q_control", V::Arr(pq))])));
    }
    }
    if let Ok(t) = std::fs::read_to_string(format!("{}/e2.json", outdir)) { res.push(("e2_product_dims".into(), parse(&t))); }

    std::fs::write(out, dump(&V::Obj(res)) + "\n").unwrap();
}

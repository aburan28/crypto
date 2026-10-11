//! Point-level validation of isogeny-loop endomorphisms on GLS twists.
//!
//! For each discriminant D: find a GLS instance (`p = (t^2 + |D| v^2)/4`, 127-bit by default,
//! `p = 3 mod 4`, `#E'(F_{p^2}) = h r`), build `E/F_p` from a root of the natively computed
//! Hilbert class polynomial `H_D`, walk every closing non-backtracking horizontal loop with the
//! cheapest loop's degree sequence, and check on points:
//!
//! 1. on `E(F_p)[r_E]` (largest prime factor `r_E` of `p + 1 - t`, when it is large) the loop
//!    acts as one of `+-lambda_E(+-s_E)` with `s_E = (2 - t)/v`, `lambda = (A + B s)/2` (pi = 1);
//! 2. on the quadratic twist `E'/F_{p^2}` (by a non-square `u`), the GLS map
//!    `psi(x, y) = (u^{1-p} x^p, u^{3(1-p)/2} y^p)` acts as `mu` or `-mu`, `mu = (p - 1)/t mod r`;
//! 3. on `E'(F_{p^2})[r]` the twisted loop acts as one of `+-lambda(+-s)` with
//!    `s = (2 mu - t)/v mod r` — the eigenvalue the lattice study
//!    (`research/isogeny_loop_glv_20261011/`) derived from the ring map `pi -> mu` but never
//!    checked on points.
//!
//! PASS requires, per closing loop, the same candidate on every sample point, and the two
//! closing loops on conjugate embeddings (`s` and `-s`).
//!
//! ```bash
//! cargo run --release --example gls_loop_validate -- [--pbits 127] [--points 24] \
//!     [--seed 20261012] [--discs -619,-23,-71,-431,-3019] [--json out.json]
//! ```
use isogeny_algos::api::{self, PrimeFieldInt};
use isogeny_algos::bigint::Big;
use isogeny_algos::curve::{pmul, pmul_big, random_point_f, Curve, Pt};
use isogeny_algos::field::{Field, Rng};
use isogeny_algos::fp2::Fp2;
use isogeny_algos::fpm::FpM;
use isogeny_algos::gls_loop::{build_curve, cheapest_loop, closing_loops, find_instance};
use isogeny_algos::int::Int;
use isogeny_algos::json::Json;

type F = FpM<2>;
type E2 = <Fp2<F> as Field>::E;

fn pow_big<G: Field>(g: &G, a: G::E, e: &Big) -> G::E {
    let mut r = g.one();
    for i in (0..e.bits()).rev() {
        r = g.mul(r, r);
        if e.bit(i) {
            r = g.mul(r, a);
        }
    }
    r
}

fn modr(x: &Int, r: &Int) -> Int {
    x.modulo(r)
}

/// The four candidate eigenvalues +-lambda(s), +-lambda(-s) for alpha = (A + B sqrt D)/2.
fn candidates(a: i64, b: i64, s: &Int, r: &Int) -> Vec<(String, Int)> {
    let inv2 = Int::from(2i64).inv_mod(r).unwrap();
    let mut out = vec![];
    for (sname, ss) in [("s", s.clone()), ("-s", (-s).modulo(r))] {
        let lam = modr(&(&(&Int::from(a) + &(&Int::from(b) * &ss)) * &inv2), r);
        out.push((format!("+lambda({sname})"), lam.clone()));
        out.push((format!("-lambda({sname})"), (-&lam).modulo(r)));
    }
    out
}

fn arg<T: std::str::FromStr>(args: &[String], name: &str, default: T) -> T {
    args.iter().position(|a| a == name).and_then(|i| args.get(i + 1)).and_then(|v| v.parse().ok()).unwrap_or(default)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let pbits: usize = arg(&args, "--pbits", 127);
    let npts: usize = arg(&args, "--points", 24);
    let seed: u64 = arg(&args, "--seed", 20_261_012);
    let discs: Vec<i64> = arg(&args, "--discs", "-619,-23,-71,-431,-3019".to_string())
        .split(',')
        .map(|s| s.trim().parse().unwrap())
        .collect();
    let json_out: Option<String> = args.iter().position(|a| a == "--json").and_then(|i| args.get(i + 1)).cloned();
    assert!(pbits <= 127, "this example uses FpM<2>; pbits must be <= 127");

    let mut all_pass = true;
    let mut rows = vec![];
    println!("# Point-level validation of loop endomorphisms on GLS twists\n");
    println!("p bits {pbits}, {npts} sample points per check, seed {seed}\n");
    println!("| D | loop | p | h | log2 r | H_D degree / precision bits / worst rounding log2 | closing loops | F_p check (r_E bits) | psi acts as | twisted loop acts as (per closing loop) | verdict |");
    println!("|--:|:--|:--|--:|--:|:--|--:|:--|:--|:--|:--|");
    for (di, &d) in discs.iter().enumerate() {
        let s_d = seed ^ ((di as u64 + 1) << 32);
        let (lp, _, _, _) = cheapest_loop(d);
        let degs = lp.degree_sequence();
        let avoid: Vec<u64> = lp.steps.iter().map(|&(l, _)| l).filter(|&l| l > 2).collect();
        let inst = find_instance(d, pbits, &avoid, 64, s_d, 20_000_000).expect("no instance");
        let cm = build_curve(&inst, s_d).expect("curve");
        let chains = closing_loops(&cm.curve, &cm.crater, &degs, s_d).expect("walk");
        let p = &inst.p;
        let f = F::from_dec(&p.to_string());
        let r = &inst.r;

        // ---- 1. F_p check on E(F_p)[r_E]
        let ne = &(p + &Int::one()) - &inst.t;
        let mut rng_f = isogeny_algos::field::Rng::new(s_d ^ 0xF00D);
        let (fac, rest) = api::factor_partial(&ne, 50_000_000, &mut rng_f);
        let mut r_e = rest.clone();
        for (q, _) in &fac {
            if *q > r_e {
                r_e = q.clone();
            }
        }
        let mut fp_note = String::from("skipped (no large prime factor)");
        let mut fp_ok = true;
        if r_e.bits() >= 40 && r_e.is_probable_prime() {
            let cof = ne.div_floor(&r_e);
            let c = Curve::new(f.elem(&cm.curve.a), f.elem(&cm.curve.b));
            let s_e = (&(&Int::from(2i64) - &inst.t) * &inst.v.inv_mod(&r_e).unwrap()).modulo(&r_e);
            let cands = candidates(lp.a, lp.b, &s_e, &r_e);
            let mut seen: Vec<Vec<String>> = vec![vec![]; chains.len()];
            let mut rng = Rng::new(s_d ^ 0xBEEF);
            for _ in 0..npts {
                let pt = pmul_big(&f, &c, &random_point_f(&f, &c, &mut rng), cof.mag());
                let Pt::Aff(x, y) = pt else { continue };
                for (ci, ch) in chains.iter().enumerate() {
                    let fc = ch.to_field(&f, &|v: &Int| f.elem(v));
                    let img = fc.eval_twisted(&f, f.one(), x, y).map(|(a, b)| Pt::Aff(a, b)).unwrap_or(Pt::Inf);
                    let hit: Vec<String> = cands
                        .iter()
                        .filter(|(_, l)| pmul_big(&f, &c, &pt, l.mag()) == img)
                        .map(|(n, _)| n.clone())
                        .collect();
                    seen[ci].push(hit.join("|"));
                }
            }
            let summary: Vec<String> = seen
                .iter()
                .map(|v| {
                    let first = v.first().cloned().unwrap_or_default();
                    if !first.is_empty() && v.iter().all(|x| *x == first) {
                        first
                    } else {
                        fp_ok = false;
                        format!("INCONSISTENT {v:?}")
                    }
                })
                .collect();
            fp_note = format!("{} ({} bits)", summary.join(", "), r_e.bits());
        }

        // ---- 2./3. GLS twist over F_{p^2}
        let g = Fp2::new(f.clone());
        let mut rng = Rng::new(s_d ^ 0xCAFE);
        let p2m1_half = {
            let pp = p.mag().mul(p.mag());
            pp.sub_small(1).shr(1)
        };
        let u: E2 = loop {
            let cand = g.random(&mut rng);
            if !g.is_zero(cand) && pow_big(&g, cand, &p2m1_half) != g.one() {
                break cand;
            }
        };
        let a0 = g.embed(f.elem(&cm.curve.a));
        let b0 = g.embed(f.elem(&cm.curve.b));
        let u2 = g.mul(u, u);
        let u3 = g.mul(u2, u);
        let ct = Curve::new(g.mul(u2, a0), g.mul(u3, b0));
        // psi constants: c1 = u / conj(u), c2 = (u^((p-1)/2))^-3
        let c1 = g.mul(u, g.inv(g.conj(u)));
        let pm1_half = p.mag().sub_small(1).shr(1);
        let uh = pow_big(&g, u, &pm1_half);
        let uhi = g.inv(uh);
        let c2 = g.mul(g.mul(uhi, uhi), uhi);
        let psi = |pt: &Pt<E2>| -> Pt<E2> {
            match *pt {
                Pt::Inf => Pt::Inf,
                Pt::Aff(x, y) => Pt::Aff(g.mul(c1, g.conj(x)), g.mul(c2, g.conj(y))),
            }
        };
        let mu = (&(p - &Int::one()) * &inst.t.inv_mod(r).unwrap()).modulo(r);
        assert!((&(&mu * &mu) + &Int::one()).modulo(r).is_zero(), "mu^2 != -1");
        let s = (&(&(&Int::from(2i64) * &mu) - &inst.t) * &inst.v.inv_mod(r).unwrap()).modulo(r);
        assert!((&(&s * &s) - &Int::from(d)).modulo(r).is_zero(), "s^2 != D");
        let cands = candidates(lp.a, lp.b, &s, r);
        let mut psi_seen = vec![];
        let mut seen: Vec<Vec<String>> = vec![vec![]; chains.len()];
        let fcs: Vec<_> = chains.iter().map(|ch| ch.to_field(&g, &|v: &Int| g.embed(f.elem(v)))).collect();
        let mut killed = true;
        for _ in 0..npts {
            let q = pmul(&g, &ct, &random_point_f(&g, &ct, &mut rng), inst.h as u128);
            let Pt::Aff(x, y) = q else { continue };
            killed &= pmul_big(&g, &ct, &q, r.mag()) == Pt::Inf;
            let pq = psi(&q);
            let ps = if pq == pmul_big(&g, &ct, &q, mu.mag()) {
                "+mu"
            } else if pq == pmul_big(&g, &ct, &q, (-&mu).modulo(r).mag()) {
                "-mu"
            } else {
                "neither"
            };
            psi_seen.push(ps.to_string());
            for (ci, fc) in fcs.iter().enumerate() {
                let img = fc.eval_twisted(&g, u, x, y).map(|(a, b)| Pt::Aff(a, b)).unwrap_or(Pt::Inf);
                let hit: Vec<String> = cands
                    .iter()
                    .filter(|(_, l)| pmul_big(&g, &ct, &q, l.mag()) == img)
                    .map(|(n, _)| n.clone())
                    .collect();
                seen[ci].push(hit.join("|"));
            }
        }
        let psi_first = psi_seen.first().cloned().unwrap_or_default();
        let psi_ok = psi_first != "neither" && psi_seen.iter().all(|x| *x == psi_first);
        let mut tw_ok = !chains.is_empty();
        let tw: Vec<String> = seen
            .iter()
            .map(|v| {
                let first = v.first().cloned().unwrap_or_default();
                if !first.is_empty() && v.iter().all(|x| *x == first) {
                    first
                } else {
                    tw_ok = false;
                    format!("INCONSISTENT {v:?}")
                }
            })
            .collect();
        // the closing loops must sit on conjugate embeddings
        let embeds: Vec<&str> = tw.iter().map(|x| if x.contains("(-s)") { "-s" } else { "s" }).collect();
        let conj_ok = embeds.len() == 2 && embeds[0] != embeds[1];
        let pass = killed && psi_ok && tw_ok && fp_ok && conj_ok;
        all_pass &= pass;
        let loop_s = format!(
            "({} + {} sqrt({d}))/2, N = {}, steps {:?}",
            lp.a, lp.b, lp.norm, degs
        );
        println!(
            "| {d} | {loop_s} | {}… ({} bits) | {} | {} | {} / {} / {:.0} | {} | {} | {} ({}/{}) | {} | {} |",
            &p.to_string()[..12],
            p.bits(),
            inst.h,
            r.bits(),
            cm.hilbert_degree,
            cm.hilbert_precision,
            cm.hilbert_error_log2,
            chains.len(),
            fp_note,
            psi_first,
            psi_seen.iter().filter(|x| **x == psi_first).count(),
            psi_seen.len(),
            tw.join("; "),
            if pass { "PASS" } else { "FAIL" }
        );
        rows.push(Json::obj(vec![
            ("D", Json::Num(d)),
            ("loop", Json::obj(vec![
                ("A", Json::Num(lp.a)),
                ("B", Json::Num(lp.b)),
                ("norm", Json::str(lp.norm)),
                ("cost_M", Json::Float(lp.cost)),
                ("degrees", Json::Arr(degs.iter().map(|&l| Json::Num(l as i64)).collect())),
            ])),
            ("instance", Json::obj(vec![
                ("p", Json::str(p)),
                ("t", Json::str(&inst.t)),
                ("v", Json::str(&inst.v)),
                ("h", Json::Num(inst.h as i64)),
                ("r", Json::str(r)),
                ("tries", Json::Num(inst.tries as i64)),
            ])),
            ("curve", Json::obj(vec![("a", Json::str(&cm.curve.a)), ("b", Json::str(&cm.curve.b)), ("j", Json::str(cm.curve.j()))])),
            ("hilbert", Json::obj(vec![
                ("degree", Json::Num(cm.hilbert_degree as i64)),
                ("precision_bits", Json::Num(cm.hilbert_precision as i64)),
                ("worst_rounding_log2", Json::Float(cm.hilbert_error_log2)),
                ("roots_in_Fp", Json::Num(cm.crater.len() as i64)),
            ])),
            ("closing_loops", Json::Arr(chains.iter().map(|ch| Json::Arr(ch.js.iter().map(Json::str).collect())).collect())),
            ("fp_check", Json::str(&fp_note)),
            ("mu", Json::str(&mu)),
            ("psi", Json::Arr(psi_seen.iter().map(Json::str).collect())),
            ("twisted_loops", Json::Arr(seen.iter().map(|v| Json::Arr(v.iter().map(Json::str).collect())).collect())),
            ("order_r_killed_all_points", Json::Bool(killed)),
            ("pass", Json::Bool(pass)),
        ]));
    }
    println!(
        "\nverdict: {}",
        if all_pass { "ALL PASS" } else { "FAIL (see rows)" }
    );
    if let Some(path) = json_out {
        let doc = Json::obj(vec![
            ("schema", Json::str("gls_loop_validate/v1")),
            ("command", Json::str(args.join(" "))),
            ("pbits", Json::Num(pbits as i64)),
            ("points", Json::Num(npts as i64)),
            ("seed", Json::str(seed)),
            ("all_pass", Json::Bool(all_pass)),
            ("rows", Json::Arr(rows)),
        ]);
        std::fs::write(&path, doc.dump()).expect("write");
        eprintln!("wrote {path}");
    }
    if !all_pass {
        std::process::exit(1);
    }
}

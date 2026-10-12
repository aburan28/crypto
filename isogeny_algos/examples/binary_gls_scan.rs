//! Binary-field GLS with a CM loop endomorphism: existence scan, Frobenius-assisted loops, and
//! 4-dimensional decomposition lattices.
//!
//! Setting.  An ordinary `E/F_q`, `q = 2^m`, with `End(E) = O_K` of discriminant `D` needs 2 to
//! split (`D = 1 mod 8`), and its Frobenius `pi` generates `p_F^m` for the prime `p_F | 2`
//! containing `pi`.  So `[p_2]^m = 1` in `Cl(O_K)`; for `m` prime and `D != -7` the class of `p_2`
//! has order exactly `m`, hence `m | h(D)`.  The GLS twist `E'/F_{q^2}` has
//! `#E'(F_{q^2}) = (q - 1)^2 + t^2`, always `2 mod 4`, and `psi` with `psi^2 = -1`.
//!
//! Frobenius-assisted loops.  In characteristic 2 the steps along `p_F` are the Frobenius
//! `F: (x, y) -> (x^2, y^2)`, two squarings, and `F^k` is a GF(2)-linear map of the coordinates.
//! If an odd split prime `l` has `[l] = [p_F]^{-k}` (always the case when `Cl = <[p_2]>`), the walk
//! `l` then `F^k` is a loop: an endomorphism of degree `l 2^k` that costs one odd `l`-isogeny plus
//! `k` Frobenius steps, whatever `k` is.  Its eigenvalue on `E'(F_{q^2})[r]` follows from the ring
//! map `pi -> mu = (q - 1)/t` exactly as over prime fields.
//!
//! This example (number theory and lattices only; no curve is built):
//! 1. scans fundamental `D = 1 mod 8` with `|D| <= --dmax` for `[p_2]^m = 1`, `[p_2] != 1`;
//! 2. recovers `pi = (t + v sqrt D)/2` from `t^2 + |D| v^2 = 2^{m+2}` (Cornacchia, verified);
//! 3. keeps `D` where `((q-1)^2 + t^2)/2 = r` is prime;
//! 4. for each odd split `l <= --lmax` finds `k` with `[l][p_2]^k = 1`, the generator `alpha` of
//!    norm `l 2^min(k, m-k)` lying in `p_F^min(k, m-k)` (2-adic valuation check), and its cost;
//! 5. measures the GLS-2 and GLV-GLS-4 lattices (exact LLL, Babai, seeded scalars) for the
//!    cheapest loop, both ring embeddings.
//!
//! ```bash
//! cargo run --release --example binary_gls_scan -- [--m 127] [--dmax 4000000] [--lmax 61] \
//!     [--keep 6] [--samples 5000] [--seed 20261012] [--json out.json]
//! ```
use isogeny_algos::cm::reduced_forms;
use isogeny_algos::field::Rng;
use isogeny_algos::gls_loop::{
    cornacchia, gauss2, is_prime_u64, lll4, prime_form, sqrt_mod_2k, squarefree, v2, Babai,
};
use isogeny_algos::int::Int;
use isogeny_algos::json::Json;
use isogeny_algos::path::relation::Form;

fn arg<T: std::str::FromStr>(args: &[String], name: &str, default: T) -> T {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let m: usize = arg(&args, "--m", 127);
    let dmax: i64 = arg(&args, "--dmax", 4_000_000);
    let lmax: u64 = arg(&args, "--lmax", 61);
    let keep: usize = arg(&args, "--keep", 6);
    let samples: usize = arg(&args, "--samples", 5000);
    let seed: u64 = arg(&args, "--seed", 20_261_012);
    let json_out: Option<String> = args
        .iter()
        .position(|a| a == "--json")
        .and_then(|i| args.get(i + 1))
        .cloned();
    let q = Int::one().shl(m);
    let t0 = std::time::Instant::now();

    let mut scanned = 0u64;
    let mut ord_hits = 0u64;
    let mut prime_r_hits = 0u64;
    let mut kept: Vec<Json> = vec![];
    let mut lines: Vec<String> = vec![];
    let mut first_hits: Vec<(i64, usize)> = vec![];
    let mut d = -15i64; // D = 1 mod 8 and D != -7
    while -d <= dmax && kept.len() < keep {
        d -= 8;
        if !squarefree(-d) {
            continue;
        }
        scanned += 1;
        let Some(f2) = prime_form(d, 2) else { continue };
        let id = Form::identity(d as i128);
        if f2 == id || f2.pow(m as i128) != id {
            continue;
        }
        ord_hits += 1;
        // Frobenius pi from t^2 + |D| v^2 = 2^(m+2)
        let sols = cornacchia(d, 1, m);
        let Some((t, v)) = sols
            .iter()
            .find(|(a, b)| !a.is_even() && !b.is_zero())
            .cloned()
        else {
            lines.push(format!(
                "| {d} | — | Cornacchia found no odd-t solution | | | | | | |"
            ));
            continue;
        };
        let h = reduced_forms(d).len();
        if first_hits.len() < 12 {
            first_hits.push((d, h));
        }
        let qm1 = &q - &Int::one();
        let nt = &(&qm1 * &qm1) + &(&t * &t);
        let r = nt.shr(1);
        if r.shl(1) != nt || !r.is_probable_prime() {
            continue;
        }
        prime_r_hits += 1;
        // orientation: the 2-adic root s2 with (t + v s2)/2 = 0 mod 2^m
        let kbits = m + 8;
        let s2 = sqrt_mod_2k(&Int::from(d), kbits)
            .into_iter()
            .find(|s| v2(&(&t + &(&v * s))) > m)
            .expect("pi lies in p_F^m");
        // Frobenius-assisted loops for odd split l
        let mut loops = vec![];
        let mut l = 3u64;
        while l <= lmax {
            if is_prime_u64(l) {
                if let Some(fl) = prime_form(d, l as i64) {
                    {
                        // k with [l][p_2]^k = 1
                        let mut g = fl;
                        let mut kfound = None;
                        for k in 0..m {
                            if g == id {
                                kfound = Some(k);
                                break;
                            }
                            g = g.compose(&f2).reduce();
                        }
                        if let Some(k) = kfound {
                            let kk = k.min(m - k);
                            // generators of norm l 2^kk in p_F^kk
                            let cands = cornacchia(d, l, kk);
                            let good: Vec<(Int, Int)> = cands
                                .iter()
                                .flat_map(|(a, b)| [(a.clone(), b.clone()), (a.clone(), -b)])
                                .filter(|(a, b)| v2(&(a + &(b * &s2))) > kk)
                                .collect();
                            if let Some((a, b)) = good.first().cloned() {
                                loops.push((l, k, kk, a, b));
                            }
                        }
                    }
                }
            }
            l += 2;
        }
        // cheapest: 7.5 l for the odd step + 2 per F step (two squarings counted as 1 M each)
        let cost = |lp: &(u64, usize, usize, Int, Int)| 7.5 * lp.0 as f64 + 2.0 * lp.2 as f64;
        loops.sort_by(|x, y| cost(x).partial_cmp(&cost(y)).unwrap());
        let Some(best) = loops.first().cloned() else {
            lines.push(format!(
                "| {d} | {h} | no Frobenius-assisted loop with l <= {lmax} | | | | | | |"
            ));
            continue;
        };
        let (bl, bk, bkk, ba, bb) = best.clone();
        // lattices
        let mu = (&qm1 * &t.inv_mod(&r).unwrap()).modulo(&r);
        assert!((&(&mu * &mu) + &Int::one()).modulo(&r).is_zero());
        let s = (&(&(&Int::from(2i64) * &mu) - &t) * &v.inv_mod(&r).unwrap()).modulo(&r);
        assert!((&(&s * &s) - &Int::from(d)).modulo(&r).is_zero());
        let inv2 = Int::from(2i64).inv_mod(&r).unwrap();
        let mut rng = Rng::new(seed ^ (-d) as u64);
        let rand_scalar = |rng: &mut Rng| {
            let mut x = Int::zero();
            for _ in 0..((m * 2) / 64 + 2) {
                x = &x.shl(64) + &Int::from(rng.next());
            }
            x.modulo(&r)
        };
        let b2 = Babai::new(gauss2(&r, &mu), r.bits());
        let mut g2max = 0;
        let mut g2sum = 0.0;
        let ks: Vec<Int> = (0..samples).map(|_| rand_scalar(&mut rng)).collect();
        for k in &ks {
            let c = b2.decompose(k);
            assert!((&(&c[0] + &(&c[1] * &mu)) - k).modulo(&r).is_zero());
            let b = c.iter().map(|x| x.bits()).max().unwrap();
            g2max = g2max.max(b);
            g2sum += b as f64;
        }
        let mut emb_stats = vec![];
        for ss in [s.clone(), (-&s).modulo(&r)] {
            let lam = (&(&ba + &(&bb * &ss)) * &inv2).modulo(&r);
            let nal = &(&(&ba * &ba) - &(&Int::from(d) * &(&bb * &bb))).div_floor(&Int::from(4i64));
            assert!(
                (&(&(&lam * &lam) - &(&ba * &lam)) + nal)
                    .modulo(&r)
                    .is_zero(),
                "lambda char poly"
            );
            let b4 = Babai::new(lll4(&r, &mu, &lam), r.bits());
            let lm = (&lam * &mu).modulo(&r);
            let mut mx = 0;
            let mut sm = 0.0;
            for k in &ks {
                let c = b4.decompose(k);
                let sum = &(&(&c[0] + &(&c[1] * &mu)) + &(&c[2] * &lam)) + &(&c[3] * &lm);
                assert!((&sum - k).modulo(&r).is_zero());
                let b = c.iter().map(|x| x.bits()).max().unwrap();
                mx = mx.max(b);
                sm += b as f64;
            }
            emb_stats.push((mx, sm / samples as f64));
        }
        let loop_cost = cost(&best);
        let line = format!(
            "| {d} | {h} | {} | {} | {} | l = {bl}, k = {bk}, F steps {bkk}, N(alpha) = {bl}*2^{bkk} | {loop_cost:.1} | {g2max} / {:.1} | {} / {:.2} ; {} / {:.2} |",
            t.bits(),
            v.bits(),
            r.bits(),
            g2sum / samples as f64,
            emb_stats[0].0,
            emb_stats[0].1,
            emb_stats[1].0,
            emb_stats[1].1
        );
        eprintln!("{line}");
        lines.push(line);
        kept.push(Json::obj(vec![
            ("D", Json::Num(d)),
            ("h", Json::Num(h as i64)),
            ("t", Json::str(&t)),
            ("v", Json::str(&v)),
            ("r", Json::str(&r)),
            ("mu", Json::str(&mu)),
            (
                "best_loop",
                Json::obj(vec![
                    ("l", Json::Num(bl as i64)),
                    ("k", Json::Num(bk as i64)),
                    ("frobenius_steps", Json::Num(bkk as i64)),
                    ("A", Json::str(&ba)),
                    ("B", Json::str(&bb)),
                    ("model_cost_M", Json::Float(loop_cost)),
                ]),
            ),
            (
                "loops",
                Json::Arr(
                    loops
                        .iter()
                        .map(|lp| {
                            Json::obj(vec![
                                ("l", Json::Num(lp.0 as i64)),
                                ("k", Json::Num(lp.1 as i64)),
                                ("frobenius_steps", Json::Num(lp.2 as i64)),
                                ("A", Json::str(&lp.3)),
                                ("B", Json::str(&lp.4)),
                                ("model_cost_M", Json::Float(cost(lp))),
                            ])
                        })
                        .collect(),
                ),
            ),
            (
                "gls2",
                Json::obj(vec![
                    ("max_bits", Json::Num(g2max as i64)),
                    ("mean_bits", Json::Float(g2sum / samples as f64)),
                ]),
            ),
            (
                "glvgls4_s",
                Json::obj(vec![
                    ("max_bits", Json::Num(emb_stats[0].0 as i64)),
                    ("mean_bits", Json::Float(emb_stats[0].1)),
                ]),
            ),
            (
                "glvgls4_minus_s",
                Json::obj(vec![
                    ("max_bits", Json::Num(emb_stats[1].0 as i64)),
                    ("mean_bits", Json::Float(emb_stats[1].1)),
                ]),
            ),
        ]));
    }
    let secs = t0.elapsed().as_secs_f64();
    println!("# Binary-field GLS with CM loop endomorphisms: scan for m = {m}\n");
    println!(
        "Scanned {scanned} fundamental D = 1 mod 8 with |D| <= {}; {ord_hits} have [p_2] of order {m}; {prime_r_hits} of those have prime r = ((2^{m} - 1)^2 + t^2)/2; {} of those have a Frobenius-assisted loop with l <= {lmax}. {secs:.1} s.\n",
        -d,
        kept.len()
    );
    println!("First discriminants with [p_2] of order {m} (D, h): {first_hits:?}\n");
    println!("| D | h | t bits | v bits | r bits | cheapest Frobenius-assisted loop | model cost M | GLS-2 max/mean bits | GLV-GLS-4 max/mean bits (s ; -s) |");
    println!("|--:|--:|--:|--:|--:|:--|--:|--:|:--|");
    for l in &lines {
        println!("{l}");
    }
    if let Some(path) = json_out {
        let doc = Json::obj(vec![
            ("schema", Json::str("binary_gls_scan/v1")),
            ("command", Json::str(args.join(" "))),
            ("m", Json::Num(m as i64)),
            ("scanned", Json::Num(scanned as i64)),
            ("order_m_hits", Json::Num(ord_hits as i64)),
            ("prime_r_hits", Json::Num(prime_r_hits as i64)),
            ("d_reached", Json::Num(d)),
            (
                "first_hits",
                Json::Arr(
                    first_hits
                        .iter()
                        .map(|(d, h)| Json::Arr(vec![Json::Num(*d), Json::Num(*h as i64)]))
                        .collect(),
                ),
            ),
            ("seconds", Json::Float(secs)),
            ("kept", Json::Arr(kept)),
        ]);
        std::fs::write(&path, doc.dump()).expect("write");
        eprintln!("wrote {path}");
    }
}

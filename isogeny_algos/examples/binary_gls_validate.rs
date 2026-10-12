//! Point-level validation of Frobenius-assisted loop endomorphisms on binary GLS curves, at a
//! reduced size (default m = 31: E/F_{2^31}, GLS twist over F_{2^62}, r about 61 bits).
//!
//! The mechanism (see `binary_gls_scan.rs`): for `D = 1 mod 8` with `[p_2]` of order `m`, an odd
//! split prime `l` with `[l][p_2]^k = 1` gives the loop `l` then `F^kk` (`kk = min(k, m - k)`
//! Frobenius steps `x -> x^2`), an endomorphism `alpha` of degree `l 2^kk`.  This example builds a
//! real curve and checks, on points of `E'(F_{q^2})[r]`:
//!
//! * `x(psi(P)) = x(P)^q` equals `x([mu] P)`, `mu = (q - 1)/t mod r`;
//! * `x(alpha'(P))` (the odd isogeny's x-map, then `kk` squarings) equals `x([lambda] P)` for
//!   `lambda = (A + B s)/2`, `s = +-(2 mu - t)/v`, for one embedding consistently.
//!
//! In characteristic 2 the twisting isomorphism `(x, y) -> (x, y + s x)` and the closing
//! isomorphism between models with equal `j` leave `x` unchanged, so x-coordinates suffice; they
//! identify `lambda` up to sign, which is all a decomposition needs.
//!
//! Curve: `j` a root of `H_D mod 2` (complex-analytic `H_D` over Z from `cm.rs`, reduced mod 2,
//! roots in `GF(2^{2m})`), `E: y^2 + xy = x^3 + 1/j`, twist `a_2 = delta` with `Tr(delta) = 1`.
//! Field: `GF2n` with `n = 2m <= 63`.
//!
//! ```bash
//! cargo run --release --example binary_gls_validate -- [--m 31] [--count 3] [--hmax 160] \
//!     [--lmax 41] [--points 24] [--seed 20261012] [--json out.json]
//! ```
use isogeny_algos::binary::{kernel_polys, kohel_codomain, x_map, BinCurve};
use isogeny_algos::cm::{hilbert_class_poly, reduced_forms};
use isogeny_algos::curve::Pt;
use isogeny_algos::field::{Field, Rng};
use isogeny_algos::gf2n::GF2n;
use isogeny_algos::gls_loop::{cornacchia, is_prime_u64, prime_form, sqrt_mod_2k, squarefree, v2};
use isogeny_algos::int::Int;
use isogeny_algos::json::Json;
use isogeny_algos::path::relation::Form;
use isogeny_algos::poly;

fn arg<T: std::str::FromStr>(args: &[String], name: &str, default: T) -> T {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .and_then(|v| v.parse().ok())
        .unwrap_or(default)
}

fn frob_pow(f: &GF2n, mut x: u64, k: usize) -> u64 {
    for _ in 0..k {
        x = f.sq(x);
    }
    x
}

fn eval_ratio(f: &GF2n, num: &[u64], den: &[u64], x: u64) -> Option<u64> {
    let d = poly::eval(f, &den.to_vec(), x);
    if d == 0 {
        return None;
    }
    Some(f.mul(poly::eval(f, &num.to_vec(), x), f.inv(d)))
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let m: usize = arg(&args, "--m", 31);
    let count: usize = arg(&args, "--count", 3);
    let hmax: usize = arg(&args, "--hmax", 160);
    let lmax: u64 = arg(&args, "--lmax", 41);
    let npts: usize = arg(&args, "--points", 24);
    let seed: u64 = arg(&args, "--seed", 20_261_012);
    let json_out: Option<String> = args
        .iter()
        .position(|a| a == "--json")
        .and_then(|i| args.get(i + 1))
        .cloned();
    assert!(2 * m <= 63, "GF2n holds n <= 63");
    let f = GF2n::new((2 * m) as u32);
    let q = Int::one().shl(m);
    let qm1 = &q - &Int::one();

    let mut rows = vec![];
    let mut all_pass = true;
    println!("# Binary GLS: point-level validation of Frobenius-assisted loops (m = {m}, field GF(2^{}))\n", 2 * m);
    println!("| D | h | r bits | loop (l, F steps) | alpha | H_D precision / worst rounding | crater roots | closing walks | psi acts as [mu] | twisted loop acts as | verdict |");
    println!("|--:|--:|--:|:--|:--|:--|--:|--:|:--|:--|:--|");
    let mut d = -15i64;
    while rows.len() < count && -d < 10_000_000 {
        d -= 8;
        if !squarefree(-d) {
            continue;
        }
        let Some(f2) = prime_form(d, 2) else { continue };
        let id = Form::identity(d as i128);
        if f2 == id || f2.pow(m as i128) != id {
            continue;
        }
        let h = reduced_forms(d).len();
        if h > hmax {
            continue;
        }
        let Some((t, v)) = cornacchia(d, 1, m)
            .into_iter()
            .find(|(a, b)| !a.is_even() && !b.is_zero())
        else {
            continue;
        };
        let nt = &(&qm1 * &qm1) + &(&t * &t);
        let r = nt.shr(1);
        if r.shl(1) != nt || !r.is_probable_prime() {
            continue;
        }
        let s2 = sqrt_mod_2k(&Int::from(d), m + 8)
            .into_iter()
            .find(|s| v2(&(&t + &(&v * s))) > m)
            .expect("orientation");
        // cheapest Frobenius-assisted loop
        let mut best: Option<(u64, usize, Int, Int)> = None;
        let mut l = 3u64;
        while l <= lmax && best.is_none() {
            if is_prime_u64(l) {
                if let Some(fl) = prime_form(d, l as i64) {
                    let mut g = fl;
                    for k in 0..m {
                        if g == id {
                            let kk = k.min(m - k);
                            let found = cornacchia(d, l, kk)
                                .into_iter()
                                .flat_map(|(a, b)| [(a.clone(), b.clone()), (a, -b)])
                                .find(|(a, b)| v2(&(a + &(b * &s2))) > kk);
                            if let Some((a, b)) = found {
                                best = Some((l, kk, a, b));
                            }
                            break;
                        }
                        g = g.compose(&f2).reduce();
                    }
                }
            }
            l += 2;
        }
        let Some((l, kk, a_al, b_al)) = best else {
            continue;
        };

        // ---- build the curve from H_D mod 2
        let hp = match hilbert_class_poly(d) {
            Ok(hp) => hp,
            Err(e) => {
                println!("| {d} | {h} | | | | H_D failed: {e} | | | | | FAIL |");
                all_pass = false;
                continue;
            }
        };
        let hpoly: Vec<u64> = hp.coeffs.iter().map(|c| c.mod_u64(2)).collect();
        let mut rng = Rng::new(seed ^ (-d) as u64);
        let crater = poly::roots(&f, &hpoly, &mut rng);
        let j0 = crater[0];
        let e0 = BinCurve::from_j(&f, j0);
        // walks: horizontal l-isogeny, then F^kk back to j0
        let mut walks: Vec<(Vec<u64>, Vec<u64>)> = vec![];
        for kp in kernel_polys(&f, &e0, l, &mut rng) {
            let cod = kohel_codomain(&f, &e0, &kp);
            let j1 = cod.j(&f);
            if !crater.contains(&j1) {
                continue;
            }
            if frob_pow(&f, j1, kk) == j0 {
                let (num, den) = x_map(&f, &kp);
                walks.push((num, den));
            }
        }
        // ---- GLS twist and eigenvalues
        let delta = loop {
            let c = f.random(&mut rng);
            if f.trace(c) == 1 {
                break c;
            }
        };
        let et = BinCurve::new(e0.a2 ^ delta, e0.a4, e0.a6);
        let mu = (&qm1 * &t.inv_mod(&r).unwrap()).modulo(&r);
        let s = (&(&(&Int::from(2i64) * &mu) - &t) * &v.inv_mod(&r).unwrap()).modulo(&r);
        assert!((&(&s * &s) - &Int::from(d)).modulo(&r).is_zero());
        let inv2 = Int::from(2i64).inv_mod(&r).unwrap();
        let lam = |ss: &Int| (&(&a_al + &(&b_al * ss)) * &inv2).modulo(&r);
        let lam_s = lam(&s);
        let lam_ms = lam(&(-&s).modulo(&r));
        let to_u128 = |x: &Int| x.to_i128().unwrap() as u128;
        let (mu_u, ls_u, lms_u, r_u) =
            (to_u128(&mu), to_u128(&lam_s), to_u128(&lam_ms), to_u128(&r));
        let mut psi_ok = 0usize;
        let mut walk_hits: Vec<Vec<&str>> = vec![vec![]; walks.len()];
        let mut tested = 0usize;
        while tested < npts {
            let p0 = et.random_point(&f, &mut rng);
            let pt = et.mul(&f, &p0, 2);
            let Pt::Aff(x, _) = pt else { continue };
            if et.mul(&f, &pt, r_u) != Pt::Inf {
                continue;
            }
            tested += 1;
            let xq = frob_pow(&f, x, m);
            if let Pt::Aff(xm, _) = et.mul(&f, &pt, mu_u) {
                if xm == xq {
                    psi_ok += 1;
                }
            }
            let x_s = if let Pt::Aff(xx, _) = et.mul(&f, &pt, ls_u) {
                Some(xx)
            } else {
                None
            };
            let x_ms = if let Pt::Aff(xx, _) = et.mul(&f, &pt, lms_u) {
                Some(xx)
            } else {
                None
            };
            for (wi, (num, den)) in walks.iter().enumerate() {
                let img = eval_ratio(&f, num, den, x).map(|y| frob_pow(&f, y, kk));
                let tag = if img.is_some() && img == x_s {
                    "lambda(s)"
                } else if img.is_some() && img == x_ms {
                    "lambda(-s)"
                } else {
                    "none"
                };
                walk_hits[wi].push(tag);
            }
        }
        let summary: Vec<String> = walk_hits
            .iter()
            .map(|v| {
                let first = v.first().copied().unwrap_or("none");
                if first != "none" && v.iter().all(|x| *x == first) {
                    format!("+-{first} ({}/{})", v.len(), v.len())
                } else {
                    format!("INCONSISTENT {v:?}")
                }
            })
            .collect();
        let pass = !walks.is_empty()
            && psi_ok == tested
            && summary.iter().all(|s| !s.starts_with("INCONSISTENT"));
        all_pass &= pass;
        println!(
            "| {d} | {h} | {} | ({l}, {kk}) | ({a_al} + {b_al} sqrt({d}))/2, N = {l}*2^{kk} | {} / {:.0} | {} | {} | {psi_ok}/{tested} | {} | {} |",
            r.bits(),
            hp.precision,
            hp.worst_error_log2,
            crater.len(),
            walks.len(),
            summary.join("; "),
            if pass { "PASS" } else { "FAIL" }
        );
        rows.push(Json::obj(vec![
            ("D", Json::Num(d)),
            ("h", Json::Num(h as i64)),
            ("t", Json::str(&t)),
            ("v", Json::str(&v)),
            ("r", Json::str(&r)),
            (
                "loop",
                Json::obj(vec![
                    ("l", Json::Num(l as i64)),
                    ("frobenius_steps", Json::Num(kk as i64)),
                    ("A", Json::str(&a_al)),
                    ("B", Json::str(&b_al)),
                ]),
            ),
            (
                "hilbert",
                Json::obj(vec![
                    ("degree", Json::Num((hp.coeffs.len() - 1) as i64)),
                    ("precision_bits", Json::Num(hp.precision as i64)),
                    ("worst_rounding_log2", Json::Float(hp.worst_error_log2)),
                ]),
            ),
            ("crater_roots", Json::Num(crater.len() as i64)),
            ("closing_walks", Json::Num(walks.len() as i64)),
            ("points", Json::Num(tested as i64)),
            ("psi_matches_mu", Json::Num(psi_ok as i64)),
            ("walks", Json::Arr(summary.iter().map(Json::str).collect())),
            ("pass", Json::Bool(pass)),
        ]));
    }
    println!(
        "\nverdict: {}",
        if all_pass && !rows.is_empty() {
            "ALL PASS"
        } else {
            "FAIL"
        }
    );
    if let Some(path) = json_out {
        let doc = Json::obj(vec![
            ("schema", Json::str("binary_gls_validate/v1")),
            ("command", Json::str(args.join(" "))),
            ("m", Json::Num(m as i64)),
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

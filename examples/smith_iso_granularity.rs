//! Isomorphism-class granularity of the Smith class-level census
//! (`research/smith_iso_granularity_20261015/PROTOCOL.md`).
//!
//! ```text
//! cargo run --release --example smith_iso_granularity -- --p 7 --out DIR
//! ```
//!
//! Analyses every census model (monic octic, zero `x^7` coefficient, and its
//! twist by the least non-residue) exactly as `smith_prop24 census` does,
//! partitions the models into `F_p`-isomorphism classes of curves, and writes
//! one row per isomorphism class to `iso_p<p>.csv.gz`
//! (`rep_code,rep_twist,a1,a2,a3,stratum,full,size`) plus the controls of
//! the protocol to `controls_p<p>.json`. Orbits are enumerated by sending
//! each non-Weierstrass point of `P^1(F_p)` to infinity, normalising the
//! leading coefficient by a square, depressing, and scaling `x -> u x`.

use crypto_lib::cryptanalysis::smith_kernel_coverage::Fpk;
use crypto_lib::cryptanalysis::smith_prop24 as sp;
use flate2::{write::GzEncoder, Compression};
use rayon::prelude::*;
use serde_json::json;
use std::collections::HashMap;
use std::io::Write;
use std::path::PathBuf;
use std::{fs, time::Instant};

fn arg(args: &[String], k: &str) -> Option<String> {
    args.iter()
        .position(|a| a == k)
        .and_then(|i| args.get(i + 1).cloned())
}

#[derive(Clone, Copy, Default)]
struct Rec {
    ok: bool,
    key: [i16; 3],
    stratum: u8,
    full: [bool; 2],
}

fn twisted_key(k: [i16; 3], tw: bool) -> [i16; 3] {
    if tw {
        [-k[0], k[1], -k[2]]
    } else {
        k
    }
}

/// `f(x + a)` for a degree-8 coefficient vector (low degree first).
fn shift(f: &[u64; 9], a: u64, p: u64) -> [u64; 9] {
    let mut g = *f;
    for i in 0..9 {
        for j in (i..8).rev() {
            g[j] = (g[j] + a * g[j + 1]) % p;
        }
    }
    g
}

/// Orbit of the model `(code, tw)` as a list of census models.
fn orbit(code: u64, tw: bool, p: u64) -> Vec<(u64, bool)> {
    let mut f = [0u64; 9];
    let mut v = code;
    for c in f.iter_mut().take(7) {
        *c = v % p;
        v /= p;
    }
    f[8] = 1;
    let inv8 = sp::inv_mod(8 % p, p);
    let mut out = Vec::new();
    // P = infinity, then P = a for each a with f(a) != 0.
    let mut images = vec![f];
    for a in 0..p {
        let g = shift(&f, a, p);
        if g[0] != 0 {
            let mut h = g;
            h.reverse();
            images.push(h);
        }
    }
    for h in images {
        let lc = h[8];
        let tw2 = tw ^ (sp::chi(lc, p) == -1);
        let il = sp::inv_mod(lc, p);
        let mut m = h.map(|c| c * il % p);
        let a = (p - m[7]) * inv8 % p;
        m = shift(&m, a, p);
        debug_assert_eq!(m[7], 0);
        for u in 1..p {
            // coefficient i scaled by u^(i-8)
            let ui = sp::inv_mod(u, p);
            let mut s = 1u64; // u^(i-8) for i = 8 downwards: ui^(8-i)
            let mut code2 = 0u64;
            let mut coeffs = [0u64; 7];
            for i in (0..8).rev() {
                s = s * ui % p;
                if i < 7 {
                    coeffs[i] = m[i] * s % p;
                }
            }
            for i in (0..7).rev() {
                code2 = code2 * p + coeffs[i];
            }
            out.push((code2, tw2));
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let p: u64 = arg(&args, "--p").expect("--p").parse().unwrap();
    let out = PathBuf::from(arg(&args, "--out").expect("--out DIR"));
    fs::create_dir_all(&out).unwrap();
    let t0 = Instant::now();
    let pu = p as u32;
    let fields = [Fpk::new(pu, 1), Fpk::new(pu, 2), Fpk::new(pu, 3)];
    let total = p.pow(7);
    let recs: Vec<Rec> = (0..total)
        .into_par_iter()
        .map(|code| {
            let mut c = [0u32; 8];
            let mut v = code;
            for cd in c.iter_mut().take(7) {
                *cd = (v % p) as u32;
                v /= p;
            }
            let mut f: sp::Poly = c.iter().map(|&x| x as u64).collect();
            f.push(1);
            let mut rng = sp::SplitMix64(code);
            let Some(o) = sp::analyze_curve(&f, p, &mut rng) else {
                return Rec::default();
            };
            let s: Vec<i64> = fields.iter().map(|fk| fk.char_sum(&c)).collect();
            let pk = [-s[0] - 1, -s[1] - 1, -s[2] - 1];
            let e1 = pk[0];
            let e2 = (e1 * pk[0] - pk[1]) / 2;
            let e3 = (e2 * pk[0] - e1 * pk[1] + pk[2]) / 3;
            let ns = o.n_tractable() as u32;
            Rec {
                ok: true,
                key: [(-e1) as i16, e2 as i16, (-e3) as i16],
                stratum: sp::STRATA.iter().position(|&x| x == ns).unwrap() as u8,
                full: [o.full_success(false), o.full_success(true)],
            }
        })
        .collect();
    let t_an = t0.elapsed().as_secs_f64();
    eprintln!("analysed {total} octics in {t_an:.1}s");

    let mut visited = vec![false; 2 * total as usize];
    let idx = |c: u64, t: bool| (2 * c + t as u64) as usize;
    let mut rows = String::from("rep_code,rep_twist,a1,a2,a3,stratum,full,size\n");
    let (mut n_iso, mut n_models, mut c2_fail) = (0u64, 0u64, 0u64);
    let mut c2_nonsqfree = 0u64;
    // per L-key: models, tractable, full (control C1)
    let mut agg: HashMap<[i16; 3], [u64; 3]> = HashMap::new();
    for code in 0..total {
        if !recs[code as usize].ok {
            continue;
        }
        for tw in [false, true] {
            if visited[idx(code, tw)] {
                continue;
            }
            let r = recs[code as usize];
            let key = twisted_key(r.key, tw);
            let full = r.full[tw as usize];
            let orb = orbit(code, tw, p);
            for &(c2, t2) in &orb {
                let r2 = recs[c2 as usize];
                if !r2.ok {
                    c2_nonsqfree += 1;
                    continue;
                }
                debug_assert!(!visited[idx(c2, t2)]);
                visited[idx(c2, t2)] = true;
                if twisted_key(r2.key, t2) != key
                    || r2.stratum != r.stratum
                    || r2.full[t2 as usize] != full
                {
                    c2_fail += 1;
                }
            }
            let size = orb.len() as u64;
            n_iso += 1;
            n_models += size;
            let e = agg.entry(key).or_default();
            e[0] += size;
            e[1] += size * (r.stratum != 0) as u64;
            e[2] += size * full as u64;
            rows.push_str(&format!(
                "{code},{},{},{},{},{},{},{size}\n",
                tw as u8, key[0], key[1], key[2], r.stratum, full as u8
            ));
        }
    }
    let unvisited = (0..total)
        .filter(|&c| recs[c as usize].ok)
        .map(|c| (!visited[idx(c, false)]) as u64 + (!visited[idx(c, true)]) as u64)
        .sum::<u64>();

    // Control C1 against the 2026-10-14 class rows.
    let census = format!(
        "research/smith_prop24_native_20261014/data/census_p{p}_classes.csv"
    );
    let text = fs::read_to_string(&census).expect("census class rows");
    let (mut c1_rows, mut c1_mismatch) = (0u64, 0u64);
    for line in text.lines().skip(1) {
        let v: Vec<i64> = line.split(',').map(|x| x.parse().unwrap()).collect();
        let key = [v[0] as i16, v[1] as i16, v[2] as i16];
        let want = [v[3] as u64, v[4] as u64, v[6] as u64];
        c1_rows += 1;
        if agg.get(&key).copied().unwrap_or_default() != want {
            c1_mismatch += 1;
        }
    }
    let c1_pass = c1_mismatch == 0 && c1_rows == agg.len() as u64;
    let secs = t0.elapsed().as_secs_f64();

    let mut gz = GzEncoder::new(Vec::new(), Compression::best());
    gz.write_all(rows.as_bytes()).unwrap();
    fs::write(out.join(format!("iso_p{p}.csv.gz")), gz.finish().unwrap()).unwrap();
    let git = std::process::Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_default();
    let ctl = json!({
        "p": p,
        "git_head": git,
        "octics": total,
        "iso_classes": n_iso,
        "models": n_models,
        "l_classes": agg.len(),
        "c1_census_rows": c1_rows,
        "c1_mismatched_rows": c1_mismatch,
        "c1_pass": c1_pass,
        "c2_inconsistent_models": c2_fail,
        "c2_orbit_hits_nonsquarefree": c2_nonsqfree,
        "c2_pass": c2_fail == 0 && c2_nonsqfree == 0,
        "c3_unvisited_models": unvisited,
        "c3_pass": unvisited == 0,
        "seconds_analysis": t_an,
        "seconds_total": secs,
        "threads": rayon::current_num_threads(),
    });
    fs::write(
        out.join(format!("controls_p{p}.json")),
        serde_json::to_string_pretty(&ctl).unwrap() + "\n",
    )
    .unwrap();
    println!("{}", serde_json::to_string_pretty(&ctl).unwrap());
}

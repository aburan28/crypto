//! Isomorphism-class granularity of the Smith class-level census
//! (`research/smith_iso_granularity_20261015/PROTOCOL.md`).
//!
//! ```text
//! cargo run --release --example smith_iso_granularity -- --p 7 --out DIR
//! cargo run --release --example smith_iso_granularity -- --p 7 --out DIR --subgroups
//! cargo run --release --example smith_iso_granularity -- summarize --dir DIR 7 11
//! cargo run --release --example smith_iso_granularity -- summarize-serre --dir DIR 7 11
//! ```
//!
//! `--subgroups` also writes per-isomorphism-class Prop. 2 / Prop. 4 subgroup
//! counts (`iso_sub_p<p>.csv.gz`) for `summarize-serre`, which runs H1 and H2
//! of `research/smith_serre_sign_20261016/PROTOCOL.md`.
//!
//! `summarize` reads `iso_p<p>.csv.gz` and `controls_p<p>.json`, runs the
//! protocol's nulls (999 within-stratum permutations of the `full_success`
//! labels over isomorphism classes, seed 20261015; and the same stratified by
//! stratum and orbit size), and writes `summary_p<p>.json` and `verdict.json`.
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
use serde_json::Value;
use std::collections::{BTreeMap, HashMap};
use std::io::Read;
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
    /// Rational tractable subgroups passing Prop. 2 (twist-invariant).
    n_trig: u8,
    /// Of those, passing Prop. 4 on the curve and on its twist.
    n_iso: [u8; 2],
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
    if args.get(1).map(String::as_str) == Some("summarize") {
        summarize_cmd(&args);
        return;
    }
    if args.get(1).map(String::as_str) == Some("summarize-serre") {
        serre_cmd(&args);
        return;
    }
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
                n_trig: o.n_trig() as u8,
                n_iso: [o.n_iso(false) as u8, o.n_iso(true) as u8],
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
    // `--subgroups`: per-iso-class Prop. 2 / Prop. 4 subgroup counts, with
    // their orbit constancy (c4) and twist complementarity (c5).
    let want_sub = args.iter().any(|a| a == "--subgroups");
    let mut sub = String::from("rep_code,rep_twist,n_trig,n_iso\n");
    let (mut c4_fail, mut c5_fail) = (0u64, 0u64);
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
                if r2.n_trig != r.n_trig || r2.n_iso[t2 as usize] != r.n_iso[tw as usize] {
                    c4_fail += 1;
                }
            }
            c5_fail += (r.n_iso[0] + r.n_iso[1] != r.n_trig) as u64;
            if want_sub {
                sub.push_str(&format!(
                    "{code},{},{},{}\n",
                    tw as u8, r.n_trig, r.n_iso[tw as usize]
                ));
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
    let census = format!("research/smith_prop24_native_20261014/data/census_p{p}_classes.csv");
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
    if want_sub {
        let mut gz = GzEncoder::new(Vec::new(), Compression::best());
        gz.write_all(sub.as_bytes()).unwrap();
        fs::write(
            out.join(format!("iso_sub_p{p}.csv.gz")),
            gz.finish().unwrap(),
        )
        .unwrap();
    }
    let mut ctl = json!({
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
    if want_sub {
        ctl["c4_subgroup_counts_inconsistent"] = json!(c4_fail);
        ctl["c4_pass"] = json!(c4_fail == 0);
        ctl["c5_twist_complementarity_violations"] = json!(c5_fail);
        ctl["c5_pass"] = json!(c5_fail == 0);
    }
    fs::write(
        out.join(format!("controls_p{p}.json")),
        serde_json::to_string_pretty(&ctl).unwrap() + "\n",
    )
    .unwrap();
    println!("{}", serde_json::to_string_pretty(&ctl).unwrap());
}

const PERMS: usize = 999;
const SEED: u64 = 20261015;

struct Iso {
    lid: Vec<usize>,
    key: Vec<[i64; 3]>,
    stratum: Vec<u8>,
    full: Vec<bool>,
    size: Vec<u64>,
    n_l: usize,
}

fn load(dir: &std::path::Path, p: u64) -> Iso {
    let mut text = String::new();
    flate2::read::GzDecoder::new(fs::File::open(dir.join(format!("iso_p{p}.csv.gz"))).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    let mut ids: HashMap<[i64; 3], usize> = HashMap::new();
    let mut iso = Iso {
        lid: vec![],
        key: vec![],
        stratum: vec![],
        full: vec![],
        size: vec![],
        n_l: 0,
    };
    for line in text.lines().skip(1) {
        let v: Vec<i64> = line.split(',').map(|x| x.parse().unwrap()).collect();
        let k = [v[2], v[3], v[4]];
        let n = ids.len();
        iso.lid.push(*ids.entry(k).or_insert(n));
        iso.key.push(k);
        iso.stratum.push(v[5] as u8);
        iso.full.push(v[6] == 1);
        iso.size.push(v[7] as u64);
    }
    iso.n_l = ids.len();
    iso
}

/// Share of models lying in an `L`-class with a `full_success` member.
fn class_share(iso: &Iso, full: &[bool], total: u64) -> f64 {
    let mut hit = vec![false; iso.n_l];
    for (i, &f) in full.iter().enumerate() {
        if f {
            hit[iso.lid[i]] = true;
        }
    }
    let s: u64 = (0..full.len())
        .filter(|&i| hit[iso.lid[i]])
        .map(|i| iso.size[i])
        .sum();
    s as f64 / total as f64
}

/// Permute `full` within groups; quantiles of the class share.
fn null(iso: &Iso, group: &[u64], total: u64, rng: &mut sp::SplitMix64) -> Vec<f64> {
    let mut groups: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
    for (i, &g) in group.iter().enumerate() {
        groups.entry(g).or_default().push(i);
    }
    let mut out = Vec::with_capacity(PERMS);
    let mut f = iso.full.clone();
    for _ in 0..PERMS {
        for idx in groups.values() {
            let mut lab: Vec<bool> = idx.iter().map(|&i| iso.full[i]).collect();
            for a in (1..lab.len()).rev() {
                let b = rng.below(a as u64 + 1) as usize;
                lab.swap(a, b);
            }
            for (j, &i) in idx.iter().enumerate() {
                f[i] = lab[j];
            }
        }
        out.push(class_share(iso, &f, total));
    }
    out.sort_by(|a, b| a.partial_cmp(b).unwrap());
    out
}

fn quant(v: &[f64], x: f64) -> f64 {
    // linear interpolation between order statistics
    let h = x * (v.len() - 1) as f64;
    let (lo, hi) = (h.floor() as usize, h.ceil() as usize);
    v[lo] + (h - lo as f64) * (v[hi] - v[lo])
}

fn band(c: f64, ns: &[f64]) -> &'static str {
    if c < quant(ns, 0.025) {
        "residual deficit beyond granularity"
    } else if c > quant(ns, 0.975) {
        "excess over granularity"
    } else {
        "granularity explains the deficit"
    }
}

fn null_json(ns: &[f64], c: f64) -> Value {
    json!({
        "perms": PERMS, "seed": SEED,
        "median": quant(ns, 0.5), "lo": quant(ns, 0.025), "hi": quant(ns, 0.975),
        "frac_le_observed": ns.iter().filter(|&&x| x <= c).count() as f64 / ns.len() as f64,
    })
}

fn summarize(dir: &std::path::Path, p: u64) -> Value {
    let ctl: Value =
        serde_json::from_str(&fs::read_to_string(dir.join(format!("controls_p{p}.json"))).unwrap())
            .unwrap();
    let valid = ["c1_pass", "c2_pass", "c3_pass"]
        .iter()
        .all(|k| ctl[k].as_bool() == Some(true));
    let iso = load(dir, p);
    let n = iso.size.len();
    let total: u64 = iso.size.iter().sum();
    let mut tract_l = vec![false; iso.n_l];
    for i in 0..n {
        if iso.stratum[i] != 0 {
            tract_l[iso.lid[i]] = true;
        }
    }
    let share = |pred: &dyn Fn(usize) -> bool| -> f64 {
        (0..n)
            .filter(|&i| pred(i))
            .map(|i| iso.size[i])
            .sum::<u64>() as f64
            / total as f64
    };
    let u = share(&|i| tract_l[iso.lid[i]]);
    let d = share(&|i| iso.full[i]);
    let c = class_share(&iso, &iso.full, total);
    let mut rng = sp::SplitMix64(SEED);
    let g1: Vec<u64> = iso.stratum.iter().map(|&s| s as u64).collect();
    let prim = null(&iso, &g1, total, &mut rng);
    let g2: Vec<u64> = (0..n)
        .map(|i| iso.stratum[i] as u64 * 1_000_000 + iso.size[i])
        .collect();
    let sec = null(&iso, &g2, total, &mut rng);

    // Per-stratum iso-level success rate q_s.
    let mut qs = [[0u64; 2]; 9];
    for i in 0..n {
        qs[iso.stratum[i] as usize][0] += 1;
        qs[iso.stratum[i] as usize][1] += iso.full[i] as u64;
    }
    let q = |s: u8| qs[s as usize][1] as f64 / qs[s as usize][0] as f64;
    // By k = number of tractable iso classes in a tractable L-class.
    let mut members: Vec<Vec<usize>> = vec![vec![]; iso.n_l];
    for i in 0..n {
        if iso.stratum[i] != 0 {
            members[iso.lid[i]].push(i);
        }
    }
    let mut kb: BTreeMap<usize, (u64, u64, f64)> = BTreeMap::new();
    let mut ks = vec![];
    for m in members.iter().filter(|m| !m.is_empty()) {
        ks.push(m.len());
        let e = kb.entry(m.len().min(10)).or_default();
        e.0 += 1;
        e.1 += !m.iter().any(|&i| iso.full[i]) as u64;
        e.2 += m.iter().map(|&i| 1.0 - q(iso.stratum[i])).product::<f64>();
    }
    ks.sort_unstable();
    let by_k: BTreeMap<String, Value> = kb
        .iter()
        .map(|(k, e)| {
            let lab = if *k == 10 {
                "10+".into()
            } else {
                k.to_string()
            };
            (
                format!("{lab:>3}"),
                json!({"classes": e.0, "no_full": e.1, "expected_no_full": e.2}),
            )
        })
        .collect();
    // Post hoc probe: success rate by stratum and v2(#J(F_p)), #J = L(1).
    let mut probe: BTreeMap<String, [u64; 2]> = BTreeMap::new();
    let pi = p as i64;
    for i in 0..n {
        if iso.stratum[i] == 0 {
            continue;
        }
        let [a1, a2, a3] = iso.key[i];
        let nj = 1 + a1 + a2 + a3 + pi * a2 + pi * pi * a1 + pi * pi * pi;
        let v2 = (nj.trailing_zeros()).min(6);
        let e = probe
            .entry(format!(
                "s{}_v2_{}",
                sp::STRATA[iso.stratum[i] as usize],
                v2
            ))
            .or_default();
        e[0] += 1;
        e[1] += iso.full[i] as u64;
    }
    let s_rates: BTreeMap<String, Value> = (1..9u8)
        .filter(|&s| qs[s as usize][0] > 0)
        .map(|s| {
            (
                format!("#S(H)={}", sp::STRATA[s as usize]),
                json!({"iso_classes": qs[s as usize][0], "full": qs[s as usize][1], "rate": q(s)}),
            )
        })
        .collect();
    json!({
        "p": p,
        "controls": ctl,
        "instrument_valid": valid,
        "models": total,
        "iso_classes": n,
        "tractable_iso_classes": (0..n).filter(|&i| iso.stratum[i] != 0).count(),
        "full_iso_classes": iso.full.iter().filter(|&&f| f).count(),
        "l_classes": iso.n_l,
        "tractable_l_classes": ks.len(),
        "median_tractable_iso_classes_per_tractable_l_class": ks[ks.len() / 2],
        "U": u, "D": d, "C": c, "r": (c - d) / (u - d),
        "iso_level_success_by_stratum": s_rates,
        "null_primary": null_json(&prim, c),
        "null_secondary": null_json(&sec, c),
        "verdict": if valid { band(c, &prim) } else { "control failed - invalid" },
        "verdict_secondary": band(c, &sec),
        "tractable_l_classes_without_full": {
            "observed": kb.values().map(|e| e.1).sum::<u64>(),
            "expected_independent_iso": kb.values().map(|e| e.2).sum::<f64>(),
        },
        "by_k": by_k,
        "post_hoc_by_stratum_and_v2_of_J": probe,
    })
}

fn summarize_cmd(args: &[String]) {
    let dir = PathBuf::from(arg(args, "--dir").expect("--dir DIR"));
    let ps: Vec<u64> = args[4..].iter().map(|s| s.parse().unwrap()).collect();
    let mut verdicts = BTreeMap::new();
    for p in ps {
        let s = summarize(&dir, p);
        fs::write(
            dir.join(format!("summary_p{p}.json")),
            serde_json::to_string_pretty(&s).unwrap() + "\n",
        )
        .unwrap();
        verdicts.insert(p.to_string(), s["verdict"].clone());
        eprintln!(
            "p={p}: C={:.4} null [{:.4}, {:.4}] -> {}",
            s["C"], s["null_primary"]["lo"], s["null_primary"]["hi"], s["verdict"]
        );
    }
    let vs: std::collections::BTreeSet<String> = verdicts.values().map(|v| v.to_string()).collect();
    let overall = if vs.len() == 1 {
        verdicts.values().next().unwrap().clone()
    } else {
        json!("prime-dependent")
    };
    fs::write(
        dir.join("verdict.json"),
        serde_json::to_string_pretty(&json!({"per_prime": verdicts, "overall": overall})).unwrap()
            + "\n",
    )
    .unwrap();
}

const SEED_SERRE: u64 = 20261016;

/// H1 and H2 of `research/smith_serre_sign_20261016/PROTOCOL.md`.
fn serre(dir: &std::path::Path, p: u64) -> Value {
    let ctl: Value =
        serde_json::from_str(&fs::read_to_string(dir.join(format!("controls_p{p}.json"))).unwrap())
            .unwrap();
    let valid = ["c1_pass", "c2_pass", "c3_pass", "c4_pass"]
        .iter()
        .all(|k| ctl[k].as_bool() == Some(true));
    let iso = load(dir, p);
    let mut text = String::new();
    flate2::read::GzDecoder::new(fs::File::open(dir.join(format!("iso_sub_p{p}.csv.gz"))).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    let sub: Vec<[u64; 2]> = text
        .lines()
        .skip(1)
        .map(|l| {
            let v: Vec<u64> = l.split(',').map(|x| x.parse().unwrap()).collect();
            [v[2], v[3]]
        })
        .collect();
    let n = iso.size.len();
    assert_eq!(sub.len(), n, "iso and iso_sub rows align");
    let total: u64 = iso.size.iter().sum();
    let n1: Vec<u64> = iso
        .key
        .iter()
        .map(|k| (p as i64 + 1 + k[0]) as u64)
        .collect();
    let nb = 2 * p as usize + 3;

    // C5': curve-side + twist-side Prop. 4 passes <= Prop. 2 passes, per N1
    // pair (models weighted; the twist of a curve with N1 has 2p+2-N1).
    let (mut t, mut ip) = (vec![0u64; nb], vec![0u64; nb]);
    for i in 0..n {
        t[n1[i] as usize] += iso.size[i] * sub[i][0];
        ip[n1[i] as usize] += iso.size[i] * sub[i][1];
    }
    let mut c5p = true;
    let all_trig: u64 = t.iter().sum();
    let mut bins = vec![];
    let (mut mono, mut maxdev, mut prev) = (true, 0f64, -1f64);
    for b in 0..nb {
        let tw = 2 * p as usize + 2 - b.min(2 * p as usize + 2);
        if b > 2 * p as usize + 2 || t[b] == 0 {
            continue;
        }
        c5p &= ip[b] + ip[tw] <= t[b];
        let nondeg = ip[b] + ip[tw];
        let r = ip[b] as f64 / nondeg as f64;
        let model = b as f64 / (2 * p + 2) as f64;
        let big = t[b] as f64 >= 0.01 * all_trig as f64;
        if big {
            mono &= r >= prev;
            prev = r;
            maxdev = maxdev.max((r - model).abs());
        }
        bins.push(json!({
            "N1": b, "prop2_passes": t[b], "prop4_passes": ip[b],
            "degenerate_share": 1.0 - nondeg as f64 / t[b] as f64,
            "r": r, "model_N1_over_2p2": model, "in_test": big,
        }));
    }
    let h1 = if mono && maxdev <= 0.10 {
        "sign law holds"
    } else {
        "sign law fails"
    };

    // H2: nulls conditioned on (stratum, N1, any Prop. 2 pass).
    let mut tract_l = vec![false; iso.n_l];
    for i in 0..n {
        if iso.stratum[i] != 0 {
            tract_l[iso.lid[i]] = true;
        }
    }
    let u = (0..n)
        .filter(|&i| tract_l[iso.lid[i]])
        .map(|i| iso.size[i])
        .sum::<u64>() as f64
        / total as f64;
    let c = class_share(&iso, &iso.full, total);
    let mut rng = sp::SplitMix64(SEED_SERRE);
    let g0: Vec<u64> = iso.stratum.iter().map(|&s| s as u64).collect();
    let null0 = null(&iso, &g0, total, &mut rng);
    let cell = |lab: &[u64]| -> Vec<u64> {
        (0..n)
            .map(|i| (iso.stratum[i] as u64 * 1000 + lab[i]) * 2 + (sub[i][0] > 0) as u64)
            .collect()
    };
    let null_n1 = null(&iso, &cell(&n1), total, &mut rng);
    let n_vals = n1.iter().collect::<std::collections::BTreeSet<_>>().len() as u64;
    let mut lrng = sp::SplitMix64(SEED_SERRE ^ 0x5eed);
    let lab_l: Vec<u64> = (0..iso.n_l).map(|_| lrng.below(n_vals)).collect();
    let rand_lab: Vec<u64> = (0..n).map(|i| lab_l[iso.lid[i]]).collect();
    let null_rand = null(&iso, &cell(&rand_lab), total, &mut rng);
    let closed = |ns: &[f64]| (u - quant(ns, 0.5)) / (u - c);
    let (cl0, cl1, clr) = (closed(&null0), closed(&null_n1), closed(&null_rand));
    let gain_n1 = cl1 - cl0;
    let gain_rand = clr - cl0;
    let h2 = if !valid {
        "control failed - invalid"
    } else if gain_rand >= 0.75 * gain_n1 {
        "inconclusive: coarse conditioning alone"
    } else if c < quant(&null_n1, 0.025) {
        "partial: residual clustering remains"
    } else if c > quant(&null_n1, 0.975) {
        "excess over the N1 null"
    } else {
        "N1 and Prop. 2 explain the clustering"
    };
    let nj = |ns: &[f64]| {
        let mut v = null_json(ns, c);
        v["seed"] = json!(SEED_SERRE);
        v
    };
    json!({
        "p": p, "controls": ctl, "instrument_valid": valid,
        "c5prime_pass": c5p,
        "h1": {"bins": bins, "monotone": mono, "max_abs_dev": maxdev, "verdict": h1},
        "h2": {
            "U": u, "C": c,
            "null_stratum": nj(&null0), "null_stratum_n1_trig": nj(&null_n1),
            "null_stratum_random_label_trig": nj(&null_rand),
            "random_label_values": n_vals,
            "gap_closed_stratum": cl0, "gap_closed_n1": cl1, "gap_closed_random": clr,
            "verdict": h2,
        },
    })
}

fn serre_cmd(args: &[String]) {
    let dir = PathBuf::from(arg(args, "--dir").expect("--dir DIR"));
    let mut verdict = BTreeMap::new();
    for p in args[4..].iter().map(|s| s.parse::<u64>().unwrap()) {
        let s = serre(&dir, p);
        fs::write(
            dir.join(format!("serre_p{p}.json")),
            serde_json::to_string_pretty(&s).unwrap() + "\n",
        )
        .unwrap();
        eprintln!(
            "p={p}: H1 {} (max dev {:.3}); H2 {} (closed {:.2} / {:.2}, random {:.2})",
            s["h1"]["verdict"],
            s["h1"]["max_abs_dev"].as_f64().unwrap(),
            s["h2"]["verdict"],
            s["h2"]["gap_closed_stratum"].as_f64().unwrap(),
            s["h2"]["gap_closed_n1"].as_f64().unwrap(),
            s["h2"]["gap_closed_random"].as_f64().unwrap()
        );
        verdict.insert(
            p.to_string(),
            json!({"h1": s["h1"]["verdict"], "h2": s["h2"]["verdict"]}),
        );
    }
    fs::write(
        dir.join("serre_verdict.json"),
        serde_json::to_string_pretty(&json!({
            "per_prime": verdict,
            "confirmatory_prime": 11,
            "note": "p = 7 was used to find the law (exploratory)",
        }))
        .unwrap()
            + "\n",
    )
    .unwrap();
}

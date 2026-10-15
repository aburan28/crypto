//! `pkm-aux-search` — the wide auxiliary-curve search of the PKM isogeny
//! tower over the NIST P-192 prime, by 2-Sylow halving
//! (`crypto_lib::cryptanalysis::two_sylow_search`).
//!
//! ```text
//! pkm-aux-search --prime nist-p192 --t T --family z2xz2|z2xz8
//!     --seed-lo LO --seed-hi HI [--arm fast|fast-chi|baseline]
//!     [--out ARTIFACT.json] [--stop first_success|none]
//!     [--summary SUMMARY.json] [--record-min D]
//!     [--checkpoint FILE [--resume]]
//! ```
//!
//! Draw `s` for `s` in `[LO, HI)` is one curve of the family, derived from
//! `s` alone. For each draw the search computes `b` (the 2-Sylow is
//! `Z/2 × Z/2^b`). The first draw with `b ≥ T` gives a point `Q` of order
//! exactly `2^T`, certified over `BigUint` and written to `--out` in the
//! format of `backlog_check_aux2t` (`{p, a, b, Qx, Qy, t}` plus the draw's
//! provenance). With `--stop none` the whole range runs (a measurement).
//! The summary JSON holds the depth histogram, every draw with
//! `b ≥ --record-min` as `[seed, b]`, the operation counts and the time.
//! Exit 0 found (or `--stop none` finished), 3 range exhausted without a
//! curve, 2 bad usage.
//!
//! `--family supersingular-j1728` skips the search: over the P-192 prime
//! (`p ≡ 3 mod 4`) `y² = x³ + x` is supersingular with `#E = p + 1` and a
//! single rational 2-torsion point, so `E(F_p)` is cyclic of order `p + 1`
//! and `v₂(p + 1) = 64`. For draw `s` it takes `P` with abscissa from `s`,
//! `Q = [(p + 1)/2^t]P`, and accepts when `Q` has order exactly `2^t`
//! (probability 1/2 per draw), all over `BigUint`.

use std::process::ExitCode;
use std::time::Instant;

use crypto_lib::cryptanalysis::fp_backend::{FpP192, PrimeField};
use crypto_lib::cryptanalysis::two_sylow_search::{
    big_double, certify, p192, Arm, Family, Halver, SplitMix,
};
use num_bigint::BigUint;
use num_traits::{One, Zero};
use serde_json::{json, Value};

struct Args {
    t: u32,
    family: Family,
    arm: Arm,
    lo: u64,
    hi: u64,
    out: Option<String>,
    stop_first: bool,
    summary: Option<String>,
    record_min: u32,
    checkpoint: Option<String>,
    resume: bool,
}

fn parse() -> Result<Args, String> {
    parse_from(std::env::args().collect())
}

fn parse_from(all: Vec<String>) -> Result<Args, String> {
    let v: Vec<String> = all.into_iter().skip(1).collect();
    let mut it = v.iter();
    let mut a = Args {
        t: 0,
        family: Family::Z2xZ8,
        arm: Arm::Fast,
        lo: 0,
        hi: 0,
        out: None,
        stop_first: true,
        summary: None,
        record_min: 16,
        checkpoint: None,
        resume: false,
    };
    let mut prime_ok = false;
    let mut have = (false, false, false);
    while let Some(k) = it.next() {
        let mut val = || it.next().cloned().ok_or(format!("{k} needs a value"));
        match k.as_str() {
            "--prime" => {
                let p = val()?;
                if p != "nist-p192" {
                    return Err(format!("unsupported prime {p}"));
                }
                prime_ok = true;
            }
            "--t" => {
                a.t = val()?.parse().map_err(|e| format!("--t: {e}"))?;
                have.0 = true;
            }
            "--family" => {
                let f = val()?;
                a.family = Family::parse(&f).ok_or(format!("unknown family {f}"))?;
            }
            "--arm" => {
                a.arm = match val()?.as_str() {
                    "fast" => Arm::Fast,
                    "fast-chi" => Arm::FastChi,
                    "baseline" => Arm::Baseline,
                    o => return Err(format!("unknown arm {o}")),
                }
            }
            "--seed-lo" => {
                a.lo = val()?.parse().map_err(|e| format!("--seed-lo: {e}"))?;
                have.1 = true;
            }
            "--seed-hi" => {
                a.hi = val()?.parse().map_err(|e| format!("--seed-hi: {e}"))?;
                have.2 = true;
            }
            "--out" => a.out = Some(val()?),
            "--stop" => {
                a.stop_first = match val()?.as_str() {
                    "first_success" => true,
                    "none" => false,
                    o => return Err(format!("unknown stop rule {o}")),
                }
            }
            "--summary" => a.summary = Some(val()?),
            "--record-min" => {
                a.record_min = val()?.parse().map_err(|e| format!("--record-min: {e}"))?
            }
            "--checkpoint" => a.checkpoint = Some(val()?),
            "--resume" => a.resume = true,
            o => return Err(format!("unknown flag {o}")),
        }
    }
    if !prime_ok || !have.0 || !have.1 || !have.2 {
        return Err("--prime, --t, --seed-lo and --seed-hi are required".into());
    }
    if !(1..=64).contains(&a.t) || a.lo >= a.hi {
        return Err("need 1 ≤ t ≤ 64 and seed-lo < seed-hi".into());
    }
    Ok(a)
}

fn write_json(path: &str, v: &Value) -> std::io::Result<()> {
    let tmp = format!("{path}.tmp");
    std::fs::write(&tmp, serde_json::to_string_pretty(v).unwrap() + "\n")?;
    std::fs::rename(tmp, path)
}

/// Affine addition over `BigUint` on `y² = x³ + ax + b` (`None` = O).
fn big_add(
    p: &BigUint,
    a: &BigUint,
    p1: &Option<(BigUint, BigUint)>,
    p2: &Option<(BigUint, BigUint)>,
) -> Option<(BigUint, BigUint)> {
    let (Some((x1, y1)), Some((x2, y2))) = (p1, p2) else {
        return if p1.is_none() { p2.clone() } else { p1.clone() };
    };
    if x1 == x2 {
        if ((y1 + y2) % p).is_zero() {
            return None;
        }
        return big_double(p, a, p1);
    }
    let inv = |v: &BigUint| v.modpow(&(p - 2u32), p);
    let l = (((y2 + p - y1) % p) * inv(&((x2 + p - x1) % p))) % p;
    let x3 = ((&l * &l) % p + 2u32 * p - x1 - x2) % p;
    let y3 = ((&l * ((x1 + p - &x3) % p)) % p + p - y1) % p;
    Some((x3, y3))
}

fn supersingular(a: &Args) -> ExitCode {
    let p = p192();
    let one = BigUint::one();
    let (ca, cb) = (one.clone(), BigUint::zero());
    let cof = (&p + 1u32) >> a.t;
    let start = Instant::now();
    for seed in a.lo..a.hi {
        let mut rng = SplitMix(seed ^ 0x5355_5045_5253_494E);
        let w = [rng.next_u64(), rng.next_u64(), rng.next_u64()];
        let x = ((BigUint::from(w[2]) << 128u32) + (BigUint::from(w[1]) << 64u32) + w[0]) % &p;
        let rhs = (&x * &x * &x + &x) % &p;
        let y = rhs.modpow(&((&p + 1u32) >> 2u32), &p);
        if (&y * &y) % &p != rhs {
            continue;
        }
        // Q = [cof]P by double-and-add.
        let mut q: Option<(BigUint, BigUint)> = None;
        let base = Some((x, y));
        for i in (0..cof.bits()).rev() {
            q = big_double(&p, &ca, &q);
            if cof.bit(i) {
                q = big_add(&p, &ca, &q, &base);
            }
        }
        let Some((qx, qy)) = q else { continue };
        // Order exactly 2^t, checked by doublings.
        let mut r = Some((qx.clone(), qy.clone()));
        let mut ok = true;
        for _ in 1..a.t {
            r = big_double(&p, &ca, &r);
            if r.is_none() {
                ok = false;
                break;
            }
        }
        if !ok || big_double(&p, &ca, &r).is_some() {
            continue;
        }
        let art = json!({
            "p": p.to_string(), "a": ca.to_string(), "b": cb.to_string(),
            "Qx": qx.to_string(), "Qy": qy.to_string(), "t": a.t,
            "family": "supersingular-j1728", "seed": seed,
            "draws": seed - a.lo + 1, "cpu_seconds": start.elapsed().as_secs_f64(),
        });
        if let Some(out) = &a.out {
            if let Some(dir) = std::path::Path::new(out).parent() {
                let _ = std::fs::create_dir_all(dir);
            }
            write_json(out, &art).expect("write artifact");
        }
        println!(
            "{}",
            json!({"found": true, "t": a.t, "seed": seed, "draws": seed - a.lo + 1,
                              "elapsed_s": start.elapsed().as_secs_f64()})
        );
        return ExitCode::SUCCESS;
    }
    println!("{}", json!({"found": false, "t": a.t}));
    ExitCode::from(3)
}

fn main() -> ExitCode {
    if std::env::args().any(|a| a == "supersingular-j1728") {
        let v: Vec<String> = std::env::args()
            .map(|a| {
                if a == "supersingular-j1728" {
                    "z2xz2".into()
                } else {
                    a
                }
            })
            .collect();
        let saved: Vec<String> = v.clone();
        return match parse_from(saved) {
            Ok(a) => supersingular(&a),
            Err(e) => {
                eprintln!("pkm-aux-search: {e}");
                ExitCode::from(2)
            }
        };
    }
    let a = match parse() {
        Ok(a) => a,
        Err(e) => {
            eprintln!("pkm-aux-search: {e}");
            return ExitCode::from(2);
        }
    };
    let f = FpP192::new();
    let p = p192();
    let mut h = Halver::new(f.clone(), a.arm);
    let mut hist = vec![0u64; 66];
    let mut records: Vec<(u64, u32)> = Vec::new();
    let mut next = a.lo;
    let mut prior_s = 0.0f64;
    if a.resume {
        if let Some(cp) = &a.checkpoint {
            if let Ok(text) = std::fs::read_to_string(cp) {
                let v: Value = serde_json::from_str(&text).expect("checkpoint JSON");
                next = v["next_seed"].as_u64().expect("next_seed");
                prior_s = v["elapsed_s"].as_f64().unwrap_or(0.0);
                for (i, c) in v["hist"].as_array().expect("hist").iter().enumerate() {
                    hist[i] = c.as_u64().unwrap_or(0);
                }
                for r in v["records"].as_array().expect("records") {
                    records.push((r[0].as_u64().unwrap(), r[1].as_u64().unwrap() as u32));
                }
            }
        }
    }
    let start = Instant::now();
    let mut found: Option<Value> = None;
    let checkpoint_every = 1u64 << 22;
    let elapsed = |s: &Instant| prior_s + s.elapsed().as_secs_f64();
    let state = |next: u64, hist: &[u64], records: &[(u64, u32)], el: f64| {
        json!({"next_seed": next, "elapsed_s": el, "hist": hist,
               "records": records.iter().map(|&(s, d)| json!([s, d])).collect::<Vec<_>>()})
    };
    while next < a.hi {
        let seed = next;
        next += 1;
        let d = h.run_draw(a.family, seed, a.t);
        hist[d.depth as usize] += 1;
        if d.depth >= a.record_min {
            records.push((seed, d.depth));
        }
        if let Some(pt) = d.point.as_ref().filter(|_| found.is_none()) {
            let e = d.curve.e.map(|v| f.to_biguint(v));
            match certify(
                &p,
                [&e[0], &e[1], &e[2]],
                &f.to_biguint(pt.x),
                &f.to_biguint(pt.y),
                a.t,
            ) {
                Ok((ca, cb, qx, qy)) => {
                    let art = json!({
                        "p": p.to_string(), "a": ca.to_string(), "b": cb.to_string(),
                        "Qx": qx.to_string(), "Qy": qy.to_string(), "t": a.t,
                        "family": a.family.name(), "seed": seed, "depth_b": d.depth,
                        "e": e.iter().map(|v| v.to_string()).collect::<Vec<_>>(),
                        "draws": next - a.lo, "cpu_seconds": elapsed(&start),
                    });
                    if let Some(out) = &a.out {
                        if let Some(dir) = std::path::Path::new(out).parent() {
                            let _ = std::fs::create_dir_all(dir);
                        }
                        write_json(out, &art).expect("write artifact");
                    }
                    found = Some(art);
                    if a.stop_first {
                        break;
                    }
                }
                Err(e) => {
                    eprintln!("pkm-aux-search: seed {seed}: certification failed: {e}");
                    return ExitCode::from(1);
                }
            }
        }
        if let Some(cp) = &a.checkpoint {
            if (next - a.lo).is_multiple_of(checkpoint_every) {
                write_json(cp, &state(next, &hist, &records, elapsed(&start)))
                    .expect("write checkpoint");
            }
        }
    }
    let el = elapsed(&start);
    if let Some(cp) = &a.checkpoint {
        write_json(cp, &state(next, &hist, &records, el)).expect("write checkpoint");
    }
    let draws = next - a.lo;
    let c = &h.cost;
    let summary = json!({
        "schema": "pkm-aux-search.summary/v1",
        "prime": "nist-p192", "t": a.t, "family": a.family.name(),
        "arm": match a.arm { Arm::Fast => "fast", Arm::FastChi => "fast-chi", Arm::Baseline => "baseline" },
        "seed_lo": a.lo, "seed_hi": a.hi, "next_seed": next, "draws": draws,
        "elapsed_s": el, "us_per_draw": if draws > 0 { el * 1e6 / draws as f64 } else { 0.0 },
        "hist_depth": hist, "record_min": a.record_min,
        "records": records.iter().map(|&(s, d)| json!([s, d])).collect::<Vec<_>>(),
        "cost_this_process": {"draws": c.draws, "singular": c.singular, "legendre": c.legendre,
                 "exps": c.exps, "halvings": c.halvings},
        "found": found,
    });
    if let Some(sp) = &a.summary {
        write_json(sp, &summary).expect("write summary");
    }
    println!(
        "{}",
        json!({"draws": draws, "elapsed_s": el, "us_per_draw": summary["us_per_draw"],
               "found": found.is_some(), "t": a.t})
    );
    if found.is_some() || !a.stop_first {
        ExitCode::SUCCESS
    } else {
        ExitCode::from(3)
    }
}

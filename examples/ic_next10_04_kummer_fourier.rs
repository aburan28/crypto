//! IC next-10 item 4 (research/ic_next10_20261014/04-kummer-fourier):
//! exhaustive enumeration of every `m`-term decomposition over a factor base
//! of a prime curve, grouped by target abscissa, for the full Fourier
//! spectrum of the decomposition indicator over the exponent group `(Z/N)^m`.
//!
//! Base `kummer`: abscissas `ζ^e`, `ζ` a primitive `N = 2^t`-th root of
//! unity, label `e` (ecbench `prime-kummer:t=T`).  Base `random`: `N`
//! distinct uniform abscissas of `F_p^*`, label = draw index (the
//! random-factor-base null).  Both signs of each on-curve abscissa.
//!
//! Every multiset of `m` signed base points with no pair `P, −P` is summed
//! once; the x-only decomposition set of target `x(R)` is the set of
//! unordered label multisets reaching it.  Output (little-endian binary):
//! `u32 m, u32 N, u32 n_labels, n_labels × u16 on-curve labels, u32 T`,
//! then `T × u16 K` (orbit counts of the targets with at least two orbits)
//! and their `ΣK·m × u16` labels, target by target.
//! A JSON line with counts goes to stdout.
//!
//! ```text
//! cargo run --release --example ic_next10_04_kummer_fourier -- \
//!   --p 1073153 --a A --b B --t 11 --m 2 --base kummer --seed 1 --out FILE
//! ```

use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};

use crypto_lib::cryptanalysis::ic_boundary::PrimeCurve;
use crypto_lib::cryptanalysis::ic_framework::prime_tower::primitive_two_power_root;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

fn mulm(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}

fn inv(a: u64, p: u64) -> u64 {
    let (mut r0, mut r1) = (p as i128, a as i128);
    let (mut s0, mut s1) = (0i128, 1i128);
    while r1 != 0 {
        let q = r0 / r1;
        (r0, r1) = (r1, r0 - q * r1);
        (s0, s1) = (s1, s0 - q * s1);
    }
    debug_assert_eq!(r0, 1);
    s0.rem_euclid(p as i128) as u64
}

/// Affine point, `None` = infinity.
type Pt = Option<(u64, u64)>;

fn add(p: u64, a: u64, u: Pt, v: Pt) -> Pt {
    let (Some((x1, y1)), Some((x2, y2))) = (u, v) else {
        return u.or(v);
    };
    let lam = if x1 == x2 {
        if (y1 + y2) % p == 0 {
            return None;
        }
        let num = (3 * mulm(x1, x1, p) + a) % p;
        mulm(num, inv(2 * y1 % p, p), p)
    } else {
        mulm((y2 + p - y1) % p, inv((x2 + p - x1) % p, p), p)
    };
    let x3 = (mulm(lam, lam, p) + 2 * p - x1 - x2) % p;
    let y3 = (mulm(lam, (x1 + p - x3) % p, p) + p - y1) % p;
    Some((x3, y3))
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let get = |k: &str| {
        args.iter()
            .position(|a| a == k)
            .map(|i| args[i + 1].clone())
    };
    let num = |k: &str| get(k).map(|s| s.parse::<u64>().expect(k));
    let p = num("--p").expect("--p");
    let a = num("--a").expect("--a");
    let b = num("--b").expect("--b");
    let t = num("--t").expect("--t") as u32;
    let m = num("--m").expect("--m") as usize;
    let seed = num("--seed").unwrap_or(1);
    let base = get("--base").unwrap_or_else(|| "kummer".into());
    let out = get("--out").expect("--out");
    assert!(m == 2 || m == 3, "m must be 2 or 3");
    let curve = PrimeCurve { p, a, b };
    let n = 1u64 << t;
    // Abscissas with labels.
    let xs: Vec<u64> = match base.as_str() {
        "kummer" => {
            let zeta = primitive_two_power_root(p, t).expect("2^t | p - 1");
            let mut z = 1u64;
            (0..n)
                .map(|_| {
                    let r = z;
                    z = mulm(z, zeta, p);
                    r
                })
                .collect()
        }
        "random" => {
            let mut rng = StdRng::seed_from_u64(seed);
            let mut seen = std::collections::HashSet::new();
            let mut v = Vec::new();
            while (v.len() as u64) < n {
                let x = rng.gen_range(1..p);
                if seen.insert(x) {
                    v.push(x);
                }
            }
            v
        }
        _ => panic!("--base kummer|random"),
    };
    // Signed points: (point, label, sign-class id = label).
    let mut pts: Vec<((u64, u64), u16)> = Vec::new();
    let mut labels: Vec<u16> = Vec::new();
    for (e, &x) in xs.iter().enumerate() {
        if let Some(y) = curve.sqrt(curve.rhs(x)) {
            labels.push(e as u16);
            pts.push(((x, y), e as u16));
            if y != 0 {
                pts.push(((x, (p - y) % p), e as u16));
            }
        }
    }
    let ns = pts.len();
    let is_neg = |i: usize, j: usize| {
        let ((x1, y1), _) = pts[i];
        let ((x2, y2), _) = pts[j];
        x1 == x2 && y1 != y2
    };
    // (target abscissa, packed sorted labels).
    let mut recs: Vec<(u64, u64)> = Vec::new();
    let mut sums = 0u64;
    for i in 0..ns {
        for j in i..ns {
            if is_neg(i, j) {
                continue;
            }
            let s = add(p, a, Some(pts[i].0), Some(pts[j].0));
            if m == 2 {
                sums += 1;
                if let Some((x, _)) = s {
                    let (u, v) = (pts[i].1.min(pts[j].1), pts[i].1.max(pts[j].1));
                    recs.push((x, ((u as u64) << 16) | v as u64));
                }
                continue;
            }
            for k in j..ns {
                if is_neg(i, k) || is_neg(j, k) {
                    continue;
                }
                sums += 1;
                if let Some((x, _)) = add(p, a, s, Some(pts[k].0)) {
                    let mut l = [pts[i].1, pts[j].1, pts[k].1];
                    l.sort_unstable();
                    recs.push((
                        x,
                        ((l[0] as u64) << 32) | ((l[1] as u64) << 16) | l[2] as u64,
                    ));
                }
            }
        }
    }
    recs.sort_unstable();
    recs.dedup();
    let mut w = BufWriter::new(File::create(&out).expect("out"));
    let wu32 = |w: &mut BufWriter<File>, v: u32| w.write_all(&v.to_le_bytes()).unwrap();
    wu32(&mut w, m as u32);
    wu32(&mut w, n as u32);
    wu32(&mut w, labels.len() as u32);
    for &l in &labels {
        w.write_all(&l.to_le_bytes()).unwrap();
    }
    // Groups.
    let mut groups: Vec<(usize, usize)> = Vec::new();
    let mut targets = 0usize;
    let mut i = 0;
    while i < recs.len() {
        let mut j = i;
        while j < recs.len() && recs[j].0 == recs[i].0 {
            j += 1;
        }
        targets += 1;
        if j - i >= 2 {
            groups.push((i, j));
        }
        i = j;
    }
    wu32(&mut w, groups.len() as u32);
    let mut pairs = 0u64;
    let mut kmax = 0usize;
    for &(i, j) in &groups {
        let k = j - i;
        kmax = kmax.max(k);
        pairs += (k * (k - 1) / 2) as u64;
        w.write_all(&(k as u16).to_le_bytes()).unwrap();
    }
    for &(i, j) in &groups {
        for r in &recs[i..j] {
            for s in (0..m).rev() {
                w.write_all(&(((r.1 >> (16 * s)) & 0xffff) as u16).to_le_bytes())
                    .unwrap();
            }
        }
    }
    w.flush().unwrap();
    println!(
        "{{\"p\":{p},\"a\":{a},\"b\":{b},\"t\":{t},\"m\":{m},\"base\":\"{base}\",\"seed\":{seed},\
         \"onc_labels\":{},\"signed_points\":{ns},\"sums\":{sums},\"orbits\":{},\"targets\":{targets},\
         \"targets_ge2\":{},\"cross_pairs\":{pairs},\"kmax\":{kmax}}}",
        labels.len(),
        recs.len(),
        groups.len()
    );
}

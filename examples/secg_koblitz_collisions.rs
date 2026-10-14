//! Which pair sums of the small-x base collide on a Koblitz curve, and why.
//!
//!     cargo run --release --example secg_koblitz_collisions -- sect163k1 256
use crypto_lib::cryptanalysis::binary_small_x::{small_x_base, BinCurve, Pt};
use crypto_lib::cryptanalysis::secg_standards::secg;
use std::collections::HashMap;

fn run<const L: usize>(name: &str, size: usize) {
    let c = BinCurve::<L>::from_secg(secg(name).unwrap()).unwrap();
    // Source abscissae: the raw small-x points, then [h] of each.
    let mut raw_curve = c.clone();
    raw_curve.cofactor = 1;
    let mut raw: Vec<Pt<L>> = Vec::new();
    let mut base: Vec<(u64, Pt<L>)> = Vec::new();
    let mut seen = std::collections::HashSet::new();
    let mut want = size;
    let mut take = size;
    while base.len() < size {
        raw = small_x_base(&raw_curve, take).unwrap().0;
        base.clear();
        seen.clear();
        for p in &raw {
            if let Some(q) = c.mul_small(p, c.cofactor) {
                if seen.insert(q.x) {
                    base.push((p.x[0], q));
                }
            }
        }
        take += want / 2 + 1;
        want = size;
    }
    base.truncate(size);
    let real = small_x_base(&c, size).unwrap().0;
    assert!(base.iter().map(|b| b.1).eq(real.iter().copied()));

    let mut by_key: HashMap<[u64; L], Vec<String>> = HashMap::new();
    for i in 0..base.len() {
        let (xi, a) = &base[i];
        if let Some(d) = c.add_points(a, a) {
            by_key.entry(d.x).or_default().push(format!("2·[{xi}]"));
        }
        for j in (i + 1)..base.len() {
            let (xj, b) = &base[j];
            for neg in [false, true] {
                let bb = if neg { c.negate(b) } else { *b };
                if let Some(s) = c.add_points(a, &bb) {
                    let sign = if neg { '−' } else { '+' };
                    by_key
                        .entry(s.x)
                        .or_default()
                        .push(format!("[{xi}]{sign}[{xj}]"));
                }
            }
        }
    }
    let mut n = 0;
    for v in by_key.values().filter(|v| v.len() > 1) {
        n += 1;
        println!("{}", v.join("  =  ±"));
    }
    println!(
        "{name}: {n} colliding abscissae among {} sums (B = {size}, h = {})",
        base.len() * base.len(),
        c.cofactor
    );
    let _ = raw;
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let name = a.get(1).map_or("sect163k1", |s| s.as_str());
    let size: usize = a.get(2).and_then(|s| s.parse().ok()).unwrap_or(256);
    match secg(name).unwrap().field_bits() {
        0..=128 => run::<2>(name, size),
        129..=192 => run::<3>(name, size),
        193..=256 => run::<4>(name, size),
        _ => panic!("pick a curve of at most 256 bits"),
    }
}

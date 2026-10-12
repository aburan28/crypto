//! The targets of `examples/hamming_nb_sat_panel`, regenerated with the same
//! generator, the same draw order and the library's own curve code, so the
//! x-coordinates are the panel's.  The pieces copied from the panel (it is a
//! standalone example, so nothing can be imported) are `SplitMix64`,
//! `enumerate_weight`, the base construction order and `random_point`; the
//! cross-check against the panel's own output is in `run_ps4.sh`.

use crypto_lib::cryptanalysis::koblitz_rotated_chain::{OnbCurve, OnbPoint, Type2Onb};

/// SplitMix64, as in the panel.
pub struct Rng(pub u64);

impl Rng {
    pub fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^ (z >> 31)
    }
}

/// Every word of weight `≤ w`, in the panel's order.
pub fn enumerate_weight(n: usize, w: usize) -> Vec<u128> {
    let mut out = vec![0u128];
    let mut cur: Vec<u128> = vec![0];
    for _ in 0..w {
        let mut next = Vec::new();
        for &v in &cur {
            let top = if v == 0 {
                0
            } else {
                128 - v.leading_zeros() as usize
            };
            for b in top..n {
                next.push(v | (1u128 << b));
            }
        }
        out.extend(next.iter().copied());
        cur = next;
    }
    out
}

fn random_point(curve: OnbCurve, n: usize, rng: &mut Rng) -> OnbPoint {
    loop {
        let x = (rng.next() as u128) & ((1u128 << n) - 1);
        let pts = curve.points_with_x(x);
        if !pts.is_empty() {
            return pts[(rng.next() & 1) as usize % pts.len()];
        }
    }
}

pub struct PanelTarget {
    pub stratum: &'static str,
    pub x: u64,
    pub y: u64,
}

/// The panel's targets for `--n n --w w --m m --targets t --planted p --seed s`,
/// and the panel's base size (points, both lifts).
pub fn panel_targets(
    n: usize,
    w: usize,
    m: usize,
    targets: usize,
    planted: usize,
    seed: u64,
) -> (Vec<PanelTarget>, usize) {
    let field = Type2Onb::new(n).expect("type-II ONB");
    let curve = OnbCurve::new(field);
    let mut base: Vec<OnbPoint> = Vec::new();
    for x in enumerate_weight(n, w) {
        for p in curve.points_with_x(x) {
            base.push(p);
        }
    }
    let base_size = base.len();
    let mut rng = Rng(seed);
    let mut out = Vec::new();
    for _ in 0..targets {
        let p = random_point(curve, n, &mut rng);
        out.push(PanelTarget {
            stratum: "uniform",
            x: p.x as u64,
            y: p.y as u64,
        });
    }
    for _ in 0..planted {
        let mut sum = OnbPoint::INFINITY;
        for _ in 0..m {
            let p = base[(rng.next() as usize) % base_size];
            sum = curve.add(sum, p);
        }
        if sum.infinity {
            continue;
        }
        out.push(PanelTarget {
            stratum: "planted",
            x: sum.x as u64,
            y: sum.y as u64,
        });
    }
    (out, base_size)
}

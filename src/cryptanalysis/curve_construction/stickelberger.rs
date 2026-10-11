//! Route 2, second family: Fermat quotients `y^m = x^a (1 − x)^b`, cyclic covers
//! of `P¹` branched at three points.
//!
//! The Jacobian factor attached to `(a, b, c)` (`a + b + c ≡ 0 mod m`) has CM
//! type `Φ = {s ∈ (Z/m)^* : ⟨sa/m⟩ + ⟨sb/m⟩ + ⟨sc/m⟩ = 1}` (Stickelberger), and it
//! is isogenous to a power of a CM elliptic curve with CM by `Q(√−7)` — the
//! only way it can carry `A(E_a)` — exactly when `Φ` is induced from
//! `Q(√−7)`, i.e. `Φ = {s : (s | 7) = ε}` for a fixed sign `ε`.  The Klein
//! quartic is `(a, b, c) = (1, 2, 4)` at `m = 7`.  At `m = 7·263 = 1841` the scan
//! decides the whole family.

use num_integer::Integer;
use rayon::prelude::*;
use serde::Serialize;
use std::collections::BTreeMap;

#[derive(Clone, Debug, Serialize)]
pub struct StickelbergerScan {
    pub m: u64,
    pub ordered_pairs_scanned: u64,
    /// Matches by the level `m / gcd(a, b, c, m)` of the factor.
    pub matches_by_level: BTreeMap<u64, u64>,
    pub examples: Vec<(u64, u64, u64, i8)>,
    /// Genus of `y^m = x^a(1−x)^b` when `a, b, a + b` are units mod `m`.
    pub genus_full_level_curve: u64,
}

pub fn scan(m: u64) -> StickelbergerScan {
    let units: Vec<u64> = (1..m).filter(|s| s.gcd(&m) == 1).collect();
    let qr7: Vec<bool> = (0..7u64)
        .map(|r| r != 0 && (1..7u64).any(|x| x * x % 7 == r))
        .collect();
    let per_a: Vec<(BTreeMap<u64, u64>, Vec<(u64, u64, u64, i8)>)> = (1..m)
        .into_par_iter()
        .map(|a| {
            let mut by_level = BTreeMap::new();
            let mut ex = Vec::new();
            for b in 1..m {
                let c = (2 * m - a - b) % m;
                if c == 0 {
                    continue;
                }
                for eps in [1i8, -1] {
                    let ok = units.iter().all(|&s| {
                        let in_phi = s * a % m + s * b % m + s * c % m == m;
                        let r = (s % 7) as usize;
                        in_phi == (if eps == 1 { qr7[r] } else { !qr7[r] })
                    });
                    if ok {
                        let level = m / a.gcd(&b).gcd(&c).gcd(&m);
                        *by_level.entry(level).or_insert(0) += 1;
                        if ex.len() < 4 {
                            ex.push((a, b, c, eps));
                        }
                    }
                }
            }
            (by_level, ex)
        })
        .collect();
    let mut matches_by_level = BTreeMap::new();
    let mut examples = Vec::new();
    for (lv, ex) in per_a {
        for (k, v) in lv {
            *matches_by_level.entry(k).or_insert(0) += v;
        }
        for e in ex {
            if examples.len() < 8 {
                examples.push(e);
            }
        }
    }
    StickelbergerScan {
        m,
        ordered_pairs_scanned: (m - 1) * (m - 1),
        matches_by_level,
        examples,
        genus_full_level_curve: (m - 1) / 2,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn klein_is_found_at_m7() {
        let s = scan(7);
        assert!(s.matches_by_level.get(&7).copied().unwrap_or(0) > 0);
        assert!(s.examples.iter().any(|&(a, b, c, _)| {
            let mut v = [a, b, c];
            v.sort();
            v == [1, 2, 4] || v == [3, 5, 6]
        }));
    }
}

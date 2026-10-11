//! # The matrix-multiplication exponent `ω`, with its provenance.
//!
//! Every Gröbner-basis and linear-algebra cost model in this repository
//! raises a matrix dimension to `ω`.  A bare float hides three facts that
//! decide whether a given `ω` may be used for a given cost:
//!
//! - **construction**: whether the bound comes with an algorithm anyone
//!   runs, one that exists in principle with astronomical constants, or
//!   only a proof that a fast decomposition exists;
//! - **characteristic**: `ω` depends only on the characteristic of the
//!   field (Schönhage), so a binary-curve cost needs the bound in
//!   characteristic 2 and a prime-field cost needs it in that `p`;
//! - **proof status**: whether the proof is published or is a manuscript
//!   nobody here has checked.
//!
//! [`BOUNDS`] lists every `ω` this repository quotes with those facts, and
//! [`OmegaBound::in_characteristic`] says what each is worth in a given
//! characteristic.
//!
//! ## The two corpus bounds (review of 2026-10)
//!
//! The ECDLP review of a mathematics corpus (`ECDLP_REVIEW.md` §6 in
//! `openai/math`, not vendored here) checked two new upper bounds on `ω`
//! for characteristic 2, the case binary-curve index calculus needs:
//!
//! - **`ω < 2.258`** transfers from characteristic 0 by a Nullstellensatz
//!   descent.  The rank schemes live in an unknown number field and the
//!   excluded primes divide a denominator `D` the argument neither
//!   computes nor bounds, so whether 2 (or any specific prime) is excluded
//!   cannot be decided from the paper: [`Applicability::Undetermined`] in
//!   every positive characteristic.
//! - **`ω ≤ 9/4`** (corpus family 107).  Its only characteristic-dependent
//!   step separates the homogeneous components of a polynomial in `ε`
//!   with an `L`-th root of unity, which needs `L` invertible; any `L`
//!   coprime to `p` works, and an extension field supplies the root (`ω`
//!   does not change under field extension).  [`separate_components`] is
//!   that step in its generic form over `GF(2^k)`, and the tests run it on
//!   every field from `GF(2^3)` to `GF(2^18)`.  The bound therefore holds
//!   in every characteristic **if** its characteristic-0 proof is right,
//!   which nobody here has verified: [`Applicability::Conditional`].
//!
//! Neither bound yields an algorithm: both prove that a fast decomposition
//! exists.  So neither changes what any solver in this repository costs;
//! they change only the constant in a heuristic exponent (see
//! [`BinaryIcHeuristic`]), and only conditionally.
//!
//! ## The binary index-calculus heuristic
//!
//! Petit–Quisquater (ASIACRYPT 2012) price Weil-descent index calculus on
//! `E(F_{2^n})` as `2^{(2ω/3)·log(n)·D}` with `D` a first-fall-degree bound
//! at decomposition arity `m ≈ n^{1/3}`; Kousidis–Wiemers (J. Math.
//! Cryptol. 2019) sharpen `D` from `m² + 1` to `m² − m + 1`.  Both rest on
//! the first-fall-degree assumption (`D_reg ≈ D_ff`), which Kosters–Yeo and
//! Huang–Kosters–Yeo (2015) give evidence against.  Reading `log` as
//! `log₁₀` is the only convention that reproduces both papers' turning
//! points (`research/notes/ecc2k130/RESEARCH_ECC2K130_IC_LITERATURE.md`
//! §2), and the tests pin those figures.  The exponent is linear in `ω`,
//! so swapping `ω` scales it by the ratio of the two values.

/// What a bound gives you besides the number.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Construction {
    /// An algorithm that is implemented and used.
    Practical,
    /// An explicit algorithm in principle, never run: its constants are
    /// astronomically large (the laser-method bounds).
    Galactic,
    /// A proof that fast decompositions exist; no algorithm can be written
    /// down from it.
    ExistenceOnly,
}

impl Construction {
    pub fn tag(self) -> &'static str {
        match self {
            Construction::Practical => "practical",
            Construction::Galactic => "galactic",
            Construction::ExistenceOnly => "existence-only",
        }
    }
}

/// Whether the proof behind a bound has been checked.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ProofStatus {
    /// Published and accepted.
    Published,
    /// A manuscript whose proof nobody here has verified.
    Unverified,
}

impl ProofStatus {
    pub fn tag(self) -> &'static str {
        match self {
            ProofStatus::Published => "published",
            ProofStatus::Unverified => "unverified",
        }
    }
}

/// The characteristics a proof covers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CharacteristicScope {
    /// Every characteristic.
    All,
    /// Characteristic 0, and every prime except the divisors of a
    /// denominator the proof does not compute.
    ZeroExceptUncomputedPrimes,
}

/// What a bound is worth in one characteristic.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Applicability {
    /// Proved there by a published proof.
    Established,
    /// Holds there if an unverified proof is right.
    Conditional,
    /// The proof cannot decide this characteristic.
    Undetermined,
}

impl Applicability {
    pub fn tag(self) -> &'static str {
        match self {
            Applicability::Established => "established",
            Applicability::Conditional => "conditional",
            Applicability::Undetermined => "undetermined",
        }
    }
}

/// One upper bound on `ω` and its provenance.
#[derive(Clone, Copy, Debug)]
pub struct OmegaBound {
    /// A stable slug, used in JSON reports.
    pub id: &'static str,
    pub value: f64,
    /// `true` when the bound is `ω < value` rather than `ω ≤ value`.
    pub strict: bool,
    pub construction: Construction,
    pub proof: ProofStatus,
    pub scope: CharacteristicScope,
    pub source: &'static str,
}

impl OmegaBound {
    /// What this bound is worth in characteristic `p` (`0` for
    /// characteristic 0).
    pub fn in_characteristic(&self, p: u64) -> Applicability {
        let covered = match self.scope {
            CharacteristicScope::All => true,
            CharacteristicScope::ZeroExceptUncomputedPrimes => p == 0,
        };
        match (covered, self.proof) {
            (false, _) => Applicability::Undetermined,
            (true, ProofStatus::Published) => Applicability::Established,
            (true, ProofStatus::Unverified) => Applicability::Conditional,
        }
    }

    pub fn relation(&self) -> &'static str {
        if self.strict {
            "<"
        } else {
            "<="
        }
    }
}

/// `ω ≤ 3`: schoolbook multiplication.
pub const SCHOOLBOOK: OmegaBound = OmegaBound {
    id: "schoolbook",
    value: 3.0,
    strict: false,
    construction: Construction::Practical,
    proof: ProofStatus::Published,
    scope: CharacteristicScope::All,
    source: "schoolbook matrix multiplication",
};

/// `ω ≤ log₂ 7`: Strassen.  The repository's Macaulay-cost models use it
/// rounded to `2.807`.
pub const STRASSEN: OmegaBound = OmegaBound {
    id: "strassen-1969",
    value: 2.807_354_922_057_604,
    strict: false,
    construction: Construction::Practical,
    proof: ProofStatus::Published,
    scope: CharacteristicScope::All,
    source: "V. Strassen, Gaussian elimination is not optimal, Numer. Math. 13 (1969)",
};

/// `ω < 2.371339`: the best published bound, a laser-method analysis of
/// the Coppersmith–Winograd tensor.  The tensor's border-rank identity has
/// integer coefficients and the rest of the method is combinatorial, so it
/// holds in every characteristic.
pub const ALMAN_ET_AL_2024: OmegaBound = OmegaBound {
    id: "advwxxz-2024",
    value: 2.371_339,
    strict: true,
    construction: Construction::Galactic,
    proof: ProofStatus::Published,
    scope: CharacteristicScope::All,
    source: "J. Alman, R. Duan, V. Vassilevska Williams, Y. Xu, Z. Xu, R. Zhou, \
             More asymmetry yields faster matrix multiplication, SODA 2025 \
             (arXiv:2404.16349)",
};

/// `ω < 2.258`: a corpus manuscript, transferred from characteristic 0 by
/// a Nullstellensatz descent that leaves the excluded primes uncomputed.
pub const CORPUS_2_258: OmegaBound = OmegaBound {
    id: "corpus-2.258",
    value: 2.258,
    strict: true,
    construction: Construction::ExistenceOnly,
    proof: ProofStatus::Unverified,
    scope: CharacteristicScope::ZeroExceptUncomputedPrimes,
    source: "openai/math corpus; characteristic transfer reviewed in ECDLP_REVIEW.md §6",
};

/// `ω ≤ 9/4`: corpus family 107.  Characteristic-free apart from a
/// root-of-unity separation that works in every characteristic
/// ([`separate_components`]); its characteristic-0 proof is unverified.
pub const CORPUS_FAMILY_107: OmegaBound = OmegaBound {
    id: "corpus-family-107-9/4",
    value: 2.25,
    strict: false,
    construction: Construction::ExistenceOnly,
    proof: ProofStatus::Unverified,
    scope: CharacteristicScope::All,
    source: "openai/math corpus family 107; characteristic-2 step checked in \
             ECDLP_REVIEW.md §6 and in this module's tests",
};

/// Every bound the repository quotes, loosest first.
pub const BOUNDS: [OmegaBound; 5] = [
    SCHOOLBOOK,
    STRASSEN,
    ALMAN_ET_AL_2024,
    CORPUS_2_258,
    CORPUS_FAMILY_107,
];

/// The tightest bound that is [`Applicability::Established`] in
/// characteristic `p`.
pub fn best_established(p: u64) -> OmegaBound {
    BOUNDS
        .iter()
        .filter(|b| b.in_characteristic(p) == Applicability::Established)
        .min_by(|a, b| a.value.total_cmp(&b.value))
        .copied()
        .expect("schoolbook is established everywhere")
}

/// The Petit–Quisquater-type heuristic cost of index calculus on
/// `E(F_{2^n})`, conditional on the first-fall-degree assumption.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BinaryIcHeuristic {
    /// `D = m² + 1` (Petit–Quisquater, ASIACRYPT 2012).
    PetitQuisquater2012,
    /// `D = m² − m + 1` (Kousidis–Wiemers, J. Math. Cryptol. 2019, Thm 3.2).
    KousidisWiemers2019,
}

/// Largest field degree [`BinaryIcHeuristic::turning_point`] searches.
pub const TURNING_POINT_SEARCH_LIMIT: u32 = 1 << 20;

impl BinaryIcHeuristic {
    pub fn id(self) -> &'static str {
        match self {
            BinaryIcHeuristic::PetitQuisquater2012 => "petit-quisquater-2012",
            BinaryIcHeuristic::KousidisWiemers2019 => "kousidis-wiemers-2019",
        }
    }

    pub fn formula(self) -> &'static str {
        match self {
            BinaryIcHeuristic::PetitQuisquater2012 => "log2 T = (2w/3) * log10(n) * (n^(2/3) + 1)",
            BinaryIcHeuristic::KousidisWiemers2019 => {
                "log2 T = (2w/3) * log10(n) * (n^(2/3) - n^(1/3) + 1)"
            }
        }
    }

    /// `log₂` of the heuristic cost at field degree `n` and exponent `omega`.
    pub fn log2_cost(self, n: f64, omega: f64) -> f64 {
        let m = n.cbrt();
        let degree = match self {
            BinaryIcHeuristic::PetitQuisquater2012 => m * m + 1.0,
            BinaryIcHeuristic::KousidisWiemers2019 => m * m - m + 1.0,
        };
        2.0 * omega / 3.0 * n.log10() * degree
    }

    /// The field degree from which the heuristic stays below the generic
    /// `2^{n/2}`: one more than the largest `n` at which it does not.  (At
    /// tiny `n` the formula also dips below `n/2`; that region is not a
    /// crossover and is skipped.)
    pub fn turning_point(self, omega: f64) -> u32 {
        let last_above = (2..=TURNING_POINT_SEARCH_LIMIT)
            .rev()
            .find(|&n| self.log2_cost(n as f64, omega) >= n as f64 / 2.0)
            .expect("the heuristic exceeds n/2 somewhere below the limit");
        assert!(
            last_above < TURNING_POINT_SEARCH_LIMIT,
            "no turning point below {TURNING_POINT_SEARCH_LIMIT} at omega = {omega}"
        );
        last_above + 1
    }
}

/// The root-of-unity separation step in its generic form: from the values
/// `P(ζ^j)`, `j = 0, …, L−1`, of a polynomial `P(ε) = Σ_{i ≤ M} c_i ε^i`
/// over `GF(2^k)`, recover `c_h = L⁻¹ Σ_j ζ^{−jh} P(ζ^j)` for each `h` in
/// `hs`.  `zeta` must have multiplicative order exactly `L` and `L` must be
/// invertible; in characteristic 2 that is `L` odd, `L | 2^k − 1`, and then
/// `L⁻¹ = 1`.  The answer is `c_h` when `M < L`; when `M ≥ L` the
/// coefficients alias (`c_h + c_{h+L} + …`).
pub fn separate_components(
    field: &crate::cryptanalysis::wide_gf2m::WideGf2,
    zeta: u128,
    values: &[u128],
    hs: impl IntoIterator<Item = usize>,
) -> Vec<u128> {
    let l = values.len();
    assert!(
        l % 2 == 1,
        "L must be invertible, so odd, in characteristic 2"
    );
    let zeta_inv = field.inv(zeta);
    hs.into_iter()
        .map(|h| {
            // step = ζ^{−h}; w runs through ζ^{−jh}.
            let mut step = 1u128;
            for _ in 0..h % l {
                step = field.mul(step, zeta_inv);
            }
            let mut w = 1u128;
            let mut acc = 0u128;
            for &v in values {
                acc ^= field.mul(w, v);
                w = field.mul(w, step);
            }
            acc
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::descent_expansion::enumerate_irreducibles;
    use crate::cryptanalysis::wide_gf2m::WideGf2;

    #[test]
    fn the_registry_applies_each_bound_where_its_proof_reaches() {
        for b in [SCHOOLBOOK, STRASSEN, ALMAN_ET_AL_2024] {
            assert_eq!(b.in_characteristic(2), Applicability::Established);
            assert_eq!(b.in_characteristic(0), Applicability::Established);
        }
        // 2.258: characteristic 0 only, and that conditionally.
        assert_eq!(
            CORPUS_2_258.in_characteristic(0),
            Applicability::Conditional
        );
        for p in [2, 3, 5, 7, (1 << 61) - 1] {
            assert_eq!(
                CORPUS_2_258.in_characteristic(p),
                Applicability::Undetermined
            );
        }
        // 9/4: every characteristic, conditionally.
        for p in [0, 2, 3, (1 << 61) - 1] {
            assert_eq!(
                CORPUS_FAMILY_107.in_characteristic(p),
                Applicability::Conditional
            );
        }
        assert_eq!(best_established(2).id, ALMAN_ET_AL_2024.id);
        assert!((STRASSEN.value - 7f64.log2()).abs() < 1e-15);
        assert!(BOUNDS.windows(2).all(|w| w[0].value > w[1].value));
    }

    /// The figures `RESEARCH_ECC2K130_IC_LITERATURE.md` §2 derives at
    /// `ω = log₂ 7`: `2^86.0` at `n = 131`, and the values at each paper's
    /// turning point (`615.9` at `n = 1250`, `986.9` at `n = 2000`).  The
    /// exact crossings sit below the papers' rounded `≈ 1250` and `≈ 2000`.
    #[test]
    fn the_heuristic_reproduces_the_literature_note() {
        let w = STRASSEN.value;
        let kw = BinaryIcHeuristic::KousidisWiemers2019;
        let pq = BinaryIcHeuristic::PetitQuisquater2012;
        assert!((kw.log2_cost(131.0, w) - 86.0).abs() < 0.06);
        assert!((kw.log2_cost(1250.0, w) - 615.9).abs() < 0.1);
        assert!((pq.log2_cost(2000.0, w) - 986.9).abs() < 0.1);
        assert_eq!(kw.turning_point(w), 1144);
        assert_eq!(pq.turning_point(w), 1876);
    }

    /// The turning point is far more sensitive to `ω` than the exponent:
    /// the generic `n/2` and the heuristic are nearly parallel where they
    /// cross, so a 5% smaller constant moves the crossing by a third.
    #[test]
    fn turning_points_for_every_bound() {
        let kw = BinaryIcHeuristic::KousidisWiemers2019;
        let pq = BinaryIcHeuristic::PetitQuisquater2012;
        let got: Vec<(u32, u32)> = BOUNDS
            .iter()
            .map(|b| (kw.turning_point(b.value), pq.turning_point(b.value)))
            .collect();
        assert_eq!(
            got,
            vec![
                (1696, 2584),
                (1144, 1876),
                (357, 802),
                (235, 619),
                (227, 607)
            ]
        );
    }

    /// The exponent is linear in `ω`: `9/4` against the best established
    /// bound is the ~5% the review reported, against Strassen ~20%.
    #[test]
    fn swapping_omega_scales_the_exponent_by_the_ratio() {
        let kw = BinaryIcHeuristic::KousidisWiemers2019;
        for n in [131.0, 163.0, 571.0] {
            let base = kw.log2_cost(n, ALMAN_ET_AL_2024.value);
            let new = kw.log2_cost(n, CORPUS_FAMILY_107.value);
            assert!((new / base - 2.25 / 2.371_339).abs() < 1e-12);
        }
        let drop = 1.0 - CORPUS_FAMILY_107.value / ALMAN_ET_AL_2024.value;
        assert!((0.05..0.055).contains(&drop), "{drop}");
        // A smaller ω moves the turning point down.
        assert!(kw.turning_point(2.25) < kw.turning_point(ALMAN_ET_AL_2024.value));
        assert!(kw.turning_point(ALMAN_ET_AL_2024.value) < kw.turning_point(STRASSEN.value));
    }

    fn pow(f: &WideGf2, mut a: u128, mut e: u128) -> u128 {
        let mut r = 1u128;
        while e > 0 {
            if e & 1 == 1 {
                r = f.mul(r, a);
            }
            a = f.mul(a, a);
            e >>= 1;
        }
        r
    }

    /// An element of order exactly `l` in `GF(2^k)^*`, `l | 2^k − 1`.
    fn root_of_unity(f: &WideGf2, k: u32, l: u128) -> u128 {
        let group = (1u128 << k) - 1;
        let (mut primes, mut r, mut q) = (Vec::new(), l, 2u128);
        while q * q <= r {
            if r.is_multiple_of(q) {
                primes.push(q);
                while r.is_multiple_of(q) {
                    r /= q;
                }
            }
            q += 1;
        }
        if r > 1 {
            primes.push(r);
        }
        (2..group)
            .map(|x| pow(f, x, group / l))
            .find(|&z| primes.iter().all(|q| pow(f, z, l / q) != 1))
            .expect("GF(2^k)^* is cyclic")
    }

    /// `P(ζ^j)` for `j = 0, …, l−1`, by Horner.
    fn evaluate(f: &WideGf2, coeffs: &[u128], zeta: u128, l: u128) -> Vec<u128> {
        let mut e = 1u128;
        (0..l)
            .map(|_| {
                let v = coeffs.iter().rev().fold(0, |acc, &c| f.mul(acc, e) ^ c);
                e = f.mul(e, zeta);
                v
            })
            .collect()
    }

    /// The review's characteristic-2 check: on every `GF(2^k)`,
    /// `3 ≤ k ≤ 18`, for every odd `L | 2^k − 1` with `L ≤ 1023` (and
    /// `L = 2^k − 1` for `k = 13, 17`, where that is prime) and every degree
    /// `M ≤ 7` below `L`, separation recovers each coefficient of a
    /// pseudo-random `P` exactly.  Every `c_h` is checked when `L ≤ 127`;
    /// above that, `h ≤ M + 8`: the coefficients a proof extracts and the
    /// first zeros past them.
    #[test]
    fn separation_recovers_coefficients_over_every_gf2k() {
        let mut state = 0x9E37_79B9_7F4A_7C15u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state as u128
        };
        let mut checked = std::collections::BTreeMap::new();
        for k in 3..=18u32 {
            let f = WideGf2::new(&enumerate_irreducibles(k, 1)[0]);
            let group = (1u128 << k) - 1;
            let mersenne_prime = k == 13 || k == 17;
            for l in (3..=group).step_by(2).filter(|l| group.is_multiple_of(*l)) {
                if l > 1023 && !(mersenne_prime && l == group) {
                    continue;
                }
                let zeta = root_of_unity(&f, k, l);
                for m in 1..=7usize.min(l as usize - 1) {
                    let coeffs: Vec<u128> = (0..=m).map(|_| next() & group).collect();
                    let values = evaluate(&f, &coeffs, zeta, l);
                    let top = if l <= 127 {
                        l as usize
                    } else {
                        (m + 9).min(l as usize)
                    };
                    let got = separate_components(&f, zeta, &values, 0..top);
                    for (h, &c) in got.iter().enumerate() {
                        assert_eq!(
                            c,
                            coeffs.get(h).copied().unwrap_or(0),
                            "k={k} L={l} M={m} h={h}"
                        );
                    }
                    *checked.entry(k).or_insert(0) += 1;
                }
            }
        }
        assert_eq!(
            checked.len(),
            16,
            "every k in 3..=18 had a case: {checked:?}"
        );
    }

    /// `L > M` is needed: at `M = L` the top coefficient aliases onto `c_0`.
    #[test]
    fn separation_aliases_when_the_degree_reaches_l() {
        let k = 6;
        let f = WideGf2::new(&enumerate_irreducibles(k, 1)[0]);
        let l = 7u128;
        let zeta = root_of_unity(&f, k, l);
        let coeffs = [5u128, 9, 0, 0, 0, 0, 0, 33];
        let values = evaluate(&f, &coeffs, zeta, l);
        let got = separate_components(&f, zeta, &values, [0usize, 1]);
        assert_eq!(got, vec![5 ^ 33, 9]);
    }
}

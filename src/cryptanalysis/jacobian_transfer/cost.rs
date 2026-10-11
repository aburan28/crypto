//! The frozen formula-level cost model every transfer route is priced with.
//!
//! Everything is a base-2 logarithm of an estimated bit-operation count, so
//! the source reference (Pollard rho on the elliptic subgroup) and each
//! destination (a Jacobian of genus `g` over a field of `log2 q` bits) are in
//! one unit and the ratio column is a difference of logarithms.  The constants
//! are parameters, written into every report, and are *estimates*: a row is a
//! formula, never a measurement (`AGENTS.md` §2, §5).

use serde::{Deserialize, Serialize};

/// `log2(sqrt(pi / 4))`: Pollard rho's expected iteration count is
/// `sqrt(pi r / 4)` on a group of order `r`.
const LOG2_SQRT_PI_OVER_4: f64 = -0.173_614_7;

/// Parameters of the cost model, frozen into every report.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CostModel {
    /// Exponent `log2 3` of Karatsuba multiplication: a multiplication in a
    /// field of `w` bits is priced at `w^karatsuba_exponent` bit operations.
    pub karatsuba_exponent: f64,
    /// Field multiplications per elliptic group operation (affine, inversion
    /// included; the cover-route protocol's "roughly ten").
    pub elliptic_group_op_mults: f64,
    /// Field multiplications per genus-`g` Jacobian group operation are
    /// `jacobian_group_op_mults_per_genus_squared * g^2` (Cantor's algorithm
    /// is quadratic in the genus; `g = 1` recovers the elliptic constant).
    pub jacobian_group_op_mults_per_genus_squared: f64,
    /// Bits a route must save over the reference to be classed a formula
    /// advantage rather than a boundary row.
    pub advantage_margin_bits: f64,
    /// Hyperelliptic index calculus is priced `g^2 * q^(2 - 2/g)` group
    /// operations (double large primes, Gaudry–Thomé–Thériault–Diem); this is
    /// the exponent of the polynomial factor in `g`.
    pub hyperelliptic_ic_genus_polynomial_exponent: f64,
    /// Enge–Gaudry's `L_N(1/2, c)` constant for a Jacobian of order `N`.
    pub enge_gaudry_c: f64,
    /// Decomposition of a point of `E(F_{p^k})` over a factor base from
    /// `F_p` solves the Semaev system whose Bézout bound is `2^(k(k-1))`; the
    /// solve is priced at that bound times this exponent.
    pub extension_decomposition_bezout_exponent: f64,
}

impl Default for CostModel {
    fn default() -> Self {
        Self {
            karatsuba_exponent: 1.585,
            elliptic_group_op_mults: 10.0,
            jacobian_group_op_mults_per_genus_squared: 10.0,
            advantage_margin_bits: 8.0,
            hyperelliptic_ic_genus_polynomial_exponent: 2.0,
            enge_gaudry_c: std::f64::consts::SQRT_2,
            extension_decomposition_bezout_exponent: 1.0,
        }
    }
}

/// Genus up to which [`CostModel::log2_smoothness_ic`] counts smooth
/// polynomials explicitly (the count is `O(g^2 log g)`).
pub const SMOOTHNESS_COUNT_MAX_GENUS: f64 = 1024.0;

/// `log2` of the number of monic irreducible polynomials of degree `d` over
/// `F_q`, `log2 q = q_bits`: Möbius inversion where it fits a double,
/// `q^d / d` above.
pub fn log2_irreducible_count(q_bits: f64, d: u32) -> f64 {
    if q_bits * f64::from(d) > 900.0 {
        return q_bits * f64::from(d) - f64::from(d).log2();
    }
    let q = q_bits.exp2();
    let mut total = 0.0;
    for e in 1..=d {
        if d % e != 0 {
            continue;
        }
        let mu = mobius(e);
        if mu != 0 {
            total += f64::from(mu) * q.powf(f64::from(d / e));
        }
    }
    (total / f64::from(d)).max(0.0).log2()
}

fn mobius(mut n: u32) -> i32 {
    let mut result = 1;
    let mut p = 2;
    while p * p <= n {
        if n % p == 0 {
            n /= p;
            if n % p == 0 {
                return 0;
            }
            result = -result;
        }
        p += 1;
    }
    if n > 1 {
        result = -result;
    }
    result
}

/// `log2` of the number of monic polynomials of degree `n` over `F_q`
/// whose irreducible factors all have degree `<= b`, for every `b` up to
/// `b_max`: `out[b-1][n]`, `n = 0..=g`.
pub fn log2_smooth_counts(q_bits: f64, g: usize, b_max: u32) -> Vec<Vec<f64>> {
    let mut cur = vec![f64::NEG_INFINITY; g + 1];
    cur[0] = 0.0;
    let mut out = Vec::new();
    for d in 1..=b_max {
        let li = log2_irreducible_count(q_bits, d);
        let i_d = li.exp2();
        let kmax = g / d as usize;
        // log2 C(I + k - 1, k), multisets of size k from I irreducibles.
        let mut cb = vec![0.0f64; kmax + 1];
        for k in 1..=kmax {
            let kf = k as f64;
            let num = if li.is_finite() && li < 50.0 {
                (i_d + kf - 1.0).log2()
            } else {
                li + ((kf - 1.0) / i_d).ln_1p() / std::f64::consts::LN_2
            };
            cb[k] = cb[k - 1] + num - kf.log2();
        }
        let mut next = vec![f64::NEG_INFINITY; g + 1];
        for n in 0..=g {
            let mut acc = f64::NEG_INFINITY;
            for k in 0..=(n / d as usize) {
                let prev = cur[n - k * d as usize];
                if prev.is_finite() && (k == 0 || li.is_finite()) {
                    acc = log2_add(acc, prev + cb[k]);
                }
            }
            next[n] = acc;
        }
        cur = next;
        out.push(cur.clone());
    }
    out
}

/// `log2(g!)`: exact product below 64, Stirling's series above.
pub fn log2_factorial(g: f64) -> f64 {
    if g < 64.0 {
        let n = g.round() as u64;
        return (2..=n).map(|i| (i as f64).log2()).sum();
    }
    let ln = g * g.ln() - g + 0.5 * (2.0 * std::f64::consts::PI * g).ln() + 1.0 / (12.0 * g);
    ln / std::f64::consts::LN_2
}

/// `log2(a + b)` from `log2 a` and `log2 b`.
pub fn log2_add(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let (hi, lo) = if a >= b { (a, b) } else { (b, a) };
    hi + (1.0 + (lo - hi).exp2()).log2()
}

impl CostModel {
    /// `log2` of the bit cost of one multiplication in a field of `field_bits` bits.
    pub fn log2_field_mult(&self, field_bits: f64) -> f64 {
        self.karatsuba_exponent * field_bits.max(1.0).log2()
    }

    /// `log2` of one elliptic group operation over a field of `field_bits` bits.
    pub fn log2_elliptic_op(&self, field_bits: f64) -> f64 {
        self.elliptic_group_op_mults.log2() + self.log2_field_mult(field_bits)
    }

    /// `log2` of one genus-`g` Jacobian group operation over a field of
    /// `field_bits` bits; `log2_genus = log2 g`.
    pub fn log2_jacobian_op(&self, log2_genus: f64, field_bits: f64) -> f64 {
        self.jacobian_group_op_mults_per_genus_squared.log2()
            + 2.0 * log2_genus
            + self.log2_field_mult(field_bits)
    }

    /// Pollard rho on a subgroup of `subgroup_bits` bits of `E(F_Q)`,
    /// `log2 Q = field_bits`: the reference every route is measured against.
    pub fn log2_rho_elliptic(&self, subgroup_bits: f64, field_bits: f64) -> f64 {
        0.5 * subgroup_bits + LOG2_SQRT_PI_OVER_4 + self.log2_elliptic_op(field_bits)
    }

    /// Pollard rho on the image of the subgroup inside a genus-`g` Jacobian
    /// over a field of `field_bits` bits.
    pub fn log2_rho_jacobian(&self, subgroup_bits: f64, log2_genus: f64, field_bits: f64) -> f64 {
        0.5 * subgroup_bits + LOG2_SQRT_PI_OVER_4 + self.log2_jacobian_op(log2_genus, field_bits)
    }

    /// Hyperelliptic index calculus on a genus-`g` Jacobian over `F_q`,
    /// `log2 q = field_bits`, in Jacobian operations: Gaudry's two terms,
    /// relation collection `g! * q` (a random reduced divisor splits over the
    /// degree-one factor base with probability about `1/g!`) and the
    /// double-large-prime linear algebra `g^e * q^(2 - 2/g)`
    /// (Gaudry–Thomé–Thériault–Diem).  The `g!` term is what confines the
    /// method to small genus beside `q`; omitting it prices a genus-64 cover
    /// over `GF(16)` at `2^30` when relation collection alone is `2^300`.
    /// `None` for genus 1 (the formula is for `g >= 2`).
    pub fn log2_hyperelliptic_ic(&self, log2_genus: f64, field_bits: f64) -> Option<f64> {
        if log2_genus < 1.0 - 1e-12 {
            return None;
        }
        let genus = log2_genus.exp2();
        let linear_algebra = (2.0 - 2.0 / genus) * field_bits
            + self.hyperelliptic_ic_genus_polynomial_exponent * log2_genus;
        let relations = log2_factorial(genus) + field_bits;
        Some(log2_add(linear_algebra, relations) + self.log2_jacobian_op(log2_genus, field_bits))
    }

    /// Enge–Gaudry subexponential index calculus, `L_N(1/2, c)` group
    /// operations for a Jacobian of `log2 N = log2_jacobian_order` bits.  Its
    /// hypothesis (genus large beside `log q`) is reported, not enforced.
    pub fn log2_enge_gaudry(&self, log2_jacobian_order: f64, log2_genus: f64, field_bits: f64) -> f64 {
        let ln_n = log2_jacobian_order * std::f64::consts::LN_2;
        if ln_n <= std::f64::consts::E {
            return self.log2_jacobian_op(log2_genus, field_bits);
        }
        let nats = self.enge_gaudry_c * (ln_n * ln_n.ln()).sqrt();
        nats / std::f64::consts::LN_2 + self.log2_jacobian_op(log2_genus, field_bits)
    }

    /// Enge–Gaudry index calculus on a genus-`g` Jacobian over `F_q` priced
    /// by an explicit smoothness count, minimised over the factor-base degree
    /// `b`: factor base `F = (1/2) Σ_(d<=b) I_d` places (about half the
    /// primes of each degree split), relations `F / P_b` Jacobian operations
    /// with `P_b` the fraction of `b`-smooth monic polynomials of degree `g`,
    /// and sparse linear algebra `g F²` operations modulo the subgroup order.
    /// Returns `(total, leading term without the operation constants, b)`;
    /// `None` when `g` exceeds [`SMOOTHNESS_COUNT_MAX_GENUS`].
    pub fn log2_smoothness_ic(
        &self,
        log2_genus: f64,
        field_bits: f64,
        subgroup_bits: f64,
    ) -> Option<(f64, f64, u32)> {
        let genus = log2_genus.exp2();
        if !(2.0 - 1e-9..=SMOOTHNESS_COUNT_MAX_GENUS).contains(&genus) {
            return None;
        }
        let g = genus.round() as usize;
        // Factor bases past 2^160 places are never optimal against these
        // subgroups; stop growing b there.
        let mut b_max = 1u32;
        let mut log2_fb = log2_irreducible_count(field_bits, 1) - 1.0;
        while (b_max as usize) < g && log2_fb < 160.0 {
            b_max += 1;
            log2_fb = log2_add(log2_fb, log2_irreducible_count(field_bits, b_max) - 1.0);
        }
        let counts = log2_smooth_counts(field_bits, g, b_max);
        let total_polys = genus * field_bits;
        let op = self.log2_jacobian_op(log2_genus, field_bits);
        let la_op = self.log2_field_mult(subgroup_bits);
        let mut fb = f64::NEG_INFINITY;
        let mut best: Option<(f64, f64, u32)> = None;
        for b in 1..=b_max {
            fb = log2_add(fb, log2_irreducible_count(field_bits, b) - 1.0);
            let log2_p = counts[(b - 1) as usize][g] - total_polys;
            let relations = fb - log2_p;
            let linear = 2.0 * fb + log2_genus;
            let leading = log2_add(relations, linear);
            let total = log2_add(relations + op, linear + la_op);
            if best.is_none_or(|(t, _, _)| total < t) {
                best = Some((total, leading, b));
            }
        }
        best
    }

    /// Gaudry–Diem index calculus on `E(F_{p^k})` with a factor base from
    /// `F_p`: `p^(2 - 2/k)` relations, each decomposition priced at the Bézout
    /// bound `2^(k(k-1))` of the Semaev system, in elliptic operations over
    /// `F_{p^k}`.
    pub fn log2_extension_ic(&self, p_bits: f64, k: u32) -> f64 {
        let k_f = f64::from(k);
        (2.0 - 2.0 / k_f) * p_bits
            + self.extension_decomposition_bezout_exponent * k_f * (k_f - 1.0)
            + self.log2_elliptic_op(p_bits * k_f)
    }
}

/// The destination solver that attains the minimum of the formulas.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DestinationSolver {
    /// Pollard rho on the carried subgroup inside the Jacobian.
    RhoOnJacobian,
    /// Double-large-prime hyperelliptic index calculus, `q^(2-2/g)`.
    HyperellipticIndexCalculus,
    /// Enge–Gaudry `L(1/2, √2)`, asymptotic formula; used only above
    /// [`SMOOTHNESS_COUNT_MAX_GENUS`] where the explicit count is not run.
    EngeGaudry,
    /// Enge–Gaudry index calculus priced by an explicit smoothness count:
    /// factor base of all places of degree `<= b`, relation probability the
    /// exact fraction of `b`-smooth monic polynomials of degree `g`, sparse
    /// linear algebra of weight `g`; minimised over `b`.
    EngeGaudrySmoothness,
    /// Gaudry–Diem decomposition index calculus on `E(F_{p^k})` itself.
    ExtensionIndexCalculus,
}

/// The cheapest formula for solving on a Jacobian destination.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DestinationCost {
    pub solver: DestinationSolver,
    pub log2_cost: f64,
    /// The leading term of the chosen solver with every constant stripped:
    /// `0.5 log2 r` for rho, `(2 - 2/g) log2 q` for hyperelliptic index
    /// calculus, the `L(1/2)` value for Enge–Gaudry.
    pub log2_leading: f64,
    /// `log2` cost of each formula considered, by solver name.
    pub considered: Vec<(DestinationSolver, f64)>,
    /// Whether the Enge–Gaudry hypothesis `g >= log2 q` held; the formula is
    /// considered only when it does.
    pub enge_gaudry_hypothesis_holds: bool,
}

impl CostModel {
    /// Price the Jacobian solvers that apply and keep the cheapest.
    pub fn jacobian_destination(
        &self,
        subgroup_bits: f64,
        log2_genus: f64,
        field_bits: f64,
    ) -> DestinationCost {
        // (solver, total, leading term)
        let mut considered = vec![(
            DestinationSolver::RhoOnJacobian,
            self.log2_rho_jacobian(subgroup_bits, log2_genus, field_bits),
            0.5 * subgroup_bits,
        )];
        if let Some(ic) = self.log2_hyperelliptic_ic(log2_genus, field_bits) {
            let g = log2_genus.exp2();
            considered.push((
                DestinationSolver::HyperellipticIndexCalculus,
                ic,
                ((2.0 - 2.0 / g) * field_bits).max(log2_factorial(g) + field_bits),
            ));
        }
        let genus = log2_genus.exp2();
        let enge_gaudry_hypothesis_holds = genus >= field_bits;
        if let Some((total, leading, _b)) =
            self.log2_smoothness_ic(log2_genus, field_bits, subgroup_bits)
        {
            considered.push((DestinationSolver::EngeGaudrySmoothness, total, leading));
        } else if enge_gaudry_hypothesis_holds {
            let log2_jacobian_order = genus * field_bits;
            let total = self.log2_enge_gaudry(log2_jacobian_order, log2_genus, field_bits);
            considered.push((
                DestinationSolver::EngeGaudry,
                total,
                total - self.log2_jacobian_op(log2_genus, field_bits),
            ));
        }
        let (solver, log2_cost, log2_leading) = considered
            .iter()
            .copied()
            .fold(None, |best: Option<(DestinationSolver, f64, f64)>, x| match best {
                Some(b) if b.1 <= x.1 => Some(b),
                _ => Some(x),
            })
            .expect("at least rho");
        DestinationCost {
            solver,
            log2_cost,
            log2_leading,
            considered: considered.iter().map(|(s, c, _)| (*s, *c)).collect(),
            enge_gaudry_hypothesis_holds,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rho_reference_scales_as_half_the_subgroup_bits() {
        let m = CostModel::default();
        let a = m.log2_rho_elliptic(160.0, 163.0);
        let b = m.log2_rho_elliptic(200.0, 163.0);
        assert!((b - a - 20.0).abs() < 1e-9);
    }

    #[test]
    fn genus_one_jacobian_is_the_elliptic_curve() {
        let m = CostModel::default();
        assert!((m.log2_jacobian_op(0.0, 100.0) - m.log2_elliptic_op(100.0)).abs() < 1e-12);
        assert!(m.log2_hyperelliptic_ic(0.0, 100.0).is_none());
    }

    #[test]
    fn oakley_group_three_relation_collection_dominates_at_genus_sixteen() {
        // research/cover_route_prime_binary_20261008 priced g = 16 over GF(2^31)
        // with the q^(2-2/g) term alone: 2^58 against rho's 2^76.  Gaudry's
        // relation term g! * q = 2^44.25 * 2^31 = 2^75.25 is larger, so the
        // exponent gain is under one bit, not eighteen.
        let m = CostModel::default();
        let dest = m.jacobian_destination(152.0, 4.0, 31.0);
        let gaudry = dest
            .considered
            .iter()
            .find(|(s, _)| *s == DestinationSolver::HyperellipticIndexCalculus)
            .unwrap()
            .1;
        let relations = log2_factorial(16.0) + 31.0;
        assert!(relations > (2.0 - 2.0 / 16.0) * 31.0 + 17.0);
        assert!((gaudry - relations - m.log2_jacobian_op(4.0, 31.0)).abs() < 0.5);
        // Whichever solver wins, the formulas leave less than the margin.
        let rho = m.log2_rho_elliptic(152.0, 155.0);
        assert!(rho - dest.log2_cost < m.advantage_margin_bits, "rho {rho} vs dest {}", dest.log2_cost);
    }

    #[test]
    fn same_field_genus_two_cover_never_beats_rho() {
        let m = CostModel::default();
        for bits in [64.0, 128.0, 256.0, 521.0] {
            let dest = m.jacobian_destination(bits, 1.0, bits);
            assert!(!dest.enge_gaudry_hypothesis_holds);
            assert!(dest.log2_cost > m.log2_rho_elliptic(bits, bits));
        }
    }

    #[test]
    fn factorial_term_confines_index_calculus_to_small_genus() {
        assert!((log2_factorial(5.0) - 120f64.log2()).abs() < 1e-12);
        let stirling = log2_factorial(64.0);
        let exact: f64 = (2..=64u32).map(|i| f64::from(i).log2()).sum();
        assert!((stirling - exact).abs() < 1e-6, "{stirling} {exact}");
        let m = CostModel::default();
        // Genus 64 over GF(16): relation collection is 64! * 16 ~ 2^300.
        let ic = m.log2_hyperelliptic_ic(6.0, 4.0).unwrap();
        assert!(ic > 290.0, "{ic}");
        // Genus 16 over GF(2^31): 16! * q = 2^75.25 now dominates q^(2-2/g) * g^2 = 2^66.
        let ic16 = m.log2_hyperelliptic_ic(4.0, 31.0).unwrap();
        let relations = log2_factorial(16.0) + 31.0;
        assert!((ic16 - (relations + m.log2_jacobian_op(4.0, 31.0))).abs() < 0.5, "{ic16}");
    }

    #[test]
    fn irreducible_counts_match_necklace_formula() {
        // Over F_2: 2, 1, 2, 3, 6 irreducibles of degree 1..5.
        for (d, n) in [(1u32, 2.0f64), (2, 1.0), (3, 2.0), (4, 3.0), (5, 6.0)] {
            assert!((log2_irreducible_count(1.0, d) - n.log2()).abs() < 1e-9, "{d}");
        }
        assert_eq!(mobius(30), -1);
        assert_eq!(mobius(12), 0);
    }

    #[test]
    fn smooth_counts_are_exact_on_small_fields() {
        // Over F_2, every monic polynomial of degree n is n-smooth: 2^n of them.
        let c = log2_smooth_counts(1.0, 6, 6);
        for n in 0..=6 {
            assert!((c[5][n] - n as f64).abs() < 1e-9, "{n}");
        }
        // 1-smooth monic polynomials of degree n over F_q: C(q + n - 1, n).
        // q = 4, n = 3: C(6, 3) = 20.
        let c = log2_smooth_counts(2.0, 3, 1);
        assert!((c[0][3] - 20f64.log2()).abs() < 1e-9);
        // Degree 4 over F_2 with factors of degree <= 2: x^4 has 3 linear
        // choices... count = all 16 minus the 3 irreducible quartics
        // minus the 2*2 (cubic irreducible)*(linear) products = 9.
        let c = log2_smooth_counts(1.0, 4, 2);
        assert!((c[1][4] - 9f64.log2()).abs() < 1e-9, "{}", c[1][4].exp2());
    }

    #[test]
    fn smoothness_ic_on_genus_sixty_four_over_gf16() {
        // Jacobian of about 2^256 elements.  The explicit count lands below the
        // asymptotic L(1/2, sqrt 2) (2^61.8) and far below rho (2^128), with an
        // interior optimum for the factor-base degree.
        let m = CostModel::default();
        let (total, leading, b) = m.log2_smoothness_ic(6.0, 4.0, 160.0).unwrap();
        let asymptotic = m.log2_enge_gaudry(256.0, 6.0, 4.0) - m.log2_jacobian_op(6.0, 4.0);
        assert!(b >= 2 && b < 20, "b = {b}");
        assert!(leading > 40.0 && leading < asymptotic + 1.0, "explicit {leading} vs asymptotic {asymptotic}");
        assert!(leading < 128.0);
        assert!(total > leading);
    }

    #[test]
    fn log2_add_is_log_sum_exp() {
        assert!((log2_add(3.0, 3.0) - 4.0).abs() < 1e-12);
        assert_eq!(log2_add(f64::NEG_INFINITY, 2.0), 2.0);
        assert!((log2_add(10.0, 0.0) - 10.0) < 0.002);
    }
}

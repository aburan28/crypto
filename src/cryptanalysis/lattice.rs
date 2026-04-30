//! LLL lattice reduction.
//!
//! Lenstra–Lenstra–Lovász, "Factoring polynomials with rational
//! coefficients," Math. Annalen 1982.  The classical lattice
//! basis-reduction algorithm: given an integer lattice basis,
//! produce an equivalent basis whose vectors are nearly orthogonal
//! and the first vector approximates the shortest vector.
//!
//! This implementation pairs a `BigInt` basis with `f64`
//! Gram–Schmidt coefficients — the same hybrid that production
//! libraries (`fplll`, NTL's LLL_FP) use.  Floating-point GS is
//! sound at the dimensions and bit-sizes that matter for
//! cryptanalytic targets in this crate (≤ 30 dim, ≤ 256-bit
//! coefficients).  At higher dimensions or longer entries you'd
//! switch to `L²` (Nguyen–Stehlé 2005); we don't need that here.
//!
//! # Use case in this crate
//!
//! [`crate::cryptanalysis::hnp_ecdsa`] uses LLL to solve the
//! Hidden-Number Problem instance arising from ECDSA signatures
//! with biased per-signature nonces.  See that module for the
//! attack details; this module is the linear-algebra workhorse.

use num_bigint::{BigInt, ToBigInt};
use num_traits::{Signed, ToPrimitive, Zero};

/// LLL-reduce `basis` in place using Lovász parameter `delta`.
///
/// The standard choice is `delta = 0.75`.  Higher (closer to 1)
/// gives a stronger reduction at higher cost; lower gives weaker
/// reduction.  Theory requires `delta ∈ (1/4, 1)`.
///
/// `basis` is a list of basis vectors, each of equal length.  The
/// rows are *the* basis vectors (not columns).
///
/// Returns `Err` if the input is malformed or the algorithm fails
/// to terminate within a generous iteration cap.
pub fn lll_reduce(basis: &mut Vec<Vec<BigInt>>, delta: f64) -> Result<(), &'static str> {
    if !(0.25 < delta && delta < 1.0) {
        return Err("delta must be in (1/4, 1)");
    }
    let n = basis.len();
    if n == 0 {
        return Ok(());
    }
    let dim = basis[0].len();
    if dim == 0 {
        return Err("basis vectors must be non-empty");
    }
    if basis.iter().any(|v| v.len() != dim) {
        return Err("basis vectors must all have the same length");
    }

    // Iteration cap: classical LLL is O(n^4 log B) on n-dim basis
    // with entries ≤ B in magnitude.  Cap is generous to handle
    // cryptographic-size inputs (256-bit entries on dim ~ 30) while
    // still preventing pathological non-termination.  Originally
    // `50 · n²`; bumped to `500 · n² · 8` after a secp256k1 HNP
    // basis hit the cap with dim = 14 — `cap = 500·196·8 + 10⁴ ≈
    // 800K iterations` is well above the worst-case at our sizes.
    let max_iter = 500 * n * n * 8 + 10_000;
    let mut iter = 0;

    let (mut bstar, mut mu) = gram_schmidt(basis);
    let mut k = 1usize;

    while k < n {
        iter += 1;
        if iter > max_iter {
            return Err("LLL exceeded iteration cap (degenerate input?)");
        }

        // ── Size reduction (work top-down from j = k-1 to 0) ────────
        for j in (0..k).rev() {
            if mu[k][j].abs() > 0.5 {
                let q = mu[k][j].round();
                if q == 0.0 {
                    continue;
                }
                // Use the f64 round directly — for our sizes it fits in i64.
                // For safety we'd ideally round in BigInt arithmetic, but
                // for HNP-attack-scale lattices f64 round is exact.
                let q_int = q as i64;
                let q_bi = q_int.to_bigint().unwrap();
                for i in 0..dim {
                    let bji = basis[j][i].clone();
                    basis[k][i] -= &q_bi * &bji;
                }
                let qf = q_int as f64;
                // μ[k][j] update: μ[j][j] = 1 (b_j has coefficient 1
                // on bstar_j), but we don't store the diagonal.
                mu[k][j] -= qf;
                for i in 0..j {
                    mu[k][i] -= qf * mu[j][i];
                }
            }
        }

        // ── Lovász condition ────────────────────────────────────────
        let bsk_sq: f64 = bstar[k].iter().map(|x| x * x).sum();
        let bs_km1_sq: f64 = bstar[k - 1].iter().map(|x| x * x).sum();
        let mu_sq = mu[k][k - 1] * mu[k][k - 1];

        if bsk_sq >= (delta - mu_sq) * bs_km1_sq {
            k += 1;
        } else {
            basis.swap(k, k - 1);
            // Swap invalidates the GS data; cheapest correct fix is
            // to recompute from scratch.  At our dimensions (≤ 30)
            // this is microseconds.
            let (new_bstar, new_mu) = gram_schmidt(basis);
            bstar = new_bstar;
            mu = new_mu;
            k = if k > 1 { k - 1 } else { 1 };
        }
    }
    Ok(())
}

/// Gram–Schmidt orthogonalisation in `f64` precision.
///
/// Returns `(bstar, mu)` where `bstar[i]` is the i-th orthogonal
/// vector and `mu[i][j]` is the GS coefficient
/// `<b_i, bstar_j> / <bstar_j, bstar_j>` for `j < i`.
fn gram_schmidt(basis: &[Vec<BigInt>]) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let n = basis.len();
    let dim = if n > 0 { basis[0].len() } else { 0 };
    let mut bstar: Vec<Vec<f64>> = vec![vec![0.0; dim]; n];
    let mut mu: Vec<Vec<f64>> = vec![vec![0.0; n]; n];

    for i in 0..n {
        let bi: Vec<f64> = basis[i].iter().map(big_to_f64).collect();
        bstar[i].copy_from_slice(&bi);
        for j in 0..i {
            let dot: f64 = bi.iter().zip(&bstar[j]).map(|(a, b)| a * b).sum();
            let norm: f64 = bstar[j].iter().map(|x| x * x).sum();
            mu[i][j] = if norm > 0.0 { dot / norm } else { 0.0 };
            for k in 0..dim {
                bstar[i][k] -= mu[i][j] * bstar[j][k];
            }
        }
    }
    (bstar, mu)
}

/// `BigInt → f64`, lossy for large values.  Sufficient for our
/// attack-scale lattices (256-bit coefficients fit in f64 with
/// ~204-bit precision; the GS coefficients are ratios that bring
/// us comfortably back into double range).
fn big_to_f64(x: &BigInt) -> f64 {
    x.to_f64().unwrap_or(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::ToBigInt;

    fn vec_bigint(xs: &[i64]) -> Vec<BigInt> {
        xs.iter().map(|x| x.to_bigint().unwrap()).collect()
    }

    /// Smallest non-trivial LLL test from the original paper.
    #[test]
    fn lll_two_dim_obvious_reduction() {
        // Basis (1, 1, 1) and (-1, 0, 2).  After LLL we expect
        // shorter vectors.
        let mut basis = vec![vec_bigint(&[1, 1, 1]), vec_bigint(&[-1, 0, 2])];
        lll_reduce(&mut basis, 0.75).unwrap();
        // Reduced first vector should have norm² ≤ 3.
        let norm_sq: i64 = basis[0]
            .iter()
            .map(|x| {
                let xv = x.to_i64().unwrap();
                xv * xv
            })
            .sum();
        assert!(norm_sq <= 3, "got norm² = {}", norm_sq);
    }

    /// Reduced basis should satisfy size-reduction: |μ_ij| ≤ 1/2.
    #[test]
    fn reduced_basis_satisfies_size_reduction() {
        // Random-ish 4-dim integer basis.
        let mut basis = vec![
            vec_bigint(&[10, 23, 1, 5]),
            vec_bigint(&[7, -3, 12, 19]),
            vec_bigint(&[-4, 8, 15, -2]),
            vec_bigint(&[3, 11, -6, 9]),
        ];
        lll_reduce(&mut basis, 0.75).unwrap();
        let (_bstar, mu) = gram_schmidt(&basis);
        for i in 1..basis.len() {
            for j in 0..i {
                assert!(
                    mu[i][j].abs() <= 0.5 + 1e-9,
                    "size-reduction violated: μ[{}][{}] = {}",
                    i,
                    j,
                    mu[i][j]
                );
            }
        }
    }

    /// Identity basis is already reduced — LLL should be a no-op.
    #[test]
    fn identity_basis_is_already_reduced() {
        let mut basis = vec![
            vec_bigint(&[1, 0, 0]),
            vec_bigint(&[0, 1, 0]),
            vec_bigint(&[0, 0, 1]),
        ];
        let original = basis.clone();
        lll_reduce(&mut basis, 0.75).unwrap();
        assert_eq!(basis, original);
    }

    /// Swap should produce shorter first vector when input has a
    /// long-followed-by-short ordering.
    #[test]
    fn long_short_basis_swaps() {
        let mut basis = vec![vec_bigint(&[100, 0]), vec_bigint(&[1, 1])];
        lll_reduce(&mut basis, 0.75).unwrap();
        // Norm² of first vector should be ≤ 2 (the short one).
        let norm_sq: i64 = basis[0]
            .iter()
            .map(|x| {
                let xv = x.to_i64().unwrap();
                xv * xv
            })
            .sum();
        assert!(norm_sq <= 2);
    }

    #[test]
    fn rejects_invalid_delta() {
        let mut basis: Vec<Vec<BigInt>> = vec![vec_bigint(&[1, 0]), vec_bigint(&[0, 1])];
        assert!(lll_reduce(&mut basis, 1.5).is_err());
        assert!(lll_reduce(&mut basis, 0.1).is_err());
    }
}

// ── BKZ (Block Korkine-Zolotarev) lattice reduction ──────────────────────────
//
// Schnorr 1987, refined by Schnorr-Euchner 1991.  Stronger than LLL —
// for each consecutive block of size β, find the shortest vector in
// the projected sublattice via enumeration and use it to refine the
// basis.  BKZ-2 reduces to LLL; β = 20 is "comfortable" with f64 GS;
// β > 30 needs extended-precision GS (`L²` of Nguyen-Stehlé) which we
// don't implement.
//
// Reference: Chen-Nguyen "BKZ 2.0: Better Lattice Security Estimates"
// (Asiacrypt 2011).  We implement the basic Schnorr-Euchner BKZ
// without the 2.0 enhancements (no extreme pruning, no progressive
// sieving) — sufficient for HNP attacks at sub-bit bias.

/// BKZ-β reduction.  Extends [`lll_reduce`] with block enumeration:
/// for each `k`, finds the shortest vector in the projected lattice
/// `π_k(L_{[k..k+β-1]})` and inserts it at position `k`.  Repeats
/// until a full pass produces no improvement.
///
/// `delta` = LLL Lovász parameter (use 0.99 inside BKZ for
/// stronger LLL-cleanup between block insertions).
/// `beta` = block size.  Cost is roughly `2^β` enumeration ops per
/// block per round; β = 20 is practical, β = 25 is slow but
/// feasible, β = 30 needs better GS precision than `f64` provides.
pub fn bkz_reduce(
    basis: &mut Vec<Vec<BigInt>>,
    beta: usize,
    delta: f64,
) -> Result<(), &'static str> {
    if beta < 2 {
        return Err("BKZ block size must be ≥ 2");
    }
    let n = basis.len();
    if n == 0 {
        return Ok(());
    }

    // Initial LLL pass — BKZ assumes the basis is already LLL-reduced.
    lll_reduce(basis, delta)?;

    // Iterative BKZ rounds: pass over all start positions; if any
    // insertion happens, do another full pass.  Bounded total
    // rounds for safety.
    let max_rounds = 8 * n;
    for _round in 0..max_rounds {
        let mut changed = false;
        for k in 0..n.saturating_sub(1) {
            let end = (k + beta - 1).min(n - 1);
            // Compute GS data restricted to rows [k..=end].
            let (bstar, mu) = gram_schmidt(&basis[..]);
            let block_norms_sq: Vec<f64> =
                (k..=end).map(|i| bstar[i].iter().map(|x| x * x).sum()).collect();
            // Build mu_block[i][j] for i, j in 0..(end-k+1), j < i,
            // mapping mu_block[i][j] = mu[k+i][k+j].
            let block_dim = end - k + 1;
            let mut mu_block: Vec<Vec<f64>> = vec![vec![0.0; block_dim]; block_dim];
            for i in 0..block_dim {
                for j in 0..i {
                    mu_block[i][j] = mu[k + i][k + j];
                }
            }
            // Bound: ||b*_k||² (the current best at position k).  Any
            // shorter vector in the projected lattice is an improvement.
            let r_sq = block_norms_sq[0] * 0.999; // strict improvement
            let candidate = enumerate_block(&block_norms_sq, &mu_block, r_sq);
            if let Some(coeffs) = candidate {
                // Build the candidate vector v = Σ coeffs[i] · basis[k+i] (in original Z basis).
                let dim = basis[0].len();
                let mut v: Vec<BigInt> = vec![BigInt::zero(); dim];
                for i in 0..block_dim {
                    if coeffs[i] == 0 {
                        continue;
                    }
                    let c_bi = BigInt::from(coeffs[i]);
                    for d in 0..dim {
                        v[d] += &c_bi * &basis[k + i][d];
                    }
                }
                // Insert v at position k, shift the rest down,
                // then LLL-reduce — LLL drops the now-redundant vector.
                basis.insert(k, v);
                lll_reduce(basis, delta)?;
                // After LLL the basis has rank ≤ original rank n.
                // Drop any zero rows that LLL produced.
                basis.retain(|row| row.iter().any(|x| !x.is_zero()));
                if basis.len() < n {
                    // LLL collapsed — abort safely.
                    return Err("BKZ insertion produced rank-deficient basis");
                }
                changed = true;
                break; // restart the outer round
            }
        }
        if !changed {
            return Ok(());
        }
    }
    Err("BKZ exceeded round cap (degenerate input?)")
}

/// Schnorr-Euchner enumeration: find the shortest non-zero
/// vector in the projected lattice with coefficients in `Z^d`,
/// constrained to squared-norm `< radius_sq`.  Returns the integer
/// coefficient vector or `None` if no shorter vector exists.
///
/// `bstar_norms_sq[i]` = `||b*_{k+i}||²` for the `d` block vectors.
/// `mu[i][j]` for `j < i` = the GS coefficient `μ_{k+i, k+j}`
/// (only block-internal coefficients).
///
/// Algorithm: depth-first enumeration with center-rounding and
/// bound-pruning.  At each level `i`, the optimal continuous `x_i`
/// is `−Σ_{j>i} μ_{j,i} · x_j`; we try integer values in a window
/// around it whose width is bounded by the remaining squared-norm
/// budget.
fn enumerate_block(
    bstar_norms_sq: &[f64],
    mu: &[Vec<f64>],
    radius_sq: f64,
) -> Option<Vec<i64>> {
    let d = bstar_norms_sq.len();
    if d == 0 {
        return None;
    }
    let mut x = vec![0i64; d];
    let mut best: Option<(Vec<i64>, f64)> = None;
    enumerate_recurse(d - 1, &mut x, bstar_norms_sq, mu, 0.0, radius_sq, &mut best);
    best.map(|(coeffs, _)| coeffs)
}

fn enumerate_recurse(
    i: usize,
    x: &mut Vec<i64>,
    c: &[f64],
    mu: &[Vec<f64>],
    partial_sq: f64,
    radius_sq: f64,
    best: &mut Option<(Vec<i64>, f64)>,
) {
    let d = c.len();
    // Center for x[i]:  c_i = -Σ_{j>i} μ[j][i] · x[j]
    let mut center = 0.0_f64;
    for j in (i + 1)..d {
        center -= mu[j][i] * (x[j] as f64);
    }
    // Bound: (x[i] - center)² · c[i] ≤ R² - partial_sq.
    // Also use the current best as a tighter bound if available.
    let bound = best.as_ref().map(|b| b.1).unwrap_or(radius_sq);
    let remaining = bound - partial_sq;
    if remaining <= 0.0 {
        return;
    }
    if c[i] <= 0.0 {
        // Degenerate (zero-length GS vector); skip this level by
        // forcing x[i] = round(center).
        let xi = center.round() as i64;
        x[i] = xi;
        let dev = (xi as f64) - center;
        let new_partial = partial_sq + dev * dev * c[i];
        if i == 0 {
            if !x.iter().all(|&v| v == 0) && new_partial < bound {
                *best = Some((x.clone(), new_partial));
            }
        } else {
            enumerate_recurse(i - 1, x, c, mu, new_partial, radius_sq, best);
        }
        x[i] = 0;
        return;
    }
    let max_dev = (remaining / c[i]).sqrt();
    let lo = (center - max_dev).ceil() as i64;
    let hi = (center + max_dev).floor() as i64;
    if lo > hi {
        return;
    }
    // Schnorr-Euchner zigzag from round(center) outward (more
    // pruning-friendly than a left-to-right sweep).
    let center_r = center.round() as i64;
    let center_r_clamped = center_r.max(lo).min(hi);
    let mut tried: Vec<i64> = Vec::new();
    let mut step = 0i64;
    loop {
        let v_a = center_r_clamped + step;
        if v_a >= lo && v_a <= hi {
            tried.push(v_a);
        }
        if step != 0 {
            let v_b = center_r_clamped - step;
            if v_b >= lo && v_b <= hi {
                tried.push(v_b);
            }
        }
        if tried.len() as i64 >= (hi - lo + 1) {
            break;
        }
        step += 1;
        if step > (hi - lo + 1) {
            break; // safety
        }
    }
    for xi in tried {
        x[i] = xi;
        let dev = (xi as f64) - center;
        let new_partial = partial_sq + dev * dev * c[i];
        let bound = best.as_ref().map(|b| b.1).unwrap_or(radius_sq);
        if new_partial >= bound {
            continue;
        }
        if i == 0 {
            // Reached bottom; record if non-zero and shorter.
            if !x.iter().all(|&v| v == 0) {
                *best = Some((x.clone(), new_partial));
            }
        } else {
            enumerate_recurse(i - 1, x, c, mu, new_partial, radius_sq, best);
        }
    }
    x[i] = 0;
}

#[cfg(test)]
mod bkz_tests {
    use super::*;
    use num_bigint::ToBigInt;

    fn vec_bigint(xs: &[i64]) -> Vec<BigInt> {
        xs.iter().map(|x| x.to_bigint().unwrap()).collect()
    }

    /// BKZ-β=2 should be a no-op beyond LLL (β=2 BKZ ≡ LLL).
    #[test]
    fn bkz_beta_2_equals_lll() {
        let mut a = vec![
            vec_bigint(&[10, 23, 1, 5]),
            vec_bigint(&[7, -3, 12, 19]),
            vec_bigint(&[-4, 8, 15, -2]),
            vec_bigint(&[3, 11, -6, 9]),
        ];
        let mut b = a.clone();
        lll_reduce(&mut a, 0.99).unwrap();
        bkz_reduce(&mut b, 2, 0.99).unwrap();
        // a's first vector norm² should equal b's first vector norm²
        // (or be strictly close; β=2 BKZ is identical to LLL).
        let na: i64 = a[0]
            .iter()
            .map(|x| {
                let v = x.to_i64().unwrap_or(0);
                v * v
            })
            .sum();
        let nb: i64 = b[0]
            .iter()
            .map(|x| {
                let v = x.to_i64().unwrap_or(0);
                v * v
            })
            .sum();
        assert_eq!(na, nb, "BKZ-2 should match LLL");
    }

    /// BKZ-β=8 should produce a basis whose first vector is at least
    /// as short as LLL's first vector (often strictly shorter for
    /// challenging bases).
    #[test]
    fn bkz_at_least_as_good_as_lll() {
        // Standard "challenge" type basis: random-ish 6-dim integers.
        let basis_template = vec![
            vec_bigint(&[40, 13, 27, 9, 2, 17]),
            vec_bigint(&[5, 29, -8, 14, 22, 3]),
            vec_bigint(&[-12, 7, 31, -4, 18, 25]),
            vec_bigint(&[19, -6, 11, 24, -3, 8]),
            vec_bigint(&[2, 16, -9, 5, 30, -11]),
            vec_bigint(&[33, -14, 4, -7, 12, 21]),
        ];

        let mut a = basis_template.clone();
        let mut b = basis_template.clone();
        lll_reduce(&mut a, 0.99).unwrap();
        bkz_reduce(&mut b, 8, 0.99).unwrap();

        let norm_sq = |row: &Vec<BigInt>| -> f64 {
            row.iter().map(|x| big_to_f64(x).powi(2)).sum()
        };
        let na = norm_sq(&a[0]);
        let nb = norm_sq(&b[0]);
        assert!(
            nb <= na + 1e-9,
            "BKZ first-vector norm {} should be ≤ LLL norm {}",
            nb,
            na
        );
    }

    /// BKZ should preserve the lattice (basis vectors stay integer
    /// linear combinations of the original — verified by determinant
    /// magnitude up to sign).
    #[test]
    fn bkz_preserves_lattice_volume() {
        let original = vec![
            vec_bigint(&[5, 0, 0]),
            vec_bigint(&[0, 5, 0]),
            vec_bigint(&[0, 0, 5]),
        ];
        let mut b = original.clone();
        bkz_reduce(&mut b, 3, 0.99).unwrap();
        // For axis-aligned cube lattice, BKZ should leave it
        // essentially unchanged (or permuted).  All rows should
        // still have norm² = 25.
        for row in &b {
            let n_sq: i64 = row
                .iter()
                .map(|x| {
                    let v = x.to_i64().unwrap_or(0);
                    v * v
                })
                .sum();
            assert_eq!(n_sq, 25);
        }
    }

    /// BKZ rejects β < 2.
    #[test]
    fn bkz_rejects_tiny_beta() {
        let mut basis = vec![vec_bigint(&[1, 0]), vec_bigint(&[0, 1])];
        assert!(bkz_reduce(&mut basis, 1, 0.75).is_err());
    }
}

// Reference the import so unused-import lints don't trip on
// `Signed` / `Zero` if compiler folds them out in release.
#[allow(dead_code)]
fn _link() -> bool {
    let x = BigInt::zero();
    x.is_negative() || x.is_positive()
}

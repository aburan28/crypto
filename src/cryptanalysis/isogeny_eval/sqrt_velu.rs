//! √élu: `Õ(√ℓ)` isogeny evaluation on x-only Montgomery curves
//! (Bernstein, De Feo, Leroux, Smith, "Faster computation of isogenies of
//! large prime degree", ANTS 2020).
//!
//! The odd multiples `S = {1, 3, …, ℓ−2}` of the kernel generator `P`
//! represent the kernel up to sign.  They are split as `S = (I ± J) ∪ K`
//! with
//!
//! * `b = ⌊√(ℓ−1)/2⌋`, `b' = ⌊(ℓ−1)/(4b)⌋`;
//! * `I = {2b(2i+1) : 0 ≤ i < b'}`, `J = {2j+1 : 0 ≤ j < b}`;
//! * `K = {4bb'+1, …, ℓ−2}` (odd).
//!
//! For each `α` the product `Π_{s∈S}(α − x([s]P))` is, up to the common
//! nonzero factor `Δ = Π_{i,j}(x_i − x_j)²`, the resultant of `Π_{i∈I}(Z −
//! x_i)` and the biquadratic `Π_{j∈J} F_α(Z, x_j)`, computed by a product
//! tree on `I` and a remainder tree; `K` is handled directly.  The codomain
//! uses `α = ±1` in [`super::velu`]'s Edwards formula; images use
//! `x(φ(Q)) = x · (Π(1 − x x_s))² / (Π(x − x_s))²`, where `Δ` cancels.
//!
//! With Karatsuba multiplication the measured cost exponent was 0.82–0.86
//! in `ℓ` over `2^10 ≤ ℓ ≤ 2^24` (round 4, lead 16.3); the asymptotic
//! `Õ(√ℓ)` needs fast (NTT/FFT) multiplication, not implemented here.

use super::counted_fp::{OpCounter, PrimeField};
use super::montgomery::{batch_to_affine, ladder, xadd, xdbl, MontgomeryCurve, XZ};
use super::poly::{product, ProductTree};
use super::velu::{check_degree, edwards_codomain, velu, IsogenyEvaluation};
use super::IsogenyEvalError;

/// The index-set sizes √élu used.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SqrtVeluShape {
    /// `b = ⌊√(ℓ−1)/2⌋` (0 when `ℓ = 3` and plain Vélu was used).
    pub b: u64,
    /// `b' = ⌊(ℓ−1)/(4b)⌋`.
    pub b_prime: u64,
    /// `|I| = b'`.
    pub i_len: usize,
    /// `|J| = b`.
    pub j_len: usize,
    /// `|K|`, the leftover odd multiples.
    pub k_len: usize,
}

/// `⌊√n⌋`.
fn isqrt(n: u64) -> u64 {
    if n < 2 {
        return n;
    }
    let mut x = (n as f64).sqrt() as u64;
    while x.checked_mul(x).is_none_or(|v| v > n) {
        x -= 1;
    }
    while (x + 1).checked_mul(x + 1).is_some_and(|v| v <= n) {
        x += 1;
    }
    x
}

/// `x([start + 2i]P)` for `i < count` (`start` odd ≥ 1).
fn odd_multiples<C: OpCounter>(
    f: &PrimeField<C>,
    curve: &MontgomeryCurve,
    p1: &XZ<u64>,
    start: u64,
    count: usize,
) -> Vec<XZ<u64>> {
    if count == 0 {
        return Vec::new();
    }
    let p2 = xdbl(f, p1, &curve.a24);
    let mut out = Vec::with_capacity(count);
    out.push(ladder(f, p1, start as u128, &curve.a24));
    if count > 1 {
        out.push(ladder(f, p1, (start + 2) as u128, &curve.a24));
    }
    for i in 2..count {
        let next = xadd(f, &out[i - 1], &p2, &out[i - 2]);
        out.push(next);
    }
    out
}

/// The precomputed √élu data for one kernel.
struct Plan {
    tree: ProductTree,
    /// For each `x_j`, `(x_j, x_j² + 1 + 2A x_j)`.
    jc: Vec<(u64, u64)>,
    xk: Vec<u64>,
}

impl Plan {
    /// `Δ·Π_{s∈S}(α − x_s)` and, if `need_num`, `Δ·Π_{s∈S}(1 − α x_s)`.
    fn eval<C: OpCounter>(&mut self, f: &PrimeField<C>, alpha: u64, need_num: bool) -> (u64, u64) {
        let a2 = f.sqr(alpha);
        let mut quads_den = Vec::with_capacity(self.jc.len());
        let mut quads_num = Vec::with_capacity(self.jc.len());
        for &(x, mid) in &self.jc {
            // c2 = (α − x)², c1 = −2(α²x + α·mid + x), c0 = (αx − 1)².
            let amx = f.sub(alpha, x);
            let c2 = f.sqr(amx);
            let axm1 = f.sub(f.mul(alpha, x), 1);
            let c0 = f.sqr(axm1);
            let inner = f.add(f.add(f.mul(a2, x), f.mul(alpha, mid)), x);
            let c1 = f.neg(f.mul_small(2, inner));
            quads_den.push(vec![c0, c1, c2]);
            if need_num {
                quads_num.push(vec![c2, c1, c0]);
            }
        }
        let ej_den = product(f, quads_den);
        let den = self
            .tree
            .evaluate(f, &ej_den)
            .into_iter()
            .fold(1u64, |acc, v| f.mul(acc, v));
        let den = self
            .xk
            .iter()
            .fold(den, |acc, &x| f.mul(acc, f.sub(alpha, x)));
        let mut num = 0u64;
        if need_num {
            let ej_num = product(f, quads_num);
            num = self
                .tree
                .evaluate(f, &ej_num)
                .into_iter()
                .fold(1u64, |acc, v| f.mul(acc, v));
            num = self
                .xk
                .iter()
                .fold(num, |acc, &x| f.mul(acc, f.sub(1, f.mul(alpha, x))));
        }
        (den, num)
    }
}

/// Evaluate the `ℓ`-isogeny with kernel `⟨P⟩`, `x(P) = kernel_x`, at the
/// x-coordinates `points`, by √élu.  For `ℓ = 3` (no useful split) it runs
/// [`velu`] and reports an all-zero shape apart from `k_len = 1`.
///
/// Preconditions and errors are those of [`velu`]; in addition
/// [`IsogenyEvalError::NotAKernelPoint`] is returned if a multiple
/// `[s]P`, `s ≤ ℓ − 2`, is at infinity (so `P` does not have order `ℓ`).
/// Agrees with [`velu`] on the codomain and on every image.
pub fn sqrt_velu<C: OpCounter>(
    f: &PrimeField<C>,
    curve: &MontgomeryCurve,
    ell: u64,
    kernel_x: u64,
    points: &[u64],
) -> Result<(IsogenyEvaluation, SqrtVeluShape), IsogenyEvalError> {
    check_degree(ell)?;
    let b = isqrt(ell - 1) / 2;
    if b == 0 {
        let out = velu(f, curve, ell, kernel_x, points)?;
        let shape = SqrtVeluShape {
            k_len: 1,
            ..SqrtVeluShape::default()
        };
        return Ok((out, shape));
    }
    let bp = (ell - 1) / (4 * b);
    let p1 = XZ::affine(kernel_x);
    // I = {2b(2i+1)}, via Q_{i+1} = Q_i + [4b]P with difference Q_{i−1}.
    let q0 = ladder(f, &p1, (2 * b) as u128, &curve.a24);
    let q4 = xdbl(f, &q0, &curve.a24);
    let mut ipts = Vec::with_capacity(bp as usize);
    ipts.push(q0.clone());
    if bp > 1 {
        ipts.push(xadd(f, &q0, &q4, &q0));
    }
    for i in 2..bp as usize {
        let next = xadd(f, &ipts[i - 1], &q4, &ipts[i - 2]);
        ipts.push(next);
    }
    let jpts = odd_multiples(f, curve, &p1, 1, b as usize);
    let kstart = 4 * b * bp + 1;
    let kcount = if kstart <= ell - 2 {
        ((ell - 2 - kstart) / 2 + 1) as usize
    } else {
        0
    };
    let kpts = odd_multiples(f, curve, &p1, kstart, kcount);
    let not_kernel = IsogenyEvalError::NotAKernelPoint;
    let xi = batch_to_affine(f, &ipts).ok_or(not_kernel.clone())?;
    let xj = batch_to_affine(f, &jpts).ok_or(not_kernel.clone())?;
    let xk = batch_to_affine(f, &kpts).ok_or(not_kernel)?;
    let a = curve.a;
    let jc: Vec<(u64, u64)> = xj
        .iter()
        .map(|&x| {
            let x2 = f.sqr(x);
            let mid = f.add(f.add(x2, 1), f.mul(f.mul_small(2, a), x));
            (x, mid)
        })
        .collect();
    let shape = SqrtVeluShape {
        b,
        b_prime: bp,
        i_len: xi.len(),
        j_len: xj.len(),
        k_len: xk.len(),
    };
    let mut plan = Plan {
        tree: ProductTree::new(f, &xi),
        jc,
        xk,
    };
    // den(−1) = Δ Π(−1 − x_s), den(1) = Δ Π(1 − x_s): equal to Π(x_s + 1)
    // and Π(x_s − 1) up to Δ and signs, which the eighth powers remove.
    let (dm1, _) = plan.eval(f, f.neg(1), false);
    let (d1, _) = plan.eval(f, 1, false);
    let codomain = edwards_codomain(f, curve, ell, dm1, d1)?;
    let images = points
        .iter()
        .map(|&x| {
            let (den, num) = plan.eval(f, x, true);
            let xn = f.mul(x, f.sqr(num));
            let zn = f.sqr(den);
            XZ { x: xn, z: zn }.to_affine(f)
        })
        .collect();
    Ok((IsogenyEvaluation { codomain, images }, shape))
}

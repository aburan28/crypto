//! lit2026 W11 (a): counted field operations per visited point for the three
//! BSGS step kernels ([`crypto_lib::cryptanalysis::bsgs_fast::StepKernel`]),
//! on the solver's own chain layout, for a curve given on the command line.
//!
//! ```sh
//! cargo run --release --example lit2026_w11_bsgs_opcount -- \
//!     <p> <a> <b> <gx> <gy> <n> <label>
//! ```
//!
//! For each kernel and window `K` it runs `W = 8` chains of `len` points of
//! the baby phase (addend `G`) and of the giant phase (addend `−(2m−1)G`
//! from a target), single-threaded on the counting field, and prints the
//! multiplications, squarings and inversions per visited point, the
//! inversion's own cost in multiplications (Fermat: counted once by running
//! it on the counting field), and whether every kernel visited the same
//! `(index, x)` pairs as the original chain.

use num_bigint::BigUint;

use crypto_lib::cryptanalysis::bsgs_fast::{
    count_kernel_ops, BsgsFastPlan, FastCurve, FastPoint, StepKernel,
};
use crypto_lib::ecc::curve::CurveParams;

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    assert!(a.len() >= 7, "usage: p a b gx gy n label");
    let u = |i: usize| -> u64 { a[i].parse().expect("integer") };
    let params = CurveParams {
        name: "lit2026-w11",
        p: BigUint::from(u(0)),
        a: BigUint::from(u(1)),
        b: BigUint::from(u(2)),
        gx: BigUint::from(u(3)),
        gy: BigUint::from(u(4)),
        n: BigUint::from(u(5)),
        h: 1,
    };
    let label = &a[6];
    let curve = FastCurve::from_params(&params).expect("curve");
    let plan = BsgsFastPlan::new(curve.n, 0, true, 1, 8);
    let s = curve.scalar_mul(curve.g, plan.stride);
    let neg_s = curve.neg(s);
    let target = curve.scalar_mul(curve.g, curve.n / 3 + 12345);

    // Cost of one inversion in multiplications (Fermat a^{p−2}).
    let fermat = {
        // `inv` is counted as one I; its internal pow is not.  Count the
        // pow's squarings and multiplications from the exponent.
        let e = curve.f.modulus() - 2;
        let sq = 63 - e.leading_zeros() as u64;
        let ml = e.count_ones() as u64;
        (sq, ml)
    };

    let len = 1u64 << 14;
    let chains = 8usize;
    for (phase, base, addend) in [("baby", None, curve.g), ("giant", Some(target), neg_s)] {
        let base = base.unwrap_or(FastPoint::INFINITY);
        let reference =
            count_kernel_ops(&curve, StepKernel::Chain, 1, chains, len, base, addend, 1);
        for (kernel, ks) in [
            (StepKernel::Chain, vec![1usize]),
            (StepKernel::Paired, vec![1, 2, 4, 8, 16, 32, 64]),
            (StepKernel::CoZ, vec![1, 2, 4, 8, 16, 32, 64]),
        ] {
            for k in ks {
                let c = count_kernel_ops(&curve, kernel, k, chains, len, base, addend, 1);
                let v = c.visited as f64;
                let m = c.ops.mul as f64 / v;
                let sq = c.ops.sqr as f64 / v;
                let inv = c.ops.inv as f64 / v;
                // M-equivalent: S = M, I = Fermat squarings + multiplications.
                let m_eq = m + sq + inv * (fermat.0 + fermat.1) as f64;
                println!(
                    "{{\"label\":\"{label}\",\"phase\":\"{phase}\",\"kernel\":\"{}\",\"pair_k\":{k},\"chains\":{chains},\"len\":{len},\"visited\":{},\"slow_rounds\":{},\"mul\":{},\"sqr\":{},\"inv\":{},\"add\":{},\"mul_per_point\":{m:.4},\"sqr_per_point\":{sq:.4},\"inv_per_point\":{inv:.6},\"fermat_sqr\":{},\"fermat_mul\":{},\"m_equiv_per_point\":{m_eq:.4},\"same_points_as_chain\":{}}}",
                    kernel.name(),
                    c.visited,
                    c.slow_rounds,
                    c.ops.mul,
                    c.ops.sqr,
                    c.ops.inv,
                    c.ops.add,
                    fermat.0,
                    fermat.1,
                    c.digest == reference.digest && c.visited == reference.visited
                );
            }
        }
    }
}

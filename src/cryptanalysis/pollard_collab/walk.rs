//! The shared pseudo-random walk and the distinguished-point records
//! it produces.
//!
//! Every peer runs exactly this function on exactly the branch table
//! derived in [`super::job`], so trails that meet anywhere in the
//! group stay merged until the next distinguished point — which is
//! what lets a collision between *any* two peers' walkers solve the
//! instance.

use num_bigint::BigUint;
use serde::{Deserialize, Serialize};

use super::job::{hex_of, parse_hex, JobContext};
use crate::cryptanalysis::ecdlp_variants::{add_mod, neg_mod};
use crate::ecc::point::Point;

/// A distinguished point reached by one walker, as shipped in a
/// check-in.  Self-certifying: anyone can check
/// `a·P + b·Q == (x, y)` and that `(x, y)` is distinguished.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct DpRecord {
    /// Walker index whose trail ended here (audit handle: re-run it).
    pub walker: u64,
    /// Steps from the walker's start to this point.
    pub steps: u64,
    /// DP coordinates (hex).
    pub x: String,
    pub y: String,
    /// Coefficients (hex, mod `n`) with `a·P + b·Q = (x, y)`.
    pub a: String,
    pub b: String,
}

impl DpRecord {
    pub fn point(&self, ctx: &JobContext) -> Result<Point, String> {
        Ok(Point::Affine {
            x: crate::ecc::field::FieldElement::new(parse_hex(&self.x)?, ctx.p.clone()),
            y: crate::ecc::field::FieldElement::new(parse_hex(&self.y)?, ctx.p.clone()),
        })
    }

    pub fn coefficients(&self) -> Result<(BigUint, BigUint), String> {
        Ok((parse_hex(&self.a)?, parse_hex(&self.b)?))
    }

    /// Full verification: parse, on-curve, distinguished, and
    /// `a·P + b·Q` equals the point.  Cost: two scalar
    /// multiplications — cheap next to the `2^dp_bits` steps the
    /// record represents, so peers verify every record they accept.
    pub fn verify(&self, ctx: &JobContext) -> Result<Point, String> {
        let pt = self.point(ctx)?;
        if !ctx.on_curve(&pt) {
            return Err("DP not on curve".into());
        }
        if !ctx.is_dp(&pt) {
            return Err("point is not distinguished".into());
        }
        let (a, b) = self.coefficients()?;
        if a >= ctx.n || b >= ctx.n {
            return Err("coefficients not reduced mod n".into());
        }
        if ctx.combine(&a, &b) != pt {
            return Err("a·P + b·Q ≠ claimed point".into());
        }
        Ok(pt)
    }
}

/// How a single walker's trail ended.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum WalkerOutcome {
    /// Reached a distinguished point.
    Dp(DpRecord),
    /// Hit the step cap (fruitless cycle or unlucky sparse region).
    DeadTrail { steps: u64 },
}

/// Branch selector from the low 64 bits of `x`.
fn branch_index(pt: &Point, r: usize) -> usize {
    match pt {
        Point::Infinity => 0,
        Point::Affine { x, .. } => {
            let lo = x.value.iter_u64_digits().next().unwrap_or(0);
            (lo % r as u64) as usize
        }
    }
}

/// Is `pt` the upper `±` representative, `2y > p`?  Compared as
/// `y > ⌊p/2⌋` against the job's precomputed [`JobContext::half_p`].
fn needs_flip(pt: &Point, half_p: &BigUint) -> bool {
    match pt {
        Point::Affine { y, .. } => y.value > *half_p,
        Point::Infinity => false,
    }
}

/// One step of the walk.  Pure in `(R, a, b)` except for the 2-cycle
/// escape under the negation map, which also looks at the previous
/// point.
pub fn step(
    ctx: &JobContext,
    r: &Point,
    a: &BigUint,
    b: &BigUint,
    last: Option<&Point>,
) -> (Point, BigUint, BigUint) {
    let n = &ctx.n;
    let idx = branch_index(r, ctx.branches.len());
    let br = &ctx.branches[idx];
    let mut np = r.add_vartime(&br.point, &ctx.a);
    let mut na = add_mod(a, &br.u, n);
    let mut nb = add_mod(b, &br.v, n);
    if ctx.spec.negation_map {
        if last == Some(&np) {
            // Fruitless 2-cycle A→B→A: break symmetry by doubling.
            np = r.double_vartime(&ctx.a);
            na = (a << 1) % n;
            nb = (b << 1) % n;
        }
        if needs_flip(&np, &ctx.half_p) {
            np = np.neg_vartime();
            na = neg_mod(&na, n);
            nb = neg_mod(&nb, n);
        }
    }
    (np, na, nb)
}

/// Run walker `i` from its derived start until a DP or the step cap.
pub fn run_walker(ctx: &JobContext, i: u64) -> WalkerOutcome {
    let (mut a, mut b, mut r) = ctx.walker_start(i);
    if ctx.spec.negation_map && needs_flip(&r, &ctx.half_p) {
        r = r.neg_vartime();
        a = neg_mod(&a, &ctx.n);
        b = neg_mod(&b, &ctx.n);
    }
    let mut last: Option<Point> = None;
    let mut steps = 0u64;
    loop {
        if ctx.is_dp(&r) {
            let (x, y) = match &r {
                Point::Affine { x, y } => (hex_of(&x.value), hex_of(&y.value)),
                Point::Infinity => unreachable!("∞ is never distinguished"),
            };
            return WalkerOutcome::Dp(DpRecord {
                walker: i,
                steps,
                x,
                y,
                a: hex_of(&a),
                b: hex_of(&b),
            });
        }
        if steps >= ctx.step_cap {
            return WalkerOutcome::DeadTrail { steps };
        }
        let (np, na, nb) = step(ctx, &r, &a, &b, last.as_ref());
        last = Some(std::mem::replace(&mut r, np));
        a = na;
        b = nb;
        steps += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::super::job::{demo_curve, JobSpec};
    use super::*;
    use num_traits::Zero;

    fn ctx(negation: bool) -> JobContext {
        let curve = demo_curve("demo-mid").unwrap();
        let q = curve
            .generator()
            .scalar_mul(&BigUint::from(4242u32), &curve.a_fe());
        let mut spec = JobSpec::new(&curve, &q, "walk", 1).unwrap();
        spec.dp_bits = 4;
        spec.negation_map = negation;
        spec.build().unwrap()
    }

    #[test]
    fn walker_records_verify_in_both_modes() {
        for negation in [false, true] {
            let c = ctx(negation);
            let mut dps = 0;
            for i in 0..64u64 {
                if let WalkerOutcome::Dp(rec) = run_walker(&c, i) {
                    dps += 1;
                    let pt = rec.verify(&c).expect("own record verifies");
                    assert!(c.is_dp(&pt));
                    assert_eq!(rec.walker, i);
                }
            }
            assert!(
                dps > 40,
                "negation={negation}: only {dps}/64 walkers reached a DP"
            );
        }
    }

    #[test]
    fn tampered_records_are_rejected() {
        let c = ctx(false);
        let rec = (0..64u64)
            .find_map(|i| match run_walker(&c, i) {
                WalkerOutcome::Dp(r) => Some(r),
                _ => None,
            })
            .unwrap();
        let mut bad = rec.clone();
        bad.a = hex_of(&((parse_hex(&rec.a).unwrap() + 1u32) % &c.n));
        assert!(bad.verify(&c).unwrap_err().contains("≠"));
        let mut bad = rec.clone();
        bad.x = "1".into();
        assert!(bad.verify(&c).is_err());
        let mut bad = rec.clone();
        bad.a = hex_of(&(&c.n + parse_hex(&rec.a).unwrap()));
        assert!(bad.verify(&c).unwrap_err().contains("reduced"));
    }

    #[test]
    fn same_walker_same_trail() {
        let c1 = ctx(true);
        let c2 = ctx(true);
        for i in 0..16u64 {
            assert_eq!(run_walker(&c1, i), run_walker(&c2, i));
        }
    }

    #[test]
    fn step_keeps_the_invariant() {
        for negation in [false, true] {
            let c = ctx(negation);
            let (mut a, mut b, mut r) = c.walker_start(9);
            if negation && needs_flip(&r, &c.half_p) {
                r = r.neg();
                a = neg_mod(&a, &c.n);
                b = neg_mod(&b, &c.n);
            }
            let mut last = None;
            for _ in 0..200 {
                let (np, na, nb) = step(&c, &r, &a, &b, last.as_ref());
                assert_eq!(c.combine(&na, &nb), np, "a·P + b·Q drifted");
                last = Some(std::mem::replace(&mut r, np));
                a = na;
                b = nb;
            }
        }
    }

    /// `add_mod`, `neg_mod`, `needs_flip` and [`JobContext::is_dp`] return
    /// exactly what the expressions they replaced did, on one-word values
    /// and on multi-word ones (their `BigUint` fallbacks).
    #[test]
    fn step_helpers_match_the_expressions_they_replace() {
        use num_bigint::RandBigInt;
        use rand::SeedableRng;
        let mut rng = rand::rngs::StdRng::seed_from_u64(0xc011ab);
        let zero = BigUint::zero();
        for n in [
            BigUint::from(2u32),
            BigUint::from(7u32),
            BigUint::from(u64::MAX - 58),
            BigUint::from(1u32) << 64,
            BigUint::from(5u32).pow(60),
        ] {
            let mut vals = vec![zero.clone(), BigUint::from(1u32) % &n, &n - 1u32];
            vals.extend((0..40).map(|_| rng.gen_biguint_below(&n)));
            let half = &n >> 1;
            for v in &vals {
                let old_neg = if v.is_zero() { zero.clone() } else { &n - v };
                assert_eq!(neg_mod(v, &n), old_neg, "n={n} v={v}");
                let y = Point::Affine {
                    x: crate::ecc::field::FieldElement::zero(n.clone()),
                    y: crate::ecc::field::FieldElement {
                        value: v.clone(),
                        modulus: n.clone(),
                    },
                };
                assert_eq!(needs_flip(&y, &half), (v << 1) > n, "n={n} y={v}");
                for w in vals.iter().step_by(7) {
                    assert_eq!(add_mod(v, w, &n), (v + w) % &n, "n={n} {v}+{w}");
                }
            }
            let wide = BigUint::from(u64::MAX);
            assert_eq!(add_mod(&wide, &wide, &n), (&wide + &wide) % &n);
        }
        let mut c = ctx(true);
        let big = BigUint::from(5u32).pow(60);
        for bits in [0u32, 1, 10, 63, 64, 65, 130] {
            c.dp_mask = (BigUint::from(1u32) << bits) - 1u32;
            for i in 0..60u32 {
                let x = (rng.gen_biguint_below(&big) >> (i % 70)) << (i % 70);
                let p = Point::Affine {
                    x: crate::ecc::field::FieldElement {
                        value: x.clone(),
                        modulus: big.clone(),
                    },
                    y: crate::ecc::field::FieldElement::zero(big.clone()),
                };
                assert_eq!(c.is_dp(&p), (&x & &c.dp_mask).is_zero(), "{bits} {x}");
            }
        }
    }
}

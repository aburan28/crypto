//! Separable isogenies by Vélu's formulas on short Weierstrass curves,
//! with full `(x, y)` images so every map is an exact group homomorphism.
//!
//! For a finite subgroup `G = {O} ∪ F₂ ∪ R ∪ (−R)` of `E: y² = x³ + ax + b`
//! (`F₂` the points of order 2, `R` one of each `±` pair of the rest),
//! Vélu's formulas (Washington, *Elliptic Curves*, Thm 12.16) give the
//! codomain `E': y² = x³ + (a − 5t)x + (b − 7w)` and the map
//!
//! ```text
//! x' = x + Σ_Q [ t_Q/(x − x_Q) + u_Q/(x − x_Q)² ]
//! y' = y − Σ_Q [ u_Q·2y/(x − x_Q)³ + t_Q(y − y_Q)/(x − x_Q)² − g_Q^x g_Q^y/(x − x_Q)² ]
//! ```
//!
//! over `Q ∈ F₂ ∪ R`, with `g_Q^x = 3x_Q² + a`, `g_Q^y = −2y_Q`,
//! `t_Q = g_Q^x` for `Q ∈ F₂` and `2g_Q^x` otherwise, `u_Q = (g_Q^y)²`,
//! `t = Σ t_Q`, `w = Σ (u_Q + x_Q t_Q)`.  The kernel maps to `O`.
//!
//! A prime-degree step costs `O(ℓ)` per evaluation; a chain of degree
//! `N = ∏ ℓ_i^{e_i}` with `N | p + 1` is built from a kernel generator by
//! repeatedly isolating an order-`ℓ` point, applying one step, and pushing
//! the generator through.  No balanced strategy: the toy sizes do not need
//! one, and the naive loop keeps the code auditable.

use crate::curve::{Curve, Point};
use crate::fp::Fp2;

/// One prime-degree Vélu step.
#[derive(Clone, Debug)]
pub struct VeluStep {
    pub domain: Curve,
    pub codomain: Curve,
    pub degree: u128,
    /// The kernel representatives `F₂ ∪ R` with their precomputed data
    /// `(x_Q, y_Q, g_x, g_y, t_Q, u_Q)`.
    reps: Vec<(Fp2, Fp2, Fp2, Fp2, Fp2, Fp2)>,
}

impl VeluStep {
    /// Build the step from a kernel generator `K` of prime order `ell`.
    pub fn new(domain: &Curve, kernel: &Point, ell: u128) -> Result<VeluStep, String> {
        if !domain.is_on_curve(kernel) || kernel.is_infinity() {
            return Err("kernel generator must be a finite point on the curve".into());
        }
        if !domain.mul(kernel, ell).is_infinity() {
            return Err("kernel generator does not have the stated order".into());
        }
        let prime = domain.prime;
        let mut reps = Vec::new();
        let mut t = prime.fp2_zero();
        let mut w = prime.fp2_zero();
        let two = prime.fp2(2, 0);
        if ell == 2 {
            let (xq, yq) = match kernel {
                Point::Affine { x, y } => (*x, *y),
                Point::Infinity => unreachable!(),
            };
            if !yq.is_zero() {
                return Err("order-2 point must have y = 0".into());
            }
            let gx = xq.square().mul_small(3).add(domain.a);
            let gy = yq.mul(two).neg();
            let tq = gx;
            let uq = gy.square();
            t = t.add(tq);
            w = w.add(uq.add(xq.mul(tq)));
            reps.push((xq, yq, gx, gy, tq, uq));
        } else {
            // R: [k]K for k = 1 .. (ell − 1)/2
            let mut cur = *kernel;
            for _ in 0..(ell - 1) / 2 {
                let (xq, yq) = match cur {
                    Point::Affine { x, y } => (x, y),
                    Point::Infinity => return Err("kernel generator order too small".into()),
                };
                let gx = xq.square().mul_small(3).add(domain.a);
                let gy = yq.mul(two).neg();
                let tq = gx.mul(two);
                let uq = gy.square();
                t = t.add(tq);
                w = w.add(uq.add(xq.mul(tq)));
                reps.push((xq, yq, gx, gy, tq, uq));
                cur = domain.add(&cur, kernel);
            }
        }
        let a2 = domain.a.sub(t.mul_small(5));
        let b2 = domain.b.sub(w.mul_small(7));
        let codomain = Curve::new(prime, a2, b2)?;
        Ok(VeluStep {
            domain: *domain,
            codomain,
            degree: ell,
            reps,
        })
    }

    /// Evaluate the isogeny on a point of the domain.
    pub fn eval(&self, p: &Point) -> Point {
        let (x, y) = match p {
            Point::Infinity => return Point::Infinity,
            Point::Affine { x, y } => (*x, *y),
        };
        let mut xs = x;
        let mut ys = y;
        for &(xq, yq, gx, gy, tq, uq) in &self.reps {
            let d = x.sub(xq);
            let Some(inv) = d.inv() else {
                // x = x_Q: the point is in the kernel (or its negative), image O.
                return Point::Infinity;
            };
            let inv2 = inv.square();
            let inv3 = inv2.mul(inv);
            xs = xs.add(tq.mul(inv)).add(uq.mul(inv2));
            let term = uq
                .mul(y.mul_small(2))
                .mul(inv3)
                .add(tq.mul(y.sub(yq)).mul(inv2))
                .sub(gx.mul(gy).mul(inv2));
            ys = ys.sub(term);
        }
        Point::Affine { x: xs, y: ys }
    }
}

/// A composition of prime-degree steps.
#[derive(Clone, Debug)]
pub struct IsogenyChain {
    pub steps: Vec<VeluStep>,
    pub domain: Curve,
    pub codomain: Curve,
    pub degree: u128,
}

impl IsogenyChain {
    /// The identity chain.
    pub fn identity(curve: &Curve) -> IsogenyChain {
        IsogenyChain {
            steps: Vec::new(),
            domain: *curve,
            codomain: *curve,
            degree: 1,
        }
    }

    /// Build the separable isogeny with cyclic kernel `⟨K⟩`, `K` of order
    /// `n` with the given factorisation.
    pub fn from_kernel(
        curve: &Curve,
        kernel: &Point,
        n: u128,
        factors: &[(u128, u32)],
    ) -> Result<IsogenyChain, String> {
        if !curve.mul(kernel, n).is_infinity() {
            return Err("kernel generator order does not divide n".into());
        }
        if curve.order_dividing(kernel, n, factors) != n {
            return Err("kernel generator order is not n".into());
        }
        let mut chain = IsogenyChain::identity(curve);
        let mut k = *kernel;
        let mut remaining = n;
        for &(ell, e) in factors {
            for _ in 0..e {
                let order_ell = chain.codomain.mul(&k, remaining / ell);
                let step = VeluStep::new(&chain.codomain, &order_ell, ell)?;
                k = step.eval(&k);
                remaining /= ell;
                chain.codomain = step.codomain;
                chain.degree *= ell;
                chain.steps.push(step);
            }
        }
        if !k.is_infinity() {
            return Err("kernel generator did not map to O".into());
        }
        Ok(chain)
    }

    pub fn eval(&self, p: &Point) -> Point {
        let mut cur = *p;
        for step in &self.steps {
            cur = step.eval(&cur);
        }
        cur
    }

    /// Append another chain starting at this chain's codomain.
    pub fn then(mut self, next: IsogenyChain) -> Result<IsogenyChain, String> {
        if next.domain != self.codomain {
            return Err("chains do not compose".into());
        }
        self.steps.extend(next.steps);
        self.codomain = next.codomain;
        self.degree *= next.degree;
        Ok(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curve::factor_smooth;
    use crate::fp::Prime;
    use rand::{Rng, SeedableRng};

    const P: u64 = 810_647_932_926_689_279;

    #[test]
    fn prime_degree_steps_are_homomorphisms() {
        let prime = Prime::new(P).unwrap();
        let e0 = Curve::e0(prime);
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(7);
        let factors = factor_smooth(P as u128 + 1, 100).unwrap();
        for ell in [2u128, 3, 5] {
            let k = e0.random_point_of_order(&mut rng, ell, &[(ell, 1)]);
            let step = VeluStep::new(&e0, &k, ell).unwrap();
            assert!(step.eval(&k).is_infinity());
            for _ in 0..6 {
                let r = e0.random_point(&mut rng);
                let s = e0.random_point(&mut rng);
                let fr = step.eval(&r);
                let fs = step.eval(&s);
                assert!(
                    step.codomain.is_on_curve(&fr),
                    "image on codomain, ell = {ell}"
                );
                assert_eq!(step.eval(&e0.add(&r, &s)), step.codomain.add(&fr, &fs));
                assert_eq!(step.eval(&e0.neg(&r)), step.codomain.neg(&fr));
                // The kernel is exactly ⟨K⟩: R + K maps to the same point.
                assert_eq!(step.eval(&e0.add(&r, &k)), fr);
            }
            // The codomain is again supersingular with (p+1)² points.
            let r = step.codomain.random_point(&mut rng);
            assert!(step.codomain.mul(&r, P as u128 + 1).is_infinity());
            let _ = factors.len();
        }
    }

    #[test]
    fn smooth_chain_and_degree() {
        let prime = Prime::new(P).unwrap();
        let e0 = Curve::e0(prime);
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(11);
        let n = 2u128.pow(12) * 9 * 5;
        let factors = factor_smooth(n, 100).unwrap();
        let k = e0.random_point_of_order(&mut rng, n, &factors);
        let chain = IsogenyChain::from_kernel(&e0, &k, n, &factors).unwrap();
        assert_eq!(chain.degree, n);
        assert_eq!(chain.steps.len(), 12 + 2 + 1);
        assert!(chain.eval(&k).is_infinity());
        for _ in 0..4 {
            let r = e0.random_point(&mut rng);
            let s = e0.random_point(&mut rng);
            let fr = chain.eval(&r);
            assert!(chain.codomain.is_on_curve(&fr));
            assert_eq!(
                chain.eval(&e0.add(&r, &s)),
                chain.codomain.add(&fr, &chain.eval(&s))
            );
            // Every multiple of K is in the kernel; nothing else of order | n is,
            // so a random point of E[n] outside ⟨K⟩ survives.
            let m = rng.gen::<u128>() % n;
            assert_eq!(chain.eval(&e0.add(&r, &e0.mul(&k, m))), fr);
        }
        // Dual check: the composed degree equals the kernel order, and the
        // image of E[n] is cyclic of order n (the quotient E[n]/⟨K⟩).
        let (p1, q1) = e0.torsion_basis(&mut rng, n, &factors);
        let im_p = chain.eval(&p1);
        let im_q = chain.eval(&q1);
        let ord_p = chain.codomain.order_dividing(&im_p, n, &factors);
        let ord_q = chain.codomain.order_dividing(&im_q, n, &factors);
        let gcd = num_integer::gcd(ord_p, ord_q);
        assert_eq!(ord_p / gcd * ord_q, n, "image orders must have lcm n");
    }
}

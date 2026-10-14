//! `prime-kummer:t=T[,g=G]`: Petit–Kosters–Messeng's factor base for prime
//! curves (PKC 2016, §3) in its Kummer form, as an `ic_framework` plug-in.
//!
//! `L(x) = (x/g)^{2^t} − 1` is the composition of `t` squarings, so its
//! zero set `V = g·μ_{2^t}` has exactly `2^t` elements in `F_p` when
//! `2^t | p − 1` (condition 1 of PKM §3.1) and is the zero set of
//! low-degree maps (condition 2).  The base is every point of the curve
//! whose abscissa is in `V`, both signs, one column an abscissa.  `g`
//! defaults to 1 (`V = μ_{2^t}`).
//!
//! Building it costs `2^t` Legendre symbols and square roots; that work is
//! counted in the factor-base phase (`legendre_symbols`, `sqrt_solves`),
//! as `prime-abscissa`'s is.

use super::stages::{FactorBaseBuilder, InstanceCtx, Params};
use crate::cryptanalysis::ic_boundary::{
    CountedGroup, FactorBase, GroupOps, PrimeCurve, PrimeInstance, PrimePoint,
};
use crate::cryptanalysis::residual_walk::pow_mod;

fn mul(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}

/// A primitive `2^t`-th root of unity in `F_p`, from the first small
/// base that gives one; `None` when `2^t ∤ p − 1`.
pub fn primitive_two_power_root(p: u64, t: u32) -> Option<u64> {
    if t == 0 || t >= 63 || !(p - 1).is_multiple_of(1u64 << t) {
        return None;
    }
    let e = (p - 1) >> t;
    (2..p.min(1 << 20)).find_map(|h| {
        let z = pow_mod(h, e, p);
        (pow_mod(z, 1u64 << (t - 1), p) != 1).then_some(z)
    })
}

/// The Kummer tower base.
pub struct PrimeKummerBase<'i> {
    pub instance: &'i PrimeInstance,
}

impl FactorBaseBuilder<PrimeCurve> for PrimeKummerBase<'_> {
    fn name(&self) -> &str {
        "prime-kummer"
    }

    fn describe(&self, params: &Params) -> String {
        format!(
            "PKM's Kummer tower: every point with x ∈ g·μ_{{2^t}}, the zeros of (x/g)^(2^t) − 1, t = {}, \
             g = {}, one column an abscissa",
            params.get("t").unwrap_or("?"),
            params.get("g").unwrap_or("1")
        )
    }

    fn parameters(&self) -> &[(&str, &str)] {
        &[
            ("t", "tower height: |V| = 2^t; needs 2^t | p − 1"),
            ("g", "the coset representative, default 1"),
        ]
    }

    fn build(
        &self,
        _ctx: &InstanceCtx<PrimeCurve>,
        params: &Params,
        _ops: &mut GroupOps,
    ) -> Result<FactorBase<PrimePoint>, String> {
        let curve = &self.instance.curve;
        let p = curve.p;
        let t = params.u64("t")? as u32;
        let g = params.u64_or("g", 1)? % p;
        if g == 0 {
            return Err("g must be non-zero".into());
        }
        let zeta = primitive_two_power_root(p, t)
            .ok_or_else(|| format!("2^{t} does not divide p − 1 = {}", p - 1))?;
        let mut xs: Vec<u64> = Vec::with_capacity(1 << t);
        let mut z = g;
        for _ in 0..(1u64 << t) {
            xs.push(z);
            z = mul(z, zeta, p);
        }
        xs.sort_unstable();
        let (mut points, mut col_of, mut coef_of) = (Vec::new(), Vec::new(), Vec::new());
        let (mut legendre, mut sqrts) = (0u64, 0u64);
        let mut columns = 0usize;
        for &x in &xs {
            let rhs = curve.rhs(x);
            legendre += 1;
            if curve.legendre(rhs) < 0 {
                continue;
            }
            sqrts += 1;
            let Some(y) = curve.sqrt(rhs) else { continue };
            let pt = PrimePoint::affine(x, y.min(if y == 0 { 0 } else { p - y }));
            let neg = curve.neg(pt);
            points.push(pt);
            col_of.push(columns);
            coef_of.push(1);
            if neg != pt {
                points.push(neg);
                col_of.push(columns);
                coef_of.push(self.instance.r - 1);
            }
            columns += 1;
        }
        if columns == 0 {
            return Err("no abscissa of the tower carries a point".into());
        }
        let mut fb = FactorBase::from_column_map(
            format!(
                "Kummer tower t = {t}, g = {g}, {} of 2^{t} abscissae on the curve",
                columns
            ),
            points,
            col_of,
            coef_of,
            columns,
            |q| curve.key(q),
            |q| curve.key(&curve.neg(*q)),
            |q| q.x,
        )?;
        fb.cost.count("legendre_symbols", legendre);
        fb.cost.count("sqrt_solves", sqrts);
        Ok(fb)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cryptanalysis::ic_boundary::prime_instance_for;

    #[test]
    fn the_tower_is_the_two_power_roots_of_unity() {
        // 2013265921 = 15·2^27 + 1, the PKM rounds' p₁.
        let p = 2_013_265_921u64;
        let z = primitive_two_power_root(p, 10).unwrap();
        assert_eq!(pow_mod(z, 1 << 10, p), 1);
        assert_ne!(pow_mod(z, 1 << 9, p), 1);
        assert!(primitive_two_power_root(p, 28).is_none());
    }

    #[test]
    fn the_base_holds_exactly_the_tower_points() {
        // Any small prime curve with 2^t | p − 1 for some t ≥ 2.
        let mut seed = 1;
        let (inst, t) = loop {
            let inst = prime_instance_for(16, seed);
            let t = (inst.curve.p - 1).trailing_zeros().min(6);
            if t >= 3 {
                break (inst, t);
            }
            seed += 1;
        };
        let ctx = InstanceCtx {
            group: &inst.curve,
            generator: inst.generator_point(),
            target: inst.generator_point(),
            r: inst.r,
            cofactor: inst.cofactor,
            group_order: inst.group_order,
            name: inst.name.clone(),
            field_degree: None,
        };
        let (_, params) = Params::parse_spec(&format!("prime-kummer:t={t}")).unwrap();
        let fb = PrimeKummerBase { instance: &inst }
            .build(&ctx, &params, &mut GroupOps::default())
            .unwrap();
        let p = inst.curve.p;
        for q in &fb.points {
            assert!(inst.curve.is_on_curve(*q));
            assert_eq!(pow_mod(q.x, 1 << t, p), 1, "abscissa outside μ_2^t");
        }
        assert!(fb.abscissae <= 1 << t);
    }
}

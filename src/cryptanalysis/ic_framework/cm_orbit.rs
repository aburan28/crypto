//! Prime-field scalar-orbit controls and CM self-isogeny orbit bases.
//!
//! The scalar arm works on both conductor strata. The Velu arm is a
//! low-degree implementation of the same scalar action, not an extra
//! generic-group automorphism credited to the DLP search.
use std::collections::HashSet;
use std::time::Instant;

use super::stages::{FactorBaseBuilder, InstanceCtx, Params};
use crate::cryptanalysis::ic_boundary::{
    fast_map, CountedGroup, FactorBase, FieldOps, GroupOps, PhaseCost, PrimeCurve, PrimeInstance,
    PrimePoint,
};

fn mul(a: u64, b: u64, p: u64) -> u64 {
    (a as u128 * b as u128 % p as u128) as u64
}
fn add(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 + b as u128) % p as u128) as u64
}
fn sub(a: u64, b: u64, p: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        p - (b - a)
    }
}
fn pow(mut a: u64, mut k: u64, p: u64) -> u64 {
    let mut x = 1;
    while k > 0 {
        if k & 1 != 0 {
            x = mul(x, a, p);
        }
        a = mul(a, a, p);
        k >>= 1;
    }
    x
}
fn isqrt(x: u128) -> u64 {
    let mut n = (x as f64).sqrt() as u64;
    while (n as u128 + 1).pow(2) <= x {
        n += 1;
    }
    while (n as u128).pow(2) > x {
        n -= 1;
    }
    n
}

#[derive(Clone, Debug)]
pub struct CmMap {
    pub lambda: u64,
    pub degree: u64,
    pub kernel_x: Vec<u64>,
    pub kernel_t: Vec<u64>,
    pub kernel_rhs: Vec<u64>,
    pub inv_c: u64,
    pub inv_cu: u64,
}

fn counted_mul(a: u64, b: u64, p: u64, f: &mut FieldOps) -> u64 {
    f.muls += 1;
    mul(a, b, p)
}
fn counted_pow(mut a: u64, mut k: u64, p: u64, f: &mut FieldOps) -> u64 {
    let mut x = 1;
    while k > 0 {
        if k & 1 != 0 {
            x = counted_mul(x, a, p, f);
        }
        f.sqrs += 1;
        a = mul(a, a, p);
        k >>= 1;
    }
    x
}
impl CmMap {
    /// Normalized Velu x map and its derivative determine both coordinates.
    /// Kernel abscissae lie in F_p, even when kernel ordinates do not.
    pub fn apply(
        &self,
        c: &PrimeCurve,
        point: PrimePoint,
        f: &mut FieldOps,
    ) -> Result<PrimePoint, String> {
        if point.infinity {
            return Ok(point);
        }
        let p = c.p;
        let mut ds = Vec::with_capacity(self.kernel_x.len());
        let mut pref = Vec::with_capacity(ds.capacity());
        let mut product = 1;
        for &k in &self.kernel_x {
            let d = sub(point.x, k, p);
            if d == 0 {
                return Err("CM map kernel meets input point".into());
            }
            pref.push(product);
            product = counted_mul(product, d, p, f);
            ds.push(d);
        }
        f.invs += 1;
        let mut inverse = counted_pow(product, p - 2, p, f);
        for i in (0..ds.len()).rev() {
            let d = ds[i];
            ds[i] = counted_mul(inverse, pref[i], p, f);
            inverse = counted_mul(inverse, d, p, f);
        }
        let mut x = point.x;
        let mut derivative = 1;
        for (i, &d) in ds.iter().enumerate() {
            f.sqrs += 1;
            let d2 = mul(d, d, p);
            let d3 = counted_mul(d2, d, p, f);
            let t2 = add(self.kernel_t[i], self.kernel_t[i], p);
            let y4 = mul(4, self.kernel_rhs[i], p);
            let y8 = add(y4, y4, p);
            f.muls += 1;
            x = add(
                x,
                add(counted_mul(t2, d, p, f), counted_mul(y4, d2, p, f), p),
                p,
            );
            derivative = sub(
                derivative,
                add(counted_mul(t2, d2, p, f), counted_mul(y8, d3, p, f), p),
                p,
            );
        }
        let out = PrimePoint::affine(
            counted_mul(x, self.inv_c, p, f),
            counted_mul(counted_mul(point.y, derivative, p, f), self.inv_cu, p, f),
        );
        Ok(out)
    }
}

/// Recover the public scalar action from pi=(t+f sqrt(-d))/2 and pi=1 on E[r].
pub fn action(inst: &PrimeInstance, abs_d: u64) -> Result<(u64, u64, u64), String> {
    if !matches!(abs_d, 19 | 43) {
        return Err("CM orbit currently supports fundamental discriminants -19 and -43".into());
    }
    let p = inst.curve.p;
    let t = (p as i128 + 1) - (inst.group_order as i128);
    let disc = 4i128 * p as i128 - t * t;
    if disc <= 0 || disc % abs_d as i128 != 0 {
        return Err("Frobenius discriminant does not match the declared CM field".into());
    }
    let square = (disc / abs_d as i128) as u128;
    let ell = isqrt(square);
    if (ell as u128).pow(2) != square || ell == 0 || ell % inst.r == 0 {
        return Err("CM Frobenius conductor is not an invertible integer".into());
    }
    let residue = (2 - t).rem_euclid(inst.r as i128) as u64;
    let lambda = mul(
        add(
            1,
            mul(residue, pow(ell, inst.r - 2, inst.r), inst.r),
            inst.r,
        ),
        (inst.r + 1) / 2,
        inst.r,
    );
    let q = (abs_d + 1) / 4;
    if add(
        sub(mul(lambda, lambda, inst.r), lambda, inst.r),
        q % inst.r,
        inst.r,
    ) != 0
    {
        return Err("CM scalar action fails its characteristic polynomial".into());
    }
    Ok((ell, lambda, q))
}

/// Build a degree-q self-isogeny from rational q-torsion on the quadratic twist.
pub fn derive_map(inst: &PrimeInstance, abs_d: u64, ops: &mut GroupOps) -> Result<CmMap, String> {
    let (_, lambda, q) = action(inst, abs_d)?;
    let c = &inst.curve;
    let p = c.p;
    let jneg = match abs_d {
        19 => 884736,
        43 => 884736000,
        _ => unreachable!(),
    };
    let den = add(
        mul(4, mul(mul(c.a, c.a, p), c.a, p), p),
        mul(27, mul(c.b, c.b, p), p),
        p,
    );
    if den == 0
        || mul(
            mul(1728, mul(4, mul(mul(c.a, c.a, p), c.a, p), p), p),
            pow(den, p - 2, p),
            p,
        ) != sub(0, jneg % p, p)
    {
        return Err("Velu CM orbit requires the maximal-order CM j-invariant".into());
    }
    let t = p as i128 + 1 - inst.group_order as i128;
    let twist_order = (p as i128 + 1 + t) as u64;
    if !twist_order.is_multiple_of(q) {
        return Err("this trace has no required rational twist kernel".into());
    }
    let d = (2..1000)
        .find(|&d| c.legendre(d) == -1)
        .ok_or("no twist multiplier found")?;
    let twist = PrimeCurve {
        p,
        a: mul(c.a, mul(d, d, p), p),
        b: mul(c.b, mul(mul(d, d, p), d, p), p),
    };
    let mut kernel = None;
    for x in 0..10000 {
        if let Some(y) = twist.sqrt(twist.rhs(x)) {
            let kernel_point = twist.mul(ops, PrimePoint::affine(x, y), twist_order / q);
            if !kernel_point.infinity && twist.mul(ops, kernel_point, q).infinity {
                kernel = Some(kernel_point);
                break;
            }
        }
    }
    let kernel_point = kernel.ok_or("rational twist kernel search exhausted")?;
    let inv_d = pow(d, p - 2, p);
    let mut xs = Vec::new();
    let mut ts = Vec::new();
    let mut ys = Vec::new();
    let mut sum_t = 0;
    let mut sum_w = 0;
    for i in 1..=(q - 1) / 2 {
        let image = twist.mul(ops, kernel_point, i);
        let x = mul(image.x, inv_d, p);
        let x2 = mul(x, x, p);
        let t_term = add(mul(3, x2, p), c.a, p);
        let w_term = add(
            add(mul(5, mul(x2, x, p), p), mul(3, mul(c.a, x, p), p), p),
            mul(2, c.b, p),
            p,
        );
        xs.push(x);
        ts.push(t_term);
        ys.push(c.rhs(x));
        sum_t = add(sum_t, t_term, p);
        sum_w = add(sum_w, w_term, p);
    }
    let a1 = sub(c.a, mul(10, sum_t, p), p);
    let b1 = sub(c.b, mul(14, sum_w, p), p);
    if a1 == 0 || b1 == 0 || c.a == 0 || c.b == 0 {
        return Err("exceptional CM model not supported".into());
    }
    let scale = mul(mul(b1, c.a, p), pow(mul(a1, c.b, p), p - 2, p), p);
    if mul(c.a, mul(scale, scale, p), p) != a1
        || mul(c.b, mul(mul(scale, scale, p), scale, p), p) != b1
    {
        return Err("CM Velu codomain is not isomorphic to source".into());
    }
    let u = c.sqrt(scale).ok_or("CM isomorphism is not rational")?;
    let mut map = CmMap {
        lambda,
        degree: q,
        kernel_x: xs,
        kernel_t: ts,
        kernel_rhs: ys,
        inv_c: pow(scale, p - 2, p),
        inv_cu: pow(mul(scale, u, p), p - 2, p),
    };
    let expected = c.mul(ops, inst.generator_point(), lambda);
    let value = map.apply(c, inst.generator_point(), &mut ops.field_ops)?;
    if value == c.neg(expected) {
        map.inv_cu = sub(0, map.inv_cu, p);
    } else if value != expected {
        return Err("actual CM map action disagrees with the public eigenvalue".into());
    }
    Ok(map)
}

pub struct CmOrbitBase<'a> {
    pub instance: &'a PrimeInstance,
}
impl FactorBaseBuilder<PrimeCurve> for CmOrbitBase<'_> {
    fn name(&self) -> &str {
        "cm-scalar-orbit"
    }
    fn describe(&self, p: &Params) -> String {
        format!(
            "{} unknown seeds, {} signed scalar-orbit levels, mode {}",
            p.get("columns").unwrap_or("?"),
            p.get("length").unwrap_or("?"),
            p.get("mode").unwrap_or("?")
        )
    }
    fn parameters(&self) -> &[(&str, &str)] {
        &[
            ("columns", "unknown seed logs"),
            ("length", "distinct abscissae per seed"),
            ("abs_d", "19 or 43"),
            ("mode", "velu or scalar"),
        ]
    }
    fn build(
        &self,
        _ctx: &InstanceCtx<PrimeCurve>,
        params: &Params,
        ops: &mut GroupOps,
    ) -> Result<FactorBase<PrimePoint>, String> {
        if params
            .0
            .keys()
            .any(|k| !matches!(k.as_str(), "columns" | "length" | "abs_d" | "mode"))
        {
            return Err("unknown CM orbit parameter".into());
        }
        let columns = params.u64("columns")? as usize;
        let length = params.u64("length")? as usize;
        let abs_d = params.u64("abs_d")?;
        let mode = params.get("mode").ok_or("CM orbit requires mode")?;
        if columns == 0 || length == 0 || columns.checked_mul(length).is_none_or(|n| n > 16384) {
            return Err("CM orbit requires 1..16384 abscissae".into());
        }
        if mode != "scalar" && mode != "velu" {
            return Err("CM orbit mode must be scalar or velu".into());
        }
        let started = Instant::now();
        let c = &self.instance.curve;
        let p = c.p;
        let r = self.instance.r;
        let (_, lambda, degree) = action(self.instance, abs_d)?;
        let map = if mode == "velu" {
            Some(derive_map(self.instance, abs_d, ops)?)
        } else {
            None
        };
        let mut fb = FactorBase {
            points: Vec::new(),
            neg_index: Vec::new(),
            col_of: Vec::new(),
            coef_of: Vec::new(),
            columns,
            abscissae: 0,
            dimension: None,
            description: self.describe(params),
            abscissa_list: Vec::new(),
            x_index: fast_map(2 * columns * length),
            point_index: fast_map(4 * columns * length),
            cost: PhaseCost::default(),
            subspace_basis: None,
        };
        let mut seen = HashSet::new();
        let mut x = 1;
        let mut extra = FieldOps::default();
        for column in 0..columns {
            let mut point = loop {
                let cand = x;
                x += 1;
                fb.cost.count("abscissae_scanned", 1);
                if cand == self.instance.generator.0 {
                    continue;
                }
                fb.cost.count("sqrt_solves", 1);
                if let Some(y) = c.sqrt(c.rhs(cand)) {
                    if y != 0 {
                        break PrimePoint::affine(cand, y.min(p - y));
                    }
                }
            };
            let mut coef = 1;
            for level in 0..length {
                if !seen.insert(point.x) {
                    return Err("CM orbit repeats an abscissa across seeds or levels".into());
                }
                if !c.is_on_curve(point) {
                    return Err("CM orbit image is not on the curve".into());
                }
                let i = fb.points.len();
                fb.points.extend([point, c.neg(point)]);
                fb.col_of.extend([column, column]);
                fb.coef_of.extend([coef, r - coef]);
                fb.neg_index.extend([i + 1, i]);
                fb.x_index.insert(point.x, vec![i, i + 1]);
                fb.point_index.insert(c.key(&point), i);
                fb.point_index.insert(c.key(&c.neg(point)), i + 1);
                fb.abscissa_list.push(point.x);
                if level + 1 < length {
                    point = if let Some(m) = &map {
                        m.apply(c, point, &mut extra)?
                    } else {
                        c.mul(ops, point, lambda)
                    };
                    coef = mul(coef, lambda, r);
                    if point.infinity {
                        return Err("CM orbit reached infinity".into());
                    }
                    if point.y > p - point.y {
                        point = c.neg(point);
                        coef = r - coef;
                    }
                    fb.cost.count("orbit_steps", 1);
                }
            }
        }
        fb.abscissa_list.sort_unstable();
        fb.abscissae = columns * length;
        fb.cost.count("cm_degree", degree);
        fb.cost.count("cm_formula_muls_uncharged", extra.muls);
        fb.cost.count("cm_formula_sqrs_uncharged", extra.sqrs);
        fb.cost.count("cm_formula_invs_uncharged", extra.invs);
        fb.cost.wall_ns = started.elapsed().as_nanos() as u64;
        Ok(fb)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn degree_eleven_second_family_map_matches_scalar_action() {
        let inst=PrimeInstance{name:"D=-43 independent PARI control".into(),
            curve:PrimeCurve{p:5325000452995313,a:130015040327813,b:260030080655626},
            group_order:5325000430993249,r:5325000430993249,cofactor:1,
            generator:(4,437435961147748)};
        let mut ops=GroupOps::default();
        let map=derive_map(&inst,43,&mut ops).unwrap();
        assert_eq!(map.degree,11);
        for k in [1,2,7,19,12345] {
            let p=inst.curve.mul(&mut ops,inst.generator_point(),k);
            assert_eq!(map.apply(&inst.curve,p,&mut FieldOps::default()).unwrap(),
                inst.curve.mul(&mut ops,p,map.lambda));
        }
    }
    fn instance() -> PrimeInstance {
        PrimeInstance {
            name: "public control".into(),
            curve: PrimeCurve {
                p: 343073,
                a: 36086,
                b: 72172,
            },
            group_order: 342833,
            r: 342833,
            cofactor: 1,
            generator: (0, 2882),
        }
    }
    #[test]
    fn velu_and_scalar_orbits_are_exactly_identical() {
        let i = instance();
        let ctx = InstanceCtx {
            group: &i.curve,
            generator: i.generator_point(),
            target: PrimePoint::INFINITY,
            r: i.r,
            cofactor: 1,
            group_order: i.r,
            name: i.name.clone(),
            field_degree: None,
        };
        let mut p = Params::default();
        p.set("columns", "4")
            .set("length", "16")
            .set("abs_d", "19")
            .set("mode", "velu");
        let a = CmOrbitBase { instance: &i }
            .build(&ctx, &p, &mut GroupOps::default())
            .unwrap();
        p.set("mode", "scalar");
        let b = CmOrbitBase { instance: &i }
            .build(&ctx, &p, &mut GroupOps::default())
            .unwrap();
        assert_eq!(a.points, b.points);
        assert_eq!(a.col_of, b.col_of);
        assert_eq!(a.coef_of, b.coef_of);
        assert_eq!(a.points.len(), 128);
        assert_eq!(a.columns, 4);
    }
    #[test]
    fn map_action_and_polynomial_on_public_generator() {
        let i = instance();
        let mut ops = GroupOps::default();
        let m = derive_map(&i, 19, &mut ops).unwrap();
        assert_eq!(m.degree, 5);
        assert_eq!(
            m.apply(&i.curve, i.generator_point(), &mut FieldOps::default())
                .unwrap(),
            i.curve.mul(&mut ops, i.generator_point(), m.lambda)
        );
        assert_eq!(
            add(sub(mul(m.lambda, m.lambda, i.r), m.lambda, i.r), 5, i.r),
            0
        );
    }
    #[test]
    fn floor_refuses_the_low_degree_map() {
        let mut i = instance();
        i.curve.a = 43182;
        i.curve.b = 57354;
        i.generator = (2, 20047);
        assert!(derive_map(&i, 19, &mut GroupOps::default()).is_err());
        assert!(action(&i, 19).is_ok());
    }
    #[test]
    fn malformed_parameters_fail_closed() {
        let i = instance();
        assert!(action(&i, 11).is_err());
        let mut wrong = i.clone();
        wrong.group_order += 2;
        assert!(action(&wrong, 19).is_err());
    }
}

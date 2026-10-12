//! Exact GLV lattice splitting and complete scalar-arithmetic diagnostics.
use crate::{
    curves::{Counts, Curve, Map, Point},
    Result,
};
use serde_json::{json, Value};

pub type Vector = (i64, i64);
pub type Lattice = (Vector, Vector);
fn dot(u: Vector, v: Vector) -> i64 {
    u.0 * v.0 + u.1 * v.1
}
fn nearest(mut a: i64, mut b: i64) -> i64 {
    assert_ne!(b, 0);
    if b < 0 {
        a = -a;
        b = -b;
    }
    a.signum() * (a.abs() / b + i64::from(2 * (a.abs() % b) >= b))
}
pub fn lattice(n: i64, value: i64, cost: &mut Counts) -> Result<Lattice> {
    // This backend is intentionally tied to the bounded laboratory subgroups.
    if !(2..=32762).contains(&n) {
        return Err("lattice subgroup outside laboratory bound".into());
    }
    let v = value.rem_euclid(n);
    let (mut u, mut w) = ((n, 0), ((-v).rem_euclid(n), 1));
    let mut done = false;
    for _ in 0..1000 {
        cost.integer_steps += 1;
        if dot(w, w) < dot(u, u) {
            std::mem::swap(&mut u, &mut w);
        }
        if 2 * dot(u, w).abs() <= dot(u, u) {
            done = true;
            break;
        }
        let q = nearest(dot(u, w), dot(u, u));
        w = (w.0 - q * u.0, w.1 - q * u.1);
    }
    let det = u.0 * w.1 - u.1 * w.0;
    if !done
        || det.abs() != n
        || (u.0 + v * u.1).rem_euclid(n) != 0
        || (w.0 + v * w.1).rem_euclid(n) != 0
    {
        return Err("invalid Gauss reduction".into());
    }
    Ok((u, w))
}
pub fn split(k: i64, n: i64, value: i64, basis: Lattice, cost: &mut Counts) -> Result<Vector> {
    if !(2..=32762).contains(&n)
        || [basis.0 .0, basis.0 .1, basis.1 .0, basis.1 .1]
            .iter()
            .any(|v| v.unsigned_abs() > n as u64)
    {
        return Err("invalid bounded decomposition lattice".into());
    }
    let k = k.rem_euclid(n);
    let (u, v) = basis;
    let det = u.0 * v.1 - u.1 * v.0;
    let lambda = value.rem_euclid(n);
    if det.abs() != n
        || (u.0 + lambda * u.1).rem_euclid(n) != 0
        || (v.0 + lambda * v.1).rem_euclid(n) != 0
    {
        return Err("invalid decomposition lattice congruences".into());
    }
    let a = nearest(k * v.1, det);
    let b = nearest(-k * u.1, det);
    let mut choices = Vec::new();
    for dx in -1..=1 {
        for dy in -1..=1 {
            cost.integer_steps += 1;
            choices.push((
                k - (a + dx) * u.0 - (b + dy) * v.0,
                -(a + dx) * u.1 - (b + dy) * v.1,
            ));
        }
    }
    choices.sort_by_key(|(x, y)| (x.abs().max(y.abs()), x.abs() + y.abs(), *x, *y));
    let (x, y) = choices[0];
    if (x + value.rem_euclid(n) * y - k).rem_euclid(n) != 0 {
        return Err("scalar decomposition failed".into());
    }
    Ok((x, y))
}
pub fn joint(c: Curve, a: i64, p: Point, b: i64, q: Point, cost: &mut Counts) -> Point {
    let p = if a < 0 { c.neg_count(p, cost) } else { p };
    let q = if b < 0 { c.neg_count(q, cost) } else { q };
    let a = a.unsigned_abs();
    let b = b.unsigned_abs();
    let table = [None, q, p, c.add_count(p, q, cost)];
    let mut r = None;
    for i in (0..64 - (a | b).leading_zeros()).rev() {
        r = c.add_count(r, r, cost);
        let d = (2 * ((a >> i) & 1) + ((b >> i) & 1)) as usize;
        if d != 0 {
            r = c.add_count(r, table[d], cost);
        }
    }
    r
}
pub fn accelerated(
    c: Curve,
    k: i64,
    p: Point,
    map: &Map,
    n: i64,
    basis: Lattice,
    cost: &mut Counts,
) -> Result<Point> {
    let (a, b) = split(
        k,
        n,
        map.eigenvalue.ok_or("missing eigenvalue")?,
        basis,
        cost,
    )?;
    let image = map.evaluate(c, p, cost);
    Ok(joint(c, a, p, b, image, cost))
}
pub fn diagnostic(c: Curve, maps: &[Map], n: i64, g: Point) -> Result<Value> {
    let expected: Vec<_> = (0..n).map(|k| c.mul(k, g)).collect();
    let mut baseline = Vec::new();
    for width in [0, 3, 4] {
        let mut cost = Counts::default();
        for k in 0..n {
            let image = if width == 0 {
                c.mul_count(k, g, &mut cost)
            } else {
                c.wnaf(k, g, width, &mut cost)
            };
            if image != expected[k as usize] {
                return Err(format!("baseline failed for k={k}"));
            }
        }
        baseline.push(
            json!({"method":if width==0{"binary"}else if width==3{"wnaf3"}else{"wnaf4"},
            "counts":cost.record(),"exhaustive_scalar_checks":n,"table_precomputation":"per_call"}),
        );
    }
    let mut results = Vec::new();
    for map in maps {
        let v = map.eigenvalue.ok_or("unverified map")?;
        if [0, 1, n - 1].contains(&v) {
            results
                .push(json!({"map_id":map.id,"status":"trivial_or_annihilating_subgroup_action"}));
            continue;
        }
        let mut setup = Counts::default();
        let basis = lattice(n, v, &mut setup)?;
        let mut cost = Counts::default();
        let mut maximum = (0, 0);
        for k in 0..n {
            let (a, b) = split(k, n, v, basis, &mut Counts::default())?;
            maximum.0 = maximum.0.max(a.abs());
            maximum.1 = maximum.1.max(b.abs());
            if accelerated(c, k, g, map, n, basis, &mut cost)? != expected[k as usize] {
                return Err(format!("joint multiplication failed for {}, k={k}", map.id));
            }
        }
        results.push(json!({"map_id":map.id,"status":"verified_scalar_stage_diagnostic","lattice":basis,
            "exhaustive_scalar_checks":n,"maximum_absolute_components":maximum,"online_counts":cost.record(),
            "one_batch_lattice_setup_counts":setup.record(),"amortization_operations":n,
            "wall_time":null,"speedup":null,"end_to_end_cost":null}));
    }
    Ok(
        json!({"kind":"scalar_arithmetic_stage_diagnostic","workload":{"scalars":"all integers 0 <= k < subgroup_order",
        "scalar_count":n,"base_points":1,"precomputation":"per_call","lattice_setup":"once per full batch"},
        "baseline_candidates":baseline,"maps":results,"operation_units":"separate uncalibrated counters; no weighted total",
        "integer_steps_definition":"Gauss-loop iterations, split neighbor checks, and wNAF digit iterations",
        "included_arithmetic":["map evaluation","per-call joint/wNAF tables","complete group arithmetic"],
        "excluded_setup":["map discovery","map verification","independent correctness replay"],
        "native_wall_timing":"unmeasured","end_to_end_cost":null,"speedup":null,"constant_time_claim":false}),
    )
}

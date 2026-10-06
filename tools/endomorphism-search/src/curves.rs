//! Small prime fields, independent affine group arithmetic and Vélú maps.
use crate::{digest, orders, scalar, Result, SCHEMA};
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

pub type Point = Option<(i64, i64)>;
pub const MAX_PRIME: i64 = 16381;

#[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
pub struct Counts {
    pub field_add: u64,
    pub field_mul: u64,
    pub field_inverse: u64,
    pub group_calls: u64,
    pub integer_steps: u64,
}
impl Counts {
    pub fn record(self) -> Value {
        json!({"field_additions":self.field_add,"field_multiplications":self.field_mul,
        "field_inversions":self.field_inverse,"group_calls":self.group_calls,"integer_steps":self.integer_steps})
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Curve {
    pub p: i64,
    pub a: i64,
    pub b: i64,
}
impl Curve {
    pub fn new(p: i64, a: i64, b: i64) -> Result<Self> {
        if p <= 3 || p > MAX_PRIME || orders::factor(p)? != vec![(p, 1)] {
            return Err(format!(
                "laboratory curves require prime 3 < p <= {MAX_PRIME}"
            ));
        }
        let c = Self {
            p,
            a: a.rem_euclid(p),
            b: b.rem_euclid(p),
        };
        if (4 * c.a * c.a * c.a + 27 * c.b * c.b) % p == 0 {
            return Err("singular curve".into());
        }
        Ok(c)
    }
    fn fa(self, a: i64, b: i64, cost: &mut Counts) -> i64 {
        cost.field_add += 1;
        (a + b).rem_euclid(self.p)
    }
    fn fs(self, a: i64, b: i64, cost: &mut Counts) -> i64 {
        cost.field_add += 1;
        (a - b).rem_euclid(self.p)
    }
    fn fm(self, a: i64, b: i64, cost: &mut Counts) -> i64 {
        cost.field_mul += 1;
        (a * b).rem_euclid(self.p)
    }
    fn inv(self, a: i64, cost: &mut Counts) -> i64 {
        cost.field_inverse += 1;
        let (mut r0, mut r1, mut s0, mut s1) = (self.p, a.rem_euclid(self.p), 0, 1);
        assert_ne!(r1, 0, "inverse of zero");
        while r1 != 0 {
            let q = r0 / r1;
            (r0, r1) = (r1, r0 - q * r1);
            (s0, s1) = (s1, s0 - q * s1);
        }
        assert_eq!(r0, 1);
        s0.rem_euclid(self.p)
    }
    pub fn pow(self, mut a: i64, mut k: i64) -> i64 {
        let mut v = 1;
        a = a.rem_euclid(self.p);
        while k > 0 {
            if k % 2 != 0 {
                v = v * a % self.p;
            }
            a = a * a % self.p;
            k /= 2;
        }
        v
    }
    pub fn contains(self, v: Point) -> bool {
        v.is_none_or(|(x, y)| {
            (0..self.p).contains(&x)
                && (0..self.p).contains(&y)
                && (y * y - x * x * x - self.a * x - self.b).rem_euclid(self.p) == 0
        })
    }
    pub fn neg(self, p: Point) -> Point {
        p.map(|(x, y)| (x, (-y).rem_euclid(self.p)))
    }
    pub fn neg_count(self, p: Point, cost: &mut Counts) -> Point {
        p.map(|(x, y)| (x, self.fs(0, y, cost)))
    }
    pub fn add(self, p: Point, q: Point) -> Point {
        self.add_count(p, q, &mut Counts::default())
    }
    pub fn add_count(self, p: Point, q: Point, cost: &mut Counts) -> Point {
        cost.group_calls += 1;
        let (Some((x1, y1)), Some((x2, y2))) = (p, q) else {
            return p.or(q);
        };
        let slope = if x1 == x2 {
            if self.fa(y1, y2, cost) == 0 {
                return None;
            }
            let sq = self.fm(x1, x1, cost);
            let num = self.fm(3, sq, cost);
            let num = self.fa(num, self.a, cost);
            let den = self.fm(2, y1, cost);
            let inv = self.inv(den, cost);
            self.fm(num, inv, cost)
        } else {
            let num = self.fs(y2, y1, cost);
            let den = self.fs(x2, x1, cost);
            let inv = self.inv(den, cost);
            self.fm(num, inv, cost)
        };
        let x3 = self.fm(slope, slope, cost);
        let x3 = self.fs(x3, x1, cost);
        let x3 = self.fs(x3, x2, cost);
        let y3 = self.fs(x1, x3, cost);
        let y3 = self.fm(slope, y3, cost);
        let y3 = self.fs(y3, y1, cost);
        Some((x3, y3))
    }
    pub fn mul(self, k: i64, p: Point) -> Point {
        self.mul_count(k, p, &mut Counts::default())
    }
    pub fn mul_count(self, k: i64, p: Point, cost: &mut Counts) -> Point {
        let p = if k < 0 { self.neg_count(p, cost) } else { p };
        let k = k.unsigned_abs();
        let mut r = None;
        for i in (0..64 - k.leading_zeros()).rev() {
            r = self.add_count(r, r, cost);
            if k & (1u64 << i) != 0 {
                r = self.add_count(r, p, cost);
            }
        }
        r
    }
    pub fn wnaf(self, k: i64, p: Point, width: u32, cost: &mut Counts) -> Point {
        assert!((3..=4).contains(&width));
        let p = if k < 0 { self.neg_count(p, cost) } else { p };
        let mut k = k.unsigned_abs() as i128;
        if k == 0 || p.is_none() {
            return None;
        }
        let twice = self.add_count(p, p, cost);
        let mut table = vec![p];
        for i in 1..(1 << (width - 2)) {
            table.push(self.add_count(table[i - 1], twice, cost));
        }
        let mut digits = Vec::new();
        while k > 0 {
            cost.integer_steps += 1;
            let mut d = 0;
            if k % 2 != 0 {
                d = k % (1 << width);
                if d >= 1 << (width - 1) {
                    d -= 1 << width;
                }
                k -= d;
            }
            digits.push(d);
            k /= 2;
        }
        let mut r = None;
        for d in digits.into_iter().rev() {
            r = self.add_count(r, r, cost);
            if d != 0 {
                let q = table[(d.unsigned_abs() as usize - 1) / 2];
                let q = if d < 0 { self.neg_count(q, cost) } else { q };
                r = self.add_count(r, q, cost);
            }
        }
        r
    }
    pub fn points(self) -> Vec<Point> {
        let mut squares = vec![Vec::new(); self.p as usize];
        for y in 0..self.p {
            squares[(y * y % self.p) as usize].push(y);
        }
        let mut out = vec![None];
        for x in 0..self.p {
            for y in &squares[((x * x * x + self.a * x + self.b) % self.p) as usize] {
                out.push(Some((x, *y)));
            }
        }
        out
    }
    pub fn model(self) -> Value {
        json!({"p":self.p,"a":self.a,"b":self.b,"model":"short_weierstrass"})
    }
    pub fn identity(self, count: i64, ring: Option<i64>, prime: i64, g: Point) -> Value {
        let model = json!({"a":self.a.to_string(),"b":self.b.to_string(),"field":format!("fp-{}",self.p),
            "form":"y^2=x^3+a*x+b","p":self.p.to_string(),"v":"1"});
        let model_json = model.to_string();
        let hash = digest(model_json.as_bytes());
        let t = self.p + 1 - count;
        let tag = if t < 0 {
            format!("tm{}", -t)
        } else {
            format!("t{t}")
        };
        let disc = (4 * self.a * self.a * self.a + 27 * self.b * self.b) % self.p;
        let j = 1728 * 4 * self.a * self.a * self.a % self.p * self.pow(disc, self.p - 2) % self.p;
        let slug = format!(
            "icv1-fp{}-{tag}-{}",
            64 - self.p.leading_zeros(),
            &hash[..8]
        );
        let field = json!({"characteristic":self.p,"degree":1,"representation":"prime",
            "modulus":self.p,"element_encoding":"canonical decimal integers"});
        let curve = json!({"model":"short_weierstrass","a":self.a,"b":self.b,"group_order":count,
            "subgroup_order":prime,"cofactor":count/prime,"generator":g,"target_group":"prime_order_subgroup"});
        let uid = digest(json!({"field":field,"curve":curve}).to_string().as_bytes());
        json!({"name":slug,"icv1":format!("ICV1:fp-{}:{t}:{count}:{j}:{}:unk:r:{}",self.p,
            ring.map_or_else(||"unk".into(),|d|d.to_string()),&hash[..12]),"model_json":model_json,
            "ec1":format!("EC1P{}Cswh{}",64-self.p.leading_zeros(),&uid[..12]),
            "curve_uid":format!("urn:ec-record:1:sha256:{uid}"),"field_sha256":digest(field.to_string().as_bytes()),
            "field":field,"curve":curve,"registration_status":"not_registered_by_this_tool"})
    }
}
pub fn subgroup(c: Curve, points: &[Point]) -> Result<(i64, Point)> {
    let total = points.len() as i64;
    let factors = orders::factor(total)?;
    for (n, _) in factors.iter().rev() {
        if *n < 5 {
            continue;
        }
        for p in &points[1..] {
            let mut order = total;
            for (q, _) in &factors {
                while order % q == 0 && c.mul(order / q, *p).is_none() {
                    order /= q;
                }
            }
            if order % n == 0 {
                let g = c.mul(order / n, *p);
                if g.is_some() && c.mul(*n, g).is_none() {
                    return Ok((*n, g));
                }
            }
        }
    }
    Err("no laboratory prime subgroup of order >= 5".into())
}

#[derive(Clone, Debug)]
pub enum Operation {
    Scale {
        scale: i64,
        x: i64,
        y: i64,
        order: i64,
    },
    Velu {
        kernel: Vec<Point>,
        target: Curve,
        scale: i64,
    },
    Frobenius {
        c: i64,
    },
}
#[derive(Clone, Debug)]
pub struct Map {
    pub id: String,
    pub degree: i64,
    pub op: Operation,
    pub traces: Vec<i64>,
    pub eigenvalue: Option<i64>,
    pub verification: Value,
}
impl Map {
    pub fn family(&self) -> &'static str {
        match self.op {
            Operation::Scale { .. } => "scaling_automorphism",
            Operation::Velu { .. } => "rational_velu_self_isogeny",
            Operation::Frobenius { .. } => "frobenius_minus_scalar",
        }
    }
    pub fn evaluate(&self, c: Curve, p: Point, cost: &mut Counts) -> Point {
        match &self.op {
            Operation::Scale { x, y, .. } => {
                p.map(|(px, py)| (c.fm(*x, px, cost), c.fm(*y, py, cost)))
            }
            Operation::Frobenius { c: k } => {
                let q = c.mul_count(*k, p, cost);
                let q = c.neg_count(q, cost);
                c.add_count(p, q, cost)
            }
            Operation::Velu { kernel, scale, .. } => {
                if kernel.contains(&p) {
                    return None;
                }
                let (mut x, mut y) = p.unwrap();
                for q in kernel.iter().flatten() {
                    let image = c.add_count(p, Some(*q), cost).expect("pole only in kernel");
                    let dx = c.fs(image.0, q.0, cost);
                    let dy = c.fs(image.1, q.1, cost);
                    x = c.fa(x, dx, cost);
                    y = c.fa(y, dy, cost);
                }
                // These constants are computed here, so charge both products.
                let xs = c.fm(*scale, *scale, cost);
                let ys = c.fm(xs, *scale, cost);
                Some((c.fm(xs, x, cost), c.fm(ys, y, cost)))
            }
        }
    }
    pub fn record(&self) -> Value {
        let construction = match &self.op {
            Operation::Scale { scale, x, y, order } => {
                json!({"scale":scale,"x_scale":x,"y_scale":y,"unit_order":order})
            }
            Operation::Velu {
                kernel,
                target,
                scale,
            } => {
                json!({"kernel":kernel,"intermediate_curve":target.model(),"return_isomorphism_scale":scale})
            }
            Operation::Frobenius { c } => {
                json!({"formula":"omega=pi-[c]","c":c,"rational_restriction":"[1-c]"})
            }
        };
        json!({"map_id":self.id,"family":self.family(),"geometric_degree":self.degree,"construction":construction,
            "trace_candidates":self.traces,"subgroup_eigenvalue":self.eigenvalue,"verification":self.verification})
    }
}
pub fn scales(c: Curve) -> Vec<Map> {
    let mut maps = Vec::new();
    for s in 2..c.p - 1 {
        if (c.pow(s, 4) - 1) * c.a % c.p != 0 || (c.pow(s, 6) - 1) * c.b % c.p != 0 {
            continue;
        }
        let (mut v, mut n) = (s, 1);
        while v != 1 {
            v = v * s % c.p;
            n += 1;
        }
        let t = match n {
            3 => -1,
            4 => 0,
            6 => 1,
            _ => unreachable!(),
        };
        maps.push(Map {
            id: format!("scale:{s}"),
            degree: 1,
            op: Operation::Scale {
                scale: s,
                x: c.pow(s, 2),
                y: c.pow(s, 3),
                order: n,
            },
            traces: vec![t],
            eigenvalue: None,
            verification: Value::Null,
        });
    }
    maps
}
pub fn velu_target(c: Curve, kernel: &[Point]) -> Result<Curve> {
    let set: BTreeSet<_> = kernel.iter().copied().collect();
    if kernel.len() < 2
        || set.len() != kernel.len()
        || !set.contains(&None)
        || kernel.iter().any(|p| !c.contains(*p))
        || kernel
            .iter()
            .any(|p| kernel.iter().any(|q| !set.contains(&c.add(*p, *q))))
    {
        return Err("kernel must be a distinct closed subgroup containing identity".into());
    }
    let mut t = 0;
    let mut w = 0;
    for (x, _) in kernel.iter().flatten() {
        t = (t + 3 * x * x + c.a) % c.p;
        w = (w + 5 * x * x * x + 3 * c.a * x + 2 * c.b) % c.p;
    }
    Curve::new(c.p, c.a - 5 * t, c.b - 7 * w)
}
pub fn rational_maps(
    c: Curve,
    points: &[Point],
    bound: i64,
    max_kernels: usize,
    seconds: f64,
) -> Result<(Vec<Map>, Value)> {
    if !(1..=1024).contains(&max_kernels) || !seconds.is_finite() || seconds <= 0.0 {
        return Err("finite positive kernel limits required; max-kernels <= 1024".into());
    }
    let start = Instant::now();
    let mut attempted = 0;
    let mut stop = None;
    let mut maps = Vec::new();
    let mut seen = BTreeSet::new();
    let degrees: Vec<_> = [2, 3, 5, 7, 11, 13]
        .into_iter()
        .filter(|d| *d <= bound)
        .collect();
    'outer: for d in &degrees {
        if points.len() as i64 % d != 0 {
            continue;
        }
        for p in &points[1..] {
            if start.elapsed().as_secs_f64() >= seconds {
                stop = Some("time_limit");
                break 'outer;
            }
            if c.mul(*d, *p).is_some() {
                continue;
            }
            let kernel: Vec<_> = (0..*d).map(|i| c.mul(i, *p)).collect();
            let mut key = kernel.clone();
            key.sort();
            if seen.contains(&key) {
                continue;
            }
            if attempted >= max_kernels {
                stop = Some("kernel_limit");
                break 'outer;
            }
            seen.insert(key);
            attempted += 1;
            let target = velu_target(c, &kernel)?;
            for s in 1..c.p {
                if start.elapsed().as_secs_f64() >= seconds {
                    stop = Some("time_limit");
                    break 'outer;
                }
                if c.pow(s, 4) * target.a % c.p != c.a || c.pow(s, 6) * target.b % c.p != c.b {
                    continue;
                }
                let radius = orders::isqrt((4 * d) as u64) as i64;
                maps.push(Map {
                    id: format!("velu:{d}:{attempted}:{s}"),
                    degree: *d,
                    op: Operation::Velu {
                        kernel: kernel.clone(),
                        target,
                        scale: s,
                    },
                    traces: (-radius..=radius).collect(),
                    eigenvalue: None,
                    verification: Value::Null,
                });
            }
        }
    }
    Ok((
        maps,
        json!({"status":if stop.is_some(){"incomplete"}else{"complete_for_declared_family"},"stop_reason":stop,
        "prime_degrees":degrees,"rational_kernels_only":true,"path_depth":1,"kernels_examined":attempted,
        "unsearched_families":["non-rational kernels","longer isogeny paths","extension-field maps"],
        "wall_seconds":null,"verification_outside_kernel_time_limit":true}),
    ))
}
pub fn verify_map(c: Curve, points: &[Point], n: i64, g: Point, map: &mut Map) -> Result<()> {
    let mut cost = Counts::default();
    let image = map.evaluate(c, g, &mut cost);
    let mut matches = BTreeMap::<i64, BTreeSet<i64>>::new();
    for t in &map.traces {
        for v in 0..n {
            if (v * v - t * v + map.degree).rem_euclid(n) == 0 && c.mul(v, g) == image {
                matches.entry(v).or_default().insert(*t);
            }
        }
    }
    if matches.len() == 1 {
        let (v, ts) = matches.into_iter().next().unwrap();
        map.eigenvalue = Some(v);
        map.traces = ts.into_iter().collect();
    }
    let on_curve = points
        .iter()
        .all(|p| c.contains(map.evaluate(c, *p, &mut cost)));
    let identity = map.evaluate(c, None, &mut cost).is_none();
    let mut seed = 203u64;
    let mut pairs = true;
    for _ in 0..64 {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        let p = points[(seed % points.len() as u64) as usize];
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        let q = points[(seed % points.len() as u64) as usize];
        pairs &= map.evaluate(c, c.add(p, q), &mut cost)
            == c.add(map.evaluate(c, p, &mut cost), map.evaluate(c, q, &mut cost));
    }
    let subgroup = map.eigenvalue.is_some_and(|v| {
        (0..n).all(|k| map.evaluate(c, c.mul(k, g), &mut cost) == c.mul((k * v) % n, g))
    });
    map.verification = json!({"algebraic_construction_conditions_checked":true,"all_rational_images_on_curve":on_curve,
        "rational_points_checked":points.len(),"identity_preserved":identity,"sample_homomorphism_checks":64,
        "sample_homomorphism_checks_passed":pairs,"sample_seed":203,"sample_rng":"lcg64/v1",
        "subgroup_action_exhaustive":subgroup,"subgroup_elements_checked":n,
        "verification_scope":"construction_conditions_and_rational_laboratory_points; samples_are_not_a_geometric_proof"});
    if !(on_curve && identity && pairs && subgroup) {
        return Err(format!("map verification failed: {}", map.id));
    }
    Ok(())
}
pub fn symmetries(values: &[i64], n: i64) -> Value {
    let mut generators = BTreeSet::from([n - 1]);
    for v in values {
        if v.rem_euclid(n) != 0 {
            generators.insert(v.rem_euclid(n));
        }
    }
    let mut group = BTreeSet::from([1]);
    let mut frontier = vec![1];
    while let Some(v) = frontier.pop() {
        for g in &generators {
            let q = v * g % n;
            if group.insert(q) {
                frontier.push(q);
            }
        }
    }
    json!({"distinct_action_count_including_negation":group.len(),"baseline_action_count":2,
        "generators":generators,"measured_rho_speedup":null,"orbit_canonicalization_cost":null})
}
pub fn probe(
    c: Curve,
    bound: i64,
    max_kernels: usize,
    seconds: f64,
    count_ops: bool,
) -> Result<Value> {
    if !(1..=1_000_000_000_000).contains(&bound) {
        return Err("degree bound must be in 1..10^12".into());
    }
    let points = c.points();
    let total = points.len() as i64;
    let t = c.p + 1 - total;
    let ordinary = t % c.p != 0;
    let (n, g) = subgroup(c, &points)?;
    let mut maps = scales(c);
    let mut full_ring = None;
    let mut binding = if ordinary {
        let (dk, f) = orders::fundamental(t * t - 4 * c.p)?;
        let reason = if c.a == 0 {
            full_ring = Some(-3);
            "ordinary_j0"
        } else if c.b == 0 {
            full_ring = Some(-4);
            "ordinary_j1728"
        } else if f == 1 {
            full_ring = Some(dk);
            "fundamental_frobenius_discriminant_forces_maximal_order"
        } else {
            "conductor_unresolved"
        };
        let discriminants: Vec<_> = (1..=f).filter(|d| f % d == 0).map(|d| dk * d * d).collect();
        if f == 1 {
            let s = dk.rem_euclid(2);
            let k = (t - s) / 2;
            let degree = (s * s - dk) / 4;
            maps.push(Map {
                id: "frobenius-ring-generator".into(),
                degree,
                op: Operation::Frobenius { c: k },
                traces: vec![s],
                eigenvalue: None,
                verification: Value::Null,
            });
        }
        json!({"status":if full_ring.is_some(){"certified_by_exact_elementary_conditions"}else{"unresolved"},
            "frobenius_discriminant":t*t-4*c.p,"fundamental_discriminant":dk,"frobenius_conductor":f,
            "full_ring_discriminant":full_ring,"reason":reason,"candidate_discriminants_before_map_search":discriminants})
    } else {
        json!({"status":"supersingular","full_geometric_ring":"quaternion_order"})
    };
    let (isogenies, coverage) = rational_maps(c, &points, bound, max_kernels, seconds)?;
    maps.extend(isogenies);
    maps.retain(|m| m.degree <= bound);
    for m in &mut maps {
        verify_map(c, &points, n, g, m)?;
    }
    if ordinary && full_ring.is_none() {
        let mut compatible = Vec::new();
        let evidence: Vec<_> = maps
            .iter()
            .filter(|m| matches!(m.op, Operation::Velu { .. }))
            .collect();
        for d in binding["candidate_discriminants_before_map_search"]
            .as_array()
            .unwrap()
        {
            let d = d.as_i64().unwrap();
            let mut okay = true;
            for m in &evidence {
                let e = orders::enumerate(d, m.degree, true, orders::Limits::default())?;
                if e["status"] == "complete"
                    && !e["candidates"]
                        .as_array()
                        .unwrap()
                        .iter()
                        .any(|x| x["degree"] == m.degree)
                {
                    okay = false;
                    break;
                }
            }
            if okay {
                compatible.push(d);
            }
        }
        if compatible.is_empty() {
            return Err("constructed maps contradict all candidate ordinary orders".into());
        }
        binding["candidate_discriminants_after_map_search"] = json!(compatible);
        binding["degree_evidence"] = json!(evidence
            .iter()
            .map(|m| json!({"map_id":m.id,"degree":m.degree}))
            .collect::<Vec<_>>());
        if compatible.len() == 1 {
            full_ring = Some(compatible[0]);
            binding["status"] = json!("certified_by_exact_elementary_conditions");
            binding["full_ring_discriminant"] = json!(full_ring);
            binding["reason"] = json!("non_scalar_map_degrees_exclude_other_candidate_conductors");
        }
    }
    let values: Vec<_> = maps.iter().filter_map(|m| m.eigenvalue).collect();
    let sym = symmetries(&values, n);
    let mut report = json!({"schema":SCHEMA,"curve":c.model(),"curve_identity":c.identity(total,full_ring,n,g),
        "point_count":total,"trace":t,"ordinary":ordinary,"subgroup":{"order":n,"generator":g,"cofactor":total/n},
        "full_field_frobenius_action":{"eigenvalue":1,"additional_orbit_gain":false},"ring_binding":binding,
        "explicit_map_search":{"coverage":coverage,"maps_constructed_and_verified":maps.len(),"maps":maps.iter().map(Map::record).collect::<Vec<_>>()},
        "rho_structural":sym,"index_calculus":{"status":"unresolved","relation_model":null,"verified_independent_relations":null,"total_cost":null},
        "end_to_end_cost":null,"speedup":null,"implementation":"native_rust_variable_time_laboratory"});
    let units = if let Some(d) = full_ring {
        let mut o = orders::scan(d, bound, orders::Limits::default())?;
        o["curve_binding"] = binding;
        let units = o["order"]["geometric_unit_count"].as_u64().unwrap();
        report["order_search"] = o;
        Some(units)
    } else {
        None
    };
    if count_ops {
        report["scalar_operation_diagnostic"] = scalar::diagnostic(c, &maps, n, g)?;
    }
    report["decision_tree"] = json!({"ring_binding":report["ring_binding"]["status"],
        "scalar_multiplication":{"status":if maps.is_empty(){"no_match_in_declared_map_families"}else{"verified_map_candidates"},
            "verified_map_count":maps.len(),"measured_gain_scope":null,"native_timing":"unmeasured"},
        "rho":{"extra_geometric_units":match units {None=>"unresolved",Some(v) if v>2=>"present",_=>"excluded_for_certified_ring"},
            "distinct_subgroup_actions":report["rho_structural"]["distinct_action_count_including_negation"],
            "cheap_orbit_handling":"unresolved","measured_gain":"unresolved"},"index_calculus":report["index_calculus"]});
    Ok(report)
}

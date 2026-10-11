//! 4N = (2a + sb)^2 + |D|b^2, with s = D mod 2.
use crate::{Result, SCHEMA};
use serde_json::{json, Value};
use std::time::Instant;

#[derive(Clone, Copy, Debug)]
pub struct Limits {
    pub work: usize,
    pub candidates: usize,
    pub seconds: f64,
}
impl Default for Limits {
    fn default() -> Self {
        Self {
            work: 1_000_000,
            candidates: 100_000,
            seconds: 10.0,
        }
    }
}
impl Limits {
    pub fn check(self) -> Result<()> {
        if self.work == 0
            || self.work > 10_000_000
            || self.candidates == 0
            || self.candidates > 1_000_000
            || !self.seconds.is_finite()
            || self.seconds <= 0.0
        {
            return Err("positive finite limits required; max-work <= 10000000 and max-candidates <= 1000000".into());
        }
        Ok(())
    }
}

pub fn isqrt(n: u64) -> u64 {
    if n < 2 {
        return n;
    }
    let mut x = 1u64 << (64 - n.leading_zeros()).div_ceil(2);
    loop {
        let y = (x + n / x) / 2;
        if y >= x {
            return x;
        }
        x = y;
    }
}
pub fn gcd(mut a: i64, mut b: i64) -> i64 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a.abs()
}
pub fn valid_d(d: i64) -> Result<()> {
    if !(-1_000_000_000_000..=-3).contains(&d) || ![0, 1].contains(&d.rem_euclid(4)) {
        return Err("D must be negative, 0 or 1 modulo 4, and |D| <= 10^12".into());
    }
    Ok(())
}
pub fn factor(mut n: i64) -> Result<Vec<(i64, u32)>> {
    if !(1..=1_000_000_000_000).contains(&n) {
        return Err("factorization input outside 1..10^12".into());
    }
    let mut factors = Vec::new();
    let mut d = 2;
    let mut trials = 0;
    while d <= n / d {
        if trials >= 1_000_000 {
            return Err("factorization budget exceeded".into());
        }
        trials += 1;
        let mut e = 0;
        while n % d == 0 {
            n /= d;
            e += 1;
        }
        if e != 0 {
            factors.push((d, e));
        }
        d = if d == 2 { 3 } else { d + 2 };
    }
    if n > 1 {
        factors.push((n, 1));
    }
    Ok(factors)
}
pub fn fundamental(d: i64) -> Result<(i64, i64)> {
    valid_d(d)?;
    let mut s = -1;
    for (p, e) in factor(-d)? {
        if e % 2 != 0 {
            s *= p;
        }
    }
    let dk = if s.rem_euclid(4) == 1 { s } else { 4 * s };
    let f = isqrt((d / dk) as u64) as i64;
    if f * f * dk != d {
        return Err("invalid conductor decomposition".into());
    }
    Ok((dk, f))
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Element {
    pub a: i64,
    pub b: i64,
    pub degree: i64,
    pub trace: i64,
}
impl Element {
    pub fn record(self) -> Value {
        json!({"a":self.a,"b":self.b,"degree":self.degree,"trace":self.trace,
        "characteristic_polynomial_coefficients":[1,-self.trace,self.degree],
        "geometrically_non_scalar":self.b!=0,"map_status":"abstract_ring_element"})
    }
}
pub fn norm(d: i64, a: i64, b: i64) -> i128 {
    let (a, b, s) = (a as i128, b as i128, d.rem_euclid(2) as i128);
    a * a + s * a * b + ((s * s - d as i128) / 4) * b * b
}
pub fn enumerate(d: i64, bound: i64, non_scalar: bool, limits: Limits) -> Result<Value> {
    valid_d(d)?;
    limits.check()?;
    if !(1..=1_000_000_000_000).contains(&bound) {
        return Err("degree bound must be in 1..10^12".into());
    }
    let b_bound = isqrt((4 * bound / -d) as u64) as i64;
    let start = Instant::now();
    let mut work = 0;
    let mut stop = None;
    let mut items = Vec::new();
    'outer: for b in -b_bound..=b_bound {
        if work >= limits.work {
            stop = Some("work_limit");
            break;
        }
        if start.elapsed().as_secs_f64() >= limits.seconds {
            stop = Some("time_limit");
            break;
        }
        if non_scalar && b == 0 {
            work += 1;
            continue;
        }
        let t_bound = isqrt((4 * bound as i128 + d as i128 * b as i128 * b as i128) as u64) as i64;
        let first = -t_bound + (-t_bound - d.rem_euclid(2) * b).rem_euclid(2);
        for t in (first..=t_bound).step_by(2) {
            if work >= limits.work {
                stop = Some("work_limit");
                break 'outer;
            }
            if start.elapsed().as_secs_f64() >= limits.seconds {
                stop = Some("time_limit");
                break 'outer;
            }
            work += 1;
            let a = (t - d.rem_euclid(2) * b) / 2;
            let n = norm(d, a, b) as i64;
            if n == 0 {
                continue;
            }
            if items.len() >= limits.candidates {
                stop = Some("candidate_limit");
                break 'outer;
            }
            items.push(Element {
                a,
                b,
                degree: n,
                trace: t,
            });
        }
    }
    items.sort_by_key(|v| (v.degree, v.trace, v.a, v.b));
    Ok(
        json!({"status":if stop.is_none(){"complete"}else{"incomplete"},"stop_reason":stop,
        "degree_bound":bound,"limits":{"max_work":limits.work,"max_candidates":limits.candidates,"seconds":limits.seconds},
        "non_scalar_only":non_scalar,"signs_and_conjugates":"both_retained","b_bound":b_bound,"work_used":work,
        "candidate_count":items.len(),"candidates":items.into_iter().map(Element::record).collect::<Vec<_>>(),
        "completeness_certificate":{"identity":"4*N=(2*a+s*b)^2+abs(D)*b^2","enumeration_complete":stop.is_none()},
        "elapsed_seconds":null,"time_limit_applies_to_enumeration":true}),
    )
}
pub fn reduced_forms(d: i64, max_work: usize) -> Result<Value> {
    valid_d(d)?;
    let mut work = 0;
    let mut forms = Vec::new();
    for a in 1..=isqrt((-d / 3) as u64) as i64 {
        let first = -a + (d - a).rem_euclid(2);
        for b in (first..=a).step_by(2) {
            if work >= max_work {
                return Ok(
                    json!({"status":"incomplete","class_number":null,"forms":forms,"work_used":work}),
                );
            }
            work += 1;
            let num = b * b - d;
            if num % (4 * a) != 0 {
                continue;
            }
            let c = num / (4 * a);
            if a > c || gcd(gcd(a, b.abs()), c) != 1 || ((b.abs() == a || a == c) && b < 0) {
                continue;
            }
            forms.push(json!({"a":a,"b":b,"c":c}));
        }
    }
    Ok(json!({"status":"complete","class_number":forms.len(),"forms":forms,"work_used":work}))
}
pub fn scan(d: i64, bound: i64, limits: Limits) -> Result<Value> {
    let search = enumerate(d, bound, true, limits)?;
    let units = enumerate(d, 1, false, Limits::default())?;
    let classes = reduced_forms(d, limits.work)?;
    let (dk, f) = fundamental(d)?;
    let s = d.rem_euclid(2);
    let status = if search["status"] != "complete" {
        "search_incomplete"
    } else if search["candidate_count"] == 0 {
        "no_candidates_within_bound"
    } else {
        "abstract_candidates"
    };
    Ok(
        json!({"schema":SCHEMA,"input":{"discriminant":d,"degree_bound":bound},
        "curve_binding":{"status":"not_bound_to_a_curve","exclusions":"exact_for_O_D; conditional_for_any_named_curve"},
        "order":{"basis":{"omega":"(s+sqrt(D))/2","s":s},"norm_coefficients":[1,s,(s*s-d)/4],
            "fundamental_discriminant":dk,"conductor":f,"minimum_non_scalar_degree":(-d+3)/4,
            "minimum_degree_witness":Element{a:0,b:1,degree:(s*s-d)/4,trace:s}.record(),
            "minimum_degree_proof":"b!=0; minimize trace at |b|=1; |b|>=2 gives N>=abs(D)",
            "geometric_units":units["candidates"],"geometric_unit_count":units["candidate_count"],
            "class_forms":classes,"ideal_classes_are_not_self_endomorphisms":true},"norm_search":search,
        "decision_tree":{"scalar_multiplication":{"status":status,"map_construction":"unresolved","benchmark":"unresolved"},
            "rho":{"extra_geometric_automorphisms":if units["candidate_count"].as_u64().unwrap()>2 {"present"}else{"excluded_for_this_order"},
                "cheap_orbit_handling":"unresolved","measured_gain":"unresolved"},
            "index_calculus":{"status":"unresolved","relation_model":null,"total_cost":null}}}),
    )
}
pub fn summary(v: &Value) -> Value {
    json!({"discriminant":v["input"]["discriminant"],"minimum_non_scalar_degree":v["order"]["minimum_non_scalar_degree"],
        "geometric_unit_count":v["order"]["geometric_unit_count"],"class_number":v["order"]["class_forms"]["class_number"],
        "abstract_candidate_count":v["norm_search"]["candidate_count"],"norm_search_status":v["norm_search"]["status"],
        "curve_binding":"not_bound_to_a_curve","measured_speedup":null})
}

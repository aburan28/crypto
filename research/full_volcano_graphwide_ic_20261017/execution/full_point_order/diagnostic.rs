//! Standard smooth-order PH and reduced Tate/Miller diagnostic. No new identity.
//! Miller arithmetic copied from the separately frozen finite_field_k1.rs;
//! that earlier source and its measured binaries are not edited.
use crypto_lib::cryptanalysis::{
    ecbench::canonical::sha256_hex,
    ic_boundary::{CountedGroup, GroupOps, IncrementalGauss, PrimeCurve, PrimePoint, RowStatus},
};
use rusqlite::{Connection, OpenFlags};
use serde_json::{json, Value};
use std::{env, fs, path::Path, time::Instant};

#[derive(Default)]
struct Work {
    ec: GroupOps,
    muls: u64,
    squares: u64,
    inversions: u64,
    exponents: u64,
    adds_subs: u64,
    square_roots: u64,
    shift_candidates: u64,
    miller_steps: u64,
}
impl Work {
    fn mul(&mut self, a: u64, b: u64, p: u64) -> u64 {
        self.muls += 1;
        ((a as u128 * b as u128) % p as u128) as u64
    }
    fn sqr(&mut self, a: u64, p: u64) -> u64 {
        self.squares += 1;
        ((a as u128 * a as u128) % p as u128) as u64
    }
    fn add(&mut self, a: u64, b: u64, p: u64) -> u64 {
        self.adds_subs += 1;
        (a + b) % p
    }
    fn sub(&mut self, a: u64, b: u64, p: u64) -> u64 {
        self.adds_subs += 1;
        (a + p - b) % p
    }
    fn pow(&mut self, mut a: u64, mut n: u64, p: u64) -> u64 {
        self.exponents += 1;
        let mut b = 1;
        while n != 0 {
            if n & 1 != 0 {
                b = self.mul(b, a, p);
            }
            n >>= 1;
            if n != 0 {
                a = self.sqr(a, p);
            }
        }
        b
    }
    fn inv(&mut self, a: u64, p: u64) -> Option<u64> {
        if a == 0 {
            return None;
        }
        self.inversions += 1;
        Some(self.pow(a, p - 2, p))
    }
    fn point(&mut self, c: &PrimeCurve, x: u64) -> Option<PrimePoint> {
        let x = x % c.p;
        let x2 = self.sqr(x, c.p);
        let x3 = self.mul(x2, x, c.p);
        let ax = self.mul(c.a, x, c.p);
        let t = self.add(x3, ax, c.p);
        let rhs = self.add(t, c.b, c.p);
        self.square_roots += 1;
        c.sqrt(rhs).map(|y| PrimePoint::affine(x, y))
    }
    fn json(&self) -> Value {
        json!({"ec_group":self.ec,"ec_field":self.ec.field_ops,
            "explicit_fp":{"multiplications":self.muls,"squarings":self.squares,
                "inversions":self.inversions,"exponentiations":self.exponents,
                "additions_subtractions":self.adds_subs},
            "explicit_fp_accounting":"Inversion exponentiation is expanded into the multiplication/square counts; do not price it twice. EC field ledgers are separate.",
            "sqrt_calls_unexpanded":self.square_roots,"shift_candidates":self.shift_candidates,
            "miller_steps":self.miller_steps})
    }
}

// g_(A,B)(R) = line_(A,B)(R) / vertical_(A+B)(R), with vertical special cases.
fn line(
    c: &PrimeCurve,
    a: PrimePoint,
    b: PrimePoint,
    q: PrimePoint,
    w: &mut Work,
) -> Option<(u64, PrimePoint)> {
    if q.infinity {
        return None;
    }
    if a.infinity {
        return Some((1, b));
    }
    if b.infinity {
        return Some((1, a));
    }
    let p = c.p;
    if a.x == b.x && (a.y != b.y || a.y == 0) {
        return Some((w.sub(q.x, a.x, p), PrimePoint::INFINITY));
    }
    let (num, den, sum) = if a == b {
        let x2 = w.sqr(a.x, p);
        let three_x2 = w.mul(3, x2, p);
        (
            w.add(three_x2, c.a, p),
            w.add(a.y, a.y, p),
            c.double(&mut w.ec, a),
        )
    } else {
        (
            w.sub(b.y, a.y, p),
            w.sub(b.x, a.x, p),
            c.add(&mut w.ec, a, b),
        )
    };
    let inv = w.inv(den, p)?;
    let slope = w.mul(num, inv, p);
    let dx = w.sub(q.x, a.x, p);
    let sx = w.mul(slope, dx, p);
    let dy = w.sub(q.y, a.y, p);
    let numerator = w.sub(dy, sx, p);
    if sum.infinity {
        return Some((numerator, sum));
    }
    let vertical = w.sub(q.x, sum.x, p);
    let inverse = w.inv(vertical, p)?;
    Some((w.mul(numerator, inverse, p), sum))
}

fn miller(c: &PrimeCurve, g: PrimePoint, q: PrimePoint, r: u64, w: &mut Work) -> Option<u64> {
    let mut f = 1;
    let mut t = g;
    let top = 63 - r.leading_zeros();
    for bit in (0..top).rev() {
        w.miller_steps += 1;
        let (v, next) = line(c, t, t, q, w)?;
        let square = w.sqr(f, c.p);
        f = w.mul(square, v, c.p);
        t = next;
        if r & (1 << bit) != 0 {
            let (v, next) = line(c, t, g, q, w)?;
            f = w.mul(f, v, c.p);
            t = next;
        }
    }
    if !t.infinity || f == 0 {
        None
    } else {
        Some(f)
    }
}

fn tate(c: &PrimeCurve, g: PrimePoint, q: PrimePoint, r: u64, w: &mut Work) -> Result<u64, String> {
    if (c.p - 1) % r != 0 || q.infinity {
        return Err("k=1 requires r | p-1 and nonzero target".into());
    }
    for x in 1..=128 {
        w.shift_candidates += 1;
        let Some(s) = w.point(c, x) else {
            continue;
        };
        let qs = c.add(&mut w.ec, q, s);
        let Some(num) = miller(c, g, qs, r, w) else {
            continue;
        };
        let Some(den) = miller(c, g, s, r, w) else {
            continue;
        };
        let Some(inv) = w.inv(den, c.p) else {
            continue;
        };
        let ratio = w.mul(num, inv, c.p);
        let value = w.pow(ratio, (c.p - 1) / r, c.p);
        if value != 0 {
            return Ok(value);
        }
    }
    Err("shift search exhausted".into())
}

use std::{
    collections::BTreeMap,
    io::{BufRead, BufReader, BufWriter, Write},
};
const C_P: u64 = 1637958981148213249;
const C_N: u64 = C_P - 1;
const PREFIX_SHA: &str = "468ce8320b6b7deae37198ca3afc07534fa4fe0e950f621e45ad322f23a2c881";

#[derive(Default)]
struct Trace {
    phases: BTreeMap<String, Work>,
    counters: BTreeMap<String, u64>,
    details: Vec<Value>,
}
impl Trace {
    fn w(&mut self, phase: &str) -> &mut Work {
        self.phases.entry(phase.into()).or_default()
    }
    fn inc(&mut self, name: &str, n: u64) {
        *self.counters.entry(name.into()).or_default() += n;
    }
    fn json(&self) -> Value {
        json!({"phases":self.phases.iter().map(|(k,v)|(k.clone(),v.json())).collect::<BTreeMap<_,_>>(),"counters":self.counters,"details":self.details})
    }
}
fn factor(n: u64) -> Result<Vec<(u64, u32, u64)>, String> {
    if n < 2 {
        return Err("point order must exceed one".into());
    }
    let mut left = n;
    let mut out = vec![];
    for q in [2, 3, 7] {
        let mut e = 0;
        let mut m = 1;
        while left % q == 0 {
            left /= q;
            e += 1;
            m *= q;
        }
        if e > 0 {
            out.push((q, e, m));
        }
    }
    if left != 1 {
        return Err("unsupported factor in point order".into());
    }
    Ok(out)
}
fn point_order(c: &PrimeCurve, p: PrimePoint, n: u64, w: &mut Work) -> Result<u64, String> {
    factor(n)?;
    if p.infinity || !c.mul(&mut w.ec, p, n).infinity {
        return Err("point/order annihilation failed".into());
    }
    let mut m = n;
    for q in [2, 3, 7] {
        while m % q == 0 && c.mul(&mut w.ec, p, m / q).infinity {
            m /= q;
        }
    }
    if m < 2 || !c.mul(&mut w.ec, p, m).infinity {
        return Err("exact point order failed".into());
    }
    for (q, _, _) in factor(m)? {
        if c.mul(&mut w.ec, p, m / q).infinity {
            return Err("point order not exact".into());
        }
    }
    Ok(m)
}
struct Input {
    c: PrimeCurve,
    g: PrimePoint,
    q: PrimePoint,
    m: u64,
}
fn fixture(c: PrimeCurve, vertex: u64, target: u64) -> Result<(Input, u64, Value), String> {
    let start = Instant::now();
    let mut w = Work::default();
    let n = c.p - 1;
    let mut chosen = None;
    for offset in 1..=4096 {
        let x = (vertex.wrapping_mul(137) + offset) % c.p;
        if let Some(g) = w.point(&c, x) {
            let m = point_order(&c, g, n, &mut w)?;
            if m % 42 == 0 {
                chosen = Some((g, m));
                break;
            }
        }
    }
    let (g, m) = chosen.ok_or("generator search exhausted")?;
    let digest = sha256_hex(format!("C-full-point-v1:{vertex}:{target}").as_bytes());
    let x = u64::from_str_radix(&digest[..16], 16).unwrap() % (m - 1) + 1;
    let q = c.mul(&mut w.ec, g, x);
    let evidence = json!({"point_order":m.to_string(),"group_order":n.to_string(),"factors":factor(m)?,"generator":[g.x.to_string(),g.y.to_string()],"target":[q.x.to_string(),q.y.to_string(),q.infinity],"operations":w.json(),"elapsed_ns":start.elapsed().as_nanos().to_string(),"scope":"Exact order of selected point, not a maximality or full-group cyclicity certificate. Shared public fixture generation, outside both solver intervals."});
    Ok((Input { c, g, q, m }, x, evidence))
}

struct Base {
    q: u64,
    alpha: u64,
    values: Vec<u64>,
    logs: Vec<u64>,
    slots: Vec<(u64, Vec<u64>)>,
}
impl Base {
    fn decompose(&self, p: u64, value: u64, t: &mut Trace, tag: &str) -> Option<(Vec<u64>, Value)> {
        t.inc(&format!("{tag}_queries"), 1);
        for (i, (left, a)) in self.slots.iter().enumerate() {
            for (j, (right, b)) in self.slots.iter().enumerate() {
                t.inc(&format!("{tag}_pair_tests"), 1);
                if t.w(tag).mul(*left, *right, p) != value {
                    continue;
                }
                let row: Vec<_> = a.iter().zip(b).map(|(x, y)| (x + y) % self.q).collect();
                t.inc("scalar_modular_ops_uncharged", row.len() as u64);
                t.inc(&format!("{tag}_relations"), 1);
                return Some((
                    row.clone(),
                    json!({"left_slot":i,"right_slot":j,"row":row,"value":value.to_string()}),
                ));
            }
        }
        None
    }
    fn build(p: u64, q: u64, alpha: u64, t: &mut Trace) -> Result<Self, String> {
        if ![2, 3, 7].contains(&q) || alpha == 1 || t.w("base_setup").pow(alpha, q, p) != 1 {
            return Err("invalid prime destination generator".into());
        }
        let size = if q <= 3 { 1 } else { 2 };
        let mut values = vec![];
        for prime in [
            2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83,
            89, 97,
        ] {
            t.inc("base_prime_candidates", 1);
            let z = t.w("base_setup").pow(prime, (p - 1) / q, p);
            if z <= 1 {
                continue;
            }
            let inv = t.w("base_setup").inv(z, p).ok_or("zero base")?;
            if values.contains(&z) || values.contains(&inv) {
                continue;
            }
            values.push(z);
            if values.len() == size {
                break;
            }
        }
        if values.len() != size {
            return Err("factor-base search exhausted".into());
        }
        let mut slots = vec![(1, vec![0; size])];
        for (i, z) in values.iter().enumerate() {
            let mut row = vec![0; size];
            row[i] = 1;
            slots.push((*z, row));
            let inv = t.w("base_setup").inv(*z, p).unwrap();
            if inv != *z {
                let mut row = vec![0; size];
                row[i] = q - 1;
                slots.push((inv, row));
            }
        }
        let mut base = Self {
            q,
            alpha,
            values,
            logs: vec![],
            slots,
        };
        let mut matrix = IncrementalGauss::new(size, q);
        let mut walk = 1;
        for shift in 0..q {
            if let Some((row, certificate)) = base.decompose(p, walk, t, "base_relations") {
                if matrix.add_row(row, shift) == RowStatus::Inconsistent {
                    return Err("inconsistent anchor relation".into());
                }
                t.details.push(json!({"kind":"anchor_relation","prime":q,"shift":shift,"certificate":certificate}));
            }
            if matrix.rank() == size {
                break;
            }
            walk = t.w("base_relations").mul(walk, alpha, p);
        }
        t.inc("matrix_entry_ops_uncharged", matrix.row_ops);
        t.inc("matrix_dependent_rows", matrix.dependent);
        if matrix.rank() != size {
            return Err("anchor rank exhausted".into());
        }
        for (i, z) in base.values.iter().enumerate() {
            let log = matrix.pinned(i).ok_or("base log unpinned")?;
            if t.w("base_verification").pow(alpha, log, p) != *z {
                return Err("base log verification failed".into());
            }
            base.logs.push(log);
        }
        t.details.push(json!({"kind":"base","prime":q,"alpha":alpha.to_string(),"values":base.values.iter().map(u64::to_string).collect::<Vec<_>>(),"recovered_logs":base.logs,"rank":matrix.rank()}));
        Ok(base)
    }
    fn solve(&self, p: u64, beta: u64, t: &mut Trace) -> Result<u64, String> {
        if beta == 0 || t.w("digit_membership").pow(beta, self.q, p) != 1 {
            return Err("destination digit outside subgroup".into());
        }
        let mut walk = beta;
        for shift in 0..self.q {
            if let Some((row, certificate)) = self.decompose(p, walk, t, "target_relations") {
                let sum: u64 = row.iter().zip(&self.logs).map(|(a, b)| a * b).sum();
                let d = (sum % self.q + self.q - shift) % self.q;
                t.inc("scalar_modular_ops_uncharged", 2 * row.len() as u64 + 3);
                if t.w("digit_verification").pow(self.alpha, d, p) != beta {
                    return Err("digit relation failed verification".into());
                }
                t.details.push(json!({"kind":"target_relation","prime":self.q,"shift":shift,"digit":d,"beta":beta.to_string(),"certificate":certificate}));
                return Ok(d);
            }
            walk = t.w("target_relations").mul(walk, self.alpha, p);
        }
        Err("target relation search exhausted".into())
    }
}

fn inverse_mod(a: u64, n: u64, t: &mut Trace) -> Result<u64, String> {
    let (mut r, mut next) = (n as i128, a as i128);
    let (mut x, mut y) = (0i128, 1i128);
    while next != 0 {
        let k = r / next;
        (r, next) = (next, r - k * next);
        (x, y) = (y, x - k * y);
        t.inc("crt_integer_divisions_uncharged", 1);
    }
    if r != 1 {
        return Err("CRT moduli not coprime".into());
    }
    Ok(x.rem_euclid(n as i128) as u64)
}
fn crt(residues: &[(u64, u64)], t: &mut Trace) -> Result<u64, String> {
    let (mut x, mut product) = (0u64, 1u64);
    for &(m, r) in residues {
        let inv = inverse_mod(product % m, m, t)?;
        let k = (((r + m - x % m) % m) as u128 * inv as u128 % m as u128) as u64;
        x = (x as u128 + product as u128 * k as u128)
            .try_into()
            .map_err(|_| "CRT overflow")?;
        product = product.checked_mul(m).ok_or("CRT modulus overflow")?;
        t.inc("crt_integer_ops_uncharged", 8);
    }
    Ok(x)
}
fn ec_digit(
    c: &PrimeCurve,
    g: PrimePoint,
    baby: &[PrimePoint],
    jump: PrimePoint,
    delta: PrimePoint,
    q: u64,
    t: &mut Trace,
) -> Result<u64, String> {
    let mut walk = delta;
    let b = baby.len() as u64;
    for i in 0..b {
        for (j, p) in baby.iter().enumerate() {
            t.inc("ec_digit_table_comparisons_uncharged", 1);
            let d = i * b + j as u64;
            if walk == *p && d < q {
                if c.mul(&mut t.w("digit_verification").ec, g, d) != delta {
                    return Err("EC digit failed verification".into());
                }
                return Ok(d);
            }
        }
        if i + 1 < b {
            walk = c.add(&mut t.w("ec_digit_search").ec, walk, jump);
        }
    }
    Err("EC digit BSGS exhausted".into())
}
fn solve(i: &Input, finite_field: bool, t: &mut Trace) -> Result<u64, String> {
    let c = &i.c;
    let components = factor(i.m)?;
    if (c.p - 1) % i.m != 0
        || i.g.infinity
        || !c.mul(&mut t.w("input_validation").ec, i.g, i.m).infinity
    {
        return Err("invalid exact input order".into());
    }
    for (q, _, _) in &components {
        if c.mul(&mut t.w("input_validation").ec, i.g, i.m / q)
            .infinity
        {
            return Err("input order is not exact".into());
        }
    }
    let mut residues = vec![];
    for (q, e, m) in components {
        let g = c.mul(&mut t.w("projection").ec, i.g, i.m / m);
        let target = c.mul(&mut t.w("projection").ec, i.q, i.m / m);
        let mut x = 0u64;
        let mut power = 1u64;
        if finite_field {
            let mut found = None;
            for rx in 1..=256 {
                t.inc("pairing_auxiliary_candidates", 1);
                let Some(r) = t.w("pairing_setup").point(c, 10000 + rx + q * 1024) else {
                    continue;
                };
                let Ok(alpha) = tate(c, g, r, m, t.w("pairing_setup")) else {
                    continue;
                };
                if t.w("pairing_setup").pow(alpha, m, c.p) == 1
                    && t.w("pairing_setup").pow(alpha, m / q, c.p) != 1
                {
                    found = Some((r, alpha));
                    break;
                }
            }
            let (r, alpha) = found.ok_or("full-order pairing image search exhausted")?;
            let beta = tate(c, target, r, m, t.w("target_pairing"))?;
            if t.w("target_pairing").pow(beta, m, c.p) != 1 {
                return Err("target image outside destination".into());
            }
            t.details.push(json!({"kind":"transfer","prime":q,"power_order":m.to_string(),"auxiliary":[r.x.to_string(),r.y.to_string()],"alpha":alpha.to_string(),"beta":beta.to_string(),"exact_image_order_checked":true,"varying_argument":"first"}));
            let gamma = t.w("digit_projection").pow(alpha, m / q, c.p);
            let base = Base::build(c.p, q, gamma, t)?;
            for digit in 0..e {
                let known = t.w("digit_projection").pow(alpha, x, c.p);
                let inv = t
                    .w("digit_projection")
                    .inv(known, c.p)
                    .ok_or("zero image")?;
                let residual = t.w("digit_projection").mul(beta, inv, c.p);
                let delta = t.w("digit_projection").pow(residual, m / (power * q), c.p);
                let d = base.solve(c.p, delta, t)?;
                x += d * power;
                power *= q;
                t.details.push(json!({"kind":"digit","prime":q,"index":digit,"digit":d,"partial":x.to_string()}));
            }
            if t.w("destination_verification").pow(alpha, x, c.p) != beta {
                return Err("prime-power destination check failed".into());
            }
        } else {
            let gamma = c.mul(&mut t.w("digit_setup").ec, g, m / q);
            let b = if q == 7 { 3 } else { 2 };
            let mut baby = vec![PrimePoint::INFINITY];
            for _ in 1..b {
                baby.push(c.add(&mut t.w("digit_setup").ec, *baby.last().unwrap(), gamma));
            }
            let jump_pos = c.mul(&mut t.w("digit_setup").ec, gamma, b as u64);
            let jump = c.neg(jump_pos);
            t.inc("ec_negations_uncharged", 1);
            let mut residual = target;
            let mut place = g;
            for digit in 0..e {
                let delta = c.mul(&mut t.w("digit_projection").ec, residual, m / (power * q));
                let d = ec_digit(c, gamma, &baby, jump, delta, q, t)?;
                x += d * power;
                power *= q;
                let correction = c.mul(&mut t.w("digit_update").ec, place, d);
                let neg = c.neg(correction);
                t.inc("ec_negations_uncharged", 1);
                residual = c.add(&mut t.w("digit_update").ec, residual, neg);
                if digit + 1 < e {
                    place = c.mul(&mut t.w("digit_update").ec, place, q);
                }
                t.details.push(json!({"kind":"digit","prime":q,"index":digit,"digit":d,"partial":x.to_string()}));
            }
            if !residual.infinity {
                return Err("EC prime-power residual nonzero".into());
            }
        }
        residues.push((m, x));
    }
    let x = crt(&residues, t)?;
    if x >= i.m || c.mul(&mut t.w("ec_verification").ec, i.g, x) != i.q {
        return Err("recovered scalar failed EC membership/verification".into());
    }
    Ok(x)
}

fn input_json(i: &Input) -> Value {
    json!({"p":i.c.p.to_string(),"a":i.c.a.to_string(),"b":i.c.b.to_string(),"G":[i.g.x.to_string(),i.g.y.to_string()],"Q":[i.q.x.to_string(),i.q.y.to_string(),i.q.infinity],"point_order":i.m.to_string(),"point_order_factors":factor(i.m).unwrap()})
}
fn arm(i: &Input, ff: bool) -> Value {
    let mut trace = Trace::default();
    let start = Instant::now();
    let result = solve(i, ff, &mut trace);
    let elapsed = start.elapsed().as_nanos().to_string();
    let mut row = json!({"method":if ff {"ff.k1_ph_relations.full_point_v1"}else{"ec.ph_bsgs.full_point_v1"},"status":if result.is_ok(){"verified"}else{"failed"},"end_to_end_ns":elapsed,"trace":trace.json(),"unpriced":["square-root internals; calls counted","integer modular arithmetic and comparisons; selected counters retained","Gaussian inversions and allocation; entry operations counted","JSON evidence bookkeeping"],"timing_status":"UNQUALIFIED_DESCRIPTIVE; not an ecbench speed claim"});
    match result {
        Ok(x) => row["recovered"] = json!(x.to_string()),
        Err(e) => row["reason"] = json!(e),
    }
    row
}
#[derive(Default)]
struct Stat {
    n: u64,
    min: f64,
    max: f64,
    sum: f64,
}
impl Stat {
    fn add(&mut self, x: f64) {
        if self.n == 0 {
            self.min = x;
            self.max = x;
        } else {
            self.min = self.min.min(x);
            self.max = self.max.max(x);
        }
        self.n += 1;
        self.sum += x;
    }
    fn json(&self) -> Value {
        json!({"n":self.n,"min":self.min,"max":self.max,"sum":self.sum,"mean":self.sum/self.n as f64})
    }
}
fn flatten(prefix: &str, v: &Value, out: &mut BTreeMap<String, f64>) {
    if let Some(n) = v.as_f64() {
        out.insert(prefix.into(), n);
    } else if let Some(map) = v.as_object() {
        for (k, v) in map {
            flatten(&format!("{prefix}.{k}"), v, out);
        }
    }
}
fn strip_clocks(v: &mut Value) {
    match v {
        Value::Object(m) => {
            m.remove("end_to_end_ns");
            m.remove("elapsed_ns");
            for x in m.values_mut() {
                strip_clocks(x);
            }
        }
        Value::Array(a) => {
            for x in a {
                strip_clocks(x);
            }
        }
        _ => {}
    }
}
fn main() {
    let a: Vec<_> = env::args().collect();
    if a.get(1).map(String::as_str) == Some("compare") {
        assert_eq!(a.len(), 5);
        assert!(!Path::new(&a[4]).exists());
        let mut checked = BTreeMap::new();
        for name in ["records.jsonl", "fixtures.jsonl", "truth.jsonl"] {
            let mut left =
                BufReader::new(fs::File::open(Path::new(&a[2]).join(name)).unwrap()).lines();
            let mut right =
                BufReader::new(fs::File::open(Path::new(&a[3]).join(name)).unwrap()).lines();
            loop {
                match (left.next(), right.next()) {
                    (None, None) => break,
                    (Some(x), Some(y)) => {
                        let mut x: Value = serde_json::from_str(&x.unwrap()).unwrap();
                        let mut y: Value = serde_json::from_str(&y.unwrap()).unwrap();
                        strip_clocks(&mut x);
                        strip_clocks(&mut y);
                        assert_eq!(x, y, "non-clock replay mismatch in {name}");
                        *checked.entry(name).or_insert(0u64) += 1;
                    }
                    _ => panic!("replay row count differs in {name}"),
                }
            }
        }
        let doc = json!({"ok":true,"checked":checked,"scope":"Every non-clock row/fixture/truth field equal; same-implementation replay, not independent arithmetic."});
        fs::write(&a[4], serde_json::to_vec_pretty(&doc).unwrap()).unwrap();
        println!("{doc}");
        return;
    }
    assert_eq!(
        a.len(),
        6,
        "diagnostic run PREFIX_DB C_CERT FRESH_OUT source|prefix"
    );
    assert_eq!(a[1], "run");
    let prefix = Path::new(&a[2]);
    let cert_path = Path::new(&a[3]);
    let out = Path::new(&a[4]);
    assert!(!out.exists());
    assert!(["source", "prefix"].contains(&a[5].as_str()));
    assert_eq!(sha256_hex(&fs::read(prefix).unwrap()), PREFIX_SHA);
    let cert: Value = serde_json::from_slice(&fs::read(cert_path).unwrap()).unwrap();
    assert_eq!(cert["db_sha256"], PREFIX_SHA);
    let uri = format!(
        "file:{}?immutable=1",
        prefix.canonicalize().unwrap().display()
    );
    let db = Connection::open_with_flags(
        uri,
        OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_URI,
    )
    .unwrap();
    let receipt: String = db
        .query_row(
            "SELECT value FROM metadata WHERE key='verification'",
            [],
            |r| r.get(0),
        )
        .unwrap();
    let receipt: Value = serde_json::from_str(&receipt).unwrap();
    assert_eq!(receipt["status"], "REPLAYED");
    assert_eq!(receipt["vertex_class_complete"], false);
    let plan: String = db
        .query_row("SELECT value FROM metadata WHERE key='plan'", [], |r| {
            r.get(0)
        })
        .unwrap();
    let plan: Value = serde_json::from_str(&plan).unwrap();
    assert_eq!(plan["order"], C_N);
    assert_eq!(plan["input"]["p"], C_P);
    let mut statement = db
        .prepare("SELECT id,a,b,expanded FROM nodes ORDER BY id")
        .unwrap();
    let all = statement
        .query_map([], |r| {
            Ok((
                r.get::<_, u64>(0)?,
                r.get::<_, u64>(1)?,
                r.get::<_, u64>(2)?,
                r.get::<_, u64>(3)?,
            ))
        })
        .unwrap()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert_eq!(all.len(), 1869);
    let cases: Vec<_> = all
        .into_iter()
        .filter(|v| a[5] == "prefix" || v.0 == 0)
        .collect();
    fs::create_dir(out).unwrap();
    let mut records = BufWriter::new(fs::File::create(out.join("records.jsonl")).unwrap());
    let mut fixtures = BufWriter::new(fs::File::create(out.join("fixtures.jsonl")).unwrap());
    let mut truth = BufWriter::new(fs::File::create(out.join("truth.jsonl")).unwrap());
    let mut groups = BTreeMap::<(String, String, String), BTreeMap<String, Stat>>::new();
    let mut populations = BTreeMap::<(String, String, String), u64>::new();
    let mut outcomes = BTreeMap::<String, u64>::new();
    let mut orders = BTreeMap::<String, u64>::new();
    let meta = json!({"schema":"C-full-selected-point-order/v1","mode":a[5],"vertices_requested":cases.len(),"method_records_requested":cases.len()*4,"prefix_sha256":PREFIX_SHA,"certificate_sha256":sha256_hex(&fs::read(cert_path).unwrap()),"source_sha256":sha256_hex(include_bytes!("diagnostic.rs")),"protocol_sha256":sha256_hex(include_bytes!("PROTOCOL.md")),"binary_sha256":sha256_hex(&fs::read(env::current_exe().unwrap()).unwrap()),"scope":"Every stored vertex in requested population; selected point exact order, not maximality. Original C full graph and larger finite-field backends remain open.","timing_status":"UNQUALIFIED_DESCRIPTIVE","known_order_factorization":[[2,26],[3,20],[7,1]],"fixed_library_source":"2bf7700b73b3a84810aa86cb29837d1359b28322"});
    fs::write(
        out.join("meta.json"),
        serde_json::to_vec_pretty(&meta).unwrap(),
    )
    .unwrap();
    for (vertex, ca, cb, expanded) in &cases {
        let annotation = &cert["rows"][*vertex as usize];
        assert_eq!(annotation["vertex"], *vertex);
        let conductor = annotation["conductor"].as_u64();
        for target in 0..2 {
            let generated = fixture(
                PrimeCurve {
                    p: C_P,
                    a: *ca,
                    b: *cb,
                },
                *vertex,
                target,
            );
            let mut rows = vec![];
            match generated {
                Ok((input, expected, evidence)) => {
                    let public = input_json(&input);
                    let binding = sha256_hex(&serde_json::to_vec(&public).unwrap());
                    *orders.entry(input.m.to_string()).or_default() += 1;
                    writeln!(fixtures,"{}",json!({"vertex":vertex,"target_index":target,"input":public,"workload_sha256":binding,"evidence":evidence})).unwrap();
                    writeln!(truth,"{}",json!({"vertex":vertex,"target_index":target,"workload_sha256":binding,"scalar":expected.to_string(),"scope":"Public fixture truth; never an argument to either solver."})).unwrap();
                    for ff in [false, true] {
                        let mut row = arm(&input, ff);
                        if row["status"] == "verified" {
                            let ok = row["recovered"] == expected.to_string();
                            row["known_answer_checked"] = json!(ok);
                            if !ok {
                                row["status"] = json!("failed");
                                row["reason"] = json!("runner known-answer check failed");
                            }
                        }
                        row["workload_sha256"] = json!(binding);
                        row["point_order"] = json!(input.m.to_string());
                        rows.push(row);
                    }
                }
                Err(error) => {
                    for method in [
                        "ec.ph_bsgs.full_point_v1",
                        "ff.k1_ph_relations.full_point_v1",
                    ] {
                        rows.push(json!({"method":method,"status":"failed","reason":error,"scope":"Fixture failure; counters/clock missing, not zero."}));
                    }
                }
            }
            for mut row in rows {
                row["vertex"] = json!(vertex);
                row["target_index"] = json!(target);
                row["prefix_expanded"] = json!(expanded == &1);
                row["conductor"] = json!(conductor);
                let method = row["method"].as_str().unwrap().to_owned();
                let status = row["status"].as_str().unwrap().to_owned();
                *outcomes.entry(format!("{method}/{status}")).or_default() += 1;
                let mut metrics = BTreeMap::new();
                if let Some(ns) = row["end_to_end_ns"].as_str() {
                    metrics.insert("end_to_end_ns".into(), ns.parse().unwrap());
                }
                flatten("phases", &row["trace"]["phases"], &mut metrics);
                flatten("counters", &row["trace"]["counters"], &mut metrics);
                for tag in ["base_relations", "target_relations"] {
                    if let Some(queries) = row["trace"]["counters"][format!("{tag}_queries")]
                        .as_u64()
                        .filter(|n| *n > 0)
                    {
                        let relations = row["trace"]["counters"][format!("{tag}_relations")]
                            .as_u64()
                            .unwrap_or(0);
                        metrics.insert(
                            format!("{tag}.equations_per_query"),
                            relations as f64 / queries as f64,
                        );
                    }
                }
                for stratum in [
                    "all_sampled".to_string(),
                    conductor.map(|n| n.to_string()).unwrap_or("unknown".into()),
                ] {
                    let key = (method.clone(), stratum, status.clone());
                    *populations.entry(key.clone()).or_default() += 1;
                    let bucket = groups.entry(key).or_default();
                    for (name, value) in &metrics {
                        bucket.entry(name.clone()).or_default().add(*value);
                    }
                }
                writeln!(records, "{row}").unwrap();
            }
        }
        if vertex % 100 == 0 {
            eprintln!("completed vertex {vertex}");
        }
    }
    records.flush().unwrap();
    fixtures.flush().unwrap();
    truth.flush().unwrap();
    let mut table = vec![];
    let mut csv = BufWriter::new(fs::File::create(out.join("statistics.csv")).unwrap());
    writeln!(csv, "method,stratum,status,metric,n,min,max,mean").unwrap();
    for ((method, stratum, status), metrics) in groups {
        let mut data = serde_json::Map::new();
        for (name, s) in metrics {
            writeln!(
                csv,
                "{method},{stratum},{status},{name},{},{},{},{}",
                s.n,
                s.min,
                s.max,
                s.sum / s.n as f64
            )
            .unwrap();
            data.insert(name, s.json());
        }
        let n = populations[&(method.clone(), stratum.clone(), status.clone())];
        table.push(
            json!({"method":method,"stratum":stratum,"status":status,"records":n,"metrics":data}),
        );
    }
    csv.flush().unwrap();
    let summary = json!({"meta":meta,"outcomes":outcomes,"point_order_histogram_targets_not_vertices":orders,"groups":table,"scope":"Native observed min/max/mean with present-value n. Separate method/status/stratum groups; unknown conductors retained. No qualified speed or population claim."});
    fs::write(
        out.join("summary.json"),
        serde_json::to_vec_pretty(&summary).unwrap(),
    )
    .unwrap();
    println!(
        "{}",
        json!({"outcomes":outcomes,"requested_records":cases.len()*4})
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    fn points(c: &PrimeCurve) -> Vec<PrimePoint> {
        let mut out = vec![];
        for x in 0..c.p {
            let rhs = ((x as u128 * x as u128 % c.p as u128 * x as u128
                + c.a as u128 * x as u128
                + c.b as u128)
                % c.p as u128) as u64;
            if let Some(y) = c.sqrt(rhs) {
                out.push(PrimePoint::affine(x, y));
                if y != 0 {
                    out.push(PrimePoint::affine(x, c.p - y));
                }
            }
        }
        out
    }
    fn control(p: u64) -> (PrimeCurve, PrimePoint, u64) {
        for a in 0..p {
            for b in 0..p {
                if (4 * a * a % p * a + 27 * b * b) % p == 0 {
                    continue;
                }
                let c = PrimeCurve { p, a, b };
                let ps = points(&c);
                if ps.len() as u64 + 1 != p - 1 {
                    continue;
                }
                let (g, m) = ps
                    .iter()
                    .map(|g| {
                        (
                            *g,
                            point_order(&c, *g, p - 1, &mut Work::default()).unwrap(),
                        )
                    })
                    .max_by_key(|(_, m)| *m)
                    .unwrap();
                return (c, g, m);
            }
        }
        panic!("no control");
    }
    #[test]
    fn every_control_scalar_recovers_in_both_arms() {
        for p in [109, 113, 127] {
            let (c, g, m) = control(p);
            assert_eq!(points(&c).len() as u64 + 1, p - 1);
            for x in 0..m {
                let q = c.mul(&mut GroupOps::default(), g, x);
                let input = Input {
                    c: c.clone(),
                    g,
                    q,
                    m,
                };
                for ff in [false, true] {
                    assert_eq!(
                        solve(&input, ff, &mut Trace::default()).unwrap(),
                        x,
                        "p={p}, m={m}, ff={ff}, x={x}"
                    );
                }
            }
        }
    }
    #[test]
    fn all_destination_generators_and_digits() {
        for (p, q) in [(109, 2), (109, 3), (113, 7)] {
            let mut w = Work::default();
            let z = (2..p)
                .map(|x| w.pow(x, (p - 1) / q, p))
                .find(|z| *z != 1)
                .unwrap();
            for u in 1..q {
                let alpha = w.pow(z, u, p);
                let base = Base::build(p, q, alpha, &mut Trace::default()).unwrap();
                for x in 0..q {
                    let beta = w.pow(alpha, x, p);
                    assert_eq!(base.solve(p, beta, &mut Trace::default()).unwrap(), x);
                }
            }
        }
    }
    #[test]
    fn original_source_full_point_orders() {
        let c = PrimeCurve {
            p: C_P,
            a: 571607911220311176,
            b: 29329524900995129,
        };
        for target in 0..2 {
            let (i, x, _) = fixture(c.clone(), 0, target).unwrap();
            assert!(i.m > 7);
            assert_eq!(i.m % 42, 0);
            for ff in [false, true] {
                assert_eq!(solve(&i, ff, &mut Trace::default()).unwrap(), x);
            }
        }
    }
    #[test]
    fn reject_wrong_order_and_degenerate_generator() {
        assert!(Base::build(113, 7, 1, &mut Trace::default()).is_err());
        let (c, g, m) = control(127);
        let q = c.mul(&mut GroupOps::default(), g, 3);
        let i = Input { c, g, q, m: m / 2 };
        for ff in [false, true] {
            assert!(solve(&i, ff, &mut Trace::default()).is_err());
        }
    }
    #[test]
    fn reject_rank_two_target_outside_actual_cyclic_subgroup() {
        let (c, g, m) = control(109);
        assert!(m < c.p - 1);
        let multiples: Vec<_> = (0..m)
            .map(|x| c.mul(&mut GroupOps::default(), g, x))
            .collect();
        let q = points(&c)
            .into_iter()
            .find(|q| !multiples.contains(q) && c.mul(&mut GroupOps::default(), *q, m).infinity)
            .unwrap();
        let i = Input { c, g, q, m };
        for ff in [false, true] {
            assert!(solve(&i, ff, &mut Trace::default()).is_err());
        }
    }
}

//! B2b's measurement 6, the subfield sweep, natively
//! (`rounds/B2b-subfield-kic-certificates/sweep.py`, ported, with the
//! conformance generator's arithmetic it used, `conformance/v2/`,
//! `v2-b2/` and `v2-b2b/make_cases.py`).
//!
//! On every pair `k` in 2–8, odd `e ≥ 3`, `16 ≤ n = k·e ≤ 40`, the first
//! curve `y² + xy = x³ + ax² + b` over GF(2^k) by the script's rule (`a`
//! in 0, 1, then the rest ascending; `b` ascending; least subfield
//! GF(2^k); `#E`'s largest prime factor `r > h` with `r² ∤ #E`) runs as
//! `kic` paired with `rho-negation` on two known-answer targets.  The
//! documents are built here, and every answer is replayed here as
//! `[d]G = Q`, in arithmetic of this module's own that shares nothing with
//! the tool.  Untimed: it runs under the benchmark lock.
//!
//! The port is checked by writing B2b's frozen C059 and C063 documents
//! again, byte for byte: the same field, subfield count, order, prime
//! factor, hashed generator, hashed logarithm and recipe rule.

use std::io::Read;
use std::path::Path;
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

use num_bigint::BigUint;
use num_traits::ToPrimitive;

use super::bround::{clmul_mod, find_irreducible_sparse, is_prime, largest_prime_factor};
use super::json::{self, obj, J};

const LABEL: &str = "ic-b2b-sweep";
const TARGETS: u32 = 2;
const TIMEOUT_S: u64 = 900;
const PANIC_STATUS: i32 = 101;
/// B2b's `RHO_SEED`: suite v1's rho seed base plus one.
pub const RHO_SEED: i128 = 0x230000 + 1;

// ── GF(2^n), n ≤ 63 ────────────────────────────────────────────────

#[derive(Clone, Copy)]
pub struct Field {
    pub n: u32,
    pub f: u64,
}

fn deg(x: u64) -> i32 {
    63 - x.leading_zeros() as i32
}

impl Field {
    pub fn new(n: u32) -> Result<Field, String> {
        Ok(Field {
            n,
            f: find_irreducible_sparse(n)?,
        })
    }

    pub fn mul(&self, a: u64, b: u64) -> u64 {
        clmul_mod(a, b, self.f, self.n)
    }

    pub fn sqr(&self, a: u64) -> u64 {
        self.mul(a, a)
    }

    /// Algorithm 2.48 of Hankerson, Menezes and Vanstone.
    pub fn inv(&self, a: u64) -> u64 {
        assert!(a != 0, "zero has no inverse");
        let (mut u, mut v, mut g1, mut g2) = (a, self.f, 1u64, 0u64);
        while u != 1 {
            let mut j = deg(u) - deg(v);
            if j < 0 {
                std::mem::swap(&mut u, &mut v);
                std::mem::swap(&mut g1, &mut g2);
                j = -j;
            }
            u ^= v << j;
            g1 ^= g2 << j;
        }
        g1
    }

    pub fn trace(&self, c: u64) -> u64 {
        let (mut t, mut s) = (0, c);
        for _ in 0..self.n {
            t ^= s;
            s = self.sqr(s);
        }
        t
    }

    pub fn sqrt(&self, mut a: u64) -> u64 {
        for _ in 0..self.n - 1 {
            a = self.sqr(a);
        }
        a
    }

    /// The root of `z² + z = c` whose bit 0 is 0, by linear algebra over
    /// GF(2), as the generator finds it.
    pub fn solve_quadratic(&self, c: u64) -> Option<u64> {
        if self.trace(c) != 0 {
            return None;
        }
        let mut basis: Vec<(u64, u64)> = Vec::new();
        for i in 0..self.n {
            let (mut v, mut combo) = (self.sqr(1 << i) ^ (1 << i), 1u64 << i);
            for &(bv, bc) in &basis {
                if v ^ bv < v {
                    v ^= bv;
                    combo ^= bc;
                }
            }
            if v != 0 {
                basis.push((v, combo));
                basis.sort_unstable_by(|x, y| y.cmp(x));
            }
        }
        let (mut rest, mut z) = (c, 0u64);
        for &(bv, bc) in &basis {
            if rest ^ bv < rest {
                rest ^= bv;
                z ^= bc;
            }
        }
        if rest != 0 {
            return None;
        }
        z ^= z & 1;
        Some(z)
    }

    fn frobenius(&self, mut x: u64, k: u32) -> u64 {
        for _ in 0..k {
            x = self.sqr(x);
        }
        x
    }

    /// Every element of GF(2^k) inside GF(2^n), ascending: the kernel of
    /// `x ↦ x^(2^k) + x`.
    pub fn subfield(&self, k: u32) -> Vec<u64> {
        let (mut pivots, mut kernel): (Vec<(u64, u64)>, Vec<u64>) = (Vec::new(), Vec::new());
        for i in 0..self.n {
            let (mut v, mut c) = (self.frobenius(1 << i, k) ^ (1 << i), 1u64 << i);
            for &(pv, pc) in &pivots {
                if v ^ pv < v {
                    v ^= pv;
                    c ^= pc;
                }
            }
            if v != 0 {
                pivots.push((v, c));
                pivots.sort_unstable_by(|x, y| y.cmp(x));
            } else {
                kernel.push(c);
            }
        }
        let mut out: Vec<u64> = (0..1u64 << kernel.len())
            .map(|m| {
                kernel
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| m >> i & 1 == 1)
                    .fold(0, |x, (_, e)| x ^ e)
            })
            .collect();
        out.sort_unstable();
        out.dedup();
        out
    }

    /// The least `d` with `a` and `b` in GF(2^d).
    pub fn least_subfield(&self, a: u64, b: u64) -> u32 {
        (1..=self.n)
            .find(|d| {
                self.n.is_multiple_of(*d)
                    && self.frobenius(a, *d) == a
                    && self.frobenius(b, *d) == b
            })
            .expect("n itself qualifies")
    }

    /// The trace from GF(2^k) to GF(2) of `c` in GF(2^k).
    fn trace_over(&self, c: u64, k: u32) -> u64 {
        let (mut t, mut s) = (0, c);
        for _ in 0..k {
            t ^= s;
            s = self.sqr(s);
        }
        t
    }
}

// ── the curve ──────────────────────────────────────────────────────

pub type Point = Option<(u64, u64)>;

#[derive(Clone, Copy)]
pub struct Curve {
    pub k: Field,
    pub a: u64,
    pub b: u64,
}

fn shake(label: &str, nbytes: usize) -> BigUint {
    BigUint::from_bytes_be(&crypto_lib::hash::sha3::shake256(label.as_bytes(), nbytes))
}

impl Curve {
    pub fn on_curve(&self, p: &Point) -> bool {
        let Some((x, y)) = *p else { return true };
        let k = &self.k;
        let x2 = k.sqr(x);
        k.sqr(y) ^ k.mul(x, y) == k.mul(x2, x) ^ k.mul(self.a, x2) ^ self.b
    }

    pub fn add(&self, p: &Point, q: &Point) -> Point {
        let (Some((x1, y1)), Some((x2, y2))) = (*p, *q) else {
            return if p.is_none() { *q } else { *p };
        };
        let k = &self.k;
        if x1 == x2 {
            if y1 ^ y2 == x1 {
                return None;
            }
            return self.dbl(p);
        }
        let lam = k.mul(y1 ^ y2, k.inv(x1 ^ x2));
        let x3 = k.sqr(lam) ^ lam ^ x1 ^ x2 ^ self.a;
        Some((x3, k.mul(lam, x1 ^ x3) ^ x3 ^ y1))
    }

    pub fn dbl(&self, p: &Point) -> Point {
        let (x, y) = (*p)?;
        if x == 0 {
            return None;
        }
        let k = &self.k;
        let lam = x ^ k.mul(y, k.inv(x));
        let x3 = k.sqr(lam) ^ lam ^ self.a;
        Some((x3, k.sqr(x) ^ k.mul(lam ^ 1, x3)))
    }

    pub fn mul(&self, e: u128, p: &Point) -> Point {
        let mut r: Point = None;
        for i in (0..128 - e.leading_zeros()).rev() {
            r = self.dbl(&r);
            if e >> i & 1 == 1 {
                r = self.add(&r, p);
            }
        }
        r
    }

    fn lift(&self, x: u64) -> Point {
        let k = &self.k;
        if x == 0 {
            return Some((0, k.sqrt(self.b)));
        }
        let c = x ^ self.a ^ k.mul(self.b, k.inv(k.sqr(x)));
        let z = k.solve_quadratic(c)?;
        let p = Some((x, k.mul(x, z)));
        debug_assert!(self.on_curve(&p));
        p
    }

    fn element(&self, label: &str) -> u64 {
        let bytes = (self.k.n as usize).div_ceil(8);
        let v = shake(label, bytes);
        let mask = (BigUint::from(1u8) << self.k.n) - 1u8;
        (v & mask).to_u64().expect("below 2^n")
    }

    /// The generator's `point`: the first labelled abscissa that lifts.
    pub fn point(&self, label: &str) -> Point {
        for i in 0..10_000 {
            let l = if i == 0 {
                label.to_string()
            } else {
                format!("{label}:{i}")
            };
            if let Some(p) = self.lift(self.element(&l)) {
                return Some(p);
            }
        }
        None
    }

    /// The generator's `subgroup_point`: `h` times the first labelled point
    /// that it does not send to the identity.
    pub fn subgroup_point(&self, label: &str, h: u128, r: u128) -> Result<Point, String> {
        for i in 0..10_000 {
            let l = if i == 0 {
                label.to_string()
            } else {
                format!("{label}#{i}")
            };
            let g = self.mul(h, &self.point(&l));
            if g.is_some() {
                if self.mul(r, &g).is_some() {
                    return Err(format!("{label}: h·P has no order r"));
                }
                return Ok(g);
            }
        }
        Err(format!("{label}: no subgroup point"))
    }
}

/// `#E(GF(2^k))` for a curve with `a`, `b` in GF(2^k).
fn count_over_subfield(c: &Curve, k: u32, elements: &[u64]) -> i128 {
    let f = &c.k;
    let mut count = 2i128;
    for &x in elements {
        if x != 0 {
            let t = x ^ c.a ^ f.mul(c.b, f.inv(f.sqr(x)));
            if f.trace_over(t, k) == 0 {
                count += 2;
            }
        }
    }
    count
}

/// `#E(GF(q^e)) = q^e + 1 − s_e`, `s_0 = 2`, `s_1 = t`, `s_i = t s_{i−1} − q s_{i−2}`.
fn order_over_extension(t: i128, q: i128, e: u32) -> i128 {
    let (mut prev, mut cur) = (2i128, t);
    for _ in 1..e {
        (prev, cur) = (cur, t * cur - q * prev);
    }
    q.pow(e) + 1 - cur
}

/// The subfield curve's `#E`, its largest prime factor `r` and cofactor
/// `h`, with `r` prime, `r² ∤ #E`, and the order checked on labelled
/// points.
pub fn subfield_curve(
    n: u32,
    k: u32,
    a: u64,
    b: u64,
    label: &str,
) -> Result<(Curve, u128, u128, u128), String> {
    let field = Field::new(n)?;
    let c = Curve { k: field, a, b };
    if field.least_subfield(a, b) != k {
        return Err(format!("{label}: the least subfield is not GF(2^{k})"));
    }
    let elements = field.subfield(k);
    let q = 1i128 << k;
    let order = order_over_extension(q + 1 - count_over_subfield(&c, k, &elements), q, n / k);
    let order = u128::try_from(order).map_err(|_| "a negative order".to_string())?;
    let r = largest_prime_factor(order)?;
    let h = order / r;
    let prime = u64::try_from(r).is_ok_and(is_prime);
    if !prime || h % r == 0 {
        return Err(format!("{label}: r is not prime, or r² divides #E"));
    }
    for i in 0..8 {
        let p = c.point(&format!("{label}/order/{i}"));
        if c.mul(order, &p).is_some() {
            return Err(format!("{label}: #E does not annihilate a point"));
        }
    }
    Ok((c, order, r, h))
}

/// The first `(a, b)` over GF(2^k), in the generator's order, with least
/// subfield GF(2^k) whose `#E` passes `accept(r, h)` with `r` prime and
/// `r² ∤ #E`.
pub fn first_subfield_curve(
    n: u32,
    k: u32,
    accept: impl Fn(u128, u128) -> bool,
) -> Result<Option<(u64, u64)>, String> {
    let field = Field::new(n)?;
    let elements = field.subfield(k);
    let mut a_order = vec![0u64, 1];
    a_order.extend(elements.iter().copied().filter(|&x| x > 1));
    let q = 1i128 << k;
    for &a in &a_order {
        for &b in &elements {
            if b == 0 || field.least_subfield(a, b) != k {
                continue;
            }
            let c = Curve { k: field, a, b };
            let order =
                order_over_extension(q + 1 - count_over_subfield(&c, k, &elements), q, n / k);
            let order = u128::try_from(order).map_err(|_| "a negative order".to_string())?;
            let r = largest_prime_factor(order)?;
            let h = order / r;
            if u64::try_from(r).is_ok_and(is_prime) && h % r != 0 && accept(r, h) {
                return Ok(Some((a, b)));
            }
        }
    }
    Ok(None)
}

/// The generator's `known_log`: a labelled SHAKE-256 value in `[1, r)`.
pub fn known_log(label: &str, r: u128) -> u128 {
    let v = shake(label, 32) % BigUint::from(r - 1);
    v.to_u128().expect("below r") + 1
}

/// Ledger §20's recipe rules with the orbit's length `e` in the place of
/// `n`, as the generator writes them (Python's `round`, ties to even).
pub fn recipe(e: u32, r: u128, columns: u64, m: u32, seed: i128) -> J {
    let points = 2 * u64::from(e) * columns;
    let pp = (points * points) as f64;
    let max_trials = (64.0 * (2.02 * r as f64 / pp).ceil()).max(100_000.0);
    let unit = (2.0 * r as f64 / (f64::from(e) * pp)).round_ties_even();
    obj([
        ("summands", J::Int(3)),
        ("descent_summands", J::Int(m.into())),
        (
            "collection_window",
            J::Int(((points as f64 / 32.0).round_ties_even() as i128).max(1)),
        ),
        ("collection_aim", J::Bool(true)),
        ("solver", J::Str("pair_table".into())),
        ("seed", J::Int(seed)),
        (
            "max_trials",
            J::Int((max_trials as i128).min((1i128 << 32) - 1)),
        ),
        (
            "collection",
            obj([
                ("unit_trials", J::Int((unit as i128).clamp(16, 65536))),
                ("units", J::Int(1)),
                ("max_units", J::Int(100_000)),
            ]),
        ),
        (
            "factor_base",
            obj([
                ("mode", J::Str("spec".into())),
                (
                    "spec",
                    obj([
                        ("kind", J::Str("subgroup_orbits".into())),
                        ("points", J::Int(points.into())),
                        ("seed", J::Int(seed)),
                    ]),
                ),
            ]),
        ),
    ])
}

fn hx(v: u64) -> J {
    J::Str(format!("0x{v:x}"))
}

/// The generator's `binary_doc`, paired with `rho` routed.
pub fn binary_doc(
    name: &str,
    c: &Curve,
    r: u128,
    h: u128,
    g: &Point,
    target: J,
    ic_pipeline: &str,
    ic_recipe: J,
) -> Result<J, String> {
    let (gx, gy) = g.ok_or("the generator is the identity")?;
    Ok(obj([
        ("schema_version", J::Int(2)),
        ("name", J::Str(name.into())),
        (
            "field",
            obj([
                ("kind", J::Str("binary".into())),
                ("degree", J::Int(c.k.n.into())),
                ("modulus", hx(c.k.f)),
            ]),
        ),
        (
            "curve",
            obj([
                ("form", J::Str("binary_weierstrass".into())),
                ("a", hx(c.a)),
                ("b", hx(c.b)),
            ]),
        ),
        (
            "subgroup",
            obj([
                ("order", J::Str(r.to_string())),
                ("cofactor", J::Str(h.to_string())),
                ("generator", obj([("x", hx(gx)), ("y", hx(gy))])),
            ]),
        ),
        ("target", target),
        (
            "method",
            obj([
                ("solve", J::Str("paired".into())),
                (
                    "index_calculus",
                    obj([
                        ("pipeline", J::Str(ic_pipeline.into())),
                        ("recipe", ic_recipe),
                    ]),
                ),
                (
                    "rho",
                    obj([
                        ("pipeline", J::Str("auto".into())),
                        ("seed", J::Int(RHO_SEED)),
                    ]),
                ),
            ]),
        ),
    ]))
}

// ── the sweep ──────────────────────────────────────────────────────

pub fn pairs() -> Vec<(u32, u32)> {
    let mut out = Vec::new();
    for k in 2..9 {
        for e in (3..41).step_by(2) {
            if (16..=40).contains(&(k * e)) {
                out.push((k, e));
            }
        }
    }
    out
}

/// The sweep's curve on `(k, e)`, or none when no curve has `r > h`.
pub fn sweep_curve(k: u32, e: u32) -> Result<Option<(Curve, u128, u128, u128)>, String> {
    let n = k * e;
    match first_subfield_curve(n, k, |r, h| r > h)? {
        Some((a, b)) => subfield_curve(n, k, a, b, &format!("{LABEL}/{k}/{e}")).map(Some),
        None => Ok(None),
    }
}

pub struct Run {
    pub k: u32,
    pub e: u32,
    pub r: u128,
    pub h: u128,
    pub columns: u64,
    pub target: u32,
    pub known_log: u128,
    pub curve: Curve,
    pub generator: Point,
    pub doc: J,
}

/// Python's `round(x, 6)`: the float nearest the correctly rounded
/// six-decimal value.
fn round6(x: f64) -> f64 {
    format!("{x:.6}").parse().expect("a formatted float parses")
}

/// The runs, and the pairs on which no curve has `r > h`.
pub fn documents() -> Result<(Vec<Run>, Vec<String>), String> {
    let (mut out, mut none) = (Vec::new(), Vec::new());
    for (k, e) in pairs() {
        let Some((curve, _order, r, h)) = sweep_curve(k, e)? else {
            none.push(format!("k{k}e{e}"));
            continue;
        };
        let generator = curve.subgroup_point(&format!("{LABEL}/{k}/{e}/generator"), h, r)?;
        let cube = round6((r as f64).powf(1.0 / 3.0));
        let columns = (cube / f64::from(2 * e)).ceil().max(8.0) as u64;
        for t in 0..TARGETS {
            let d = known_log(&format!("{LABEL}/{k}/{e}/known_log/{t}"), r);
            let doc = binary_doc(
                &format!("B2b sweep: k = {k}, e = {e}, n = {}, target {t}", k * e),
                &curve,
                r,
                h,
                &generator,
                obj([("known_log", J::Str(d.to_string()))]),
                "kic",
                recipe(e, r, columns, 2, 201),
            )?;
            out.push(Run {
                k,
                e,
                r,
                h,
                columns,
                target: t,
                known_log: d,
                curve,
                generator,
                doc,
            });
        }
    }
    Ok((out, none))
}

fn replay(run: &Run, scalar: Option<&str>) -> bool {
    let Some(d) = scalar.and_then(|s| s.parse::<u128>().ok()) else {
        return false;
    };
    0 < d
        && d < run.r
        && run.curve.mul(d, &run.generator) == run.curve.mul(run.known_log, &run.generator)
}

fn dig<'a>(doc: Option<&'a J>, keys: &[&str]) -> Option<&'a J> {
    let mut cur = doc?;
    for k in keys {
        cur = cur.get(k)?;
    }
    Some(cur)
}

/// Every run, once, into `out` (which must not exist yet): each document,
/// each report, and `results.json`.  The second value is whether no run
/// had a finding.
pub fn sweep(ic: &Path, out: &Path) -> Result<(J, bool), String> {
    if out.exists() {
        return Err(format!(
            "{} exists; the sweep writes a fresh directory",
            out.display()
        ));
    }
    std::fs::create_dir_all(out.join("params")).map_err(|e| format!("{}: {e}", out.display()))?;
    let (runs_list, no_curve) = documents()?;
    let mut rows = Vec::new();
    for run in &runs_list {
        let stem = format!("k{}e{}-t{}", run.k, run.e, run.target);
        let params = out.join("params").join(format!("{stem}.json"));
        let text = json::dumps(&run.doc, 1) + "\n";
        std::fs::write(&params, &text).map_err(|e| format!("{}: {e}", params.display()))?;
        let report_path = out.join(format!("{stem}.report.json"));
        let started = Instant::now();
        let mut child = Command::new(ic)
            .arg("price")
            .arg("--params")
            .arg(&params)
            .args(["--json", "--out"])
            .arg(&report_path)
            .args(["--repeats", "1", "--repeats-fast", "1"])
            .env_clear()
            .env("RAYON_NUM_THREADS", "1")
            .env("PATH", "/usr/bin:/bin")
            .stdin(Stdio::null())
            .stdout(Stdio::null())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| format!("{}: {e}", ic.display()))?;
        let mut err = child.stderr.take().expect("piped stderr");
        let reader = std::thread::spawn(move || {
            let mut s = Vec::new();
            let _ = err.read_to_end(&mut s);
            String::from_utf8_lossy(&s).into_owned()
        });
        let deadline = started + Duration::from_secs(TIMEOUT_S);
        let status = loop {
            if let Some(s) = child.try_wait().map_err(|e| e.to_string())? {
                break Some(s);
            }
            if Instant::now() >= deadline {
                let _ = child.kill();
                let _ = child.wait();
                break None;
            }
            std::thread::sleep(Duration::from_millis(20));
        };
        let stderr = reader.join().unwrap_or_default();
        let wall = started.elapsed().as_secs_f64();
        let exit: J = match status {
            None => J::Str("timeout".into()),
            Some(s) => J::Int(s.code().map_or(-1, i128::from)),
        };
        let code = status.and_then(|s| s.code());
        let report = std::fs::read_to_string(&report_path)
            .ok()
            .and_then(|t| json::parse(&t).ok());
        let get = |keys: &[&str]| dig(report.as_ref(), keys).cloned();
        let mut findings = Vec::new();
        if code == Some(PANIC_STATUS) || stderr.contains("panicked") {
            findings.push("panic".to_string());
        }
        if !matches!(code, Some(0) | Some(1)) {
            findings.push(format!(
                "exit {}",
                json::dumps_line(&exit, false).trim_matches('"')
            ));
        }
        if report.is_none() && matches!(code, Some(0) | Some(1)) {
            findings.push("no report".into());
        }
        let status_s = get(&["status"]);
        let scalar = get(&["result", "scalar"]);
        let scalar_s = scalar.as_ref().and_then(J::as_str);
        if code == Some(0) && status_s.as_ref().and_then(J::as_str) == Some("complete") {
            if get(&["result", "verified"]) != Some(J::Bool(true)) {
                findings.push("complete but not verified".into());
            }
            if scalar_s != Some(run.known_log.to_string().as_str()) || !replay(run, scalar_s) {
                findings.push("wrong answer".into());
            }
            let route = (
                get(&["route", "ic", "pipeline"]),
                get(&["route", "rho", "pipeline"]),
            );
            if route
                != (
                    Some(J::Str("kic".into())),
                    Some(J::Str("rho-negation".into())),
                )
            {
                findings.push("route".into());
            }
        }
        let tail = if findings.is_empty() {
            String::new()
        } else {
            stderr
                .chars()
                .rev()
                .take(400)
                .collect::<Vec<_>>()
                .into_iter()
                .rev()
                .collect()
        };
        println!(
            "k={} e={} n={} t={}: exit {}, {}, {wall:.1} s{}",
            run.k,
            run.e,
            run.k * run.e,
            run.target,
            json::dumps_line(&exit, false),
            status_s
                .as_ref()
                .map_or("None".into(), |s| json::dumps_line(s, false)),
            if findings.is_empty() {
                String::new()
            } else {
                format!(", FINDINGS {findings:?}")
            }
        );
        rows.push(obj([
            ("k", J::Int(run.k.into())),
            ("e", J::Int(run.e.into())),
            ("n", J::Int((run.k * run.e).into())),
            ("r", J::Int(run.r as i128)),
            ("h", J::Int(run.h as i128)),
            ("columns", J::Int(run.columns.into())),
            ("target", J::Int(run.target.into())),
            ("known_log", J::Int(run.known_log as i128)),
            (
                "params_sha256",
                J::Str(hex::encode(crypto_lib::hash::sha256::sha256(
                    text.as_bytes(),
                ))),
            ),
            ("exit", exit),
            ("status", status_s.unwrap_or(J::Null)),
            ("scalar", scalar.unwrap_or(J::Null)),
            ("verified", get(&["result", "verified"]).unwrap_or(J::Null)),
            ("wall_s", J::Float((wall * 1000.0).round() / 1000.0)),
            (
                "ic_online_ms",
                get(&["median", "ic_online_wall_ms"]).unwrap_or(J::Null),
            ),
            (
                "rho_online_ms",
                get(&["median", "rho_online_wall_ms"]).unwrap_or(J::Null),
            ),
            (
                "findings",
                J::Arr(findings.into_iter().map(J::Str).collect()),
            ),
            ("stderr_tail", J::Str(tail)),
        ]));
    }
    let count = |pred: &dyn Fn(&J) -> bool| rows.iter().filter(|r| pred(r)).count() as i128;
    let complete = count(&|r| r.get("status").and_then(J::as_str) == Some("complete"));
    let with_findings = count(&|r| {
        r.get("findings")
            .and_then(J::as_arr)
            .is_some_and(|f| !f.is_empty())
    });
    let not_complete: Vec<J> = rows
        .iter()
        .filter(|r| r.get("status").and_then(J::as_str) != Some("complete"))
        .map(|r| {
            J::Str(format!(
                "k{}e{}-t{}",
                r.get("k")
                    .map_or(String::new(), |v| json::dumps_line(v, false)),
                r.get("e")
                    .map_or(String::new(), |v| json::dumps_line(v, false)),
                r.get("target")
                    .map_or(String::new(), |v| json::dumps_line(v, false)),
            ))
        })
        .collect();
    let summary = obj([
        ("runs", J::Int(rows.len() as i128)),
        (
            "no_curve",
            J::Arr(no_curve.into_iter().map(J::Str).collect()),
        ),
        ("complete", J::Int(complete)),
        ("not_complete", J::Arr(not_complete)),
        ("findings", J::Int(with_findings)),
        (
            "harness",
            J::Str("icprog b2b sweep (sweep.py, ported)".into()),
        ),
    ]);
    let doc = obj([("summary", summary.clone()), ("rows", J::Arr(rows))]);
    std::fs::write(out.join("results.json"), json::dumps(&doc, 1) + "\n")
        .map_err(|e| format!("{}: {e}", out.display()))?;
    Ok((summary, with_findings == 0))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn frozen(name: &str) -> String {
        std::fs::read_to_string(
            Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("research/ic_tool_program/conformance/v2-b2b/params")
                .join(name),
        )
        .unwrap()
    }

    /// C059: a curve over GF(4) taken over GF(2^22), its generator, its
    /// known logarithm and its recipe, as the generator wrote them.
    #[test]
    fn c059_is_written_again_byte_for_byte() {
        let field = Field::new(22).unwrap();
        let w = field.solve_quadratic(1).unwrap();
        let (c, _, r, h) = subfield_curve(22, 2, w, 1, "ic-conformance-b2b/curve-b").unwrap();
        assert_eq!((r, h), (2_097_349, 2));
        let g = c
            .subgroup_point("ic-conformance-b2b/curve-b/generator", h, r)
            .unwrap();
        let d = known_log("ic-conformance-b2b/curve-b/known_log", r);
        let doc = binary_doc(
            "C059: a curve over GF(4) taken over GF(2^22), kic paired with rho",
            &c,
            r,
            h,
            &g,
            obj([("known_log", J::Str(d.to_string()))]),
            "kic",
            recipe(11, r, 24, 2, 201),
        )
        .unwrap();
        assert_eq!(
            json::dumps(&doc, 1) + "\n",
            frozen("C059-subfield-curve-paired.json")
        );
    }

    /// C063: the first curve over GF(8) at n = 33 with r ≥ 2^29, found by
    /// the same search the sweep runs.
    #[test]
    fn c063_is_written_again_byte_for_byte() {
        let (a, b) = first_subfield_curve(33, 3, |r, h| r > h && r >= 1 << 29)
            .unwrap()
            .unwrap();
        let (c, _, r, h) = subfield_curve(33, 3, a, b, "ic-conformance-b2b/curve-c").unwrap();
        assert_eq!((r, h), (715_829_951, 12));
        let g = c
            .subgroup_point("ic-conformance-b2b/curve-c/generator", h, r)
            .unwrap();
        let d = known_log("ic-conformance-b2b/curve-c/known_log", r);
        let doc = binary_doc(
            "C063: a curve over GF(8) taken over GF(2^33), both pipelines routed",
            &c,
            r,
            h,
            &g,
            obj([("known_log", J::Str(d.to_string()))]),
            "auto",
            recipe(11, r, 48, 2, 201),
        )
        .unwrap();
        assert_eq!(
            json::dumps(&doc, 1) + "\n",
            frozen("C063-subfield-curve-routed.json")
        );
    }

    /// The protocol's count: 21 pairs, and only `k = 2`, `e = 15` without a
    /// curve, so 20 curves and 40 runs.
    #[test]
    fn the_sweep_has_twenty_curves_and_forty_runs() {
        assert_eq!(pairs().len(), 21);
        let (runs, none) = documents().unwrap();
        assert_eq!(none, vec!["k2e15".to_string()]);
        assert_eq!(runs.len(), 40);
        for run in runs.iter().step_by(2) {
            assert!(run.curve.on_curve(&run.generator));
            assert!(run.r > run.h);
            assert!(replay(run, Some(&run.known_log.to_string())));
            assert!(!replay(
                run,
                Some(&(run.known_log % (run.r - 1) + 1).to_string())
            ));
        }
    }

    #[test]
    fn order_over_extension_agrees_with_counting() {
        // y² + xy = x³ + 1 over GF(2^6) (k = 2, e = 3, a = 0, b = 1): count.
        let field = Field::new(6).unwrap();
        let c = Curve {
            k: field,
            a: 0,
            b: 1,
        };
        let mut count = 1i128;
        for x in 0..64u64 {
            for y in 0..64u64 {
                if c.on_curve(&Some((x, y))) {
                    count += 1;
                }
            }
        }
        let elements = field.subfield(2);
        assert_eq!(elements.len(), 4);
        let t = 5 - count_over_subfield(&c, 2, &elements);
        assert_eq!(order_over_extension(t, 4, 3), count);
    }
}

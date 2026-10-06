//! **Native checker for the PKM tower oracle's rows.**
//!
//! `research/notes/index-calculus/RESEARCH_PKM_TOWER_ORACLE.md` §16.3. It
//! replaces the Python scripts `verify.py` and `analyze.py` of
//! `research/pkm_tower_pilot_20260924/` (`AGENTS.md`, "no Python"). It keeps
//! their contracts and prints what they printed, line for line, and adds the
//! identity checks and the size rule of round 7.
//!
//! - **`verify FILE…`** is the independent cross-check of `AGENTS.md` §6.
//!   Without F4, it recomputes how many tuples of `V^m` solve each finished
//!   tower row's summation system:
//!   - `m = 2`: ordered pairs `(x1, x2)` with `S3(x1, x2, xR) = 0`;
//!   - `m = 3`: ordered triples with `S4(x1, x2, x3, xR) = 0`, where
//!     `S4 = Res_U(S3(x1, x2, U), S3(x3, xR, U))`. This is the chain system
//!     with its free unknown eliminated, so it counts solutions whose partial
//!     sum is defined over `F_{p²}` only, as the `F_p` ideal does;
//!   - `m = 4`: ordered quadruples with
//!     `S5 = Res_U(S4(x1, x2, x3, U), S3(U, x4, xR)) = 0`, as a Sylvester
//!     determinant.
//!
//!   It then checks the one thing that can be checked without trusting F4: a
//!   refuted system has no solution on the grid, a system F4 did not refute
//!   has one, and a planted one has one. `V` is the row's own for isogeny
//!   towers, and rebuilt from `g, zeta` (Kummer) or `u, zeta` (Dickson)
//!   otherwise. The null, naive and ladder controls are skipped, since their
//!   systems are not the summation system.
//! - **`analyze FILE…`** prints, per (engine, kind, `m`, control, `p`):
//!   - one row per `N`;
//!   - the fit `D = α + β·N` of note §5.2 with its bootstrap interval, drawn
//!     from Python's Mersenne Twister as the script drew it;
//!   - the leave-one-`N`-out slopes, the upper-half slope, the final plateau
//!     and A1's reading (§§10.8, 11.6);
//!   - the slope of `log₂` width.
//!
//!   Then it prints the two engines side by side, the confirmation runs, the
//!   generator ladder and the planted-solution count. It exits non-zero if a
//!   repeat disagrees with itself, a planted solution was lost, the engines'
//!   verdicts differ, or a confirmation run contradicts its full run. Each of
//!   those is a bug, not data. Its numbers follow CPython 3.12's semantics:
//!   `sum` with Neumaier compensation, `statistics.median`, true division,
//!   `x ** 2` through C `pow`, and `format(x, "g")`.
//! - **`compare-rows REF NEW`** checks the rows of `NEW` against the rows of
//!   `REF` for the same systems. Every field must be equal but:
//!   - the wall clock;
//!   - `--may-fall` fields, which may only fall;
//!   - `--ignore` fields.
//! - **`compare-trace REF NEW`** checks the `tower step` lines of two logs,
//!   step for step. Every count must be equal: the times and the memory are
//!   the only fields left out.
//! - **`size-cap LOG`** applies note §16.3's rule for stage 2's cap to the
//!   `tower stop in step` line of a stage-1 log.
//!
//! ```bash
//! cargo run --release --example pkm_tower_check -- verify research/pkm_tower_round2_20260925/runs/*.jsonl
//! cargo run --release --example pkm_tower_check -- analyze research/pkm_tower_pilot_20260924/runs/*.jsonl
//! cargo test --release --example pkm_tower_check
//! ```

use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::process::ExitCode;

use serde_json::{json, Value};

// ── Rows ───────────────────────────────────────────────────────────

/// Every JSON object of a JSONL file, blank lines skipped.
fn read_rows(path: &str) -> Vec<Value> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    text.lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .map(|l| serde_json::from_str(l).unwrap_or_else(|e| panic!("{path}: {e}")))
        .collect()
}

fn field<'a>(r: &'a Value, k: &str) -> &'a Value {
    r.get(k).unwrap_or(&Value::Null)
}

fn int(r: &Value, k: &str) -> i64 {
    field(r, k)
        .as_i64()
        .unwrap_or_else(|| panic!("field `{k}` is not an integer"))
}

fn text<'a>(r: &'a Value, k: &str) -> &'a str {
    field(r, k)
        .as_str()
        .unwrap_or_else(|| panic!("field `{k}` is not a string"))
}

/// Python's truth value of a field.
fn truthy(v: &Value) -> bool {
    match v {
        Value::Null => false,
        Value::Bool(b) => *b,
        Value::Number(n) => n.as_f64() != Some(0.0),
        Value::String(s) => !s.is_empty(),
        Value::Array(a) => !a.is_empty(),
        Value::Object(o) => !o.is_empty(),
    }
}

/// A value with every object's keys sorted, for keys
/// (`json.dumps(…, sort_keys=True)`).
fn canonical(v: &Value) -> String {
    fn walk(v: &Value, out: &mut String) {
        match v {
            Value::Object(o) => {
                let mut keys: Vec<&String> = o.keys().collect();
                keys.sort();
                out.push('{');
                for (i, k) in keys.into_iter().enumerate() {
                    if i > 0 {
                        out.push(',');
                    }
                    out.push_str(&Value::String(k.clone()).to_string());
                    out.push(':');
                    walk(&o[k], out);
                }
                out.push('}');
            }
            Value::Array(a) => {
                out.push('[');
                for (i, x) in a.iter().enumerate() {
                    if i > 0 {
                        out.push(',');
                    }
                    walk(x, out);
                }
                out.push(']');
            }
            other => out.push_str(&other.to_string()),
        }
    }
    let mut out = String::new();
    walk(v, &mut out);
    out
}

/// Python's `==` on two JSON values: numbers by value, the rest as they are.
fn py_eq(a: &Value, b: &Value) -> bool {
    match (a, b) {
        (Value::Number(x), Value::Number(y)) => num_of(a)
            .zip(num_of(b))
            .map_or(x == y, |(x, y)| num_cmp(x, y) == Ordering::Equal),
        (Value::Bool(x), Value::Number(_)) | (Value::Number(_), Value::Bool(x)) => {
            let n = if let Value::Number(_) = a { a } else { b };
            num_of(n).is_some_and(|n| num_cmp(n, Num::Int(i64::from(*x))) == Ordering::Equal)
        }
        (Value::Array(x), Value::Array(y)) => {
            x.len() == y.len() && x.iter().zip(y).all(|(p, q)| py_eq(p, q))
        }
        (Value::Object(x), Value::Object(y)) => {
            x.len() == y.len() && x.iter().all(|(k, v)| y.get(k).is_some_and(|w| py_eq(v, w)))
        }
        _ => a == b,
    }
}

/// The engine, and for the signature engine its variant: two variants of
/// one system are two measurements. Round 4's rows predate the steps field
/// and ran the default steps, by signature degree.
fn engine(r: &Value) -> String {
    let e = r.get("engine").and_then(Value::as_str).unwrap_or("f4_fp");
    if e != "f4_fp_tower_sig" {
        return e.to_string();
    }
    let part = |k: &str, default: &str| -> String {
        r.get(k)
            .and_then(Value::as_str)
            .filter(|s| !s.is_empty())
            .unwrap_or(default)
            .to_string()
    };
    format!(
        "{e}/{}/{}/{}",
        part("sig_order", "PositionFirst"),
        part("sig_rewrite", "Ratio"),
        part("sig_steps", "SignatureDegree")
    )
}

/// The system itself, whichever engine measured it.
fn instance_key(r: &Value) -> String {
    let curve = field(r, "curve");
    canonical(&json!([
        field(r, "p"),
        field(r, "kind"),
        field(r, "m"),
        field(r, "control"),
        field(r, "t"),
        field(r, "g"),
        field(r, "target"),
        field(r, "target_index"),
        field(r, "x_r"),
        field(curve, "a"),
        field(curve, "b"),
        field(r, "tower"),
    ]))
}

/// One measurement: a system, the engine that measured it, and the degree
/// bound it ran under.
fn system_key(r: &Value) -> String {
    format!(
        "{}|{}|{}",
        instance_key(r),
        engine(r),
        field(r, "max_degree_bound")
    )
}

/// Python's `repr` of a string.
fn py_repr(s: &str) -> String {
    let quote = if s.contains('\'') && !s.contains('"') {
        '"'
    } else {
        '\''
    };
    let mut out = String::new();
    out.push(quote);
    for c in s.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\t' => out.push_str("\\t"),
            c if c == quote => {
                out.push('\\');
                out.push(c);
            }
            c => out.push(c),
        }
    }
    out.push(quote);
    out
}

fn py_bool(b: bool) -> &'static str {
    if b {
        "True"
    } else {
        "False"
    }
}

// ── Python's numbers ─────────────────────────────────────────────────

/// A Python `int` or `float`, as the analysis meets them: JSON integers,
/// JSON floats, and what `statistics.median` returns.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Num {
    Int(i64),
    Float(f64),
}

impl Num {
    fn f(self) -> f64 {
        match self {
            Num::Int(i) => i as f64,
            Num::Float(x) => x,
        }
    }
}

fn num_of(v: &Value) -> Option<Num> {
    if let Some(i) = v.as_i64() {
        Some(Num::Int(i))
    } else if let Some(u) = v.as_u64() {
        Some(Num::Int(i64::try_from(u).expect("an integer below 2^63")))
    } else {
        v.as_f64().map(Num::Float)
    }
}

fn num(r: &Value, k: &str) -> Num {
    num_of(field(r, k)).unwrap_or_else(|| panic!("field `{k}` is not a number"))
}

/// Python's comparison of two numbers. The integers here are far below
/// 2^53, where a float holds them exactly.
fn num_cmp(a: Num, b: Num) -> Ordering {
    match (a, b) {
        (Num::Int(x), Num::Int(y)) => x.cmp(&y),
        _ => a.f().partial_cmp(&b.f()).expect("no NaN"),
    }
}

fn py_add(a: Num, b: Num) -> Num {
    match (a, b) {
        (Num::Int(x), Num::Int(y)) => Num::Int(x + y),
        _ => Num::Float(a.f() + b.f()),
    }
}

/// `a / n`, Python's true division. For an integer `a` it is correctly
/// rounded, as Python's is: both operands are exact floats.
fn py_div(a: Num, n: usize) -> f64 {
    a.f() / n as f64
}

/// CPython 3.12's `sum` with the default start `0`: integers exactly until
/// the first float, then floats with Neumaier's compensation (integers
/// after it added plainly), and the compensation added at the end.
fn py_sum(items: impl IntoIterator<Item = Num>) -> Num {
    let mut it = items.into_iter();
    let mut i_result: i64 = 0;
    let mut f_result = loop {
        match it.next() {
            None => return Num::Int(i_result),
            Some(Num::Int(b)) => i_result += b,
            Some(Num::Float(x)) => break i_result as f64 + x,
        }
    };
    let mut c = 0.0f64;
    for item in it {
        match item {
            Num::Float(x) => {
                let t = f_result + x;
                if f_result.abs() >= x.abs() {
                    c += (f_result - t) + x;
                } else {
                    c += (x - t) + f_result;
                }
                f_result = t;
            }
            Num::Int(v) => f_result += v as f64,
        }
    }
    if c != 0.0 && c.is_finite() {
        f_result += c;
    }
    Num::Float(f_result)
}

extern "C" {
    /// The C library's `pow`, which CPython's `float ** int` calls.
    fn pow(x: f64, y: f64) -> f64;
}

/// `x ** 2` for a Python float: C `pow`, which the compiler must not turn
/// into `x * x`.
fn py_square(x: f64) -> f64 {
    // SAFETY: `pow` is a pure function of two doubles.
    unsafe { pow(x, std::hint::black_box(2.0)) }
}

/// `statistics.median`: the middle value as it is for an odd count, else
/// the mean of the middle two, `(a + b) / 2`, a float.
fn median(vals: impl IntoIterator<Item = Num>) -> Num {
    let mut v: Vec<Num> = vals.into_iter().collect();
    assert!(!v.is_empty(), "median of no data");
    v.sort_by(|a, b| num_cmp(*a, *b));
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        Num::Float(py_div(py_add(v[n / 2 - 1], v[n / 2]), 2))
    }
}

/// `format(x, "g")`: six significant digits, fixed notation for decimal
/// exponents from −4 to 5, trailing zeros removed.
fn fmt_g(x: Num) -> String {
    let x = x.f();
    if x == 0.0 {
        return if x.is_sign_negative() { "-0" } else { "0" }.to_string();
    }
    let strip = |s: String| -> String {
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            s
        }
    };
    let sci = format!("{x:.5e}");
    let (mantissa, exponent) = sci.split_once('e').expect("an exponent");
    let exponent: i32 = exponent.parse().expect("an exponent");
    if !(-4..6).contains(&exponent) {
        let sign = if exponent < 0 { '-' } else { '+' };
        format!(
            "{}e{sign}{:02}",
            strip(mantissa.to_string()),
            exponent.unsigned_abs()
        )
    } else {
        strip(format!("{x:.*}", (5 - exponent) as usize))
    }
}

/// CPython's `random.Random`: MT19937 seeded by `init_by_array` with the
/// seed's 32-bit words, and `choice` by rejection on `getrandbits`.
struct PyRandom {
    mt: [u32; 624],
    index: usize,
}

impl PyRandom {
    fn new(seed: u64) -> Self {
        let mut key: Vec<u32> = Vec::new();
        let mut s = seed;
        while s > 0 {
            key.push(s as u32);
            s >>= 32;
        }
        if key.is_empty() {
            key.push(0);
        }
        let mut r = PyRandom {
            mt: [0; 624],
            index: 624,
        };
        r.init_genrand(19_650_218);
        let n = 624usize;
        let (mut i, mut j) = (1usize, 0usize);
        for _ in 0..n.max(key.len()) {
            let prev = r.mt[i - 1] ^ (r.mt[i - 1] >> 30);
            r.mt[i] = (r.mt[i] ^ prev.wrapping_mul(1_664_525))
                .wrapping_add(key[j])
                .wrapping_add(j as u32);
            i += 1;
            j += 1;
            if i >= n {
                r.mt[0] = r.mt[n - 1];
                i = 1;
            }
            if j >= key.len() {
                j = 0;
            }
        }
        for _ in 0..n - 1 {
            let prev = r.mt[i - 1] ^ (r.mt[i - 1] >> 30);
            r.mt[i] = (r.mt[i] ^ prev.wrapping_mul(1_566_083_941)).wrapping_sub(i as u32);
            i += 1;
            if i >= n {
                r.mt[0] = r.mt[n - 1];
                i = 1;
            }
        }
        r.mt[0] = 0x8000_0000;
        r.index = n;
        r
    }

    fn init_genrand(&mut self, s: u32) {
        self.mt[0] = s;
        for i in 1..624 {
            let prev = self.mt[i - 1] ^ (self.mt[i - 1] >> 30);
            self.mt[i] = prev.wrapping_mul(1_812_433_253).wrapping_add(i as u32);
        }
        self.index = 624;
    }

    fn genrand_u32(&mut self) -> u32 {
        const UPPER: u32 = 0x8000_0000;
        const LOWER: u32 = 0x7fff_ffff;
        const MATRIX_A: u32 = 0x9908_b0df;
        if self.index >= 624 {
            for k in 0..624 {
                let y = (self.mt[k] & UPPER) | (self.mt[(k + 1) % 624] & LOWER);
                let mag = if y & 1 == 1 { MATRIX_A } else { 0 };
                self.mt[k] = self.mt[(k + 397) % 624] ^ (y >> 1) ^ mag;
            }
            self.index = 0;
        }
        let mut y = self.mt[self.index];
        self.index += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^= y >> 18;
        y
    }

    /// `random.random()`: 53 bits from two draws.
    #[cfg(test)]
    fn random(&mut self) -> f64 {
        let a = f64::from(self.genrand_u32() >> 5);
        let b = f64::from(self.genrand_u32() >> 6);
        (a * 67_108_864.0 + b) * (1.0 / 9_007_199_254_740_992.0)
    }

    /// `getrandbits(k)` for `1 ≤ k ≤ 32`.
    fn getrandbits(&mut self, k: u32) -> u64 {
        assert!((1..=32).contains(&k));
        u64::from(self.genrand_u32() >> (32 - k))
    }

    /// `_randbelow(n)`: `k = n.bit_length()` bits, redrawn while `≥ n`.
    fn randbelow(&mut self, n: usize) -> usize {
        assert!(n > 0);
        let k = usize::BITS - n.leading_zeros();
        loop {
            let r = self.getrandbits(k) as usize;
            if r < n {
                return r;
            }
        }
    }
}

// ── verify ─────────────────────────────────────────────────────────

/// Arithmetic modulo a prime below 2^32.
#[derive(Clone, Copy)]
struct Fp(u64);

impl Fp {
    fn add(self, a: u64, b: u64) -> u64 {
        (a + b) % self.0
    }
    fn sub(self, a: u64, b: u64) -> u64 {
        (a + self.0 - b % self.0) % self.0
    }
    fn mul(self, a: u64, b: u64) -> u64 {
        (u128::from(a) * u128::from(b) % u128::from(self.0)) as u64
    }
    fn neg(self, a: u64) -> u64 {
        self.sub(0, a)
    }
    fn pow(self, mut a: u64, mut e: u64) -> u64 {
        let mut r = 1 % self.0;
        a %= self.0;
        while e > 0 {
            if e & 1 == 1 {
                r = self.mul(r, a);
            }
            a = self.mul(a, a);
            e >>= 1;
        }
        r
    }
    fn inv(self, a: u64) -> u64 {
        self.pow(a, self.0 - 2)
    }
}

/// `S3(x1, x2, x3)` on `y² = x³ + a x + b`.
fn s3(f: Fp, x1: u64, x2: u64, x3: u64, a: u64, b: u64) -> u64 {
    let (aa, bb, cc) = s3_coeffs(f, x1, x2, a, b);
    f.add(f.add(f.mul(f.mul(aa, x3), x3), f.mul(bb, x3)), cc)
}

/// `S3(x1, x2, U)` as `A U² + B U + C`.
fn s3_coeffs(f: Fp, x1: u64, x2: u64, a: u64, b: u64) -> (u64, u64, u64) {
    let d = f.sub(x1, x2);
    let s = f.add(x1, x2);
    let pr = f.mul(x1, x2);
    let aa = f.mul(d, d);
    let bb = f.neg(f.mul(2, f.add(f.mul(s, f.add(pr, a)), f.mul(2, b))));
    let q = f.sub(pr, a);
    let cc = f.sub(f.mul(q, q), f.mul(f.mul(4, b), s));
    (aa, bb, cc)
}

/// The resultant of two quadratics `A U² + B U + C`.
fn res2(f: Fp, g: (u64, u64, u64), h: (u64, u64, u64)) -> u64 {
    let (a, b, c) = g;
    let (a2, b2, c2) = h;
    let u = f.sub(f.mul(a, c2), f.mul(a2, c));
    let v = f.sub(f.mul(a, b2), f.mul(a2, b));
    let w = f.sub(f.mul(b, c2), f.mul(b2, c));
    f.sub(f.mul(u, u), f.mul(v, w))
}

/// Product of two polynomials in `U` (coefficients, low to high).
fn pmul(f: Fp, g: &[u64], h: &[u64]) -> Vec<u64> {
    let mut out = vec![0; g.len() + h.len() - 1];
    for (i, &x) in g.iter().enumerate() {
        for (j, &y) in h.iter().enumerate() {
            out[i + j] = f.add(out[i + j], f.mul(x, y));
        }
    }
    out
}

fn psub(f: Fp, g: &[u64], h: &[u64]) -> Vec<u64> {
    (0..g.len().max(h.len()))
        .map(|i| {
            f.sub(
                g.get(i).copied().unwrap_or(0),
                h.get(i).copied().unwrap_or(0),
            )
        })
        .collect()
}

fn scale(f: Fp, c: u64, g: &[u64]) -> Vec<u64> {
    g.iter().map(|&x| f.mul(c, x)).collect()
}

/// `S3(x3, U, w)` as a quadratic in `w`, its coefficients polynomials in `U`
/// (low to high).
fn s3_in_w(f: Fp, x3: u64, a: u64, b: u64) -> [Vec<u64>; 3] {
    let m2 = f.neg(2);
    let ax_2b = f.add(f.mul(a, x3), f.mul(2, b));
    let a2 = vec![f.mul(x3, x3), f.mul(m2, x3), 1];
    let b2 = vec![
        f.mul(m2, ax_2b),
        f.mul(m2, f.add(f.mul(x3, x3), a)),
        f.mul(m2, x3),
    ];
    let c2 = vec![
        f.sub(f.mul(a, a), f.mul(f.mul(4, b), x3)),
        f.mul(m2, ax_2b),
        f.mul(x3, x3),
    ];
    [a2, b2, c2]
}

/// `S4(x1, x2, x3, U) = Res_w(S3(x1, x2, w), S3(w, x3, U))`, a polynomial in
/// `U`: the two-quadratic resultant with polynomial coefficients.
fn s4_in_u(f: Fp, x1: u64, x2: u64, x3: u64, a: u64, b: u64) -> Vec<u64> {
    let (aa, bb, cc) = s3_coeffs(f, x1, x2, a, b);
    let [a2, b2, c2] = s3_in_w(f, x3, a, b);
    let u = psub(f, &scale(f, aa, &c2), &scale(f, cc, &a2));
    let v = psub(f, &scale(f, aa, &b2), &scale(f, bb, &a2));
    let w = psub(f, &scale(f, bb, &c2), &scale(f, cc, &b2));
    psub(f, &pmul(f, &u, &u), &pmul(f, &v, &w))
}

/// Whether `Res(g, h) = 0` for polynomials in `U`, as the Sylvester
/// determinant at their formal degrees.
fn res_zero(f: Fp, g: &[u64], h: &[u64]) -> bool {
    let (dg, dh) = (g.len() - 1, h.len() - 1);
    let n = dg + dh;
    let mut rows: Vec<Vec<u64>> = Vec::with_capacity(n);
    for i in 0..dh {
        let mut row = vec![0; n];
        for (k, &c) in g.iter().rev().enumerate() {
            row[i + k] = c;
        }
        rows.push(row);
    }
    for i in 0..dg {
        let mut row = vec![0; n];
        for (k, &c) in h.iter().rev().enumerate() {
            row[i + k] = c;
        }
        rows.push(row);
    }
    // The determinant is zero iff some column has no pivot.
    for c in 0..n {
        let Some(piv) = (c..n).find(|&r| rows[r][c] != 0) else {
            return true;
        };
        rows.swap(c, piv);
        let inv = f.inv(rows[c][c]);
        for r in c + 1..n {
            if rows[r][c] != 0 {
                let k = f.mul(rows[r][c], inv);
                let pivot = rows[c].clone();
                for (x, &y) in rows[r].iter_mut().zip(&pivot) {
                    *x = f.sub(*x, f.mul(k, y));
                }
            }
        }
    }
    false
}

/// The tower's `V`: the row's own, or rebuilt from `g, zeta` (Kummer) or
/// `u, zeta` (Dickson).
fn tower_v(row: &Value) -> Option<Vec<u64>> {
    let info = field(row, "tower");
    let p = int(row, "p") as u64;
    let f = Fp(p);
    let t = int(row, "t") as u32;
    if let Some(v) = info.get("v").and_then(Value::as_array) {
        return Some(v.iter().map(|x| x.as_u64().expect("V")).collect());
    }
    let get = |k: &str| info.get(k).and_then(Value::as_u64).expect("a tower field");
    match text(row, "kind") {
        "kummer" => {
            let (g, z) = (get("g"), get("zeta"));
            Some((0..1u64 << t).map(|k| f.mul(g, f.pow(z, k))).collect())
        }
        "dickson" => {
            let (u, z) = (get("u"), get("zeta"));
            Some(
                (0..1u64 << t)
                    .map(|k| {
                        let w = f.mul(u, f.pow(z, k));
                        f.add(w, f.pow(w, p - 2))
                    })
                    .collect(),
            )
        }
        _ => None,
    }
}

/// How many tuples of `V^m` solve the row's summation system.
fn count(row: &Value, v: &[u64]) -> Option<u64> {
    let f = Fp(int(row, "p") as u64);
    let curve = field(row, "curve");
    let (a, b) = (
        curve["a"].as_u64().expect("a"),
        curve["b"].as_u64().expect("b"),
    );
    let xr = int(row, "x_r") as u64;
    match int(row, "m") {
        2 => Some(
            v.iter()
                .flat_map(|&x1| v.iter().map(move |&x2| (x1, x2)))
                .filter(|&(x1, x2)| s3(f, x1, x2, xr, a, b) == 0)
                .count() as u64,
        ),
        3 => {
            let tails: Vec<(u64, u64, u64)> =
                v.iter().map(|&x3| s3_coeffs(f, x3, xr, a, b)).collect();
            let mut n = 0;
            for &x1 in v {
                for &x2 in v {
                    let g = s3_coeffs(f, x1, x2, a, b);
                    n += tails.iter().filter(|&&h| res2(f, g, h) == 0).count() as u64;
                }
            }
            Some(n)
        }
        4 => {
            // S5 = Res_U(S4(x1, x2, x3, U), S3(U, x4, x_R)), and
            // S3(U, x4, x_R) = S3(x4, x_R, U).
            let tails: Vec<Vec<u64>> = v
                .iter()
                .map(|&x4| {
                    let (aa, bb, cc) = s3_coeffs(f, x4, xr, a, b);
                    vec![cc, bb, aa]
                })
                .collect();
            let mut n = 0;
            for &x1 in v {
                for &x2 in v {
                    for &x3 in v {
                        let s4 = s4_in_u(f, x1, x2, x3, a, b);
                        n += tails.iter().filter(|h| res_zero(f, &s4, h)).count() as u64;
                    }
                }
            }
            Some(n)
        }
        _ => None,
    }
}

fn verify(paths: &[String]) -> ExitCode {
    let (mut checked, mut agreed, mut skipped, mut bounded) = (0, 0, 0, 0);
    let mut bad: Vec<String> = Vec::new();
    let mut seen: HashSet<String> = HashSet::new();
    for path in paths {
        for row in read_rows(path) {
            if row.get("summary").is_some()
                || text(&row, "control") != "tower"
                || truthy(field(&row, "timed_out"))
            {
                continue;
            }
            // A run that stopped at its degree bound with pairs left (a
            // confirmation run, note §11.4) claims nothing about solutions:
            // it neither refuted nor finished.
            let pairs_above = row
                .get("pairs_above_bound")
                .and_then(Value::as_i64)
                .unwrap_or(0);
            let inconsistent = truthy(field(&row, "inconsistent"));
            if pairs_above > 0 && !inconsistent {
                bounded += 1;
                continue;
            }
            // A system measured by two runs is checked once.
            let curve = field(&row, "curve");
            let key = canonical(&json!([
                field(&row, "p"),
                field(&row, "kind"),
                field(&row, "m"),
                field(&row, "t"),
                field(&row, "target"),
                field(&row, "target_index"),
                field(&row, "x_r"),
                field(curve, "a"),
                field(curve, "b"),
                field(&row, "tower"),
                engine(&row),
            ]));
            if !seen.insert(key) {
                continue;
            }
            let Some(v) = tower_v(&row) else {
                skipped += 1;
                continue;
            };
            let n = count(&row, &v).expect("m is 2, 3 or 4");
            checked += 1;
            let mut ok = (n == 0) == inconsistent;
            if text(&row, "target") == "planted" {
                ok = ok && n >= 1;
            }
            if ok {
                agreed += 1;
            } else {
                bad.push(format!(
                    "({}, {}, {}, {}, {}, {}, {n}, {})",
                    py_repr(path),
                    py_repr(&engine(&row)),
                    py_repr(text(&row, "kind")),
                    int(&row, "m"),
                    int(&row, "N"),
                    py_repr(text(&row, "target")),
                    py_bool(inconsistent)
                ));
            }
        }
    }
    println!(
        "checked {checked} tower rows against exhaustive search; {agreed} agree; {skipped} skipped (no V); {bounded} stopped at their degree bound, with no verdict to check."
    );
    for b in &bad {
        println!("DISAGREE {b}");
    }
    if bad.is_empty() {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}

// ── analyze ────────────────────────────────────────────────────────

/// Fields that depend only on the system and the algorithm, not on the
/// machine.
const DETERMINISTIC: [&str; 7] = [
    "solving_degree_max",
    "last_productive_degree",
    "degree_reached",
    "max_cols_to_solution",
    "steps_to_solution",
    "inconsistent",
    "basis_len",
];
/// The fields a run halted by the staircase stop shares with a full run of
/// the same system.
const STOP_COMPARABLE: [&str; 2] = ["solving_degree_max", "max_cols_to_solution"];

fn timed_out(r: &Value) -> bool {
    truthy(field(r, "timed_out"))
}

/// The example's degree bound without `--cap`: `n + d + 6`, or `2d + 8` for
/// the naive control.
fn default_bound(r: &Value) -> i64 {
    let d = field(r, "input_degrees")
        .as_array()
        .expect("input_degrees")
        .iter()
        .map(|x| x.as_i64().expect("a degree"))
        .max()
        .expect("an input degree");
    if text(r, "control") == "naive" {
        2 * d + 8
    } else {
        int(r, "n_vars") + d + 6
    }
}

/// A run with a degree bound below the default: a confirmation run (note
/// §11.4), which is expected to stop at the bound and is not a measurement
/// of `D`.
fn capped(r: &Value) -> bool {
    field(r, "max_degree_bound")
        .as_i64()
        .is_some_and(|b| b < default_bound(r))
}

/// Every row, a system measured twice counted once (see `analyze`'s doc).
fn load(paths: &[String]) -> (Vec<Value>, usize) {
    let mut rows: Vec<Value> = Vec::new();
    let mut seen: HashMap<String, usize> = HashMap::new();
    let (mut duplicates, mut old_format) = (0, 0);
    let mut mismatches: Vec<String> = Vec::new();
    for path in paths {
        for r in read_rows(path) {
            if r.get("summary").is_some() {
                continue;
            }
            if r.get("solving_degree_max").is_none() {
                old_format += 1;
                continue;
            }
            let k = system_key(&r);
            let Some(&i) = seen.get(&k) else {
                seen.insert(k, rows.len());
                rows.push(r);
                continue;
            };
            duplicates += 1;
            let first = &rows[i];
            let mismatch = |diff: &[&str]| {
                format!(
                    "({}, {}, {}, {}, {}, [{}])",
                    py_repr(path),
                    py_repr(text(&r, "kind")),
                    int(&r, "m"),
                    py_repr(text(&r, "control")),
                    int(&r, "N"),
                    diff.iter()
                        .map(|d| py_repr(d))
                        .collect::<Vec<_>>()
                        .join(", ")
                )
            };
            let (cut_r, cut_first) = (timed_out(&r), timed_out(first));
            if !(cut_r || cut_first) {
                // A copy halted by the staircase stop ends before the basis
                // is certified: it shares only the degree and the width.
                let stopped = !field(&r, "staircase_at_stop").is_null()
                    || !field(first, "staircase_at_stop").is_null();
                let fields: &[&str] = if stopped {
                    &STOP_COMPARABLE
                } else {
                    &DETERMINISTIC
                };
                let diff: Vec<&str> = fields
                    .iter()
                    .copied()
                    .filter(|f| !py_eq(field(&r, f), field(first, f)))
                    .collect();
                if !diff.is_empty() {
                    mismatches.push(mismatch(&diff));
                }
            } else if cut_r != cut_first {
                // One copy finished and one stopped early, whose degree is
                // only a lower bound: the finished copy stands for the
                // system, and its degree must not fall below that bound.
                let (done, cut) = if cut_r { (first, &r) } else { (&r, first) };
                if int(done, "solving_degree_max") < int(cut, "solving_degree_max") {
                    mismatches.push(mismatch(&[
                        "solving_degree_max below a stopped copy's lower bound",
                    ]));
                }
                if !cut_r {
                    rows[i] = r;
                }
            }
        }
    }
    println!(
        "{} distinct systems; {duplicates} repeated measurements ({} disagreeing on a deterministic field); {old_format} rows in the format before `solving_degree_max`, skipped.",
        rows.len(),
        mismatches.len()
    );
    for m in &mismatches {
        println!("MISMATCH {m}");
    }
    (rows, mismatches.len())
}

/// The least-squares line through `(xs, ys)`, as `(intercept, slope)`.
fn lsq(xs: &[Num], ys: &[Num]) -> Option<(f64, f64)> {
    let n = xs.len();
    let mx = py_div(py_sum(xs.iter().copied()), n);
    let my = py_div(py_sum(ys.iter().copied()), n);
    let sxx = py_sum(xs.iter().map(|x| Num::Float(py_square(x.f() - mx)))).f();
    if sxx == 0.0 {
        return None;
    }
    let sxy = py_sum(
        xs.iter()
            .zip(ys)
            .map(|(x, y)| Num::Float((x.f() - mx) * (y.f() - my))),
    )
    .f();
    let b = sxy / sxx;
    Some((my - b * mx, b))
}

/// The rows of one cell, finished, by `N` (ascending, rows in file order).
type ByN<'a> = BTreeMap<i64, Vec<&'a Value>>;

/// Resample targets within each `N`; refit the slope each time.
fn bootstrap_slope(by_n: &ByN, key: impl Fn(&Value) -> Num) -> Option<(f64, f64)> {
    let mut rng = PyRandom::new(20_260_924);
    if by_n.len() < 2 {
        return None;
    }
    let mut slopes: Vec<f64> = Vec::new();
    for _ in 0..4000 {
        let (mut xs, mut ys) = (Vec::new(), Vec::new());
        for (&n, vals) in by_n {
            for _ in 0..vals.len() {
                xs.push(Num::Int(n));
                ys.push(key(vals[rng.randbelow(vals.len())]));
            }
        }
        if let Some((_, b)) = lsq(&xs, &ys) {
            slopes.push(b);
        }
    }
    slopes.sort_by(|a, b| a.partial_cmp(b).expect("no NaN"));
    let len = slopes.len() as f64;
    let lo = slopes[(0.025 * len) as usize];
    let hi = slopes[(0.975 * len) as usize - 1];
    Some((lo, hi))
}

fn medians(by_n: &ByN, key: &impl Fn(&Value) -> Num) -> Vec<(i64, Num)> {
    by_n.iter()
        .map(|(&n, rs)| (n, median(rs.iter().map(|r| key(r)))))
        .collect()
}

/// Least-squares slopes over the per-`N` medians with each `N` left out in
/// turn: the pre-registered bootstrap collapses when every target at an `N`
/// gives the same degree.
fn jackknife_slopes(by_n: &ByN, key: impl Fn(&Value) -> Num) -> Option<(f64, f64)> {
    if by_n.len() < 3 {
        return None;
    }
    let med = medians(by_n, &key);
    let out: Vec<f64> = med
        .iter()
        .filter_map(|&(drop, _)| {
            let keep: Vec<&(i64, Num)> = med.iter().filter(|(n, _)| *n != drop).collect();
            let xs: Vec<Num> = keep.iter().map(|(n, _)| Num::Int(*n)).collect();
            let ys: Vec<Num> = keep.iter().map(|(_, y)| *y).collect();
            lsq(&xs, &ys).map(|fit| fit.1)
        })
        .collect();
    // Python's `min` and `max`: the first of equal values stays.
    let (mut min, mut max) = (out[0], out[0]);
    for &x in &out[1..] {
        if x < min {
            min = x;
        }
        if x > max {
            max = x;
        }
    }
    Some((min, max))
}

fn upper_half_slope(by_n: &ByN, key: impl Fn(&Value) -> Num) -> Option<(Option<f64>, i64, i64)> {
    if by_n.len() < 4 {
        return None;
    }
    let med = medians(by_n, &key);
    let upper = &med[med.len() / 2..];
    let xs: Vec<Num> = upper.iter().map(|(n, _)| Num::Int(*n)).collect();
    let ys: Vec<Num> = upper.iter().map(|(_, y)| *y).collect();
    Some((
        lsq(&xs, &ys).map(|fit| fit.1),
        upper[0].0,
        upper[upper.len() - 1].0,
    ))
}

/// The `N` range since `D` last rose, and `D` there (note §10.8).
fn final_plateau(by_n: &ByN, key: impl Fn(&Value) -> Num) -> (i64, i64, Num) {
    let med = medians(by_n, &key);
    let mut start = med[0].0;
    for i in 1..med.len() {
        if num_cmp(med[i].1, med[i - 1].1) == Ordering::Greater {
            start = med[i].0;
        }
    }
    let last = med[med.len() - 1];
    (start, last.0, last.1)
}

/// Amendment A1 of the note (§10.8), as §11.6 adopted it: H1a if the final
/// plateau is longer than 10 in `N`; H0 if the least-squares slope over the
/// upper half of the range is at least 1/4; inconclusive otherwise.
fn a1_reading(by_n: &ByN, key: impl Fn(&Value) -> Num + Copy) -> String {
    let (lo_n, hi_n, _) = final_plateau(by_n, key);
    if hi_n - lo_n > 10 {
        return format!("H1a (final plateau L = {} > 10)", hi_n - lo_n);
    }
    match upper_half_slope(by_n, key) {
        Some((Some(slope), first, last)) => {
            let reading = if slope >= 0.25 { "H0" } else { "inconclusive" };
            format!(
                "{reading} (upper-half slope {slope:.3} over N = {first}…{last}, final plateau L = {})",
                hi_n - lo_n
            )
        }
        _ => format!("inconclusive (only {} values of N)", by_n.len()),
    }
}

/// The decision rule of note §5.2, on `β` alone.
fn verdict(lo: f64, hi: f64) -> &'static str {
    if lo > 0.25 {
        "H0 (beta lower bound > 0.25)"
    } else if hi < 0.10 {
        "H1 candidate (beta upper bound < 0.10; also needs D below the null)"
    } else {
        "inconclusive"
    }
}

/// Systems both engines finished: do their solving degrees agree? Returns
/// the number on which the verdicts (refuted or not) disagree: a bug in one
/// engine, where a different degree is data.
fn compare_engines(rows: &[&Value]) -> usize {
    let mut order: Vec<String> = Vec::new();
    let mut by_instance: HashMap<String, HashMap<String, &Value>> = HashMap::new();
    for r in rows.iter().copied().filter(|r| !timed_out(r)) {
        let k = instance_key(r);
        if !by_instance.contains_key(&k) {
            order.push(k.clone());
        }
        by_instance.entry(k).or_default().insert(engine(r), r);
    }
    let pairs: Vec<(&Value, &Value)> = order
        .iter()
        .filter_map(|k| {
            let v = &by_instance[k];
            Some((*v.get("f4_fp")?, *v.get("f4_fp_tower")?))
        })
        .collect();
    if pairs.is_empty() {
        return 0;
    }
    let d = |r: &Value| int(r, "solving_degree_max");
    let refuted = |r: &Value| truthy(field(r, "inconsistent"));
    let agree = pairs.iter().filter(|(a, b)| d(a) == d(b)).count();
    let verdicts = pairs
        .iter()
        .filter(|(a, b)| refuted(a) == refuted(b))
        .count();
    println!("\n### The two engines on the same systems\n");
    println!(
        "{} systems finished by both: the solving degree agrees on {agree}; the verdict (refuted or not) agrees on {verdicts}.",
        pairs.len()
    );
    let mut diff: BTreeMap<(String, i64, String, i64, String, i64, i64), usize> = BTreeMap::new();
    for (a, b) in &pairs {
        if d(a) != d(b) {
            *diff
                .entry((
                    text(a, "kind").to_string(),
                    int(a, "m"),
                    text(a, "control").to_string(),
                    int(a, "N"),
                    text(a, "target").to_string(),
                    d(a),
                    d(b),
                ))
                .or_default() += 1;
        }
    }
    if !diff.is_empty() {
        println!("\n| kind | m | control | N | target | D f4_fp | D f4_fp_tower | systems |");
        println!("|:--|--:|:--|--:|:--|--:|--:|--:|");
        for ((kind, m, control, n, target, d1, d2), c) in &diff {
            println!("| {kind} | {m} | {control} | {n} | {target} | {d1} | {d2} | {c} |");
        }
    }
    let split: Vec<&(&Value, &Value)> = pairs
        .iter()
        .filter(|(a, b)| refuted(a) != refuted(b))
        .collect();
    for (a, b) in &split {
        println!(
            "VERDICT MISMATCH {} {} {} {} {} {} f4_fp: {} f4_fp_tower: {}",
            text(a, "kind"),
            int(a, "m"),
            text(a, "control"),
            int(a, "N"),
            text(a, "target"),
            int(a, "target_index"),
            py_bool(refuted(a)),
            py_bool(refuted(b))
        );
    }
    split.len()
}

/// The confirmation runs of note §11.4 beside the full runs of the same
/// systems. Returns how many contradict their full run: a capped run that
/// refutes although the full run needed a step above the bound.
fn confirmations(capped_rows: &[&Value], rows: &[&Value]) -> usize {
    if capped_rows.is_empty() {
        return 0;
    }
    let full: HashMap<(String, String), &Value> = rows
        .iter()
        .map(|r| ((instance_key(r), engine(r)), *r))
        .collect();
    println!("\n### Confirmation runs (degree bound below the default, section 11.4)\n");
    println!("| kind | m | control | p | N | target | bound | D, full run | pairs above the bound | refuted | confirms |");
    println!("|:--|--:|:--|--:|--:|:--|--:|--:|--:|:--|:--|");
    let mut sorted: Vec<&Value> = capped_rows.to_vec();
    sorted.sort_by(|a, b| {
        let k = |r: &Value| {
            (
                text(r, "kind").to_string(),
                int(r, "m"),
                text(r, "control").to_string(),
                int(r, "p"),
                int(r, "N"),
                int(r, "target_index"),
            )
        };
        k(a).cmp(&k(b))
    });
    let mut bad = 0;
    for r in sorted {
        let f = full.get(&(instance_key(r), engine(r)));
        let d_full = f
            .filter(|f| !timed_out(f))
            .map(|f| int(f, "solving_degree_max"));
        let bound = int(r, "max_degree_bound");
        let inconsistent = truthy(field(r, "inconsistent"));
        let status = if timed_out(r) {
            "no (timed out)".to_string()
        } else if inconsistent {
            let mut s = "no: refuted under the bound".to_string();
            if d_full.is_some_and(|d| d > bound) {
                bad += 1;
                s.push_str(" (CONTRADICTS the full run)");
            }
            s
        } else if !field(r, "staircase_at_stop").is_null() {
            "no (staircase stop)".to_string()
        } else if int(r, "pairs_above_bound") > 0 {
            "yes".to_string()
        } else {
            "no (finished below the bound)".to_string()
        };
        println!(
            "| {} | {} | {} | {} | {} | {} {} | {bound} | {} | {} | {} | {status} |",
            text(r, "kind"),
            int(r, "m"),
            text(r, "control"),
            int(r, "p"),
            int(r, "N"),
            text(r, "target"),
            int(r, "target_index"),
            d_full.map_or("—".to_string(), |d| d.to_string()),
            int(r, "pairs_above_bound"),
            if inconsistent { "yes" } else { "no" }
        );
    }
    bad
}

fn degree(r: &Value) -> Num {
    num(r, "solving_degree_max")
}

fn log2_width(r: &Value) -> Num {
    Num::Float((int(r, "max_cols_to_solution").max(1) as f64).log2())
}

/// Tables and fits per cell, and the cross-checks, as `analyze.py`.
///
/// Files may overlap. Every cell draws its tower, curve and targets from a
/// seed fixed by (seed, kind, `m`, `t`, `g`), so a cell that two runs share
/// is the same system measured twice. Such a system counts once, as the copy
/// read first, unless that copy timed out and a later one finished: the
/// finished copy then stands for the system, and its degree must not fall
/// below the lower bound of the copy that stopped. Finished copies must
/// agree on every deterministic field. Rows written before
/// `solving_degree_max` existed are skipped and counted. Runs that hit the
/// budget are reported and kept out of the fits: their degree is a lower
/// bound, never a value.
fn analyze(paths: &[String]) -> ExitCode {
    let (all, mismatches) = load(paths);
    let capped_rows: Vec<&Value> = all.iter().filter(|r| capped(r)).collect();
    let rows: Vec<&Value> = all.iter().filter(|r| !capped(r)).collect();
    let mut groups: BTreeMap<(String, String, i64, String, i64), Vec<&Value>> = BTreeMap::new();
    for r in rows.iter().copied() {
        if text(r, "control") == "ladder" {
            continue;
        }
        groups
            .entry((
                engine(r),
                text(r, "kind").to_string(),
                int(r, "m"),
                text(r, "control").to_string(),
                int(r, "p"),
            ))
            .or_default()
            .push(r);
    }
    for ((eng, kind, m, control, prime), rs) in &groups {
        println!("\n### {kind}, m = {m}, {control}, p = {prime}, {eng}\n");
        println!("| N | runs | timeouts | D (median, range) | width to solution (median) | F4 steps to solution (median) | F4 ms (median) |");
        println!("|--:|--:|--:|:--|--:|--:|--:|");
        let mut by_n: ByN = BTreeMap::new();
        for r in rs {
            by_n.entry(int(r, "N")).or_default().push(r);
        }
        let mut fit_by_n: ByN = BTreeMap::new();
        for (&n, vals) in &by_n {
            let done: Vec<&Value> = vals.iter().copied().filter(|r| !timed_out(r)).collect();
            let to = vals.len() - done.len();
            let (d_cell, w_cell, s_cell, t_cell) = if done.is_empty() {
                let lower = vals
                    .iter()
                    .map(|r| int(r, "solving_degree_max"))
                    .max()
                    .unwrap();
                (
                    format!("≥ {lower} (all timed out)"),
                    "—".to_string(),
                    "—".to_string(),
                    "—".to_string(),
                )
            } else {
                let ds: Vec<i64> = done.iter().map(|r| int(r, "solving_degree_max")).collect();
                let cells = (
                    format!(
                        "{} ({}–{})",
                        fmt_g(median(ds.iter().map(|&d| Num::Int(d)))),
                        ds.iter().min().unwrap(),
                        ds.iter().max().unwrap()
                    ),
                    fmt_g(median(done.iter().map(|r| num(r, "max_cols_to_solution")))),
                    fmt_g(median(done.iter().map(|r| num(r, "steps_to_solution")))),
                    format!("{:.1}", median(done.iter().map(|r| num(r, "ms"))).f()),
                );
                fit_by_n.insert(n, done);
                cells
            };
            println!(
                "| {n} | {} | {to} | {d_cell} | {w_cell} | {s_cell} | {t_cell} |",
                vals.len()
            );
        }
        if fit_by_n.len() < 2 {
            continue;
        }
        let xs: Vec<Num> = fit_by_n
            .iter()
            .flat_map(|(&n, vs)| vs.iter().map(move |_| Num::Int(n)))
            .collect();
        let ys: Vec<Num> = fit_by_n.values().flatten().map(|r| degree(r)).collect();
        let (a, b) = lsq(&xs, &ys).expect("two values of N");
        let ci = bootstrap_slope(&fit_by_n, degree).expect("two values of N");
        let wys: Vec<Num> = fit_by_n.values().flatten().map(|r| log2_width(r)).collect();
        let (_, wb) = lsq(&xs, &wys).expect("two values of N");
        let degenerate = ci.0 == ci.1;
        println!(
            "\nD = {a:.2} + {b:.3}·N over N = {}…{} ({} values). Pre-registered bootstrap (targets within N): beta 95% [{:.3}, {:.3}]{}; rule on it: {}.",
            fit_by_n.keys().next().unwrap(),
            fit_by_n.keys().next_back().unwrap(),
            fit_by_n.len(),
            ci.0,
            ci.1,
            if degenerate {
                " (degenerate: no variance within any N)"
            } else {
                ""
            },
            verdict(ci.0, ci.1)
        );
        if let Some((lo, hi)) = jackknife_slopes(&fit_by_n, degree) {
            println!("Leave-one-N-out slopes of D: {lo:.3} to {hi:.3}.");
        }
        let up = upper_half_slope(&fit_by_n, degree);
        if let Some((Some(slope), first, last)) = up {
            println!("Slope of D over the upper half, N = {first}…{last}: {slope:.3}.");
        }
        let (lo_n, hi_n, d_last) = final_plateau(&fit_by_n, degree);
        println!(
            "Final plateau: D = {} over N = {lo_n}…{hi_n}, length L = {}.",
            fmt_g(d_last),
            hi_n - lo_n
        );
        println!("A1 (section 11.6): {}.", a1_reading(&fit_by_n, degree));
        let upper_width = if up.is_some() {
            let (slope, _, _) = upper_half_slope(&fit_by_n, log2_width).unwrap();
            format!("; {:.3} over the upper half.", slope.expect("a slope"))
        } else {
            ".".to_string()
        };
        println!("log2(width) slope {wb:.3} per unit N over the whole range{upper_width}");
    }
    let split_verdicts = compare_engines(&rows);
    let contradicted = confirmations(&capped_rows, &rows);
    let ladder: Vec<&Value> = rows
        .iter()
        .copied()
        .filter(|r| text(r, "control") == "ladder")
        .collect();
    if !ladder.is_empty() {
        println!("\n### Generator-count ladder (planted targets)\n");
        println!("| kind | N | g | runs | D (median, range) | width to solution (median) |");
        println!("|:--|--:|--:|--:|:--|--:|");
        let mut cells: BTreeMap<(String, i64, i64), Vec<&Value>> = BTreeMap::new();
        for r in ladder {
            cells
                .entry((text(r, "kind").to_string(), int(r, "N"), int(r, "g")))
                .or_default()
                .push(r);
        }
        for ((kind, n, g), rs) in &cells {
            let done: Vec<&Value> = rs.iter().copied().filter(|r| !timed_out(r)).collect();
            if done.is_empty() {
                println!("| {kind} | {n} | {g} | {} | all timed out | — |", rs.len());
            } else {
                let ds: Vec<i64> = done.iter().map(|r| int(r, "solving_degree_max")).collect();
                println!(
                    "| {kind} | {n} | {g} | {} | {} ({}–{}) | {} |",
                    rs.len(),
                    fmt_g(median(ds.iter().map(|&d| Num::Int(d)))),
                    ds.iter().min().unwrap(),
                    ds.iter().max().unwrap(),
                    fmt_g(median(done.iter().map(|r| num(r, "max_cols_to_solution"))))
                );
            }
        }
    }
    let bad = rows
        .iter()
        .filter(|r| field(r, "planted_ok") == &Value::Bool(false))
        .count();
    println!(
        "\nPlanted-solution violations: {bad} of {} rows.",
        rows.len()
    );
    // A repeat that disagrees, a lost planted solution, two engines that
    // disagree on whether a system has a solution, or a confirmation run
    // that contradicts its full run: each is a bug, not data.
    if mismatches > 0 || bad > 0 || split_verdicts > 0 || contradicted > 0 {
        ExitCode::FAILURE
    } else {
        ExitCode::SUCCESS
    }
}

// ── compare-rows, compare-trace and size-cap (round 7) ────────────

/// Flags of the form `--name a,b,c` after the positional arguments.
fn list_flag(args: &[String], name: &str) -> Vec<String> {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .map(|v| v.split(',').map(str::to_string).collect())
        .unwrap_or_default()
}

/// The rows of `NEW` against the rows of `REF` for the same systems: every
/// field equal but the wall clock, the `--may-fall` fields (which may only
/// fall) and the `--ignore` fields.
fn compare_rows(reference: &str, new: &str, args: &[String]) -> ExitCode {
    let may_fall: HashSet<String> = list_flag(args, "--may-fall").into_iter().collect();
    let ignore: HashSet<String> = list_flag(args, "--ignore")
        .into_iter()
        .chain(["ms".to_string(), "wall_s".to_string()])
        .collect();
    let refs: Vec<Value> = read_rows(reference)
        .into_iter()
        .filter(|r| r.get("N").is_some())
        .collect();
    let mut diffs = 0;
    let mut compared = 0;
    for r in read_rows(new).iter().filter(|r| r.get("N").is_some()) {
        let key = (instance_key(r), engine(r));
        let matches: Vec<&Value> = refs
            .iter()
            .filter(|o| (instance_key(o), engine(o)) == key)
            .collect();
        let label = format!(
            "{} m={} {} N={} target {}",
            text(r, "kind"),
            int(r, "m"),
            text(r, "control"),
            int(r, "N"),
            int(r, "target_index")
        );
        let Some(o) = matches
            .iter()
            .find(|o| timed_out(o) == timed_out(r))
            .or(matches.first())
        else {
            diffs += 1;
            println!("{label}: no reference row");
            continue;
        };
        compared += 1;
        let mut fields: Vec<&String> = r
            .as_object()
            .unwrap()
            .keys()
            .chain(o.as_object().unwrap().keys())
            .filter(|f| !ignore.contains(*f))
            .collect();
        fields.sort();
        fields.dedup();
        let mut bad: Vec<String> = Vec::new();
        for f in fields {
            let (x, y) = (field(r, f), field(o, f));
            let ok = if may_fall.contains(f) {
                let as_num = |v: &Value| num_of(v).unwrap_or(Num::Int(0));
                num_cmp(as_num(x), as_num(y)) != Ordering::Greater
            } else {
                py_eq(x, y)
            };
            if !ok {
                bad.push(format!("    {f}: new {x} / reference {y}"));
            }
        }
        if bad.is_empty() {
            let fell: Vec<String> = may_fall
                .iter()
                .filter(|f| !py_eq(field(r, f), field(o, f)))
                .map(|f| format!("{f} {} -> {}", field(o, f), field(r, f)))
                .collect();
            println!(
                "{label}: as required{}",
                if fell.is_empty() {
                    String::new()
                } else {
                    format!(" ({})", fell.join(", "))
                }
            );
        } else {
            diffs += 1;
            println!("DIFF {label}:");
            for b in bad {
                println!("{b}");
            }
        }
    }
    println!("{compared} rows compared, {diffs} differences");
    if diffs == 0 && compared > 0 {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}

/// A `tower step` line without its times and its memory: what is left is
/// determined by the system and the algorithm.
fn step_counts(line: &str) -> String {
    let body = line.split_once(": ").map_or(line, |(_, b)| b);
    let body = body.split(", memory ").next().unwrap_or(body);
    // Drop "X ms (rows …, A …, B …, update …), " after the pairs left.
    match (body.find(" ms (rows "), body.find("), B' ")) {
        (Some(ms), Some(end)) if ms < end => {
            let start = body[..ms].rfind(", ").map_or(ms, |i| i + 2);
            format!("{}{}", &body[..start], &body[end + 3..])
        }
        _ => body.to_string(),
    }
}

/// The `tower step` lines of a log, one list per system.
fn trace_systems(path: &str) -> Vec<Vec<String>> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut systems: Vec<Vec<String>> = Vec::new();
    let mut cur: Vec<String> = Vec::new();
    for line in text.lines() {
        if line.starts_with("tower step ") {
            cur.push(step_counts(line));
        } else if line.starts_with("tower stop in step ") {
            continue;
        } else if !cur.is_empty() {
            systems.push(std::mem::take(&mut cur));
        }
    }
    if !cur.is_empty() {
        systems.push(cur);
    }
    systems
}

/// The steps of `NEW`'s systems against `REF`'s, system by system, up to
/// `--steps K` steps each or as many as both printed.
fn compare_trace(reference: &str, new: &str, args: &[String]) -> ExitCode {
    let limit: Option<usize> = list_flag(args, "--steps")
        .first()
        .map(|k| k.parse().expect("--steps K"));
    let (a, b) = (trace_systems(reference), trace_systems(new));
    let (mut steps, mut same, mut short) = (0, 0, 0);
    if a.len() != b.len() {
        println!(
            "{} systems traced in the reference, {} in the new log",
            a.len(),
            b.len()
        );
    }
    for (s, (x, y)) in a.iter().zip(&b).enumerate() {
        let n = limit.unwrap_or(x.len().min(y.len()));
        if x.len() < n || y.len() < n {
            short += 1;
            println!(
                "system {}: {} reference steps and {} new, {n} required",
                s + 1,
                x.len(),
                y.len()
            );
        }
        for (k, (p, q)) in x.iter().zip(y).take(n).enumerate() {
            steps += 1;
            if p == q {
                same += 1;
            } else {
                println!(
                    "system {} step {}:\n    reference {p}\n    new       {q}",
                    s + 1,
                    k + 1
                );
            }
        }
    }
    println!("{steps} steps compared, {same} identical");
    if same == steps && short == 0 && steps > 0 && a.len() == b.len() {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}

/// The number after `prefix` in `line`, up to the next space or comma.
fn number_after(line: &str, prefix: &str) -> Option<u64> {
    let i = line.find(prefix)? + prefix.len();
    let digits: String = line[i..].chars().take_while(char::is_ascii_digit).collect();
    digits.parse().ok()
}

/// A log's `tower step` lines, and its `tower stop in step` line, as a
/// Markdown table: what each step reduced, what it found, what it cost.
fn trace_table(log: &str) -> ExitCode {
    let text = std::fs::read_to_string(log).unwrap_or_else(|e| panic!("{log}: {e}"));
    println!("| step | degree | S-rows | pivot rows | columns | without a divisor | residues | new elements (lowest degree) | `B'` entries | kept entries | multiply-adds | s | memory MB |");
    println!("|--:|--:|--:|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|");
    let mut rows = 0;
    for line in text.lines() {
        let stop = line.starts_with("tower stop in step ");
        if !stop && !line.starts_with("tower step ") {
            continue;
        }
        let n = |prefix: &str| number_after(line, prefix);
        let show = |v: Option<u64>| v.map_or("—".to_string(), |x| x.to_string());
        let step = if stop {
            n("tower stop in step ")
        } else {
            n("tower step ")
        };
        let residues = line
            .split_once("residue ")
            .map(|(_, rest)| rest.split_whitespace().collect::<Vec<_>>());
        let (res, q) = match residues.as_deref() {
            Some([r, "x", q, ..]) => (r.to_string(), q.trim_end_matches(',').to_string()),
            _ => ("—".to_string(), "—".to_string()),
        };
        let seconds = line
            .split_once("pairs left ")
            .and_then(|(_, rest)| rest.split_once(", "))
            .and_then(|(_, rest)| rest.split_once(" ms"))
            .and_then(|(ms, _)| ms.parse::<f64>().ok())
            .map_or("—".to_string(), |ms| format!("{:.1}", ms / 1e3));
        let pivots = n("S-rows + ")
            .zip(n("reducers + "))
            .map(|(r, p)| r + p);
        println!(
            "| {}{} | {} | {} | {} | {} | {q} | {} | {} | {} | {} | {} | {seconds} | {} |",
            show(step),
            if stop { " (stopped)" } else { "" },
            show(n(": degree ")),
            show(n("pairs, ")),
            show(pivots),
            show(n("promoted x ")),
            if stop { "—".to_string() } else { res },
            if stop {
                "—".to_string()
            } else {
                format!("{} ({})", show(n("fresh ")), show(n("lowest degree ")))
            },
            show(n("B' ")),
            show(n("elements with ")),
            if stop {
                "—".to_string()
            } else {
                show(n("muladds "))
            },
            show(n("memory ")),
        );
        rows += 1;
    }
    if rows > 0 {
        ExitCode::SUCCESS
    } else {
        println!("{log}: no tower step lines");
        ExitCode::FAILURE
    }
}

/// Note §16.3's rule for stage 2's cap, from the `tower stop in step` line
/// of a stage-1 log: `C` is the largest multiple of `10⁸` with
/// `H + 4C ≤ A − 2.5·10⁹` bytes, where `A` is the address-space cap
/// (`ulimit -v 14000000`, in KiB) and `H` the live heap at the stop (the
/// line's MB are 2^20 bytes).
fn size_cap(log: &str) -> ExitCode {
    let text = std::fs::read_to_string(log).unwrap_or_else(|e| panic!("{log}: {e}"));
    let Some(line) = text.lines().find(|l| l.starts_with("tower stop in step ")) else {
        println!("{log}: no `tower stop in step` line");
        return ExitCode::FAILURE;
    };
    let step = number_after(line, "tower stop in step ").expect("a step");
    let e = number_after(line, "B' ").expect("B' entries");
    let q = number_after(line, "residue 0 x ").expect("columns without a divisor");
    let h_mb = number_after(line, "live heap ").expect("a live heap");
    let a: u64 = 14_000_000 * 1024;
    let reserve: u64 = 2_500_000_000;
    let h = h_mb << 20;
    let c = a.saturating_sub(reserve).saturating_sub(h) / 4 / 100_000_000 * 100_000_000;
    println!("stopped step {step}: E = {e} entries of B', q = {q} columns without a divisor, live heap H = {h_mb} MB = {h} bytes");
    println!(
        "A = {a} bytes, reserve 2.5e9 bytes: C = {c} ({})",
        if c >= e {
            "C >= E: the step fits by the rule"
        } else {
            "C < E: the rule does not expect the step to fit; stage 2 runs once with --max-dense E"
        }
    );
    ExitCode::SUCCESS
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let usage = "usage: pkm_tower_check verify FILE… | analyze FILE… | compare-rows REF NEW [--may-fall F,…] [--ignore F,…] | compare-trace REF NEW [--steps K] | size-cap LOG | trace-table LOG";
    let files = |from: usize| -> Vec<String> {
        args[from..]
            .iter()
            .take_while(|a| !a.starts_with("--"))
            .cloned()
            .collect()
    };
    match args.first().map(String::as_str) {
        Some("verify") => verify(&files(1)),
        Some("analyze") => analyze(&files(1)),
        Some("compare-rows") if args.len() >= 3 => compare_rows(&args[1], &args[2], &args[3..]),
        Some("compare-trace") if args.len() >= 3 => compare_trace(&args[1], &args[2], &args[3..]),
        Some("size-cap") if args.len() == 2 => size_cap(&args[1]),
        Some("trace-table") if args.len() == 2 => trace_table(&args[1]),
        _ => {
            eprintln!("{usage}");
            ExitCode::from(2)
        }
    }
}

#[cfg(test)]
mod tests {
    //! The fixtures of `test_verify.py`, the legacy tests of `verify.py`:
    //! planted decompositions on random curves over `F_1009`, and the two
    //! orders of eliminating the `m = 4` chain's free unknowns. Then the
    //! Python semantics `analyze` relies on.
    use super::*;
    use rand::rngs::StdRng;
    use rand::seq::SliceRandom;
    use rand::{Rng, SeedableRng};

    const P: u64 = 1009;

    type Point = Option<(u64, u64)>;

    fn add(f: Fp, p: Point, q: Point, a: u64) -> Point {
        let (Some((x1, y1)), Some((x2, y2))) = (p, q) else {
            return p.or(q);
        };
        if x1 == x2 && f.add(y1, y2) == 0 {
            return None;
        }
        let slope = if p == q {
            f.mul(f.add(f.mul(3, f.mul(x1, x1)), a), f.inv(f.mul(2, y1)))
        } else {
            f.mul(f.sub(y2, y1), f.inv(f.sub(x2, x1)))
        };
        let x3 = f.sub(f.sub(f.mul(slope, slope), x1), x2);
        Some((x3, f.sub(f.mul(slope, f.sub(x1, x3)), y1)))
    }

    /// Random curves over `F_P` with their points.
    fn curves(rng: &mut StdRng, count: usize) -> Vec<(u64, u64, Vec<(u64, u64)>)> {
        let f = Fp(P);
        let mut out = Vec::new();
        while out.len() < count {
            let (a, b) = (rng.gen_range(0..P), rng.gen_range(1..P));
            let disc = f.add(f.mul(4, f.pow(a, 3)), f.mul(27, f.mul(b, b)));
            if disc == 0 {
                continue;
            }
            let pts: Vec<(u64, u64)> = (0..P)
                .flat_map(|x| (0..P).map(move |y| (x, y)))
                .filter(|&(x, y)| f.mul(y, y) == f.add(f.add(f.pow(x, 3), f.mul(a, x)), b))
                .collect();
            out.push((a, b, pts));
        }
        out
    }

    /// `m` points with distinct abscissae and the abscissa of their sum.
    fn planted(rng: &mut StdRng, a: u64, pts: &[(u64, u64)], m: usize) -> Option<(Vec<u64>, u64)> {
        let f = Fp(P);
        let chosen: Vec<(u64, u64)> = loop {
            let c: Vec<(u64, u64)> = pts.choose_multiple(rng, m).copied().collect();
            let xs: HashSet<u64> = c.iter().map(|q| q.0).collect();
            if xs.len() == m {
                break c;
            }
        };
        let total = chosen
            .iter()
            .fold(None, |acc, &q| add(f, acc, Some(q), a))?;
        Some((chosen.iter().map(|q| q.0).collect(), total.0))
    }

    #[test]
    fn planted_decompositions_are_zeros_of_the_summation_polynomials() {
        let f = Fp(P);
        let mut rng = StdRng::seed_from_u64(20_260_924);
        let cs = curves(&mut rng, 12);
        let mut checked = [0; 3];
        for (a, b, pts) in &cs {
            if let Some((x, xr)) = planted(&mut rng, *a, pts, 2) {
                assert_eq!(s3(f, x[0], x[1], xr, *a, *b), 0);
                checked[0] += 1;
            }
            if let Some((x, xr)) = planted(&mut rng, *a, pts, 3) {
                let g = s3_coeffs(f, x[0], x[1], *a, *b);
                let h = s3_coeffs(f, x[2], xr, *a, *b);
                assert_eq!(res2(f, g, h), 0);
                checked[1] += 1;
            }
            if let Some((x, xr)) = planted(&mut rng, *a, pts, 4) {
                let (aa, bb, cc) = s3_coeffs(f, x[3], xr, *a, *b);
                let s4 = s4_in_u(f, x[0], x[1], x[2], *a, *b);
                assert!(res_zero(f, &s4, &[cc, bb, aa]));
                checked[2] += 1;
            }
        }
        assert!(checked.iter().all(|&c| c >= 8), "{checked:?}");
    }

    /// `Res_w(S3(U, x3, w), S3(w, x4, xR))` as a polynomial in `U`: the
    /// chain's other end eliminated first.
    fn s4_other_order(f: Fp, x3: u64, x4: u64, xr: u64, a: u64, b: u64) -> Vec<u64> {
        let (aa, bb, cc) = s3_coeffs(f, x4, xr, a, b);
        let [a2, b2, c2] = s3_in_w(f, x3, a, b);
        let u = psub(f, &scale(f, aa, &c2), &scale(f, cc, &a2));
        let v = psub(f, &scale(f, aa, &b2), &scale(f, bb, &a2));
        let w = psub(f, &scale(f, bb, &c2), &scale(f, cc, &b2));
        psub(f, &pmul(f, &u, &u), &pmul(f, &v, &w))
    }

    #[test]
    fn both_elimination_orders_of_m4_agree() {
        let f = Fp(P);
        let mut rng = StdRng::seed_from_u64(20_260_925);
        let mut zeros = 0;
        for (a, b, pts) in curves(&mut rng, 2) {
            let mut v: Vec<u64> = pts.choose_multiple(&mut rng, 7).map(|q| q.0).collect();
            v.sort_unstable();
            v.dedup();
            let xr = pts.choose(&mut rng).unwrap().0;
            for &x1 in &v {
                for &x2 in &v {
                    let g = s3_coeffs(f, x1, x2, a, b);
                    for &x3 in &v {
                        let s4 = s4_in_u(f, x1, x2, x3, a, b);
                        for &x4 in &v {
                            let (aa, bb, cc) = s3_coeffs(f, x4, xr, a, b);
                            let one = res_zero(f, &s4, &[cc, bb, aa]);
                            let other =
                                res_zero(f, &[g.2, g.1, g.0], &s4_other_order(f, x3, x4, xr, a, b));
                            assert_eq!(one, other, "{x1} {x2} {x3} {x4}");
                            zeros += usize::from(one);
                        }
                    }
                }
            }
        }
        // The grid is small against P, but not so small that no tuple
        // solves it.
        assert!(zeros > 0);
    }

    /// The `m = 2` counter against the group law: for distinct abscissae of
    /// points over `F_p`, `S3(x1, x2, xR) = 0` exactly when some choice of
    /// the points' signs sums to a point with abscissa `xR`.
    #[test]
    fn the_m2_counter_matches_point_sums() {
        let f = Fp(P);
        let mut rng = StdRng::seed_from_u64(7);
        for (a, b, pts) in curves(&mut rng, 6) {
            let xs: Vec<u64> = {
                let mut v: Vec<u64> = pts.iter().map(|q| q.0).collect();
                v.sort_unstable();
                v.dedup();
                v.into_iter().take(24).collect()
            };
            let xr = pts.choose(&mut rng).unwrap().0;
            let by_x = |x: u64| pts.iter().filter(move |q| q.0 == x).copied();
            for &x1 in &xs {
                for &x2 in &xs {
                    if x1 == x2 {
                        continue;
                    }
                    let sums = by_x(x1).flat_map(|p| by_x(x2).map(move |q| (p, q)));
                    let hit = sums
                        .filter_map(|(p, q)| add(f, Some(p), Some(q), a))
                        .any(|s| s.0 == xr);
                    assert_eq!(s3(f, x1, x2, xr, a, b) == 0, hit, "{x1} {x2} {xr}");
                }
            }
        }
    }

    #[test]
    fn the_twister_draws_what_python_draws() {
        // random.seed(0); random.random() and random.seed(42); random.random().
        assert_eq!(PyRandom::new(0).random(), 0.844_421_851_525_048_1);
        assert_eq!(PyRandom::new(42).random(), 0.639_426_798_457_883_7);
        // _randbelow stays below its bound and reaches every value.
        let mut r = PyRandom::new(20_260_924);
        let mut seen = [false; 5];
        for _ in 0..200 {
            seen[r.randbelow(5)] = true;
        }
        assert!(seen.iter().all(|&s| s));
    }

    #[test]
    fn numbers_follow_python() {
        // sum() compensates since 3.12: ten 0.1s make 1.0 exactly.
        assert_eq!(py_sum(vec![Num::Float(0.1); 10]), Num::Float(1.0));
        assert_eq!(py_sum(vec![Num::Int(2), Num::Int(3)]), Num::Int(5));
        assert_eq!(median(vec![Num::Int(6), Num::Int(7)]), Num::Float(6.5));
        assert_eq!(
            median(vec![Num::Int(7), Num::Int(6), Num::Int(9)]),
            Num::Int(7)
        );
        assert_eq!(py_square(1.5), 2.25);
        for (x, s) in [
            (6.0, "6"),
            (6.5, "6.5"),
            (238.5, "238.5"),
            (31397.0, "31397"),
            (104600.0, "104600"),
            (1234567.0, "1.23457e+06"),
            (0.0001, "0.0001"),
            (0.00001234, "1.234e-05"),
            (999999.5, "1e+06"),
            (2.5, "2.5"),
        ] {
            assert_eq!(fmt_g(Num::Float(x)), s);
        }
        assert_eq!(py_repr("copy's"), "\"copy's\"");
        assert_eq!(py_repr("kummer"), "'kummer'");
        let line = "tower step 3: degree 5, 4 critical + 5 tower pairs, 9 S-rows + 13 reducers + 4 promoted x 99 cols, nnz 326, residue 5 x 82, fresh 5 (lowest degree 4), basis 10, pairs left 31, 0.7 ms (rows 0, A 0, B 0, update 0), B' 1203 entries, kept 10 elements with 452 entries, muladds 2915, memory 26 MB (peak 64 MB)";
        assert_eq!(
            step_counts(line),
            "degree 5, 4 critical + 5 tower pairs, 9 S-rows + 13 reducers + 4 promoted x 99 cols, nnz 326, residue 5 x 82, fresh 5 (lowest degree 4), basis 10, pairs left 31, B' 1203 entries, kept 10 elements with 452 entries, muladds 2915"
        );
    }
}

//! JSON and numbers the way CPython reads and prints them, for the
//! ECC2K-130 tools ported from `ecc2k130/aws/*.py`.
//!
//! A port has to write what the Python wrote, byte for byte, so this is
//! `json.loads` and `json.dumps` with Python's dict semantics (insertion
//! order, in-place update, `True == 1`), `repr` and `%g` of a float,
//! `int()` and `float()` of a JSON value, `repr` and `str` of one, and the
//! builtin `sum` of floats.

// ── JSON as Python reads and writes it ─────────────────────────────

/// A JSON value with Python's semantics: objects keep insertion order and
/// update a key in place, and `True == 1`.
#[derive(Clone, Debug)]
pub enum Json {
    Null,
    Bool(bool),
    Int(i128),
    /// Written only by [`dumps_floats`]: no merge file carries one.
    Float(f64),
    Str(String),
    List(Vec<Json>),
    Obj(Obj),
}

#[derive(Clone, Debug, Default)]
pub struct Obj(Vec<(String, Json)>);

impl Obj {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    pub fn get(&self, key: &str) -> Option<&Json> {
        self.0.iter().find(|(k, _)| k == key).map(|(_, v)| v)
    }

    pub fn get_mut(&mut self, key: &str) -> Option<&mut Json> {
        self.0.iter_mut().find(|(k, _)| k == key).map(|(_, v)| v)
    }

    /// `d[key] = value`: an existing key keeps its position.
    pub fn set(&mut self, key: &str, value: Json) {
        match self.get_mut(key) {
            Some(slot) => *slot = value,
            None => self.0.push((key.to_string(), value)),
        }
    }

    pub fn contains(&self, key: &str) -> bool {
        self.get(key).is_some()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = &(String, Json)> {
        self.0.iter()
    }

    pub fn with(mut self, key: &str, value: Json) -> Self {
        self.set(key, value);
        self
    }
}

impl PartialEq for Obj {
    fn eq(&self, other: &Self) -> bool {
        self.len() == other.len() && self.iter().all(|(k, v)| other.get(k) == Some(v))
    }
}

impl PartialEq for Json {
    fn eq(&self, other: &Self) -> bool {
        use Json::*;
        match (self, other) {
            (Null, Null) => true,
            (Bool(a), Bool(b)) => a == b,
            (Bool(a), Int(b)) | (Int(b), Bool(a)) => i128::from(*a) == *b,
            (Bool(a), Float(b)) | (Float(b), Bool(a)) => f64::from(u8::from(*a)) == *b,
            (Int(a), Int(b)) => a == b,
            (Int(a), Float(b)) | (Float(b), Int(a)) => (*a as f64) == *b,
            (Float(a), Float(b)) => a == b,
            (Str(a), Str(b)) => a == b,
            (List(a), List(b)) => a == b,
            (Obj(a), Obj(b)) => a == b,
            _ => false,
        }
    }
}

impl Json {
    /// Python truthiness.
    pub fn truthy(&self) -> bool {
        match self {
            Json::Null => false,
            Json::Bool(b) => *b,
            Json::Int(i) => *i != 0,
            Json::Float(f) => *f != 0.0,
            Json::Str(s) => !s.is_empty(),
            Json::List(l) => !l.is_empty(),
            Json::Obj(o) => !o.is_empty(),
        }
    }

    /// `type(x) is int`: a bool is not one.
    pub fn as_int(&self) -> Option<i128> {
        match self {
            Json::Int(i) => Some(*i),
            _ => None,
        }
    }

    pub fn as_str(&self) -> Option<&str> {
        match self {
            Json::Str(s) => Some(s),
            _ => None,
        }
    }

    pub fn as_obj(&self) -> Option<&Obj> {
        match self {
            Json::Obj(o) => Some(o),
            _ => None,
        }
    }

    pub fn str(s: &str) -> Json {
        Json::Str(s.to_string())
    }
}

/// `json.loads`: keys in document order, a repeated key updating the first.
pub fn parse_json(text: &str) -> Result<Json, String> {
    let mut p = Parser {
        s: text.as_bytes(),
        i: 0,
    };
    p.ws();
    let v = p.value()?;
    p.ws();
    if p.i != p.s.len() {
        return Err(format!("extra data at byte {}", p.i));
    }
    Ok(v)
}

struct Parser<'a> {
    s: &'a [u8],
    i: usize,
}

impl Parser<'_> {
    fn peek(&self) -> Option<u8> {
        self.s.get(self.i).copied()
    }

    fn ws(&mut self) {
        while matches!(self.peek(), Some(b' ' | b'\t' | b'\n' | b'\r')) {
            self.i += 1;
        }
    }

    fn expect(&mut self, b: u8) -> Result<(), String> {
        if self.peek() == Some(b) {
            self.i += 1;
            Ok(())
        } else {
            Err(format!("expected {:?} at byte {}", b as char, self.i))
        }
    }

    fn literal(&mut self, word: &str, v: Json) -> Result<Json, String> {
        if self.s[self.i..].starts_with(word.as_bytes()) {
            self.i += word.len();
            Ok(v)
        } else {
            Err(format!("invalid literal at byte {}", self.i))
        }
    }

    fn value(&mut self) -> Result<Json, String> {
        match self.peek() {
            Some(b'{') => self.object(),
            Some(b'[') => self.array(),
            Some(b'"') => Ok(Json::Str(self.string()?)),
            Some(b't') => self.literal("true", Json::Bool(true)),
            Some(b'f') => self.literal("false", Json::Bool(false)),
            Some(b'n') => self.literal("null", Json::Null),
            Some(b'N') => self.literal("NaN", Json::Float(f64::NAN)),
            Some(b'I') => self.literal("Infinity", Json::Float(f64::INFINITY)),
            Some(b'-') if self.s[self.i..].starts_with(b"-Infinity") => {
                self.literal("-Infinity", Json::Float(f64::NEG_INFINITY))
            }
            Some(b'-' | b'0'..=b'9') => self.number(),
            _ => Err(format!("expecting value at byte {}", self.i)),
        }
    }

    fn object(&mut self) -> Result<Json, String> {
        self.expect(b'{')?;
        let mut obj = Obj::new();
        self.ws();
        if self.peek() == Some(b'}') {
            self.i += 1;
            return Ok(Json::Obj(obj));
        }
        loop {
            self.ws();
            let key = self.string()?;
            self.ws();
            self.expect(b':')?;
            self.ws();
            let v = self.value()?;
            obj.set(&key, v);
            self.ws();
            match self.peek() {
                Some(b',') => self.i += 1,
                Some(b'}') => {
                    self.i += 1;
                    return Ok(Json::Obj(obj));
                }
                _ => return Err(format!("expecting ',' or '}}' at byte {}", self.i)),
            }
        }
    }

    fn array(&mut self) -> Result<Json, String> {
        self.expect(b'[')?;
        let mut out = Vec::new();
        self.ws();
        if self.peek() == Some(b']') {
            self.i += 1;
            return Ok(Json::List(out));
        }
        loop {
            self.ws();
            out.push(self.value()?);
            self.ws();
            match self.peek() {
                Some(b',') => self.i += 1,
                Some(b']') => {
                    self.i += 1;
                    return Ok(Json::List(out));
                }
                _ => return Err(format!("expecting ',' or ']' at byte {}", self.i)),
            }
        }
    }

    fn hex4(&mut self) -> Result<u32, String> {
        let digits = self
            .s
            .get(self.i..self.i + 4)
            .ok_or("truncated \\u escape")?;
        let text = std::str::from_utf8(digits).map_err(|_| "invalid \\u escape")?;
        let v = u32::from_str_radix(text, 16).map_err(|_| "invalid \\u escape")?;
        self.i += 4;
        Ok(v)
    }

    fn string(&mut self) -> Result<String, String> {
        self.expect(b'"')?;
        let mut out: Vec<u8> = Vec::new();
        loop {
            let c = self.peek().ok_or("unterminated string")?;
            self.i += 1;
            match c {
                b'"' => break,
                b'\\' => {
                    let e = self.peek().ok_or("unterminated escape")?;
                    self.i += 1;
                    let ch = match e {
                        b'"' => '"',
                        b'\\' => '\\',
                        b'/' => '/',
                        b'b' => '\u{8}',
                        b'f' => '\u{c}',
                        b'n' => '\n',
                        b'r' => '\r',
                        b't' => '\t',
                        b'u' => {
                            let hi = self.hex4()?;
                            let code = if (0xD800..0xDC00).contains(&hi)
                                && self.s[self.i..].starts_with(b"\\u")
                            {
                                self.i += 2;
                                let lo = self.hex4()?;
                                if !(0xDC00..0xE000).contains(&lo) {
                                    return Err("unpaired surrogate".into());
                                }
                                0x10000 + ((hi - 0xD800) << 10) + (lo - 0xDC00)
                            } else {
                                hi
                            };
                            char::from_u32(code).ok_or("unpaired surrogate")?
                        }
                        _ => return Err(format!("invalid escape at byte {}", self.i)),
                    };
                    let mut buf = [0u8; 4];
                    out.extend_from_slice(ch.encode_utf8(&mut buf).as_bytes());
                }
                0x00..=0x1f => return Err(format!("control character at byte {}", self.i)),
                _ => out.push(c),
            }
        }
        String::from_utf8(out).map_err(|_| "invalid UTF-8 in string".to_string())
    }

    fn number(&mut self) -> Result<Json, String> {
        let start = self.i;
        if self.peek() == Some(b'-') {
            self.i += 1;
        }
        match self.peek() {
            Some(b'0') => self.i += 1,
            Some(b'1'..=b'9') => {
                while matches!(self.peek(), Some(b'0'..=b'9')) {
                    self.i += 1;
                }
            }
            _ => return Err(format!("invalid number at byte {start}")),
        }
        let mut float = false;
        if self.peek() == Some(b'.') && matches!(self.s.get(self.i + 1), Some(b'0'..=b'9')) {
            float = true;
            self.i += 1;
            while matches!(self.peek(), Some(b'0'..=b'9')) {
                self.i += 1;
            }
        }
        if matches!(self.peek(), Some(b'e' | b'E')) {
            let mut j = self.i + 1;
            if matches!(self.s.get(j), Some(b'+' | b'-')) {
                j += 1;
            }
            if matches!(self.s.get(j), Some(b'0'..=b'9')) {
                float = true;
                self.i = j;
                while matches!(self.peek(), Some(b'0'..=b'9')) {
                    self.i += 1;
                }
            }
        }
        let text = std::str::from_utf8(&self.s[start..self.i]).expect("ASCII digits");
        if float {
            text.parse::<f64>()
                .map(Json::Float)
                .map_err(|e| e.to_string())
        } else {
            text.parse::<i128>()
                .map(Json::Int)
                .map_err(|e| e.to_string())
        }
    }
}

/// How `json.dumps` was called.
#[derive(Clone, Copy)]
pub enum Style {
    /// `separators=(",", ":")`.
    Compact,
    /// `indent=n`.
    Indent(usize),
}

/// `json.dumps(v, ...)` with `ensure_ascii`, for files that carry no float:
/// one is refused rather than written.
pub fn dumps(v: &Json, style: Style, sort_keys: bool) -> Result<String, String> {
    let mut out = String::new();
    write_json(&mut out, v, style, sort_keys, false, 0)?;
    Ok(out)
}

/// `json.dumps(v, ...)` with `ensure_ascii` and the default `allow_nan`:
/// floats as `repr` prints them, and the non-finite ones as `NaN`,
/// `Infinity` and `-Infinity`.
pub fn dumps_floats(v: &Json, style: Style, sort_keys: bool) -> String {
    let mut out = String::new();
    write_json(&mut out, v, style, sort_keys, true, 0).expect("floats are written");
    out
}

fn write_json(
    out: &mut String,
    v: &Json,
    style: Style,
    sort_keys: bool,
    floats: bool,
    level: usize,
) -> Result<(), String> {
    match v {
        Json::Null => out.push_str("null"),
        Json::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
        Json::Int(i) => out.push_str(&i.to_string()),
        Json::Float(_) if !floats => {
            return Err("refusing to write a float: no merge file carries one".into())
        }
        Json::Float(f) if f.is_nan() => out.push_str("NaN"),
        Json::Float(f) if f.is_infinite() => {
            out.push_str(if *f > 0.0 { "Infinity" } else { "-Infinity" })
        }
        Json::Float(f) => out.push_str(&float_repr(*f)),
        Json::Str(s) => escape_into(out, s),
        Json::List(items) => {
            if items.is_empty() {
                out.push_str("[]");
                return Ok(());
            }
            out.push('[');
            for (n, item) in items.iter().enumerate() {
                if n > 0 {
                    out.push(',');
                }
                newline(out, style, level + 1);
                write_json(out, item, style, sort_keys, floats, level + 1)?;
            }
            newline(out, style, level);
            out.push(']');
        }
        Json::Obj(obj) => {
            if obj.is_empty() {
                out.push_str("{}");
                return Ok(());
            }
            let mut entries: Vec<&(String, Json)> = obj.iter().collect();
            if sort_keys {
                entries.sort_by(|a, b| a.0.cmp(&b.0));
            }
            out.push('{');
            for (n, (k, item)) in entries.into_iter().enumerate() {
                if n > 0 {
                    out.push(',');
                }
                newline(out, style, level + 1);
                escape_into(out, k);
                out.push_str(match style {
                    Style::Compact => ":",
                    Style::Indent(_) => ": ",
                });
                write_json(out, item, style, sort_keys, floats, level + 1)?;
            }
            newline(out, style, level);
            out.push('}');
        }
    }
    Ok(())
}

fn newline(out: &mut String, style: Style, level: usize) {
    if let Style::Indent(n) = style {
        out.push('\n');
        out.extend(std::iter::repeat_n(' ', n * level));
    }
}

/// `py_encode_basestring_ascii`: printable ASCII stays, the rest is escaped.
fn escape_into(out: &mut String, s: &str) {
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{8}' => out.push_str("\\b"),
            '\u{c}' => out.push_str("\\f"),
            ' '..='~' => out.push(c),
            _ => {
                let mut units = [0u16; 2];
                for unit in c.encode_utf16(&mut units) {
                    out.push_str(&format!("\\u{unit:04x}"));
                }
            }
        }
    }
    out.push('"');
}

// ── Python numbers and text ────────────────────────────────────────

/// `repr(x)`: the shortest digits that read back as `x`, laid out fixed
/// while the decimal point falls within 16 digits of them and in
/// exponent form otherwise.
pub fn float_repr(x: f64) -> String {
    if !x.is_finite() {
        return non_finite(x);
    }
    let (digits, point) = shortest_digits(x.abs());
    let n = digits.len() as i32;
    let body = if -4 < point && point <= 16 {
        if point <= 0 {
            format!("0.{}{digits}", "0".repeat(-point as usize))
        } else if point >= n {
            format!("{digits}{}.0", "0".repeat((point - n) as usize))
        } else {
            let (whole, frac) = digits.split_at(point as usize);
            format!("{whole}.{frac}")
        }
    } else {
        let (first, rest) = digits.split_at(1);
        let dot = if rest.is_empty() { "" } else { "." };
        format!("{first}{dot}{rest}{}", exponent(point - 1))
    };
    let sign = if x.is_sign_negative() { "-" } else { "" };
    format!("{sign}{body}")
}

/// The shortest digits that read back as `x >= 0`, and where the decimal
/// point falls in them. When `x` sits exactly halfway between two such
/// strings, Rust's shortest formatting takes the upper and CPython's
/// `dtoa` the one ending in an even digit, so that one is taken here.
fn shortest_digits(x: f64) -> (String, i32) {
    let split = |s: String| -> (String, i32) {
        let (mantissa, exp) = s.split_once('e').expect("LowerExp has an exponent");
        let digits = mantissa.chars().filter(|c| *c != '.').collect();
        (digits, exp.parse::<i32>().expect("LowerExp exponent") + 1)
    };
    // LowerExp without a precision is the shortest round-trip digits; with
    // 800 it is every digit of the exact binary value.
    let (digits, point) = split(format!("{x:e}"));
    let (exact, exact_point) = split(format!("{x:.800e}"));
    let exact = exact.trim_end_matches('0');
    let n = digits.len();
    if exact.len() != n + 1 || !exact.ends_with('5') || exact_point != point {
        return (digits, point);
    }
    let lower = exact[..n].to_string();
    let Some(upper) = increment(&lower) else {
        return (digits, point);
    };
    let even = if lower.ends_with(['0', '2', '4', '6', '8']) {
        lower
    } else {
        upper
    };
    if even != digits && format!("0.{even}e{point}").parse::<f64>() == Ok(x) {
        (even, point)
    } else {
        (digits, point)
    }
}

/// A digit string plus one in its last place, if that keeps its length.
fn increment(digits: &str) -> Option<String> {
    let mut d = digits.as_bytes().to_vec();
    for i in (0..d.len()).rev() {
        if d[i] == b'9' {
            d[i] = b'0';
        } else {
            d[i] += 1;
            return String::from_utf8(d).ok();
        }
    }
    None
}

/// `'%.{precision}g' % x`.
pub fn format_g(x: f64, precision: usize) -> String {
    if !x.is_finite() {
        return non_finite(x);
    }
    if x == 0.0 {
        return if x.is_sign_negative() { "-0" } else { "0" }.into();
    }
    let precision = precision.max(1);
    let sci = format!("{:.*e}", precision - 1, x);
    let (mantissa, exp) = sci.split_once('e').expect("LowerExp has an exponent");
    let exp = exp.parse::<i32>().expect("LowerExp exponent");
    if exp < -4 || exp >= precision as i32 {
        format!("{}{}", strip_fraction_zeros(mantissa), exponent(exp))
    } else {
        let decimals = (precision as i32 - 1 - exp) as usize;
        strip_fraction_zeros(&format!("{x:.decimals$}")).to_string()
    }
}

fn non_finite(x: f64) -> String {
    if x.is_nan() {
        "nan".into()
    } else if x > 0.0 {
        "inf".into()
    } else {
        "-inf".into()
    }
}

/// `e+16`, `e-05`: a sign and at least two digits.
fn exponent(e: i32) -> String {
    format!("e{}{:02}", if e < 0 { '-' } else { '+' }, e.abs())
}

fn strip_fraction_zeros(s: &str) -> &str {
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.')
    } else {
        s
    }
}

/// `repr(v)`.
pub fn py_repr(v: &Json) -> String {
    match v {
        Json::Null => "None".into(),
        Json::Bool(b) => if *b { "True" } else { "False" }.into(),
        Json::Int(i) => i.to_string(),
        Json::Float(f) => float_repr(*f),
        Json::Str(s) => str_repr(s),
        Json::List(items) => {
            let parts: Vec<String> = items.iter().map(py_repr).collect();
            format!("[{}]", parts.join(", "))
        }
        Json::Obj(obj) => {
            let parts: Vec<String> = obj
                .iter()
                .map(|(k, v)| format!("{}: {}", str_repr(k), py_repr(v)))
                .collect();
            format!("{{{}}}", parts.join(", "))
        }
    }
}

/// `str(v)`: a string is itself, anything else its `repr`.
pub fn py_str(v: &Json) -> String {
    match v {
        Json::Str(s) => s.clone(),
        _ => py_repr(v),
    }
}

/// `repr(s)` of a str: single quotes unless only double quotes avoid an
/// escape, and anything unprintable escaped.
fn str_repr(s: &str) -> String {
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
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c == quote => {
                out.push('\\');
                out.push(c);
            }
            c if printable(c) => out.push(c),
            c if (c as u32) < 0x100 => out.push_str(&format!("\\x{:02x}", c as u32)),
            c if (c as u32) < 0x10000 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push_str(&format!("\\U{:08x}", c as u32)),
        }
    }
    out.push(quote);
    out
}

/// `str.isprintable` for one character: false for the general categories
/// Cc, Cf, Co, Zs (but the space), Zl and Zp. Unassigned code points (Cn)
/// need the Unicode database and count as printable here.
fn printable(c: char) -> bool {
    if c == ' ' {
        return true;
    }
    if c.is_control() || c.is_whitespace() {
        return false;
    }
    let format = matches!(
        c as u32,
        0xAD | 0x600..=0x605
            | 0x61C
            | 0x6DD
            | 0x70F
            | 0x890..=0x891
            | 0x8E2
            | 0x180E
            | 0x200B..=0x200F
            | 0x202A..=0x202E
            | 0x2060..=0x2064
            | 0x2066..=0x206F
            | 0xFEFF
            | 0xFFF9..=0xFFFB
            | 0x110BD
            | 0x110CD
            | 0x13430..=0x1343F
            | 0x1BCA0..=0x1BCA3
            | 0x1D173..=0x1D17A
            | 0xE0001
            | 0xE0020..=0xE007F
    );
    let private = matches!(c as u32, 0xE000..=0xF8FF | 0xF0000..=0xFFFFD | 0x100000..=0x10FFFD);
    !(format || private)
}

fn type_name(v: &Json) -> &'static str {
    match v {
        Json::Null => "NoneType",
        Json::Bool(_) => "bool",
        Json::Int(_) => "int",
        Json::Float(_) => "float",
        Json::Str(_) => "str",
        Json::List(_) => "list",
        Json::Obj(_) => "dict",
    }
}

/// `int(v)`. Strings take ASCII digits only, where Python also takes other
/// scripts' decimal digits.
pub fn py_int(v: &Json) -> Result<i128, String> {
    match v {
        Json::Bool(b) => Ok(i128::from(*b)),
        Json::Int(i) => Ok(*i),
        Json::Float(f) if f.is_nan() => Err("cannot convert float NaN to integer".into()),
        Json::Float(f) if f.is_infinite() => Err("cannot convert float infinity to integer".into()),
        Json::Float(f) => {
            let t = f.trunc();
            if t.abs() < 2f64.powi(127) {
                Ok(t as i128)
            } else {
                Err(format!("{} does not fit in 128 bits", float_repr(*f)))
            }
        }
        Json::Str(s) => int_literal(s)
            .ok_or_else(|| format!("invalid literal for int() with base 10: {}", str_repr(s))),
        _ => Err(format!(
            "int() argument must be a string, a bytes-like object or a real number, not '{}'",
            type_name(v)
        )),
    }
}

/// `float(v)`.
pub fn py_float(v: &Json) -> Result<f64, String> {
    match v {
        Json::Bool(b) => Ok(f64::from(u8::from(*b))),
        Json::Int(i) => Ok(*i as f64),
        Json::Float(f) => Ok(*f),
        Json::Str(s) => float_literal(s)
            .ok_or_else(|| format!("could not convert string to float: {}", str_repr(s))),
        _ => Err(format!(
            "float() argument must be a string or a real number, not '{}'",
            type_name(v)
        )),
    }
}

/// `str.isspace`.
pub fn py_space(c: char) -> bool {
    c.is_whitespace() || ('\u{1c}'..='\u{1f}').contains(&c)
}

/// `str.strip()`.
pub fn py_strip(s: &str) -> &str {
    s.trim_matches(py_space)
}

/// `str.splitlines()`.
pub fn py_splitlines(s: &str) -> Vec<String> {
    let mut lines = Vec::new();
    let mut cur = String::new();
    let mut chars = s.chars().peekable();
    while let Some(c) = chars.next() {
        match c {
            '\r' => {
                if chars.peek() == Some(&'\n') {
                    chars.next();
                }
                lines.push(std::mem::take(&mut cur));
            }
            '\n' | '\u{b}' | '\u{c}' | '\u{1c}' | '\u{1d}' | '\u{1e}' | '\u{85}' | '\u{2028}'
            | '\u{2029}' => lines.push(std::mem::take(&mut cur)),
            _ => cur.push(c),
        }
    }
    if !cur.is_empty() {
        lines.push(cur);
    }
    lines
}

/// `int(s)`: surrounding whitespace, a sign, and digits that underscores
/// may separate one at a time. The whitespace is `trim`'s: unlike
/// `str.strip`, CPython's number parsers keep `\x1c`..`\x1f`.
pub fn int_literal(s: &str) -> Option<i128> {
    let t = s.trim();
    let (negative, body) = match t.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, t.strip_prefix('+').unwrap_or(t)),
    };
    let digits = without_digit_underscores(body)?;
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let magnitude = digits.parse::<i128>().ok()?;
    Some(if negative { -magnitude } else { magnitude })
}

/// `float(s)`.
pub fn float_literal(s: &str) -> Option<f64> {
    let t = s.trim();
    let plain = without_digit_underscores(t)?;
    let body = plain.trim_start_matches(['+', '-']);
    let named = ["inf", "infinity", "nan"]
        .iter()
        .any(|w| body.eq_ignore_ascii_case(w));
    let numeric = !body.is_empty()
        && body
            .bytes()
            .all(|b| b.is_ascii_digit() || matches!(b, b'.' | b'e' | b'E' | b'+' | b'-'));
    if plain.len() - body.len() > 1 || !(named || numeric) {
        return None;
    }
    plain.parse::<f64>().ok()
}

/// `s` without its underscores, if every one sits between two digits.
fn without_digit_underscores(s: &str) -> Option<String> {
    let b = s.as_bytes();
    for (i, &c) in b.iter().enumerate() {
        if c == b'_' {
            let before = i > 0 && b[i - 1].is_ascii_digit();
            let after = b.get(i + 1).is_some_and(u8::is_ascii_digit);
            if !(before && after) {
                return None;
            }
        }
    }
    Some(s.replace('_', ""))
}

/// The builtin `sum` of floats. It starts from the int 0, so no terms is
/// the int `0`, and since CPython 3.12 it carries a Neumaier compensation.
pub fn py_sum_floats(terms: &[f64]) -> Json {
    let Some((&first, rest)) = terms.split_first() else {
        return Json::Int(0);
    };
    let mut total = 0.0 + first;
    let mut compensation = 0.0;
    for &x in rest {
        let t = total + x;
        compensation += if total.abs() >= x.abs() {
            (total - t) + x
        } else {
            (x - t) + total
        };
        total = t;
    }
    if compensation != 0.0 && compensation.is_finite() {
        total += compensation;
    }
    Json::Float(total)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn doc() -> Json {
        parse_json(
            r#"{"z": [1, {"b": null, "a": true}], "a": "quote\" back\\ nl\n tab\t del\u007f e\u00e9 \u20ac astral\ud83d\ude00", "m": {}, "l": [], "n": -5, "big": 18446744073709551615}"#,
        )
        .unwrap()
    }

    /// `json.dumps(doc, ...)` as CPython 3.12 printed it on 2026-10-06.
    #[test]
    fn dumps_matches_python_in_every_style_the_merge_uses() {
        let ascii = r#""quote\" back\\ nl\n tab\t del\u007f e\u00e9 \u20ac astral\ud83d\ude00""#;
        assert_eq!(
            dumps(&doc(), Style::Indent(2), true).unwrap(),
            format!("{{\n  \"a\": {ascii},\n  \"big\": 18446744073709551615,\n  \"l\": [],\n  \"m\": {{}},\n  \"n\": -5,\n  \"z\": [\n    1,\n    {{\n      \"a\": true,\n      \"b\": null\n    }}\n  ]\n}}")
        );
        assert_eq!(
            dumps(&doc(), Style::Indent(1), false).unwrap(),
            format!("{{\n \"z\": [\n  1,\n  {{\n   \"b\": null,\n   \"a\": true\n  }}\n ],\n \"a\": {ascii},\n \"m\": {{}},\n \"l\": [],\n \"n\": -5,\n \"big\": 18446744073709551615\n}}")
        );
        assert_eq!(
            dumps(&doc(), Style::Compact, true).unwrap(),
            format!("{{\"a\":{ascii},\"big\":18446744073709551615,\"l\":[],\"m\":{{}},\"n\":-5,\"z\":[1,{{\"a\":true,\"b\":null}}]}}")
        );
    }

    #[test]
    fn objects_behave_like_python_dicts() {
        let mut o = parse_json(r#"{"b": 1, "a": 2, "b": 3}"#).unwrap();
        let Json::Obj(ref mut obj) = o else { panic!() };
        assert_eq!(
            dumps(&Json::Obj(obj.clone()), Style::Compact, false).unwrap(),
            r#"{"b":3,"a":2}"#
        );
        obj.set("a", Json::Int(4));
        obj.set("c", Json::Null);
        assert_eq!(
            dumps(&o, Style::Compact, false).unwrap(),
            r#"{"b":3,"a":4,"c":null}"#
        );
        assert_eq!(
            parse_json(r#"{"x": 1, "y": [true]}"#).unwrap(),
            parse_json(r#"{"y": [1], "x": true}"#).unwrap()
        );
        assert!(dumps(&parse_json("1.5").unwrap(), Style::Compact, false).is_err());
        for bad in [
            "01",
            "[1,]",
            "{\"a\" 1}",
            "\"\u{1}\"",
            "\"\\ud800\"",
            "nul",
            "1 2",
        ] {
            assert!(parse_json(bad).is_err(), "{bad}");
        }
    }

    #[test]
    fn python_text_helpers() {
        assert_eq!(
            py_splitlines("a\r\nb\rc\nd\u{b}e\n"),
            ["a", "b", "c", "d", "e"]
        );
        assert_eq!(py_splitlines("x\n\ny"), ["x", "", "y"]);
        assert_eq!(py_strip("\u{1f} k \n"), "k");
    }

    /// `repr(x)` and `'%.4g' % x` as CPython 3.12 printed them on 2026-10-06.
    #[test]
    fn floats_print_like_cpython() {
        let two_to_60_9 = f64::from_bits(0x43bd_db68_0117_ab0a);
        for (x, want) in [
            (0.0, "0.0"),
            (-0.0, "-0.0"),
            (1.0, "1.0"),
            (0.1, "0.1"),
            (0.30000000000000004, "0.30000000000000004"),
            (1e16, "1e+16"),
            (1e15, "1000000000000000.0"),
            (9999999999999998.0, "9999999999999998.0"),
            (0.0001, "0.0001"),
            (0.00001, "1e-05"),
            (1.5e-7, "1.5e-07"),
            (123.456, "123.456"),
            (6.160384e15, "6160384000000000.0"),
            (two_to_60_9, "2.151427600900885e+18"),
            (14.123e9, "14123000000.0"),
            (1e22, "1e+22"),
            (5e-324, "5e-324"),
            (f64::MAX, "1.7976931348623157e+308"),
            (-2.5, "-2.5"),
            (1.2345678901234568e18, "1.2345678901234568e+18"),
            // Exactly halfway between two shortest strings: the even one.
            (f64::from_bits(0x4309_b670_d8d2_c15a), "904683776006187.2"),
            (f64::from_bits(0xc2d0_a186_0217_1168), "-73143696251973.62"),
            (f64::INFINITY, "inf"),
            (f64::NEG_INFINITY, "-inf"),
            (f64::NAN, "nan"),
        ] {
            assert_eq!(float_repr(x), want, "{x:e}");
        }
        for (x, want) in [
            (0.0, "0"),
            (-0.0, "-0"),
            (6.160384e15, "6.16e+15"),
            (12345678.0, "1.235e+07"),
            (0.00012345, "0.0001234"),
            (0.000012345, "1.234e-05"),
            (9.9995, "9.999"),
            (1000.0, "1000"),
            (1e-5, "1e-05"),
            (123.0, "123"),
            (1.5, "1.5"),
            (two_to_60_9, "2.151e+18"),
            (99995.0, "1e+05"),
        ] {
            assert_eq!(format_g(x, 4), want, "{x:e}");
        }
    }

    /// `json.dumps(doc, indent=1)` as CPython 3.12 printed it.
    #[test]
    fn floats_are_written_like_json_dumps() {
        let doc = Json::Obj(
            Obj::new()
                .with("a", Json::Float(f64::INFINITY))
                .with("b", Json::Float(1e16))
                .with("c", Json::List(vec![Json::Float(0.5), Json::Float(-0.0)]))
                .with("d", Json::Float(f64::NAN))
                .with("e", Json::Float(f64::NEG_INFINITY))
                .with("f", Json::Int(3)),
        );
        assert_eq!(
            dumps_floats(&doc, Style::Indent(1), false),
            "{\n \"a\": Infinity,\n \"b\": 1e+16,\n \"c\": [\n  0.5,\n  -0.0\n ],\n \"d\": NaN,\n \"e\": -Infinity,\n \"f\": 3\n}"
        );
        assert!(dumps(&doc, Style::Indent(1), false).is_err());
        let back = parse_json(&dumps_floats(&doc, Style::Compact, false)).unwrap();
        let Json::Obj(back) = back else { panic!() };
        assert_eq!(back.get("a"), Some(&Json::Float(f64::INFINITY)));
        assert_eq!(back.get("e"), Some(&Json::Float(f64::NEG_INFINITY)));
        assert!(matches!(back.get("d"), Some(Json::Float(x)) if x.is_nan()));
        assert!(parse_json("-Inf").is_err());
    }

    #[test]
    fn int_and_float_read_values_like_python() {
        for (text, want) in [
            ("12", Some(12)),
            (" 7 ", Some(7)),
            ("+3", Some(3)),
            ("-0", Some(0)),
            ("1_000", Some(1000)),
            ("007", Some(7)),
            ("\u{a0}9\u{2003}", Some(9)),
            ("\u{85}9", Some(9)),
            ("\u{1c}1", None),
            ("1__0", None),
            ("_1", None),
            ("1_", None),
            ("", None),
            ("1.5", None),
            ("+-3", None),
            ("abc", None),
        ] {
            assert_eq!(py_int(&Json::str(text)).ok(), want, "{text:?}");
        }
        for (text, want) in [
            ("1.5", Some(1.5)),
            (" 2 ", Some(2.0)),
            ("\u{b}2.5\u{c}", Some(2.5)),
            ("1\u{1d}", None),
            ("1e3", Some(1000.0)),
            ("1e+3", Some(1000.0)),
            ("inf", Some(f64::INFINITY)),
            ("-Infinity", Some(f64::NEG_INFINITY)),
            ("1_0.5", Some(10.5)),
            (".5", Some(0.5)),
            ("1.", Some(1.0)),
            ("1._5", None),
            ("--1", None),
            ("abc", None),
            ("", None),
            ("0x10", None),
        ] {
            assert_eq!(py_float(&Json::str(text)).ok(), want, "{text:?}");
        }
        assert!(py_float(&Json::str("nan")).unwrap().is_nan());
        assert_eq!(py_int(&Json::Float(1.9)), Ok(1));
        assert_eq!(py_int(&Json::Float(-1.9)), Ok(-1));
        assert_eq!(
            py_int(&Json::Float(1e30)),
            Ok(1_000_000_000_000_000_019_884_624_838_656)
        );
        assert_eq!(py_int(&Json::Bool(true)), Ok(1));
        for (v, error) in [
            (
                Json::str("abc"),
                "invalid literal for int() with base 10: 'abc'",
            ),
            (
                Json::Float(f64::INFINITY),
                "cannot convert float infinity to integer",
            ),
            (Json::Float(f64::NAN), "cannot convert float NaN to integer"),
            (
                Json::List(vec![Json::Int(1)]),
                "int() argument must be a string, a bytes-like object or a real number, not 'list'",
            ),
        ] {
            assert_eq!(py_int(&v).unwrap_err(), error);
        }
        assert_eq!(
            py_float(&Json::str("abc")).unwrap_err(),
            "could not convert string to float: 'abc'"
        );
        assert_eq!(
            py_float(&Json::Null).unwrap_err(),
            "float() argument must be a string or a real number, not 'NoneType'"
        );
    }

    /// `sum(terms)` as CPython 3.12 returned it.
    #[test]
    fn sum_is_the_compensated_builtin() {
        assert!(matches!(py_sum_floats(&[]), Json::Int(0)));
        let total = |terms: &[f64]| match py_sum_floats(terms) {
            Json::Float(x) => float_repr(x),
            other => panic!("{other:?}"),
        };
        assert_eq!(total(&[0.1; 10]), "1.0");
        assert_eq!(total(&[1e100, 1.0, -1e100]), "1.0");
        assert_eq!(total(&[14.1e9, 13.9e9, 0.1, 0.2]), "28000000000.3");
        assert_eq!(total(&[-0.0]), "0.0");
    }

    #[test]
    fn repr_and_str_follow_python() {
        let v = Json::List(vec![
            Json::str("a"),
            Json::str("it's"),
            Json::str("say \"hi\""),
            Json::str("both ' and \""),
            Json::str("tab\tnl\n\u{7f}\u{0}"),
            Json::str("é€😀"),
            Json::str("\u{a0}\u{ad}\u{200b}\u{2028}"),
            Json::Null,
            Json::Bool(true),
            Json::Float(1.5),
            Json::Obj(Obj::new().with("k", Json::List(vec![Json::Int(1), Json::str("v")]))),
        ]);
        assert_eq!(
            py_repr(&v),
            r#"['a', "it's", 'say "hi"', 'both \' and "', 'tab\tnl\n\x7f\x00', 'é€😀', '\xa0\xad\u200b\u2028', None, True, 1.5, {'k': [1, 'v']}]"#
        );
        assert_eq!(py_str(&Json::str("x'y")), "x'y");
        assert_eq!(py_str(&Json::Null), "None");
        assert_eq!(py_str(&Json::Float(2.0)), "2.0");
    }
}

//! JSON as the programme's records hold it: objects keep their key order,
//! integers stay integers, and output is formatted as the earlier rounds'
//! analysis files were (one-space indent, floats in their shortest
//! round-trip form with an exponent below `1e-4` and from `1e16`, non-ASCII
//! escaped).  A record read and written back is therefore unchanged, and a
//! frozen analysis can be compared byte for byte with a native one.

use std::fmt::Write as _;
use std::path::Path;

#[derive(Clone, Debug, PartialEq)]
pub enum J {
    Null,
    Bool(bool),
    Int(i128),
    Float(f64),
    Str(String),
    Arr(Vec<J>),
    Obj(Vec<(String, J)>),
}

impl J {
    pub fn get(&self, key: &str) -> Option<&J> {
        match self {
            J::Obj(kv) => kv.iter().find(|(k, _)| k == key).map(|(_, v)| v),
            _ => None,
        }
    }

    /// `self[key]`, or an error naming the key.
    pub fn at(&self, key: &str) -> Result<&J, String> {
        self.get(key).ok_or_else(|| format!("no key `{key}`"))
    }

    pub fn as_f64(&self) -> Option<f64> {
        match self {
            J::Int(i) => Some(*i as f64),
            J::Float(x) => Some(*x),
            _ => None,
        }
    }

    pub fn as_i128(&self) -> Option<i128> {
        match self {
            J::Int(i) => Some(*i),
            _ => None,
        }
    }

    pub fn as_str(&self) -> Option<&str> {
        match self {
            J::Str(s) => Some(s),
            _ => None,
        }
    }

    pub fn as_bool(&self) -> Option<bool> {
        match self {
            J::Bool(b) => Some(*b),
            _ => None,
        }
    }

    pub fn as_arr(&self) -> Option<&[J]> {
        match self {
            J::Arr(v) => Some(v),
            _ => None,
        }
    }

    pub fn as_obj(&self) -> Option<&[(String, J)]> {
        match self {
            J::Obj(kv) => Some(kv),
            _ => None,
        }
    }

    /// Python's truth value, as `bool(x)` reads a record's field.
    pub fn truthy(&self) -> bool {
        match self {
            J::Null => false,
            J::Bool(b) => *b,
            J::Int(i) => *i != 0,
            J::Float(x) => *x != 0.0,
            J::Str(s) => !s.is_empty(),
            J::Arr(v) => !v.is_empty(),
            J::Obj(kv) => !kv.is_empty(),
        }
    }
}

/// An object built in order: `obj([("a", J::Int(1)), ...])`.
pub fn obj<const N: usize>(kv: [(&str, J); N]) -> J {
    J::Obj(kv.into_iter().map(|(k, v)| (k.to_string(), v)).collect())
}

pub fn opt_f64(x: Option<f64>) -> J {
    x.map_or(J::Null, J::Float)
}

// ── reading ─────────────────────────────────────────────────────────

pub fn parse(text: &str) -> Result<J, String> {
    let mut p = Parser {
        s: text.as_bytes(),
        i: 0,
    };
    p.ws();
    let v = p.value()?;
    p.ws();
    if p.i != p.s.len() {
        return Err(format!("trailing characters at byte {}", p.i));
    }
    Ok(v)
}

pub fn read(path: &Path) -> Result<J, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    parse(&text).map_err(|e| format!("{}: {e}", path.display()))
}

/// The record at `path`, or `None` when there is no file.
pub fn read_opt(path: &Path) -> Result<Option<J>, String> {
    if path.exists() {
        read(path).map(Some)
    } else {
        Ok(None)
    }
}

struct Parser<'a> {
    s: &'a [u8],
    i: usize,
}

impl Parser<'_> {
    fn ws(&mut self) {
        while self.i < self.s.len() && matches!(self.s[self.i], b' ' | b'\t' | b'\n' | b'\r') {
            self.i += 1;
        }
    }

    fn err<T>(&self, what: &str) -> Result<T, String> {
        Err(format!("{what} at byte {}", self.i))
    }

    fn eat(&mut self, lit: &str) -> bool {
        if self.s[self.i..].starts_with(lit.as_bytes()) {
            self.i += lit.len();
            true
        } else {
            false
        }
    }

    fn value(&mut self) -> Result<J, String> {
        match self.s.get(self.i) {
            None => self.err("unexpected end"),
            Some(b'{') => {
                self.i += 1;
                let mut kv = Vec::new();
                self.ws();
                if self.eat("}") {
                    return Ok(J::Obj(kv));
                }
                loop {
                    self.ws();
                    if self.s.get(self.i) != Some(&b'"') {
                        return self.err("expected a key");
                    }
                    let k = self.string()?;
                    self.ws();
                    if !self.eat(":") {
                        return self.err("expected `:`");
                    }
                    self.ws();
                    let v = self.value()?;
                    kv.push((k, v));
                    self.ws();
                    if self.eat(",") {
                        continue;
                    }
                    if self.eat("}") {
                        return Ok(J::Obj(kv));
                    }
                    return self.err("expected `,` or `}`");
                }
            }
            Some(b'[') => {
                self.i += 1;
                let mut v = Vec::new();
                self.ws();
                if self.eat("]") {
                    return Ok(J::Arr(v));
                }
                loop {
                    self.ws();
                    v.push(self.value()?);
                    self.ws();
                    if self.eat(",") {
                        continue;
                    }
                    if self.eat("]") {
                        return Ok(J::Arr(v));
                    }
                    return self.err("expected `,` or `]`");
                }
            }
            Some(b'"') => self.string().map(J::Str),
            Some(b't') if self.eat("true") => Ok(J::Bool(true)),
            Some(b'f') if self.eat("false") => Ok(J::Bool(false)),
            Some(b'n') if self.eat("null") => Ok(J::Null),
            Some(b'N') if self.eat("NaN") => Ok(J::Float(f64::NAN)),
            Some(b'I') if self.eat("Infinity") => Ok(J::Float(f64::INFINITY)),
            Some(b'-') if self.eat("-Infinity") => Ok(J::Float(f64::NEG_INFINITY)),
            Some(_) => self.number(),
        }
    }

    fn number(&mut self) -> Result<J, String> {
        let start = self.i;
        let mut float = false;
        while let Some(&c) = self.s.get(self.i) {
            match c {
                b'0'..=b'9' | b'-' | b'+' => {}
                b'.' | b'e' | b'E' => float = true,
                _ => break,
            }
            self.i += 1;
        }
        let text = std::str::from_utf8(&self.s[start..self.i]).expect("ASCII");
        if text.is_empty() {
            return self.err("unexpected character");
        }
        if float {
            text.parse::<f64>()
                .map(J::Float)
                .map_err(|e| format!("number `{text}`: {e}"))
        } else {
            text.parse::<i128>()
                .map(J::Int)
                .map_err(|e| format!("integer `{text}`: {e}"))
        }
    }

    fn hex4(&mut self) -> Result<u32, String> {
        let h = self
            .s
            .get(self.i..self.i + 4)
            .and_then(|b| std::str::from_utf8(b).ok())
            .and_then(|t| u32::from_str_radix(t, 16).ok());
        match h {
            Some(v) => {
                self.i += 4;
                Ok(v)
            }
            None => self.err("bad \\u escape"),
        }
    }

    fn string(&mut self) -> Result<String, String> {
        self.i += 1; // the opening quote
        let mut out = String::new();
        loop {
            let start = self.i;
            while let Some(&c) = self.s.get(self.i) {
                if c == b'"' || c == b'\\' {
                    break;
                }
                self.i += 1;
            }
            out.push_str(std::str::from_utf8(&self.s[start..self.i]).map_err(|e| e.to_string())?);
            match self.s.get(self.i) {
                None => return self.err("unterminated string"),
                Some(b'"') => {
                    self.i += 1;
                    return Ok(out);
                }
                Some(_) => {
                    self.i += 1;
                    let c = *self.s.get(self.i).ok_or("unterminated escape")?;
                    self.i += 1;
                    match c {
                        b'"' => out.push('"'),
                        b'\\' => out.push('\\'),
                        b'/' => out.push('/'),
                        b'b' => out.push('\u{8}'),
                        b'f' => out.push('\u{c}'),
                        b'n' => out.push('\n'),
                        b'r' => out.push('\r'),
                        b't' => out.push('\t'),
                        b'u' => {
                            let hi = self.hex4()?;
                            let cp = if (0xD800..0xDC00).contains(&hi) && self.eat("\\u") {
                                let lo = self.hex4()?;
                                0x10000 + ((hi - 0xD800) << 10) + (lo - 0xDC00)
                            } else {
                                hi
                            };
                            out.push(char::from_u32(cp).ok_or("bad code point")?);
                        }
                        _ => return self.err("bad escape"),
                    }
                }
            }
        }
    }
}

// ── writing ─────────────────────────────────────────────────────────

/// A float as Python's `repr` writes it: the shortest digits that read
/// back to the same value, in positional form from `1e-4` up to `1e16`
/// and with an exponent outside it.
pub fn py_float(x: f64) -> String {
    if x.is_nan() {
        return "NaN".into();
    }
    if x.is_infinite() {
        return if x > 0.0 { "Infinity" } else { "-Infinity" }.into();
    }
    if x == 0.0 {
        return if x.is_sign_negative() { "-0.0" } else { "0.0" }.into();
    }
    // `{:e}` gives the shortest round-trip digits as d.ddd…e±x.
    let sci = format!("{:e}", x.abs());
    let (mant, exp) = sci.split_once('e').expect("exponent");
    let exp: i32 = exp.parse().expect("exponent digits");
    let digits: String = mant.chars().filter(|c| *c != '.').collect();
    let decpt = exp + 1; // x = 0.d1d2… × 10^decpt
    let sign = if x < 0.0 { "-" } else { "" };
    let nd = digits.len() as i32;
    if decpt <= -4 || decpt > 16 {
        let m = if digits.len() == 1 {
            digits.clone()
        } else {
            format!("{}.{}", &digits[..1], &digits[1..])
        };
        let e = decpt - 1;
        format!(
            "{sign}{m}e{}{:02}",
            if e < 0 { '-' } else { '+' },
            e.unsigned_abs()
        )
    } else if decpt <= 0 {
        format!("{sign}0.{}{digits}", "0".repeat((-decpt) as usize))
    } else if decpt < nd {
        format!(
            "{sign}{}.{}",
            &digits[..decpt as usize],
            &digits[decpt as usize..]
        )
    } else {
        format!("{sign}{digits}{}.0", "0".repeat((decpt - nd) as usize))
    }
}

fn write_str(out: &mut String, s: &str) {
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
                let mut buf = [0u16; 2];
                for unit in c.encode_utf16(&mut buf) {
                    let _ = write!(out, "\\u{unit:04x}");
                }
            }
        }
    }
    out.push('"');
}

fn write(out: &mut String, v: &J, indent: usize, depth: usize) {
    let pad = |out: &mut String, d: usize| {
        out.push('\n');
        out.push_str(&" ".repeat(indent * d));
    };
    match v {
        J::Null => out.push_str("null"),
        J::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
        J::Int(i) => {
            let _ = write!(out, "{i}");
        }
        J::Float(x) => out.push_str(&py_float(*x)),
        J::Str(s) => write_str(out, s),
        J::Arr(items) if items.is_empty() => out.push_str("[]"),
        J::Arr(items) => {
            out.push('[');
            for (k, item) in items.iter().enumerate() {
                if k > 0 {
                    out.push(',');
                }
                pad(out, depth + 1);
                write(out, item, indent, depth + 1);
            }
            pad(out, depth);
            out.push(']');
        }
        J::Obj(kv) if kv.is_empty() => out.push_str("{}"),
        J::Obj(kv) => {
            out.push('{');
            for (k, (key, item)) in kv.iter().enumerate() {
                if k > 0 {
                    out.push(',');
                }
                pad(out, depth + 1);
                write_str(out, key);
                out.push_str(": ");
                write(out, item, indent, depth + 1);
            }
            pad(out, depth);
            out.push('}');
        }
    }
}

/// `v` as the earlier rounds wrote their records: `indent` spaces a level.
pub fn dumps(v: &J, indent: usize) -> String {
    let mut out = String::new();
    write(&mut out, v, indent, 0);
    out
}

/// Every object's keys in sorted order, at every depth.
pub fn sorted(v: &J) -> J {
    match v {
        J::Arr(items) => J::Arr(items.iter().map(sorted).collect()),
        J::Obj(kv) => {
            let mut kv: Vec<(String, J)> = kv.iter().map(|(k, v)| (k.clone(), sorted(v))).collect();
            kv.sort_by(|a, b| a.0.cmp(&b.0));
            J::Obj(kv)
        }
        other => other.clone(),
    }
}

/// `v` on one line with `", "` and `": "` between items, as a JSON lines
/// record is written; with `sort_keys`, every object's keys sorted.
pub fn dumps_line(v: &J, sort_keys: bool) -> String {
    fn line(out: &mut String, v: &J) {
        match v {
            J::Arr(items) => {
                out.push('[');
                for (k, item) in items.iter().enumerate() {
                    if k > 0 {
                        out.push_str(", ");
                    }
                    line(out, item);
                }
                out.push(']');
            }
            J::Obj(kv) => {
                out.push('{');
                for (k, (key, item)) in kv.iter().enumerate() {
                    if k > 0 {
                        out.push_str(", ");
                    }
                    write_str(out, key);
                    out.push_str(": ");
                    line(out, item);
                }
                out.push('}');
            }
            scalar => write(out, scalar, 0, 0),
        }
    }
    let mut out = String::new();
    if sort_keys {
        line(&mut out, &sorted(v));
    } else {
        line(&mut out, v);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn floats_print_as_python_repr_does() {
        for (x, want) in [
            (1.0, "1.0"),
            (0.1, "0.1"),
            (120.561, "120.561"),
            (1e-5, "1e-05"),
            (0.0001, "0.0001"),
            (1.5e-7, "1.5e-07"),
            (1e16, "1e+16"),
            (1234567890123456.0, "1234567890123456.0"),
            (12345678901234567.0, "1.2345678901234568e+16"),
            (-2.5, "-2.5"),
            (1.0195642272162666, "1.0195642272162666"),
            (1e22, "1e+22"),
            (123.0, "123.0"),
            (0.001, "0.001"),
        ] {
            assert_eq!(py_float(x), want, "{x:e}");
        }
    }

    #[test]
    fn a_record_reads_and_writes_back_unchanged() {
        let text = "{\n \"b\": 1,\n \"a\": [\n  1.5,\n  \"x\\u00e9\\n\",\n  {}\n ],\n \"c\": [],\n \"d\": null,\n \"e\": true,\n \"f\": 1e-05\n}";
        let v = parse(text).unwrap();
        assert_eq!(dumps(&v, 1), text);
        assert_eq!(v.get("b"), Some(&J::Int(1)));
    }

    #[test]
    fn non_bmp_characters_escape_as_surrogate_pairs() {
        let v = J::Str("\u{1F600}".into());
        assert_eq!(dumps(&v, 1), "\"\\ud83d\\ude00\"");
        assert_eq!(parse("\"\\ud83d\\ude00\"").unwrap(), v);
    }
}

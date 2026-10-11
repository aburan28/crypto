//! Minimal JSON values, serialisation and parsing for the command-line records and the
//! benchmark files (the crate has no dependencies). Object keys keep insertion order; integers
//! that may exceed 2^53 are written by the callers as decimal strings. Parsed numbers keep their
//! source text (`Json::Raw`), so re-serialising reproduces them exactly.
use std::fmt::Write;

#[derive(Clone, Debug, PartialEq)]
pub enum Json {
    Null,
    Bool(bool),
    Num(i64),
    Float(f64),
    Str(String),
    /// a number as written in parsed input
    Raw(String),
    Arr(Vec<Json>),
    Obj(Vec<(String, Json)>),
}

impl Json {
    /// An object from (key, value) pairs.
    pub fn obj<K: Into<String>>(pairs: Vec<(K, Json)>) -> Json {
        Json::Obj(pairs.into_iter().map(|(k, v)| (k.into(), v)).collect())
    }
    pub fn str<S: ToString>(s: S) -> Json {
        Json::Str(s.to_string())
    }
    /// Compact serialisation.
    pub fn dump(&self) -> String {
        let mut out = String::new();
        self.write(&mut out);
        out
    }
    fn write(&self, out: &mut String) {
        match self {
            Json::Null => out.push_str("null"),
            Json::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
            Json::Num(n) => {
                let _ = write!(out, "{n}");
            }
            Json::Float(x) => {
                if x.is_finite() {
                    let _ = write!(out, "{x}");
                } else {
                    out.push_str("null");
                }
            }
            Json::Str(s) => escape(s, out),
            Json::Raw(r) => out.push_str(r),
            Json::Arr(v) => {
                out.push('[');
                for (i, x) in v.iter().enumerate() {
                    if i > 0 {
                        out.push(',');
                    }
                    x.write(out);
                }
                out.push(']');
            }
            Json::Obj(kv) => {
                out.push('{');
                for (i, (k, v)) in kv.iter().enumerate() {
                    if i > 0 {
                        out.push(',');
                    }
                    escape(k, out);
                    out.push(':');
                    v.write(out);
                }
                out.push('}');
            }
        }
    }
}

fn escape(s: &str, out: &mut String) {
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out.push('"');
}

impl Json {
    /// Parse one JSON value (the whole input).
    pub fn parse(s: &str) -> Result<Json, String> {
        let b = s.as_bytes();
        let mut i = 0;
        let v = value(b, &mut i)?;
        ws(b, &mut i);
        if i != b.len() {
            return Err(format!("trailing data at byte {i}"));
        }
        Ok(v)
    }
    /// Field of an object.
    pub fn get(&self, k: &str) -> Option<&Json> {
        match self {
            Json::Obj(kv) => kv.iter().find(|(key, _)| key == k).map(|(_, v)| v),
            _ => None,
        }
    }
    /// Numeric value of a number (parsed or built).
    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Json::Num(n) => Some(*n as f64),
            Json::Float(x) => Some(*x),
            Json::Raw(r) => r.parse().ok(),
            _ => None,
        }
    }
    pub fn as_str(&self) -> Option<&str> {
        match self {
            Json::Str(s) => Some(s),
            _ => None,
        }
    }
}

fn ws(b: &[u8], i: &mut usize) {
    while *i < b.len() && matches!(b[*i], b' ' | b'\t' | b'\n' | b'\r') {
        *i += 1;
    }
}

fn value(b: &[u8], i: &mut usize) -> Result<Json, String> {
    ws(b, i);
    match b.get(*i) {
        Some(b'{') => {
            *i += 1;
            let mut kv = vec![];
            ws(b, i);
            if b.get(*i) == Some(&b'}') {
                *i += 1;
                return Ok(Json::Obj(kv));
            }
            loop {
                ws(b, i);
                let Json::Str(k) = value(b, i)? else {
                    return Err("object key must be a string".into());
                };
                ws(b, i);
                if b.get(*i) != Some(&b':') {
                    return Err(format!("expected ':' at byte {i}"));
                }
                *i += 1;
                let v = value(b, i)?;
                kv.push((k, v));
                ws(b, i);
                match b.get(*i) {
                    Some(b',') => *i += 1,
                    Some(b'}') => {
                        *i += 1;
                        return Ok(Json::Obj(kv));
                    }
                    _ => return Err(format!("expected ',' or '}}' at byte {i}")),
                }
            }
        }
        Some(b'[') => {
            *i += 1;
            let mut v = vec![];
            ws(b, i);
            if b.get(*i) == Some(&b']') {
                *i += 1;
                return Ok(Json::Arr(v));
            }
            loop {
                v.push(value(b, i)?);
                ws(b, i);
                match b.get(*i) {
                    Some(b',') => *i += 1,
                    Some(b']') => {
                        *i += 1;
                        return Ok(Json::Arr(v));
                    }
                    _ => return Err(format!("expected ',' or ']' at byte {i}")),
                }
            }
        }
        Some(b'"') => {
            *i += 1;
            let mut s = String::new();
            loop {
                match b.get(*i) {
                    None => return Err("unterminated string".into()),
                    Some(b'"') => {
                        *i += 1;
                        return Ok(Json::Str(s));
                    }
                    Some(b'\\') => {
                        let c = *b.get(*i + 1).ok_or("bad escape")?;
                        *i += 2;
                        match c {
                            b'"' => s.push('"'),
                            b'\\' => s.push('\\'),
                            b'/' => s.push('/'),
                            b'n' => s.push('\n'),
                            b't' => s.push('\t'),
                            b'r' => s.push('\r'),
                            b'b' => s.push('\u{8}'),
                            b'f' => s.push('\u{c}'),
                            b'u' => {
                                let h =
                                    std::str::from_utf8(b.get(*i..*i + 4).ok_or("bad \\u escape")?)
                                        .map_err(|e| e.to_string())?;
                                let cp = u32::from_str_radix(h, 16).map_err(|e| e.to_string())?;
                                s.push(char::from_u32(cp).unwrap_or('\u{fffd}'));
                                *i += 4;
                            }
                            _ => return Err("bad escape".into()),
                        }
                    }
                    Some(&lead) => {
                        // copy one UTF-8 character (the input is a &str, so it is well formed)
                        let len = match lead {
                            0..=0x7f => 1,
                            0xc0..=0xdf => 2,
                            0xe0..=0xef => 3,
                            _ => 4,
                        };
                        let ch = std::str::from_utf8(b.get(*i..*i + len).ok_or("truncated UTF-8")?)
                            .map_err(|e| e.to_string())?;
                        s.push_str(ch);
                        *i += len;
                    }
                }
            }
        }
        Some(b't') if b[*i..].starts_with(b"true") => {
            *i += 4;
            Ok(Json::Bool(true))
        }
        Some(b'f') if b[*i..].starts_with(b"false") => {
            *i += 5;
            Ok(Json::Bool(false))
        }
        Some(b'n') if b[*i..].starts_with(b"null") => {
            *i += 4;
            Ok(Json::Null)
        }
        Some(c) if *c == b'-' || c.is_ascii_digit() => {
            let start = *i;
            while *i < b.len() && matches!(b[*i], b'-' | b'+' | b'.' | b'e' | b'E' | b'0'..=b'9') {
                *i += 1;
            }
            Ok(Json::Raw(
                std::str::from_utf8(&b[start..*i]).unwrap().to_string(),
            ))
        }
        _ => Err(format!("unexpected input at byte {i}")),
    }
}

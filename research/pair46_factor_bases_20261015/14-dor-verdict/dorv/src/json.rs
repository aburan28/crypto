//! Minimal JSON reader/writer (no dependencies).

#[derive(Clone, Debug)]
pub enum V {
    Null,
    Bool(bool),
    Num(f64),
    Str(String),
    Arr(Vec<V>),
    Obj(Vec<(String, V)>),
}

impl V {
    pub fn get(&self, k: &str) -> Option<&V> {
        match self {
            V::Obj(kv) => kv.iter().find(|(a, _)| a == k).map(|(_, v)| v),
            _ => None,
        }
    }
    pub fn num(&self) -> Option<f64> {
        if let V::Num(x) = self { Some(*x) } else { None }
    }
    pub fn str(&self) -> Option<&str> {
        if let V::Str(s) = self { Some(s) } else { None }
    }
    pub fn arr(&self) -> Option<&Vec<V>> {
        if let V::Arr(a) = self { Some(a) } else { None }
    }
    pub fn obj(&self) -> Option<&Vec<(String, V)>> {
        if let V::Obj(a) = self { Some(a) } else { None }
    }
}

pub fn parse(s: &str) -> V {
    let b = s.as_bytes();
    let mut i = 0;
    let v = val(b, &mut i);
    v
}

fn ws(b: &[u8], i: &mut usize) {
    while *i < b.len() && (b[*i] as char).is_whitespace() {
        *i += 1;
    }
}

fn val(b: &[u8], i: &mut usize) -> V {
    ws(b, i);
    match b[*i] {
        b'{' => {
            *i += 1;
            let mut kv = Vec::new();
            loop {
                ws(b, i);
                if b[*i] == b'}' { *i += 1; break; }
                let k = match val(b, i) { V::Str(s) => s, _ => panic!("key") };
                ws(b, i);
                assert_eq!(b[*i], b':');
                *i += 1;
                let v = val(b, i);
                kv.push((k, v));
                ws(b, i);
                if b[*i] == b',' { *i += 1; }
            }
            V::Obj(kv)
        }
        b'[' => {
            *i += 1;
            let mut a = Vec::new();
            loop {
                ws(b, i);
                if b[*i] == b']' { *i += 1; break; }
                a.push(val(b, i));
                ws(b, i);
                if b[*i] == b',' { *i += 1; }
            }
            V::Arr(a)
        }
        b'"' => {
            *i += 1;
            let mut s = String::new();
            while b[*i] != b'"' {
                if b[*i] == b'\\' {
                    *i += 1;
                    match b[*i] {
                        b'n' => s.push('\n'),
                        b't' => s.push('\t'),
                        b'u' => {
                            let h = std::str::from_utf8(&b[*i + 1..*i + 5]).unwrap();
                            s.push(char::from_u32(u32::from_str_radix(h, 16).unwrap()).unwrap_or('?'));
                            *i += 4;
                        }
                        c => s.push(c as char),
                    }
                    *i += 1;
                } else {
                    // utf8 passthrough
                    let start = *i;
                    while b[*i] != b'"' && b[*i] != b'\\' { *i += 1; }
                    s.push_str(std::str::from_utf8(&b[start..*i]).unwrap());
                }
            }
            *i += 1;
            V::Str(s)
        }
        b't' => { *i += 4; V::Bool(true) }
        b'f' => { *i += 5; V::Bool(false) }
        b'n' => { *i += 4; V::Null }
        _ => {
            let start = *i;
            while *i < b.len() && (b"+-.eE0123456789".contains(&b[*i]) || b[*i].is_ascii_alphabetic()) { *i += 1; }
            let t = std::str::from_utf8(&b[start..*i]).unwrap();
            V::Num(t.parse::<f64>().unwrap_or(f64::NAN))
        }
    }
}

pub fn dump(v: &V) -> String {
    let mut s = String::new();
    w(v, &mut s, 0);
    s
}

fn esc(x: &str) -> String {
    let mut o = String::from("\"");
    for c in x.chars() {
        match c {
            '"' => o.push_str("\\\""),
            '\\' => o.push_str("\\\\"),
            '\n' => o.push_str("\\n"),
            c => o.push(c),
        }
    }
    o.push('"');
    o
}

fn w(v: &V, s: &mut String, ind: usize) {
    let pad = |n: usize| " ".repeat(n);
    match v {
        V::Null => s.push_str("null"),
        V::Bool(b) => s.push_str(if *b { "true" } else { "false" }),
        V::Num(x) => {
            if !x.is_finite() { s.push_str("null") }
            else if x.fract() == 0.0 && x.abs() < 1e15 { s.push_str(&format!("{}", *x as i64)) }
            else { s.push_str(&format!("{}", x)) }
        }
        V::Str(x) => s.push_str(&esc(x)),
        V::Arr(a) => {
            let simple = a.iter().all(|e| matches!(e, V::Num(_) | V::Str(_) | V::Null | V::Bool(_)));
            if simple {
                s.push('[');
                for (k, e) in a.iter().enumerate() { if k > 0 { s.push_str(", "); } w(e, s, ind); }
                s.push(']');
            } else {
                s.push_str("[\n");
                for (k, e) in a.iter().enumerate() {
                    s.push_str(&pad(ind + 1));
                    w(e, s, ind + 1);
                    if k + 1 < a.len() { s.push(','); }
                    s.push('\n');
                }
                s.push_str(&pad(ind)); s.push(']');
            }
        }
        V::Obj(kv) => {
            s.push_str("{\n");
            for (k, (key, e)) in kv.iter().enumerate() {
                s.push_str(&pad(ind + 1));
                s.push_str(&esc(key)); s.push_str(": ");
                w(e, s, ind + 1);
                if k + 1 < kv.len() { s.push(','); }
                s.push('\n');
            }
            s.push_str(&pad(ind)); s.push('}');
        }
    }
}

pub fn n(x: f64) -> V { V::Num(x) }
pub fn st(x: &str) -> V { V::Str(x.to_string()) }
pub fn obj(kv: Vec<(&str, V)>) -> V { V::Obj(kv.into_iter().map(|(k, v)| (k.to_string(), v)).collect()) }

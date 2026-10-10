enum Op {
    Text(f64, f64, f64, String),
    Rect(f64, f64, f64, f64, String),
    Line(f64, f64, f64, f64, String),
    Dot(f64, f64, f64, String),
}
fn xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}
fn svg(ops: &[Op]) -> String {
    let mut s=String::from("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1100\" height=\"760\" viewBox=\"0 0 1100 760\"><rect width=\"1100\" height=\"760\" fill=\"white\"/><g font-family=\"Arial,Helvetica,sans-serif\" fill=\"#172637\">");
    for op in ops {
        match op {
        Op::Text(x,y,f,t)=>write!(s,"<text x=\"{x}\" y=\"{y}\" font-size=\"{f}\">{}</text>",xml(t)).unwrap(),
        Op::Rect(x,y,w,h,c)=>write!(s,"<rect x=\"{x}\" y=\"{y}\" width=\"{w}\" height=\"{h}\" rx=\"9\" fill=\"{c}\"/>").unwrap(),
        Op::Line(x,y,a,b,c)=>write!(s,"<line x1=\"{x}\" y1=\"{y}\" x2=\"{a}\" y2=\"{b}\" stroke=\"{c}\" stroke-width=\"2\"/>").unwrap(),
        Op::Dot(x,y,r,c)=>write!(s,"<circle cx=\"{x}\" cy=\"{y}\" r=\"{r}\" fill=\"{c}\"/>").unwrap(),
    }
    }
    s.push_str("</g></svg>\n");
    s
}
fn escape(s: &str) -> String {
    s.replace('\\', "\\\\")
        .replace('(', "\\(")
        .replace(')', "\\)")
}
fn text(s: &mut String, x: f64, y: f64, size: f64, t: &str) {
    writeln!(
        s,
        "0.08 0.14 0.2 rg BT /F1 {size} Tf {x} {y} Td ({}) Tj ET",
        escape(t)
    )
    .unwrap();
}
fn wrapped(s: &mut String, y: &mut f64, t: &str, size: f64) {
    let mut line = String::new();
    for word in t.split_whitespace() {
        if line.len() + word.len() > 96 {
            text(s, 42., *y, size, &line);
            *y -= 16.;
            line.clear();
        }
        if !line.is_empty() {
            line.push(' ');
        }
        line.push_str(word);
    }
    if !line.is_empty() {
        text(s, 42., *y, size, &line);
        *y -= 16.;
    }
    *y -= 8.;
}
fn rgb(c: &str) -> (f64, f64, f64) {
    let n = u32::from_str_radix(&c[1..], 16).unwrap();
    (
        ((n >> 16) & 255) as f64 / 255.,
        ((n >> 8) & 255) as f64 / 255.,
        (n & 255) as f64 / 255.,
    )
}

fn emit_pdf(path: &Path, pages: Vec<(u32, u32, String)>) {
    let mut objects: Vec<String> = vec![
        String::new(),
        "<< /Type /Catalog /Pages 2 0 R >>".into(),
        String::new(),
        "<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>".into(),
    ];
    let mut kids = vec![];
    for (w, h, body) in pages {
        let n = objects.len();
        kids.push(format!("{n} 0 R"));
        objects.push(format!("<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {w} {h}] /Resources << /Font << /F1 3 0 R >> >> /Contents {} 0 R >>",n+1));
        objects.push(format!(
            "<< /Length {} >>\nstream\n{}endstream",
            body.len(),
            body
        ));
    }
    objects[2] = format!(
        "<< /Type /Pages /Count {} /Kids [{}] >>",
        kids.len(),
        kids.join(" ")
    );
    let mut bytes = String::from("%PDF-1.4\n");
    let mut offsets = vec![0];
    for (i, o) in objects.iter().enumerate().skip(1) {
        offsets.push(bytes.len());
        writeln!(bytes, "{i} 0 obj\n{o}\nendobj").unwrap();
    }
    let xref = bytes.len();
    writeln!(bytes, "xref\n0 {}\n0000000000 65535 f ", objects.len()).unwrap();
    for offset in offsets.iter().skip(1) {
        writeln!(bytes, "{offset:010} 00000 n ").unwrap();
    }
    writeln!(
        bytes,
        "trailer\n<< /Size {} /Root 1 0 R >>\nstartxref\n{xref}\n%%EOF",
        objects.len()
    )
    .unwrap();
    fs::write(path, bytes).unwrap();
}
fn vector_pdf(ops: &[Op]) -> String {
    let mut s = String::from("q 0.71 0 0 -0.71 26 568 cm\n");
    for op in ops {
        match op {
            Op::Text(x, y, size, t) => writeln!(
                s,
                "0.08 0.14 0.2 rg BT /F1 {size} Tf 1 0 0 -1 {x} {y} Tm ({}) Tj ET",
                escape(t)
            )
            .unwrap(),
            Op::Rect(x, y, w, h, c) => {
                let (r, g, b) = rgb(c);
                writeln!(s, "{r} {g} {b} rg {x} {y} {w} {h} re f").unwrap();
            }
            Op::Line(x, y, a, b, c) => {
                let (r, g, bb) = rgb(c);
                writeln!(s, "{r} {g} {bb} RG 2 w {x} {y} m {a} {b} l S").unwrap();
            }
            Op::Dot(x, y, r, c) => {
                let (rr, g, b) = rgb(c);
                writeln!(
                    s,
                    "{rr} {g} {b} rg {} {} {} {} re f",
                    x - r,
                    y - r,
                    2. * r,
                    2. * r
                )
                .unwrap();
            }
        }
    }
    s.push_str("Q\n");
    s
}

//! Human-readable and graphical renderings of a [`super::TransferFile`]:
//! a Markdown summary table, and the ratio-versus-subgroup-size chart as
//! SVG (for the docs) and as TikZ (for the PDF report), both from the same
//! data so neither is hand-edited.

use super::{CurveTransfers, RouteClass, TransferFile};

fn best(c: &CurveTransfers) -> Option<&super::TransferRoute> {
    c.best_route.map(|i| &c.routes[i])
}

/// One Markdown table: every curve with a standard name plus every curve
/// whose verdict is a formula advantage, best carried route and its ratio.
pub fn markdown_table(file: &TransferFile) -> String {
    let mut rows: Vec<&CurveTransfers> = file
        .curves
        .iter()
        .filter(|c| {
            !c.standard_names.is_empty()
                || matches!(
                    c.verdict,
                    RouteClass::FormulaAdvantage
                        | RouteClass::FormulaMarginal
                        | RouteClass::Unresolved
                )
        })
        .collect();
    let key = |c: &CurveTransfers| best(c).map_or(f64::NEG_INFINITY, |r| r.cost.ratio_log2);
    rows.sort_by(|a, b| {
        key(b)
            .partial_cmp(&key(a))
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.slug.cmp(&b.slug))
    });
    let mut s = String::new();
    s.push_str("| curve | field | r bits | rho log2 | best carried route | genus | total log2 | exponent ratio | ratio log2 | class |\n");
    s.push_str("|:--|:--|--:|--:|:--|--:|--:|--:|--:|:--|\n");
    for c in rows {
        let name = c
            .standard_names
            .first()
            .map(|n| format!("{n} (`{}`)", c.slug))
            .unwrap_or_else(|| format!("`{}`", c.slug));
        let field = match c.field.kind.as_str() {
            "binary" => format!("GF(2^{})", c.field.degree),
            "extension" => format!(
                "GF(p^{}), p {} bits",
                c.field.degree,
                c.field.log2_base.round()
            ),
            _ => format!("GF(p), {} bits", c.field.log2_size.round()),
        };
        match best(c) {
            Some(r) => s.push_str(&format!(
                "| {name} | {field} | {:.0} | {:.1} | {} | {} | {:.1} | {:+.1} | {:+.1} | {} |\n",
                c.subgroup.bits,
                c.reference_log2,
                r.label,
                if r.destination.genus.len() > 12 {
                    format!("2^{:.0}", r.destination.log2_genus)
                } else {
                    r.destination.genus.clone()
                },
                r.cost.log2_total,
                r.cost.exponent_ratio_log2,
                r.cost.ratio_log2,
                super::class_name(c.verdict)
            )),
            None => s.push_str(&format!(
                "| {name} | {field} | {:.0} | {:.1} | no carried route | — | — | — | — | {} |\n",
                c.subgroup.bits,
                c.reference_log2,
                super::class_name(c.verdict)
            )),
        }
    }
    s
}

/// Chart data: one point per curve, `(subgroup bits, best ratio over carried
/// routes, family, label, class)`; the label is the standard name, else the slug.
pub fn points(file: &TransferFile) -> Vec<(f64, f64, String, String, RouteClass)> {
    file.curves
        .iter()
        .filter_map(|c| {
            best(c).map(|r| {
                (
                    c.subgroup.bits,
                    r.cost.ratio_log2,
                    c.family.clone(),
                    // Standard name when there is one: the chart labels it.
                    c.standard_names
                        .first()
                        .cloned()
                        .unwrap_or_else(|| c.slug.clone()),
                    c.verdict,
                )
            })
        })
        .collect()
}

const FAMILIES: [(&str, &str); 5] = [
    ("prime", "#1f77b4"),
    ("binary", "#d62728"),
    ("koblitz", "#2ca02c"),
    ("subfield", "#9467bd"),
    ("extension", "#ff7f0e"),
];

fn colour(family: &str) -> &'static str {
    FAMILIES
        .iter()
        .find(|(f, _)| *f == family)
        .map(|(_, c)| *c)
        .unwrap_or("#7f7f7f")
}

struct Axes {
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
}

fn axes(pts: &[(f64, f64, String, String, RouteClass)]) -> Axes {
    let x_max = pts.iter().map(|p| p.0).fold(0.0, f64::max).max(1.0);
    let y_min = pts.iter().map(|p| p.1).fold(0.0, f64::min);
    let y_max = pts.iter().map(|p| p.1).fold(0.0, f64::max);
    let pad = ((y_max - y_min) * 0.08).max(2.0);
    Axes {
        x_min: 0.0,
        x_max: (x_max * 1.05).ceil(),
        y_min: (y_min - pad).floor(),
        y_max: (y_max + pad).ceil(),
    }
}

/// The ratio chart as a self-contained SVG.
pub fn svg_chart(file: &TransferFile) -> String {
    let pts = points(file);
    let ax = axes(&pts);
    let (w, h, ml, mr, mt, mb) = (900.0, 520.0, 70.0, 30.0, 40.0, 60.0);
    let pw = w - ml - mr;
    let ph = h - mt - mb;
    let sx = |x: f64| ml + (x - ax.x_min) / (ax.x_max - ax.x_min) * pw;
    let sy = |y: f64| mt + (ax.y_max - y) / (ax.y_max - ax.y_min) * ph;
    let mut s = String::new();
    s.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{w}\" height=\"{h}\" viewBox=\"0 0 {w} {h}\" font-family=\"Helvetica, Arial, sans-serif\" font-size=\"12\">\n"
    ));
    s.push_str("<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n");
    s.push_str(&format!(
        "<text x=\"{}\" y=\"22\" text-anchor=\"middle\" font-size=\"15\">Jacobian transfer routes: formula ratio to Pollard rho, every registered curve (schema {})</text>\n",
        w / 2.0, file.schema
    ));
    // Axes and grid.
    s.push_str(&format!(
        "<line x1=\"{ml}\" y1=\"{mt}\" x2=\"{ml}\" y2=\"{}\" stroke=\"black\"/>\n<line x1=\"{ml}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"black\"/>\n",
        mt + ph, mt + ph, ml + pw, mt + ph
    ));
    let x_step = if ax.x_max > 400.0 { 100.0 } else { 50.0 };
    let mut x = 0.0;
    while x <= ax.x_max {
        s.push_str(&format!(
            "<line x1=\"{0}\" y1=\"{1}\" x2=\"{0}\" y2=\"{2}\" stroke=\"#ddd\"/><text x=\"{0}\" y=\"{3}\" text-anchor=\"middle\">{4}</text>\n",
            sx(x), mt, mt + ph, mt + ph + 16.0, x
        ));
        x += x_step;
    }
    let y_step = ((ax.y_max - ax.y_min) / 8.0).max(1.0).round();
    let mut y = (ax.y_min / y_step).ceil() * y_step;
    while y <= ax.y_max {
        s.push_str(&format!(
            "<line x1=\"{0}\" y1=\"{1}\" x2=\"{2}\" y2=\"{1}\" stroke=\"#ddd\"/><text x=\"{3}\" y=\"{4}\" text-anchor=\"end\">{5}</text>\n",
            ml, sy(y), ml + pw, ml - 6.0, sy(y) + 4.0, y
        ));
        y += y_step;
    }
    // Zero line and margin line.
    s.push_str(&format!(
        "<line x1=\"{ml}\" y1=\"{0}\" x2=\"{1}\" y2=\"{0}\" stroke=\"black\" stroke-dasharray=\"6 3\"/>\n<text x=\"{2}\" y=\"{3}\" text-anchor=\"end\" fill=\"#333\">rho boundary (ratio 0)</text>\n",
        sy(0.0), ml + pw, ml + pw - 4.0, sy(0.0) - 5.0
    ));
    let margin = file.cost_model.advantage_margin_bits;
    s.push_str(&format!(
        "<line x1=\"{ml}\" y1=\"{0}\" x2=\"{1}\" y2=\"{0}\" stroke=\"#555\" stroke-dasharray=\"2 3\"/>\n<text x=\"{2}\" y=\"{3}\" text-anchor=\"end\" fill=\"#555\">advantage margin (+{margin} bits)</text>\n",
        sy(margin), ml + pw, ml + pw - 4.0, sy(margin) - 5.0
    ));
    s.push_str(&format!(
        "<text x=\"{}\" y=\"{}\" text-anchor=\"middle\">subgroup order r, bits</text>\n",
        ml + pw / 2.0,
        h - 16.0
    ));
    s.push_str(&format!(
        "<text transform=\"translate(18,{}) rotate(-90)\" text-anchor=\"middle\">log2(rho on E) − log2(best carried route), bits</text>\n",
        mt + ph / 2.0
    ));
    for (xv, yv, family, slug, class) in &pts {
        let r = if matches!(
            class,
            RouteClass::FormulaAdvantage | RouteClass::FormulaMarginal | RouteClass::Unresolved
        ) {
            6.0
        } else {
            3.5
        };
        s.push_str(&format!(
            "<circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"{r}\" fill=\"{}\" fill-opacity=\"0.8\" stroke=\"{}\"><title>{slug}: ratio {yv:+.1} bits</title></circle>\n",
            sx(*xv),
            sy(*yv),
            colour(family),
            match class { RouteClass::FormulaAdvantage => "black", RouteClass::FormulaMarginal => "#555", RouteClass::Unresolved => "#e6a700", _ => "none" }
        ));
    }
    // Label the advantage points.
    let mut dy = 0.0;
    for (xv, yv, _, slug, class) in &pts {
        if matches!(
            class,
            RouteClass::FormulaAdvantage | RouteClass::FormulaMarginal | RouteClass::Unresolved
        ) {
            s.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{:.1}\" font-size=\"11\">{slug}</text>\n",
                sx(*xv) + 8.0,
                sy(*yv) + 4.0 + dy
            ));
            dy += 13.0;
        }
    }
    // Legend.
    let mut ly = mt + 10.0;
    for (f, c) in FAMILIES {
        s.push_str(&format!(
            "<circle cx=\"{}\" cy=\"{ly}\" r=\"4\" fill=\"{c}\"/><text x=\"{}\" y=\"{}\">{f}</text>\n",
            ml + 14.0, ml + 24.0, ly + 4.0
        ));
        ly += 16.0;
    }
    s.push_str(&format!(
        "<text x=\"{}\" y=\"{}\" font-size=\"11\" fill=\"#333\">{} curves; large rimmed points: formula_advantage (black), formula_marginal (grey), unresolved (amber; plotted at the favourable bound); every value is a formula under the frozen cost model, not a measurement</text>\n",
        ml, h - 36.0, pts.len()
    ));
    s.push_str("</svg>\n");
    s
}

/// The same chart as a TikZ picture (pgfplots is not assumed).
pub fn tikz_chart(file: &TransferFile) -> String {
    let pts = points(file);
    let ax = axes(&pts);
    let (pw, ph) = (14.0, 8.0);
    let sx = |x: f64| (x - ax.x_min) / (ax.x_max - ax.x_min) * pw;
    let sy = |y: f64| (y - ax.y_min) / (ax.y_max - ax.y_min) * ph;
    let mut s = String::new();
    s.push_str("\\begin{tikzpicture}[x=1cm,y=1cm]\n");
    s.push_str(&format!(
        "\\draw[->] (0,0) -- ({:.2},0) node[below right] {{$r$ bits}};\n",
        pw + 0.3
    ));
    s.push_str(&format!(
        "\\draw[->] (0,0) -- (0,{:.2}) node[above] {{ratio, bits}};\n",
        ph + 0.3
    ));
    let x_step = if ax.x_max > 400.0 { 100.0 } else { 50.0 };
    let mut x = 0.0;
    while x <= ax.x_max {
        s.push_str(&format!(
            "\\draw ({:.2},0) -- ++(0,-0.1) node[below] {{\\tiny {x}}};\n",
            sx(x)
        ));
        x += x_step;
    }
    let y_step = ((ax.y_max - ax.y_min) / 8.0).max(1.0).round();
    let mut y = (ax.y_min / y_step).ceil() * y_step;
    while y <= ax.y_max {
        s.push_str(&format!(
            "\\draw (0,{:.2}) -- ++(-0.1,0) node[left] {{\\tiny {y}}};\n",
            sy(y)
        ));
        y += y_step;
    }
    s.push_str(&format!(
        "\\draw[dashed] (0,{0:.2}) -- ({1:.2},{0:.2}) node[right] {{\\tiny rho}};\n",
        sy(0.0),
        pw
    ));
    s.push_str(&format!(
        "\\draw[dotted] (0,{0:.2}) -- ({1:.2},{0:.2}) node[right] {{\\tiny +{2} margin}};\n",
        sy(file.cost_model.advantage_margin_bits),
        pw,
        file.cost_model.advantage_margin_bits
    ));
    for (xv, yv, family, slug, class) in &pts {
        let (r, draw) = match class {
            RouteClass::FormulaAdvantage => ("4pt", "draw=black,"),
            RouteClass::FormulaMarginal => ("4pt", "draw=gray,"),
            RouteClass::Unresolved => ("4pt", "draw=orange,"),
            _ => ("1.8pt", ""),
        };
        let rgb = colour(family).trim_start_matches('#');
        let (cr, cg, cb) = (
            u8::from_str_radix(&rgb[0..2], 16).unwrap_or(0),
            u8::from_str_radix(&rgb[2..4], 16).unwrap_or(0),
            u8::from_str_radix(&rgb[4..6], 16).unwrap_or(0),
        );
        s.push_str(&format!(
            "\\fill[{draw}fill={{rgb,255:red,{cr};green,{cg};blue,{cb}}},fill opacity=0.8] ({:.2},{:.2}) circle ({r});\n",
            sx(*xv),
            sy(*yv)
        ));
        if matches!(
            class,
            RouteClass::FormulaAdvantage | RouteClass::FormulaMarginal | RouteClass::Unresolved
        ) {
            s.push_str(&format!(
                "\\node[right,font=\\tiny] at ({:.2},{:.2}) {{\\texttt{{{}}}}};\n",
                sx(*xv) + 0.15,
                sy(*yv) + 0.12,
                slug.replace('_', "\\_")
            ));
        }
    }
    let mut ly = ph - 0.3;
    for (f, c) in FAMILIES {
        let rgb = c.trim_start_matches('#');
        s.push_str(&format!(
            "\\fill[fill={{rgb,255:red,{};green,{};blue,{}}}] (0.4,{ly:.2}) circle (2pt); \\node[right,font=\\tiny] at (0.5,{ly:.2}) {{{f}}};\n",
            u8::from_str_radix(&rgb[0..2], 16).unwrap_or(0),
            u8::from_str_radix(&rgb[2..4], 16).unwrap_or(0),
            u8::from_str_radix(&rgb[4..6], 16).unwrap_or(0)
        ));
        ly -= 0.35;
    }
    s.push_str("\\end{tikzpicture}\n");
    s
}

//! Editable source for peak_ratios_all_patterns.svg and .csv.
//! Run: cargo run --locked --release --example exp7_render_all_patterns -- research/prime_fourier_flatness_exp7_20261007

use serde_json::Value;
use std::fmt::Write as _;
use std::fs;
use std::path::Path;

pub fn run() -> Result<(), Box<dyn std::error::Error>> {
    let directory = std::env::args()
        .nth(1)
        .ok_or("usage: exp7_render_all_patterns STUDY_DIRECTORY")?;
    let directory = Path::new(&directory);
    let results: Value =
        serde_json::from_slice(&fs::read(directory.join("results_all_patterns.json"))?)?;
    let curves = results["curves"].as_array().ok_or("missing curves")?;
    let cases = [
        ("smooth", "B-smooth x"),
        ("cf16", "CF quotients ≤16"),
        ("cantor3", "Cantor digits"),
        ("hamming", "Hamming weight"),
        ("farey", "Farey height"),
        ("legendre+++", "Legendre +++"),
        ("legendre++-", "Legendre ++−"),
        ("legendre+-+", "Legendre +−+"),
        ("legendre+--", "Legendre +−−"),
        ("legendre-++", "Legendre −++"),
        ("legendre-+-", "Legendre −+−"),
        ("legendre--+", "Legendre −−+"),
        ("legendre---", "Legendre −−−"),
        ("sha-x", "SHA(x) control"),
        ("random-pairs", "Random-pair control"),
        ("log-interval", "Log-interval control"),
    ];
    let mut csv = String::from("case,kind,p,n,m,observed_peak_z,null_q95_z,ratio,empirical_p\n");
    let mut svg = String::from(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1080\" height=\"1160\" viewBox=\"0 0 1080 1160\">\n",
    );
    svg.push_str("<rect width=\"1080\" height=\"1160\" fill=\"#fffdf9\"/>\n");
    svg.push_str("<text x=\"36\" y=\"48\" font-family=\"sans-serif\" font-size=\"24\" font-weight=\"bold\">EXP7 all-pattern Fourier peaks against matched nulls</text>\n");
    svg.push_str("<text x=\"36\" y=\"77\" font-family=\"sans-serif\" font-size=\"16\">Cell = observed maximum / 95th percentile of 2048 same-size null maxima</text>\n");
    for (c, curve) in curves.iter().enumerate() {
        let p = curve["p"].as_u64().ok_or("missing p")?;
        let n = curve["n"].as_u64().ok_or("missing n")?;
        let x = 311 + c * 185;
        writeln!(svg, "<text x=\"{x}\" y=\"123\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"17\" font-weight=\"bold\">p={p}</text>")?;
        writeln!(svg, "<text x=\"{x}\" y=\"144\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"14\">n={n}</text>")?;
    }
    for (r, &(name, label)) in cases.iter().enumerate() {
        let y = 173 + r * 55;
        let kind = if r < 13 { "candidate" } else { "control" };
        writeln!(
            svg,
            "<text x=\"40\" y=\"{}\" font-family=\"sans-serif\" font-size=\"16\">{label}</text>",
            y + 31
        )?;
        for (c, curve) in curves.iter().enumerate() {
            let row = curve["cases"]
                .as_array()
                .ok_or("missing cases")?
                .iter()
                .find(|row| row["name"] == name)
                .ok_or("missing case")?;
            let p = curve["p"].as_u64().ok_or("missing p")?;
            let n = curve["n"].as_u64().ok_or("missing n")?;
            let m = row["members"].as_u64().ok_or("missing m")?;
            let observed = row["peak_z"].as_f64().ok_or("missing peak")?;
            let q95 = row["null_q95_z"].as_f64().ok_or("missing null q95")?;
            let empirical_p = row["empirical_p"].as_f64().ok_or("missing p-value")?;
            let ratio = observed / q95;
            let x = 228 + c * 185;
            let color = if ratio > 2.0 {
                "#ee9a82"
            } else if ratio > 1.0 {
                "#f6dfa6"
            } else {
                "#c8e4eb"
            };
            writeln!(svg, "<rect x=\"{x}\" y=\"{y}\" width=\"166\" height=\"48\" rx=\"6\" fill=\"{color}\" stroke=\"#536b75\"/>")?;
            writeln!(svg, "<text x=\"{}\" y=\"{}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"19\" font-weight=\"bold\">{ratio:.2}×</text>", x + 83, y + 33)?;
            writeln!(
                csv,
                "{name},{kind},{p},{n},{m},{observed:.9},{q95:.9},{ratio:.9},{empirical_p:.9}"
            )?;
        }
    }
    svg.push_str("<rect x=\"40\" y=\"1075\" width=\"20\" height=\"20\" fill=\"#c8e4eb\" stroke=\"#536b75\"/><text x=\"68\" y=\"1091\" font-family=\"sans-serif\" font-size=\"14\">≤ null 95th percentile</text>\n");
    svg.push_str("<rect x=\"320\" y=\"1075\" width=\"20\" height=\"20\" fill=\"#f6dfa6\" stroke=\"#536b75\"/><text x=\"348\" y=\"1091\" font-family=\"sans-serif\" font-size=\"14\">&gt; null 95th percentile</text>\n");
    svg.push_str("<rect x=\"627\" y=\"1075\" width=\"20\" height=\"20\" fill=\"#ee9a82\" stroke=\"#536b75\"/><text x=\"655\" y=\"1091\" font-family=\"sans-serif\" font-size=\"14\">&gt; 2× null 95th percentile</text>\n");
    svg.push_str("<text x=\"40\" y=\"1135\" font-family=\"sans-serif\" font-size=\"13\">52 candidate cells; threshold p ≤ 0.05/52. A colored cell alone is not a lead.</text>\n</svg>\n");
    fs::write(directory.join("peak_ratios_all_patterns.svg"), svg)?;
    fs::write(directory.join("peak_ratios_all_patterns.csv"), csv)?;
    Ok(())
}

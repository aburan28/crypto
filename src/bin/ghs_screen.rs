//! CLI for the native classical-binary GHS structural screener.

use clap::Parser;
use crypto_lib::cryptanalysis::ghs_screen::{screen_curve, GhsCurveInput, GhsScreenReport};
use num_bigint::BigUint;
use serde::Serialize;
use std::{fs, path::PathBuf};

#[derive(Parser)]
#[command(about = "Screen every GHS field factorisation and emit exact cover towers")]
struct Args {
    /// Absolute degree N of the top field F_(2^N).
    #[arg(long)]
    degree: u32,
    /// Binary field polynomial bitset, including z^N (decimal or 0x hex).
    #[arg(long)]
    modulus: String,
    /// Curve coefficient a (decimal or 0x hex).
    #[arg(long, default_value = "0")]
    a: String,
    /// Nonzero curve coefficient b (decimal or 0x hex).
    #[arg(long)]
    b: String,
    /// Structural candidate bound for the descended genus.
    #[arg(long, default_value = "4")]
    genus_bound: String,
    /// Write JSON to this path instead of stdout.
    #[arg(long)]
    output: Option<PathBuf>,
}

#[derive(Serialize)]
struct JsonCurve {
    model: &'static str,
    absolute_degree: u32,
    modulus: String,
    a: String,
    b: String,
}

#[derive(Serialize)]
struct JsonSummary {
    factorisations: usize,
    candidates_within_bound: usize,
    genus_bound: String,
    best_relative_degree: Option<u32>,
    best_base_degree: Option<u32>,
    best_magic_number: Option<u32>,
    best_genus: Option<String>,
    interpretation: &'static str,
}

#[derive(Serialize)]
struct JsonTower {
    relative_degree: u32,
    base_degree: u32,
    field_tower: String,
    cover_tower: String,
    magic_number: u32,
    genus: String,
    ghs_type: &'static str,
    artin_schreier_cover_degree: String,
    within_genus_bound: bool,
}

#[derive(Serialize)]
struct JsonReport {
    schema: &'static str,
    curve: JsonCurve,
    summary: JsonSummary,
    towers: Vec<JsonTower>,
}

fn parse_number(value: &str, name: &str) -> Result<BigUint, String> {
    let cleaned: String = value
        .chars()
        .filter(|character| *character != '_')
        .collect();
    if cleaned.len() > 1300 {
        return Err(format!("{name} exceeds the 4096-bit input limit"));
    }
    let (digits, radix) = cleaned
        .strip_prefix("0x")
        .or_else(|| cleaned.strip_prefix("0X"))
        .map_or((cleaned.as_str(), 10), |digits| (digits, 16));
    if digits.is_empty() {
        return Err(format!("invalid {name}"));
    }
    BigUint::parse_bytes(digits.as_bytes(), radix).ok_or_else(|| format!("invalid {name}"))
}

fn json_report(input: &GhsCurveInput, report: &GhsScreenReport) -> JsonReport {
    let best = report.best();
    JsonReport {
        schema: "ghs-screen/v1",
        curve: JsonCurve {
            model: "y^2+xy=x^3+a*x^2+b",
            absolute_degree: input.absolute_degree,
            modulus: format!("0x{:x}", input.modulus),
            a: format!("0x{:x}", input.a),
            b: format!("0x{:x}", input.b),
        },
        summary: JsonSummary {
            factorisations: report.rows.len(),
            candidates_within_bound: report
                .rows
                .iter()
                .filter(|row| row.within_genus_bound)
                .count(),
            genus_bound: report.genus_bound.to_string(),
            best_relative_degree: best.map(|row| row.relative_degree),
            best_base_degree: best.map(|row| row.base_degree),
            best_magic_number: best.map(|row| row.magic_number),
            best_genus: best.map(|row| row.genus.to_string()),
            interpretation: "structural screen only; no DLP cost or vulnerability claim",
        },
        towers: report
            .rows
            .iter()
            .map(|row| JsonTower {
                relative_degree: row.relative_degree,
                base_degree: row.base_degree,
                field_tower: row.field_tower(report.absolute_degree),
                cover_tower: row.cover_tower(report.absolute_degree),
                magic_number: row.magic_number,
                genus: row.genus.to_string(),
                ghs_type: if row.type_i { "I" } else { "II" },
                artin_schreier_cover_degree: row.cover_degree.to_string(),
                within_genus_bound: row.within_genus_bound,
            })
            .collect(),
    }
}

fn run() -> Result<(), String> {
    let args = Args::parse();
    let input = GhsCurveInput {
        absolute_degree: args.degree,
        modulus: parse_number(&args.modulus, "modulus")?,
        a: parse_number(&args.a, "a")?,
        b: parse_number(&args.b, "b")?,
    };
    let genus_bound = parse_number(&args.genus_bound, "genus bound")?;
    let report = screen_curve(&input, &genus_bound).map_err(|error| error.to_string())?;
    let text = serde_json::to_string_pretty(&json_report(&input, &report))
        .map_err(|error| error.to_string())?
        + "\n";
    if let Some(output) = args.output {
        fs::write(output, text).map_err(|error| error.to_string())?;
    } else {
        print!("{text}");
    }
    Ok(())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("ghs_screen: {error}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::Zero;

    #[test]
    fn parses_decimal_hex_and_separators() {
        assert_eq!(parse_number("27", "n").unwrap(), BigUint::from(27u32));
        assert_eq!(
            parse_number("0x1_1b", "n").unwrap(),
            BigUint::from(0x11bu32)
        );
        assert!(parse_number("0x", "n").is_err());
        assert!(parse_number("-1", "n").is_err());
    }

    #[test]
    fn json_preserves_exact_large_genus_strings() {
        let input = GhsCurveInput {
            absolute_degree: 33,
            modulus: BigUint::zero(),
            a: BigUint::zero(),
            b: BigUint::from(1u32),
        };
        let genus = BigUint::from(1u32) << 32usize;
        let report = GhsScreenReport {
            absolute_degree: 33,
            genus_bound: BigUint::from(4u32),
            rows: vec![crypto_lib::cryptanalysis::ghs_screen::GhsCoverTower {
                relative_degree: 33,
                base_degree: 1,
                magic_number: 33,
                genus: genus.clone(),
                type_i: false,
                cover_degree: BigUint::from(1u32) << 33usize,
                within_genus_bound: false,
            }],
        };
        let value = serde_json::to_value(json_report(&input, &report)).unwrap();
        assert_eq!(value["summary"]["best_genus"], genus.to_string());
        assert_eq!(value["towers"][0]["genus"], "4294967296");
    }
}

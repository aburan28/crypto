//! Execute the checked trace branch of binary GHS descent.

use clap::Parser;
use crypto_lib::cryptanalysis::ghs_descent::Pt;
use crypto_lib::cryptanalysis::ghs_screen::GhsCurveInput;
use crypto_lib::cryptanalysis::ghs_transport::{trace_transport, BinaryPointInput, TraceStatus};
use num_bigint::BigUint;
use serde::Serialize;

#[derive(Parser)]
#[command(about = "Check binary composite-field trace transport on a supplied subgroup")]
struct Args {
    #[arg(long)]
    degree: u32,
    #[arg(long)]
    modulus: String,
    #[arg(long)]
    a: String,
    #[arg(long)]
    b: String,
    #[arg(long)]
    relative_degree: u32,
    #[arg(long)]
    p_x: String,
    #[arg(long)]
    p_y: String,
    #[arg(long)]
    q_x: String,
    #[arg(long)]
    q_y: String,
    /// Claimed subgroup annihilator. This command checks annihilation, not primality.
    #[arg(long)]
    order: String,
}

#[derive(Serialize)]
struct JsonPoint {
    infinity: bool,
    x: Option<String>,
    y: Option<String>,
}

#[derive(Serialize)]
struct JsonReport {
    schema: &'static str,
    characteristic: u8,
    absolute_degree: u32,
    relative_degree: u32,
    base_degree: u32,
    status: &'static str,
    generator_image: Option<JsonPoint>,
    target_image: Option<JsonPoint>,
    prime_order_certified: bool,
    target_relation_verified: bool,
}

fn parse(value: &str, name: &str) -> Result<BigUint, String> {
    let cleaned: String = value.chars().filter(|c| *c != '_').collect();
    if cleaned.len() > 1300 {
        return Err(format!("{name} exceeds the input limit"));
    }
    let (digits, radix) = cleaned
        .strip_prefix("0x")
        .or_else(|| cleaned.strip_prefix("0X"))
        .map_or((cleaned.as_str(), 10), |digits| (digits, 16));
    BigUint::parse_bytes(digits.as_bytes(), radix).ok_or_else(|| format!("invalid {name}"))
}

fn json_point(point: Pt) -> JsonPoint {
    match point {
        Pt::Inf => JsonPoint {
            infinity: true,
            x: None,
            y: None,
        },
        Pt::Aff { x, y } => JsonPoint {
            infinity: false,
            x: Some(format!("0x{:x}", x.to_biguint())),
            y: Some(format!("0x{:x}", y.to_biguint())),
        },
    }
}

fn run() -> Result<(), String> {
    let args = Args::parse();
    let curve = GhsCurveInput {
        absolute_degree: args.degree,
        modulus: parse(&args.modulus, "modulus")?,
        a: parse(&args.a, "a")?,
        b: parse(&args.b, "b")?,
    };
    let p = BinaryPointInput {
        x: parse(&args.p_x, "p-x")?,
        y: parse(&args.p_y, "p-y")?,
    };
    let q = BinaryPointInput {
        x: parse(&args.q_x, "q-x")?,
        y: parse(&args.q_y, "q-y")?,
    };
    let order = parse(&args.order, "order")?;
    let result = trace_transport(&curve, args.relative_degree, &p, &q, &order)?;
    let status = match result.status {
        TraceStatus::CurveNotDefinedOverSubfield => "curve_not_defined_over_subfield",
        TraceStatus::GeneratorKilledByTrace => "generator_killed_by_trace",
        TraceStatus::NonzeroSubgroupImage => "nonzero_subgroup_image",
    };
    let report = JsonReport {
        schema: "ghs-transport/v1",
        characteristic: 2,
        absolute_degree: args.degree,
        relative_degree: result.relative_degree,
        base_degree: result.base_degree,
        status,
        generator_image: result.generator_image.map(json_point),
        target_image: result.target_image.map(json_point),
        prime_order_certified: false,
        target_relation_verified: false,
    };
    println!(
        "{}",
        serde_json::to_string_pretty(&report).map_err(|error| error.to_string())?
    );
    Ok(())
}

fn main() {
    if let Err(error) = run() {
        eprintln!("ghs_transport: {error}");
        std::process::exit(1);
    }
}

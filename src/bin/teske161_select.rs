//! Reproduce the coefficient-selection step of Teske's F_(2^161) construction.
//!
//! Edlyn Teske, "An Elliptic Curve Trapdoor System", J. Cryptology 19
//! (2006), Section 2.1 and Theorem 1, DOI 10.1007/s00145-004-0328-3.
//! This selects structural candidates only. Group order, class group,
//! isogeny route, explicit GHS map and DLP cost are separate obligations.

use clap::{Parser, ValueEnum};
use crypto_lib::binary_ecc::{F2mElement, IrreduciblePoly};
use crypto_lib::cryptanalysis::ec_trapdoor::FieldTower;
use num_bigint::BigUint;
use num_traits::{One, Zero};
use rand::{rngs::SmallRng, RngCore, SeedableRng};
use serde::Serialize;

const N: u32 = 161;
const L: u32 = 23;

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Component {
    W1,
    W2,
}

impl Component {
    fn label(self) -> &'static str {
        match self {
            Self::W1 => "W1: ker(sigma^3 + sigma^2 + 1)",
            Self::W2 => "W2: ker(sigma^3 + sigma + 1)",
        }
    }
}

#[derive(Parser)]
#[command(about = "Select and certify a Teske-style genus-7/8 structural candidate")]
struct Args {
    /// Frobenius component of b.
    #[arg(long, value_enum, default_value = "w1")]
    component: Component,
    /// Desired genus (7: pure W1/W2; 8: add nonzero W0).
    #[arg(long, default_value_t = 7)]
    genus: u8,
    /// Reproducible research seed; this is not a secret generator.
    #[arg(long, default_value_t = 20261006)]
    seed: u64,
    /// Bounded rejection budget.
    #[arg(long, default_value_t = 100)]
    max_attempts: u32,
    /// Choose the twist representative a=0 or a=1.
    #[arg(long, default_value_t = 0)]
    a: u8,
}

#[derive(Serialize)]
struct Report {
    schema: &'static str,
    source: &'static str,
    field: &'static str,
    modulus: String,
    curve_model: &'static str,
    a: u8,
    b: String,
    component: &'static str,
    seed: u64,
    attempt: u32,
    magic_number: u32,
    relative_trace_zero: bool,
    genus: u8,
    verified_conditions: Vec<&'static str>,
    outstanding: Vec<&'static str>,
}

fn tower() -> FieldTower {
    FieldTower::new(
        N,
        7,
        L,
        IrreduciblePoly {
            degree: N,
            low_terms: vec![0, 18],
        },
    )
}

fn random_element(rng: &mut SmallRng) -> F2mElement {
    let mut bytes = [0u8; 21];
    rng.fill_bytes(&mut bytes);
    bytes[20] &= 1;
    F2mElement::from_biguint(&BigUint::from_bytes_le(&bytes), N)
}

/// For sigma of order seven, F_2[sigma] factors as
/// (sigma+1)(sigma^3+sigma^2+1)(sigma^3+sigma+1).
/// The complementary-factor products map onto W1 or W2.
fn project(t: &FieldTower, x: &F2mElement, component: Component) -> F2mElement {
    let mut powers = vec![x.clone()];
    for _ in 0..4 {
        let next = t.frobenius(powers.last().unwrap());
        powers.push(next);
    }
    let mut result = powers[0].add(&powers[2]).add(&powers[4]);
    result.add_assign(match component {
        Component::W1 => &powers[3],
        Component::W2 => &powers[1],
    });
    result
}

fn trace(t: &FieldTower, x: &F2mElement) -> F2mElement {
    let mut sum = F2mElement::zero(N);
    let mut cur = x.clone();
    for _ in 0..7 {
        sum.add_assign(&cur);
        cur = t.frobenius(&cur);
    }
    sum
}

/// Teske's equation (1) is the F_2 rank of (1, sqrt(sigma^i(b))).
/// The leading marker is a separate coordinate, not a field element.
fn magic_number(t: &FieldTower, b: &F2mElement) -> u32 {
    let mut cur = t.sqrt(b);
    let mut basis: Vec<BigUint> = Vec::new();
    for _ in 0..7 {
        let mut v = (cur.to_biguint() << 1usize) | BigUint::one();
        for pivot in &basis {
            if v.bit(pivot.bits() - 1) {
                v ^= pivot;
            }
        }
        if !v.is_zero() {
            basis.push(v);
        }
        cur = t.frobenius(&cur);
    }
    basis.len() as u32
}

fn in_component(t: &FieldTower, b: &F2mElement, component: Component) -> bool {
    let s1 = t.frobenius(b);
    let s2 = t.frobenius(&s1);
    let s3 = t.frobenius(&s2);
    let value = match component {
        Component::W1 => b.add(&s2).add(&s3),
        Component::W2 => b.add(&s1).add(&s3),
    };
    value.is_zero()
}

fn select(args: &Args) -> Result<Report, String> {
    if !matches!(args.genus, 7 | 8) {
        return Err("--genus must be 7 or 8".into());
    }
    if args.a > 1 {
        return Err("--a must be 0 or 1".into());
    }
    if args.max_attempts == 0 || args.max_attempts > 100_000 {
        return Err("--max-attempts must be in 1..=100000".into());
    }
    let t = tower();
    let mut rng = SmallRng::seed_from_u64(args.seed);
    for attempt in 1..=args.max_attempts {
        let w = project(&t, &random_element(&mut rng), args.component);
        if w.is_zero() || !in_component(&t, &w, args.component) {
            continue;
        }
        let w0 = if args.genus == 8 {
            trace(&t, &random_element(&mut rng))
        } else {
            F2mElement::zero(N)
        };
        if args.genus == 8 && w0.is_zero() {
            continue;
        }
        let b = w.add(&w0);
        let tr = trace(&t, &b);
        let m = magic_number(&t, &b);
        if m != 4 || tr.is_zero() != (args.genus == 7) {
            continue;
        }
        // In characteristic two, trace has seven identical W0 terms,
        // so the trace is exactly the selected W0 component.
        if tr != w0 {
            continue;
        }
        return Ok(Report {
            schema: "teske161-select/v1",
            source: "Teske 2006, Section 2.1 and Theorem 1; DOI 10.1007/s00145-004-0328-3",
            field: "F_(2^161) = F_((2^23)^7)",
            modulus: format!("0x{:x}", (BigUint::one() << 161usize) | (BigUint::one() << 18usize) | BigUint::one()),
            curve_model: "y^2 + x*y = x^3 + a*x^2 + b",
            a: args.a,
            b: format!("0x{:x}", b.to_biguint()),
            component: args.component.label(),
            seed: args.seed,
            attempt,
            magic_number: m,
            relative_trace_zero: tr.is_zero(),
            genus: args.genus,
            verified_conditions: vec![
                "b in W0 plus one nonzero W1/W2 component",
                "augmented Frobenius orbit has F2 rank four",
                "relative trace determines genus seven or eight",
            ],
            outstanding: vec![
                "group order, prime subgroup, cofactor and embedding degree",
                "squarefree Frobenius discriminant and class-group criteria",
                "certified binary-field isogeny path and signed point maps",
                "explicit smooth GHS curve, subgroup-preserving map and full-cost DLP comparison",
            ],
        });
    }
    Err("selection budget exhausted".into())
}

fn main() {
    let args = Args::parse();
    match select(&args) {
        Ok(report) => println!("{}", serde_json::to_string_pretty(&report).unwrap()),
        Err(error) => {
            eprintln!("teske161_select: {error}");
            std::process::exit(1);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn both_components_and_both_genus_branches() {
        for component in [Component::W1, Component::W2] {
            for genus in [7, 8] {
                let report = select(&Args {
                    component,
                    genus,
                    seed: 20261006,
                    max_attempts: 10,
                    a: 0,
                })
                .unwrap();
                assert_eq!(report.magic_number, 4);
                assert_eq!(report.genus, genus);
                assert_eq!(report.relative_trace_zero, genus == 7);
                let b = F2mElement::from_biguint(
                    &BigUint::parse_bytes(&report.b.as_bytes()[2..], 16).unwrap(),
                    N,
                );
                let t = tower();
                let w = b.add(&trace(&t, &b));
                assert!(in_component(&t, &w, component));
                assert!(!w.is_zero());
            }
        }
    }

    #[test]
    fn rejects_invalid_inputs() {
        let args = Args {
            component: Component::W1,
            genus: 16,
            seed: 0,
            max_attempts: 1,
            a: 0,
        };
        assert!(select(&args).is_err());
    }
}

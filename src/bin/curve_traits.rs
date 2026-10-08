//! `curve_traits`: size-independent traits of every registered curve, and
//! queries that group curves by them or rank curves by likeness to one.
//!
//! ```text
//! curve_traits build                       # docs/curves/registry.json → docs/curves/traits.json
//! curve_traits check                       # the committed file is current
//! curve_traits show ECC2K-130              # one curve's record
//! curve_traits group --by cm,descent       # curves sharing those keys, across sizes
//! curve_traits similar sect163k1 --other-sizes
//! curve_traits keys                        # what each key means
//! ```
//!
//! See `docs/curves/TRAITS.md`.

use std::process::ExitCode;

use clap::{Parser, Subcommand};
use crypto_lib::cryptanalysis::curve_traits::query::{self, DEFAULT_WEIGHTS};
use crypto_lib::cryptanalysis::curve_traits::{
    signature, CurveTraits, TraitFile, DEFAULT_BUDGET, DEFAULT_SIGNATURE, KEYS,
};

const REGISTRY: &str = "docs/curves/registry.json";
const TRAITS: &str = "docs/curves/traits.json";

#[derive(Parser)]
#[command(
    about = "Size-independent traits of registered curves: build, check, group, find similar"
)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Derive every curve's traits from the registry and write the file.
    Build {
        #[arg(long, default_value = REGISTRY)]
        registry: String,
        #[arg(long, default_value = TRAITS)]
        out: String,
        /// Pollard-rho iterations per composite.
        #[arg(long, default_value_t = DEFAULT_BUDGET)]
        budget: u64,
    },
    /// Rebuild with the committed file's budget and fail if it differs.
    Check {
        #[arg(long, default_value = REGISTRY)]
        registry: String,
        #[arg(long, default_value = TRAITS)]
        traits: String,
    },
    /// Print one curve's record (slug, standard name or legacy spelling).
    Show {
        name: String,
        #[arg(long, default_value = TRAITS)]
        traits: String,
        /// The raw JSON record.
        #[arg(long)]
        json: bool,
    },
    /// Group curves that share every key in --by.
    Group {
        /// Comma-separated keys (`curve_traits keys` lists them).
        #[arg(long, value_delimiter = ',', default_values_t = DEFAULT_SIGNATURE.iter().map(|s| s.to_string()))]
        by: Vec<String>,
        /// Hide groups with fewer members.
        #[arg(long, default_value_t = 1)]
        min: usize,
        /// Only groups spanning at least this many field sizes.
        #[arg(long, default_value_t = 1)]
        min_sizes: usize,
        #[arg(long, default_value = TRAITS)]
        traits: String,
        #[arg(long)]
        json: bool,
    },
    /// Rank every other curve by likeness to one.
    Similar {
        name: String,
        #[arg(long, default_value_t = 15)]
        top: usize,
        /// Leave out curves over a field of the same size.
        #[arg(long)]
        other_sizes: bool,
        /// Comma-separated keys to compare, each weighted 1, instead of
        /// the default weights.
        #[arg(long, value_delimiter = ',')]
        keys: Vec<String>,
        #[arg(long, default_value = TRAITS)]
        traits: String,
        #[arg(long)]
        json: bool,
    },
    /// List the grouping keys and what they mean.
    Keys,
}

fn load(path: &str) -> Result<TraitFile, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
    serde_json::from_str(&text).map_err(|e| format!("{path}: {e}"))
}

fn find<'a>(file: &'a TraitFile, name: &str) -> Result<&'a CurveTraits, String> {
    file.find(name)
        .ok_or_else(|| format!("{name}: not a registered curve name"))
}

fn field_label(c: &CurveTraits) -> String {
    match c.field.kind.as_str() {
        "binary" => format!("GF(2^{})", c.field.degree),
        _ => format!("GF(p:{}b)", c.field.size_bits),
    }
}

fn status(s: impl serde::Serialize) -> String {
    serde_json::to_value(s)
        .ok()
        .and_then(|v| v.as_str().map(str::to_string))
        .unwrap_or_default()
}

fn show(c: &CurveTraits) {
    let f = &c.frobenius;
    println!("{}  ({}, {})", c.slug, c.family, field_label(c));
    println!("  signature   {}", signature(c, DEFAULT_SIGNATURE));
    println!(
        "  order       {}  [{}: {}]",
        c.order,
        status(c.order_check.status),
        c.order_check.method
    );
    println!("  trace       {}  (t/2√q = {})", c.trace, c.trace_ratio);
    println!("  Δ = t²−4q   {}  [{}]", f.disc, status(f.disc_status));
    println!(
        "  d_K         {}  h(d_K) = {}",
        f.cm_disc.as_deref().unwrap_or("unknown"),
        f.class_number.map_or("—".to_string(), |h| h.to_string())
    );
    let primes: Vec<String> = f
        .conductor_primes
        .iter()
        .map(|(p, e)| {
            if *e == 1 {
                p.clone()
            } else {
                format!("{p}^{e}")
            }
        })
        .chain(
            f.conductor_unfactored
                .iter()
                .map(|(c, e)| format!("[{c}]^{e}")),
        )
        .collect();
    println!(
        "  conductor   v {} {} ({} bits, fraction {}) = {}  [{}]",
        if f.cm_disc.is_some() { "=" } else { "≥" },
        f.conductor,
        f.conductor_bits,
        f.conductor_fraction
            .map_or("—".to_string(), |x| x.to_string()),
        if primes.is_empty() {
            "1".into()
        } else {
            primes.join(" · ")
        },
        status(f.status)
    );
    let s = &c.subfield;
    if c.field.kind == "binary" {
        println!(
            "  subfields   j ∈ GF(2^{}), defined over GF(2^{}), t_k = {}{}  [{}]",
            s.j_field_degree.unwrap_or(0),
            s.definition_degree.unwrap_or(0),
            if s.base_trace_sign_free { "±" } else { "" },
            s.base_trace.as_deref().unwrap_or("—"),
            status(s.status)
        );
    }
    println!(
        "  subgroup    r: {} bits, cofactor {}  [{}]",
        c.subgroup.bits.map_or("?".into(), |b| b.to_string()),
        c.subgroup.cofactor.as_deref().unwrap_or("?"),
        status(c.subgroup.status)
    );
    let e = &c.embedding;
    println!(
        "  embedding   k = {}  [{}]",
        e.degree
            .clone()
            .or(e.lower_bound.map(|b| format!("> {b}")))
            .unwrap_or("?".into()),
        status(e.status)
    );
    println!(
        "  twist       cofactor {}, largest prime {} bits  [{}]",
        c.twist.cofactor.as_deref().unwrap_or("?"),
        c.twist
            .largest_prime_bits
            .map_or("?".into(), |b| b.to_string()),
        status(c.twist.status)
    );
    let ells: Vec<String> = c
        .small_primes
        .iter()
        .map(|s| format!("{}{}{}", s.ell, &s.splitting[..1], s.depth))
        .collect();
    println!(
        "  ℓ ≤ 31      {}  (s split, i inert, r ramified; then the ℓ-volcano depth)",
        ells.join(" ")
    );
}

fn run(cli: Cli) -> Result<bool, String> {
    match cli.cmd {
        Cmd::Build {
            registry,
            out,
            budget,
        } => {
            let text =
                std::fs::read_to_string(&registry).map_err(|e| format!("{registry}: {e}"))?;
            let file = TraitFile::build(&text, budget)?;
            std::fs::write(&out, file.to_json()).map_err(|e| format!("{out}: {e}"))?;
            eprintln!("wrote {} curves to {out}", file.curves.len());
            Ok(true)
        }
        Cmd::Check { registry, traits } => {
            let committed =
                std::fs::read_to_string(&traits).map_err(|e| format!("{traits}: {e}"))?;
            let budget = load(&traits)?.factor_budget;
            let text =
                std::fs::read_to_string(&registry).map_err(|e| format!("{registry}: {e}"))?;
            let fresh = TraitFile::build(&text, budget)?.to_json();
            if fresh == committed {
                eprintln!("{traits} is current");
                return Ok(true);
            }
            eprintln!(
                "{traits} is stale: run `cargo run --release --bin curve_traits -- build --budget {budget}`"
            );
            Ok(false)
        }
        Cmd::Show { name, traits, json } => {
            let file = load(&traits)?;
            let c = find(&file, &name)?;
            if json {
                println!("{}", serde_json::to_string_pretty(c).expect("serialisable"));
            } else {
                show(c);
            }
            Ok(true)
        }
        Cmd::Group {
            by,
            min,
            min_sizes,
            traits,
            json,
        } => {
            let file = load(&traits)?;
            check_keys(&by)?;
            let keys: Vec<&str> = by.iter().map(String::as_str).collect();
            let groups: Vec<_> = query::group(&file.curves, &keys)
                .into_iter()
                .filter(|g| g.members.len() >= min && g.sizes().len() >= min_sizes)
                .collect();
            if json {
                let v: Vec<_> = groups
                    .iter()
                    .map(|g| {
                        serde_json::json!({
                            "key": keys.iter().zip(&g.values).map(|(k, v)| (k.to_string(), v.clone())).collect::<std::collections::BTreeMap<_, _>>(),
                            "sizes": g.sizes(),
                            "members": g.members.iter().map(|c| c.slug.clone()).collect::<Vec<_>>(),
                        })
                    })
                    .collect();
                println!(
                    "{}",
                    serde_json::to_string_pretty(&v).expect("serialisable")
                );
                return Ok(true);
            }
            for g in &groups {
                let label: Vec<String> = keys
                    .iter()
                    .zip(&g.values)
                    .map(|(k, v)| format!("{k}={v}"))
                    .collect();
                let sizes: Vec<String> = g.sizes().iter().map(u64::to_string).collect();
                println!(
                    "{}  — {} curves, sizes {}",
                    label.join(" "),
                    g.members.len(),
                    sizes.join(",")
                );
                for c in &g.members {
                    println!("    {:<14} {:<10} {}", field_label(c), c.family, c.slug);
                }
            }
            eprintln!("{} groups", groups.len());
            Ok(true)
        }
        Cmd::Similar {
            name,
            top,
            other_sizes,
            keys,
            traits,
            json,
        } => {
            let file = load(&traits)?;
            let target = find(&file, &name)?;
            check_keys(&keys)?;
            let custom: Vec<(&str, u32)> = keys.iter().map(|k| (k.as_str(), 1)).collect();
            let weights = if custom.is_empty() {
                DEFAULT_WEIGHTS
            } else {
                &custom[..]
            };
            let ranked = query::similar(&file.curves, target, weights, other_sizes);
            if json {
                let v: Vec<_> = ranked
                    .iter()
                    .take(top)
                    .map(|m| {
                        serde_json::json!({
                            "slug": m.curve.slug,
                            "field": field_label(m.curve),
                            "mismatch": m.mismatch,
                            "distance": m.distance,
                            "differs": m.differs,
                        })
                    })
                    .collect();
                println!(
                    "{}",
                    serde_json::to_string_pretty(&v).expect("serialisable")
                );
                return Ok(true);
            }
            println!(
                "{}  ({}, {})",
                target.slug,
                target.family,
                field_label(target)
            );
            println!("  {}", signature(target, DEFAULT_SIGNATURE));
            println!("  rank  miss  dist     field          family     slug  [differs in]");
            for (i, m) in ranked.iter().take(top).enumerate() {
                println!(
                    "  {:>4}  {:>4}  {:<7.4}  {:<14} {:<10} {}  [{}]",
                    i + 1,
                    m.mismatch,
                    m.distance,
                    field_label(m.curve),
                    m.curve.family,
                    m.curve.slug,
                    m.differs.join(" ")
                );
            }
            Ok(true)
        }
        Cmd::Keys => {
            for (k, d) in KEYS {
                println!("{k:<15} {d}");
            }
            println!("\ndefault signature: {}", DEFAULT_SIGNATURE.join(","));
            Ok(true)
        }
    }
}

fn check_keys(keys: &[String]) -> Result<(), String> {
    for k in keys {
        if !KEYS.iter().any(|(name, _)| name == k) {
            return Err(format!("unknown key {k}; `curve_traits keys` lists them"));
        }
    }
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(true) => ExitCode::SUCCESS,
        Ok(false) => ExitCode::FAILURE,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

//! Producer-only implementation.
//!
//! This module is deliberately private to `p192_weighted_cm`.  The independent
//! verifier is a separate Cargo package with its own arithmetic.

mod arithmetic;
mod artifacts;
mod certificate;
mod controls;
mod digest;
mod evidence;
mod factor_base;
mod pocklington;
mod primality;
mod provenance;
mod ref0;
mod schema;

use std::ffi::OsString;
use std::path::PathBuf;

pub type Result<T> = std::result::Result<T, String>;

const USAGE: &str = "p192_weighted_cm preflight --curve p192 --sieve-bound 65521 --map-bound 113 --large-prime-bound 2147483647 --positive-control D=-23 --reference-shell REF-0 --reference-v-max 2 --reference-x-max-inclusive 1024 --out RUN_DIR\n";

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct PreflightArgs {
    pub(super) curve: String,
    pub(super) sieve_bound: u64,
    pub(super) map_bound: u64,
    pub(super) large_prime_bound: u64,
    pub(super) positive_control: String,
    pub(super) reference_shell: String,
    pub(super) reference_v_max: u64,
    pub(super) reference_x_max_inclusive: u64,
    pub(super) out: PathBuf,
}

fn parse_u64(value: &str, option: &str) -> Result<u64> {
    value
        .parse::<u64>()
        .map_err(|_| format!("invalid {option}: {value}"))
}

fn parse_preflight(args: &[String]) -> Result<PreflightArgs> {
    if args.len() != 19 || args.first().map(String::as_str) != Some("preflight") {
        return Err(format!(
            "only the preflight command is implemented\n{USAGE}"
        ));
    }
    let expected_options = [
        "--curve",
        "--sieve-bound",
        "--map-bound",
        "--large-prime-bound",
        "--positive-control",
        "--reference-shell",
        "--reference-v-max",
        "--reference-x-max-inclusive",
        "--out",
    ];
    for (offset, expected) in expected_options.iter().enumerate() {
        let actual = &args[1 + offset * 2];
        if actual != expected {
            return Err(format!(
                "frozen argv requires {expected} at position {}, found {actual}",
                2 + offset * 2
            ));
        }
    }
    let parsed = PreflightArgs {
        curve: args[2].clone(),
        sieve_bound: parse_u64(&args[4], "--sieve-bound")?,
        map_bound: parse_u64(&args[6], "--map-bound")?,
        large_prime_bound: parse_u64(&args[8], "--large-prime-bound")?,
        positive_control: args[10].clone(),
        reference_shell: args[12].clone(),
        reference_v_max: parse_u64(&args[14], "--reference-v-max")?,
        reference_x_max_inclusive: parse_u64(&args[16], "--reference-x-max-inclusive")?,
        out: PathBuf::from(&args[18]),
    };
    if parsed.curve != "p192"
        || parsed.sieve_bound != 65_521
        || parsed.map_bound != 113
        || parsed.large_prime_bound != 2_147_483_647
        || parsed.positive_control != "D=-23"
        || parsed.reference_shell != ref0::REFERENCE_SHELL
        || parsed.reference_v_max != ref0::REFERENCE_V_MAX
        || parsed.reference_x_max_inclusive != ref0::REFERENCE_X_MAX_INCLUSIVE
        || parsed.out.as_os_str().is_empty()
    {
        return Err("preflight arguments differ from the frozen protocol".to_owned());
    }
    Ok(parsed)
}

pub fn run<I, S>(args: I) -> Result<String>
where
    I: IntoIterator<Item = S>,
    S: Into<OsString>,
{
    let args = args
        .into_iter()
        .map(|value| {
            value
                .into()
                .into_string()
                .map_err(|_| "arguments must be valid UTF-8".to_owned())
        })
        .collect::<Result<Vec<_>>>()?;
    let parsed = parse_preflight(&args)?;
    let provenance = provenance::clean()?;
    artifacts::write_preflight(&parsed, &provenance)?;
    Ok(String::new())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn frozen_argv() -> Vec<&'static str> {
        vec![
            "preflight",
            "--curve",
            "p192",
            "--sieve-bound",
            "65521",
            "--map-bound",
            "113",
            "--large-prime-bound",
            "2147483647",
            "--positive-control",
            "D=-23",
            "--reference-shell",
            "REF-0",
            "--reference-v-max",
            "2",
            "--reference-x-max-inclusive",
            "1024",
            "--out",
            "RUN_DIR",
        ]
    }

    #[test]
    fn only_the_exact_ref0_preflight_is_admitted() {
        assert!(parse_preflight(
            &frozen_argv()
                .into_iter()
                .map(str::to_owned)
                .collect::<Vec<_>>()
        )
        .is_ok());
        for forbidden in [
            vec!["sieve"],
            vec!["preflight", "BOX-0"],
            vec!["preflight", "BOX-1"],
            vec!["preflight", "BOX-2"],
        ] {
            let owned = forbidden.into_iter().map(str::to_owned).collect::<Vec<_>>();
            assert!(parse_preflight(&owned).is_err(), "admitted {owned:?}");
        }
    }

    #[test]
    fn emitted_argv_contains_no_search_dispatch_token() {
        let argv = frozen_argv();
        assert!(!argv
            .iter()
            .any(|value| matches!(*value, "sieve" | "BOX-0" | "BOX-1" | "BOX-2")));
    }
}

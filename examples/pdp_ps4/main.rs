//! PS-4 of the PDP oracle speedup program
//! (`research/pdp_oracle_speedup_20261011/protocols/PS-4.md`): candidate
//! reduction by Frobenius-orbit representatives and by the 2-torsion
//! symmetry, on a Hamming-weight base in a type-II optimal normal basis,
//! with a scalar group-law pair kernel.
//!
//! ```bash
//! cargo build --release --example pdp_ps4
//! pdp_ps4 check-field   --out raw/field_checks.jsonl
//! pdp_ps4 check-curve   --out raw/curve_checks.jsonl
//! pdp_ps4 check-theorem --out raw/theorem_checks.jsonl
//! pdp_ps4 run     --n 11 --w 3 --out raw/runs.jsonl [--panel raw/panel_n11_w3.jsonl]
//! pdp_ps4 torsion --n 11 --w 3 --out raw/torsion.jsonl
//! pdp_ps4 control --l 11 --corpus ../baseline_v1/corpus/corpus_l11.txt --out raw/control.jsonl   # --features pdp-opcount
//! pdp_ps4 summary --raw raw
//! ```
//!
//! This is a stage measurement (`AGENTS.md` §8): it counts the work of the
//! pairs-and-solve decision of one abscissa on toy fields, never an end-to-end
//! ECDLP.  Wall times in the records are labelled `ms_unisolated` and are not
//! results.

mod base;
mod checks;
#[cfg(feature = "pdp-opcount")]
mod control;
mod curve;
mod field;
mod kernel;
mod model;
mod oracle;
mod run;
mod summary;
mod targets;
mod torsion;

pub fn arg<T: std::str::FromStr>(args: &[String], name: &str) -> Option<T> {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .and_then(|v| v.parse().ok())
}

pub fn arg_s(args: &[String], name: &str) -> Option<String> {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .cloned()
}

pub fn hex(x: u64) -> String {
    format!("0x{x:x}")
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let Some(cmd) = args.first().cloned() else {
        eprintln!("usage: pdp_ps4 (check-field|check-curve|check-theorem|run|torsion|control|summary) ...");
        std::process::exit(2);
    };
    let rest = &args[1..];
    match cmd.as_str() {
        "check-field" => checks::cmd_check_field(rest),
        "check-curve" => checks::cmd_check_curve(rest),
        "check-theorem" => checks::cmd_check_theorem(rest),
        "run" => run::cmd_run(rest),
        "torsion" => torsion::cmd_torsion(rest),
        "control" => model::cmd_control(rest),
        "summary" => summary::cmd_summary(rest),
        other => {
            eprintln!("unknown command {other}");
            std::process::exit(2);
        }
    }
}

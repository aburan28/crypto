//! Emit the worked identity certificate for Schoenhage's ten-multiplication
//! identity (Lemma 6 of `arXiv:2610.06783v1`) as JSON on stdout.
//!
//! ```sh
//! cargo run --example identity_certificate_schoenhage > docs/identity-certificates/schoenhage-10mult.json
//! ```
//!
//! The committed file is pinned by
//! `cryptanalysis::identity_certificate::tests::committed_schoenhage_certificate_still_checks`.
//! Regenerating it with the same seed reproduces it byte for byte; a change
//! to the statement changes the id, which is the point.

use crypto_lib::cryptanalysis::identity_certificate::{check, issue, schoenhage_identity, Verdict};

fn main() {
    let cert = issue(
        schoenhage_identity(),
        0x2610_0678_3000_0001,
        32,
        "issued by examples/identity_certificate_schoenhage.rs from the module constructor; \
         checked independently by cryptanalysis::identity_certificate::check replaying from the sealed id and seed",
        None,
    )
    .expect("Lemma 6 holds");
    match check(&cert).expect("well-formed") {
        Verdict::Accept { .. } => {}
        other => panic!("self-check failed: {other:?}"),
    }
    println!(
        "{}",
        serde_json::to_string_pretty(&cert).expect("serialises")
    );
}

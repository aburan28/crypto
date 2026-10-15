# Post-run Rust lint repair and source custody

CI on stacked PR #1920 reported `manual_is_multiple_of` warnings under
Rust 1.98 in the χ and pairing-batch examples. After the V3 χ and V7
pairing measurements, the changed lines were limited to integer
divisibility predicates and their `is_multiple_of` spelling. The
historical binaries, specs and raw records were not replaced or rerun.

| Example | Historical source at `cb5b9eb7a` | Post-lint source | Scope |
| --- | --- | --- | --- |
| `examples/chi_polynomial_decomposition.rs` | SHA-256 `98e2c2588f70f4b6f50f9d28523c403c3692a39440878e14a8e6e264bea08993` | SHA-256 `51a57fca0cb50bda5fd9279877c98688e74e68da9b7b04274e0e7d2c7bb4cd1b` | Three divisibility predicates; rustfmt wrapped one assertion. |
| `examples/torsion_pairing_batch.rs` | SHA-256 `55e23e6c7c361ce295ea7c3313f82ad4734ef900eb7a71bdafc9ac78dbbcea6e` | SHA-256 `5d63fb62d38c6be05275e99aafe59976457a83c68d74f0034a8c79d3bbdd44f2` | Three divisibility predicates. |

The historical sources remain obtainable with `git show
cb5b9eb7a:examples/<name>.rs`; the frozen `RUN_SPEC_V3.json` and
`RUN_SPEC_V7.json` identify those exact source hashes. The post-lint
sources receive correctness/lint validation only. Earlier timing and
operation receipts refer exclusively to their hash-pinned historical
binaries, so they must not be relabeled as measurements of the
post-lint build. Any future performance comparison needs a newly
frozen source and binary spec for every arm.

The same CI run also reported an inherited curve registry mismatch and
an older isogeny example's `RhoClasses::canon` signature mismatch. The
registry files, their checker and builder, the older example, and the
trait definition have no diff between PR #1920 and its PR #1749 base.
Those failures are tracked as upstream CI gates rather than attributed
to the χ and pairing measurements. The cancelled integration job also
remains non-green; its next head run must complete before merging.

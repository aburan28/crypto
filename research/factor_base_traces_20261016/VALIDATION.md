# Validation receipt — historical factor-base trace import

This is artifact validation, not a performance benchmark or a research run.
Frozen #1898 inputs were only read. No draft experiment was dispatched.

## Native trace tool

- Rust 1.98.0; `cargo test --release --bin factor_base_trace`: **4 passed**.
- `cargo clippy --release --bin factor_base_trace -- -D warnings`: **passed**.
  The first lint run rejected a manual range test; the source was corrected
  and both tests and lint were rerun. The failed lint log is retained.
- `cargo build --release --bin factor_base_trace`: built the production tool.
- `factor_base_trace verify FROZEN_GENERIC_JSON EXPORTED_DIRECTORY`:
  **all 384 source-backed traces and caches verified**, including signed
  lookup of every archived representative.
- The catalog records **52,876 representative rows across cases** and full
  SHA-256 commitments. It does not claim global uniqueness across bases.
- Changed Rust formatting and diff-whitespace checks pass.

Raw trace-test and lint receipts are in [validation/](validation/).

## Common library validation

The root release-library suite was run during the companion #1963 CI fix:
4,271 tests, initially **4,146 passed, 8 failed, 117 ignored**. Six failures
were missing committed runtime fixtures in the sparse checkout:

- `research/notes/ecc2k130/n37_four_policy_support_20261003/SOURCE42.jsonl`
- `research/cryptopro_b_chain_sweep_20261006/frozen_inputs/chains.constants.json`

Those files were materialized unchanged. Two Jacobian-cover tests hit their
solver budgets during the parallel run. All eight failed tests were rerun
serially using the same release test binary: **8 passed**, with unchanged
budgets, assertions and ignore settings. The initial full-suite command is
not represented as a single passing run. Together these runs cover
**4,154 passing library tests and 117 existing ignored tests**.

The [initial full log](validation/library-initial.log.gz) and
[eight-test rerun](validation/library-eight-retry.log) retain the failures
and their resolution. The common suite was built from the #1963 working tree
at `d50daae602a2314bdb6fe0d9aa14c4969fe1eb4f`. The trace PR's starting main
commit is `78723c49b8bc947a9cce9013eb03f3e51cee05ec`.
The GitHub comparison and the CI scope tool showed **no changed root library
sources, Cargo configuration or embedded inputs** between these bases.
This PR adds a separate binary and 48 data-only curve registrations and
alias entries. Every previous registry record is unchanged.

After those additions, the library suite was rebuilt and rerun serially:
**4,145 passed, 9 failed, 117 ignored**. All nine failures were local mock
TCP tests denied `bind` by the execution sandbox. The same nine tests passed
with networking enabled, using the same rebuilt release test binary and
unchanged assertions. The [rebuilt full run](validation/library-registered.log.gz)
and [nine-test rerun](validation/library-network-retry.log) preserve both
results: **4,154 unique library passes and 117 existing ignores** for the
updated registry, without representing the restricted full run as a single
passing command. The registry naming check reports **0 problems**, and the
native importer verifies all 48 registered models across the 384 cases.

## Storage and documents

- Frozen source SHA-256:
  `07ee0c22dcb48969581739bdc91e4ff6c78a4ad86d52ba6792d773fa3778e248`.
- Store name: `factor-base-traces/historical-20261015/registered-20261016.zip`.
- Bundle: **2,462,662 bytes**, SHA-256
  `923579b932103882f85b4475118a5edc7404fa4138823c1bdd59a633ca6f9c9c`.
- Native store push succeeded. Native `store verify --remote` reported
  **cache ok, remote ok; 1 entries verified**. The initial unauthenticated
  attempt failed with `NoCredentials`; its receipt is retained. Credentials
  are not included in source, artifacts, index metadata or these logs.
- The initial export remains at its original indexed name and digest.
  Every point, order and derivation field, and all 384 cache files were
  compared unchanged; only curve-registration status metadata changed.
  The final archive was written as a new immutable object, without replacing
  the existing archive. Its public-source Git blob matches the frozen
  repository archive, and destination ownership was verified against the
  authenticated AWS account.
- Markdown, editable DOT, SVG, PNG and PDF were rendered. Every page of the
  three-page PDF was visually inspected; text bounding boxes do not overrun
  the printable page width.
- Conductor is unavailable; the pinned `.conductor/WORKFLOW.md` returns 404.
  These are reported unavailable, not passed.

No native timing, relation yield, logarithm-recovery or end-to-end speed claim
is promoted from these validation diagnostics. Historical Sage/C evidence and
the separate unrun, zero-funded draft card remain unchanged.

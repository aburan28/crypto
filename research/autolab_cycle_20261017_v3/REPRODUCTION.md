# Reproduction and scope

All scientific computation, exact small oracles, field arithmetic and worker
harnesses are Rust. The user-required existing Python scheduler is used through
the already-reviewed bounded-retention driver lab/runner_stream_20261017.py;
its original strict codec verifies the entire immutable hash chain and the
worker process supervisor is unchanged. Each attempt includes its driver and
source/input identity receipts. The full multi-gigabyte journal is not bundled.

Native control command per crate:
CARGO_HOME=/workspace/tools/cargo RUSTUP_HOME=/workspace/tools/rustup/state
CARGO_TARGET_DIR=/workspace/crypto-rootfinding-pr/target
/workspace/tools/cargo/bin/cargo test --release --locked --offline
--manifest-path research/<variant>/Cargo.toml -- --test-threads=1
The current host has cached dependencies; offline build elsewhere requires the
same Cargo.lock dependencies. Frozen Linux binaries are also in the bundle.
Do not execute comparison workers outside their declared budgets.

The pair and tangency workers accept experiment id, fixture JSON path and
four-profile JSON path. The fixed original fixture is
research/advance/instances/n83-d18-m5-T001-offset0.json, with input SHA pinned
by every SOURCE_MANIFEST/preregistration. Presets/field size are not a hardness
certificate. Source tuples, auxiliary fibers and whole original all-roots are
not output by these norm-family component workers.

Code and native tests are added to the existing source-only PR2049. Reports,
PDFs, plots and raw evidence remain local pending explicit publication approval
after automatic approval review rejected the benchmark payload upload. No
alternative format is uploaded to bypass that rejection. Conductor is not
installed and .conductor/WORKFLOW.md is absent in the available checkout; no
reservation is claimed. Repository-wide CI is not represented as passed.

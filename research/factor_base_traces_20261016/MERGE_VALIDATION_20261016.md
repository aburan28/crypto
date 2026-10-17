# Merge maintenance receipt

This is PR #1975 integration validation, not a research run or a performance
measurement. The published report, catalog, trace/cache archive and original
validation logs remain unchanged.

The creation head is `0b4a7451c7bb080642f03fccb641066dd11f795d`.
Main merged into it is `1e03f1d0b9c77d72ee281561fd7f58405c4ae0b0`.
The Copilot request received no response or pushed resolution after about
4 hours 51 minutes; the takeover was documented [before local resolution][takeover].
There was no rebase or force push.

Only `store/index.json` conflicted: the two branches inserted different names.
The [native merge helper][helper] uses the repository's `Index::load`, `put`
with replacement disabled, and atomic `save`. It refuses any shared name whose
complete entry differs. Both sides' 10 shared entries agree; their union has
17 entries, retaining all five main-only judge entries and both trace entries.
No remote object was uploaded, replaced or removed during this maintenance.

The existing metadata-only `harvest_native_factor_base_history` checks all
48 native identity records. The existing `curve_id.alias_map_text` regenerates
the alias map from the combined registry without point generation or counting.
Exact-integer JSON comparisons confirm that every registry record and alias
from both parents is unchanged: 2,160 models and 5,289 alias entries in total.
All 24 published non-store Git blobs retain their creation hashes. The only
published file whose hash changes is the additive store index.

The local final archive still has SHA-256
`923579b932103882f85b4475118a5edc7404fa4138823c1bdd59a633ca6f9c9c`.
The [original receipt](VALIDATION.md) retains the prior owner-store verification;
this check does not claim a new authenticated S3 verification.

No historical #1898 evidence or #1963 draft card is modified. No historical
CPython generation, experiment dispatch, point-log recovery or new timing was
run. The card remains draft, unrun, zero-funded and unleased. The qualitative
[diagram](flow.svg) and published PDF need no regeneration because their
contents, data and scope are unchanged.

Validation results and exact command receipts are recorded separately below.
Conductor remains unavailable (`command not found`), and its workflow file
at the merged main returned 404. Neither check is reported as passed.

| Check on the combined working tree | Result |
| --- | --- |
| `cargo check --all-targets` | Passed |
| Release trace tests | 4 passed, none failed or ignored |
| Strict release Clippy for the trace binary | Passed, warnings denied |
| Production trace build and source-backed verification | All 384 traces/caches verified |
| Unchanged store integration tests against rebuilt CLIs | 5 passed, none failed or ignored |
| Existing metadata suites | 5 identity, 3 normalization, 3 archived-label tests passed |
| Naming and metadata preservation checks | Passed; 0 naming problems |
| `cargo test --release --lib -- --test-threads=1` | 4,154 passed, 0 failed, 117 existing ignores |

The complete library invocation passed with networking enabled; no library
retry, assertion change or new ignore was needed. This new maintenance run
does not erase the initial authoring failures and retries in the original
receipt. Raw output is retained in [the maintenance validation directory][logs],
including the [full library log][library]. Commands, toolchain and environment
are in [commands.txt][commands]. GitHub CI and required reviews must still pass
on the published exact head before merge.

The optional `cargo test --release --test store_cli` command was interrupted
with status 130 during compilation because it also rebuilt unrelated binary
targets. It ran no tests and is not called a passing command. Its partial log
is retained. After it rebuilt the `crypto` and `ic` production binaries, the
unchanged `tests/store_cli.rs` was compiled with the same Rust 1.98 toolchain,
`--test -C opt-level=3`, and `CARGO_BIN_EXE_crypto`/`CARGO_BIN_EXE_ic` pointing
to those rebuilt binaries. Every test is run, with unchanged assertions.
All five passed. Their S3 transport uses the existing local stand-in executable;
this is not an additional live S3 check. Separately, native `ecbench store verify`
rehashes both previously cached trace bundles successfully against the combined
index, with no `--remote` request.

The existing curve-identity suite passed all five tests, and the archived-label
suite passed all three. Prime-model normalization initially reported two
failures and one error because the tracked historical `out_seed7.json` fixture
was outside the sparse checkout. That failure log is retained; the committed
fixture was materialized unchanged before retrying. No test source, assertion
or input bytes were modified. All three normalization tests passed on retry.

[takeover]: https://github.com/aburan28/crypto/pull/1975#issuecomment-6214922026
[helper]: validation/merge-20261016/merge_store_index.rs
[logs]: validation/merge-20261016/
[library]: validation/merge-20261016/library.log.gz
[commands]: validation/merge-20261016/commands.txt

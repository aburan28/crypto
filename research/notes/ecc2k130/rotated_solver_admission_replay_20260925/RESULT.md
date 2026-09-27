# Solver-admission CI replay after legitimate source drift

The 2026-09-25 static admission result remains **zero solver arms admitted for
its hash-pinned interfaces**. The original runner, freeze, input, result and
receipt have not changed. Its `ci_replay.py` failed on main at `8b1e2d3`
because the historical `koblitz_groebner.rs` hash was
`bd95ee98ad095f93aec5a2b79e7f4840dcf4473db174418aef9b77ebae79c015`
while the current source hash was
`429be4bcfb520e12481d1f95a7c54862cab3f691ab585adbb60be1ba7d033518`.
All ten reference bytes matched the original `8b640f3` commit when checked
against that commit. The builder changed later to use checked layout arithmetic
and to report unsupported inputs as inconclusive. This was a replay-environment
failure, not a new PDP or changed historical blocker result.

The new replay reconstructs a temporary historical tree from the nine
still-matching reference files and the committed 52,043-byte gzip snapshot of
the original builder. It verifies every decompressed/reference SHA-256 against
the untouched original `INPUT.json`, then runs the **original** `ci_replay.py`
both without and with the original final receipt. Both pass. The separate
current-tree check reports the direct builder, memoized template builder,
reuse wrapper and Gröbner frontend still have one shared basis, checked
arithmetic and a 64-variable limit; the raw n13/m5/d2 layout is 49 bits and
n19/m6/d2 is 88 bits. The latter remains unsupported by this generic path,
which is inconclusive rather than a refutation. This narrow result does not
assess newer O-aware exporters or solver runtime.

| Check | Recorded outcome |
|:--|:--|
| Historical frozen replay and original receipt | `HISTORICAL_REPLAY_PASS` |
| Current direct, reuse/template and frontend contract | `CURRENT_GENERIC_CONTRACT_PASS` |
| Python fail-closed controls | 3 passed; corrupt old hash, removed width guard and path traversal rejected |
| Hosted n19/m6/d2 unsupported regression | 1 passed, 0 failed |
| Hosted existing unsupported-input regression | 1 passed, 0 failed |

The first full PR head `9388d3d` preserved a **rustfmt failure** in two call
layouts in the new regression; formatting-only commit `a297bb4` fixed it.
The [hosted admission job on `a297bb4`](https://github.com/aburan28/crypto/actions/runs/36165271175/job/108171368325)
passed and its log confirms both Rust tests actually ran. On the same head,
rustfmt, both clippy jobs and workflow syntax passed. A local `cargo test`
attempt exited 101 before compiling because the sandbox could not resolve
`index.crates.io`; this is an environment failure, not a failed Rust test.
The hosted CI is the Rust validation. The evidence commit that adds this note
will receive its own exact-head CI check.

Commands: `python3 replay.py --mode historical`, `python3 replay.py --mode
current`, `python3 -m unittest discover -s . -p test_replay.py -v`, and the
two `cargo test --lib pdp_admission_... -- --nocapture` commands in the
path-filtered workflow. Compact stdout and the exact current source digests
are in [`evidence/`](evidence/). No solver process, ECDLP timing, full-DLP
cost, matched-rho ratio or ECC2K-130 transfer result is claimed.


## Round 2: a second reference file drifted, generalized to a snapshotted set

Two days after round 1, the original, untouched
`ci_replay.py` failed again on `main` at `6c0324c9`: **two** reference
files no longer matched `INPUT.json`, not the one round 1 fixed. Reproduced
directly (no CI needed):

```
python3 research/notes/ecc2k130/rotated_solver_admission_20260925/ci_replay.py
```

| reference key | path | frozen sha256 | `6c0324c9` sha256 |
|:--|:--|:--|:--|
| `koblitz_builder` | `src/cryptanalysis/koblitz_groebner.rs` | `bd95ee98…79c015` | `1590fd29…881a63c0` |
| `sat_example` | `examples/koblitz_s5_sat_instance.rs` | `c2bc8b05…2115f09` | `4ce1a95c…d550c8eaf` |

`koblitz_builder` had already moved once, past round 1's snapshotted
`429be4bc…` value, to a third hash; `sat_example` drifted for the first
time. Both files sit under `src/cryptanalysis/` and `examples/`, ordinary
locations unrelated cryptanalysis work keeps editing, so patching this one
file again would only defer the next break. Round 2 (`PROTOCOL.md`) instead
snapshots all five reference-file keys whose paths live under `src/` or
`examples/` (`binary_semaev`, `koblitz_builder`, `sat_example`,
`solver_adapter`, `wdsat_adapter`), each independently verified against
`INPUT.json`'s original hash and stored as its own
`historical_snapshots/<key>.rs.gz`; the round-1 single-file
`historical_koblitz_groebner.rs.gz` is superseded by
`historical_snapshots/koblitz_builder.rs.gz` (identical bytes). The five
evidence-file keys inside this thread's own `rotated_*_20260925/`
directories remain verified live, unchanged from round 1.

Merging round 1 onto current `main` also exposed a second, independent
problem: `main` had added a required `derived` field to
`FrobeniusFactorBase` since round 1's regression test was written, so the
sentinel factor-base literal in
`pdp_admission_n19_m6_d2_layout_is_unsupported_not_refuted` no longer
compiled (`E0063`). Fixed with `derived: Default::default()`, matching the
pattern used elsewhere in the same file; this is a struct-literal update
for an added field, not a change to the test's assertions or to
`groebner_decompose`'s admission logic.

| Check | Recorded outcome |
|:--|:--|
| `replay.py --mode freeze` | `NEW_AND_ORIGINAL_FREEZE_PASS` |
| `replay.py --mode historical` | `HISTORICAL_REPLAY_PASS` (all five snapshotted keys, all five live keys) |
| `replay.py --mode current` | `CURRENT_GENERIC_CONTRACT_PASS` |
| `python3 -m unittest discover -s . -p test_replay.py -v` | 4 passed (mismatch rejected for every snapshotted key and for a live-verified key; width-guard and path-traversal controls unchanged) |
| `cargo test --lib pdp_admission_n19_m6_d2_layout_is_unsupported_not_refuted -- --nocapture` | 1 passed, 0 failed |
| `cargo test --lib pdp_admission_unsupported_is_not_a_refutation -- --nocapture` | 1 passed, 0 failed |
| `rustfmt --check --edition 2021 src/cryptanalysis/koblitz_index_calculus.rs` | clean |
| `cargo clippy --lib --tests -- -D warnings` | clean (one pre-existing unrelated `unknown_lints` warning on an unrecognized clippy lint name, not from this change) |

All commands ran locally to completion; this round's local sandbox could
reach `crates.io` (round 1's could not), so no result here depends on
hosted CI to have actually executed. No solver process, ECDLP timing,
full-DLP cost, matched-rho ratio, or ECC2K-130 transfer result is claimed,
and the 2026-09-25 admission verdict is unchanged: zero solver arms
admitted for its hash-pinned interfaces.

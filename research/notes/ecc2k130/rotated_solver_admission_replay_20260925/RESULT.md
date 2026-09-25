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

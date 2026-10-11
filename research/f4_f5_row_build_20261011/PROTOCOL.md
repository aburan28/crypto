# Boolean F4/F5 row construction audit

Requested scope: optimize the existing F4/F5 implementations in
`aburan28/crypto` or `aburan28/cryptanalysis`. This iteration targets the shared
Boolean Macaulay row builder and the Rust Boolean F5 criterion in `crypto`.

| Requirement | Evidence to retain |
| --- | --- |
| Improve F4/F5 construction overhead | Matched native baseline/candidate allocation records |
| Preserve correctness | Exact reduced rows, pruning reports, word-XOR counts, and library tests |
| Reproducible comparison | Frozen seeds, driver, input/output hashes, source hashes, compiler and host |
| Deliver reviewable changes | Feature-branch commits and PR, or exact access blocker |

Hypothesis: compacting product parity in place removes a second F4 scratch
buffer; deduplicating columns during collection avoids materializing and
sorting repeated monomials; XOR-packing F5 criterion products directly and retaining scratch rows
removes repeated temporary allocations and sorting. A monotone leading-word
cursor avoids rescanning known-zero prefix words during criterion reduction.

The reference is repository revision
`c95a14de9bfa2711643a9318f3728d5d1109367b`, with the pre-existing duplicated
assertion line in `koblitz_factor_base_search.rs`, duplicated constructor and
field in `koblitz_index_calculus.rs`, and duplicated assertion fragment in
`prime_fast.rs` repaired identically for both arms. The surviving constructor,
helper and assertions retain the existing implementations. Further duplicated
imports, coordinate helpers, legacy return fragments, an Artin-Schreier return
fragment and two `SolveStats` fields also block the baseline; remove the
duplicates while preserving the existing fast implementations and field
defaults. These prerequisite repairs are separate from the optimization.
Record any further
prerequisites before performance edits. Save the
baseline binary before modifying the kernels. Use the same diagnostic source
and toolchain in both arms.

The test build additionally exposed an undefined current-parser fixture and
a removed sorting utility still referenced by the historical stable-sort
regression. Restore the fixture declaration and a test-only radix utility;
retain the original assertions. These repairs are outside production calls.
After formatting the diagnostic driver, rebuild it against the saved baseline
library and retain an A/A audit; require exact allocation and output agreement
with the first baseline before comparing the candidate.

Frozen workloads: eight planted quadratic Boolean systems, `(n, m, degree)`
in `(8, 8, 3)`, `(8, 8, 4)`, `(10, 10, 4)`, `(10, 10, 5)`, with seeds
`0x2545f4914f6cdd1d` and `0x9e3779b97f4a7c15`. The second seed is a fresh
holdout relative to the existing F4 benchmark's seed. Run three repetitions,
one Rayon worker, a warmup per case, and retain every row. No adaptive seed or
workload selection. The native driver generates inputs outside the counted
interval and audits complete F4/F5 calls and criterion construction separately.

Primary units are allocation/reallocation requests and their requested bytes;
requested bytes are cumulative traffic, **not peak memory**. Both include
packing, elimination and result construction for a complete call. Fixture
generation, warmup, output hashing and verification are outside the interval.
This is a solver construction diagnostic, not an ECDLP measurement. The WDSat
regression and matched-rho workload do not measure this standalone algebra
contract. No full-DLP, scaling, general-field, hardware throughput or CPU
wall-time gain is asserted. Full-pipeline cost and speedup remain unset.

Success requires exact output/pruning/count agreement on all frozen inputs
and holdouts, passing relevant correctness tests, and fewer F4 scratch
allocations plus at least 20% fewer criterion allocation requests. Retain any
allocation regressions and qualify the result by workload. Stop and correct
any correctness mismatch; do not relax the acceptance conditions. No wall
clock is collected by the audit driver because this host has no demonstrated
exclusive CPU isolation. Retain five additional alternating baseline/candidate
pairs as a reproducibility control, with the same three repetitions per pair;
compare each against the primary JSONL files byte for byte. A controlled
wall-time comparison remains a separate measurement obligation.

Commands: `cargo build --locked --release --example f4_f5_row_build_bench -j 2`;
save the binary; run baseline and candidate with argument `3`. Run
`cargo test --locked --release --lib` and the touched example's tests after
the optimization, plus focused F4/F5 suites. Preserve failure logs. Full-library
verification has a 20-minute execution budget. If a pre-existing test stalls,
continue the remaining tests independently, preserve the full incomplete log
and timeout, and do not claim full-suite success or weaken the stalled test.

Workflow limits at setup: `conductor check` could not run because the command
is not installed and this checkout has no `.conductor/WORKFLOW.md`. The work
uses a dedicated feature branch and explicitly named paths. `gh api user`
returned HTTP 401; PR publication requires functioning repository credentials.

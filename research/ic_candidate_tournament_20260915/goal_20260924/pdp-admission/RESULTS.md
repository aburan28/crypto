# PDP outcome correctness results

The [registered controls](PROTOCOL.md) passed locally: 22 release Rust tests,
121 Python harness tests with two Linux-specific skips, and compiler checks for
`koblitz_unknown_scalar_panel` and `ic_tournament_worker`. These are toy
correctness/accounting controls, not a new improvement round or a speed claim.
The 39 site-build tests also pass. The first CI head exposed two relative
scoreboard evidence links that did not survive publication; both now use the
repository's absolute GitHub evidence-link convention.
Linux integration remains the PR acceptance gate; its receipt is linked in the PR.

The defect was in reporting an unattempted solve. The narrow Gröbner and Crossbred
frontends used default statistics when their encoder rejected infinity, fewer
than two summands, or an oversized Boolean layout. The relation driver interpreted
the absence of both an answer and exhaustion as a mathematical refutation. SAT
already reported unknown but omitted these rejections from its unknown counter.

| Input/outcome | Corrected classification | Required evidence |
|---|---|---|
| Input cannot be encoded | `Unsupported`, legacy incomplete flag set | No algebraic refutation or solver work credited |
| Supported input exhausts its budget | `Unknown`, unsupported flag clear | Incomplete search remains in the ledger |
| Completed solve has no lifting root | `Refuted` | Existing enumeration/UNSAT controls still pass |
| Completed solve produces a witness | `RelationFound` | Summands independently checked in the group |

The new tests exercise known decompositions of infinity, a one-summand target,
and a nonidentity 16-summand target whose narrow layout is too wide. The F4, F5,
inherited F4, Crossbred, native-XOR SAT and CNF SAT frontends all preserve the
unsupported outcome. Four pipeline configurations retain eight attempts each
without creating a false refutation. A supported zero-node solve remains budget
exhaustion. Both the direct and memoized encoder paths reject integer overflow.

The existing 39-target comparison now executes F5 as well as the default Gröbner
and SAT paths, comparing existence with enumeration and verifying returned sums.
Existing genuine refutations, three-summand chains, inherited reduction,
Crossbred search and template-equivalence controls also pass. This finite corpus
does not establish completeness or performance on every field or encoding.

## Reproduction and retained evidence

The local run used macOS arm64, Rust 1.93.1, the checked-in dependency lock and one
test thread. [source-control.json](source-control.json) binds the changed source
files and lockfile to their SHA-256 hashes and the accepted base commit.
[local-regressions.log](local-regressions.log) contains every command and result;
[python-tests.log](python-tests.log) records the harness outcome. Compilation
emitted existing platform/dead-code warnings; this is not an all-warning lint pass.
An initial test compilation used a nonexistent negation method; it was corrected
to the existing `point_neg` helper before these passing runs. That failed build
is retained in [initial-build.log](initial-build.log).

```sh
cp research/ic_candidate_tournament_20260915/ci/Cargo.lock Cargo.lock
cargo test --locked --release --lib pdp_admission_ -- --test-threads=1
cargo test --locked --release --lib all_three_oracles_answer_every_target_identically -- --test-threads=1
cargo test --locked --release --lib undecomposable_targets -- --test-threads=1
cargo test --locked --release --lib chained_s3_handles_three_summands -- --test-threads=1
cargo test --locked --release --lib the_two_algebraic_engines_agree -- --test-threads=1
cargo test --locked --release --lib cryptanalysis::crossbred::tests -- --test-threads=1
cargo test --locked --release --lib cryptanalysis::polynomial_reuse::tests -- --test-threads=1
cargo check --locked --release --example koblitz_unknown_scalar_panel --example ic_tournament_worker
python3.12 -m unittest discover -s research/ic_candidate_tournament_20260915 -p 'test_*.py'
```

The Linux workflow pins Rust 1.94.1 and additionally runs existing descent/rho
regressions plus complete-solve/profiler integration. Raw controls are not
performance samples: online time, cold time, instructions, S and speedup stay null
in the source-control record. No archived measured source was modified.

Generic `RelationCollector` and individual-log dispatch still discard per-attempt
solver statistics. Public-point timing, exclusive phases, actual engine/base/LA
identities, ordinary-query yield and independent stage certificates must be wired
before these engines enter the ranked panel. The goal remains active: round one
retained the incumbent, with two bounded improvement attempts remaining.

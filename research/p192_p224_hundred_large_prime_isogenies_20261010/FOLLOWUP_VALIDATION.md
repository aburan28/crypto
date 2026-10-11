# Required library gate follow-up, 2026-10-10

The scientific campaign is unchanged: 100 independently verified maps across
two prime degrees, with 50 linked edges in each public P192/P224 class. The
frozen producer is source 49a280dc6; the standalone 449-model archive is source
3c39c5f23. A fresh native audit again matched all 891 files in MANIFEST.json.
The archived executable again passed a two-step chain replay with complete=true.

The previous full-library run stopped before a final summary. Its retained
root-complete-prefix.log records at least ten failures; it is not a passing
suite. The old temporary log and process were unavailable when this follow-up
started. The committed checkpoint and every original map receipt remain intact.

## Reproductions and repairs

The saved release test executable reproduced three degree-71 construction
failures and the isogenous-model factor-base failure in isolation. The
linearised-kernel elimination packed field images and preimages into one u64,
discarding coordinates above bit 63. Packing both words into u128 restores the
admitted wide-field range. A regression checks X^2+X, the identity and the zero
map at degree 71, including all high-coordinate basis vectors.

The packed orbit walk applied Frobenius even when the model's coefficients
were not fixed by it. Models without that endomorphism now use the existing
general orbit walk, whose step is the model's declared identity operation.
The existing fast point lifts remain available, and ordinary Koblitz models
retain their packed orbit path. Both existing isogenous-model tests pass.

The j=0 stage-accounting fixture failed on the second isolated repetition.
It demanded full relation rank with only two extra random rows. The fixture
now uses 32 redundant rows and checks the corresponding exact raw-row count;
the solver and its failure behavior are unchanged. This remains a stochastic
fixture, rather than a deterministic rank guarantee.

The three cover-solver failures, the Edwards four-point failure and the gate
benchmark failure passed in isolated runs. Their full-run failures may depend
on host contention, timing or shared process state; the isolated passes do not
prove which cause applied. Their arithmetic and acceptance tests were retained.
The full rerun uses two concurrent test threads with localhost access.

## Checks completed so far

- 57 wide-arithmetic tests pass, including the new high-coordinate regression.
- The degree-71 curve constructor and subgroup-order check pass.
- Both isogenous-model tests pass.
- The revised j=0 fixture passes 50 consecutive isolated repetitions; all
  individual test outputs are retained in one concatenated log.
- All 12 active cover-module tests pass in isolation; its existing ignored
  exhaustive test remains ignored.
- The Edwards four-point test and gate benchmark pass in isolation.
- The original artifact audit matches all 891 recorded files.
- The archived 449-model executable replays two saved P224/1117 steps.
- Formatting and git diff --check pass for the follow-up changes.

The full command is:

```sh
cargo test --release --lib --target-dir target/isogeny-root-gate-20261010 -- --test-threads=2
```

Its durable log is evidence/validation/root-repaired-full-20261010.log and its
terminal status is written to root-repaired-full-20261010.exit. Until that
command finishes successfully, the required library gate remains pending and
the new commits are not pushed. Hosted CI and ReleaseMe publication follow
their existing gates. These follow-up logs are additional integration evidence
bound by their Git commit; the frozen scientific manifest is not regenerated.

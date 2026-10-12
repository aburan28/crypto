# Strict control environment repair

The first PR #1093 head
`ddbfe57fe1a9577f477a5cac620783b2c606f7ef` passed Linux integration, Python
harness, candidate controls, both/pairinv archival transport and formatting/
Clippy. Its [scaled control job](https://github.com/aburan28/crypto/actions/runs/36782205723/job/110114956233)
failed in the newly added isolated native tests, after the controlled generic
build succeeded. Three tests rejected before any target attempt:

```text
undeclared algorithm environment override: IC_SOURCE_MANIFEST_SHA256
```

The workflow's earlier archived-producer step intentionally writes that legacy
producer provenance into `GITHUB_ENV`. The generic worker admits only its own
frozen default runtime environment and correctly rejects inherited `IC_*`
variables. The mode-rejection test passed because that guard precedes setup;
the success/capped/log-corruption controls stopped at the environment guard.
This is a test-command integration failure, not an observed F5 search outcome,
zero-yield cell or source-bound production invocation. The failed run and its
uploaded control artifact remain available; it was not redispatched or erased.

The repair invokes only the new native test command through
`env -u IC_SOURCE_MANIFEST_SHA256 RAYON_NUM_THREADS=1 ...`. It leaves the earlier
archived provenance intact for subsequent driver checks. No other environment
variable is removed or allowed, no native/source guard is relaxed, and no
historical manifest, evaluator or measured registration is changed. A new PR
head must pass every applicable check before merge.

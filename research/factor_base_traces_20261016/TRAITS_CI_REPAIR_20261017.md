# Current-head native catalog repair

PR #1975's maintained head `1a83e379ade9c2cd4c52e426cf4f4777d863e6c8`
merged main `b59187e5549adb02a43c977f4bea9979a46bf3a3`. This repair is a
direct descendant of that shared head. It does not overwrite either older
archive branch, merge main again, rebase, or force-push.

The current-head [traits check][failure] failed because the generated file
was stale. Normal Rust CodeQL succeeded on that head, so the earlier scan
timeout is historical, not its current blocker. The neutral CodeQL summary
warned about missing baseline configurations; it is not itself a passing
scan. The actual five language-analysis jobs succeeded.

## Regeneration and exact scope

A fresh Rust 1.98 production `curve_traits` build regenerated all 2,336
records with the unchanged **262,144** factoring budget and all six exact
tracked derivation inputs. Comparing that output with the committed file
reproduced the freshness discrepancy. The output was generated anew, not
copied from an older archive.

Exactly one field changes, for `icv1-fp24-t925-8af643e2`:
`endomorphism.levels[0].method` changes from
“Frobenius is a scalar on E[ℓ]: above the floor” to
“Frobenius is a scalar on E[ℓ]: the crater”. Numerical values, array order,
headers, operation budget, and every one of the 2,288 main-parent records
remain unchanged. This is generated description consistency, not a new
mathematical or performance result.

The final fixture was committed as
`a4dd1db93fefc35ba5210a10288c83b69762f695` before test compilation.
Its SHA-256 is
`45e928f3290f1e3f76eb887e233658492d91fff29c6212c658b68f2b0fe198c2`.
A frozen manifest pins all 1,374 tracked native build inputs, including
compile-time fixtures and configuration. Earlier snapshot test counts do
not validate this source; its validation is recorded separately below.

## Validation

Validation passed on source commit `a4dd1db9`, with networking enabled:

- Full `cargo test --release --lib -- --test-threads=1`: **4,239 passed,
  zero failed, 117 existing ignores**, one complete command (4,356 tests).
- Independent production `curve_traits check`: the final file is current.
- Production trace build, all four native trace tests and verification of
  all 384 source-backed traces/caches, including signed lookups.
- `cargo check --all-targets` and strict release trace Clippy with warnings
  denied.
- Final frozen-input and preservation audits.

The unit binary SHA-256 is
`5874319da87e0ff14df739153636ea7ce53489cd85882c2392d2169077bab47a`;
it contains the exact final fixture bytes. Test compilation completed with
the source and final fixture already pinned. The actual serial test run
started after catalog checking and compilation finished, without competing
catalog computations. Two Rayon workers and two Cargo build jobs were
configured; command details and resource limits are retained. No counts
are borrowed from an earlier run.

Exact commands, actual exits, complete logs and source/input manifests are
retained in `validation/traits-20261017/runs.tar.xz`. Its per-member checksum
manifest and fixed-metadata repacking instructions are beside the receipt.

## Preservation and limits

The read-only audit verifies all 55 prior non-store payload blobs and 27
frozen historical-evaluation blobs unchanged, both parents' complete 2,160
model / 5,289 alias / 23 store entries, browser identities and frozen
leaderboard fields. The #1963 planning card remains identical to main:
draft, unrun, zero-funded, unleased, without a frozen protocol, checkpoint
or runnable command. All previous validation receipts and failed attempts
remain unchanged.

Sparse materialization was followed by restoration of the unchanged,
tracked off-cone isogeny-route input JSONs; they were not regenerated.
A transient command-connection interruption lost the original session
handle, but process inspection confirmed compilation continued. No duplicate
native generation was started. Conductor was unavailable (127), and the
pinned `.conductor/WORKFLOW.md` returned 404. Neither is a pass.

Historical CPython generation was not rerun. No point logarithms, new
candidate traces, native timings or end-to-end gains are claimed. There was
no S3 operation or renewed live verification of the existing owner-store
bundle. Its historical 384 bases, 52,876 representative rows, 48 registered
ICV1 models and FB1 identities remain frozen. Future measured comparisons
retain symmetric native-operation pricing and matched external references;
verified engineering improvements need not beat rho to count as successes.

Applicable checks and review gates must complete on the published final
head before an exact-head-guarded merge. No test, assertion, ignore, limit,
source coverage or workflow was weakened.

[failure]: https://github.com/aburan28/crypto/actions/runs/39006728222/job/117320603054

## Concurrent shared-branch update

Before publication, another session merged newer main
`a95fcd95305af6ac98773011e7730203e927ae70` into the shared feature branch as
`b6f55641aa9650c19abe78c9dd18f96dd92f09db`. That head already contains
the diagnosed “the crater” description, and its catalog now has 2,344
records at the same budget. It also changes point-group native code and
other sources. This completed local snapshot is therefore archived
separately, without overwriting or merging it over the newer shared head.
**The 4,239-pass result does not validate the newer shared-head sources.**

A read-only comparison of that head finds all 53 immutable prior payload
blobs and all 27 frozen historical-evaluation blobs unchanged. The two
additive registry/alias containers have upstream additions: 2,168 models
and 5,312 aliases. All 48 historical models remain identical, every
main-parent model matches main, all previous alias mappings survive, and
all 23 store entries from both parents survive unchanged. One inherited
73-bit model has upstream provenance/representation enrichment matching
main; its existing mathematical field/curve identity is unchanged. This
is metadata inspection, not a new native replay or live S3 verification.
Current-head native freshness, CI and review gates remain pending.

The first naming command failed when an off-cone Git diff needed a
promisor fetch. Repeating it with normal authorized network access and
no-op fetch negotiation passed with zero problems. Packaging was first
attempted before that interrupted step had produced the run directory;
the resulting incomplete unpublished archive was rejected and retained
separately as a failed attempt. Only the subsequently hash-verified,
round-trip-tested archive is accepted. No failed or partial outcome is
described as a pass.

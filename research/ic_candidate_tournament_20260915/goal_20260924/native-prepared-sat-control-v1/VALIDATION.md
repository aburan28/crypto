# Native controller implementation controls

The implementation is tested on the already retained synthetic n17 preparation,
ANF/CNF exports and disclosed model. No test launches the SAT exporter or
CryptoMiniSat. Test processes exercise transport deadlines with ordinary shell
processes only. Model substitution is a correctness control, not solver search.

Local controls on macOS ARM64 with Homebrew Rust 1.93.1:

- Six producer controls pass: complete retained-model recovery without a scalar
  input, preservation of the earlier inconclusive attempt, all 767 individual
  expanded-model bit corruptions rejected, invalid model retained as a failure,
  partial/duplicate/unterminated/oversized model rejection, native status/exit
  consistency, preparation/config/source-layout gates and exact exclusive
  phase summation.
- The independent query sampler matches native rand 0.8 on 128 scalars for
  each of five seeds, including zero and `u64::MAX`. The retained trial-one
  coefficient pair is `(32326,42888)`.
- The native custody/transport controls check create-only output, immutable
  tree mutation rejection, retained source metadata, exit preservation,
  process-group timeout cleanup and accepted archive extraction of all fourteen
  files without executing those binaries.
- The complete independent-auditor control reconstructs the producer's
  query sequence, checks source model and group relation, derives the scalar
  and rejects changed query coefficients, scalar, witness, model validity/hash,
  timing, outcome, source file and attempted promotion. It uses explicitly
  labelled mock process receipts and does not enter the source-bound admission
  wrapper. A passing mock control cannot admit scientific runtime execution.

The first archive extraction control failed on its checksum-header parsing;
[the retained failure](first-custody-test-failure.txt) records the cause and
correction. The corrected five custody/full-auditor controls pass. Neither that
failure nor its rerun consumed a scientific invocation.
The first full local harness run passed 47 harness and four worker unit tests,
then hit the unchanged R05 regression's Linux `/bin/true` placeholder on macOS.
That test now selects macOS `/usr/bin/true` for its preflight; frozen outputs
still prevent any binary invocation, and its expected pin remains unchanged.
The original failure is retained alongside the archive-parser failure.
Source-freeze validation also stopped before registration when offline Cargo
vendoring tried to unpack locked `openssl-probe 0.2.1` into the restricted user
cache. The [complete log](first-freeze-vendor.log) and
[build-step receipt](first-freeze-vendor.receipt.json) retain exit 101 and the
original inputs. The failed capsule remains unregistered and unexecuted.
Prepare the locked package cache before offline freeze; use a new output
directory. This is a build/setup failure, not a new solver outcome or retry.
The next source-freeze validation completed from a clean native snapshot:
5,957 immutable files, including 152 MiB of offline vendor sources, producer
and checker binaries, accepted native assets and all build receipts. That
capsule is unexecuted and remains build validation; its source revision predates
the final portable test and launch/admission guards. Freeze and publish the final accepted
registration separately before a scientific dispatch.

The first Linux CI pass reached the new watchdog test, where `/bin/sh` was a
symlink and the regular-file source gate correctly rejected it. The test now
canonicalizes that ordinary shell fixture; the scientific executable gate
still rejects symlinks. The complete first failing Linux/harness logs are
retained. This is a test portability correction, not a solver result.
Final native launch controls require the consumed execution, live producer
parent, worker-start marker and registered query cap before a helper can exec
a role binary. The source-bound audit wrapper explicitly rejects labelled mock
receipts; mock controls still exercise its arithmetic verifier without admitting
execution. These launch-context and complete-auditor controls pass locally.
The changed preparation helper signatures were also formatted after the CI
formatter flagged them; unrelated recursively visited files remain unchanged.

The CI workflow runs these native controls on Linux x86-64 and macOS ARM64,
plus the existing independent preparation and source replay tests and harness
output regressions. Actual solver hardware admission remains macOS ARM64 only.
CI tests no solver search, native yield or timing claim. Exact-head CI results
must be checked before merging; the protocol's actual frozen build/registration,
sole dispatch and independent admission are still pending.

The full F4/F5-plus-SAT goal remains active. This implementation control does
not complete either the fresh paired comparison or the globally best pipeline
search, and it makes no cross-method speedup claim.

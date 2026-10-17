# Clean native retake at crypto015fa7bb504e

This is the clean source/control checkpoint for TASK-20261016-de577e.
Source: 015fa7bb504ef63953e0e04e464bb266f4f207bf.
PR: https://github.com/aburan28/crypto/pull/1991.

Reference and candidate were built from exactly this clean commit with the
same rustc1.93.1, Cargo.lock, aarch64-apple-darwin target and release profile.
Only the explicit cached-fixed-modulus-reduction feature differs. Both
installed version records say git_dirty=false and name the expected policy.
The four native executables remain outside Git at the exact paths/hashes in
native-manifest.json. No native executable is inside native-controls.tar.gz.

| Arm | isogeny SHA256 | Native suite |
|---|---|---|
| Reference | 2ff7d57174cda2d0e90088eb6ce217cbd1826e7697ed12ec7770d3d8f395d379 | 128/0; one parent-invoked helper ignored |
| Candidate | 305bdff7db7bc25a2ccdce3f44e71c3d47269e76afa0ecd1703803365ab9b207 | 128/0; one parent-invoked helper ignored |

Compiled and installed bytes match for both isogeny and the research companion.
Six complete mathematical verification JSON pairs are byte-identical, with no
excluded fields: public P192/13, P521/7, P224/1471, and three P192 trace
controls. All 20 transports per map are still required. Nine adversarial trace
pairs and both rehashed P224 cache negatives have identical status/reason/score.
Each arm retains 16 lost/zero-RSS injection receipts (two worker shapes,
before/after acknowledgement, two independently compiled parent binaries)
and two deadline receipts. The injected4096-byte RSS is fixture data, not a
real cryptographic memory measurement. Actual public workers retain their
own sampled RSS receipts.

native-controls.tar.gz is 2,049,977 bytes,
SHA25609347071c8857216ae1597898d345307fc2a4086930c4252150cc0c6e2081f76.
All 1261 regular files replayed against the manifest. It includes the exact
compiled Lean object, build/test/proof logs, controls, and the failed
default-parallel root gate plus focused serial diagnostic. The later complete
serial root retake is retained separately in root-gates/. The original first
metadata path error and its source/log were preserved in the sealed archive;
it caused no native rerun.

The final-source default-parallel root gate failed4146/2/117 ignored in two
unchanged JV admission-count tests (observed21 and22 versus thresholds25 and30).
The complete JV group then passed14/0/1 ignored with one test thread.
The full root suite then passed4148/0/117 ignored with one test thread,
without omissions, threshold changes or source repairs. Source/root executable
e6b0366e26794fe703741ba309cd3c7cb168150c827fa2a23ba068038aec7399 is the same
binary used by the earlier default-parallel pass and failure. This establishes
observed variation, not a diagnosed cause. No global parallel-stability claim
is made.

The serialized full-library body reported1047.40s. The enclosing OS receipt
reported1047.87s real,799.80s user,101.98s system. These intervals overlap and
are not added as disjoint costs. Build/test timing and memory receipts are
operational engineering data on a shared Apple M4 Pro host; there is no
performance comparison. Earlier development CPU costs were not measured
completely. compiler-cfg.txt is a default rustc inventory query, not a claim
about release debug assertions; actual release commands, environment and
feature arguments are retained.

A separate native generate control used only the already accepted P192/13
map, SHA2567a1809d56cb774f681e06ecf1dad9e04d97abe70a3eb2738bff9d54affc861d1.
Both worker receipts are PASS with zero unavailable RSS samples; eight
dataset cases were produced. The9.8MiB datasets are byte-identical,
SHA2562b5bb3f2cdd04ddb3ab068e99c7a81e22e3801bcfd0471759ca993ec9e5fd871.
Exact files are retained locally and bound by generation-control.json;
they are not Git payloads. This checks input preparation; generation computes
reference expectations in both builds. It is not a cached-reduction timing
comparison. No measure/analyze campaign or sealed50341 extraction/replay ran.

Native identity IDC1h41e9973f839c85bf and both exact-source Lean lemmas
rechecked successfully with pinned Lean4.35.0. The full formal-series and
quotient-reversal Lean bridges remain open as the proof map states.

Remote custody is LOCAL_ONLY_PENDING_AUTHORIZATION; zero custody uploads
executed. Independent source/resource review has not started, and the source
gate is incomplete. The source claim remains active. Nothing here promotes
a degree, a speed ratio, a curve board or a scientific run. The public release
fix PR1976 merged atbd7eb117f094a336d3e3e8f31facb0c6b86f95bc; its newer main
parent leaves the frozen native closure unchanged. Its official asset release
and the unrelated reports-s3 credential failure are followed separately.


The first two report-rendering attempts failed at Typst's project-root boundary; both sources/logs are retained. The successful command passes --root research/isogeny_world_frontier_20261016/replay-optimization-v1. The resulting two-page report was visually inspected. No native rerun was caused by this presentation repair.

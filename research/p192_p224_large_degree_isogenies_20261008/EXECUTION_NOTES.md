# Execution notes

The first construction launch is retained in `evidence/` and `SEARCH.log`.
It stopped at the P-192 order preflight: the new CLI panicked with
`use fewer limbs for this p`. Its dispatcher chose four Montgomery limbs
for a modulus requiring three. The native dispatcher fix selects the exact
limb count and adds P-192 order/map and all-width boundary regression checks.
The rerun uses a new evidence directory; the first launch is not overwritten.

The required root `cargo test --release --lib` was attempted on the unchanged
source revision c70c32d486a3ac7531fe27f7193d9f09caa58344. It failed to compile
with 643 errors, including duplicate module definitions/reexports and stale
Koblitz/SAT structure fields. The full log is in `validation/lib-test.log`.
These errors precede this search's changes. The independent replay executable
therefore imports the five existing walker source modules directly, with
their original implementations, rather than building the unrelated library
modules. The imported source paths and hashes are recorded in the manifest.
This keeps the independent field, polynomial, torsion, subgroup-closure,
and codomain verifier; it does not make the root library check pass.

Conductor allowed the initial research scope and created task T-65. Its
service later returned HTTP 500 for progress and the dispatcher scope
expansion. The scope expansion was reported in chat before editing the two
dispatcher/regression files. No private conversation was published as metadata.

The initial full worktree checkout hit the disk limit and Git removed its
incomplete tracked checkout. The protocol and runner were preserved, then a
sparse worktree was initialized for source, documentation and required fixtures.
No existing worktree or historical evidence was deleted.

The session interruption stopped the search during P-192 degree 149, before a
completed receipt. No search child survived the interruption. The resumed
native supervisor checks every sealed command and stdout/stderr digest before
reusing its receipt. Unsealed output is moved to a numbered `interrupted/`
directory with its byte hashes and an explicit unknown elapsed time before
retrying the same degree. The frozen degree plan and resource limits are
unchanged. The resume change affects orchestration, not the construction CLI.

The first independent degree-73 replay failed to parse the base registry:
the E-382 entry's final `trace` line was immediately followed by the next
extension-curve entry's `slug`, without an object boundary. The original
registry is frozen as `validation/registry-before-repair.json`; the repair
adds the missing closing brace, comma, and opening brace, preserving both
entries. The failed replay log is retained separately before rerunning it.

The actual construction CLI uses the default optimized release build (level
3). An alternate level-1 build was not used for any recorded construction.
macOS rejected the trial virtual-memory `ulimit`; the actual native supervisor
uses the predeclared 8 GiB sampled resident-memory limit instead.

The feature branch was rebased onto upstream 2fe5cec8a4a9d8d55c8e0abec9fe22852f3a3726,
which includes the registry/dashboard repairs. The original registry input and
failed replay remain frozen. The mandatory root check was rerun after rebasing
and again failed with 643 compiler errors (`validation/lib-test-after-rebase.log`).
No isogeny construction implementation changed in this rebase.

The complete standalone release suite passed 110 tests
(`validation/algorithms-tests-resume.log`), including the new API regressions.
The current search example's two tests passed independently against the actual
construction library (`validation/example-tests-resume.log`). A direct six-test
API run also passed (`validation/api-tests.log`); its first manual build lacked
Cargo's CLI/version compile-time environment, and that build failure is retained.
The earlier incomplete test/build logs ended with the interrupted session; they
are not counted as passing checks.

The native `catalogue-covers` executable imports the unchanged root cover checker,
model conversions, linkage graph, field arithmetic, and prime screen directly.
Its baseline `--check` reproduced all 321 existing cover and graph records, with
zero unsupported or invalid models (`validation/catalogue-covers-baseline.log`).
The native `catalogue-views` executable performs incremental joins only; its
rehearsal reproduced the existing alias map, leaderboard JSON/Markdown/HTML,
and browser JSON byte for byte, preserving measurement bytes. New prime model
rows use the exact registered identity, subgroup, mapped generator, and native
cover certificates. The frozen IC measurement rows are not recalculated.

Conductor recovered during this continuation. The existing T-65 task was
reattached and the API, native study, and affected canonical view scopes were
granted before the catalogue refresh. Earlier HTTP/service failures remain
part of the execution record.

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

Final replay adds an exact rational-map curve-equation substitution check,
using the separate native polynomial implementation. It also checks numerator
degree/normalization and coprimality with the verified kernel denominator.
The known degree-three map test passes and mutations of its numerator and
target coefficient are rejected (`validation/exact-map-tests.log`). This is
an implementation verification of the supplied maps, not a proposed new
algebraic rule. The degree-199 replay passes these stronger checks; the final
whole-search replay applies them to every completed construction. Earlier
partial replay records remain frozen as the checks that were performed then.

The first complete 40-trial summary is retained in `evidence-v2/search.json`.
P-192 produced six maps at degrees 73, 193 and 199, and 21 bounded timeouts.
All 16 P-224 construction trials failed immediately with
`characteristic must exceed l + 1`. The source-order preflight had passed.
The modular-polynomial precondition used the low-word `Field::char()` accessor;
P-224's low word is 1 although its full characteristic is 224 bits. The native
fix adds an exact small-integer characteristic comparison, overridden by
multiword prime fields, and uses it for the Hecke/Newton precondition. The
regression checks a P-224 modular polynomial against the independent linear
algebra route and requires two verified maps at a structurally split degree.

The corrected P-224 trials use `evidence-v3/`. P-192's sealed screen and trial
files were copied byte for byte and reused after command/digest checks; they
are not represented as newly executed constructions. The complete earlier
failure run is preserved. The original CLI executable was saved as
`isogeny-algos-pre-p224-fix` before rebuilding, so both construction versions
can be identified by executable hashes. Resource limits, seeds, source models,
candidate degrees and the frozen protocol remain unchanged.

The corrected CLI built in 7m52s and its complete release suite passed 111
tests (`validation/algorithms-tests-p224-fix.log`). The first rerun, retained
in `evidence-v3/`, stopped at an RSS-monitor shutdown race: P-224's preflight
emitted `PASS` and exited with code 0, but `proc_taskinfo` disappeared just
before `waitpid` reported the exit. The supervisor now permits up to 100 ms
for a transient unavailable sample to resolve; persistent unavailable
monitoring still fails closed, and the original 180-second deadline continues
to apply. Missing sample counts and the grace interval are recorded.
The final rerun uses `evidence-v4/`, with the same sealed P-192 trial files
copied from v2 and fresh P-224 trials. The earlier v3 monitor receipt is
preserved without reclassifying it.

The exact characteristic comparison also forwards through the existing generic
quadratic-field wrapper. Its regression uses the public Mersenne modulus
2^127 - 1 and compares against u64::MAX. This consistency fix does not change
the prime-field construction route used by any recorded search invocation.
The final search CLI was built from source revision
631d857e8414aed95c1b0254f0365bc0ff56f3c5; the test build for this final wrapper
change uses a separate target directory so it cannot replace that executable
while the search is still running.

The wrapper scope's Conductor CLI conflict check was clear. Expanding the
existing task lease then hit its retry-budget limit; the MCP scope-expansion
endpoint failed with an HTTPS/HTTP protocol mismatch. The existing T-65 task
and reserved study/catalogue paths are retained; these coordination failures
are not interpreted as validation passes.

The final full standalone release suite passed 111 tests in 43 result groups
(`validation/algorithms-tests-final.log`), including the wrapper check.
Catalogue roster rendering retains the existing builder's large-prime
coefficient-label convention so its generated contract remains compatible.
New targets have empty standard-name lists and exact parameters/generators in
the registry; the label does not assign them a named curve standard.

The final whole-search replay passed for all 22 maps and checked stdout/stderr
digests for all 40 degree outcomes (`validation/replay-final.log`). P-192 has
three completed degrees and 21 timeouts; corrected P-224 has eight completed
degrees and eight timeouts. All six probes at or above 1009 timed out.

Before the canonical update, Conductor reported that T-145 held the registry
while repairing the root library gate. The report was therefore rehearsed in
an isolated temporary copy. Its first registration attempt exposed a parser
assumption that `curves` was the final top-level registry member; the current
registry also has a trailing `standards_source`. The failed generation log is
retained (`validation/report-first-generation-failure.log`). The corrected
stream parser replaces only the curve array and preserves the complete
top-level suffix. The canonical registry was not edited during this failure.

T-145's registry reservation cleared before the actual target registration.
The construction task T-65 had completed its search and exhausted its four
lease-claim attempts; catalogue closeout was recorded separately as T-148,
with the study and all eight canonical paths reserved. No additional research
agent or construction campaign was dispatched.

The final pre-catalogue root release-library check again failed with the same
643 compiler errors (`validation/lib-test-final.log`). The report's three-page
PDF was rendered with Poppler at 130 dpi and visually checked in an isolated
rehearsal; all identifiers, arrows, coverage labels and evidence paths were
legible with no clipping. Final source and canonical-output hashes are frozen
after the actual catalogue refresh.

The first native catalogue pass verified all 343 models, with zero invalid or
unsupported inputs. Diff review showed that globally sorting the expanded
registry moved existing standard records. Registration now preserves all
existing raw entries and their order and appends only the 22 replayed models.
The derived views were regenerated from this preserved-order registry.

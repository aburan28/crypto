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

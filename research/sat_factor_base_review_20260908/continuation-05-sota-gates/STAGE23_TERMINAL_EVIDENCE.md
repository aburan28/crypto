# Stage 23 compact terminal evidence

The Stage-23 terminal-evidence tools turn one already terminal unknown-scalar
run into a smaller directory that can be checked without the original absolute
paths. Packaging never launches a producer and refuses a run whose
`run-seal.json` is missing, nonterminal, self-inconsistent, or different from
the files on disk.

The compact archive preserves the original run seal and every sealed run,
task, input, and copied-binary byte outside `build-target/`. It also preserves
the complete outer-process directory and the complete project-verification
directory. The omitted set is derived from the original run-seal inventory and
contains every and only `build-target/` path. The two retained Stage-23
executables remain bound to their omitted
`build-target/release/examples/<name>` outputs by the original seal's byte
counts and SHA-256 values.

The retained run grammar permits only the three inputs, the run summary, the
two bound binaries after a successful build, and the exact files implied by
each reconstructed task receipt. The outer tree is exactly
`driver.metrics.json`, `driver.stdout`, and `driver.stderr`; project
verification is exactly `verification.json` and `verification-seal.json`.
Correctly re-sealing an extra file in any of these trees is still rejected. A
regular file named exactly `build-target` is invalid; only descendants of the
real build directory form the omitted partition.

Paths are interpreted through typed roots for the run, outer receipts, project
verification, and source repository. The verifier compares only schema-defined
path fields with the corresponding typed original path. It never rewrites
arbitrary strings. Archive paths must be normalized relative POSIX paths; `..`,
absolute paths, backslashes, duplicate paths, duplicate JSON keys, links,
non-regular files, and even correctly re-sealed extra files are rejected.
Reads and inventories walk path components through directory descriptors with
`O_NOFOLLOW`; each regular file is read and checked through one descriptor.
New files use `O_EXCL|O_NOFOLLOW` and receive their mode with `fchmod` before
that same descriptor closes.
File reads require mode, link count, device, inode, size, nanosecond mtime,
nanosecond ctime, and bytes-read length to remain identical across the read.
Every directory must retain the same metadata and entry-name snapshot across
all child traversal. Immediately before success, the verifier inventories the
entire bundle again and requires exact equality with its initial authenticated
snapshot, including the bundle seal.

## Trusted verification entrypoint

Obtain `verify_koblitz_stage23_terminal_evidence.py` and
`koblitz_stage23_terminal_evidence.py` independently from a reviewed repository
checkout. Pin that checkout by an expected Git commit and compare both tool
SHA-256 values with values received through a separate trusted channel before
opening a bundle. Run this independently obtained verifier as the first and
only verification entrypoint:

```text
python3 /trusted/reviewed-checkout/scripts/verify_koblitz_stage23_terminal_evidence.py \
  /absolute/path/to/terminal-evidence
```

The copies under `terminal-evidence/verifier/` are evidence bytes retained to
describe the packaging implementation. They are not a trust anchor, must never
be executed to authenticate their own bundle, and do not replace the
independently pinned verifier. The trusted verifier treats every bundled source
and verifier file as inert data.

The first independently reviewed publication of this trusted verifier is
commit `f977ab014a5ba00aa6084542e2b09cc766a20092`. Its wrapper SHA-256 is
`c3a2431c7d6c99d0e9efcf13b441d00e86937220bd9a883faa42c47f12026da9`,
and its core SHA-256 is
`307c4933fc18fdbd2385224a6eab78e6283509c765dece45fde6a07ca91310bb`.
The reviewed successor-03 bundle seal is
`d2b096d39b9158e8da844d89c58732b83dfb9a38ebe08df4f1229a88a581ca5c`;
its compressed publication SHA-256 is
`7886a767b0fe0f3372bd034429e3616632c40ea1781de835602fa08b72e9dd44`.

The typed external-tool map retains Python, Cargo, and rustc identities for
command reconstruction. Each mapped path, size, and SHA-256 must equal the
corresponding immutable entry in the archived source binding; the tool map
cannot substitute a coherently re-sealed executable identity.

## Packaging

Package a terminal run into a new directory:

```text
python3 scripts/package_koblitz_stage23_terminal_evidence.py \
  --run-root /absolute/path/to/terminal-run \
  --outer-root /absolute/path/to/outer-directory \
  --outer-metrics /absolute/path/to/outer-directory/driver.metrics.json \
  --project-verification-root /absolute/path/to/project-verification \
  --output /absolute/path/to/new-terminal-evidence
```

The packager checks the finished archive before returning success, but that
self-check is not independent authentication. The trusted verifier is read-only
and disables Python bytecode writes. It validates the internal consistency of
the bundle and original seals, derives the retained and omitted partitions
again, checks the binary/build cross-bindings, and reopens every raw intent,
metric, receipt, stdout, stderr, and parsed result. It reconstructs the exact
task prefix, commands, watchdogs, environment, inputs, rows, aggregate
resources, ratios, critical-failure state, and panel status. It then checks the
original outer receipt and project-verification record against that
reconstruction.

## Source boundary

The archive contains every regular file in the bound Git commit, the raw commit
object, every blob identity, and each Git mode. Verification recomputes every
blob ID, the recursive root-tree ID, the commit ID, and the exact
`git ls-tree -r --full-tree` SHA-256 recorded by the run.

The repository-root `Cargo.lock` is ignored and is not claimed as tracked. It
is archived separately and must be byte-for-byte equal to the tracked frozen
lock at
`research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-20-rust-build/Cargo.lock`.
This gives a complete project-source closure and a frozen dependency-version
binding. Registry source trees and the external Python/Cargo/rustc executable
bytes are not archived, so this is not a hermetic rebuild packet.

## Evidence boundary

This verifier performs project-authored custody, schema, transcript-hash,
modular relation-matrix, accounting, and derivability checks. It does not
independently replay all finite-field point generation, factor-base
materialization, decomposition witnesses, curve validity of relation rows,
elliptic-curve scalar checks, or the signed-Frobenius rho walk. Those checks are
listed explicitly as pending in every successful verification result.

Accordingly, `scientific_measurement_admitted`,
`external_portable_verification_satisfied`,
`independent_external_reproduction_satisfied`, `full_cost_gate_passed`, and
`koblitz_index_calculus_sota` remain false. A structurally verified archive is
custody evidence for the finite public-synthetic panel; it is not an external
reproduction, a cryptographic-size result, an asymptotic conclusion, or a SOTA
claim.

The machine-readable protocol is
`stage-23-terminal-evidence-protocol.json`. Adversarial controls are in
`scripts/test_koblitz_stage23_terminal_evidence.py`.

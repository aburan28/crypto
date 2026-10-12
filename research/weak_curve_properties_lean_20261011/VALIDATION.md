# Checked proof sources and retained validation

The Lean package builds, passes the all-declaration transitive axiom audit,
and passes replay of its compiled proofs in a fresh Lean kernel. It contains
48 explicit public theorem declarations. The audit covers 71 declarations,
including generated helpers, with only `propext`, `Classical.choice` and
`Quot.sound` as transitive axioms. No incomplete terms or custom research
axioms are admitted.

## Final local checks

Commands run in this directory unless a different directory is named.
The final build completed before the final audit and kernel replay were run.

| Command | Exit/status | Evidence |
| --- | --- | --- |
| `lake build` | 0; 1,638 build jobs complete | `build_final.stdout`, `build_final.stderr` |
| `lake env lean AxiomAudit.lean` | 0; `AXIOM_AUDIT_COMPLETE\|declarations=71` | `axioms_final.stdout`, `axioms_final.stderr` |
| `lake env leanchecker --fresh --verbose WeakCurveProperties` | 0; explicit root module replayed with `--fresh` | `kernel_final.stdout`, `kernel_final.stderr` |
| `rustc --edition=2021 proof_check.rs -o .lake/build/bin/proof_check` and `.lake/build/bin/proof_check --check` | Native source and receipt seal verified | `EVIDENCE_SHA256.tsv`; verification result reported by the command |
| `cargo test --release --lib` from the repository root | 101; 3,880 passed, 14 failed, 110 ignored | `root_lib.stdout`, `root_lib.stderr` |

The bundled `leanchecker` prints `replaying WeakCurveProperties with --fresh`
and returns zero after checking the imported and local declarations. Its
final stderr is empty. The explicit root module is necessary because the
Lake package name and root Lean module name differ.

The source/receipt checker hashes every file in this proof directory except
`.lake/` and the seal itself, plus the repository SHA-256 implementation,
CI workflow. The AGENTS instruction is reviewed as a separate policy addition,
so unrelated future policy changes do not invalidate the mathematical proof
seal. The checker also requires the exact frozen input hash
`9738f6bd6c016732f88e4e7979ac09e9aa6fd16eba96220ec8781c0017b6055e`.
The checker verifies source identity; the Lean build and kernel replay supply
the proof checks. CI repeats both rather than relying on archived stdout.

## Broader Rust baseline

The mandatory repository library suite was run in the same worktree. This
addition changes no pre-existing Rust library or example source. Five tests
failed assertions in `jv_cover` or `jv_quintic_edwards`; nine network tests
failed while binding local sockets with `PermissionDenied: Operation not
permitted`. Their complete names and assertion output are preserved in the
raw receipts. They remain failed checks, not passing validation, and are
separate from the successful proof-package checks. No timing ratio is derived
from the shared-host run.

## Earlier attempts and setup

Raw setup, cache and build attempts are retained with their original output.
Their whitespace is preserved and exempted from Git's whitespace checks;
proof sources and documentation retain normal whitespace validation.
The initial package update was interrupted before bulk cache acquisition;
the successful update used `MATHLIB_NO_CACHE_ON_UPDATE=1`. Cache acquisition
then used an external-volume cache directory. The first cache request used
file paths relative to the wrong root; the module-name request succeeded.
Early proof builds retain parser, type and normalization errors corrected in
the final sources. They are development attempts, not final pass receipts.

The files named `*_before_build_completion` record checks interrupted because
the aggregate build had not yet finished. `axioms_run2.stdout` includes audit
output but also a parser error and exited unsuccessfully; it is not the final
audit. Only the explicitly named final receipts above establish the pass.

## Formalization scope

The finite-field/subfield image criterion, actual short-model twist
isomorphism classifications, normalized generic-model descent, point-group
torsion factoring and order arithmetic are checked. The ordinary-curve
conductor proposition remains partial: the actual CM-order/Frobenius basis,
geometric rationality bridge and algebraic-endomorphism factorization are
the concrete remaining obligations in [PROOF_MAP.md](PROOF_MAP.md).

The 814 numerical records, independent PARI checks, and earlier census reports
keep their original verification status. This is a versioned proof addition;
it does not relabel the numerical rows as individual Lean certificates.

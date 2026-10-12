# Validation and remaining research work

The conjecture/methodology deliverable is implemented. Its eight-entry
ledger contains six mathematical conjectures, one statistical hypothesis
and one refuted finite screen. Exact domains and logical falsifiers are
documented; all new geometric witness searches remain pending.

## Completed checks

| Check | Result | Receipt |
| --- | --- | --- |
| Standalone native review of seven frozen censuses | Exit 0; 87,486 rows, 2,074 union zeros, 1,184 boundary-strip violations | `review_final.stdout`, `review_final.stderr`, `REVIEW.tsv`, `FALSIFIERS.tsv` |
| Recompile and replay comparison | All three output tables/hashes byte-identical to the first run | `REVIEW_run1.tsv`, `FALSIFIERS_run1.tsv`, `INPUTS_SHA256_run1.tsv` and final counterparts |
| Independent input hashing | System `shasum -a 256` agrees with the native SHA-256 implementation on all seven inputs | `independent-input-hashes.tsv`, `INPUTS_SHA256.tsv` |
| Native parser/depth controls | 3 passed; inconsistent union and unknown Boolean label rejected | `tests_final.stdout`, `tests_final.stderr` |
| `jq -e -f validate.jq conjectures.json` | Exit 0; unique IDs, 7 open/1 refuted, required fields and valid source references | Validator source and ledger |
| `lake build` | Exit 0; actual three-theorem module built | `build_run1.stdout`, `build_run1.stderr` |
| `lake env lean AxiomAudit.lean` | Exit 0; all three public theorems use only `propext`, `Quot.sound` | `axioms_final.stdout`, `axioms_final.stderr` |
| `lake env leanchecker --fresh --verbose CensusCounting` | Exit 0; fresh kernel replay completed | `kernel_final.stdout`, `kernel_final.stderr` |

The final source/receipt seal is `EVIDENCE_SHA256SUMS`; verify it with
`shasum -a 256 -c EVIDENCE_SHA256SUMS` locally or `sha256sum -c` in Linux CI.
The workflow then reruns the native review, compares outputs, checks the
ledger and rebuilds/replays the Lean proofs. Hash verification does not
prove a conjecture. The native review proves consistency/completeness of
the imported labelled table; correctness of geometric labels retains the
parent native/PARI evidence.

The existing Rust library and examples are unchanged. Applicable checks
target the new standalone review and proof package. This addition introduces
no arithmetic solver, selected subgroup benchmark or measured speedup.

## Retained development failures

The first standalone test compilation also compiled the repository SHA-256
module's tests, which require the Cargo `hex` dependency; `rustc --test`
without that dependency failed. `tests_compile_run1.stderr` retains the
error. The corrected standalone tests exercise the review parser/depth
functions; the production hashing path is independently compared against
system SHA-256 instead. The final test build's dead-field warning is retained.
Both subsequent test runs pass all three controls. A byte comparison of
their stdout differed because the parallel test harness printed tests in
different orders; both raw outputs are retained. Test assertions and exit
status establish the control result. Deterministic scientific tables and
input hashes are compared byte-for-byte separately.

The first axiom command failed because the standalone Std-only module did
not import `Lean.Elab.Command`. `axioms_run1.stdout/stderr` retain that error.
The final audit imports the required command module and succeeds. Only the
named final receipts establish its pass. Raw output whitespace is preserved
through local attributes.

## Research and novelty status

- Fixed-census exploration and the boundary-screen refutation are complete.
- The counting theorems are Lean-checked. They verify aggregation semantics,
  not geometric membership or any open conjecture.
- Full Type I/II support outside the current union, primitive two-tower
  witnesses, smaller-genus primitive quotients, the general torus-degree
  theorem and polarized target profiles remain open as mapped in
  `PROOF_MAP.md`.
- The primary-source audit found a published-class backlog. Detailed
  comparison with every classification table and later citation chain is
  still required before a novelty claim.
- The availability audit found old p53/p59/p61/p67 cubic artifacts without
  reading their outcome labels. Proposed prospective primes p71/p73/p79
  had no matching filenames in the inspected ISO-1 archive; execution
  registry confirmation and freezing remain required.
- New 192–252-bit structural panels and construction-cost experiments are
  planned, not launched or completed here. Existing 512-source class labels
  retain their previous unresolved/bounded status.

No parent census, verified family record, 814-model property record or prior
proof source is modified. The catalogue README links this separate
conjecture ledger, preserving the existing family count and labels.

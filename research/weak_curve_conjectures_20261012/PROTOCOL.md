# Conjectures and methodology: requested scope and first exploration

Requested: assemble conjectures about additional weak-curve families and
distinctive properties, investigate prior analysis, and build a methodology
for discovering and testing them. Preserve the 192–252-bit field agenda and
the standing requirement to attempt faithful Lean proofs when feasible.

| Requirement | Deliverable |
| --- | --- |
| Candidate families and properties | Individually stated conjectures, domains, exclusions and falsifiers |
| Prior analysis and novelty | Primary-source audit; published families retain attribution; exact new statements have bounded novelty status |
| Discovery methodology | Model/class/cover/property schema, frozen pilots, independent replay and promotion criteria |
| Exploration now | Native review of seven fixed, complete, previously verified censuses |
| Larger fields | Gated structural-validation protocol for 192–252-bit fields; new construction and performance runs remain pending |
| Lean | Checked universal counting lemmas and an explicit formalization map for open geometric conjectures |

Before the fixed-census review, register two exploratory tests:

1. Count branch-only, shared and union-zero classes with one row per signed
   ordinary trace `t=2 mod 4`. Verify the CSV union labels, the complete trace
   domain and the existing cubic conductor-depth restriction independently.
2. Test the simple boundary-strip statement: every union-zero row in the
   existing seven censuses satisfies `2*p^3-|t| <= 64*p`. Preserve each
   violating class if this statement fails. This finite statement does not
   constitute an asymptotic density conjecture.

Inputs are ONLY `p{7,11,13,17,19,29,37}_classes.csv` in the frozen
`iso1_weak_classes_20261007/two_branch_20261009` study. The native program
accepts no externally supplied curves, public keys or input path. It checks
already certified labels and computes descriptive arithmetic; it generates
no new curve-family support certificates. Save source/input hashes, raw
stdout/stderr and exact commands. Preserve zeros and exception cases.

Use a standalone Rust program and the installed pinned Lean 4.34.1 toolchain.
The Lean lemmas check the universal two-predicate counting identities, not
the truth of imported curve-family labels or any open conjecture. A separate
proof map identifies the actual geometric obligations. No paid compute or
new large-field arithmetic campaign is part of this initial exploration.

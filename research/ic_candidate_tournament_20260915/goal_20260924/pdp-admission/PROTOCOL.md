# PDP outcome admission controls

This is a correctness/accounting change before admitting algebraic challengers
to the remaining two bounded improvement rounds. It does not consume an
improvement attempt or change the sealed round-one source, protocol or results.
Round-one evidence merged in [PR 782](https://github.com/aburan28/crypto/pull/782)
at `b338522537f316f2379ed8cc11bb317212e573ec`; no challenger qualified.

## Defect and hypothesis

The narrow Gröbner and Crossbred frontends returned default statistics on an
unencodable input. The relation-attempt reporter classified an empty result
without `exhausted` as `Refuted`. Narrow SAT already classified the result as
unknown but omitted it from its unknown counter. These are measurement errors:
infinity, a one-summand request, or a layout wider than 64 variables says
nothing about whether a decomposition exists. The uncached layout arithmetic
also lacked the template builder's checked overflow handling.

Hypothesis: explicit unsupported-input metadata, propagated to the attempt
ledger and accompanied by the legacy incomplete flag, prevents false negative
claims while preserving completed solver results and actual budget exhaustion.
The status distinguishes an unsupported encoding from an exhausted search.

## Frozen controls and stop conditions

Use the checked-in `ci/Cargo.lock`, one test thread, release builds and the
existing Linux integration workflow (Rust 1.94.1). Local macOS results are
correctness controls only. All inputs below are synthetic toy curves.

- On `K_1/F_(2^9)`, the existing factor-base index zero supplies known witnesses
  for infinity (`P + -P`), one summand (`P`), and a 16-summand target (`[16]P`).
  Test F4, F5, inherited F4, Crossbred, native-XOR SAT and CNF SAT. Require no
  answer, explicit `unsupported`, incomplete status, zero solver calls/reductions
  and no refutation. Check layout overflow through both builders.
- On `K_0/F_(2^9)`, give the three Gröbner engines a zero-node budget on a
  supported two-summand input. Require incomplete without unsupported/refutation.
- Run eight relation attempts per engine on `K_1/F_(2^9)` with 16 summands,
  default frozen seed and direct recovery disabled. Require every attempted
  nonidentity query to remain `Unsupported`, never `Refuted`, zero relations,
  and SAT unknown counts consistent with the retained attempt records.
- Extend the existing 39-target enumeration/Gröbner/SAT comparison to F5,
  checking every returned witness in the group. Preserve genuine UNSAT
  controls, three-summand chain controls, Crossbred tests and template tests.
  Any incorrect witness, false refutation, crash or lost receipt fails admission.

The regression commands and observed outcomes are retained in `RESULTS.md`.
These controls contain no performance selection, new held-out confirmation or
natural-yield estimate. Stop on a failing correctness check and repair it before
the next gate; do not turn unsupported inputs into favorable timing samples.

## Accounting and remaining gates

Classification: **accounting/correctness**. Complete IC online time, cold time,
instructions, S, boundary ratios and speedup are **unknown** for this change.
No benchmark speedup is claimed. Planted witnesses are correctness controls.

This does not yet qualify the generic engines for the full tournament. Actual
dispatch identities, public-point inputs, exclusive stage costs, ordinary-query
status/yield, subgroup base census, matrix/descent certificates and pinned
resource accounting remain required. Macaulay GF(2) work stays within PDP;
relation LA uses the subgroup scalar modulus. Round two must retain the qualified
optimized incumbent, strong matched rho, all prior point exclusions and the
predeclared confirmation rule. Do not tune on round-one confirmation or replay.

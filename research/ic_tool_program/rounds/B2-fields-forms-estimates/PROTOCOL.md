# B2: prime and extension fields, the other curve forms, the importers, and estimates

**Declared 2026-10-01, before any B2 code.** This is Track B's third step
in `research/notes/index-calculus/IC_TOOL_PROGRAM.md` §9 ("validation"),
built to the design in
[`../../design/schema-v2.md`](../../design/schema-v2.md), §10's B2 row.
B2's cases, C032–C051, are frozen with this declaration in
[`../../conformance/v2-b2/`](../../conformance/v2-b2/cases.json), as the
design's §9 requires. Nothing below changes after the first measurement,
except by a dated amendment at the end.

## What B2 does

- **Prime fields** (design §3.1, §4.2):
  - `p ≥ 5` up to 1,024 bits;
  - `p-range` and `p-composite`, under the design's primality test
    (exact below `ψ₁₃`, a disclosed screen above it).
- **Extension fields `GF(p^k)`.**
  - Parsed and validated: `degree-range`; `modulus-reducible`, by
    Rabin's test over `GF(p)`; and the curve and subgroup checks in
    `GF(p^k)` arithmetic.
  - Routed only from B5. Until then the code is `no-pipeline-for-field`.
- **The other curve forms** (design §3.2):
  - `general_weierstrass` over a binary or a prime field, `montgomery`
    and `twisted_edwards`.
  - Each is converted, with its points, to `binary_weierstrass` or
    `short_weierstrass` by the standard isomorphisms and birational
    equivalences.
  - A binary conversion with `n` odd also takes `a` to 0 or 1, so that a
    Koblitz curve in disguise reaches `kic` (C037).
  - A binary `general_weierstrass` curve with `a₁ = 0` is refused as
    `supersingular-binary`.
- **The group order** (design §4.5): method 3, the unique multiple of
  `r` in the Hasse interval when `r > 4√q`, and method 4, enumeration
  for `q ≤ 2^32`. Neither applying is `cardinality-unknown`.
- **The disclosures for prime fields**: `embedding-degree` and
  `anomalous`.
- **The importers**:
  - `rho-negation`: `ic_boundary`'s negation-map rho, binary `n ≤ 62`,
    prime `p < 2^63`.
  - `ic-binary-s4`: `ic boundary`'s generic binary pipeline,
    `n` in 5..=32.
  - `ic-prime-s3`: `ic boundary`'s prime pipeline, `p < 2^63`.
  - `kic` on curves over a subfield with `k > 1`, paired with
    `rho-negation`.
  - **`rho-bignum`**: new. It is rho with the negation map on
    big-integer arithmetic, and it takes any valid instance whose
    estimate fits the budget, whatever the field's width.
- **Named curves** (design §3.7): `named` stands for a published curve's
  field, curve and subgroup, from the curve registry and
  `params::named`. A name without a recorded generator is refused as
  `named-incomplete`.
- **Estimates and budgets** (design §5.4, §3.6):
  - Every route carries an F2 estimate for each arm.
  - A run whose estimate exceeds its budget is refused before it starts,
    with exit status 4.
  - `fidelity: F2` returns the estimates and runs nothing.
- **`recipe: auto`** (design §5.3), and **`solve: index_calculus`**
  alone.

B2 leaves two refusals in place: F1 (B7), and a budget of more than one
thread, which A6 takes up.

### What the cases read that the design left open

C032–C051 fix their expectations as report keys. The design states
these keys; this declaration makes them exact, and its amendment to the
design records them:
- **`conversion`** (§6): `{from, to, map}` for a converted curve. `from`
  is the input's form; `to` is `short_weierstrass` or
  `binary_weierstrass`; `map` states the substitution.
- **`estimate.<arm>`** (§5.4):
  - `{pipeline, seconds, model}`, and for rho also `steps` and
    `steps_bits`;
  - `steps_bits` is `⌈log₂ √(πr/2A)⌉`, with `A = 2n` for signed
    Frobenius and `A = 2` for the negation map;
  - `seconds` uses the pipeline's step cost on the reference host, and
    `model` names the source.
- **`status: estimated`** for F2.
- **The disclosure `study-pipeline`** (§4.6), whose value says why. It
  is set when the route uses `ic-prime-s3`, since no subexponential
  index calculus is known for prime fields.

## Class

**Robustness.** B2 widens what the tool accepts and runs. It changes
nothing on any input B1 already ran, and it claims no speedup.

## Arms

- **The base**: the newest accepted baseline when B2 runs, recorded in
  B2's manifest. B1 comes before B2. B3, if it is accepted first, comes
  before B2 too.
- **B2**: that commit plus B2's change, built with
  `IC_BUILD_COMMIT=$(git rev-parse HEAD)`.

## Measurements

1. **Tests.** `cargo test --release --bin ic`, and the library tests of
   every module B2 touches.
   - `rho-bignum` is tested against `rho-negation` on instances both
     admit: the same answer on every target tried.
2. **Conformance.**
   `conformance/run.py --steps <the steps accepted>,B2`. For example,
   `B0,B1,B3,B2` if B3 is accepted first. That runs C001–C051, the
   later steps' cases, and the `until` rule. The same runner also runs
   on the base, and its failures are recorded as what the base could
   not do.
3. **The pin, untimed.** All 90 suite rows from their v1 files. Every
   output must equal the base's.
4. **The translation, untimed, on all 90 rows.** This is B1's check
   repeated: `ic check --translate`, then `ic price` on each
   translation, with outputs identical to the v1 file's.
5. **No slowdown, timed.** The base against B2 on `M1`'s 22 rows: five
   rounds ABAB, isolated, 220 processes. The figure is the paired
   cold-time ratio per size.
6. **The estimate's error.**
   - B2's F2 estimate of the cold cost at the eleven suite sizes,
     against v0's measured cold times (R01).
   - Reported as a table of ratios with their geometric mean and range.
     This is a model's error, recorded, not gated.
7. **`rho-bignum`'s step cost.** One-thread steps at three field widths
   (binary `n = 100` and `127`, prime `p ≈ 2^127`), isolated. These are
   the step costs its estimates use from then on.

## Acceptance

B2 is **accepted** when all of the following hold:
- every test passes;
- every case the runner selects passes;
- the pin holds on every row;
- both parts of the translation check hold on every row;
- no size regresses beyond its A/A band, judged as in B1;
- the estimate's error and `rho-bignum`'s step costs are reported.

It is **stopped** if any v1 output differs. It is **rejected** if a test
or a case fails, the translation check fails, or a size regresses. A
rejected B2 is fixed and declared again, by amendment.

## Inadmissible

- Loosening a case after a run. The only exception is the `until` rule,
  applied by a later step.
- Changing a file under `conformance/v2-b2/` after this declaration.
  `SHA256SUMS` pins B2's files.
- Changing `conformance/run.py`'s rules other than by a step's
  declaration.
- Tuning an estimate on the suite rows its error is measured on, after
  measuring it.
- Counting a contended run.

## Cost

- The tests and the conformance suite: about ten minutes.
- The pin and the translation: 270 untimed processes.
- The timing check: 220 processes, about 90 minutes.
- The estimate's error: no processes, since it reads R01's figures.
- `rho-bignum`'s step costs: three short runs.

## Amendment 1 (2026-10-01, before any B2 measurement)

Written once B2's code passed every selected case (C001–C058 under
`--steps B0,B1,B3,B2`, 57 of 57), and before the first measurement. It
records where the code differs from the text above or makes it exact.
None of it loosens a case: no case file changes.

1. **A step of its own: B2b.** Three items move out of B2:
   - `kic` on a curve over a subfield with `k > 1`;
   - `solve: index_calculus` on `kic`;
   - order certificates.

   `ic price` prices Koblitz curves only. Pairing `kic` with
   `rho-negation`, or running `kic` alone, needs a second rho arm and a
   new report shape in the single-target pricer, whose outputs the pin
   holds fixed. Until B2b:
   - `kic`'s gate for `k > 1` stays `subfield-curve-unsupported`;
   - `kic` alone is `not-yet-supported`, naming B2b;
   - a certificate is recorded as `not_checked`.

   The imported index calculus pipelines, `ic-binary-s4` and
   `ic-prime-s3`, run alone from B2.
2. **How method 4 counts.**
   - Up to `q = 2^16` it enumerates.
   - Above that it uses Mestre's method: the orders of points on the
     curve and on its quadratic twist, intersected in the Hasse
     interval. That is exact as well.
   - If 64 points leave more than one value, the result is
     `cardinality-unknown`.

   Enumeration alone would take minutes near `2^32`. The tests check the
   counts against enumeration over `F_p`, `GF(17^3)` and `GF(2^13)`.
3. **A recipe object belongs to `kic`.** A recipe object holds v1's
   knobs, and only `kic` reads them. When a document gives one, every
   other index calculus pipeline is gated `recipe-not-taken`, a new gate
   code. So `auto` never routes the arm to a pipeline that would ignore
   the user's recipe.

   B1's frozen case C030 is the instance this protects. It is a Koblitz
   curve over `GF(2^32)` with an explicit `kic` recipe; without this
   rule, `ic-binary-s4` (`n ≤ 32`) would take it.
4. **The other field's pipelines.** A pipeline with no implementation
   for the instance's field kind is gated `no-pipeline-for-field`:
   - `kic`, `rho-koblitz` and `ic-binary-s4` on a prime field;
   - `ic-prime-s3` on a binary one;
   - every pipeline on `GF(p^k)`. There the refusal is
     `no-pipeline-for-field` itself, before any arm is chosen, as C050
     reads.
5. **The gates F2 skips** are the five word-width codes:
   - `field-wider-than-one-word` and `field-wider-than-two-words`;
   - `prime-wider-than-one-word`;
   - `scalar-wider-than-63-bits` and `scalar-wider-than-127-bits`.

   `enumeration-bound` is not one of them. `ic-binary-s4`'s enumeration
   of `GF(2^n)` is the pipeline's method, not a word.
6. **The budget compares the arms' sum.** A paired run's arms share one
   process and one wall budget, so the router refuses a run whose arms'
   estimates together exceed it. That implies the per-arm test above.
   Under `max_iterations`, rho's estimate counts that many steps.
7. **The estimates' models.** Every constant beyond `baselines.json` is
   in [`estimates.json`](estimates.json), with its source.
   - **`kic`.** §20's `ic_phases`, evaluated for one target (`K = 1`,
     the single-target rule's count, where `predict.py` used 32).
     - `K = 1` was chosen before any comparison with v0.
     - A development check against v0's cold times at the eleven sizes
       was then computed: geometric mean 0.91, range 0.34–1.98.
     - Nothing in the model changed after it. Measurement 6 reports the
       binary's own figures.
   - **`recipe: auto`.**
     - At the suite's sizes it is the suite's row. There are twelve with
       the smoke size, keyed by `(a, n, r)`.
     - Elsewhere it is `make_params.py` at §20's optimum for one target.
     - It rounds as Python's `round` does. A test reproduces every suite
       file's recipe.
   - **`rho-koblitz`.** The newest baseline's row nearest in `log₂ r`:
     `s_rho_online × unit_ns × √(4n/π)` nanoseconds a step.
   - **`rho-negation`.** No baseline row times it, so it uses the
     boundary ledger's matched rho at the largest rung of its regime.
   - **`ic-prime-s3` and `ic-binary-s4`.** Each pipeline's measured cost
     at the ledger's largest rung, carried along its fitted exponent.
   - **`rho-bignum`.** Provisional development-host step costs, until
     measurement 7 replaces them in `estimates.json`.
8. **What the imported pipelines run.**
   - **`ic-prime-s3`** is the ledger's `mitm_m2_negfold_walk_balanced`
     row: two summands (`S₃`) by meet in the middle, on the base the
     family's shape law asks for.
   - **`ic-binary-s4`** is `mitm_m3_negfold_walk`: three summands
     (`S₄`) on the `⌈n/3⌉` subspace base. Where a census finds three
     cannot reach the subgroup, it uses two, and says so.
   - Each times its set-up apart from its target-dependent work. Neither
     splits that work into AGENTS.md's five exclusive phases, so their
     reports mark the pair `speedup_eligible: false`.
9. **Targets v1's rules do not cover.** On an instance outside v1's
   rules (not a Koblitz curve on `kic` or `rho-koblitz`),
   `random_seed` and `public_hash_seed` follow v2 rules. Each is SHA-256
   in counter mode under a fixed label, and the report states it.
10. **Characteristic 3.**
    - A `prime_extension` with `p = 3` is refused with
      `no-pipeline-for-field` once its field is validated (design §11).
      B2 validates no curve over it.
    - `p < 3` is `p-range`.
11. **Named curves.**
    - `named` takes the names `params::NAMES` lists.
    - A name with no recorded generator or polynomial basis is
      `named-incomplete`, class `unsupported`, exit 3. Today that is
      `ecc2k-130`.
12. **Measurement 7** runs `examples/rho_bignum_rate.rs`:
    - stand-in curves of the declared widths;
    - a fixed step count, isolated, on one thread.

## Amendment 2 (2026-10-01, before any measurement)

**Measurement 5, the timing check, runs in the chain of Track B steps.**
- **What it was:** the base against B2 on `M1`'s 22 rows, five rounds
  ABAB, 220 processes, for this step alone.
- **What it is now:** one interleave over the newest accepted baseline
  and every Track B step's arm, in the queue's order: the baseline, B0,
  B1, B3, B2, B2b, B7a and B3b (`bround.py chain`). It runs on `M1`'s 22
  rows, five rounds, isolated, and the order reverses every other
  round.
- **B2's figure** is the paired cold-time ratio of the arm before it
  in the chain over B2's own, per size. It is read against R01's A/A
  bands, as before.
- **Why the chain is a valid base.** Each step's arm is built on the one
  before, so the arm before B2 in the chain is its base. The two run
  back to back in every round, with ten pairs a size, as in its own
  ABAB.
- **What it saves.** The seven steps' checks take 880 processes
  together, against 1,540 one by one.
- **If an earlier step is rejected,** the steps after it wait. The chain
  runs again from that step, once it is fixed or removed, since the
  later arms carry its change.
- **Nothing else changes:** the acceptance rule, the A/A bands and the
  other measurements.

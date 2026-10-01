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

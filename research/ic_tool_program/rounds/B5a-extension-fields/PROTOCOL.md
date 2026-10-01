# B5a: extension fields `GF(p^k)`, at F0

**Declared 2026-10-01, before any B5a code or measurement.** This is
Track B in `research/notes/index-calculus/IC_TOOL_PROGRAM.md` (§9, B5).
- **Its design:** [`../../design/extension-fields.md`](../../design/extension-fields.md).
- **Its cases:** [`../../conformance/v2-b5a/`](../../conformance/v2-b5a/cases.json).

Nothing below changes after the first measurement, except by a dated
amendment appended at the end.

## What B5a adds

**Every valid instance over `GF(p^k)` gets a route**, where today every one
is refused as `no-pipeline-for-field` (design §2):
- **`rho-negation` on `GF(p^k)`, for `q ≤ 2^62` and `r < 2^63`:** the
  matched rho's walk over the document's own field.
- **`rho-bignum` on `GF(p^k)`:** any `p` and any `k`.
- **`ic-gaudry-cubic`:** Gaudry's index calculus on `E(GF(p³))`, imported
  from the residual-walk thread's `gaudry_cubic`. It admits:
  - the modulus `t³ − c`;
  - `h = 1`;
  - `r < 2^63`.

  Elsewhere it is refused with a gate code that names the condition.
- **ICV1's extension part:**
  - the specification and the reference implementation, in this
    declaration;
  - the Rust port, with the code. The report's `curve_id` names every
    extension instance.

**B5a changes no binary or prime path.** The suite's rows, and every
earlier step's instances, run exactly what they ran.

**The plan's B5 is split.** B5a is extension fields. Prime fields past one
word are B5b, declared on their own later. B5's row in the plan is done
when both are accepted.

## Class

**Robustness.** B5a adds instances the tool can solve. It claims no
speedup. Its figures are new rows, not gains.

## Arms

- **The base:** the newest accepted baseline when B5a runs, with Track
  B's stack through B4, as a commit recorded in the manifest.
- **B5a:** that commit plus B5a's change.

## Measurements

1. **Tests.** `cargo test --release --bin ic`, and the library tests of
   every module B5a touches. They must include design §5's tests:
   - the field against `Fpk`;
   - the group law and the keys;
   - both rho pipelines on small instances;
   - sameness with `gaudry_cubic` on its own instances;
   - `curve_id::extension` against the reference's vectors.
2. **Conformance.** `conformance/run.py --steps <accepted>,B5a`, with B5a's
   cases (below). The same runner also runs on the base, and the base's
   failures are recorded as what the base cannot do.
3. **The pin, untimed.** All 90 suite rows from their v1 files. Every
   output must equal the base's.
4. **No slowdown, timed.** By the chain's rule (`bround.py chain`), with
   B5a the next arm after B4.
5. **F0 on extension fields.** Each run is one process, isolated, on one
   core, with an hour (`timeout 3600`). Each instance below (design §4)
   is run on two targets: its known-answer target and its public point
   `T001`. `H1` is C108's refusal, and is not run here.
   - `G1`–`G3`: `ic price`, paired, `ic-gaudry-cubic` against
     `rho-negation`.
   - `E2`, `E5` and `E11`: `solve: rho` with `rho-negation`.
   - C050's document, on its own known-answer target only:
     `rho-negation`.
   - `B2`: `solve: rho` with `rho-bignum`.

   Recorded for each run:
   - the logarithm, checked in the run and replayed outside it in the
     generator's arithmetic;
   - each arm's `S`;
   - every phase's cost and counts;
   - the host manifest and the isolation record.

   [`run.py`](run.py) runs them (`manifest`, then `f0`) through the
   programme's runner. [`analyse.py`](analyse.py) reads them, replaying
   every arm's certificate in the generator's arithmetic.
6. **Calibration of the estimates.** Five isolated runs each:
   - `rho-negation` on `E2`, `E5` and `E11`;
   - `rho-bignum` on `B2`;
   - `ic-gaudry-cubic` on `G1`–`G3`.

   Their medians replace the provisional constants in
   [`estimates.json`](estimates.json): the step costs by width, and
   `ic-gaudry-cubic`'s anchor and exponent, the exponent refitted over
   `G1`–`G3`. No case's expectation depends on these constants.

## B5a's cases

`conformance/v2-b5a/` is written by its `make_cases.py`. Its arithmetic
for `GF(p^k)` is its own, and shares nothing with the tool:
- **Group orders:** found exactly, as B2's generator does.
- **Instances:** read from `instances.json`, every property re-checked.
- **ICV1 slugs:** the reference implementation's.

| case | what it checks |
|:--|:--|
| C103 | C050's successor (it `supersedes` C050, whose `until` is B5): C050's document, `GF(1009³)` with a general modulus and `h = 1524`, recovers its known logarithm with `rho-negation`, verified |
| C104–C106 | `ic-gaudry-cubic` paired with `rho-negation` at F0 on `G1`–`G3`. Both arms recover the known logarithm and are verified, and `curve_id.slug` is the registry's |
| C107 | C050's document under `paired`: `no-ic-route`, with `ic-gaudry-cubic` refused as `modulus-not-binomial` |
| C108 | `H1` under `paired`: `no-ic-route`, with `cofactor-not-one` |
| C109 | `E2` under `paired`: `no-ic-route`, with `extension-degree-not-three` |
| C110 | `E2` under `solve: rho`: `rho-negation` at `q ≈ 2^62`, verified |
| C111 | `E5` under `solve: rho`: verified |
| C112 | `E11` under `solve: rho`: verified, with 11 coefficients |
| C113 | `B2` under `solve: rho`: `rho-negation` refused as `field-wider-than-one-word` and `rho-bignum` admitted, verified |
| C114 | `G1`'s curve in `general_weierstrass` form: converted, the conversion recorded, and the logarithm verified by rho |
| C115 | `G1` under `check`: exit 0, with `ic-gaudry-cubic` and `rho-negation` admitted |
| C116 | `G1` at `fidelity: F2`: `status: estimated`, with an estimate for each arm |
| C117 | `ic-gaudry-cubic` named on a binary document: refused as `no-pipeline-for-field` |
| C118 | `G2` with `recipe: {"oracle": "groebner"}`: the logarithm recovered and verified with Gaudry's `S₄` solve |

C050 names B5 as its `until`, and B5 is now two steps, so C103 retires it
by the `supersedes` rule. The runner's list of steps gains `B5a` and
`B5b`; this declaration makes that change to `conformance/run.py`.

## Acceptance

B5a is **accepted** when:
- every test passes, sameness included;
- every case the steps select passes;
- the pin holds on every row;
- no size regresses beyond its A/A band;
- every run in measurement 5 recovers its logarithm, verified in the run
  and by the replay.

It is **stopped** if a binary or prime output differs. It is **rejected**
if a test or a case fails, if a size regresses, or if a run in
measurement 5 gives a wrong answer, or none within its hour.

**What an accepted B5a may claim:** the tool solves the discrete
logarithm on curves over `GF(p^k)` at F0 where the subgroup allows:
- with rho, for any `p` and any `k`;
- with Gaudry's index calculus, on `E(GF(p³))` in the module's basis and
  of prime order;

with every answer verified. It claims nothing about the index calculus's
cost against rho beyond the rows it measures, which the residual-walk
thread priced (design §2.4).

## What B5a does not do

- Prime fields past one word (B5b).
- A basis change for other cubic moduli, `p ≡ 2 (mod 3)`, or cofactors
  in the index calculus (design §6).
- `k = 4` (`gaudry_quartic.rs`), or any other degree, for the index
  calculus.
- Characteristic 3.

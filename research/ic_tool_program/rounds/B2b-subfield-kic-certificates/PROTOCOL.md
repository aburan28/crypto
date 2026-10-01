# B2b: `kic` on curves over a subfield, `kic` alone, and order certificates

**Declared 2026-10-01, before any B2b code.** B2's amendment 1
([`../B2-fields-forms-estimates/PROTOCOL.md`](../B2-fields-forms-estimates/PROTOCOL.md))
moved three items out of B2 into this step:
- `kic` on curves defined over a subfield `GF(2^k)` with `k > 1`,
  paired with `rho-negation`;
- `kic` alone;
- the verification of order certificates.

The design is [`../../design/schema-v2.md`](../../design/schema-v2.md),
§3.3, §4.2, §4.4, §5.1–§5.3 and §10's B2b bullet. B2b's cases,
C059–C070, are frozen with this declaration in
[`../../conformance/v2-b2b/`](../../conformance/v2-b2b/cases.json).
Nothing below changes after the first measurement, except by a dated
amendment at the end.

## What B2b does

### 1. `kic` on a curve over a subfield

The instance is a binary curve `y² + xy = x³ + ax² + b` over `GF(2^n)`.
Its subfield degree `k` is the least `d` with `a, b ∈ GF(2^d)`, as the
`curve-subfield` disclosure already reports. Let `e = n/k`.

- **The gates.** `kic` admits `1 < k ≤ 8` when `e` is odd and at least
  3. Its gates apply in this order, and the first that applies is
  reported:
  1. `subfield-too-large` (`k > 8`), or `not-a-subfield-curve`
     (`k = n`);
  2. `field-wider-than-one-word` (`n > 62`);
  3. `even-extension-degree` (`e` even, or `e < 3`);
  4. `scalar-wider-than-63-bits`;
  5. `subgroup-smaller-than-cofactor` (`r ≤ h`).

  For `k = 1`, `e = n`, so the gates are B1's. `rho-koblitz` keeps
  `subfield-curve-unsupported` for every `k > 1`.
- **The importer** takes the curve as given: the modulus, `a`, `b`, `r`,
  `h` and `G`. The `q`-power Frobenius, `q = 2^k`, acts on it, and its
  trace comes from `E(GF(2^k))` counted by enumeration. The importer
  finds `λ` on `⟨G⟩` and requires what `KoblitzCurve::import` requires:
  - `h·r = #E`;
  - `r > h` and `r ∤ h`;
  - `G ≠ O` on the curve, with `[r]G = O`.
- **The pair.** The rho arm is `rho-negation`. The report:
  - carries the disclosure `rho-frobenius-unused`, valued
    `{"subfield_degree": k, "extension_degree": e}`;
  - sets `speedup_eligible: false`, with `speedup_ineligible_reason`
    saying the same thing. The reference leaves the `q`-power Frobenius
    unused (plan §8, A5). A walk on signed `q`-Frobenius classes would
    take about `√e` times fewer steps, and the tool has none for
    `k > 1`.
- **The recipe.**
  - `recipe: auto` is gated `no-recipe` (design §5.3).
  - A recipe the document gives is taken as it is.
  - A factor base that is not a `spec` stays outside the single-target
    pricer, as it is for `k = 1`.
- **The estimate** (F2, and the budget) is §20's phase model with `e` in
  the place of `n`: a signed orbit has `2e` points, not `2n`. Its
  `model` says the model was fitted for `k = 1` only. `rho-negation`'s
  estimate is unchanged.
- **The intervals.**
  - The index calculus's online interval and its five phases are
    `k = 1`'s.
  - The rho arm's interval is the library call
    `rho_reference_negation`. That call builds its jump table inside
    itself, so the interval includes the table. The table's `setup_ops`
    are reported apart, in operations.
  - The replay `[d]G = Q`, in the general arithmetic, follows each arm
    and is timed apart.

### 2. `kic` alone

`solve: index_calculus` on `kic` runs the paired price's index calculus
arm alone, for `k = 1` and for `k > 1`. It has the same reusable set-up,
online interval, phases, counts and certificate.

The report has:
- `operation: ic_single_target` and `ic.pipeline: kic`;
- `result.scalar`, `result.verified` and `result.known_answer`;
- `counts` and `certificates.ic`;
- `speedup_eligible: false`, with the reason "no rho arm".

The rho arm's keys are absent. The imported index calculus pipelines
run alone report `operation: ic_single_target` too. Under B2 they said
`price_single_target`.

### 3. Order certificates

- **Where.**
  - `subgroup.order_certificate` certifies `r`: the given order, or the
    one derived when `order` is omitted.
  - `field.p_certificate` certifies `p`, on fields of kind `prime` and
    `prime_extension`.
  - On a binary field, `p_certificate` is `unknown-key`.
- **The form.**
  - A certificate is `{"kind": "pocklington", "factors": [...]}`, with
    no other keys.
  - Each factor is `{"prime", "exponent", "witness"}`, and optionally
    `"certificate"`: a certificate for that prime, in the same form.
  - `prime` and `witness` follow the schema's integer syntax. `exponent`
    is a JSON integer, 1–1,024.
- **The check**, for a number `N ≥ 3`:
  1. The primes are distinct.
  2. With `F = ∏ prime^exponent`, `F` divides `N − 1` and
     `(F + 1)² > N`.
  3. For each factor `q` with witness `a`:
     - `1 < a < N`;
     - `a^{N−1} ≡ 1 (mod N)`;
     - `gcd(a^{(N−1)/q} − 1, N) = 1`.
  4. Each `q` is prime: below `ψ₁₃` by the exact Miller–Rabin test,
     otherwise by its own certificate, checked the same way.

  Then every prime factor of `N` is `1 (mod F)`, so it exceeds `√N`, and
  `N` is prime (Pocklington; Brillhart, Lehmer and Selfridge, 1975).
- **Limits**: 64 levels, and 1,024 factor entries in all.
- **Failure.** Anything that fails refuses the input with
  `certificate-invalid`, class `invalid`, exit 2. That covers the form,
  the syntax, a limit and an arithmetic condition. The message names the
  condition and the factor it failed at.
- **Success.**
  - The check `certificate-invalid` passes, exact, and names the number
    certified.
  - The primality check it replaces (`order-composite` or `p-composite`)
    passes exact, naming the certificate as its method.
  - That number gets no `primality-screen` disclosure.

  The certificate is checked before that primality check.

### 4. The runner

The step B2b is added to `../../conformance/run.py`'s `STEPS`, after
B2. The runner's rules change only by a step's declaration, and this is
B2b's.

## The cases

C059–C070 are in
[`../../conformance/v2-b2b/cases.json`](../../conformance/v2-b2b/cases.json).
`make_cases.py` built them, and `SHA256SUMS` pins them. Their arithmetic
is the generator's own. Each curve's order is found twice: by the
subfield count with the trace recurrence, and by the orders of points.

| case | input | expected |
|:--|:--|:--|
| C059 | a curve over `GF(4)` taken over `GF(2^22)` (`e = 11`, `r = 2097349`, `h = 2`); `kic` with a recipe, paired with rho `auto` | exit 0, verified; `rho-negation`; `rho-frobenius-unused`; `speedup_eligible: false` |
| C060 | C059 with `recipe: auto` | exit 3, `no-recipe` |
| C061 | C059 with `solve: index_calculus` | exit 0, `ic_single_target`; `counts` and `certificates.ic` equal C059's |
| C062 | B1's C010 (curve A, `k = 1`) with `solve: index_calculus` | exit 0; `counts` and `certificates.ic` equal C010's paired run's |
| C063 | a curve over `GF(8)` taken over `GF(2^33)` (`r ≈ 2^29.4`); both pipelines `auto`, with a recipe | exit 0; `kic` and `rho-negation`; `ic-binary-s4` gated `recipe-not-taken` |
| C064 | a curve over `GF(2^9)` taken over `GF(2^27)`; `kic` named | exit 3, `subfield-too-large` |
| C065 | a curve over `GF(4)` taken over `GF(2^20)` (`e = 10`); `kic` named | exit 3, `even-extension-degree` |
| C066 | the challenge file with a certificate for its `r` (`2^129`, one level) | `ic check` passes; `order-composite` exact |
| C067 | C066 without the certificate's largest factor | exit 2, `certificate-invalid` |
| C068 | `y² = x³ + x` over a 201-bit `p`, with `p` certified two levels deep | `ic check` passes; `p-composite` exact |
| C069 | C068 with a witness failing at the second level | exit 2, `certificate-invalid` |
| C070 | C059 at `fidelity: F2` | `status: estimated`; `kic` and `rho-negation` |

**Disclosure.** The recipes of C059 and C063 follow §20's rules with
`e` for `n` (`make_cases.recipe`). Before freezing, each was tried on
v1's `ic workflow`, on the curves v1's constructor builds with the same
subfield, field degree and group order. Every run completed. Those
runs build the same `kic` pipeline, so they check the recipes, not
B2b's code.

## Class

**Robustness.** B2b widens what the tool accepts and runs. It changes
nothing on any input B2 already ran, except the `operation` label of an
imported index calculus run alone. It claims no speedup, and its pair
on a subfield curve is not speedup-eligible.

## Arms

- **The base**: the newest accepted baseline when B2b runs, recorded in
  B2b's manifest. B2 comes before B2b.
- **B2b**: that commit plus B2b's change, built with
  `IC_BUILD_COMMIT=$(git rev-parse HEAD)`.

## Measurements

The runs use [`../../harness/bround.py`](../../harness/bround.py) where
it applies.

1. **Tests.** `cargo test --release --bin ic`, and the library tests of
   every module B2b touches. The certificate check is tested on:
   - primes, certified one level and two levels deep;
   - composites, including Carmichael numbers and a strong pseudoprime;
   - every way a certificate can be malformed.
2. **Conformance.** `conformance/run.py --steps <the steps accepted>,B2b`
   (`bround.py conformance`).
3. **The pin, untimed.** All 90 suite rows from their v1 files. Every
   output must equal the base's (`bround.py pin`).
4. **The translation, untimed, on all 90 rows**, as in B1 and B2
   (`bround.py translate`).
5. **No slowdown, timed.** The base against B2b on `M1`'s 22 rows: five
   rounds ABAB, isolated, 220 processes (`bround.py timing`). B2b
   changes the single-target pricer, which every Koblitz row runs, so
   this is checked as in B1.
6. **The subfield sweep, untimed**
   ([`sweep.py`](sweep.py), pinned by [`SHA256SUMS`](SHA256SUMS)).
   - The pairs are every `k` in 2–8 and odd `e ≥ 3` with
     `16 ≤ n = k·e ≤ 40`: 21 pairs.
   - On each, the first curve by the script's rule.
   - Two known-answer targets on each, `kic` paired with `rho-negation`.
   - Every answer is replayed in the script's own arithmetic.
   - The pair `k = 2`, `e = 15` has no curve with `r > h`. For a
     composite `e`, `#E(GF(q^e))` is a multiple of `#E(GF(q^d))` for
     every `d` dividing `e`. The script lists that pair, so the sweep is
     20 curves and 40 runs.
   - It runs under the benchmark lock.

## Acceptance

B2b is **accepted** when all of the following hold:
- every test passes;
- every case the runner selects passes;
- the pin holds on every row;
- both parts of the translation check hold on every row;
- no size regresses beyond its A/A band, judged as in B1;
- the sweep has no panic, no exit other than 0 or 1, no missing report,
  no wrong answer and no complete run with an unverified answer.

A sweep run that ends `not_recovered` within its recipe is listed with
its counts. It is not a wrong answer and does not block acceptance.
The results say how many there were.

B2b is **stopped** if any v1 output differs. It is **rejected** if a
test or a case fails, the translation check fails, a size regresses, or
the sweep finds a panic or a wrong answer. A rejected B2b is fixed and
declared again, by amendment.

## Inadmissible

- Loosening a case after a run. The only exception is the `until` rule,
  applied by a later step.
- Changing a file under `conformance/v2-b2b/`, or `sweep.py`, after this
  declaration.
- Changing `conformance/run.py`'s rules, other than by a step's
  declaration.
- Weakening a rule of the certificate check to pass a case.
- Reporting the `kic` and `rho-negation` pair on a subfield curve as a
  speedup.
- Counting a contended run.

## Cost

- The tests and the conformance suite: about ten minutes.
- The pin and the translation: 270 untimed processes.
- The timing check: 220 processes, about 90 minutes.
- The sweep: 40 processes, a few minutes.

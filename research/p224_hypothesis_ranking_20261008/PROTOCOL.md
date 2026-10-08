# Protocol: P-224 checks behind a ranked hidden-weakness hypothesis list

Status: **frozen before the run** (2026-10-08). The report is
[`README.md`](README.md), written after the run.

## 1. What is being checked

A ranked list of four hypotheses was supplied for assessment. It describes how
a 1999-era designer could have weakened the NIST prime curves, ranked by how
well each fits the public evidence and the compute of the time:

1. a weak class recognisable from `(p, a, b)` without point counting, found by
   grinding seeds;
2. a weakness in the special-form prime rather than in the curve;
3. hidden endomorphism or isogeny structure;
4. a new ECDLP algorithm known only to the designer.

The list cites this repository as evidence: `b_seed_profile`, the
Solinas-correlation experiment, the cyclotomic `p − 1` finding, the
isogeny-walk thread and `p256_speculation`. Most of that evidence was gathered
on P-256. This protocol fixes the checks that give the P-224 counterparts and
the cross-curve facts the list compares P-224 with. It is an assessment under
`AGENTS.md` "Intermediate-cost and asymmetric-access hypotheses". It is not an
ECDLP measurement: no attack is run, no cost is priced, and no scoreboard row
is owed (`AGENTS.md` §7a lists no panel for structural checks).

## 2. Frozen inputs

- Curve parameters: `CurveParams::{p192, p224, p256, p384, p521,
  brainpool_p224r1}` in `src/ecc/curve.rs` and `src/ecc/curve_zoo.rs`, at
  the commit that adds this file. brainpoolP224r1 is the control: a 224-bit
  prime chosen pseudo-randomly, not in Solinas form.
- Seeds and printed `c` values: FIPS 186-4 Appendix D, transcribed in
  `nist_seeded_curves()` in `src/cryptanalysis/p224_hypotheses.rs`.
- External cross-checks: SafeCurves' factorisations of `4p − t²` (P-224 and
  P-256) and of the P-224 twist order, from <https://safecurves.cr.yp.to/>
  (`disc.html` and `twist.html`, fetched 2026-10-08). They are transcribed in
  `published_checks()` in `src/bin/p224_hypotheses.rs` and verified by
  product and primality, independently of how they were found.
- Null label `p224-b-null/v1`. Null draw `i` is
  `SHA-1(label ‖ i ‖ 0) ‖ SHA-1(label ‖ i ‖ 1)` reduced mod the P-224 prime.

## 3. Checks

| id | check | hypothesis | outcome rule |
|:--|:--|:--|:--|
| C1 | X9.62 seed → `c`, then `b²c ≡ −27 (mod p)`, all five NIST curves | 1 | pass iff all five hold |
| C2 | `v₂(p ∓ 1)` and factorisation of `p ∓ 1`, five NIST primes and the control | 2 | descriptive |
| C3 | factorisation of `4p − t²` for P-224 and P-256; conductor and `D_K` when complete | 3 | descriptive |
| C4 | factorisation of the P-224 and P-256 twist orders | implementation track | descriptive |
| C5 | embedding degree of P-224 and P-256, exact when `n − 1` factors | 3 (transfer) | descriptive |
| C6 | KS profile of `(b mod q)/q`, `5 ≤ q ≤ Q`, P-224 and P-256, `Q ∈ {10⁵, 10⁶, 10⁷, 10⁸}`; null panel of 1,000 draws at `Q = 10⁶` | 1 | see §4 |
| C7 | SafeCurves factorisations: product and primality | 3, implementation track | pass iff every listed factorisation checks |

Factoring budget, fixed: trial division below `2²⁰`, then Pollard–Brent with at
most `2²⁶` iterations per integer. A cofactor that does not split within the
budget is reported as unfactored. It is not rerun with a larger budget without
an amendment to this file.

Prime certainty: a factor is *proved* below `2⁴⁰` by trial division and below
`3.317·10²⁴` by Miller–Rabin on the first 13 prime bases (deterministic
there). Above that it is a *probable* prime: Miller–Rabin on twenty bases.

## 4. The b-residue test (C6): threshold, calibration, power

- **Flag.** P-224's asymptotic KS `p < 0.001` at any `Q`, or its statistic
  above all 1,000 null statistics at `Q = 10⁶` (empirical `p ≤ 1/1001`). A
  flag is reported as a flag. It needs its own protocol before anyone reads
  it as structure. The eight `(curve, Q)` cells are nested, so they are not
  independent tests.
- **Calibration gate.** The null panel's count of asymptotic `p < 0.05` should
  fall in `50 ± 3·√(1000·0.05·0.95)`, that is 30–70. If it falls outside, the
  asymptotic p-values are reported as miscalibrated and only the empirical
  p-value is read.
- **Power.** The test can detect only a class that biases `b mod q` for primes
  `5 ≤ q ≤ Q`. A class defined by any other property of `b` (a Legendre
  pattern of a derived quantity, a property of `j`, of the trace or of an
  isogenous curve) is invisible to it. A null result is therefore a statement
  about this family of classes only. It does not exclude hypothesis 1.

## 5. Disclosure of prior computation

These were run during development, before this file was frozen:

- the unit tests in `p224_hypotheses.rs`: C1 for all five curves, `v₂` for the
  five primes, the complete P-224 `p − 1` factorisation, the P-224 conductor
  `f = 3` and `D_K`, and the P-256 published factorisation;
- the existing P-256 b-profile test at `Q = 10⁵`. Its output reproduced the
  recorded `D = 0.011490`, `p = 0.1589`
  (`research/notes/cm-isogeny/RESEARCH_P256.md`).

These are deterministic algebraic facts, or a reproduction of a recorded
statistic. **The P-224 b-residue statistic (C6) had not been computed** when
this file was frozen.

## 6. Command

```bash
cargo build --release --bin p224_hypotheses
target/release/p224_hypotheses --out research/p224_hypothesis_ranking_20261008/run
```

Defaults: `--trial-bound 1048576 --rho-budget 67108864
--bounds 100000,1000000,10000000,100000000 --null-bound 1000000
--null-draws 1000 --null-label p224-b-null/v1`. `run/run.json` records the
arguments, revision, compiler, thread count and wall time per check. Wall time
is a practicality note, not a metric.

## 7. Stop conditions

A crash, a failed C1 or C7, or a reproduction of the P-256 `Q = 10⁵` row that
differs from the recorded value stops the analysis. The failure is reported as
it stands.

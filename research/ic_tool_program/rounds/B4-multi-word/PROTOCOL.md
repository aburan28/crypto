# B4: binary fields of three words and more, at F0

**Declared 2026-10-01, before any B4 measurement.** This is Track B in
`research/notes/index-calculus/IC_TOOL_PROGRAM.md` (§9, B4). Its design
is [`../../design/multi-word.md`](../../design/multi-word.md), and its
cases are [`../../conformance/v2-b4/`](../../conformance/v2-b4/cases.json).
Nothing below changes after the first measurement, except by a dated
amendment appended at the end.

## What B4 adds

**`rho-koblitz` and `kic` run on binary fields of 3 to 9 words**,
odd `n` with `127 ≤ n ≤ 574`, end to end at F0 where the subgroup is
small enough:
- multi-word kernels for the field, the curve, the class key, the walk,
  and `kic`'s select, build, collect, logs and descent (design §2);
- ECC2K-130's published points read at three words, in the polynomial
  basis they are published in (design §4).

**B4 changes no one- or two-word kernel.** The suite's rows and B3b's
instances run exactly what they ran. The two-word pipelines gain three of
v1's rules that B3 and B3b refused past one word: the hashed target, the
random target and the generator rule (design §2). Each is v1's own rule
at one word.

**The plan's B4 is three words** (plan §9; a class key at three words
needs `n ≤ 190`). B4 is generic in the word count, so it takes nine,
`n ≤ 574`. The standard curves
`sect233k1` to `sect571k1` then pass the width gate. Their F0 cost is far
over any budget, so they gain checks, routing and F1's estimates, not
answers.

**`--kic-multi`, a test switch on `ic price`, runs the multi-word
pipeline where a narrower one would run,** at its narrowest width,
`W = 3`, so that the pipelines can be compared end to end. It changes
nothing else. Sameness at `W = 1` and `W = 2` is a unit test (design §3).

**The report says which pipeline ran.** `ic.words` is 1, 2, or `W`.

## Class

**Robustness.** B4 adds instances the tool can solve. It claims no
speedup, and its multi-word figures are new rows, not gains.

## Arms

- **The base**: the newest accepted baseline when B4 runs, with B3b's
  stack, as a commit recorded in the manifest. B4 builds on B3b's field
  and curve.
- **B4**: that commit plus B4's change.

## Measurements

1. **Tests.** `cargo test --release --bin ic`, and the library tests of
   every module B4 touches or adds. They must include the design's §5:
   - the field against `binary_ecc` at every degree from 127 to 574, and
     on every registered binary modulus;
   - sameness at `W = 1` and `W = 2`, counter for counter, on the smoke
     size, `M1`'s rows at `icv1-f2m53-tm56619371-dac20a85`, and B3b's
     C078 and C080 instances;
   - v1's three rules past one word, each equal to v1's own at every
     one-word degree;
   - ECC2K-130's published points from B1's v2 document.
2. **Conformance.** `conformance/run.py --steps <accepted>,B4`, with
   B4's cases (below). The same runner also runs on the base, and its
   failures are recorded as what the base cannot do.
3. **The pin, untimed.** All 90 suite rows from their v1 files. Every
   output must equal the base's.
4. **No slowdown, timed.** By the chain's rule (`bround.py chain`), with
   B4 the next arm after B3b.
5. **F0 at three words.** Each of the six instances runs as `ic price` at
   F0 on two targets, its known-answer target and a public point `T001`.
   Each run is one process, isolated, on one core.

   | `n` | `a` | `log₂ r` | curve |
   |--:|--:|--:|:--|
   | 127 | 0 | 30.06 | `icv1-f2m127-t24589614856193402413-a00e5890` |
   | 137 | 0 | 37.48 | `icv1-f2m137-t574016314927011818501-11ceb4d6` |
   | 151 | 0 | 34.71 | `icv1-f2m151-tm97937651335354183120307-ee31df7e` |
   | 157 | 1 | 47.14 | `icv1-f2m157-t158060339213695215877259-816eafb2` |
   | 173 | 1 | 27.63 | `icv1-f2m173-tm67870783603944754053042229-ddf2de63` |
   | 179 | 0 | 39.58 | `icv1-f2m179-t1681527843948629186391379613-d1b73c24` |

   Recorded for each run: the logarithm, checked in the run and replayed
   outside it; `S` against rho's on the same point; every phase's cost
   and counts; the host manifest and the isolation record. A run has one
   hour.

## B4's cases

`conformance/v2-b4/`, written by its `make_cases.py` with B1's
generator code, which shares no arithmetic with the tool:

| case | what it checks |
|:--|:--|
| C088–C093 | `kic` at F0 on each three-word instance (design §4), paired with `rho-koblitz`: the known logarithm recovered and verified by both, `ic.words` 3 |
| C094 | `--kic-multi` on the smoke row's translation: the one-word pipeline's counts, logarithm and certificates |
| C095 | `--kic-multi` on C080's document (`n = 79`): the two-word pipeline's |
| C096 | C055's and C086's successor: on the challenge's document, `n = 131` passes the field gate, and `kic` and `rho-koblitz` refuse `r = 2^129` as `scalar-wider-than-127-bits` |
| C097 | C056's successor: v1's hashed target at `n = 83`, derived by design §2's rule; the derived point is the one the generator computes with its own BLAKE3, and rho alone stops at its step cap |
| C098 | the width gate: at `n = 577`, `kic` and `rho-koblitz` refuse the field as `field-wider-than-nine-words` |
| C099 | `sect163k1`: `n = 163` passes the field gate at three words, and both pipelines refuse `r = 2^162` by the same scalar gate |
| C100 | v1's random target past one word, at `n = 67`: drawn by design §2's rule, recovered and verified by both arms |
| C101 | v1's generator rule past one word, at `n = 67`: the generator found by design §2's rule, and the known logarithm recovered and verified by both arms |
| C102 | a 95-bit subgroup at `n = 127`: both pipelines admit it by field and scalar width, and the paired price at F0 is refused as over budget |

C031, C055, C056 and C086 name B4 as their `until`, so the runner retires
them when B4 is among the steps.

## Acceptance

B4 is **accepted** when:
- every test passes, sameness included;
- every case the steps select passes;
- the pin holds on every row;
- no size regresses beyond its A/A band;
- every run in measurement 5 recovers its logarithm, verified in the run
  and by the replay.

It is **stopped** if a one- or two-word output differs. It is
**rejected** if a test or a case fails, if a size regresses, or if a run
in measurement 5 gives a wrong answer or none within its hour.

**What an accepted B4 may claim:** the tool runs `rho-koblitz` and `kic`
end to end on Koblitz curves over fields of up to nine words, at F0
where the subgroup allows, with the answers verified. Subgroups stay
below `2^127`, as at two words: a larger `r` is refused as
`scalar-wider-than-127-bits`, which B7b lifts for F1 at `n = 131`.
Nothing about the cost of ECC2K-130's subgroup, which is B7b's.

## What B4 does not do

- F1 at `n = 131` (B7b), and scalars of 127 bits or more, which F1 there
  needs (B7b).
- Prime fields beyond one word for `kic` (B5).
- Any method using a proper intermediate subfield.

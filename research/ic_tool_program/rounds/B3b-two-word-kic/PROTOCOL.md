# B3b: the index calculus on two-word binary fields, at F0

**Declared 2026-10-01, before any B3b code.** This is Track B in
`research/notes/index-calculus/IC_TOOL_PROGRAM.md` (§9, B3b). Its design
is [`../../design/two-word-kic.md`](../../design/two-word-kic.md), and its
cases are [`../../conformance/v2-b3b/`](../../conformance/v2-b3b/cases.json).
Nothing below changes after the first measurement, except by a dated
amendment appended at the end.

## What B3b adds

**`kic` runs on Koblitz curves over two-word binary fields,**
`63 ≤ n ≤ 126`, end to end at F0:
- two-word kernels for select, build, collect, logs and descent
  (design §2);
- any subgroup order `r < 2^127`. The sparse solver's residues go to two
  limbs from `2^63` up.

**B3b changes no one-word code.** The suite's rows, all at `n ≤ 61`, run
exactly what they ran.

**`--kic-wide`, a test switch on `ic price`, runs the two-word pipeline
where the one-word pipeline would run.** It exists so that the two can be
compared end to end, and it changes nothing else.

**The report says which pipeline ran.** `ic.words` is 1 or 2.

## Class

**Robustness.** B3b adds instances the tool can solve. It claims no
speedup, and its two-word figures are new rows, not gains.

## Arms

- **The base**: the newest accepted baseline when B3b runs, as a commit
  recorded in the manifest. B3 comes before B3b, since B3b builds on
  B3's field and curve.
- **B3b**: that commit plus B3b's change.

## Measurements

1. **Tests.** `cargo test --release --bin ic`, and the library tests of
   every module B3b touches or adds. They must include:
   - **the two-word kernels against the general arithmetic**
     (`binary_ecc`), at every field degree from 63 to 127 and on the
     gate's modulus;
   - **sameness at one word** (design §3), each with its suite recipe:
     the smoke size `icv1-f2m31-tm90707-c95f16f5`, and `M1`'s rows at
     `icv1-f2m41-tm2308219-7f48b14a`,
     `icv1-f2m53-tm56619371-dac20a85` and
     `icv1-f2m59-tm943548413-98844ecc`. Sameness means the same base, the
     same table, the same relations in the same order, the same
     logarithms and the same counters;
   - **the two-limb residues against BigUint**. A synthetic sparse
     system is solved modulo a 70-bit and an 81-bit prime, and the
     solution is checked.
2. **Conformance.** `conformance/run.py --steps <accepted>,B3b`. B3b's
   cases are C078–C085:
   - **C078–C082**: two-word `kic` at F0, paired with rho on the same
     point, with a known answer, on the five instances below.
   - **C083**: C054's successor (B3 marked C054 `until: B3b`). The gate
     curve's paired price is refused as over budget, both pipelines
     admitting it by width.
   - **C084 and C085**: `--kic-wide` against the one-word pipeline, on
     the smoke row's translation and on the `n = 61` row's. They must
     give the same counts, logarithms and certificates.
   The same runner also runs on the base, and its failures are recorded
   as what the base cannot do.
3. **The pin, untimed.** All 90 suite rows from their v1 files. Every
   output must equal the base's.
4. **No slowdown, timed.** The base against B3b on `M1`'s 22 rows: five
   rounds ABAB, isolated, 220 processes, with the programme's runner.
   The figure is the paired cold-time ratio per size.
5. **F0 at two words.** Each of the five instances runs as `ic price` at
   F0 on two targets:
   - its known-answer target, from C078–C082's document;
   - a public point `T001`, hashed from a public label as the gate's
     `T001` is, whose logarithm nobody knows before the run.

   v2's hashed-seed targets wait for B4 past one word (B3's C056), so
   both targets are points. Each run is one process, isolated, on one
   core.

   | case | curve | `n` | `log₂ r` | note |
   |:--|:--|--:|--:|:--|
   | C078 | `icv1-f2m67-tm19346764963-82c84cca` | 67 | 26.19 | |
   | C079 | `icv1-f2m67-t19346764963-e760df38` | 67 | 35.99 | |
   | C080 | `icv1-f2m79-tm420247971347-2a24b892` | 79 | 33.74 | |
   | C081 | `icv1-f2m71-tm48653080717-f25c4638` | 71 | 49.06 | |
   | C082 | `icv1-f2m83-t6151469093347-cdcc5432` | 83 | 52.93 | the gate's field, not the gate's curve |

   Recorded for each run:
   - the logarithm, checked in the run and replayed outside it;
   - `S`, against rho's on the same point;
   - every phase's cost and counts;
   - the host manifest and the isolation record.
6. **The two-word premium, a stage diagnostic.** At `n = 47, 53` and `61`,
   on `M1`'s rows, three rounds, `--kic-wide` against the default:
   - the build per stored pair;
   - the scan per summand;
   - the descent per probe.

   Each is a ratio of the two pipelines' own reported costs. B7b reads
   F1 at `n = 83` with it. It is not a speed figure.

## Acceptance

B3b is **accepted** when all of the following hold:
- every test passes, sameness included;
- every case the steps select passes;
- the pin holds on every row;
- no size regresses beyond its A/A band;
- every run in measurement 5 recovers its logarithm, verified in the run
  and by the replay.

It is **stopped** if a one-word output differs. It is **rejected** if a
test or a case fails, if a size regresses, or if a run in measurement 5
gives a wrong answer or none within the budget.

**What an accepted B3b may claim:**
- `kic` solves the listed two-word instances at F0, with their costs and
  their ratios to rho on the same points;
- `C082` is on the gate's field.

**What it may not claim:**
- anything about the gate's curve (`E_0`, `r ≈ 2^81`), which needs F1
  (B7b);
- a speedup. Measurement 6 is a stage diagnostic.

## Inadmissible

- Changing a one-word kernel to make sameness hold.
- Dropping an instance or a target, or loosening a case, after a run.
- Counting a contended run.
- Reading an `n = 83` figure here as the gate's.

## Cost

- Tests and conformance: minutes. C082 runs `kic` at `r ≈ 2^53`, for a
  few minutes.
- The pin: 90 untimed processes.
- The timing check: 220 processes, about 40 minutes.
- Measurement 5: ten processes, minutes each at the largest.
- Measurement 6: 36 runs, about 20 minutes.

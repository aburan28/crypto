# B3: two-word binary fields, and rho alone at `m = 83`

**Declared 2026-10-01, before any B3 measurement.** This is the first
half of Track B's B3 in `research/notes/index-calculus/IC_TOOL_PROGRAM.md`
§9: two-word fields and the rho reference on them. The index calculus
on two-word fields is the second half, **B3b**, declared separately.

Unlike B1, B3's code was written before this declaration, on a local
branch, so the declaration states what that code does. The rates quoted
under "What the code does" are development diagnostics taken while
building it, not B3's measurements. Nothing below changes after the
first measurement, except by a dated amendment at the end.

## What B3 does

- **`gf2_wide`: GF(2^n) for n ≤ 127, an element in a `u128`.**
  - A product is four 64-bit carry-less multiplications.
  - Reduction folds when the modulus's low terms fit one word, as for
    every trinomial and pentanomial the standards use. That is three
    folds on the gate's modulus `z^83 + z^45 + z^2 + z + 1`. Any other
    modulus takes polynomial Barrett.
  - Inversion is Itoh–Tsujii.
  - On x86-64 an SSE kernel keeps elements in vector registers through
    the product and the reduction. The modulus's shape is fixed at
    compile time, with a run-time detected PCLMULQDQ, SSSE3 and POPCNT
    and a portable fallback. Arm64 uses PMULL.
  - Tests check every operation against the general big-integer field
    at degrees 31 to 127, and the kernels against each other.
- **`koblitz_wide`: Koblitz curves over those fields, and `WideRho`.**
  `WideRho` is `ParallelRho`'s signed-Frobenius r-adding walk, ported
  step for step: the same seed draws, jumps, starts, jump index,
  distinguished points, normal element, least rotation and tie rule,
  look-ahead, cycle window and escape, with coefficients in parts.
  - A lane keeps the point the addition left and the automorphism that
    names its class representative, and adds the jump's Frobenius image
    rather than the jump. A step therefore reads the representative's
    abscissa and one bit of its ordinate, not both coordinates.
  - It charges its group operations exactly as `ParallelRho` charges
    them.
  - A test checks, at `n = 23`, 31 and 37, that the two walks agree on:
    the steps, walks, distinguished points, rounds and answer; the setup
    and target operations; and every shared counter.
  - Another test checks that the vector kernel's walk equals the
    portable one's at `n = 83`.
- **Routing (design §5).** `rho-koblitz` gets its own gate. It takes
  Koblitz curves with odd `n` up to 126 and `r < 2^127`. Two new gate
  codes, defined in the design's amendment below, refuse instances past
  that. `kic` keeps its one-word gates until B3b.
- **`solve: rho`.** Rho alone runs at F0, on one thread, with no
  recipe, at every degree `rho-koblitz` admits.
  - At `n ≤ 62` it gives the paired price's rho arm's counts and answer
    (C057, C058).
  - Its timing follows the single-target rule. The jump table is the
    reusable set-up and is timed apart. The online interval opens at the
    first target-dependent walk and closes on the checked logarithm,
    with exclusive phases `rho_solve` and `recovery_check`.
  - After the interval comes the replay, `[d]G = Q` in the general
    big-integer arithmetic, timed separately.
  - The report carries the walk's ledger, `S = target operations / √r`
    against the floor `√(π/2A)` with `A = 2n`, and a replay certificate
    with its SHA-256.
- **Targets under `solve: rho`.** A point or a known logarithm is taken
  at any degree. v1's hashed and random targets are taken where v1's
  curve is defined (`n ≤ 63`) and refused for B4 beyond that.

**Why 126, not 127.** A class key packs to `2(x + 1) + s`, as
`FastPoint::pack` does, and at `n = 127` that needs 129 bits. `n = 127`
therefore waits for B4's three-word fields.

**What B3 does not do.**
- **The index calculus at `n > 62`.** That is B3b (the plan's "index
  calculus at F1").
- **`rho-negation`.** B2 imports it first.
- **The `m = 83` gate itself.** AGENTS.md §8a asks for the index
  calculus's baseline and candidate there. B3 supplies the rho side on
  one target, and the index-calculus side stays unestablished.

**Development diagnostics.** These are rates of the walk alone at
`n = 83`, from `examples/wide_rho_rate.rs` on one pinned core of this
host. They size the gate run and are not results:
- 352 ns a walk operation at the first port;
- about 140–160 ns after the class-only canonicalisation and the vector
  kernel;
- 1,270 instructions a step by callgrind.

## Class

**Robustness**, with one **measurement**. B3 widens what the tool runs.
It changes nothing on any input the tool already ran, and it claims no
speedup. The gate's rho at F0 is a measurement of the reference at
`n = 83`, not a comparison.

## Arms

- **The base**: the newest accepted baseline when B3 runs, recorded in
  B3's manifest. B1 comes before B3, because B3 routes through B1's
  schema.
- **B3**: that commit plus B3's change, built with
  `IC_BUILD_COMMIT=$(git rev-parse HEAD)`.

## Measurements

1. **Tests.** `cargo test --release --lib -- gf2_wide koblitz_wide` and
   `cargo test --release --bin ic`.
2. **Conformance.** `conformance/run.py --steps B0,B1,B3` on B3, which
   runs:
   - v1's cases;
   - v2's cases with the `until` rule applied, which retires C027;
   - B3's C052–C058.

   The same runner also runs on the base, and its failures are recorded
   as what the base cannot do.
3. **The pin, untimed.** B3 runs all 90 suite rows from their v1 files.
   Every output must equal the base's: the counts, both arms' scalars,
   rho's counts and the verification flags.
4. **No slowdown, timed.** The base against B3 on `M1`'s 22 rows from
   their v1 files: five rounds ABAB, isolated, 220 processes, with the
   programme's runner. The figure is the paired cold-time ratio per
   size. B3 touches no code those rows run, so this is a check, not a
   hope.
5. **The gate's rho at F0.** `ic price` on the gate document
   (`conformance/v2/params/gate-m83-T001.json`) with
   `method = {"solve": "rho", "fidelity": "F0", "rho": {"pipeline": "auto", "seed": 2293761}}`.
   That seed is the document's own.
   - One run, isolated, on one core, under the default budget of
     86,400 s.
   - The floor is `(πr/4n)^{1/2} = 1.513 × 10^{11}` steps. At the
     development rate that is about 6.3 hours, and the budget is about
     3.8 times it.
   - The chance that a rho run needs more than 3.8 times its expectation
     is about `exp(−π·3.8²/4) ≈ 10^{−5}`.
   - Recorded: the steps and walk operations against the floor; `S`
     against `√(π/4n)`; the online interval and its phases; the replay;
     the certificate and its SHA-256; and the host manifest and isolation
     record.

## Acceptance

B3 is **accepted** when all of the following hold:
- every test passes;
- every case `conformance/run.py --steps B0,B1,B3` selects passes;
- the pin holds on every row;
- no size regresses beyond its A/A band, judged as in B1;
- the gate's rho recovers the logarithm, checked in the walk and
  replayed in the general arithmetic, with its certificate.

It is **stopped** if any v1 output differs. It is **rejected** if a test
or a case fails, or if a size regresses.

If the gate run exhausts its budget, B3's code may still be accepted on
the other conditions. The run is then recorded as a timeout and the
gate's rho side stays unestablished. The walk is deterministic in its
seed, so it is not repeated under the same seed; a new seed would be a
new, declared run.

**What an accepted B3 may claim.** A rho measurement at F0 at `n = 83`
on the gate curve's target T001, on this host, with its rate and `S`.
It may not claim that the gate is discharged, nor anything about the
index calculus at `n = 83`.

## Inadmissible

- Loosening a case after a run. The only exception is the `until` rule,
  applied by a later step.
- Changing a file under `conformance/v2-b3/` after this declaration.
  `SHA256SUMS` pins B3's files.
- Changing `conformance/run.py`'s rules other than by a later step's
  declaration. A later step may add its own case set beside B3's.
- Running the gate under another seed and reporting the better run.
- Counting a contended run.

## Cost

- The tests and the conformance suite: a few minutes.
- The pin: 90 untimed processes.
- The timing check: 220 processes, about 90 minutes.
- The gate's rho: about 6.3 hours, and at most 24.

That is about 8 to 9 hours of the benchmark lock in all.

## Amendment 1 (2026-10-01, before any measurement)

**Measurement 4, the timing check, runs in the chain of Track B steps.**
- **What it was:** the base against B3 on `M1`'s 22 rows, five rounds
  ABAB, 220 processes, for this step alone.
- **What it is now:** one interleave over the newest accepted baseline
  and every Track B step's arm, in the queue's order: the baseline, B0,
  B1, B3, B2, B2b, B7a and B3b (`bround.py chain`). It runs on `M1`'s 22
  rows, five rounds, isolated, and the order reverses every other
  round.
- **B3's figure** is the paired cold-time ratio of the arm before it
  in the chain over B3's own, per size. It is read against R01's A/A
  bands, as before.
- **Why the chain is a valid base.** Each step's arm is built on the one
  before, so the arm before B3 in the chain is its base. The two run
  back to back in every round, with ten pairs a size, as in its own
  ABAB.
- **What it saves.** The seven steps' checks take 880 processes
  together, against 1,540 one by one.
- **If an earlier step is rejected,** the steps after it wait. The chain
  runs again from that step, once it is fixed or removed, since the
  later arms carry its change.
- **Nothing else changes:** the acceptance rule, the A/A bands and the
  other measurements.

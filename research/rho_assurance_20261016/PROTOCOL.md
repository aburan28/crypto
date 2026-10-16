# Strong-rho assurance round: protocol

Declared on 2026-10-16, before any evidence in this directory was produced.

Every Koblitz `vs_rho` comparison prices index calculus against the strong rho
(`src/cryptanalysis/koblitz_strong_rho.rs`, `koblitz_rho_fixture … strong`).
If that rho is wrong, or needs more steps or more work per step than it should,
every comparison against it is unfair. This round asks three questions about
it. It changes no committed result.

## Q1. Is it correct? (formal and exhaustive checks)

- **TLA+ (`tla/`).** A model of `solve` over `Z_r`: lockstep lanes, walk starts,
  distinguished points, the table, collisions, fruitless-cycle abandon and
  escape, and caps. TLC checks exhaustively:
  - the representation invariant `P = a·G + b·Q` on every lane and table entry;
  - canonical lane points;
  - no failed collision;
  - a correct answer whenever it answers.

  Mutated specs with one planted bug each must fail. A model that cannot fail
  checks nothing.
- **Lean (`lean/`).** Core Lean 4, no `sorry`:
  - the least-rotation lemma behind `NormalBasis::least_rotation`;
  - Montgomery batch inversion as the Rust implements it, zero entries
    included;
  - the walk algebra: steps, the orbit multiplier, escapes and the collision
    formula, as preservation of the representation invariant.

  `PROOF_MAP.md` records each theorem's hypotheses and the gap to the Rust
  code.
- **Rust.**
  - An exhaustive test of `least_rotation` over every word of every width
    n ≤ 16.
  - A debug-build assertion of the representation invariant at every
    distinguished point `solve` stores or matches. It ties the code to the
    TLA+ invariant on every test solve. Release builds are unaffected.

Pass: every expected-pass model, proof and test passes, and every mutant fails.
Anything else is reported as found.

## Q2. Does it take the number of steps it should? (counted units)

- **Command.** `koblitz_rho_fixture n a signed_frobenius 400 strong`, with
  `RAYON_NUM_THREADS=1`, for (n, a) = (41, 0), (47, 1) and (53, 0): 400
  published synthetic fixtures each. Each fixture has its own target and jump
  seed, as the example derives them.
- **Recorded.** Per fixture: `walk_steps / ideal_steps`, where
  `ideal_steps = √(πr/2A)` and A = 2n is the expected first-collision time of
  an ideal random walk on the orbit classes. Also recorded: walks, fruitless
  and capped walks, and the share of steps spent in abandoned walks.
- **Reference.** An ideal random mapping gives a mean ratio of 1 with a
  standard deviation of about 0.52 per fixture. Known inflations:
  - a 32-jump additive walk, about 1/√(1 − 1/32) ≈ 1.016;
  - the distinguished-point lag, about `lanes·2^dp_bits` = 512 steps per
    solve.

  At these sizes both are a few percent at most.
- **Decision, per size, on the 95% interval of the mean ratio:**
  - upper end below 1.10: consistent with an efficient walk;
  - lower end above 1.10: inflated, report the factor and look for the cause;
  - otherwise: not established.

Counted units do not depend on contention, so Q2 may run beside other work.

## Q3. What does a step cost? (instructions, then native time)

- **Callgrind** on `koblitz_rho_fixture 53 0 signed_frobenius 4 strong`.
  - Instructions per walk step, overall and by inclusive function:
    canonicalization (basis conversions and least rotation), batched
    inversion, the point addition, the two hashes, the table and walk starts.
  - **Field-multiplication unit M.** Instructions per `Gf2::mul` at n = 53,
    from two runs of the field-operation example that differ only in the
    number of multiplications (difference ÷ count).
  - **The step's cost in M.** Set against the floor of an affine batched step:
    1 slope multiplication, 1 squaring and 1 multiplication for y, 3
    multiplications for the batch inversion, and one inversion shared by 32
    lanes. That is about 5M + 1S + I/32, derived, not measured.
- **Native time.** ns per step = `walk_ms / walk_steps` from 5 runs under
  `isolated_bench run --cpus 3`, with `RAYON_NUM_THREADS=1`. Report the median
  and range.
- **IPC.**
  - Cycles, instructions and IPC come from `perf stat` through the new
    `perf-ipc` skill, but only on a host whose PMU exposes hardware counters.
  - This container's VM exposes none: `/sys/bus/event_source/devices/` has no
    `cpu` PMU, and `perf stat -e cycles:u` reports `<not supported>`.
  - IPC is therefore recorded as unavailable here, not estimated.

Q3 characterizes the cost. A step above about twice the floor in M is reported
as a candidate inefficiency of the reference, with its largest components.

## Limits declared in advance

- **Scope of the formal checks.** TLC checks small instances of a model, not
  the Rust, and not every r. The Lean proofs model the subgroup as a field and
  words as bit lists. The link to the Rust is the exhaustive and debug-assert
  tests, not a refinement proof.
- **Fixtures.** Q2 uses the published synthetic fixtures of three sizes. It
  says nothing about m = 83 beyond the shared code.
- **Hosts.** One x86-64 cloud VM with no hardware counters. No other host is
  measured.

## Addendum, declared before the change was written: Q4, a faster basis conversion

Q3's line-level profile (`runs/q3/`) was taken from a line-table build of the
same source on 1 fixture. It shows:

- **Where the step goes.** Canonicalization is about 585 instructions of a
  1,196-instruction step. About 425 of those are the four byte-table basis
  conversions (`Linear::apply`), at about 14 instructions per table lookup:
  byte extraction by a variable shift, the generic XOR and the slice
  iterator.
- **The candidate.** `Linear::apply` shifts the word down one byte per table.
  A new `apply2` converts x and y in one pass over the tables. `canonicalize`
  uses `apply2` both ways.
- **Unchanged.** The tables, their values and every result are unchanged, so
  the walk must be bit-identical.
- **Identity gate.**
  - Every output line of `koblitz_rho_fixture 53 0 signed_frobenius 4 strong`
    and `… 61 0 … 4 strong` is identical once the `*_ms` and `peak_rss_bytes`
    timing fields are set aside.
  - The strong-rho unit tests pass, including the wide-field instantiation's
    step-for-step test.
- **Measure.**
  - Callgrind instructions per walk step on the n53 4-fixture command, one run
    per arm, deterministic. This is the primary metric.
  - Native time: `walk_ms / walk_steps` under `isolated_bench`, an A/A pass
    (baseline × 6) and 5 interleaved A/B rounds, once no other job runs on the
    machine. A gain inside the A/A spread is not established.
- **Baseline.** The library at this round's parent commit plus the Rust test
  and debug-assertion additions above. Those compile out of release builds.

### Q4, second part (declared after the first Callgrind pair, before this code was written)

The first pair (`runs/q4/callgrind/`) cut instructions per step from 1,196 to
1,014, with identical outputs. Its profile shows two more costs that change no
result:

- **`add_mod`.** The one-word scalar `add_mod` reduces with a 128-bit `%`
  (`__umodti3`), about 15 instructions per step. Every caller passes operands
  already below the modulus: the walk's `a + c`, the start cursor and the
  escape's doubling. One conditional subtraction gives the same value. A
  `debug_assert` pins that precondition.
- **The table's hasher.** The distinguished-point table hashes its keys with
  std's randomly keyed SipHash, about 6% of instructions. A fixed 64-bit mixer
  computes the same map, because lookups and inserts do not depend on the
  hash, only its speed does.

Both join the candidate. The identity gate and the measurements are the ones
declared above. The final candidate is measured against the same baseline
binary.

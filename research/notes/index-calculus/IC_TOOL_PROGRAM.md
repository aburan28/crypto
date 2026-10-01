# The `ic` tool programme: fidelity, speed and generality, baseline by baseline

**Status (2026-10-01): declared. Nothing below is measured unless it
says so.** Execution is tracked round by round in
[`research/ic_tool_program/`](../../ic_tool_program/README.md), one
pull request per round.

**User direction (2026-10-01).** Make the `ic` binary the repository's
go-to, highest-fidelity index-calculus emulation tool. Make it as fast
as it can be made, round after round: establish a baseline, iterate,
measure, and establish the next baseline. Make it robust for anything a
user can pass as parameters, across field families.

This note turns that into a programme under AGENTS.md:
- §1 says what the programme may claim.
- §2 states the goals, each with a test.
- §3 defines fidelity levels.
- §4 freezes the suite.
- §5 fixes the measurement.
- §6 is the loop and §7 the ledger.
- §8 and §9 are the two backlogs.
- §10 says what does not count, §11 when to stop, and §12 what runs first.

## 1. What the programme may claim

- **"Faster" is a comparison with this programme's previous baseline.**
  It is measured end to end on the frozen suite (§4), cold, with
  identical outputs.
- **"Fastest worldwide" cannot be certified.** No published
  index-calculus implementation reports a matched instance under this
  repository's accounting. A claim against an outside implementation
  needs one of two things:
  - that implementation run on the same instance, in the same unit; or
  - its published figure on an identical instance, with its accounting
    stated and converted at a measured factor.

  §8's item A8 collects the candidates. Until one is matched, the
  programme says "faster than its own previous baseline", and nothing
  more.
- **Engineering, not advance, by default.** Making the same algorithm
  faster lowers `S` and leaves the ratio to the floor flat: AGENTS.md
  §3's engineering class. A round is an advance only if its ratio to the
  floor falls, and it is labelled by that test.
- **The method's verdict stands until a round changes it.** Ledger §23
  (2026-10-01), on one unseen point at `2^36.6`–`2^47.2`:
  - online, the index calculus is 8.8–17.5× faster than rho on the same
    point;
  - cold, it is 4.9–31.6× slower;
  - online, it is 2.0–8.0× above a generic walk given the same
    precomputation (a model).

  By AGENTS.md §2 the method is not faster than rho. The programme makes
  the tool faster. It does not presuppose that the method crosses rho,
  and it keeps rho as strong as it can make it (A5).

## 2. Goals

| goal | statement | test |
|:--|:--|:--|
| **G1 fidelity** | The tool runs the problem named in its input: the exact field, curve, subgroup, generator and Frobenius action. Every phase is priced. One unseen target is solved against rho on the same point. Answers are verified in the process and replayed outside it. Where a full run does not fit, the tool says so and runs what does fit, at a labelled level (§3). | every suite row carries its identities, its level and a replayed answer; no row's level is overstated |
| **G2 speed** | The end-to-end cold cost of one verified target falls, round by round, with identical outputs. | each accepted round's paired cold-time ratio over the previous baseline has its 95% interval above 1 at the sizes it declares, and no size regresses beyond its A/A band (§5) |
| **G3 generality and robustness** | Every well-formed parameter set either runs correctly or is refused with a precise reason. No input produces a panic or a wrong answer. | the conformance suite (§4) passes, and the fuzzing budget (B6) finds no panic and no wrong answer |

## 3. Fidelity levels

Every output names its level.

| level | what ran | what it may claim |
|:--|:--|:--|
| **F0, full** | Every phase, cold, to a verified scalar on one unseen target, with rho on the same point and both answers replayed | a measurement |
| **F1, sampled** | The real field and curve, with each phase's kernel run on a sample: relation yield per trial, build cost per row, linear algebra per non-zero, descent per trial. The totals are extrapolated. | an extrapolation, labelled, with its samples, its exponents and its formula |
| **F2, model** | formulas only | a model |

The `m = 83` gate (AGENTS.md §8a) needs F0 at `n = 83` to count as
evidence there. F1 at `n = 83` or `131` is an extrapolation, and it
never discharges the gate.

## 4. The frozen programme suite, v1

The suite is frozen with baseline v0 (round R01) in
`research/ic_tool_program/suite/v1/`. Its `SUITE.json` holds every
file's SHA-256. A later version adds a directory and keeps v1.

- **The speed suite (S).** The eleven Koblitz curves the tool runs
  today: §20's nine and §23's two, `r = 2^18`–`2^47.2`.
  - Each uses §20's frozen recipe at its swept sizes and §20's model
    optimum at the two it never ran, as in §23.
  - **Eight processes per size.** The set-up depends on the recipe
    seed, not on the target, so the suite varies both:
    - the four seeds of §20's measurement sets `M1`–`M4` (201–204);
    - two public hash-to-curve targets per seed, §23's `T01`–`T08` in
      order.

    `M1`'s two rows at §23's six sizes are §23's own files, so v0 can
    be checked against §23 directly.
  - **Holdouts.** Each round also draws two holdout processes per size,
    on seed 205 with fresh targets, from a target seed declared before
    the round's first run.
  - **Each row** is `ic price --single-target` with its default three
    in-process repetitions, one process per row, through
    `tools/isolated_bench.py`, as in §23.
- **The smoke tier.** It is for development and for pinning outputs
  before any timed run. It is not evidence. It has two parts:
  - `E_0/GF(2^31)` (`r = 2^20.5`), AGENTS.md §8b's primary exploratory
    size, from `docs/ic/params/k0n31`, which CI's
    `ic-e2e-benchmark.yml` already runs;
  - the S suite's three smallest sizes.
- **The conformance suite (C).** Parameter files across field families,
  valid and invalid, each with its expected outcome:
  - a verified scalar, at a stated level; or
  - a refusal with a named reason code.

  It grows with Track B and never shrinks. A file once added stays.
- **The gate files.** Two curves, run at F1 until F0 fits:
  - `E_0/GF(2^83)` with `z^83 + z^45 + z^2 + z + 1`, and its subgroup,
    generator and cofactor clearing frozen as AGENTS.md §8a requires;
  - `E_0/GF(2^131)`, the challenge field.
- **Equivalence with the WDSat regression suite.** AGENTS.md §8 asks
  every index-calculus performance iteration to run the frozen WDSat
  suite. The alternative is to freeze an equivalent matched suite under
  the parent accounting contract and say why WDSat's protocol is
  inapplicable.
  - It is inapplicable here. The Koblitz collection pipeline decomposes
    points by table lookup, not with a SAT or Gröbner solver, so that
    suite has no stage of this pipeline to measure.
  - The S suite is the equivalent suite, frozen under
    `research/index_calculus_baseline_20260914/ec_index_calculus_contract.json`'s
    accounting.
  - A round that touches a solver the WDSat suite covers also runs that
    suite, as AGENTS.md §8 requires. Such solvers are
    `ic descent --solver` and the Weil-descent paths.

## 5. Measurement, per round

AGENTS.md §2, §5, §8 and §10 and the single-target rule apply
unchanged, as in §23. On top of them:

- **The baseline.** The baseline binary is built from the unmodified
  head and kept outside the tree, with its SHA-256 and a host manifest.
- **Pinned outputs.** Baseline and candidate must give identical results
  on every S-suite row and every smoke row:
  - counts, relations and logarithms;
  - rho's walks and the recovered scalars.

  A round that changes any of them changes the algorithm, and it is
  declared as such.
- **The primary metric of an engineering round** is the paired ratio of
  the index calculus's cold time, baseline over candidate, per size.
  - Cold time is the set-up plus the online interval, at the median
    in-process repetition.
  - Each row runs in ABAB order over five rounds, isolated.
  - The ratio is the geometric mean of the paired ratios over the size's
    eight rows and five rounds. Its 95% interval is a `t` interval on
    the logarithms, §22's method.
  - Beside it, the round reports:
    - the online interval's ratio and each set-up phase's ratio;
    - an instruction-count ratio from `valgrind --tool=callgrind` at two
      sizes, as the deterministic cross-check.
- **`S` in a pinned unit.** Across baselines, `S` is quoted in v0's
  unit: v0's median batched addition (`unit_ns`) at each size, from R01,
  on the reference host. §21.4 showed that the in-process unit moves
  between binaries. Each process's own unit is reported beside it.
- **The rule's comparison**, at every new baseline. The comparison is
  one unseen target, the index calculus against rho on the same point,
  online and cold.
  - It runs at §23's six sizes with 64 targets, as §23 did.
  - It moves the scoreboard's §23 row.
- **The A/A.** Every round runs the baseline against a copy of itself,
  on `M1`'s rows for five rounds. That gives the noise floor, and a
  ratio inside the A/A band is not a gain. It is needed every round
  because this container is rebuilt between sessions, so the host can
  change.
- **Threads.** One thread is primary (`RAYON_NUM_THREADS=1`). A parallel
  round also reports all four cores, and it must not regress one thread.
- **The hardware class.** The reference host is x86-64 with AVX-512,
  PCLMULQDQ, VPCLMULQDQ and GFNI: this container, whose manifest is §23's
  `runs/host.json`.
  - Results hold for that class. Arm64 and GPU are named as unmeasured.
  - Kernels use runtime feature detection with a portable fallback, and
    CI tests the fallback.

## 6. The loop

Each round is one branch and one pull request:

1. **Declare** before any candidate code. The declaration goes in
   `research/ic_tool_program/rounds/R<k>-<slug>/PROTOCOL.md` and states:
   - the hypothesis;
   - the phase it targets, and that phase's share in the current
     profile;
   - the predicted ratio, the success condition and the stop condition;
   - the classes the round could fall in.
2. **Implement**, with tests. Outputs are pinned on the smoke tier
   before any timed run.
3. **Measure** the declared comparison: the A/A first, then ABAB.
4. **Decide.** A candidate is accepted if:
   - it meets the declared success condition;
   - its outputs are identical;
   - no size regresses beyond its A/A band.

   A rejected round keeps its code on its branch, linked from the
   ledger, and its numbers in the ledger.
5. **Re-baseline.** An accepted candidate becomes baseline `v<k>`, with
   its binary hash, host and suite results. The previous row stays as
   the "before".
6. **Update** the ledger (§7), the profile, and the scoreboard's
   programme panel, in the same PR (AGENTS.md §7).

A Track B round may change no speed at all. It must then show:
- the S suite unchanged, with identical outputs and time inside the A/A
  band;
- the conformance suite grown and passing.

## 7. The baseline ledger

Rows are baselines and columns are one unit. R01 filled v0 (2026-10-01),
and R03 v1 (2026-10-01). A baseline's `S` comes from the round that made
it, beside its base's `S` in the same runs: runs hours apart on this host
differ by several per cent, so only a paired ratio compares across
rounds.

| baseline | commit | class | `S` cold at the six top sizes (v0 unit) | cold-time ratio over the previous [95%] | online / rho | cold / rho | C suite | largest `n` at F0 | PR |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| v0 | `46ae2014` (`src/` tree `003badc2`) | accounting | 4.43, 7.65, 5.42, 7.53, 8.77, 14.49 | — | IC 8.8–17.5× faster (§23) | IC 4.9–31.6× slower (§23) | none yet | 61 | #1104, R01 results |
| v1 | `30f6c153` on v0′ `c1a2e5f8` | engineering | 4.39, 4.45, 5.83, 8.16, 9.43, 15.61 (R03's runs; v0′ there: 4.45, 8.00, 5.90, 8.04, 9.17, 15.82) | 1.771 [1.694, 1.851] at `2^38.0`, holdouts 1.922 [1.756, 2.104]; 0.970–1.049 at the nine prime sizes | not re-measured | not re-measured | none yet | 61 | #1119, #1166 |

**Rounds that did not become baselines.** A rejected round keeps its
numbers here and its code on record (§6, step 4). Its ratio is paired
within its own runs: its arms' `S` come from those runs, and a
different run of the same binary can read several per cent apart on
this host.

| round | candidate | would-be class | decision | cold-time ratio over its base [95%], at the sizes it targets | elsewhere | record |
|:--|:--|:--|:--|:--|:--|:--|
| R02, the AVX-512 batched addition for `n + deg t = 66` | `a1645ac6` on v0′ `c1a2e5f8` | engineering | **rejected**: the holdouts at `2^44.5` read 1.161 [1.089, 1.239], and the lower end is not above 1.10; the callgrind control also failed, by 0.1% of instructions in the compiler's code for two scan functions | suite 1.182 [1.151, 1.213] at `2^44.5`, 1.267 [1.242, 1.292] at `2^47.2`; holdouts 1.161 and 1.262 [1.216, 1.310] | 0.973–1.027, inside every A/A band | [`rounds/R02-wide-tail-kernel/`](../../ic_tool_program/rounds/R02-wide-tail-kernel/README.md), `candidate.patch` |

## 8. Track A: speed

**What §23 already measured.** Under the single-target rule the cold
cost is the set-up: it is 62–292 times the online interval (§23). §23's
reports split the set-up by phase. These figures are at `0bf67f16`,
one thread, isolated, with medians over 64 targets and three
repetitions:

| curve | log₂ r | set-up | collection | build | other | ns per summand scanned | stored pairs |
|:--|--:|--:|--:|--:|:--|--:|--:|
| `icv1-f2m47-t22705043-f4e44623` | 36.6 | 37 ms | 42% | 25% | selection 12% | 39 | 0.15 M |
| `icv1-f2m57-tm747311035-c1f545af` | 38.0 | 129 ms | 31% | 14% | curve construction 41% | 44 | 0.24 M |
| `icv1-f2m41-tm2308219-7f48b14a` | 39.0 | 100 ms | 67% | 18% | selection 6% | 33 | 0.27 M |
| `icv1-f2m53-tm56619371-dac20a85` | 44.3 | 1.03 s | 65% | 28% | selection 3% | 49 | 3.7 M |
| `icv1-f2m59-tm943548413-98844ecc` | 44.5 | 1.37 s | 71% | 22% | selection 3% | 71 | 3.9 M |
| `icv1-f2m61-t158598901-ab42b6c5` | 47.2 | 5.74 s | 89% | 9% | selection 1% | 69 | 6.3 M |

**The collection scan dominates at the top.**
- **Its cost per summand is not set by the table's size alone.** A
  scanned summand costs 33–49 ns (1.4–1.6 units) at the four sizes up
  to 3.7 M stored pairs, and 69–71 ns (2.3–2.4 units) at the two
  largest. 3.7 M pairs cost 49 ns, but 3.9 M pairs cost 71 ns.
- **It already prefetches.** The scan keys a block of 1,024 rests at a
  time. It prefetches the presence filter 32 keys ahead, and the bucket
  of each admitted key before probing it (KIC:5244-5260).
- **Cloud memory is slow.** On the 4-vCPU Xeon VM of an earlier study,
  which was not this host, a random access beyond 2–4 MiB cost
  90–130 ns, and 220–245 ns over 256 MiB with 4 KiB pages
  (`RESEARCH_CACHEGRIND_N41_20260929.md`). R01 measures this host's.
- **The profile names no cause yet.** It is a profile of `0bf67f16`, not
  a measurement of the cause. v0's profile repeats it on main and adds
  instruction and cache-miss counts from callgrind. Each round takes the
  largest phase it can move.

Paths below are to `src/cryptanalysis/koblitz_index_calculus.rs`
(`KIC`) at `0bf67f16`.

| id | lever | why it is worth trying | what could go wrong |
|:--|:--|:--|:--|
| **A2** | **the collection scan** (`collect_walked`, KIC:9635; its scan `witnesses_fast_scan`, KIC:5172). Which lever depends on v0's profile: a longer prefetch distance or a smaller filter if the filter's misses are exposed; huge pages for the table and filter if address translation is (this host's THP is `madvise`, and nothing calls it); the batched addition or the key if the scan is compute-bound | 65–89% of the set-up at `2^39`–`2^47.2`, at 2.3–2.4 units a summand at the two largest sizes against 1.4–1.6 below | huge pages cost first-touch kernel time, which the cold time includes: in that study they cut user time by 14.6% and raised system time sixfold |
| **A1** | **the pair-table build** | 9–28% of the set-up | §19.4's two passes per stored pair survive only in the compact tier (KIC:3822-3828). The folded tier, which every suite size uses, builds in one pass up to 2 GiB of scratch (KIC:4027-4030). So the lever is the one pass's own cost, which needs v0's profile to locate. |
| **A9** | **the curve's construction** | 41% of the set-up at `2^38.0`. `factorise_u64` (KIC:540) tests primality at every trial divisor. | it is shared by both arms, so it moves cold time, not the ratio |
| **A10** | **relation verification and the log system** | each check builds a fresh `FastCurve` inside `kc.mul` (KIC:9791, 828-838). `LogSystem::push` makes a dense big-integer row per relation (KIC:6301-6322). Each extra solve attempt refilters from scratch. | 1–7% of the set-up on the suite, so the gain is mostly at small sizes |
| **A11** | **the scan's canonical key** (`PairSumTable::keys_of`) | 46% of the set-up's instructions at `icv1-f2m61-t158598901-ab42b6c5`, about 565 a scanned summand, and no cache misses (R01's callgrind, added 2026-10-01). Valgrind hides AVX-512, so that is the portable key's cost. The timed runs use the AVX-512 key, about 70 instructions a key by its code (R01's correction). Natively, outside the scan, it costs 5.6–10.2 ns a key, mostly the rotation loop, and the 8-lane batched subtraction 8.3–8.4 ns a summand at `n = 53`: together about half of a scanned summand's ~40 ns ([exploration, 2026-10-01](../../ic_tool_program/explorations/A11-run-kernel-20261001/README.md)). | a vectorised run-based kernel gave the same keys but was slower at `n = 41` and 53 and equal at 61, so it was rejected before declaration (same exploration). R04 then priced the scan's other half inside the scan, with probes compiled into a separate binary (below) |
| **A12** | **re-optimising the recipe for one target** | §20 chose its recipes for `K = 32` targets sharing one set-up; the programme's metric is one target, cold. §20's own frozen model (`predict.py`), re-run at `K = 1` (2026-10-01), predicts 1.00–1.09× at the top six sizes, where descent is under 2% of one target's cost and set-up is collection against build, balanced whatever `K` is. It predicts 2.3–4.7× at `r < 2^25`, where selection dominates. | it changes every output, so it is an algorithm change and is declared as one; it moves only the small sizes, so it waits behind the scan |
| **A3** | online lookups | `target_PDP` is 53–96% of the online interval (§23) | the online interval is under 2% of cold, so this moves the rule's online ratio, not the cold one |
| **A4** | linear algebra: block Wiedemann over `u64` residues (`koblitz_sparse_la.rs`) | under 1% of the set-up at these sizes | it grows faster than the other phases (AGENTS.md §5); price it at every baseline |
| **A5** | rho parity | the batched walk's step costs 1.86–2.25× the canonical step (§23) | it moves the ratio against the index calculus, which is the point: a reference must be as strong as it can be made |
| **A6** | multi-core set-up | four cores here | it must not regress one thread |
| **A7** | field kernels per hardware class | the portable path is unmeasured | each class is reported separately |
| **A8** | an external reference audit | needed before any outside comparison (§1) | published figures are usually stage-only and on other hardware |

**R04 priced the scan's stages (2026-10-01).** It is a stage
diagnostic, class accounting: no code on the timed path changed, and
nothing here is a speedup
([results](../../ic_tool_program/rounds/R04-scan-probes/README.md)).
Time-stamp-counter probes, compiled in only by a `scan-probes` feature,
split the scan into the batched subtraction, the key, the presence
filter and the admitted keys. Both arms were built from `9a48b389` (v0′
plus the probes) and run on `M1`'s 22 rows, three rounds, isolated.
- **The subtraction leads at the top.** It holds 37–41% of the scan at
  three of the four largest sizes. At `2^44.3` the admitted keys lead
  instead, with 42%.
- **The two wide-tail sizes pay for a scalar subtraction:** 32.5 and
  33.3 ns a summand, against 13.1 at `2^44.3`, where the 8-lane kernel
  runs. That is R02's kernel's target, and R02b (#1157) re-tests it.
- **The key and the filter are steady:** 9–16 ns and 4–9 ns a summand
  at the top six sizes. The filter admits 12–21% of keys there.
- **The shares are reported, not trusted.** The protocol trusts a size
  where the probes' overhead interval lies inside R01's A/A band. That
  held at no size: with six pairs a size, each interval is as wide as
  its band. Ten of eleven intervals contain 1, and the ranking rests on
  gaps of 7–18 points of share.
- **A check after the run, not declared:** the filter is one bit per
  stored pair under one hash, at 4–8 bits a pair. From `2^36.6` up, its
  predicted pass rate for absent keys matches the admitted fraction to
  within 2%, and true hits are one summand in 7,200 to 221,000. So the
  admitted stage is spent almost entirely on false positives, at
  118–159 ns an admitted key at the three largest sizes. A sharper
  filter is the next lever: R05, below.

**R03 made v1 (2026-10-01).** `factorise_u64` now tests primality only
when a division changes what is left, and skips even divisors after 2
([results](../../ic_tool_program/rounds/R03-curve-construction/README.md)).
- **At `2^38.0` the construction fell from 56.0 ms to 1.28 ms**, and the
  cold time is 1.771 [1.694, 1.851] times faster on the suite and 1.922
  [1.756, 2.104] on the holdouts. Every pinned output is v0's.
- **At `2^24.8`, the other composite size, nothing moved** (1.007
  [0.985, 1.030]). Its order's cofactor ends at 211, so the old loop
  tested primality only about two hundred times; the prediction of 1.3×
  had assumed `n = 57`'s profile.
- **The prime sizes read 0.970–1.049,** every interval overlapping its
  A/A band.
- **The class is engineering:** the factors are identical, so the ratio
  to the floor is flat. The construction's median is now 2.2 ms or less
  at every suite size.

**R05 tests a sharper presence filter, before R02b (declared
2026-10-01).** Three bits a key in one 64-bit word, at eight to sixteen
bits a stored pair, keep the probe at one word read. At the three
largest sizes they cut the predicted false-positive rate from 0.17–0.21
to 0.019–0.031 ([protocol](../../ic_tool_program/rounds/R05-presence-filter/PROTOCOL.md)).
- **An exploration came first, and is disclosed.** On Track B's stack,
  the prototype ran 1.27× and 1.37× faster in cold time at `2^44.3` and
  `2^47.2`, with identical outputs
  ([record](../../ic_tool_program/explorations/R05-filter-20261001/README.md)).
  It is not a round's measurement.
- **The order changes.** R05 runs after R03's decision and before R02b
  (R02b's amendment 1). Both are on the scan, after R02's rejection
  there, so under §11 a failed R05 sets the scan aside, and R02b does
  not run.

## 9. Track B: generality and robustness

**Where the tool stands (at `0bf67f16`).**
- **The runnable pipelines** (`ic run`, `workflow`, `price`) take Koblitz
  curves `K_a`, and curves defined over `F_{2^k}` (`k ≤ 8`) taken over
  `F_{2^n}` with `n/k` odd. Their field elements are one `u64`.
- **The fast arithmetic stops at `n = 62`** (`FastCurve::MAX_DEGREE`).
  `MAX_N = 63` reaches only the slow paths.
- **Scalars are `u64`**, and the subgroup order is read as its low word
  (KIC:9358-9364, 10569). An `r ≥ 2^64` would be truncated silently.
  No runnable curve reaches that today, since `n ≤ 63` keeps `r` below
  `2^62`. The sparse linear algebra needs `r < 2^63`.
- **Other entry points.**
  - `ic boundary`, `bench` and `rho` add generic binary curves (with
    point counting by enumeration) and prime fields with `p < 2^63`.
  - `ic inspect` checks far larger parameters (binary degree up to 571,
    primes up to 512 bits) but runs nothing.
  - `ic fixed` runs an index calculus on `K_0` up to `n = 131`, in
    Python (`ecc2k130/codegen/indexcalc_fixed.py`).
  - `icx` holds a 52-curve catalog. It runs a small same-family
    substitute, never the named curve.

**Defects the survey found, each a C-suite file once fixed.**
- `r ≥ 2^64` would be truncated rather than refused. It is latent until
  B3.
- `experiment.rs:64-67`'s error text says `n ≤ 62` while 63 passes.
- A foreign `state.json` with a short digest panics on a slice
  (`workflow.rs:1133`).
- A huge `points` overflows in debug builds (KIC:2311).
- Large generic-binary degrees hang in point counting
  (`ic_boundary.rs:3370-3401`), and `ic bench --char2-degree` has no
  range check.
- A too-small `pair_table_bytes` is reported as "field too wide".
- `ic search`'s validation child loses `--subfield`, `--curve-b` and
  `--wdsat-binary` (`experiment.rs:1201-1231`), so a subfield search
  fails its own curve check.
- A report's `commit` is `git rev-parse HEAD`, run in the working
  directory when the report is written (`src/bin/ic.rs:176-188`). A
  binary run from another checkout therefore reports that checkout's
  commit. The harness records the build commit beside it, as §23's host
  manifest did, until B0 embeds it at build time.

| id | step | done when |
|:--|:--|:--|
| **B0** | **Refuse what is wrong today.** Each defect above becomes a precise refusal or a fix, with a C-suite file. `r ≥ 2^63` is refused until B3 lands. | the C suite holds a file per defect and passes |
| **B1** | **Parameter schema v2.** One schema for binary, prime and extension fields, with strict parsing that refuses unknown keys. It has four parts: the field (kind, degree, and modulus or `p`); the curve (general Weierstrass coefficients, or a named family); the subgroup (order, cofactor, generator); and the target and method. | every C-suite file parses or is refused with its expected code |
| **B2** | **Validation.** The modulus is irreducible, the curve is non-singular, the generator is on the curve, `[r]G = O`, `r` is prime by a deterministic test, and the cofactor is consistent. Each failure is a refusal with a code. | every invalid C-suite file is refused with its code; no panic is reachable from input (B6) |
| **B3** | **Two-word binary fields** (`u128`), with multi-word scalars and residues. This lifts the fast path from `n ≤ 62` to `n ≤ 127` and `r` past `2^63`, keeping the one-word kernels' speed at `n ≤ 62` by dispatch. | rho runs at F0 at `n = 83` and the index calculus at F1; the S suite is unchanged |
| **B4** | **Three-word fields** (`n ≤ 191`), for `n = 131` | `E_0/GF(2^131)` at F1 |
| **B5** | **Prime fields beyond 64 bits, and extension fields** `F_{p^k}` | the C suite's files for each, and verified small instances |
| **B6** | **Fuzzing and differential checks.** A seeded generator of parameter files, valid and corrupted, runs in CI with a time budget. Answers are checked against `oracle.py` for `n ≤ 61`. Larger fields are checked against a slow reference implementation; `ic fixed`'s Python is a candidate if an audit shows it shares no arithmetic with the Rust. | no panic and no wrong answer within the budget |
| **B7** | **The F1 sampled mode** | `n = 83` and `131` reported as labelled extrapolations, with their samples |

**The design for B1 and B2** (2026-10-01) is
[`research/ic_tool_program/design/schema-v2.md`](../../ic_tool_program/design/schema-v2.md).
It sets out:
- one schema for every field kind;
- checks with stable codes, each exact or marked as a screen;
- routing to the pipelines that exist, naming the gate that refuses
  each instance and the step that lifts it;
- the conformance cases, of which B1's are frozen in
  `research/ic_tool_program/conformance/v2/`.

**B3 is split in two** (2026-10-01).
- **B3** lifts the rho reference to two-word fields (`n ≤ 126`) and
  adds `solve: rho`. Its measurement is the gate's rho at F0 at
  `n = 83`.
- **B3b** lifts the index calculus to two-word fields, at F0 where the
  subgroup is small enough to afford, declared 2026-10-01
  ([`research/ic_tool_program/rounds/B3b-two-word-kic/PROTOCOL.md`](../../ic_tool_program/rounds/B3b-two-word-kic/PROTOCOL.md)).
  F1 at `n = 83` moves to B7b, which needs its kernels.

B3's declaration is
[`research/ic_tool_program/rounds/B3-two-word-rho/PROTOCOL.md`](../../ic_tool_program/rounds/B3-two-word-rho/PROTOCOL.md).
The row above still states the done-when for both halves together.

**B7 is split in two** (2026-10-01), by its design,
[`research/ic_tool_program/design/f1-sampled.md`](../../ic_tool_program/design/f1-sampled.md).
- **B7a** builds the F1 sampled level where F0 also runs, at one-word
  sizes. It measures F1's error there against F0, size by size, with a
  falsification target declared in the design. It also measures the
  yield constant B7b has to carry.
- **B7b** runs F1 at `n = 83` and `131`, after B3b and B4 give `kic`'s
  kernels two and three words.

**Track B's code is on record before its measurements** (2026-10-01).
- B0 to B3b are implemented, each on the one before, and B4 on B3b;
  none is measured yet. Each step's code reaches `src/` only through its
  results pull request, in the queue's order.
- Until then the stack is kept as a verified git bundle in
  [`research/ic_tool_program/track-b/`](../../ic_tool_program/track-b/README.md),
  so that no local branch is its only record.
- B3b's amendments 1–2, made before any measurement, fix three defects
  in its declaration that the implementation found:
  - the `r ≤ h` refusal is lifted past one word;
  - C031 and C053, whose `kic` width gate B3b moves, are superseded by
    C086 and C087 under a new `supersedes` rule;
  - C082's document name is shortened to the schema's limit.

**The steps' timing checks run as one chain** (an amendment to each B
protocol, 2026-10-01, before any B measurement).
- One interleave covers the newest baseline and each step's arm, in the
  queue's order, on `M1`'s 22 rows, five rounds (`bround.py chain`).
- Each step is read against the arm before it, which is its base, with
  ten pairs a size as before.
- The seven checks take 880 processes, not 1,540.

**B4 is declared, generic to nine words** (2026-10-01, before any B4
measurement; [protocol](../../ic_tool_program/rounds/B4-multi-word/PROTOCOL.md),
[design](../../ic_tool_program/design/multi-word.md)).
- **One more pipeline, generic in the word count:** `rho-koblitz` and
  `kic` on fields of 3 to 9 words, odd `n ≤ 574`. The row above says
  three words; the generic pipeline takes nine at no extra code, so the
  standard curves up to `sect571k1` pass the field gate. Subgroups stay
  below `2^127`, as at two words: ECC2K-130's and the standard curves'
  meet `scalar-wider-than-127-bits`, which B7b lifts for F1 at
  `n = 131`.
- **v1's rules past one word,** which B3 and B3b refused: the hashed and
  random targets and the generator rule, each v1's own at one word.
- **Six F0 instances at three words,** Koblitz curves of prime degree
  127–179 with a 28- to 47-bit subgroup, each solved and verified on a
  known answer and a public point.
- **A correction to the record:** ECC2K-130's published points are in
  the polynomial basis of `x^131 + x^13 + x^2 + x + 1`, in which they lie
  on the curve. Schema v2 §10–11 had deferred a normal-basis import to
  B4; none is needed.
- **Its code is on record** in a second bundle,
  [`track-b/stack-20261001-b4.bundle`](../../ic_tool_program/track-b/README.md),
  on B3b's branch. As development checks, not B4's measurement, all
  fifteen of its frozen cases pass on it. One case's expectation was
  corrected before any measurement: C098's document, at `n = 577`, is
  refused by schema v2's degree range (2..=571) before any gate.

## 10. What does not count

- A speedup with any output changed, unless the round declared an
  algorithm change.
- A ratio inside its A/A band.
- A size dropped from the suite, or a target re-drawn after a run.
- Timings from contended runs.
- A gain measured at one thread count, reported as a gain at another.
- An F1 extrapolation quoted as a measurement.
- "Fastest worldwide" without a matched external comparison.

## 11. Stopping and reporting

- **A phase is set aside** after two consecutive declared rounds on it
  fail their success condition. Track A then moves to the next phase in
  the profile.
- **Every baseline is reported** to the user with its ledger row and
  PR.
- **The programme is open-ended.** It ends when the user ends it.

## 12. The first rounds

- **R01: baseline v0.**
  - Freeze suite v1.
  - Pin main's outputs against §23's at its six sizes.
  - Run the A/A, and the S suite once at one thread.
  - Profile every phase: in-process phase times at all sizes, and
    callgrind at two sizes.
  - Its rule comparison is §23 if main's binary reproduces §23's
    outputs, and a re-run otherwise.
- **R02:** the largest set-up phase in v0's profile. That is A2, the
  collection scan, if v0 repeats §23's profile. Its lever is the one
  callgrind and a huge-page probe point to.
- **R03:** the next phase: A1 or A9, by v0's profile.
- **B0 first, then B1–B2:** the defects and refusals, then the schema
  and its validation, as robustness rounds alongside Track A. B0 comes
  first because it is cheap, and because lifting the field limit (B3)
  would turn the latent truncation of `r` into wrong answers.

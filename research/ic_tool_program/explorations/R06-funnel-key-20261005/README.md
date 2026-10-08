# R06 exploration: the scan's key by funnel shifts, and an invariant filter set aside

**What this is.** These are stage timings and diagnostics for the plan's A11, the scan's canonical key
(`research/notes/index-calculus/IC_TOOL_PROGRAM.md` §8). They ran on
2026-10-02 and 2026-10-05, before R06 was declared, to decide whether to
declare it and what to predict. It is a record, not a round's
measurement. R06's protocol discloses it, and none of its runs is pooled
with R06's.

## The candidate

`least_rotation16_avx512` takes the least rotation of each key's normal
coordinates by `n − 1` chained steps: two shifts, a mask and a minimum
on eight keys, each rotation made from the one before.

The candidate, `least_rotation16_funnel_avx512` (`funnel.patch`), takes
every rotation from the same two words:
- **The doubled word.** It writes the word twice, top-aligned, as
  `hi:lo`.
- **One shift a rotation.** One `vpshldq` (AVX-512 VBMI2) forms
  rotation `t`'s window. The minimum runs over windows compared
  top-aligned.
- **The bits below the window** are never masked: two windows that tie
  on their top `n` bits are equal rotations, whichever wins.
- **The cost.** One shift and one minimum a rotation, where the chained
  kernel takes two shifts, a mask and a minimum, in a dependent chain.
- **The control.** `KIC_CANON_KERNEL=chain` keeps the chained kernel as
  a same-binary control.

**The keys are identical.** The patch's test,
`the_funnel_kernel_matches_the_chained_one`, compares the two kernels
and the scalar search for every `n` from 2 to 63. It covers random,
sparse, dense and periodic words, mixed within blocks of sixteen. It
passes on main (`995ea207`) with the other 29 `koblitz_fast` tests
(`r06m-tests.log`, 30 passed).

## 1. The kernels alone (`kbench/`)

**What ran.** `kbench 41 53 59 61` times each kernel on 2^20 random
coordinate words, best of seven passes, in ns a key. It ran twice,
isolated (`isolated_bench run --cpus 2`, both uncontended). The
checksums of every key agree between the chained and funnel kernels.

| `n` | chained | funnel | invariant, 10 terms | 12 terms | 14 terms |
|--:|--:|--:|--:|--:|--:|
| 41 | 5.43, 7.05 | 2.73, 2.38 | 1.74, 1.63 | 2.06, 1.96 | 2.00, 2.03 |
| 53 | 5.75, 5.87 | 2.52, 2.45 | 1.58, 1.51 | 1.95, 1.83 | 1.96, 1.85 |
| 59 | 7.11, 6.12 | 4.42, 2.65 | 1.51, 1.45 | 1.83, 1.98 | 1.86, 1.98 |
| 61 | 6.12, 6.58 | 2.72, 2.83 | 1.48, 1.60 | 1.84, 1.78 | 1.86, 1.77 |

**The funnel kernel is 2.2–2.5× the chained one,** except one 4.42 at
`n = 59` in the first pass.

## 2. An invariant filter hash, considered and set aside

**The idea.** The presence filter needs only a hash that is the same for
every rotation of the coordinates. The exact least rotation would then
be computed only for the keys the filter admits, 2–3% of summands on v2.
- **The invariant.** The popcount, with the popcounts of
  `c AND rotl^k(c)` for `k = 1..K`, is such a hash.
- **Its cost.** In AVX-512 (VBMI2 and VPOPCNTDQ) it costs about 1.5–2.0
  ns a key at `K = 10–14` (the table above).

**The catch is collisions.** Two orbits with the same invariant pass the
filter together. `invsim` (`kbench/src/main.rs`, `invsim.out`) stores 6.3
M random orbits, the stored-pair count at
`icv1-f2m61-t158598901-ab42b6c5`. It then probes 4 M random words and
counts those whose invariant matches a stored one while their orbit is
not stored:

| `n` | 4 terms | 6 | 8 | 10 | 8, spread `k` |
|--:|--:|--:|--:|--:|--:|
| 61 | 0.991 | 0.761 | 0.105 | 0.0031 | 0.098 |
| 53 | 0.994 | 0.834 | 0.185 | 0.0074 | 0.172 |

**So 12 or more terms are needed** to keep false matches well under the
filter's own 2–3%.

**That leaves the invariant only 0.8–1 ns a key cheaper than the funnel
kernel.** It would still change what the filter means:
- every builder;
- the GPU fold kernel;
- the external table format;
- a second lookup for admitted keys.

**Set aside.** The funnel kernel takes most of the gain and changes
nothing but the arithmetic.

## 3. Whole processes, on v2 (stopped)

**What ran.** `M1`'s two rows at `icv1-f2m53-tm56619371-dac20a85`,
`icv1-f2m59-tm943548413-98844ecc` and `icv1-f2m61-t158598901-ab42b6c5`,
ABAB, isolated, base v2 against v2 plus the patch (`377fcde7`).
- It ran on 2026-10-05 between R02b's timed processes.
- It stopped after 16 of 36 processes, when R02b stopped.
- The runs are in `runs-v2.tar.xz`. Every recovered logarithm matched
  between the arms.

**Cold time, base over candidate, round 1** (one pair a row):

| curve | T01 | T02 |
|:--|--:|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 1.085 | 1.037 |
| `icv1-f2m59-tm943548413-98844ecc` | 0.984 | 1.067 |
| `icv1-f2m61-t158598901-ab42b6c5` | 1.048 | 1.075 |

**Earlier diagnostics** ran on 2026-10-02 on v2, unisolated, with
`taskset` only: three repetitions at `icv1-f2m61-t158598901-ab42b6c5`,
same binary, `KIC_CANON_KERNEL` switching the kernel.
- Collection fell from 3.74–4.12 s to 3.58–3.75 s.
- Set-up fell from 4.61 to 4.26 s at the median.

## 4. Whole processes, on main (`995ea207`)

**What ran.**
- The same rows, three rounds, ABAB, isolated.
- Base main's head, against main's head plus the patch (`878775fc`).
- 36 processes, all uncontended, interleaved with R07's A/A under the
  benchmark lock, on 2026-10-05.
- Every logarithm matched.
- The runs, `explore.sh` and the pairs are in `runs-main.tar.xz`.

**Cold time, base over candidate, geometric mean of six pairs:**

| curve | cold | collection | per-pair spread (sd of log) |
|:--|--:|--:|--:|
| `icv1-f2m53-tm56619371-dac20a85` | 1.015 | 1.036 | 0.048 |
| `icv1-f2m59-tm943548413-98844ecc` | 1.086 | 1.099 | 0.070 |
| `icv1-f2m61-t158598901-ab42b6c5` | 1.068 | 1.088 | 0.052 |

**What it sets in R06's protocol:**
- the prediction, about 1.07–1.09× at the two wide sizes and about
  1.02× at `2^44.3`;
- the power estimate: a spread of 0.05–0.07 a pair gives a 95%
  half-width of about 1.6–2.2% at forty pairs.

## Files

- `funnel.patch`: the candidate, `git am` on main at `995ea207`. It
  applies unchanged to v2's `edcb0bec`.
- `kbench/`: the two benches, a crate of their own
  (`cargo build --release`), and their outputs and isolation records.
- `runs-v2.tar.xz`, `runs-main.tar.xz`: the process outputs, each with
  its isolation record, and the scripts that ran them.
- `SHA256SUMS`, `binaries.sha256`: the files and the binaries.

The host is the one in R07's manifest: `6.18.44-fc-v70`, the 4-vCPU
Xeon with AVX-512, VBMI2, VPOPCNTDQ and GFNI.

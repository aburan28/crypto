# Protocol: the P-192 auxiliary-curve search by 2-Sylow halving (pre-registered)

Written and committed before any measured cell below ran. Code at the commit
that adds this file (branch `claude/isogeny-tower-p192`); the solver is
`pkm-aux-search` (`src/bin/pkm_aux_search.rs`), the library is
`src/cryptanalysis/two_sylow_search.rs`, the component timer
`examples/two_sylow_bench.rs`, the tower check `examples/p192_isogeny_tower.rs`,
the independent checker `backlog_check_aux2t`.

## Question

The PKM isogeny-tower factor base at the NIST P-192 prime
(`p = 2^192 − 2^64 − 1`, `v₂(p − 1) = 1`, `v₂(p + 1) = 64`) needs an auxiliary
curve `E'/F_p` with a rational point of order `2^t`. How many curves must a
wide search draw, what does a draw cost at 192 bits, and what does that make
the `t = 48` search of `backlog/XP-20261015-p192-isogeny-aux-search-t48`
cost?

## Method (fixed)

- **Families.** `z2xz2`: `y² = x(x − α)(x − β)`, `α, β` uniform. `z2xz8`: the
  `X₁(2, 8)` family `Y² = X(X + B⁴)(X + A⁴)`, `A = m² − 1`, `B = 2m`, `m`
  uniform, which has `b ≥ 3` for every member.
- **Seeds.** Draw `s` of a family uses SplitMix64 seeded with
  `s ⊕ tag` (`tag = 0x5A325A3200000001` for `z2xz2`, `0x5A325A3800000001` for
  `z2xz8`), one output discarded, then three outputs per field element as a
  little-endian 192-bit integer reduced mod `p`.
- **Instrument.** For every draw the search computes `b`, where the 2-Sylow
  subgroup is `Z/2 × Z/2^b` (`p ≡ 3 mod 4` rules out `Z/4 × Z/4`), by x-only
  point halving. A point of order `2^t` exists iff `b ≥ t`. Correctness gate
  (unit tests, run before this protocol): `b` equals the exhaustive group
  structure on about 300 curves over primes `103…359`, in all three arms;
  the `z2xz8` family forces `b ≥ 3`; the supersingular `y² = x³ − x` at P-192
  gives `b = 63 = v₂(p + 1) − 1` in all arms; found points certify over
  `BigUint`.
- **Arms.** `baseline` (three `PrimeField::sqrt` per halving, the four halves
  tested one by one with two Legendre symbols each), `fast-chi` (one Legendre
  symbol per level test, two `(p − 3)/4` addition-chain exponentiations per
  halving), `fast` (the next level's exponentiation is the test). All arms
  give the same `b` on every draw (tested).

## Cells

| cell | what | command (all `taskset -c 2`, `nice -n 0`, `prlimit --as=4294967296`) |
|:--|:--|:--|
| C1 | µs per draw, 3 arms × 2 families | `pkm-aux-search --prime nist-p192 --t 64 --family F --arm A --seed-lo 1000000000 --seed-hi 1000065536 --stop none --summary …`; 7 repetitions, arms interleaved within each repetition, order rotated per repetition |
| C2 | primitive timings | `two_sylow_bench --reps 11 --json …` (mul, sqr, isqrt chain, `sqrt`, binary Jacobi, posdivsteps Jacobi, inversion, draw) |
| C3 | depth distribution, `z2xz2` | `--arm fast --t 64 --seed-lo 0 --seed-hi 33554432 --stop none --record-min 16` (`2^25` draws) |
| C4 | depth distribution, `z2xz8` (the ~1 core-hour run) | `--arm fast --t 64 --seed-lo 0 --seed-hi 268435456 --stop none --record-min 16` (`2^28` draws) |
| C5 | artifacts | for `t ∈ {16, 24, 28}` and the largest `b` C4 finds: rerun the recorded seed with `--t t --seed-lo s --seed-hi s+1 --out …`, then `backlog_check_aux2t … --expect-t t` |
| C6 | supersingular control | `--family supersingular-j1728 --t 48` and `--t 64`, checker on both |
| C7 | the tower | `p192_isogeny_tower ARTIFACT T` for `T ∈ {8, 12, 16}` on the C4 curve and `T = 8, 16` on the C6 curve: `|V| = 2^T`, `L(x) = c` on all of `V` (field backend and `BigUint`), no hit on 4096 random abscissae outside `V` |

The host is shared with another lane (4 cores); timings are labelled
noisy. Development timings (40,000 draws, one run each, before this
protocol) were: `z2xz8` baseline 32.5, fast-chi 12.8, fast 9.9 µs/draw;
`z2xz2` 25.5, 17.2, 14.2 µs/draw. The thresholds below were set after
seeing those, which is disclosed.

## Predictions (Haar model of the 2-adic Frobenius, written before data)

With `π = 1 + 2M` (full 2-torsion), `det π = p ≡ −1 mod 2^64`, Haar measure on
the fibre gives `P(b ≥ t) = 3·2^(−t)` for `t ≥ 2` on random full-2-torsion
curves, and `P(b = 1) = 1/4`. The `z2xz8` family conditions on `b ≥ 3`:
`P(b ≥ t) = 2^(3 − t)`.

- **H1** (`z2xz2`, C3): for every `t ≥ 2` whose predicted count is at least 20,
  the observed count of draws with `b ≥ t` lies in the two-sided 99.9 %
  Poisson interval around `N·3·2^(−t)`; and `P(b = 1)` within that interval of
  `N/4`.
- **H2** (`z2xz8`, C4): the same with `N·2^(3−t)`, `t ≥ 4`.
- **Rejection.** A `t` outside its interval rejects the model at that `t`; the
  report then uses the measured fraction, not the model, for the extrapolation.

These are expectations, not certainties: no theorem covers `P-192`'s single
prime (the model is the equidistribution heuristic). The earlier card
hypothesis was `C·2^(t − k)` with `k = 2` (`z2xz2`) and `k = 4` (`Z/2 × Z/8`);
H1 and H2 replace it with `2^t/3` and `2^(t−3)`.

## Speed claims (pre-registered thresholds, C1 medians over 7 repetitions)

- **S1** `z2xz8`: `fast ≤ 0.40 × baseline`.
- **S2** `z2xz8`: `fast ≤ 0.85 × fast-chi`.
- **S3** `z2xz8`: `fast ≤ 12 µs/draw`.
- **S4** `z2xz2`: `fast ≤ 0.65 × baseline`.
- A claim passes only if it holds for the medians; the report also gives the
  min–max spread and says when spreads overlap.

## Derived quantity (labelled derived in the report)

`t = 48` cost on `z2xz8` = `2^45` draws (H2, or `1/P̂` extrapolated if H2 is
rejected) × the C1 `fast` median µs/draw, in core-hours, with the range from
the C1 min and max. Draws to the first success at `t = 16` and `24` are read
off C4 (mean gap between consecutive draws with `b ≥ t`, with its standard
error) and compared with `2^(t − 3)`.

## What this does not decide

Whether a long tower helps PKM at all (lead 1, the solving degree). Whether a
supersingular `E'` (C6) is acceptable to PKM's relation phase is not tested
here; C6 only shows a point of order `2^48` and `2^64` exists on one and the
checker accepts it.

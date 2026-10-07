# A strong-rho reference on the m = 83 gate curve

Registered 2026-10-05, before either part below ran. Depends on
aburan28/crypto#1360 (wide Koblitz curves in `ecbench`).

## Question

What does the strong single-target rho cost, in `S = group-addition
equivalents / √r`, on `icv1-f2m83-tm6151469093347-debefd74`? That is the m = 83
confidence-gate curve of AGENTS.md §8a: `E_0` over `GF(2^83)`, with `r ≈ 2^81`
and `A = 166`. Every IC candidate that claims to transfer toward ECC2K-130
must be read against this figure. None has been measured there.

## Why the reference needs a variant first

`rho.signed_frobenius_strong` cannot run at m = 83 as committed:

- **Memory.** At its `dp_bits = 4` it stores about 6.7 % of all steps (the
  measured `table_entries / walk_steps` on every curve in
  `ecbench_pair_claw_20261003`). An m = 83 solve takes `√(πr/4n) ≈ 2^37.1`
  steps, which puts about `9·10^9` entries in the table.
- **Fruitless cycles.** Raising `dp_bits` to shrink the table lengthens
  the walks. The reference abandons a walk that meets a fruitless cycle,
  and the same records show about one such cycle every 3,300 steps (0.4 to
  0.5 % of 15-step walks). At `dp_bits = 11`, about 38 % of walks would die
  mid-way. That inflates the reference and flatters every candidate read
  against it.

`rho.signed_frobenius_strong_escape` keeps the reference unchanged except
for one thing. On a fruitless cycle the state becomes the canonical form of
`[2]P_next`, with `a` and `b` doubled, and the walk continues. The escape is
a function of the walk's own state, so walks that have merged escape alike,
and it is charged as one group operation. With the escape off, the code is
the reference bit for bit. 160 of 160 committed reference replays were
re-checked identical after the change.

## Boundaries

- **Floor:** `√(π/2A) = √(π/332) = 0.0973`.
- **Reference behaviour so far:** the reference measured 0.81 to 1.35 × the
  floor on `icv1-f2m43` through `icv1-f2m61` in
  `ecbench_pair_claw_20261003`.

## Part A: calibration (local, before the cloud run)

[`spec-calibration.json`](spec-calibration.json). Arms:

- the reference (`dp_bits 4`);
- the escape variant at `dp_bits` 4, 8 and 10.

Curves: `icv1-f2m53-tm56619371-dac20a85`, `icv1-f2m59-tm943548413-98844ecc`
and `icv1-f2m61-t158598901-ab42b6c5`. Eight public targets each, three
interleaved rounds.

- **A1.** Every run verifies.
- **A2.** `escape-dp4`: the mean `S` per curve is within `[0.9, 1.1]` × the
  reference. At 15-step walks the two differ only on the rare fruitless
  cycle.
- **A3.** `escape-dp8` and `escape-dp10`: the mean `S` per curve is within
  `[0.85, 1.25]` × the reference. The in-flight overhead
  `lanes · 2^dp_bits` is at most 3 % of `√(πr/4n)` at m = 61.
- **A4.** `table_entries / walk_steps` is within 25 % of `2^−dp_bits` for
  every escape arm, and `fruitless_walks` is 0.

**Stop condition:** if A1 or A3 fails, the m = 83 run does not launch.

## Part B: the m = 83 run (AWS)

- **Specs:** [`specs/target-1.json`](specs) … `target-8.json`. Each is one
  public target, with target seeds 2026100501 to 2026100508 and the method
  `rho.signed_frobenius_strong_escape` at `dp_bits = 12`. That puts about
  `2^25` table entries per process, with an in-flight overhead of
  `32 · 2^12 ≈ 2^17` steps, about `10^-4` of a solve.
- **Host:** one AWS `r7i.4xlarge` (8 physical cores, 128 GB, us-west-2, key
  pair `meow34` per AGENTS.md §9). Eight sessions run side by side, each
  pinned to its own physical core with its own lock. Operation counts do
  not depend on contention. Wall time is descriptive.
- **Binary:** `ecbench` is built on the host from this branch's commit, and
  the commit and binary hash are recorded.
- **Unattended safety:** the instance's shutdown behaviour is `stop`. A
  watcher stops it when all eight sessions end, and it stops itself after
  48 hours regardless.

Predictions:

- **B1.** Every target verifies as `[d]G = Q`.
- **B2.** The pooled mean `S` over the eight targets is within
  `[0.6, 1.5]` × the floor. A single rho run is roughly Rayleigh-distributed
  (coefficient of variation ≈ 0.52), so a mean of eight has a CV of about
  0.18.
- **B3.** No process exceeds 12 GB of resident memory.

## Inadmissible

- Changing `dp_bits`, the targets or the method after launch.
- Dropping a slow or failed target.
- Quoting wall time as the result.
- Reading part B without part A having passed.

Failures, timeouts and stops are committed as they are.

## Amendment 1 (2026-10-05, before part B launched)

AWS refused the launch: the account is on an account-verification hold
(`RunInstances`: "This account is currently blocked"). Part B therefore runs
on the development host instead. That host is an Apple M4 Pro with 10
performance cores, 4 efficiency cores and 48 GB of RAM, running macOS.

- **Host.** Eight sessions run side by side under `caffeinate`, launched by
  [`run-local.sh`](run-local.sh). macOS has no CPU pinning, so every run is
  isolation L0. Operation counts are the result, as before. Wall time is
  descriptive and contended.
- **`dp_bits` 13, not 12**, so that eight tables fit in 48 GB. That is
  about `2^24` entries per process. The in-flight overhead,
  `32 · 2^13 = 2^18` steps, is still about `2·10^-6` of a solve. Part A
  calibrated `dp_bits` 4, 8 and 10. The table scales as measured
  (`2^−dp_bits` to within 7 %), and nothing in the escape depends on walk
  length.
- **B3 becomes:** no process exceeds 6 GB of resident memory.
- **No unattended stop.** There is no instance to stop. The sessions end
  when their targets are solved, or at the spec's 96-hour timeout.

Nothing else changes: the targets, the method, B1 and B2.

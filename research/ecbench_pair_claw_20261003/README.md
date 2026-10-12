# The pair claw of cryptanalysis#175 against strong rho, on the same targets

**Result.** Ported natively into `ecbench` and run cold against
`rho.signed_frobenius_strong` on the same 48 one-target workloads, the
four-point signed-Frobenius pair claw of
[aburan28/cryptanalysis#175](https://github.com/aburan28/cryptanalysis/pull/175)
costs **2.20× [1.90, 2.55] strong rho at its cost-minimising shape** and
**32.6× [25.3, 41.0] at the shape the PR ran at n = 83**, in
`S = group-addition equivalents / √r`, pooled over six Koblitz curves
from `r = 2^32.1` to `2^47.2`. The ratio does not shrink with `r`. All
960 runs verified, and 12 of 12 audit replays were identical. All five
preregistered predictions pass.

That is what a generic algorithm should do. Every factor-base log is
known, so the claw is a randomised baby-step giant-step on signed
Frobenius classes. Its cost is `(c + 1/c)/√n` against rho's `√(π/4n)`.
It is not index calculus in the sense that could beat rho, and nothing
here is a speedup.

| | |
|---|---|
| protocol | [`PROTOCOL.md`](PROTOCOL.md), committed before the session ran |
| spec | [`spec.json`](spec.json) |
| method | `claw.pair_table`, [`src/cryptanalysis/ecbench/claw.rs`](../../src/cryptanalysis/ecbench/claw.rs) |
| session | [`sessions/koblitz`](sessions/koblitz) `ECBS1ha960fa84a3c2`, 960 records (240 warm-up + 720 measured), status `complete` |
| audit | [`audit-koblitz.json`](audit-koblitz.json): `ok`, 12 of 12 replays `identical` |
| table | [`table-koblitz.md`](table-koblitz.md), written by `ecbench table --reference rho-strong` |
| comparisons | [`sessions/koblitz/comparisons/`](sessions/koblitz/comparisons), one `ecbench compare --save` per arm against `rho-strong` |
| log | [`run-koblitz.log`](run-koblitz.log); [`launch-refused-lock.log`](launch-refused-lock.log) is a first launch the bench lock refused before anything ran |
| host | macOS development host (Apple silicon), `--cpus none`: isolation L0 on every run, so **operation counts are the result and wall time is descriptive** |

## The one table

Mean `S` over 24 verified runs per cell (8 targets × 3 rounds, a fresh
algorithm seed each round), with the ratio to the session's
`rho-strong`. The floor is `√(π/4n)`.

| arm | f2m43 (`2^32.1`) | f2m47 (`2^36.6`) | f2m41 (`2^39.0`) | f2m53 (`2^44.3`) | f2m59 (`2^44.5`) | f2m61 (`2^47.2`) |
|---|---:|---:|---:|---:|---:|---:|
| `rho-strong` (reference) | S 0.183, floor×1.35 | S 0.164, floor×1.27 | S 0.138, floor×1.00 | S 0.099, floor×0.81 | S 0.114, floor×0.99 | S 0.127, floor×1.12 |
| `claw-c1` (`c = 1`) | S 0.338, **ref 1.85** | S 0.338, **ref 2.07** | S 0.347, **ref 2.51** | S 0.293, **ref 2.96** | S 0.236, **ref 2.07** | S 0.258, **ref 2.04** |
| `claw-pr` (PR's n=83 shape) | S 4.79, **ref 26.2** | S 4.52, **ref 27.6** | S 5.77, **ref 41.7** | S 3.38, **ref 34.2** | S 4.80, **ref 42.2** | S 3.64, **ref 28.7** |
| `bsgs-neg` (no Frobenius fold) | S 0.910, ref 4.98 | S 1.021, ref 6.24 | S 0.915, ref 6.61 | S 1.056, ref 10.68 | S 0.937, ref 8.24 | S 0.904, ref 7.15 |
| `rho-strong-aa` (control) | ref 1.000 | ref 1.000 | ref 1.000 | ref 1.000 | ref 1.000 | ref 1.000 |

Curve slugs: `icv1-f2m43-tm998717-e2e742b0`, `icv1-f2m47-t22705043-f4e44623`,
`icv1-f2m41-tm2308219-7f48b14a`, `icv1-f2m53-tm56619371-dac20a85`,
`icv1-f2m59-tm943548413-98844ecc`, `icv1-f2m61-t158598901-ab42b6c5`.
Unlike the r ≤ 2^21 sessions of `ecbench_all_candidates_20261003`, the
reference sits at 0.81–1.35× its floor here, so the ratios do not
flatter the candidate.

Pooled over six curves, 144 pairs each, from the saved comparisons:

| B / `rho-strong` | ratio | 95 % (cluster bootstrap) |
|---|---:|---|
| `claw-c1` | 2.195 | [1.897, 2.546] |
| `claw-pr` | 32.65 | [25.29, 41.04] |
| `bsgs-neg` | 6.969 | [6.063, 8.024] |
| `rho-strong-aa` | 1.000 | [1.000, 1.000] |

Every row is a lower bound by the README §7 rule: canonicalisations,
table inserts and probes are counted, not charged, on both sides.

## Predictions, scored

- **P1, correctness: pass.** 960 of 960 runs verified, including all 384
  claw runs (warm-up included); no wrong answer and no exhaustion.
- **P2, `claw-c1` within `[0.8, 1.25] × 2/√n`: pass.** Measured over
  expected: 1.108 (n = 43), 1.159 (47), 1.111 (41), 1.067 (53), 0.906
  (59), 1.007 (61).
- **P3, `claw-c1` at least 1.5× rho on every curve: pass.** Per-curve
  ratios 1.85 to 2.96. The falsifier did not fire: on every curve the
  claw's 95 % interval of mean `S` lies wholly above rho's.
- **P4, `claw-pr` at least 10× rho on every curve: pass.** Ratios 26.2
  to 42.2. Its means over `(32 + 1/32)/√n` are 0.98, 0.97, 1.15, 0.77,
  1.15 and 0.89. No run exhausted the query domain in this session; one
  scratch smoke run did, before the protocol (see `PROTOCOL.md`).
- **P5, flat in `r`: pass.** A least-squares fit of `log₂(S·√n)` for
  `claw-c1` against `log₂ r` over the six curves has slope **−0.0153**,
  inside `[−0.05, 0.05]`. Total operations grow as `r^{0.48}`, with rho's
  exponent. Reproduce the fit from the table:

  ```sh
  awk -F'|' '$3 ~ /claw-c1/ {gsub(/ /,"",$6); gsub(/ /,"",$8); n=$5; sub(/.*f2m/,"",n); sub(/-.*/,"",n); x=$6; y=log($8*sqrt(n))/log(2); sx+=x; sy+=y; sxx+=x*x; sxy+=x*y; k++} END {printf "slope %.4f\n", (k*sxy-sx*sy)/(k*sxx-sx*sx)}' research/ecbench_pair_claw_20261003/table-koblitz.md
  ```

## Class

**Accounting** (AGENTS.md §3). The algorithm is the PR's. This session
prices it, for the first time, in the repository's unit against a matched
reference on the same targets. Its ratio to the floor is a constant,
2.0–2.6× at `c = 1` and 28–42× at the PR's shape, so it is not an advance,
and no tuning of `c` can make it one: `(c + 1/c)` has its minimum, 2, at
`c = 1`.

## What this says about n = 83 (extrapolation)

`ecbench` handles only `m ≤ 62`, so n = 83 is **not measured**. With P5's
flat ratio and the protocol's `S` laws, strong rho on
`EC1N83Ckb1h876c2921cb64` (`r ≈ 2^81`) costs about `√(πr/332) ≈ 2^37.1`
additions. The claw costs about `2^38.3` at `c = 1` and about `2^42.3` at
the PR's shape. The PR's own modelled `2^45.5` field calls for completed
receipts is consistent with the latter at a few field operations per
addition. These are extrapolations from the six measured curves, resting
on a slope of −0.015 in `log₂(S·√n)`.

## Wall time (descriptive only)

All runs are L0. `ecbench compare` reports the wall ratio as descriptive:
for `claw-c1` the median is 10.5× rho. On f2m61 the mean solve time was
341 ms for `rho-strong`, 3.78 s for `claw-c1`, 46.2 s for `claw-pr` and
12.0 s for `bsgs-neg`, with peak RSS of 47, 42, 8 and 205 MiB. The claw's
additions are unbatched affine additions with a hash probe each, against
rho's batched-inversion lanes. A batched claw would close part of the
wall gap and none of the `S` gap.

## Reproduce

```sh
cargo build --release --bin ecbench
./target/release/ecbench run --spec research/ecbench_pair_claw_20261003/spec.json --out /tmp/claw-rerun --cpus none --wait
./target/release/ecbench verify --dir research/ecbench_pair_claw_20261003/sessions/koblitz --replay 12 --exit-code
./target/release/ecbench table --dir research/ecbench_pair_claw_20261003/sessions/koblitz --reference rho-strong
./target/release/ecbench compare --dir research/ecbench_pair_claw_20261003/sessions/koblitz --a rho-strong --b claw-c1
```

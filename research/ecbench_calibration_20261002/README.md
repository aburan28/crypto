# ecbench calibration: results

**Decision: the harness's accounting is calibrated for the seven generic
methods tested.** Each method's mean `S` reaches the constant its analysis
gives, on fresh data, at `r ≈ 2^25.4` (prime) and `r = 2^39` (Koblitz).
Two of round 1's preregistered predictions failed as written. Both
failures stand, and a preregistered follow-up on fresh targets traced
them to sample size, not to the harness.

This is a validation of the measurement instrument, not an attack
result. No method was changed and no ratio below is a speedup.
AGENTS.md §3's classes do not apply, and the IC scoreboard is untouched
because no IC variant was measured.

| | |
|---|---|
| protocols | [`PROTOCOL.md`](PROTOCOL.md) (round 1), [`PROTOCOL-2.md`](PROTOCOL-2.md) (follow-up), each committed before its sessions ran |
| binary | `ecbench`, release, SHA-256 `32390c09f730971ad83eeff8ff97262be1150dc5e4c90657784f798827ef32af`, built from `78a2b631f`; the follow-up ran it from `92a6a219c` (only test code changed) |
| host | Apple M4 Pro, macOS 26 (Darwin 25.6.0), arm64, class `ECBENV1h9391b7a712c6`; no affinity control, so every run is **L0**; operation counts only |
| sessions | [`sessions/prime`](sessions/prime) `ECBS1h162089ae4d82` (1 008 runs), [`sessions/koblitz`](sessions/koblitz) `ECBS1h73210d525d59` (864), [`sessions/followup-prime`](sessions/followup-prime) `ECBS1he0789e5685e0` (256), [`sessions/followup-koblitz`](sessions/followup-koblitz) `ECBS1h6f9ddee85cae` (256) |
| correctness | 2 384 of 2 384 executions verified (`[k]G = Q`, checked by the runner and again by the audit) |
| replay certificates | audit receipts, each `ok` with 12 of 12 replays identical: `audit-prime.json` `f8d850c7…0842`, `audit-koblitz.json` `6eba1d5c…08a4`, `audit-followup-prime.json` `ab52f59d…7467`, `audit-followup-koblitz.json` `8fbb9b97…0048` |
| full tables | [`table-round1.md`](table-round1.md), [`table-followup.md`](table-followup.md), both written by `ecbench table` |

## Round 1: the predictions, scored

At each family's largest curve: `icv1-fp26-tm1775-7e8fb6df`
(`r = 2^25.35`) and `icv1-f2m41-tm2308219-7f48b14a` (`r = 2^39`). There
were 8 targets × 3 rounds per cell. Intervals are 95 % two-stage
bootstrap intervals of the mean `S` (workloads, then runs).

| prediction | method | curve | mean `S` | 95 % interval | theory | `S/theory` | verdict |
|---|---|---|---:|---|---:|---:|---|
| P1 | `bsgs.textbook` | fp26 | 1.466 | [1.350, 1.572] | 1.500 | 0.977 | pass |
| P1 | `bsgs.interleaved` | fp26 | 1.396 | [1.197, 1.587] | 1.333 | 1.047 | pass |
| P1 | `bsgs.negation` | fp26 | 0.966 | [0.850, 1.073] | 1.000 | 0.966 | pass |
| P1 | `bsgs.interleaved` | f2m41 | 1.219 | [0.974, 1.466] | 1.333 | 0.914 | pass |
| P1 | `bsgs.negation` | f2m41 | 0.980 | [0.806, 1.163] | 1.000 | 0.980 | pass |
| P2 | `rho.plain` | fp26 | 1.306 | [1.075, 1.561] | 1.253 | 1.042 | pass |
| P2 | `rho.negation` | fp26 | 0.935 | [0.681, 1.198] | 0.886 | 1.055 | pass |
| P2 | `rho.signed_frobenius` | f2m41 | 0.134 | [0.106, 0.166] | 0.138 | 0.970 | pass |
| P2 | `rho.negation` | f2m41 | 0.875 | [0.647, 1.128] | 0.886 | 0.987 | pass |
| P2 | `rho.plain` | f2m41 | 0.993 | [0.688, 1.307] | 1.253 | **0.792** | **fail** (band `[0.85, 1.35]`; interval contains theory) |
| P3 | `kangaroo.vow` | f2m41 | 1.934 | [1.549, 2.323] | 2.000 | 0.967 | pass |
| P3 | `kangaroo.vow` | fp26 | 1.203 | [0.684, 1.761] | 2.000 | **0.601** | **fail** (band `[0.8, 1.5]`; interval excludes theory) |
| P4 | all | all | | | | | pass: every execution verified, every replay identical |

`rho.frozen_reference` was reported, not tested. It runs at
`S/theory = 1.31` on fp26 and 5.1 on fp16, the historical "before" walk
the tuned walks replaced.

## Why the two failed, and the follow-up

The investigation is recorded in full in [`PROTOCOL-2.md`](PROTOCOL-2.md).
In short:

- **Rho:** a rho run's cost varies with a coefficient of variation near
  0.5, so a 24-run mean moves about ±20 % at 95 %. The band was set on a
  point estimate without that spread. The interval contained the theory
  constant throughout.
- **Kangaroo:** `kangaroo.vow` is deterministic given its target. The
  three rounds of each target in `sessions/prime` are identical, so the
  effective sample was 8 targets. Its cost tracks the distance
  `|k − r/2|`: `S` ran from 0.13 at distance `0.017 r` to 2.48 at
  `0.347 r`. Those 8 targets averaged `0.119 r` against a uniform
  `0.25 r`, about −2.6σ. The derivation is unbiased: a native test over
  20 000 draws passes a 1 % Kolmogorov–Smirnov test with mean distance
  `0.2499` (`workload::tests::planted_scalars_are_uniform`).

The follow-up reran exactly those cells on **64 fresh targets** × 2
rounds (`target_seed` 20261003):

| prediction | method | curve | mean `S` | 95 % interval | theory | `S/theory` | verdict |
|---|---|---|---:|---|---:|---:|---|
| F1 | `rho.plain` | f2m41 | 1.240 | [1.109, 1.389] | 1.253 | 0.989 | pass |
| F2 | `kangaroo.vow` | fp26 | 1.972 | [1.759, 2.178] | 2.000 | 0.986 | pass |
| F3 | all | both | | | | | pass: 512 of 512 verified, 24 of 24 replays identical |
| (reference) | `rho.signed_frobenius` | f2m41 | 0.139 | [0.125, 0.154] | 0.138 | 1.006 | |
| (reference) | `rho.negation` | fp26 | 0.991 | [0.884, 1.109] | 0.886 | 1.119 | |

By the interpretation fixed in `PROTOCOL-2.md`, round 1's failures are
attributed to sample size and the calibration is accepted.

## What this establishes, and what it does not

- **Established:** on these curves, `ecbench` charges plain, negation and
  signed-Frobenius rho, the three BSGS forms and the kangaroo what their
  analyses predict, to within the stated intervals. A comparison between
  any two of them in `ecbench` measures the methods, not a bookkeeping
  difference. Every count reproduces bit for bit from its session files.
- **Not established:** wall time (every run is L0); index calculus (it
  has no closed-form constant to calibrate against); behaviour above
  `r = 2^39` or on curves the word-size groups cannot hold; the NUMA
  path (single-node host).
- **Lessons for later protocols.** A tolerance on a mean is set from that
  mean's sampling spread, or stated as "the constant lies in the
  interval". A method that is deterministic per target needs targets,
  not rounds. `kangaroo.vow` could be made seed-dependent by randomising
  its starting offsets; that would be a method change and would need its
  own calibration.

## Reproduce

```bash
cargo build --release --bin ecbench
```

```bash
./target/release/ecbench verify --dir research/ecbench_calibration_20261002/sessions/prime --replay 24 --exit-code
```

```bash
./target/release/ecbench table --dir research/ecbench_calibration_20261002/sessions/prime research/ecbench_calibration_20261002/sessions/koblitz
```

To rerun from scratch, run each spec in `specs/` into a new directory.
Operation counts will match these sessions exactly on any host built from
the same commit, and CI checks that on Linux x86-64 for every committed
session (`.github/workflows/ecbench.yml`).

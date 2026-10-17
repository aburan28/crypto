# Protocol: `two_sylow_search` on `fp_backend::Fp64`, before and after the branch-free selections

Frozen before any run below. The protocol commit is the last commit that
changes this file or the code it tests before the first run; it is recorded
in `runs/protocol_commit.txt`, and both binaries are built from it.

Follows `research/fp64_branchfree_20261016/` (PR #2034). That study made
`fp_backend::Fp64`'s reduction, `add` and `sub` branch-free
(`select_unpredictable`). It measured the rho kernel 2.3–2.5× faster and
dependent `mul`/`sqr` chains unchanged. A dependent chain of `sub`s with a
predictable wrap ran 18–33% slower. Its "Next" item 1 asks whether a
latency-bound production user of `fp_backend::Fp64` relies on such chains.
This study times the first one named: `two_sylow_search`.

## Workload

`two_sylow_search::Halver::run_draw`: one draw of the 2-Sylow halving
search. It draws a curve, then descends levels by square roots or
`(p − 3)/4` powers and Legendre symbols, recording a point of order `2⁶`
when the depth reaches 6 (the module's own tests use `t = 6`).

- **Field and prime.** `fp_backend::Fp64` at the largest prime
  `p ≡ 3 (mod 4)` below 2⁶³ (9223372036854775783). `Fp64` takes the generic
  `HalvingField` defaults: `pow` for `a^((p−3)/4)` and `jacobi_u64` for the
  Legendre symbol.
- **Arms.** `tsy.<family>.<arm>` for both families (`z2xz2`, `z2xz8`) and all
  three search arms (`baseline`: three square roots per level; `fast-chi`;
  `fast`). Each sample runs consecutive draws from seed 0, continuing across
  calls, and reports wall ns per draw.

Production `pkm_aux_search` runs the same search on `FpP192`. This study
measures the `Fp64` instantiation only, which is the one the change
touches.

**Gate (before timing, every run).** Over 2,000 draws per family, the three
search arms must agree on every draw's depth. The run logs the depth sum,
the number of recorded points and an XOR of each arm's recorded `x`. The
baseline and candidate runs must log identical lines, because the
arithmetic is bit for bit the same.

## Arms and runs

- **Binaries.**
  - *Candidate (B):* `cargo build --release --bin hwcost` at the protocol
    commit (branch-free `fp_backend::Fp64`).
  - *Baseline (A):* the same tree with `src/cryptanalysis/fp_backend.rs`
    restored to `main` `f25d5e520`. Only `fp_backend::Fp64` differs.

  Both use the default x86-64 target, and each run records its binary's
  SHA-256.
- **Host.** This session's cloud container: Intel Xeon @ 2.10 GHz, family 6
  model 207, 4 vCPUs, KVM. Isolation is L0. There is no PMU. No other build
  or test runs during the runs.
- **Runs.** Four runs back to back: A1 baseline, B1 candidate, A2 baseline,
  B2 candidate. Each is `hwcost measure --fields fp64,w1 --two-sylow` with
  round-1 parameters (7 repetitions, 100 ms per arm), writing `runs/base-1`,
  `runs/cand-1`, `runs/base-2` and `runs/cand-2`.
- **Controls.**
  - *Positive:* the `fp64` rho kernel, whose 2.3–2.5× gain PR #2034
    measured. It shows the comparison detects a known effect.
  - *Negative:* the `w1` (`FpWide<1>`) kernel, whose code does not change.

`hwcost compare` computes, per arm, `r_k = B_k / A_k` and the drift
`d = |A2/A1 − 1|`. The claims are computed by
`hwcost compare --max-ratio/--max-dev` (Rust).

| id | role | claim | rule |
| --- | --- | --- | --- |
| S1 | the question | no `two_sylow_search` arm gets slower | every `tsy.*` arm: max(r1, r2) ≤ 1.10 |
| K1 | positive control | the known effect is seen | `kernel.fp64` at every B ≥ 4096: max(r1, r2) ≤ 0.85 |
| C1 | negative control | unchanged code does not move | `kernel.w1`, B ≥ 4096: max \|r_k − 1\| ≤ 0.10 |
| G | gate | both binaries compute the same search | the gate lines in `stderr.txt` are identical across the four runs |
| D1 | descriptive | each `tsy` arm's ratio and class (faster / slower / no detected effect beyond max(2d, 2%)) | no threshold |

How to read the result:

- **S1 fails:** the change slows this latency-bound user beyond 10%. The
  result is reported with the arm and size, as a cost of PR #2034.
- **S1 passes:** D1 says whether the search got faster, stayed the same, or
  got slower by less than 10%.
- **K1 fails or C1 fails:** the comparison is reported as unable to see the
  effect, or as confounded, respectively.

## Outputs

- `runs/{base,cand}-{1,2}/`: `measurements.json`, `samples.csv`,
  `medians.csv`, `stderr.txt`, `binary.txt` and timestamps.
- `runs/compare/`, `runs/protocol_commit.txt`.
- A run report and a results report, with graphs and PDFs.
- `benchmarks/20261017-two-sylow-fp64/`.

## Commands

```sh
S=<scratch dir>; D=research/two_sylow_fp64_20261017/runs; F=src/cryptanalysis/fp_backend.rs
cargo build --release --bin hwcost && cp target/release/hwcost $S/hwcost-cand
git show f25d5e520:$F > $F && cargo build --release --bin hwcost && cp target/release/hwcost $S/hwcost-base
git checkout HEAD -- $F
for r in base-1 cand-1 base-2 cand-2; do
  mkdir -p $D/$r; $S/hwcost-${r%-*} measure --fields fp64,w1 --two-sylow --out $D/$r 2> $D/$r/stderr.txt
done
mkdir -p $D/compare
$S/hwcost-cand compare --base $D/base-1,$D/base-2 --cand $D/cand-1,$D/cand-2 --out $D/compare \
  --max-ratio S1:tsy.:0:1.10 --max-ratio K1:kernel.fp64.:4096:0.85 \
  --max-dev C1:kernel.w1.:4096:0.10 > $D/compare/summary.txt
for r in base-1 cand-1 base-2 cand-2; do grep '^two_sylow ' $D/$r/stderr.txt > $D/$r/gate.txt; done
cmp $D/base-1/gate.txt $D/cand-1/gate.txt && cmp $D/base-1/gate.txt $D/base-2/gate.txt \
  && cmp $D/base-1/gate.txt $D/cand-2/gate.txt && echo "G PASS" > $D/compare/gate.txt
```

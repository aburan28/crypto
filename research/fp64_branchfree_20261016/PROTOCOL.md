# Protocol: branch-free selections in production `fp_backend::Fp64`

Frozen before any run below; the commit that adds this file (with the code
change it tests) is the protocol commit. Follows `hwcost` round 2
(`research/hwcost_round2_20261016/`, PR #2003), whose addendum 1 found that
writing `Fp64`'s data-dependent selections with
`std::hint::select_unpredictable` (`hwcost::kernel::Fp64Cm`) cut the rho
kernel to 0.40–0.50× of `Fp64`'s time at every B ≥ 4096 on this host's
baseline-ISA build, bit for bit the same results.

## The change

In `src/cryptanalysis/fp_backend.rs`, `Fp64` (odd `p < 2⁶³`, Montgomery `u64`):

- `redc`: the final conditional subtraction `u ≥ p ? u − p : u`;
- `PrimeField::add`: the final conditional subtraction of `a + b`;
- `PrimeField::sub`: the wrap-around `a ≥ b ? a − b : a + p − b`;

each becomes one `overflowing_sub` and one `select_unpredictable`, exactly
as `Fp64Cm`. `neg` keeps its branch (`a == 0` is almost never true on data,
so the branch predicts). Every `Fp64` user inherits the change (`to_mont`,
`from_mont`, `mul`, `sqr`, `wide_reduce`, `inv`'s final product, the `Fp52`
IFMA backend's scalar operations, F4 rows' reductions).

**Correctness gate (before any timing).** A new unit test
(`fp64_selects_match_branchy_formulas`) compares `redc`, `mul`, `sqr`,
`add`, `double` and `sub` with the previous branchy formulas, bit for bit,
on every pair of an edge set (0, 1, 2, p − 2, p − 1, ⌊p/2⌋, ⌊p/2⌋ + 1),
on 20,000 random pairs, and `redc` on 20,000 random inputs across its whole
range `t < p·2⁶⁴`, for every prime of the existing `SMALL`, `P32` and `P64`
lists; the existing `fp64_matches_biguint`, `f4_tower_generic`
`Fp64`, `Fp52` and `two_sylow_search` tests must pass unchanged; and every
`hwcost measure` run checks its walks against an unbatched reference and
`Fp64Cm` against `Fp64` (round 2's gate). Results are bit for bit the
previous ones, so no recorded output of any other study changes.

## Arms and runs

- **Baseline (A):** `cargo build --release --bin hwcost` at the main commit
  this branch starts from (`f25d5e520`), built before any edit.
  **Candidate (B):** the same command at the protocol commit. Each run
  records its binary's commit and SHA-256 (`binary.txt`). Both
  default target (x86-64 baseline), the repository's production build.
- Host: this session's cloud container (Intel Xeon @ 2.10 GHz, family 6
  model 207, 4 vCPUs, KVM; isolation L0; no PMU, so no branch-miss counts).
- Runs, back to back: **A1** baseline, **B1** candidate, **A2** baseline,
  **B2** candidate, each `hwcost measure --fields fp64,fp64cm,w1` (round-1
  parameters: 7 repetitions, 100 ms, batches 2⁴ … 2²², all arm groups) →
  `runs/base-1`, `runs/cand-1`, `runs/base-2`, `runs/cand-2`, each fitted
  with `hwcost fit`.
- Fields: `fp64` (production `Fp64`, the arm under test); `fp64cm`
  (`hwcost::kernel::Fp64Cm`, whose selections were already branch-free and
  whose code does not change: a control and the reference P3 compares to);
  `w1` (`FpWide<1>`, untouched: a control).

`hwcost compare` computes per arm `r_k = B_k / A_k` and the drift
`d = |A2/A1 − 1|`; the claims below are computed by `hwcost compare
--max-ratio/--max-dev` (added in the protocol commit, with a unit test) and
`hwcost pair` (Rust), not by the report.

| id | role | claim | rule |
| --- | --- | --- | --- |
| P1 | gate | the rho kernel gets faster | `kernel.fp64` at every B ≥ 4096 (6 sizes): max(r1, r2) ≤ 0.85 |
| P2a | guard | a dependent multiply chain does not get slower | `fp64.mul.latency`: max(r1, r2) ≤ 1.10 |
| P2b | guard | a dependent squaring chain does not get slower | `fp64.sqr.latency`: max(r1, r2) ≤ 1.10 |
| P3 | match | production `Fp64` matches the study's `Fp64Cm` | in cand-1 and cand-2, `fp64 / fp64cm` kernel ratio ≤ 1.10 at every B ≥ 4096 |
| C1a | control | unchanged `Fp64Cm` kernel does not move | `kernel.fp64cm`, B ≥ 4096: max \|r_k − 1\| ≤ 0.10 |
| C1b | control | unchanged `FpWide<1>` kernel does not move | `kernel.w1`, B ≥ 4096: max \|r_k − 1\| ≤ 0.10 |
| D1 | descriptive | every arm's class (faster / slower / no detected effect beyond max(2d, 2%)) | counts per group, no threshold |
| D2 | descriptive | the machine code changed as intended | instructions, conditional jumps and `cmov` of `step_fp64` and `step_fp64cm` in both binaries (`objdump -d`, round 2's `f3_disasm.sh`) |

P2 guards the latency-bound users of `Fp64` (exponentiation chains in
`pow`, `sqrt`, Legendre): on a chain a `cmov` adds latency where a
well-predicted branch does not. The change lands only if P1, P2a and P2b
pass; P3 and C1 qualify the reading (a control outside ±10% means the
builds differ in more than the change, and P1 is reported as confounded).
Any failure is reported as measured.

## Outputs

`runs/{base,cand}-{1,2}/` (`measurements.json`, `samples.csv`,
`medians.csv`, `validation.{json,csv}`, `profile.json`, `fit.txt`,
`stderr.txt`, `commit.txt`, `binary.txt`, timestamps), `runs/compare/`,
`runs/cm-pair/`, `runs/disasm/`, a run report and a results report with graphs and PDFs, and
`benchmarks/20261016-fp64-branchfree/`.

## Commands

```sh
S=<scratch dir>; D=research/fp64_branchfree_20261016/runs
git checkout f25d5e520 && cargo build --release --bin hwcost && cp target/release/hwcost $S/hwcost-base
git checkout <protocol commit> && cargo build --release --bin hwcost && cp target/release/hwcost $S/hwcost-cand
for r in base-1 cand-1 base-2 cand-2; do
  mkdir -p $D/$r; $S/hwcost-${r%-*} measure --fields fp64,fp64cm,w1 --out $D/$r 2> $D/$r/stderr.txt
done
for r in base-1 cand-1 base-2 cand-2; do $S/hwcost-cand fit --run $D/$r > $D/$r/fit.txt; done
mkdir -p $D/compare $D/cm-pair
$S/hwcost-cand compare --base $D/base-1,$D/base-2 --cand $D/cand-1,$D/cand-2 --out $D/compare \
  --max-ratio P1:kernel.fp64.:4096:0.85 --max-ratio P2a:fp64.mul.latency:0:1.10 \
  --max-ratio P2b:fp64.sqr.latency:0:1.10 --max-dev C1a:kernel.fp64cm.:4096:0.10 \
  --max-dev C1b:kernel.w1.:4096:0.10 > $D/compare/summary.txt
$S/hwcost-cand pair --runs $D/cand-1,$D/cand-2 --a fp64cm --b fp64 --min-batch 4096 \
  --threshold 1.10 --id P3 --out $D/cm-pair > $D/cm-pair/summary.txt
F3=research/hwcost_round2_20261016/f3_disasm.sh
$F3 $S/hwcost-base baseline $D/disasm; $F3 $S/hwcost-cand candidate $D/disasm
```

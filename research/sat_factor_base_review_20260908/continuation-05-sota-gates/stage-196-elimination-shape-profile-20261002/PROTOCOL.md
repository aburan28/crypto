# Stage 196 protocol: matrix-shape elimination profile

## Purpose and claim boundary

The selected current F4 engine performs 147,794,583,858 table-assisted word
XORs on the public `n=59, ell=9, m=3` true-negative. The existing opt-in
full-matrix M4RI control reduces that count substantially but previously
regressed whole-process wall time. Stages 194 and 195 rejected scheduler and
allocator explanations.

Stage 196 measures where elimination time and work occur by matrix row-count
bin in one selected-`BlockTables` run and one full-M4RI-control run. It does
not select a runtime implementation or report a speedup. Its only decision is
whether the data justify a separately preregistered shape-selective hybrid
screen.

## Frozen source, target, and arms

- Branch base after the inherited Stage 195 and GLV E15 merges:
  `60fb2ae47fcca9352d8bae2b280fce77371e57f0`.
- Input manifest:
  `../stage-175-current-f4-single-target-20261001/input/manifest.json`.
- Manifest SHA-256:
  `188dfdcb7e9398b0ab03193265f7b3035776d4a8a52de7dd4a51475dfc1a572c`.
- Public source-instance id:
  `954e10f8bf0280094fed195280b203d7cd613150b17339716eb469b88ffa9ac7`.
- Equation fingerprint:
  `02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb`.
- Algebraic factor base: `span_F2(1,z,...,z^8)`, with target-subgroup
  enumeration false and discrete-log-label use false.
- Both arms use fixed-X1 direct symmetrised S4, all 512 masks, all 242 rational
  systems, dense exact critical-pair selection, build-serial inner policy, and
  one twelve-worker outer batch.
- Current arm: `F4_F2_FULL_M4RI=0`.
- Full control: `F4_F2_FULL_M4RI=1`.
- Both arms: `F4_F2_MATRIX_PROFILE=1`.

No matrix selector or elimination implementation changes in this stage.

## Profile bins and counters

Every `echelon` call is assigned by matrix row count to exactly one bin:

```text
lt256
r256_511
r512_1023
r1024_2047
r2048_4095
ge4096
```

For each bin, aggregate:

- matrix count;
- row and column sums;
- elimination wall nanoseconds summed across F4 calls;
- row-equivalent logical word XORs;
- actually performed word XORs; and
- number of matrices routed through full M4RI.

The sums across bins must exactly reproduce the run's matrix/elimination
counters. Profiling adds one timer and fixed counter updates per matrix. These
profiled process times are diagnostics and are not compared to unprofiled
Stages 192, 194, or 195 as speed measurements.

## Native execution and accounting

Use the merged Rust `koblitz_f4_stage188 meter` and a Rust Stage 196
composer/verifier. No Python process may execute, verify, compose, or hash new
evidence. Both arms set:

```text
RAYON_NUM_THREADS=12
PQ_F4_X1_BATCH=512
F4_F2_DENSE_PAIR_SELECT=1
F4_F2_DISABLE_INNER_BUILD_PARALLEL=1
F4_F2_MATRIX_PROFILE=1
VECLIB_MAXIMUM_THREADS=1
OPENBLAS_NUM_THREADS=1
OMP_NUM_THREADS=1
MKL_NUM_THREADS=1
BLIS_NUM_THREADS=1
NUMEXPR_NUM_THREADS=1
```

Run exactly:

```text
current-profile, full-m4ri-profile
```

Charge the fresh build, mode-specific tests, both profiled processes, all
workers, memory, failures, composition, and replay. Report process wall,
total core-seconds, peak RSS, conflicts `null`, and single-core `null`.

## Correctness gates

Both arms must authenticate the source and equation fingerprint, use no target
subgroup enumeration or discrete-log labels, visit all 512 masks, complete all
242 rational systems, find zero roots, and return exhaustive `UNSAT`.
Timeout, resource refusal, a malformed profile, or any terminal mismatch makes
the profile unusable.

Full M4RI may choose a different valid pivot basis, so downstream reducer and
logical-XOR counters may differ. Source equations, terminal classification,
factor-base contract, mask/system counts, degree bound, and exact curve result
must agree.

## Frozen selector rule

For each candidate suffix beginning at row thresholds `256`, `512`,
`1024`, `2048`, or `4096`, sum the corresponding bins in both arms.

A threshold may be recommended for a separate hybrid screen only when:

1. both arms have nonzero matrices in the suffix;
2. the full/current profiled elimination-nanosecond ratio is below `0.90`;
3. the full/current performed-XOR ratio is below `0.90`; and
4. the current suffix accounts for at least `0.05` of current profiled
   elimination nanoseconds.

Choose the largest qualifying threshold, which minimizes the scope of the
experimental full-M4RI route. If no threshold qualifies, record
`NO_HYBRID_THRESHOLD`. This stage never changes the repository default.

## Boundaries

Direct MITM, the named SAT controls, unknown-scalar evidence, full rho
comparison, and independent-review gates remain unchanged. Complete campaign
cost remains `null`; this is a one-target implementation profile.

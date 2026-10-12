# Stage 189 protocol: outer-only parallelism for fixed-X1 F4

## Hypothesis

The selected `n=59, ell=9, m=3` Phase B schedule places all 242 rational
fixed-X1 systems in one outer Rayon batch with twelve workers. Each concurrent
F4 solve can also submit parallel product construction, symbolic preprocessing,
row packing, block-table construction, and row-reduction work to the same
global pool. Rayon preserves correctness under nesting, but nested work stealing
may enlarge simultaneously live matrix working sets, disrupt cache locality,
and add scheduling overhead.

The pre-registered candidate retains the twelve-way outer batch but executes
each individual F4 solve serially internally. The hypothesis is that independent
system-level parallelism is sufficient to occupy the pool and that removing
nested parallel sections lowers both wall time and total core-seconds.

This is a scheduling experiment on one opened decomposition target. It is not
a new algorithm, factor base, relation-yield result, DLP, rho comparison,
external reproduction, novelty result, or SOTA claim.

## Frozen target and selected algorithm

- Branch base and inherited Stage 188 merge:
  `66a02ac25ba14355f67046f2084cfd9621a1bade`.
- Input manifest:
  `../stage-175-current-f4-single-target-20261001/input/manifest.json`.
- Input manifest SHA-256:
  `188dfdcb7e9398b0ab03193265f7b3035776d4a8a52de7dd4a51475dfc1a572c`.
- Source-instance id:
  `954e10f8bf0280094fed195280b203d7cd613150b17339716eb469b88ffa9ac7`.
- Equation fingerprint:
  `02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb`.
- Algebraic factor base: `span_F2(1,z,...,z^8)`, without target-subgroup
  enumeration or known discrete-log labels.
- Selected solver: dense exact critical-pair updates, deterministic linear
  symbolic-reducer scan, five-column `BlockTables`, and full M4RI disabled.

No equation, mask order, batch size, pair rule, reducer choice, matrix shape,
elimination policy, extraction rule, or terminal verification may change.

## Native execution and arms

Use the merged Rust-native `koblitz_f4_stage188 meter` command for fresh-process
`wait4` wall/CPU/RSS accounting, output hashing, process-group custody, and
watchdogs. A Stage 189 Rust verifier authenticates the meter receipts and
backend records, checks exact cross-arm structure, computes paired ratios, and
composes the decision. No Python process may execute, verify, compose, or hash
new Stage 189 evidence.

Both arms explicitly set:

```text
RAYON_NUM_THREADS=12
PQ_F4_X1_BATCH=512
F4_F2_DENSE_PAIR_SELECT=1
F4_F2_FULL_M4RI=0
VECLIB_MAXIMUM_THREADS=1
OPENBLAS_NUM_THREADS=1
OMP_NUM_THREADS=1
MKL_NUM_THREADS=1
BLIS_NUM_THREADS=1
NUMEXPR_NUM_THREADS=1
```

- Control: `F4_F2_DISABLE_INNER_PARALLEL=0`.
- Candidate: `F4_F2_DISABLE_INNER_PARALLEL=1`.

Unset remains equivalent to control unless confirmation selects the candidate.
The candidate must report one disabled-inner-parallel policy count for every F4
call while retaining outer parallel execution. Charge builds, tests, every
child process, all twelve workers, failures, verification, composition, wall,
total core-seconds, and peak RSS. `single_core_seconds` remains null because
the outer batch is parallel.

## Correctness and identity gates

1. Existing Boolean-F4, fixed-X1 specialization, and backend tests pass before
   measurement.
2. Every arm authenticates the same source, visits all 512 masks, skips the
   same 270 non-rational masks, constructs and completes all 242 systems, finds
   zero roots, and returns exhaustive `UNSAT`.
3. Target-subgroup enumeration and discrete-log-label use remain false;
   conflicts remain `null` because F4 has no SAT-conflict counter.
4. Arms agree exactly on the equation fingerprint, equations and terms, pair
   and field-pair counters, dense-pair counters, divisor tests, matrix rows and
   columns, basis, extraction, degree, logical XORs, actually performed XORs,
   full-M4RI counters, and peak matrix/table/scratch accounting.
5. The only permitted backend differences are scheduling-policy fields,
   disabled-policy call counts, timing, process RSS, and scheduling-dependent
   internal nanosecond diagnostics.

Any correctness or structural mismatch rejects the candidate without reading
timing. Timeout and resource refusal are censored rather than `UNSAT`.

## Frozen schedule and decision

Run one screen pair:

```text
current, outer-only
```

Continue only if both records are correct and outer-only/current wall and
total-core ratios are below `1.00`.

If the screen continues, run three confirmation pairs in fixed order:

```text
current, outer-only, outer-only, current, current, outer-only
```

Select outer-only scheduling as the repository default only if every gate
passes and the median paired wall and total-core ratios are both strictly below
`0.97`. Report RSS without adding a post-hoc RSS gate. A candidate that moves
work between outer and inner scheduling without reducing total process cost is
rejected as relabelling.

## Boundary

The same-binary direct-MITM decomposition and full automorphism-aware Pollard
rho references inherited through Stage 188 remain unchanged. Stage 189 cannot
alter any seven-gate row from a solver scheduling result. Complete campaign
cost remains null unless every inherited and new component is measured.

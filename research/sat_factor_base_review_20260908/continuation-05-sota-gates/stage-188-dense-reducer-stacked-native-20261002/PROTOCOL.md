# Stage 188 protocol: native stacked dense-reducer replay

## Question and hypothesis

Stage 178 tested a reusable exact leading-monomial array before dense critical-
pair selection became the repository default.  That candidate reduced symbolic
divisor probes from `4,190,633,182` to `103,532,494` and improved the median
paired wall ratio to `0.867978`, but its `0.976622` total-core ratio missed the
frozen `0.97` selection threshold and the implementation was reverted.

Stage 184 subsequently selected exact dense critical-pair updates.  The
pre-registered hypothesis is that, with pair UPDATE work reduced, symbolic
reducer lookup is now a larger fraction of the selected F4 path and the exact
Stage 178 dense reducer index may cross the unchanged wall and total-core gate
when stacked with the selected pair implementation.

This is a transfer test of a previously rejected engineering mechanism.  It is
not a new factor-base, relation-yield, DLP, rho, asymptotic, or SOTA claim.

## Frozen source and target

- Branch base: `9f941e2934c6730a2b1a3baca39f2f074ac50abe`.
- Stage 178 rejected patch SHA-256:
  `42ae7103f5810c51bb8e480a13c7a39e495bd0279d100f3ef779b54e0ba948b2`.
- Input manifest:
  `../stage-175-current-f4-single-target-20261001/input/manifest.json`.
- Input manifest SHA-256:
  `188dfdcb7e9398b0ab03193265f7b3035776d4a8a52de7dd4a51475dfc1a572c`.
- Opened public instance: `n=59, ell=9, m=3`, source-instance id
  `954e10f8bf0280094fed195280b203d7cd613150b17339716eb469b88ffa9ac7`.
- Fixed-X1 equation fingerprint:
  `02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb`.
- Algebraic factor base: `span_F2(1,z,...,z^8)`.  Target-subgroup enumeration
  and known discrete-log labels remain forbidden.

The candidate may port the exact Stage 178 reducer semantics to current
`pq_f4_f2.rs`; it may not change Semaev equations, the factor base, X1 order,
pair selection, matrix elimination, extraction, or verification.

## Native execution contract

New measurements use the Rust example `koblitz_f4_stage188`; no Python code or
Python process may construct, run, verify, compose, or hash the experiment.
The runner creates a fresh process group for every backend invocation, meters
wall/user/system time and peak RSS with `wait4`, enforces a 360-second watchdog,
hashes inputs and outputs with the repository SHA-256 implementation, validates
the backend terminal record, and writes immutable per-run receipts plus a phase
result and verification record.  A retry must use a new output directory.

Build both Rust examples from one clean candidate commit with the exact
evidence-scoped Stage 185 `Cargo.lock`.  Record the candidate commit and hashes
of the lock, backend binary, runner binary, manifest, protocol, and source files.
Compilation, tests, failed attempts, and every measured child process are
charged.  The twelve Rayon workers are included in total core-seconds; the
multi-worker experiment leaves `single_core_seconds` null.

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

The control adds `F4_F2_INDEXED_REDUCERS=0`.  The candidate adds
`F4_F2_INDEXED_REDUCERS=1`.  The repository default remains unchanged until a
confirmation passes.

## Correctness gates

1. A randomized exhaustive unit test must compare indexed and linear reducer
   answers over every monomial in small Boolean domains, including duplicate
   leading monomials, shortest-row ties, index ties, and both adaptive paths.
2. Current Boolean-F4, fixed-X1 specialization, and Phase B backend tests must
   pass before measurement.
3. Every measured run must authenticate the frozen source, visit all 512 X1
   masks, skip the same 270 non-rational masks, construct and complete all 242
   rational systems, find zero algebraic roots, and return exhaustive `UNSAT`.
4. The factor-base contract must report both target-subgroup enumeration and
   discrete-log labels as false.  Conflicts remain `null`, because F4 does not
   expose a SAT-conflict counter.
5. Control and candidate must agree exactly on equation fingerprint, pair and
   field-pair counts, matrix shapes, basis, extraction, degree, logical XORs,
   actually performed XORs, and every dense-pair-selector counter.  Only
   reducer-index counters, timing, RSS, and the charged peak-memory upper bound
   may differ.
6. The candidate must account for index allocation, and its reported divisor
   total must equal submask lookups plus linear tests.  The control must make
   zero submask lookups and allocate zero reducer-index bytes.

Any failed correctness gate rejects the candidate without interpreting timing.
A timeout or resource refusal is censored, not `UNSAT` and not evidence against
the mathematical method.

## Frozen schedule and decision

Run one screen pair in this order:

```text
linear, indexed
```

Continue only if both records are correct and the candidate/control wall and
total-core ratios are below `1.00`.  Otherwise retain the attempt and reject the
stacked transfer without confirmation.

If the screen continues, run three interleaved confirmation pairs in this
fixed order:

```text
linear, indexed, indexed, linear, linear, indexed
```

Select the indexed reducer as the repository default only if every correctness
gate passes and the median of the three paired candidate/control ratios is
strictly below `0.97` for both wall time and total core-seconds.  The threshold
is inherited unchanged from Stage 178.  Report peak RSS and absolute process
metrics, but do not add an unregistered RSS selection threshold.

If selected, make unset `F4_F2_INDEXED_REDUCERS` choose the indexed path and
retain `F4_F2_INDEXED_REDUCERS=0` as the exact linear control.  Then run a clean
default-mode replay and both explicit modes from the selected commit before
updating the scoreboard or default claim.

## Boundary and charging

The reference boundaries remain the same-binary direct MITM decomposition and
the full-cost automorphism-aware rho result already recorded by Stage 187.
Stage 188 is a solver-stage diagnostic and cannot change either boundary.  Its
table must label the current selected F4 control, the stacked candidate, direct
MITM, and rho in their existing units; unavailable full-pipeline quantities
remain null.

Charge native-runner and backend builds, tests, source/predicate validation,
every control and candidate process, total user plus system seconds, wall time,
peak RSS, all parallel work, verification, composition, and failed attempts.
Append the unique Stage 188 charge to the campaign lower bound without
rewriting Stages 175-187.

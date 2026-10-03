# Prepared one-target F5 versus inherited F4: frozen comparison protocol

Status: registration before any solver run. The existing source-bound n17
MatrixF5 complete solve and the prepared-state controls are historical inputs,
not measurements in this comparison. No old one-use registration is resumed.

## Question and boundary

The latest Boolean matrix-F5 complete-call optimization is in
`matrix_f5_f2.rs`; the IC worker calls `koblitz_groebner.rs` instead. This
experiment measures a production-path alternative that is already implemented
in the IC library: `InheritedF4`, which carries a Macaulay basis down the
Boolean splitting tree, against the worker's bounded `MatrixF5` engine. The
F5 criterion pruned zero rows in the documented degree-three ladder, whereas
inherited F4 greatly reduced counted work there. This is a solver-family
comparison, **not** an attribution of the separate F5 kernel optimization.

The primary boundary is a verified one-target signed-Frobenius rho solve on
the identical public point and one-worker resource envelope. For each IC
candidate report `rho_online_ns / ic_online_ns`. Also pair the two IC engines
on that point and report `f5_online_ns / inherited_f4_online_ns`. A 2.0×
or greater median paired online improvement, with all three repetitions
verified, is the predeclared development success gate. Any failure, timeout,
different point, incomplete phase ledger, or target-dependent work moved
outside the online interval fails that gate. This n17 toy is not a security
claim or ECC2K-130 extrapolation.

The measured objective is one target. Preparation, including the base,
certified logarithms and symbolic template, is reusable and excluded from
online timing. The five exclusive IC phases are target query, target PDP
including all failed attempts, relation check, descent, and scalar recovery
check. Cold/setup costs are supplementary. The rho interval starts with its
first target-dependent walk and ends after scalar replay. Retain attempts,
status mix, rank, base size, folded columns, wall time, operations and RSS
(null if unavailable). A missing phase or unverified scalar yields no ratio.

## Frozen inputs and arms

- Curve: Koblitz `a=1` over the existing n17 field, subgroup and generator
  from `KoblitzCurve::new(1,17)`; both IC arms import the exact
  `prepared-f5-runtime-v1/mathematics.json` state with 63 geometric points,
  62 usable projected points, 29 folded columns and verified rank 29.
- Public target: the worker's existing public hash-to-curve fixture policy
  with **seed `20261003017`**. The fixture process emits the point before any
  measured solver run. No scalar is constructed or supplied. All arms receive
  that same point as `public_targets` and an empty `target_seeds` list.
- IC algorithm seed: `20261003032`. Configuration: `summands=3`,
  `groebner_degree=3`, `node_budget=8192`, `batch_trials=1`,
  `max_trials=32`, `linear_algebra=dense`, standard subspace dimension 6,
  `exclusive_phases=true`, one Rayon worker. The only arm difference is
  `config.solver`: `f5` versus `inherited_f4`.
- Rho uses the same point and one worker, signed-Frobenius policy, seed
  `20261003032`, `rho_parallel_walks=1`, `max_trials=1000000`; the other
  config fields remain valid defaults. One rho observation accompanies each
  IC repetition.
- Three fresh processes per arm, in order F5, inherited F4, rho;
  inherited F4, F5, rho; F5, inherited F4, rho. This order is fixed before
  reading outcomes. Retain every attempted process and failure. Do not
  substitute another seed or target after seeing a result.

The patch to `ic_tournament_worker.rs` only admits the already existing
`inherited_f4` solver into the prepared mode guard. It changes no solver
implementation, preparation mathematics or existing `f5` dispatch. Compile
the single release worker once and run both IC arms from that exact binary.
No `KIC_`, `F4_`, `SOLVER_` or undeclared `IC_` environment override is
allowed. Record source SHA-256, binary SHA-256, OS/CPU, compiler, exact input,
stdout/stderr/exit status, and online phase ledgers. Use native Rust or
simple shell build/run commands; do not use Python for research execution or
analysis in this repository.

## Identities, verification and stop conditions

After the source and public point are frozen, create distinct canonical IC
candidate manifests for `PDP3f5` and `PDP3f4`, plus a canonical one-target
workload manifest and `R1`–`R3` rows. Preserve the prepared state as an
explicit input; actual usable `fb62` and folded columns 29 are separate
fields. The source digest and solver selection distinguish candidate IDs.
The worker must replay `[d]G=Q` for every successful solve. The separate rho
engine must recover that same scalar and replay it; the result record retains
both certificates. If an external Sage replay is added, launch it through
`/Volumes/SSD990/cryptanalysis/sage` and save its `--runtime-info` first.

Stop after the fixed three repetitions. If the fixture is malformed, the
prepared state fails verification, or the binary cannot be built within the
host's available disk/resource envelope, record a no-run failure and stop.
Do not infer a speed from a partial solver row. A successful development
measurement motivates a fresh unseen-target panel with new candidate IDs and
paired rho before promotion.

## Algorithm literature screen

There is no established Faugère “F6” successor in the checked algorithm
literature. [F4/5](https://arxiv.org/abs/1006.4933) combines F4-style matrix
reduction with F5 signatures; [F5C](https://arxiv.org/abs/0906.2967) reduces
intermediate bases; [M5GB](https://arxiv.org/abs/2208.00844) combines
signature criteria with cached tail-reduced reductors. These are concrete
alternatives to screen, not presumed improvements on Boolean Semaev systems.
The zero-pruning degree-three observation makes a new signature criterion
low priority on this frozen n17 cell. A follow-on M5GB-style experiment must
first show fewer complete-call operations and identical Boolean roots on
the exact production systems, then satisfy the one-target IC gate above.

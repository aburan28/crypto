# Stage 175 protocol: current repository F4 on one frozen target

## Question

Does the current repository Boolean F4 engine, with adaptive leading-block
`BlockTables`, improve the complete native-F4 decomposition cost on the one
already-opened Phase B target without changing the target, factor base,
fixed-X1 formulation, schedule, or correctness boundary?

This is a single-target implementation experiment. It is not a target-yield
experiment, a fresh holdout, an index-calculus relation campaign, a completed
unknown-scalar DLP, a rho crossover, or a SOTA claim.

## Frozen input

- Cell: `n=59, ell=9, m=3`, binary Koblitz curve
  `y^2 + xy = x^3 + x^2 + 1`.
- Blind target: `b-421e22a9c1c3b9d56396c8bbd0e46185bbebc32de0306f2bee6e9585703c2be4`.
- Source instance id:
  `954e10f8bf0280094fed195280b203d7cd613150b17339716eb469b88ffa9ac7`.
- Manifest SHA-256:
  `188dfdcb7e9398b0ab03193265f7b3035776d4a8a52de7dd4a51475dfc1a572c`.
- Expected fixed-X1 equation BLAKE3:
  `02341a5f51fd237b6a3fab8a82517047b974cd75664e6d9b02e4895e33252beb`.
- The factor base is the algebraic predicate
  `x = sum(c_i z^i), c_i in F_2, 0 <= i < ell`; the manifest states and the
  backend verifies that it neither enumerates the target subgroup nor uses
  discrete-log labels.

The four source artifacts are copied byte-for-byte under `input/`. The
manifest authenticates the three exported solver representations by BLAKE3.

## Candidate and reference

The candidate is the exact Git commit containing this protocol, built in
release mode from a clean checkout. It uses the current repository
`src/cryptanalysis/pq_f4_f2.rs` engine and the Phase B fixed-X1 adapter:

```text
RAYON_NUM_THREADS=12
PQ_F4_X1_BATCH=512
VECLIB_MAXIMUM_THREADS=1
OPENBLAS_NUM_THREADS=1
OMP_NUM_THREADS=1
MKL_NUM_THREADS=1
BLIS_NUM_THREADS=1
NUMEXPR_NUM_THREADS=1
koblitz_pdp_backend native-f4 input/manifest.json 300
```

No Stage 174 `PQ_F4_*` solver switches are valid in this candidate. The live
engine controls are `F4_F2_BITMAP_SEEN=0` and `F4_F2_BATCH_INSERTS=0`; neither
is set in the candidate.

The immutable Stage 174 selected reference used the same target and fixed-X1
schedule. Its three-run median was `26.358217916989815` wall seconds,
`147.841771` total core-seconds, and `2756362240` bytes peak RSS. Its selected
representative was `26.358217916989815` wall seconds, `143.569806` total
core-seconds, and `2676736000` bytes peak RSS. Those figures remain historical
measurements at source commit `e51efb219edd4511a08d545de21184928e5e771c`.

The candidate binary also runs direct MITM three times on the identical input.
This provides a same-binary, same-host decomposition reference. It does not
stand in for full index-calculus cost or automorphism-optimized rho.

## Repeats and accounting

- Three fresh candidate processes and three fresh direct-MITM processes.
- `scripts/process_meter.py` records wall time, user plus system core-seconds,
  and peak RSS for every process.
- The F4 report charges source authentication, algebraic factor-base
  membership checks, construction of every rational fixed-X1 system, all F4
  calls including failed branches, solution extraction, and exact curve lift.
- `cost.ops` is the stable row-by-row-equivalent elimination XOR count.
  `cost.extra.word_xors_performed` separately records actual table-assisted
  XORs including table construction. Matrix and table memory are both exposed.
- The twelve requested Rayon workers are charged through total core-seconds;
  there is no per-worker wall-time discount.
- Build cost is recorded separately and is not called query time.

## Correctness gates

Every candidate repeat must satisfy all of these or the candidate is rejected:

1. exit code zero, no watchdog timeout, and backend status `unsat`;
2. authenticated source and regenerated source both exact;
3. factor-base contract remains `false` for target-subgroup enumeration and
   discrete-log labels;
4. all 512 X1 masks visited, all 242 rational systems constructed and
   completed, and `exhaustive=true`;
5. zero algebraic roots, no witness, and `conflicts=null`;
6. equation fingerprint equals the frozen value above;
7. all repeated operation counts and structural counts agree exactly.

Direct MITM must return the same exhaustive `unsat` result for the same source
instance id on every repeat.

## Decision rule

This stage is an engineering improvement only if the candidate's three-run
median wall time and total core-seconds are both below the Stage 174 medians.
Peak RSS is reported as a separate resource ratio and may not be hidden. A
failure of either timing condition, any correctness gate, or any structural
disagreement rejects the candidate.

Even a passing result changes none of the seven SOTA gates: it is one opened
public target and one decomposition backend. Full factor-base discovery,
natural independent-relation yield, relation collection, linear algebra,
unknown-scalar recovery, all same-cell external backends, matched
automorphism-rho cost, independent reproduction, and novelty review remain
separate requirements.

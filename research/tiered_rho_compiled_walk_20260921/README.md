# Compiled reservoir-guided tiered rho study

This package specifies and checks a theoretical current-plus-archive collision
walk with a target-independent compiled schedule, bounded RAM locality, Bloom
coverage guidance, and an exact NVMe archive model.

The study is public-synthetic. It has no external target interface and makes no
private-key, practical ECDLP, 256-bit feasibility, asymptotic, or measured
wall-clock speedup claim.

## Main result

The online recurrence is

```text
X[t+1] = X[t] + X[t+1-d[t]],  1 <= d[t] <= 64.
```

A bottom-k/Bloom coefficient compiler chooses the six-bit offsets `d[t]` once
per public group order. Online execution keeps the current record and 64 hot
records in RAM, performs one group addition, appends one exact record, and
probes an exact collision directory per step.

The pinned 321-byte schedule is exhaustively checked across all 65,537 target
logs. Its exact capped mean is `19027631/65537` additions, versus
`19114086/65537` for the matched unguided schedule. The finite online reduction
is `86455/65537`, or 0.4523%, after compilation.

The schedule SHA-256 is:

```text
3a07aeaf38483b9facabe21ef8415be95617534af03a6328d05d846a2efbccd3
```

## Start here

- `COMPILED_RESERVOIR_GUIDED_WALK.md` — mechanism, controls, scaling, capacity,
  amortization, and novelty boundary.
- `TIERED_RHO_COMPLETION_AUDIT.md` — requirement-by-requirement claim audit.
- `certified_compiled_walk.json` — exact all-target stopping histogram and
  source/schedule hashes.
- `certified_compiled_walk_schedule.bin` — the pinned six-bit schedule.

`COVERAGE_GUIDED_CHAIN.md` records the full-archive upper-bound mechanism that
motivated the bounded compiler. `MULTI_TARGET_LINEAR_COLLISIONS.md` records the
known multi-user `sqrt(KN)` positive control and its rank-deficiency fallback.

## Reproduce

From this directory, using Python 3:

```bash
python3 -m py_compile check_*.py
python3 check_certified_compiled_walk.py
python3 check_walk_reservoir_guided.py
python3 check_walk_reservoir_bloom.py
python3 check_compiled_walk_scaling.py
python3 check_compiled_walk_capacity.py
python3 check_coverage_guided_chain.py
python3 check_multi_target_linear_collisions.py
```

Each checker uses fixed public parameters and seeds. Generated JSON records the
SHA-256 of the exact checker source. The certificate checker also regenerates
the binary schedule and checks its hash.

## Claim layers

Established in this package:

- exact finite coefficient and collision correctness;
- exhaustive all-target stopping improvement for the pinned schedule;
- paired guided/unguided finite controls;
- exact/Bloom transcript equality at the tested load;
- three-order finite scaling screens; and
- exact decimal 1 TB RAM / 100 TB peak-disk envelopes.

Not established:

- nonvanishing asymptotic improvement;
- positive measured wall-clock amortization;
- classical fixed-state rho coalescence;
- random-probe service for a 100 TB external compiler sketch;
- cryptographic-scale record feasibility; or
- exhaustive independent publication or patent novelty.

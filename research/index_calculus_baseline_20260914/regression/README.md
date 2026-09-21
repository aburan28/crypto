# Frozen WDSat regression target

This directory turns the six-case pilot into a repeatable **60-input, two-configuration,
three-repetition** benchmark: 360 solver processes per batch. It covers the full
pinned S4 corpus at `(n,l) = (15,5), (17,6), (19,6)` on
`y² + xy = x³ + x² + 1`. The two configurations use the same WDSat binary with
symmetry breaking on and off, with Gaussian elimination disabled in both.

The accepted target is [results/baseline_v2/summary.json](results/baseline_v2/summary.json).
[RESULTS.md](RESULTS.md) carries its measurements and the identical scoreboard table.
[contract.json](contract.json) fixes inputs, repetitions, run-order seed, 15-second
per-process timeout and 8 GiB address-space cap. Three repetitions of the same
input measure repeatability; they are not three independent problem instances.

## Run every iteration against this target

From this directory, with Python 3.10+, Linux, GCC and libm:

```bash
# Build the unmodified reference once. Choose a fresh build directory.
python3 ../pilot/build_pilot.py --cache-dir /tmp/ec-baseline-cache --build-dir /tmp/ec-reference-build

# Baseline binary replay. Omit --offline on the first run to fetch missing inputs.
python3 run.py --cache-dir /tmp/ec-baseline-cache \
  --solver /tmp/ec-reference-build/wdsat_solver \
  --source-revision 61c6ff3f49445af2729274819edb27b85a89efc1 \
  --output /tmp/ec-iteration-001

python3 compare.py results/baseline_v2 /tmp/ec-iteration-001 \
  --output /tmp/ec-iteration-001-comparison.json
```

For a candidate, supply its actual executable and source revision to the same
`run.py` command, with a new output directory. The executable must support the
WDSat command/output protocol, identical ANF semantics and the same conflict
counter definition. This script is not an adapter for Magma or a different SAT
counter. Broader solver/family experiments use the parent accounting contract.

The comparison refuses missing/duplicate trials, changed input fingerprints,
changed outcomes, incompatible contracts/checkers/counter units, incomplete runs,
and corrupted evidence. It reports **candidate / frozen baseline** per input and
configuration, including every per-case conflict regression above 10%.
The optional `--require-counter-target` exits unsuccessfully unless both
configurations cut total per-case median conflicts by at least 20%, with no case
more than 10% worse. That is a diagnostic solver-counter target, not an end-to-end
operation or complexity improvement. Without that option, successful exit means
the comparison is complete and compatible; read its ratios to assess regressions.

These 60 inputs are a fixed regression set. Independent holdout curves and uniform
targets are required before claiming a general research improvement.

Keep `results/baseline_v2` immutable. Save each iteration under a new directory,
commit its evidence and comparison, and append its row to the scoreboard. If the
corpus or benchmark protocol changes, create a new version and retain the old
baseline. The scripts refuse existing output paths. The benchmark is manually
invoked; this PR does not install a scheduled job.

## What is measured

Each raw JSONL record retains exact stdout/stderr, command, exit status, input hash,
conflicts, process wall time, child CPU time, and verification evidence. Timers
separate the solver process from SAT checking and the exhaustive UNSAT oracle.
The oracle runs once per distinct UNSAT input and is reused only within that batch.
Three-repetition nearest-rank p90 is the maximum observed repetition, not a tail
latency estimate from a large sample. No failed runs are dropped.

The loop wall time includes orchestration and independent checking but excludes
upstream compilation/download and initial input preparation. It is not the cost
of a relation collector or complete ECDLP. Peak RSS is unmeasured; the address-space
cap is enforced. The host is shared and CPU affinity is recorded, not pinned.
For a timing claim, rerun both reference and candidate binaries on the same host.

The primary repository metric remains **total calibrated common operations /
sqrt(subgroup order)**, with ratios to measured rho and an applicable derived
floor. Those quantities are null for this solver-only corpus. Conflicts are a
named solver counter, not curve additions or total operations. No family-wide
success rate, scaling exponent, or rho crossover follows from these inputs.

## Algebraic SAT is not always a point relation

The original strict point-only batch completed all 360 processes but failed its
admission gate on six repetitions of `n19l6-19-U`. Its unmodified raw data and
contract are retained under [results/strict_point_v1](results/strict_point_v1), with
the historical runner in [history/strict_point_v1_run.py](history/strict_point_v1_run.py).
That archived runner expects its historical contract and directory layout;
use the maintained commands above for current replays.

The existing repository audit already identifies that filename as algebraically
SAT. This run further checks the point-lifting distinction: its S4 root
`(20,34,41)` satisfies the ANF equations, but neither those abscissas nor target
`x=171112` lifts to this curve over `F_(2^19)`. The separate trace certificate is
recorded in [negative_control.json](negative_control.json).

Before the second batch, v2 explicitly designated this input as an
**algebraic-only rejection control**. It must satisfy ANF/S4, fail target lifting,
and contribute zero point relations. The full v2 batch was then rerun from
scratch. The useful target is 30 point-relation inputs, 29 certified polynomial
UNSAT inputs, and one algebraic-only input, each exercised in both configurations.
Changing that last case into a claimed point relation fails comparison.

## Verify saved evidence

```bash
python3 -m unittest discover -s . -p 'test_*.py' -v
python3 compare.py results/baseline_v2 results/baseline_v2 \
  --output /tmp/ec-self-comparison.json
```

Self-comparison must give ratio 1 for each configuration and must not meet the
20% improvement target. The parent's artifact manifest hashes all committed
results, contracts, scripts and notes. Upstream GPL source and generated input
files stay in the external, hash-verified cache; they are not vendored here.

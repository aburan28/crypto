# EC index-calculus baseline, 2026-09-14

An accounting baseline for Koblitz, other binary, and generic prime-field research:
[report](EC_Index_Calculus_Performance_Baseline.md), [contract](ec_index_calculus_contract.json),
[budget calculator](ec_index_calculus_budget.py), and a certified six-instance WDSat pilot.
The solver is unchanged. No full-DLP speedup or complexity improvement is claimed.

For the fixed target used by subsequent iterations, use the
[full-corpus regression suite](regression/README.md). It expands the pilot to
60 inputs, two configurations and three repetitions, with saved raw results,
an explicit algebraic-only rejection control, and a comparison command.

The repository already has the [corpus audit](../../RESEARCH_TRIMOSKA_BENCHMARKS.md),
[F4 notes](../../RESEARCH_GROEBNER_F4.md), [Koblitz experiments](../../RESEARCH_KOBLITZ_SCALING_TARGET.md),
and [Nagao comparison](../nagao_relations/solver_09/README.md). Use those runners for
subsequent matched experiments; this directory freezes an independent control and
an accounting contract. The [scoreboard](../../docs/index-calculus-scoreboard.html#ec-baseline-20260914)
records its admission status. Common-operation conversions, full-DLP S and rho
ratios are unmeasured, represented by null in [summary.json](summary.json).

## Reproduce

Requires Python 3.10+, Linux, GCC and libm. From this directory:

```bash
python3 ec_index_calculus_budget.py
python3 -m unittest discover -s . -p 'test_*.py' -v
python3 pilot/build_pilot.py --cache-dir /tmp/ec-baseline-cache --build-dir /tmp/ec-baseline-build
python3 pilot/run_pilot.py --cache-dir /tmp/ec-baseline-cache --solver /tmp/ec-baseline-build/wdsat_solver --output /tmp/ec-baseline-replay
```

Choose new build/output directories for each replay: existing directories are refused.
The build fetches the exact files in [source_manifest.json](pilot/source_manifest.json),
including upstream licence notices, and verifies SHA-256 before compilation. An
already populated, hash-verified cache can be used with `--offline`. Source and
input files are downloaded separately, following the repository's existing
non-vendoring convention. They are not linked into the Rust crate.

WDSat is pinned to `61c6ff3f49445af2729274819edb27b85a89efc1` and the corpus to
`a053971cb155fbab30eaed3e2d11df850b4f3fbb`. Upstream static configuration is unchanged:
GCC `-O3 -Wall`, no OpenMP, symmetry on, GE off, allocation configured for l=6.
Compiler warnings are retained. The l=5 cases are untuned allocation controls.
Do not generalize the small pilot runner into a production or large-instance solver.

The runner removes blank lines from ANF inputs into the new output directory,
checks every SAT assignment against all equations, then S4, point lifting and
independent group addition. Every UNSAT input is exhaustively checked over all
unordered support triples with repetition: 97,504 S4 evaluations across three
cases. First-solution process time excludes the separately recorded verification
time and is not an all-phase decomposition price. The pilot does not verify a
subgroup logarithm. A timeout stops the replay with incomplete evidence retained;
exceptions or nonzero exits invalidate the run. Do not run Python with `-O`.

## Frozen evidence and model

`pilot/results/` preserves the original 18 solver outputs, compiler diagnostics,
measurements and historical runner/build scripts. Those original scripts wrote
into their own directory; use the maintained commands above for safe replays.
The new runner strengthens certificate checks and records build/binary provenance.
[artifact_manifest.json](artifact_manifest.json) protects all committed files in
this directory except itself. Neither build products nor downloaded inputs belong
in a commit.

The default budget is an **illustrative model**, not an actual curve measurement:
N=2^60, M=2^20, B=2^19, m=3, ideal negation quotient and an assumed rho rate.
Use actual subgroup and support sizes, measured yield, calibrated rho rate and
all other-stage costs for a machine-specific planning budget. Zero success yields
no finite forecast. The attempt-count floor includes repeated summands and is
limited to uniform independent targets in first-relation mode. Matrix work is
estimated in the subgroup scalar field, separately from Boolean elimination.
The calculator never invents common-operation conversions or measured rho ratios.

The contract defines phase-exclusive timers, typed counters, status/certificate
requirements and proposed parameter sweeps. Its 20% paired operation-cost gate is
an engineering threshold. Full-DLP admission also requires every phase priced and
the recovered scalar verified; an exponent claim needs at least four sizes.

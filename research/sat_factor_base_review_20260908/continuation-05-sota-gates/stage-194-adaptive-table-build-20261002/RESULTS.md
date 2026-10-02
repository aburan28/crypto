# Stage 194: adaptive BlockTables construction scheduling

## Decision

`REJECTED_SCREEN`. The candidate does not clear the preregistered joint
`0.98` wall-and-total-core gate, so the three-pair confirmation is prohibited
and the selected runtime remains unchanged.

| arm | wall seconds | total core-seconds | peak RSS |
|:---|---:|---:|---:|
| current unconditional Rayon table adds | 16.869204 | 180.921311 | 3,939,205,120 B |
| adaptive table-build threshold | 16.840656 | 183.553080 | 3,835,707,392 B |
| **candidate / current** | **0.998308** | **1.014546** | **0.973726** |

The 0.17 percent wall change is below the frozen useful-effect threshold, and
total CPU regresses 1.45 percent. The 2.63 percent RSS reduction is reported
but cannot rescue a failed wall-and-CPU gate.

## Mechanism and correctness

The candidate moves 1,647 of 5,500 `BlockTables::add` calls from a Rayon
parallel iterator to serial iteration. It retains 3,853 parallel calls. Both
arms exactly agree on:

- 728,503 table runs and 1,700,553,592 scheduled table words;
- 147,794,583,858 actually performed table/row XORs and 319,313,687,585
  row-equivalent logical XORs;
- all 512 fixed-X1 masks, 270 non-rational skips, and 242 completed rational
  systems;
- the source instance, equation fingerprint, equation/term counts, pair and
  field-pair counters, reducer rows, matrices, degrees, basis, and zero roots;
  and
- exhaustive `UNSAT`, with target-subgroup enumeration and known
  discrete-log-label use both false.

The aggregate in-call diagnostics move from 42.351720 to 43.972401 seconds for
matrix build work and from 127.463648 to 120.368083 seconds for elimination.
Those are sums across concurrently executing F4 calls, not process wall or CPU
times. The charged process totals decide the gate.

The candidate implementation is preserved in
`candidate-adaptive-table-build.patch` and reverted from runtime source. The
first launch failure is also preserved: the meter stopped before spawning the
backend because the screen parent directory was absent. The metered layout
correction then ran the unchanged frozen order exactly once.

## Verification and accounting

The Rust verifier authenticates commands, explicit arm modes, meter receipts,
source and equation identities, algebraic factor-base claims, all invariant
work counters, the scheduling split, ratios, decision, candidate patch, and
runtime reversion. Every one of the 16 charged metrics files has exactly one
authenticated receipt. Final replay passes `19/19`; result SHA-256 is
`918687ed0358356d1309c866ede15b003eec9d3b94a09ce9d5b5b65ece4abaad`.

Stage 194 contributes a measured lower bound of 16 components,
`216.945372` wall-seconds, `1,203.588335` total core-seconds, and
`5,264,408,576` bytes peak RSS. The cumulative measured campaign lower bound
is 630 components, `24,380.108566` wall-seconds, `63,014.492177`
core-seconds, and `6,310,576,128` bytes maximum RSS.

Complete campaign cost remains `null`. This is a rejected one-target
scheduling experiment. It changes no natural relation-yield, unknown-scalar
recovery, full automorphism-aware rho, independent reproduction, novelty, or
Koblitz index-calculus SOTA gate.

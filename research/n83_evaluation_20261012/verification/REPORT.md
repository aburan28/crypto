# Native independent replay of the retained degree-83 online observation

Dated 2026-10-12. Source snapshot: `aburan28/crypto` commit
`a8194d616ed8c04489ec46177d369f1b1670360a`.

The standalone C++ validator independently reproduced the archived curve,
subgroup, factor-base action, three four-point witnesses and six recovered
scalar answers. It used ordinary polynomial-basis GF(2^83) arithmetic and
imported no research solver. This is a new archive correctness check; the IC
and rho algorithms were not rerun or retimed by this validator.

![Scope of the independent arithmetic replay; the fresh algorithm comparison remains separate.](verification.pdf)

## Exact mathematical instance

The model is `icv1-f2m83-t6151469093347-cdcc5432`:
`y²+xy=x³+x²+1`, modulo `x^83+x^7+x^4+x^2+1`.
Independent Rabin irreducibility passed. The Koblitz trace recurrence gives
`#E = 9671406556910881928556062`, equal to subgroup prime
`r = 8569786107849059` times cofactor `1128547018`.
The subgroup prime was checked by the deterministic seven-base
Miller–Rabin test valid over 64-bit integers.

The generator and public target have these polynomial-basis coordinates:

```
G_x = 6614244651398184860292015
G_y = 8868260410029249391539529
Q_x = 2655898345537166180532960
Q_y = 3881526307991673265292615
```
Both are nonidentity points on the curve, annihilated by r.
The validation-only sidecar scalar `3141592653589793` replays to this target.
No new unknown-scalar solution was computed.

The scalar `lambda = 2998004823950457` satisfies
`lambda²-lambda+2=0 (mod r)`, has exact order 83 and obeys
`tau(G)=[lambda]G`. Every one of the 600 factor-base representatives lies
in the subgroup. All 99,600 retained points were checked to be distinct,
on curve, exactly the 166 signed Frobenius images of their representatives,
with correct `(column,±lambda^k)` labels. Each orbit closes after 83 shifts.
All three four-point target witnesses re-add to Q; their x codes and two
pinned intermediate sum abscissae match the retained receipts.

## Recomputed archived interval arithmetic

These are historical recorded milliseconds on the original contended ARM64
host, not measurements on this Linux checker. Three repeats share one public
point. Ratios use internal rho recorded total / IC target interval.

| Run | IC target ms | Rho internal total ms | Rho/IC | Unassigned IC ms |
|---|---:|---:|---:|---:|
| R1 | 7939.137708 | 64952.483042 | 8.181301979 | 0.000458 |
| R2 | 7746.014250 | 42987.579541 | 5.549638582 | 0.000084 |
| R3 | 5623.275125 | 71253.689916 | 12.671208207 | 0.000541 |

Each five-phase sum matches the recorded `target_phase_sum_ms` but differs
slightly from `target_ms`. The gap is less than one microsecond; an exact
exclusive-phase contract must assign that elapsed bookkeeping cost rather
than state exact equality. Rho totals equal setup+walk+validation; setup is
7.556, 7.816083 and 7.516250 ms. Its target dependence remains unresolved.
The source README rounded rho figures differ from the retained internal
totals; the CSV and this check use the latter.

Each IC repeat records 8,845,441 probes. Rho walk steps are 6,608,900,
5,958,775 and 12,179,440. Those counters count different work and do not
form an algorithm speedup metric. Reusable IC construction/rank setup is
excluded from this online comparison.

## What remains unresolved

The 600 rank equations and recovered factor-log database are not retained
in the checked summary rows; their independent algebraic replay remains
unverified. The final recovered scalars and retained target witnesses pass
independent arithmetic replay. This does not prove the implementation's
rank construction history or absence of scalar input beyond the retained
metadata/launch description.

The archived observations have no qualifying isolation receipt, enforced
matched memory cap, complete canonical candidate/workload manifests or
independent execution replay of both algorithms. The IC process used 12
threads for reusable rank setup; the current source defaults target
extraction to one thread unless `KIC_TARGET_THREADS` is set. No target-thread
envelope is frozen in the historical producer receipts.

The requested fresh same-point benchmark remains a separate obligation,
being handled by the parent task. This validator does not discharge it.
No repository, ledger or canonical graph was edited. The archival
observation retains its exploratory status; existing scoreboard conclusions
do not change from this arithmetic recheck.

## Reproduction and evidence

```
g++ -std=c++17 -O2 -Wall -Wextra -Werror check.cpp -o check
CHECK_SOURCE_DIR=/workspace/repos/crypto/experiments
CHECK_SOURCE_DIR=$CHECK_SOURCE_DIR/koblitz-single-target-n83-20261006
./check "$CHECK_SOURCE_DIR" receipt.json
```

Compiler: Debian g++ 14.2.0; Linux x86-64. The run completed successfully.
`source_hashes.sha256` binds the fixture, full base and all three producer
and summary records. `receipt.json` contains exact arithmetic outputs.
`check.cpp` is the editable independent arithmetic implementation.

Receipt SHA-256:
`7d4e290524ea9b83ff75d48d20b17449fb60a6991bae787fdfee53ef28d461d5`.

The [pinned experiment evidence](https://github.com/aburan28/crypto/tree/a8194d616ed8c04489ec46177d369f1b1670360a/experiments/koblitz-single-target-n83-20261006)
contains `frozen/fixture.json`, `base_n83_K600.jsonl`, and the R1, R2 and R3
receipt sets under `runs/`. Exact source-file hashes are retained beside this report.

The explanatory [verification diagram](verification.svg) shows the scope of
this check. Its editable SVG is included with this report.

# Fresh ordinary-query natural-yield sample for static wide-S4 SAT

The [earlier static controls](../static-cms-s4-controls/RESULTS.md) established
that the source-receipted wide-S4 CryptoMiniSat stage can return valid point
decompositions. Those six points were selected after their exact outcomes
were disclosed. This [panel](panel.json) instead freezes the first **32**
queries of a new trial-keyed `StdRng08` stream at algorithm seed `2026092935`.
The stream's independent Python replay is already tested against the pinned
Rust query law. It yields 32 distinct n17a1 public points, with zero overlap
with the 104 previously retained F5 ordinary queries. The schedule was fixed
before any SAT instance or new exact-group classification. It contains no
feasibility or prior-solver outcomes. The solver subprocess receives one
public point and public nonce at a time, never an exact label.

Run the same pinned-source Rust wide-S4 `--export-only` frontend and the
source-receipted static CryptoMiniSat 5.14.7 executable as the six-control
panel. Verify the archived Phase-B build bundle, parent worker/base source,
curve/subgroup, all 32 scalar-to-point calculations and query-law values,
then preflight the copied executable's linkage and `--version` at its final
run path. Freeze one thread, one returned model, one million conflicts,
60-second exporter and 120-second SAT watchdogs. Execute exactly 32 in trial
order with no retries, replacements or outcome-dependent stopping. Keep all
export and solver stdout/stderr, instance bytes, metrics, statuses, invalid
models, nonlifting models, timeouts, errors and zero-yield rows.

For each SAT answer, check the full XOR-DIMACS source assignment, lift the
three abscissae to the verified geometric factor base and group-replay a
witness against its public point. A source SAT assignment without a group
witness is nonlifting, not a relation. A timeout or unknown is not UNSAT.
After the one-shot schedule finishes, use a **separate** exhaustive group
three-sum oracle to classify each query's mathematical feasibility, including
every solver failure. Cross-check measured point witnesses against those
labels. Report verified witness count per 32 ordinary queries, mathematical
feasibility, missed feasible points, attempt-status mix and a Wilson 95%
interval as a small-sample uncertainty description. Retain process wall,
CPU and peak RSS per query, but do not compare these costs directly with
prior F5 rows that lack per-query wall intervals.

This is stage evidence on a fresh ordinary-query law, **not** a target DLP
workload. It does not collect a full-rank matrix, solve factor-base logs,
descend an unseen target or run matched rho. Thus it cannot produce an `IC1`
candidate result, one-target online time or speedup. A subsequent protocol
must implement and measure the complete SAT pipeline, F5/incumbent and rho
on new paired one-target workloads. No sealed confirmation set is reopened.

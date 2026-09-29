# Disclosed-point F4/F5/SAT complete-recovery pilot

This registration extends the [source-bound encoder and base audit](../generic-backend-qualification-v2/STATIC-FEASIBILITY.md)
and the [one-query dispatch pilot](../generic-backend-disclosed-pilot/RESULTS.md).
It executes the [frozen panel](panel.json) on the five previously disclosed
public points in `standard-subspace-d6-inventory-control.json`. The latter's
exact file digest, point coordinates, field/curve/subgroup records and base
inventories are fixed before measurement. No new target or sealed
confirmation/replay point is used. The separate F4/F5-only [recovery plan](../generic-f4-subspace-pilot/PROTOCOL.md)
is not dispatched by this runner; this panel also fixes both in-process SAT
encodings and its resource limits.
The registered panel SHA-256 is
`bf586fa75f9b5a00c6de95793c73502a2fb064898f1d4eb09c84685082b7b03d`.

**Hypothesis.** On n17a1, at least one of the F4/F5 engines and at least one
of native-XOR/CNF SAT can reach full relation rank and recover the supplied
target through a verified IC descent under the frozen limits. On all five
cells, report the natural ordinary-query outcome mix, usable base, folded
columns, rank and bounded completion status for every arm. A failure on any
cell stays in the table; it cannot be replaced with a new seed or limit.

All twenty cell-by-solver jobs run once in panel order. The pinned worker
source is `765c3c5f19032bd852163805f257c56babef2040`, built by
`generic_build.py` with a retained source/dependency manifest and executable
hash. The polynomial subspace has dimension six; this is an explicit
factor-base-policy choice, not an isolated solver change from the orbit-base
incumbent. Each job uses three summands, dense final relation LA, one Rayon
thread, the same point and algorithm seed across solvers, and the per-cell
collection/target trial caps in the panel. F4/F5 use degree 3 and 4096 search
nodes per PDP attempt. SAT uses the worker's native-XOR or Tseitin-CNF CDCL
backend and 100000 conflicts per attempt. The n17a1 128-trial budget tests
complete recovery; n19a0 has 24 trials, n23 cells four, and n31a0 one. Each
attempted job has its own 300- or 180-second cap and an 8 GiB child-RSS
threshold sampled every 100 ms. Short memory peaks may be missed, so this is
an approximate diagnostic resource guard, not a strict memory comparison.

Keep raw stdin job, stdout/stderr, process exit, wall/RSS diagnostic and the
controlled build receipt. `status: incomplete` at CLI exit code 2 is an
auditable bounded report, not a process failure. Independently reconstruct
the actual factor base and its cofactor image, replay every collection and
descent query and group relation, check the observed engine, matrix rows/rank
and factor logs, verify target descent/scalar replay when complete, and check
exclusive phase closure with `generic_admission.admit`. A timeout, OOM,
missing report, unsupported encoding, failed audit or unverified scalar is
not a completed solve. Charge failed attempts to the run and preserve their
outcomes. Missing phases and complete costs remain null.

The pass condition for this **disclosed-point admission pilot** is at least
one verified complete F4/F5 run and one verified complete SAT run on the same
n17a1 point. This does not qualify either for the tournament or establish a
speedup. No incumbent or rho worker runs here. Local macOS wall and RSS values
are feasibility diagnostics only; the competitive unit requires a separately
frozen Linux/Valgrind profile with a strong source-bound incumbent and matched
one-target rho, fresh paired points, no prior exposures, and confidence
intervals. Failed cells and solvers are retained rather than filtered out.

The v2 seed `2026092902` remains a single live dispatch
([Actions 36580669479](https://github.com/aburan28/crypto/actions/runs/36580669479))
and is never retried. Its artifact, when terminal, must be independently
audited before any fresh competitive registration. This pilot neither
relabels nor reopens that campaign.

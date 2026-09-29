# Fresh paired target allocation

This [target panel](target-panel.json) fixes one new n17a1 public point before
the static-SAT correctness pilot has completed. It is a target allocation,
not yet an executable comparison registration. All four listed arms must run
on `[52411,72106]` regardless of the first arm's outcome; no early winner or
survivor selection changes the lineup. The SAT and F5 arms use the same
standard-subspace dimension-six geometric base. The qualified pair-table
incumbent may use its own factor-base policy, so its comparison explicitly
changes both algorithm and base policy.

Construct the point with SHA-256 over the ASCII domain
`ic-paired-target-v1`, a zero byte, seed `2026092948` as little-endian u64,
and a little-endian u64 counter. Take the low 17 bits of the first eight
digest bytes (little-endian) for the abscissa, choose lift number
`digest[8] mod lift_count`, multiply by the cofactor, and skip points sent
to the identity. Counter one is the first valid point. No scalar is supplied
to any solver. Its subgroup membership and lack of prior exposure must be
checked before execution.

Before any arm starts, seal exact candidate and workload manifests, source
receipts, per-arm collection and target limits, memory/wall policy, process
order and a physical-host calibration record. Run each arm once under that
policy, retaining all failures. Charge every target-dependent attempt to the
single target and stop the online clock at independent scalar replay.
Report setup and online costs separately. An unverified result, a missing
phase, or a failed same-point rho reference leaves speedup unknown. This
one-point comparison can diagnose engineering differences; it cannot prove
a global optimum or a population-level speedup.

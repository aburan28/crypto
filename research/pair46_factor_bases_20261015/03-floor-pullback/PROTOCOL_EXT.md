# PROTOCOL_EXT — 03-floor-pullback extension E1 (frozen before E1 data)

Written 2026-10-15 after the main protocol's population turned out to hold
only 3 instances (v in {359, 409, 443}; main results in results.json) and the
main run finished in ~13 s of the 40-min budget. Not edited after.

Change from PROTOCOL.md: population only. Toy family, arms, metric, base size,
seeds, controls, statistics and decision rule are identical (same script,
same functions). Population E1: every prime v in [920, 2000] with
p = 5v^2 - v + 1 and N = 5v^2 - 2v + 4 both prime: v in {1213, 1283, 1409,
1549, 1613} (the list was computed before any E1 metric data; it carries no
outcome information). Instances run in increasing v; stop rule as before.
Output: results_ext.json, run_ext.log.

Deviation D1 (applies to both runs, recorded here because PROTOCOL.md is
frozen): PROTOCOL.md says order-v twist points come from [5v]*random points.
The crater twist has exponent 5v (structure Z/5v x Z/v, checked), so [5v]
kills every point; the script uses [5]*random point instead. This changes no
arm, metric or decision rule; it was found when the first run hung.
Implementation note: log tables use int32 arrays (memory), same values.

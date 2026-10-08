# Pilot R1 by a concurrent session (2026-10-03 ~00:58-01:01 UTC)

A concurrent agent session froze the n=71 fixture (target + scalar sidecar)
and attempted a first paired observation in this directory before the
decimal-string wire-format fix landed in `koblitz_orbit_dlp_fast.rs`:

- `freeze_receipt.json`, `public_target.json`, `target_points.jsonl`,
  `target_scalar_validation_only.txt` — the frozen fixture (scalar
  314159265358979, since solved by an independent smoke run; kept here
  for provenance, NOT used as evidence).
- `rho.jsonl` — a VALID rho observation with the wide backend
  (`packed_u128_polynomial_basis`, n=71, A=142): same public point,
  scalar recovered and verified, walk 28,496.7 ms over 5,315,897 steps.
  It ran against a binary built from the in-tree wide rho code.
- `ic.jsonl`, `ic.stdout.txt`, `ic.stderr.txt`, `ic_external_time.txt` —
  EMPTY: the IC arm ran against a binary predating the wire-format fix
  and panicked before emitting rows. No IC evidence from this attempt.

The n=71 record below (R1-R3) was produced independently with fresh
scalars after the fix, and does not depend on this pilot.

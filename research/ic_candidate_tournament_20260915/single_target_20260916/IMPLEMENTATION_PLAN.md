# Round 0008: fixed two-orbit implementation costs

Predeclared before new fixtures or measurements. Parent: audited round 0007,
conditional on promotion of `orbits2`. If that condition fails, do not run.
Only development phase costs from round 0007 motivate these mechanisms.

One cold target, unchanged signed base and m=3, same public targets and rho API.
Coverage is at most C(B+2,3); the full-rank collector requires at least K
relation-producing trials. No change to either boundary is proposed. These are
engineering changes, not a new ECDLP exponent or a generic-group advance.

Candidates against the promoted parent:

* `orbit_projection`: compute cofactor multiplication once per signed orbit;
  reconstruct conjugates and signs in original point order. Setup and descent
  should fall; O(B) temporary memory remains. Exact pointwise equivalence tests
  are required before admission.
* `serial_collection`: use serial iteration for the ordinary collector's batch.
  Batch size remains one; preserve probe scalars, relations and ordering. Removes
  dispatch overhead on the declared single CPU, without a multicore claim.
* `combined`: both mechanisms, measured as a whole executable.
* `combined_dense`: combined implementation with dense scalar elimination for
  the two-column system; certificates and rank requirements remain unchanged.

Each mechanism is falsified as a gain if complete cost does not fall. Promotion
retains the sealed >=20% reduction in BOTH instructions and native time, upper
paired 95% intervals below one and no cell >10% worse, confirmation and replay.
Strict rho beating additionally requires BOTH metric upper intervals <1 and
EVERY cell ratio <1 on both stages. Report this separately from IC promotion
and the existing <=1.10 parity rule. Never relax a gate after measurement.

Fresh seed 2026091608; pilot five-cell confirmation, 60 targets x 3 repetitions;
CPU 7, 8 GiB per child, 60 seconds per child, at most 1800 paired jobs. Build on
CPUs 0-1, and start measurement only after parent completion and preflight.
All failed trials, raw profiles, receipts and frozen sources are retained.
Stop after this bounded round for audit and publication, whatever the verdict.

# Separately declared cold single-target factor-base policy panel

Reason for the new panel: round-0006b development measures ~21% lower instruction
cost but ~11% lower native time for `combined_batch1`. The original support still
comes from eight sampled orbits at each batch boundary. Small one-target jobs
pay to construct and solve more columns than their requested point count needs.
The candidate changes the factor base deliberately; it is not admitted as a
fixed-support implementation comparison and does not revise round-0006b.

Reference workload: the same complete cold one-target public ECDLPs, five small
Koblitz cells, matched signed-Frobenius rho, all setup/verification/failed work
charged. Incumbent remains the last promoted single-target source/configuration.
Execution waits for the preceding full decision. Fresh seed 2026091607, 60 unseen
confirmation fixtures, three repetitions and replay, CPU 7, 8 GiB and 60 seconds
per child, <=1,800 profiled/native job pairs.

For each arm separately, B signed points and m=3 allow at most binomial(B+2,3)
target images. Changing B changes that coverage bound and the column/rank floor.
Report actual B and rank per arm, and complete cost / sqrt(r) against the same
rho reference. A smaller ratio to a changed floor is not an algorithmic advance.
Classification: engineering and parameter policy; no new exponent or family-wide
claim.

Candidates, fixed before inspecting any new selection/confirmation data:

- `cold_context`: reuse one fast field context during all factor-base lifts;
  preserve the original support. Its parent is round-0006b's measured combined
  implementation, which is retained as a candidate even if not promoted.
- `orbits1`, `orbits2`, `orbits3`, `orbits4`: sample one orbit at a time and stop
  at 2*n times the specified orbit count. Same seed/domain and exact point checks.
- `cube_root`: the same one-orbit sampler, requesting at least max(2*n,
  ceil(cuberoot(r/2))) points. This balances the rough folded setup term B^2/(2n)
  with the collection term r/(2n*B); the model is a proposal, not a measured bound.

The policy worker supports explicit orbit count and cube-root policy knobs and
records the effective base recipe in its complete report. A new sealed contract
declares `comparison_kind=factor-base-policy`. Fixed-support comparisons still
reject differing base fingerprints. Policy comparisons require stable support
within every arm/case, identical public targets across arms, valid independently
checked relations/logs/scalars and complete equal workloads.

Promotion: >=20% lower instructions AND native time, both paired 95% upper limits
below one and no cell worse by >10%, in confirmation and replay. Beating rho is
reported only if both upper 95% ratios and every cell ratio are below one on both
stages. Keep all intermediate and losing policies. Do not select a new policy
from final holdouts or replace an unsuccessful original panel with this result.

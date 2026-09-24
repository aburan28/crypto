# Frozen PDP feature experiment, 2026-09-24

Question: On the saved natural two-summand S3 instances, can the affine row-span contradiction screen save a complete root solve while retaining every useful group relation? Does the lower quadratic rank (equivalently target x in the product span) enrich useful hits among surviving instances?

This is a toy structural and diagnostic-cost experiment. It makes no ECC2K-130 solving-degree, full-DLP, or rho speed claim. The full producer and its original contract are preserved as profile_probe.py and pdp_profile_contract.md. The source dependency scripts/ecc2k130_point_decomposition.py is tracked on main at SHA-256 5259a42b3613e315835546cf3936a4e1dc77c3a7560f3da36dde4e1b233c26d5.

Inputs are the producer's 640 natural cases: n=7,9,13,17; five bases per field; 32 targets shared across bases. All 76 planted cases remain in the producer's audit. No natural case is dropped, and repeated target abscissae are clustered by their binary Frobenius orbit.

Frozen split: in each field, sort the five base IDs by SHA-256 of the UTF-8 string "pdp-base-holdout-v1|" plus base ID; hold out the first. Canonicalize each target x by the minimum of its repeated squares in F_(2^n), sort distinct canonical x values by SHA-256 of "pdp-orbit-holdout-v1|" plus n plus "|" plus canonical x, and hold out ceil(one third). The confirmation set is the intersection of held-out bases and held-out target orbits. The training set is the intersection of the complements. The two mixed cells are reported but cannot choose a rule. This deterministic split was specified before computing split-specific outcomes. The dataset already existed and its aggregate outcomes had been inspected, so this is a retrospective grouped holdout, not an untouched prospective test.

Policies are fixed without fitting:
- Baseline: solve all two-summand S3 systems by the linear-fiber solver.
- Safe screen: compute product-span, quadratic rank, and complete affine consequences in the original equation row span; skip the fiber solver only if the affine rows contain 1=0. A missed group hit or an algebraic root in a skipped case fails the experiment.
- Rank priority diagnostic: among nonrefuted cases, compare group-hit rates for x in the product span versus x outside it. This is a priority hypothesis, never a hard UNSAT decision.

For every natural case, the analysis recomputes the profile and complete fiber roots from the tracked code, verifies feature and algebraic-root counts against the saved producer, and verifies group labels on the double-held-out confirmation cases. It records field multiplications, squarings, affine elimination row XORs, fiber elimination row XORs, and per-case process CPU nanoseconds. It runs three paired repetitions in fixed order and retains all repetitions. Field setup, base generation, target generation, and group verification are common or audit work and are reported separately, not hidden in a claimed solver speedup. Native counters are not converted to curve additions.

Acceptance: zero correctness discrepancies; report both train and double-held-out hit retention, certified refutations, and total screen-plus-remaining-fiber cost against fiber-only in each native unit and CPU. A positive cost result is only a diagnostic on this saved corpus, because the run is not a paired fresh full-DLP benchmark. Reject the rank priority hypothesis if it does not have a consistent held-out enrichment across field degrees with hits. Preserve negative outcomes.

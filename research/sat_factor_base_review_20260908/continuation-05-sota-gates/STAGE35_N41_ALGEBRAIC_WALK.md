# Stage 35: walked two-summand descent on the algebraic n=41 base

Stage 35 tests whether the walked individual-logarithm descent from PR #197 improves the same algebraically discovered factor base and five public targets used by Stage 33. It does not substitute the pseudo-random subgroup-orbit factor base used by the newer crossover experiments.

The collection factor base remains the target-independent search over the six-mask Frobenius union `[1,2,4,8,16,32]` and its two-torsion saturation. Search receives 1,024 separate sample points, no future Stage-35 target, no subgroup enumeration, and no factor-base log labels. Relation collection uses three summands. Only the individual descent changes: it walks probes by `+G` and asks the existing pair table for two-summand decompositions.

The previous workflow validator capped `max_trials` at one million even though walked probes are deliberately cheap and the merged optimization documents a 100-million cap. This integration mismatch censored the exact public targets: the old Stage-33 cap of 200,000 solved one of five; one million solved four of five; the fifth completed after 1,054,786 probes. Stage 35 raises only the workflow validation ceiling and freezes the intended 100-million cap. It also emits the effective `descent_summands` in the workflow report so the two-summand result cannot be mistaken for the three-summand collection setting.

A local unbound feasibility run completed all five exact targets with 2,417,575 total walked probes. It used the same 4,759-point, 29-column algebraic base and 39 relation rows as Stage 33. Its online IC/rho wall ratio was 1.2069 in rho's favor, an improvement over Stage 33's 2.2494 ratio but not a crossover. This local number is diagnostic only; the hosted singleton-CPU run supplies the claimable result.

The hosted workflow charges a fresh dependency acquisition, four-job clean build, true one-CPU algebraic search, relation collection, sparse log solve, five public descents, and signed-Frobenius rho on the same targets. It retains process and process-tree memory, core-seconds, wall time, and the complete walked-probe counts. Pair-table descent exposes no SAT conflicts; Phase B remains the same-instance SAT conflict comparison.

Even if the online ratio crosses on a different host, this is a finite constant-factor result over a 40-bit subgroup. It does not supply licensed Magma, unaffiliated reproduction, an exponent change, or a Koblitz index-calculus SOTA result.

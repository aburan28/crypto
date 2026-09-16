# Cold single-target continuation

Objective: a complete cold one-target IC solve cheaper than matched signed-Frobenius
rho in both Valgrind amd64 instructions and native process time. All setup,
verification and failed work remain charged. No amortization across targets.

Incumbent: round-0004 `folded_lift_batch4`, with the exact frozen source/config.
Measured confirmation costs were 3.0429 times rho's instructions and 1.5213 times
rho's time. Factor-base/table setup costs 9–22 million instructions per cold job;
descent adds 3–9 million. These are the first bottlenecks to attack.

Boundaries remain fixed: signed base B and m=3 give at most binomial(B+2,3)
target images. The weak rank/instruction floor and matched rho reference retain
the existing accounting contract. All proposed changes are engineering, with
identical factor-base points, ordering, relations and complete certificates.

Round 0006 uses fresh seed 2026091606, one target, the existing five-cell panel,
60 fresh confirmation fixtures, three repetitions, replay, CPU 7, 8 GiB and
60 seconds per child, at most 1,800 paired profiled/native jobs.

Candidates:

1. `descent`: test round-0005's lazy first probe and fast internal scalar checks
   on the single-target workload for the first time.
2. `lazy_field`: avoid constructing polynomial field structure for the pair-table
   and enumeration strategies, which do not use it. Algebraic strategies retain it.
3. `fast_orbits`: construct the same ordinary and signed orbit tables using exact
   single-word points and keys; retain the general implementation as fallback and
   an equivalence oracle. Expected memory decreases by avoiding large-integer keys.
4. `combined`: measure these mechanisms together in a complete solve.
5. `combined_batch1`: test reduced surplus collection on the combined source.

Promotion keeps the existing requirement: >=20% reductions in instructions and
native time, paired 95% upper limits below one, no cell worse by >10%, in both
confirmation and replay. Separately, **beating rho** requires both paired 95%
upper ratios below one and every cell ratio below one in both final stages.
The existing <=1.10 parity flag alone does not establish beating rho.

Reject any changed base support, incorrect orbit mapping, bad scalar or missing
cost. Preserve losing candidates, failures and timeouts. Stop this bounded round
if these mechanisms fail; any further round requires a new measured hypothesis
and fresh confirmation seed. No mathematical or cryptographic-size claim.

## Preparation correction

The first preparation completed all five builds but rejected CPU 7 because the
controller was restricted to CPUs 0 and 1. It produced no fixtures or measured
trials. Preserve that attempt as `round-0006-single`; `round-0006b-single` uses
the declared CPU 7 with an allowed affinity mask. Reuse the completed binaries
only after rehashing all 450 source files against each requested candidate and
verifying exact binary copies. `build-reuse.json` records each provenance link.
The seed, candidates, workload, costs and gates are unchanged.

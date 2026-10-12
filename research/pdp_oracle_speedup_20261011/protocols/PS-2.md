# PS-2: scalar reference kernel with incremental updates

Part of [`../PLAN.md`](../PLAN.md); frozen when the plan merges.  Repo: crypto.
Needs PS-1.

**Hypothesis.**  Updating the quartic's coefficients incrementally along a
Gray code over one summand, and testing roots with PS-1's cheaper
formulation, is at least twice as fast as the baseline on one core before any
vectorisation, and its gate count `G` per pair is small enough for
bit-slicing to reach 5× (`G ≤ 1.8·10^5`).

**Frozen inputs.**  The PS-0 corpus and baseline; PS-1's formulas.

**Method.**  A new module `semaev_gray_pairs` in `src/cryptanalysis/` with a
scalar kernel and a `cfg`-free gate-count export (the straight-line program
as data, hashed), tested against the corpus and `semaev_decomp`.

**Unit.**  ns per pair on an isolated core, counted bit operations per pair,
`G`, and the SHA-256 of the verdict and root-set stream.

**Success.**  Digests equal PS-0's on every rung; scalar ns per pair at most
half the baseline's; `G` reported with its derivation.

**Stop.**  `G > 1.8·10^5` per pair (PLAN §5 model: bit-slicing cannot reach
5×): PS-3 does not start and the finding is recorded.  Scalar slower than
the baseline alone does not stop the program, because bit-slicing changes
the cost model, but it must be reported.

**Inadmissible.**  Changing the corpus; counting setup (table build) outside
the timed region without reporting it; comparing against a different
compiler profile than the baseline's.

**Class.**  Engineering.  **Effort.**  3 days.

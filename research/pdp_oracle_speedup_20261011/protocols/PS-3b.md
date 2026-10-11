# PS-3b: port the CPU kernel into the cryptanalysis suite

Part of [`../PLAN.md`](../PLAN.md); frozen when the plan merges.  Repo:
cryptanalysis (paired with PS-3).  Needs PS-3 merged in `crypto`.  The
companion protocol in that repository is
`experiments/pdp-oracle-speedup-20261011/README.md`.

**Hypothesis.**  The kernel can be ported to `suite/src/cryptanalysis/` with
the corpus digests unchanged, even though `semaev_decomp.rs` has diverged
between the two repositories (about 1,500 diff lines).

**Frozen inputs.**  PS-3's merged commit, PS-0's corpus and digests, the
`suite`'s own test suite and CI.

**Unit and success.**  The PS-0 digests reproduce from the port; the suite's
tests and CI pass; per-core pairs per second within the A/A spread of PS-3's
on the same host class.

**Stop.**  A faithful port needs the other repository's oracle module
refactored first: stop, report the minimal refactor as its own task, do not
fork the algorithm.

**Class.**  Engineering (no new claim).  **Effort.**  2–3 days.

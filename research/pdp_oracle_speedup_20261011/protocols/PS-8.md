# PS-8: the meet-in-the-middle pair table against the pairs-and-solve kernels

Added after PS-0…PS-7 were frozen, as a **new lead** produced by their results; it changes none of them
([`../LEADS.md`](../LEADS.md) records where it came from).  Frozen when this file is committed, before any PS-8
measurement.  Repo: crypto.  Needs PS-2 (the group-law helpers and the harness); uses PS-3's binary only as a context arm.

**Why.**  PS-2 found that the pairs-and-solve test for `m = 3` is a group-law computation: the candidate third
abscissa of the pair `(P₁, P₂)` is `x(P₁ ± R ± P₂)`.  `src/cryptanalysis/ic_boundary.rs` already carries the oracle that
computes the same thing the other way round: `Oracle::Mitm` builds the table of `x(Pᵢ ± Pⱼ)` once and probes it once per
`R ∓ P_k` (`decompose_mitm`).  The table's keys do not depend on the target; the kernels recompute the equivalent
quantities for every target.  No measurement in this program puts the two on one corpus: PS-6 compares kernels with one
another.  PS-8 measures the table in the same harness, so that a speedup of the pairs-and-solve stage is read against the
repository's own alternative and not only against the loop it replaces.

**Observations, not results (from reading the code, before measuring).**  The table holds `|B|²` keys for `|B|` base
abscissae (`|B|(|B|+1)/2` sums and `|B|(|B|−1)/2` differences) and a target costs `2|B|` probes, against
`|B|(|B|+1)/2` pairs and four candidates per pair for the kernels.

**Hypotheses (each can fail).**
H1. On the whole PS-0 corpus the table arm returns the verdict digest and every root set of PS-0, exactly as `group-scalar`.
H2. At `ℓ ≥ 10` the per-target probe time is at most 1/10 of `group-scalar`'s.
H3. The table's build time divided by `(t_group − t_probe)` (the break-even number of targets `N*`) is at most 100 at every rung `ℓ = 8…12`.
H4. Table memory is at most 40 bytes per key at every rung measured.

**Frozen inputs.**  PS-0's corpus and digests (`ℓ = 8…12`); `ℓ = 13` timing sample of PS-3 (3 refuted targets, corpus in
`ps3_bitsliced/raw/corpus_l13.txt`) if its verdicts are first replayed by the table; the PS-0 targets for timing (the
first 16 refuted targets, 6 at `ℓ = 12`, 3 at `ℓ = 13`); rounds ≥ 5; `tools/isolated_bench.py run`, one pinned core.

**Method.**  Arm `mitm` in `examples/pdp_speedup` (`mitm.rs`): open-addressing hash table keyed by the abscissa, built once
per run from the liftable base points of `V` (the same `B` as the group kernels); probe = `x(R + P_k)` and `x(R − P_k)` for
every `k`, denominators inverted by Montgomery's trick, every table hit confirmed with the exact reference (`group::thirds`).
Targets with `x_R ∈ V ∪ V⁻¹ ∪ {0}` go to the baseline, as in PS-2 (same exactness theorem).  Arms of the ABAB comparison:
`baseline`, `baseline2` (A/A), `group-scalar`, `bs512` (PS-3's widest vector arm, as context), `mitm`.

**Unit.**  Wall nanoseconds per refuted target for the probe stage; **separately and always beside it**: table build seconds
(`setup_secs` of the same runs), keys, bytes; `N* = build / (t_group − t_mitm)` and the same for `bs512`, with the arithmetic
shown.  `ns per own pair` is not reported for `mitm` (its unit is a probe; `pairs_per_target` is `2|B|`).

**Success.**  H1 holds at every rung (identical digests and root sets).  The remaining hypotheses are reported as held or
failed with their measured values; a failed hypothesis is a result, not a reason to retune.

**Stop.**  Any H1 discrepancy: report the first differing target, keep the arm out of every timing table, and treat the arm as
unverified.

**Inadmissible.**  A probe time without the build time and the memory beside it; a break-even that amortises over an
unstated target count; comparing the table with itself (a target whose own sums are in the table is not a probe of the
table: targets are PS-0's, not base-point sums); pooling contended and uncontended rounds; any statement about an end-to-end
or ECC2K-130 cost.

**Not claimed.**  The table is `Θ(|B|²)` memory.  At the sizes where the program's open question lives (`ℓ ≈ 44` for
ECC2K-130) it is `2^{86}` keys and is out of reach; PS-8 says nothing about that regime.  Relation yield, the number of
relations needed, linear algebra and the charge in the counted unit `S` belong to PS-6; PS-8 is a stage measurement on one
host class.

**Class.**  Engineering measurement.  **Effort.**  1 day.

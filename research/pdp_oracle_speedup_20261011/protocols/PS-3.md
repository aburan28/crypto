# PS-3: bit-sliced CPU kernel

Part of [`../PLAN.md`](../PLAN.md); frozen when the plan merges.  Repo: crypto.
Needs PS-2.

**Hypothesis.**  Evaluating PS-2's straight-line program over 64, 256 and 512
pairs per register, with the 512-lane form using `vpternlogd`, reaches at
least 10× the baseline per core at `ℓ = 11` on an AVX-512F/BW host, and
scales to the host's four cores at 85% efficiency or better.

**Frozen inputs.**  PS-0 corpus and baseline; PS-2's program (hash recorded);
rungs `ℓ = 8…12` plus `ℓ = 13` for timing only (decided against the
unreduced loop on a sample).

**Method.**  Portable `u64` path; AVX2 and AVX-512 paths behind
`is_x86_feature_detected!` and `#[target_feature]`, with a runtime dispatch
and no global `target-cpu`; threads by `std::thread`; one-time table builds
reported separately.

**Unit.**  Pairs per second per core and counted bit operations per pair;
instructions retired under `callgrind`; wall figures only from
`tools/isolated_bench.py` runs, ABAB with at least five rounds, with the A/A
spread; thread counts 1 and 4.

**Success.**  Identical digests on the whole corpus at every width and thread
count; ≥ 10× the PS-0 baseline per core at `ℓ = 11` outside the A/A spread;
four-thread efficiency ≥ 85%; no regression of the single-thread figure when
threads are enabled.

**Stop.**  Below 5× after tuning at all three widths; record the best width
and why.  Hardware classes not measured (Arm64 NEON, Metal) are named, not
implied.

**Inadmissible.**  Reporting the 512-lane figure without the AVX2 and
portable ones; pooling contended and uncontended runs; a global `target-cpu`
build as the headline.

**Class.**  Engineering.  **Effort.**  5–8 days.

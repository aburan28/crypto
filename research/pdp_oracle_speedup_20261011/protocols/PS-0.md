# PS-0: freeze the baseline, host manifest and verdict corpus

Part of [`../PLAN.md`](../PLAN.md); frozen when the plan merges.  Repo: crypto.

**Hypothesis.**  None.  This is accounting: it fixes what every later task is
compared against.

**Frozen inputs.**

- Reference: `examples/semaev_decomp_bench` pairs-and-solve (one-word field,
  `gcd(q, L_V mod q)`), built from the `main` commit recorded at run time,
  binary kept.  Any instrumentation needed to count field operations is a
  separate build; if it changes the timed code, both builds are recorded.
- Rungs `ℓ ∈ {8, 9, 10, 11, 12}` with the bench's own moduli (`n = 3ℓ`).  The
  curve is named by its ICV1 slug once PS-0 records it; no slug is cited
  before it is registered.
- Corpus per rung: 64 uniform targets and 16 planted decomposable targets,
  SplitMix64 seeds `20261011 + ℓ`; every target decided exhaustively and
  every decomposition re-added in the group.  The corpus stores targets,
  verdicts and full root sets.

**Unit and accounting.**  Per pair: `callgrind` instructions retired at
`ℓ = 8, 9`; counted field multiplications, squarings and inversions; ns per
pair single-threaded on a core reserved with `tools/isolated_bench.py`,
reported with the number of contended runs.  A/A (baseline against a copy of
itself) and ABAB interleaving, at least five rounds.

**Success.**  The A/A spread is recorded (a spread above 10% is reported, not
hidden, and widens every later gate by that spread); the SHA-256 of the
verdict and root-set stream per rung is identical in three runs.

**Stop.**  None: PS-0 cannot fail, only record.

**Inadmissible.**  Reusing the unpinned numbers of PLAN §2 as the baseline;
dropping contended runs; changing the corpus after a candidate has run.

**Outputs.**  `research/pdp_oracle_speedup_20261011/baseline_v1/` (host
manifest, binary hash, commands, raw runs, corpus, digests, A/A summary).

**Class.**  Accounting.  **Effort.**  1 day.

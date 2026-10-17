# Protocol: PKM's isogeny-tower factor base in ecbench, with its auxiliary-curve search charged

Registered 2026-10-15, before any measured session of this directory ran.
The code (`prime-isogeny:t=T[,seed=S]` in
`src/cryptanalysis/ic_framework/prime_tower.rs`), the curve search
(`examples/prime_tower_curves.rs`), the factor-base sizing reads
(`ecbench fb`, below) and `ecbench plan` on the specs precede it; every
measured run follows the commit of this file and the specs.

## Question

Lead 2 of the factor-base encyclopedia (`docs/ic/factor-base-encyclopedia/README.md`
§8): the isogeny form of Petit–Kosters–Messeng's tower base (PKM note
`research/notes/index-calculus/RESEARCH_PKM_TOWER_ORACLE.md` §2.1),
`V = x(R' + ⟨Q⟩)` for a point `Q` of order `2^t` on an auxiliary curve
`E'/F_p`, the zeros of a chain of `t` Vélu 2-isogeny x-maps, is the only
tower shape that needs no condition on `p`; it is the shape PKM attach to
P-192 (`v₂(p − 1) = 1`).  Until now it existed only in
`examples/pkm_tower_pilot.rs`, with naive `O(p)` point counting, and the
pilot reported but never charged its auxiliary-curve search.  This round
asks what index calculus with that base costs end to end in ecbench's unit
`S` when the search is charged, against the same index calculus on a plain
base of equal size, on the Kummer tower (where it exists), and against
Pollard rho with the negation map.

## What is charged (the instrument)

`prime-isogeny` draws `(a', b')` from a SplitMix64 stream seeded by
`(p, t, seed)` and, on every curve,

1. **filter**: one random point `P` and a baby-step giant-step over the
   multiples of `2^t` in the Hasse interval `[p + 1 − ⌊2√p⌋, p + 1 + ⌊2√p⌋]`
   for `[k]P = O`; no hit proves `2^t ∤ #E'` and drops the curve;
2. **count** (survivors only): `#E'` by Mestre's baby-step giant-step
   (each point narrows the class of `#E'` until one multiple is left in the
   Hasse interval);
3. `Q = [#E'/2^t]P'` until `[2^{t−1}]Q ≠ O`;

then `⟨Q⟩` (`2^t − 1` additions), `R'` with `2R' ∉ ⟨Q⟩` (one doubling a
draw), `V` (`2^t` additions) and the chain (`K_j` by doublings, Vélu
steps).  Every addition and doubling on `E'` (a curve over the same `F_p`,
same affine formulas as `E`) is on the factor-base phase's ledger and so in
`S`.  The record's `factor_base` phase also carries, as counters,
`aux_curves_tried`, `aux_search_gae` (filter + counts, every curve),
`aux_count_gae` (counts of the survivors), `aux_tower_build_gae` (`Q`,
`⟨Q⟩`, `R'`, `V`, chain), the rejection reasons, and the field work outside
the group law: `aux_legendre_symbols`/`aux_sqrt_solves` (also added to the
phase's `legendre_symbols`/`sqrt_solves`, unpriced as for every prime base),
`aux_bsgs_table_entries` (memory), and `velu_evaluations_uncharged`
(`2t` inversions).  The factor base is every point of `E` whose abscissa is
in `V`, both signs, one column an abscissa, exactly as `prime-kummer` lifts
its `V` (shared code).

## Curves

**Panel 1** (two-adicity 12): the five curves of
`research/pkm_tower_ecbench_20261013` (`curves.jsonl` there, copied here as
`curves-panel1.jsonl`), registered slugs `icv1-fp20-tm728-d215a2a1`,
`icv1-fp22-tm2947-e15fed69`, `icv1-fp24-t2705-27cbb01b`,
`icv1-fp26-t2838-8e779ca8`, `icv1-fp28-t4155-50f492b1`; the same `t`
(`8, 9, 9, 10, 11`).

**Panel 2** (the P-192 shape, `v₂(p − 1) = 1`, so no Kummer tower):
`prime_tower_curves --bits 20,22,24 --two-adicity 1 --seed 20261015
--max-cofactor 4` (`curves-panel2.jsonl`, deterministic), registered in
`docs/curves/registry.json` before they are cited; `t` by the same rule
`t = ⌈(bits + log₂ 6)/3⌉` (`8, 9, 9`).

## Arms

`seed = 0` for every `prime-isogeny` arm.  `S_plain` = the abscissa count
of the `prime-isogeny:t=T` base on that curve, read with `ecbench fb`
before the runs (a sizing read, not a measurement; listed in `sizing.txt`,
which also shows each curve's `seed = 0` search outcome: 3 to 744 auxiliary
curves tried, the build being deterministic, this is the value every run
of that curve will charge).  The search outcome is one draw per curve from
a wide distribution; panel 3 measures that distribution, and the
per-curve bootstrap intervals (over targets) do not include it.

| arm | role | method | panels |
|:--|:--|:--|:--|
| `rho-neg` | reference | `rho.negation` | 1, 2 |
| `ic-iso-rootfind` | candidate | `ic.pipeline`, `prime-isogeny:t=T`, `rootfind-mitm:m=3` | 1, 2 |
| `ic-iso-mitm` | candidate | `ic.pipeline`, `prime-isogeny:t=T`, `mitm:m=3,negation_folded=1` | 1, 2 |
| `ic-kummer-rootfind` | comparison | `ic.pipeline`, `prime-kummer:t=T`, `rootfind-mitm:m=3` | 1 |
| `ic-plain` | control (baseline) | `ic.pipeline`, `prime-abscissa:size=S_plain`, `rootfind-mitm:m=3` | 1, 2 |
| `rho-neg-aa` | negative control (A/A) | `rho.negation` | 1, 2 |

`max_trials = 200000`, `linalg = incremental-gauss` (default), eight
one-target workloads a curve (panel 1: the 2026-10-13 study's
`target_seed` and measurement `seed`, so the targets are its targets;
panel 2: `202610150100 + bits` and `202610153100 + bits`), **two** measured rounds, `warmup = 1`,
`order = alternate`, `isolation_required = L0`, 3,600 s a run.

**Scale-down, stated in advance.** The 2026-10-13 study ran three measured
rounds at `L2`.  Here: two rounds (to keep each sealed `records.jsonl` under
1 MB, the lane's file cap, so sessions can be committed and CI can replay
them; CPU is not the constraint, the earlier study took four minutes), and
`--cpus none --allow-busy` (L0): the host is shared with four other lanes,
so **no wall-clock figure is admitted**; every claim is in counted `S`.

**Panel 3** (search-cost calibration, not an ecbench session): `ecbench fb
--factor-base prime-isogeny:t=T,seed=S` on the 24-bit curve of each panel,
`t = 4 … 14`, `S = 0 … 15`, recording `aux_curves_tried` and the charged
group operations; fitted as `curves_tried ≈ C·2^t` for the projection
below.

## Hypotheses and predictions

P1. Every run verifies and the replay audit (`--replay-all`) is `OK` with
    every replay identical, on every curve.
P2. The isogeny base's x-set is `V` (tests: `|V| = 2^t`, every element a
    root of the composed Vélu x-map, and on small `p` nothing else is), and
    the run's `fb_id` equals the `ecbench fb` dump's.
P3. `aux_curves_tried / 2^t` is of order 1 (between 0.25 and 8) at every
    `t`; the charged search per curve grows like `√(4√p / 2^t) + log p`.
P4. `ic-iso-rootfind ÷ ic-plain` in `S` is **above 1** on every curve with a
    95 % paired-bootstrap interval excluding 1, and equals
    `1 + (charged factor-base GAE)/(plain total)` within the interval: the
    plain base has no search; the tower structure buys a table oracle nothing
    (2026-10-13 P3: Kummer ÷ plain within 1.2 %).
P5. The aux-search share of the isogeny arm's total `S` (mean of
    factor-base GAE ÷ total GAE over its runs) falls with `p` at the
    PKM `t`: about `p^{1/3}·(p^{1/12} + log p)` against the pair table's
    `p^{2/3}`.
P6. `ic-iso-rootfind ÷ ic-kummer-rootfind` (panel 1) is above 1 by about the
    aux share (bases of the same `t`; their lifting abscissa counts differ
    by sampling).
P7. Every IC arm is above rho (`S_IC/S_rho > 1`), as in 2026-10-13.
Negative control: `rho-neg-aa ÷ rho-neg` = 1.000 in `S`.
Positive control: `ic-iso-mitm` vs `ic-iso-rootfind` reproduces the
2026-10-13 rootfind ÷ mitm ratio (1.03–1.21, intervals excluding 1),
showing the instrument resolves differences of a few per cent on these
curves.

## Decision rule

- Audit `OK` on a curve, else no result on that curve.
- Per curve, from `ecbench compare` (`Σ S_B / Σ S_A`, 95 % cluster
  bootstrap): **isogeny ÷ plain** and **isogeny ÷ kummer**.  Verdict
  "the charged search is a measured cost" on a curve when the isogeny ÷
  plain interval lies above 1; "not resolved" when it contains 1;
  "isogeny cheaper" when it lies below 1 (which would mean an accounting
  error to explain, since the plain base has no search).
- The **aux-search share** of total `S` per curve, with its range over
  runs; P5 holds when the share falls from the smallest to the largest
  panel-1 curve.
- The P-192 projection is derived, never measured: `curves_tried ≈ C·2^t`
  with `C` from panel 3, times the per-curve cost of (a) this BSGS filter
  at `log₂ p = 192`, `t = 48`, against rho's `√(πr/4)`; SEA, CM
  construction and 2-Sylow halving are named and not implemented.
- Class per AGENTS.md §3: expected **accounting** (a cost the pilot left
  out is now in `S`) for the isogeny arm against the Kummer arm, and
  **engineering** for nothing; no exponent claim is made (sizes 20–28 bits).

## Inadmissible

Changing curves, `t`, seeds, arms, rounds or accounting after seeing
results; quoting wall time (L0); pooling sizes into a scaling claim;
reporting the isogeny-vs-plain ratio without the aux share beside it.

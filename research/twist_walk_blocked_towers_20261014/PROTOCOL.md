# Protocol: can a short 3-isogeny walk leave a blocked twist tower?

Version 1, 2026-10-14. Frozen before the catalog and null runs; the sha256 of
this file is recorded in `manifest.json` next to the results. Skill:
`.claude/skills/controlled-experiment/SKILL.md`. Instrument:
`examples/twist_walk_blocked_towers.rs` (native Rust). Conductor is not
installed in this environment, so no task was attached.

## Question and claims

MMT Lemma 6 condition (8) blocks the GHS descent of a tower `(N, l, n = N/l)`
when `Tr_{K/F_2}(a) = 1`, `n` is even and `Tr_{K/k}(b) = 0`. On the 91 binary
catalog curves no tower is blocked as recorded. On the **quadratic twist** (same
`b`, flipped `Tr(a)`), 13 towers on 5 curves are blocked (the pilot in the experiment-proposals PR,
the pilot; the frozen list is `FROZEN_TOWERS` in the example).

Question: from each of those 5 curves, does a walk of at most 4 steps in the
`Φ_3 mod 2` isogeny graph (Frobenius `j ↦ j²` free) reach a curve whose twist
tower is **not** blocked, and with what magic number `m'`?

| claim | class (§1 of the skill) | licensed by |
|:--|:--|:--|
| "tower T of curve C has an unblocked 3-isogenous twist within depth `d`" | exact computation | the BFS result for T, reproducible from the registry |
| "the first such depth is `d*` and the least `m'` there is `m*`" | exact computation | same |
| "random blocked starts in the same `(N, l)` escape within depth 4 at rate `p0`" | exact computation on a sample, with a Clopper–Pearson interval | the null rows |
| "catalog curves are unusual / typical against random starts" | comparative | **not made**: see Decision rule |

## Instrument

- Edges: roots in `F_{2^N}` of `Φ_3(X, j) mod 2`, from the repository's tabled
  `Φ_3` (`binary_isogeny::phi_l_mod2_in_x`) and its root finder. Nodes: Frobenius
  orbits of `j`; the start is depth 0. Degree-2 steps are Frobenius in
  characteristic 2 and are not walked; they preserve `m` and whether
  `Tr_{K/k}(b) = 0`.
- Per node `b = 1/j`: `m' = ec_trapdoor::magic_number` (MMT value), and the twist is
  unblocked iff `ec_trapdoor::mmt_condition_8(tower, a_tw, b)` with `a_tw` a
  `z^i` of absolute trace 1.
- Depth cap 4. For each tower: nodes per depth, unblocked nodes per depth,
  least `m'` among unblocked nodes within depth `≤ d`, first unblocked depth.
- Same instrument, depth cap and metrics on the catalog and on the null.

## Controls

- **Positive control A (root finder, run before anything else).** At `N = 20`
  (the Cantor–Zassenhaus path) the roots of `Φ_3(X, j)` equal an exhaustive
  scan of all `2^20` elements that evaluates `Φ_3` from its tabled monomials, for
  8 values of `j` that have at least two roots. **Result: PASS** (46 roots over
  23 draws, 8 draws with 4 roots), `control_out.txt`.
- **Positive control B (walk and evaluation).** At `N = 18, l = 3`, `walk` +
  `eval_tower` equals the same BFS fed by exhaustive neighbour scans, on 10 random
  blocked starts with at least one neighbour (2 to 17 nodes; first unblocked depth
  1, 2 or none). **Result: PASS**, `control_out.txt`.
- **Null population.** For each distinct `(N, l)` among the 13 towers (12
  configurations): `M = 2000` random blocked starts. `b` uniform, then forced to
  `Tr_{K/k}(b) = 0` by subtracting `c·w0` with `Tr(c·w0) = Tr(b)` (the
  construction asserts the result is blocked). Same field modulus as the catalog
  curve of that degree. The null matches everything except the property under
  test (`b` has no special structure). Seed `20261014`.
- The catalog `b` are structured (for example `c2pnb208w1` has magic number 7 at
  `n = 104`); a walk neighbour is a different, unstructured `b'`.

## Decision rule

None of the catalog outcomes is tested against the null. Reason (power check
for the comparative claim, `2/(N+1) < α/m`): the catalog supplies 5 independent
`j`-graphs (the 13 towers share graphs by curve) and a binary outcome per tower,
so with `m = 13` corrected tests (`α = 0.05`) a catalog outcome can reject only if
the null rate `p0` is below `α/m = 0.0038` or above `0.996`; for any moderate
`p0` the design cannot report anything but "within the null". We therefore report
the catalog outcomes as exact computations and the null rates `p0` with 95%
Clopper–Pearson intervals beside them, and **write no sentence that ranks the
catalog against the null.** Verdict wording: "unblocked at depth d (m' = x)",
"no unblocked twist within depth 4", "no rational 3-isogeny from this curve".

## Pre-flight runtime estimate

Dry run on scratch output, null population only (the catalog was not run), 3
starts per configuration: 0.4 s in total, so about 11 ms per start, dominated by
the `N = 368` configurations. Estimate for `M = 2000`: `12 · 2000 · 11 ms ≈ 4.5`
minutes; catalog (13 towers on 5 graphs): under 1 s. Record the actual time in
`manifest.json`.

## What does not count

Changing `M`, the depth cap, the null construction or the tower list after seeing
results. Any such change is a labelled new protocol version, and both are reported.

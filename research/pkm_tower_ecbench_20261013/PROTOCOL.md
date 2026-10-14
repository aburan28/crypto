# Protocol: Petit–Kosters–Messeng's prime-field index calculus in ecbench, with the root-finding oracles, against rho

Registered 2026-10-13, before any session of this directory ran.  The
curve search (`examples/prime_tower_curves.rs`) and `ecbench plan` on the
specs precede it; every measured run follows the commit of this file and
the specs.

## Question

PKM (PKC 2016, §3) replace the binary subspace by the zero set of a
composed low-degree map over `F_p`; in its Kummer form, `V = μ_{2^t}`, the
zeros of `x^{2^t} − 1`, and the factor base is the points above `V`.  They
solve the decomposition system with Gröbner bases and leave open whether
anything better exists.  `research/pdp_root_finder_20261012` measured the
root-finding oracles on PKM's own cells at 10⁵–10¹² times fewer field
operations than `f4_fp_tower`.  This round asks what the whole index
calculus costs end to end on prime curves with that base and those
oracles — factor base, oracle set-up, relations, linear algebra,
verification, in ecbench's unit `S` — against Pollard rho with the
negation map on the same targets, over five field sizes.

## Curves

`y² = x³ + ax + b` over `F_p`, `p = c·2^{12} + 1` of 20, 22, 24, 26 and 28
bits, `#E = h·r` with `r` prime and `h ≤ 4`, found by
`prime_tower_curves --bits 20,22,24,26,28 --two-adicity 12 --seed 20261013
--max-cofactor 4` (deterministic), entered as `prime_explicit` curves and
registered (AGENTS.md §11) before they are cited.

## Base and arms

Tower height `t = ⌈(bits + log₂ 6)/3⌉` (`8, 9, 9, 10, 11`): the smallest
`t` at which `|F|³/6` sums of three signed base points can cover the group,
the PKM choice `|V| ≈ p^{1/m}` at `m = 3`.

| arm | role | method |
|:--|:--|:--|
| `rho-neg` | reference | `rho.negation` |
| `ic-mitm` | baseline | `ic.pipeline`, `prime-kummer:t=T`, `mitm:m=3,negation_folded=1` |
| `ic-rootfind` | candidate | `ic.pipeline`, `prime-kummer:t=T`, `rootfind-mitm:m=3` |
| `ic-rootfind-pairs` | candidate, `t ≤ 9` | `ic.pipeline`, `prime-kummer:t=T`, `rootfind-pairs:m=3` |
| `ic-abscissa` | control | `ic.pipeline`, `prime-abscissa:size=S` with `S` = the tower base's abscissa count, `rootfind-mitm:m=3` |
| `rho-neg-aa` | control | `rho.negation` (A/A) |

`max_trials = 200000`, `linalg = incremental-gauss`, eight one-target
workloads a curve, three measured rounds, `warmup = 1`, `order =
alternate`, `isolation_required = L2`, 3,600 s a run.

No Gröbner arm: `f4_fp_tower` is not an `ic.pipeline` oracle, and its per-target
cost on the PKM cells is already frozen (rounds 2–7, re-read in
`research/pdp_root_finder_20261012`, where the root-finder costs
`10⁵`–`10¹²` times fewer operations).  The `ic-abscissa` control asks
whether the tower structure itself buys anything for an oracle that does
not use it: PKM's structure helps an algebraic solver, not a table.

## Predictions

P1. Every IC arm verifies every run; the audit replays every run
    identically.
P2. `S(ic) / S(rho)` is above 1 on every size and grows with `p`: the
    MITM index calculus costs about `p^{2/3}` group operations at `m = 3`
    against rho's `p^{1/2}` (PKM's own family floor is at or above rho for
    `m ≤ 3`).  A fitted slope of `log₂ S(ic)` against `log₂ p` near
    `2/3 − 1/2 = 1/6` per bit-doubling unit, that is `S_ic/S_rho ∝ p^{1/6}`.
P3. `ic-abscissa` and `ic-rootfind` agree in `S` within 20 % on every size:
    the tower is irrelevant to a table oracle.
P4. `ic-rootfind` is below `ic-mitm` in wall time and within 1.3× of it in
    `S` (the binary ladder's reading).

## Decision rule

Audit `OK` on every size, else no result on that size.  The per-size table
reads the paired ratios; the growth claim needs the fit over at least four
sizes with its interval.  Class per AGENTS.md §3: expected
**engineering**; an exponent signal only if `S(ic)/S(rho)` falls with `p`
with an interval excluding no change.

## Inadmissible

Changing curves, `t`, seeds, arms or accounting after seeing results;
pooling contended or below-level wall runs; quoting a lower-bound ratio as
a speedup.

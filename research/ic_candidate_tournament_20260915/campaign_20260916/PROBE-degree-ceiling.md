# Probe: can the crossover be measured at all?

Exploratory, not a round. No pre-registration, no promotion, nothing measured
here is admissible as evidence about index calculus. It is recorded because it
answers the question rounds 0019 and 0020 left open — and answers it *no*.

## The question

Rounds 0019 and 0020 hold `beats_rho_strict` on the eight-cell panel, whose
largest subgroup is `n23a1` at `r = 4,196,903`. The boundary in
ROUND19-single-target.md §3 puts this collector at `Θ(r^{2/3})` against rho's
`Θ(r^{1/2})`, and the same-degree contrast measures the ratio rising 16.6% per
doubling of `r` at degree 23. That is an extrapolation from two points, and
the next-proposal named a ninth cell beyond `r = 8·10⁶` as the measurement
that would settle it.

## Three new cells exist

Censusing degrees 23–41 with the frozen worker returns "worker accepts bounded
odd-degree Koblitz fixtures (5..31)" for everything above 31. Lifting that one
dispatch bound to 41 in a probe build and re-censusing:

| cell | r | cofactor | vs `n23a1` |
|:--|--:|--:|--:|
| `n39a0` | 68,616,367 | 8,012 | 16× |
| `n37a0` | 230,603,167 | 596 | 55× |
| `n41a0` | 549,756,390,943 | 4 | 131,000× |

(`n33a0/1`, `n35a0/1`, `n37a1`, `n39a1`, `n41a1` return "curve has no usable
subgroup" — mathematics, not a guard.)

## But the promoted candidate cannot run there

`koblitz_tiny_ic.rs` declares `MAX_DEGREE = 31`, and its `supports()` requires
`kc.n <= MAX_DEGREE` **and** `is_prime(kc.n)`. The worker's `(5..=31)` dispatch
bound is not an arbitrary guard sitting on top of a capable implementation; it
matches the tiny single-word module's own ceiling. Above degree 31 the worker
silently falls back to a general sparse path.

That fallback is visible in the report shape. At `n23a1` a solution carries a
`relation` — the descent certificate — and the base is reported as
`factor_base_orbits` (8 entries). At `n39a0` the solution carries only
`recovered` and `trials`, the base is a flat `factor_base` of 624 points, and
a `sparse_report` appears. `oracle.py`, with its own degree bound lifted to
match, accepts the degree and then rejects the certificate:

```
missing descent relation: logarithm not certified as index calculus
```

**This is the finding.** The promoted candidate — the executable rounds 0017
through 0020 measured, `e9f263b8…` — is structurally capped at degree 31. The
crossover cannot be measured with it, at any `r` above the panel's, because it
does not run there.

## What the fallback numbers are, and are not

Native wall, one fixture a cell, probe build:

| cell | r | IC path | IC wall | rho wall | ratio |
|:--|--:|:--|--:|--:|--:|
| `n23a1` | 4,196,903 | tiny, certified | 0.001s | 0.001s | 0.84 |
| `n39a0` | 68,616,367 | general, **uncertified** | 0.013s | 0.002s | 5.46 |
| `n37a0` | 230,603,167 | general, **uncertified** | 0.015s | 0.003s | 5.81 |
| `n41a0` | 549,756,390,943 | general, **incomplete** | 0.398s | 0.044s | — |

The 5.46 and 5.81 are **not** a measurement of the promoted candidate's
crossover and must not be read as one. They compare a different solver, on a
different code path, producing no certificate this campaign accepts. They are
recorded only to say what the fallback does.

At `n41a0` neither arm completes within the frozen `max_trials`, so that row
carries no comparison at all — an incomplete workload is not a cost, and by
the campaign's own rule a timeout is never negative mathematical evidence.

## What this changes

The next-proposal's "ninth cell beyond r = 8·10⁶" is **not reachable** as
written. Measuring the crossover needs the tiny implementation extended past
one-word arithmetic first: `MAX_DEGREE = 31` exists because the packed pair
table stores coordinates as `u32` and the field is held in one `u64`. Degree 37
and 41 are prime and would qualify on every other count; degree 39 is not
prime and would still be refused.

So the honest state of the campaign is:

- The strict win over rho is established and replicated, at subgroup orders up
  to `4·10⁶`, for an implementation that **cannot be run above degree 31**.
- The `Θ(r^{2/3})` against `Θ(r^{1/2})` boundary remains derived, and its only
  measured support remains the same-degree contrast at degrees 19 and 23.
- Extending `koblitz_tiny_ic` past 31 is now the prerequisite for any further
  scaling claim, and is a larger piece of work than a tournament round.

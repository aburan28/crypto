# Amendment 1 to PROTOCOL.md: sample sizes that can reject

**Status: frozen 2026-10-11, before any run at protocol N, v1 or v2.**
`PROTOCOL.md` is not edited. Both versions are run and reported:

- **v1** is `PROTOCOL.md` as frozen.
- **v2** is `PROTOCOL.md` with only the changes below.

## Why

A dry run at small N (`--protocol quick`, scratch output, not committed)
showed that v1's decision rule cannot reject in any comparative family.
With the add-one empirical tail the smallest attainable two-sided p is
`2/(N + 1)`:

| family | tests m | α/m | v1 N | v1 min p | v2 N | v2 min p |
|:--|--:|--:|--:|--:|--:|--:|
| E1 | 18 | 0.00278 | 400 | 0.00499 (**cannot reject**) | 2000 | 0.00100 |
| E3 | 6 | 0.00833 | 200 | 0.00995 (**cannot reject**) | 1000 | 0.00200 |
| E4 | 6 | 0.00833 | 100 | 0.01980 (**cannot reject**) | 500 | 0.00399 |

v1 can therefore only ever report "within the control distribution", and
its nulls are vacuous. That is a design error in the protocol, not a
result. v2's N are set by the rule now in the skill, `N + 1 ≥ 4m/α`:
1440, 480 and 480, rounded up. No result entered the choice.

## What changes in v2

- E1: `N = 2000` controls per curve (was 400).
- E3: `N = 1000` random points per curve (was 200).
- E4: `N = 500` random pairs per curve (was 100).
- E2 is not re-run. It has no rank test, and v1's E2 stands as v2's.

Nothing else changes: instruments, budgets, metrics, positive controls,
seeds (`20261011 + offset`), α, the Bonferroni families and the verdict
wording are all as in `PROTOCOL.md`. Because every control sample is drawn
first from its seeded stream, v1's control samples are exact prefixes of
v2's for E1, E3 and E4. The positive-control draws come after the controls
and differ between versions.

## Disclosure: what the dry run showed before this amendment

The dry run used N = 40 / 20 / 10 for E1 / E3 / E4, 5 planted values per
curve for the E1 and E3 positive controls, and 3 curves per kind for E2.
Its output was looked at. In particular:

- **E1.** Target M1 values: P-256 `n − 1` 125.0, against a dry-run control
  maximum of 91.1 (N = 40). All other targets were inside the dry-run
  control range. Planted positive controls scored M1 128–137.
- **E2.** Every planted point was found (27/27). No random-`b` curve gave a
  find (0/9).
- **E3.** Small-`x` positive controls had `ĥ/ln p` 0.65–0.67 and
  `h_x − ĥ` from −65 to −97 nats, against a generator `ĥ/ln p` of
  0.99–1.00.
- **E4.** P-192 and P-224 dry-run rows were produced. Their summary
  statistics were not examined before this amendment.

None of this changes v2's N, which the power rule fixes. It is recorded so
that a reader can discount v2 for having been designed after seeing it.

## Implementation notes (clarifications, not changes)

- E2 planted rational points draw `v` uniformly in `[0, √(p·d⁶))`. The
  protocol's "near `√(p d⁶)`" means the same bit size range, so `b` is at
  most `p` in absolute value, as for the NIST `b`.
- Instrument calls run in parallel across threads. Samples are drawn
  sequentially from the seeded stream first, and results are collected in
  sample order, so the output does not depend on the thread count.

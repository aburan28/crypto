# XP-20261015-p192-isogeny-aux-search-t48: an auxiliary curve with a point of order 2^48 over the P-192 prime

**Status: abandoned (superseded, 2026-10-16 claim audit).** The card is kept
as the record; do not run it. Three independent reasons, each sufficient:

1. **`t = 48` is not the attack-optimal parameter.** The PKM note's
   linear-algebra floor
   ([`RESEARCH_PKM_TOWER_ORACLE.md`](../../research/notes/index-calculus/RESEARCH_PKM_TOWER_ORACLE.md)
   §3.7) puts the optimal base at `|V| ≈ (m!·r·T)^{1/(m+1)}`, about `2^39`
   at `m = 4` for a 192-bit `r` (derived). `t = 48` is `192 / m`, the size at
   which every point decomposes, not the cost optimum, so the card priced a
   base about `2^9` larger than the one an attack would use. At `t ≈ 39` the
   ordinary `z2xz8` search is about `2^36` draws, roughly 190 core-hours at
   the measured 10.06 µs (derived), and the 45.6-minute run below already
   found `b = 38`.
2. **The argv contradicted the estimate.** `run.argv` passed
   `--family z2xz2` while `est_core_hours` (98,352) is the `z2xz8` figure;
   `z2xz2` at `t = 48` is 362,005 core-hours. The argv now says `z2xz8`, for
   the record only.
3. **The search is unnecessary.** The supersingular `y² = x³ + x` gives a
   point of order `2^t`, every `t ≤ 64`, in one draw, with F4 degree and
   relation yield indistinguishable from the ordinary tower (Lane H, below);
   the Dickson torus tower (`2^t | p + 1`) also exists on P-192 for `t ≤ 64`
   ([`research/pkm_dickson_torus_20261015/`](../../research/pkm_dickson_torus_20261015/README.md), A5c).

The gate the card named, lead 1, has also moved: at `m = 4, N = 20` the F4
solving degree is **`D ≥ 8`, NEGATIVE** for the `m = 4` plateau
([`research/ic_next10_20261014/03-pkm-degree/`](../../research/ic_next10_20261014/03-pkm-degree/README.md)).
Any further P-192 tower work is gated on the `m = 3` cell of
[XP-20261015-pkm-m3-n21-m4-n20](../XP-20261015-pkm-m3-n21-m4-n20/README.md).
The text below is the card as written on 2026-10-15.

*Original status line:* blocked on the wide-field isogeny-tower builder, the
192-bit field backend and PR #1889; recommended to run after lead 1 came back
bounded. Funding 131,072, the lowest of the four seed cards.

## Question

PKM's isogeny-tower factor base for a P-192-shaped prime (`v₂(p − 1) = 1`, so
there is no Kummer tower; encyclopedia §0.2) needs an auxiliary curve `E'/F_p`
with a rational point of order `2^t`, `t ≈ 48` at `m = 4`. How many curves must
be drawn before one has such a point, and what does each draw cost at 192 bits?

## Why it matters, and what it does not decide

This run prices *construction*. A curve found here is useful only if the
tower's F4 solving degree stays bounded as the tower grows (encyclopedia §8
lead 1). That question is unresolved, and the evidence so far points the
wrong way: `D` rises at `m = 2, N = 20`. So:

- If lead 1 comes back bounded, the curve found here is the input the relation
  phase needs, and the measured draws-per-success fixes construction's share of
  the cost.
- If lead 1 stays unbounded, this run only measures a feasibility constant.
- Either way it says nothing about P-192's security. Not finding a curve within
  the seed space is a bounded negative result about this family and this filter.

## Prerequisites

| needed | where | state |
|:--|:--|:--|
| solver binary `pkm-aux-search`: deterministic seed → curve draw, 2-Sylow halving filter, checkpoint/resume | `claude/isogeny-tower-wide` | in progress, not pushed |
| 192-bit `F_p` Legendre symbol and square root, measured | `claude/fp-backend` (Stage 0) | in progress, not pushed |
| the measured small-field search (`C = 1.05`) and its projection | PR #1889 | open |
| the checker `backlog_check_aux2t` (`src/bin/backlog_check_aux2t.rs`) | this repository | built |
| `PROTOCOL.md`: seed derivation, family, filter, success check, stop rule | this directory | not written |

The solver's flags in `card.json` are provisional; the builder lane fixes them.

## Measured cost (2026-10-15)

Study: [`research/pkm_isogeny_tower_p192_20261015/`](../../research/pkm_isogeny_tower_p192_20261015/README.md)
(pre-registered protocol, branch `claude/isogeny-tower-p192`). The solver
`pkm-aux-search` is built there (`src/bin/pkm_aux_search.rs`; library
`src/cryptanalysis/two_sylow_search.rs`).

- **Per draw (measured).** 10.06 µs on the `z2xz8` family (`X₁(2, 8)`,
  `b ≥ 3` forced), 13.89 µs on `z2xz2`. Median of 7 × 65,536 draws, one core
  of a shared 4-core host. The baseline implementation took 30.83 µs.
- **Success probability (measured, consistent with the model).** `2^28`
  `z2xz8` draws and `2^25` `z2xz2` draws give `P(b ≥ t)` within the
  pre-registered 99.9 % bands of `2^(3−t)` and `3·2^(−t)` at every tested `t`.
  (`b` is the exponent of the cyclic factor of the 2-Sylow `Z/2 × Z/2^b`.)
  Because `v₂(p − 1) = 1`, `Z/2 × Z/8` forces
  `k = 3` bits toward a point of order `2^t`, not 4.
- **`t = 48` (derived).** About `2^45` `z2xz8` draws × 10.06 µs gives
  **98,352 core-hours ≈ 11.2 core-years**, or 95,617–99,624 core-hours over
  the timing spread. Each `2^32`-draw shard takes about 12 core-hours. With
  `z2xz2`, `2^48/3` draws take 362,005 core-hours. The table below (the
  placeholder 2–20 µs) is superseded.
- **Found.** In 45.6 minutes of one core: an ordinary curve with `b = 38`
  (seed 196253408). The checker accepts it at `t = 38`.
- **Shortcut (measured).** Over this prime, `y² = x³ + x` is supersingular,
  with `#E = p + 1` and a cyclic 2-Sylow of order `2^64`.
  `pkm-aux-search --family supersingular-j1728 --t 48` returns a point of
  order `2^48` in one draw (7 ms), and `backlog_check_aux2t --expect-t 48`
  accepts it (also `t = 64`). That answers the card's construction question
  without the search. Whether a supersingular `E'` is acceptable to PKM's
  relation phase is not tested.

## Expected cost (derived before the measurement; superseded above)

**Draws.** PR #1889 measured the number of curves tried before one had a point
of order `2^t` at `C·2^t`, pooled `C = 1.05` over 224 builds at `t = 8…14`
(16 seeds, two curves; `research/pkm_isogeny_tower_ecbench_20261015/README.md`,
panel 3, on `claude/leads-isogeny-tower`). A family whose every member already
has `k` bits of 2-power torsion forces part of the condition, so the
*hypothesis* (to be measured at small `t` before the run) is

```
expected draws ≈ C · 2^(t − k)
```

with `k = 2` for full rational 2-torsion (`Z/2 × Z/2`, the card's default
`--family z2xz2`) and up to `k = 4` for a `Z/2 × Z/8` family from `X₁(2, 8)`.
At `t = 48`: `2^46 ≈ 7.0·10¹³` draws for `k = 2`, `2^44 ≈ 1.8·10¹³` for
`k = 4` (times `C ≈ 1.05`).

**Per draw.** With a halving filter (a point of order `2^j` is halved by
taking square roots; a draw that cannot be halved past a few levels is
rejected), most draws are rejected after `O(1)` Legendre symbols and square
roots at 192 bits. Stage 0 has not measured that cost yet, so the card carries
a placeholder range of **2–20 µs per draw**:

```
core-hours = C · 2^(48 − k) · (µs per draw) / 3.6·10⁹
```

| family | draws | 2 µs/draw | 20 µs/draw |
|:--|--:|--:|--:|
| `k = 2` (`Z/2 × Z/2`) | `2^46` | 39,000 core-h ≈ 4.5 core-years | 391,000 core-h ≈ 45 core-years |
| `k = 4` (`Z/2 × Z/8`) | `2^44` | 9,800 core-h ≈ 1.1 core-years | 98,000 core-h ≈ 11 core-years |

(`C = 1` in the table; multiply by 1.05. A core-year is 8,766 h.)

**The filter is not optional.** The search as measured in PR #1889 charges a
point draw, `[2^t]P` and a baby-step giant-step over `4√p/2^t` indices per
curve: about `2^26` group operations per curve at `b = 192, t = 48`, and
`≈ 2^74` in all (#1889's projection). That is `2^26`–`2^28` times the table
above; the card is only worth running with the halving filter.

**Parallelism.** Draws are independent, so wall time is core-hours divided by
cores. The card splits the seed space `[0, 2^48)` into 65,536 shards of `2^32`
draws (2.4–23.9 core-hours each). `2^48` seeds are 4 times the `k = 2`
expectation, so the probability that every shard fails is about `e^{−4} ≈ 2 %`
if the hypothesis holds. `stop = first_success`: the first verified curve ends
the search.

## What "done" means

- **Success:** a curve `(a, b)` and point `Q` with `[2^48]Q = O`,
  `[2^47]Q ≠ O`, checked by the independent checker (`outputs.verify`), plus
  the total draws, CPU time and per-draw cost actually measured, landed in
  `research/pkm_aux_search_p192_t48_<date>/` by a PR with its report.
- **Bounded failure:** every shard ran with no curve; reported with the exact
  seed space, family and filter. It does not show such curves are rare
  elsewhere.
- **Execution (today):** cairn lab tasks `<ID>/shard-<k>`, one seed range of
  `2^32` draws each, under the first-success stop rule of
  [backlog/README.md](../README.md): before each claim a runner checks
  `backlog cairn-found <ID> --lab DIR` (exit 3 means a curve is already
  recorded: stop); a shard that finds a curve writes
  `{out}/found/artifact.json`, published by `cairn lab exec` under
  `runs/<ID>/shard-<k>/found/`, and the runner, still holding the lease,
  records it with `backlog cairn-found … --lab-path … --checker-bin …`, which
  runs `backlog_check_aux2t --expect-t 48` from the pinned build on the
  lab's copy and sets `found` only on accept. One shard takes 2.4–23.9 h of one core at 2–20 µs per draw,
  so `ttl_s` is 6 h with a re-claim every 2 h; the checkpoint makes a shard
  resumable and so preemptible.
- **Upgrade path:** once cairn can pin a native checker
  ([problem record](../../problems/cairn-certificate-checker-python-only.md)),
  post the search as an objective (`backlog cairn-post`), so any solver can
  contribute and the pinned verifier decides each claim. As piecework it
  would pay `floor(reward / units)` per distinct accepted curve and never stop
  at the first; to pay once for the first curve, post it without piecework,
  so it settles on its first accepted claim. The validator's
  `funding ≥ piecework.units` rule applies to the piecework form.

## Is the search needed? A supersingular `E'` tested (2026-10-15, Lane H)

Study: [`research/pkm_supersingular_aux_20261015/`](../../research/pkm_supersingular_aux_20261015/README.md)
(pre-registered protocol, branch `claude/pkm-supersingular-aux`). It answers
the blocker this card records: may the PKM relation phase use the
supersingular `y² = x³ + x`, which over the P-192 prime gives a point of order
`2^48` in one draw?

- **Structure (exact, small primes).** The supersingular tower set
  `x(R' + ⟨Q⟩)` on `y² = x³ + x` is *not* the Dickson torus set under any
  Möbius change of `x` (no match in 6,507 sets at `t = 3…5` over 103 primes;
  one isolated match for the twist `y² = x³ − x`). It is a Vélu tower like the
  searched curve's, so it does not inherit the torus family's behaviour.
- **F4 solving degree (measured, `p < 2^31`, `N ≤ 20`).** The supersingular
  tower gives exactly the ordinary isogeny tower's `D` (and width and step
  count to within one column) in every cell, for `y² = x³ ± x`, `y² = x³ + 1`
  and its twist; both differ from the Dickson torus (5 against 4 at
  `m = 2, N = 4…10`).
- **Relation phase (measured, ecbench, 20–24-bit P-192-shape curves).** With
  `prime-isogeny:t=T,aux=supersingular` the auxiliary cost falls from 10–27 %
  of `S` (the charged search) to 0.9–1.9 % (the build alone), the relation
  yield is unchanged, and `S` is within 2.2 % of a plain base of the same size
  with unpriced work uncharged; under the AGENTS.md §5a pricing band (0, 0.1, 1
  GAE per unpriced operation) that comparison is unresolved on two curves and
  1–7 % dearer on the third, while the supersingular arm stays cheaper than the
  searched ordinary arm on the two larger curves at every weight.

**What it means for this card (derived).** Nothing measured distinguishes a
supersingular `E'` from an ordinary one except the cost of finding it, so the
~11-core-year ordinary search is not needed to supply the relation phase's
input: `pkm-aux-search --family supersingular-j1728 --t 48` already does, in
one draw. The card stays `blocked` (its status is not changed here): running
it would now only measure a construction constant for ordinary curves, and
the decision whether a reason outside these experiments (none is known) should
require an ordinary curve rests with the coordinator. Lead 1 (whether `D`
stays bounded as `N` grows) is untouched and remains the gate for any use of
the tower at P-192.

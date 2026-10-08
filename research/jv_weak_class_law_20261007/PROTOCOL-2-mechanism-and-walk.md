# Protocol 2: the mechanism behind the depth-1 exclusion, and the refuse-ascend-walk policy against a random walk

Frozen 2026-10-08, before the instrument's new modes are built or run.
**Stage diagnostic, toy sizes.**  `S`, end-to-end cost and speedup are
**unset**.  Nothing here concerns a prime-field or deployed curve.  Builds
on [`PROTOCOL.md`](PROTOCOL.md) and [`README.md`](README.md) (PR #1556).

## Part A — mechanism (toward a proof of A1)

For a full-2-torsion curve `y² = (x − e₁)(x − e₂)(x − e₃)`, the point
`(eᵢ, 0)` is halvable (lies in `2E(F)`) iff `eᵢ − eⱼ` and `eᵢ − eₖ` are
both squares.  The **pattern** of a curve is its number of halvable
2-torsion points, in `{0, 1, 3}`.  For a weak curve the protocol's
derivation gives: `0` iff `α` is a non-square, `1` iff `α` is a square
and `α − σα` is not, `3` iff both are squares.  The 2-isogeny with
kernel `⟨(eᵢ, 0)⟩` has a **direction**: up, down or level, read from the
census heights.

Mode `patterns` tabulates, per size, (i) nodes and weak nodes by
`(height, pattern)`, (ii) for nodes with exactly one halvable point,
whether the isogeny with that kernel goes up or down, and (iii) for
nodes with pattern 3, the directions of the three isogenies.

- **M1 (ascent rule).**  For nodes with exactly one halvable point and
  height below the crater, one of the two rules "the halvable point's
  isogeny ascends" / "the halvable point's isogeny descends" holds on at
  least 99% of nodes at every size.  The rule is then the intrinsic
  ascent test the walk uses; it is read from the data, not assumed.
- **M2 (pattern at height 1).**  At height 1 in a class of depth ≥ 2,
  weak curves have pattern 1 or 3 only; pattern 0 never occurs for a
  weak curve at height 1.  *Falsified by one weak node of pattern 0 at
  height 1.*  If M2 holds, A1 reduces to: a depth-1 crater curve of
  pattern 1 or 3 is never weak, which is the statement to prove.
- **M3 (pattern 3 means height ≥ 2).**  Every node of pattern 3 has
  height ≥ 2.  This is the textbook fact `E[4] ⊂ E(F) ⟺ height ≥ 2`
  and is a correctness check on the instrument.

## Part B — the walk policy, measured on ground truth

Both policies use the same move set on the same exhaustive node tables:
rational 2-isogenies whose codomain keeps full 2-torsion, and rational
3-isogenies (kernel a root of the 3-division polynomial in `F_{p⁶}`,
codomain in Legendre form from the images of the 2-torsion points).
Degrees 5 and 7, which the ledger's walk also used, are **not** in this
instrument; both arms lack them equally, and the comparison is between
policies, not against the ledger's absolute numbers.

- **Policy R (random, §17-like).**  Start at a random full-2-torsion
  curve of the class; at each step choose uniformly among the available
  moves; stop at a weak curve, after `3q` steps, or when `50 × (distinct
  j seen)` steps pass without a new `j`.
- **Policy N (refuse, ascend, walk level).**  Compute `v₂(f)` from the
  public trace; if it is 1, refuse at zero steps.  Otherwise ascend by
  the M1 rule until height `min(3, v₂(f))`, then move only by
  3-isogenies and by 2-isogenies whose codomain has the same height;
  same stop rules.

Twenty starts per class, every class with a full-2-torsion curve, sizes
`p ∈ {7, 11, 13}`, seed `20261008`.  A refusal of a class that holds a
weak curve counts as a failure of N.

- **W-1 (refusal exactness).**  N refuses no class that holds a weak
  curve.
- **W-2 (success).**  On classes that hold a weak curve, N succeeds in
  at least 90% of starts at every size; R in fewer than 75%.
- **W-3 (steps).**  On classes where both succeed, the median curves
  visited by N is at most `0.4×` R's at `p = 11, 13`.
- **W-4 (ascent cost).**  N's ascent takes at most `v₂(f)` steps on
  every start, as the M1 rule predicts.

## Decision rule

M1–M3 passing fixes the intrinsic ascent test and reduces A1 to one
statement about depth-1 crater curves.  W-1 to W-4 passing makes N the
registered replacement for the cover walk's §17 policy, to be wired into
`jv_isogeny_walk.rs` on branch `research/jv-cover-end-to-end-20261007`
as the next PR, with degrees 5 and 7 restored there.  Any failure is
reported with its table.  Class: **engineering** for the walk, **boundary**
for the mechanism statements.

## Inadmissible

Reading the ascent rule off the heights at walk time (the walk must use
the halvability test); giving N moves R does not have; counting a refusal
as a success; comparing against the ledger's success rates, which used
more degrees.


## Status, 2026-10-08

Run.  M1, M3 pass; M2 falsified; W-1, W-4 pass; W-2 fails for N (three
disclosed versions) and its R clause was mis-set; W-3 met by N where it
succeeds, not by A at `p = 11`.  Policy A (refuse, ascend, then walk
freely) was added **after** N failed W-2 and is labelled post hoc in the
results.  See `README-2-mechanism-and-walk.md`.

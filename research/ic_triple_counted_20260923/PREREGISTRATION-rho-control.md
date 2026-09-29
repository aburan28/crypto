# Pre-registration: rho with the IC arm's own canonical form

Written and committed **before the control exists as code and before it has
run on anything**, per `AGENTS.md` §4.  Written *after* the counted-sizing
confirmation (`PREREGISTRATION.md`, item 3: counted below the worker's rho at
`n37a0`, `n43a1`, `n59a0` and `n61a1`), and prompted by what profiling found
while that run was going.  Nothing below is a result.

## What the profile found

Every IC/rho ratio in this lane is taken against the worker's matched rho. That
means `koblitz_signed_frobenius_rho` with automorphism order `2n`, the rho of
every tournament round since round 0006.

**Its canonical form walks the Frobenius orbit by squaring.**
`FastRhoWalk::canonicalize` (`koblitz_index_calculus.rs`) runs `x` through all
`n` Frobenius images with `sqr_k`, then squares `y` up to the winning power.
That is about `1.5n` field squarings every step.

The IC arm names the same orbits with a normal-basis rotation.
`NormalBasis::canon` in `koblitz_tiny_ic.rs` is nibble tables plus a zero-run
search, and it arrived in round 0007 (`RESULTS.md`, "nibble-table normal
basis"). The rho baseline never got an equivalent.

**What it costs**, from `rho-lines-n43a1.txt`: a line-table build, and
callgrind's self cost by source line over one exposed probe fixture's
`rho_solve`.

- 54.7% is field squaring and reduction (`semaev_decomp.rs`).
- 15.0% is the bounds-checked reduction-table index.
- 21% is the orbit loop's compare and range.
- The point addition is under 1%.
- Per step, that is 9,930 instructions, of which about 9,100 are the canonical
  form.
- An IC probe costs about 1,100, and it is also one addition, one
  canonicalisation and one lookup.

`rho-phases-probe.json` has the phase split at every cell, on the same exposed
stream.

**So the registered comparison was against a rho carrying a canonicalisation
the IC arm had engineered away.** This control removes that one difference.

## The control

`rho-normal-basis.patch`, on top of `triple-table.patch`, changes only
`FastRhoWalk::canonicalize`:

- **Name the orbit with the IC arm's `NormalBasis::canon`.** This is the same
  code, made `pub(crate)` and constructed once from the rho field's own
  squaring. It gives the least rotation of `x`'s normal coordinates and the
  power `k`.
- **Recover `φ^k(x)` and `φ^k(y)` with the inverse basis change.** That is
  nibble tables over the basis columns, plus a rotation of `y`'s normal
  coordinates. Only this is new code.
- **Choose the sign exactly as now** (the lesser of `y` and `x ⊕ y`), and scale
  the coefficients by `±λ^k` exactly as now.
- **Fall back to the squaring loop** when `x`'s normal coordinates are
  rotation-invariant (all zeros or all ones). No least rotation names the
  point there. This happens at two points of the curve at most.

The representative changes: least normal-basis rotation instead of least
polynomial-basis abscissa. It is still a function of the class
`{±φ^j(P)}`, so the walk remains a deterministic map on classes, and
collisions, distinguished points and fruitless-cycle escapes keep their
meaning. Walk lengths differ fixture by fixture, and their expectation does
not.

## Prediction

The canonical form costs about `212·n` instructions a step (fitted at
`n43a1`), and the replacement about 400. Applied to each cell's measured mean
rho cost from the confirmation run, with its fixed phases from
`rho-phases-probe.json`:

| cell | counted/rho (measured) | control/rho (predicted) | **counted/control** | band |
|---|--:|--:|--:|--:|
| `n23a1` | 1.210 | 0.66 | **≈ 1.85** | 1.2–3 |
| `n37a0` | 0.795 | 0.25 | **≈ 3.1** | 2–5 |
| `n43a1` | 0.824 | 0.16 | **≈ 5.2** | 3–8 |
| `n59a0` | 0.823 | 0.66 | **≈ 1.25** | 1.05–1.6 |
| `n61a1` | 0.714 | 0.18 | **≈ 3.9** | 2.5–6 |

- **`n59a0` is dampened, not different.** Its fixture set-up costs 259M
  instructions, and IC and rho both pay it. Take that out and the cell looks
  like `n61a1`.
- **The bands are wide on purpose.** The constant per step is fitted at one
  cell and one fixture, and walk lengths vary from fixture to fixture.

## Measurement

**Fixtures:** the confirmation run's own fixtures (`confirm.json`): 32 a cell,
and 128 at `n43a1`, with the same seeds. The control runs once on each, in rho
mode, counted by callgrind. `counted` and the current rho are the committed
counts from that run, so every ratio here is paired on the fixture.

**Statistic:** as registered for the confirmation. Per cell, the geometric
mean of paired ratios with a 95% percentile bootstrap interval (10,000
resamples, seed 0).

## Decision rules, fixed now

**Gate, absolute.** Every control report verifies in `oracle.py`'s rho mode
and recovers the same logarithm as the other three arms on that fixture. One
failure is a defect in the control, and no ratio is reported. An unfinished
run removes its cell, as before.

**1. The control is a rho, not a different walk.** Its mean walk additions
over the current rho's lie in `[0.85, 1.15]` at every cell. This is reported
as a check. Outside it, the control's walk is suspect, and item 2 is reported
with that caveat rather than scored.

**2. The e2e claim under a matched canonical form.**

- **Supported (the IC arm loses to a matched rho):** `counted`/control is
  *above one* at all five cells, meaning the interval's lower end is above one.
- **Refuted:** `counted`/control is *below one* at any cell, meaning the
  interval's upper end is below one.
- **Otherwise:** partial, cell by cell.

**What each outcome means for the lane.**

- **If supported,** every IC-below-rho result that compared against this rho
  is a comparison with a handicapped baseline. That covers this study's item
  3, the triple study's secondary table, and the tournament's
  `beats_rho_strict` rounds 0007–0020. The IC arm has not beaten rho; it has
  out-engineered one step of it. That would be an **accounting** correction
  in `AGENTS.md` §3's terms. It says nothing against index calculus in
  general, and nothing about ECC2K-130, where the verdict already stood
  against IC.
- **If refuted,** the canonical form was not the gap, and the note will say
  what was.

## Scope, fixed now

- **Instructions only, on one machine.**
- **The control changes one function**, and it is not offered as the best
  possible rho. Other asymmetries may remain, in either direction.
- **This is not a tournament round.** Whether the tournament's matched rho
  should carry this change is the tournament lane's call, and a pointer to
  this note goes to that lane's open PR only with the user's agreement.

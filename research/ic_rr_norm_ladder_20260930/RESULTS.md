# The Riemann–Roch norm form at `m = 3`: results

This was written after the run. [PREREGISTRATION.md](PREREGISTRATION.md) is unchanged since
its registration commit `9f0dc489`, apart from its two dated, additive amendments (resumable
cells; order of measurement). The readout is [runs/registered/readout.txt](runs/registered/readout.txt).

## Run

| | |
|:--|:--|
| instrument | `examples/rr_degree_ladder.rs`; the `rr`, `x4` and `ctrl` builders are byte-identical to `9f0dc489` (the later commits added `--resume` and the two-pass order) |
| binaries | `b77b246a…` (first resume, amendment 1), `f6c5d1d9…` (second resume, amendment 2); recorded in `runs/registered/resumes.txt` |
| wall | 2026-10-01: 00:47Z (cut within a minute), 01:52–03:52Z, 15:10–17:10Z; three lanes on CPUs 1–3 under the benchmark lock |
| cells | 14; all nine `ℓ ≤ 4` cells complete (exit 0, four rootless draws each); the five `ℓ ≥ 5` cells cut at their CPU limit with attempts left (see below); no cell hit the matrix caps |

## Registered verdict: constant lever

**The norm form does not flatten the degree in `ℓ`.** `s̄_rr = 1.167`, far above the 0.35
threshold, with all three curves fitted:

| `ℓ` | `rr` unknowns | `K₀/2¹³` | `K₁/2¹⁷` | `K₁/2¹⁹` |
|--:|--:|:--|:--|:--|
| 2 | 7 | 4 5 4 4 | 4 4 4 4 | 4 4 4 4 |
| 3 | 10 | 6 5 6 6 | 5 5 5 5 | 5 5 5 5 |
| 4 | 13 | 7 7 7 7 | 6 6 6 6 | 6 6 6 6 |
| 5 | 16 | ≥8 ×4 | ≥8 4 ≥8 ≥8 | 4 7 7 7 |
| 6 | 19 | — | ≥7 ×3 (scan to 6) | ≥7 ×2 |
| **`s_n`** | | **1.500** (ℓ 2–4) | **1.000** (ℓ 2–4) | **1.000** (ℓ 2–5) |

Resolved cells give medians 4, 6, 7 / 4, 5, 6 / 4, 5, 6, 7: about `ℓ + 2`, the same law the
chained `S₃` ladder (5, 6, 6, ≥7 at `ℓ = 2…5`) and the torsion-symmetrised ladder (`s̄ =
1.033`) read. The full-scan lower bounds at `ℓ = 5` (`K₀/2¹³`, `K₁/2¹⁷`) and `ℓ = 6`
(`K₁/2¹⁷`) are consistent with it and were not needed for the verdict.

## Predictions

1. **Constant lever: holds.**
2. **`rr` below `x4` on the same draw, by 1 or 2: holds on the direction, wrong on the
   size.** On all 36 draws resolved on both arms, `rr` is below `x4`: mean −2.44, range −3
   to −1, never equal or above. The gap is mostly 2–3, not 1–2. The direct `S₄`
   presentation reads 6, 8, 9 at `ℓ = 2, 3, 4` on every curve (`s̄_x4 = 1.5`); the norm form
   reads 4–6 there. Its degree-4 generating set of the same ideal is caught 2–3 degrees
   earlier by the Macaulay scan, and that advantage is a constant: the two slopes are 1.0–1.5
   and 1.5.
3. **The control is not below `rr`: holds.** A random system of `rr`'s shape resolves at the
   same degree or higher in every cell (medians 5, 6, 7 against 4, 6, 7 on `K₀/2¹³`; 5, 6, ≥8
   against 4, 5, 6 on `K₁/2¹⁷`; ≥8, ≥8, 7 against 4, 5, 6 on `K₁/2¹⁹`), often by pinning
   rather than refutation. `rr`'s structure makes it a little easier than random, by a
   constant.

## What the heavy cells show, and what they cost

Reported, not decisive. Per-draw wall seconds (medians), from the logged `secs`:

| cell | `rr` (`d_max`) | `x4` (`d_max` 9) | `ctrl` |
|:--|--:|--:|--:|
| `K₀/2¹³ ℓ = 4` | 50 (7) | 25 | 127 |
| `K₁/2¹⁷ ℓ = 4` | 4 (7) | 39 | 163 |
| `K₁/2¹⁹ ℓ = 4` | 5 (7) | 45 | 182 |
| `K₀/2¹³ ℓ = 5` | 579 (7) | 147 | 2,433 (one) |
| `K₁/2¹⁷ ℓ = 5` | 1,020 (7) | 279 | — |
| `K₁/2¹⁹ ℓ = 5` | 1,264 (7) | 383 | 5,486 (one) |
| `K₁/2¹⁷ ℓ = 6` | 106 (6) | 1,584 | 238 |
| `K₁/2¹⁹ ℓ = 6` | 131 (6) | 2,300 | 330 |

- At `ℓ = 5` `rr` resolves at 7 on `K₁/2¹⁹` (three of four draws, one at 4) and does not
  resolve by 7 on the other two curves: `≥ 8`. `x4` reads `≥ 10` on every `ℓ = 5` draw.
- At `ℓ = 6`, scanned to 6, `rr` reads `≥ 7` on every draw, in about two minutes each; `x4`
  to degree 9 is the expensive arm there, 25–40 minutes a draw, and reads `≥ 10` or `≥ 9`.
- The random control is the costliest arm at `ℓ = 5`: 40–90 minutes a draw. Its `ℓ = 5`
  medians could not be formed.
- Three `rr` draws refute at the system degree 4 at `ℓ = 5` (one on `K₁/2¹⁷`, one on
  `K₁/2¹⁹`) and one at `ℓ = 2` reads 5; the rest of a cell's draws agree with each other.
  A degree-4 refutation means the 2n + 1 equations alone are inconsistent at degree 4,
  the Kosters–Yeo regime the RR panel §9 describes for `3ℓ ≪ n`.

**Attempts.** Amendment 1 allows three attempts of 6,000 CPU-s per heavy cell. After the
second window every curve had at least three resolved `ℓ`, so the fit, and the verdict, were
complete; the remaining attempts (two for four cells, one for `K₁/2¹⁹ ℓ = 5`) were **not
run**, by the author's decision, because they can only add lower bounds at `ℓ = 5, 6` and the
session environment does not survive unattended multi-hour runs (three cuts in this run).
The files resume with `run.sh … --resume` unchanged if anyone wants them; a cut cell's
missing draws are censored, never negative evidence.

## What this says

- **The Riemann–Roch norm form is a cheaper presentation, not a slope lever.** At `m = 3` it
  needs `3ℓ + 1` unknowns against the chain's `3ℓ + n`, refutes 2–3 degrees below the
  direct `S₄` and at about the chain's own degree, and its degree grows with `ℓ` at the
  chain's rate. That is a constant, and the registration's "constant lever" reading
  applies: the `m = 4` norm form (`4ℓ + n + 1` unknowns, cubic) is **not** built for the
  exponent audit. It stays available as an engineering arm: a cheaper way to run the
  audit, not a way to change its slope.
- **With the search form (survey §3.3, RR panel §8: `|F|²`, Shoup) and the support form
  ([support-degree.txt](support-degree.txt): system degree rising 2 per unit `ℓ`), every
  algebraic reading of the Nagao / Riemann–Roch encoding is now closed for the exponent at
  these sizes.** What remains of §3.3 is Nagao's first-fall-degree claim on the
  incidence/auxiliary-variable form, which is the SEMBIN lane's question and is not
  duplicated here.
- **The symmetric-group lever of survey §3.2 is measured.** The norm form is fully symmetric
  in the summands; its slope is 1.0–1.5. With the torsion lever (`s̄ = 1.033`) and the
  Frobenius item
  ([closed by structure](../ic_candidate_tournament_20260915/campaign_20260916/NOTE-20260930-frobenius-orbit-coordinates.md)),
  every symmetry in §3.2 has now been slope-tested or ruled out, and none flattens the
  degree.

## Scope

Two Koblitz curves, `n` 13–19, `ℓ` 2–6, `m = 3`, four rootless draws per cell (three at
`K₁/2¹⁷ ℓ = 6`, two at `K₁/2¹⁹ ℓ = 6`), one engine's sparse Macaulay scan to degree 7 (6 at
`ℓ = 6`; 9 for `x4`). A refutation degree is a stage diagnostic (AGENTS.md §2, §5): no
end-to-end cost, no yield, nothing at `n ≈ 83` or 131. The constant it finds is real at
these sizes and says nothing about any other presentation of the ideal.

## Erratum 1 (2026-10-03): the instrument's floor, and the "degree-4 refutations"

Found while calibrating an external engine on the same systems
([ic_gb_ladder_20261003](../ic_gb_ladder_20261003/PREREGISTRATION.md), §2). Additive;
nothing above is edited.

- **The scan has a floor at the system degree.** `measure()` in
  `examples/rr_degree_ladder.rs` starts at `system_degree(polys)`, so no reading below 4
  (`rr`) or 6 (`x4`) was possible. A system refuted below its own degree reads as refuted
  *at* it. Singular's truncated homogeneous basis, which has no such floor, reads the
  direct `S₄` at `ℓ = 2` on `K₁/2¹⁷` as **5** on three of four draws (6 on the fourth), not
  the 6 reported above. The `x4` slope quoted above (1.5 over `ℓ = 2…4`) is therefore an
  underestimate; it does not change the verdict, which rests on `rr`'s slope and on `rr`
  never falling below the control.
- **The two "degree-4 refutations at `ℓ = 5`" on `K₁` are not Kosters–Yeo refutations.**
  `K₁/2¹⁷ ℓ = 5` draw 1 and `K₁/2¹⁹ ℓ = 5` draw 0 contain the **constant equation `1`**:
  the `rr` form's trace condition `Tr(e₁) = 0` is constant when `V ⊂ ker Tr` and
  `Tr(x_R) = 1` (probability about `2^{−ℓ−1}` per draw), and such a draw is inconsistent at
  degree 1 for a reason that has nothing to do with the summation ideal. The floor
  reported it as 4. The sentence above that reads these as "the 2n + 1 equations alone
  are inconsistent at degree 4, the Kosters–Yeo regime" is withdrawn. (`K₁/2¹⁷ ℓ = 2`
  draw 1 is the same kind of draw and read 4 at the floor; the cell's other three draws
  read 4 genuinely.) The cell medians are unchanged: 4 was never a median. The
  external-engine ladder reports such draws as `triv` and excludes them from medians.
- The direct `S₄` descent does not contain that linear condition explicitly; on the same
  draws `x4` read `≥ 10`. That the norm form exposes the trace obstruction as one linear
  equation is a (constant-size) feature of the presentation, consistent with the verdict.

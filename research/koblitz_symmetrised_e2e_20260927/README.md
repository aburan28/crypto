# Torsion-symmetrised index calculus on Koblitz curves, end to end: frozen contract

Companion to
[`research/notes/index-calculus/RESEARCH_EXOTIC_COORDINATES.md`](../notes/index-calculus/RESEARCH_EXOTIC_COORDINATES.md)
(§8, §11, §17, §18: the symmetrised system, oracle-only measurements) and
[`research/notes/ecdlp-general/RESEARCH_TORSION_AUXILIARY_INPUTS.md`](../notes/ecdlp-general/RESEARCH_TORSION_AUXILIARY_INPUTS.md)
(the torsion thread this round closes).  This file is the contract:
what is measured, in what unit, against which boundaries, and what
would count as success.  It was written before the runs in `results/`
were made; `run.sh` is the exact command list.

## The question

Every Koblitz curve `K_a : y² + xy = x³ + a x² + 1` over `F_{2^n}` has
one rational 2-torsion point `T = (0, 1)`.  Translation by `T` acts on
abscissae as `x ↦ 1/x`, and in the frame `u = 1/(x + 1)` (which the
`F_2`-Frobenius preserves) as `u ↦ u + 1`.  A factor base
`F_u = {P : u(P) ∈ V}` for a Frobenius-stable `V ∋ 1` is closed under
`+T`, and the summation polynomial can be rewritten in the
translation invariants `w = u² + u` (Artin–Schreier) and `s = Σ u`,
which lowers its degrees from `[4,4,4,4]` to `[2,2,2,2,1]` for `S₄`.
Earlier rounds measured the *oracle alone* (one decomposition at a
time) and found the symmetrised system refutes non-decomposable targets
×300–6000 faster at `n = 15` but loses ×1.6 at `n = 17, m = 3`.

**Does the symmetrised system, run as a complete index calculus
(factor base, relation collection, linear algebra, verification) with
every phase priced, beat the plain `x`-frame decomposition on the same
curve, and does either beat Pollard rho, in the repository's unit?**

## Boundaries, stated first

- **Reference.**  The matched Pollard rho `ic bench` runs on the same
  instance and the same planted targets: on a Koblitz curve the lower
  mean `S` of the signed-Frobenius walk (`A = 2n`) and the negation walk
  (`A = 2`), each counted, distinguished points, verified.  Ledger
  values on these instances: `S ≈ 0.9–3.6` (`RESEARCH_IC_BOUNDARY_LEDGER.md`).
- **Generic floor.**  `Ω(√(r / A))` group operations for any generic
  algorithm on a group with an automorphism group of order `A`.  The
  2-torsion translation is not a group automorphism (it is a
  translation, `P ↦ P + T`), and `u ↦ u + 1` acts on the *factor base*,
  not on the DLP group; nothing in this directory touches the floor.
- **Counting floor for the factor base** (torsion note §4): an
  `m`-fold decomposition over a base of `2^ℓ` points folded by a
  symmetry of order `γ` needs `2^{ℓ m} / (m! · γ) ≳ #E / P_hit` targets
  to be representable; the symmetry lowers the column count and the
  degrees, never the trial count per relation.  So the *expected*
  classification is **engineering** (a constant factor on the solve),
  and any claim of more needs the ratio to grow with `n`.

## Unit

`S = total group-addition equivalents / √r`, taken over the whole
pipeline by `ic bench`: factor-base construction, the target walk,
every decomposition attempt (hit or miss) at the solver's own operation
count (word XORs for both algebraic arms, pinned at the repository's
`ns_per_word_xor`; lookups for the combinatorial arms), the relation
matrix, and the final verification.  Wall time rides along and is
not the metric.  The `vs rho` column is `S / S_rho` on the matched
reference.

## Instances

Koblitz instances with a prime-order subgroup in the roster (`K₀/2¹⁷`
and `K₁/2³¹` have none):

| curve | `r` | `log₂ r` | `V` (indices into the factors of `xⁿ − 1`, index 0 = `x + 1`) | `dim V` |
|:--|--:|--:|:--|--:|
| `icv1-f2m17-tm101-00378d4e` | 65587 | 16.0 | `0;1` | 9 |
| `icv1-f2m23-t5197-69e76b73` | 2095853 | 21.0 | `0;1` | 12 |
| `icv1-f2m23-tm5197-1f85e9e1` | 4196903 | 22.0 | `0;1` | 12 |
| `icv1-f2m31-tm90707-c95f16f5` | 1439393 | 20.5 | `0;1;2` | 11 |

`n = 13, 19, 29, 37` have no proper Frobenius-stable `V ∋ 1` (2 is
primitive mod `n`); composite `n` is degenerate (subfield points
collapse the columns); `n = 41` has only `dim V = 21` and is out of
range.  `n = 31, dim 16` and every `m = 3` system at `n = 23` (34
boolean variables) are out of range by §17's oracle-only measurements
and are not run.

## Arms

Same seed, same planted targets and same rho reference for every arm on
an instance.  `D` is the divisor above.

| arm | factor base | decomposition | solver | `m` |
|:--|:--|:--|:--|--:|
| `sym` | `koblitz-symmetrised:divisor=D` | `symmetrised` (`w, s` system) | inherited-F4, Macaulay cap 3, 4096 splits | 2 |
| `sym4` | same | same, Macaulay cap 4 | | 2 |
| `x` | `koblitz-orbit:divisor=D` | `descent-algebraic` (plain `x`, Weil descent) | inherited-F4 | 2 |
| `mitm-x` | `koblitz-orbit:divisor=D` | `mitm-frobenius` (combinatorial control) | — | 2 |
| `mitm-u` | `koblitz-symmetrised:divisor=D` | `mitm-frobenius` on the `u`-frame base | — | 2 |
| `sym-m3`, `x-m3` | as `sym`, `x` | | | 3 (`n = 17` only) |

Three repeats per arm (three planted logarithms), eight counted rho
runs per instance.  Every run has a 3-hour wall cap; a run that hits it
is kept as a failure in `results/` with its log.

## What is frozen

For every arm, `results/<instance>__<arm>.json` (the full `ic bench`
report: calibration, rho reference per run, every phase's cost, the
system shape, the recovered logarithm and its verification) and
`results/<instance>__<arm>.fb.json` (the factor base point by point:
curve, `V` basis, each point's `x`, `y`, `u = 1/(x + 1)`, its column and
its negative, and a BLAKE3 hash of the point list).  `results/host.txt`
is the machine and binary; `run.sh` the commands; `report.py` the
tabulation that produces `RESULTS.md`.

## Falsification target, declared in advance

The symmetrised system counts as an **engineering result** if, on
every instance, with the logarithm verified on every repeat:

- `S_sym / S_x < 0.8` (the symmetrised solve is cheaper than the plain
  one at the same `m`, same `V` dimension and same seeds), with the
  relation counts and the columns showing the same structure (the fold
  does not change what a relation is).

It counts as an **advance** only if `S_sym / S_x` falls with `n` across
the three degrees (a growing gain, which the counting floor above says
should not happen), or if any arm reaches `S / S_rho < 1` at fixed
accounting.

It is **refuted** (the symmetry is not worth its rewriting) if
`S_sym / S_x ≥ 1` on the instances where both finish, or if the
symmetrised arm fails to finish where the plain one does.

Inadmissible: changing `V`, `m`, the seed or the pricing between arms;
reporting the oracle's cost without the trials that missed; a rho
column from a different instance or unmatched seeds.

# Auxiliary inputs and torsion points against Pollard rho: frozen contract

Companion to [`research/notes/ecdlp-general/RESEARCH_TORSION_AUXILIARY_INPUTS.md`](../notes/ecdlp-general/RESEARCH_TORSION_AUXILIARY_INPUTS.md),
which carries the tables, the reading and the verdict.  This file is the
contract: what is measured, in what unit, against which boundaries, and
what would count as success.  It was written before the sweeps in
`results/` were run; the harness docstrings state the same boundaries.

## The question

Do *auxiliary inputs* (`[α^i]G` leaked by a protocol) or *torsion points*
(small-order points on the curve, or an auxiliary curve's torsion
subgroup used as the embedding group) let a collision search or an index
calculus beat Pollard rho on a prime-order elliptic-curve subgroup, and
by how much, in operations?

## Boundaries, stated first

- **Reference.**  Pollard rho with the negation map on the same curve,
  subgroup, target and accounting (`rho.py`): `S ≈ 1.3` in this
  repository's unit, historically measured 1.2–1.8.
- **Floor for the plain ECDLP.**  `Ω(√p)` group operations for any
  generic algorithm given only `G` and `[α]G` (Shoup 1997).  No auxiliary
  input can be manufactured from those two points without solving CDH,
  and a torsion point of the curve is just another group element, so this
  floor is untouched by anything in this directory.
- **Floor for the DLP with auxiliary inputs.**  `Ω(√(p/d))` group
  operations for a generic algorithm given `G, [α]G, …, [α^d]G`
  (Cheon 2010, Theorem on the generic bilinear model; Kim–Cheon–Song
  survey §3).  Cheon's `p ± 1` algorithms attain it up to the cost of an
  exponentiation, and only when `d` divides `p − 1` or `p + 1`.

## Unit

`S = total group operations / √p`, with every affine addition or
doubling counted as one operation (`ec.Counter`).  Every phase is
charged: fixed-base tables, the polynomial multi-scalar sums of the
`p + 1` case, the auxiliary-curve arithmetic of the torsion route, the
kangaroo restarts, and the final verification `[α]G = [α]G`.  The
auxiliary inputs themselves are *given* by the problem and are not
charged; the note says so wherever it matters.

Two ratio columns: `S / S_rho` (the multiple of the reference) and
`ops / √(p/d)` (the multiple of the DLPwAI floor).

## Rows

| variant | inputs | what it is |
|:--|:--|:--|
| `rho_reference` | `G, [α]G` | the boundary |
| `plain_bsgs_as_in_rust` | `G, [α]G` | what `cheon_attack.rs` computed before this round: BSGS with `d` baby steps, auxiliary input unused |
| `cheon_p-1_bsgs_naive` | `G, [α]G, [α^d]G`, `d \| p−1` | Cheon Theorem 1, every step a full double-and-add |
| `cheon_p-1_bsgs_comb` | same | same, fixed-base comb tables (Kozaki–Kutsuma–Matsuo) |
| `cheon_p-1_kangaroo_comb` | same | Cheon §3.1, memoryless, distinguished points |
| `cheon_p+1_bsgs_comb` | `G, …, [α^{2d}]G`, `d \| p+1` | Cheon Theorem 2 through the norm-one torus of `F_{p²}` |
| `torsion_embedding_honest` | `G, …, [α^{d²}]G`, `d \| #E'(F_p)` | Kim–Cheon elliptic embedding with the only computable comparison (pairwise) |
| `torsion_embedding_oracle_quotient` | same + `α` | the same route handed the exponent quotient for free: shows the `√δ` birthday count is real and the division is the wall |

Sizes: `p` of 24–48 bits for `p − 1` (d ≈ p^{1/2}), 24–40 bits for
`p + 1` (d ≈ p^{1/3}; a 44-bit run was started three times and lost to
container restarts each time, so it is not in the tables), 14–20 bits
for the torsion route (d ≈ p^{1/4}, the
`d²` auxiliary inputs and the degree-`d²` division polynomials cap it).
Three planted secrets per size; the reference runs on the same three.

## Falsification target, declared in advance

The thread counts as a **result for the DLPwAI** if, on every seed, with
the answer verified against the planted `α`:

- `ops / √(p/d) ≤ c · log₂ p` for a constant `c` that does not grow with
  `p` across at least four sizes (the algorithm sits on its floor up to
  the exponentiation cost), and
- `S / S_rho` falls with `p` at a fitted exponent of `−1/4 ± 0.05` for the
  `p − 1` case at `d ≈ √p` and `−1/6 ± 0.05` for the `p + 1` case at
  `d ≈ p^{1/3}`.

It counts as a **result for the plain ECDLP** only if some row that sees
nothing but `G` and `[α]G` has `S / S_rho < 0.9` at fixed accounting,
which the floor above says cannot happen; a row that appears to do so is
an accounting error to be found.

The torsion route counts as **closed** if the honest row costs
`Θ(p/d)` operations, i.e. its `S` grows like `p^{1/2 − ε}` rather than
falling, while the oracle row finds its collision in `O(√(p/d))` samples.

Inadmissible: changing `d` after seeing a result, dropping the table or
polynomial phases from the count, counting an unverified `α`, choosing
seeds, or reporting exponentiations as if they were group operations.

## Reproduce

```
python3 run.py --bits 24 28 32 36 40 44 48 --cases p-1 --seeds 3 --rho-seeds 3 --rho-cap-bits 48
python3 run.py --bits 24 28 32 36 40       --cases p+1 --seeds 3 --rho-seeds 3 --rho-cap-bits 40
python3 torsion_embedding.py --bits 14 16 18 20 --seeds 3
python3 report.py > RESULTS.md
```

Pure Python 3, no dependencies.  Curves are generated deterministically
from the size (`find_curve(bits, seed=1000·bits + case)`), so a rerun
lands on the same curves, secrets and walks.

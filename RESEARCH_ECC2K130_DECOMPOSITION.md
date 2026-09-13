# Decomposing an ECC2K-130 point into a sum of factor-base points

**Experiment:** `scripts/ecc2k130_point_decomposition.py`
**Frozen artefact:** `experiments/ecc2k130_point_decomposition.json`
**Related:** [`RESEARCH_SEMAEV_DECOMPOSITION.md`](RESEARCH_SEMAEV_DECOMPOSITION.md)
(the oracle, built and measured at toy sizes),
[`RESEARCH_KOBLITZ_INDEX_CALCULUS.md`](RESEARCH_KOBLITZ_INDEX_CALCULUS.md)
(the pipeline the oracle feeds),
[`RESEARCH_ECC2K130_EXTENSION.md`](RESEARCH_ECC2K130_EXTENSION.md)
(the same barrier reached by enlarging the field),
[`RESEARCH_QUASI_SUBFIELD.md`](RESEARCH_QUASI_SUBFIELD.md)
(the same barrier reached from the factor-base side).

**The question.**  Index calculus on an elliptic curve needs one thing that
Pollard rho does not: a way to write a target point as a sum of several
points drawn from a small, fixed set.  For ECC2K-130 — `K_0 : y² + xy = x³ + 1`
over `F_2^131`, `#E = 4r` with `r` a 129-bit prime — that means asking, for a
target `R ∈ ⟨G⟩` and a factor base

```text
    F = { P ∈ E(F_2^131) : x(P) ∈ V },      V ⊆ F_2^131 an F_2-subspace, dim V = l,
```

whether `R = P₁ + … + P_m` with every `P_i ∈ F`.

**Bottom line.  Yes, such decompositions exist; yes, they are admissible at
these parameters; and no, none of them can be found for less than `2^131`
operations by any oracle this repository has built or measured.**  Those are
three separate questions with three different answers, and conflating them is
how this route gets over-sold.  Concretely:

- **Existence** is not the obstacle.  The yield law
  `E[#decompositions] = C(|F|, m)/#E` is confirmed on twelve toy cells to
  within a factor `1.38`, and at `n = 131` it says a typical target has a
  decomposition as soon as `m·l ≥ 131 + log₂ m!` — `l = 45` at `m = 3` (§2).
- **Admissibility** is not the obstacle either.  The cofactor class that makes
  odd `m` decompose *nothing* on `K_1/F_2^7` does not bite here: measured on the
  challenge curve, a subspace meeting both trace values gives all four classes
  and one inside `ker Tr` gives the two even ones, and **every `m ≥ 2` is
  admissible either way** (§3).  §4 exhibits certified `m = 2, 3, 4`
  decompositions at full size, with Semaev's `S₄` vanishing on the `m = 3` one.
- **Finding one** is the whole difficulty, and it costs **`m·2^131`** group
  operations for the entire attack — `2^71.77` times the rho reference at its
  best `m` — *independently of `l`*, because a larger factor base needs
  proportionally fewer targets and makes each target proportionally dearer
  (§5).
- Stated as a demand on the oracle rather than as a cost: the decomposition
  oracle would have to beat exhaustive search over **its own candidate set**
  by a factor `2^{70.19 + log₂ m}`, and that number is the same for every `m`
  and every `l` (§5.3).  The two algebraic routes that might have done it are
  already measured and do not (§5.4).

## 0. The boundary, stated before anything is measured

```text
#E(F_2^131) = 4·r,   r = 680564733841876926932320129493409985129   (prime, 2^129.0000)
```

recomputed here from the Koblitz recurrence and matching
`experiments/ecc2k130_extension_field_boundary.json` to the digit.  The
reference is Pollard rho on `⟨G⟩` with the `⟨−1⟩ × ⟨π⟩` speed-up, the same
accounting the rest of the repository uses:

| quantity | value |
|---|---|
| automorphisms on `⟨G⟩` | `2 · 131 = 262` |
| **rho, reference** | **`2^60.8090`**, `S = 0.077430` |

Two boundaries, both derived, per §1 of `AGENTS.md`:

- **The reference:** rho at `2^60.8090`, above.
- **The floor:** the *free-oracle floor*.  Charge one operation per target
  tried and `m·2^{2l}` for the linear algebra, and charge **nothing at all**
  for deciding whether a target decomposes.  Nothing built on an `m`-summand
  subspace factor base can go below that line however good the algebra gets.
  It is `2^56.40` at `m = 4` and falls with `m` (§5.3) — *below* rho, which is
  the point: the floor is not what protects ECC2K-130.  The oracle is.

**Falsification target.**  This thread is a success if a decomposition oracle
is exhibited that answers "is `R` a sum of `m` points of `F`?" in fewer than
`2^{70.19 + log₂ m}`-th of the `C(|F|, m−1)` candidate sub-tuples' cost, at any
`(m, l)` with `m·l ≥ 131`, carrying the whole run — collection, table setup,
verification and linear algebra — below `2^60.81`.  It is abandoned if the
product `relations × targets × oracle` can be shown constant at `2^131` for
every oracle of the fix-some-summands form, which is what §5.1 does.

Inadmissible, per §6 of `AGENTS.md`: dropping the linear algebra from the
budget, quoting a per-target oracle cost as if it were an attack cost,
counting a memory-bought speed-up without pricing the memory against a
generic algorithm given the same memory, and treating the toy yield
measurements as if they were runs at `n = 131`.

## 1. What is measured and what is derived

Worth separating up front, because this note mixes both.

**Measured**, on the real curve over `F_2^131` with its challenge reduction
polynomial `x^131 + x^13 + x² + x + 1` (checked irreducible at startup):

- the cofactor-class histogram of the factor base, and which `m` it admits;
- the correspondence `class parity = Tr(x(P))`, on every sampled point;
- certified `m = 2, 3, 4` decompositions at full size, `S₄` verified.

**Measured at toy sizes** (`K_0/F_2^n`, `n = 11, 13, 17, 19`, every point of
the curve enumerated and the count checked against the Koblitz recurrence):

- the yield law, over twelve `(n, m, l)` cells and 512 targets each.

**Derived**, not measured — and marked as such everywhere below:

- every cost at `n = 131`.  Nothing in §5 is a run; they are operation counts
  from the model in §5.1, whose ingredients (`|F| ≈ 2^l`, yield `C(|F|,m)/#E`,
  oracle `C(|F|, m−split)`, Wiedemann `m·2^{2l}`) are each either measured
  above or standard.

## 2. Decompositions exist — the yield law

The `m`-subsets of `F` number `C(|F|, m)` and their sums are spread over the
whole group, so a fixed target is hit

```text
    λ = C(|F|, m) / #E  ≈  2^{ml} / (m! · 2^n)
```

times on average, and decomposes with probability `1 − e^{−λ}`.  That is the
entire content of "can we decompose": it is a counting statement, and it has
nothing to do with how hard the search is.

Measured on toy Koblitz curves, 512 targets per cell, every target drawn from
`⟨G⟩` as `[h]P` and every reported witness re-summed and checked against its
target:

| `n` | `m` | `l` | `\|F\|` | rate measured | rate predicted | ratio |
|---:|---:|---:|---:|---:|---:|---:|
| 11 | 2 | 6 | 71 | 0.664 | 0.691 | 0.96 |
| 11 | 3 | 4 | 21 | 0.562 | 0.467 | 1.20 |
| 11 | 4 | 4 | 21 | 1.000 | 0.941 | 1.06 |
| 13 | 2 | 7 | 129 | 0.712 | 0.643 | 1.11 |
| 13 | 3 | 5 | 33 | 0.521 | 0.494 | 1.06 |
| 13 | 4 | 4 | 19 | 0.530 | 0.384 | 1.38 |
| 17 | 2 | 9 | 523 | 0.617 | 0.647 | 0.95 |
| 17 | 3 | 6 | 65 | 0.336 | 0.284 | 1.19 |
| 17 | 4 | 5 | 35 | 0.363 | 0.330 | 1.10 |
| 19 | 2 | 10 | 1 005 | 0.643 | 0.619 | 1.04 |
| 19 | 3 | 7 | 139 | 0.602 | 0.567 | 1.06 |
| 19 | 4 | 5 | 29 | 0.049 | 0.044 | 1.10 |

Measured over predicted stays in `[0.95, 1.38]` across rates spanning
`0.049` to `1.000`.  The bias is upward and largest where the base is
smallest (`|F| = 19`), which is where the Poisson approximation is worst; the
law is right to about twenty per cent and that is all §5 needs of it.

Carried to the challenge parameters, `λ ≥ 1` needs `m·l ≥ 131 + log₂ m!`:

| `m` | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---:|---:|---:|---:|---:|---:|---:|
| saturating `dim V` | 66.00 | 44.53 | 33.90 | 27.58 | 23.42 | 20.47 | 18.29 |

So a 45-dimensional subspace of `F_2^131` already gives most targets a
3-decomposition, and a 34-dimensional one gives most targets a
4-decomposition.  **There is no shortage of decompositions at ECC2K-130.**

## 3. They are admissible — the cofactor class, measured on the challenge curve

Existence in the *whole* group is not enough.  `E(F_2^131)` is cyclic of order
`4r` — the only rational 2-torsion point is `(0, √b)` — so every point has a
class in `E/⟨G⟩ ≅ Z/4`, computed as `[r]P`, and a decomposition of a target in
`⟨G⟩` needs the summands' classes to sum to zero mod 4.  This is not a
formality: on `K_1/F_2^7` the whole base sits off `⟨G⟩` and **odd `m`
decomposes nothing at any `|F|`** (`RESEARCH_KOBLITZ_INDEX_CALCULUS.md`, "the
cofactor class decides which `m` can work").

Measured here by computing `[r]P` for 256 sampled base points, in lockstep
with Montgomery batch inversion, on two 45-dimensional subspaces.  The
abscissae are sampled rather than enumerated — `2^45` inversions is not a
thing to do in this language — and the measured point density in the subspace,
`0.494` and `0.452`, is the check that the sample is a factor base and that
`|F| ≈ 2^l`:

| `V` | class 0 | 1 | 2 | 3 | admissible `m` |
|---|---:|---:|---:|---:|---|
| `span{1, z, …, z^44}` | 60 | 64 | 68 | 64 | every `m ≥ 2` |
| `span{z, z², …, z^45}` ⊂ `ker Tr` | 128 | 0 | 128 | 0 | every `m ≥ 2` |

Two things fall out, and the second explains the first.

**The class parity is the trace of the abscissa.**  `P ∈ 2E` iff
`Tr(x(P)) = Tr(a) = 0`, so the class of `P` is even exactly when `Tr(x(P)) = 0`
— checked point by point on every sample, no exceptions.  For this field's
reduction polynomial `Tr(z^i) = 1` only at `i ∈ {0, 129}`, so `span{1, …, z^44}`
meets both parities and `span{z, …, z^45}` lies inside `ker Tr` and meets only
the even classes.  The zeros in the odd columns are that prediction; the even
split landing at exactly 128/128 is the draw.

**Neither choice is an obstruction at `n = 131`.**  With all four classes
present every `m ≥ 2` works; with only `{0, 2}` present, `2 + 2 ≡ 0`, so every
`m ≥ 2` still works.  The degeneracy that kills odd `m` at `n = 7` needs a
base whose classes generate a proper subgroup *avoiding* zero, and no
subspace of `F_2^131` produced one.

## 4. A witness at full size

Existence arguments are cheap; here are the points.  Random abscissae drawn
from `V = span{1, z, …, z^44}` (dimension 45), the summands' classes made to
cancel, everything re-checked on the curve:

```text
m = 3,  V = span{1, z, …, z^44} ⊂ F_2^131

  x₁ = 0x17e0fcba1baf
  x₂ = 0x1ef434ab8c49
  x₃ = 0x1bf1faf66546
  R  = (0x5f51ebb34b6499bab2900fa7c032d9a06,
        0x2894e97000896584123c8ff637b9d1b81)

  P₁ + P₂ + P₃ = R,   R on the curve,   [r]R = O  (so R ∈ ⟨G⟩)
  S₄(x₁, x₂, x₃, x(R)) = 0
```

`m = 2` and `m = 4` witnesses are in the frozen artefact alongside it.  The
`S₄` check is the one that matters for index calculus: the polynomial whose
roots the oracle would have to find really does vanish on a genuine
decomposition at the challenge parameters, so the algebraic formulation is
sound here and the difficulty is not a modelling artefact.

This is a *planted* decomposition — the summands were chosen and the target
computed from them, not the other way round.  That is precisely the asymmetry
the rest of the note prices.

## 5. Finding one — the cost surface

### 5.1 The product law

The model, with every phase priced (`AGENTS.md` §5):

| phase | cost |
|---|---|
| relations needed | `\|F\| = 2^l` — no Frobenius orbit collapse is available (§6) |
| targets per relation | `2^n / C(\|F\|, m)` |
| oracle per target | `C(\|F\|, m−1)`: walk the sub-tuples, root-find the last summand inside `V` |
| linear algebra | `m·2^{2l}`, sparse block Wiedemann on an `\|F\| × \|F\|` matrix mod `r` |

The oracle row is the one this repository has actually built and measured:
`RESEARCH_SEMAEV_DECOMPOSITION.md` turns `S₄` into a quartic in the last
variable and intersects its roots with the subspace via
`gcd(q, L_V mod q)`, so a candidate sub-tuple costs `O(l)` field operations
and no search over `F` at all.  It is `926×` faster than enumeration at
`l = 8`.  Multiply the column out anyway:

```text
    2^l  ·  2^n / C(|F|, m)  ·  C(|F|, m−1)  =  m · 2^n
```

**independent of `l`** — because `C(|F|,m−1)/C(|F|,m) ≈ m/|F|`, and the `|F|`
cancels the factor base.  Measured as flatness across the usable range, at
`m = 3`:

| `dim V` | 12 | 20 | 28 | 36 | 44 | 52 |
|---|---:|---:|---:|---:|---:|---:|
| targets per relation | `2^97.58` | `2^73.58` | `2^49.58` | `2^25.58` | `2^1.58` | `2^0` |
| oracle per target | `2^23.00` | `2^39.00` | `2^55.00` | `2^71.00` | `2^87.00` | `2^103.00` |
| **total** | `2^132.58` | `2^132.58` | `2^132.58` | `2^132.58` | `2^132.58` | `2^155.00` |

Flat to the last digit until the base outgrows what the method needs, then
rising.  This is the same law `RESEARCH_SEMAEV_DECOMPOSITION.md` states and
measures at toy sizes — its `projected collection / 2^n` column is flat within
`1.5×` while `2^n` grows four-thousand-fold — instantiated at `n = 131`.

### 5.2 The table, one unit, one ratio column

`log₂` group operations, everything inside: targets tried, oracle work, table
setup, verification and linear algebra.  Every row is **derived**; none is a
run.  `dim V` is the operating point, the saturating dimension from §2 for the
memory-free rows and the cost-minimising one for the rest.

| variant | `dim V` | table entries | `log₂` ops | ratio to rho | class |
|---|---:|---:|---:|---:|---|
| **Pollard rho with `⟨−1⟩ × ⟨π⟩`** *(reference)* | — | `O(1)` | **60.81** | `2^0` | reference |
| free-oracle floor, `m = 4` *(floor)* | 26.83 | — | 56.40 | `2^−4.41` | floor |
| `m = 2`, built oracle | 66.00 | — | 133.58 | `2^+72.77` | derived |
| **`m = 3`, built oracle** | 44.53 | — | **132.58** | **`2^+71.77`** | derived |
| `m = 4`, built oracle | 33.90 | — | 133.00 | `2^+72.19` | derived |
| `m = 5`, built oracle | 27.58 | — | 133.32 | `2^+72.51` | derived |
| `m = 6`, built oracle | 23.42 | — | 133.58 | `2^+72.77` | derived |
| `m = 2`, pair table, unbounded memory | 43.25 | `2^85.50` | 89.36 | `2^+28.55` | derived |
| `m = 3`, pair table, unbounded memory | 43.50 | `2^86.00` | 90.58 | `2^+29.77` | derived |
| `m = 4`, triple table, unbounded memory | 27.50 | `2^79.92` | 81.29 | `2^+20.48` | derived |
| `m = 4`, all 4-subset sums tabulated | 20.00 | `2^75.42` | 76.50 | `2^+15.69` | generic in costume |
| best cell, `≤ 2^40` entries | 8.00 | `2^38.51` | 100.49 | `2^+39.68` | generic in costume |
| best cell, `≤ 2^50` entries | 8.00 | `2^48.70` | 90.30 | `2^+29.49` | generic in costume |
| best cell, `≤ 2^60` entries | 9.25 | `2^58.70` | 81.55 | `2^+20.74` | generic in costume |
| best cell, `≤ 2^70` entries | 10.50 | `2^68.70` | 72.88 | `2^+12.07` | generic in costume |

The best memory-free cell is `2^132.58`, `2^71.77×` rho.  That is within a bit
and a half of `2^70.19×`, the cheapest cell in the extension-field note's large
horn — which is the same fact from the other side, and the reason it is worth
saying plainly:

> **Index calculus is priced by the ambient field; rho is priced by the
> subgroup.**  ECC2K-130 has `r ≈ 2^129` filling a `2^131` field, so there is
> no slack between them to trade.

**Memory does not rescue it, and the rows that look like it does are not doing
index calculus.**  Tabulating every `m`-subset sum makes the table a baby-step
table, so the honest comparison is against a generic algorithm handed the same
memory.  Baby-step giant-step on `⟨G⟩` with `2^μ` stored steps costs
`max(2^μ, r/2^μ)`:

| table entries | best decomposition cell | BSGS, same memory | rho, no memory |
|---:|---:|---:|---:|
| `2^30` | `2^110.08` | `2^99.00` | `2^60.81` |
| `2^40` | `2^100.49` | `2^89.00` | `2^60.81` |
| `2^50` | `2^90.30` | `2^79.00` | `2^60.81` |
| `2^60` | `2^81.55` | `2^69.00` | `2^60.81` |
| `2^70` | `2^72.88` | `2^70.00` | `2^60.81` |

At **every** budget the generic algorithm is cheaper than the decomposition
method given the same table, and rho is cheaper than both while storing
nothing.  There is no memory budget at which buying a table turns
decomposition into a win.

### 5.3 What the oracle would have to cost

Turn the product law around.  At the `(m, l)` that minimises the free-oracle
floor, an oracle costing `2^w` per target puts the whole run at
`2^{targets + w}`; solve for the `w` that reaches rho, and compare it with the
cost of simply searching the oracle's own candidate set:

| `m` | `dim V` | floor | search space per target | oracle budget | required speed-up over search |
|---:|---:|---:|---:|---:|---:|
| 2 | 43.33 | `2^89.25` | `2^43.33` | `2^−27.86` | `2^−71.19` |
| 3 | 33.00 | `2^68.58` | `2^65.00` | `2^−6.78` | `2^−71.78` |
| 4 | 26.83 | `2^56.40` | `2^77.91` | `2^+5.71` | `2^−72.19` |
| 5 | 22.76 | `2^48.44` | `2^86.46` | `2^+13.94` | `2^−72.51` |
| 6 | 19.89 | `2^42.85` | `2^92.54` | `2^+19.77` | `2^−72.78` |
| 8 | 16.12 | `2^35.61` | `2^100.54` | `2^+27.35` | `2^−73.19` |

Three readings, in increasing order of usefulness.

**The floor is not what protects the curve.**  From `m = 4` up it sits below
rho — `2^56.40` against `2^60.81` — so a free decomposition oracle *would*
break ECC2K-130 by index calculus.  Anyone claiming the method is
structurally impossible here is claiming more than the counting argument
supports.

**The budget is absurd in absolute terms.**  At `m = 4` the oracle gets
`2^5.71 ≈ 52` operations to decide whether a target is a sum of four points
whose abscissae lie in a 27-dimensional subspace.

**And the last column is the real invariant.**  The required speed-up over
exhaustive search is `2^{−(70.19 + log₂ m)}` — it barely moves with `m`, and it
does not move with `l` at all, because it *is* the product law:

```text
    required speed-up  =  rho / (m · 2^131)  =  2^{-(70.19 + log₂ m)} .
```

Choosing a different `m`, a different `l`, a bigger factor base or a smaller
one moves work between the columns and never changes that number.  **Every
tuning knob on this attack is inside the identity, not outside it.**

### 5.4 The two routes that might have moved it are already measured

Both are in `RESEARCH_SEMAEV_DECOMPOSITION.md`, and both were measured
against exactly the `2^{(m−1)l}` boundary the table above charges:

- **Gröbner (F4/F5) on the descended system with one summand fixed.**  The
  deciding degree is *not* bounded: its median runs `5, 6, 6, 7, 7` over
  `l = 3…7`, a least-squares slope of `0.500`, so `D ≈ l/2`.  The Macaulay
  matrix then has `≈ 2^{1.81l}` columns and its elimination alone is
  `≈ 2^{3.6l}` at `ω = 2` — past the `2^{2l}` enumeration boundary for a
  *single* fixed summand, before the sweep is paid at all.  The measured total
  is `2^{6.26·l}` — an overestimate of the asymptote at those sizes, as that
  note says, but the deciding-degree slope is what settles it and `0.5` settles
  it against.  This is the wrong direction by a factor exponential in `l`.
- **SAT (CDCL) on the Weil-descended system.**  One conflict per candidate
  tuple, no pruning; a strictly worse enumerator, and the gap widens with `l`
  — `5×` at `l = 5`, `62×` at `l = 8`.

So the two candidates for the `2^{−72}` are both measured, and both go the
wrong way.  What is *not* ruled out is an oracle that is not of the
fix-some-summands form at all; this note prices the family, not the universe.

## 6. The Frobenius saving is not available at `n = 131`

GGMP's headline is that a **Frobenius-stable** factor base collapses `|F|`
unknowns into `|F|/n` orbit unknowns: `≈ n` off relation collection and `≈ n²`
off the linear algebra.  None of that is available here.

A Frobenius-stable `F_2`-subspace of `F_2^n` is a binary cyclic code of length
`n`, so its dimension is a subset sum of the 2-cyclotomic coset sizes mod `n`.
At `n = 131`, `2` is a primitive root, `ord_131(2) = 130`, there are exactly
two cosets of sizes `1` and `130`, and therefore

```text
    available invariant dimensions at n = 131:   0, 1, 130, 131
```

— so the only invariant factor bases are `E(F_2)`, four points, and one of
size `2^130`: nothing in between.  This is the
same dichotomy `RESEARCH_ECC2K130_EXTENSION.md` proves survives every base
change, and the census in `RESEARCH_QUASI_SUBFIELD.md` finds no non-subfield
quasi-subfield polynomial at `n = 131` either.

It would not have been enough in any case.  Grant both savings for free on top
of the memory-free row:

```text
    m · 2^131 / 131²  =  2^{118.52}  =  2^{57.71} × rho .
```

A factor of `n²` — `2^14.07` — against a gap of `2^71.77`.

## 7. What this does not settle

- **The oracle family, not every oracle.**  §5 prices methods that fix some
  summands and solve or look up the rest.  A decomposition method of a
  different shape is outside it, and §5.3's invariant is a restatement of the
  product law, not a lower bound on all algorithms.
- **The linear-algebra model is standard, not measured here.**  `m·2^{2l}`
  sparse Wiedemann; it is never the binding term in any row above, so the
  verdict does not rest on it.
- **The yield law is measured at `n ≤ 19` and applied at `n = 131`.**  That is
  an extrapolation of a counting argument, not of a timing, and §4's witnesses
  are a direct check of its conclusion at full size — but the *rate* at
  `n = 131` is not measured and cannot be.
- **Non-subspace factor bases are untouched.**  Every base here is
  `{P : x(P) ∈ V}` for an `F_2`-subspace `V`, which is what makes the subspace
  polynomial `L_V` and the whole root-finding shortcut work.  Diem's and
  Gaudry's frameworks both use them; a base of another shape is a different
  question.
- **The cofactor measurement samples, it does not exhaust.**  256 points per
  subspace at `n = 131`; the parity-equals-trace correspondence is a theorem
  and holds on every sample, but "all four classes occur" is a measurement on
  a sample, not a proof about the whole base.
- **The `2^131` is an upper bound on this family's cost and a statement about
  today's oracles**, not a security proof.  §5.3 says exactly what would
  falsify it, and the number is written down.

## References

- **I. Semaev**, *Summation polynomials and the discrete logarithm problem on
  elliptic curves*, ePrint 2004/031 — the `S_{m+1}` this is all built on.
- **P. Gaudry**, *Index calculus for abelian varieties of small dimension and
  the elliptic curve discrete logarithm problem*, J. Symbolic Comput. 44
  (2009) — the pairs-and-solve structure.
- **C. Diem**, *On the discrete logarithm problem in elliptic curves*,
  Compositio Math. 147 (2011) — factor bases as subspaces.
- **S. Galbraith, R. Granger, S.-P. Merz, C. Petit**, *On index calculus
  algorithms for subfield curves*, SAC 2020 (ePrint 2020/1315) — the
  Frobenius-invariant factor bases §6 shows are unavailable at 131.
- **C. Petit, J.-J. Quisquater**, *On polynomial systems arising from a Weil
  descent*, ASIACRYPT 2012 — the first-fall-degree route §5.4 measures.
- **D. Bailey et al.**, *Breaking ECC2K-130*, ePrint 2009/541 — the challenge,
  and the rho attack the reference column is.
- **P. L. Montgomery**, *Speeding the Pollard and elliptic curve methods of
  factorization*, Math. Comp. 48 (1987) — the batch inversion §3 runs on.

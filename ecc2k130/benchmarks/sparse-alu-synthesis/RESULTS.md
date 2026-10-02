# Sparse-basis ALU linear-map synthesis: negative static result

Decision: **none of the bounded ALU/synthesized map families survives; do not
implement or run a GPU/full-walk candidate.**  Even an impossible oracle that
assigns the cheapest map `L_3` to every update, makes every word shift and
dynamic selection free, and charges only one fused `LOP3` per nonzero masked
word contribution needs **1,826.125 ALU operations/update** for the required
maps.  The entire 26 B/s walk has room for only **1,124.529 ALU
lane-instructions/update** on the favourable 188-SM, 2.430-GHz, 64-lane model.

This is a native static circuit result.  No CUDA compilation, GPU run, sparse
walk, search, distinguished-point collection, collision recovery, or key
recovery was performed.

## Exact checks and boundaries

The checker reconstructs the selected sparse field and all eleven required
linear maps from the same beta root and repository transforms as the measured
table experiment.  It passes:

- all 131 basis vectors;
- 4,096 deterministic dense vectors;
- 32,768 dense `L_j` comparisons against both sparse Frobenius and repository
  beta/normal arithmetic;
- ranks 131 for all conversion maps and 130 for every `L_j`; and
- 59 exact composition identities covering `L_(a+b)`, `L_(2a)`, and
  `L_(j+1)` for powers through ten.

The complete-walk ALU budgets are:

| boundary | ALU lane-instructions/update |
|:--|--:|
| confirmed fused 15.436677 B/s | 1,894.044942 |
| 26 B/s objective | **1,124.529231** |

## Generated masked-word circuits

Each nonzero `(diagonal, output-word)` term can be accumulated with one
`LOP3(out, shifted, mask)`.  The exact circuit also needs one shifted word per
term except on the zero diagonal.  “Free shift” below is deliberately
impossible and gives the candidate the most favourable count.

| map | matrix ones | diagonals | word terms / free-shift `LOP3` | exact shift + `LOP3` |
|:--|--:|--:|--:|--:|
| sparse-to-normal | 8,541 | 261 | **754** | 1,503 |
| sparse-to-beta | 8,559 | 259 | 751 | 1,497 |
| beta-to-sparse | 8,551 | 259 | 755 | 1,505 |
| `L_3` | 1,474 | 223 | **489** | 973 |
| `L_4` | 2,567 | 233 | 569 | 1,133 |
| `L_5` | 3,894 | 243 | 643 | 1,281 |
| `L_6` | 5,302 | 254 | 687 | 1,369 |
| `L_7` | 6,511 | 258 | 724 | 1,443 |
| `L_8` | 7,261 | 257 | 729 | 1,453 |
| `L_9` | 7,790 | 256 | 730 | 1,455 |
| `L_10` | 8,266 | 257 | 746 | 1,487 |

The impossible oracle map ledger is therefore:

```
normal + 2*L3 + (toBeta + fromBeta)/16
= 754 + 2*489 + (751+755)/16
= 1,826.125 ALU/update.
```

That map-only count is 1.624 times the total 26 B/s ALU budget.  Its isolated
ALU ceiling is 16.010821 B/s even though all shifts, selection, arithmetic,
state, control, and reporting are free.  The exact masked-word circuit costs
3,636.625 operations/update.

## Dynamic selection and common structure

For uniform random field coordinates, the actual Hamming selector is close to
uniform but not assumed exactly uniform.  Its exact binomial model predicts
**7.874659 distinct `j` branches per 32-lane warp**.  A switch over fixed
circuits therefore serializes essentially all eight arms.

Across the eight `L_j` matrices:

- 769 `(diagonal,word)` positions occur in at least one map;
- there are 5,317 map memberships;
- those positions carry 5,010 distinct nonzero masks;
- zero masks are identical across all eight maps; and
- all 1,048 output rows are distinct.

Sharing every shifted word, materializing all eight outputs, and selecting
branchlessly still costs 13,922.625 operations/update.  The warp-divergent
fixed-circuit model costs 22,532.377.  Perfect block compaction with free
compaction, free shifts, and the real selector distribution costs 2,162.867.

## Greedy XOR CSE

A deterministic Paar-style common-subexpression heuristic was run with every
bit extraction, placement, and intermediate-liveness cost omitted and every
bit-level XOR optimistically priced as one GPU ALU instruction:

| synthesized outputs | temporaries | total XORs |
|:--|--:|--:|
| sparse-to-normal | 953 | 2,999 |
| fixed `L_3` | 200 | 735 |
| normal and `L_3` jointly on X | 1,114 | 3,702 |
| joint normal/`L_3` on X plus `L_3` on Y | | **4,437** |

The normal circuit alone exceeds both complete-walk budgets and exposes 1,084
signals; the joint circuit exposes 1,245.  This heuristic is not a universal
XOR lower bound, but it supplies no implementable candidate and confirms that
the matrices do not contain the common structure the diagonal density might
have hidden.

## Frobenius composition

The checker verifies the identities

```
L_(a+b) = L_a + sigma^a L_b
L_(2a)  = L_a composed with L_a
L_(j+1) = sigma L_j + L_1.
```

They preserve correctness but do not reduce complete work in the sparse
basis.  One sparse ALU square is five 15-operation spreads plus the 82-op
two-fold reducer: 157 source operations.  Even perfect `j` compaction averages
13 squares/update across X and Y, for 3,638.125 map operations after the
normal and amortized basis maps.  Keeping those map squares on `CLMAD` costs
65 `CLMAD`s/update and has a 14.056615 B/s isolated ceiling before the walk's
ordinary field products.  The branchless `3+1+2+4` chain is worse.

## Complete comparable B16 source ledger

| row | maps | sparse reductions | ALU-square spread | comparable total |
|:--|--:|--:|--:|--:|
| reference beta selector/direct-reduction stages | 1,230.000 | 1,046.250 | — | **2,276.250** |
| impossible oracle sparse candidate | 1,826.125 | 476.625 | 75.000 | **2,377.750** |

The candidate is already **4.459% worse** in the comparable source unit while
receiving free shifts, free selection, and `L_3` for every update.  Its
source-ledger ALU ceiling is 12.296398 B/s.  Unchanged products, inverse
transforms, state traffic, control, and reporting are outside both rows, so
adding them cannot reverse the decision.  Source operations are not SASS or a
measured rate; the stronger target-budget failure above does not require this
translation.

## Scope of the stop

This closes the generated masked-diagonal, all-output shared-shift, divergent
fixed-map, greedy XOR-CSE, and repeated-Frobenius families for this sparse
basis.  It does not prove a lower bound for every possible GF(2) circuit.  A
future route must present a materially different circuit whose complete map
cost is below 1,124 ALU/update before any implementation or GPU request is
justified.

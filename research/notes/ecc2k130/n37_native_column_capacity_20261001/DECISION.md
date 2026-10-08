# n37 fixed-leaf native base: six summands are the first capacity-admissible arity at 42 log columns

**Decision: `NATIVE_K42_UP_TO_FIVE_SUMMANDS_CAPACITY_NO_GO`.** This is an exact
necessary support bound, not a natural-target PDP measurement or an attack
speed result. It applies to a single fixed descendant-native base with at
most 42 distinct nonzero projected sign classes in the prime subgroup and no
additional cheap same-curve log action. It permits repeated factors, either
sign, the identity, and decompositions with *at most* the stated number of
summands. The matching source control is the verified n37/L1,024 compact
batch's `K=42` log columns and subgroup order `r=230603167` in the
[frozen input](../disjoint_cold_v2_20261001/FROZEN.json). The degree-73 map,
equal source/leaf group order `137439487532`, and prime subgroup are pinned
by the [repaired archive](../../../koblitz_isogeny_descent_37_results_20260925/RESULTS.md).

Every sum of at most `m` signed factors from 42 classes determines an integer
coefficient vector `c` in `Z^42` with `||c||_1 ≤ m`. Choosing `s` nonzero
coordinates, their signs, and positive magnitudes with total at most `m`
gives exactly

```text
V(K,m) = sum_{s=0}^{min(K,m)} 2^s * C(K,s) * C(m,s).
```

Different vectors may still produce the same group point, so the number of
supported targets is **at most** `min(V(K,m),r)`. For a uniformly selected
prime-subgroup target, its support probability is at most that count divided
by `r`. This is a rigorous ceiling for *any* choice of 42 native points under
the stated action and arity. It is not a claim that the sums are uniformly
distributed or that an algebraic solver can find them.

| Maximum summands `m` | Formal sums `V(42,m)` | Capacity ratio `V/r` | Uniform-target support ceiling | First `K` whose counting ratio reaches 1 |
| ---: | ---: | ---: | ---: | ---: |
| 2 | 3,613 | 0.000015667608 | 0.00157% | 10,738 |
| 3 | 102,425 | 0.000444161289 | 0.04442% | 557 |
| 4 | 2,179,241 | 0.009450178106 | 0.94502% | 136 |
| 5 | 37,129,037 | 0.161008356837 | **16.10084%** | 61 |
| 6 | 527,810,725 | 2.288826870275 | ≤100%; yield unknown | **37** |
| 7 | 6,440,955,121 | 27.930904873480 | ≤100%; yield unknown | 26 |

Thus at `K=42`, at least **83.89916%** of uniformly sampled subgroup points
cannot be expressed using five or fewer native factors, even with a perfect
PDP solver. Six summands are merely the first arity that escapes this counting
no-go. Increasing the native base to at least 61 columns also escapes the
five-summand count, but that changes the rank and index budget. The source's
signed-Frobenius orbit can represent 74 physical points per complete log
column, or 3,108 physical points for 42 full orbits; a fixed leaf has no
corresponding known cheap Frobenius action. A transported source orbit keeps
its source log relations by mapping *each* member; it is a different base
policy from choosing 42 native leaf sign classes.

The archive also certifies a useful pullback fact. Since `73` does not divide
the common group order (`137439487532 mod 73 = 67`), the degree-73 isogeny's
rational kernel is trivial. Equal source/target cardinalities make its
restriction to rational points bijective, so every rational leaf point has a
unique rational source preimage. This proves an exact pullback policy exists;
it does **not** price or implement one. An explicit inverse-by-rational-map
implementation, with costs charged, is the next pullback prerequisite.

**Next gate.** Freeze new, orbit-disjoint natural Q and compare the same
`K=42` useful log-column budget across original source, transported source,
native leaf, and exact pullback bases. The native arm must start at `m=6`
under this fixed-column model; source/transported and pullback may use their
own fixed arities, with every extra physical factor, map evaluation, inverse,
rank row, failed PDP attempt, and verified scalar charged. First measure
natural support and full-rank admission, then run cold single-target and
L≥1,024 batches against the same-Q strong signed-Frobenius rho. Keep proved
UNSAT separate from UNKNOWN/resource exits. A ratio from this counting note
cannot be promoted to a solver or end-to-end speedup, and n131 transfer
remains unestablished.

[BOUND.json](BOUND.json) contains every integer and input hash. Its producer
[bound.py](bound.py) uses the closed form above; [verify.py](verify.py)
recomputes the coefficients by independent dynamic programming and checks
the source hashes, first-admissible thresholds, and decision. Reproduce with:

```sh
python3 research/notes/ecc2k130/n37_native_column_capacity_20261001/bound.py \
  --out /tmp/n37-native-column-bound.json
cmp /tmp/n37-native-column-bound.json \
  research/notes/ecc2k130/n37_native_column_capacity_20261001/BOUND.json
python3 research/notes/ecc2k130/n37_native_column_capacity_20261001/verify.py \
  --out /tmp/n37-native-column-verify.json
cmp /tmp/n37-native-column-verify.json \
  research/notes/ecc2k130/n37_native_column_capacity_20261001/VERIFY.json
```

Class: **accounting / capacity**. Full-DLP total cost, operation-normalized
`S`, rho/floor ratios, natural-target yield, and speedup remain unset.

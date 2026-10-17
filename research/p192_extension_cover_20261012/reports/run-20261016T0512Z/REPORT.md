# The complete P-192 Type I quotient base as reduced genus-three classes

**Result.** All 9,784 new quotient coordinates in the frozen P-192 Type I
list have fixed-origin, nonspecial degree-three representatives on the
descended smooth plane quartic. The earlier 141 representatives and these
9,784 give **9,925 distinct reduced Jacobian classes**; their quotient
points are also distinct after sign reversal. The checked Sage producer
constructed both cubic cuts for every new class. An independent PARI/GP
replay verified their residual resultants, point coordinates, squarefree
supports, and numerator interpolation. The [receipt audit](../../type_i_line_new_all_receipt_verification.json)
checks contiguous chunk custody, every quotient-coordinate match and
nonspeciality flag, and spaced subgroup and character controls.
The new reduced classes occupy the four \((Q,H)\) character cells in counts
`++` 2,409, `+-` 2,435, `-+` 2,424, and `--` 2,516.

The complete source-line geometry gives **4,891 independent three-sum
collision rows**. The 4,892 original split quartic lines form a connected
graph on the 141 older points. Their 9,784 three-sum presentations have 141
distinct targets and 9,643 repeated decompositions. The
[graph certificate](../../type_i_line_full_pencil_graph.json) includes a
4,891-edge spanning tree of the line-adjacency graph, which supplies an
explicit triangular collision basis. An [independent PARI replay](../../type_i_line_full_pencil_gp_verification.log)
checks all 4,892 quotient line-section sums and the 141 star targets.

| Exact source-pencil quantity | Count |
| --- | ---: |
| Older reduced classes | 141 |
| Newly reduced classes | 9,784 |
| Distinct signed quotient classes in the resulting base | 9,925 |
| Original split source lines | 4,892 |
| Three-sum decompositions supplied by those lines | 9,784 |
| Distinct source-pencil targets | 141 |
| Repeated pencil decompositions | 9,643 |
| Independent pencil collision rows | 4,891 |

The 141 pencil targets reach all four geometric character cells: `++` 34,
`+-` 39, `-+` 29, and `--` 39. Counting decompositions with multiplicity
gives 2,306, 2,742, 2,028, and 2,708 respectively. Each older point is
incident to between 51 and 86 source lines.

## Fixed-origin reduction and native replay

For each new quotient point \(Q\), the checked Sage wrapper computes its
prime-order half \(H=[2^{-1}\bmod r]Q\) and applies the frozen inverse-cover
constructor to both fibers. It obtains two degree-three residual divisors
\(A=R_Q\) and \(B=R_H\). With the same effective origin \(T\) used for the
older 141 classes, it finds cubic sections cutting
\(G=B+T+U\) and \(F=A+U+D\). Their equal degree-twelve line-bundle
classes give \([D-T]=[B-A]\). A squarefree degree-three support with a
nonzero quadratic term in its interpolated \(y_D(x)\) is noncollinear and
has \(h^0(D)=1\) on a smooth plane quartic; thus \(D\) is the unique
effective representative of that fixed-origin class.

The quotient-coordinate identification uses the previously proved generic
three-quotient pullback relation: the two fiber sums at \(Q\) and its half
\(H\) push forward with difference \((Q,Q,Q)\). The original 141 classes
have a complete independent quotient-image replay. For the new list, the
[spaced image control](../../type_i_line_new64_spaced_quotient_images.json)
recomputes all three pushforwards on 64 indices spanning the frozen order,
and a [last-index control](../../type_i_line_new_last_quotient_image.json)
covers the short final chunk;
the remaining new coordinates use the same algebraic cover identity and
checked constructor. This is a sampled coordinate replay, distinct from
the exhaustive cubic-cut replay.

The [resumable driver](../../line_new_quotient_all_driver.py) ran the checked
repository Sage launcher in 20 disjoint chunks, retaining the accepted
runtime receipt and each producer log, stderr, and exact divisor record.
The [native serializer](../../prepare_line_new_all_fixed_origin_gp.py) audits
every chunk hash and input index before writing the
[GP case list](../../type_i_line_new_all_fixed_origin_input.gp). The
[native verifier](../../verify_line_new_all_fixed_origin_native.gp) checks
both cubic cuts per class: the degree-twelve resultant, known-factor
division, residual \(x\)- and \(y\)-coordinates, squarefreeness, disjoint
support, and rank-nine numerator interpolation. Its pass marker is
`NEW_ALL_LINE_FIXED_ORIGIN_NATIVE_PASS cases=9784 prime_bits=192`.

## Why the pencil rank does not yet constrain older logarithms

Let a source line \(e=(i,j)\) have private new residual classes
\(N_{e,1},N_{e,2}\). Its quotient line-section identity is

\[
Q_i+Q_j+N_{e,1}+N_{e,2}=S,
\]

with the same \(S\) for every source line. Fixing a vertex \(i\) makes
\(Q_j+N_{e,1}+N_{e,2}=S-Q_i\) a three-sum presentation, one for each
incident line. The source graph is connected, so comparisons of incident
lines connect all 4,892 line-section vectors. Each vector has a private
new-point column, hence the vectors are independent. Differences along
the recorded spanning tree are triangular in discovery order and have
rank \(4892-1=4891\) over the P-192 subgroup order.

The same private columns explain the useful-rank boundary. For arbitrary
values on the 141 older classes, an arbitrary common \(S\), and one
arbitrary member of each new pair, the other member is uniquely determined
by its line identity. Thus the pencil rows leave all 141 older logarithms
free. The 9,925-column pencil system has nullity \(9925-4891=5034\);
after fixing the old columns, 4,893 new/common-value degrees remain.
Connecting private residuals through independent line cuts or another
divisor relation is the next rank-bearing experiment.

## Known-scalar collector and scale

The [205-class collector run](../run-20261016T0432Z/REPORT.md) exactly
enumerated 1,456,935 three-sums, found 1,456,904 distinct targets, and
compared a `++` filter with a balanced four-cell query law. Every accepted
query retained the equation \(T_a=-2[a]G\). The four-cell law made all
1,456,904 targets eligible versus 364,036 for `++`, and filled its frozen
64-query quota in 76 rather than 241 raw scalar trials. Both arms had zero
matches and zero query-row rank. The four-cell supports are nearly equal,
so the change improves coverage per *raw scalar trial*; it does not create
a comparable gain per already accepted query.

For the complete 9,925-class base, the exact unordered three-sum capacity
is \(\binom{9927}{3}=162993977475\). This is a cap, not a full support
enumeration. The [Lean integer certificate](../../lean/LineSupport.lean),
combined with the all-cell accepted-scalar lower bound, places an
ideal-uniform conditional hit probability below \(2^{-152}\) per accepted
query even at this capacity. The observed pencil collisions are exact
algebraic relations, but their private columns keep the older-log rank at
zero. Collector wall time, relation yield outside these structured pencil
identities, final linear algebra, and unseen-target descent require their
own measured run.

## Proof and replay map

| Result | Verification |
| --- | --- |
| Two cubic sections preserve \([D-T]=[B-A]\) | `fixed_origin_cubic_class_identity` in [CrossRatio.lean](../../lean/CrossRatio.lean), plus independent native checks of all concrete sections. |
| Quotient image of the fiber difference | `three_quotients_inverse_on_odd_torsion` in [CrossRatio.lean](../../lean/CrossRatio.lean) under explicit push-pull hypotheses; complete prior 141-class replay, 64 spaced new-class controls, and the final index. A Lean model of the concrete quotient maps remains open. |
| Nonspecial degree-three representatives | Native squarefreeness and curve-incidence replay for all 9,784; receipt checks nonzero quadratic \(y_D\)-coefficient for each. The Riemann--Roch bridge to the concrete quartic remains a written proof. |
| Shared-point pencil three-sum equality | `pencil_three_sum_eq_of_line_sections` in [PencilThreeSum.lean](../../lean/PencilThreeSum.lean), with all 4,892 quotient identities checked by Python and PARI. |
| Older assignments extend across private new pairs | `line_private_pairs_extend` in [LineSupport.lean](../../lean/LineSupport.lean). |
| Full pencil rank 4,891 | Connected source graph, explicit 4,891-row spanning-tree certificate, and private pivots. The graph-to-rank bridge is a written linear-algebra proof; a Lean statement for the concrete finite graph remains open. |
| Four-cell known-scalar density bound | `four_cell_character_lower_bound` in [FourCellDensity.lean](../../lean/FourCellDensity.lean), with its concrete character-sum hypotheses supplied by the frozen density receipt. |

The cited Lean modules built locally, were audited with `#print axioms`,
and passed fresh `leanchecker` replay. Their structural hypotheses are
explicit; the concrete finite-field quartic, its Riemann--Roch reduction,
and the 4,892-line graph are retained in coordinate certificates rather
than encoded as one Lean theorem.

Replay from the repository root:

```sh
/Volumes/SSD990/cryptanalysis/sage --runtime-info
python3 research/p192_extension_cover_20261012/line_new_quotient_all_driver.py --workers 4 --chunk-size 512
/Volumes/SSD990/cryptanalysis/sage -python research/p192_extension_cover_20261012/verify_line_new_quotient_images.py --start 0 --count 64 --stride 153 --out research/p192_extension_cover_20261012/type_i_line_new64_spaced_quotient_images.json
/Volumes/SSD990/cryptanalysis/sage -python research/p192_extension_cover_20261012/verify_line_new_quotient_images.py --start 9783 --count 1 --out research/p192_extension_cover_20261012/type_i_line_new_last_quotient_image.json
python3 research/p192_extension_cover_20261012/prepare_line_new_all_fixed_origin_gp.py
gp -s 512000000 -q research/p192_extension_cover_20261012/verify_line_new_all_fixed_origin_native.gp
python3 research/p192_extension_cover_20261012/verify_line_new_all_receipt.py
python3 research/p192_extension_cover_20261012/line_full_pencil_graph.py
python3 research/p192_extension_cover_20261012/prepare_line_full_pencil_gp.py
gp -s 128000000 -q research/p192_extension_cover_20261012/verify_line_full_pencil.gp
```

The [run manifest](RUN.json) hashes these inputs, sources, receipts, logs,
proof modules, and replay outputs. Local elapsed durations in the receipts
are unisolated correctness diagnostics.

# The common P-192 split-line quotient constant and its centered form

**Result.** All 4,892 frozen split line cuts have the same nonidentity
first-quotient subgroup sum. The complete 9,784-point coordinate replay
matches the earlier independent projection digest. Translating each
quotient point by the fourth part of the common sum makes every four-point
line equation homogeneous, while preserving distinct point labels. The
previous private-column result still applies: the 4,892 equations have
rank 4,892 on their own line-row system and no nonzero combination
eliminates every new-point column to constrain only the original 141.

The field prime and subgroup order are

\[
\begin{aligned}
p&=6277101735386680763835789423207666416083908700390324961279,\\
r&=6277101735386680763835789423176059013767194773182842284081.
\end{aligned}
\]

For the fixed rational plane-quartic anchor \(A\), set
\(q(P)=h\phi([P-A])\), where \(\phi\) is the first elliptic quotient
map and \(h\) is the frozen cofactor projection. A line section with
old points \(P_i,P_j\) and residual points \(R_k,S_k\) is canonical,
so every split cut satisfies

\[
 q(P_i)+q(P_j)+q(R_k)+q(S_k)
   =h\phi([K-4A])=:C.
\]

The complete [checked-Sage replay](../../line_full_constant_replay.py)
projects all 9,784 new rational quartic points, checks their \(r\)-torsion
and signed distinctness, checks the first 64 in all three quotient
components, and verifies this equation for all 4,892 cuts. The ordered
point-key digest is the prior full-projection digest
`c139b81dfc72e3d7a71b37d4512c2a1b22082fd4e1721b84e912576e5764005a`.
An [independent prime-field replay](../../verify_line_full_constant.py)
reads the stored coordinates and recomputes every sum with a separate
elliptic-curve group law. It also checks \(C\) against the earlier 32-cut
sample, verifies \(rC=0\), and checks the center is the unique fourth
part of \(C\) in the \(r\)-subgroup.

| Exact quotient-group value | \(x\) | \(y\) |
| --- | ---: | ---: |
| Common line sum \(C\) | 5084569872798834823362581587648610749445923590216715218881 | 744970661357237284008179939626879199549785034920678237478 |
| Center \(D=[4^{-1}\bmod r]C\) | 2036922914201919277430838227202717277511694987407646586794 | 1745346030196626306372530920185700559379931648537267977553 |

For every line, \(\sum(q(P)-D)=C-4D=0\). The [Lean module](../../lean/LineCenter.lean)
proves this four-point identity for any commutative group and proves
translation by \(D\) is injective. The concrete \(C\), \(D\), quartic,
and quotient map are checked by the Sage and independent integer-field
receipts; the full finite-field curve model has not been formalized in
Lean. Neither \(C\) nor \(D\) is among the frozen 9,925 quotient points.

Centering removes the common-constant column from the homogeneous line
matrix, but it does not remove the private fresh-point columns. It also
changes a known-scalar target equation: a target \([a]G\) becomes
\([a]G-D\), whose scalar relative to \(G\) requires the scalar of \(D\)
or an explicit additional unknown. The danger criterion therefore tracks
offset handling, target anchoring, and core rank together.

## Cached projection candidate

A one-second profile of the established inverse found repeated large
finite-field powers in the point checks. The
[cached inverse](../../plane_point_inverse_cached.py) computes Frobenius
conjugates of model constants once per batch and retains the same inverse,
norm, round-trip, and quotient equations. A [frozen paired replay](../../verify_cached_plane_inverse.py)
matched exact cover coordinates and all three quotient images on 16 old
plus 64 new points, alternating implementation order. The local
unisolated elapsed totals were 8.362 seconds for the established path
and 0.332 seconds for the cached path. These are implementation
diagnostics on a heavily contended host; a Tier 0 or Tier 1 paired
benchmark with A/A noise controls is the next performance gate.

## Proof map and next experiment

| Result | Lean declaration | Concrete bridge |
| --- | --- | --- |
| Four centered points sum to zero when their original sum is \(4D\) | `P192ExtensionCover.four_point_line_centered_zero` | Full 4,892-cut Sage and independent prime-field replays. |
| Translation preserves distinct point labels | `P192ExtensionCover.translate_point_labels_injective` | Exact stored coordinates and signed-distinctness replay. |
| Every frozen split line has the same \(C\) | Open concrete elliptic-curve formalization | Canonical-section proof plus complete independent coordinate replay. |
| Private-column rank persists after centering | `P192ExtensionCover.private_column_coefficients_zero` and translation theorem | Full line-incidence and distinctness receipts. |

The checked Sage runtime, source and coordinate hashes, Lean build, axiom
audit, fresh kernel replay, and independent verifier are sealed in the
run manifest. The next attack-relevant experiment should first seek
relations whose fresh-point columns recur across independent cover
families, then measure known-scalar target yield and useful core rank.

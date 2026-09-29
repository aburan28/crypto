# Producer correctness argument

This argument covers the bounded standalone Boolean worker. It is not an external
audit or a machine-checked proof. The executable accepts the declared n12/16/20/24
fixtures. Tests also exercise small exhaustive domains and retained fallbacks.

## The affine tail is the complete intersection with the original row span

Write the original coefficient rows as M=[Q | L | c], with every quadratic
coordinate before the affine coordinates. Echelon reduction preserves the row
span. Rows whose leading coordinates are quadratic have distinct ordered pivots.
In any nonzero combination of those rows, the earliest selected quadratic pivot
cannot cancel. Consequently every affine row in span(M) belongs to the span of
the affine-tail rows. Conversely, every tail row belongs to span(M). Reducing the
tail to RREF therefore gives exactly span(M) intersected with the affine space.

No multiplication closure is performed. An affine member of the whole ideal
that requires multiplication may remain undiscovered. Every constraint actually
used here is still an exact linear combination of original equations.

## A projected-rank certificate can only skip an empty affine tail

Let T map each quadratic monomial coordinate to its declared one-hot coordinate
in F_2^64, extended linearly over GF(2). The implementation's monomial hash fixes T;
it is not a probabilistic assertion about the current rows.

If vQ=0, then vQT=0. Thus if QT has full row rank, its left kernel is zero and
the left kernel of Q is also zero. There is no nonzero combination of input rows
whose quadratic part vanishes, so the affine intersection is {0}. Skipping the
full reduction in this case is exact. This does not assume that T is injective.

A collision, or any other failure of the projected rows to have full rank, proves
nothing about the original rank. The worker then performs the complete reduction
above. It never substitutes projected equality for original coefficient equality.
Repeated input monomials cancel by XOR before rank is assessed. Tests include a
deliberate collision between distinct quadratic monomials, a dependent pair that
derives 1, and the exhaustive small row-space corpus.

## Restriction and recovery preserve precisely the original solutions

An affine RREF containing the constant row 1 certifies inconsistency. Otherwise,
its r pivot coordinates are affine expressions in the n-r free coordinates.
Keep the free coordinates in increasing original-label order and call the resulting
map x=A y+b. Its free-coordinate rows form an identity submatrix, so the map is
injective and its image is exactly the affine solution space.

Every original solution satisfies the affine consequences and therefore has a
unique preimage y. Substitution is performed in the Boolean quotient, including
y_i squared=y_i and GF(2) parity. Hence for every original equation f,
the transformed polynomial evaluates to f(A y+b). The recovered transformed
solutions are exactly the original solutions, even when terms cancel, degrees
drop, or new quadratic terms appear. No equation is discarded to fit a word.

The packed implementation applies the same transformation to 32 coefficient
planes simultaneously. For affine images a and b, the constant-linear terms
and diagonal products are added separately with XOR; off-diagonal contributions
toggle both symmetric locations. The direct list implementation independently
multiplies and canonicalizes monomial images. Exhaustive maps include noninjective
maps as a stronger test of the coefficient-transform primitive; solver recovery
uses only the injective RREF maps.

## Factored block transport is exact at each Gray step

Fix four low Boolean coordinates y and let h be the remaining coordinates. A
block stores all sixteen words F(y,h), one bit per original equation. For reflected
Gray order g(t)=t XOR (t>>1), step t flips j=ctz(t). Just before such a flip, the
lower high-coordinate bits are zero except bit j-1, which is one when j>0.

The change in every block entry is therefore

    D_j + W_j(y),
    D_j = linear[j] + quadratic[j,j-1] + sum(i>j) quadratic[j,i] h_i,
    W_j(y) = sum(i<4) quadratic[low_i,j] y_i.

The neighbor term is omitted for j=0. W_j is fixed for the entire enumeration.
Flipping high coordinate j changes only the scheduled D_i with i<j, each by
quadratic[i,j]. Keeping scalar D_i and fixed vector W_i therefore maintains all
sixteen equation values without reconstructing the block. An induction from the
directly evaluated h=0 block proves pointwise equality with direct evaluation.

The enumeration order and block accounting are unchanged from the retained
sixteen-point method. Scalar and native paths choose the same first zero lane,
charge the entire block, and preserve the wrapping checksum. The checksum is
only diagnostic; direct pointwise tests, original-model checks and the completed
retained search reference supply separate correctness evidence. NEON uses unsigned
zero detection; SSE2 uses exact equality masks. All loads and stores address full
four-u32 arrays. Architecture-specific correctness is checked on local AArch64
and, through CI, x86_64.

## Limits and costs

An enumeration cap reached before the next whole block returns UNKNOWN. Source
domains outside the specialized implementations take their retained bounded
fallbacks. A tree leaf scatters free coordinates back through its original labels;
initial restriction evaluates the RREF recovery map. Every returned model is
independently evaluated on the original equations by the worker.

The full cost includes certification, failed certificates, affine extraction,
map/transform construction, block setup, all updates and queries, recovery and
verification. A domain reduction from 2^n to 2^(n-r) is an algebraic fact; it is
not a solve-time ratio because setup, throughput and early-SAT ordering also change.
Tests and CI do not establish a production or cryptanalytic result.

# Exact CRT source-product equality and strong reference

Scope: representation of the equality between a monic quintic H over
GF(2^n) and five source roots whose polynomial-basis coordinates have degree
less than d. This component does not generate norm coefficients or solve
an original multivariate system. Common norm/input-solving costs are omitted
from both arms. No source root, speed gain or novelty is certified by counts.

Assume 4(d-1)<n. The first four elementary coefficients of the unreduced
product have degrees at most17j in the n83,d18 case, or j(d-1) generally,
below n. Only the constant coefficient can wrap: its unique lift is H0+f*q,
deg(q)<=5(d-1)-n, with a single shared q across all places. Higher coefficients
use canonical n-bit representatives; equality forces their extra bits to0.
Every source assignment satisfying the original field coefficient identity
therefore has exactly one admitted q, and every lifted identity reduces to
that original identity. Multiplicities and ordered source assignments remain.

Every difference coefficient has degree at most M=5(d-1), provided n<=M+1;
the selected worker panels have n<=M except the tiny n5,d2 control. That
control has M=n, also satisfies this bound. A panel with n>M+1 would instead
require degree max(M,n-1)+1 and is not admitted by this implementation.
Choose pairwise coprime moduli with total degree at least M+1. Vanishing in
each quotient ring is equivalent to divisibility by their product. As its
degree exceeds the difference degree, this iff the difference polynomial is0.
Rings may have nilpotents: T^4,(T+1)^4 and (T^2+T+1)^2 are valid and
pairwise coprime; irreducible places of degrees3,4,5,... supply the rest.
No independent q choices or marginal rejection-rate multiplication is used.

The comparator computes the same unreduced coefficient product directly,
using true per-coefficient degree bounds and sequential sparse polynomial
products. It does not pad every multiplication to n bits. Both use identical
Boolean AND/XOR hash-consing and constant/idempotence simplification. The
reported standard Tseitin counts use3clauses per AND,4per XOR, one unit per
nontrivial zero constraint. These are deterministic representation counts,
not counts from a SAT solve or evidence of end-to-end speed.

Independent controls enumerate1024 source tuples in n5,d2, form literal
carryless products, assign canonical H and the true quotient bit, then flip
every single H/q bit. Both circuits agree with all27648 literal cases. A
truncated CRT using T^4 alone accepts a changed T^4 coefficient that the
exact polynomial identity rejects; this is retained as an explicit false
positive for insufficient total degree. Generic proof covers the declared
six wide binary panels; no prime/odd arithmetic transfer is claimed.

Prediction: CRT reduces standard clause count by at least25percent on every
panel versus the strong low-degree reference. Any loss rejects blanket
replacement; retain all counts. Native90seconds/4GiB/onecore, source/input
hashes frozen before profiles. Late full serialization demotes all profile
flags; supervisor charges complete process lifetime. No wall-speed ratio.
Review is one-model plus independently formulated literal oracles, not Lean.

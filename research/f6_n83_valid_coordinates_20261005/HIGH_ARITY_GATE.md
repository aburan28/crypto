# Why a wider free-input suffix does not solve n83 eight-summand F6

This is an algebraic audit of the existing `wide_groebner::decompose_rec`
and `WideSystem::build_suffix` paths, not a performance measurement or a
lower bound on a different algorithm.

The n83 dimension-12 source base has 2,029 valid abscissa coordinates
representing 4,057 curve points. A direct eight-summand chained system
requires `8·12 + (8−2)·83 = 594` Boolean variables, while the current
monomial mask admits 128. The current fallback therefore constructs a
one-link suffix (`1·(12+83)=95` variables) and recursively branches on
its selected point, leaving an m=7 problem. It repeats until an m=3
head fits in 119 variables. Merely widening the mask to fit two suffix
links (`2·(12+83)=190` variables) does not change the following issue.

For a target group point `R` and any valid trailing points
`P₁,…,P_k`, set `W₀ = R − Σ_i P_i` and
`W_j = W₀ + Σ_{i≤j} P_i`. Then `W_k = R`, and each link
`S₃(x(W_{j−1}), x(P_j), x(W_j)) = 0` holds whenever those points have
affine x-coordinates. This constructs a root of the **free-input suffix**
for every ordinary valid trailing tuple. An intermediate at infinity
needs an exceptional-case representation, but does not create useful
selectivity. The suffix equations alone cannot reject the generic
trailing choices; the head constraint is where a relation can be
decided.

With nondecreasing coordinate order, the number of five-coordinate
tails before a three-summand head is

`C(2029+5−1,5) = C(2033,5) = 287,983,658,457,236`.

This is the exhaustion size of that enumerative recursion shape, not
a measured runtime, an expected time to a witness, or a complexity
lower bound for F6/IC. Signs, exceptional points, head restrictions,
and early witnesses alter actual work. It shows why the measured fast
three-summand block does not yet compose into a usable ordinary
eight-summand solver. A next method must couple a head constraint to
tail choices before enumerating them, or derive a smaller supported
set of intermediate group points. Giving a free incoming point more
Boolean mask words does not supply that constraint.

The simple group identity above also limits a generalized-birthday
shortcut: elliptic-curve point coordinates do not provide a known
public additive bit projection. A candidate that buckets point sums
must prove its bucket rule and charge the point additions used to
form those sums. No such method or speedup is established here.

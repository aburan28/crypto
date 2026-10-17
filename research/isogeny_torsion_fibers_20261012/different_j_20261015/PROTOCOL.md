# Frozen different-j degree-three attribution gate (2026-10-15)

Status: selector and checks frozen before running the search. This is distinct
from the scaling-isomorphic `x^3+1 -> x^3-27` control.

## Candidate and selector

For ascending primes `p` in `[61,271]` with `p=1 (mod 3)`, start from
`T=(3,1)` on `E:y^2=x^3-21x+37` over `F_p`, and its rational kernel
`{O,T,-T}`. Check nonsingularity, `[3]T=O`, and construct the degree-three
Vélu map by adding the two kernel translates and subtracting the kernel
coordinates. The destination coefficients to check independently are
`a'=-9a-30x_T^2`, `b'=-27b-70x_T^3-42a*x_T`. Require the source and
destination `j` values to differ in `F_(p^2)`, exact equal group orders,
and a prime factor `r>=19`, `r != 3`, with multiplicity one in the group
order. Reject a candidate if source and destination are isomorphic over
`F_(p^2)`; differing `j` is a certificate of this for nonsingular curves.
Keep every rejected candidate with its exact failed predicate. Select the
first passing `p` and the largest qualifying `r`; cap the screen at `271`.

## Paired ECDLP comparison

Freeze ICV1/EC1 identities, field modulus, `r`, cofactor, public generator and
target from seed `SHA256("different-j/2026-10-15/p=<p>/r=<r>")`, a line
chosen in order `s,1+s,2+s,...` with eight distinct projected columns, and
three query seeds `0x20261015 XOR p + round*0x9e3779b9`. The answer scalar is
excluded from solver input. Check the map on every point for the first cell
when feasible, plus sampled homomorphism and kernel checks, and independently
replay each recovered source scalar. Run source projected-pair, isogeny-image
projected-pair, destination direct relation solve, same-target negation rho,
and A/A controls on identical target/seed/trial cap `min(10^6,20r)`.

Charge map construction, source/destination base construction, map evaluation,
every failed query, relation validation, rank, target recovery, scalar replay,
memory, thread CPU and wall. Report cold `S=total counted source-addition
equivalents/sqrt(r)` only where conversion covers all native work, plus
separate phase counts and L0 timings. The generic two-class floor is
`sqrt(pi/4)=0.886226925`; the measured same-target rho is the reference.
Missing pair, map check, full solve, or isolated timing remains an open gate,
not evidence of a speed change.

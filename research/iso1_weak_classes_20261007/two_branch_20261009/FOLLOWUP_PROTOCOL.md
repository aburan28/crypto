# Adaptive small-field construction follow-up

Registered after the degree-2/3/5/7 validation. All 64 source fixtures remain
in the denominator. The original degree-2/3 run constructs 55 endpoints; the
broader validation leaves three additional exact-positive classes without a
route and one exact zero, and reaches a vertex cap on a source already solved
by the first run. This follow-up selects seeds 202610090033, 202610090042 and
202610090052 by that recorded status, retaining their original inputs.

Search rational cyclic isogenies of degrees 2,3,5,7,11,13,17,19,23 with at most
1024 distinct invariants and 90 seconds per source. Enumerate kernel factors
whose degrees divide `(ell-1)/2`, generate each cyclic subgroup over its
minimal point field, and project its kernel polynomial to the original field.
Check subgroup order, all coefficient projections, isogeny point images,
every retained route point count, and the final literal-family conversion.

This is an adaptive construction check, not a new population. Report its
cost alongside all prior attempts on the same source. Caps leave the exact
positive membership intact and the construction unresolved. The 192–252-bit
run retains its independently frozen degree-2/3 policy and caps.

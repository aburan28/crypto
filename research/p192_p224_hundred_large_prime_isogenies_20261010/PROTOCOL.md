# Hundred large prime-degree maps: frozen protocol

Requested outcome: 100 new distinct independently verified prime-degree maps
within the public P192 and P224 isogeny classes, with degree greater than 1009.
The number of distinct degrees is reported separately. Replays, duplicate maps,
and the previously verified P192/10453 and P224/1471 maps do not count.

Screen every prime in [1010, 1000000] on both published source curves. Retain
both Frobenius eigenvalues and their multiplicative orders. The initial kernel
constructor dispatch supports extension degrees 2 through 16, using the smallest
eigenvalue order (ties resolved by eigenvalue). Order one uses extension degree
two. Every extension modulus must pass the existing Rabin irreducibility checks.
The initial new walks use P192 degree 2341 and P224 degree 1117, both order nine.
Repeated degrees are allowed; a map counts only when its source model and
normalized kernel define a previously uncounted map. Aim for 50 on each class.

Construct one eigenline per step. Preserve coverage PARTIAL: one of two lines.
Retain the existing exact torsion, Frobenius, kernel descent and producer checks.
Independent replay uses separate limb field arithmetic and polynomial routines,
kernel validation, exact rational curve-equation substitution, subgroup order,
and all 20 existing public scalar-transport checks. No sampled identity replaces
the exact polynomial checks. A walk source must be the verified codomain of its
preceding map, anchored in the published source and generator. Record every
source and destination ICV1 and EC1 representation, and link raw map hashes.

Run native dedicated binaries locally. Constructor budget: 1800 seconds and
8 GiB RSS per invocation; retain timeouts and failures. Use at most two concurrent
walks, one per class. Record construction and independent replay elapsed times
separately and complete invocation cost. Failed or interrupted attempts remain
in the evidence; resume only from fully sealed verified steps. Stop on actual
resource exhaustion and retain the remaining requested count as unresolved.

Freeze code by commit before measured execution; retain binary digest, source
revision, dirty state, input parameters and host receipts. This campaign makes
map construction and correctness claims, with same-host independent arithmetic
replay. It does not establish a performance comparison or security conclusion.

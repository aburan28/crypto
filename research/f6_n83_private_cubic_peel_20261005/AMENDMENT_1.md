# Amendment 1: count every cubic monomial in the residual core

Registered after the sampled result was frozen in commit `538c647c4`
and before changing the probe or running this follow-on. Sampling 32
degree-three monomials per residual row certified only 719 of the 15,822
quartic-unresolved products, leaving 15,103. Because sampling can miss
a private cubic term, that count is a lower bound on certified rows.
This amendment measures the exact ceiling of the same single-pass
degree-three certificate by selecting **every** degree-three monomial
from each residual product.

Keep the original K0 curve, standard dimension-18 base, public T001
at torsion offset zero, 332 original equations, all 90 source
multipliers in the frozen order, and the same 32-candidate quartic
certificate. Require the frozen 15,822 quartic-unresolved rows and
BLAKE3 witness digest
`153342de74d4e40a2ea5f96ba35d24b7bdf2726700ed9b04714c9c89ca49dbdd`.
The exact cubic pass indexes every degree-three monomial in those rows
and counts it in **all** residual products plus all original equations.
A row is certified only if an exact monomial count is one. The old
sampled result remains untouched. Repeat the planted `[0,2,4,6,8]`
control and full group replay.

If at most 2,000 rows remain, construct the same original-plus-residual
core and run the deterministic reducer with the 6,500,000-column cap;
success requires a nonconstant source-only affine row or contradiction.
If more survive, stop after reporting the exact certified/unresolved
counts and ordered witness digest. Then the remaining count is an exact
bound for this *one-pass cubic-private* criterion, but says nothing about
other elimination strategies. A completed reduction without a positive
row rejects this two-level private-column screen for T001 offset zero.

Build one new native Rust release binary. Run planted then ordinary
once each, one thread, 300-second process timeout, and a live 7-GiB
observed RSS kill threshold. Preserve raw stdout/stderr, status, RSS
samples, build log, source/binary hashes and SHA-256 manifest. Time is
an exploratory structural diagnostic on this contended host, not an
F6 or end-to-end IC speedup.

# Protocol: exact degree-263 dual composition

Registered 2026-09-25 before the dual-construction run. This follows merged [#717](https://github.com/aburan28/crypto/pull/717) and the kernel-direction certificate [#743](https://github.com/aburan28/crypto/pull/743).

## Frozen input and reference

Use both preflight seeds 20260924 and 20260925 and all eight saved #717 kernel representatives: two descending and two horizontal lines per seed. Frozen SHA-256: `twist_torsion_results.json=6f1e7b22f3471214d38ec0d3196c88fe9edaf45791e9c05ea67764831fd75368`; `oriented_velu.py=a1a656cc32efd612250b1ecf50dff18af79a8f62f970ee79841bfd390603423b`; `redteam_velu_replay.py=c08fbb8de9d54906d783be967a36a6073b2a1d1a3583292badd3adcdd7126790`; `fastfield.py=b5990fda51700bbfba363251c41f02febfd0fc92edfea4f5ddd647bab472789e`; `relations.py=0180501ef8fe00f5c54f910822a28cd86dc3c2c1af164348976c1c2e25cc763f`. The independent reference is the #703 bit-polynomial full-point group law and direct extension-field Vélu sum, which imports no production arithmetic.

## Construction and exact checks

Let G be the saved twist generator of a 263-kernel line and let V be the second saved twist-basis point. Each saved representative line has first basis coefficient 1, so V is complementary to G. Construct the full oriented map phi:E→E', and the same-kernel map phi_tw on the quadratic twist E_tw→E'_tw. Set H'=phi_tw(V). Require exact order 263, require H' not infinity, and build psi:E'→E'' using H' on E'_tw as a generator of the reverse kernel. Require the normalized E'' model to equal E. Test psi(phi(P)) and psi(phi(Q)) against both +[263] and −[263] using full coordinates, record the sign, and require the sign to agree on every independently selected public scalar control. Check phi, phi_tw, psi and psi_tw on infinity; all 262 nonzero forward and reverse twist-kernel points map to infinity; off-curve points with a kernel abscissa reject before the exceptional path. Compare independently computed direct extension-field Vélu coordinates for forward and reverse images on P, Q, P+Q and 2P. Check reverse composition at full coordinates and challenge subgroup order.

## Cost and failure accounting

Charge field initialization, construction of all four maps from supplied generators, forward and reverse evaluations, independent reference calls, kernel enumeration and all failed attempts separately. The saved torsion basis is supplied; obtaining it from scratch is outside this stage and must be charged in a later cold experiment. Record operation counts (field multiplication, squaring and inversion, including nested inversion work), median wall time over 11 repeated production compositions after warmup, host, Python version, source/input hashes, and a deterministic digest excluding timing. Stop at 300 seconds and preserve a FAIL receipt with the exact failed stage. A failed construction is useful evidence: retain its counts and identify the next missing primitive rather than silently switching kernels.

## Claim boundary and decision

A verified dual identity proves correctness of transport and fixes its sign convention. It does not prove cheaper PDP, independent relation yield, or full ECDLP speedup. This package will contain a focused CI replay and a scoreboard update once the in-flight #746 scoreboard PR has merged. The next attack experiment remains matched native/transported bases and held-out targets with cold setup, rank and comparison to rho charged.
